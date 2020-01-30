
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "meryl.H"
#include "strings.H"
#include "system.H"

#include "sweatShop.H"

//  The number of KB to use for a merylCountArray segment.
#define SEGMENT_SIZE       64
#define SEGMENT_SIZE_BITS  (SEGMENT_SIZE * 1024 * 8)



#if 0
class mcTheadData {
public:
  mcThreadData() {
  };
  ~mcThreadData() {
  };
};
#endif




class mcGlobalData {
public:
  mcGlobalData(vector<merylInput *>   &inputs,
               merylOp                 op,
               uint64                  nPrefix,
               uint32                  wData,
               kmdata                  wDataMask,
               uint64                  maxMemory,
               merylFileWriter        *output) : _inputs(inputs) {
    _operation      = op;
    _nPrefix        = nPrefix;
    _wData          = wData;
    _wDataMask      = wDataMask;

    _data           = new merylCountArray [_nPrefix];
    _output         = output;
    _writer         = output->getBlockWriter();

    _maxMemory      = maxMemory;
    _memBase        = getProcessSize();
    _memUsed        = _memBase;
    _memReported    = 0;

    _kmersAdded     = 0;

    _inputPos       = 0;
    //_inputs         = inputs;

    _lastBuffer[0]  = 0;

    _computeWait    = 0;
    _numComputing   = 0;

    for (uint32 pp=0; pp<nPrefix; pp++)       //  Initialize each bucket.
      _memUsed += _data[pp].initialize(pp, wData, SEGMENT_SIZE);
  };

  ~mcGlobalData() {
    delete [] _data;
    delete [] _writer;
  };

  merylOp                _operation;        //  Parameters.
  uint64                 _nPrefix;
  uint32                 _wData;
  kmdata                 _wDataMask;

  merylCountArray       *_data;             //  Data for counting.
  merylFileWriter       *_output;
  merylBlockWriter      *_writer;           //  Data for writing.

  uint64                 _maxMemory;        //  Maximum memory we can use.
  uint64                 _memBase;          //  Overhead memory.
  uint64                 _memUsed;          //  Sum of actual memory used.
  uint64                 _memReported;      //  Memory usage at last report.

  uint64                 _kmersAdded;       //  Boring statistics for the user.

  uint32                 _inputPos;         //  Input files.
  vector<merylInput *>  &_inputs;

  char                   _lastBuffer[64];   //  Wrap-around from the last buffer.

  uint32                 _computeWait;
  uint32                 _numComputing;
};



class mcComputation {
public:
  mcComputation() {
    _bufferMax  = 1024 * 1024;
    _bufferLen  = 0;
    _buffer     = new char [_bufferMax];

    memset(_buffer, 0, sizeof(char) * _bufferMax);

    _memUsed    = 0;
    _kmersAdded = 0;
  };

  ~mcComputation() {
    delete [] _buffer;
  };

  //  Data for input sequences.
  uint64        _bufferMax;
  uint64        _bufferLen;
  char         *_buffer;

  //  Data for converting sequence to kmers.
  kmerIterator  _kiter;

  //  Data for debugging.
  //char          _fstr[65];      //  For debugging
  //char          _rstr[65];

  uint64        _memUsed;
  uint64        _kmersAdded;
};




void *
loadBases(void *G) {
  mcGlobalData     *g  = (mcGlobalData  *)G;
  mcComputation    *s  = new mcComputation();
  uint32            kl = kmerTiny::merSize() - 1;

  //  Where to put this?  Is it even useful?
  //fprintf(stderr, "Loading kmers from '%s' into buckets.\n", _inputs[ii]->_name);

  //  If no more inputs, we're done.
 loadBasesAgain:
  if (g->_inputPos >= g->_inputs.size())
    return(NULL);

  //  Copy the end of the last block into our buffer.

  if (g->_lastBuffer[0] != 0) {
    memcpy(s->_buffer, g->_lastBuffer, sizeof(char) * kl);
    s->_bufferLen += kl;
  }

  //  Try to load bases

  uint64  bMax     = s->_bufferMax - s->_bufferLen;
  uint64  bLen     = 0;
  bool    endOfSeq = false;

  g->_inputs[g->_inputPos]->loadBases(s->_buffer + s->_bufferLen,
                                      bMax,
                                      bLen, endOfSeq);

  //  If nothing loaded, move to the next file and try again.  If
  //  no more files, we're done.

  if (bLen == 0) {
    delete g->_inputs[g->_inputPos]->_sequence;
    g->_inputs[g->_inputPos]->_sequence = NULL;

    g->_inputPos++;

    goto loadBasesAgain;
  }

  //fprintf(stderr, "read " F_U64 " bases from '%s'\n", bufferLen, _inputs[ii]->_name);

  //  With bases loaded, we need to save the last few bases for the next buffer,
  //  and tell the kmerIterator about the bases we loaded.

  if (endOfSeq == false) {
    memcpy(g->_lastBuffer, s->_buffer + s->_bufferLen - kl, sizeof(char) * kl);
  }

  else {
    g->_lastBuffer[0] = 0;
  }

  //  Now just tell the iterator about the buffer.

  s->_kiter.addSequence(s->_buffer, s->_bufferLen);

  //  That's it.

  return(s);
}



void
insertKmers(void *G, void *T, void *S) {
  mcGlobalData     *g = (mcGlobalData  *)G;
  //mcThreadData     *t = (mcThreadData  *)T;
  mcComputation    *s = (mcComputation *)S;

  uint64  memUsed    = 0;
  uint64  kmersAdded = 0;

  while (s->_kiter.nextMer()) {
    bool    useF = (g->_operation == opCountForward);
    kmdata  pp   = 0;
    kmdata  mm   = 0;

    if (g->_operation == opCount)
      useF = (s->_kiter.fmer() < s->_kiter.rmer());

    if (useF == true) {
      pp = (kmdata)s->_kiter.fmer() >> g->_wData;
      mm = (kmdata)s->_kiter.fmer()  & g->_wDataMask;
      //fprintf(stderr, "useF F=%s R=%s ms=%u pp %lu mm %lu\n", s->_kiter.fmer().toString(fstr), s->_kiter.rmer().toString(rstr), s->_kiter.fmer().merSize(), pp, mm);
    }

    else {
      pp = (kmdata)s->_kiter.rmer() >> g->_wData;
      mm = (kmdata)s->_kiter.rmer()  & g->_wDataMask;
      //fprintf(stderr, "useR F=%s R=%s ms=%u pp %lu mm %lu\n", s->_kiter.fmer().toString(fstr), s->_kiter.rmer().toString(rstr), s->_kiter.rmer().merSize(), pp, mm);
    }

    assert(pp < g->_nPrefix);

    s->_memUsed    += g->_data[pp].add(mm);
    s->_kmersAdded += 1;
  }
}




void
writeBatch(void *G, void *S) {
  mcGlobalData     *g = (mcGlobalData  *)G;
  mcComputation    *s = (mcComputation *)S;

  //  Udpate memory used and kmers added.  There's only one writer thread,
  //  so this is thread safe!

  g->_memUsed    += s->_memUsed;
  g->_kmersAdded += s->_kmersAdded;

  //  Do some accounting.

  if (g->_memUsed - g->_memReported > (uint64)128 * 1024 * 1024) {
    g->_memReported = g->_memUsed;

    fprintf(stderr, "Used %3.3f GB out of %3.3f GB to store %12lu kmers.\n",
            g->_memUsed   / 1024.0 / 1024.0 / 1024.0,
            g->_maxMemory / 1024.0 / 1024.0 / 1024.0,
            g->_kmersAdded);
  }

  //  Free the input buffer.

  delete s;

  //  If we haven't hit the memory limit yet, just return.

  if (g->_memUsed < g->_maxMemory)
    return;

  //  Otherwise, we're out of space and need to dump the data.  But first,
  //  we need to wait for all compute threads to finish.

  //  Do we need someway of knowing how many compute threads are actually computing?
  //  Or can we just wait until all the expectd number of threads finishes?


  //  If not the first one to wait, signal that we're waiting too,
  //  then wait for the first one to write the data out.

  if (g->_computeWait > 0) {
    g->_computeWait++;

    while (g->_computeWait > 0)
      usleep(100);

    return;
  }

  //  We're the first one to wait.  Wait for all threads to finish.

  while (g->_computeWait < g->_numComputing)
    usleep(100);

  //  Write data!

  fprintf(stderr, "Memory full.  Writing results to '%s', using " F_S32 " threads.\n",
          g->_output->filename(), omp_get_max_threads());
  fprintf(stderr, "\n");

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<g->_output->numberOfFiles(); ff++) {
    //fprintf(stderr, "thread %2u writes file %2u with prefixes 0x%016lx to 0x%016lx\n",
    //        omp_get_thread_num(), ff, g->_output->firstPrefixInFile(ff), g->_output->lastPrefixInFile(ff));

    for (uint64 pp=g->_output->firstPrefixInFile(ff); pp <= g->_output->lastPrefixInFile(ff); pp++) {
      g->_data[pp].countKmers();                   //  Convert the list of kmers into a list of (kmer, count).
      g->_data[pp].dumpCountedKmers(g->_writer);   //  Write that list to disk.
      g->_data[pp].removeCountedKmers();           //  And remove the in-core data.
    }
  }

  g->_writer->finishBatch();

  //  Reset accounting.

  g->_memUsed    = g->_memBase;

  for (uint32 pp=0; pp<g->_nPrefix; pp++)
    g->_memUsed += g->_data[pp].usedSize();

  g->_kmersAdded = 0;

  //  Signal that threads can proceeed.

  g->_computeWait = 0;
}



void
merylOperation::countThreads(uint32  wPrefix,
                             uint64  nPrefix,
                             uint32  wData,
                             kmdata  wDataMask) {

  //  If we're only configuring, stop now.

  if (_onlyConfig)
    return;

  //  Configure the writer for the prefix bits we're counting with.
  //
  //  We split the kmer into wPrefix and wData (bits) pieces.
  //  The prefix is used by the filewriter to decide which file to write to, by simply
  //  shifting it to the right to keep the correct number of bits in the file.
  //
  //           kmer -- [ wPrefix (18) = prefixSize               | wData (36) ]
  //           file -- [ numFileBits  | prefixSize - numFileBits ]

  _output->initialize(wPrefix);

  //  Initialize the counter.

  mcGlobalData  *g = new mcGlobalData(_inputs, _operation, nPrefix, wData, wDataMask, _maxMemory, _output);

  //  Set up a sweatShop and run it.

  sweatShop    *ss = new sweatShop(loadBases, insertKmers, writeBatch);
  uint32        nt = omp_get_max_threads();

  ss->setLoaderQueueSize(4 * nt);
  ss->setWriterQueueSize(1 * nt);

  ss->setNumberOfWorkers(1 * nt);

  //for (uint32 w=0; w<nt; w++)
  //  ss->setThreadData(w, td[w] = new mcThreadData(g, w));

  ss->run(g, true);

  delete ss;

  //for (uint32 w=0; w<nt; w++)
  //  delete td[w];
  //delete [] td;

  //  All data loaded.  Write the output.

  fprintf(stderr, "\n");
  fprintf(stderr, "Writing results to '%s', using " F_S32 " threads.\n",
          _output->filename(), omp_get_max_threads());

  //for (uint64 pp=0; pp<nPrefix; pp++)
  //  fprintf(stderr, "Prefix 0x%016lx writes to file %u\n", pp, _output->fileNumber(pp));

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<_output->numberOfFiles(); ff++) {
    //fprintf(stderr, "thread %2u writes file %2u with prefixes 0x%016lx to 0x%016lx\n",
    //        omp_get_thread_num(), ff, _output->firstPrefixInFile(ff), _output->lastPrefixInFile(ff));

    for (uint64 pp=_output->firstPrefixInFile(ff); pp <= _output->lastPrefixInFile(ff); pp++) {
      g->_data[pp].countKmers();                   //  Convert the list of kmers into a list of (kmer, count).
      g->_data[pp].dumpCountedKmers(g->_writer);   //  Write that list to disk.
      g->_data[pp].removeCountedKmers();           //  And remove the in-core data.
    }
  }

  //  Merge any iterations into a single file, or just rename
  //  the single file to the final name.

  g->_writer->finish();

  delete g->_writer;   //  Explicitly delete the writer.
  g->_writer = NULL;

  //  Cleanup.
  delete g;

  fprintf(stderr, "\n");
  fprintf(stderr, "Finished counting.\n");
}

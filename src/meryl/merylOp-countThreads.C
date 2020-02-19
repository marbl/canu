
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

#include <atomic>

//  The number of KB to use for a merylCountArray segment.
#define SEGMENT_SIZE       64
#define SEGMENT_SIZE_BITS  (SEGMENT_SIZE * 1024 * 8)



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

    _dumping        = false;

    _lock           = new std::atomic_flag [_nPrefix];
    _data           = new merylCountArray  [_nPrefix];
    _output         = output;
    _writer         = output->getBlockWriter();

    _maxMemory      = maxMemory;
    _memBase        = getProcessSize();
    _memUsed        = _memBase;
    _memReported    = 0;

    _kmersAdded     = 0;

    _inputPos       = 0;
    //_inputs         = inputs;

    for (uint32 ii=0; ii<65; ii++)
      _lastBuffer[ii] = 0;

    _computeWait    = 0;
    _numComputing   = 0;

    for (uint32 pp=0; pp<_nPrefix; pp++) {      //  Initialize each bucket.
      _lock[pp].clear();
      _memUsed += _data[pp].initialize(pp, wData, SEGMENT_SIZE);
    }
  };

  ~mcGlobalData() {
    delete [] _lock;
    delete [] _data;
    delete [] _writer;
  };

  merylOp                _operation;        //  Parameters.
  uint64                 _nPrefix;
  uint32                 _wData;
  kmdata                 _wDataMask;

  bool                   _dumping;

  std::atomic_flag      *_lock;
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

  char                   _lastBuffer[65];   //  Wrap-around from the last buffer.

  uint32                 _computeWait;
  uint32                 _numComputing;
};



#define BUF_MAX  (1 * 1024 * 1024)

class mcComputation {
public:
  mcComputation() {
    _bufferMax  = BUF_MAX;
    _bufferLen  = 0;

    _memUsed    = 0;
    _kmersAdded = 0;
  };

  ~mcComputation() {
  };

  //  Data for input sequences.
  uint64        _bufferMax;
  uint64        _bufferLen;
  char          _buffer[BUF_MAX];

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

  //  Copy the end of the last block into our buffer.

  assert(s->_bufferLen == 0);

  if (g->_lastBuffer[0] != 0) {
    memcpy(s->_buffer, g->_lastBuffer, sizeof(char) * kl);

    s->_bufferLen += kl;

    g->_lastBuffer[0] = 0;
  }

  //  If no more inputs, we're done.

  if (g->_inputPos >= g->_inputs.size())
    return(NULL);

  //  Try to load bases.  Keep loading until the buffer is filled
  //  or we exhaust the file.

  while (1) {
    uint64  bMax     = s->_bufferMax - s->_bufferLen;
    uint64  bLen     = 0;
    bool    endOfSeq = false;

    if (bMax < 512)   //  If the buffer is full enough,
      break;          //  stop loading.

    //  Load bases, but reserve 2 characters in the buffer for a
    //  sequence terminating ','.

    bool success = g->_inputs[g->_inputPos]->loadBases(s->_buffer + s->_bufferLen,
                                                       bMax - 2,
                                                       bLen, endOfSeq);

    //  If no bases loaded, we've exhausted the file.  Close it,
    //  and move to the next one.  Bail on loading more bases;
    //  let the test above figure out if we're all done on
    //  the next call to loadBases().

    if (success == false) {
      assert(bLen == 0);

      s->_buffer[s->_bufferLen++] = '.';   //  Insert a mer-breaker, just to be safe.

      delete g->_inputs[g->_inputPos]->_sequence;
      g->_inputs[g->_inputPos]->_sequence = NULL;

      g->_inputPos++;

      break;
    }

    //  Account for whatever we just loaded.

    s->_bufferLen += bLen;

    assert(s->_bufferLen+1 <= s->_bufferMax);

    //  If we're at the end of a sequence, append a mer-breaker.

    if (endOfSeq == true)
      s->_buffer[s->_bufferLen++] = '.';
  }

  //  With bases loaded, we need to save the last few bases for the next buffer,
  //  and tell the kmerIterator about the bases we loaded.

  if (s->_buffer[s->_bufferLen-1] != '.')
    memcpy(g->_lastBuffer, s->_buffer + s->_bufferLen - kl, sizeof(char) * kl);

  //  Now just tell the iterator about the buffer.

  s->_kiter.addSequence(s->_buffer, s->_bufferLen);

  //  That's it.

  return(s);
}



void
insertKmers(void *G, void *T, void *S) {
  mcGlobalData     *g = (mcGlobalData  *)G;
  mcComputation    *s = (mcComputation *)S;

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

    //  If we're dumping data, stop immediately and sleep until dumping is
    //  finished.

    while (g->_dumping == true)
      usleep(1000);

    //  We need exclusive access to this specific merylCountArray, so busy
    //  wait on a lock until we get it.

    while (g->_lock[pp].test_and_set(std::memory_order_relaxed) == true)
      ;

    s->_memUsed    += g->_data[pp].add(mm);
    s->_kmersAdded += 1;

    g->_lock[pp].clear(std::memory_order_relaxed);
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
  //  Otherwise, dump data.

  if (g->_memUsed < g->_maxMemory)
    return;

  //  Tell all the threads to pause, then grab all the locks to ensure nobody
  //  is still writing.

  g->_dumping = true;

  for (uint32 pp=0; pp<g->_nPrefix; pp++)
    while (g->_lock[pp].test_and_set(std::memory_order_relaxed) == true)
      ;

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

  for (uint32 pp=0; pp<g->_nPrefix; pp++)
    g->_lock[pp].clear(std::memory_order_relaxed);

  g->_dumping = false;
}



void
merylOperation::countThreads(uint32  wPrefix,
                             uint64  nPrefix,
                             uint32  wData,
                             kmdata  wDataMask) {

  //  If we're only configuring, stop now.

  if (_onlyConfig)
    return;

  fprintf(stderr, "THREADED VERSION\n");

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

  ss->setLoaderBatchSize(1 * nt);
  ss->setLoaderQueueSize(2 * nt);
  ss->setWriterQueueSize(1 * nt);
  ss->setNumberOfWorkers(1 * nt);

  ss->run(g, false);

  delete ss;

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

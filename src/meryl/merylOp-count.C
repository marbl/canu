
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

//  The number of KB to use for a merylCountArray segment.
#define SEGMENT_SIZE       64
#define SEGMENT_SIZE_BITS  (SEGMENT_SIZE * 1024 * 8)


//
//  mcaSize       = sizeof(merylCountArray)  == 80
//  ptrSize       = sizeof(uint64 *)         == 8
//  segSize       = (as above)               ~= 64 KB
//
//  mersPerPrefix = nKmers / nPrefix+1
//  mersPerSeg    = 8 * segSize / 2*merSize-prefixSize
//
//  nSegPerPrefix = mersPerPrefix / mersPerSeg
//
//             basic structure             pointers to data                         data
//           |-----------------|   |------------------------------|   |------------------------------|
//           |                 |   |                              |   |                              |
//  memory = (mcaSize * nPrefix) + (ptrSize * nPrefix * nSegPrefix) + (segSize * nPrefix * nSegPrefix)
//
//  Solving for nKmers
//
//  memory = (mcaSize * nPrefix) + nSegPerPrefix * (ptrSize * nPrefix + segSize * nPrefix)
//
//  memory - (mcaSize * nPrefix) = nSegPerPrefix * (ptrSize * nPrefix + segSize * nPrefix)
//
//  nSegPerPrefix = (memory - mcaSize * nPrefix) / (ptrSize * nPrefix + nPrefix * segSize)
//
//  (nKmers / nPrefix+1) / mersPerSeg = (memory - mcaSize * nPrefix) / (ptrSize * nPrefix + segSize * nPrefix)
//   nKmers / nPrefix+1  = mersPerSeg * (memory - mcaSize * nPrefix) / (ptrSize * nPrefix + segSize * nPrefix)

uint64
findMaxInputSizeForMemorySize(uint32 merSize, uint64 memSize) {
  uint64  mcaSize = sizeof(merylCountArray<uint32>);
  uint64  ptrSize = sizeof(uint64 *);
  uint64  segSize = SEGMENT_SIZE * 1024;

  //  Free variable - prefixSize (wp)

  fprintf(stderr, "For memory size %lu bytes\n", memSize);
  fprintf(stderr, "                                                  |---------memory-breakdown--------\n");
  fprintf(stderr, "pBits    nKmers     nPrefix   nSegment   mers/seg  structure   pointers         data\n");
  fprintf(stderr, "----- ---------- ---------- ---------- ---------- ---------- ---------- ------------\n");

  for (uint32 wp=1; wp < 2*merSize; wp++) {
    uint64  nPrefix       = (uint64)1 << wp;

    if (mcaSize * nPrefix > memSize)
      break;

    //  Compute how many mer suffixes we can fit in a segment of memory.

    uint64  mersPerSeg = (segSize * 8) / (2 * merSize - wp);

    //  Compute the number of segments we can fit in memeory.  Each prefix must have a segment.  We
    //  don't actually know how many pointers we'll need -- it depends on the number of prefixes and
    //  number of segments -- so we assume there are 64 per prefix.
    //
    uint64  nSeg = (memSize - mcaSize * nPrefix - ptrSize * 64 * nPrefix) / segSize;

    if (nSeg < nPrefix)
      nSeg = nPrefix;

    //  Assuming some segment loading average, compute the number of kmers we can fit.  Each prefix
    //  has a 2/3 full segment, then any left over segments are 100% full.

    uint64  nKmers = nPrefix * mersPerSeg * 0.666;

    if (nPrefix < nSeg)
      nKmers += (nSeg - nPrefix) * mersPerSeg;

    //  For flavoring, show the memory breakdown.  nSegPrefix underflows, so needs to be included directly.

    uint64  ptrPerPrefix = 64; //32 * (nSeg / nPrefix + 1);  //  Somewhat incorrect.  The basic allocation is 32 pointers, then it doubles.

    uint64  basicMem     = mcaSize * nPrefix;                                    //  Size of the merCountArray structure itself.
    uint64  ptrMem       = ptrSize * nPrefix * ptrPerPrefix;                     //  Additional allocation for pointers to segments.
    uint64  dataMem      = segSize * nSeg;                                       //  Additional allocation of segments.

    if (basicMem + ptrMem + dataMem > 4 * memSize)
      break;

    fprintf(stderr, "%5u %10lu %10lu %10lu %10lu %10lu %10lu %12lu%s\n",
            wp, nKmers,
            nPrefix, nSeg,
            mersPerSeg,
            basicMem, ptrMem, dataMem,
            (basicMem + ptrMem + dataMem < memSize) ? "" : " INVALID");
  }

  exit(0);
}




//  Memory used by the simple algorithm depends only on the kmer size and the count
//  of the most frequent kmer (or highest count allowed).
//
void
findExpectedSimpleSize(uint64  nKmerEstimate,
                       uint64 &memoryUsed_) {
  uint32   lowBitsSize     = sizeof(lowBits_t) * 8;
  uint64   nEntries        = (uint64)1 << (2 * kmerTiny::merSize());

  uint64   expMaxCount     = 0.004 * nKmerEstimate;
  uint64   expMaxCountBits = countNumberOfBits64(expMaxCount) + 1;
  uint64   extraBits       = (expMaxCountBits < lowBitsSize) ? (0) : (expMaxCountBits - lowBitsSize);

  uint64   lowMem          = nEntries * lowBitsSize;
  uint64   highMem         = nEntries * extraBits;
  uint64   totMem          = (lowMem + highMem) / 8;

  memoryUsed_ = UINT64_MAX;

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "SIMPLE MODE\n");
  fprintf(stderr, "-----------\n");
  fprintf(stderr, "\n");

  if (kmerTiny::merSize() > 20) {
    fprintf(stderr, "  Disabled for mers larger than 20.\n");
    return;
  }

  fprintf(stderr, "  %2u-mers\n", kmerTiny::merSize());
  fprintf(stderr, "    -> %lu entries for counts up to %u.\n", nEntries, ((uint32)1 << lowBitsSize) - 1);
  fprintf(stderr, "    -> %lu %cbits memory used\n", scaledNumber(lowMem),  scaledUnit(lowMem));
  fprintf(stderr, "\n");
  fprintf(stderr, "  %lu input bases\n", nKmerEstimate);
  fprintf(stderr, "    -> expected max count of %lu, needing %lu extra bits.\n", expMaxCount, extraBits);
  if (extraBits > 0)
    fprintf(stderr, "    -> %lu %cbits memory used\n", scaledNumber(highMem), scaledUnit(highMem));
  else
    fprintf(stderr, "    -> no memory used\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  %lu %cB memory needed\n", scaledNumber(totMem),  scaledUnit(totMem));

  memoryUsed_ = totMem;
}



void
findBestPrefixSize(uint64  nKmerEstimate,
                   uint32 &bestPrefix_,
                   uint64 &memoryUsed_) {
  uint32  merSize    = kmerTiny::merSize();

  bestPrefix_  = 0;
  memoryUsed_  = UINT64_MAX;

  //  IMPORTANT!  Prefixes must be at least six bits - the number of bits
  //  we use for deciding on a file - and probably need to then leave one
  //  bit for a block id.  So only save a bestPrefix if it is at least 7.

  //  IMPORTANT!  Smaller prefixes mean bigger blocks in the output file,
  //  and any operation on those outputs needs to load all kmers in
  //  the block into memory, as full kmers.  Only save a bestPrefix
  //  if it is at least 10 - giving us at least 16 blocks per file.

  //  IMPORTANT!  Start searching at 1 and stop at merSize-1, otherwise
  //  we end up with a prefix or a suffix of size zero.

  for (uint32 wp=1; wp < 2 * merSize - 1; wp++) {
    uint64  nPrefix          = (uint64)1 << wp;                          //  Number of prefix == number of blocks of data
    uint64  kmersPerPrefix   = nKmerEstimate / nPrefix + 1;              //  Expected number of kmers we need to store per prefix
    uint64  kmersPerSeg      = SEGMENT_SIZE_BITS / (2 * merSize - wp);   //  Kmers per segment
    uint64  segsPerPrefix    = kmersPerPrefix / kmersPerSeg + 1;         //

    if (wp + countNumberOfBits64(segsPerPrefix) + countNumberOfBits64(SEGMENT_SIZE) + 10 >= 64)
      break;   //  Otherwise, dataMemory overflows.

    uint64  structMemory     = ((sizeof(merylCountArray<uint32>) * nPrefix) +          //  Basic structs
                                (sizeof(uint64 *)        * nPrefix * segsPerPrefix));  //  Pointers to segments
    uint64  dataMemory       = nPrefix * segsPerPrefix * SEGMENT_SIZE * 1024;
    uint64  totalMemory      = structMemory + dataMemory;

    //  Pick a larger prefix if it is dramatically smaller than what we have.
    //  More prefixes seem to run a bit slower, but also have smaller buckets
    //  for sorting at the end.

    if ((wp > 9) && (totalMemory + 16 * wp * 1024 * 1024 < memoryUsed_)) {
      memoryUsed_ = totalMemory;
      bestPrefix_ = wp;
    }

    if (totalMemory > 16 * memoryUsed_)
      break;
  }
}



void
findBestValues(uint64  nKmerEstimate,
               uint32  bestPrefix,
               uint64  memoryUsed,
               uint32 &wPrefix_,
               uint64 &nPrefix_,
               uint32 &wData_,
               uint64 &wDataMask_) {
  uint32  merSize = kmerTiny::merSize();

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "COMPLEX MODE\n");
  fprintf(stderr, "------------\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "prefix     # of   struct   kmers/    segs/     data    total\n");
  fprintf(stderr, "  bits   prefix   memory   prefix   prefix   memory   memory\n");
  fprintf(stderr, "------  -------  -------  -------  -------  -------  -------\n");

  for (uint32 wp=1; wp < 2 * merSize - 1; wp++) {
    uint64  nPrefix          = (uint64)1 << wp;                          //  Number of prefix == number of blocks of data
    uint64  kmersPerPrefix   = nKmerEstimate / nPrefix + 1;              //  Expected number of kmers we need to store per prefix
    uint64  kmersPerSeg      = SEGMENT_SIZE_BITS / (2 * merSize - wp);   //  Kmers per segment
    uint64  segsPerPrefix    = kmersPerPrefix / kmersPerSeg + 1;         //

    if (wp + countNumberOfBits64(segsPerPrefix) + countNumberOfBits64(SEGMENT_SIZE) + 10 >= 64)
      break;   //  Otherwise, dataMemory overflows.

    uint64  structMemory     = ((sizeof(merylCountArray<uint32>) * nPrefix) +          //  Basic structs
                                (sizeof(uint64 *)        * nPrefix * segsPerPrefix));  //  Pointers to segments
    uint64  dataMemory       = nPrefix * segsPerPrefix * SEGMENT_SIZE * 1024;
    uint64  totalMemory      = structMemory + dataMemory;

    fprintf(stderr, "%6" F_U32P "  %4" F_U64P " %cP  %4" F_U64P " %cB  %4" F_U64P " %cM  %4" F_U64P " %cS  %4" F_U64P " %cB  %4" F_U64P " %cB",
            wp,
            scaledNumber(nPrefix),        scaledUnit(nPrefix),
            scaledNumber(structMemory),   scaledUnit(structMemory),
            scaledNumber(kmersPerPrefix), scaledUnit(kmersPerPrefix),
            scaledNumber(segsPerPrefix),  scaledUnit(segsPerPrefix),
            scaledNumber(dataMemory),     scaledUnit(dataMemory),
            scaledNumber(totalMemory),    scaledUnit(totalMemory));

    if (wp == bestPrefix) {
      fprintf(stderr, "  Best Value!\n");

      wPrefix_   = wp;
      nPrefix_   = nPrefix;
      wData_     = 2 * merSize - wp;
      wDataMask_ = uint64MASK(wData_);

    } else {
      fprintf(stderr, "\n");
    }

    if (totalMemory > 16 * memoryUsed)
      break;
  }
}



void
reportNumberOfOutputs(uint64   nKmerEstimate,
                      uint64   memoryUsed,        //  expected memory needed for counting in one block
                      uint64   memoryAllowed,     //  memory the user said we can use
                      bool     useSimple) {
  uint32  nOutputsI      = memoryUsed / memoryAllowed + 1;
  double  nOutputsD      = (double)memoryUsed / memoryAllowed - (nOutputsI - 1);


  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "FINAL CONFIGURATION\n");
  fprintf(stderr, "-------------------\n");

  if (useSimple == true) {
    assert(nOutputsI == 1);
  }

  else {
    char    batchString[64] = { 0 };

    if      (nOutputsD < 0.2) {
      nOutputsI += 0;
      snprintf(batchString, 42, "split into up to %u (possibly %u)", nOutputsI-1, nOutputsI);
    }

    else if (nOutputsD < 0.8) {
      nOutputsI += 0;
      snprintf(batchString, 42, "split into up to %u", nOutputsI);
    }

    else {
      nOutputsI += 1;
      snprintf(batchString, 42, "split into up to %u (possibly %u)", nOutputsI, nOutputsI+1);
    }


    if (nOutputsI > 1) {
      fprintf(stderr, "\n");
      fprintf(stderr, "WARNING:\n");
      fprintf(stderr, "WARNING: Cannot fit into " F_U64 " %cB memory limit.\n", scaledNumber(memoryAllowed), scaledUnit(memoryAllowed));
      fprintf(stderr, "WARNING: Will %s batches, and merge them at the end.\n", batchString);
      fprintf(stderr, "WARNING:\n");
    }

    if (nOutputsI > 32) {
      fprintf(stderr, "WARNING: Large number of batches.  Increase memory for better performance.\n");
      fprintf(stderr, "WARNING:\n");
    }
  }

  //  This is parsed by Canu.  Do not change.

  fprintf(stderr, "\n");
  fprintf(stderr, "Configured %s mode for %.3f GB memory per batch, and up to %u batch%s.\n",
          (useSimple == true) ? "simple" : "complex",
          ((memoryUsed < memoryAllowed) ? memoryUsed : memoryAllowed) / 1024.0 / 1024.0 / 1024.0,
          nOutputsI,
          (nOutputsI == 1) ? "" : "es");
  fprintf(stderr, "\n");
}



void
merylOperation::configureCounting(uint64   memoryAllowed,      //  Input:  Maximum allowed memory in bytes
                                  bool    &useSimple_,         //  Output: algorithm to use
                                  uint32  &wPrefix_,           //  Output: Number of bits in the prefix (== bucket address)
                                  uint64  &nPrefix_,           //  Output: Number of prefixes there are (== number of buckets)
                                  uint32  &wData_,             //  Output: Number of bits in kmer data
                                  uint64  &wDataMask_) {       //  Output: A mask to return just the data of the mer

  //
  //  Check kmer size, presence of output, and guess how many bases are in the inputs.
  //

  uint32  merSize = kmerTiny::merSize();

  if (kmerTiny::merSize() == 0)
    fprintf(stderr, "ERROR: Kmer size not supplied with modifier k=<kmer-size>.\n"), exit(1);

  if ((_output == NULL) && (_onlyConfig == false))
    fprintf(stderr, "ERROR: No output specified for count operation.\n"), exit(1);

  if (_expNumKmers == 0)
    _expNumKmers = guesstimateNumberOfkmersInInput();

  //
  //  Report what we're trying.
  //

  fprintf(stderr, "\n");
  fprintf(stderr, "Counting %lu (estimated) %s %s%s%s " F_U32 "-mers from " F_SIZE_T " input file%s:\n",
          scaledNumber(_expNumKmers), scaledName(_expNumKmers),
          (_operation == opCount)        ? "canonical" : "",
          (_operation == opCountForward) ? "forward" : "",
          (_operation == opCountReverse) ? "reverse" : "",
          kmerTiny::merSize(), _inputs.size(), (_inputs.size() == 1) ? "" : "s");

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    fprintf(stderr, "  %15s: %s\n", _inputs[ii]->inputType(), _inputs[ii]->_name);

  //
  //  Set up to use the simple algorithm.
  //

  uint64  memoryUsedSimple = UINT64_MAX;

  findExpectedSimpleSize(_expNumKmers, memoryUsedSimple);

  //
  //  Set up to use the complex algorithm.
  //

  uint64   memoryUsedComplex = UINT64_MAX;
  uint32   bestPrefix        = 0;

  findBestPrefixSize(_expNumKmers, bestPrefix, memoryUsedComplex);
  findBestValues(_expNumKmers, bestPrefix, memoryUsedComplex, wPrefix_, nPrefix_, wData_, wDataMask_);

  //
  //  Decide simple or complex.
  //

  bool    useSimple  = false;
  uint64  memoryUsed = 0;


  if ((memoryUsedSimple < memoryUsedComplex) &&
      (memoryUsedSimple < memoryAllowed)) {
    useSimple  = true;
    memoryUsed = memoryUsedSimple;
  }

  else {
    useSimple  = false;
    memoryUsed = memoryUsedComplex;
  }

  reportNumberOfOutputs(_expNumKmers, memoryUsed, memoryAllowed, useSimple);
}



//  Return a complete guess at the number of kmers in the input files.  No
//  rigorous went into the multipliers, just looked at a few sets of lambda reads.
uint64
merylOperation::guesstimateNumberOfkmersInInput_dnaSeqFile(dnaSeqFile *sequence) {
  uint64  numMers = 0;
  char   *name    = sequence->filename();
  uint32  len     = strlen(name);

  if ((name[0] == '-') && (len == 1))
    return(0);

  uint64  size = AS_UTL_sizeOfFile(name);

  if      ((len > 3) && (name[len-3] == '.') && (name[len-2] == 'g') && (name[len-1] == 'z'))
    numMers += size * 3;

  else if ((len > 4) && (name[len-4] == '.') && (name[len-3] == 'b') && (name[len-2] == 'z') && (name[len-1] == '2'))
    numMers += size * 3.5;

  else if ((len > 3) && (name[len-3] == '.') && (name[len-2] == 'x') && (name[len-1] == 'z'))
    numMers += size * 4.0;

  else
    numMers += size;

  return(numMers);
}


#ifdef CANU
uint64
merylOperation::guesstimateNumberOfkmersInInput_sqStore(sqStore *store, uint32 bgnID, uint32 endID) {
  uint64  numMers = 0;

  for (uint32 ii=bgnID; ii<endID; ii++)
    numMers += store->sqStore_getRead(ii)->sqRead_sequenceLength();

  return(numMers);
}
#endif


uint64
merylOperation::guesstimateNumberOfkmersInInput(void) {
  uint64  guess = 0;

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    if (_inputs[ii]->isFromSequence())
      guess += guesstimateNumberOfkmersInInput_dnaSeqFile(_inputs[ii]->_sequence);

#ifdef CANU
    if (_inputs[ii]->isFromStore())
      guess += guesstimateNumberOfkmersInInput_sqStore(_inputs[ii]->_store, _inputs[ii]->_sqBgn, _inputs[ii]->_sqEnd);
#endif
  }

  return(guess);
}



void
merylOperation::count(uint32  wPrefix,
                      uint64  nPrefix,
                      uint32  wData,
                      uint64  wDataMask) {

  //configureCounting(_maxMemory, useSimple, wPrefix, nPrefix, wData, wDataMask);

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

  kmerCountBlockWriter  *_writer = _output->getBlockWriter();

  //  Allocate buckets.  The buckets don't allocate space for mers until they're added,
  //  and allocate space for these mers in blocks of 8192 * 64 bits.
  //
  //  Need someway of balancing the number of prefixes we have and the size of each
  //  initial allocation.

  merylCountArray<uint32>  *data = new merylCountArray<uint32> [nPrefix];

  //  Load bases, count!

  uint64          bufferMax  = 1024 * 1024;
  uint64          bufferLen  = 0;
  char           *buffer     = new char [bufferMax];
  bool            endOfSeq   = false;

  memset(buffer, 0, sizeof(char) * bufferMax);

  //char            fstr[65];
  //char            rstr[65];

  uint64          memBase     = getProcessSize();   //  Overhead memory.
  uint64          memUsed     = 0;                  //  Sum of actual memory used.
  uint64          memReported = 0;                  //  Memory usage at last report.

  memUsed = memBase;

  for (uint32 pp=0; pp<nPrefix; pp++)
    memUsed += data[pp].initialize(pp, wData, SEGMENT_SIZE);

  uint64          kmersAdded  = 0;

  kmerIterator    kiter;

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    fprintf(stderr, "Loading kmers from '%s' into buckets.\n", _inputs[ii]->_name);

    while (_inputs[ii]->loadBases(buffer, bufferMax, bufferLen, endOfSeq)) {
      if (bufferLen == 0)
        continue;

      //fprintf(stderr, "read " F_U64 " bases from '%s'\n", bufferLen, _inputs[ii]->_name);

      kiter.addSequence(buffer, bufferLen);

      while (kiter.nextMer()) {
        bool    useF = (_operation == opCountForward);
        uint64  pp   = 0;
        uint64  mm   = 0;

        if (_operation == opCount)
          useF = (kiter.fmer() < kiter.rmer());

        if (useF == true) {
          pp = (uint64)kiter.fmer() >> wData;
          mm = (uint64)kiter.fmer()  & wDataMask;
          //fprintf(stderr, "F %s %s %u pp %lu mm %lu\n", kiter.fmer().toString(fstr), kiter.rmer().toString(rstr), kiter.fmer().merSize(), pp, mm);
        }

        else {
          pp = (uint64)kiter.rmer() >> wData;
          mm = (uint64)kiter.rmer()  & wDataMask;
          //fprintf(stderr, "R %s %s %u pp %lu mm %lu\n", kiter.fmer().toString(fstr), kiter.rmer().toString(rstr), kiter.rmer().merSize(), pp, mm);
        }

        assert(pp < nPrefix);

        memUsed += data[pp].add(mm);

        kmersAdded++;
      }

      if (endOfSeq)      //  If the end of the sequence, clear
        kiter.reset();   //  the running kmer.

      //  Report that we're actually doing something.

      if (memUsed - memReported > (uint64)128 * 1024 * 1024) {
        memReported = memUsed;

        fprintf(stderr, "Used %3.3f GB out of %3.3f GB to store %12lu kmers.\n",
                memUsed    / 1024.0 / 1024.0 / 1024.0,
                _maxMemory / 1024.0 / 1024.0 / 1024.0,
                kmersAdded);
      }

      //  If we're out of space, process the data and dump.

      if (memUsed > _maxMemory) {
        fprintf(stderr, "Memory full.  Writing results to '%s', using " F_S32 " threads.\n",
                _output->filename(), omp_get_max_threads());
        fprintf(stderr, "\n");

#pragma omp parallel for schedule(dynamic, 1)
        for (uint32 ff=0; ff<_output->numberOfFiles(); ff++) {
          //fprintf(stderr, "thread %2u writes file %2u with prefixes 0x%016lx to 0x%016lx\n",
          //        omp_get_thread_num(), ff, _output->firstPrefixInFile(ff), _output->lastPrefixInFile(ff));

          for (uint64 pp=_output->firstPrefixInFile(ff); pp <= _output->lastPrefixInFile(ff); pp++) {
            data[pp].countKmers();                //  Convert the list of kmers into a list of (kmer, count).
            data[pp].dumpCountedKmers(_writer);   //  Write that list to disk.
            data[pp].removeCountedKmers();        //  And remove the in-core data.
          }
        }

        _writer->finishBatch();

        kmersAdded = 0;

        memUsed = memBase;                        //  Reinitialize or memory used.
        for (uint32 pp=0; pp<nPrefix; pp++)
          memUsed += data[pp].usedSize();
      }

    }

    //  Would like some kind of report here on the kmers loaded from this file.

    delete _inputs[ii]->_sequence;
    _inputs[ii]->_sequence = NULL;
  }

  //  Finished loading kmers.  Free up some space.

  //delete [] kmers;
  delete [] buffer;

  //  Sort, dump and erase each block.
  //
  //  A minor complication is that within each output file, the blocks must be in order.

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
      data[pp].countKmers();                //  Convert the list of kmers into a list of (kmer, count).
      data[pp].dumpCountedKmers(_writer);   //  Write that list to disk.
      data[pp].removeCountedKmers();        //  And remove the in-core data.
    }
  }

  //  Merge any iterations into a single file, or just rename
  //  the single file to the final name.

  _writer->finish();

  delete _writer;
  _writer = NULL;

  //  Cleanup.

  delete [] data;

  fprintf(stderr, "\n");
  fprintf(stderr, "Finished counting.\n");
}

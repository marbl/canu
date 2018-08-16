
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
 */

#include "meryl.H"
#include "strings.H"



void
merylOperation::countSimple(void) {
  uint64          bufferMax  = 1300000;
  uint64          bufferLen  = 0;
  char           *buffer     = new char     [bufferMax];
  bool            endOfSeq   = false;

  uint64          kmersLen   = 0;
  kmerTiny       *kmers      = new kmerTiny [bufferMax];

  kmerTiny        fmer;
  kmerTiny        rmer;

  uint32          kmerLoad   = 0;
  uint32          kmerValid  = fmer.merSize() - 1;
  uint32          kmerSize   = fmer.merSize();

  uint64          maxKmer    = (uint64)1 << (2 * kmerSize);

  char            str[32];
  char            strf[32];
  char            strr[32];

  if (kmerSize == 0)
    fprintf(stderr, "ERROR: Kmer size not supplied with modifier k=<kmer-size>.\n"), exit(1);

  if ((_output == NULL) && (_onlyConfig == false))
    fprintf(stderr, "ERROR: No output specified for count operation.\n"), exit(1);

  omp_set_num_threads(_maxThreads);

  fprintf(stderr, "\n");
  fprintf(stderr, "Counting %s%s%s " F_U32 "-mers from " F_SIZE_T " input file%s:\n",
          (_operation == opCount)        ? "canonical" : "",
          (_operation == opCountForward) ? "forward" : "",
          (_operation == opCountReverse) ? "reverse" : "",
          kmerSize, _inputs.size(), (_inputs.size() == 1) ? "" : "s");

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    fprintf(stderr, "  %s\n", _inputs[ii]->_name);

  //  Allocate memory.

  //  Tiny counter uses a fixed-width array (either 8- or 16-bit counts) and a
  //  variable-width array to store counts.
  //
  //  Estimating how much memory is needed is the same as estimating the highest
  //  count in the dataset.
  //  For 40x corrected human,  with 115 Gbp, the largest count is  25 million, 25 bits, 0.02% of the input.
  //  For 68x raw pacbio human, with 192 Gbp, the largest count is 152 million, 28 bits, 0.08% of the input.

  typedef uint16 lowBits_t;

  uint32        lowBitsSize = sizeof(lowBits_t) * 8;
  uint32        lowBitsMax  = ((uint32)1 << lowBitsSize) - 1;
  uint32        highBitMax  = 0;

  uint64        nKmersGuess = guesstimateNumberOfkmersInInput();
  uint64        expMaxCount = 0.004 * nKmersGuess;

  uint64        expMemory   = (maxKmer * sizeof(lowBits_t) +                                   //  Fixed data,   16-bits wide
                               maxKmer * (logBaseTwo64(expMaxCount) + 1 - lowBitsSize) / 8);   //  Variable data, 1-bit  wide

  fprintf(stderr, "\n");
  fprintf(stderr, "lowBitsSize %u\n", lowBitsSize);
  fprintf(stderr, "lowBitsMax  %u\n", lowBitsMax);
  fprintf(stderr, "highBitMax  %u\n", highBitMax);
  fprintf(stderr, "nKmersGuess %lu\n", nKmersGuess);
  fprintf(stderr, "expMaxCount %lu\n", expMaxCount);
  fprintf(stderr, "expMemory   %lu\n", expMemory);
  fprintf(stderr, "maxKmer     %lu\n", maxKmer);

  fprintf(stderr, "\n");
  fprintf(stderr, "Expecting to use " F_U64 " %cB memory to count " F_U64 "%s " F_U32 "-mers.\n",
          scaledNumber(expMemory),         scaledUnit(expMemory),
          scaledNumber(nKmersGuess, 1000), scaledName(nKmersGuess, 1000),
          kmerSize);
  fprintf(stderr, "\n");

  //  If we're only configuring, stop now.

  if (_onlyConfig)
    return;

  lowBits_t    *lowBits     = new lowBits_t    [maxKmer];
  bitArray     *highBits    = new bitArray [64];

  memset(lowBits,  0, sizeof(lowBits_t) * maxKmer);

  //  Load bases, count!

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    fprintf(stderr, "Loading kmers from '%s' into buckets.\n", _inputs[ii]->_name);

    while (_inputs[ii]->loadBases(buffer, bufferMax, bufferLen, endOfSeq)) {
      if (bufferLen == 0)
        continue;

      //fprintf(stderr, "read %lu bases from '%s'\n", bufferLen, _inputs[ii]->_name);

      kmersLen = 0;

      for (uint64 bb=0; bb<bufferLen; bb++) {
        if ((buffer[bb] != 'A') && (buffer[bb] != 'a') &&   //  If not valid DNA, don't
            (buffer[bb] != 'C') && (buffer[bb] != 'c') &&   //  make a kmer, and reset
            (buffer[bb] != 'G') && (buffer[bb] != 'g') &&   //  the count until the next
            (buffer[bb] != 'T') && (buffer[bb] != 't')) {   //  valid kmer is available.
          kmerLoad = 0;
          continue;
        }

        fmer.addR(buffer[bb]);
        rmer.addL(buffer[bb]);

        if (kmerLoad < kmerValid) {   //  If not a full kmer, increase the length we've
          kmerLoad++;                 //  got loaded, and keep going.
          continue;
        }

        if      (_operation == opCount)
          kmers[kmersLen++] = (fmer < rmer) ? fmer : rmer;

        else if (_operation == opCountForward)
          kmers[kmersLen++] = fmer;

        else
          kmers[kmersLen++] = rmer;
      }

      if (endOfSeq)                   //  If the end of the sequence, clear
        kmerLoad = 0;                 //  the running kmer.

      //  Now, just pass our list of kmers to the counting engine.

      for (uint64 kk=0; kk<kmersLen; kk++) {
        uint64  kidx = (uint64)kmers[kk];
        uint32  hib  = 0;

        assert(kidx < maxKmer);

        //  If we can add one to the low bits, do it and get outta here.

        if (lowBits[kidx] < lowBitsMax) {
          lowBits[kidx]++;
          continue;
        }

        //  Otherwise, we need to do some manual addition.

        lowBits[kidx] = 0;

        for (uint32 hib=0; hib < 64; hib++) {
          if (highBits[hib].isAllocated() == false) {
            fprintf(stderr, "Increasing to %u-bit storage (for kmer 0x%016lx).\n",
                    lowBitsSize + hib + 1, kidx);
            highBits[hib].allocate(maxKmer);
          }

          if (highBits[hib].flipBit(kidx) == 0) {    //  If not set, set it,
            highBitMax = max(highBitMax, hib + 1);   //  remember the possible maximum bit set,
            break;                                   //  and stop.
          }
        }
      }
    }

    //  Would like some kind of report here on the kmers loaded from this file.

    delete _inputs[ii]->_sequence;
    _inputs[ii]->_sequence = NULL;
  }

  //  Finished loading kmers.  Free up some space before dumping.

  delete [] kmers;
  delete [] buffer;

  //  Then dump.

  uint32                 wPrefix    = 10;
  uint32                 wSuffix    = fmer.merSize() * 2 - wPrefix;

  uint64                 nPrefix    = ((uint64)1 << wPrefix);
  uint64                 nSuffix    = ((uint64)1 << wSuffix);

  uint64                 sMask      = ((uint64)1 << wSuffix) - 1;

  _output->initialize(wPrefix);

  fprintf(stderr, "\n");
  fprintf(stderr, "Writing results to '%s', using " F_S32 " threads.\n",
          _output->filename(), omp_get_max_threads());
  fprintf(stderr, "  wPrefix  %u\n",  wPrefix);
  fprintf(stderr, "  wSuffix  %u\n",  wSuffix);
  fprintf(stderr, "  nPrefix  %lu\n", nPrefix);
  fprintf(stderr, "  nSuffix  %lu\n", nSuffix);
  fprintf(stderr, "  sMask    0x%016lx\n", sMask);

#if 0
  for (uint32 ff=0; ff<_output->numberOfFiles(); ff++) {
    uint64  kStart   = (_output->firstPrefixInFile(ff) << wSuffix);
    uint64  kEnd     = (_output->lastPrefixInFile(ff)  << wSuffix) | sMask;

    fprintf(stderr, "file %2u with prefixes 0x%016lx to 0x%016lx for kmers 0x%016lx to 0x%016lx\n",
            ff,
            _output->firstPrefixInFile(ff),
            _output->lastPrefixInFile(ff),
            kStart,
            kEnd);
  }
#endif

  fprintf(stderr, "\n");


#pragma omp parallel for
  for (uint32 ff=0; ff<_output->numberOfFiles(); ff++) {
    uint64  bStart = (_output->firstPrefixInFile(ff));
    uint64  bEnd   = (_output->lastPrefixInFile(ff));

    //fprintf(stderr, "thread %2u writes file %2u with prefixes 0x%016lx to 0x%016lx for kmers 0x%016lx to 0x%016lx\n",
    //        omp_get_thread_num(), ff, bStart, bEnd,
    //        (_output->firstPrefixInFile(ff) << wSuffix),
    //        (_output->lastPrefixInFile(ff)  << wSuffix) | sMask);

    uint64  *sBlock  = new uint64 [nSuffix];
    uint32  *cBlock  = new uint32 [nSuffix];

    //  Iterate over kmers that belong in this data file.  For each kmer, reconstruct the count
    //  from our bit-sliced array, adding it to the list of kmers to output if it exists.

    for (uint64 bp=bStart; bp<=bEnd; bp++) {
      uint64  kStart = (bp << wSuffix);
      uint64  kEnd   = (bp << wSuffix) | sMask;
      uint64  nKmers = 0;

      //fprintf(stderr, "thread %2u writes file %2u - prefix 0x%016lx for kmers 0x%016lx to 0x%016lx\n",
      //        omp_get_thread_num(), ff, bp, kStart, kEnd);

      for (uint64 kk=kStart; kk<=kEnd; kk++) {
        uint32  count = 0;

        for (uint32 aa=highBitMax; aa-- > 0; ) {    //  Reconstruct the count.
          count <<= 1;
          count  |= highBits[aa].getBit(kk);
        }

        count <<= lowBitsSize;
        count  |= lowBits[kk];

        if (count > 0) {
          sBlock[nKmers] = kk & sMask;
          cBlock[nKmers] = count;
          nKmers++;
        }
      }

      //  With the kmers reconstructed, write this block of data to the file.
      _output->addBlock(bp, nKmers, sBlock, cBlock);
    }

    delete [] sBlock;
    delete [] cBlock;
  }

  //  Even though there are no iterations, we still need to finish.

  _output->finishIteration();

  //  Cleanup.

  delete [] lowBits;
  delete [] highBits;
}

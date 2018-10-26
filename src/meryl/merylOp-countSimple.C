
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



void
merylOperation::countSimple(void) {
  uint64          bufferMax  = 1300000;
  uint64          bufferLen  = 0;
  char           *buffer     = new char     [bufferMax];
  bool            endOfSeq   = false;

  uint64          kmersLen   = 0;
  kmerTiny       *kmers      = new kmerTiny [bufferMax];

  uint64          nEntries   = (uint64)1 << (2 * kmerTiny::merSize());

  //  If we're only configuring, stop now.

  if (_onlyConfig)
    return;

  uint32          lowBitsSize     = sizeof(lowBits_t) * 8;
  uint32          lowBitsMax      = ((uint32)1 << lowBitsSize) - 1;   //  Largest value that can be stored in lowBits.
  lowBits_t      *lowBits         = new lowBits_t    [nEntries];

  uint32          highBitMax      = 0;
  bitArray       *highBits        = new bitArray [64];

  memset(lowBits,  0, sizeof(lowBits_t) * nEntries);

  //  Load bases, count!

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    fprintf(stderr, "Loading kmers from '%s' into buckets.\n", _inputs[ii]->_name);

    while (_inputs[ii]->loadBases(buffer, bufferMax, bufferLen, endOfSeq)) {
      if (bufferLen == 0)
        continue;

      //fprintf(stderr, "read %lu bases from '%s'\n", bufferLen, _inputs[ii]->_name);

      kmersLen = 0;

      kmerIterator  kiter(buffer, bufferLen);

      while (kiter.nextMer()) {
        if      (_operation == opCount)
          kmers[kmersLen++] = (kiter.fmer() < kiter.rmer()) ? kiter.fmer() : kiter.rmer();

        else if (_operation == opCountForward)
          kmers[kmersLen++] = kiter.fmer();

        else
          kmers[kmersLen++] = kiter.rmer();
      }

      if (endOfSeq)                   //  If the end of the sequence, clear
        kiter.reset();                //  the running kmer.

      //  Now, just pass our list of kmers to the counting engine.

      for (uint64 kk=0; kk<kmersLen; kk++) {
        uint64  kidx = (uint64)kmers[kk];
        uint32  hib  = 0;

        assert(kidx < nEntries);

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
            highBits[hib].allocate(nEntries);
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
  uint32                 wSuffix    = kmerTiny::merSize() * 2 - wPrefix;

  uint64                 nPrefix    = ((uint64)1 << wPrefix);
  uint64                 nSuffix    = ((uint64)1 << wSuffix);

  uint64                 sMask      = ((uint64)1 << wSuffix) - 1;

  _output->initialize(wPrefix);

  kmerCountBlockWriter  *_writer = _output->getBlockWriter();

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
      _writer->addBlock(bp, nKmers, sBlock, cBlock);
    }

    delete [] sBlock;
    delete [] cBlock;
  }

  //  Even though there are no iterations, we still need to finish.

  _writer->finish();

  delete _writer;
  _writer = NULL;

  //  Cleanup.

  delete [] lowBits;
  delete [] highBits;
}

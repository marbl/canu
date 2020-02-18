
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

  uint64          nEntries   = (uint64)1 << (2 * kmerTiny::merSize() - 2 * _countSuffixLength);

  //  If we're only configuring, stop now.

  if (_onlyConfig)
    return;

  uint32          lowBitsSize     = sizeof(lowBits_t) * 8;
  uint32          lowBitsMax      = ((uint32)1 << lowBitsSize) - 1;   //  Largest value that can be stored in lowBits.
  lowBits_t      *lowBits         = new lowBits_t    [nEntries];

  uint32          highBitMax      = 0;
  bitArray       *highBits        = new bitArray [64];

  memset(lowBits,  0, sizeof(lowBits_t) * nEntries);

  //  Generate a mask for the count-suffix.

  kmdata  suffixMask = 0;
  kmdata  suffixTest = _countSuffix;

  for (uint32 ii=0; ii<_countSuffixLength; ii++) {
    suffixMask <<= 2;
    suffixMask  |= 0x3;
  }

  //  Load bases, count!

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    fprintf(stderr, "Loading kmers from '%s' into buckets.\n", _inputs[ii]->_name);

    while (_inputs[ii]->loadBases(buffer, bufferMax, bufferLen, endOfSeq)) {
      if (bufferLen == 0)
        continue;

      kmerIterator  kiter(buffer, bufferLen);

      while (kiter.nextMer()) {
        kmdata  kidx;
        uint32  hib  = 0;

        //  Figure out if we want the forward or reverse kmer.

        if      (_operation == opCount)
          kidx = (kiter.fmer() < kiter.rmer()) ? kiter.fmer() : kiter.rmer();

        else if (_operation == opCountForward)
          kidx = kiter.fmer();

        else
          kidx = kiter.rmer();

        //  If there is a suffix, test the kidx against it and strip it off.

        if (_countSuffixLength > 0) {
          if ((kidx & suffixMask) != suffixTest)
            continue;

          kidx >>= (2 * _countSuffixLength);
        }

        //  Check that we fit into the array.

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
            highBits[hib].allocate(nEntries);
          }

          if (highBits[hib].flipBit(kidx) == 0) {    //  If not set, set it,
            highBitMax = max(highBitMax, hib + 1);   //  remember the possible maximum bit set,
            break;                                   //  and stop.
          }
        }
      }

      if (endOfSeq)                   //  If the end of the sequence, clear
        kiter.reset();                //  the running kmer.
    }

    //  Would like some kind of report here on the kmers loaded from this file.

    delete _inputs[ii]->_sequence;
    _inputs[ii]->_sequence = NULL;
  }

  //  Finished loading kmers.  Free up some space before dumping.

  delete [] buffer;

  //  Then dump.
  //
  //  A kmer is [ file ][ blockPrefix ][ suffix ][ count-suffix ].
  //
  //  The kmer counting array does not include [count-suffix].
  //  It is indexed by [ file ][ blockPrefix ][ suffix ].
  //
  //  We parallelize over [file], which is hard coded to be 6 bits.
  //
  //  Within each file, we retrieve all the data for a single [blockPrefix]
  //  at a time, convert this into an array of actual kmers (by adding any
  //  count-suffix, but removing the [file][blockPrefix]) and values, then
  //  writing this block to disk.
  //
  //  So we can give the file writer something to compress, and to make our
  //  memory usage constant, we want each block to have a reasonably large
  //  number of kmers.
  //
  //    wSuffix  wSuffix    nSuffix      memory   (kmdata = 128 bits = 16 bytes)
  //    (bases)   (bits)    (kmers)        (MB)   (kmvalu =  32 bits =  4 bytes)
  //         10       20    1048576          20   (per thread)
  //          4        8        256           -
  //          2        4         16           -
  //
  //  Our constraints are:
  //     wPrefix >=  6 -- one block per file
  //     wSuffix >=  5 -- 32 suffixes per block
  //     wSuffix <= 20 -- 1048576 suffixes per block
  //
  //  psbits is the number of bits we can use for [blockPrefix][suffix].
  //  We'll use at most 20 bits of that for [suffix] and the rest for
  //  [blockPrefix].  Plus 6 because wPrefix is both [file][blockPrefix].
  //

  uint32                 psbits     = kmerTiny::merSize() * 2 - _countSuffixLength * 2 - 6;

  uint32                 wSuffix    = (psbits > 20) ? 20 : psbits;
  uint32                 wPrefix    = 6 + psbits - wSuffix;

  uint64                 nPrefix    = ((uint64)1 << wPrefix);
  uint64                 nSuffix    = ((uint64)1 << wSuffix);

  uint64                 sMask      = ((uint64)1 << wSuffix) - 1;
  uint64                 cSuffix    = ((kmdata)_countSuffix);

  _output->initialize(wPrefix);

  merylBlockWriter  *_writer = _output->getBlockWriter();

  fprintf(stderr, "\n");
  fprintf(stderr, "Writing results to '%s', using " F_S32 " threads.\n",
          _output->filename(), omp_get_max_threads());
  fprintf(stderr, "             [ file ][  prefix ][  suffix ][ count-suffix ]\n");
  fprintf(stderr, "   widths    [    6 ][ %7u ][ %7u ][ %12u ]\n", wPrefix - 6, wSuffix, 2 * _countSuffixLength);
  fprintf(stderr, "   number    [   64 ][ %7lu ][ %7lu ][ %12s ]\n", nPrefix / 64, nSuffix, _countSuffixString);
  fprintf(stderr, "\n");


#pragma omp parallel for
  for (uint32 ff=0; ff<_output->numberOfFiles(); ff++) {
    uint64  bStart = (_output->firstPrefixInFile(ff));
    uint64  bEnd   = (_output->lastPrefixInFile(ff));

    kmdata  *sBlock  = new kmdata [nSuffix];   //  Suffixes  -- kk and sMask should properly be kmdata too.
    kmvalu  *cBlock  = new kmvalu [nSuffix];   //  Counts       but we're guaranteed to have small mers here.

    //  Iterate over kmers that belong in this data file.  For each kmer,
    //  reconstruct the count from our bit-sliced array, adding it to the
    //  list of kmers to output if it exists.

    for (uint64 bp=bStart; bp<=bEnd; bp++) {      //  Get kmer data from kmer 'kStart' to 'kEnd'.
      uint64  kStart = (bp << wSuffix);
      uint64  kEnd   = (bp << wSuffix) | sMask;
      uint64  nKmers = 0;

#if 1
      fprintf(stderr, "thread %2d wokring on block 0x%08lx<0x%08lx<0x%08lx %8lu kmers between 0x%016lx and 0x%016lx\n",
              omp_get_thread_num(),
              bStart, bp, bEnd,
              nSuffix,
              kStart, kEnd);
#endif

      for (uint64 kk=kStart; kk<=kEnd; kk++) {
        uint32  count = 0;

        for (uint32 aa=highBitMax; aa-- > 0; ) {         //  Reconstruct the count.
          count <<= 1;
          count  |= highBits[aa].getBit(kk);
        }

        count <<= lowBitsSize;
        count  |= lowBits[kk];

        if (count > 0) {                                 //  If the count is non-zero, add the kmer
          sBlock[nKmers]   = kk & sMask;                 //  to the block to output.
          sBlock[nKmers] <<= (2 * _countSuffixLength);
          sBlock[nKmers]  |= cSuffix;

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

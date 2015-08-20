
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
 *    Brian P. Walenz from 2005-FEB-07 to 2014-APR-11
 *      are Copyright 2005-2007,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <unistd.h>
#include <time.h>
#include <math.h>

#include "util++.H"

uint32 wordSize = 41;
uint32 testSize =  1 * 1024 * 1024;
uint32 arrySize =  1 * 1024 * 1024;

int
uint64compare(const void *a, const void *b) {
  const uint64   A = *(const uint64 *)a;
  const uint64   B = *(const uint64 *)b;
  if (A<B) return(-1);
  if (A>B) return(1);
  return(0);
}

int
main(int argc, char **argv) {

  mt_s *mtctx = mtInit(time(NULL));

  //  Test the bitPackedArray by writing a bunch of random gibberish
  //  to it, and see if it's the same.

  uint32  *pos = new uint32 [testSize];
  uint64  *val = new uint64 [testSize];
  uint64  *ans = new uint64 [arrySize];

  bitPackedArray *ARR  = new bitPackedArray(wordSize, 16);
  uint32          fail = uint32ZERO;

#if 1
  fprintf(stderr, "Touching the end of the array and clearing.\n");
  //ARR->set(arrySize, 0);
  //ARR->clear();

  fprintf(stderr, "Generating random test data.\n");

  //  Hit every element first, just to do it
  for (uint32 i=0; i<arrySize; i++) {
    pos[i]       = i;
    val[i]       = mtRandom64(mtctx);
    val[i]      &= uint64MASK(wordSize);
    ans[pos[i]]  = val[i];
  }

  //  Then hit random elements, with replacement, looking for bugs
  for (uint32 i=arrySize; i<testSize; i++) {
    pos[i]       = mtRandom32(mtctx) % arrySize;
    val[i]       = mtRandom64(mtctx);
    val[i]      &= uint64MASK(wordSize);
    ans[pos[i]]  = val[i];
  }

  fprintf(stderr, "Filling array.\n");

  for (uint32 i=0; i<testSize; i++)
    ARR->set(pos[i], val[i]);

  fprintf(stderr, "Validating array.\n");

  for (uint32 i=0; i<arrySize; i++)
    if (ARR->get(i) != ans[i]) {
      fprintf(stderr, "FAIL at i="uint32FMT"\n", i);
      fail++;

      if (fail > 1024) {
        fprintf(stderr, "bitPackedArray has errors, aborting!\n");
        return(1);
      }
    }

  if (fail) {
    fprintf(stderr, "bitPackedArray had "uint32FMT" errors.\n", fail);
    return(1);
  }

  fprintf(stderr, "OK!\n");
#endif

  delete    ARR;
  delete [] pos;
  delete [] val;
  delete [] ans;

  //
  //
  //

  for (uint32 testNum=0; testNum<32; testNum++) {
    uint32  thisTestSize = 0;
    uint32  thisWordSize = 0;

    //  Test a BIG heap the first iteration.
    if (testNum == 0) {
      thisTestSize = 857353; //23987153;
      thisWordSize = 63;

      fprintf(stderr, "Building heap "uint32FMT" (wordsize="uint32FMT" testsize="uint32FMT").\n",
              testNum, thisWordSize, thisTestSize);
    } else {
      thisTestSize = (mtRandom64(mtctx) % (2 * testNum)) * 1024 + 1024;
      thisWordSize = (mtRandom64(mtctx) % 63) + 1;
    }

    uint32  blockSize = mtRandom64(mtctx) % 32 + 1;
    bitPackedHeap  *HEAP = new bitPackedHeap(thisWordSize, blockSize);

    val = new uint64 [thisTestSize];
    for (uint32 i=0; i<thisTestSize; i++) {
      val[i]       = mtRandom64(mtctx);
      val[i]      &= uint64MASK(thisWordSize);
      HEAP->add(val[i]);
    }

    fprintf(stderr, "Testing heap "uint32FMT" (wordsize="uint32FMT" testsize="uint32FMT").\n",
            testNum, thisWordSize, thisTestSize);

    qsort(val, thisTestSize, sizeof(uint64), uint64compare);

    for (uint32 i=0; i<thisTestSize; i++) {
      uint64  h = HEAP->get();

      //fprintf(stderr, "val["uint32FMT"]="uint64FMT" -- HEAP="uint64FMT"\n", i, val[i], h);

      if (val[i] != h) {
        fprintf(stderr, "val["uint32FMT"]="uint64FMT" !! HEAP="uint64FMT"\n", i, val[i], h);
        fail++;
        if (fail > 25) {
          fprintf(stderr, "bitPackedHeap has errors, aborting!\n");
          return(1);
        }
      }
    }

    if (fail) {
      fprintf(stderr, "bitPackedHeap had "uint32FMT" errors.!\n", fail);
      return(1);
    }

    delete    HEAP;
    delete [] val;
  }

  fprintf(stderr, "OK!\n");

  return(fail);
}


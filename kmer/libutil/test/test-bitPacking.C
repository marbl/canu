
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
 *    Brian P. Walenz on 2004-MAY-06
 *      are Copyright 2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-MAY-24 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2007-DEC-05 to 2014-APR-11
 *      are Copyright 2007,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "util++.H"

//  An integer multiplier on the test length.  1 is pretty quick, but
//  10 is the default.
//
#define TEST_LENGTH (10 * 1024 * 1024)

//  We test
//
//  1)  binary encoding/decoding
//
//  2)  pre/post increment of binary encoding
//
//  3) Perform some testing on the fibonacci encoded bit-packed stream
//      -- encode a bunch of random 64-bit numbers, make sure we can
//      decode back to the same number.
//
//  NOTES: pre/post increment/decrement work modulo whatever size they
//  are.  So, if you have a 6-bit value of zero, and you decrement,
//  you end up with a 6-bit value of all 1's, or 63.


void
testBinaryEncoding(void) {
  time_t     mtseed     = time(0L);
  mt_s      *mtctx      = 0L;

  uint32     iterations = TEST_LENGTH;

  uint64    *bits        = new uint64 [iterations + 2];
  uint64     bpos        = uint64ZERO;

  uint64    *V           = new uint64 [iterations];
  uint64    *C           = new uint64 [iterations];
  uint64    *S           = new uint64 [iterations];

  uint32     failed      = 0;
  uint32     errors      = 0;

  fprintf(stderr, "Starting test of binary encoding\n");

  bpos   = uint64ZERO;
  mtctx = mtInit(mtseed);

  //  Build some values to stuff into the bits

  for (uint32 j=0; j < iterations; j++) {
    S[j] = (mtRandom32(mtctx) % 63) + 1;
    V[j] = mtRandom64(mtctx) & uint64MASK(S[j]);
    //fprintf(stderr, "[%2d] S="uint64FMT" V="uint64HEX"\n", j, S[j], V[j]);
  }

  //  Stuff them in, in blocks of some size.  At the same time, decode
  //  (this has found bugs in the past).

  failed = 0;
  for (uint32 j=0; j < iterations; ) {
    uint64 num = (mtRandom32(mtctx) % 8);

    if (j + num > iterations)
      num = iterations - j;

    if (num == 0) {
      setDecodedValue(bits, bpos, S[j], V[j]);
      C[j] = getDecodedValue(bits, bpos, S[j]);
      //fprintf(stderr, "[%2d] V="uint64HEX" C="uint64HEX" single\n", j, V[j], C[j]);
      bpos += S[j];
    } else {
      uint64 newp1 = setDecodedValues(bits, bpos, num, S+j, V+j);
      uint64 newp2 = getDecodedValues(bits, bpos, num, S+j, C+j);

      if (newp1 != newp2) {
        //  not perfect; we should be checking the values too, but we do that later.
        for (uint32 x=0; x<num; x++)
          fprintf(stderr, "[%2d] #1 V="uint64HEX" C="uint64HEX" multiple "uint32FMT" %s\n",
                  j+x, V[j+x], C[j+x], num, (V[j+x] == C[j+x]) ? "" : "FAILED");
        failed++;
      }

      bpos = newp2;
    }

    j += num;
    if (num == 0)
      j++;
  }
  if (failed) {
    fprintf(stderr, "binEncoding #1 failed encoding "uint32FMT" times.\n", failed);
    errors++;
  }

  //  Check that V == C

  failed = 0;
  for (uint32 j=0; j<iterations; j++) {
    if (V[j] != C[j]) {
      fprintf(stderr, "[%2d] #2 V="uint64HEX" C="uint64HEX" S="uint32FMT"\n",
              j, V[j], C[j], S[j]);
      failed++;
    }
  }
  if (failed) {
    fprintf(stderr, "binEncoding #2 failed encode/decode "uint32FMT" times.\n", failed);
    errors++;
  }


  //  Decode independently, with different nums

  bpos = 0;  //  reset to start of bits

  for (uint32 j=0; j < iterations; ) {
    uint64 num = (mtRandom32(mtctx) % 8);

    if (j + num > iterations)
      num = iterations - j;

    if (num == 0) {
      C[j] = getDecodedValue(bits, bpos, S[j]);
      bpos += S[j];
    } else {
      bpos = getDecodedValues(bits, bpos, num, S+j, C+j);
    }

    j += num;
    if (num == 0)
      j++;
  }

  //  Check that V == C

  failed = 0;
  for (uint32 j=0; j<iterations; j++) {
    if (V[j] != C[j]) {
      fprintf(stderr, "[%2d] #3 V="uint64HEX" C="uint64HEX" S="uint32FMT"\n",
              j, V[j], C[j], S[j]);
      failed++;
    }
  }
  if (failed) {
    fprintf(stderr, "binEncoding #3 failed decoding "uint32FMT" times.\n", failed);
    errors++;
  }

  //  Clean.

  delete [] bits;
  delete [] V;
  delete [] C;
  delete [] S;

  if (errors)
    exit(1);
}




void
testBinaryEncodingPrePost(void) {
  time_t     mtseed     = time(0L);
  mt_s      *mtctx      = 0L;

  uint32     iterations = TEST_LENGTH;

  uint64    *bits        = new uint64 [2 * iterations];
  uint64     bpos        = uint64ZERO;
  uint32     siz1       = uint64ZERO;
  uint64     val1       = uint64ZERO;
  uint64     val2       = uint64ZERO;

  fprintf(stderr, "Starting test of binary encoding pre/post increment\n");

  bpos   = uint64ZERO;
  mtctx = mtInit(mtseed);

  for (uint32 j=0; j < iterations; j++) {
    siz1  = (mtRandom32(mtctx) % 63) + 1;
    val1  = mtRandom64(mtctx) & uint64MASK(siz1);

    setDecodedValue(bits, bpos, siz1, val1);

    val2 = postDecrementDecodedValue(bits, bpos, siz1);
    if (val2 != val1) {
      fprintf(stderr, "postDec1 failed: got "uint64FMT" expected "uint64FMT" siz="uint32FMT"\n",
              val2, val1, siz1);
      exit(1);
    }
    val2  = getDecodedValue(bits, bpos, siz1) + 1;
    val2 &= uint64MASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "postDec2 failed: got "uint64FMT" expected "uint64FMT" siz="uint32FMT"\n",
              val2, val1, siz1);
      exit(1);
    }

    val2 = preDecrementDecodedValue(bits, bpos, siz1) + 2;
    val2 &= uint64MASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "preDec failed: got "uint64FMT" expected "uint64FMT" siz="uint32FMT"\n",
              val2, val1, siz1);
      exit(1);
    }

    val2 = postIncrementDecodedValue(bits, bpos, siz1) + 2;
    val2 &= uint64MASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "postInc failed: got "uint64FMT" expected "uint64FMT"\n", val2+2, val1-2);
      exit(1);
    }
    val2  = getDecodedValue(bits, bpos, siz1) + 1;
    val2 &= uint64MASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "postInc2 failed: got "uint64FMT" expected "uint64FMT" siz="uint32FMT"\n",
              val2, val1, siz1);
      exit(1);
    }

    val2 = preIncrementDecodedValue(bits, bpos, siz1);
    //  Should be back to original value, so no mask
    if (val2 != val1) {
      fprintf(stderr, "preInc failed: got "uint64FMT" expected "uint64FMT"\n", val2, val1);
      exit(1);
    }

    switch (j % 4) {
      case 0:
        val2 = postDecrementDecodedValue(bits, bpos, siz1);
        break;
      case 1:
        val2 = preDecrementDecodedValue(bits, bpos, siz1);
        break;
      case 2:
        val2 = postIncrementDecodedValue(bits, bpos, siz1);
        break;
      case 3:
        val2 = preIncrementDecodedValue(bits, bpos, siz1);
        break;
    }

    bpos += siz1;
  }

  bpos   = uint64ZERO;
  mtctx = mtInit(mtseed);

  //for (j=0; j < iterations; j++) {
  //}

  delete [] bits;
}






void
testFibonacciEncoding(void) {
  time_t     mtseed     = time(0L);
  mt_s      *mtctx      = 0L;

  uint32     iterations = TEST_LENGTH / 4;

  uint64    *bits        = new uint64 [3 * iterations];
  uint64     bpos        = uint64ZERO;

  uint32     failed = 0;
  uint32     errors = 0;

  fprintf(stderr, "Starting test of fibonacci encoding\n");

  bpos   = uint64ZERO;
  mtctx  = mtInit(mtseed);
  failed = 0;

  for (uint32 j=0; j < iterations; j++) {
    uint64 siz1  = (mtRandom32(mtctx) % 63) + 1;
    uint64 val1  = mtRandom64(mtctx) & uint64MASK(siz1);
    uint64 siz2  = siz1;

    setFibonacciEncodedNumber(bits, bpos, &siz1, val1);

    uint64 val2 = getFibonacciEncodedNumber(bits, bpos, &siz2);

    if ((val1 != val2) || (siz1 != siz2)) {
      fprintf(stderr, "fibEnc #1 failed on "uint32FMT": got "uint64FMT" expected "uint64FMT"\n", j, val2, val1);
      failed++;
    }
    bpos += siz1;
  }
  if (failed) {
    fprintf(stderr, "fibEnc #1 failed "uint32FMT" times.\n", failed);
    errors++;
  }

  bpos   = uint64ZERO;
  mtctx  = mtInit(mtseed);
  failed = 0;

  for (uint32 j=0; j < iterations; j++) {
    uint64 siz1  = (mtRandom32(mtctx) % 63) + 1;
    uint64 val1  = mtRandom64(mtctx) & uint64MASK(siz1);
    uint64 val2  = getFibonacciEncodedNumber(bits, bpos, &siz1);

    if (val1 != val2) {
      fprintf(stderr, "fibEnc #2 failed on "uint32FMT": got "uint64FMT" expected "uint64FMT"\n", j, val2, val1);
      failed++;
    }
    bpos += siz1;
  }
  if (failed) {
    fprintf(stderr, "fibEnc #2 failed "uint32FMT" times.\n", failed);
    errors++;
  }

  delete [] bits;

  if (errors)
    exit(1);
}





int
main(int argc, char **argv) {
  testBinaryEncoding();
  testBinaryEncodingPrePost();
  testFibonacciEncoding();
  return(0);
}

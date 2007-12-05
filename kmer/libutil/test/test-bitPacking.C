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

  u32bit     iterations = TEST_LENGTH;

  u64bit    *bits        = new u64bit [iterations + 2];
  u64bit     bpos        = u64bitZERO;

  u64bit    *V           = new u64bit [iterations];
  u64bit    *C           = new u64bit [iterations];
  u64bit    *S           = new u64bit [iterations];

  u32bit     failed      = 0;
  u32bit     errors      = 0;

  fprintf(stderr, "Starting test of binary encoding\n");

  bpos   = u64bitZERO;
  mtctx = mtInit(mtseed);

  //  Build some values to stuff into the bits

  for (u32bit j=0; j < iterations; j++) {
    S[j] = (mtRandom32(mtctx) % 63) + 1;
    V[j] = mtRandom64(mtctx) & u64bitMASK(S[j]);
    //fprintf(stderr, "[%2d] S="u64bitFMT" V="u64bitHEX"\n", j, S[j], V[j]);
  }

  //  Stuff them in, in blocks of some size.  At the same time, decode
  //  (this has found bugs in the past).

  failed = 0;
  for (u32bit j=0; j < iterations; ) {
    u64bit num = (mtRandom32(mtctx) % 8);

    if (j + num > iterations)
      num = iterations - j;

    if (num == 0) {
      setDecodedValue(bits, bpos, S[j], V[j]);
      C[j] = getDecodedValue(bits, bpos, S[j]);
      //fprintf(stderr, "[%2d] V="u64bitHEX" C="u64bitHEX" single\n", j, V[j], C[j]);
      bpos += S[j];
    } else {
      u64bit newp1 = setDecodedValues(bits, bpos, num, S+j, V+j);
      u64bit newp2 = getDecodedValues(bits, bpos, num, S+j, C+j);

      if (newp1 != newp2) {
        //  not perfect; we should be checking the values too, but we do that later.
        for (u32bit x=0; x<num; x++)
          fprintf(stderr, "[%2d] #1 V="u64bitHEX" C="u64bitHEX" multiple "u32bitFMT" %s\n",
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
    fprintf(stderr, "binEncoding #1 failed encoding "u32bitFMT" times.\n", failed);
    errors++;
  }

  //  Check that V == C

  failed = 0;
  for (u32bit j=0; j<iterations; j++) {
    if (V[j] != C[j]) {
      fprintf(stderr, "[%2d] #2 V="u64bitHEX" C="u64bitHEX" S="u32bitFMT"\n",
              j, V[j], C[j], S[j]);
      failed++;
    }
  }
  if (failed) {
    fprintf(stderr, "binEncoding #2 failed encode/decode "u32bitFMT" times.\n", failed);
    errors++;
  }


  //  Decode independently, with different nums

  bpos = 0;  //  reset to start of bits

  for (u32bit j=0; j < iterations; ) {
    u64bit num = (mtRandom32(mtctx) % 8);

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
  for (u32bit j=0; j<iterations; j++) {
    if (V[j] != C[j]) {
      fprintf(stderr, "[%2d] #3 V="u64bitHEX" C="u64bitHEX" S="u32bitFMT"\n",
              j, V[j], C[j], S[j]);
      failed++;
    }
  }
  if (failed) {
    fprintf(stderr, "binEncoding #3 failed decoding "u32bitFMT" times.\n", failed);
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

  u32bit     iterations = TEST_LENGTH;

  u64bit    *bits        = new u64bit [2 * iterations];
  u64bit     bpos        = u64bitZERO;
  u32bit     siz1       = u64bitZERO;
  u64bit     val1       = u64bitZERO;
  u64bit     val2       = u64bitZERO;

  fprintf(stderr, "Starting test of binary encoding pre/post increment\n");

  bpos   = u64bitZERO;
  mtctx = mtInit(mtseed);

  for (u32bit j=0; j < iterations; j++) {
    siz1  = (mtRandom32(mtctx) % 63) + 1;
    val1  = mtRandom64(mtctx) & u64bitMASK(siz1);

    setDecodedValue(bits, bpos, siz1, val1);

    val2 = postDecrementDecodedValue(bits, bpos, siz1);
    if (val2 != val1) {
      fprintf(stderr, "postDec1 failed: got "u64bitFMT" expected "u64bitFMT" siz="u32bitFMT"\n",
              val2, val1, siz1);
      exit(1);
    }
    val2  = getDecodedValue(bits, bpos, siz1) + 1;
    val2 &= u64bitMASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "postDec2 failed: got "u64bitFMT" expected "u64bitFMT" siz="u32bitFMT"\n",
              val2, val1, siz1);
      exit(1);
    }

    val2 = preDecrementDecodedValue(bits, bpos, siz1) + 2;
    val2 &= u64bitMASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "preDec failed: got "u64bitFMT" expected "u64bitFMT" siz="u32bitFMT"\n",
              val2, val1, siz1);
      exit(1);
    }

    val2 = postIncrementDecodedValue(bits, bpos, siz1) + 2;
    val2 &= u64bitMASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "postInc failed: got "u64bitFMT" expected "u64bitFMT"\n", val2+2, val1-2);
      exit(1);
    }
    val2  = getDecodedValue(bits, bpos, siz1) + 1;
    val2 &= u64bitMASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "postInc2 failed: got "u64bitFMT" expected "u64bitFMT" siz="u32bitFMT"\n",
              val2, val1, siz1);
      exit(1);
    }

    val2 = preIncrementDecodedValue(bits, bpos, siz1);
    //  Should be back to original value, so no mask
    if (val2 != val1) {
      fprintf(stderr, "preInc failed: got "u64bitFMT" expected "u64bitFMT"\n", val2, val1);
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

  bpos   = u64bitZERO;
  mtctx = mtInit(mtseed);

  //for (j=0; j < iterations; j++) {
  //}

  delete [] bits;
}






void
testFibonacciEncoding(void) {
  time_t     mtseed     = time(0L);
  mt_s      *mtctx      = 0L;

  u32bit     iterations = TEST_LENGTH / 4;

  u64bit    *bits        = new u64bit [3 * iterations];
  u64bit     bpos        = u64bitZERO;

  u32bit     failed = 0;
  u32bit     errors = 0;

  fprintf(stderr, "Starting test of fibonacci encoding\n");

  bpos   = u64bitZERO;
  mtctx  = mtInit(mtseed);
  failed = 0;

  for (u32bit j=0; j < iterations; j++) {
    u64bit siz1  = (mtRandom32(mtctx) % 63) + 1;
    u64bit val1  = mtRandom64(mtctx) & u64bitMASK(siz1);
    u64bit siz2  = siz1;

    setFibonacciEncodedNumber(bits, bpos, &siz1, val1);

    u64bit val2 = getFibonacciEncodedNumber(bits, bpos, &siz2);

    if ((val1 != val2) || (siz1 != siz2)) {
      fprintf(stderr, "fibEnc #1 failed on "u32bitFMT": got "u64bitFMT" expected "u64bitFMT"\n", j, val2, val1);
      failed++;
    }
    bpos += siz1;
  }
  if (failed) {
    fprintf(stderr, "fibEnc #1 failed "u32bitFMT" times.\n", failed);
    errors++;
  }

  bpos   = u64bitZERO;
  mtctx  = mtInit(mtseed);
  failed = 0;

  for (u32bit j=0; j < iterations; j++) {
    u64bit siz1  = (mtRandom32(mtctx) % 63) + 1;
    u64bit val1  = mtRandom64(mtctx) & u64bitMASK(siz1);
    u64bit val2  = getFibonacciEncodedNumber(bits, bpos, &siz1);

    if (val1 != val2) {
      fprintf(stderr, "fibEnc #2 failed on "u32bitFMT": got "u64bitFMT" expected "u64bitFMT"\n", j, val2, val1);
      failed++;
    }
    bpos += siz1;
  }
  if (failed) {
    fprintf(stderr, "fibEnc #2 failed "u32bitFMT" times.\n", failed);
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

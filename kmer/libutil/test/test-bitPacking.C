#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "bri++.H"

//  An integer multiplier on the test length.  1 is pretty quick, but
//  10 is the default.
//
#define TEST_LENGTH 10



//  Perform some testing on the fibonacci encoded bit-packed stream
//
//  2)  Encode a bunch of random 64-bit numbers, make sure we can
//      decode back to the same number.
//
//  We test
//  1)  binary encoding/decoding
//  2)  pre/post increment of binary encoding
//
//  NOTES: pre/post increment/decrement work modulo whatever size they
//  are.  So, if you have a 6-bit value of zero, and you decrement,
//  you end up with a 6-bit value of all 1's, or 63.
//
//  3)  Fibonacci encoding/decoding
//
//




void
testBinaryEncoding(void) {
  u32bit     iterations = TEST_LENGTH * 1024 * 1024;
  u32bit     j          = u32bitZERO;
  u64bit    *ptr        = new u64bit [2 * iterations];
  u64bit     pos        = u64bitZERO;
  u64bit     siz1       = u64bitZERO;
  u64bit     siz2       = u64bitZERO;
  u64bit     val1       = u64bitZERO;
  u64bit     val2       = u64bitZERO;
  time_t     mtseed     = time(0L);
  mt_s      *mtctx      = mtInit(mtseed);

  fprintf(stderr, "Starting test of binary encoding\n");

  pos   = u64bitZERO;
  mtctx = mtInit(mtseed);

  for (j=0; j < iterations; j++) {
    siz1  = (mtRandom32(mtctx) % 63) + 1;
    val1  = mtRandom64(mtctx) & u64bitMASK(siz1);
    siz2  = siz1;
    setDecodedValue(ptr, pos, siz1, val1);
    val2 = getDecodedValue(ptr, pos, siz1);
    if ((val1 != val2) || (siz1 != siz2)) {
      fprintf(stderr, "binEnc1 failed on "u32bitFMT": got "u64bitFMT" expected "u64bitFMT"\n", j, val2, val1);
      exit(1);
    }
    pos += siz1;
  }

  pos   = u64bitZERO;
  mtctx = mtInit(mtseed);

  for (j=0; j < iterations; j++) {
    siz1  = (mtRandom32(mtctx) % 63) + 1;
    val1  = mtRandom64(mtctx) & u64bitMASK(siz1);
    val2  = getDecodedValue(ptr, pos, siz1);
    if (val1 != val2) {
      fprintf(stderr, "binEnc2 failed on "u32bitFMT": got "u64bitFMT" expected "u64bitFMT"\n", j, val2, val1);
      exit(1);
    }
    pos += siz1;
  }

  delete [] ptr;
}




void
testBinaryEncodingPrePost(void) {
  u32bit     iterations = TEST_LENGTH * 1024 * 1024;
  u32bit     j          = u32bitZERO;
  u64bit    *ptr        = new u64bit [2 * iterations];
  u64bit     pos        = u64bitZERO;
  u32bit     siz1       = u64bitZERO;
  u64bit     val1       = u64bitZERO;
  u64bit     val2       = u64bitZERO;
  time_t     mtseed     = time(0L);
  mt_s      *mtctx      = mtInit(mtseed);

  fprintf(stderr, "Starting test of binary encoding pre/post increment\n");

  pos   = u64bitZERO;
  mtctx = mtInit(mtseed);

  for (j=0; j < iterations; j++) {
    siz1  = (mtRandom32(mtctx) % 63) + 1;
    val1  = mtRandom64(mtctx) & u64bitMASK(siz1);

    setDecodedValue(ptr, pos, siz1, val1);

    val2 = postDecrementDecodedValue(ptr, pos, siz1);
    if (val2 != val1) {
      fprintf(stderr, "postDec1 failed: got "u64bitFMT" expected "u64bitFMT" siz="u32bitFMT"\n",
              val2, val1, siz1);
      exit(1);
    }
    val2  = getDecodedValue(ptr, pos, siz1) + 1;
    val2 &= u64bitMASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "postDec2 failed: got "u64bitFMT" expected "u64bitFMT" siz="u32bitFMT"\n",
              val2, val1, siz1);
      exit(1);
    }

    val2 = preDecrementDecodedValue(ptr, pos, siz1) + 2;
    val2 &= u64bitMASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "preDec failed: got "u64bitFMT" expected "u64bitFMT" siz="u32bitFMT"\n",
              val2, val1, siz1);
      exit(1);
    }

    val2 = postIncrementDecodedValue(ptr, pos, siz1) + 2;
    val2 &= u64bitMASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "postInc failed: got "u64bitFMT" expected "u64bitFMT"\n", val2+2, val1-2);
      exit(1);
    }
    val2  = getDecodedValue(ptr, pos, siz1) + 1;
    val2 &= u64bitMASK(siz1);
    if (val2 != val1) {
      fprintf(stderr, "postInc2 failed: got "u64bitFMT" expected "u64bitFMT" siz="u32bitFMT"\n",
              val2, val1, siz1);
      exit(1);
    }

    val2 = preIncrementDecodedValue(ptr, pos, siz1);
    //  Should be back to original value, so no mask
    if (val2 != val1) {
      fprintf(stderr, "preInc failed: got "u64bitFMT" expected "u64bitFMT"\n", val2, val1);
      exit(1);
    }

    switch (j % 4) {
      case 0:
        val2 = postDecrementDecodedValue(ptr, pos, siz1);
        break;
      case 1:
        val2 = preDecrementDecodedValue(ptr, pos, siz1);
        break;
      case 2:
        val2 = postIncrementDecodedValue(ptr, pos, siz1);
        break;
      case 3:
        val2 = preIncrementDecodedValue(ptr, pos, siz1);
        break;
    }

    pos += siz1;
  }

  pos   = u64bitZERO;
  mtctx = mtInit(mtseed);

  for (j=0; j < iterations; j++) {
  }

  delete [] ptr;
}






void
testFibonacciEncoding(void) {
  u32bit     iterations = TEST_LENGTH * 256 * 1024;
  u32bit     j          = u32bitZERO;
  u64bit    *ptr        = new u64bit [3 * iterations];
  u64bit     pos        = u64bitZERO;
  u64bit     siz1       = u64bitZERO;
  u64bit     siz2       = u64bitZERO;
  u64bit     val1       = u64bitZERO;
  u64bit     val2       = u64bitZERO;
  time_t     mtseed     = time(0L);
  mt_s      *mtctx      = mtInit(mtseed);

  fprintf(stderr, "Starting test of fibonacci encoding\n");

  pos   = u64bitZERO;
  mtctx = mtInit(mtseed);

  for (j=0; j < iterations; j++) {
    siz1  = (mtRandom32(mtctx) % 63) + 1;
    val1  = mtRandom64(mtctx) & u64bitMASK(siz1);
    siz2  = siz1;
    setFibonacciEncodedNumber(ptr, pos, &siz1, val1);
    val2 = getFibonacciEncodedNumber(ptr, pos, &siz2);
    if ((val1 != val2) || (siz1 != siz2)) {
      fprintf(stderr, "fibEnc1 failed on "u32bitFMT": got "u64bitFMT" expected "u64bitFMT"\n", j, val2, val1);
      exit(1);
    }
    pos += siz1;
  }

  pos   = u64bitZERO;
  mtctx = mtInit(mtseed);

  for (j=0; j < iterations; j++) {
    siz1  = (mtRandom32(mtctx) % 63) + 1;
    val1  = mtRandom64(mtctx) & u64bitMASK(siz1);
    val2  = getFibonacciEncodedNumber(ptr, pos, &siz1);
    if (val1 != val2) {
      fprintf(stderr, "fibEnc2 failed on "u32bitFMT": got "u64bitFMT" expected "u64bitFMT"\n", j, val2, val1);
      exit(1);
    }
    pos += siz1;
  }

  delete [] ptr;
}





int
main(int argc, char **argv) {

  testBinaryEncoding();
  testBinaryEncodingPrePost();
  testFibonacciEncoding();

  return(0);
}

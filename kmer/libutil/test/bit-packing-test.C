#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "bit-packing.H"

//  Perform some testing on the fibonacci encoded bit-packed stream
//
//  1)  Encode a bunch of consecutive (small) numbers, decoding as
//      we encode, then decoding after all encoding is done.
//
//  2)  Same as #1, but with large numbers that should span
//      mutiple 64-bit words when they're encoded, and not
//      consecutive.
//
//  3)  Same as #1, but with a bunch of random (kind of) values.
//

int
main(int argc, char **argv) {

  u32bit     iterations = 16 * 1024 * 1024;

  if (argc > 1) {
    iterations = atoi(argv[1]);
  }

  u64bit    *ptr        = new u64bit [2 * iterations];
  u64bit     pos        = 0;
  u64bit     siz1       = 0;
  u64bit     siz2       = 0;
  u64bit     rnd        = 0;
  u64bit     val        = 0;
  u64bit     max        = fibonacciValues[90];



  //  Test #1
  //
  fprintf(stderr, "Starting test #1 of fibonacci encoded values.\n");
  pos = 0;
  for (u64bit i=0; i<iterations; i++) {
    setFibonacciEncodedNumber(ptr, pos, &siz1, i);
    val = getFibonacciEncodedNumber(ptr, pos, &siz2);
    if (val != i) {
      fprintf(stderr, "Fibonacci encoding test FAILED!\n");
      fprintf(stderr, "val=%lu != correct=%lu  (siz1=%lu siz2=%lu)\n", val, i, siz1, siz2);
      exit(1);
    }
    if (siz1 != siz2) {
      fprintf(stderr, "Fibonacci encoding test FAILED!\n");
      fprintf(stderr, "val=%lu == correct=%lu, but siz1=%lu != siz2=%lu\n", val, i, siz1, siz2);
      exit(1);
    }

    pos += siz1;
  }

  pos = 0;

  for (u64bit i=0; i<iterations; i++) {
    val = getFibonacciEncodedNumber(ptr, pos, &siz1);
    if (val != i) {
      fprintf(stderr, "Fibonacci encoding test FAILED!\n");
      fprintf(stderr, "val=%lu != correct=%lu\n", val, i);
      exit(1);
    }
    pos += siz1;
  }


  //  Test #2
  //
  fprintf(stderr, "Starting test #2 of fibonacci encoded values.\n");
  pos = 0;
  for (u64bit i=20609890381824llu, j=0; j < iterations; i += 13573457, j++) {
    setFibonacciEncodedNumber(ptr, pos, &siz1, i);
    val = getFibonacciEncodedNumber(ptr, pos, &siz2);
    if (val != i) {
      fprintf(stderr, "Fibonacci encoding test FAILED!\n");
      fprintf(stderr, "val=%lu != correct=%lu  (siz1=%lu siz2=%lu)\n", val, i, siz1, siz2);
      exit(1);
    }
    if (siz1 != siz2) {
      fprintf(stderr, "Fibonacci encoding test FAILED!\n");
      fprintf(stderr, "val=%lu == correct=%lu, but siz1=%lu != siz2=%lu\n", val, i, siz1, siz2);
      exit(1);
    }

    pos += siz1;
  }

  pos = 0;

  for (u64bit i=20609890381824llu, j=0; j < iterations; i += 13573457, j++) {
    val = getFibonacciEncodedNumber(ptr, pos, &siz1);
    if (val != i) {
      fprintf(stderr, "Fibonacci encoding test FAILED!\n");
      fprintf(stderr, "val=%lu != correct=%lu\n", val, i);
      exit(1);
    }
    pos += siz1;
  }



  //  Test #3
  //
  pos = 0;
  fprintf(stderr, "Starting test #3 of fibonacci encoded values.\n");

  long   seed = time(NULL);
  srand48(seed);

  for (u64bit i=0; i<iterations; i++) {
    u64bit rnd;

    rnd   = lrand48();
    rnd <<= 32;
    rnd  |= lrand48();

    setFibonacciEncodedNumber(ptr, pos, &siz1, rnd);
    val = getFibonacciEncodedNumber(ptr, pos, &siz2);
    if (val != rnd) {
      fprintf(stderr, "Fibonacci encoding test FAILED!\n");
      fprintf(stderr, "val=%lu != correct=%lu  (siz1=%lu siz2=%lu)\n", val, rnd, siz1, siz2);
      exit(1);
    }
    if (siz1 != siz2) {
      fprintf(stderr, "Fibonacci encoding test FAILED!\n");
      fprintf(stderr, "val=%lu == correct=%lu, but siz1=%lu != siz2=%lu\n", val, rnd, siz1, siz2);
      exit(1);
    }

    pos += siz1;
  }

  srand48(seed);
  pos = 0;

  for (u64bit i=0; i<iterations; i++) {
    u64bit rnd;

    rnd   = lrand48();
    rnd <<= 32;
    rnd  |= lrand48();

    val = getFibonacciEncodedNumber(ptr, pos, &siz1);
    if (val != rnd) {
      fprintf(stderr, "Fibonacci encoding test FAILED!\n");
      fprintf(stderr, "val=%lu != correct=%lu\n", val, rnd);
      exit(1);
    }
    pos += siz1;
  }

  fprintf(stderr, "All tests passed!\n");

  return(0);
}

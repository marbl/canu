#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "bitPackedFile.H"
#include "bit-packing.H"


//  Two tests are defined (OK, three): one that tests fixed-width bit
//  packing, one that tests the fibonacci numbers, and one that tests
//  both.
//
#define TESTSIZE         4929482
#define TEST_BITS
#define TEST_NUMBERS

void
main(void) {
  u32bit    i;
  u32bit   *siz = new u32bit [TESTSIZE];
  u64bit   *val = new u64bit [TESTSIZE];
  u32bit    errs = 0;

  srand48(time(NULL));
  for (i=0; i<TESTSIZE; i++) {
    siz[i] =  (u32bit)floor(drand48() * 64) + 1;
    val[i] = ((u64bit)lrand48()) << 32 | ((u64bit)lrand48());
    val[i] &= u64bitMASK(siz[i]);
  }

  bitPackedFileWriter *W = new bitPackedFileWriter("bittest.junk");
  for (i=0; i<TESTSIZE; i++) {
#ifdef TEST_BITS
    W->putBits(val[i], siz[i]);
#endif

#ifdef TEST_NUMBERS
    W->putNumber(val[i]);
#endif
  }
  delete W;

  bitPackedFileReader *R = new bitPackedFileReader("bittest.junk");
  for (i=0; i<TESTSIZE; i++) {
#ifdef TEST_BITS
    u64bit r = R->getBits(siz[i]);
    if (r != val[i]) {
      fprintf(stderr, "%6d] ERROR in getBits()   -- retrieved 0x%016lx != expected 0x%016lx (%2d bits).\n", i, r, val[i], siz[i]);
      errs++;
#if 0
    } else {
      fprintf(stderr, "%6d]       in getBits()   -- retrieved 0x%016lx != expected 0x%016lx (%2d bits).\n", i, r, val[i], siz[i]);
#endif
    }
#endif

#ifdef TEST_NUMBERS
    u64bit v = R->getNumber();
    if (v != val[i]) {
      fprintf(stderr, "%6d] ERROR in getNumber() -- retrieved 0x%016lx != expected 0x%016lx.\n", i, v, val[i]);
      errs++;
#if 0
    } else {
      fprintf(stderr, "%6d]       in getNumber() -- retrieved 0x%016lx != expected 0x%016lx.\n", i, v, val[i]);
#endif
    }
#endif
  }
  delete R;

  fprintf(stderr, "There are %u errors.\n", errs);
}

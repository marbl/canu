#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "bit-packing.H"

#define TESTSIZE 4929482

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
  for (i=0; i<TESTSIZE; i++)
    W->putBits(val[i], siz[i]);
  delete W;

  bitPackedFileReader *R = new bitPackedFileReader("bittest.junk");
  for (i=0; i<TESTSIZE; i++) {
    u64bit r = R->getBits(siz[i]);

    if (r != val[i]) {
      fprintf(stderr, "%6d] ERROR -- 0x%016lx != 0x%016lx (%2d).\n", i, r, val[i], siz[i]);
      errs++;
    }
  }
  delete R;

  fprintf(stderr, "There are %u errors.\n", errs);
  for (u32bit i=0; i<=64; i++) {
    u32bit c=0;
    for (u32bit j=0; j<TESTSIZE; j++) {
      if (siz[j] == i)
        c++;
    }
    fprintf(stderr, "%8d at size = %d\n", c, i);
  }
}

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>

#include "bp.H"
#include "bitPackedFile.H"


int
main(int argc, char **argv) {
  u64bit   *A = new u64bit [1024];
  u64bit    w = u64bitZERO;

  long seed = time(NULL);
  srand48(seed);
  fprintf(stderr, "Using seed = %ld\n", seed);

  fprintf(stderr, "\nMASKS\n");
  fprintf(stderr, "64 -- 0x%016lx\n", u64bitMASK(64));
  fprintf(stderr, "63 -- 0x%016lx\n", u64bitMASK(63));
  fprintf(stderr, "33 -- 0x%016lx\n", u64bitMASK(33));
  fprintf(stderr, "32 -- 0x%016lx\n", u64bitMASK(32));
  fprintf(stderr, "31 -- 0x%016lx\n", u64bitMASK(31));
  fprintf(stderr, "01 -- 0x%016lx\n", u64bitMASK(1));
  fprintf(stderr, "00 -- 0x%016lx\n", u64bitMASK(0));   //  this should generate a warning

  fprintf(stderr, "\nARRAY\n");
  setDecodedValue(A, w, 32, 0x00000000ffffffff);  w += 32;
  setDecodedValue(A, w, 32, 0x00000000eeeeeeee);  w += 32;
  setDecodedValue(A, w, 32, 0x00000000dddddddd);  w += 32;
  setDecodedValue(A, w, 32, 0x00000000cccccccc);  w += 32;
  setDecodedValue(A, w, 32, 0x00000000bbbbbbbb);  w += 32;
  setDecodedValue(A, w, 32, 0x00000000aaaaaaaa);  w += 32;

  w = 0;

  fprintf(stderr, "0x%016lx\n", getDecodedValue(A, w, 32));  w += 32;
  fprintf(stderr, "0x%016lx\n", getDecodedValue(A, w, 32));  w += 32;
  fprintf(stderr, "0x%016lx\n", getDecodedValue(A, w, 32));  w += 32;
  fprintf(stderr, "0x%016lx\n", getDecodedValue(A, w, 32));  w += 32;
  fprintf(stderr, "0x%016lx\n", getDecodedValue(A, w, 32));  w += 32;
  fprintf(stderr, "0x%016lx\n", getDecodedValue(A, w, 32));  w += 32;


  u64bit   val[0x003fffff];
  u64bit   siz[0x003fffff];

  for (u64bit i=0; i<0x00000000003fffff; i++) {
    val[i] =  (lrand48() >> 6);
    siz[i] = ((lrand48() >> 6) & 0x3e) + 1;
  }

  fprintf(stderr, "\nFILE\n");
  bitPackedFileWriter  *W = new bitPackedFileWriter("test-bp-file");
  W->putBits(0xaaaaaaaaaaaaaaaa, 32);
  W->putBits(0xbbbbbbbbbbbbbbbb, 64);
  W->putBits(0xcccccccccccccccc, 32);
  W->putBits(0xdddddddddddddddd, 64);
  W->putBits(0xeeeeeeeeeeeeeeee, 32);
  W->putBits(0xffffffffffffffff, 64);
  for (u64bit i=0; i<0x00000000000fffff; i++) {
    W->putBits(i, 32);
  }
  for (u64bit i=0; i<0x00000000003fffff; i++) {
    W->putBits(val[i], siz[i]);
  }
  delete W;


  bitPackedFileReader  *R = new bitPackedFileReader("test-bp-file");
  fprintf(stderr, "0x%016lx\n", R->getBits(32));
  fprintf(stderr, "0x%016lx\n", R->getBits(64));
  fprintf(stderr, "0x%016lx\n", R->getBits(32));
  fprintf(stderr, "0x%016lx\n", R->getBits(64));
  fprintf(stderr, "0x%016lx\n", R->getBits(32));
  fprintf(stderr, "0x%016lx\n", R->getBits(64));
  for (u64bit i=0; i<0x00000000000fffff; i++) {
    if (i != R->getBits(32)) {
      fprintf(stderr, "Error at i=%lu\n", i);
      exit(1);
    }
  }
  for (u64bit i=0; i<0x00000000003fffff; i++) {
    if ((val[i] & u64bitMASK(siz[i])) != R->getBits(siz[i])) {
      fprintf(stderr, "Error at i=%lu (val=0x%016lx  siz=%lu)\n", i, val[i], siz[i]);
      exit(1);
    }
  }
  delete R;

  unlink("test-bp-file");
}

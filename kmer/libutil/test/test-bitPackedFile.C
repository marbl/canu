#include <unistd.h>
#include <time.h>
#include <math.h>

#include "util++.H"

u32bit   testSize = 5000000;
u32bit   testIter = 200;

mt_s  *mtctx;


//  Generate a list of random 64-bit numbers, remember the number and the size
//
void
generateRandom(u32bit *siz, u64bit *val) {

  for (u32bit i=0; i<testSize; i++) {
#if 0
    //  For debugging
    siz[i] =  13;
    val[i] =  (i % 2) ? u64bitZERO : ~u64bitZERO;
#else
    siz[i] = (mtRandom32(mtctx) % 63) + 1;
    val[i] = mtRandom64(mtctx);
#endif
    val[i] &= u64bitMASK(siz[i]);
  }
}



void
testStreaming(void) {
  bitPackedFile *F = 0L;
  u32bit         i;
  u32bit        *siz = new u32bit [testSize];
  u64bit        *val = new u64bit [testSize];
  u32bit         errs = 0;

  generateRandom(siz, val);

  //  Write those numbers to a bitPackedFile, both binary encoded and
  //  fibonacci encoded.
  //
  F = new bitPackedFile("bittest.junk");
  for (i=0; i<testSize; i++) {
    F->putBits(val[i], siz[i]);
    F->putNumber(val[i]);
  }
  delete F;

  //  Open the file and check what we just wrote.
  //
  F = new bitPackedFile("bittest.junk");
  for (i=0; i<testSize; i++) {
    u64bit v;

    v = F->getBits(siz[i]);
    if (v != val[i]) {
      fprintf(stderr, u32bitFMT"] ERROR in getBits()   -- retrieved "u64bitHEX" != expected "u64bitHEX" ("u32bitFMT" bits).\n", i, v, val[i], siz[i]);
      errs++;
    }

    v = F->getNumber();
    if (v != val[i]) {
      fprintf(stderr, u32bitFMT"] ERROR in getNumber() -- retrieved "u64bitHEX" != expected "u64bitHEX".\n", i, v, val[i]);
      errs++;
    }
  }
  delete F;

  delete [] val;
  delete [] siz;

  if (errs > 0) {
    fprintf(stderr, "There are "u32bitFMT" errors in the stream test.\n", errs);
    exit(1);
  } else {
    fprintf(stderr, "The stream test PASSED.\n");
  }

  unlink("bittest.junk");
}


void
testRandomReading(void) {
  bitPackedFile *F = 0L;
  u32bit         i;
  u32bit        *siz = new u32bit [testSize];
  u64bit        *val = new u64bit [testSize];
  u32bit         errs = 0;

  generateRandom(siz, val);

  //  Create a new bitpacked file, writing just numbers as binary encoded.
  //
  F = new bitPackedFile("bittest.junk");
  for (i=0; i<testSize; i++)
    F->putBits(val[i], siz[i]);
  delete F;


  //  Do several seek tests.  Seek to a random element, and read it.
  //
  F = new bitPackedFile("bittest.junk");
  for (i=0; i<testIter; i++) {
    u32bit idx = (u32bit)lrand48() % testSize;
    u64bit pos = 0;

    for (u32bit j=0; j<idx; j++)
      pos += siz[j];

    F->seek(pos);
    u64bit r = F->getBits(siz[idx]);

    if (r != val[idx]) {
      fprintf(stderr, u32bitFMT"] ERROR in seek()/getBits()   -- retrieved "u64bitHEX" != expected "u64bitHEX" ("u32bitFMT" bits).\n", i, r, val[i], siz[i]);
      errs++;
    }
  }
  delete F;

  delete [] val;
  delete [] siz;

  if (errs > 0) {
    fprintf(stderr, "There are "u32bitFMT" errors in the random access.\n", errs);
    exit(1);
  } else {
    fprintf(stderr, "The seek test PASSED.\n");
  }

  unlink("bittest.junk");
}





void
testReWrite(void) {
  bitPackedFile *F = 0L;
  u32bit         i;
  u32bit        *siz = new u32bit [testSize];
  u64bit        *val = new u64bit [testSize];
  u32bit         errs = 0;
  u64bit         pos = u64bitZERO;

  generateRandom(siz, val);

  //  First, write zeros to the file
  //
  F = new bitPackedFile("bittest.junk");
  for (i=0; i<testSize; i++)
    F->putBits(u64bitZERO, siz[i]);
  delete F;

  fprintf(stderr, "WRITING FORWARDS!\n");

  //  Now, write every other number to the file
  //
  F = new bitPackedFile("bittest.junk");
  for (i=0; i<testSize; i++) {
    if ((i % 2) == 1) {
      F->seek(pos);
      F->putBits(val[i], siz[i]);
    }
    pos += siz[i];
  }
  F->showStats(stderr);
  delete F;

  fprintf(stderr, "WRITING BACKWARDS!\n");

  //  And go backwards and write the other set of numbers to the file
  //
  F = new bitPackedFile("bittest.junk");
  for (i=testSize; i--; ) {
    pos -= siz[i];
    if ((i % 2) == 0) {
      F->seek(pos);
      F->putBits(val[i], siz[i]);
    }
  }
  F->showStats(stderr);
  delete F;

  //  Now, stream through the file and see if we wrote what we should have
  //
  F = new bitPackedFile("bittest.junk");
  for (i=0; i<testSize; i++) {
    u64bit v;

    v = F->getBits(siz[i]);
    if (v != val[i]) {
      fprintf(stderr, u32bitFMT"] ERROR in seekstream/getBits()   -- retrieved "u64bitHEX" != expected "u64bitHEX" ("u32bitFMT" bits).\n", i, v, val[i], siz[i]);
      errs++;
    }
  }
  F->showStats(stderr);
  delete F;

  delete [] val;
  delete [] siz;

  if (errs > 0) {
    fprintf(stderr, "There are "u32bitFMT" errors in the rewrite test.\n", errs);
    exit(1);
  } else {
    fprintf(stderr, "The rewrite test PASSED.\n");
  }

  unlink("bittest.junk");
}



int
main(int argc, char **argv) {

  if (argc != 3) {
    fprintf(stderr, "usage: %s testSize testIter\n", argv[0]);
    fprintf(stderr, "  This will perform various tests on the bitPackedFile class,\n");
    fprintf(stderr, "  returning 0 if OK and 1 if error.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  testSize -- the number of words to use in a write then read test\n");
    fprintf(stderr, "  testIter -- the number of random access tests to do\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "I'll assume reasonable values and continue\n");
    fprintf(stderr, "  testSize = "u32bitFMT"\n", testSize);
    fprintf(stderr, "  testIter = "u32bitFMT"\n", testIter);
  } else {
    testSize = strtou32bit(argv[1], 0L);
    testIter = strtou32bit(argv[2], 0L);
  }

  mtctx = mtInit(time(NULL));

  testStreaming();
  testRandomReading();
  testReWrite();
}


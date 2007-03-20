#include <unistd.h>
#include <time.h>
#include <math.h>

#include "util++.H"

//  This will perform various tests on the bitPackedFile class,
//  returning 0 if OK and 1 if error.
//
//  testSize -- the number of words to use in a write then read test
//  testIter -- the number of random access tests to do

u32bit   testSize = 2000000;
u32bit   testIter = 50;

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
testRandomReading(bool inCore) {
  bitPackedFile *F = 0L;
  u32bit         i;
  u32bit        *siz = new u32bit [testSize + 1];
  u64bit        *val = new u64bit [testSize];
  u32bit         errs = 0;

  fprintf(stderr, "BUILDING random test set.\n");
  generateRandom(siz, val);

  //  Create a new bitpacked file, writing just numbers as binary encoded.
  //
  fprintf(stderr, "SAVING random test set.\n");
  F = new bitPackedFile("bittest.junk");
  for (i=0; i<testSize; i++)
    F->putBits(val[i], siz[i]);
  delete F;

  //  Covert the siz[] into offsets
  //
  u32bit t = siz[0];
  siz[0] = 0;
  for (u32bit i=1; i<testSize; i++) {
    u32bit x = siz[i];
    siz[i] = t;
    t += x;
  }
  siz[testSize] = t;

  //  Attempt to flush memory
  //
  {
    u32bit  ll = 400 * 1024 * 1024 / 8;
    u64bit *xx = new u64bit [ll];
    xx[0] = 1;
    xx[1] = 1;
    for (u32bit i=2; i<ll; i++)
      xx[i] = xx[i-1] + xx[i-2];
    fprintf(stdout, "FLUSHED: "u32bitFMT"\n", xx[ll-1]);
    delete [] xx;
  }

  //  Do several seek tests.  Seek to a random element, and read it.
  //
  F = new bitPackedFile("bittest.junk");

  if (inCore) {
    F->loadInCore();
    fprintf(stderr, "Begin INCORE seek test!\n");
  } else {
    fprintf(stderr, "Begin DISKBASED seek test!\n");
  }

  double  startTime = getTime();

  for (i=0; i<testIter; i++) {
    u32bit idx = (u32bit)lrand48() % testSize;

    F->seek(siz[idx]);
    u64bit r = F->getBits(siz[idx+1] - siz[idx]);

    if (r != val[idx]) {
      fprintf(stderr, u32bitFMT"] ERROR in seek()/getBits()   -- retrieved "u64bitHEX" != expected "u64bitHEX" ("u32bitFMT" bits).\n", i, r, val[i], siz[i]);
      errs++;
    }
  }
  delete F;

  if (errs > 0) {
    fprintf(stderr, "There are "u32bitFMT" errors in the %s random access.\n", errs, (inCore) ? "inCore" : "disk");
    exit(1);
  } else {
    fprintf(stderr, "The %s seek test PASSED (%f seconds).\n",
            (inCore) ? "inCore" : "disk",
            getTime() - startTime);
  }

  delete [] val;
  delete [] siz;

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

  mtctx = mtInit(time(NULL));

  testSize = 30000000;
  testIter = 2000;
  //testStreaming();
  //testReWrite();

  testSize = 40000000;
  testIter = 10000;
  testRandomReading(false);
  testRandomReading(true);
}

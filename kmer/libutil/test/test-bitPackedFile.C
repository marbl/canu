
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
 *    Brian P. Walenz on 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2006-OCT-22 to 2014-APR-11
 *      are Copyright 2006-2007,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <unistd.h>
#include <time.h>
#include <math.h>

#include "util++.H"

//  This will perform various tests on the bitPackedFile class,
//  returning 0 if OK and 1 if error.
//
//  testSize -- the number of words to use in a write then read test
//  testIter -- the number of random access tests to do

uint32   testSize = 2000000;
uint32   testIter = 50;

mt_s  *mtctx;

//  Generate a list of random 64-bit numbers, remember the number and the size
//
void
generateRandom(uint32 *siz, uint64 *val) {

  for (uint32 i=0; i<testSize; i++) {
#if 0
    //  For debugging
    siz[i] =  13;
    val[i] =  (i % 2) ? uint64ZERO : ~uint64ZERO;
#else
    siz[i] = (mtRandom32(mtctx) % 63) + 1;
    val[i] = mtRandom64(mtctx);
#endif
    val[i] &= uint64MASK(siz[i]);
  }
}


void
testStreaming(void) {
  bitPackedFile *F = 0L;
  uint32         i;
  uint32        *siz = new uint32 [testSize];
  uint64        *val = new uint64 [testSize];
  uint32         errs = 0;

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
    uint64 v;

    v = F->getBits(siz[i]);
    if (v != val[i]) {
      fprintf(stderr, uint32FMT"] ERROR in getBits()   -- retrieved "uint64HEX" != expected "uint64HEX" ("uint32FMT" bits).\n", i, v, val[i], siz[i]);
      errs++;
    }

    v = F->getNumber();
    if (v != val[i]) {
      fprintf(stderr, uint32FMT"] ERROR in getNumber() -- retrieved "uint64HEX" != expected "uint64HEX".\n", i, v, val[i]);
      errs++;
    }
  }
  delete F;

  delete [] val;
  delete [] siz;

  if (errs > 0) {
    fprintf(stderr, "There are "uint32FMT" errors in the stream test.\n", errs);
    exit(1);
  } else {
    fprintf(stderr, "The stream test PASSED.\n");
  }

  unlink("bittest.junk");
}


void
testRandomReading(bool inCore) {
  bitPackedFile *F = 0L;
  uint32         i;
  uint32        *siz = new uint32 [testSize + 1];
  uint64        *val = new uint64 [testSize];
  uint32         errs = 0;

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
  uint32 t = siz[0];
  siz[0] = 0;
  for (uint32 i=1; i<testSize; i++) {
    uint32 x = siz[i];
    siz[i] = t;
    t += x;
  }
  siz[testSize] = t;

  //  Attempt to flush memory
  //
  {
    uint32  ll = 400 * 1024 * 1024 / 8;
    uint64 *xx = new uint64 [ll];
    xx[0] = 1;
    xx[1] = 1;
    for (uint32 i=2; i<ll; i++)
      xx[i] = xx[i-1] + xx[i-2];
    fprintf(stdout, "FLUSHED: "uint32FMT"\n", xx[ll-1]);
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
    uint32 idx = (uint32)lrand48() % testSize;

    F->seek(siz[idx]);
    uint64 r = F->getBits(siz[idx+1] - siz[idx]);

    if (r != val[idx]) {
      fprintf(stderr, uint32FMT"] ERROR in seek()/getBits()   -- retrieved "uint64HEX" != expected "uint64HEX" ("uint32FMT" bits).\n", i, r, val[i], siz[i]);
      errs++;
    }
  }
  delete F;

  if (errs > 0) {
    fprintf(stderr, "There are "uint32FMT" errors in the %s random access.\n", errs, (inCore) ? "inCore" : "disk");
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
  uint32         i;
  uint32        *siz = new uint32 [testSize];
  uint64        *val = new uint64 [testSize];
  uint32         errs = 0;
  uint64         pos = uint64ZERO;

  generateRandom(siz, val);

  //  First, write zeros to the file
  //
  F = new bitPackedFile("bittest.junk");
  for (i=0; i<testSize; i++)
    F->putBits(uint64ZERO, siz[i]);
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
    uint64 v;

    v = F->getBits(siz[i]);
    if (v != val[i]) {
      fprintf(stderr, uint32FMT"] ERROR in seekstream/getBits()   -- retrieved "uint64HEX" != expected "uint64HEX" ("uint32FMT" bits).\n", i, v, val[i], siz[i]);
      errs++;
    }
  }
  F->showStats(stderr);
  delete F;

  delete [] val;
  delete [] siz;

  if (errs > 0) {
    fprintf(stderr, "There are "uint32FMT" errors in the rewrite test.\n", errs);
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

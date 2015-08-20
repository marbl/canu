
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
 *    Brian P. Walenz from 2004-APR-27 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2008-JAN-29 to 2014-APR-11
 *      are Copyright 2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "util++.H"


uint64  numLoops = 1;
uint64  numNums  = 4000000;
uint64  numSize  = 300;

//  The space in bits that we can play with, and the pointer to said space.
//
uint64  spa = 128 * 1024 * 1024 * 8;
uint64 *ptr = 0L;
uint64 *rnd = 0L;

void
testUnary(void) {
  uint64  pos = uint64ZERO;
  uint64  siz = uint64ZERO;
  uint64  val = uint64ZERO;
  uint64  i   = uint64ZERO;

  for (i=0; i<numNums; i++) {
    setUnaryEncodedNumber(ptr, pos, &siz, rnd[i]);
    pos += siz;
    if (pos + 1000 >= spa) {
      fprintf(stderr, "ERROR:  Ran out of space in testUnary at number "uint64FMT" out of "uint64FMT"\n", i, numNums);
      exit(1);
    }
  }

  //fprintf(stderr, "unaryEncodedNumbers used "uint64FMT"MB of storage out of "uint64FMT"MB.\n", pos >> 23, spa >> 23);

  pos = uint64ZERO;

  for (i=0; i<numNums; i++) {
    val = getUnaryEncodedNumber(ptr, pos, &siz);
    if (val != rnd[i]) {
      fprintf(stderr, "Number "uint64FMT" at bitpos "uint64FMT" failed.  Desired "uint64FMT" got "uint64FMT"\n", i, pos, rnd[i], val);
      exit(1);
    }
    pos += siz;
  }

  fprintf(stderr, "unary encoded numbers OK!\n");
}



void
testGeneralizedUnary(void) {
  uint64  pos = uint64ZERO;
  uint64  siz = uint64ZERO;
  uint64  val = uint64ZERO;
  uint64  i   = uint64ZERO;

  for (i=0; i<numNums; i++) {
    setGeneralizedUnaryEncodedNumber(ptr, pos, &siz, rnd[i]);
    pos += siz;
    if (pos + 1000 >= spa) {
      fprintf(stderr, "ERROR:  Ran out of space in testGeneralizedUnary at number "uint64FMT" out of "uint64FMT"\n", i, numNums);
      exit(1);
    }
  }

  //fprintf(stderr, "generalizedUnaryEncodedNumbers used "uint64FMT"MB of storage out of "uint64FMT"MB.\n", pos >> 23, spa >> 23);

  pos = uint64ZERO;

  for (i=0; i<numNums; i++) {
    val = getGeneralizedUnaryEncodedNumber(ptr, pos, &siz);
    if (val != rnd[i]) {
      fprintf(stderr, "Number "uint64FMT" at bitpos "uint64FMT" failed.  Desired "uint64FMT" got "uint64FMT"\n", i, pos, rnd[i], val);
      exit(1);
    }
    pos += siz;
  }

  fprintf(stderr, "generalized unary encoded numbers OK!\n");
}




void
testEliasGamma(void) {
  uint64  pos = uint64ZERO;
  uint64  siz = uint64ZERO;
  uint64  val = uint64ZERO;
  uint64  i   = uint64ZERO;

  for (i=0; i<numNums; i++) {
    setEliasGammaEncodedNumber(ptr, pos, &siz, rnd[i]);
    pos += siz;
    if (pos + 1000 >= spa) {
      fprintf(stderr, "ERROR:  Ran out of space in testGeneralizedUnary at number "uint64FMT" out of "uint64FMT"\n", i, numNums);
      exit(1);
    }
  }

  //fprintf(stderr, "eliasGammaEncodedNumbers used "uint64FMT"MB of storage out of "uint64FMT"MB.\n", pos >> 23, spa >> 23);

  pos = uint64ZERO;

  for (i=0; i<numNums; i++) {
    val = getEliasGammaEncodedNumber(ptr, pos, &siz);
    if (val != rnd[i]) {
      fprintf(stderr, "Number "uint64FMT" at bitpos "uint64FMT" failed.  Desired "uint64FMT" got "uint64FMT"\n", i, pos, rnd[i], val);
      exit(1);
    }
    pos += siz;
  }

  fprintf(stderr, "Elias gamma encoded numbers OK!\n");
}



void
testEliasDelta(void) {
  uint64  pos = uint64ZERO;
  uint64  siz = uint64ZERO;
  uint64  val = uint64ZERO;
  uint64  i   = uint64ZERO;

  for (i=0; i<numNums; i++) {
    setEliasDeltaEncodedNumber(ptr, pos, &siz, rnd[i]);
    pos += siz;
    if (pos + 1000 >= spa) {
      fprintf(stderr, "ERROR:  Ran out of space in testGeneralizedUnary at number "uint64FMT" out of "uint64FMT"\n", i, numNums);
      exit(1);
    }
  }

  //fprintf(stderr, "eliasDeltaEncodedNumbers used "uint64FMT"MB of storage out of "uint64FMT"MB.\n", pos >> 23, spa >> 23);

  pos = uint64ZERO;

  for (i=0; i<numNums; i++) {
    val = getEliasDeltaEncodedNumber(ptr, pos, &siz);
    if (val != rnd[i]) {
      fprintf(stderr, "Number "uint64FMT" at bitpos "uint64FMT" failed.  Desired "uint64FMT" got "uint64FMT"\n", i, pos, rnd[i], val);
      exit(1);
    }
    pos += siz;
  }

  fprintf(stderr, "Elias delta encoded numbers OK!\n");
}





int
main(int argc, char **argv) {

  if (argc != 3) {
    fprintf(stderr, "usage: %s <num-loops> <num-nums-per-loop>\n", argv[0]);
    fprintf(stderr, "  -> DEFAULTS USED <-\n");
  } else {
    numLoops = strtouint32(argv[1], 0L);
    numNums  = strtouint32(argv[2], 0L);
  }

  rnd = new uint64 [numNums];
  ptr = new uint64 [spa >> 6];

  mt_s *ctx = mtInit(time(NULL));

  //  Generate some random numbers to store
  //
  while (numLoops--) {

    //  Test out unary encodings on small numbers
    //
    for (uint64 i=0; i<numNums; i++)
      rnd[i] = mtRandom32(ctx) % numSize;
    testUnary();

    //  Generalized unary encoding can handle larger numbers
    //
    for (uint64 i=0; i<numNums; i++)
      rnd[i] = mtRandom32(ctx);
    testGeneralizedUnary();

    //  Elias Gamma and Delta codes are probably pretty good
    //
    for (uint64 i=0; i<numNums; i++)
      rnd[i] = mtRandom64(ctx);
    testEliasGamma();
    testEliasDelta();
  }

  delete [] rnd;
  delete [] ptr;

  exit(0);
}



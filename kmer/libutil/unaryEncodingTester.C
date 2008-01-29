#include "util++.H"


u64bit  numLoops = 1;
u64bit  numNums  = 4000000;
u64bit  numSize  = 300;

//  The space in bits that we can play with, and the pointer to said space.
//
u64bit  spa = 128 * 1024 * 1024 * 8;
u64bit *ptr = 0L;
u64bit *rnd = 0L;

void
testUnary(void) {
  u64bit  pos = u64bitZERO;
  u64bit  siz = u64bitZERO;
  u64bit  val = u64bitZERO;
  u64bit  i   = u64bitZERO;

  for (i=0; i<numNums; i++) {
    setUnaryEncodedNumber(ptr, pos, &siz, rnd[i]);
    pos += siz;
    if (pos + 1000 >= spa) {
      fprintf(stderr, "ERROR:  Ran out of space in testUnary at number "u64bitFMT" out of "u64bitFMT"\n", i, numNums);
      exit(1);
    }
  }

  //fprintf(stderr, "unaryEncodedNumbers used "u64bitFMT"MB of storage out of "u64bitFMT"MB.\n", pos >> 23, spa >> 23);

  pos = u64bitZERO;

  for (i=0; i<numNums; i++) {
    val = getUnaryEncodedNumber(ptr, pos, &siz);
    if (val != rnd[i]) {
      fprintf(stderr, "Number "u64bitFMT" at bitpos "u64bitFMT" failed.  Desired "u64bitFMT" got "u64bitFMT"\n", i, pos, rnd[i], val);
      exit(1);
    }
    pos += siz;
  }

  fprintf(stderr, "unary encoded numbers OK!\n");
}



void
testGeneralizedUnary(void) {
  u64bit  pos = u64bitZERO;
  u64bit  siz = u64bitZERO;
  u64bit  val = u64bitZERO;
  u64bit  i   = u64bitZERO;

  for (i=0; i<numNums; i++) {
    setGeneralizedUnaryEncodedNumber(ptr, pos, &siz, rnd[i]);
    pos += siz;
    if (pos + 1000 >= spa) {
      fprintf(stderr, "ERROR:  Ran out of space in testGeneralizedUnary at number "u64bitFMT" out of "u64bitFMT"\n", i, numNums);
      exit(1);
    }
  }

  //fprintf(stderr, "generalizedUnaryEncodedNumbers used "u64bitFMT"MB of storage out of "u64bitFMT"MB.\n", pos >> 23, spa >> 23);

  pos = u64bitZERO;

  for (i=0; i<numNums; i++) {
    val = getGeneralizedUnaryEncodedNumber(ptr, pos, &siz);
    if (val != rnd[i]) {
      fprintf(stderr, "Number "u64bitFMT" at bitpos "u64bitFMT" failed.  Desired "u64bitFMT" got "u64bitFMT"\n", i, pos, rnd[i], val);
      exit(1);
    }
    pos += siz;
  }

  fprintf(stderr, "generalized unary encoded numbers OK!\n");
}




void
testEliasGamma(void) {
  u64bit  pos = u64bitZERO;
  u64bit  siz = u64bitZERO;
  u64bit  val = u64bitZERO;
  u64bit  i   = u64bitZERO;

  for (i=0; i<numNums; i++) {
    setEliasGammaEncodedNumber(ptr, pos, &siz, rnd[i]);
    pos += siz;
    if (pos + 1000 >= spa) {
      fprintf(stderr, "ERROR:  Ran out of space in testGeneralizedUnary at number "u64bitFMT" out of "u64bitFMT"\n", i, numNums);
      exit(1);
    }
  }

  //fprintf(stderr, "eliasGammaEncodedNumbers used "u64bitFMT"MB of storage out of "u64bitFMT"MB.\n", pos >> 23, spa >> 23);

  pos = u64bitZERO;

  for (i=0; i<numNums; i++) {
    val = getEliasGammaEncodedNumber(ptr, pos, &siz);
    if (val != rnd[i]) {
      fprintf(stderr, "Number "u64bitFMT" at bitpos "u64bitFMT" failed.  Desired "u64bitFMT" got "u64bitFMT"\n", i, pos, rnd[i], val);
      exit(1);
    }
    pos += siz;
  }

  fprintf(stderr, "Elias gamma encoded numbers OK!\n");
}



void
testEliasDelta(void) {
  u64bit  pos = u64bitZERO;
  u64bit  siz = u64bitZERO;
  u64bit  val = u64bitZERO;
  u64bit  i   = u64bitZERO;

  for (i=0; i<numNums; i++) {
    setEliasDeltaEncodedNumber(ptr, pos, &siz, rnd[i]);
    pos += siz;
    if (pos + 1000 >= spa) {
      fprintf(stderr, "ERROR:  Ran out of space in testGeneralizedUnary at number "u64bitFMT" out of "u64bitFMT"\n", i, numNums);
      exit(1);
    }
  }

  //fprintf(stderr, "eliasDeltaEncodedNumbers used "u64bitFMT"MB of storage out of "u64bitFMT"MB.\n", pos >> 23, spa >> 23);

  pos = u64bitZERO;

  for (i=0; i<numNums; i++) {
    val = getEliasDeltaEncodedNumber(ptr, pos, &siz);
    if (val != rnd[i]) {
      fprintf(stderr, "Number "u64bitFMT" at bitpos "u64bitFMT" failed.  Desired "u64bitFMT" got "u64bitFMT"\n", i, pos, rnd[i], val);
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
    numLoops = strtou32bit(argv[1], 0L);
    numNums  = strtou32bit(argv[2], 0L);
  }

  rnd = new u64bit [numNums];
  ptr = new u64bit [spa >> 6];

  mt_s *ctx = mtInit(time(NULL));

  //  Generate some random numbers to store
  //
  while (numLoops--) {

    //  Test out unary encodings on small numbers
    //
    for (u64bit i=0; i<numNums; i++)
      rnd[i] = mtRandom32(ctx) % numSize;
    testUnary();

    //  Generalized unary encoding can handle larger numbers
    //
    for (u64bit i=0; i<numNums; i++)
      rnd[i] = mtRandom32(ctx);
    testGeneralizedUnary();

    //  Elias Gamma and Delta codes are probably pretty good
    //
    for (u64bit i=0; i<numNums; i++)
      rnd[i] = mtRandom64(ctx);
    testEliasGamma();
    testEliasDelta();
  }

  delete [] rnd;
  delete [] ptr;

  exit(0);
}



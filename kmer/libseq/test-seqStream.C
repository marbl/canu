#include "util.h"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

#include "test-correctSequence.H"


u32bit
testIndexing(u32bit numSeq, char sep, u32bit sepLen) {
  u32bit     err = 0;
  seqStream *ST  = new seqStream("test-correctSequence.fasta");

  ST->setSeparator(sep, sepLen);

  u32bit  maxLen = ST->startOf(numSeq-1) + ST->lengthOf(numSeq-1);

  //  Basic checks on the reverse lookup

  if (ST->sequenceNumberOfPosition(maxLen) != ~u32bitZERO) {
    fprintf(stderr, "maxLen 1 wrong.\n");
    err++;
  }
  if (ST->sequenceNumberOfPosition(maxLen - 1) == ~u32bitZERO) {
    fprintf(stderr, "maxLen 2 wrong.\n");
    err++;
  }

  //  Check all lookups

  u32bit p = 0;
  for (u32bit sid=0; sid<numSeq; sid++) {
    for (u32bit ppp=0; ppp<correctSequence[sid].sequenceLength; ppp++, p++) {
      if (ST->sequenceNumberOfPosition(p) != sid) {
        fprintf(stderr, "sequenceNumberOfPosition "u32bitFMT" returned "u32bitFMT", not correct "u32bitFMT".\n",
                p, ST->sequenceNumberOfPosition(p), sid);
        err++;
      }
    }
  }
  if (p != maxLen) {
    fprintf(stderr, "maxLen wrong.\n");
    err++;
  }

  //  Check the separator.  Seek to a spot right before one, and count
  //  that we have the correct length.

  ST->setRange(ST->startOf(numSeq/3) + ST->lengthOf(numSeq/3)-1, ~u32bitZERO);
  ST->get();
  for (u32bit x=0; x<sepLen; x++) {
    char s = 'x';
    if (ST->get() != sep) {
      fprintf(stderr, "wrong separator at sep "u32bitFMT" got %d expected %d\n", x, s, sep);
      err++;
    }
  }

  


  delete ST;

  return(err);
}






int
main(int argc, char **argv) {
  u32bit     minLen = 100;
  u32bit     maxLen = 20000;
  u32bit     numSeq = 1000;
  seqStream *ST     = 0L;
  u32bit     err    = 0;

  generateCorrectSequence(minLen, maxLen, numSeq);

  err += testIndexing(numSeq, '.', 1);
  err += testIndexing(numSeq, ':', 10);
  err += testIndexing(numSeq, 'z', 100);
  err += testIndexing(numSeq, '-', 1000);

  //  1 - string
  //  2 - file, separator length 1
  //  3 - file, separator length 100
  //  4 - file, lots of short sequences
  //  5 - file, few long sequences
  //  6 - file, one long sequence

  //  a - limit (0, middle)
  //  b - limit (middle, end)
  //  c - limit (middle, middle)
  //  d - rewind & random limits in a loop

  //  operationally, build the expected sequence chain, then compare
  //  the stream (get() & eof()) against that chain.  While doing
  //  that, check seqPos(), seqIID(), strPos().


  removeCorrectSequence(numSeq);

  exit(err > 0);
}

#include "util.h"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

#include "test-correctSequence.H"

#define FAIL() { err++; assert(0); }

#warning HOW DO WE TEST IF WE GET ALL THE MERS?


u32bit
testMerStreamSimple(merStream *MS, u32bit merSize, char *seq, u32bit *SP) {
  u32bit   err = 0;
  u32bit   pos = 0;
  char     testmer[32];
  bool     verbose = true;
  bool     nm      = false;

  if (verbose)
    fprintf(stdout, "testMSsimple() begins.\n");

  //  Until we have no more mers in the input

  while (seq[pos + merSize - 1] != 0) {
    nm = MS->nextMer();

    MS->theFMer().merToString(testmer);

    if (verbose) {
      fprintf(stdout, "MS pos="u32bitFMT" posInSeq="u64bitFMT" posInStr="u64bitFMT" seqNum="u64bitFMT"\n",
              pos,
              MS->thePositionInSequence(),
              MS->thePositionInStream(),
              MS->theSequenceNumber());
      if (strncmp(testmer, seq + pos, merSize))
        fprintf(stdout, "MS pos="u32bitFMT" failed '%s' != '%s'.\n", pos, testmer, seq + pos);
    }

    assert(nm == true);
    assert(MS->thePositionInSequence() == SP[pos]);
    assert(MS->thePositionInStream()   == SP[pos]);
    assert(MS->theSequenceNumber()     == 0);
    assert(strncmp(testmer, seq + pos, merSize) == 0);

    pos++;
  }

  //  Should have no more mers

  nm = MS->nextMer();
  assert(nm == false);

  return(err);
}



u32bit
testMerStreamOperation(merStream *MS, u32bit beg, u32bit end, u32bit sepLen) {
  u32bit  err = 0;

  char fmerstr[256];
  char rmerstr[256];
  char cmerstr[256];
  char tmerstr[256];

  while (MS->nextMer()) {
    MS->theFMer().merToString(fmerstr);
    MS->theRMer().merToString(rmerstr);
    MS->theCMer().merToString(cmerstr);

    if ((strcmp(fmerstr, cmerstr) != 0) && ((strcmp(rmerstr, cmerstr) != 0))) {
      fprintf(stderr, "mer strings disagree; F:%s R:%s C:%s\n", fmerstr, rmerstr, cmerstr);
      FAIL();
    }

    reverseComplementSequence(rmerstr, strlen(rmerstr));

    if (strcmp(fmerstr, rmerstr) != 0) {
      fprintf(stderr, "mer strings disagree after reverse; F:%s R:%s\n", fmerstr, rmerstr);
      FAIL();
    }

    u32bit  pseq = MS->thePositionInSequence();
    u32bit  pstr = MS->thePositionInStream();
    u32bit  piid = MS->theSequenceNumber();

    u32bit  mersize = MS->theFMer().getMerSize();
    u32bit  merspan = MS->theFMer().getMerSpan();

#if 0
    if (beg > 10) {
      u32bit  pp = pstr + piid * sepLen - 10;
      u32bit  xx = 0;

      fprintf(stderr, "beg="u32bitFMT" pstr="u32bitFMT" '", beg, pstr);

      for (xx=0; xx<10; xx++, pp++)
        fprintf(stderr, "%c", chainSeq[pp]);
      fprintf(stderr, ":");
      for (xx=0; xx<merspan; xx++, pp++)
        fprintf(stderr, "%c", chainSeq[pp]);
      fprintf(stderr, ":");
      for (xx=0; xx<10; xx++, pp++)
        fprintf(stderr, "%c", chainSeq[pp]);

      fprintf(stderr, "'\n");
    }
#endif

    if (mersize == merspan) {
      strncpy(tmerstr, correctSequence[piid].sequence + pseq, mersize);
      tmerstr[mersize] = 0;

      if (strcmp(fmerstr, tmerstr) != 0) {
        fprintf(stderr, "mer string doesn't agree with sequence; '%s' vs '%s'.\n", fmerstr, tmerstr);
        FAIL();
      }

      if ((pstr < beg) || (end < pstr)) {
        fprintf(stderr, "mer stream position out of range; at "u32bitFMT", range "u32bitFMT"-"u32bitFMT"\n",
                pstr, beg, end);
        FAIL();
      }

      //  The pstr returned above is the ACGT position, not the
      //  chainSeq position we expect.  Trusting that the IID is
      //  correct (if not, the previous strcmp() would have failed) we
      //  can add in the missing separators to get a chainSeq
      //  position.

      strncpy(tmerstr, chainSeq + pstr + piid * sepLen, mersize);
      tmerstr[mersize] = 0;

      if (strcmp(fmerstr, tmerstr) != 0) {
        fprintf(stderr, "mer string doesn't agree with stream; '%s' vs '%s'.\n", fmerstr, tmerstr);
        FAIL();
      }
    }
  }

  return(err);
}




u32bit
testMerStream(kMerBuilder *KB, u32bit numSeq, char sep, u32bit sepLen) {
  u32bit      err = 0;
  seqStream  *ST  = 0L;
  merStream  *MS  = 0L;

  generateChainedAnswer(numSeq, sep, sepLen);

  if (numSeq > 1) {
    ST = new seqStream("test-correctSequence.fasta");
    ST->setSeparator(sep, sepLen);
  } else {
    ST = new seqStream(correctSequence[0].sequence, correctSequence[0].sequenceLength);
  }

  MS = new merStream(KB, ST, true, true);

  u32bit  maxLen = ST->startOf(numSeq-1) + ST->lengthOf(numSeq-1);

  //  Whole thing, rewind, whole thing

  fprintf(stderr, "whole thing.\n");

  err += testMerStreamOperation(MS, 0, maxLen, sepLen);
  MS->rewind();
  err += testMerStreamOperation(MS, 0, maxLen, sepLen);


  //  Random subsets - we're not terribly interested in streaming,
  //  just getting the start/end correct.

  fprintf(stderr, "subsets.\n");

  for (u32bit iter=0; iter<500; iter++) {
    u32bit beg = mtRandom32(mtctx) % maxLen;
    u32bit end = (beg + 10000 < maxLen) ? (beg + 10000) : maxLen;

    //fprintf(stderr, "subsets - "u32bitFMT"-"u32bitFMT"\n", beg, end);

    MS->setBaseRange(beg, end);

    err += testMerStreamOperation(MS, beg, end, sepLen);
    MS->rewind();
    err += testMerStreamOperation(MS, beg, end, sepLen);
  }

  delete MS;

  return(err);
}





int
main(int argc, char **argv) {
  u32bit     minLen = 1000;
  u32bit     maxLen = 200000;
  u32bit     numSeq = 1000;
  u32bit     err    = 0;

  //  Very simple merStream test

  {
    fprintf(stdout, "merStream(kMerBuilder(20), ...)\n");

    merStream *MS = new merStream(new kMerBuilder(20),
                                  new seqStream("GGGTCAACTCCGCCCGCACTCTAGC", 25),
                                  true, true);
    u32bit     SP[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    testMerStreamSimple(MS, 20, "GGGTCAACTCCGCCCGCACTCTAGC", SP);
    MS->rewind();
    testMerStreamSimple(MS, 20, "GGGTCAACTCCGCCCGCACTCTAGC", SP);
    MS->rewind();
    MS->rewind();
    testMerStreamSimple(MS, 20, "GGGTCAACTCCGCCCGCACTCTAGC", SP);

    delete MS;

    fprintf(stdout, "merStream(kMerBuilder(20), ...) - PASSED\n");
  }

  {
    fprintf(stdout, "merStream(kMerBuilder(20, 1), ...)\n");

    merStream *MS = new merStream(new kMerBuilder(20, 1),
                                  new seqStream("GGGAATTTTCAACTCCGCCCGCACTCTAGCCCAAA", 35),
                                  true, true);
    u32bit     SP[10] = { 0, 3, 5, 9, 10, 12 };

    testMerStreamSimple(MS, 20, "GATCACTCGCGCACTCTAGCA", SP);
    MS->rewind();
    testMerStreamSimple(MS, 20, "GATCACTCGCGCACTCTAGCA", SP);
    MS->rewind();
    MS->rewind();
    testMerStreamSimple(MS, 20, "GATCACTCGCGCACTCTAGCA", SP);

    delete MS;

    fprintf(stdout, "merStream(kMerBuilder(20, 1), ...) - PASSED\n");
  }

  //  Move on to harder tests

  generateCorrectSequence(minLen, maxLen, numSeq);

  //  Tests seqStream(string, strlen) construction method

  fprintf(stderr, "err += testMerStream(new kMerBuilder(20, 0, 0L),    1, '.', 1);\n");
  err += testMerStream(new kMerBuilder(20, 0, 0L),    1, '.', 1);

  fprintf(stderr, "err += testMerStream(new kMerBuilder(22, 1, 0L),    1, '.', 1);\n");
  err += testMerStream(new kMerBuilder(22, 1, 0L),    1, '.', 1);

  //  Tests seqStream(filename) construction method

  fprintf(stderr, "err += testMerStream(new kMerBuilder(20, 0, 0L), numSeq, '.',   1);\n");
  err += testMerStream(new kMerBuilder(20, 0, 0L), numSeq, '.',   1);

  fprintf(stderr, "err += testMerStream(new kMerBuilder(28, 0, 0L), numSeq, '.', 100);\n");
  err += testMerStream(new kMerBuilder(28, 0, 0L), numSeq, '.', 100);

  fprintf(stderr, "err += testMerStream(new kMerBuilder(24, 4, 0L), numSeq, '.', 100);\n");
  err += testMerStream(new kMerBuilder(24, 4, 0L), numSeq, '.', 100);

  removeCorrectSequence(numSeq);

  if (err == 0)
    fprintf(stderr, "Success!\n");

  exit(err > 0);
}

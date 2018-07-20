
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "util.h"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

#include "test-correctSequence.H"

#define FAIL() { err++; assert(0); }

#warning HOW DO WE TEST IF WE GET ALL THE MERS?


uint32
testMerStreamSimple(merStream *MS, uint32 merSize, char *seq, uint32 *SP) {
  uint32   err = 0;
  uint32   pos = 0;
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
      fprintf(stdout, "MS pos="uint32FMT" posInSeq="uint64FMT" posInStr="uint64FMT" seqNum="uint64FMT"\n",
              pos,
              MS->thePositionInSequence(),
              MS->thePositionInStream(),
              MS->theSequenceNumber());
      if (strncmp(testmer, seq + pos, merSize))
        fprintf(stdout, "MS pos="uint32FMT" failed '%s' != '%s'.\n", pos, testmer, seq + pos);
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



uint32
testMerStreamOperation(merStream *MS, uint32 beg, uint32 end, uint32 sepLen) {
  uint32  err = 0;

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

    uint32  pseq = MS->thePositionInSequence();
    uint32  pstr = MS->thePositionInStream();
    uint32  piid = MS->theSequenceNumber();

    uint32  mersize = MS->theFMer().getMerSize();
    uint32  merspan = MS->theFMer().getMerSpan();

#if 0
    if (beg > 10) {
      uint32  pp = pstr + piid * sepLen - 10;
      uint32  xx = 0;

      fprintf(stderr, "beg="uint32FMT" pstr="uint32FMT" '", beg, pstr);

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
        fprintf(stderr, "mer stream position out of range; at "uint32FMT", range "uint32FMT"-"uint32FMT"\n",
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




uint32
testMerStream(kMerBuilder *KB, uint32 numSeq, char sep, uint32 sepLen) {
  uint32      err = 0;
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

  uint32  maxLen = ST->startOf(numSeq-1) + ST->lengthOf(numSeq-1);

  //  Whole thing, rewind, whole thing

  fprintf(stderr, "whole thing.\n");

  err += testMerStreamOperation(MS, 0, maxLen, sepLen);
  MS->rewind();
  err += testMerStreamOperation(MS, 0, maxLen, sepLen);


  //  Random subsets - we're not terribly interested in streaming,
  //  just getting the start/end correct.

  fprintf(stderr, "subsets.\n");

  for (uint32 iter=0; iter<500; iter++) {
    uint32 beg = mtRandom32(mtctx) % maxLen;
    uint32 end = (beg + 10000 < maxLen) ? (beg + 10000) : maxLen;

    //fprintf(stderr, "subsets - "uint32FMT"-"uint32FMT"\n", beg, end);

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
  uint32     minLen = 1000;
  uint32     maxLen = 200000;
  uint32     numSeq = 1000;
  uint32     err    = 0;

  //  Very simple merStream test

  {
    fprintf(stdout, "merStream(kMerBuilder(20), ...)\n");

    merStream *MS = new merStream(new kMerBuilder(20),
                                  new seqStream("GGGTCAACTCCGCCCGCACTCTAGC", 25),
                                  true, true);
    uint32     SP[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

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
    uint32     SP[10] = { 0, 3, 5, 9, 10, 12 };

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

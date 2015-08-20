
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


uint32
testIndexing(uint32 numSeq, char sep, uint32 sepLen) {
  uint32     err = 0;
  seqStream *ST  = 0L;

  fprintf(stderr, "testIndexing()-- numSeq="uint32FMT" sep=%c sepLen="uint32FMT"\n", numSeq, sep, sepLen);

  generateChainedAnswer(numSeq, sep, sepLen);

  if (numSeq > 1) {
    ST = new seqStream("test-correctSequence.fasta");
    ST->setSeparator(sep, sepLen);
  } else {
    ST = new seqStream(correctSequence[0].sequence, correctSequence[0].sequenceLength);
  }

  uint32  maxLen = ST->startOf(numSeq-1) + ST->lengthOf(numSeq-1);

  //  Basic checks on the reverse lookup - this is state independent;
  //  it changes only based on the separator length.  In other words,
  //  there is no need to check this while iterating through the
  //  seqStream.

  fprintf(stderr, "IGNORE THIS WARNING:  ");
  if (ST->sequenceNumberOfPosition(maxLen) != ~uint32ZERO) {
    fprintf(stderr, "maxLen too small.\n");
    FAIL();
  }
  if (ST->sequenceNumberOfPosition(maxLen - 1) == ~uint32ZERO) {
    fprintf(stderr, "maxLen too big.\n");
    FAIL();
  }

  //  Check all lookups - lengthOf() and IIDOf() are implicitly
  //  checked by the operation of seqStream (get() mostly).  startOf()
  //  isn't, but inserting errors in setRange() led to
  //  infinite-looking loops.

  uint64 pos = 0;
  uint64 sta = 0;

  for (uint32 sid=0; sid<numSeq; sid++) {
    if (ST->lengthOf(sid) != correctSequence[sid].sequenceLength) {
      fprintf(stderr, "lengthOf "uint32FMT" returned "uint32FMT", not correct "uint32FMT"\n",
              sid, ST->lengthOf(sid), correctSequence[sid].sequenceLength);
      FAIL();
    }
    if (ST->startOf(sid)  != sta) {
      fprintf(stderr, "startOf "uint32FMT" returned "uint64FMT", not correct "uint64FMT"\n",
              sid, ST->startOf(sid), sta);
      FAIL();
    }
    if (ST->IIDOf(sid)    != sid) {
      fprintf(stderr, "IIDOf "uint32FMT" returned "uint32FMT", not correct "uint32FMT"\n",
              sid, ST->IIDOf(sid), sid);
      FAIL();
    }

    sta += correctSequence[sid].sequenceLength;

    for (uint32 ppp=0; ppp<correctSequence[sid].sequenceLength; ppp++, pos++) {
      if (ST->sequenceNumberOfPosition(pos) != sid) {
        fprintf(stderr, "sequenceNumberOfPosition "uint64FMT" returned "uint32FMT", not correct "uint32FMT".\n",
                pos, ST->sequenceNumberOfPosition(pos), sid);
        FAIL();
      }
    }
  }
  if (pos != maxLen) {
    fprintf(stderr, "maxLen wrong.\n");
    FAIL();
  }

  //  Check the separator.  Seek to a spot right before one, and count
  //  that we have the correct length.  More rigorously tested in
  //  testChaining().

  for (uint32 sid=0; sid<numSeq-1; sid++) {
    ST->setRange(ST->startOf(sid) + ST->lengthOf(sid)-1, ~uint64ZERO);
    ST->get();
    for (uint32 x=0; x<sepLen; x++) {
      char s = ST->get();
      if (s != sep) {
        fprintf(stderr, "wrong separator at sep "uint32FMT" got %d expected %d\n", x, s, sep);
        FAIL();
      }
    }
    if (ST->get() == sep) {
      fprintf(stderr, "too many separators!\n");
      FAIL();
    }
  }

  delete ST;

  return(err);
}



uint32
testSeqStream(seqStream *ST, uint32 sib, uint32 sie, char sep) {
  uint32  err = 0;

  while (ST->eof() == false) {
    uint32   sp = ST->seqPos();
    uint32   si = ST->seqIID();
    uint64   st = ST->strPos();
    char     ch = ST->get();

    if (ch != 0) {
      if (ch != chainSeq[sib]) {
        fprintf(stderr, "sp="uint32FMT" si="uint32FMT" st="uint64FMT" ch=%c -- letter wrong got'%c'\n", sp, si, st, ch, chainSeq[sib]);
        FAIL();
      }
      if ((ch != sep) && (sp != chainSeqPos[sib])) {
        fprintf(stderr, "sp="uint32FMT" si="uint32FMT" st="uint64FMT" ch=%c -- seqPos wrong got "uint32FMT"\n", sp, si, st, ch, chainSeqPos[sib]);
        FAIL();
      }
      if ((ch != sep) && (si != chainSeqIID[sib])) {
        fprintf(stderr, "sp="uint32FMT" si="uint32FMT" st="uint64FMT" ch=%c -- seqIID wrong got"uint32FMT"\n", sp, si, st, ch, chainSeqIID[sib]);
        FAIL();
      }
      if ((ch != sep) && (st != chainStrPos[sib])) {
        fprintf(stderr, "sp="uint32FMT" si="uint32FMT" st="uint64FMT" ch=%c -- strPos wrong got "uint64FMT"\n", sp, si, st, ch, chainStrPos[sib]);
        FAIL();
      }

      sib++;
    }
  }

  if (sib != sie) {
    fprintf(stderr, "iterated length wrong; sib="uint32FMT" sie="uint32FMT"\n", sib, sie);
    FAIL();
  }

  return(err);
}



uint32
testChaining(uint32 numSeq, char sep, uint32 sepLen) {
  uint32     err = 0;
  seqStream *ST  = 0L;

  fprintf(stderr, "testChaining()-- numSeq="uint32FMT" sep=%c sepLen="uint32FMT"\n", numSeq, sep, sepLen);

  generateChainedAnswer(numSeq, sep, sepLen);

  if (numSeq > 1) {
    ST = new seqStream("test-correctSequence.fasta");
    ST->setSeparator(sep, sepLen);
  } else {
    ST = new seqStream(correctSequence[0].sequence, correctSequence[0].sequenceLength);
  }

  //  Do a test on the whole thing.

  {
    uint32  sib = 0;
    uint32  sie = strlen(chainSeq);

    fprintf(stderr, "initial test with full range\n");
    testSeqStream(ST, sib, sie, sep);

    fprintf(stderr, "initial test with full range (rewind)\n");
    ST->rewind();
    testSeqStream(ST, sib, sie, sep);
  }


  //  Set the range to random values, and check all the results.
  //  We've already verified the index works, so we're free to use
  //  that (but we currently don't).

  uint32  maxLen = ST->startOf(numSeq-1) + ST->lengthOf(numSeq-1);

  fprintf(stderr, "test on subranges\n");

  for (uint32 iter=0; iter<500; iter++) {
    uint32 beg = mtRandom32(mtctx) % maxLen;
    uint32 end = mtRandom32(mtctx) % maxLen;
    if (beg > end) {
      uint32 t = end;
      end = beg;
      beg = t;
    }

    ST->setRange(beg, end);

    //  Compute the position in our stream for the ACGT based beg and
    //  end.  The quirk here is that our stream includes the
    //  separator.

    uint32 sib = 0;  //  chainSeq position
    uint32 sie = 0;

    for (uint32 ppp=0, sid=0; sid<numSeq; sid++) {
      uint32 len = correctSequence[sid].sequenceLength;

      if ((ppp <= beg) && (beg < ppp + len)) {
        sib += beg - ppp;
        break;
      }

      ppp += len;
      sib += len + sepLen;
    }

    for (uint32 ppp=0, sid=0; sid<numSeq; sid++) {
      uint32 len = correctSequence[sid].sequenceLength;

      if ((ppp <= end) && (end < ppp + len)) {
        sie += end - ppp;
        break;
      }

      ppp += len;
      sie += len + sepLen;
    }

    //  Optionally do a rewind in the middle

    if (iter % 2) {
      //fprintf(stderr, "Random iter "uint32FMT" (with rewind)\n", iter);
      while (ST->eof() == false)
        ST->get();
      ST->rewind();
    } else {
      //fprintf(stderr, "Random iter "uint32FMT"\n", iter);
    }


    testSeqStream(ST, sib, sie, sep);
  }

  return(err > 0);
}



int
main(int argc, char **argv) {
  uint32     minLen = 100;
  uint32     maxLen = 20000;
  uint32     numSeq = 1000;
  uint32     err    = 0;

  generateCorrectSequence(minLen, maxLen, numSeq);

  //  Tests seqStream(string, strlen) construction method

  err += testIndexing(1, '.', 1);
  err += testChaining(1, '.', 1);

  //  Tests seqStream(filename) construction method

  err += testIndexing(numSeq, '.', 1);
  err += testIndexing(numSeq, ':', 10);
  err += testIndexing(numSeq, 'z', 100);
  err += testIndexing(numSeq, '-', 1000);

  err += testChaining(numSeq, '.', 1);
  err += testChaining(numSeq, ':', 10);
  err += testChaining(numSeq, 'z', 100);
  err += testChaining(numSeq, '-', 1000);

  removeCorrectSequence(numSeq);

  if (err == 0)
    fprintf(stderr, "Success!\n");

  exit(err > 0);
}

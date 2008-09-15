#include "util.h"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

#include "test-correctSequence.H"


char     *chainSeq;
u32bit   *chainSeqPos;
u32bit   *chainSeqIID;
u32bit   *chainStrPos;


u32bit
testIndexing(u32bit numSeq, char sep, u32bit sepLen) {
  u32bit     err = 0;
  seqStream *ST  = 0L;

  fprintf(stderr, "testIndexing()-- numSeq="u32bitFMT" sep=%c sepLen="u32bitFMT"\n", numSeq, sep, sepLen);

  if (numSeq > 1) {
    ST = new seqStream("test-correctSequence.fasta");
    ST->setSeparator(sep, sepLen);
  } else {
    ST = new seqStream(correctSequence[0].sequence, correctSequence[0].sequenceLength);
  }

  u32bit  maxLen = ST->startOf(numSeq-1) + ST->lengthOf(numSeq-1);

  //  Basic checks on the reverse lookup - this is state independent;
  //  it changes only based on the separator length.  In other words,
  //  there is no need to check this while iterating through the
  //  seqStream.

  fprintf(stderr, "IGNORE THIS WARNING:  ");
  if (ST->sequenceNumberOfPosition(maxLen) != ~u32bitZERO) {
    fprintf(stderr, "maxLen 1 wrong.\n");
    err++;
  }
  if (ST->sequenceNumberOfPosition(maxLen - 1) == ~u32bitZERO) {
    fprintf(stderr, "maxLen 2 wrong.\n");
    err++;
  }

  //  Check all lookups - lengthOf() and IIDOf() are implicitly
  //  checked by the operation of seqStream (get() mostly).  startOf()
  //  isn't, but inserting errors in setRange() led to
  //  infinite-looking loops.

  u32bit pos = 0;
  u32bit sta = 0;

  for (u32bit sid=0; sid<numSeq; sid++) {
    if (ST->lengthOf(sid) != correctSequence[sid].sequenceLength) {
      fprintf(stderr, "lengthOf "u32bitFMT" returned "u32bitFMT", not correct "u32bitFMT"\n",
              sid, ST->lengthOf(sid), correctSequence[sid].sequenceLength);
      err++;
    }
    if (ST->startOf(sid)  != sta) {
      fprintf(stderr, "lengthOf "u32bitFMT" returned "u32bitFMT", not correct "u32bitFMT"\n",
              sid, ST->startOf(sid), sta);
      err++;
    }
    if (ST->IIDOf(sid)    != sid) {
      fprintf(stderr, "IIDOf "u32bitFMT" returned "u32bitFMT", not correct "u32bitFMT"\n",
              sid, ST->IIDOf(sid), sid);
      err++;
    }

    sta += correctSequence[sid].sequenceLength;

    for (u32bit ppp=0; ppp<correctSequence[sid].sequenceLength; ppp++, pos++) {
      if (ST->sequenceNumberOfPosition(pos) != sid) {
        fprintf(stderr, "sequenceNumberOfPosition "u32bitFMT" returned "u32bitFMT", not correct "u32bitFMT".\n",
                pos, ST->sequenceNumberOfPosition(pos), sid);
        err++;
      }
    }
  }
  if (pos != maxLen) {
    fprintf(stderr, "maxLen wrong.\n");
    err++;
  }

  //  Check the separator.  Seek to a spot right before one, and count
  //  that we have the correct length.  More rigorously tested in
  //  testChaining().

  for (u32bit sid=0; sid<numSeq-1; sid++) {
    ST->setRange(ST->startOf(sid) + ST->lengthOf(sid)-1, ~u32bitZERO);
    ST->get();
    for (u32bit x=0; x<sepLen; x++) {
      char s = ST->get();
      if (s != sep) {
        fprintf(stderr, "wrong separator at sep "u32bitFMT" got %d expected %d\n", x, s, sep);
        err++;
      }
    }
    if (ST->get() == sep) {
      fprintf(stderr, "too many separators!\n");
      err++;
    }
  }

  delete ST;

  return(err);
}



u32bit
testChaining(u32bit numSeq, char sep, u32bit sepLen) {
  u32bit     err = 0;
  seqStream *ST  = 0L;

  fprintf(stderr, "testChaining()-- numSeq="u32bitFMT" sep=%c sepLen="u32bitFMT"\n", numSeq, sep, sepLen);

  if (numSeq > 1) {
    ST = new seqStream("test-correctSequence.fasta");
    ST->setSeparator(sep, sepLen);
  } else {
    ST = new seqStream(correctSequence[0].sequence, correctSequence[0].sequenceLength);
  }

  //  Build the answer

  //  Purposely allocate more space than real, because we _ALWAYS_
  //  tack on an extra separator at the end of the chain we build
  //  below.

  {
    u32bit  maxLen = ST->startOf(numSeq-1) + ST->lengthOf(numSeq-1) + numSeq * sepLen + sepLen;

    chainSeq    = new char   [maxLen];
    chainSeqPos = new u32bit [maxLen];
    chainSeqIID = new u32bit [maxLen];
    chainStrPos = new u32bit [maxLen];
  }

  u32bit p      = 0;
  u32bit strpos = 0;

  for (u32bit sid=0; sid<numSeq; sid++) {
    for (u32bit ppp=0; ppp<correctSequence[sid].sequenceLength; ppp++, p++, strpos++) {
      chainSeq[p]    = correctSequence[sid].sequence[ppp];
      chainSeqPos[p] = ppp;
      chainSeqIID[p] = sid;
      chainStrPos[p] = strpos;
    }
    if (sid+1 < numSeq)
      for (u32bit ppp=0; ppp<sepLen; ppp++, p++) {
        chainSeq[p]    = sep;
        chainSeqPos[p] = ~u32bitZERO;
        chainSeqIID[p] = ~u32bitZERO;
        chainStrPos[p] = ~u32bitZERO;
      }
  }

  //  Do a test on the whole thing.

  {
    u32bit  sib = 0;
    u32bit  sie = p;

    fprintf(stderr, "initial test with full range\n");

    while (ST->eof() == false) {
      u64bit   sp = ST->seqPos();
      u64bit   si = ST->seqIID();
      u64bit   st = ST->strPos();
      char     ch = ST->get();

      if (ch != 0) {
        if (ch != chainSeq[sib]) {
          fprintf(stderr, "letter wrong '%c' vs '%c'\n", ch, chainSeq[sib]);
          err++;
        }
        if ((ch != sep) && (sp != chainSeqPos[sib])) {
          fprintf(stderr, "seqPos wrong "u32bitFMT" vs "u32bitFMT"\n", sp, chainSeqPos[sib]);
          err++;
        }
        if ((ch != sep) && (si != chainSeqIID[sib])) {
          fprintf(stderr, "seqIID wrong "u32bitFMT" vs "u32bitFMT"\n", si, chainSeqIID[sib]);
          err++;
        }
        if ((ch != sep) && (st != chainStrPos[sib])) {
          fprintf(stderr, "strPos wrong "u32bitFMT" vs "u32bitFMT"\n", st, chainStrPos[sib]);
          err++;
        }

        sib++;
      }
    }

    if (sib != sie) {
      fprintf(stderr, "iterated length wrong; sib="u32bitFMT" sie="u32bitFMT"\n", sib, sie);
      err++;
    }

    ST->rewind();

    sib = 0;

    fprintf(stderr, "initial test with full range (rewind)\n");

    while (ST->eof() == false) {
      u64bit   sp = ST->seqPos();
      u64bit   si = ST->seqIID();
      u64bit   st = ST->strPos();
      char     ch = ST->get();

      if (ch != 0) {
        if (ch != chainSeq[sib]) {
          fprintf(stderr, "letter wrong '%c' vs '%c'\n", ch, chainSeq[sib]);
          err++;
        }
        if ((ch != sep) && (sp != chainSeqPos[sib])) {
          fprintf(stderr, "seqPos wrong "u32bitFMT" vs "u32bitFMT"\n", sp, chainSeqPos[sib]);
          err++;
        }
        if ((ch != sep) && (si != chainSeqIID[sib])) {
          fprintf(stderr, "seqIID wrong "u32bitFMT" vs "u32bitFMT"\n", si, chainSeqIID[sib]);
          err++;
        }
        if ((ch != sep) && (st != chainStrPos[sib])) {
          fprintf(stderr, "strPos wrong "u32bitFMT" vs "u32bitFMT"\n", st, chainStrPos[sib]);
          err++;
        }

        sib++;
      }
    }

    if (sib != sie) {
      fprintf(stderr, "iterated length wrong; sib="u32bitFMT" sie="u32bitFMT"\n", sib, sie);
      err++;
    }
  }


  //  Set the range to random values, and check all the results.
  //  We've already verified the index works, so we're free to use
  //  that (but we currently don't).

  u32bit  maxLen = ST->startOf(numSeq-1) + ST->lengthOf(numSeq-1);

  fprintf(stderr, "test on subranges\n");

  for (u32bit iter=0; iter<100; iter++) {
    u32bit beg = mtRandom32(mtctx) % maxLen;
    u32bit end = mtRandom32(mtctx) % maxLen;
    if (beg > end) {
      u32bit t = end;
      end = beg;
      beg = t;
    }

    ST->setRange(beg, end);

    //  Compute the position in our stream for the ACGT based beg and
    //  end.  The quirk here is that our stream includes the
    //  separator.

    u32bit sib = 0;  //  chainSeq position
    u32bit sie = 0;

    for (u32bit ppp=0, sid=0; sid<numSeq; sid++) {
      u32bit len = correctSequence[sid].sequenceLength;

      if ((ppp < beg) && (beg < ppp + len)) {
        sib += beg - ppp;
        break;
      }

      ppp += len;
      sib += len + sepLen;
    }

    for (u32bit ppp=0, sid=0; sid<numSeq; sid++) {
      u32bit len = correctSequence[sid].sequenceLength;

      if ((ppp < end) && (end < ppp + len)) {
        sie += end - ppp;
        break;
      }

      ppp += len;
      sie += len + sepLen;
    }

    //  Optionally do a rewind in the middle

    if (iter % 2) {
      //fprintf(stderr, "Random iter "u32bitFMT" (with rewind)\n", iter);
      while (ST->eof() == false)
        ST->get();
      ST->rewind();
    } else {
      //fprintf(stderr, "Random iter "u32bitFMT"\n", iter);
    }

    //  Run the seqStream until it stops, checking all the interesting
    //  values, and counting the number of non-separator we see.

    while (ST->eof() == false) {
      u64bit   sp = ST->seqPos();
      u64bit   si = ST->seqIID();
      u64bit   st = ST->strPos();
      char     ch = ST->get();

      if (ch != 0) {
        if (ch != chainSeq[sib]) {
          fprintf(stderr, "letter wrong '%c' vs '%c'\n", ch, chainSeq[sib]);
          err++;
        }
        if ((ch != sep) && (sp != chainSeqPos[sib])) {
          fprintf(stderr, "seqPos wrong "u32bitFMT" vs "u32bitFMT"\n", sp, chainSeqPos[sib]);
          err++;
        }
        if ((ch != sep) && (si != chainSeqIID[sib])) {
          fprintf(stderr, "seqIID wrong "u32bitFMT" vs "u32bitFMT"\n", si, chainSeqIID[sib]);
          err++;
        }
        if ((ch != sep) && (st != chainStrPos[sib])) {
          fprintf(stderr, "strPos wrong "u32bitFMT" vs "u32bitFMT"\n", st, chainStrPos[sib]);
          err++;
        }

        sib++;
      }
    }

    if (sib != sie) {
      fprintf(stderr, "iterated length wrong; sib="u32bitFMT" sie="u32bitFMT" beg="u32bitFMT" end="u32bitFMT"\n", sib, sie, beg, end);
      err++;
    }
  }

  return(err > 0);
}



int
main(int argc, char **argv) {
  u32bit     minLen = 100;
  u32bit     maxLen = 20000;
  u32bit     numSeq = 1000;
  u32bit     err    = 0;

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

  exit(err > 0);
}

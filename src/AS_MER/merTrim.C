
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2010, J. Craig Venter Institute
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id: merTrim.C,v 1.15 2011-08-22 16:44:19 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>

#include "AS_global.h"
#include "AS_UTL_reverseComplement.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlapStore.h"

#include "AS_MER_gkpStore_to_FastABase.H"

#include "bio++.H"
#include "sweatShop.H"
#include "existDB.H"
#include "positionDB.H"
#include "libmeryl.H"

#include "merTrimResult.H"

uint32  VERBOSE = 0;

//  Correction parameters.  Bases with less than MIN_CORRECT evidence are subject to correction.
//  They will be corrected to something with at least MIN_VERIFIED evidence, or left alone if
//  the correction is weaker.

#define MIN_CORRECT         6
#define MIN_VERIFIED        2

#define ALLGOOD 1
#define ALLCRAP 2
#define ATTEMPTCORRECTION 3

//  Doesn't work, gets different results.
#undef USE_MERSTREAM_REBUILD

#undef TEST_TESTBASE

class mertrimGlobalData {
public:
  mertrimGlobalData() {
    gkpPath         = 0L;
    merCountsFile   = 0L;
    merSize         = 22;
    compression     = 0;
    numThreads      = 4;
    beVerbose       = false;
    forceCorrection = false;

    gkRead          = NULL;
    edb             = NULL;

    resFile         = NULL;

    tBgn = 0;
    tEnd = 0;
    tCur = 0;
  };

  ~mertrimGlobalData() {
    delete edb;
    delete gkRead;
    if (resFile != NULL)
      fclose(resFile);
  };

  void              initialize(void) {

    fprintf(stderr, "opening gkStore '%s'\n", gkpPath);
    gkRead  = new gkStore(gkpPath, FALSE, FALSE);

    if (tBgn == 0) {
      tBgn = 1;
      tEnd = gkRead->gkStore_getNumFragments();
    }

    tCur = tBgn;

    if (tBgn > tEnd)
      fprintf(stderr, "ERROR: invalid range:  -b ("F_U32") >= -e ("F_U32").\n",
              tBgn, tEnd), exit(1);
    if (tEnd > gkRead->gkStore_getNumFragments())
      fprintf(stderr, "ERROR: invalid range:  -e ("F_U32") > num frags ("F_U32").\n",
              tEnd, gkRead->gkStore_getNumFragments()), exit(1);

    fprintf(stderr, "loading mer database.\n");
    edb    = new existDB(merCountsFile, merSize, existDBcounts, MIN_VERIFIED, ~0);
  };

public:

  //  Command line parameters
  //
  char         *gkpPath;

  char         *merCountsFile;

  uint32        merSize;
  uint32        compression;
  uint32        numThreads;
  bool          beVerbose;
  bool          forceCorrection;

  FILE         *resFile;

  //  Global data
  //
  gkStore      *gkRead;
  existDB      *edb;

  //  Input State
  //
  uint32        tBgn;
  uint32        tCur;
  uint32        tEnd;
};




class mertrimThreadData {
public:
  mertrimThreadData(mertrimGlobalData *g) {
    kb     = new kMerBuilder(g->merSize, g->compression, 0L);
  };
  ~mertrimThreadData() {
    delete kb;
  };

public:
  kMerBuilder  *kb;
};




class mertrimComputation {
public:
  mertrimComputation() {
    origSeq    = NULL;
    origQlt    = NULL;
    corrSeq    = NULL;
    corrQlt    = NULL;
    seqMap     = NULL;
    ms         = NULL;

    disconnect = NULL;
    deletion   = NULL;
    coverage   = NULL;
    corrected  = NULL;
  }
  ~mertrimComputation() {
    delete [] origSeq;
    delete [] origQlt;
    delete [] corrSeq;
    delete [] corrQlt;
    delete [] seqMap;
    delete    ms;

    delete [] disconnect;
    delete [] deletion;
    delete [] coverage;
    delete [] corrected;
  }

  void   initialize(mertrimGlobalData *g_) {
    g  = g_;

    readIID = fr.gkFragment_getReadIID();
    seqLen  = fr.gkFragment_getSequenceLength();

    origSeq = new char   [AS_READ_MAX_NORMAL_LEN];
    origQlt = new char   [AS_READ_MAX_NORMAL_LEN];
    corrSeq = new char   [AS_READ_MAX_NORMAL_LEN];
    corrQlt = new char   [AS_READ_MAX_NORMAL_LEN];

    seqMap  = new uint32 [AS_READ_MAX_NORMAL_LEN];

    nMersExpected = 0;
    nMersTested   = 0;
    nMersFound    = 0;

    ms = NULL;

    disconnect = NULL;
    deletion   = NULL;
    coverage   = NULL;

    corrected  = NULL;

    strcpy(origSeq, fr.gkFragment_getSequence());
    strcpy(origQlt, fr.gkFragment_getQuality());
    strcpy(corrSeq, fr.gkFragment_getSequence());
    strcpy(corrQlt, fr.gkFragment_getQuality());

    //  Replace Ns with a random low-quality base.  This is necessary, since the mer routines
    //  will not make a mer for N, and we never see it to correct it.

    char  letters[4] = { 'A', 'C', 'G', 'T' };

    //  Not really replacing N with random ACGT, but good enough for us.
    for (uint32 i=0; i<seqLen; i++)
      if (corrSeq[i] == 'N') {
        corrSeq[i] = letters[i & 0x03];
        corrQlt[i] = '0';
      }

    clrBgn = 0;
    clrEnd = seqLen;

    suspectedChimer    = false;
    suspectedChimerBgn = 0;
    suspectedChimerEnd = 0;

    for (uint32 i=0; i<AS_READ_MAX_NORMAL_LEN; i++)
      seqMap[i] = i;
  };

  uint32     evaluate(bool postCorrection);

  void       reverse(void);
  void       analyze(void);

  uint32     testBases(char *bases, uint32 basesLen);
  uint32     testBaseChange(uint32 pos, char replacement);
  uint32     testBaseIndel(uint32 pos, char replacement);

  uint32     attemptCorrection(bool isReversed);

  uint32     BADBASE(uint32 i);
  uint32     attemptTrimming(void);

  void       analyzeChimer(void);

  uint32     getNumCorrected(void) {
    uint32 nCorr = 0;

    if (corrected)
      for (uint32 pos=clrBgn; pos<clrEnd; pos++)
        if ((corrected[pos] == 'C') || (corrected[pos] == 'D') || (corrected[pos] == 'I'))
          nCorr++;

    return(nCorr);
  };

  uint32     getClrBgn(void) { return(seqMap[clrBgn]); };
  uint32     getClrEnd(void) { return(seqMap[clrEnd]); };
  uint32     getSeqLen(void) { return(seqMap[seqLen]); };

  void       dump(FILE *F, char *label);

  //  Public for the writer.
  gkFragment           fr;

  mertrimGlobalData   *g;
  mertrimThreadData   *t;

  AS_IID     readIID;
  uint32     seqLen;

  char      *origSeq;
  char      *origQlt;
  char      *corrSeq;
  char      *corrQlt;

  uint32    *seqMap;

  uint32     nMersExpected;
  uint32     nMersTested;
  uint32     nMersFound;

  merStream *ms;

  uint32     clrBgn;
  uint32     clrEnd;

  bool       suspectedChimer;
  uint32     suspectedChimerBgn;
  uint32     suspectedChimerEnd;

  uint32    *disconnect;  //  per base - a hole before this base
  uint32    *deletion;    //  
  uint32    *coverage;    //  per base - mer coverage

  //  per base: C - this base corrected
  //            A - this base failed to be corrected, no answer
  //            X - this base failed to be corrected, conflicts
  uint32    *corrected;

  uint32     nHole;  //  Number of spaces (between bases) with no mer coverage
  uint32     nCorr;  //  Number of bases corrected
  uint32     nFail;  //  Number of bases uncorrected because no answer found
  uint32     nConf;  //  Number of bases uncorrected because multiple answers found

  char       merstring[256];
};




//  Scan the sequence, counting the number of kmers verified.  If we find all of them, we're done.
//  
uint32
mertrimComputation::evaluate(bool postCorrection) {

  if (ms == NULL)
    ms = new merStream(t->kb, new seqStream(corrSeq, seqLen), false, true);

  ms->rewind();

  nMersExpected = clrEnd - clrBgn - g->merSize + 1;
  nMersTested   = 0;
  nMersFound    = 0;

  while ((ms->nextMer()) &&
         (ms->thePositionInSequence() + g->merSize - 1 < clrEnd)) {
    if (ms->thePositionInSequence() < clrBgn)
      //  Mer before the clear range begins.
      continue;

    if (clrEnd <= ms->thePositionInSequence() + g->merSize - 1)
      //  Mer after the clear range ends
      continue;

    nMersTested++;

    if (g->edb->exists(ms->theCMer()))
      //  kmer exists in the database, assumed to be at least MIN_VERIFIED
      nMersFound++;
  }

  if (nMersFound == nMersExpected)
    //  All mers confirmed, read is 100% verified!
    return(ALLGOOD);

  if (nMersFound == 0)
    //  No kMers confirmed, read is 100% garbage (or 100% unique).
    return(ALLCRAP);

  //  Attempt correction.
  return(ATTEMPTCORRECTION);
}



void
mertrimComputation::reverse(void) {
  uint32  c = 0;
  uint32 *s = NULL;
  uint32 *S = NULL;

  reverseComplement(corrSeq, corrQlt, seqLen);

  delete ms;
  ms = new merStream(t->kb, new seqStream(corrSeq, seqLen), false, true);

  if (corrected) {
    s = corrected;
    S = corrected + seqLen - 1;

    while (s < S) {
      if (*s == 'X')  *s = 0;
      if (*S == 'X')  *S = 0;

      c    = *s;
      *s++ =  *S;
      *S-- =  c;
    }
  }

  s = seqMap;
  S = seqMap + seqLen - 1;

  while (s < S) {
    c    = *s;
    *s++ =  *S;
    *S-- =  c;
  }
}


//  Analyze the read.  Annotate with the depth of mer coverage.  This lets us identify
//  mismatches/insertions (though we cannot distinguish thses) and deletions.
//
void
mertrimComputation::analyze(void) {

  if (coverage == 0L)
    coverage   = new uint32 [AS_READ_MAX_NORMAL_LEN + 1];

  if (disconnect == 0L)
    disconnect = new uint32 [AS_READ_MAX_NORMAL_LEN + 1];

  if (corrected == 0L) {
    corrected  = new uint32 [AS_READ_MAX_NORMAL_LEN + 1];
    memset(corrected, 0, sizeof(uint32) * (AS_READ_MAX_NORMAL_LEN + 1));
  }

  memset(coverage,   0, sizeof(uint32) * (AS_READ_MAX_NORMAL_LEN + 1));
  memset(disconnect, 0, sizeof(uint32) * (AS_READ_MAX_NORMAL_LEN + 1));

  assert(ms != NULL);

  ms->rewind();

  while (ms->nextMer()) {
    u32bit  posBgn = ms->thePositionInSequence();
    u32bit  posEnd = ms->thePositionInSequence() + g->merSize;

    //fprintf(stderr, "posBgn %u %d\n", posBgn, posEnd);
 
    assert(posEnd <= seqLen);

    if (g->edb->exists(ms->theCMer()) == false)
      //  This mer is too weak for us.  SKip it.
      continue;

    //  If we aren't the first mer, then there should be coverage for our first base.  If not,
    //  we have found a correctable error, an uncorrectable error, or a chimeric read.
    if ((posBgn > 0) &&
        (coverage[posBgn-1] > 0) && (coverage[posBgn] == 0))
      disconnect[posBgn-1] = disconnect[posBgn] = 'D';

    //  Add coverage for the good mer.
    for (u32bit add=posBgn; add<posEnd; add++)
      coverage[add]++;

  }  //  Over all mers

  ms->rewind();

  if (VERBOSE > 1) {
    dump(stderr, "ANALYZE");
  }  //  VERBOSE
}



uint32
mertrimComputation::attemptCorrection(bool isReversed) {

  assert(coverage);
  assert(disconnect);
  assert(corrected);

  assert(ms != NULL);

  ms->rewind();

  while (ms->nextMer()) {
    uint32 pos = ms->thePositionInSequence() + g->merSize - 1;

    //  Sometimes, we see long bogus 454 reads.  Gatekeeper will trim these to 2043 bases, leaving us
    //  no space to insert bases.  If this happens, well, there isn't a whole lot we can really do,
    //  and we'll just return ALLCRAP.

    if (seqLen >= AS_READ_MAX_NORMAL_LEN)
      return(ALLCRAP);
    assert(seqLen < AS_READ_MAX_NORMAL_LEN);

    uint32  count = g->edb->count(ms->theCMer());

    //fprintf(stderr, "MER at %d is %s has count %d %s\n",
    //        pos,
    //        ms->theFMer().merToString(merstring),
    //        (count >= MIN_CORRECT) ? "CORRECT" : "ERROR",
    //        count);

    if (count >= MIN_CORRECT)
      //  Mer exists, no need to correct.
      continue;

    //  These asserts are no longer true.  The mer can exist in the table,
    //  but just at below MIN_CORRECT count.
    //
    //assert(g->edb->exists(ms->theFMer()) == false);
    //assert(g->edb->exists(ms->theRMer()) == false);

    //  State the minimum number of mers we'd accept as evidence any change we make is correct.  The
    //  penalty for OVER correcting (too low a threshold) is potentially severe -- we could
    //  insert/delete over and over and over eventually blowing up.  The penalty for UNDER
    //  correcting, however, is that we trim a read too aggressively.
    //
    //  Being strictly greater than before works for mismatches and deletions.
    //
    //  For insertions, especially insertions in single nucleotide runs, this doesn't work so well.
    //  There are two cases, an insertion before a run, and an insertion after a run.  Below, X
    //  represents ACGT, and the A's are the run.
    //    XXXXXAAAAA: inserting an A after the X's, returns the same sequence, and the count of
    //                good mers doesn't change.
    //    AAAAAXXXXX: inserting an A before the X's changes the sequence, returning AAAAAAXXXX,
    //                and it is possible (likely) that this new sequence will have a higher count.
    //
    //  We therefore require that we find at least TWO more good mers before accepting a change.
    //
    //  The drawback of this is that we cannot correct two adjacent errors.  The first error
    //  (the one we're currently working on) is corrected and adds one to the coverage count,
    //  but then we hit that second error and do not find any more mers.

    //  A solution would be to retry any base we cannot correct and allow a positive change of
    //  one mer to accept the change.  (in other words, change +1 below to +0).

    uint32 mCur = testBaseChange(pos, corrSeq[pos]) + 1;
    uint32 mNum = mCur;

    //  This is the number of replacements we found.  If zero, we mark this mer as not confirmed
    //  ('Z').  If more than one, we mark this mer as conflicting ('X').
    //
    uint32 nR = 0;

    //  Test if we can repair the sequence with a single base change.
    if (1) {
      uint32 nA = (corrSeq[pos] != 'A') ? testBaseChange(pos, 'A') : 0;
      uint32 nC = (corrSeq[pos] != 'C') ? testBaseChange(pos, 'C') : 0;
      uint32 nG = (corrSeq[pos] != 'G') ? testBaseChange(pos, 'G') : 0;
      uint32 nT = (corrSeq[pos] != 'T') ? testBaseChange(pos, 'T') : 0;
      uint32 rB = 0;

      if (VERBOSE > 2) {
        if (nA > mNum)  fprintf(stderr, "testA at %d -- %d req=%d\n", pos, nA, mNum);
        if (nC > mNum)  fprintf(stderr, "testC at %d -- %d req=%d\n", pos, nC, mNum);
        if (nG > mNum)  fprintf(stderr, "testG at %d -- %d req=%d\n", pos, nG, mNum);
        if (nT > mNum)  fprintf(stderr, "testT at %d -- %d req=%d\n", pos, nT, mNum);
      }  //  VERBOSE

      //  If we found a single perfectly correct choice, ignore all the other solutions.

      nR = 0;

      if (nA == g->merSize)  nR++;
      if (nC == g->merSize)  nR++;
      if (nG == g->merSize)  nR++;
      if (nT == g->merSize)  nR++;

      if (nR == 1) {
        if (nA != g->merSize)  nA = 0;
        if (nC != g->merSize)  nC = 0;
        if (nG != g->merSize)  nG = 0;
        if (nT != g->merSize)  nT = 0;
      }

      //  Count the number of viable solutions.

      nR = 0;

      if (nA > mNum)  { nR++;  rB = 'A'; }
      if (nC > mNum)  { nR++;  rB = 'C'; }
      if (nG > mNum)  { nR++;  rB = 'G'; }
      if (nT > mNum)  { nR++;  rB = 'T'; }

      //  If we found a single choice, correct it.

      if (nR == 1) {
        if (VERBOSE > 0) {
          fprintf(stderr, "Correct read %d at position %d from %c to %c (QV %d) (%s)\n",
                  readIID,
                  (isReversed == false) ? pos : seqLen - pos,
                  corrSeq[pos],
                  rB,
                  corrQlt[pos],
                  (isReversed == false) ? "fwd" : "rev");
        }  //  VERBOSE

        corrSeq[pos] = rB;

        corrected[pos] = 'C';

        //  Rebuild the merStream to use the corrected sequence, then move to the same position in
        //  the sequence.  This is done because we cannot simply change the string -- we need to
        //  change the state of the kMerBuilder associated with the merStream, and we can't do that.

#ifdef USE_MERSTREAM_REBUILD
        ms->rebuild();
#else
        pos = ms->thePositionInSequence();

        delete ms;
        ms = new merStream(t->kb, new seqStream(corrSeq, seqLen), false, true);
        ms->nextMer();

        while (pos != ms->thePositionInSequence())
          ms->nextMer();
#endif

        continue;

      } else if (nR > 1) {
        corrected[pos] = 'X';
      }
    }  //  End of mismatch change test

    //  Since we couldn't correct the single-base mismatch, see if it is an insertion or a deletion
    //  in the read.

    nR = 0;

    if (1) {
      uint32 nD = testBaseIndel(pos, '-');
      uint32 nA = testBaseIndel(pos, 'A');
      uint32 nC = testBaseIndel(pos, 'C');
      uint32 nG = testBaseIndel(pos, 'G');
      uint32 nT = testBaseIndel(pos, 'T');
      char   rB = 0;

      if (nD > mNum)  { nR++;  rB = '-'; }
      if (nA > mNum)  { nR++;  rB = 'A'; }
      if (nC > mNum)  { nR++;  rB = 'C'; }
      if (nG > mNum)  { nR++;  rB = 'G'; }
      if (nT > mNum)  { nR++;  rB = 'T'; }

      if (VERBOSE > 2) {
        if (nD > mNum)  fprintf(stderr, "test-- %d -- %d req=%d\n", pos, nD, mNum);
        if (nA > mNum)  fprintf(stderr, "test+A %d -- %d req=%d\n", pos, nA, mNum);
        if (nC > mNum)  fprintf(stderr, "test+C %d -- %d req=%d\n", pos, nC, mNum);
        if (nG > mNum)  fprintf(stderr, "test+G %d -- %d req=%d\n", pos, nG, mNum);
        if (nT > mNum)  fprintf(stderr, "test+T %d -- %d req=%d\n", pos, nT, mNum);
      }  //  VERBOSE

      if (nR == 0) {
        //corrected[pos] = 'Z';  (don't mark this, use coverage[] instead
        continue;
      }

      if (nR > 1) {
        corrected[pos] = 'X';
        continue;
      }

      if (nD > mNum) {
        if (VERBOSE > 0) {
          fprintf(stderr, "Correct read %d at position %d from %c to DELETE (QV %d) (%s)\n",
                  readIID,
                  (isReversed == false) ? pos : seqLen - pos,
                  corrSeq[pos],
                  corrQlt[pos] - '0',
                  (isReversed == false) ? "fwd" : "rev");

        }  //  VERBOSE
        for (uint32 i=pos; i<seqLen; i++) {
          corrSeq[i] = corrSeq[i+1];
          corrQlt[i] = corrQlt[i+1];
          seqMap[i]  = seqMap[i+1];
        }

        seqLen--;
        clrEnd--;

        corrected[pos] = 'D';

      } else {
        if (VERBOSE > 0) {
          fprintf(stderr, "Correct read %d at position %d INSERT %c (%s)\n",
                  readIID,
                  (isReversed == false) ? pos : seqLen - pos,
                  rB,
                  (isReversed == false) ? "fwd" : "rev");

        }  //  VERBOSE
        for (uint32 i=seqLen+1; i>pos; i--) {
          corrSeq[i] = corrSeq[i-1];
          corrQlt[i] = corrQlt[i-1];
          seqMap[i]  = seqMap[i-1];
        }

        corrSeq[pos] = rB;
        corrQlt[pos] = '5';
        seqMap[pos]  = seqMap[pos-1];

        seqLen++;
        clrEnd++;

        corrected[pos] = 'I';
      }

      //  Rebuild the merstream.  When we call analyze() the stream gets rewound.  If we do not
      //  restore our position, it is possible to get stuck in an infinite loop, inserting and
      //  deleting the same base over and over.  This happens because we don't explicitly require
      //  that the mer we are at be found, just that we find enough mers to make the change.
      //
      //  So, on the next pass through, we'd encounter the same mer we didn't find before, attempt
      //  to change it again, and possibly delete the base we inserted.
      //
      //  test+A 57 -- 6 req=2
      //  Correct read 6 at position 57 INSERT A
      //
      //  test-- 58 -- 3 req=2
      //  Correct read 6 at position 58 from A to DELETE (QV 5)
      //
      //  The first time through, we insert an A (with 6 mers agreeing).  The second time through,
      //  since we didn't fix the mer we were at, our choice is to delete the base we inserted (3
      //  mers tell us to do so).

#ifdef USE_MERSTREAM_REBUILD
      ms->rebuild();
#else
      pos = ms->thePositionInSequence();

      delete ms;
      ms = new merStream(t->kb, new seqStream(corrSeq, seqLen), false, true);

      analyze();

      ms->nextMer();
      while (pos != ms->thePositionInSequence())
        ms->nextMer();
#endif
    }  //  End of indel change test

  }  //  Over all mers

  if (VERBOSE > 1) {
    dump(stderr, "POSTCORRECT");
  }  //  VERBOSE

  analyze();

  return(evaluate(true));
}



uint32
mertrimComputation::testBases(char *bases, uint32 basesLen) {
  uint32  offset       = 0;
  uint32  numConfirmed = 0;

  //
  //  UNTESTED with KMER_WORDS != 1
  //

  kMer F(g->merSize);
  kMer R(g->merSize);

  for (uint32 i=1; i<g->merSize && offset<basesLen; i++, offset++) {
    F += letterToBits[bases[offset]];
    R -= letterToBits[complementSymbol[bases[offset]]];
  }

  for (uint32 i=0; i<g->merSize && offset<basesLen; i++, offset++) {
    F += letterToBits[bases[offset]];
    R -= letterToBits[complementSymbol[bases[offset]]];

    F.mask(true);
    R.mask(false);

    if (F < R) {
      if (g->edb->exists(F))
        numConfirmed++;
    } else {
      if (g->edb->exists(R))
        numConfirmed++;
    }
  }

  return(numConfirmed);
}



//  Attempt to change the base at pos to make the kmers spanning it agree.
//  Returns the number of kmers validated, and the letter to change to.
//
uint32
mertrimComputation::testBaseChange(uint32 pos, char replacement) {
  uint32   numConfirmed = 0;
  char     originalBase = corrSeq[pos];
  uint32   offset       = pos + 1 - g->merSize;

  corrSeq[pos] = replacement;

  numConfirmed = testBases(corrSeq + offset, MIN(seqLen - offset, 2 * g->merSize - 1));

#ifdef TEST_TESTBASE
  {
    uint32 oldConfirmed = 0;

    merStream  *localms = new merStream(new kMerBuilder(merSize, compression, 0L),
                                        new seqStream(corrSeq + offset, seqLen - offset),
                                        true,
                                        true);

    //  Test
    for (uint32 i=0; i<g->merSize && localms->nextMer(); i++)
      if (g->edb->exists(localms->theCMer()))
        oldConfirmed++;

    delete localms;

    assert(oldConfirmed == numConfirmed);
  }
#endif

  corrSeq[pos] = originalBase;

  //if (numConfirmed > 0)
  //  fprintf(stderr, "testBaseChange() pos=%d replacement=%c confirmed=%d\n",
  //          pos, replacement, numConfirmed);

  return(numConfirmed);
}




uint32
mertrimComputation::testBaseIndel(uint32 pos, char replacement) {
  uint32   numConfirmed = 0;
  char     testStr[128] = {0};
  uint32   len          = 0;
  uint32   offset       = pos + 1 - g->merSize;
  uint32   limit        = g->merSize * 2 - 1;

  assert(2 * g->merSize < 120);  //  Overly pessimistic

  //  Copy the first merSize bases.

  while (len < g->merSize - 1)
    testStr[len++] = corrSeq[offset++];

  //  Copy the second merSize bases, but overwrite the last base in the first copy (if we're testing
  //  a de;etion) or insert a replacement base.

  if (replacement == '-') {
    offset++;
  } else {
    testStr[len++] = replacement;
  }

  //  Copy the rest of the bases.

  while ((len < limit) && (corrSeq[offset]))
    testStr[len++] = corrSeq[offset++];

  numConfirmed = testBases(testStr, len);

#ifdef TEST_TESTBASE
  {
    uint32 oldConfirmed = 0;

    merStream  *localms = new merStream(new kMerBuilder(merSize, compression, 0L),
                                        new seqStream(testStr, len),
                                        true,
                                        true);

    //  Test
    for (uint32 i=0; i<g->merSize && localms->nextMer(); i++)
      if (g->edb->exists(localms->theCMer()))
        oldConfirmed++;

    delete localms;

    assert(oldConfirmed == numConfirmed);
  }
#endif

  //if (numConfirmed > 0)
  //  fprintf(stderr, "testBaseIndel() pos=%d replacement=%c confirmed=%d\n",
  //          pos, replacement, numConfirmed);

  return(numConfirmed);
}




//  A BADBASE is one that was:
//
//  o Corrected by an base change, insert or delete.
//  o Has very very low quality, at or below QV 6 == 25% chance of error.
//  o Has lower quality (QV of 16 == 2.5% chance of error), and is marked as being corrected.
//    The marking here should be either "mer wasn't found" or "correction was conflicting", but we
//    don't actually care what the mark says.
//
//  We don't want to trust many of these bases close together.
//
uint32
mertrimComputation::BADBASE(uint32 i) {

  if ((corrected) && ((corrected[i] == 'C') || (corrected[i] == 'I') || (corrected[i] == 'D')))
    //  We made a correction.
    return(1);

#if 0
  if (corrQlt[i] - '0' <= 6)
    //  Low quality base - too aggressive?
    return(1);
#endif

  if ((coverage) && ((coverage[i] == 'A') && (corrQlt[i] - '0' < 16)))
    //  Lower quality base, and we couldn't correct it or find it in the mers.
    return(1);

  return(0);
}

uint32
mertrimComputation::attemptTrimming(void) {

  //  Any 5-base window with more two or more errors or low quality.  Trim to the inner-most.

  if (corrected == NULL) {
    corrected  = new uint32 [AS_READ_MAX_NORMAL_LEN + 1];
    memset(corrected,  0, sizeof(uint32) * (AS_READ_MAX_NORMAL_LEN + 1));
  }

  //  Compute the number of errors in each 5-base window.
  //
  //  In general, wErr[i] = B[i-2] + B[i-1] + B[i] + B[i+1] + B[i+2].
  //
  uint32  wErr[AS_READ_MAX_NORMAL_LEN] = {0};

  wErr[0] = BADBASE(0) + BADBASE(1) + BADBASE(2);
  wErr[1] = wErr[0] + BADBASE(3);
  wErr[2] = wErr[1] + BADBASE(4);

  for (uint32 i=3; i<seqLen-2; i++)
    wErr[i] = wErr[i-1] - BADBASE(i-3) + BADBASE(i+2);

  wErr[seqLen-2] = wErr[seqLen-3] - BADBASE(seqLen-4);
  wErr[seqLen-1] = wErr[seqLen-2] - BADBASE(seqLen-3);

#if 0
  fprintf(stderr, "B=");
  for (uint32 i=0; i<seqLen; i++) {
    fprintf(stderr, "%d", BADBASE(i));
  }
  fprintf(stderr, "\nW=");
  for (uint32 i=0; i<seqLen; i++) {
    fprintf(stderr, "%d", wErr[i]);
  }
  fprintf(stderr, "\n");
#endif

  //  Find the largest region with < 2 errors.  Whenever we hit a region with
  //  two errors, clear the old region (saving it if it was the biggest).

  uint32  lBgn=0, lEnd=0;  //  Largest found so far
  uint32  tBgn=2, tEnd=2;  //  The one we're looking at

  for (uint32 i=2; i<seqLen-2; i++) {

    if (wErr[i] >= 2) {
      //  Lots of errors in this window, save any good region found so far.
      if (tEnd - tBgn > lEnd - lBgn) {
        lBgn = tBgn;
        lEnd = tEnd;
      }
      //  Any good region must start at least at the next position.
      tBgn = i+1;
      tEnd = i+1;

    } else {
      //  Good window, extend any existing region.
      tEnd = i;
    }
  }

  //  Save any remaining window.
  if (tEnd - tBgn > lBgn - lEnd) {
    lBgn = tBgn;
    lEnd = tEnd;
  }

  //  The (pre) clear range is jsut the largest window found, but trimmed back (into potentially
  //  good sequence) to the end/start of the window.
  //  limit to valid sequences
  clrBgn = (lBgn <= seqLen - 2 ? lBgn + 2 : seqLen);
  clrEnd = (lEnd >= 2 ? lEnd - 2 : 0);

  if (VERBOSE > 2) {
    fprintf(stderr, "TRIM: %d,%d (pre)\n", clrBgn, clrEnd);
  }  //  VERBOSE

  //  Adjust borders to the first bad base.  This undoes the (pre) window trimming.  Remember,
  //  clrBgn includes the base (and so we must conditionally move ahead one), clrEnd does not.

  while ((clrBgn > 0) && (BADBASE(clrBgn) == 0))
    clrBgn--;
  if (BADBASE(clrBgn))
    clrBgn++;

  while ((clrEnd < seqLen) && (BADBASE(clrEnd) == 0))
    clrEnd++;

  if (clrBgn >= clrEnd) {
    clrBgn = 0;
    clrEnd = 0;
  }

  if (VERBOSE > 1) {
    fprintf(stderr, "TRIM: %d,%d (post)\n", clrBgn, clrEnd);
    dump(stderr, "TRIM");
  }  //  VERBOSE

  return(evaluate(true));
}



void
mertrimComputation::analyzeChimer(void) {

  //  Examine the coverage for a specific pattern that indicates a chimeric read (or a read with an
  //  uncorrected error):  a valley in the coverage:  ..../---\./---\...

  assert(suspectedChimer    == false);
  assert(suspectedChimerBgn == 0);
  assert(suspectedChimerEnd == 0);

  if (coverage == NULL)
    //  No coverage?  Must have been a perfect read.
    return;

  int32   floc = 0, rloc = seqLen-1;
  int32   fcov = coverage[floc];
  int32   rcov = coverage[rloc];

  //  Search in from the ends while the coverage monotonically increases.

  while ((floc < seqLen) && (fcov <= coverage[floc])) {
    fcov = coverage[floc];
    floc++;
  }

  while ((rloc > 0) && (rcov <= coverage[rloc])) {
    rcov = coverage[rloc];
    rloc--;
  }

  if (rloc <= floc)
    //  If ranges are flipped, we can stop.  No chimer found.
    return;

  //  Otherwise, there is a dip in the coverage starting at floc until rloc.  Continue searching in
  //  while the coverage monotonically drops.

  while ((floc < seqLen) && (fcov >= coverage[floc])) {
    fcov = coverage[floc];
    floc++;
  }

  while ((rloc > 0) && (rcov >= coverage[rloc])) {
    rcov = coverage[rloc];
    rloc--;
  }

  //  The searches always go one too far.  This backs up the point to be the last monotonically
  //  decreasing value.

  floc--;
  rloc++;

  //  If the coverages are different, or if floc < rloc, then we're not at a junction.  Something
  //  bizarre happened in this read, and there are two valleys instead of one.

  if ((floc < rloc) || (fcov != rcov))
    return;

  assert(fcov == rcov);
  assert(floc >= rloc);

  //  By chance, mers will cross the junction.  This will inflate the coverage count.
  //
  //  Ideally, there is a single base added between the genomic sequences, and so our
  //  coverage should drop to zero.  Or the chimera is from an abutment and the coverage will be one.
  //
  //  The probability that we're extending past the chimeric region by X bases is 0.25^X:
  //      0 - 0.25^0 = 100%
  //      1 - 0.25^1 =  25%
  //      2 - 0.25^2 =   6.25%
  //      3 - 0.25^3 =   1.5625%
  //      4 - 0.25^4 =   0.3906%
  //      5 - 0.25^5 =   0.0977%
  //      6 - 0.25^6 =   0.0244%
  //      7 - 0.25^7 =   0.0061%
  //      8 - 0.25^8 =   0.0015%
  //
  //  Unfortunately, we cannot reliably tell how far we've passed the junction from coverage alone.
  //  We blindly declare that 7 is too high.

  if (fcov >= 6)
    return;

  //  For a true chimeric junction, the pattern we should see is:
  //
  //    fcov == rcov == X
  //    floc == rloc + X     (X is also the number of bases the junction is crossed by)
  //
  //  In general, uncorrected errors in the read should not have this pattern; only errors near SNPs
  //  or in diverged repeats should be spuriously spanned.  The uncorrected error is, by definition,
  //  different from the real sequence -- where in the chimeric junction case, the base after the
  //  junction has a 25% chance of being the same as the true next base.
  //
#warning NOT CORRECT
  if (floc == rloc + fcov) {
    suspectedChimer    = true;
    suspectedChimerBgn = rloc;
    suspectedChimerEnd = floc;
    return;
  }

  //  If the 'loc' pattern isn't met, there must be something else going on.  Do we err on the side
  //  of caution and label this as chimeric read??
  //
  //  Examples:
  //   * two uncorrected errors next to each other.  We cannot correct these.
  //   * a pile of bases in the middle of a read with lots of low quality on the end.
  //     the bases were composed of T's and A's only.

  //fprintf(stderr, "CHIMER?  floc=%d rloc=%d  cov=%d\n",
  //        floc, rloc, fcov);

  return;
}




void
mertrimComputation::dump(FILE *F, char *label) {
  fprintf(F, "%s read %d len %d (trim %d-%d)\n", label, readIID, seqLen, clrBgn, clrEnd);
  for (uint32 i=0; origSeq[i]; i++) {
    if (i == clrBgn)
      fprintf(F, "-[");
    fprintf(F, "%c", origSeq[i]);
    if (i+1 == clrEnd)
      fprintf(F, "]-");
  }
  fprintf(F, " (ORI)\n");
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn)
      fprintf(F, "-[");
    fprintf(F, "%c", corrSeq[i]);
    if (i+1 == clrEnd)
      fprintf(F, "]-");
  }
  fprintf(F, " (SEQ)\n");
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn)
      fprintf(F, "-[");
    fprintf(F, "%c", corrQlt[i]);
    if (i+1 == clrEnd)
      fprintf(F, "]-");
  }
  fprintf(F, " (QLT)\n");
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn)
      fprintf(F, "-[");
    fprintf(F, "%c", (coverage) ? coverage[i] + 'A' : 'A');
    if (i+1 == clrEnd)
      fprintf(F, "]-");
  }
  fprintf(F, " (COVERAGE)\n");
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn)
      fprintf(F, "-[");
    fprintf(F, "%c", (corrected && corrected[i]) ? corrected[i] : ' ');
    if (i+1 == clrEnd)
      fprintf(F, "]-");
  }
  fprintf(F, " (CORRECTIONS)\n");
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn)
      fprintf(F, "-[");
    fprintf(F, "%c", (disconnect && disconnect[i]) ? disconnect[i] : ' ');
    if (i+1 == clrEnd)
      fprintf(F, "]-");
  }
  fprintf(F, " (DISCONNECTION)\n");
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn)
      fprintf(F, "-[");
    fprintf(F, "%c", (deletion && deletion[i]) ? deletion[i] : ' ');
    if (i+1 == clrEnd)
      fprintf(F, "]-");
  }
  fprintf(F, " (DELETIONS)\n");
}





void
mertrimWorker(void *G, void *T, void *S) {
  mertrimGlobalData    *g = (mertrimGlobalData  *)G;
  mertrimThreadData    *t = (mertrimThreadData  *)T;
  mertrimComputation   *s = (mertrimComputation *)S;

  if (VERBOSE) {
    fprintf(stderr, "\nPROCESS\n");
  }  //  VERBOSE

  s->t = t;

  uint32  eval = s->evaluate(false);

  //  Read is perfect!
  //
  if      (eval == ALLGOOD) {
    goto finished;
  }

  //  Read had NO mers found.  We'll let it through, but it might just end up a singleton.  There
  //  is a slight chance it has a 2-copy mer that will pull it in, or maybe the mate will.
  //
  if (eval == ALLCRAP) {
    s->attemptTrimming();
    goto finished;
  }

  s->reverse();
  s->analyze();
  eval = s->attemptCorrection(true);

  s->reverse();
  s->analyze();
  eval = s->attemptCorrection(false);

  //  Correction worked perfectly, we're done.
  //
  if ((eval == ALLGOOD) && (s->getNumCorrected() < 4)) {
    goto finished;
  }

  eval = s->attemptTrimming();

  if ((eval == ALLGOOD) && (s->getNumCorrected() < 4)) {
    goto finished;
  }


 finished:

  //  SKIPPING until the heuristics are worked out.
  //s->analyzeChimer();

  if (VERBOSE) {
    s->dump(stderr, "FINAL");
  }  //  VERBOSE
}


void *
mertrimReader(void *G) {
  mertrimGlobalData    *g = (mertrimGlobalData  *)G;
  mertrimComputation   *s = NULL;

  while ((g->tCur <= g->tEnd) &&
         (s == NULL)) {
    s = new mertrimComputation();

    g->gkRead->gkStore_getFragment(g->tCur, &s->fr, GKFRAGMENT_QLT);
    g->tCur++;

    if ((g->forceCorrection) ||
        (g->gkRead->gkStore_getLibrary(s->fr.gkFragment_getLibraryIID())->doTrim_initialMerBased)) {
      s->initialize(g);
    } else {
      delete s;
      s = NULL;
    }
  }

  return(s);
}



void
mertrimWriter(void *G, void *S) {
  mertrimGlobalData    *g = (mertrimGlobalData  *)G;
  mertrimComputation   *s = (mertrimComputation *)S;
  mertrimResult         res;

  //  The get*() functions return positions in the original uncorrected sequence.  They
  //  map from positions in the corrected sequence (which has inserts and deletes) back
  //  to the original sequence.
  //
  assert(s->getClrBgn() <= s->getClrEnd());
  assert(s->getClrBgn() <  s->getSeqLen());
  assert(s->getClrEnd() <= s->getSeqLen());

  res.readIID = s->fr.gkFragment_getReadIID();

  if ((s->getClrEnd() <= s->getClrBgn()) ||
      (s->getClrEnd() - s->getClrBgn() < AS_READ_MIN_LEN))
    res.deleted = true;
  else
    res.deleted = false;

  res.clrBgn  = s->getClrBgn();
  res.clrEnd  = s->getClrEnd();

  res.chimer  = s->suspectedChimer;
  res.chmBgn  = s->suspectedChimerBgn;
  res.chmEnd  = s->suspectedChimerEnd;

  res.writeResult(g->resFile);

  delete s;
}










int
main(int argc, char **argv) {
  mertrimGlobalData  *g        = new mertrimGlobalData;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      g->gkpPath = argv[++arg];

    } else if (strcmp(argv[arg], "-m") == 0) {
      g->merSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      g->compression = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-mc") == 0) {
      g->merCountsFile = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      g->numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-b") == 0) {
      g->tBgn = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      g->tEnd = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-v") == 0) {
      g->beVerbose = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      VERBOSE++;
 
    } else if (strcmp(argv[arg], "-f") == 0) {
      g->forceCorrection = true;
 
    } else if (strcmp(argv[arg], "-o") == 0) {
      errno = 0;
      g->resFile = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open output file '%s': %s\n", argv[arg], strerror(errno)), exit(1);

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((g->gkpPath == 0L) || (err)) {
    fprintf(stderr, "usage: %s -g gkpStore -m merSize -mc merCountsFile [-v]\n", argv[0]);
    exit(1);
  }

  gkpStoreFile::registerFile();

  g->initialize();

  gkFragment   fr;

#if 0
  //  DEBUG, non-threaded version.
  speedCounter SC(" Trimming: %11.0f reads -- %7.5f reads/second\r", 1.0, 0x1fff, true);

  mertrimThreadData *t = new mertrimThreadData(g);

  g->tBgn = 140222;
  g->tCur = 140222;
  g->tEnd = 140222;

  mertrimComputation *s = (mertrimComputation *)mertrimReader(g);
  while (s) {
    mertrimWorker(g, t, s);
    mertrimWriter(g, s);
    SC.tick();
    s = (mertrimComputation *)mertrimReader(g);
  }

  delete t;
#else
  //  PRODUCTION, threaded version
  sweatShop *ss = new sweatShop(mertrimReader, mertrimWorker, mertrimWriter);

  ss->setLoaderQueueSize(16384);
  ss->setWriterQueueSize(1024);

  ss->setNumberOfWorkers(g->numThreads);

  for (u32bit w=0; w<g->numThreads; w++)
    ss->setThreadData(w, new mertrimThreadData(g));  //  these leak

  ss->run(g, g->beVerbose);  //  true == verbose
#endif

  delete g;

  fprintf(stderr, "\nSuccess!  Bye.\n");

  return(0);
}


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

const char *mainid = "$Id: merTrim.C,v 1.19 2012-01-23 06:42:54 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_reverseComplement.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"  //  QUALITY_MAX, QV conversion
#include "AS_OVS_overlapStore.h"

#include <algorithm>

#include "AS_MER_gkpStore_to_FastABase.H"

#include "bio++.H"
#include "sweatShop.H"
#include "existDB.H"
#include "positionDB.H"
#include "libmeryl.H"

#include "merTrimResult.H"

uint32  VERBOSE = 0;

#define ALLGOOD 1
#define ALLCRAP 2
#define ATTEMPTCORRECTION 3

//  Doesn't work, gets different results.
#undef USE_MERSTREAM_REBUILD

#undef TEST_TESTBASE

class mertrimGlobalData {
public:
  mertrimGlobalData() {
    gkpPath                   = 0L;
    fqInputPath               = 0L;
    fqOutputPath              = 0L;

    merCountsFile             = 0L;

    merSize                   = 22;
    compression               = 0;
    numThreads                = 4;

    beVerbose                 = false;
    forceCorrection           = false;
    correctMismatch           = true;
    correctIndel              = true;

    actualCoverage            = 0;          //  Estimate of coverage
    minCorrectFraction        = 1.0 / 3.0;  //  Base can be corrected if less than 1/3 coverage
    minCorrect                = 0;          //
    minVerifiedFraction       = 1.0 / 4.0;  //  Base can be corrected to something with only 1/4 coverage
    minVerified               = 0;          //

    endTrimDefault            = true;
    endTrimNum                = 2;

    endTrimWinScale[0]        = 0.50;
    endTrimErrAllow[0]        = 2;

    endTrimWinScale[1]        = 0.25;
    endTrimErrAllow[0]        = 0;

    discardZeroCoverage       = false;
    discardImperfectCoverage  = false;

    gkRead                    = NULL;

    fqInput                   = NULL;
    fqOutput                  = NULL;
    fqLog                     = NULL;

    edb                       = NULL;

    resPath                   = NULL;
    resFile                   = NULL;

    gktBgn                    = 0;
    gktEnd                    = 0;
    gktCur                    = 0;
  };

  ~mertrimGlobalData() {
    delete gkRead;

    if (fqInput)
      fclose(fqInput);
    if (fqOutput)
      fclose(fqOutput);
    if (fqLog)
      fclose(fqLog);

    delete edb;

    if (resFile != NULL)
      fclose(resFile);
  };

  void              initializeGatekeeper(void) {

    if (gkpPath == NULL)
      return;

    fprintf(stderr, "opening gkStore '%s'\n", gkpPath);
    gkRead  = new gkStore(gkpPath, FALSE, FALSE);

    if (gktBgn == 0) {
      gktBgn = 1;
      gktEnd = gkRead->gkStore_getNumFragments();
    }

    gktCur = gktBgn;

    if (gktBgn > gktEnd)
      fprintf(stderr, "ERROR: invalid range:  -b ("F_U32") >= -e ("F_U32").\n",
              gktBgn, gktEnd), exit(1);
    if (gktEnd > gkRead->gkStore_getNumFragments())
      fprintf(stderr, "ERROR: invalid range:  -e ("F_U32") > num frags ("F_U32").\n",
              gktEnd, gkRead->gkStore_getNumFragments()), exit(1);

    errno = 0;
    resFile = fopen(resPath, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", resPath, strerror(errno)), exit(1);
  };

  void              initializeFASTQ(void) {

    if (fqInputPath == NULL)
      return;

    if (fqOutputPath == NULL)
      return;

    errno = 0;
    fqInput = fopen(fqInputPath, "r");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", fqInputPath, strerror(errno)), exit(1);

    errno = 0;
    fqOutput = fopen(fqOutputPath, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", fqOutputPath, strerror(errno)), exit(1);

    char fqName[FILENAME_MAX];

    sprintf(fqName, "%s.log", fqOutputPath);

    errno = 0;
    fqLog = fopen(fqName, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", fqName, strerror(errno)), exit(1);
  };

  void              initialize(void) {
    initializeGatekeeper();
    initializeFASTQ();

    if (actualCoverage == 0) {
      merylStreamReader  *MF = new merylStreamReader(merCountsFile);

      uint32  i  = 0;
      uint32  iX = 0;

      fprintf(stderr, "distinct: "u64bitFMT"\n", MF->numberOfDistinctMers());
      fprintf(stderr, "unique:   "u64bitFMT"\n", MF->numberOfUniqueMers());
      fprintf(stderr, "total:    "u64bitFMT"\n", MF->numberOfTotalMers());

      fprintf(stderr, "Xcoverage zero 1 0 "F_U64"\n", MF->histogram(1));

      for (i=2; (i < MF->histogramLength()) && (MF->histogram(i-1) > MF->histogram(i)); i++)
        fprintf(stderr, "Xcoverage drop "F_U32" "F_U64" "F_U64"\n", i, MF->histogram(i-1), MF->histogram(i));

      iX = i - 1;

      for (; i < MF->histogramLength(); i++) {
        if (MF->histogram(iX) < MF->histogram(i)) {
          fprintf(stderr, "Xcoverage incr "F_U32" "F_U64" "F_U64"\n", i, MF->histogram(iX), MF->histogram(i));
          iX = i;
        } else {
          //fprintf(stderr, "Xcoverage drop "F_U32" "F_U64" "F_U64"\n", i, MF->histogram(iX), MF->histogram(i));
        }
      }

      fprintf(stderr, "Guessed X coverage is "F_U32"\n", iX);

      delete MF;

      actualCoverage = iX;
    }

    if (minCorrectFraction > 0)
      minCorrect = (int)floor(minCorrectFraction * actualCoverage);

    if (minVerifiedFraction > 0)
      minVerified = (int)floor(minVerifiedFraction * actualCoverage);

    fprintf(stderr, "loading mer database.\n");
    edb    = new existDB(merCountsFile, merSize, existDBcounts, minVerified, ~0);
  };

public:

  //  Command line parameters
  //
  char         *gkpPath;
  char         *fqInputPath;
  char         *fqOutputPath;

  char         *merCountsFile;

  uint32        merSize;
  uint32        compression;
  uint32        numThreads;

  bool          beVerbose;
  bool          forceCorrection;
  bool          correctMismatch;
  bool          correctIndel;

  char         *resPath;
  FILE         *resFile;

  bool          endTrimDefault;
  uint32        endTrimNum;
  double        endTrimWinScale[16];
  uint32        endTrimErrAllow[16];

  bool          discardZeroCoverage;
  bool          discardImperfectCoverage;

  //  Global data
  //
  gkStore      *gkRead;

  uint32        actualCoverage;
  double        minCorrectFraction;
  uint32        minCorrect;
  double        minVerifiedFraction;
  uint32        minVerified;

  FILE         *fqInput;
  FILE         *fqOutput;
  FILE         *fqLog;

  existDB      *edb;

  //  Input State
  //
  uint32        gktBgn;
  uint32        gktCur;
  uint32        gktEnd;
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
    readName   = NULL;

    origSeq    = NULL;
    origQlt    = NULL;
    corrSeq    = NULL;
    corrQlt    = NULL;
    seqMap     = NULL;
    ms         = NULL;

    disconnect = NULL;
    coverage   = NULL;
    corrected  = NULL;
  }
  ~mertrimComputation() {
    delete [] readName;
    delete [] origSeq;
    delete [] origQlt;
    delete [] corrSeq;
    delete [] corrQlt;
    delete [] seqMap;
    delete    ms;

    delete [] disconnect;
    delete [] coverage;
    delete [] corrected;
  }


  void   initializeGatekeeper(mertrimGlobalData *g_) {
    g  = g_;

    readIID  = fr.gkFragment_getReadIID();
    seqLen   = fr.gkFragment_getSequenceLength();
    allocLen = seqLen + seqLen;

    readName = NULL;

    origSeq  = new char   [allocLen];
    origQlt  = new char   [allocLen];
    corrSeq  = new char   [allocLen];
    corrQlt  = new char   [allocLen];

    seqMap   = new uint32 [allocLen];

    nMersExpected = 0;
    nMersTested   = 0;
    nMersFound    = 0;

    ms = NULL;

    disconnect = NULL;
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

    gapInConfirmedKmers   = false;
    hasNoConfirmedKmers   = false;
    imperfectKmerCoverage = false;

    suspectedChimer       = false;
    suspectedChimerBgn    = 0;
    suspectedChimerEnd    = 0;

    for (uint32 i=0; i<allocLen; i++)
      seqMap[i] = i;
  };


  bool    initializeFASTQ(mertrimGlobalData *g_) {
    g  = g_;

    readIID  = g->gktCur++;
    seqLen   = 0;
    allocLen = AS_READ_MAX_NORMAL_LEN + 1;  //  Used for seq/qlt storage only

    readName = new char   [1024];

    origSeq  = new char   [allocLen];
    origQlt  = new char   [allocLen];
    corrSeq  = new char   [allocLen];
    corrQlt  = new char   [allocLen];

    seqMap   = new uint32 [allocLen];

    nMersExpected = 0;
    nMersTested   = 0;
    nMersFound    = 0;

    ms = NULL;

    disconnect = NULL;
    coverage   = NULL;
    corrected  = NULL;

    fgets(readName, 1024,     g->fqInput);
    fgets(origSeq,  allocLen, g->fqInput);
    fgets(origQlt,  allocLen, g->fqInput);  //  qv name line, ignored
    fgets(origQlt,  allocLen, g->fqInput);

    if (feof(g->fqInput))
      return(false);

    chomp(readName);
    chomp(origSeq);
    chomp(origQlt);

    //  Adjust QVs to CA encoding.

#warning ASSUMING SANGER QV ENCODING

    for (uint32 i=0; origQlt[i]; i++) {
      if (origQlt[i] < '!')
        origQlt[i] = '!';
      origQlt[i] -= '!';
      if (origQlt[i] > QUALITY_MAX)
        origQlt[i] = QUALITY_MAX;
      origQlt[i] += '0';
    }

    //  Copy to the corrected sequence

    strcpy(corrSeq, origSeq);
    strcpy(corrQlt, origQlt);

    //  Replace Ns with a random low-quality base.  This is necessary, since the mer routines
    //  will not make a mer for N, and we never see it to correct it.

    char  letters[4] = { 'A', 'C', 'G', 'T' };

    //  Not really replacing N with random ACGT, but good enough for us.
    for (uint32 i=0; corrSeq[i]; i++)
      if (corrSeq[i] == 'N') {
        corrSeq[i] = letters[i & 0x03];
        corrQlt[i] = '0';
      }

    seqLen   = strlen(origSeq);
    allocLen = seqLen + seqLen;

    clrBgn = 0;
    clrEnd = seqLen;

    gapInConfirmedKmers   = false;
    hasNoConfirmedKmers   = false;
    imperfectKmerCoverage = false;

    suspectedChimer       = false;
    suspectedChimerBgn    = 0;
    suspectedChimerEnd    = 0;

    for (uint32 i=0; i<allocLen; i++)
      seqMap[i] = i;

    return(true);
  };


  uint32     evaluate(void);

  void       reverse(void);
  void       analyze(void);

  uint32     testBases(char *bases, uint32 basesLen);
  uint32     testBaseChange(uint32 pos, char replacement);
  uint32     testBaseIndel(uint32 pos, char replacement);

  void       attemptCorrection(bool isReversed);

  void       attemptTrimming(void);
  void       attemptTrimming5End(uint32 *errorPos, uint32 endWindow, uint32 endAllowed);
  void       attemptTrimming3End(uint32 *errorPos, uint32 endWindow, uint32 endAllowed);

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
  uint32     allocLen;

  char      *readName;

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

  bool       gapInConfirmedKmers;
  bool       hasNoConfirmedKmers;
  bool       imperfectKmerCoverage;

  bool       suspectedChimer;
  uint32     suspectedChimerBgn;
  uint32     suspectedChimerEnd;

  uint32    *disconnect;  //  per base - a hole before this base
  uint32    *coverage;    //  per base - mer coverage

  uint32    *corrected;   //  per base - type of correction here

  uint32     nHole;  //  Number of spaces (between bases) with no mer coverage
  uint32     nCorr;  //  Number of bases corrected
  uint32     nFail;  //  Number of bases uncorrected because no answer found
  uint32     nConf;  //  Number of bases uncorrected because multiple answers found

  char       merstring[256];
};




//  Scan the sequence, counting the number of kmers verified.  If we find all of them, we're done.
//  
uint32
mertrimComputation::evaluate(void) {

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
      //  kmer exists in the database, assumed to be at least g->minVerified
      nMersFound++;
  }

  if (nMersFound == nMersExpected)
    //  All mers confirmed, read is 100% verified!
    return(ALLGOOD);

  if (nMersFound == 0) {
    //  No kMers confirmed, read is 100% garbage (or 100% unique).
    hasNoConfirmedKmers = true;
    return(ALLCRAP);
  }

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

  if (coverage == NULL)
    coverage = new uint32 [allocLen];
  if (disconnect == NULL)
    disconnect = new uint32 [allocLen];

  memset(coverage,   0, sizeof(uint32) * (allocLen));
  memset(disconnect, 0, sizeof(uint32) * (allocLen));

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

  if (VERBOSE > 1)
    dump(stderr, "ANALYZE");
}



void
mertrimComputation::attemptCorrection(bool isReversed) {

  if (corrected == NULL) {
    corrected = new uint32 [allocLen];
    memset(corrected, 0, sizeof(uint32) * (allocLen));
  }

  assert(coverage);
  assert(disconnect);
  assert(corrected);

  assert(ms != NULL);

  ms->rewind();

  while (ms->nextMer()) {
    uint32  pos   = ms->thePositionInSequence() + g->merSize - 1;
    uint32  count = g->edb->count(ms->theCMer());

    //fprintf(stderr, "MER at %d is %s has count %d %s\n",
    //        pos,
    //        ms->theFMer().merToString(merstring),
    //        (count >= g->minCorrect) ? "CORRECT" : "ERROR",
    //        count);

    if (count >= g->minCorrect)
      //  Mer exists, no need to correct.
      continue;

    //  These asserts are no longer true.  The mer can exist in the table,
    //  but just at below g->minCorrect count.
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
    if (g->correctMismatch) {
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

    if ((g->correctIndel) &&
        (g->merSize     < pos) &&
        (pos            < seqLen - g->merSize)) {
      uint32 nD = testBaseIndel(pos, '-');
      uint32 nA = testBaseIndel(pos, 'A');
      uint32 nC = testBaseIndel(pos, 'C');
      uint32 nG = testBaseIndel(pos, 'G');
      uint32 nT = testBaseIndel(pos, 'T');
      char   rB = 0;

      mNum += g->merSize / 3;

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

      if (nR == 0)
        continue;

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




void
mertrimComputation::attemptTrimming5End(uint32 *errorPos, uint32 endWindow, uint32 errAllow) {

  for (bool doTrim = true; doTrim; ) {
    uint32  endFound   = 0;
    uint32  endTrimPos = 0;

    uint32   bgn = clrBgn;
    uint32   end = clrBgn + endWindow;

    if (end > clrEnd)
      end = clrEnd;

    for (uint32 i=bgn; i<end; i++) {
      if (errorPos[i]) {
        endFound++;
        endTrimPos = i + 1;
      }
    }

    if (VERBOSE > 1)
      fprintf(stderr, "BGNTRIM found=%u pos=%u from %u to %u\n",
              endFound, endTrimPos,
              clrBgn, clrBgn + endWindow);

    if (endFound > errAllow)
      clrBgn = endTrimPos;

    else
      doTrim = false;
  }
}


void
mertrimComputation::attemptTrimming3End(uint32 *errorPos, uint32 endWindow, uint32 errAllow) {

  for (bool doTrim = true; doTrim; ) {
    uint32  endFound   = 0;
    uint32  endTrimPos = clrEnd;

    uint32   bgn = clrBgn;
    uint32   end = clrEnd;

    if (clrBgn + endWindow <= clrEnd)
      bgn = clrEnd - endWindow;

    assert(bgn >= clrBgn);
    assert(bgn <= clrEnd);

    for (int32 i=bgn; i<end; i++) {
      if (errorPos[i]) {
        endFound++;
        if (i < endTrimPos)
          endTrimPos = i;
      }
    }

    if (VERBOSE > 1)
      fprintf(stderr, "ENDTRIM found=%u pos=%u from %u to %u\n",
              endFound, endTrimPos,
              clrEnd - endWindow, clrEnd);

    if (endFound > errAllow)
      clrEnd = endTrimPos;

    else
      doTrim = false;
  }
}


void
mertrimComputation::attemptTrimming(void) {

  //  Just bail if the read is all junk.  Nothing to do here.
  //
  if (hasNoConfirmedKmers) {
    clrBgn = 0;
    clrEnd = 0;
    return;
  }

  assert(coverage != NULL);

  //  Lop off the ends with no confirmed kmers
  //
  while ((clrBgn < clrEnd) && (coverage[clrBgn] == 0))
    clrBgn++;

  while ((clrEnd > clrBgn) && (coverage[clrEnd-1] == 0))
    clrEnd--;


  //  True if there is an error at position i.
  //
  uint32  *errorPos = new uint32 [seqLen];

  for (uint32 i=0; i<seqLen; i++) {
    errorPos[i] = ((corrected[i] == 'C') || (corrected[i] == 'I') || (corrected[i] == 'D'));
  }


  //  If there are more than 'errAllow' corrections in the last 'winScale' (x kmerSize) bases of the
  //  read, trim all the errors off.

  for (uint32 i=0; i<g->endTrimNum; i++) {
    //fprintf(stderr, "endTrim: %f %u\n", g->endTrimWinScale[i], g->endTrimErrAllow[i]);
    attemptTrimming5End(errorPos, g->merSize * g->endTrimWinScale[i], g->endTrimErrAllow[i]);
    attemptTrimming3End(errorPos, g->merSize * g->endTrimWinScale[i], g->endTrimErrAllow[i]);
  }

  delete [] errorPos;

  //  If there are zero coverage areas in the interior, trash the whole read.

  if (g->discardZeroCoverage == true) {
    for (uint32 i=clrBgn; i<clrEnd; i++)
      if (coverage[i] == 0)
        gapInConfirmedKmers = true;
  }

  //  If the coverage isn't perfect (ramp up, constant, ramp down), trash the whole read.

  if (g->discardImperfectCoverage == true) {
    uint32  bgn = clrBgn;
    uint32  end = clrEnd - 1;

    while ((bgn + 1 < end) && (coverage[bgn] < coverage[bgn+1]))
      bgn++;

    while ((bgn + 1 < end) && (coverage[end-1] > coverage[end]))
      end--;

    uint32  bgnc = coverage[bgn];
    uint32  endc = coverage[end];

    if (VERBOSE)
      fprintf(stderr, "IMPERFECT: bgn=%u %u  end=%u %u\n",
              bgn, bgnc,
              end, endc);

    if (bgnc != endc)
      imperfectKmerCoverage = true;

    else
      for (uint32 i=bgn; i<=end; i++)
        if (coverage[i] != bgnc)
          imperfectKmerCoverage = true;
  }

  if (VERBOSE > 1) {
    fprintf(stderr, "TRIM: %d,%d (post)\n", clrBgn, clrEnd);
    dump(stderr, "TRIM");
  }  //  VERBOSE
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
}





void
mertrimWorker(void *G, void *T, void *S) {
  mertrimGlobalData    *g = (mertrimGlobalData  *)G;
  mertrimThreadData    *t = (mertrimThreadData  *)T;
  mertrimComputation   *s = (mertrimComputation *)S;

  if (VERBOSE)
    fprintf(stderr, "\nPROCESS\n");

  s->t = t;

  uint32  eval = s->evaluate();

  //  Attempt correction if there are kmers to correct from.

  if (eval == ATTEMPTCORRECTION) {
    s->reverse();
    s->analyze();
    s->attemptCorrection(true);

    s->reverse();
    s->analyze();
    s->attemptCorrection(false);
  }

  //  Attempt trimming if the read wasn't perfect

  if (eval != ALLGOOD) {
    s->analyze();
    s->attemptTrimming();
  }

  //  SKIPPING until the heuristics are worked out.
  //s->analyzeChimer();

  if (VERBOSE)
    s->dump(stderr, "FINAL");
}


mertrimComputation *
mertrimReaderGatekeeper(mertrimGlobalData *g) {
  mertrimComputation   *s = NULL;

  while ((g->gktCur <= g->gktEnd) &&
         (s == NULL)) {
    s = new mertrimComputation();

    g->gkRead->gkStore_getFragment(g->gktCur, &s->fr, GKFRAGMENT_QLT);
    g->gktCur++;

    if ((g->forceCorrection) ||
        (g->gkRead->gkStore_getLibrary(s->fr.gkFragment_getLibraryIID())->doTrim_initialMerBased)) {
      s->initializeGatekeeper(g);
    } else {
      delete s;
      s = NULL;
    }
  }

  return(s);
}


mertrimComputation *
mertrimReaderFASTQ(mertrimGlobalData *g) {
  mertrimComputation   *s = new mertrimComputation();

  if (s->initializeFASTQ(g) == false) {
    delete s;
    s = NULL;
  }

  return(s);
}


void *
mertrimReader(void *G) {
  mertrimGlobalData    *g = (mertrimGlobalData  *)G;
  mertrimComputation   *s = NULL;

  if (g->gkRead)
    s = mertrimReaderGatekeeper(g);

  if (g->fqInput)
    s = mertrimReaderFASTQ(g);

  return(s);
}



void
mertrimWriterGatekeeper(mertrimGlobalData *g, mertrimComputation *s) {
  mertrimResult         res;

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
}


void
mertrimWriterFASTQ(mertrimGlobalData *g, mertrimComputation *s) {

  //  Note that getClrBgn/getClrEnd return positions in the ORIGINAL read, not the correct read.
  //  DO NOT USE HERE!

  uint32   seqOffset = 0;
  char     label[16];

  if ((s->suspectedChimer       == false) &&
      (s->hasNoConfirmedKmers   == false) &&
      (s->gapInConfirmedKmers   == false) &&
      (s->imperfectKmerCoverage == false)) {
    strcpy(label, "CLEAN");  //  Free of chimer!

    seqOffset = s->clrBgn;

    s->corrSeq[s->clrEnd] = 0;
    s->corrQlt[s->clrEnd] = 0;

  } else if (s->hasNoConfirmedKmers == true) {
    strcpy(label, "NOKMER");

    seqOffset = 0;

    s->corrSeq[0] = 0;
    s->corrQlt[0] = 0;

  } else if (s->gapInConfirmedKmers == true) {
    strcpy(label, "GAPKMER");

    seqOffset = 0;

    s->corrSeq[0] = 0;
    s->corrQlt[0] = 0;

  } else if (s->imperfectKmerCoverage == true) {
    strcpy(label, "IMPERFECTKMER");

    seqOffset = 0;

    s->corrSeq[0] = 0;
    s->corrQlt[0] = 0;

  } else if (s->suspectedChimerBgn - s->clrBgn >= AS_READ_MIN_LEN) {
    strcpy(label, "MP");  //  Junction read, longest portion would make an MP pair

    seqOffset = s->clrBgn;

    s->corrSeq[s->suspectedChimerBgn] = 0;
    s->corrQlt[s->suspectedChimerBgn] = 0;

  } else if (s->clrEnd - s->suspectedChimerEnd >= AS_READ_MIN_LEN) {
    strcpy(label, "PE");  //  Junction read, longest portion would make a PE pair

    seqOffset = s->suspectedChimerEnd;

    s->corrSeq[s->clrEnd] = 0;
    s->corrQlt[s->clrEnd] = 0;

  } else {
    strcpy(label, "SHORT1");//  Junction read, neither side is long enough to save, so save nothing.

    seqOffset = 0;

    s->corrSeq[0] = 0;
    s->corrQlt[0] = 0;
  }

  if (strlen(s->corrSeq + seqOffset) < 64) {
    strcpy(label, "SHORT2");

    seqOffset = 0;

    s->corrSeq[0] = 0;
    s->corrQlt[0] = 0;
  }

  fprintf(g->fqLog, F_U32"\t"F_U32"\tchimer="F_U32"\t"F_U32"\t"F_U32"\t%s\t%s\n",
          s->clrBgn,
          s->clrEnd,
          s->suspectedChimer,
          s->suspectedChimerBgn,
          s->suspectedChimerEnd,
          label,
          s->readName);

  //  Convert from CA QV to Sanger QV
  for (uint32 i=0; s->corrQlt[i]; i++) {
    s->corrQlt[i] -= '0';
    s->corrQlt[i] += '!';
  }

  fprintf(g->fqOutput, "%s type=%s\n%s\n+\n%s\n",
          s->readName,
          label,
          s->corrSeq + seqOffset,
          s->corrQlt + seqOffset);
}


void
mertrimWriter(void *G, void *S) {
  mertrimGlobalData    *g = (mertrimGlobalData  *)G;
  mertrimComputation   *s = (mertrimComputation *)S;

  //  The get*() functions return positions in the original uncorrected sequence.  They
  //  map from positions in the corrected sequence (which has inserts and deletes) back
  //  to the original sequence.
  //
  assert(s->getClrBgn() <= s->getClrEnd());
  //assert(s->getClrBgn() <  s->getSeqLen());
  assert(s->getClrEnd() <= s->getSeqLen());

  if (g->resFile)
    mertrimWriterGatekeeper(g, s);

  if (g->fqOutput)
    mertrimWriterFASTQ(g, s);

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

    } else if (strcmp(argv[arg], "-F") == 0) {
      g->fqInputPath = argv[++arg];

    } else if (strcmp(argv[arg], "-m") == 0) {
      g->merSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      g->compression = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-mc") == 0) {
      g->merCountsFile = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      g->numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-b") == 0) {
      g->gktBgn = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      g->gktEnd = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-v") == 0) {
      g->beVerbose = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      VERBOSE++;
 
    } else if (strcmp(argv[arg], "-coverage") == 0) {
      g->actualCoverage = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-correct") == 0) {
      g->minCorrectFraction = atof(argv[++arg]);
      if (g->minCorrectFraction >= 1) {
        g->minCorrect         = (uint32)g->minCorrectFraction;
        g->minCorrectFraction = 0;
      }

    } else if (strcmp(argv[arg], "-evidence") == 0) {
      g->minVerifiedFraction = atof(argv[++arg]);
      if (g->minVerifiedFraction >= 1) {
        g->minVerified         = (uint32)g->minVerifiedFraction;
        g->minVerifiedFraction = 0;
      }

    } else if (strcmp(argv[arg], "-endtrim") == 0) {
      if (g->endTrimDefault == true)
        g->endTrimNum = 0;
      g->endTrimWinScale[g->endTrimNum] = atof(argv[++arg]);
      g->endTrimErrAllow[g->endTrimNum] = atoi(argv[++arg]);
      g->endTrimNum++;

    } else if (strcmp(argv[arg], "-discardzero") == 0) {
      g->discardZeroCoverage = true;

    } else if (strcmp(argv[arg], "-discardimperfect") == 0) {
      g->discardImperfectCoverage = true;

    } else if (strcmp(argv[arg], "-f") == 0) {
      g->forceCorrection = true;
 
    } else if (strcmp(argv[arg], "-NM") == 0) {
      g->correctMismatch = false;
 
    } else if (strcmp(argv[arg], "-NI") == 0) {
      g->correctIndel = false;
 
    } else if (strcmp(argv[arg], "-o") == 0) {
      g->fqOutputPath = argv[++arg];
      g->resPath      = argv[arg];

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((g->gkpPath == 0L) && (g->fqInputPath == 0L))
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -g gkpStore -m merSize -mc merCountsFile [-v]\n", argv[0]);
    exit(1);
  }

  fprintf(stderr, "zerocoverage: %d\n", g->discardZeroCoverage);
  fprintf(stderr, "imperfect:    %d\n", g->discardImperfectCoverage);

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

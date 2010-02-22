
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

const char *mainid = "$Id: merTrim.C,v 1.1 2010-02-22 06:38:40 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlapStore.h"

#include "AS_MER_gkpStore_to_FastABase.H"

#include "bio++.H"
#include "sweatShop.H"
#include "existDB.H"
#include "positionDB.H"
#include "libmeryl.H"

//
//  THIS FAILS TO CORRECT ERRORS IN THE FIRST MER.  We need to reverse the string and try from the other end.
//

#define EXISTDB_MIN_COUNT   3


class mertrimGlobalData {
public:
  mertrimGlobalData() {
    gkpPath        = 0L;
    merCountsFile  = 0L;
    merSize        = 23;
    compression    = 1;
    numThreads     = 4;
    beVerbose      = false;
    kb             = NULL;
    edb            = NULL;
  };

  ~mertrimGlobalData() {
    delete kb;
    delete edb;
  };

  void              initialize(void) {
    kb  = new kMerBuilder(merSize, compression, 0L);
    edb = new existDB(merCountsFile, merSize, existDBcounts, EXISTDB_MIN_COUNT, ~0);
  };

public:

  //  Command line parameters
  //
  char    *gkpPath;
  char    *merCountsFile;
  uint32   merSize;
  uint32   compression;
  uint32   numThreads;
  bool     beVerbose;

  //  Global data
  //
  kMerBuilder  *kb;
  existDB      *edb;
};





#define ALLGOOD 1
#define MISSINGMERS 2
#define ALLCRAP 3
#define ATTEMPTCORRECTION 4


class mertrimComputation {
public:
  mertrimComputation(mertrimGlobalData *g_, gkFragment &fr_) {
    g = g_;

    readIID = fr_.gkFragment_getReadIID();
    seqLen  = fr_.gkFragment_getSequenceLength();

    origSeq = new char [AS_READ_MAX_NORMAL_LEN];
    corrSeq = new char [AS_READ_MAX_NORMAL_LEN];
    origQlt = new char [AS_READ_MAX_NORMAL_LEN];

    nMersExpected = 0;
    nMersTested   = 0;
    nMersFound    = 0;

    ms = NULL;

    disconnect = NULL;
    coverage   = NULL;

    corrected  = NULL;

    strcpy(origSeq, fr_.gkFragment_getSequence());
    strcpy(corrSeq, fr_.gkFragment_getSequence());
    strcpy(origQlt, fr_.gkFragment_getQuality());

    clrBgn = 0;
    clrEnd = seqLen;
  };
  ~mertrimComputation() {
    delete [] origSeq;
    delete [] corrSeq;
    delete [] origQlt;
    delete    ms;
    delete [] disconnect;
    delete [] coverage;
    delete [] corrected;
  }

  uint32     evaluate(bool postCorrection);
  uint32     attemptCorrection(void);
  uint32     attemptTrimming(void);

  uint32     getNumUncoveredSpaces(void) { return(nHole); };
  uint32     getNumUncoveredBases(void)  { return(nZero); };

  uint32     getNumCorrected(void)       { return(nCorr); };

  uint32     getClrBgn(void) { return(clrBgn); };
  uint32     getClrEnd(void) { return(clrEnd); };

  void       dump(FILE *F, char *label);

private:
  mertrimGlobalData   *g;

  AS_IID     readIID;
  uint32     seqLen;

  char      *origSeq;
  char      *corrSeq;
  char      *origQlt;

  uint32     nMersExpected;
  uint32     nMersTested;
  uint32     nMersFound;

  merStream *ms;

  uint32     clrBgn;
  uint32     clrEnd;

  uint32    *disconnect;  //  per base - a hole before this base
  uint32    *coverage;    //  per base - mer coverage

  //  per base: C - this base corrected
  //            A - this base failed to be corrected, no answer
  //            X - this base failed to be corrected, conflicts
  uint32    *corrected;

  uint32     nHole;  //  Number of spaces (between bases) with no mer coverage
  uint32     nZero;  //  Number of bases with no mer coverage

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
    ms = new merStream(g->kb, new seqStream(corrSeq, seqLen), false, true);

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
      nMersFound++;
  }

  if (disconnect && coverage && corrected) {
    nHole = 0;
    nZero = 0;

    nCorr = 0;
    nFail = 0;
    nConf = 0;

    //  Count the number of holes and gaps over the clear range.
    //  Also note any corrected bases here.
    //
    for (uint32 pos=clrBgn; pos<clrEnd; pos++) {
      if (disconnect[pos] != 0)
        nHole++;
      if (coverage[pos] == 0)
        nZero++;
      if (corrected[pos] == 'C')
        nCorr++;
    }

    //  Failed and Conflicting mers do NOT count against us unless the mer is completely within the
    //  clear range.
    //
    for (uint32 pos=clrBgn + g->merSize - 1; pos<clrEnd; pos++) {
      if (corrected[pos] == 'F')
        nFail++;
      if (corrected[pos] == 'X')
        nConf++;
    }

    //  We used to check that if nMersFound == nMersExpected (that is, ALLGOOD), then the clear
    //  range sequence had no holes, was all covered by mers, no mers failed or were conflicts.
    //  None of those conditions are strictly true when the clear range is not the whole read.
  }

  //fprintf(stderr, "evaluate()-- nMersFound=%d nMersTested=%d nMersExpected=%d\n",
  //        nMersFound, nMersTested, nMersExpected);

  if (nMersFound == nMersExpected)
    //  All mers confirmed, read is 100% verified!
    return(ALLGOOD);

  if (nMersFound == nMersTested)
    //  All mers we could test were found.  We missed some of the expected mers
    //  due to N's in the sequence.
    return(MISSINGMERS);

  if (nMersFound == 0)
    //  No kMers confirmed, read is 100% garbage (or 100% unique).
    return(ALLCRAP);

  //  Attempt correction.
  return(ATTEMPTCORRECTION);
}



uint32
mertrimComputation::attemptCorrection(void) {

  coverage   = new uint32 [AS_READ_MAX_NORMAL_LEN + 1];
  disconnect = new uint32 [AS_READ_MAX_NORMAL_LEN + 1];
  corrected  = new uint32 [AS_READ_MAX_NORMAL_LEN + 1];

  memset(coverage,   0, sizeof(uint32) * (AS_READ_MAX_NORMAL_LEN + 1));
  memset(disconnect, 0, sizeof(uint32) * (AS_READ_MAX_NORMAL_LEN + 1));
  memset(corrected,  0, sizeof(uint32) * (AS_READ_MAX_NORMAL_LEN + 1));

  assert(ms != NULL);

  ms->rewind();

  while (ms->nextMer()) {

    //  If the kMer is not present in our list of good mers, search for a kmer that does exist.
    //
    if (g->edb->exists(ms->theCMer()) == false) {
      kMer    testF = ms->theFMer();
      kMer    testR = ms->theRMer();

      uint32  numReplacements = 0;
      char    replacementBase = 0;

      //  Replace the last base in the F mer with all four bases.  The mechanism here is to shift
      //  the last base out (using -= to append a base on the beginning of the mer), then to shift
      //  on a new last base (using += to append a base on the end of the mer).  Instead of
      //  reverse complementing the mer a bunch of times (somewhat expensive) we do the same
      //  operations on the R mer.

      testF -= 0x00;
      testF += letterToBits['A'];   testF.mask(true);
      testR += 0x00;                testR.mask(false);
      testR -= letterToBits['T'];
      if (g->edb->exists(testF) || g->edb->exists(testR)) {
        numReplacements++;
        replacementBase = 'A';
      }

      testF -= 0x00;
      testF += letterToBits['C'];    testF.mask(true);
      testR += 0x00;                 testR.mask(false);
      testR -= letterToBits['G'];
      if (g->edb->exists(testF) || g->edb->exists(testR)) {
        numReplacements++;
        replacementBase = 'C';
      }

      testF -= 0x00;
      testF += letterToBits['G'];    testF.mask(true);
      testR += 0x00;                 testR.mask(false);
      testR -= letterToBits['C'];
      if (g->edb->exists(testF) || g->edb->exists(testR)) {
        numReplacements++;
        replacementBase = 'G';
      }

      testF -= 0x00;
      testF += letterToBits['T'];    testF.mask(true);
      testR += 0x00;                 testR.mask(false);
      testR -= letterToBits['A'];
      if (g->edb->exists(testF) || g->edb->exists(testR)) {
        numReplacements++;
        replacementBase = 'T';
      }

      uint32 pos = ms->thePositionInSequence() + g->merSize - 1;

      if (numReplacements == 0)
        corrected[pos] = 'F';

      if (numReplacements > 1)
        corrected[pos] = 'X';

      if (numReplacements == 1) {
        corrected[pos] = 'C';

        corrSeq[pos] = replacementBase;

        //fprintf(stderr, "Correct read %d at position %d from %c to %c (QV %d)\n",
        //        readIID, pos, origSeq[pos], corrSeq[pos], origQlt[pos]);

        //  Rebuild the merStream to use the corrected sequence, then move to the same position in
        //  the sequence.  This is done because we cannot simply change the string -- we need to
        //  change the state of the kMerBuilder associated with the merStream, and we can't do that.

        pos = ms->thePositionInSequence();

        delete ms;
        ms = new merStream(g->kb, new seqStream(corrSeq, seqLen), false, true);
        ms->nextMer();

        while (pos != ms->thePositionInSequence())
          ms->nextMer();
      }
    }

    //  After any corrections for this base are done, update the analysis of the read.  We want to
    //  know if there are any bases not covered by mers (coverage[]) or if there are any discontinuities
    //  in the mer coverage (disconnect[]).

    if (g->edb->exists(ms->theCMer())) {
      u32bit  p = ms->thePositionInSequence();

      //  This is an indication that our read is missing a base; there are two mers abutting
      //  together, with no mer covering the junction.
      if ((p != 0) && (coverage[p-1] != 0) && (coverage[p] == 0))
        disconnect[p] = 1;

      for (u32bit add=0; add < g->merSize; add++)
        coverage[p + add]++;

      assert(p + g->merSize - 1 < seqLen);
    }
  }  //  Over all mers


  return(evaluate(true));
}


//  A BADBASE is one that was:
//    corrected
//    caused a mer to not exist
//    has a low QV
//
#define BADBASE(i) (((corrected[i] == 'C') || (origQlt[i] < '5')) ? 1 : 0)

uint32
mertrimComputation::attemptTrimming(void) {

  //  Any 5-base window with more two or more errors or low quality.  Trim to the inner-most.

  if (corrected == NULL) {
    corrected  = new uint32 [AS_READ_MAX_NORMAL_LEN + 1];
    memset(corrected,  0, sizeof(uint32) * (AS_READ_MAX_NORMAL_LEN + 1));
  }

  uint32  wErr[AS_READ_MAX_NORMAL_LEN] = {0};

  for (uint32 i=0; i<seqLen-4; i++) {
    wErr[i] += (BADBASE(i+0) +
                BADBASE(i+1) +
                BADBASE(i+2) +
                BADBASE(i+3) +
                BADBASE(i+4));
  }

  //  Find the largest region with < 2 errors.

  uint32  lBgn=0, lEnd=0;  //  Largest found so far
  uint32  tBgn=0, tEnd=0;  //  The one we're looking at

  for (uint32 i=0; i<seqLen-3; i++) {
    if ((i == seqLen - 4) || (wErr[i] > 1)) {
      if (tEnd - tBgn > lBgn - lEnd) {
        lBgn = tBgn;
        lEnd = tEnd;
      }
      tBgn = i;
      tEnd = i;
    } else {
      tEnd = i;
    }
  }

  if (tEnd - tBgn > lBgn - lEnd) {
    lBgn = tBgn;
    lEnd = tEnd;
  }

  clrBgn = lBgn;
  clrEnd = lEnd;

  //fprintf(stderr, "TRIM: %d,%d (pre)\n", clrBgn, clrEnd);

  //  Adjust borders to the first bad base.  Remember, clrBgn includes the base, clrEnd does not.

  while ((clrBgn > 0) && (BADBASE(clrBgn) == 0))
    clrBgn--;
  if (BADBASE(clrBgn))
    clrBgn++;

  while ((clrEnd < seqLen) && (BADBASE(clrEnd) == 0))
    clrEnd++;

  //fprintf(stderr, "TRIM: %d,%d (post)\n", clrBgn, clrEnd);

  //if (clrEnd - clrBgn < AS_READ_MIN_LEN)
  //  dump(stdout, "SHORT");

  return(evaluate(true));
}


void
mertrimComputation::dump(FILE *F, char *label) {
  fprintf(F, "%s read %d len %d (trim %d-%d) with %d uncovered bases and %d uncovered spaces.  corrected %d conflicting %d failed %d\n",
          label, readIID, seqLen, getClrBgn(), getClrEnd(), nZero, nHole, nCorr, nConf, nFail);
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
    fprintf(F, "%c", origQlt[i]);
    if (i+1 == clrEnd)
      fprintf(F, "]-");
  }
  fprintf(F, " (QLT)\n");
  if (coverage) {
    for (uint32 i=0; i<seqLen; i++) {
      if (i == clrBgn)
        fprintf(F, "-[");
      fprintf(F, "%c", coverage[i] + 'A');
      if (i+1 == clrEnd)
        fprintf(F, "]-");
    }
    fprintf(F, " (COVERAGE)\n");
  }
  if (corrected) {
    for (uint32 i=0; i<seqLen; i++) {
      if (i == clrBgn)
        fprintf(F, "-[");
      fprintf(F, "%c", (corrected[i]) ? corrected[i] : ' ');
      if (i+1 == clrEnd)
        fprintf(F, "]-");
    }
    fprintf(F, " (CORRECTIONS)\n");
  }
}



int
main(int argc, char **argv) {
  mertrimGlobalData  *g = new mertrimGlobalData;

  //assert(sizeof(kmerhit) == 8);

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

#if 0
    } else if (strcmp(argv[arg], "-t") == 0) {
      g->numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-tb") == 0) {
      g->tBeg = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-te") == 0) {
      g->tEnd = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-qb") == 0) {
      g->qBeg = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-qe") == 0) {
      g->qEnd = atoi(argv[++arg]);
#endif

    } else if (strcmp(argv[arg], "-v") == 0) {
      g->beVerbose = true;

#if 0
    } else if (strcmp(argv[arg], "-o") == 0) {
      g->outputName = argv[++arg];
#endif

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((g->gkpPath == 0L) || (err)) {
    exit(1);
  }

  gkpStoreFile::registerFile();

  //  Read mers into an existDB.

  g->initialize();

  //  Open gkpStore, stream reads, compute trim.

  fprintf(stderr, "opening gkStore '%s'\n", g->gkpPath);

  gkStore     *gk = new gkStore(g->gkpPath, FALSE, TRUE);

  gk->gkStore_enableClearRange(AS_READ_CLEAR_OBTINITIAL);

  uint32       tBeg = 1;
  uint32       tEnd = gk->gkStore_getNumFragments();

  gkFragment   fr;
  //gkStream    *gs = new gkStream(gk, tBeg, tEnd, GKFRAGMENT_QLT);

  uint32       resultCorrect      = 0;
  uint32       resultCrap         = 0;
  uint32       resultNotAMer      = 0;
  uint32       resultCorrected    = 0;
  uint32       resultTrimmed      = 0;
  uint32       resultUncovered    = 0;
  uint32       resultDeleted      = 0;

  uint32       resultCorrectedHistogram[100] = {0};

  //  while (gs->next(&fr)) {

  for (uint32 iid=1; iid<=gk->gkStore_getNumFragments(); iid++) {
    gk->gkStore_getFragment(iid, &fr, GKFRAGMENT_QLT);

    bool                delFrag = false;
    mertrimComputation *C = new mertrimComputation(g, fr);

    //assert((tBeg <= iid) && (iid <= tEnd));

    uint32  eval = C->evaluate(false);

    //  Read is perfect!
    //
    if      (eval == ALLGOOD) {
      resultCorrect++;
      goto finished;
    }

    //  Read had NO mers found.  We'll let it through, but it might just end up a singleton.  There
    //  is a slight chance it has a 2-copy mer that will pull it in, or maybe the mate will.
    //
    if (eval == ALLCRAP) {
      resultCrap++;
      C->attemptTrimming();
      //C->dump(stderr, "CRAP");
      goto finished;
    }

    eval = C->attemptCorrection();

    //  Correction worked perfectly, we're done.
    //
    if ((eval == ALLGOOD) && (C->getNumCorrected() < 4)) {
      resultCorrected++;
      resultCorrectedHistogram[C->getNumCorrected()]++;
      goto finished;
    }

    eval = C->attemptTrimming();

    if ((eval == ALLGOOD) && (C->getNumCorrected() < 4)) {
      resultTrimmed++;
      //resultCorrectedHistogram[C->getNumCorrected()]++;
      goto finished;
    }

    //  And then what??
    //C->dump(stderr, "TRIM");
    resultUncovered++;

    //  DO NOT delete.
    //delFrag = true;

  finished:

    if (C->getClrEnd() - C->getClrBgn() < AS_READ_MIN_LEN)
      delFrag = true;

    fr.gkFragment_setClearRegion(C->getClrBgn(), C->getClrEnd(), AS_READ_CLEAR_OBTINITIAL);
    gk->gkStore_setFragment(&fr);

    if (delFrag) {
      resultDeleted++;
      gk->gkStore_delFragment(fr.gkFragment_getReadIID(), false);
    }

    delete C;
  }

  //delete gs;
  delete gk;

  fprintf(stderr, "resultCorrect        %d\n", resultCorrect);
  fprintf(stderr, "resultCrap           %d\n", resultCrap);
  fprintf(stderr, "resultNotAMer        %d\n", resultNotAMer);
  fprintf(stderr, "resultCorrected      %d\n", resultCorrected);
  fprintf(stderr, "resultTrimmed        %d\n", resultTrimmed);
  fprintf(stderr, "resultUncovered      %d\n", resultUncovered);
  fprintf(stderr, "resultDeleted        %d\n", resultDeleted);

  //for (uint32 i=0; i<100; i++)
  //  fprintf(stderr, "%d\n", resultCorrectedHistogram[i]);

  fprintf(stderr, "\nSuccess!  Bye.\n");

  return(0);
}

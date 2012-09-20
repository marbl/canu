
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
static char *rcsid = "$Id: CIScaffoldT_Merge_CGW.c,v 1.85 2012-09-20 19:18:40 brianwalenz Exp $";

//
//  The ONLY exportable function here is MergeScaffoldsAggressive.
//

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Instrument_CGW.h"
#include "ChiSquareTest_CGW.h"
#include "InterleavedMerging.h"

#include "CIScaffoldT_MergeScaffolds.h"

#include "CIScaffoldT_Analysis.H"

#include <vector>
#include <algorithm>
using namespace std;


#define PREFERRED_GAP_SIZE  (-500)

#define CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD MIN_EDGES

#define EDGE_QUANTA 5.0

#define EDGE_STRENGTH_FACTOR  MIN_EDGES

#undef ALTERNATE_MATE_TEST_RULES

//  Map an edge ID to a mate pair test result.
#ifdef TRACK_MATE_PAIR_TEST
map<uint32, mergeTestResult>    matePairTestResult;
#endif


int
isBadScaffoldMergeEdge(SEdgeT * edge, ChunkOverlapperT * overlapper);


static
int
TouchesMarkedScaffolds(SEdgeT *curEdge) {
  CIScaffoldT *scaffoldA = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idA);
  CIScaffoldT *scaffoldB = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idB);

  if (curEdge->orient.isAB_AB()) {
    if ((scaffoldA->numEssentialB == 0) &&
        (scaffoldB->numEssentialA == 0)) {
      assert(scaffoldA->essentialEdgeB == NULLINDEX);
      assert(scaffoldB->essentialEdgeA == NULLINDEX);
      return FALSE;
    } else {
      return TRUE;
    }
  }

  if (curEdge->orient.isAB_BA()) {
    if ((scaffoldA->numEssentialB == 0) &&
        (scaffoldB->numEssentialB == 0)) {
      assert(scaffoldA->essentialEdgeB == NULLINDEX);
      assert(scaffoldB->essentialEdgeB == NULLINDEX);
      return FALSE;
    } else {
      return TRUE;
    }
  }

  if (curEdge->orient.isBA_AB()) {
    if ((scaffoldA->numEssentialA == 0) &&
        (scaffoldB->numEssentialA == 0)) {
      assert(scaffoldA->essentialEdgeA == NULLINDEX);
      assert(scaffoldB->essentialEdgeA == NULLINDEX);
      return FALSE;
    } else {
      return TRUE;
    }
  }

  if (curEdge->orient.isBA_BA()) {
    if ((scaffoldA->numEssentialA == 0) &&
        (scaffoldB->numEssentialB == 0)) {
      assert(scaffoldA->essentialEdgeA == NULLINDEX);
      assert(scaffoldB->essentialEdgeB == NULLINDEX);
      return FALSE;
    } else {
      return TRUE;
    }
  }

  assert(0);
  return(FALSE);
}









#define MAX_FRAC_BAD_TO_GOOD .3






static
int
isQualityScaffoldMergingEdgeNEW(SEdgeT                     *curEdge,
                                CIScaffoldT                *scaffoldA,
                                CIScaffoldT                *scaffoldB,
                                ScaffoldInstrumenter       *si,
                                VA_TYPE(MateInstrumenterP) *MIs,
                                double                       minSatisfied,
                                double                       maxDelta) {

#warning GET RID OF THIS STATIC
  static vector<instrumentLIB>   libs;
  static bool                    libsInit = false;

  if (libsInit == false) {
    libsInit = true;

    for (int32 i=0; i<GetNumDistTs(ScaffoldGraph->Dists); i++) {
      DistT *dptr = GetDistT(ScaffoldGraph->Dists, i);

      libs.push_back(instrumentLIB(i, dptr->mu, dptr->sigma, true));
    }
  }

  // NOTE, we should cache these single scaffold instrumenters.

  instrumentSCF   A(scaffoldA);
  A.analyze(libs);

  instrumentSCF   B(scaffoldB);
  B.analyze(libs);

  instrumentSCF   P(scaffoldA, curEdge, scaffoldB);
  P.analyze(libs);

#define VERBOSE_QUALITY_MERGE_EDGE

#ifdef VERBOSE_QUALITY_MERGE_EDGE
  fprintf(stderr, "isQualityScaffoldMergingEdge()--   scaffold %d instrumenter happy %.1f gap %.1f misorient close %.1f correct %.1f far %.1f oriented close %.1f far %.1f missing %.1f external %.1f\n",
          scaffoldA->id,
          A.numHappy, A.numGap, A.numMisClose, A.numMis, A.numMisFar, A.numTooClose, A.numTooFar, A.numMissing, A.numExternal);

  fprintf(stderr, "isQualityScaffoldMergingEdge()--   scaffold %d instrumenter happy %.1f gap %.1f misorient close %.1f correct %.1f far %.1f oriented close %.1f far %.1f missing %.1f external %.1f\n",
          scaffoldB->id,
          B.numHappy, B.numGap, B.numMisClose, B.numMis, B.numMisFar, B.numTooClose, B.numTooFar, B.numMissing, B.numExternal);

  fprintf(stderr, "isQualityScaffoldMergingEdge()--   scaffold (new) instrumenter happy %.1f gap %.1f misorient close %.1f correct %.1f far %.1f oriented close %.1f far %.1f missing %.1f external %.1f\n",
          P.numHappy, P.numGap, P.numMisClose, P.numMis, P.numMisFar, P.numTooClose, P.numTooFar, P.numMissing, P.numExternal);
#endif

#if 0
  //  Contig positions in the mock merged scaffold are messed up.  The greedy gap size estimator
  //  seems to over estimate gap sizes, which might be better for the few sizes that are messed up,
  //  but also over estimates the size for the correct gaps.
  //
  P.estimateGaps(libs);
  P.analyze(libs);
  fprintf(stderr, "isQualityScaffoldMergingEdge()--   scaffold (new) instrumenter happy %.1f gap %.1f misorient close %.1f correct %.1f far %.1f oriented close %.1f far %.1f missing %.1f external %.1f\n",
          P.numHappy, P.numGap, P.numMisClose, P.numMis, P.numMisFar, P.numTooClose, P.numTooFar, P.numMissing, P.numExternal);
#endif

  int32  mBeforeGood = A.numEcstatic + B.numEcstatic;
  int32  mAfterGood  = P.numEcstatic;

  int32  mBefore     = A.numMateInternal + A.numMateExternal + B.numMateInternal + B.numMateExternal;
  int32  mAfter      = P.numMateInternal + P.numMateExternal;

  //  Mates that fall into gaps are neither good nor bad.  They just don't exist as far as we are concerned.

  mBeforeGood -= A.numGap + B.numGap;
  mBefore     -= A.numGap + B.numGap;

  mAfterGood  -= P.numGap;
  mAfter      -= P.numGap;

  //  Mates before that fall into contigs, but aren't present, are still bad.

  int32  mBeforeBad = A.numDejected + B.numDejected;
  int32  mAfterBad  = P.numDejected;

  double fractMatesHappyBefore = (mAfter > 0) ? ((double)mBeforeGood / mAfter) : 0.0;
  double fractMatesHappyAfter  = (mAfter > 0) ? ((double)mAfterGood  / mAfter) : 0.0;

  //  The number of satisfied mates we expect to gain from this edge.
  //  Why -1?  Cheap way of fixing rounding issues.

  int32  mEdge = curEdge->edgesContributing * fractMatesHappyBefore - 1;

  fprintf(stderr, "isQualityScaffoldMergingEdge()--   before: %.3f satisfied (%d/%d good/bad mates)  after: %.3f satisfied (%d/%d good/bad mates)\n",
          fractMatesHappyBefore, mBeforeGood, mBeforeBad,
          fractMatesHappyAfter,  mAfterGood,  mAfterBad);


#ifdef TRACK_MATE_PAIR_TEST
  //  Why is this aborted?  Because the edge itself doesn't have an ID, and the only way
  //  to get an ID is by asing the graph.  We don't have a graph handy, and are buried
  //  FAR below where we last had it.
  //
  mergeTestResult  MTR;

  MTR.edgeID            = edgeID;  //GetVAIndex_CIEdgeT(graph->ContigGraph->edges, curEdge);
  MTR.edgeLength        = curEdge->distance.mean;
  MTR.edgeVariance      = curEdge->distance.variance;
  MTR.edgeWeight        = curEdge->edgesContributing;

  MTR.scaffoldAid       = scaffoldA->id;
  MTR.scaffoldAlength   = scaffoldA->bpLength.mean;
  MTR.scaffoldAcontigs  = scaffoldA->info.Scaffold.numElements;

  MTR.scaffoldBid       = scaffoldB->id;
  MTR.scaffoldBlength   = scaffoldB->bpLength.mean;
  MTR.scaffoldBcontigs  = scaffoldB->info.Scaffold.numElements;

  MTR.Ahappy            = A.numHappy;
  MTR.Agap              = A.numGap;
  MTR.Abad              = A.numDejected;

  MTR.Bhappy            = B.numHappy;
  MTR.Bgap              = B.numGap;
  MTR.Bbad              = B.numDejected;

  MTR.Mhappy            = P.numHappy;
  MTR.Mgap              = P.numGap;
  MTR.Mbad              = P.numDejected;

  MTR.beforeSatisfied   = fractMatesHappyBefore;
  MTR.beforeGood        = mBeforeGood;
  MTR.beforeBad         = mBeforeBad;

  MTR.afterSatisfied    = fractMatesHappyAfter;
  MTR.afterGood         = mAfterGood;
  MTR.afterBad          = mAfterBad;

  MTR.mergeAccepted     = false;
#endif



  //  NOTE, still need the metagenomic test

  bool failsOLD = false;
  bool failsNEW = false;

#if 1
  double badGoodRatio      = 1.0;

  if (mAfterGood > mBeforeGood)
    badGoodRatio = (double)(mAfterBad - mBeforeBad) / (double)(mAfterGood - mBeforeGood);

  bool failsMinimum       = (fractMatesHappyAfter < minSatisfied);
  bool failsToGetHappier1 = (fractMatesHappyAfter < fractMatesHappyBefore);

#ifndef ALTERNATE_MATE_TEST_RULES
  bool failsToGetHappier2 = (mAfterGood < mBeforeGood) || (badGoodRatio > MAX_FRAC_BAD_TO_GOOD);
  failsOLD = (failsMinimum && failsToGetHappier1 && failsToGetHappier2);
#else
#warning ALTERNATE RULES
#warning ALTERNATE RULES
#warning ALTERNATE RULES
  bool failsToGetHappier2 = (badGoodRatio > MAX_FRAC_BAD_TO_GOOD);
  failsOLD = (failsMinimum) && (failsToGetHappier1 || failsToGetHappier2);
#endif

  if (failsOLD)
    fprintf(stderr, "isQualityScaffoldMergingEdge()--   not happy enough to merge %d%d%d (%.3f < %.3f) && (%.3f < %.3f) && ((%d < %d) || (%0.3f > %.3f))\n",
            failsMinimum, failsToGetHappier1, failsToGetHappier2,
            fractMatesHappyAfter, minSatisfied,
            fractMatesHappyAfter, fractMatesHappyBefore,
            mAfterGood, mBeforeGood, badGoodRatio, MAX_FRAC_BAD_TO_GOOD);
  else
    fprintf(stderr, "isQualityScaffoldMergingEdge()--   ARE happy enough to merge %d%d%d (%.3f > %.3f) || (%.3f > %.3f) || ((%d > %d) && (%0.3f < %.3f))\n",
            failsMinimum, failsToGetHappier1, failsToGetHappier2,
            fractMatesHappyAfter, minSatisfied,
            fractMatesHappyAfter, fractMatesHappyBefore,
            mAfterGood, mBeforeGood, badGoodRatio, MAX_FRAC_BAD_TO_GOOD);
#endif

#if 0
  //  Looks good on paper, not ready for prime time.
  bool failsToGetHappierA = (fractMatesHappyAfter < fractMatesHappyBefore);           //  Satisfaction dropped
  bool failsToGetHappierB = (mAfterGood < mBeforeGood + mEdge);                       //  Good mates dropped

  failsNEW = (failsToGetHappierA || failsToGetHappierB);

  if (failsNEW)
    fprintf(stderr, "isQualityScaffoldMergingEdge()--   not happy enough to merge %d%d happiness (%.3f < %.3f) || mates (%d < %d + %d)\n",
            failsToGetHappierA, failsToGetHappierB,
            fractMatesHappyAfter, fractMatesHappyBefore,
            mAfterGood, mBeforeGood, mEdge);
  else
    fprintf(stderr, "isQualityScaffoldMergingEdge()--   ARE happy enough to merge %d%d happiness (%.3f >= %.3f) && mates (%d >= %d + %d)\n",
            failsToGetHappierA, failsToGetHappierB,
            fractMatesHappyAfter, fractMatesHappyBefore,
            mAfterGood, mBeforeGood, mEdge);

 if (failsNEW != failsOLD)
   fprintf(stderr, "isQualityScaffoldMergingEge()--  DIFFER.\n");
#endif

 if (GlobalData->mergeFilterLevel >= 5) {
   bool doTest1 = true;
   bool doTest2 = true;
   bool doTest3 = true;
   bool doTest4 = true;
   bool passTest1 = true;
   bool passTest2 = true;
   bool passTest3 = true;
   bool passTest4 = true;
   
   bool passAllTests = false;
   // We altered this test, which was probably too permissive. Using mAfter in denominator, it only required that happy mate count increased.
   //// fractMatesHappyBefore = (mAfter > 0) ? ((double)mBeforeGood / mAfter) : 0.0;
   fractMatesHappyBefore = (mBefore > 0) ? ((double)mBeforeGood / mBefore) : 0.0;
   fractMatesHappyAfter  = (mAfter > 0) ? ((double)mAfterGood  / mAfter) : 0.0;
   if (doTest1) {
     // Rule: Fraction of scored mates that are scored happy must not decrease.
     // These fractions can differ by tiny amounts. Use double precision.
     passTest1 = (fractMatesHappyAfter) >= (fractMatesHappyBefore); 
   }
   int32 observedHappyMateIncrease = (mAfterGood>mBeforeGood) ? (mAfterGood-mBeforeGood) : 0;
   int32 expectedHappyMateIncrease = (int32) ( fractMatesHappyBefore * curEdge->edgesContributing - 1);
   if (doTest2) {
     // Number of scored mates that are scored happy must increase by the expected amount.
     // Expected amount = (rate of happiness in scaffolds prior to merge) * (#mates in edge provoking the merge) - (fudge for integer conversion).
     passTest2 = (observedHappyMateIncrease >= expectedHappyMateIncrease);
   }
   int32 sadnessIncrease =   (mAfterBad >mBeforeBad ) ? (mAfterBad -mBeforeBad ) : 0;
   int32 happinessIncrease = (mAfterGood>mBeforeGood) ? (mAfterGood-mBeforeGood) : 0;
   if (doTest3) {
     // Sadness Increase must not exceed 30% of Happiness Increase.
     passTest3 = ( sadnessIncrease <= 1.0 * MAX_FRAC_BAD_TO_GOOD * happinessIncrease);
   }
   if (doTest4) {
     // Before each test, call MarkInternal to make sure all edges are properly marked. This should be a no-op.
     MarkInternalEdgeStatus(ScaffoldGraph, scaffoldA, 0, TRUE);  //  Merged
     MarkInternalEdgeStatus(ScaffoldGraph, scaffoldA, 0, FALSE);  //  Raw
     bool A_is_connected = (true==IsScaffold2EdgeConnected(ScaffoldGraph, scaffoldA));
     MarkInternalEdgeStatus(ScaffoldGraph, scaffoldB, 0, TRUE);  //  Merged
     MarkInternalEdgeStatus(ScaffoldGraph, scaffoldB, 0, FALSE);  //  Raw
     bool B_is_connected = (true==IsScaffold2EdgeConnected(ScaffoldGraph, scaffoldB));
     // Disconnected scaffolds should never happen but they do arise from interleaved scaffold merging.
     // Disconnects lead to a Least Squares failure (because one row is all zero in the coefficient matrix).
     // Least Squares failures induce black holes: scaffolds that grow to enormous size through iterations of merge scaffolds aggressive.
     // A continuing failure of Least Squares means gap sizes grow on every iteration.
     // A black hole's large gap sizes effectively suck in more contigs, which fuels further growth.
     // The test here will not fix a black hole but it can preclude any black hole growth.
     // This test says I to refuse to merge with any scaffold that shows signs of a disconnect problem.
     passTest4 = (A_is_connected) && (B_is_connected);
   }
   // Set the value that gets returned below...
   passAllTests = (passTest1) && (passTest2) && (passTest3) && (passTest4);
   failsOLD = !passAllTests;
   
   fprintf(stderr,"isQualityScaffoldMergingEdge()-- filter=5 pass=%d based on test1=%d, test2=%d, test3=%d (sad/happy=%d/%d), test4=%d.\n",
	   passAllTests,passTest1,passTest2,passTest3, sadnessIncrease, happinessIncrease, passTest4);
 }

 return(failsOLD == false);
}


static
int
isQualityScaffoldMergingEdgeOLD(SEdgeT                     *curEdge,
                                CIScaffoldT                *scaffoldA,
                                CIScaffoldT                *scaffoldB,
                                ScaffoldInstrumenter       *si,
                                VA_TYPE(MateInstrumenterP) *MIs,
                                double                       minSatisfied,
                                double                       maxDelta) {

  MateInstrumenter matesBefore;
  MateInstrumenter matesAfter;

  MateInstrumenter *sABefore = (GetNumVA_MateInstrumenterP(MIs) > scaffoldA->id) ? ((MateInstrumenter *)*GetVA_MateInstrumenterP(MIs, scaffoldA->id)) : NULL;
  MateInstrumenter *sBBefore = (GetNumVA_MateInstrumenterP(MIs) > scaffoldB->id) ? ((MateInstrumenter *)*GetVA_MateInstrumenterP(MIs, scaffoldB->id)) : NULL;

  double fractMatesHappyAfter  = 1.0;
  double fractMatesHappyBefore = 1.0;

  ResetMateInstrumenterCounts(&matesBefore);
  ResetMateInstrumenterCounts(&matesAfter);

  if (sABefore == NULL) {
    InstrumentScaffold(ScaffoldGraph, scaffoldA, si, InstrumenterVerbose2, stderr);

    sABefore = (MateInstrumenter *)safe_calloc(1, sizeof(MateInstrumenter));
    GetMateInstrumenterFromScaffoldInstrumenter(sABefore, si);
    SetVA_MateInstrumenterP(MIs, scaffoldA->id, &sABefore);
  }

  if (sBBefore == NULL) {
    InstrumentScaffold(ScaffoldGraph, scaffoldB, si, InstrumenterVerbose2, stderr);

    sBBefore = (MateInstrumenter *) safe_calloc(1, sizeof(MateInstrumenter));
    GetMateInstrumenterFromScaffoldInstrumenter(sBBefore, si);
    SetVA_MateInstrumenterP(MIs, scaffoldB->id, &sBBefore);
  }

  AddMateInstrumenterCounts(&matesBefore, sABefore);
  AddMateInstrumenterCounts(&matesBefore, sBBefore);

  InstrumentScaffoldPair(ScaffoldGraph, curEdge, si, InstrumenterVerbose2, stderr);

  GetMateInstrumenterFromScaffoldInstrumenter(&matesAfter, si);

#ifdef VERBOSE_QUALITY_MERGE_EDGE
  fprintf(stderr, "isQualityScaffoldMergingEdge()--   scaffold %d instrumenter (intra/inter): happy %d/%d bad %d/%d missing %d/%d\n",
          scaffoldA->id,
          GetMateStatsHappy(&sABefore->intra),   GetMateStatsHappy(&sABefore->inter),
          GetMateStatsBad(&sABefore->intra),     GetMateStatsBad(&sABefore->inter),
          GetMateStatsMissing(&sABefore->intra), GetMateStatsMissing(&sABefore->inter));

  fprintf(stderr, "isQualityScaffoldMergingEdge()--   scaffold %d instrumenter (intra/inter): happy %d/%d bad %d/%d missing %d/%d\n",
          scaffoldB->id,
          GetMateStatsHappy(&sBBefore->intra),   GetMateStatsHappy(&sBBefore->inter),
          GetMateStatsBad(&sBBefore->intra),     GetMateStatsBad(&sBBefore->inter),
          GetMateStatsMissing(&sBBefore->intra), GetMateStatsMissing(&sBBefore->inter));

  fprintf(stderr, "isQualityScaffoldMergingEdge()--   scaffold (new) instrumenter (intra/inter): happy %d/%d bad %d/%d missing %d/%d\n",
          GetMateStatsHappy(&matesAfter.intra),   GetMateStatsHappy(&matesAfter.inter),
          GetMateStatsBad(&matesAfter.intra),     GetMateStatsBad(&matesAfter.inter),
          GetMateStatsMissing(&matesAfter.intra), GetMateStatsMissing(&matesAfter.inter));
#endif

  int32   mBeforeGood = GetMateStatsHappy(&matesBefore.intra) + GetMateStatsHappy(&matesBefore.inter);
  int32   mBeforeBad  = GetMateStatsBad(&matesBefore.intra)   + GetMateStatsBad(&matesBefore.inter);

  int32   mAfterGood  = GetMateStatsHappy(&matesAfter.intra)  + GetMateStatsHappy(&matesAfter.inter);
  int32   mAfterBad   = GetMateStatsBad(&matesAfter.intra)    + GetMateStatsBad(&matesAfter.inter);

  //  Add in mates that should have been satisfied, but weren't.
#warning THIS IS TOO AGGRESSIVE
  mAfterBad += GetMateStatsMissing(&matesAfter.inter);

  //  This should only be counted for 'inter' (== inter-contig?) and not for 'intra'.
  assert(GetMateStatsMissing(&matesAfter.intra) == 0);

  // since we expect some set of mates to be missing due to divergence between closely related species, we don't perform this check the same way for metagenomics
  if (GetMateStatsMissing(&matesAfter.inter) > 0) {
	  if (((double)GetMateStatsMissing(&matesAfter.inter) / mAfterGood) < GlobalData->mergeScaffoldMissingMates || GlobalData->mergeScaffoldMissingMates == -1) {
		  mAfterBad -= GetMateStatsMissing(&matesAfter.inter);
	  }
  }

  int32   mBeforeSum  = mBeforeGood + mBeforeBad;
  int32   mAfterSum   = mAfterGood  + mAfterBad;

  if (mBeforeSum > 0)
    fractMatesHappyBefore = (double)mBeforeGood / mBeforeSum;

  if (mAfterSum > 0)
    fractMatesHappyAfter  = (double)mAfterGood / mAfterSum;

  fprintf(stderr, "isQualityScaffoldMergingEdge()--   before: %.3f satisfied (%d/%d good/bad mates)  after: %.3f satisfied (%d/%d good/bad mates)\n",
          fractMatesHappyBefore, mBeforeGood, mBeforeBad,
          fractMatesHappyAfter,  mAfterGood,  mAfterBad);

  if ((maxDelta > 0.0) &&
      (fractMatesHappyBefore - fractMatesHappyAfter > maxDelta)) {
    fprintf(stderr, "isQualityScaffoldMergingEdge()--   satisfied dropped by too much (%.3f > %.3f)\n",
            fractMatesHappyBefore - fractMatesHappyAfter, maxDelta);
    return(FALSE);
  }

  //  failsMinimum       -- true if the fraction happy after merging is below some arbitrary threshold.
  //
  //  failsToGetHappier1 -- true if the merged result is less happy than individually.
  //                     -- (but special case to allow exactly one bad mate in the merge -- old
  //                         option never used)
  //
  //  failsToGetHappier2 -- true if there are fewer happy mates after, or there are a whole lot more
  //                        bad mates after.  The original version of this test was screwing up the
  //                        compute of badGoodRatio, by omitting some parens.

  double badGoodRatio      = 1.0;

  if (mAfterGood > mBeforeGood)
    badGoodRatio = (double)(mAfterBad - mBeforeBad) / (double)(mAfterGood - mBeforeGood);

  bool  failsMinimum       = (fractMatesHappyAfter < minSatisfied);
  bool  failsToGetHappier1 = (fractMatesHappyAfter < fractMatesHappyBefore);
  bool  failsToGetHappier2 = (mAfterGood < mBeforeGood) || (badGoodRatio > MAX_FRAC_BAD_TO_GOOD);

  if (failsMinimum && failsToGetHappier1 && failsToGetHappier2) {
    fprintf(stderr, "isQualityScaffoldMergingEdge()--   not happy enough to merge (%.3f < %.3f) && (%.3f < %.3f) && ((%d < %d) || (%0.3f > %.3f))\n",
            fractMatesHappyAfter, minSatisfied,
            fractMatesHappyAfter, fractMatesHappyBefore,
            mAfterGood, mBeforeGood, badGoodRatio, MAX_FRAC_BAD_TO_GOOD);
    return(FALSE);
  } else {
    fprintf(stderr, "isQualityScaffoldMergingEdge()--   ARE happy enough to merge (%.3f < %.3f) && (%.3f < %.3f) && ((%d < %d) || (%0.3f > %.3f))\n",
            fractMatesHappyAfter, minSatisfied,
            fractMatesHappyAfter, fractMatesHappyBefore,
            mAfterGood, mBeforeGood, badGoodRatio, MAX_FRAC_BAD_TO_GOOD);
  }

  return(TRUE);
}




//static
int
isQualityScaffoldMergingEdge(SEdgeT                     *curEdge,
                             CIScaffoldT                *scaffoldA,
                             CIScaffoldT                *scaffoldB,
                             ScaffoldInstrumenter       *si,
                             VA_TYPE(MateInstrumenterP) *MIs,
                             double                       minSatisfied,
                             double                       maxDelta) {

  static int32  oldPASS = 0;
  static int32  oldFAIL = 0;
  static int32  newPASS = 0;
  static int32  newFAIL = 0;

  assert(si != NULL);

  if ((minSatisfied <= 0.0) &&
      (maxDelta     <= 0.0))
    return(TRUE);

  fprintf(stderr,"isQualityScaffoldMergingEdge()-- Merge scaffolds "F_CID" (%.1fbp) and "F_CID" (%.1fbp): gap %.1fbp +- %.1fbp weight %d %s edge\n",
          scaffoldA->id, scaffoldA->bpLength.mean,
          scaffoldB->id, scaffoldB->bpLength.mean,
          curEdge->distance.mean,
          sqrt(curEdge->distance.variance),
          curEdge->edgesContributing,
          ((curEdge->orient.isAB_AB()) ? "AB_AB" :
           ((curEdge->orient.isAB_BA()) ? "AB_BA" :
            ((curEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))));

#ifdef COMPARE_NEW_OLD
  bool   resnew = isQualityScaffoldMergingEdgeNEW(curEdge, scaffoldA, scaffoldB, si, MIs, minSatisfied, maxDelta);
  bool   resold = isQualityScaffoldMergingEdgeOLD(curEdge, scaffoldA, scaffoldB, si, MIs, minSatisfied, maxDelta);

  if (resnew)  newPASS++;  else  newFAIL++;
  if (resold)  oldPASS++;  else  oldFAIL++;

  fprintf(stderr, "isQualityScaffoldMergingEdge()--  NEW %s (%d/%d) OLD %s (%d/%d)\n",
          resnew ? "pass" : "fail", newPASS, newFAIL,
          resold ? "pass" : "fail", oldPASS, oldFAIL);
#else
  bool   resnew = isQualityScaffoldMergingEdgeNEW(curEdge, scaffoldA, scaffoldB, si, MIs, minSatisfied, maxDelta);

  if (resnew)  newPASS++;  else  newFAIL++;

  fprintf(stderr, "isQualityScaffoldMergingEdge()--  NEW %s (%d/%d)\n",
          resnew ? "pass" : "fail", newPASS, newFAIL);
#endif

  return(resnew);
}





//ExamineUsableSEdges
// This function examines the candidate SEdges in order of weight and decides
// which edges to use for scaffold merges.  We currently do not do interleaved merging,
// but we do try to merge the scaffolds in an order that reduces the need for interleaved
// merging.  A proposed merge can be rejected if:
//     1) It appears to have too long an overlap and is therefor probably an interleaved merge
//     2) One of the scaffold ends involved has already been dedicated to a merge
//

void
ExamineSEdgeForUsability_Interleaved(SEdgeT            *curEdge,
                                     InterleavingSpec  *iSpec,
                                     CIScaffoldT       *scaffoldA,
                                     CIScaffoldT       *scaffoldB);


static
void
ExamineSEdgeForUsability(SEdgeT            *curEdge,
                         InterleavingSpec  *iSpec) {

  if (CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD > curEdge->edgesContributing) {
    fprintf(stderr, "ExamineSEdgeForUsability()-- to few edges %d < %d\n", curEdge->edgesContributing, CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD);
    return;
  }

  CIScaffoldT *scaffoldA = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idA);
  CIScaffoldT *scaffoldB = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idB);

  //PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph, "S ", curEdge, curEdge->idA);

  // We don't want to stick a teeny tiny element in the middle of a gap between two
  // giants.  Often this will later prevent the giants from merging.  So, we
  // prevent it.  We only place stuff such that the size of the element being placed
  // is at least 20% of the size of the implied gap.  So, for a 50k gap, we need a 10k
  // scaffold.  For a 10k gap, we need a 2k scaffold, etc.
  double length_to_dist = 1.0;
  double min_scaffold_length = MIN(scaffoldA->bpLength.mean, scaffoldB->bpLength.mean);

  if (curEdge->distance.mean > 0)
    length_to_dist = min_scaffold_length / curEdge->distance.mean;

  if ((iSpec->checkForTinyScaffolds) &&
      (min_scaffold_length < 5000) &&
      (length_to_dist < 0.20)) {
    fprintf(stderr, "ExamineSEdgeForUsability()-- Scaffolds %d and %d are too short (%.0f and %.0f bp) relative to edge length (%.0f).  Skip.\n",
            curEdge->idA,
            curEdge->idB, 
            scaffoldA->bpLength.mean,
            scaffoldB->bpLength.mean,
            curEdge->distance.mean);
    return;
  }

  //  Makre sure that a weak link doesn't preempt a strong link
#warning thread unsafe
  if (scaffoldA->flags.bits.walkedAlready ||
      scaffoldB->flags.bits.walkedAlready ) {
    scaffoldA->flags.bits.walkedAlready = 1;
    scaffoldB->flags.bits.walkedAlready = 1;
    return;
  }

  if (TouchesMarkedScaffolds(curEdge)) {
    //fprintf(stderr, "* Edge ("F_CID","F_CID",%c) touches marked scaffolds\n",
    //        curEdge->idA, curEdge->idB, curEdge->orient.toLetter());
    return;
  }


  // See if we already want to merge these two scaffolds, but in an opposite orientation
  // For instance, sometimes there will be evidence suggesting a merge AB_BA and a merge
  // BA_BA between a given pair of scaffolds.
  if ((scaffoldA->essentialEdgeB != NULLINDEX &&
       (scaffoldA->essentialEdgeB == scaffoldB->essentialEdgeA ||
        scaffoldA->essentialEdgeB == scaffoldB->essentialEdgeB)) ||
      (scaffoldA->essentialEdgeA != NULLINDEX &&
       (scaffoldA->essentialEdgeA == scaffoldB->essentialEdgeA ||
        scaffoldA->essentialEdgeA == scaffoldB->essentialEdgeB))) {
    //fprintf(stderr, "* We're already trying to merge scaffold ("F_CID","F_CID") ...back off!\n",
    //        scaffoldA->id, scaffoldB->id);
    return;
  }

  ExamineSEdgeForUsability_Interleaved(curEdge, iSpec, scaffoldA, scaffoldB);
}





static
bool
ExamineUsableSEdges(vector<SEdgeT *>  &sEdges,
                    double             weightScale,
                    InterleavingSpec  *iSpec) {
  double minWeightThreshold = 0;
  int32  maxWeightEdge      = 0;

  if (sEdges.size() == 0)
    return(false);

  for (uint32 i=0; i<sEdges.size(); i++) {
    if (isBadScaffoldMergeEdge(sEdges[i], iSpec->badSEdges))
      continue;

    if (maxWeightEdge >= sEdges[i]->edgesContributing)
      continue;

    fprintf(stderr, "ExamineUsableSEdges()- maxWeightEdge from "F_S32" to "F_S32" at idx "F_S32" out of "F_SIZE_T"\n",
            maxWeightEdge, sEdges[i]->edgesContributing, i, sEdges.size());
    maxWeightEdge = sEdges[i]->edgesContributing;
  }

  //  We now recompute the min weight allowed to merge each iteration.  Previous to this commit
  //  (1.66) min was computed once at the start and slowly decremented each iteration.

  minWeightThreshold = MAX(maxWeightEdge * weightScale, GlobalData->minWeightToMerge);

  fprintf(stderr, "* Considering edges with weight >= %.2f (maxWeightEdge %d weightScale %.4f)\n",
          minWeightThreshold,
          maxWeightEdge,
          weightScale);

  //  Examine the edges

  int32    edgeListLen  = sEdges.size();
  int32    edgeListStep = edgeListLen / 100;

  //#pragma omp parallel for schedule(dynamic, edgeListStep)
  for (int i=0; i<edgeListLen; i++) {
    if (sEdges[i]->edgesContributing < minWeightThreshold)
      continue;

    ExamineSEdgeForUsability(sEdges[i], iSpec);
  }

  return(minWeightThreshold > GlobalData->minWeightToMerge);
}



////////////////////////////////////////



static
bool
CompareSEdgesContributing(const SEdgeT *s1, const SEdgeT *s2) {
  int32  n2 = (s2->edgesContributing - (isOverlapEdge(s2) ? 1 : 0));
  int32  n1 = (s1->edgesContributing - (isOverlapEdge(s1) ? 1 : 0));

  if (n2 < n1)
    //  Edge1 higher weight
    return(true);

  if (n1 < n2)
    //  Edge1 not higher weight
    return(false);

  //  Otherwise, they're the same weight.  Break ties using the insert size.

  double  d1 = fabs(PREFERRED_GAP_SIZE - s1->distance.mean) + sqrt(MAX(1., s1->distance.variance));
  double  d2 = fabs(PREFERRED_GAP_SIZE - s2->distance.mean) + sqrt(MAX(1., s2->distance.variance));

  //  Edge1 closer to truth
  return(d1 < d2);
}




//  Check the edges off the other end to see if there is a stronger edge from there.
//
static
bool
OtherEndHasStrongerEdgeToSameScaffold(CDS_CID_t scfIID, SEdgeT * curSEdge) {
  int       otherEnd;
  CDS_CID_t otherScaffoldID;

  if (scfIID == curSEdge->idA) {
    otherEnd        = (curSEdge->orient.isAB_AB() || curSEdge->orient.isAB_BA()) ? A_END : B_END;
    otherScaffoldID = curSEdge->idB;
  } else {
    otherEnd        = (curSEdge->orient.isAB_BA() || curSEdge->orient.isBA_BA()) ? A_END : B_END;
    otherScaffoldID = curSEdge->idA;
  }

  GraphEdgeIterator SEdges(ScaffoldGraph->ScaffoldGraph, scfIID, otherEnd, ALL_EDGES);
  SEdgeT           *SEdge;

  while ((SEdge = SEdges.nextMerged()) != NULL)
    if ((SEdge->idA == otherScaffoldID || SEdge->idB == otherScaffoldID) &&
        (SEdge->edgesContributing > curSEdge->edgesContributing))
      return(true);

  return(false);
}



static
void
BuildUsableSEdges(vector<SEdgeT *>   &sEdges,
                  vector<SEdgeT *>   &oEdges,
                  int32               verbose) {

  //  Reset edges
  for (int i=0; i<GetNumGraphEdges(ScaffoldGraph->ScaffoldGraph); i++) {
    EdgeCGW_T * edge = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, i);

    edge->quality                    = 1.0;
    edge->flags.bits.isBogus         = 0;
    edge->flags.bits.isProbablyBogus = 0;

    if (edge->flags.bits.MeanChangedByWalking == TRUE) {
      edge->flags.bits.MeanChangedByWalking = FALSE;
      edge->distance.mean                   = edge->minDistance;
    }
  }

  //  Reset scaffolds
  for (int i=0; i<GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); i++) {
    CIScaffoldT * scaffold = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, i);

    if (isDeadCIScaffoldT(scaffold) || (scaffold->type != REAL_SCAFFOLD))
      continue;

    scaffold->numEssentialA  = 0;
    scaffold->numEssentialB  = 0;
    scaffold->essentialEdgeA = NULLINDEX;
    scaffold->essentialEdgeB = NULLINDEX;
    scaffold->setID          = NULLINDEX;
    scaffold->flags.bits.smoothSeenAlready = 0;
    scaffold->flags.bits.walkedAlready     = 0;
  }


  for (int i=0; i<GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); i++) {
    CIScaffoldT  *thisScaffold = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, i);
    CIScaffoldT  *thatScaffold = NULL;

    if (isDeadCIScaffoldT(thisScaffold) || (thisScaffold->type != REAL_SCAFFOLD))
      continue;

    GraphEdgeIterator SEdges(ScaffoldGraph->ScaffoldGraph, thisScaffold->id, ALL_END, ALL_EDGES);
    SEdgeT           *SEdge;

    while ((SEdge = SEdges.nextMerged()) != NULL) {
      if (SEdge->flags.bits.isBogus)
        // This edge has already been visited by the recursion
        continue;

      if (SEdge->flags.bits.isDeleted)
        //  Deleted edge?  Really?  We delete edges?
        continue;

      if ((SEdge->edgesContributing - (isOverlapEdge(SEdge) ? 1 : 0)) < CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD)
        //  Edge weight to weak.
        continue;

      assert(SEdge->idA != NULLINDEX);
      assert(SEdge->idB != NULLINDEX);

      if (SEdge->idA != thisScaffold->id)
        //  Not canonical edge
        continue;

      thatScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, SEdge->idB);

      if (isDeadCIScaffoldT(thatScaffold) || (thatScaffold->type != REAL_SCAFFOLD))
        continue;

      if (thatScaffold->flags.bits.smoothSeenAlready)
        continue;

      if (TouchesMarkedScaffolds(SEdge))
        break;

      if (OtherEndHasStrongerEdgeToSameScaffold(thisScaffold->id, SEdge))
        continue;
      
      sEdges.push_back(SEdge);
    }
  }

  sort(sEdges.begin(), sEdges.end(), CompareSEdgesContributing);
}











//static
void
MergeScaffoldsAggressive(ScaffoldGraphT *graph, char *logicalcheckpointnumber, int verbose) {
  InterleavingSpec iSpec;

  ScaffoldSanity(ScaffoldGraph);

  fprintf(stderr, "MergeScaffoldsAggressive()-- begins.  minWeightToMerge="F_S32"\n", GlobalData->minWeightToMerge);

  iSpec.sai                    = CreateScaffoldAlignmentInterface();
  iSpec.contigNow              = TRUE;
  iSpec.checkForTinyScaffolds  = FALSE;
  iSpec.checkAbutting          = TRUE;
  iSpec.minSatisfied           = 0.985;  //  0.985 default
  iSpec.maxDelta               = -1;     //  0.005
  iSpec.MIs                    = CreateVA_MateInstrumenterP(GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));
  iSpec.badSEdges              = CreateChunkOverlapper();

  //  Be conservative, and rebuild the full edge set once.

  {
    vector<CDS_CID_t>  rawEdges;

    BuildSEdges(rawEdges, TRUE);
    MergeAllGraphEdges(graph->ScaffoldGraph, rawEdges, TRUE, FALSE);
  }



  {
    vector<SEdgeT *>      sEdges;
    vector<SEdgeT *>      oEdges;
    set<EdgeCGWLabel_T>   bEdges;  //  Bad edges - new scaffold is not connected if these are used

    int32             iterations         = 0;
    double            weightScaleInit    = 0.75;
    double            weightScale        = weightScaleInit;
    time_t            lastCkpTime        = time(0) - 90 * 60;


    // loop until nothing gets merged
    while (1) {
      time_t t = time(0);

      //  Checkpoint periodically - every two hours seems nice!  The
      //  first checkpoint is done after 30 minutes of work here,
      //  though.
      //
      if (t - lastCkpTime > 120 * 60) {
        char  where[1024];

        sprintf(where, "after MergeScaffoldsAggressive iteration %d", iterations);
        CheckpointScaffoldGraph(logicalcheckpointnumber, where);
        lastCkpTime = t;
      }

      if (GetNumGraphEdges(graph->ScaffoldGraph) == 0) {
        fprintf(stderr, "MergeScaffoldsAggressive()-- No additional scaffold merging is possible.\n");
        break;
      }

      //  Build sEdges for merging

      sEdges.clear();
      oEdges.clear();
      BuildUsableSEdges(sEdges, oEdges, verbose);

      bool  moreWork = ExamineUsableSEdges(sEdges, weightScale, &iSpec);

      //  Eventyally, we want to pass in matePairTestResult for scoring of which mate pair tests
      //  resulted in successful merges.  At the moment, we have no inexpensive way of traching
      //  SEdges from there to here.
      if (MergeScaffolds(&iSpec, bEdges)) {
        fprintf(stderr, "MergeScaffoldsAggressive()-- iter %d -- continue because we merged scaffolds.\n",
                iterations);

#ifdef TRACK_MATE_PAIR_TEST
        fprintf(stderr, "The following edges successfully merged.\n");
        for (map<uint32, char*>::iterator it = matePairTestResult.begin(); it != matePairTestResult.end(); it++) {
          if (id->second.mergeAccepted == true)
            it->second.print();
        }

        fprintf(stderr, "The following edges failed to merge.\n");
        for (map<uint32, char*>::iterator it = matePairTestResult.begin(); it != matePairTestResult.end(); it++) {
          if (id->second.mergeAccepted == false)
            it->second.print();
        }

        matePairTestResult.clear();
#endif

#if 0
        //  Left in, just in case.  MergeScaffolds() now rebuilds only those scaffold graph edges
        //  that have changed.  If you enable this, ALL scaffold graph edges will be rebuilt.
        {
          vector<CDS_CID_t>  rawEdges;

          BuildSEdges(rawEdges, TRUE);
          MergeAllGraphEdges(graph->ScaffoldGraph, rawEdges, TRUE, FALSE);
        }
#endif

        //  Reset to original weight scaling - this is a bit too aggressive, but seems to be safe.  We do not
        //  end up doing tons and tons of work on failed edges.
        weightScale = weightScaleInit;

      } else if (moreWork == true) {
        fprintf(stderr, "MergeScaffoldsAggressive()-- iter %d -- decrease weight scaling from %.4f to %.4f and continue.\n",
                iterations, weightScale, weightScale * 0.8);

        //  Do we need to clean up the edges/scaffolds here?  Can we jump right back to ExamineUsableSEdges()?

        weightScale *= 0.8;  //  Drop weight scaling by a bit

      } else {
        fprintf(stderr, "MergeScaffoldsAggressive()-- iter %d -- no additional scaffold merging is possible.\n",
                iterations);
        break;
      }

      iterations++;
    }
  }

  DeleteScaffoldAlignmentInterface(iSpec.sai);

  for (int i=0; i<GetNumVA_PtrT(iSpec.MIs); i++)
    safe_free(*(MateInstrumenter **)GetVA_PtrT(iSpec.MIs, i));
  DeleteVA_PtrT(iSpec.MIs);

  DestroyChunkOverlapper(iSpec.badSEdges);

  ScaffoldSanity(ScaffoldGraph);
}

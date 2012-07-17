
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
static char *rcsid = "$Id: CIScaffoldT_Merge_CGW.c,v 1.69 2012-07-17 15:19:24 brianwalenz Exp $";

//
//  The ONLY exportable function here is MergeScaffoldsAggressive.
//

#undef DEBUG_MERGE_EDGE_INVERT


#define OTHER_END_CHECK
#undef  GENERAL_STRONGER_CHECK

#undef  ALLOW_NON_INTERLEAVED_MERGING

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

#define MAX_SCAFFOLD_GAP_OVERLAP 1000

#define CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD MIN_EDGES

#define EDGE_QUANTA 5.0
#define OVERLAP_QUANTA -10000.

#define EDGE_STRENGTH_FACTOR  MIN_EDGES

#define MAX_SLOP_IN_STD 3.5

#define EDGE_WEIGHT_FACTOR  MIN_EDGES



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




static
int
isLargeOverlap(CIEdgeT *curEdge) {
  double maxGap   = curEdge->distance.mean + 3.0 * sqrt(curEdge->distance.variance);
  double minGap   = curEdge->distance.mean - 3.0 * sqrt(curEdge->distance.variance);
  int32  numEdges = curEdge->edgesContributing;

  double numOverlapQuanta = minGap   / OVERLAP_QUANTA;
  double numEdgeQuanta    = numEdges / EDGE_QUANTA;

  fprintf(stderr, "isLargeOverlapSize()--  edges %d gap %f %f ovlQuanta %f (%f) edgeQuanta %f (%f)\n",
          numEdges,
          minGap, maxGap, 
          numOverlapQuanta, OVERLAP_QUANTA,
          numEdgeQuanta,    EDGE_QUANTA);

  if (maxGap > OVERLAP_QUANTA)
    return FALSE;

  return(numEdgeQuanta < numOverlapQuanta);
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
  bool  failsMinimum;
  bool  failsToGetHappier1;
  bool  failsToGetHappier2;

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

  uint32  mBeforeGood = A.numEcstatic + B.numEcstatic;
  uint32  mAfterGood  = P.numEcstatic;

  uint32  mBefore     = A.numMateInternal + A.numMateExternal + B.numMateInternal + B.numMateExternal;
  uint32  mAfter      = P.numMateInternal + P.numMateExternal;

  //  Mates that fall into gaps are neither good nor bad.  They just don't exist as far as we are concerned.

  mBeforeGood -= A.numGap + B.numGap;
  mBefore     -= A.numGap + B.numGap;

  mAfterGood  -= P.numGap;
  mAfter      -= P.numGap;

  //  Mates before that fall into contigs, but aren't present, are still bad.

  uint32  mBeforeBad = A.numDejected + B.numDejected;
  uint32  mAfterBad  = P.numDejected;

  double fractMatesHappyBefore = (mAfter > 0) ? ((double)mBeforeGood / mAfter) : 0.0;
  double fractMatesHappyAfter  = (mAfter > 0) ? ((double)mAfterGood  / mAfter) : 0.0;

  fprintf(stderr, "isQualityScaffoldMergingEdge()--   before: %.3f satisfied (%d/%d good/bad mates)  after: %.3f satisfied (%d/%d good/bad mates)\n",
          fractMatesHappyBefore, mBeforeGood, mBeforeBad,
          fractMatesHappyAfter,  mAfterGood,  mAfterBad);

  double badGoodRatio      = 1.0;

  //  NOTE, still need the metagenomic test

  if (mAfterGood > mBeforeGood)
    badGoodRatio = (double)(mAfterBad - mBeforeBad) / (double)(mAfterGood - mBeforeGood);

  failsMinimum       = (fractMatesHappyAfter < minSatisfied);
  failsToGetHappier1 = (fractMatesHappyAfter < fractMatesHappyBefore);
  failsToGetHappier2 = (mAfterGood < mBeforeGood) || (badGoodRatio > MAX_FRAC_BAD_TO_GOOD);

  if (failsMinimum && failsToGetHappier1 && failsToGetHappier2) {
    fprintf(stderr, "isQualityScaffoldMergingEdge()--   not happy enough to merge (%.3f < %.3f) && (%.3f < %.3f) && ((%d < %d) || (%0.3f > %.3f))\n",
            fractMatesHappyAfter, minSatisfied,
            fractMatesHappyAfter, fractMatesHappyBefore,
            mAfterGood, mBeforeGood, badGoodRatio, MAX_FRAC_BAD_TO_GOOD);
    return(FALSE);
  } else {
    fprintf(stderr, "isQualityScaffoldMergingEdge()--   ARE happy enough to merge (%.3f > %.3f) || (%.3f > %.3f) || ((%d > %d) && (%0.3f < %.3f))\n",
            fractMatesHappyAfter, minSatisfied,
            fractMatesHappyAfter, fractMatesHappyBefore,
            mAfterGood, mBeforeGood, badGoodRatio, MAX_FRAC_BAD_TO_GOOD);
    return(TRUE);
  }
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

  bool   resnew = isQualityScaffoldMergingEdgeNEW(curEdge, scaffoldA, scaffoldB, si, MIs, minSatisfied, maxDelta);
  bool   resold = isQualityScaffoldMergingEdgeOLD(curEdge, scaffoldA, scaffoldB, si, MIs, minSatisfied, maxDelta);

  if (resnew)  newPASS++;  else  newFAIL++;
  if (resold)  oldPASS++;  else  oldFAIL++;

  fprintf(stderr, "isQualityScaffoldMergingEdge()--  NEW %s (%d/%d) OLD %s (%d/%d)\n",
          resnew ? "pass" : "fail", newPASS, newFAIL,
          resold ? "pass" : "fail", oldPASS, oldFAIL);

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
ExamineSEdgeForUsability_NonInterleaved(SEdgeT * curEdge, InterleavingSpec * iSpec,
                                        CIScaffoldT *scaffoldA,
                                        CIScaffoldT *scaffoldB,
                                        double minMergeDistance,
                                        double maxMergeDistance,
                                        int32 mayOverlap,
                                        int32 mustOverlap);


void
ExamineSEdgeForUsability_Interleaved(SEdgeT * curEdge, InterleavingSpec * iSpec,
                                     CIScaffoldT *scaffoldA,
                                     CIScaffoldT *scaffoldB,
                                     double minMergeDistance,
                                     double maxMergeDistance,
                                     int32 mayOverlap,
                                     int32 mustOverlap);


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

  if (scaffoldA->flags.bits.walkedAlready ||
      scaffoldB->flags.bits.walkedAlready ) {
    //  We want to make sure that a week link doesn't preempt a strong link
    scaffoldA->flags.bits.walkedAlready = scaffoldB->flags.bits.walkedAlready = 1;
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

  double mergeDistance    = curEdge->distance.mean;
  double minMergeDistance = mergeDistance - MAX_SLOP_IN_STD * sqrt(curEdge->distance.variance);
  double maxMergeDistance = mergeDistance + MAX_SLOP_IN_STD * sqrt(curEdge->distance.variance);

  //fprintf(stderr, "* curEdge mergeDistance = (%g,%g) min:%g max:%g\n",
  //        mergeDistance, curEdge->distance.variance,
  //        minMergeDistance, maxMergeDistance);

  // Look for an overlap

  int32 mayOverlap  = ((minMergeDistance < CGW_MISSED_OVERLAP) && (maxMergeDistance > CGW_MISSED_OVERLAP));
  int32 mustOverlap = ((minMergeDistance < CGW_MISSED_OVERLAP) && (maxMergeDistance < CGW_MISSED_OVERLAP));

  // If it is a really heavy edge, treat like a may overlap edge
  if ((mustOverlap) && (curEdge->edgesContributing > EDGE_QUANTA)) {
    mustOverlap = 0;
    mayOverlap  = 1;
  }

  ExamineSEdgeForUsability_Interleaved(curEdge, iSpec, scaffoldA, scaffoldB, minMergeDistance, maxMergeDistance, mayOverlap, mustOverlap);
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

  minWeightThreshold = MAX(maxWeightEdge * weightScale, EDGE_WEIGHT_FACTOR);

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

  return(minWeightThreshold > EDGE_WEIGHT_FACTOR);
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


#ifdef GENERAL_STRONGER_CHECK
/*
  it would be better to do this during scaffold edge creation
*/

static
int
ThereIsAStrongerEdgeToSameScaffoldReport(SEdgeT *curSEdge, int32 retVal) {
  if (retVal == 0)
    fprintf(stderr,
            "SCF MERGE CONFLICT: "F_CID","F_CID"  %s  %dbp  %dvar  %dec\n",
            curSEdge->idA, curSEdge->idB,
            ((curSEdge->orient.isAB_AB()) ? "AB_AB" :
             ((curSEdge->orient.isAB_BA()) ? "AB_BA" :
              ((curSEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
            (int) curSEdge->distance.mean,
            (int) curSEdge->distance.variance,
            curSEdge->edgesContributing);

  fprintf(stderr, "\t"F_CID","F_CID", %s, %dbp  %dvar  %dec\n",
          sEdge->idA, sEdge->idB,
          ((sEdge->orient.isAB_AB()) ? "AB_AB" :
           ((sEdge->orient.isAB_BA()) ? "AB_BA" :
            ((sEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
          (int) sEdge->distance.mean,
          (int) curSEdge->distance.variance,
          sEdge->edgesContributing);

  return(retVal + 1);
}

static
int
ThereIsAStrongerEdgeToSameScaffold(CDS_CID_t scfIID, SEdgeT * curSEdge) {
  SEdgeTIterator SEdges;
  SEdgeT * sEdge;
  int32 orientValue;
  CDS_CID_t otherScaffoldID;
  int32 retVal = 0;

  otherScaffoldID = (scfIID == curSEdge->idA) ? curSEdge->idB : curSEdge->idA;
  orientValue =
    (curSEdge->orient.isAB_AB() || curSEdge->orient.isBA_BA()) ? 1 : -1;

  /*
    iterate over otherEnd's merged edges to see if there is a
    stronger one to the other scaffold
  */
  InitSEdgeTIterator(ScaffoldGraph, scfIID,
                     FALSE, FALSE, ALL_END, FALSE, &SEdges);
  while ((sEdge = NextSEdgeTIterator(&SEdges)) != NULL) {

    if (sEdge->idA == otherScaffoldID ||
        (sEdge->idB == otherScaffoldID && sEdge != curSEdge)) {
      int newOrientValue =
        (sEdge->orient.isAB_AB() || sEdge->orient.isBA_BA()) ? 1 : -1;

      /*
        Does the pair of edges agree (shift) or disagree (reversal)?
        If reversal, assume greater weight one is correct
        If shift,
        if both are positive, take one with greater weight
        if one is negative, take the postive one
        if both are negative, take the shortest one (smallest overlap)
      */
      if (orientValue != newOrientValue) {
        // reversal
        if (sEdge->edgesContributing > curSEdge->edgesContributing) {
          retVal = ThereIsAStrongerEdgeToSameScaffoldReport(curSEdge, retVal);
        }
      } else {
        // shift
        if (curSEdge->distance.mean > 0) {
          if (sEdge->distance.mean > 0) {
            // both are positive, prefer stronger
            if (sEdge->edgesContributing > curSEdge->edgesContributing) {
              retVal = ThereIsAStrongerEdgeToSameScaffoldReport(curSEdge, retVal);
            }
          } else {
            /*
              curSEdge is positive, sEdge is negative
              prefer curSEdge unless sEdge is much stronger
            */
            if (sEdge->edgesContributing >
                EDGE_STRENGTH_FACTOR * curSEdge->edgesContributing) {
              retVal = ThereIsAStrongerEdgeToSameScaffoldReport(curSEdge, retVal);
            }
          }
        } else if (sEdge->distance.mean > 0) {
          /*
            sEdge is positive, curSEdge is negative
            prefer sEdge unless curSEdge is much stronger
          */
          if (curSEdge->edgesContributing <
              EDGE_STRENGTH_FACTOR * sEdge->edgesContributing) {
            retVal = ThereIsAStrongerEdgeToSameScaffoldReport(curSEdge, retVal);
          }
        } else {
          /*
            both negative
            prefer shorter overlap unless longer is much stronger
            prefer sEdge if much stronger, or shorter & strong enough
          */
          if (sEdge->edgesContributing >
              EDGE_STRENGTH_FACTOR * curSEdge->edgesContributing ||
              (sEdge->distance.mean > curSEdge->distance.mean &&
               EDGE_STRENGTH_FACTOR * sEdge->edgesContributing >
               curSEdge->edgesContributing)) {
            retVal = ThereIsAStrongerEdgeToSameScaffoldReport(curSEdge, retVal);
          }
        }
      }
    }
  }
  return retVal;
}
#endif





#ifdef OTHER_END_CHECK
static
int
OtherEndHasStrongerEdgeToSameScaffold(CDS_CID_t scfIID, SEdgeT * curSEdge) {
  SEdgeTIterator SEdges;
  SEdgeT * sEdge;
  int otherEnd;
  CDS_CID_t otherScaffoldID;

  if (scfIID == curSEdge->idA) {
    otherEnd = (curSEdge->orient.isAB_AB() || curSEdge->orient.isAB_BA()) ? A_END : B_END;
    otherScaffoldID = curSEdge->idB;
  } else {
    otherEnd = (curSEdge->orient.isAB_BA() || curSEdge->orient.isBA_BA()) ? A_END : B_END;
    otherScaffoldID = curSEdge->idA;
  }

  //  iterate over otherEnd's merged edges to see if there is a
  //  stronger one to curSEdge->idB

  InitSEdgeTIterator(ScaffoldGraph, scfIID, FALSE, FALSE, otherEnd, FALSE, &SEdges);

  while ((sEdge = NextSEdgeTIterator(&SEdges)) != NULL)
    if ((sEdge->idA == otherScaffoldID || sEdge->idB == otherScaffoldID) &&
        (sEdge->edgesContributing > curSEdge->edgesContributing))
      return 1;

  return 0;
}
#endif





// Find all merge candidates incident on scaffoldA
// Returns TRUE if marked edges are encountered
//
static
int
FindAllMergeCandidates(vector<SEdgeT *>   &sEdges,
                       vector<SEdgeT *>   &oEdges,
                       CIScaffoldT        *fromScaffold,
                       int                 fromEnd,
                       CIScaffoldT        *ignoreToScaffold,
                       int                 canonicalOnly,
                       int                 minWeight,
                       int                 verbose) {
  SEdgeTIterator  SEdges;
  SEdgeT         *curSEdge;
  CIScaffoldT    *otherScaffold;
  CDS_CID_t       otherScaffoldID;

  InitSEdgeTIterator(ScaffoldGraph, fromScaffold->id, FALSE, FALSE, fromEnd, FALSE, &SEdges);

  while ((curSEdge = NextSEdgeTIterator(&SEdges)) != NULL) {
    if (curSEdge->flags.bits.isBogus)
      // This edge has already been visited by the recursion
      continue;

    if ( curSEdge->idA != fromScaffold->id) {
      if (canonicalOnly)
        continue;
      otherScaffoldID = curSEdge->idA;
    } else {
      otherScaffoldID = curSEdge->idB;
    }

    otherScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, otherScaffoldID);

    if (otherScaffold->flags.bits.smoothSeenAlready)
      continue;

    if (otherScaffold == ignoreToScaffold)
      continue;

    if (TouchesMarkedScaffolds(curSEdge))
      return TRUE;

#ifdef OTHER_END_CHECK
    if (OtherEndHasStrongerEdgeToSameScaffold(fromScaffold->id, curSEdge))
      continue;
#endif

#ifdef GENERAL_STRONGER_CHECK
    if (ThereIsAStrongerEdgeToSameScaffold(fromScaffold->id, curSEdge))
      continue;
#endif

    if (curSEdge->flags.bits.isDeleted ||
        isDeadCIScaffoldT(otherScaffold) ||
        otherScaffold->type != REAL_SCAFFOLD)
      continue;

    assert((curSEdge->idA != NULLINDEX) && (curSEdge->idB != NULLINDEX));

    if ((curSEdge->edgesContributing - (isOverlapEdge(curSEdge) ? 1 : 0)) < minWeight)
      continue;

    sEdges.push_back(curSEdge);
  }

  return FALSE;
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
    CIScaffoldT * scaffold = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, i);

    if (isDeadCIScaffoldT(scaffold) || (scaffold->type != REAL_SCAFFOLD))
      continue;

    FindAllMergeCandidates(sEdges, oEdges,
                           scaffold,
                           ALL_END, NULL,
                           TRUE,
                           CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD,
                           verbose);
  }

  sort(sEdges.begin(), sEdges.end(), CompareSEdgesContributing);
}




////////////////////////////////////////








//static
void
MergeScaffoldsAggressive(ScaffoldGraphT *graph, char *logicalcheckpointnumber, int verbose) {
  InterleavingSpec iSpec;

  if (verbose)
    fprintf(stderr, "Starting MergeScaffoldsAggressive\n");

  CheckCIScaffoldTs(ScaffoldGraph);
  CheckCIScaffoldTLengths(ScaffoldGraph);
  fprintf(stderr, "* Successfully passed checks at beginning of scaffold merging\n");

  iSpec.sai                    = CreateScaffoldAlignmentInterface();
  iSpec.contigNow              = TRUE;
  iSpec.checkForTinyScaffolds  = FALSE;
  iSpec.checkAbutting          = TRUE;
  iSpec.minSatisfied           = 0.985;  //  0.985 default
  iSpec.maxDelta               = -1;     //  0.005
  iSpec.MIs                    = CreateVA_MateInstrumenterP(GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));
  iSpec.badSEdges              = CreateChunkOverlapper();

  {
    vector<SEdgeT *>  sEdges;
    vector<SEdgeT *>  oEdges;

    int32             iterations         = 0;
    double            weightScaleInit    = 0.75;
    double            weightScale        = weightScaleInit;
    time_t            lastCkpTime        = time(0) - 90 * 60;

    //  Create a scaffold edge for every inter-scaffold contig edge, then merge compatible ones
    BuildSEdges(graph, TRUE, TRUE);
    MergeAllGraphEdges(graph->ScaffoldGraph, TRUE, FALSE);

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

      if (MergeScaffolds(&iSpec, verbose)) {
        fprintf(stderr, "MergeScaffoldsAggressive()-- iter %d -- continue because we merged scaffolds.\n",
                iterations);

        CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);
        BuildSEdges(graph, TRUE, TRUE);
        MergeAllGraphEdges(graph->ScaffoldGraph, TRUE, FALSE);

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

  CheckCIScaffoldTs(ScaffoldGraph);
}

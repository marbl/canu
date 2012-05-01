
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
static char *rcsid = "$Id: CIScaffoldT_Merge_Interleaved.c,v 1.1 2012-05-01 04:18:19 brianwalenz Exp $";

//  These are private functions used in CIScaffoldT_Merge_CGW.c.
//  These functions are specialized for interleaved merging.

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Instrument_CGW.h"
#include "ChiSquareTest_CGW.h"
#include "InterleavedMerging.h"

#include "CIScaffoldT_MergeScaffolds.h"

#include <vector>
#include <algorithm>
using namespace std;



int
isQualityScaffoldMergingEdge(SEdgeT                     *curEdge,
                             CIScaffoldT                *scaffoldA,
                             CIScaffoldT                *scaffoldB,
                             ScaffoldInstrumenter       *si,
                             VA_TYPE(MateInstrumenterP) *MIs,
                             double                       minSatisfied,
                             double                       maxDelta);



//  After noticing with interleaved scaffold merging that there can be some rather high-weight
//  negative edges that are slightly beyond the limit of abutting by the chi-squared check, we
//  wanted to let them be abutted.
//
//  Two conditions must be satisfied:
//
//  1. The number of stddev's from a -20 edge per unit of edge weight is less than .5
//
//  2. The length of the edge is a small fraction of the shorter of the two scaffolds involved.
//
//  The specific numbers are based on an examination of edges with weight >= 25 in the rat assembly
//  12/16/2002.
//
//  After noticing that this was still not aggressive enough during the Macaque assembly we added
//  another condition. We also hope that this will allow ECR to merge the abbutted contigs.
//
//  3. If the overlap is < 2kbp and the overlap is less than 1/2 of the shorter scaffold then abut.
//
//  Granger 8/22/05.
//
#define STDDEVS_PER_WEIGHT_THRESHOLD                 0.5
#define EDGE_PER_MIN_SCAFFOLD_LENGTH_THRESHOLD       0.002

#define MAX_OVERLAP_TO_ABUT                          2000
#define MAX_PERC_SCAFFOLD_LEN                        0.5



static
void
SaveEdgeMeanForLater(SEdgeT * edge) {
  if (edge->flags.bits.MeanChangedByWalking == FALSE) {
    edge->flags.bits.MeanChangedByWalking = TRUE;
    edge->minDistance = edge->distance.mean;
  }
}


static
void
SaveBadScaffoldMergeEdge(SEdgeT * edge,
                         ChunkOverlapperT * overlapper) {
  ChunkOverlapCheckT * lookup;
  ChunkOverlapSpecT spec;
  double delta = sqrt(edge->distance.variance) * 3.;

  InitCanonicalOverlapSpec(edge->idA, edge->idB, edge->orient, &spec);
  if ((lookup = LookupCanonicalOverlap(overlapper, &spec)) == NULL) {
    ChunkOverlapCheckT olap;
    memset(&olap, 0, sizeof(ChunkOverlapCheckT));
    // create
    FillChunkOverlapWithEdge(edge, &olap);
    olap.minOverlap = (int32)(-edge->distance.mean - delta);
    olap.maxOverlap = (int32)(-edge->distance.mean + delta);
    InsertChunkOverlap(overlapper, &olap);
  } else {
    // update
    lookup->minOverlap = (int32)MIN(lookup->minOverlap, -edge->distance.mean - delta);
    lookup->maxOverlap = (int32)MAX(lookup->maxOverlap, -edge->distance.mean + delta);
  }
}


//static
int
isBadScaffoldMergeEdge(SEdgeT * edge, ChunkOverlapperT * overlapper) {
  ChunkOverlapCheckT * lookup;
  ChunkOverlapSpecT spec;
  double delta = sqrt(edge->distance.variance) * 3.;
  int32 minOverlap = (int32)(-edge->distance.mean - delta);
  int32 maxOverlap = (int32)(-edge->distance.mean + delta);

  InitCanonicalOverlapSpec(edge->idA, edge->idB, edge->orient, &spec);
  if ((lookup = LookupCanonicalOverlap(overlapper, &spec)) != NULL) {
    if (minOverlap >= lookup->minOverlap - 10 &&
        maxOverlap <= lookup->maxOverlap + 10)
      return TRUE;
  }
  return FALSE;
}



static
void
MarkScaffoldsForMerging(SEdgeT *curEdge) {
  CIScaffoldT *scaffoldA = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idA);
  CIScaffoldT *scaffoldB = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idB);

  CDS_CID_t edgeIndex = (CDS_CID_t)GetVAIndex_SEdgeT(ScaffoldGraph->ScaffoldGraph->edges, curEdge);

  //  Note the differences in what is an "A" and what is a "B" below.

  assert(curEdge->orient.isUnknown() == false);

  if (curEdge->orient.isAB_AB()) {
    assert(scaffoldA->essentialEdgeB == NULLINDEX);
    assert(scaffoldB->essentialEdgeA == NULLINDEX);
    scaffoldA->essentialEdgeB = edgeIndex;
    scaffoldB->essentialEdgeA = edgeIndex;
    scaffoldA->numEssentialB = 1;
    scaffoldB->numEssentialA = 1;
  }
  if (curEdge->orient.isAB_BA()) {
    assert(scaffoldA->essentialEdgeB == NULLINDEX);
    assert(scaffoldB->essentialEdgeB == NULLINDEX);
    scaffoldA->essentialEdgeB = edgeIndex;
    scaffoldB->essentialEdgeB = edgeIndex;
    scaffoldA->numEssentialB = 1;
    scaffoldB->numEssentialB = 1;
  }
  if (curEdge->orient.isBA_AB()) {
    assert(scaffoldA->essentialEdgeA == NULLINDEX);
    assert(scaffoldB->essentialEdgeA == NULLINDEX);
    scaffoldA->essentialEdgeA = edgeIndex;
    scaffoldB->essentialEdgeA = edgeIndex;
    scaffoldA->numEssentialA = 1;
    scaffoldB->numEssentialA = 1;
  }
  if (curEdge->orient.isBA_BA()) {
    assert(scaffoldA->essentialEdgeA == NULLINDEX);
    assert(scaffoldB->essentialEdgeB == NULLINDEX);
    scaffoldA->essentialEdgeA = edgeIndex;
    scaffoldB->essentialEdgeB = edgeIndex;
    scaffoldA->numEssentialA = 1;
    scaffoldB->numEssentialB = 1;
  }

  // ignore these two scaffolds for the rest of this iteration
  scaffoldA->flags.bits.walkedAlready  = 1;
  scaffoldB->flags.bits.walkedAlready  = 1;
}

//   UNUSED?
//static
//void
//MarkScaffoldForNotMerging(CIScaffoldT *scaffoldA) {
//  scaffoldA->numEssentialB = 1;
//  scaffoldA->numEssentialA = 1;
//}







//  The BEndCI listed in the scaffold may be contained; find the true one.
static
ChunkInstanceT *
GetTrueBEndCI(ScaffoldGraphT *graph, CIScaffoldT *scaffold) {
  CIScaffoldTIterator ci;
  ChunkInstanceT     *iti;
  ChunkInstanceT     *ret = GetGraphNode(graph->ContigGraph, scaffold->info.Scaffold.BEndCI);
  double              retmax = MAX(ret->offsetAEnd.mean, ret->offsetBEnd.mean);
  double              retmin = MIN(ret->offsetAEnd.mean, ret->offsetBEnd.mean);

  InitCIScaffoldTIterator(graph, scaffold, FALSE, FALSE, &ci);
  while ((iti = NextCIScaffoldTIterator(&ci)) != NULL) {
    double  mymax = MAX(iti->offsetAEnd.mean, iti->offsetBEnd.mean);

    if (retmax < mymax) {
      ret    = iti;
      retmax = MAX(ret->offsetAEnd.mean, ret->offsetBEnd.mean);
      retmin = MIN(ret->offsetAEnd.mean, ret->offsetBEnd.mean);
    }

    if (retmin > mymax)
      break;
  }

  return(ret);
}


// TranslateScaffoldOverlapToContigOverlap
//       Prepare the groundwork for computing an overlap between the appropriate extremal
//       contigs in a pair of scaffolds.
//
static
void
TranslateScaffoldOverlapToContigOverlap(CIScaffoldT     *scaffoldA,
                                        CIScaffoldT     *scaffoldB,
                                        PairOrient       scaffoldEdgeOrient,
                                        NodeCGW_T      *&endNodeA,
                                        NodeCGW_T      *&endNodeB,
                                        SequenceOrient  &orientEndNodeA,
                                        SequenceOrient  &orientEndNodeB,
                                        PairOrient      &edgeEndsOrient,
                                        double          &extremalGapSizeA,
                                        double          &extremalGapSizeB) {
  int AGapTowardAEnd = FALSE;
  int BGapTowardAEnd = FALSE;

#if 0
  fprintf(stderr,"* Translate (" F_CID "," F_CID ",%c)   " F_CID " has (" F_CID "," F_CID ")  " F_CID " has (" F_CID "," F_CID ")\n",
          scaffoldA->id, scaffoldB->id, scaffoldEdgeOrient,
          scaffoldA->id, scaffoldA->info.Scaffold.AEndCI, scaffoldA->info.Scaffold.BEndCI,
          scaffoldB->id, scaffoldB->info.Scaffold.AEndCI, scaffoldB->info.Scaffold.BEndCI);
#endif

  if (scaffoldEdgeOrient.isAB_AB()) {
    endNodeA       = GetTrueBEndCI(ScaffoldGraph, scaffoldA);
    endNodeB       = GetGraphNode(ScaffoldGraph->ContigGraph, scaffoldB->info.Scaffold.AEndCI);
    orientEndNodeA = GetNodeOrient(endNodeA);
    orientEndNodeB = GetNodeOrient(endNodeB);
    AGapTowardAEnd = TRUE;
    BGapTowardAEnd = FALSE;

  } else if (scaffoldEdgeOrient.isAB_BA()) {
    endNodeA = GetTrueBEndCI(ScaffoldGraph, scaffoldA);
    endNodeB = GetTrueBEndCI(ScaffoldGraph, scaffoldB);
    orientEndNodeA = GetNodeOrient(endNodeA);
    orientEndNodeB = GetNodeOrient(endNodeB);
    orientEndNodeB.flip();
    AGapTowardAEnd = TRUE;
    BGapTowardAEnd = TRUE;

  } else if (scaffoldEdgeOrient.isBA_AB()) {
    endNodeA = GetGraphNode(ScaffoldGraph->ContigGraph, scaffoldA->info.Scaffold.AEndCI);
    endNodeB = GetGraphNode(ScaffoldGraph->ContigGraph, scaffoldB->info.Scaffold.AEndCI);
    orientEndNodeA = GetNodeOrient(endNodeA);
    orientEndNodeA.flip();
    orientEndNodeB = GetNodeOrient(endNodeB);
    AGapTowardAEnd = FALSE;
    BGapTowardAEnd = FALSE;

  } else if (scaffoldEdgeOrient.isBA_BA()) {
    endNodeA = GetGraphNode(ScaffoldGraph->ContigGraph, scaffoldA->info.Scaffold.AEndCI);
    endNodeB = GetTrueBEndCI(ScaffoldGraph, scaffoldB);
    orientEndNodeA = GetNodeOrient(endNodeA);
    orientEndNodeB = GetNodeOrient(endNodeB);
    orientEndNodeA.flip();
    orientEndNodeB.flip();
    AGapTowardAEnd = FALSE;
    BGapTowardAEnd = TRUE;

  } else {
    assert(0);
  }

  if (orientEndNodeA.isForward()) {
    if (orientEndNodeB.isForward())
      edgeEndsOrient.setIsAB_AB();
    else
      edgeEndsOrient.setIsAB_BA();
  } else {
    if (orientEndNodeB.isForward())
      edgeEndsOrient.setIsBA_AB();
    else
      edgeEndsOrient.setIsBA_BA();
  }

  NodeCGW_T *nextNodeA = GetGraphNode(ScaffoldGraph->ContigGraph, AGapTowardAEnd ? endNodeA->AEndNext : endNodeA->BEndNext);
  NodeCGW_T *nextNodeB = GetGraphNode(ScaffoldGraph->ContigGraph, BGapTowardAEnd ? endNodeB->AEndNext : endNodeB->BEndNext);

  extremalGapSizeA = 1000000.0; // large number, anything will fit
  extremalGapSizeB = 1000000.0; // large number, anything will fit

  if (nextNodeA) {
    if (AGapTowardAEnd)
      extremalGapSizeA = MIN(endNodeA->offsetAEnd.mean, endNodeA->offsetBEnd.mean) - MAX(nextNodeA->offsetAEnd.mean, nextNodeA->offsetBEnd.mean);
    else
      extremalGapSizeA = MIN(nextNodeA->offsetAEnd.mean, nextNodeA->offsetBEnd.mean) - MAX(endNodeA->offsetAEnd.mean, endNodeA->offsetBEnd.mean);
  }

  if (nextNodeB) {
    if (BGapTowardAEnd)
      extremalGapSizeB = MIN(endNodeB->offsetAEnd.mean, endNodeB->offsetBEnd.mean) - MAX(nextNodeB->offsetAEnd.mean, nextNodeB->offsetBEnd.mean);
    else
      extremalGapSizeB = MIN(nextNodeB->offsetAEnd.mean, nextNodeB->offsetBEnd.mean) - MAX(endNodeB->offsetAEnd.mean, endNodeB->offsetBEnd.mean);
  }
}







static
int
LooseAbuttingCheck(SEdgeT * curEdge,
                   CIScaffoldT * scaffoldA,
                   CIScaffoldT * scaffoldB) {

  double stddevsPerWeight = -(curEdge->distance.mean + CGW_MISSED_OVERLAP) / (curEdge->edgesContributing * sqrt(MAX(0.01,curEdge->distance.variance)));
  double edgeMinScaffoldLengthRatio = -curEdge->distance.mean / MIN(scaffoldA->bpLength.mean, scaffoldB->bpLength.mean);

  //  The original abut rule.
  if ((stddevsPerWeight < STDDEVS_PER_WEIGHT_THRESHOLD) &&
      (edgeMinScaffoldLengthRatio < EDGE_PER_MIN_SCAFFOLD_LENGTH_THRESHOLD))
    return(TRUE);

  //  Be more aggressive about allowing abutments.
  if ((edgeMinScaffoldLengthRatio < MAX_PERC_SCAFFOLD_LEN) &&
      (-curEdge->distance.mean < MAX_OVERLAP_TO_ABUT))
    return(TRUE);

  return(FALSE);
}




/*
  Prior to marking scaffolds for merging, mean may be changed but not variance
  So, see if a change of mean to -CGW_MISSED_OVERLAP will allow abutting
  after which edge can be marked as trusted for recomputing offsets
*/
static
int
AbuttingWillWork(SEdgeT * curEdge,
                 CIScaffoldT * scaffoldA,
                 CIScaffoldT * scaffoldB,
                 InterleavingSpec * iSpec) {
  int isAbuttable = FALSE;
  double chiSquaredValue;

  if (curEdge->distance.mean >= -CGW_MISSED_OVERLAP)
    return TRUE;

  // the edge must be compatible with a -20 edge
  // or be strong & short relative to the scaffold lengths
  if (PairwiseChiSquare(-CGW_MISSED_OVERLAP,
                        curEdge->distance.variance,
                        curEdge->distance.mean,
                        curEdge->distance.variance,
                        (LengthT *) NULL,
                        &chiSquaredValue,
                        PAIRWISECHI2THRESHOLD_CGW) ||
      (LooseAbuttingCheck(curEdge, scaffoldA, scaffoldB))) {
    double originalMean = curEdge->distance.mean;

    curEdge->distance.mean = -CGW_MISSED_OVERLAP;

    // a -20 edge must have high satisfied mate %
    if (isQualityScaffoldMergingEdge(curEdge,
                                     scaffoldA,
                                     scaffoldB,
                                     iSpec->sai->scaffInst,
                                     iSpec->MIs,
                                     iSpec->minSatisfied,
                                     iSpec->maxDelta))
      isAbuttable = TRUE;

    curEdge->distance.mean = originalMean;
  }

  if (!isAbuttable)
    SaveBadScaffoldMergeEdge(curEdge, iSpec->badSEdges);

  return isAbuttable;
}




//static
void
ExamineSEdgeForUsability_Interleaved(SEdgeT * curEdge, InterleavingSpec * iSpec,
                                     CIScaffoldT *scaffoldA,
                                     CIScaffoldT *scaffoldB,
                                     double minMergeDistance,
                                     double maxMergeDistance,
                                     int32 mayOverlap,
                                     int32 mustOverlap) {

  EdgeCGW_T  *mergeEdge = NULL;

  ScaffoldAlignmentInterface * sai = iSpec->sai;

  if (minMergeDistance < -1000000. &&
      scaffoldA->bpLength.mean > 1000000. &&
      scaffoldB->bpLength.mean > 1000000.) {
    //fprintf(stderr, "Edge is too negative for scaffold lengths\n");
    return;
  }

  if (isBadScaffoldMergeEdge(curEdge, iSpec->badSEdges)) {
    //fprintf(stderr, "Edge previously marked as bad for merging.\n");
    return;
  }

  if (!isQualityScaffoldMergingEdge(curEdge,
                                    scaffoldA,
                                    scaffoldB,
                                    iSpec->sai->scaffInst,
                                    iSpec->MIs,
                                    iSpec->minSatisfied,
                                    iSpec->maxDelta)) {
    SaveBadScaffoldMergeEdge(curEdge, iSpec->badSEdges);
    return;
  }

  mergeEdge = curEdge;

  if ((mayOverlap  == false) &&
      (mustOverlap == false)) {
    // if here, edge was non-negative
    MarkScaffoldsForMerging(mergeEdge);
    return;
  }



  NodeCGW_T *endNodeA = NULL;
  NodeCGW_T *endNodeB = NULL;

  SequenceOrient orientEndNodeA;
  SequenceOrient orientEndNodeB;
  PairOrient     edgeEndsOrient;

  double aGapSize = 0.0;
  double bGapSize = 0.0;

  TranslateScaffoldOverlapToContigOverlap(scaffoldA,
                                          scaffoldB,
                                          mergeEdge->orient,
                                          endNodeA,
                                          endNodeB,
                                          orientEndNodeA,
                                          orientEndNodeB,
                                          edgeEndsOrient,
                                          aGapSize,
                                          bGapSize);

#if 0
  fprintf(stderr, "scaffoldA length: %d, elements: %d, end node + end gap: %d\n",
          (int) scaffoldA->bpLength.mean, scaffoldA->info.Scaffold.numElements, (int) (endNodeA->bpLength.mean + aGapSize));
  fprintf(stderr, "scaffoldB length: %d, elements: %d, end node + end gap: %d\n",
          (int) scaffoldB->bpLength.mean, scaffoldB->info.Scaffold.numElements, (int) (endNodeB->bpLength.mean + bGapSize));
  fprintf(stderr, "edge length: %.f, variance: %.f, weight: %d\n",
          mergeEdge->distance.mean, mergeEdge->distance.variance, mergeEdge->edgesContributing);
#endif


  // take the easy road if the edge is shorter than end nodes/gaps
  if (endNodeA->bpLength.mean + aGapSize > -minMergeDistance &&
      endNodeB->bpLength.mean + bGapSize > -minMergeDistance) {
    // single contig involved in each scaffold

    double chiSquaredValue = 0.0;
    int32  alternate = FALSE;

    EdgeCGW_T *overlapEdge = FindOverlapEdgeChiSquare(ScaffoldGraph,
                                                      endNodeA, endNodeB->id,
                                                      edgeEndsOrient,
                                                      mergeEdge->distance.mean,
                                                      mergeEdge->distance.variance,
                                                      &chiSquaredValue,
                                                      2 * PAIRWISECHI2THRESHOLD_CGW,
                                                      &alternate,
                                                      FALSE);

    if (overlapEdge != NULL) {
      // overlap edge found. use it
      SaveEdgeMeanForLater(mergeEdge);

      if (edgeEndsOrient!=GetEdgeOrientationWRT(overlapEdge,endNodeA->id))
        overlapEdge->distance.mean = -(endNodeA->bpLength.mean + endNodeB->bpLength.mean + overlapEdge->distance.mean);

      mergeEdge->distance.mean = overlapEdge->distance.mean;

      MarkScaffoldsForMerging(mergeEdge);
    }

    else if (!iSpec->checkAbutting || AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec)) {
      // no overlap edge, but abutting will work
      SaveEdgeMeanForLater(mergeEdge);
      mergeEdge->distance.mean = MAX(-CGW_MISSED_OVERLAP, mergeEdge->distance.mean);
      MarkScaffoldsForMerging(mergeEdge);
    }

    else {
      // else don't prevent scaffolds from merging via other edges
      SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
    }

    return;
  }



  // large negative edge - may involve interleaving

  //more contigs may be involved in one or both scaffolds
  if (PopulateScaffoldAlignmentInterface(scaffoldA, scaffoldB, mergeEdge, sai) != 0) {
    fprintf(stderr, "Failed to populate scaffold alignment interface!\n");
    SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
    return;
  }

#if 0
  // punt if the edge distance is too negative & there are no overlaps
  if ( (int) GetNumSegmentsInList(sai->segmentList) == 0 &&
       endNodeA->bpLength.mean + aGapSize < -minMergeDistance &&
       endNodeB->bpLength.mean + bGapSize < -minMergeDistance) {
    //fprintf(stderr, "Large scaffold overlap with no contig overlaps.  Edge not usable for merging.\n");
    SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
    return;
  }
#endif

  /* THERE MAY BE A MEMORY LEAK HERE:
     RETURN VALUE OF ALIGN_SCAFFOLD IS NULL IF NOTHING FOUND,
     BUT THIS DOESN'T MEAN THAT THE SEGMENT LIST GOT FREED ...
  */

  // this is still a problem but why?  (MP)
  assert(sai->numSegs==GetNumSegmentsInList(sai->segmentList));

  sai->segmentList = Align_Scaffold(sai->segmentList,
                                    sai->numSegs,
                                    sai->varWin,
                                    sai->scaffoldA->scaffold,
                                    sai->scaffoldB->scaffold,
                                    &(sai->best),
                                    sai->scaffoldA->bandBeg,
                                    sai->scaffoldA->bandEnd);

  if ((sai->segmentList == NULL) && (sai->best < 0)) {
    //  no overlap or interleaving possible
    //
    // if edge is not highly negative, abutt. Otherwise abort.
    // criterion is if -20 abutting still leaves edge as trusted

    if (!iSpec->checkAbutting || AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec)) {
      SaveEdgeMeanForLater(mergeEdge);
      mergeEdge->distance.mean = MAX(-CGW_MISSED_OVERLAP, mergeEdge->distance.mean);
      MarkScaffoldsForMerging(mergeEdge);

    } else {
      // else don't prevent scaffolds from merging via other edges
      // record that this scaffold overlap is bad
      SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
    }

    return;
  }


  if (sai->segmentList != NULL || (iSpec->checkAbutting && !AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec))) {

    // if there are overlaps or abutting isn't an option
    EdgeCGW_T *overlapEdge = MakeScaffoldAlignmentAdjustments(scaffoldA, scaffoldB, mergeEdge, sai);

    if (overlapEdge != NULL) {
      SaveEdgeMeanForLater(mergeEdge);
      mergeEdge->distance.mean = overlapEdge->distance.mean;
      MarkScaffoldsForMerging(mergeEdge);

    } else {
      // else don't prevent scaffolds from merging via other edges
      SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
    }

    return;
  }

  // no overlaps && abutting will work
  SaveEdgeMeanForLater(mergeEdge);
  mergeEdge->distance.mean = MAX(-CGW_MISSED_OVERLAP, mergeEdge->distance.mean);
  MarkScaffoldsForMerging(mergeEdge);
}


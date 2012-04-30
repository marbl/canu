
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
static char *rcsid = "$Id: CIScaffoldT_Merge_CGW.c,v 1.64 2012-04-30 22:12:11 brianwalenz Exp $";

//
//  The ONLY exportable function here is MergeScaffoldsAggressive.
//

#undef DEBUG_MERGE_EDGE_INVERT

//  Define this to check (and assert) if the graph is not internally
//  connected before recomputing offsets.  It's expensive, and if you
//  already know it's OK (debugging, maybe??) you can skip it.
#define CHECKCONNECTED

#define OTHER_END_CHECK
#undef  GENERAL_STRONGER_CHECK

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Instrument_CGW.h"
#include "ChiSquareTest_CGW.h"
#include "InterleavedMerging.h"


#define SCAFFOLD_MERGE_CHI2_THRESHHOLD (2.f*(double)PAIRWISECHI2THRESHOLD_CGW)

#define PREFERRED_GAP_SIZE  (-500)

#define MAX_SCAFFOLD_GAP_OVERLAP 1000

#define CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD MIN_EDGES

#define EDGE_QUANTA 5.0
#define OVERLAP_QUANTA -10000.

#define EDGE_STRENGTH_FACTOR  MIN_EDGES

#define MAX_SLOP_IN_STD 3.5

#define EDGE_WEIGHT_FACTOR  MIN_EDGES



static
void
SaveEdgeMeanForLater(SEdgeT * edge) {
  if (edge->flags.bits.MeanChangedByWalking == FALSE) {
    edge->flags.bits.MeanChangedByWalking = TRUE;
    edge->minDistance = edge->distance.mean;
  }
}





static
double
GetVarianceOffset(CIScaffoldT * scaffold, double meanOffset, int isAB) {
  double lengthDelta = (isAB) ? 0.0 : scaffold->bpLength.mean;
  double varDelta = (isAB) ? 0.0 : scaffold->bpLength.variance;
  CIScaffoldTIterator CIs;
  ChunkInstanceT * CI;
  LengthT lowerOffset = {-1, 0};
  LengthT higherOffset = {0, 0};

  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, isAB, FALSE, &CIs);
  while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
    double aEnd = fabs(lengthDelta - CI->offsetAEnd.mean);
    double bEnd = fabs(lengthDelta - CI->offsetBEnd.mean);
    if (aEnd > meanOffset || bEnd > meanOffset) {
      if (aEnd < meanOffset) {
        lowerOffset = CI->offsetAEnd;
        higherOffset = CI->offsetBEnd;
        break;
      } else if (bEnd < meanOffset) {
        lowerOffset = CI->offsetBEnd;
        higherOffset = CI->offsetAEnd;
        break;
      } else {
        higherOffset = (aEnd < bEnd) ? CI->offsetAEnd : CI->offsetBEnd;
        break;
      }
    }
    lowerOffset = (aEnd < bEnd) ? CI->offsetBEnd : CI->offsetAEnd;
  }
  if (higherOffset.mean != 0)
    return fabs(varDelta - lowerOffset.variance) +
      (fabs(varDelta - higherOffset.variance) - fabs(varDelta - lowerOffset.variance)) *
      (meanOffset - fabs(lengthDelta - lowerOffset.mean)) /
      (fabs(lengthDelta - higherOffset.mean) - fabs(lengthDelta - lowerOffset.mean));
  else
    return scaffold->bpLength.variance * meanOffset / scaffold->bpLength.mean;
}


static
void
InsertScaffoldContentsIntoScaffold(ScaffoldGraphT *sgraph,
                                   CDS_CID_t newScaffoldID,
                                   CDS_CID_t oldScaffoldID,
                                   SequenceOrient orient,
                                   LengthT * offset,
                                   int contigNow) {
  CIScaffoldT *oldScaffold = GetCIScaffoldT(sgraph->CIScaffolds, oldScaffoldID);
  CIScaffoldT *newScaffold = GetCIScaffoldT(sgraph->CIScaffolds, newScaffoldID);

  // CheckCIScaffoldTLength(sgraph, oldScaffold);
  // CheckCIScaffoldTLength(sgraph, newScaffold);

  fprintf(stderr,"InsertScaffoldContentsIntoScaffold()-- Insert scaffold "F_CID" (%.0fbp) into scaffold "F_CID" (%.0fbp) at offset %.3f +/- %.3f orient %c\n",
          oldScaffoldID, oldScaffold->bpLength.mean,
          newScaffoldID, newScaffold->bpLength.mean,
          offset->mean, sqrt(offset->variance),
          orient.toLetter());

  assert(offset->mean     >= 0.0);
  assert(offset->variance >= 0.0);

  //  Shouldn't be necessary, occasionally (in GOSIII) the length was wrong.
  SetCIScaffoldTLength(sgraph, newScaffold, FALSE);
  SetCIScaffoldTLength(sgraph, oldScaffold, FALSE);

  CIScaffoldTIterator CIs;
  ChunkInstanceT     *CI;

  InitCIScaffoldTIterator(sgraph, oldScaffold, orient.isForward(), FALSE, &CIs);
  while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
    LengthT offsetAEnd;
    LengthT offsetBEnd;

    //  These should be guaranteed by the SetCIScaffoldTLength() above.
    assert(CI->offsetAEnd.mean <= oldScaffold->bpLength.mean);
    assert(CI->offsetBEnd.mean <= oldScaffold->bpLength.mean);

    //  These guaranteed too.  Previous versions would detect the case and simply reset the variance
    //  of the placed contig (that's offsetBEnd) to the variance of the 'offset'.  That is:
    //    offsetBEnd.variance = offset->variance
    //    offsetAEnd.variance = offset->variance
    assert(CI->offsetAEnd.variance <= oldScaffold->bpLength.variance);
    assert(CI->offsetBEnd.variance <= oldScaffold->bpLength.variance);

    if (orient.isForward()) {
      offsetAEnd.mean     = offset->mean     + CI->offsetAEnd.mean;
      offsetAEnd.variance = offset->variance + CI->offsetAEnd.variance;
      offsetBEnd.mean     = offset->mean     + CI->offsetBEnd.mean;
      offsetBEnd.variance = offset->variance + CI->offsetBEnd.variance;
    } else {
      offsetAEnd.mean     = offset->mean     + (oldScaffold->bpLength.mean     - CI->offsetAEnd.mean);
      offsetAEnd.variance = offset->variance + (oldScaffold->bpLength.variance - CI->offsetAEnd.variance);
      offsetBEnd.mean     = offset->mean     + (oldScaffold->bpLength.mean     - CI->offsetBEnd.mean);
      offsetBEnd.variance = offset->variance + (oldScaffold->bpLength.variance - CI->offsetBEnd.variance);
    }

    //  Unfortunately, if oldScaffold has screwed up variances, 
    //assert(offsetAEnd.variance >= 0.0);
    //assert(offsetBEnd.variance >= 0.0);

    fprintf(stderr,"InsertScaffoldContentsIntoScaffold()-- Insert CI "F_CID" (%.0fbp) at offset (%.0f,%.0f); was at (%.0f,%.0f)\n",
            CI->id,
            CI->bpLength.mean,
            offsetAEnd.mean,     offsetBEnd.mean,
            CI->offsetAEnd.mean, CI->offsetBEnd.mean);

    InsertCIInScaffold(sgraph, CI->id, newScaffoldID, offsetAEnd, offsetBEnd, TRUE, contigNow);
  }

  //  Shouldn't be necessary, occasionally (in GOSIII) the length was wrong.
  SetCIScaffoldTLength(sgraph, newScaffold, FALSE);

  //  Check all the contig coordinatess in the new scaffold.
#if 1
  LengthT lastMin   = {0, 0};
  LengthT lastMax   = {0, 0};
  int     thisIdx   = 0;
  int     thisFails = 0;

  InitCIScaffoldTIterator(ScaffoldGraph, newScaffold, TRUE, FALSE, &CIs);
  while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {

    if ((CI->offsetAEnd.mean < 0.0) || (CI->offsetAEnd.variance < 0.0) ||
        (CI->offsetBEnd.mean < 0.0) || (CI->offsetBEnd.variance < 0.0)) {
      fprintf(stderr, "InsertScaffoldContentsIntoScaffold()-- Negative contig mean or variance\n");
      thisFails++;
    }

    if ((thisIdx++ != 0) && ((CI->offsetAEnd.mean == 0.0) ||
                             (CI->offsetBEnd.mean == 0.0))) {
      fprintf(stderr, "InsertScaffoldContentsIntoScaffold()-- Zero offset of internal contig\n");
      thisFails++;
    }

    LengthT thisMin = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;
    LengthT thisMax = (CI->offsetAEnd.mean > CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;

    if (thisMax.mean < lastMin.mean) {
      fprintf(stderr, "InsertScaffoldContentsIntoScaffold()-- Seriously out of order contigs\n");
      thisFails++;
    }

    if (thisMin.mean < lastMin.mean) {
      fprintf(stderr, "InsertScaffoldContentsIntoScaffold()-- Out of order contigs\n");
      thisFails++;
    }

    lastMin = thisMin;
    lastMax = thisMax;
  }

  if (thisFails)
    DumpCIScaffold(stderr, ScaffoldGraph, newScaffold, FALSE);
  //assert(thisFails == 0);
#endif
}



static
void
MarkUnderlyingRawCIEdgeTrusted(ScaffoldGraphT * sgraph, EdgeCGW_T * raw) {
  CIEdgeT *ciedge, *topCIEdge;
  ContigT *contigA, *contigB;

  ciedge = GetGraphEdge(sgraph->ContigGraph,raw->referenceEdge);
  topCIEdge = GetGraphEdge(sgraph->ContigGraph, ciedge->topLevelEdge);
  assert(topCIEdge->idA != NULLINDEX && topCIEdge->idB != NULLINDEX);
  contigA = GetGraphNode(sgraph->ContigGraph, topCIEdge->idA);
  contigB = GetGraphNode(sgraph->ContigGraph, topCIEdge->idB);
  assert(contigA->scaffoldID == contigB->scaffoldID);
  assert(contigA->flags.bits.isUnique);
  assert(contigB->flags.bits.isUnique);
  AssertPtr(topCIEdge);
  SetGraphEdgeStatus(sgraph->ContigGraph, topCIEdge, TRUSTED_EDGE_STATUS);
  if (GlobalData->debugLevel > 0) {
    fprintf(stderr,"* Marked contig edge " F_CID " (" F_CID "," F_CID ")%c as trusted(inside scaf " F_CID ")\n",
            topCIEdge->topLevelEdge, topCIEdge->idA, topCIEdge->idB, topCIEdge->orient.toLetter(),
            contigA->scaffoldID);
  }
}


static
void 
MarkUnderlyingCIEdgesTrusted(ScaffoldGraphT *sgraph, SEdgeT *edge) {
  if (GlobalData->debugLevel > 0)
    fprintf(stderr,"* MarkUnderlyingCIEdgesTrusted on SEdge (" F_CID "," F_CID ")%c nextRaw = " F_CID "\n",
            edge->idA, edge->idB, edge->orient.toLetter(), edge->nextRawEdge);

  if (edge->flags.bits.isRaw) {
    MarkUnderlyingRawCIEdgeTrusted(sgraph, edge);
  } else {
    SEdgeT *raw = edge;
    while ((raw = GetGraphEdge(sgraph->ScaffoldGraph, raw->nextRawEdge)) != NULL)
      MarkUnderlyingRawCIEdgeTrusted(sgraph, raw);
  }
}



static
double
FindMaxGapInScaffold(ScaffoldGraphT *graph, CIScaffoldT *scaffold) {
  CIScaffoldTIterator Nodes;
  NodeCGW_T *thisNode;
  LengthT *prevEnd, *thisBegin, *thisEnd;
  double maxGap;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &Nodes);
  thisNode = NextCIScaffoldTIterator(&Nodes);
  if (GetNodeOrient(thisNode).isForward()) {
    prevEnd = &(thisNode->offsetBEnd);
  } else {
    prevEnd = &(thisNode->offsetAEnd);
  }
  for(maxGap = 0.0; (thisNode = NextCIScaffoldTIterator(&Nodes)) != NULL; prevEnd = thisEnd) {
    double thisGap;
    if (GetNodeOrient(thisNode).isForward()) {
      thisBegin = &(thisNode->offsetAEnd);
      thisEnd = &(thisNode->offsetBEnd);
    } else {
      thisEnd = &(thisNode->offsetAEnd);
      thisBegin = &(thisNode->offsetBEnd);
    }
    thisGap = (thisBegin->mean - prevEnd->mean) + (3.0 * sqrt(thisBegin->variance - prevEnd->variance));
    if (thisGap > maxGap) {
      maxGap = thisGap;
    }
  }
  return(maxGap);
}


static
int
CompareSEdgesContributing(const void *c1, const void *c2) {
  SEdgeT *s1 = *(SEdgeT **)c1;
  SEdgeT *s2 = *(SEdgeT **)c2;

  int32  n2 = (s2->edgesContributing - (isOverlapEdge(s2) ? 1 : 0));
  int32  n1 = (s1->edgesContributing - (isOverlapEdge(s1) ? 1 : 0));

  if (n2 < n1)
    //  Edge1 higher weight
    return(-1);

  if (n1 < n2)
    //  Edge1 not higher weight
    return(1);

  //  Otherwise, they're the same weight.  Break ties using the insert size.

  double  d1 = fabs(PREFERRED_GAP_SIZE - s1->distance.mean) + sqrt(MAX(1., s1->distance.variance));
  double  d2 = fabs(PREFERRED_GAP_SIZE - s2->distance.mean) + sqrt(MAX(1., s2->distance.variance));

  if (d1 < d2)
    //  Edge1 closer to truth
    return(-1);

  if (d2 < d1)
    //  Edge1 not closer to truth
    return(1);

  //  Tied.
  return(0);
}




typedef struct {
  NodeCGW_T *nodeA;
  NodeCGW_T *nodeB;
  int numToMerge;
  double rating;
  double maxSpace;
  double maxGapSize;
  LengthT gapSize;
}SortGapMergeT;

typedef struct {
  CDS_CID_t nodeIndex;
  LengthT offsetAEnd;
  LengthT offsetBEnd;
  LengthT gapBeginOffset;
  LengthT gapDelta;
  int numInGap;
  int order;
}OrderCIsTmpT;


static
int
CompareGapRating(const void *c1, const void *c2) {
  SortGapMergeT *s1 = (SortGapMergeT *)c1;
  SortGapMergeT *s2 = (SortGapMergeT *)c2;

  return(((s1->rating - s2->rating) < 0.0) ? -1 : 1);
}


static
void
SetOrderInfoForGap(OrderCIsTmpT *orderAnchor, int32 index,
                   int firstTime, double oldMergeGapVariance,
                   LengthT anchorGap,
                   double newMergeGapVariance,
                   double gapBeginMean, int pathFromA,
                   LengthT insertSize,
                   double gapEndMean, int pathFromB,
                   int numNodesInGap, int endOfMerge) {
  double deltaMean;
  double gapEndVariance;
  double gapEndStdDev;
  double gapBeginStdDev;
  double ratioStdDev;

  orderAnchor[index].numInGap = numNodesInGap;
  if (firstTime)
    orderAnchor[index].gapBeginOffset.variance = anchorGap.variance;
  else
    orderAnchor[index].gapBeginOffset.variance = MAX(anchorGap.variance, oldMergeGapVariance);

  gapBeginStdDev = sqrt(orderAnchor[index].gapBeginOffset.variance);
  orderAnchor[index].gapDelta.variance = orderAnchor[index].gapBeginOffset.variance + insertSize.variance;

  if (endOfMerge)
    gapEndVariance = anchorGap.variance;
  else
    gapEndVariance = MAX(anchorGap.variance, newMergeGapVariance);

  gapEndStdDev = sqrt(gapEndVariance);
  ratioStdDev = gapBeginStdDev / (gapBeginStdDev + gapEndStdDev);
  orderAnchor[index].gapDelta.variance += gapEndVariance;
  orderAnchor[index].gapDelta.variance -= anchorGap.variance;
  if (pathFromA && pathFromB) {
    if ((gapBeginMean < - (double)CGW_MISSED_OVERLAP) &&
        (gapEndMean < - (double)CGW_MISSED_OVERLAP)) {
      orderAnchor[index].gapBeginOffset.mean = gapBeginMean;
      orderAnchor[index].gapDelta.mean = (gapBeginMean + insertSize.mean +
                                          gapEndMean) - anchorGap.mean;
    }else if (gapBeginMean < - (double)CGW_MISSED_OVERLAP) {
      orderAnchor[index].gapBeginOffset.mean = gapBeginMean;
      deltaMean = (gapBeginMean + insertSize.mean) - anchorGap.mean;
      if (deltaMean <= (double)CGW_MISSED_OVERLAP) {
        orderAnchor[index].gapDelta.mean = 0.0;
      } else {
        orderAnchor[index].gapDelta.mean = deltaMean -
          (double)CGW_MISSED_OVERLAP;
      }
    }else if (gapEndMean < - (double)CGW_MISSED_OVERLAP) {
      deltaMean = (gapEndMean + insertSize.mean) - anchorGap.mean;
      if (deltaMean <= (double)CGW_MISSED_OVERLAP) {
        orderAnchor[index].gapDelta.mean = 0.0;
        orderAnchor[index].gapBeginOffset.mean = - deltaMean;
      } else {
        orderAnchor[index].gapDelta.mean = deltaMean -
          (double)CGW_MISSED_OVERLAP;
        orderAnchor[index].gapBeginOffset.mean = - (double)CGW_MISSED_OVERLAP;
      }
    } else {
      deltaMean = (gapBeginMean + insertSize.mean + gapEndMean) -
        anchorGap.mean;
      if (deltaMean <= 0.0) {
        orderAnchor[index].gapDelta.mean = 0.0;
        orderAnchor[index].gapBeginOffset.mean = gapBeginMean -
          (deltaMean * ratioStdDev);
      } else {
        orderAnchor[index].gapBeginOffset.mean = gapBeginMean -
          (deltaMean * ratioStdDev);
        if (orderAnchor[index].gapBeginOffset.mean <
            - (double)CGW_MISSED_OVERLAP) {
          orderAnchor[index].gapBeginOffset.mean = - (double)CGW_MISSED_OVERLAP;
        }
        deltaMean = (orderAnchor[index].gapBeginOffset.mean + insertSize.mean)
          - anchorGap.mean;
        if (deltaMean <= (double)CGW_MISSED_OVERLAP) {
          orderAnchor[index].gapDelta.mean = 0.0;
        } else {
          orderAnchor[index].gapBeginOffset.mean -= deltaMean -
            (double)CGW_MISSED_OVERLAP;
          if (orderAnchor[index].gapBeginOffset.mean <
              - (double)CGW_MISSED_OVERLAP) {
            orderAnchor[index].gapBeginOffset.mean =
              - (double)CGW_MISSED_OVERLAP;
          }
          deltaMean = (orderAnchor[index].gapBeginOffset.mean +
                       insertSize.mean) - anchorGap.mean;
          if (deltaMean <= (double)CGW_MISSED_OVERLAP) {
            orderAnchor[index].gapDelta.mean = 0.0;
          } else {
            orderAnchor[index].gapDelta.mean = deltaMean -
              (double)CGW_MISSED_OVERLAP;
          }
        }
      }
    }
  }else if (pathFromA) {
    if (gapBeginMean < - (double)CGW_MISSED_OVERLAP) {
      orderAnchor[index].gapBeginOffset.mean = gapBeginMean;
      deltaMean = (gapBeginMean + insertSize.mean) - anchorGap.mean;
      if (deltaMean <= (double)CGW_MISSED_OVERLAP) {
        orderAnchor[index].gapDelta.mean = 0.0;
      } else {
        orderAnchor[index].gapDelta.mean = deltaMean -
          (double)CGW_MISSED_OVERLAP;
      }
    } else {
      deltaMean = (gapBeginMean + insertSize.mean) - anchorGap.mean;
      orderAnchor[index].gapBeginOffset.mean = gapBeginMean;
      if (deltaMean <= (double)CGW_MISSED_OVERLAP) {
        orderAnchor[index].gapDelta.mean = 0.0;
      } else {
        orderAnchor[index].gapBeginOffset.mean -= deltaMean -
          (double)CGW_MISSED_OVERLAP;
        if (orderAnchor[index].gapBeginOffset.mean <
            - (double)CGW_MISSED_OVERLAP) {
          orderAnchor[index].gapBeginOffset.mean =
            - (double)CGW_MISSED_OVERLAP;
        }
        deltaMean = (orderAnchor[index].gapBeginOffset.mean +
                     insertSize.mean) - anchorGap.mean;
        if (deltaMean <= (double)CGW_MISSED_OVERLAP) {
          orderAnchor[index].gapDelta.mean = 0.0;
        } else {
          orderAnchor[index].gapDelta.mean = deltaMean -
            (double)CGW_MISSED_OVERLAP;
        }
      }
    }
  }else if (pathFromB) {
    if (gapEndMean < - (double)CGW_MISSED_OVERLAP) {
      deltaMean = (gapEndMean + insertSize.mean) - anchorGap.mean;
      if (deltaMean <= (double)CGW_MISSED_OVERLAP) {
        orderAnchor[index].gapDelta.mean = 0.0;
        orderAnchor[index].gapBeginOffset.mean = - deltaMean;
      } else {
        orderAnchor[index].gapDelta.mean = deltaMean -
          (double)CGW_MISSED_OVERLAP;
        orderAnchor[index].gapBeginOffset.mean = - (double)CGW_MISSED_OVERLAP;
      }
    } else {
      deltaMean = (gapEndMean + insertSize.mean) - anchorGap.mean;
      if (deltaMean <= (double)CGW_MISSED_OVERLAP) {
        orderAnchor[index].gapDelta.mean = 0.0;
        orderAnchor[index].gapBeginOffset.mean = - deltaMean;
      } else {
        gapEndMean -= deltaMean - (double)CGW_MISSED_OVERLAP;
        if (gapEndMean < - (double)CGW_MISSED_OVERLAP) {
          gapEndMean = - (double)CGW_MISSED_OVERLAP;
        }
        deltaMean = (gapEndMean + insertSize.mean) - anchorGap.mean;
        if (deltaMean <= (double)CGW_MISSED_OVERLAP) {
          orderAnchor[index].gapDelta.mean = 0.0;
          orderAnchor[index].gapBeginOffset.mean = - deltaMean;
        } else {
          orderAnchor[index].gapDelta.mean = deltaMean -
            (double)CGW_MISSED_OVERLAP;
          orderAnchor[index].gapBeginOffset.mean = - (double)CGW_MISSED_OVERLAP;
        }
      }
    }
  } else {
    assert(FALSE);//!pathFromA && !pathFromB
  }

  return;
}




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
isLargeOverlapSize(LengthT *overlap, int32 numEdges) {
  double maxGap = overlap->mean + 3.0 * sqrt(overlap->variance);
  double minGap = overlap->mean - 3.0 * sqrt(overlap->variance);

  double numOverlapQuanta = minGap/OVERLAP_QUANTA;
  double numEdgeQuanta    = numEdges/EDGE_QUANTA;

  fprintf(stderr, "isLargeOverlapSize()--  edges %d gap %f %f ovlQuanta %f (%f) edgeQuanta %f (%f)\n",
          numEdges,
          minGap, maxGap, 
          numOverlapQuanta, OVERLAP_QUANTA,
          numEdgeQuanta,    EDGE_QUANTA);

  if (maxGap > OVERLAP_QUANTA)
    return FALSE;

  return(numEdgeQuanta < numOverlapQuanta);
}

static
int
isLargeOverlap(CIEdgeT *curEdge) {
  return(isLargeOverlapSize(&(curEdge->distance), curEdge->edgesContributing));
}

#ifdef GENERAL_STRONGER_CHECK
/*
  it would be better to do this during scaffold edge creation
*/
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
          if (retVal == 0) {
            fprintf(stderr,
                    "SCF MERGE CONFLICT: " F_CID "," F_CID "  %s  %dbp  %dvar  %dec\n",
                    curSEdge->idA, curSEdge->idB,
                    ((curSEdge->orient.isAB_AB()) ? "AB_AB" :
                     ((curSEdge->orient.isAB_BA()) ? "AB_BA" :
                      ((curSEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
                    (int) curSEdge->distance.mean,
                    (int) curSEdge->distance.variance,
                    curSEdge->edgesContributing);
          }
          fprintf(stderr, "\t" F_CID "," F_CID ", %s, %dbp  %dvar  %dec\n",
                  sEdge->idA, sEdge->idB,
                  ((sEdge->orient.isAB_AB()) ? "AB_AB" :
                   ((sEdge->orient.isAB_BA()) ? "AB_BA" :
                    ((sEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
                  (int) sEdge->distance.mean,
                  (int) curSEdge->distance.variance,
                  sEdge->edgesContributing);
          retVal++;
        }
      } else {
        // shift
        if (curSEdge->distance.mean > 0) {
          if (sEdge->distance.mean > 0) {
            // both are positive, prefer stronger
            if (sEdge->edgesContributing > curSEdge->edgesContributing) {
              if (retVal == 0) {
                fprintf(stderr,
                        "SCF MERGE CONFLICT: " F_CID "," F_CID "  %s  %dbp  %dvar  %dec\n",
                        curSEdge->idA, curSEdge->idB,
                        ((curSEdge->orient.isAB_AB()) ? "AB_AB" :
                         ((curSEdge->orient.isAB_BA()) ? "AB_BA" :
                          ((curSEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
                        (int) curSEdge->distance.mean,
                        (int) curSEdge->distance.variance,
                        curSEdge->edgesContributing);
              }
              fprintf(stderr, "\t" F_CID "," F_CID ", %s, %dbp  %dvar  %dec\n",
                      sEdge->idA, sEdge->idB,
                      ((sEdge->orient.isAB_AB()) ? "AB_AB" :
                       ((sEdge->orient.isAB_BA()) ? "AB_BA" :
                        ((sEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
                      (int) sEdge->distance.mean,
                      (int) curSEdge->distance.variance,
                      sEdge->edgesContributing);
              retVal++;
            }
          } else {
            /*
              curSEdge is positive, sEdge is negative
              prefer curSEdge unless sEdge is much stronger
            */
            if (sEdge->edgesContributing >
                EDGE_STRENGTH_FACTOR * curSEdge->edgesContributing) {
              if (retVal == 0) {
                fprintf(stderr,
                        "SCF MERGE CONFLICT: " F_CID "," F_CID "  %s  %dbp  %dvar  %dec\n",
                        curSEdge->idA, curSEdge->idB,
                        ((curSEdge->orient.isAB_AB()) ? "AB_AB" :
                         ((curSEdge->orient.isAB_BA()) ? "AB_BA" :
                          ((curSEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
                        (int) curSEdge->distance.mean,
                        (int) curSEdge->distance.variance,
                        curSEdge->edgesContributing);
              }
              fprintf(stderr, "\t" F_CID "," F_CID ", %s, %dbp  %dvar  %dec\n",
                      sEdge->idA, sEdge->idB,
                      ((sEdge->orient.isAB_AB()) ? "AB_AB" :
                       ((sEdge->orient.isAB_BA()) ? "AB_BA" :
                        ((sEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
                      (int) sEdge->distance.mean,
                      (int) curSEdge->distance.variance,
                      sEdge->edgesContributing);
              retVal++;
            }
          }
        } else if (sEdge->distance.mean > 0) {
          /*
            sEdge is positive, curSEdge is negative
            prefer sEdge unless curSEdge is much stronger
          */
          if (curSEdge->edgesContributing <
              EDGE_STRENGTH_FACTOR * sEdge->edgesContributing) {
            if (retVal == 0) {
              fprintf(stderr,
                      "SCF MERGE CONFLICT: " F_CID "," F_CID "  %s  %dbp  %dvar  %dec\n",
                      curSEdge->idA, curSEdge->idB,
                      ((curSEdge->orient.isAB_AB()) ? "AB_AB" :
                       ((curSEdge->orient.isAB_BA()) ? "AB_BA" :
                        ((curSEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
                      (int) curSEdge->distance.mean,
                      (int) curSEdge->distance.variance,
                      curSEdge->edgesContributing);
            }
            fprintf(stderr, "\t" F_CID "," F_CID ", %s, %dbp  %dvar  %dec\n",
                    sEdge->idA, sEdge->idB,
                    ((sEdge->orient.isAB_AB()) ? "AB_AB" :
                     ((sEdge->orient.isAB_BA()) ? "AB_BA" :
                      ((sEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
                    (int) sEdge->distance.mean,
                    (int) curSEdge->distance.variance,
                    sEdge->edgesContributing);
            retVal++;
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
            if (retVal == 0) {
              fprintf(stderr,
                      "SCF MERGE CONFLICT: " F_CID "," F_CID "  %s  %dbp  %dvar  %dec\n",
                      curSEdge->idA, curSEdge->idB,
                      ((curSEdge->orient.isAB_AB()) ? "AB_AB" :
                       ((curSEdge->orient.isAB_BA()) ? "AB_BA" :
                        ((curSEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
                      (int) curSEdge->distance.mean,
                      (int) curSEdge->distance.variance,
                      curSEdge->edgesContributing);
            }
            fprintf(stderr, "\t" F_CID "," F_CID ", %s, %dbp  %dvar  %dec\n",
                    sEdge->idA, sEdge->idB,
                    ((sEdge->orient.isAB_AB()) ? "AB_AB" :
                     ((sEdge->orient.isAB_BA()) ? "AB_BA" :
                      ((sEdge->orient.isBA_AB()) ? "BA_AB" : "BA_BA"))),
                    (int) sEdge->distance.mean,
                    (int) curSEdge->distance.variance,
                    sEdge->edgesContributing);
            retVal++;
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
FindAllMergeCandidates(VA_TYPE(PtrT) *sEdges,
                       VA_TYPE(PtrT) *overlapSEdges,
                       CIScaffoldT *fromScaffold,
                       int fromEnd,
                       CIScaffoldT *ignoreToScaffold,
                       int canonicalOnly,
                       int minWeight, int doInterleaving,
                       int verbose) {
  SEdgeTIterator SEdges;
  SEdgeT *curSEdge;
  CIScaffoldT *otherScaffold;
  CDS_CID_t otherScaffoldID;

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

    if (!doInterleaving && isLargeOverlap(curSEdge)) {
      if (overlapSEdges && curSEdge->edgesContributing >= CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD)
        AppendPtrT(overlapSEdges, (void **)&curSEdge);
      continue;
    }

    AppendPtrT(sEdges, (void **)&curSEdge);
  }

  return FALSE;
}


static
void
BuildUsableSEdges(VA_TYPE(PtrT) *sEdges,
                  VA_TYPE(PtrT) *overlapSEdges,
                  int32 verbose) {

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

    FindAllMergeCandidates(sEdges, overlapSEdges,
                           scaffold,
                           ALL_END, NULL,
                           TRUE,
                           CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD,
                           GlobalData->doInterleavedScaffoldMerging,
                           verbose);
  }

  if (GetNumPtrTs(sEdges) > 0) {
    SEdgeT **sEdge = (SEdgeT **)GetPtrT(sEdges,0);

    qsort(sEdge, GetNumPtrTs(sEdges), sizeof(SEdgeT *), CompareSEdgesContributing);
  }
}



static
ChunkInstanceT *
GetTrueBEndCI(ScaffoldGraphT * graph, CIScaffoldT * scaffold) {
  ChunkInstanceT * returnContig =
    GetGraphNode(graph->ContigGraph, scaffold->info.Scaffold.BEndCI);
  CIScaffoldTIterator contigIterator;
  ChunkInstanceT * myContig;

  InitCIScaffoldTIterator(graph, scaffold, FALSE, FALSE, &contigIterator);
  while ((myContig = NextCIScaffoldTIterator(&contigIterator)) != NULL) {
    if (MAX(returnContig->offsetAEnd.mean, returnContig->offsetBEnd.mean) <
        MAX(myContig->offsetAEnd.mean, myContig->offsetBEnd.mean))
      returnContig = myContig;
    if (MIN(returnContig->offsetAEnd.mean, returnContig->offsetBEnd.mean) >
        MAX(myContig->offsetAEnd.mean, myContig->offsetBEnd.mean))
      break;
  }
  return returnContig;
}


// TranslateScaffoldOverlapToContigOverlap
//       Prepare the groundwork for computing an overlap between the appropriate extremal
//       contigs in a pair of scaffolds.
//
static
void
TranslateScaffoldOverlapToContigOverlap(CIScaffoldT *scaffoldA, CIScaffoldT *scaffoldB,
                                        PairOrient scaffoldEdgeOrient,
                                        NodeCGW_T **endNodeA, NodeCGW_T **endNodeB,
                                        SequenceOrient *orientEndNodeA,  SequenceOrient *orientEndNodeB,
                                        PairOrient *edgeEndsOrient,
                                        double *extremalGapSizeA, double *extremalGapSizeB) {
  NodeCGW_T *nextNodeA, *nextNodeB;
  int AGapTowardAEnd, BGapTowardAEnd;
#if 0
  fprintf(stderr,"* Translate (" F_CID "," F_CID ",%c)   " F_CID " has (" F_CID "," F_CID ")  " F_CID " has (" F_CID "," F_CID ")\n",
          scaffoldA->id, scaffoldB->id, scaffoldEdgeOrient,
          scaffoldA->id, scaffoldA->info.Scaffold.AEndCI, scaffoldA->info.Scaffold.BEndCI,
          scaffoldB->id, scaffoldB->info.Scaffold.AEndCI, scaffoldB->info.Scaffold.BEndCI);
#endif
  if (scaffoldEdgeOrient.isAB_AB()) {
    // BEndCI may be contained
    *endNodeA = GetTrueBEndCI(ScaffoldGraph, scaffoldA);
    *endNodeB = GetGraphNode(ScaffoldGraph->ContigGraph,
                             scaffoldB->info.Scaffold.AEndCI);
    *orientEndNodeA = GetNodeOrient(*endNodeA);
    *orientEndNodeB = GetNodeOrient(*endNodeB);

    AGapTowardAEnd = TRUE;
    BGapTowardAEnd = FALSE;
  }else if (scaffoldEdgeOrient.isAB_BA()) {
    // BendCI may be contained
    *endNodeA = GetTrueBEndCI(ScaffoldGraph, scaffoldA);
    // BendCI may be contained
    *endNodeB = GetTrueBEndCI(ScaffoldGraph, scaffoldB);
    *orientEndNodeA = GetNodeOrient(*endNodeA);
    *orientEndNodeB = GetNodeOrient(*endNodeB);
    orientEndNodeB->flip();
    AGapTowardAEnd = TRUE;
    BGapTowardAEnd = TRUE;
  }else if (scaffoldEdgeOrient.isBA_AB()) {
    *endNodeA = GetGraphNode(ScaffoldGraph->ContigGraph, scaffoldA->info.Scaffold.AEndCI);
    *endNodeB = GetGraphNode(ScaffoldGraph->ContigGraph, scaffoldB->info.Scaffold.AEndCI);
    *orientEndNodeA = GetNodeOrient(*endNodeA);
    orientEndNodeA->flip();
    *orientEndNodeB = GetNodeOrient(*endNodeB);
    AGapTowardAEnd = FALSE;
    BGapTowardAEnd = FALSE;
  } else {//curEdge->orient.isBA_BA()
    *endNodeA = GetGraphNode(ScaffoldGraph->ContigGraph, scaffoldA->info.Scaffold.AEndCI);
    *endNodeB = GetTrueBEndCI(ScaffoldGraph, scaffoldB);
    *orientEndNodeA = GetNodeOrient(*endNodeA);
    *orientEndNodeB = GetNodeOrient(*endNodeB);
    orientEndNodeA->flip();
    orientEndNodeB->flip();
    AGapTowardAEnd = FALSE;
    BGapTowardAEnd = TRUE;
  }
  if (orientEndNodeA->isForward()) {
    if (orientEndNodeB->isForward())
      edgeEndsOrient->setIsAB_AB();
    else
      edgeEndsOrient->setIsAB_BA();
  } else {
    if (orientEndNodeB->isForward())
      edgeEndsOrient->setIsBA_AB();
    else
      edgeEndsOrient->setIsBA_BA();
  }
  nextNodeA = GetGraphNode(ScaffoldGraph->ContigGraph, (AGapTowardAEnd?(*endNodeA)->AEndNext:(*endNodeA)->BEndNext));
  nextNodeB = GetGraphNode(ScaffoldGraph->ContigGraph, (BGapTowardAEnd?(*endNodeB)->AEndNext:(*endNodeB)->BEndNext));
  if (nextNodeA) {
    if (AGapTowardAEnd) {
      *extremalGapSizeA = MIN((*endNodeA)->offsetAEnd.mean, (*endNodeA)->offsetBEnd.mean) -
        MAX(nextNodeA->offsetAEnd.mean, nextNodeA->offsetBEnd.mean);
    } else {
      *extremalGapSizeA =
        MIN(nextNodeA->offsetAEnd.mean, nextNodeA->offsetBEnd.mean)-
        MAX((*endNodeA)->offsetAEnd.mean, (*endNodeA)->offsetBEnd.mean);
    }
  } else {
    *extremalGapSizeA = 1000000.0; // large number, anything will fit
  }
  if (nextNodeB) {
    if (BGapTowardAEnd) {
      *extremalGapSizeB = MIN((*endNodeB)->offsetAEnd.mean, (*endNodeB)->offsetBEnd.mean) -
        MAX(nextNodeB->offsetAEnd.mean, nextNodeB->offsetBEnd.mean);
    } else {
      *extremalGapSizeB =
        MIN(nextNodeB->offsetAEnd.mean, nextNodeB->offsetBEnd.mean)-
        MAX((*endNodeB)->offsetAEnd.mean, (*endNodeB)->offsetBEnd.mean);
    }
  } else {
    *extremalGapSizeB = 1000000.0; // large number, anything will fit
  }

}



static
void
MarkScaffoldsForMerging(SEdgeT *curEdge, int markForMerging) {
  CIScaffoldT *scaffoldA = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idA);
  CIScaffoldT *scaffoldB = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idB);

  if (markForMerging) {
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
}


static
void
MarkScaffoldForNotMerging(CIScaffoldT *scaffoldA) {
  scaffoldA->numEssentialB = 1;
  scaffoldA->numEssentialA = 1;
}



static
int
DoesScaffoldCFit(CIScaffoldT *scaffoldA,
                 CIScaffoldT *scaffoldB,
                 CIScaffoldT *scaffoldC,
                 SEdgeT *edgeAB,
                 SEdgeT *edgeToC,
                 VA_TYPE(PtrT) *edgesFromA,
                 VA_TYPE(PtrT) *edgesFromB,
                 int32 *overlapFound,
                 int checkForTinyScaffolds,
                 int32 verbose) {
  int isEdgeFromA = (edgeToC->idA == scaffoldA->id) || (edgeToC->idB == scaffoldA->id);
  double chiSquaredValue;
  LengthT gap;
  LengthT alternateGap;
  SEdgeT *otherEdgeToC = NULL;
  CIScaffoldT *otherScaffold = NULL;
  int i;
  VA_TYPE(PtrT) *otherEdgeList = NULL;
  EdgeCGW_T *overlapEdge = NULL;
  int needNotOverlap = FALSE;
  double lengthC_to_dist;
  *overlapFound = FALSE;


  if (verbose) {
    fprintf(stderr,"* DoesScaffold " F_CID " fit of length %f fit between scaffolds (" F_CID "," F_CID ") in a gap of size %f +/- %f?\n",
            scaffoldC->id, scaffoldC->bpLength.mean, scaffoldA->id, scaffoldB->id, edgeAB->distance.mean, sqrt(edgeAB->distance.variance));
    DumpCIScaffold(stderr, ScaffoldGraph, scaffoldC, FALSE);
  }


  if (isEdgeFromA) {
    otherEdgeList = edgesFromB;
    otherScaffold = scaffoldA;
  } else {
    otherEdgeList = edgesFromA;
    otherScaffold = scaffoldB;
  }


  // We don't want to stick a teeny tiny element in the middle of a gap between two
  // giants.  Often this will later prevent the giants from merging.  So, we
  // prevent it.  We only place stuff such that the size of the element being placed
  // is at least 20% of the size of the implied gap.  So, for a 50k gap, we need a 10k
  // scaffold.  For a 10k gap, we need a 2k scaffold, etc.
  lengthC_to_dist = 1.0;
  if (edgeToC->distance.mean > 0) {
    lengthC_to_dist = (scaffoldC->bpLength.mean / edgeToC->distance.mean);
  }

  if (checkForTinyScaffolds &&
      scaffoldC->bpLength.mean < 5000 &&
      lengthC_to_dist < 0.20) {
    if (verbose) {
      fprintf(stderr,"Scaffold C is too short (%g) relative to edge length (%g)\n",
              scaffoldC->bpLength.mean, edgeToC->distance.mean);

    }
    return FALSE;
  }

  // We compute the gap size for two placements
  // Primary Placement 1:
  //      A                     C                       B
  // <===========       ================>        ==============>
  //              <-AC-->                <--gap-->
  //              <------------AB---------------->
  //
  //  gap.mean = AB.mean - C.length.mean - AC.mean

  // Primary Placement 2:
  //      A                     C                       B
  // <===========       ================>        ==============>
  //              <-gap>                <--CB -->
  //              <------------AB---------------->
  //
  //  gap.mean = AB.mean - C.length.mean - CB.mean

  gap.mean = edgeAB->distance.mean - edgeToC->distance.mean - scaffoldC->bpLength.mean;
  gap.variance = edgeAB->distance.variance + edgeToC->distance.variance + scaffoldC->bpLength.variance;

  // Alternate Placement 1:
  //      A                 B                    C
  // <===========       ================>        ==============>
  //              <-AB->                <--gap -->
  //              <------------AC---------------->
  //
  //  gap.mean = AC.mean - B.length.mean - AB.mean


  // Alternate Placement 2:
  //      C                 A                    B
  // <===========       ================>        ==============>
  //              <-gap->                <--AB -->
  //              <------------CB---------------->
  //
  //  gap.mean = AC.mean - B.length.mean - BC.mean

  alternateGap.mean = edgeToC->distance.mean - edgeAB->distance.mean - (isEdgeFromA?scaffoldB->bpLength.mean:scaffoldA->bpLength.mean);
  alternateGap.variance = edgeAB->distance.variance + edgeToC->distance.variance + (isEdgeFromA?scaffoldB->bpLength.variance:scaffoldA->bpLength.variance);

  if (verbose) {
    fprintf(stderr,"Gap = (%f +/- %f)  Alternate = (%f +/- %f)\n",
            gap.mean, sqrt(gap.variance), alternateGap.mean, sqrt(alternateGap.variance));

    fprintf(stderr,"*DoesScaffoldCFit  CLength = %g  gap estimate = %g +/- %g\n",
            scaffoldC->bpLength.mean, gap.mean, gap.variance);
  }

  if (gap.mean - 3.0 * sqrt(gap.variance)< -5000) {
    if (verbose) {
      fprintf(stderr,"* DoesScaffoldCFit fails on gap.mean < -5000  gap = %g\n", gap.mean);
    }
    return FALSE;
  }
  // Look through all of the edges in the list, and try to find one that fits
  //
  if (!edgeToC->flags.bits.isProbablyBogus) { //  we should find one if this flag isn't set
    for(i = 0; i < GetNumPtrTs(otherEdgeList); i++) {
      otherEdgeToC = (SEdgeT *) *GetPtrT(otherEdgeList,i);
      if (otherEdgeToC->idA == scaffoldC->id ||
          otherEdgeToC->idB == scaffoldC->id ) {
        break;
      } else {
        otherEdgeToC = NULL;
      }
    }
    assert(otherEdgeToC != NULL); // we should find one!
  }

  // If we found an edge from A (B) to C, then test its length
  if (otherEdgeToC) {
    if (!PairwiseChiSquare((double)otherEdgeToC->distance.mean,
                           (double)otherEdgeToC->distance.variance,
                           (double)gap.mean,
                           (double)gap.variance, (LengthT *)NULL,
                           &chiSquaredValue, SCAFFOLD_MERGE_CHI2_THRESHHOLD)) { // fails Chi-square
      if (verbose) {
        fprintf(stderr,"* 2 DoesScaffoldCFit fails pairwise chi square %f\n", chiSquaredValue);
      }
      return FALSE;
    }
  }


  // If we got here, either there is no edge between scaffoldC and scaffoldA or scaffoldB, or
  // the edge length is chi-square compatible with the gap estimate.
  // Now we refine further to handle the other cases.

  // If C can't overlap with A or B, we are done
  if ((gap.mean - 3.0 * sqrt(gap.variance)) > - CGW_MISSED_OVERLAP) {
    if (verbose) {
      fprintf(stderr,"* gap %g +/- %g ...can't overlap... returning TRUE\n", gap.mean, gap.variance);
    }
    return TRUE;
  } else {
    if (verbose)fprintf(stderr,"* Checking overlap\n");
  }

  // C need not overlap with B
  if (PairwiseChiSquare((double)gap.mean,
                        (double)gap.variance,
                        (double)-CGW_MISSED_OVERLAP,
                        (double)1.0, (LengthT *)NULL,
                        &chiSquaredValue, SCAFFOLD_MERGE_CHI2_THRESHHOLD)) { // passes Chi-square
    if (verbose)
      fprintf(stderr,"* gap.mean is chi-squ compatible with  %g +/- %g ... returning TRUE\n", (double)-CGW_MISSED_OVERLAP, 1.0);
    needNotOverlap = TRUE;
    // If this is a one sided edge, see if the alternate placement is also compatible with
    // a no overlap scenario.  If so, bail.
    if (edgeToC->flags.bits.isProbablyBogus) {
      if (PairwiseChiSquare((double)alternateGap.mean,
                            (double)alternateGap.variance,
                            (double)-CGW_MISSED_OVERLAP,
                            (double)1.0, (LengthT *)NULL,
                            &chiSquaredValue, (double)SCAFFOLD_MERGE_CHI2_THRESHHOLD)) { // passes Chi-square
        if (verbose)
          fprintf(stderr,"* One sided edge failed alternate chi square test\n");
        return FALSE;
      }
    }

  } else {
    needNotOverlap = FALSE;
    if (verbose)
      fprintf(stderr,"* gap.mean is NOT chi-squ compatible with  %g +/- %g ... \n", (double)-CGW_MISSED_OVERLAP, 1.0);
  }


  // C MUST overlap
  // look for overlap
  {
    SequenceOrient orientEndNodeA, orientEndNodeB;
    PairOrient edgeEndsOrient;
    NodeCGW_T *endNodeA, *endNodeB;
    double aGapSize, bGapSize;
    int alternate;

    TranslateScaffoldOverlapToContigOverlap(otherScaffold, scaffoldC, edgeToC->orient, &endNodeA, &endNodeB, &orientEndNodeA, &orientEndNodeB, &edgeEndsOrient,
                                            &aGapSize, &bGapSize);

    if (verbose)
      fprintf(stderr,"* Looking for overlap nodeA:" F_CID " nodeB: " F_CID ", endAOrient:%c endBOrient:%c orient:%c distance:%g\n",
              endNodeA->id, endNodeB->id, orientEndNodeA.toLetter(), orientEndNodeB.toLetter(), edgeEndsOrient.toLetter(), edgeToC->distance.mean);

    overlapEdge = FindOverlapEdgeChiSquare(ScaffoldGraph, endNodeA, endNodeB->id,
                                           edgeEndsOrient, edgeToC->distance.mean,
                                           edgeToC->distance.variance, &chiSquaredValue,
                                           (double)SCAFFOLD_MERGE_CHI2_THRESHHOLD, &alternate, verbose);

    // assert(!alternate); // shouldn't get the wrong orientation, should we? What about interleaving?
    if (alternate)
      fprintf( stderr, "Warning: got an alternate edge orientation in DoesScaffoldCFit!!!!\n");
  }
  if (!overlapEdge) {
    if (verbose)
      fprintf(stderr,"* 3 DoesScaffoldCFit fails %s overlap and doesn't\n",(needNotOverlap?"need not ":" must "));
    if (!needNotOverlap) {
      SaveEdgeMeanForLater(edgeToC);
      edgeToC->distance.mean = - CGW_MISSED_OVERLAP;
    }
    return needNotOverlap;//  C Doesn't Fit
  }
  *overlapFound = TRUE;
  if (verbose)
    PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph, " Overlap Found! ", overlapEdge, overlapEdge->idA);

  // If overlap is X-square with gap ==> it fits
  if (PairwiseChiSquare(gap.mean,
                        gap.variance,
                        overlapEdge->distance.mean,
                        overlapEdge->distance.variance,
                        (LengthT *)NULL,
                        &chiSquaredValue, (double)SCAFFOLD_MERGE_CHI2_THRESHHOLD)) { // passes Chi-square

    SaveEdgeMeanForLater(edgeToC);
    edgeToC->distance.mean = overlapEdge->distance.mean;

    if (verbose)
      fprintf(stderr,"* Overlap is chi-sq OK...returning TRUE\n");
    return TRUE;
  } else {
    if (verbose)fprintf(stderr,"* Overlap is NOT chi-sq OK.\n");
    return needNotOverlap;
  }

  return FALSE;
}


static CDS_CID_t thisID = NULLINDEX;

static
int
CompareSEdgesByOtherID(const void *c1, const void *c2) {
  SEdgeT *s1 = *(SEdgeT **)c1;
  SEdgeT *s2 = *(SEdgeT **)c2;
  CDS_CID_t other1 = (s1->idA == thisID?s1->idB:s1->idA);
  CDS_CID_t other2 = (s2->idA == thisID?s2->idB:s2->idA);
  int32 diff = other2 - other1;

  if (diff > 0)
    return (int)1;

  if (diff < 0)
    return -1;

  diff = s1->edgesContributing - s2->edgesContributing;

  if (diff > 0)
    return (int)1;

  if (diff < 0)
    return -1;

  return 0;
}


static
void
SortSEdgesByOtherID(VA_TYPE(PtrT) *edges, CDS_CID_t ThisID) {
  thisID = ThisID;
  qsort((void *)GetPtrT(edges,0),GetNumPtrTs(edges),
        sizeof(PtrT), CompareSEdgesByOtherID);
}



static
void
ZeroEdgeWeights(VA_TYPE(PtrT) *edges) {
  for(int i=0; i<GetNumPtrTs(edges); i++) {
    SEdgeT *edge = (SEdgeT *) *GetPtrT(edges,i);
    edge->quality = 0.0;
    edge->flags.bits.isBogus = FALSE;
    edge->flags.bits.isProbablyBogus = FALSE;
  }
}


static
int
AssignEdgeWeights(VA_TYPE(PtrT) *mergedEdges,
                  SEdgeT *curEdge,
                  CIScaffoldT *scaffoldA,
                  VA_TYPE(PtrT) *edgesFromA,
                  CIScaffoldT *scaffoldB,
                  VA_TYPE(PtrT) *edgesFromB,
                  int doInterleaving,
                  int verbose) {
  SEdgeT *edge;  // Edge from A->C
  int i;

  if (verbose)
    fprintf(stderr,"* Assign Edge Weights\n");

  // Iterate through edgesFromA. (A,C)  These are sorted by increasing otherID and then by decreasing weight
  if (verbose)
    fprintf(stderr,"* Looking through %d edgesFromA\n",
            (int) GetNumPtrTs(edgesFromA));

  for(i = 0; i < GetNumPtrTs(edgesFromA); i++) {
    int j;
    CDS_CID_t otherID;
    int foundBC = FALSE;
    edge = (SEdgeT *) *GetPtrT(edgesFromA,i);
    if (edge->idA == scaffoldA->id) {
      otherID = edge->idB;
    } else {
      otherID = edge->idA;
    }
    if (!doInterleaving && isLargeOverlap(edge))
      continue;

    if (verbose) {
      fprintf(stderr,"* Found an edge AC (" F_CID "," F_CID ",%c)\n",
              scaffoldA->id, otherID, GetEdgeOrientationWRT(edge, scaffoldA->id).toLetter());


      fprintf(stderr,"* Looking through %d edgesFromB\n",
              (int) GetNumPtrTs(edgesFromB));
    }
    // For each edge, look through edgesFromB for an edge to the same Contig from (B,C)
    // If present, weight = weight(A,C), add the heavier of (A,C) (B,C) to mergedEdges.
    // If absent, add (A,C) to mergedEdges
    for(j = 0; j < GetNumPtrTs(edgesFromB); j++) {
      SEdgeT *edgeB = (SEdgeT *) *GetPtrT(edgesFromB,j);

      if (verbose)
        PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph, "EdgeFromBC? ", edgeB, edgeB->idA);

      if (edgeB->quality > 0.0)
        // This edge is already spoken for
        continue;

      CDS_CID_t otherBID = (edgeB->idA == scaffoldB->id) ? edgeB->idB : edgeB->idA;

      if (otherBID != otherID)
        continue;

      // Make sure we have orientation correct

      // We have (A,B), (A,C) AND (C,B)

      PairOrient   edgeAB = GetEdgeOrientationWRT(curEdge, scaffoldA->id);  //  Orientation of AB
      PairOrient   edgeAC    = GetEdgeOrientationWRT(edge,scaffoldA->id);      //  Orientation of AC
      PairOrient   edgeBC   = GetEdgeOrientationWRT(edgeB, scaffoldB->id);    //  Orientation of BC

      assert(edgeAB.isUnknown() == false);
      assert(edgeAC.isUnknown() == false);
      assert(edgeBC.isUnknown() == false);

      //  This block used to report ("if (verbose)") the continue conditions.  The logging was
      //  mostly useless, stuff like "* 1 Orientation is %c should be BA_AB".  The code was
      //  structured as nested switch statements, each (of the primary four) was about 80 lines.
      //  Removing all that crap revealed a very simple function.

      if (edgeAB.isAB_AB()) {
        if (edgeAC.isAB_BA() && !edgeBC.isBA_AB())
          continue;
        if (edgeAC.isAB_AB() && !edgeBC.isBA_BA())
          continue;
        if (edgeAC.isBA_BA())
          continue;
        if (edgeAC.isBA_AB())
          continue;
      }

      if (edgeAB.isAB_BA()) {
        if (edgeAC.isAB_BA() && !edgeBC.isAB_AB())
          continue;
        if (edgeAC.isAB_AB() && !edgeBC.isAB_BA())
          continue;
        if (edgeAC.isBA_BA())
          continue;
        if (edgeAC.isBA_AB())
          continue;
      }

      if (edgeAB.isBA_AB()) {
        if (edgeAC.isBA_BA() && !edgeBC.isBA_AB())
          continue;
        if (edgeAC.isBA_AB() && !edgeBC.isBA_BA())
          continue;
        if (edgeAC.isAB_BA())
          continue;
        if (edgeAC.isAB_AB())
          continue;
      }

      if (edgeAB.isBA_BA()) {
        if (edgeAC.isBA_BA() && !edgeBC.isAB_AB())
          continue;
        if (edgeAC.isBA_AB() && !edgeBC.isAB_AB())
          continue;
        if (edgeAC.isAB_AB())
          continue;
        if (edgeAC.isAB_BA())
          continue;
      }

      foundBC = TRUE;

      if (verbose)
        fprintf(stderr,"* Found all three edges (" F_CID "," F_CID ",%c)%d, (" F_CID "," F_CID ",%c)%d AND (" F_CID "," F_CID ",%c)%d\n",
                scaffoldA->id, scaffoldB->id, edgeAB.toLetter(), curEdge->edgesContributing,
                scaffoldA->id, otherID,       edgeAC.toLetter(), edge->edgesContributing,
                otherBID,scaffoldB->id,       edgeBC.toLetter(), edgeB->edgesContributing);

      if (edge->edgesContributing > edgeB->edgesContributing) {
        edge->quality = edge->edgesContributing + edgeB->edgesContributing;
        edgeB->quality = -1.0;
        // If we have a singleton from each side, we accept this
        if (edge->quality +0.1 >= (double)CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD)
          AppendPtrT(mergedEdges, (void **) &edge);

      } else {
        edgeB->quality = edge->edgesContributing + edgeB->edgesContributing;
        edge->quality = -1.0;
        if (edge->quality + 0.1 >= (double)CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD)
          AppendPtrT(mergedEdges, (void **) &edge);
        
        AppendPtrT(mergedEdges, (void **)&edgeB);
      }
    }  //  End of for j loop


    if (foundBC == FALSE) {
      if (verbose)
        fprintf(stderr,"* Found %sconfirmed AC edge only\n", (edge->edgesContributing < CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD?"UN":""));
      if (edge->edgesContributing < CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD)
        continue;

      //AppendPtrT(mergedEdges, (void **)&edge);
      AppendPtrT(mergedEdges, (void **) &edge);

      edge->quality = edge->edgesContributing;
      edge->flags.bits.isProbablyBogus = TRUE; // overloaded for one-sided

    }
  }

  // Repeat, iterating over edgesFromB
  // If an edge has already been weighted (weight != 0), skip it

  if (verbose)
    fprintf(stderr,"* Iterating over %d edges from B\n",
            (int) GetNumPtrTs(edgesFromB));
  for(i = 0; i < GetNumPtrTs(edgesFromB); i++) {
    CDS_CID_t otherID;
    edge = (SEdgeT *) *GetPtrT(edgesFromB,i);
    // We already dealt with this from the other direction
    if (edge->quality != 0.0)
      continue;
    if (!doInterleaving && isLargeOverlap(edge))
      continue;
    if (edge->idA == scaffoldB->id) {
      otherID = edge->idB;
    } else {
      otherID = edge->idA;
    }
    if (edge->edgesContributing < CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD)
      continue;
    if (verbose)
      fprintf(stderr,"* Found a BC edge (" F_CID "," F_CID ",%c)\n",
              scaffoldB->id, otherID, GetEdgeOrientationWRT(edge, scaffoldB->id).toLetter());


    //AppendPtrT(mergedEdges, (void **)&edge);
    AppendPtrT(mergedEdges, (void **) &edge);

    edge->quality = edge->edgesContributing;
    edge->flags.bits.isProbablyBogus = TRUE; // overloaded for one-sided
  }


  return GetNumPtrTs(mergedEdges);
}



static
SEdgeT *
FindMoreAttractiveMergeEdge(SEdgeT *curEdge,
                            CIScaffoldT *scaffoldA,
                            CIScaffoldT *scaffoldB,
                            int32 recDepth,
                            int doInterleaving,
                            int checkForTinyScaffolds,
                            int32 verbose) {
  VA_TYPE(PtrT) *edgesFromA = CreateVA_PtrT(10);
  VA_TYPE(PtrT) *edgesFromB = CreateVA_PtrT(10);
  VA_TYPE(PtrT) *mergedEdges = CreateVA_PtrT(10);
  int32 endA, endB;
  int32 markedA, markedB;

  assert(scaffoldA->id == curEdge->idA || scaffoldA->id == curEdge->idB);
  assert(scaffoldB->id == curEdge->idA || scaffoldB->id == curEdge->idB);

  if (verbose)
    fprintf(stderr,"* FMA (" F_CID "," F_CID ") depth:%d\n", scaffoldA->id, scaffoldB->id, recDepth);

  if (scaffoldA->flags.bits.smoothSeenAlready &&
      scaffoldB->flags.bits.smoothSeenAlready) {

    if (verbose)
      fprintf(stderr,"* We've already seen both scaffolds...returning\n");
    return NULL;
  }
  // Cut off recursion back to these nodes
  scaffoldA->flags.bits.smoothSeenAlready = 1;
  scaffoldB->flags.bits.smoothSeenAlready = 1;
  curEdge->flags.bits.isBogus = TRUE; // We mark all edges visited by the recursion so we don't retrace our steps

  {
    SequenceOrient orientEndNodeA, orientEndNodeB;
    PairOrient edgeEndsOrient;
    NodeCGW_T *endNodeA, *endNodeB;
    double aGapSize, bGapSize;

    if (verbose)
      fprintf(stderr,"* FindMoreAttractiveMergeEdge %d\n", recDepth);
    TranslateScaffoldOverlapToContigOverlap(scaffoldA, scaffoldB, curEdge->orient, &endNodeA, &endNodeB, &orientEndNodeA, &orientEndNodeB, &edgeEndsOrient,
                                            &aGapSize, &bGapSize);
  }


  assert(curEdge->orient.isUnknown() == false);

  if (curEdge->orient.isAB_AB()) {
    endA = B_END;
    endB = A_END;
  }
  if (curEdge->orient.isAB_BA()) {
    endA = B_END;
    endB = B_END;
  }
  if (curEdge->orient.isBA_AB()) {
    endA = A_END;
    endB = A_END;
  }
  if (curEdge->orient.isBA_BA()) {
    endA = A_END;
    endB = B_END;
  }

  // Find all merge candidates incident on scaffoldA
  markedA = FindAllMergeCandidates(edgesFromA, NULL,
                                   scaffoldA, endA,
                                   scaffoldB,  FALSE, 1,
                                   doInterleaving, verbose);
  SortSEdgesByOtherID(edgesFromA, scaffoldA->id);
  ZeroEdgeWeights(edgesFromA);

  // Find all merge candidates incident on scaffoldB
  markedB = FindAllMergeCandidates(edgesFromB, NULL,
                                   scaffoldB, endB,
                                   scaffoldA,  FALSE, 1,
                                   doInterleaving, verbose);
  SortSEdgesByOtherID(edgesFromB, scaffoldB->id);
  ZeroEdgeWeights(edgesFromB);

  if (markedA || markedB) {
    if (verbose)
      fprintf(stderr,"* Marked Edges Encountered...returning\n");
    //            Free memory
    DeleteVA_PtrT(edgesFromA);
    DeleteVA_PtrT(edgesFromB);
    DeleteVA_PtrT(mergedEdges);

    // Clear marks
    scaffoldA->flags.bits.smoothSeenAlready = 0;
    scaffoldB->flags.bits.smoothSeenAlready = 0;
    return NULL;
  }
  // Nothing to do...no confirmed edge from either scaffold end
  if (GetNumPtrTs(edgesFromA) == 0 &&
      GetNumPtrTs(edgesFromB) == 0) {
    if (verbose) {
      fprintf(stderr,"* No edges from from A OR B!\n");
    }
    // Clear marks
    scaffoldA->flags.bits.smoothSeenAlready = 0;
    scaffoldB->flags.bits.smoothSeenAlready = 0;
    return curEdge;
  }

  // Assign a weight to all edges in the two lists
  // Store it in the overlap quality field
  // If there is a scaffold C s.t. there is a candidate edge A,C and and edge C,B
  // choose the heavier of the two and give it thes core w(A,b) + w(B,C)

  // merge list = listA + listB and sort by decreasing weight




  if (AssignEdgeWeights(mergedEdges, curEdge, scaffoldA,edgesFromA, scaffoldB, edgesFromB, doInterleaving, verbose)) {
    int i;
    //for all edges curEdge in list
    for(i = 0; i < GetNumPtrTs(mergedEdges); i++) {
      // fEdge may be an edge from EITHER scaffoldA OR ScaffoldB to
      // a scaffoldC that is in between A and B
      SEdgeT *fEdge = (SEdgeT *)*GetPtrT(mergedEdges,i);
      CDS_CID_t sourceID, targetID;
      int32 overlapFound;
      int32 isFromA;
      CIScaffoldT *targetScaffold, *sourceScaffold;

      if (verbose) {
        PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph, "AEW", fEdge, fEdge->idA);
        fprintf(stderr," Edge has quality %g\n", fEdge->quality);
      }
      if (fEdge->idA == scaffoldA->id || fEdge->idA == scaffoldB->id) {
        targetID = fEdge->idB;
        sourceID = fEdge->idA;
      } else {
        targetID = fEdge->idA;
        sourceID = fEdge->idB;
      }

      if (sourceID == scaffoldA->id) {
        isFromA = TRUE;
      } else {
        isFromA = FALSE;
      }

      targetScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,targetID);
      sourceScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,sourceID);

      if (verbose)
        PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph, "fedge ", fEdge, fEdge->idA);

      if (DoesScaffoldCFit(scaffoldA, scaffoldB, targetScaffold,curEdge, fEdge, edgesFromA, edgesFromB, &overlapFound, checkForTinyScaffolds, verbose)) {
        SEdgeT *gEdge;

        if (verbose) {
          fprintf(stderr,"*### Scaffold " F_CID " fits between (" F_CID "," F_CID ")  curEdge = (" F_CID "," F_CID ",%c)   fEdge = (" F_CID "," F_CID ",%c)\n",
                  targetID, scaffoldA->id, scaffoldB->id,
                  curEdge->idA, curEdge->idB, curEdge->orient.toLetter(),
                  fEdge->idA, fEdge->idB, fEdge->orient.toLetter());

          fprintf(stderr,"* Calling FindMoreAttractive recursively on " F_CID " " F_CID "\n", sourceScaffold->id, targetScaffold->id);
        }

        gEdge = FindMoreAttractiveMergeEdge(fEdge,sourceScaffold, targetScaffold, recDepth + 1, doInterleaving, checkForTinyScaffolds, verbose);


        if (gEdge != NULL) {
          if (verbose) {
            PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph, "FMA returned gedge ", gEdge,gEdge->idA);
          }
          MarkScaffoldsForMerging(fEdge, FALSE);
          MarkScaffoldForNotMerging(targetScaffold);
        } else {
          if (verbose)
            fprintf(stderr,"*gedge = NULL returning fEdge\n");
          gEdge = fEdge;
        }

        //            Free memory
        DeleteVA_PtrT(edgesFromA);
        DeleteVA_PtrT(edgesFromB);
        DeleteVA_PtrT(mergedEdges);
        // Clear marks
        scaffoldA->flags.bits.smoothSeenAlready = 0;
        scaffoldB->flags.bits.smoothSeenAlready = 0;

        return gEdge;
      } else {
        if (verbose)
          fprintf(stderr,"*XXX Scaffold " F_CID " DOES NOT fit between (" F_CID "," F_CID ")  curEdge = (" F_CID "," F_CID ",%c)   fEdge = (" F_CID "," F_CID ",%c)\n",
                  targetID, scaffoldA->id, scaffoldB->id,
                  curEdge->idA, curEdge->idB, curEdge->orient.toLetter(),
                  fEdge->idA, fEdge->idB, fEdge->orient.toLetter());

        //MarkScaffoldForNotMerging(targetScaffold);
      }

    }

  } else {
    if (verbose) {
      fprintf(stderr,"No edges found\n");
    }
  }

  // Clear marks
  scaffoldA->flags.bits.smoothSeenAlready = 0;
  scaffoldB->flags.bits.smoothSeenAlready = 0;

  //            Free memory
  DeleteVA_PtrT(edgesFromA);
  DeleteVA_PtrT(edgesFromB);
  DeleteVA_PtrT(mergedEdges);

  // If we didn't find anything
  return curEdge;
}



static
int
isQualityScaffoldMergingEdge(SEdgeT                     *curEdge,
                             CIScaffoldT                *scaffoldA,
                             CIScaffoldT                *scaffoldB,
                             ScaffoldInstrumenter       *si,
                             VA_TYPE(MateInstrumenterP) *MIs,
                             double                       minSatisfied,
                             double                       maxDelta) {

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

  int32   mBeforeGood = GetMateStatsHappy(&matesBefore.intra) + GetMateStatsHappy(&matesBefore.inter);
  int32   mBeforeBad  = GetMateStatsBad(&matesBefore.intra)   + GetMateStatsBad(&matesBefore.inter);

  int32   mAfterGood  = GetMateStatsHappy(&matesAfter.intra)  + GetMateStatsHappy(&matesAfter.inter);
  int32   mAfterBad   = GetMateStatsBad(&matesAfter.intra)    + GetMateStatsBad(&matesAfter.inter);

  //  Add in mates that should have been satisfied, but weren't.
  mAfterBad += GetMateStatsMissing(&matesAfter.inter);

  //  This should only be counted for 'inter' (== inter-contig?) and not for 'intra'.
  assert(GetMateStatsMissing(&matesAfter.intra) == 0);

  // since we expect some set of mates to be missing due to divergence between closely related species, we don't perform this check the same way for metagenomics
  if (GetMateStatsMissing(&matesAfter.inter) > 0) {
	  if (((double)GetMateStatsMissing(&matesAfter.inter) / mAfterGood) < GlobalData->mergeScaffoldMissingMates || GlobalData->mergeScaffoldMissingMates == -1) {
	fprintf(stderr, "DOWNCOUNDING THE MISSING MATES BY %d\n", GetMateStatsMissing(&matesAfter.inter));
		  mAfterBad -= GetMateStatsMissing(&matesAfter.inter);
	  }
  }

  int32   mBeforeSum  = mBeforeGood + mBeforeBad;
  int32   mAfterSum   = mAfterGood  + mAfterBad;

  if (mBeforeSum > 0)
    fractMatesHappyBefore = (double)mBeforeGood / mBeforeSum;

  if (mAfterSum > 0)
    fractMatesHappyAfter  = (double)mAfterGood / mAfterSum;

  fprintf(stderr, "isQualityScaffoldMergingEdge()--   before: %.3f satisfied (%d/%d good/bad mates)  after: %.3f satisfied (%d/%d good/bad mates; bad missing %d)\n",
          fractMatesHappyBefore, mBeforeGood, mBeforeBad,
          fractMatesHappyAfter,  mAfterGood,  mAfterBad,
          GetMateStatsMissing(&matesAfter.inter));

  if ((maxDelta > 0.0) &&
      (fractMatesHappyBefore - fractMatesHappyAfter > maxDelta)) {
    fprintf(stderr, "isQualityScaffoldMergingEdge()--   satisfied dropped by too much (%.3f > %.3f) - won't merge\n",
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
#define MAX_FRAC_BAD_TO_GOOD .3

  double badGoodRatio      = 1.0;

  if (mAfterGood > mBeforeGood)
    badGoodRatio = (double)(mAfterBad - mBeforeBad) / (double)(mAfterGood - mBeforeGood);

  bool  failsMinimum       = (fractMatesHappyAfter < minSatisfied);
  bool  failsToGetHappier1 = (fractMatesHappyAfter < fractMatesHappyBefore);
  bool  failsToGetHappier2 = (mAfterGood < mBeforeGood) || (badGoodRatio > MAX_FRAC_BAD_TO_GOOD);

  if (failsMinimum && failsToGetHappier1 && failsToGetHappier2) {
    fprintf(stderr, "isQualityScaffoldMergingEdge()--   not happy enough to merge (%.3f < %.3f) && (%.3f < %.3f) && ((%d < %d) || (%0.3f > %.3f)) - won't merge\n",
            fractMatesHappyAfter, minSatisfied,
            fractMatesHappyAfter, fractMatesHappyBefore,
            mAfterGood, mBeforeGood, badGoodRatio, MAX_FRAC_BAD_TO_GOOD);
    return(FALSE);
  }

  return(TRUE);
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


/*
  After noticing with interleaved scaffold merging that there can be
  some rather high-weight negative edges that are slightly beyond the
  limit of abutting by the chi-squared check, we wanted to let them
  be abutted.
  Two conditions must be satisfied:
  1. The number of stddev's from a -20 edge per unit of edge weight
  is less than .5
  2. The length of the edge is a small fraction of the shorter of the
  two scaffolds involved.
  The specific numbers are based on an examination of edges with
  weight >= 25 in the rat assembly 12/16/2002.

  After noticing that this was still not aggressive enough during the
  Macaque assembly we added another condition. We also hope that this
  will allow ECR to merge the abbutted contigs.
  3. If the overlap is < 2kbp and the overlap is less than 1/2 of
  the shorter scaffold then abut.
  Granger 8/22/05.
*/

#define STDDEVS_PER_WEIGHT_THRESHOLD                 0.5
#define EDGE_PER_MIN_SCAFFOLD_LENGTH_THRESHOLD       0.002

#define MAX_OVERLAP_TO_ABUT                          2000
#define MAX_PERC_SCAFFOLD_LEN                        0.5

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
  if (PairwiseChiSquare(-CGW_MISSED_OVERLAP, curEdge->distance.variance, curEdge->distance.mean, curEdge->distance.variance,
                        (LengthT *) NULL, &chiSquaredValue, PAIRWISECHI2THRESHOLD_CGW) ||
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
                                     iSpec->maxDelta)) {
      isAbuttable = TRUE;
    }
    curEdge->distance.mean = originalMean;
  }

  if (!isAbuttable)
    SaveBadScaffoldMergeEdge(curEdge, iSpec->badSEdges);

  return isAbuttable;
}


static
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


//ExamineUsableSEdges
// This function examines the candidate SEdges in order of weight and decides
// which edges to use for scaffold merges.  We currently do not do interleaved merging,
// but we do try to merge the scaffolds in an order that reduces the need for interleaved
// merging.  A proposed merge can be rejected if:
//     1) It appears to have too long an overlap and is therefor probably an interleaved merge
//     2) One of the scaffold ends involved has already been dedicated to a merge
//
static
void
ExamineSEdgeForUsability(VA_TYPE(PtrT) * sEdges,
                         SEdgeT * curEdge, InterleavingSpec * iSpec,
                         int verbose) {
  SEdgeT *mergeEdge = NULL;
  CIScaffoldT *scaffoldA, *scaffoldB;
  double mergeDistance,maxMergeDistance, minMergeDistance;
  int mayOverlap, mustOverlap;
  double length_to_dist;
  double min_scaffold_length;

  if (CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD > curEdge->edgesContributing) {
    fprintf(stderr, "ExamineSEdgeForUsability()-- to few edges %d < %d\n",
            curEdge->edgesContributing, CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD);
    return;
  }

  if (GlobalData->doInterleavedScaffoldMerging == FALSE && isLargeOverlap(curEdge)) {   /// too big an overlap
    MarkScaffoldsForMerging(curEdge, FALSE  /* mark for simply exclusion from other merges */);
    fprintf(stderr, "ExamineSEdgeForUsability()-- overlap too big\n");
    return;
  }

  scaffoldA = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idA);
  scaffoldB = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idB);
  if (verbose)
    PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph,
                   "S ", curEdge, curEdge->idA);

  // We don't want to stick a teeny tiny element in the middle of a gap between two
  // giants.  Often this will later prevent the giants from merging.  So, we
  // prevent it.  We only place stuff such that the size of the element being placed
  // is at least 20% of the size of the implied gap.  So, for a 50k gap, we need a 10k
  // scaffold.  For a 10k gap, we need a 2k scaffold, etc.
  length_to_dist = 1.0;
  min_scaffold_length = MIN(scaffoldA->bpLength.mean,
                            scaffoldB->bpLength.mean);
  if (curEdge->distance.mean > 0) {
    length_to_dist = ( min_scaffold_length / curEdge->distance.mean);
  }

  if (iSpec->checkForTinyScaffolds &&
      min_scaffold_length < 5000 &&
      length_to_dist < 0.20) {
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
    scaffoldA->flags.bits.walkedAlready =
      scaffoldB->flags.bits.walkedAlready = 1;
    //      fprintf(stderr,"* Skipping edge due to marking!\n");
    return;
  }

  /******	******	******	******	******	******	******	******
   ******	New Recursive Code      ******	******	******	******
   ******	******	******	******	******	******	******	******	*/
  if (TouchesMarkedScaffolds(curEdge)) {
    if (verbose)
      fprintf(stderr,
              "* Edge (" F_CID "," F_CID ",%c) touches marked scaffolds\n",
              curEdge->idA, curEdge->idB, curEdge->orient.toLetter());
    return;
  }


  // See if we already want to merge these two scaffolds, but in an opposite orientation
  // For instance, sometimes there will be evidence suggesting a merge AB_BA and a merge
  // BA_BA between a given pair of scaffolds.
  if ((scaffoldA->essentialEdgeB != NULLINDEX &&
       (scaffoldA->essentialEdgeB == scaffoldB->essentialEdgeA ||
        scaffoldA->essentialEdgeB == scaffoldB->essentialEdgeB )) ||
      (scaffoldA->essentialEdgeA != NULLINDEX &&
       (scaffoldA->essentialEdgeA == scaffoldB->essentialEdgeA ||
        scaffoldA->essentialEdgeA == scaffoldB->essentialEdgeB ))) {

    fprintf(stderr,
            "* We're already trying to merge scaffold (" F_CID "," F_CID ") ...back off!\n",
            scaffoldA->id, scaffoldB->id);
    return;
  }


  if (verbose) {
    fprintf(stderr,
            "* Top Level call to FindMoreAttractiveMergeEdge (" F_CID "," F_CID ",%c) (" F_CID "," F_CID ") gap = %g\n",
            curEdge->idA, curEdge->idB, curEdge->orient.toLetter(),
            scaffoldA->id, scaffoldB->id, curEdge->distance.mean);
  }

  mergeDistance = curEdge->distance.mean;
  minMergeDistance =
    mergeDistance - MAX_SLOP_IN_STD * sqrt(curEdge->distance.variance);
  maxMergeDistance =
    mergeDistance + MAX_SLOP_IN_STD * sqrt(curEdge->distance.variance);

  if (verbose) {
    fprintf(stderr,
            "* curEdge mergeDistance = (%g,%g) min:%g max:%g\n",
            mergeDistance, curEdge->distance.variance,
            minMergeDistance, maxMergeDistance);
  }
  // Look for an overlap

  mayOverlap = (minMergeDistance < CGW_MISSED_OVERLAP &&
                maxMergeDistance > CGW_MISSED_OVERLAP);
  mustOverlap =(minMergeDistance < CGW_MISSED_OVERLAP &&
                maxMergeDistance < CGW_MISSED_OVERLAP);

  // If it is a really heavy edge, treat like a may overlap edge
  if (mustOverlap && curEdge->edgesContributing > EDGE_QUANTA) {
    mustOverlap = 0;
    mayOverlap = 1;
  }


  if (GlobalData->doInterleavedScaffoldMerging == FALSE) {
    if (mustOverlap) {
      // All edges that got this far passed the !isLargeOverlap test, so we will
      // abut them if we can't find the overlap
      mergeEdge = curEdge;

    } else {

      mergeEdge = FindMoreAttractiveMergeEdge(curEdge, scaffoldA, scaffoldB, 0, GlobalData->doInterleavedScaffoldMerging, iSpec->checkForTinyScaffolds, verbose);
      if (verbose) {
        if (!mergeEdge) {
          fprintf(stderr,"*(NONE) Return from Top Level call to FindMoreAttractiveMergeEdge (" F_CID "," F_CID ",%c) (" F_CID "," F_CID ") ==> NONE!\n",
                  curEdge->idA, curEdge->idB, curEdge->orient.toLetter(), scaffoldA->id, scaffoldB->id);
        } else {
          fprintf(stderr,"*(%s) Return from Top Level call to FindMoreAttractiveMergeEdge (" F_CID "," F_CID ",%c) (" F_CID "," F_CID ") ==> merge (" F_CID "," F_CID ",%c)\n",
                  (curEdge == mergeEdge?"N":"Y"),
                  curEdge->idA, curEdge->idB, curEdge->orient.toLetter(), scaffoldA->id, scaffoldB->id,
                  mergeEdge->idA, mergeEdge->idB, mergeEdge->orient.toLetter());
        }
      }

      // If we backed off on the merge, continue
      if (mergeEdge == NULL) {
        // Should we propagate the marks to prevent anything from happening in this neighborhood?
        MarkScaffoldsForMerging(curEdge,  FALSE  /* mark for simply exclusion from other merges */);
        return;
      }

      mergeDistance = mergeEdge->distance.mean;
      minMergeDistance = mergeDistance - MAX_SLOP_IN_STD * sqrt(mergeEdge->distance.variance);
      maxMergeDistance = mergeDistance + MAX_SLOP_IN_STD * sqrt(mergeEdge->distance.variance);

      if (verbose) {
        fprintf(stderr,"* mergeDistance = (%g,%g) min:%g max:%g\n",
                mergeDistance, mergeEdge->distance.variance,
                minMergeDistance, maxMergeDistance);
      }
      // Look for an overlap

      mayOverlap = (minMergeDistance < CGW_MISSED_OVERLAP &&    maxMergeDistance > CGW_MISSED_OVERLAP);
      mustOverlap =(minMergeDistance < CGW_MISSED_OVERLAP &&    maxMergeDistance < CGW_MISSED_OVERLAP);
      //      fprintf(stderr,"* mayOverlap = %d mustOverlap %d\n", mayOverlap, mustOverlap);


      if (mustOverlap && mergeEdge->edgesContributing > EDGE_QUANTA) {
        mustOverlap = 0;
        mayOverlap = 1;
      }

      // Reset scaffoldA or scaffoldB so scaffoldA is to the left of scaffoldB
      if (scaffoldA->id == mergeEdge->idA || scaffoldA->id == mergeEdge->idB) {
        scaffoldB = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                 (scaffoldA->id == mergeEdge->idA) ?
                                 mergeEdge->idB : mergeEdge->idA);
      } else {
        scaffoldA = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                 (scaffoldB->id == mergeEdge->idA) ?
                                 mergeEdge->idB : mergeEdge->idA);
      }
    }

    if (isBadScaffoldMergeEdge(curEdge, iSpec->badSEdges))
      return;

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

    if (mayOverlap || mustOverlap) { // may overlap
      SequenceOrient orientEndNodeA, orientEndNodeB;
      PairOrient edgeEndsOrient;
      NodeCGW_T *endNodeA, *endNodeB;
      double aGapSize, bGapSize;
      EdgeCGW_T *overlapEdge;
      double chiSquaredValue;
      int alternate;

      TranslateScaffoldOverlapToContigOverlap(scaffoldA, scaffoldB, mergeEdge->orient, &endNodeA, &endNodeB, &orientEndNodeA, &orientEndNodeB, &edgeEndsOrient,
                                              &aGapSize, &bGapSize);

      overlapEdge = FindOverlapEdgeChiSquare(ScaffoldGraph, endNodeA, endNodeB->id,
                                             edgeEndsOrient, mergeEdge->distance.mean,
                                             mergeEdge->distance.variance, &chiSquaredValue,
                                             (double)SCAFFOLD_MERGE_CHI2_THRESHHOLD, &alternate, verbose);
      if (overlapEdge == (EdgeCGW_T *)NULL) {
        double minSizeScaffoldA = scaffoldA->bpLength.mean -
          (MAX_SLOP_IN_STD * sqrt(scaffoldA->bpLength.variance));
        double minSizeScaffoldB = scaffoldB->bpLength.mean -
          (MAX_SLOP_IN_STD * sqrt(scaffoldB->bpLength.variance));
        double maxGapScaffoldA = FindMaxGapInScaffold(ScaffoldGraph, scaffoldA);
        double maxGapScaffoldB = FindMaxGapInScaffold(ScaffoldGraph, scaffoldB);

        if (iSpec->checkAbutting &&
            !AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec))
          return;

        if (verbose) {
          fprintf(stderr, "Could not find overlap " F_CID "(" F_CID ") %c " F_CID "(" F_CID ")\n",
                  scaffoldA->id, endNodeA->id, edgeEndsOrient.toLetter(),
                  scaffoldB->id, endNodeB->id);
        }
        if (mustOverlap ) {

          if ((minSizeScaffoldA > maxGapScaffoldB) &&
              (minSizeScaffoldB > maxGapScaffoldA)) {
            SaveEdgeMeanForLater(mergeEdge);
            mergeEdge->distance.mean = MAX(- CGW_MISSED_OVERLAP, mergeEdge->distance.mean);
            fprintf(stderr,"* Must overlap edge -- can't confirm --abutting..going ahead distance is %g (%g,%g)\n", mergeEdge->distance.mean, minMergeDistance, maxMergeDistance);
            PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph, "  MustOverlap ", mergeEdge, mergeEdge->idA);
          } else {
            SaveEdgeMeanForLater(mergeEdge);
            mergeEdge->distance.mean = MAX(- CGW_MISSED_OVERLAP, mergeEdge->distance.mean);
            fprintf(stderr,"* Must overlap edge -- can't confirm, may be containment--abutting..going ahead distance is %g (%g,%g)\n", mergeEdge->distance.mean, minMergeDistance, maxMergeDistance);
            PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph, "  MustOverlap ", mergeEdge, mergeEdge->idA);
          }

        } else {
          // we don't know what to do, so
          if (mergeEdge->edgesContributing > EDGE_QUANTA ||
              ((minSizeScaffoldA > maxGapScaffoldB) &&
               (minSizeScaffoldB > maxGapScaffoldA)) ) {
            SaveEdgeMeanForLater(mergeEdge);
            mergeEdge->distance.mean = MAX(- CGW_MISSED_OVERLAP, mergeEdge->distance.mean);
            if (verbose) {
              fprintf(stderr,"* May overlap edge -- can't confirm ...going ahead distance is %g (%g,%g)\n", mergeEdge->distance.mean, minMergeDistance, maxMergeDistance);
              PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph, "  MayOverlap ", mergeEdge, curEdge->idA);
            }
          } else {
            SaveEdgeMeanForLater(mergeEdge);
            mergeEdge->distance.mean = MAX(- CGW_MISSED_OVERLAP, mergeEdge->distance.mean);
            if (verbose) {
              fprintf(stderr,"* May overlap edge -- can't confirm --abutting..going ahead distance is %g (%g,%g)\n", mergeEdge->distance.mean, minMergeDistance, maxMergeDistance);
              PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph, "  MustOverlap ", mergeEdge, mergeEdge->idA);
            }
          }
        }
      }else if (alternate) {
        // This means we were expecting the extremal contigs to overlap as follows:
        //    Expected                                      Found
        // -------A---------> 	                            -------A--------->
        //           -------B--------->           -------B--------->
        //           |-------|                    |--------------------------|


        fprintf(stderr,"*** %s:Alternate found: NEED TO FIX THIS OVERLAP (" F_CID "," F_CID ",%c) is really (" F_CID "," F_CID ",%c)...setting to CGW_MISSED_OVERLAP\n",
                __FILE__, endNodeA->id, endNodeB->id, edgeEndsOrient.toLetter(), overlapEdge->idA, overlapEdge->idB, overlapEdge->orient.toLetter());

        // To deal with this case we need to handle all of the ramifications.  once A and B have slid by each other, B may overlap with
        // the contig to the left of A, and A may overlap with the contig to the right of B.  If either of these overlaps fails to exist,
        // we need to update the scaffold positions of these neighbors so that they do not cause a merge multialigns failure.

        // As a first cut, see if the extremal contig will fit in the gap of the other scaffold, so no positions need to be updated
        if (endNodeB->bpLength.mean + overlapEdge->distance.mean < aGapSize &&
            endNodeA->bpLength.mean + overlapEdge->distance.mean < bGapSize) {
          fprintf(stderr,"* We can safely interleave them endalength %g  endblength %g  aGap %g bGap %g\n",
                  endNodeA->bpLength.mean,
                  endNodeB->bpLength.mean,
                  aGapSize, bGapSize);
          SaveEdgeMeanForLater(mergeEdge);
          mergeEdge->distance.mean = -(endNodeA->bpLength.mean + endNodeB->bpLength.mean + overlapEdge->distance.mean);

          fprintf(stderr,"*** Fixed so that overlap is length " F_CID " (%g) + length " F_CID " (%g) + overlap (%g) = %g\n",
                  endNodeA->id,endNodeA->bpLength.mean, endNodeB->id,endNodeB->bpLength.mean,overlapEdge->distance.mean, mergeEdge->distance.mean);
        } else {
          fprintf(stderr,"*** Couldn't easily fix...gaps too small\n");
          if (iSpec->checkAbutting &&
              !AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec))
            return;
          SaveEdgeMeanForLater(mergeEdge);
          mergeEdge->distance.mean = -CGW_MISSED_OVERLAP;
          fprintf(stderr,"* We CAN'T safely interleave them, gaps too small:  endalength %g  endblength %g  aGap %g bGap %g\n",
                  endNodeA->bpLength.mean,
                  endNodeB->bpLength.mean,
                  aGapSize, bGapSize);
        }
      } else {
        fprintf(stderr,"* Confirmed overlap (" F_CID "," F_CID ",%c) overlapEdge:%g   mergeEdge:%g\n",
                overlapEdge->idA, overlapEdge->idB, overlapEdge->orient.toLetter(), overlapEdge->distance.mean, mergeEdge->distance.mean);
        if (overlapEdge->orient != edgeEndsOrient && overlapEdge->idA == endNodeA->id) {
          //  There is no easy way to get this assert.
          //assert(overlapEdge->orient == InvertEdgeOrient(edgeEndsOrient));
          {
            PairOrient B = edgeEndsOrient;
            B.flip();
            assert(overlapEdge->orient == B);
          }
          overlapEdge->distance.mean = endNodeA->bpLength.mean + endNodeB->bpLength.mean + overlapEdge->distance.mean;
        } else {
          mergeEdge->distance.mean = overlapEdge->distance.mean;
        }
        if (fabs(overlapEdge->distance.mean - mergeEdge->distance.mean) > MAX_OVERLAP_SLOP_CGW) {
          SaveEdgeMeanForLater(mergeEdge);
          mergeEdge->distance.mean = overlapEdge->distance.mean;
        }
      }
    }
    MarkScaffoldsForMerging(mergeEdge,  TRUE /* mark for merging, not simply exclusion from other merges */);
  } else {
    ScaffoldAlignmentInterface * sai = iSpec->sai;

    if (minMergeDistance < -1000000. &&
        scaffoldA->bpLength.mean > 1000000. &&
        scaffoldB->bpLength.mean > 1000000.) {
      // edge is too negative given scaffold lengths
      if (verbose)
        fprintf(stderr, "Edge is too negative for scaffold lengths\n");
      return;
    }

    if (isBadScaffoldMergeEdge(curEdge, iSpec->badSEdges)) {
      if (verbose)
        fprintf(stderr, "Edge previously marked as bad for merging.\n");
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
    if (mayOverlap || mustOverlap) {
      // may overlap
      SequenceOrient orientEndNodeA, orientEndNodeB;
      PairOrient edgeEndsOrient;
      NodeCGW_T *endNodeA, *endNodeB;
      double aGapSize, bGapSize;
      EdgeCGW_T *overlapEdge = NULL;
      double chiSquaredValue;
      int alternate = FALSE;

      TranslateScaffoldOverlapToContigOverlap(scaffoldA, scaffoldB,
                                              mergeEdge->orient,
                                              &endNodeA, &endNodeB,
                                              &orientEndNodeA,
                                              &orientEndNodeB,
                                              &edgeEndsOrient,
                                              &aGapSize, &bGapSize);

      if (verbose) {
        fprintf(stderr, "scaffoldA length: %d, elements: %d, end node + end gap: %d\n",
                (int) scaffoldA->bpLength.mean, scaffoldA->info.Scaffold.numElements, (int) (endNodeA->bpLength.mean + aGapSize));
        fprintf(stderr, "scaffoldB length: %d, elements: %d, end node + end gap: %d\n",
                (int) scaffoldB->bpLength.mean, scaffoldB->info.Scaffold.numElements, (int) (endNodeB->bpLength.mean + bGapSize));
        fprintf(stderr, "edge length: %.f, variance: %.f, weight: %d\n",
                mergeEdge->distance.mean, mergeEdge->distance.variance, mergeEdge->edgesContributing);
      }

      // take the easy road if the edge is shorter than end nodes/gaps
      if (endNodeA->bpLength.mean + aGapSize > -minMergeDistance &&
          endNodeB->bpLength.mean + bGapSize > -minMergeDistance) {
        // single contig involved in each scaffold
        overlapEdge = FindOverlapEdgeChiSquare(ScaffoldGraph,
                                               endNodeA, endNodeB->id,
                                               edgeEndsOrient,
                                               mergeEdge->distance.mean,
                                               mergeEdge->distance.variance,
                                               &chiSquaredValue,
                                               SCAFFOLD_MERGE_CHI2_THRESHHOLD,
                                               &alternate,
                                               verbose);
        if (verbose)
          fprintf(stderr, "Small overlap - work with end nodes; End node overlap %sfound\n", (overlapEdge) ? "" : "not ");

        if (overlapEdge != NULL) {
          // overlap edge found. use it
          SaveEdgeMeanForLater(mergeEdge);

          if (edgeEndsOrient!=GetEdgeOrientationWRT(overlapEdge,endNodeA->id)) {
            overlapEdge->distance.mean =
              -(endNodeA->bpLength.mean + endNodeB->bpLength.mean +
                overlapEdge->distance.mean);
          }
          mergeEdge->distance.mean = overlapEdge->distance.mean;

          MarkScaffoldsForMerging(mergeEdge, TRUE);
        } else if (!iSpec->checkAbutting ||
                   AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec)) {
          // no overlap edge, but abutting will work
          SaveEdgeMeanForLater(mergeEdge);
          mergeEdge->distance.mean = MAX(-CGW_MISSED_OVERLAP,
                                         mergeEdge->distance.mean);
          MarkScaffoldsForMerging(mergeEdge, TRUE);
        } else {
          // else don't prevent scaffolds from merging via other edges
          SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
        }
      } else {
        // large negative edge - may involve interleaving
        if (verbose)
          fprintf(stderr,
                  "Large overlap - attempt interleaved merging\n");

        //more contigs may be involved in one or both scaffolds
        if (PopulateScaffoldAlignmentInterface(scaffoldA, scaffoldB,
                                               mergeEdge, sai) == 0) {

          if (verbose)
            fprintf(stderr, "Populated ScaffoldAlignmentInterface - %d segments in possible scaffold overlap.\n", (int) GetNumSegmentsInList(sai->segmentList));

#if 0
          // punt if the edge distance is too negative & there are no overlaps
          if ( (int) GetNumSegmentsInList(sai->segmentList) == 0 &&
               endNodeA->bpLength.mean + aGapSize < -minMergeDistance &&
               endNodeB->bpLength.mean + bGapSize < -minMergeDistance) {
            if (verbose) {
              fprintf(stderr,
                      "Large scaffold overlap with no contig overlaps.\n");
              fprintf(stderr, "Edge not usable for merging.\n");
            }
            SaveBadScaffoldMergeEdge(mergeEdge,
                                     iSpec->badSEdges);
            return;
          }
#endif

          if (verbose)
            fprintf(stderr, "Trying to find overlap with band %d,%d\n", sai->scaffoldA->bandBeg,sai->scaffoldA->bandEnd);

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

          if (sai->segmentList != NULL ||             sai->best >= 0) {
            if (verbose) {
              fprintf(stderr, "Align_Scaffold returned best = %d%s\n", sai->best, (sai->best == 0) ? " (this should mean pure interleaving)" : "");
              fprintf(stderr, "Adjusting scaffold contig positions\n");
            }

            if (sai->segmentList != NULL ||
                (iSpec->checkAbutting &&
                 !AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec))) {
              if (verbose)
                fprintf(stderr, "%d segments in scaffold overlap.\n", GetNumSegmentsInList(sai->segmentList));

              // if there are overlaps or abutting isn't an option
              overlapEdge = MakeScaffoldAlignmentAdjustments(scaffoldA, scaffoldB, mergeEdge, sai);

              if (overlapEdge != NULL) {
                SaveEdgeMeanForLater(mergeEdge);
                mergeEdge->distance.mean = overlapEdge->distance.mean;
                if (verbose)
                  fprintf(stderr, "Made scaffold adjustments for merging.\n");
                MarkScaffoldsForMerging(mergeEdge, TRUE);
              } else {
                // else don't prevent scaffolds from merging via other edges
                if (verbose)
                  fprintf(stderr, "Failed to make scaffold adjustments for merging.\n");
                SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
              }
            } else {
              // no overlaps && abutting will work
              SaveEdgeMeanForLater(mergeEdge);
              mergeEdge->distance.mean = MAX(-CGW_MISSED_OVERLAP, mergeEdge->distance.mean);
              MarkScaffoldsForMerging(mergeEdge, TRUE);
            }
          } else {
            // if here, no overlap or interleaving possible
            if (verbose)
              fprintf(stderr, "Align_Scaffold returned best = %d\n", sai->best);

            /* if edge is not highly negative, abutt. Otherwise abort.
               criterion is if -20 abutting still leaves edge as trusted
            */
            if (!iSpec->checkAbutting ||
                AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec)) {
              SaveEdgeMeanForLater(mergeEdge);
              mergeEdge->distance.mean = MAX(-CGW_MISSED_OVERLAP, mergeEdge->distance.mean);
              if (verbose)
                fprintf(stderr, "Abutting will work.\n");
              MarkScaffoldsForMerging(mergeEdge, TRUE);
            } else {
              // else don't prevent scaffolds from merging via other edges
              // record that this scaffold overlap is bad
              if (verbose)
                fprintf(stderr, "Edge not usable for merging.\n");
              SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
            }
          }
        } else {
          fprintf(stderr, "Failed to populate scaffold alignment interface!\n");
          SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
        }
      }
    } else {
      // if here, edge was non-negative
      MarkScaffoldsForMerging(mergeEdge, TRUE);
    }
  }
}


static
void
ExamineUsableSEdgeSet(VA_TYPE(PtrT) *sEdges,
                      int minWeightThreshold,
                      InterleavingSpec * iSpec,
                      int verbose) {

  for (int i=0; i<GetNumPtrTs(sEdges); i++) {
    SEdgeT *curEdge = *(SEdgeT **)GetPtrT(sEdges,i);

    if (curEdge->edgesContributing < minWeightThreshold)
      continue;

    ExamineSEdgeForUsability(sEdges, curEdge, iSpec, verbose);
  }
}


static
void
ExamineUsableSEdges(VA_TYPE(PtrT) *sEdges,
                    int minWeightThreshold,
                    InterleavingSpec * iSpec,
                    int verbose) {

  if (GetNumPtrTs(sEdges) > 0) {
    SEdgeT **sEdge = (SEdgeT **)GetPtrT(sEdges,0);
    int32    maxWeightEdge = 0;

    //  If not interlaved merging, find the max weight non-negative distance edge

    if (! GlobalData->doInterleavedScaffoldMerging) {
      for (int i=0; i<GetNumVA_PtrT(sEdges); i++) {
        if (sEdge[i]->distance.mean > 0) {
          if (isBadScaffoldMergeEdge(sEdge[i], iSpec->badSEdges))
            continue;

          if (maxWeightEdge < sEdge[i]->edgesContributing) {
            fprintf(stderr, "ExamineUsableSEdges()- maxWeightEdge from "F_S32" to "F_S32" at idx "F_S32" out of "F_S32"\n",
                    maxWeightEdge, sEdge[i]->edgesContributing, i, GetNumVA_PtrT(sEdges));
            maxWeightEdge = sEdge[i]->edgesContributing;
          }
        }
      }
    } else {
      for (int i=0; i<GetNumVA_PtrT(sEdges); i++) {
        if (isBadScaffoldMergeEdge(sEdge[i], iSpec->badSEdges))
          continue;

        if (maxWeightEdge < sEdge[i]->edgesContributing) {
          fprintf(stderr, "ExamineUsableSEdges()- maxWeightEdge from "F_S32" to "F_S32" at idx "F_S32" out of "F_S32"\n",
                  maxWeightEdge, sEdge[i]->edgesContributing, i, GetNumVA_PtrT(sEdges));
          maxWeightEdge = sEdge[i]->edgesContributing;
        }
      }
    }

    minWeightThreshold = MIN(minWeightThreshold, maxWeightEdge / EDGE_WEIGHT_FACTOR);
    minWeightThreshold = MAX(minWeightThreshold, EDGE_WEIGHT_FACTOR);

    fprintf(stderr, "* Considering edges with weight >= %d (maxWeightEdge/%d: %d)\n",
            minWeightThreshold,
            (int)EDGE_WEIGHT_FACTOR,
            (int)(maxWeightEdge/EDGE_WEIGHT_FACTOR));
  }

  // examine the edges
  ExamineUsableSEdgeSet(sEdges, minWeightThreshold, iSpec, verbose);
}



static
int
BuildSEdgesForMerging(ScaffoldGraphT * graph,
                      VA_TYPE(PtrT) ** sEdges,
                      VA_TYPE(PtrT) ** overlapSEdges,
                      double * minWeightThreshold,
                      int adjustThreshold,
                      InterleavingSpec * iSpec,
                      int32 verbose) {
  int firstTime;

  if (*sEdges == NULL) {
    *sEdges = CreateVA_PtrT(100);
    *overlapSEdges = CreateVA_PtrT(100);
    firstTime = TRUE;
  } else {
    ResetVA_PtrT(*sEdges);
    ResetVA_PtrT(*overlapSEdges);
    firstTime = FALSE;
  }

  //  find all merge candidates and create a list of pointers to scaffold edges that may be used for
  //  merging
  //
  BuildUsableSEdges(*sEdges, *overlapSEdges, verbose);

  if (adjustThreshold) {
    if (firstTime == TRUE && GetPtrT(*sEdges, 0) != NULL) {
      SEdgeT ** sEdge = (SEdgeT **)GetPtrT(*sEdges,0);

      if (GlobalData->doInterleavedScaffoldMerging) {
        for(int j=0; j<GetNumVA_PtrT(*sEdges); j++) {
          if (sEdge[j]->distance.mean > 0) {

            if (*minWeightThreshold < sEdge[j]->edgesContributing / EDGE_WEIGHT_FACTOR) {
              fprintf(stderr, "BuildSEdgesForMerging()-- minWeightThreshold from %f to %f at iter %d out of %d\n",
                      *minWeightThreshold,
                      (double)sEdge[j]->edgesContributing / EDGE_WEIGHT_FACTOR,
                      j, GetNumVA_PtrT(*sEdges));
              *minWeightThreshold = sEdge[j]->edgesContributing / EDGE_WEIGHT_FACTOR;
            }
          }
        }
      } else {
        *minWeightThreshold = sEdge[0]->edgesContributing / EDGE_WEIGHT_FACTOR;
      }
      fprintf(stderr, "initially setting minWeightThreshold to %f\n", *minWeightThreshold);
    } else {
      //*minWeightThreshold -= 0.01;  //  AZ suggested
      *minWeightThreshold -= 0.2;  //  original
    }

    *minWeightThreshold = MAX( *minWeightThreshold, EDGE_WEIGHT_FACTOR);
  }

  //  mark scaffolds for merging & associate scaffold ends
  //
  ExamineUsableSEdges(*sEdges, (int) *minWeightThreshold, iSpec, verbose);

  return(0);
}



static
void
RemoveDeadRefsFromSEdge(ScaffoldGraphT * graph, SEdgeT * sEdge) {
  if (!sEdge->flags.bits.isRaw) {
    SEdgeT * raw = sEdge;
    SEdgeT * lastRaw = sEdge;

    while ((raw = GetGraphEdge(graph->ScaffoldGraph, raw->nextRawEdge)) != NULL) {
      // referenceEdge references inducing contig edge
      // topLevelEdge references top contig edge
      CIEdgeT * ciEdge = GetGraphEdge(graph->ContigGraph, raw->referenceEdge);
      if (ciEdge->idA == NULLINDEX || ciEdge->idB == NULLINDEX) {
        lastRaw->nextRawEdge = raw->nextRawEdge;
        raw = lastRaw;
      }
    }
  } else {
    // do what?
  }
}


static
int
MergeScaffolds(InterleavingSpec * iSpec, int32 verbose) {
  int mergedSomething = 0;
  GraphNodeIterator scaffolds;
  CIScaffoldT *thisScaffold;
  CDS_CID_t thisScaffoldID; /* The index of the thisScaffold. We have to be careful about thisScaffold since
                               it is a pointer to a an element of the scaffolds array, which may be reallocated
                               during the loop. */
  CDS_CID_t currentSetID = 0;
  int cntScaffold = 1;
  int cnt;

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);

  while ((thisScaffold = NextGraphNodeIterator(&scaffolds)) != NULL) {
    const LengthT  nullLength = {0.0, 0.0};
    CIScaffoldT   *AendScaffold = NULL;
    CIScaffoldT   *BendScaffold = NULL;
    SEdgeT        *edgeA = NULL;
    SEdgeT        *edgeB = NULL;
    SEdgeT        *edge = NULL;
    CIScaffoldT   *neighbor = NULL;
    CDS_CID_t      neighborID = -1;
    LengthT        currentOffset = nullLength;
    CIScaffoldT    CIScaffold;
    SequenceOrient orientCI;
    CDS_CID_t      newScaffoldID = -1;
    int            numMerged = 0;

    thisScaffoldID = thisScaffold->id;

    //fprintf(stderr,"* Examining scaffold %d " F_CID "\n",
    //        cntScaffold, thisScaffoldID);

    if (thisScaffold->type != REAL_SCAFFOLD) {
      //fprintf(stderr, "Not a REAL_SCAFFOLD\n");
      continue;
    }

    if (thisScaffold->setID != NULLINDEX) {
      // This Scaffold has already been placed in a Scaffold.
      //fprintf(stderr, "Already placed in a scaffold.\n");
      continue;
    }

    if (thisScaffold->essentialEdgeA != NULLINDEX) {
      edgeA        = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeA);
      AendScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                  (edgeA->idA == thisScaffold->id) ? edgeA->idB : edgeA->idA);
    }

    if (thisScaffold->essentialEdgeB != NULLINDEX) {
      edgeB        = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeB);
      BendScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                  (edgeB->idA == thisScaffold->id) ? edgeB->idB : edgeB->idA);
    }

    if ((AendScaffold != NULL) &&
        (BendScaffold != NULL)) {
      // This CI is not a starting point for a Scaffold.
      //fprintf(stderr, "CI not a starting point.\n");
      continue;
    }

    if        (BendScaffold != NULL) {
      orientCI.setIsForward();
      edge       = edgeB;
      neighbor   = BendScaffold;
      neighborID = neighbor->id;

    } else if (AendScaffold != NULL) {
      orientCI.setIsReverse();
      edge       = edgeA;
      neighbor   = AendScaffold;
      neighborID = neighbor->id;

    } else {
      // Singleton Scaffold
      //fprintf(stderr, "singleton scaffold.\n");
      continue;
    }

    //  This is guarding against a failure in GOS III; the very first edge in a scaffold merge
    //  had negative variance, which would trigger an assert in the while() loop below.  In general,
    //  we should be checking all edges used in the merge, not just the first edge.
    //
    if (edge->distance.variance <= 0.0)
      continue;

    //  This is our last chance to abort before we create a new scaffold!

    InitializeScaffold(&CIScaffold, REAL_SCAFFOLD);

    CIScaffold.info.Scaffold.AEndCI = NULLINDEX;
    CIScaffold.info.Scaffold.BEndCI = NULLINDEX;
    CIScaffold.info.Scaffold.numElements = 0;
    CIScaffold.edgeHead = NULLINDEX;
    CIScaffold.bpLength = nullLength;
    thisScaffold->setID = currentSetID;
    newScaffoldID = CIScaffold.id = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);

    CIScaffold.flags.bits.isDead = FALSE;
    CIScaffold.numEssentialA = CIScaffold.numEssentialB = 0;
    CIScaffold.essentialEdgeB = CIScaffold.essentialEdgeA = NULLINDEX;
    CIScaffold.setID = NULLINDEX;
    thisScaffold->flags.bits.isDead = TRUE;  // Mark the old scaffold dead

    AppendGraphNode(ScaffoldGraph->ScaffoldGraph, &CIScaffold);  /* Potential realloc of ScaffoldGraph->ScaffoldGraph->nodes */

    thisScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, thisScaffoldID);
    neighbor     = GetGraphNode(ScaffoldGraph->ScaffoldGraph, neighborID);

    if (verbose) {
      fprintf(stderr,"* START: Inserting scaffold " F_CID " into scaffold " F_CID "\n", thisScaffold->id, newScaffoldID);
    }
    assert(thisScaffold->bpLength.variance >= 0);
    assert(neighbor->bpLength.variance >= 0);

    cnt = 0;

    //AppendVA_CDS_CID_t(deadScaffoldIDs, &(thisScaffold->id));

    /* Potential realloc of ScaffoldGraph->ScaffoldGraph->nodes */

    //  BPW - being the first merge, all we're doing is copying the
    //  original scaffold into the new scaffold.
    //
    InsertScaffoldContentsIntoScaffold(ScaffoldGraph,
                                       newScaffoldID, thisScaffold->id,
                                       orientCI, &currentOffset,
                                       iSpec->contigNow);

    thisScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, thisScaffoldID);
    neighbor     = GetGraphNode(ScaffoldGraph->ScaffoldGraph, neighborID);

    currentOffset = thisScaffold->bpLength;
    numMerged = 1;
    mergedSomething++;
    while (neighbor != (CIScaffoldT *)NULL) {
      PairOrient edgeOrient = GetEdgeOrientationWRT(edge, thisScaffold->id);

      cnt++;

      assert(edge->distance.variance >= 0);

      assert(currentOffset.mean >= 0);
      assert(currentOffset.variance >= 0);

      assert(thisScaffold->bpLength.variance >= 0);
      assert(neighbor->bpLength.variance >= 0);

      neighborID = neighbor->id;

      assert(orientCI.isUnknown() == false);

      if (orientCI.isForward()) {
        if (edgeOrient.isAB_AB()) {
          orientCI.setIsForward();
        }else if (edgeOrient.isAB_BA()) {
          orientCI.setIsReverse();
        } else {
          assert(0);
        }
      } else {
        if (edgeOrient.isBA_AB()) {
          orientCI.setIsForward();
        }else if (edgeOrient.isBA_BA()) {
          orientCI.setIsReverse();
        } else {
          assert(0);
        }
      }
      thisScaffold = neighbor;
      thisScaffoldID = thisScaffold->id;

      assert(thisScaffold->bpLength.variance >= 0);
      assert(edge->distance.variance >= 0);
      assert(currentOffset.mean >= 0);
      assert(currentOffset.variance >= 0);
      assert(orientCI.isUnknown() == false);

      currentOffset.mean += edge->distance.mean;
      currentOffset.variance += edge->distance.variance;

      if (edge->distance.mean < -CGW_MISSED_OVERLAP) {
        CIScaffoldT * newScaffold =
          GetCIScaffoldT(ScaffoldGraph->CIScaffolds, newScaffoldID);
        if (currentOffset.mean < 0.0) {
          /*
            0
            newScaffold:                --------------->
            thisScaffold:      ------------------
          */
          CIScaffoldTIterator CIs;
          ChunkInstanceT * CI;
          double variance = GetVarianceOffset(thisScaffold, -currentOffset.mean, orientCI.isForward());

          fprintf(stderr, "Adjusting newScaffold offsets by mean - %f variance + %f\n", currentOffset.mean, variance);

          InitCIScaffoldTIterator(ScaffoldGraph, newScaffold,
                                  TRUE, FALSE, &CIs);
          while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
            CI->offsetAEnd.mean -= currentOffset.mean;
            CI->offsetAEnd.variance += variance;
            CI->offsetBEnd.mean -= currentOffset.mean;
            CI->offsetBEnd.variance += variance;
          }
          currentOffset.mean = currentOffset.variance = 0.0;
        } else {
          /*
            0
            newScaffold:     --------------->
            thisScaffold:              ------------------
          */
          currentOffset.variance = GetVarianceOffset(newScaffold, currentOffset.mean, TRUE);
        }
      }
      assert(currentOffset.variance >= 0);

      if (verbose) {
        fprintf(stderr,"* Adding to scaffold " F_CID " scaffold " F_CID " at orient %c offset %g\n",
                newScaffoldID, thisScaffoldID, orientCI.toLetter(), currentOffset.mean);
      }
#ifdef CHECK_CONTIG_ORDERS
      AddScaffoldToContigOrientChecker(ScaffoldGraph, thisScaffold, coc);
#endif

      //AppendVA_CDS_CID_t(deadScaffoldIDs, &(thisScaffold->id));

      InsertScaffoldContentsIntoScaffold(ScaffoldGraph,
                                         newScaffoldID, thisScaffold->id,
                                         orientCI, &currentOffset,
                                         iSpec->contigNow);

      if (iSpec->contigNow == TRUE)
        // remove references to edges of dead contigs from scaffold edge
        RemoveDeadRefsFromSEdge(ScaffoldGraph, edge);
      else
        // need to set edge statuses so scaffold is clone-connected
        MarkUnderlyingCIEdgesTrusted(ScaffoldGraph, edge);

      neighbor = GetGraphNode(ScaffoldGraph->ScaffoldGraph, neighborID);

      assert(orientCI.isUnknown() == false);
      assert(edge->distance.variance >= 0);
      assert(thisScaffold->bpLength.variance >= 0);
      assert(currentOffset.variance >= 0);

      thisScaffold->setID = currentSetID;
      thisScaffold->flags.bits.isDead = TRUE; // Mark the old scaffold dead
      currentOffset.mean += thisScaffold->bpLength.mean;
      currentOffset.variance += thisScaffold->bpLength.variance;

      assert(thisScaffold->bpLength.variance >= 0);
      assert(currentOffset.mean >= 0);
      assert(currentOffset.variance >= 0);

      numMerged++;
      if (orientCI.isForward()) {
        if (thisScaffold->essentialEdgeB != NULLINDEX) {
          edge = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeB);
          neighbor = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                  (edge->idA == thisScaffold->id) ? edge->idB : edge->idA);
          neighborID = neighbor->id;
        } else {// End of Scaffold
          edge = (SEdgeT *)NULL;
          neighbor = (CIScaffoldT *)NULL;
          neighborID = NULLINDEX;
        }
      } else {// orientCI == B_A
        if (thisScaffold->essentialEdgeA != NULLINDEX) {
          edge = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeA);
          neighbor = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                  (edge->idA == thisScaffold->id) ? edge->idB : edge->idA);
          neighborID = neighbor->id;
        } else {// End of Scaffold
          edge = (SEdgeT *)NULL;
          neighbor = (CIScaffoldT *)NULL;
          neighborID = NULLINDEX;
        }
      }
    }  //  while (neighbor != (CIScaffoldT *)NULL)

    // New scaffold fully popuated now.

    if (iSpec->contigNow == TRUE &&
        GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                     newScaffoldID)->info.Scaffold.numElements > 1) {

      int32  numScaffoldsBefore = GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds);
      int32  status             = RECOMPUTE_SINGULAR;

      int32 recomputeIteration = 0;

      while (recomputeIteration < 3 &&
             (status == RECOMPUTE_SINGULAR ||
              status == RECOMPUTE_CONTIGGED_CONTAINMENTS)) {
        // need to make sure scaffold is connected with trusted raw edges

        MarkInternalEdgeStatus(ScaffoldGraph,
                               GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                            newScaffoldID),
                               PAIRWISECHI2THRESHOLD_CGW,
                               1000.0 * SLOPPY_EDGE_VARIANCE_THRESHHOLD,
                               TRUE, TRUE, 0, TRUE);

#ifdef CHECKCONNECTED
        assert(IsScaffoldInternallyConnected(ScaffoldGraph,
                                             GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                                          newScaffoldID),
                                             ALL_EDGES));
#endif

        status = RecomputeOffsetsInScaffold(ScaffoldGraph,
                                            newScaffoldID,
                                            TRUE, TRUE, FALSE);
        recomputeIteration++;
      }

      if (status != RECOMPUTE_OK)
        fprintf(stderr, "ReomputeOffsetsInScaffold failed (%d) for scaffold " F_CID " in MergeScaffolds\n", status, newScaffoldID);

      int32  numScaffoldsAfter = GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds);

      //  If numScaffoldsBefore is not the same as numScaffoldsAfter, then RecomputeOffsetsInScaffold() decided
      //  to split a scaffold we just merged.  This (once) led to an infinite loop in scaffold merging where
      //  a two contig and a one contig scaffold were merged, then ROIS() would unmerge them into the original
      //  layout.

      if (status == RECOMPUTE_FAILED_CONTIG_DELETED)
        fprintf(stderr, "MergeScaffolds()-- WARNING:  ROIS failed, contig deleted.\n");

      if         (numScaffoldsAfter - numScaffoldsBefore == numMerged - 1) {
        fprintf(stderr, "MergeScaffolds()-- WARNING:  merged %d scaffolds into 1, then unmerged into %d new scaffolds.  NO MERGING DONE?!\n",
                numMerged, numScaffoldsAfter - numScaffoldsBefore);
        mergedSomething--;

      } else  if (numScaffoldsAfter != numScaffoldsBefore) {
        fprintf(stderr, "MergeScaffolds()-- WARNING:  merged %d scaffolds into 1, then unmerged into %d new scaffolds.  Partial merging done?\n",
                numMerged, numScaffoldsAfter - numScaffoldsBefore);
      }
    }  //  if (iSpec->contigNow == TRUE && .....)

#warning numLiveScaffolds is probably broken, and is nearly unused
    ScaffoldGraph->numLiveScaffolds += (1 - numMerged);

    currentSetID++;
  }

  return mergedSomething;
}


static
void
MergeScaffoldsExhaustively(ScaffoldGraphT * graph,
                           InterleavingSpec * iSpec,
                           char *logicalcheckpointnumber,
                           int verbose) {
  static VA_TYPE(PtrT) *sEdges = NULL;
  static VA_TYPE(PtrT) *overlapSEdges = NULL;

  int32  mergedSomething    = TRUE;
  int32  iterations         = 0;
  double  minWeightThreshold = 0.0;
  time_t lastCkpTime        = time(0) - 90 * 60;

  //  Create a scaffold edge for every inter-scaffold contig edge, then merge compatible ones
  BuildSEdges(graph, TRUE, GlobalData->doInterleavedScaffoldMerging);
  MergeAllGraphEdges(graph->ScaffoldGraph, TRUE, FALSE);

  // loop until nothing gets merged
  while (mergedSomething) {
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

    BuildSEdgesForMerging(graph,
                          &sEdges, &overlapSEdges,
                          &minWeightThreshold, TRUE,
                          iSpec, verbose);

    mergedSomething = MergeScaffolds(iSpec, verbose);

    if (mergedSomething) {
      fprintf(stderr, "MergeScaffoldsAggressive()-- iter %d -- continue because we merged scaffolds.\n",
              iterations);

      //  Cleanup, build new edges, merge.
      CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);
      BuildSEdges(graph, TRUE, GlobalData->doInterleavedScaffoldMerging);
      MergeAllGraphEdges(graph->ScaffoldGraph, TRUE, FALSE);

    } else if (minWeightThreshold > 2.0) {
      fprintf(stderr, "MergeScaffoldsAggressive()-- iter %d -- continue because minWeightThreshold is %f (decrease by %f).\n",
              iterations, minWeightThreshold, MAX(1.0, minWeightThreshold / 100.0));

      //  Do we need to clean up the edges/scaffolds here?

      mergedSomething     = 1;
      minWeightThreshold -= MAX(1.0, minWeightThreshold / 100.0);

    } else {
      fprintf(stderr, "MergeScaffoldsAggressive()-- iter %d -- no additional scaffold merging is possible.\n",
              iterations);
    }

    iterations++;
  }
}



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

  if (GlobalData->doInterleavedScaffoldMerging) {
    fprintf(stderr, "** Merging scaffolds with interleaving.\n");
    fprintf(stderr, "** MinSatisfied: %f, MaxDelta: %f\n", iSpec.minSatisfied, iSpec.maxDelta);

    iSpec.checkForTinyScaffolds = FALSE;
    //LeastSquaresGapEstimates(graph, TRUE, FALSE, TRUE, TRUE, FALSE);
    MergeScaffoldsExhaustively(graph, &iSpec, logicalcheckpointnumber, verbose);

  } else {
    fprintf(stderr, "** Merging scaffolds without interleaving.\n");
    fprintf(stderr, "** MinSatisfied: %f, MaxDelta: %f\n", iSpec.minSatisfied, iSpec.maxDelta);

    iSpec.checkForTinyScaffolds = TRUE;
    MergeScaffoldsExhaustively(graph, &iSpec, logicalcheckpointnumber, verbose);
  }

  DeleteScaffoldAlignmentInterface(iSpec.sai);

  for (int i=0; i<GetNumVA_PtrT(iSpec.MIs); i++)
    safe_free(*(MateInstrumenter **)GetVA_PtrT(iSpec.MIs, i));
  DeleteVA_PtrT(iSpec.MIs);

  DestroyChunkOverlapper(iSpec.badSEdges);

  CheckCIScaffoldTs(ScaffoldGraph);
}

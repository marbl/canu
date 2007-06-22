
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
static char CM_ID[] = "$Id: CIScaffoldT_Merge_CGW.c,v 1.31 2007-06-22 20:25:48 eliv Exp $";


#undef ORIG_MERGE_EDGE_INVERT
#define MINSATISFIED_CUTOFF 0.985
#undef DEBUG_MERGE_EDGE_INVERT	  
#undef  DEBUG_BAD_MATE_RATIO


//  Define this to check (and assert) if the graph is not internally
//  connected before recomputing offsets.  It's expensive, and if you
//  already know it's OK (debugging, maybe??) you can skip it.
#define CHECKCONNECTED

//  Define this to enable more aggressive scaffold abutment rules.
#define AGGRESSIVE_ABUTTING


//  Draw .cam files for bad mates, if they have more than N
//  contributing edges.  '8' should reduce the number of cam files to
//  a manageable number.
//
//#define DRAW_BAD_MATE_CAMS 8

#ifdef DRAW_BAD_MATE_CAMS
extern int do_draw_frags_in_CelamyScaffold;
#endif


//#define INSTRUMENT_CGW
//#define CHECK_CONTIG_ORDERS
#define OTHER_END_CHECK
//#define GENERAL_STRONGER_CHECK

#define SORT_BY_EDGE_WEIGHTS
//#define DEBUG1
//#define DEBUG 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_fileIO.h"
#include "UtilsREZ.h"
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "UnionFind_AS.h"
#include "UtilsREZ.h"
#include "ChiSquareTest_CGW.h"
#include "MultiAlignment_CNS.h"
#include "DataTypesREZ.h"
#include "CommonREZ.h"
#include "Stats_CGW.h"   // for collecting scaffold merging stats
#ifdef INSTRUMENT_CGW
#include "Instrument_CGW.h"
#endif
#ifdef CHECK_CONTIG_ORDERS
#include "AS_CGW_EdgeDiagnostics.h"
#endif
#include "InterleavedMerging.h"

#undef REQUIRE_MORE_THAN_ONE_BAD_TO_REJECT
#define REQUIRE_BAD_APPROACHING_HAPPY_TO_REJECT
#define MAX_FRAC_BAD_TO_GOOD .3

void SaveEdgeMeanForLater(SEdgeT * edge)
{
  if(edge->flags.bits.MeanChangedByWalking == FALSE)
    {
      edge->flags.bits.MeanChangedByWalking = TRUE;
      edge->minDistance = edge->distance.mean;
    }
}

void PrintEdgesBetweenScaffolds(ScaffoldGraphT * graph,
                                CDS_CID_t sid1,
                                CDS_CID_t sid2)
{
  // for dros faster scaffold merging debugging
  ChunkInstanceT * scaffold;
  SEdgeTIterator SEdges;
  SEdgeT * sEdge;
  CDS_CID_t minID = MIN(sid1, sid2);
  CDS_CID_t maxID = MAX(sid1, sid2);
  
  scaffold = GetCIScaffoldT(graph->CIScaffolds, minID);
  if(scaffold != NULL && !scaffold->flags.bits.isDead)
    {
      InitSEdgeTIterator(ScaffoldGraph, scaffold->id,
                         TRUE, FALSE, ALL_END, FALSE, &SEdges);
      while((sEdge = NextSEdgeTIterator(&SEdges)) != NULL)
        {
          if(sEdge->idB == maxID)
            {
              PrintSEdgeT(GlobalData->stderrc, graph, "R ", sEdge, sEdge->idA);
              fprintf(GlobalData->stderrc, "all field: " F_U64 "\n", sEdge->flags.all);
            }
        }
      InitSEdgeTIterator(ScaffoldGraph, scaffold->id,
                         FALSE, FALSE, ALL_END, FALSE, &SEdges);
      while((sEdge = NextSEdgeTIterator(&SEdges)) != NULL)
        {
          if(sEdge->idB == maxID)
            {
              PrintSEdgeT(GlobalData->stderrc, graph, "M ", sEdge, sEdge->idA);
              fprintf(GlobalData->stderrc, "all field: " F_U64 "\n", sEdge->flags.all);
            }
        }
    }
}

#define DUMP_SCAFFOLDS

void DumpScaffoldsToFile(ScaffoldGraphT * graph,
                         char * suffix,
                         int number)
{
  char filename[1024];
  FILE * fp;
  
#ifdef DUMP_SCAFFOLDS
  sprintf(filename, "iteration%d%s", number, suffix);
  fp = fopen(filename, "w");
  assert(fp != NULL);
  DumpCIScaffolds(fp, graph, TRUE);
  fclose(fp);
#endif
}

/****************************************************************************/

#define SCAFFOLD_MERGE_CHI2_THRESHHOLD (2.f*(float)PAIRWISECHI2THRESHOLD_CGW)

int ContigCoordinatesOkay(CIScaffoldT * scaffold)
{
  CIScaffoldTIterator CIs;
  ChunkInstanceT * CI;
  LengthT lastMin = {0,0};
  LengthT lastMax = {0,0};
  int i = 0;
  
  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);
  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL)
    {
      LengthT thisMin;
      LengthT thisMax;
      if(CI->offsetAEnd.mean < 0.0 ||
         CI->offsetAEnd.variance < 0.0 ||
         CI->offsetBEnd.mean < 0.0 ||
         CI->offsetBEnd.variance < 0.0)
        {
          fprintf(GlobalData->stderrc, "Negative contig mean or variance\n");
          DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffold, FALSE);
          return FALSE;
        }
      if(i != 0 &&
         (CI->offsetAEnd.mean == 0.0 ||
          CI->offsetBEnd.mean == 0.0))
        {
          fprintf(GlobalData->stderrc, "Zero offset of internal contig\n");
          DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffold, FALSE);
          return FALSE;
        }
    
      thisMin = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ?
        CI->offsetAEnd : CI->offsetBEnd;
      thisMax = (CI->offsetAEnd.mean > CI->offsetBEnd.mean) ?
        CI->offsetAEnd : CI->offsetBEnd;
    
      if(thisMax.mean < lastMin.mean)
        {
          fprintf(GlobalData->stderrc, "Seriously out of order contigs\n");
          DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffold, FALSE);
          return FALSE;
        }
    
      if(thisMin.mean < lastMin.mean)
        {
          fprintf(GlobalData->stderrc, "Out of order contigs\n");
          DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffold, FALSE);
          return FALSE;
        }
    
      lastMin = thisMin;
      lastMax = thisMax;
    
      i++;
    }
  return TRUE;
}


double GetVarianceOffset(CIScaffoldT * scaffold, double meanOffset, int isAB)
{
  double lengthDelta = (isAB) ? 0.0 : scaffold->bpLength.mean;
  double varDelta = (isAB) ? 0.0 : scaffold->bpLength.variance;
  CIScaffoldTIterator CIs;
  ChunkInstanceT * CI;
  LengthT lowerOffset = {-1, 0};
  LengthT higherOffset = {0, 0};
  
  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, isAB, FALSE, &CIs);
  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL)
    {
      double aEnd = fabs(lengthDelta - CI->offsetAEnd.mean);
      double bEnd = fabs(lengthDelta - CI->offsetBEnd.mean);
      if(aEnd > meanOffset || bEnd > meanOffset)
        {
          if(aEnd < meanOffset)
            {
              lowerOffset = CI->offsetAEnd;
              higherOffset = CI->offsetBEnd;
              break;
            }
          else if(bEnd < meanOffset)
            {
              lowerOffset = CI->offsetBEnd;
              higherOffset = CI->offsetAEnd;
              break;
            }
          else
            {
              higherOffset = (aEnd < bEnd) ? CI->offsetAEnd : CI->offsetBEnd;
              break;
            }
        }
      lowerOffset = (aEnd < bEnd) ? CI->offsetBEnd : CI->offsetAEnd;
    }
  if(higherOffset.mean != 0)
    return fabs(varDelta - lowerOffset.variance) +
      (fabs(varDelta - higherOffset.variance) - fabs(varDelta - lowerOffset.variance)) *
      (meanOffset - fabs(lengthDelta - lowerOffset.mean)) /
      (fabs(lengthDelta - higherOffset.mean) - fabs(lengthDelta - lowerOffset.mean));
  else
    return scaffold->bpLength.variance * meanOffset / scaffold->bpLength.mean;
}


int InsertScaffoldContentsIntoScaffold(ScaffoldGraphT *sgraph,
                                       CDS_CID_t newScaffoldID,
                                       CDS_CID_t oldScaffoldID,
                                       FragOrient orient,
                                       LengthT * offset,
                                       int contigNow){
  CIScaffoldT *oldScaffold = GetCIScaffoldT(sgraph->CIScaffolds, oldScaffoldID);
  CIScaffoldT *newScaffold = GetCIScaffoldT(sgraph->CIScaffolds, newScaffoldID);
  CIScaffoldTIterator CIs;
  ChunkInstanceT *CI;
  
  // CheckCIScaffoldTLength(sgraph, oldScaffold);
  // CheckCIScaffoldTLength(sgraph, newScaffold);
  
#ifdef DEBUG_INSERTSCAFFOLDCONTENTS
  fprintf(GlobalData->stderrc,"* InsertContents of Scaffold " F_CID " into scaffold " F_CID " at offset %g +/- %g orient %c oldScaffold length %g\n",
          oldScaffoldID, newScaffoldID, offset->mean, sqrt(offset->variance), orient, oldScaffold->bpLength.mean);
#endif
  assert(offset->mean >= 0.0 && offset->variance >= 0.0);
  
  InitCIScaffoldTIterator(sgraph, oldScaffold, (orient == A_B), FALSE, &CIs);
  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
    LengthT offsetAEnd, offsetBEnd;

#ifdef DEBUG_INSERTSCAFFOLDCONTENTS
    fprintf(GlobalData->stderrc, "* CI->offsetAEnd=%d  CI->offsetBEnd=%d  oldScaffold->bpLength=%d\n",
            (int)CI->offsetAEnd.mean, (int)CI->offsetBEnd.mean, (int)oldScaffold->bpLength.mean);
#endif

    if(orient == A_B){
      offsetAEnd.mean      = offset->mean     + CI->offsetAEnd.mean;
      offsetAEnd.variance  = offset->variance + CI->offsetAEnd.variance;
      offsetBEnd.mean      = offset->mean     + CI->offsetBEnd.mean;
      offsetBEnd.variance  = offset->variance + CI->offsetBEnd.variance;
    }else{
      if ((CI->offsetAEnd.mean > oldScaffold->bpLength.mean) ||
          (CI->offsetBEnd.mean > oldScaffold->bpLength.mean)){
        fprintf(GlobalData->stderrc, "* CI extends beyond end of scaffold!\n");
        fprintf(GlobalData->stderrc, "* offsetAEnd = %d offsetBEnd = %d scaffoldLength = %d\n",
                (int)CI->offsetAEnd.mean, (int)CI->offsetBEnd.mean, (int)oldScaffold->bpLength.mean);
        assert(0);
      }
      offsetAEnd.mean     = offset->mean     + (oldScaffold->bpLength.mean     - CI->offsetAEnd.mean);
      offsetAEnd.variance = offset->variance + (oldScaffold->bpLength.variance - CI->offsetAEnd.variance);
      offsetBEnd.mean     = offset->mean     + (oldScaffold->bpLength.mean     - CI->offsetBEnd.mean);
      offsetBEnd.variance = offset->variance + (oldScaffold->bpLength.variance - CI->offsetBEnd.variance);

      if(CI->offsetBEnd.variance > oldScaffold->bpLength.variance){
        if(GlobalData->debugLevel > 0)
          fprintf(GlobalData->stderrc,"* CI " F_CID " has BEnd variance %g > scaffold bpLength.variance %g\n",
                  CI->id, CI->offsetBEnd.variance, oldScaffold->bpLength.variance);
        offsetBEnd.variance = offset->variance;
      }
      if(CI->offsetAEnd.variance > oldScaffold->bpLength.variance){
        if(GlobalData->debugLevel > 0)
          fprintf(GlobalData->stderrc,"* CI " F_CID " has AEnd variance %g > scaffold bpLength.variance %g\n",
                  CI->id, CI->offsetAEnd.variance, oldScaffold->bpLength.variance);
        offsetAEnd.variance = offset->variance;
      }
    }

    if ((offsetBEnd.variance < 0.0) || (offsetAEnd.variance < 0.0)) {
      fprintf(stderr, "oldScaffold:\n");
      DumpCIScaffold(GlobalData->stderrc, sgraph, oldScaffold, FALSE);

      fprintf(stderr, "newScaffold:\n");
      DumpCIScaffold(GlobalData->stderrc, sgraph, newScaffold, FALSE);

      fprintf(GlobalData->stderrc,"offsetAEnd mean:%d variance:%g offsetBEnd mean:%d variance:%g < 0...sigh...\n",
              (int)offsetAEnd.mean, offsetAEnd.variance,
              (int)offsetBEnd.mean, offsetBEnd.variance);
      assert(0);
    }

#ifdef DEBUG_INSERTSCAFFOLDCONTENTS
    fprintf(GlobalData->stderrc,"* Inserting CI " F_CID " at offset (%d,%d) was at (%d,%d)\n",
            CI->id,
            (int) offsetAEnd.mean, (int) offsetBEnd.mean,
            (int) CI->offsetAEnd.mean, (int) CI->offsetBEnd.mean);
#endif
    InsertCIInScaffold(sgraph, CI->id, newScaffoldID,
                       offsetAEnd, offsetBEnd,
                       TRUE, contigNow);
#if 1
    if(!ContigCoordinatesOkay(newScaffold))
      {
	fprintf(GlobalData->stderrc, "Problem with contig coordinates!\n");
      }
#endif
  }
  CheckCIScaffoldTLength(sgraph, newScaffold);
  
  return TRUE;
}



void MarkUnderlyingRawCIEdgeTrusted(ScaffoldGraphT * sgraph, EdgeCGW_T * raw)
{
  CIEdgeT *ciedge, *topCIEdge;
  ContigT *contigA, *contigB;
  
  // fprintf(GlobalData->stderrc,"* raw->referenceEdge = " F_CID "\n", raw->referenceEdge);
  ciedge = GetGraphEdge(sgraph->RezGraph,raw->referenceEdge);
  // fprintf(GlobalData->stderrc,"* ciedge->toplevelEdge = " F_CID " numCIEdges %d\n", ciedge->topLevelEdge, (int) GetNumGraphEdges(sgraph->RezGraph));
  topCIEdge = GetGraphEdge(sgraph->RezGraph, ciedge->topLevelEdge);
  assert(topCIEdge->idA != NULLINDEX && topCIEdge->idB != NULLINDEX);
  contigA = GetGraphNode(sgraph->RezGraph, topCIEdge->idA);
  contigB = GetGraphNode(sgraph->RezGraph, topCIEdge->idB);
  assert(contigA->scaffoldID == contigB->scaffoldID);
  assert(contigA->flags.bits.isUnique);
  assert(contigB->flags.bits.isUnique);
  AssertPtr(topCIEdge);
  SetGraphEdgeStatus(sgraph->RezGraph, topCIEdge, TRUSTED_EDGE_STATUS);
  if(GlobalData->debugLevel > 0){
    fprintf(GlobalData->stderrc,"* Marked contig edge " F_CID " (" F_CID "," F_CID ")%c as trusted(inside scaf " F_CID ")\n",
            topCIEdge->topLevelEdge, topCIEdge->idA, topCIEdge->idB, topCIEdge->orient,
            contigA->scaffoldID);
  }
}


/***************************************************************************/
void  MarkUnderlyingCIEdgesTrusted(ScaffoldGraphT *sgraph, SEdgeT *edge){
  if(GlobalData->debugLevel > 0)
    fprintf(GlobalData->stderrc,"* MarkUnderlyingCIEdgesTrusted on SEdge (" F_CID "," F_CID ")%c nextRaw = " F_CID "\n",
            edge->idA, edge->idB, edge->orient, edge->nextRawEdge);
  
  if(edge->flags.bits.isRaw)
    {
      MarkUnderlyingRawCIEdgeTrusted(sgraph, edge);
    }
  else
    {
      SEdgeT *raw = edge;
      while((raw = GetGraphEdge(sgraph->ScaffoldGraph,
                                raw->nextRawEdge)) != NULL)
        {
          MarkUnderlyingRawCIEdgeTrusted(sgraph, raw);
        }
    }
}


/**************************************************************************/

double FindMaxGapInScaffold(ScaffoldGraphT *graph, CIScaffoldT *scaffold){
  CIScaffoldTIterator Nodes;
  NodeCGW_T *thisNode;
  LengthT *prevEnd, *thisBegin, *thisEnd;
  double maxGap;
  
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &Nodes);
  thisNode = NextCIScaffoldTIterator(&Nodes);
  if(GetNodeOrient(thisNode) == A_B){
    prevEnd = &(thisNode->offsetBEnd);
  }else{
    prevEnd = &(thisNode->offsetAEnd);
  }
  for(maxGap = 0.0; (thisNode = NextCIScaffoldTIterator(&Nodes)) != NULL;
      prevEnd = thisEnd){
    double thisGap;
    if(GetNodeOrient(thisNode) == A_B){
      thisBegin = &(thisNode->offsetAEnd);
      thisEnd = &(thisNode->offsetBEnd);
    }else{
      thisEnd = &(thisNode->offsetAEnd);
      thisBegin = &(thisNode->offsetBEnd);
    }
    thisGap = (thisBegin->mean - prevEnd->mean) +
      (3.0 * sqrt(thisBegin->variance - prevEnd->variance));
    if(thisGap > maxGap){
      maxGap = thisGap;
    }
  }
  return(maxGap);
}

/**************************************************************************/

static int CompareSEdgesContributing(const void *c1, const void *c2){
  SEdgeT *s1 = *(SEdgeT **)c1;
  SEdgeT *s2 = *(SEdgeT **)c2;
  
  return((int)((s2->edgesContributing - (isOverlapEdge(s2) ? 1 : 0)) -
               (s1->edgesContributing - (isOverlapEdge(s1) ? 1 : 0))));
}

#define PREFERRED_GAP_SIZE  (-500)

static int CompareSEdgeGaps(const void *c1, const void *c2){
  SEdgeT *s1 = *(SEdgeT **)c1;
  SEdgeT *s2 = *(SEdgeT **)c2;
  
  // if one has 2+ weight and the other does not, favor the 2+ weight edge
  if(s1->edgesContributing > 1 && s2->edgesContributing < 2)
    {
      return -1;
    }
  else if(s2->edgesContributing > 1 && s1->edgesContributing < 2)
    {
      return 1;
    }
  
  // both or neither have 2+ weight
  // Use proximity to ideal gap size + stddev to order edges
  return((int) ((fabs(PREFERRED_GAP_SIZE - s1->distance.mean) +
                 sqrt(MAX(1.,s1->distance.variance))) -
                (fabs(PREFERRED_GAP_SIZE - s2->distance.mean) +
                 sqrt(MAX(1.,s2->distance.variance)))));
}

/***************************************************************************/

#define MAX_SCAFFOLD_GAP_OVERLAP 1000

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

/****************************************************************************/

static int CompareGapRating(const void *c1, const void *c2){
  SortGapMergeT *s1 = (SortGapMergeT *)c1;
  SortGapMergeT *s2 = (SortGapMergeT *)c2;
  
  if((s1->rating - s2->rating) < 0.0){
    return(-1);
  }
  //else
  return(1);
}
static void SetOrderInfoForGap(OrderCIsTmpT *orderAnchor, int32 index,
                               int firstTime, double oldMergeGapVariance,
                               LengthT anchorGap,
                               double newMergeGapVariance,
                               double gapBeginMean, int pathFromA,
                               LengthT insertSize,
                               double gapEndMean, int pathFromB,
                               int numNodesInGap, int endOfMerge){
  double deltaMean;
  double gapEndVariance;
  double gapEndStdDev;
  double gapBeginStdDev;
  double ratioStdDev;
  
  orderAnchor[index].numInGap = numNodesInGap;
  if(firstTime){
    orderAnchor[index].gapBeginOffset.variance = anchorGap.variance;
  }else{
    orderAnchor[index].gapBeginOffset.variance = MAX(anchorGap.variance,
                                                     oldMergeGapVariance);
  }
  gapBeginStdDev = sqrt(orderAnchor[index].gapBeginOffset.variance);
  orderAnchor[index].gapDelta.variance =
    orderAnchor[index].gapBeginOffset.variance + insertSize.variance;
  if(endOfMerge){
    gapEndVariance = anchorGap.variance;
  }else{
    gapEndVariance = MAX(anchorGap.variance, newMergeGapVariance);
  }
  gapEndStdDev = sqrt(gapEndVariance);
  ratioStdDev = gapBeginStdDev / (gapBeginStdDev + gapEndStdDev);
  orderAnchor[index].gapDelta.variance += gapEndVariance;
  orderAnchor[index].gapDelta.variance -= anchorGap.variance;
  if(pathFromA && pathFromB){
    if((gapBeginMean < - (float)CGW_MISSED_OVERLAP) &&
       (gapEndMean < - (float)CGW_MISSED_OVERLAP)){
      orderAnchor[index].gapBeginOffset.mean = gapBeginMean;
      orderAnchor[index].gapDelta.mean = (gapBeginMean + insertSize.mean +
                                          gapEndMean) - anchorGap.mean;
    }else if(gapBeginMean < - (float)CGW_MISSED_OVERLAP){
      orderAnchor[index].gapBeginOffset.mean = gapBeginMean;
      deltaMean = (gapBeginMean + insertSize.mean) - anchorGap.mean;
      if(deltaMean <= (float)CGW_MISSED_OVERLAP){
        orderAnchor[index].gapDelta.mean = 0.0;
      }else{
        orderAnchor[index].gapDelta.mean = deltaMean -
          (float)CGW_MISSED_OVERLAP;
      }
    }else if(gapEndMean < - (float)CGW_MISSED_OVERLAP){
      deltaMean = (gapEndMean + insertSize.mean) - anchorGap.mean;
      if(deltaMean <= (float)CGW_MISSED_OVERLAP){
        orderAnchor[index].gapDelta.mean = 0.0;
        orderAnchor[index].gapBeginOffset.mean = - deltaMean;
      }else{
        orderAnchor[index].gapDelta.mean = deltaMean -
          (float)CGW_MISSED_OVERLAP;
        orderAnchor[index].gapBeginOffset.mean = - (float)CGW_MISSED_OVERLAP;
      }
    }else{
      deltaMean = (gapBeginMean + insertSize.mean + gapEndMean) -
        anchorGap.mean;
      if(deltaMean <= 0.0){
        orderAnchor[index].gapDelta.mean = 0.0;
        orderAnchor[index].gapBeginOffset.mean = gapBeginMean - 
          (deltaMean * ratioStdDev);
      }else{
        orderAnchor[index].gapBeginOffset.mean = gapBeginMean - 
          (deltaMean * ratioStdDev);
        if(orderAnchor[index].gapBeginOffset.mean <
           - (float)CGW_MISSED_OVERLAP){
          orderAnchor[index].gapBeginOffset.mean = - (float)CGW_MISSED_OVERLAP;
        }
        deltaMean = (orderAnchor[index].gapBeginOffset.mean + insertSize.mean)
          - anchorGap.mean;
        if(deltaMean <= (float)CGW_MISSED_OVERLAP){
          orderAnchor[index].gapDelta.mean = 0.0;
        }else{
          orderAnchor[index].gapBeginOffset.mean -= deltaMean -
            (float)CGW_MISSED_OVERLAP;
          if(orderAnchor[index].gapBeginOffset.mean <
             - (float)CGW_MISSED_OVERLAP){
            orderAnchor[index].gapBeginOffset.mean =
              - (float)CGW_MISSED_OVERLAP;
          }
          deltaMean = (orderAnchor[index].gapBeginOffset.mean +
                       insertSize.mean) - anchorGap.mean;
          if(deltaMean <= (float)CGW_MISSED_OVERLAP){
            orderAnchor[index].gapDelta.mean = 0.0;
          }else{
            orderAnchor[index].gapDelta.mean = deltaMean -
              (float)CGW_MISSED_OVERLAP;
          }
        }
      }
    }
  }else if(pathFromA){
    if(gapBeginMean < - (float)CGW_MISSED_OVERLAP){
      orderAnchor[index].gapBeginOffset.mean = gapBeginMean;
      deltaMean = (gapBeginMean + insertSize.mean) - anchorGap.mean;
      if(deltaMean <= (float)CGW_MISSED_OVERLAP){
        orderAnchor[index].gapDelta.mean = 0.0;
      }else{
        orderAnchor[index].gapDelta.mean = deltaMean -
          (float)CGW_MISSED_OVERLAP;
      }
    }else{
      deltaMean = (gapBeginMean + insertSize.mean) - anchorGap.mean;
      orderAnchor[index].gapBeginOffset.mean = gapBeginMean;
      if(deltaMean <= (float)CGW_MISSED_OVERLAP){
        orderAnchor[index].gapDelta.mean = 0.0;
      }else{
        orderAnchor[index].gapBeginOffset.mean -= deltaMean -
          (float)CGW_MISSED_OVERLAP;
        if(orderAnchor[index].gapBeginOffset.mean <
           - (float)CGW_MISSED_OVERLAP){
          orderAnchor[index].gapBeginOffset.mean =
            - (float)CGW_MISSED_OVERLAP;
        }
        deltaMean = (orderAnchor[index].gapBeginOffset.mean +
                     insertSize.mean) - anchorGap.mean;
        if(deltaMean <= (float)CGW_MISSED_OVERLAP){
          orderAnchor[index].gapDelta.mean = 0.0;
        }else{
          orderAnchor[index].gapDelta.mean = deltaMean -
            (float)CGW_MISSED_OVERLAP;
        }
      }
    }
  }else if(pathFromB){
    if(gapEndMean < - (float)CGW_MISSED_OVERLAP){
      deltaMean = (gapEndMean + insertSize.mean) - anchorGap.mean;
      if(deltaMean <= (float)CGW_MISSED_OVERLAP){
        orderAnchor[index].gapDelta.mean = 0.0;
        orderAnchor[index].gapBeginOffset.mean = - deltaMean;
      }else{
        orderAnchor[index].gapDelta.mean = deltaMean -
          (float)CGW_MISSED_OVERLAP;
        orderAnchor[index].gapBeginOffset.mean = - (float)CGW_MISSED_OVERLAP;
      }
    }else{
      deltaMean = (gapEndMean + insertSize.mean) - anchorGap.mean;
      if(deltaMean <= (float)CGW_MISSED_OVERLAP){
        orderAnchor[index].gapDelta.mean = 0.0;
        orderAnchor[index].gapBeginOffset.mean = - deltaMean;
      }else{
        gapEndMean -= deltaMean - (float)CGW_MISSED_OVERLAP;
        if(gapEndMean < - (float)CGW_MISSED_OVERLAP){
          gapEndMean = - (float)CGW_MISSED_OVERLAP;
        }
        deltaMean = (gapEndMean + insertSize.mean) - anchorGap.mean;
        if(deltaMean <= (float)CGW_MISSED_OVERLAP){
          orderAnchor[index].gapDelta.mean = 0.0;
          orderAnchor[index].gapBeginOffset.mean = - deltaMean;
        }else{
          orderAnchor[index].gapDelta.mean = deltaMean -
            (float)CGW_MISSED_OVERLAP;
          orderAnchor[index].gapBeginOffset.mean = - (float)CGW_MISSED_OVERLAP;
        }
      }
    }
  }else{
    assert(FALSE);//!pathFromA && !pathFromB
  }
  
  return;
}

int FindScaffoldMerge(ScaffoldGraphT *graph, CIScaffoldT *scaffoldAnchor,
                      CIScaffoldT *scaffoldMerge, ChunkOrient orientAnchor,
                      ChunkOrient orientMerge, NodeCGW_T *endNodeMerge,
                      ChunkOrient orientEndNodeMerge, double intervalBegin,
                      double intervalEnd, double intervalGapBegin,
                      double intervalGapEnd, int firstTime,
                      double maxMergeGapSize, double mergeGapVariance,
                      LengthT lastAnchorOffset, OrderCIsTmpT *orderAnchor,
                      int verbose){
  CIScaffoldTIterator Nodes;
  NodeCGW_T *thisNode, *prevNode;
  double minEndNodeMergeSize;
  int foundBeginNode;
  SortGapMergeT *sortGaps, *sortGapsPtr, *sortGapsEnd;
  int numGapsRated;
  
  if(verbose){
    fprintf(GlobalData->stderrc, "FindScaffoldMerge AnchorID " F_CID " MergeID " F_CID "\norientA %c orientM %c nodeM " F_CID " orientNode %c\nBeg %f End %f GapBeg %f GapEnd %f first %d\n maxSize %f Var %f OffsetMean %f OffsetVariance %f\n",
            scaffoldAnchor->id, scaffoldMerge->id, (char)orientAnchor,
            (char)orientMerge, endNodeMerge->id, (char)orientEndNodeMerge,
            intervalBegin, intervalEnd, intervalGapBegin, intervalGapEnd,
            firstTime, maxMergeGapSize, mergeGapVariance,
            lastAnchorOffset.mean, lastAnchorOffset.variance);
  }
  sortGaps = (SortGapMergeT *)safe_malloc(
                                          (scaffoldAnchor->info.Scaffold.numElements - 1) *
                                          sizeof(*sortGaps));
  AssertPtr(sortGaps);
  InitCIScaffoldTIterator(graph, scaffoldAnchor, (orientAnchor == B_A),
                          FALSE, &Nodes);
  minEndNodeMergeSize = endNodeMerge->bpLength.mean -
    (3.0 * sqrt(endNodeMerge->bpLength.variance));
  if(verbose){
    fprintf(GlobalData->stderrc, "minEndNodeMergeSize %f\n", minEndNodeMergeSize);
  }
  for(prevNode = (NodeCGW_T *)NULL, sortGapsPtr = sortGaps,
        foundBeginNode = FALSE;
      (thisNode = NextCIScaffoldTIterator(&Nodes)) != NULL;
      prevNode = thisNode){
    if(!foundBeginNode){
      if(((orientAnchor == B_A) && // traversing from A end
          ((thisNode->offsetAEnd.mean > intervalBegin) &&
           (thisNode->offsetBEnd.mean > intervalBegin))) ||
         ((orientAnchor == A_B) && // traversing from B end
          ((thisNode->offsetAEnd.mean < intervalBegin) &&
           (thisNode->offsetBEnd.mean < intervalBegin)))){
        foundBeginNode = TRUE;
        if(verbose){
          fprintf(GlobalData->stderrc, "Found beginNode " F_CID "\n", thisNode->id);
        }
        if(firstTime && (prevNode == (NodeCGW_T *)NULL)){
          continue;
        }
      }else{
        continue;
      }
    }
    if(prevNode != (NodeCGW_T *)NULL){
      LengthT gap;
      double rating;
      double maxGapSize;
      LengthT anchorInsertSize;
      double minAnchorInsertSize = 0.;
      
      if(orientAnchor == B_A){// traversing from A end
        if(GetNodeOrient(thisNode) == A_B){
          anchorInsertSize.mean = lastAnchorOffset.mean -
            thisNode->offsetAEnd.mean;
          anchorInsertSize.variance = lastAnchorOffset.variance -
            thisNode->offsetAEnd.variance;
          gap.mean = thisNode->offsetAEnd.mean;
          gap.variance = thisNode->offsetAEnd.variance;
          if(thisNode->offsetAEnd.mean > intervalGapEnd){
            rating = 0.0;
          }else{
            rating = intervalGapEnd - thisNode->offsetAEnd.mean;
          }
        }else{
          anchorInsertSize.mean = lastAnchorOffset.mean -
            thisNode->offsetBEnd.mean;
          anchorInsertSize.variance = lastAnchorOffset.variance -
            thisNode->offsetBEnd.variance;
          gap.mean = thisNode->offsetBEnd.mean;
          gap.variance = thisNode->offsetBEnd.variance;
          if(thisNode->offsetBEnd.mean > intervalGapEnd){
            rating = 0.0;
          }else{
            rating = intervalGapEnd - thisNode->offsetBEnd.mean;
          }
        }
        if(GetNodeOrient(prevNode) == A_B){
          gap.mean -= prevNode->offsetBEnd.mean;
          gap.variance -= prevNode->offsetBEnd.variance;
          if(prevNode->offsetBEnd.mean > intervalGapBegin){
            rating += prevNode->offsetBEnd.mean - intervalGapBegin;
          }
        }else{
          gap.mean -= prevNode->offsetAEnd.mean;
          gap.variance -= prevNode->offsetAEnd.variance;
          if(prevNode->offsetAEnd.mean > intervalGapBegin){
            rating += prevNode->offsetAEnd.mean - intervalGapBegin;
          }
        }
      }else{// traversing from B end
        if(GetNodeOrient(prevNode) == A_B){
          gap.mean = prevNode->offsetAEnd.mean;
          gap.variance = prevNode->offsetAEnd.variance;
          if(prevNode->offsetAEnd.mean > intervalGapBegin){
            rating = 0.0;
          }else{
            rating = intervalGapBegin - prevNode->offsetAEnd.mean;
          }
        }else{
          gap.mean = prevNode->offsetBEnd.mean;
          gap.variance = prevNode->offsetBEnd.variance;
          if(prevNode->offsetBEnd.mean > intervalGapBegin){
            rating = 0.0;
          }else{
            rating = intervalGapBegin - prevNode->offsetBEnd.mean;
          }
        }
        if(GetNodeOrient(thisNode) == A_B){
          anchorInsertSize.mean = thisNode->offsetBEnd.mean -
            lastAnchorOffset.mean;
          anchorInsertSize.variance = thisNode->offsetBEnd.variance -
            lastAnchorOffset.variance;
          gap.mean -= thisNode->offsetBEnd.mean;
          gap.variance -= thisNode->offsetBEnd.variance;
          if(thisNode->offsetBEnd.mean > intervalGapEnd){
            rating += thisNode->offsetBEnd.mean - intervalGapEnd;
          }
        }else{
          anchorInsertSize.mean = thisNode->offsetAEnd.mean -
            lastAnchorOffset.mean;
          anchorInsertSize.variance = thisNode->offsetAEnd.variance -
            lastAnchorOffset.variance;
          gap.mean -= thisNode->offsetAEnd.mean;
          gap.variance -= thisNode->offsetAEnd.variance;
          if(thisNode->offsetAEnd.mean > intervalGapEnd){
            rating += thisNode->offsetAEnd.mean - intervalGapEnd;
          }
        }
      }
      maxGapSize = gap.mean + (3.0 * sqrt(gap.variance)) +
        MAX_SCAFFOLD_GAP_OVERLAP;
      if(!firstTime){
        minAnchorInsertSize = anchorInsertSize.mean -
          (3.0 * sqrt(anchorInsertSize.variance));
      }
      if((minEndNodeMergeSize <= maxGapSize) &&
         (firstTime || (minAnchorInsertSize <= maxMergeGapSize))){
        CIScaffoldTIterator mergeNodes;
        NodeCGW_T *mergeNode;
        LengthT beginOffset, endOffset;
        double complement;
        double minCombinedSize;
        int numToMerge;
        sortGapsPtr->nodeA = thisNode;
        sortGapsPtr->nodeB = prevNode;
        sortGapsPtr->rating = rating;
        sortGapsPtr->maxSpace = maxGapSize - minEndNodeMergeSize;
        sortGapsPtr->maxGapSize = maxGapSize;
        sortGapsPtr->gapSize = gap;
        if(orientMerge == A_B){
          complement = 1.0;
        }else{
          complement = -1.0;
        }
        if(orientEndNodeMerge == A_B){
          beginOffset = endNodeMerge->offsetAEnd;
        }else{
          beginOffset = endNodeMerge->offsetBEnd;
        }
        InitCIScaffoldTIteratorFromCI(graph, scaffoldMerge,
                                      endNodeMerge->id,
                                      (orientMerge == A_B),
                                      FALSE, &mergeNodes);
        for(mergeNode = NextCIScaffoldTIterator(&mergeNodes),
              numToMerge = 1;
            (mergeNode = NextCIScaffoldTIterator(&mergeNodes)) !=
              NULL; numToMerge++){
          ChunkOrient orientMergeNode;
          orientMergeNode = GetNodeOrient(mergeNode);
          if(orientMerge == B_A){
            orientMergeNode = FlipNodeOrient(orientMergeNode);
          }
          if(orientMergeNode == A_B){
            endOffset = mergeNode->offsetBEnd;
          }else{
            endOffset = mergeNode->offsetAEnd;
          }
          minCombinedSize = (complement * (endOffset.mean -
                                           beginOffset.mean)) -
            (3.0 * sqrt(complement * (endOffset.variance -
                                      beginOffset.variance)));
          if(minCombinedSize > maxGapSize){
            break;
          }
        }
        sortGapsPtr->numToMerge = numToMerge;
        if(verbose){
          fprintf(GlobalData->stderrc, "Gap(" F_CID "," F_CID ") rating %f maxS %f maxG %f\ngapMean %f gapStdDev %f numToMerge %d\n",
                  (sortGapsPtr->nodeA)->id, (sortGapsPtr->nodeB)->id,
                  sortGapsPtr->rating, sortGapsPtr->maxSpace,
                  sortGapsPtr->maxGapSize, sortGapsPtr->gapSize.mean,
                  sqrt(sortGapsPtr->gapSize.variance),
                  sortGapsPtr->numToMerge);
        }
        sortGapsPtr++;
      }
    }else{
      LengthT gap;
      double rating;
      double maxGapSize;
      LengthT anchorInsertSize;
      double minAnchorInsertSize = 0.;
      
      if(orientAnchor == B_A){// traversing from A end
        if(GetNodeOrient(thisNode) == A_B){
          anchorInsertSize.mean = lastAnchorOffset.mean -
            thisNode->offsetAEnd.mean;
          anchorInsertSize.variance = lastAnchorOffset.variance -
            thisNode->offsetAEnd.variance;
          gap.mean = thisNode->offsetAEnd.mean;
          if(thisNode->offsetAEnd.mean > intervalGapEnd){
            rating = 0.0;
          }else{
            rating = intervalGapEnd - thisNode->offsetAEnd.mean;
          }
        }else{
          anchorInsertSize.mean = lastAnchorOffset.mean -
            thisNode->offsetBEnd.mean;
          anchorInsertSize.variance = lastAnchorOffset.variance -
            thisNode->offsetBEnd.variance;
          gap.mean = thisNode->offsetBEnd.mean;
          if(thisNode->offsetBEnd.mean > intervalGapEnd){
            rating = 0.0;
          }else{
            rating = intervalGapEnd - thisNode->offsetBEnd.mean;
          }
        }
        gap.mean -= intervalGapBegin;
        gap.variance = (intervalGapBegin - intervalBegin) / 3.0;
        gap.variance *= gap.variance;
      }else{// traversing from B end
        gap.mean = intervalGapBegin;
        gap.variance = (intervalBegin - intervalGapBegin) / 3.0;
        gap.variance *= gap.variance;
        if(GetNodeOrient(thisNode) == A_B){
          anchorInsertSize.mean = thisNode->offsetBEnd.mean -
            lastAnchorOffset.mean;
          anchorInsertSize.variance = thisNode->offsetBEnd.variance -
            lastAnchorOffset.variance;
          gap.mean -= thisNode->offsetBEnd.mean;
          if(thisNode->offsetBEnd.mean < intervalGapEnd){
            rating = 0.0;
          }else{
            rating = thisNode->offsetBEnd.mean - intervalGapEnd;
          }
        }else{
          anchorInsertSize.mean = thisNode->offsetAEnd.mean -
            lastAnchorOffset.mean;
          anchorInsertSize.variance = thisNode->offsetAEnd.variance -
            lastAnchorOffset.variance;
          gap.mean -= thisNode->offsetAEnd.mean;
          if(thisNode->offsetAEnd.mean < intervalGapEnd){
            rating = 0.0;
          }else{
            rating = thisNode->offsetAEnd.mean - intervalGapEnd;
          }
        }
      }
      maxGapSize = gap.mean + (3.0 * sqrt(gap.variance)) +
        MAX_SCAFFOLD_GAP_OVERLAP;
      if(!firstTime){
        minAnchorInsertSize = anchorInsertSize.mean -
          (3.0 * sqrt(anchorInsertSize.variance));
      }
      if((minEndNodeMergeSize <= maxGapSize) &&
         (firstTime || (minAnchorInsertSize <= maxMergeGapSize))){
        sortGapsPtr->nodeA = thisNode;
        sortGapsPtr->nodeB = prevNode;
        sortGapsPtr->rating = rating;
        sortGapsPtr->maxSpace = maxGapSize - minEndNodeMergeSize;
        sortGapsPtr->maxGapSize = maxGapSize;
        sortGapsPtr->gapSize = gap;
        sortGapsPtr->numToMerge = 1;
        if(verbose){
          fprintf(GlobalData->stderrc, "Gap(" F_CID "," F_CID ") rating %f maxS %f maxG %f\ngapMean %f gapStdDev %f numToMerge %d\n",
                  (sortGapsPtr->nodeA)->id, NULLINDEX,
                  sortGapsPtr->rating, sortGapsPtr->maxSpace,
                  sortGapsPtr->maxGapSize, sortGapsPtr->gapSize.mean,
                  sqrt(sortGapsPtr->gapSize.variance),
                  sortGapsPtr->numToMerge);
        }
        sortGapsPtr++;
      }
    }
    if(((orientAnchor == B_A) && // traversing from A end
        ((thisNode->offsetAEnd.mean > intervalEnd) ||
         (thisNode->offsetBEnd.mean > intervalEnd))) ||
       ((orientAnchor == A_B) && // traversing from B end
        ((thisNode->offsetAEnd.mean < intervalEnd) ||
         (thisNode->offsetBEnd.mean < intervalEnd)))){
      if(verbose){
        fprintf(GlobalData->stderrc, "Found endNode " F_CID "\n", thisNode->id);
      }
      break;
    }
  }
  sortGapsEnd = sortGapsPtr;
  numGapsRated = sortGapsEnd - sortGaps;
  if(numGapsRated == 0){
    safe_free(sortGaps);
    if(verbose){
      fprintf(GlobalData->stderrc, "endNodeMerge did not fit in any of scaffoldAnchor's gaps.\n");
    }
    return(FALSE);
  }
  qsort((void *)sortGaps, numGapsRated, sizeof(*sortGaps),
        CompareGapRating);
  for(sortGapsPtr = sortGaps; sortGapsPtr < sortGapsEnd;
      sortGapsPtr++){
    Target_Info_t *mergeTargetsA, *mergeTargetPtrA, *mergeTargetEndA;
    int firstTargetA, numTargetsBestPathA, firstTargetBestPathA;
    Target_Info_t *mergeTargetsB, *mergeTargetPtrB, *mergeTargetEndB;
    int firstTargetB, numTargetsBestPathB, firstTargetBestPathB;
    int anchorNodeEnd = NO_END;
    ChunkOrient orientEndNodeAnchor, orientNewEndNodeMerge;
    CDS_COORD_t *simCoordsWhere, *simCoordsPtr;
    CIScaffoldTIterator mergeNodes;
    NodeCGW_T *thisNode, *newEndNodeMerge;
    int numToMerge;
    LengthT beginOffset;
    LengthT endOffset;
    LengthT targetPosition;
    double complement;
    int FoundPathToEndNodeMerge = FALSE;
    double endNodeMergeGapOffsetA = 0., endNodeMergeGapOffsetB = 0.;
    int pathFromA = FALSE, pathFromB = FALSE;
    int numNodesInGap = 0;
    
    if(verbose){
      fprintf(GlobalData->stderrc, "Working on Gap(" F_CID "," F_CID ") rating %f maxS %f maxG %f\ngapMean %f gapStdDev %f numToMerge %d\n",
              (sortGapsPtr->nodeA)->id, (sortGapsPtr->nodeB != NULL) ?
              (sortGapsPtr->nodeB)->id : NULLINDEX,
              sortGapsPtr->rating, sortGapsPtr->maxSpace,
              sortGapsPtr->maxGapSize, sortGapsPtr->gapSize.mean,
              sqrt(sortGapsPtr->gapSize.variance),
              sortGapsPtr->numToMerge);
    }
    mergeTargetsA = (Target_Info_t *)safe_malloc(sortGapsPtr->numToMerge *
                                                 sizeof(*mergeTargetsA));
    AssertPtr(mergeTargetsA);
    mergeTargetEndA = mergeTargetsA + sortGapsPtr->numToMerge;
    mergeTargetsB = (Target_Info_t *)safe_malloc(sortGapsPtr->numToMerge *
                                                 sizeof(*mergeTargetsB));
    AssertPtr(mergeTargetsB);
    mergeTargetEndB = mergeTargetsB + sortGapsPtr->numToMerge;
    simCoordsWhere = (CDS_COORD_t *)safe_malloc(sortGapsPtr->numToMerge *
                                                sizeof(*simCoordsWhere));
    AssertPtr(simCoordsWhere);
    /* Try to find a path from the "A" side of the anchor scaffold
       gap to the merge scaffold end node. */
    orientEndNodeAnchor = GetNodeOrient(sortGapsPtr->nodeA);
    if(orientAnchor == B_A){
      orientEndNodeAnchor = FlipNodeOrient(orientEndNodeAnchor);
    }
    if(orientMerge == A_B){
      complement = 1.0;
    }else{
      complement = -1.0;
    }
    if(orientEndNodeMerge == A_B){
      beginOffset = endNodeMerge->offsetAEnd;
    }else{
      beginOffset = endNodeMerge->offsetBEnd;
    }
    InitCIScaffoldTIteratorFromCI(graph, scaffoldMerge,
                                  endNodeMerge->id,
                                  (orientMerge == A_B),
                                  FALSE, &mergeNodes);
    for(numToMerge = sortGapsPtr->numToMerge,
          mergeTargetPtrA = mergeTargetsA,
          simCoordsPtr = simCoordsWhere;
        numToMerge > 0; numToMerge--, mergeTargetPtrA++,
          simCoordsPtr++){
      ChunkOrient orientThisNode;
      LengthT bpLength;
      if((thisNode = NextCIScaffoldTIterator(&mergeNodes)) == NULL)
        assert(0);
      assert(mergeTargetPtrA < mergeTargetEndA);
      orientThisNode = GetNodeOrient(thisNode);
      if(orientMerge == B_A){
        orientThisNode = FlipNodeOrient(orientThisNode);
      }
      if(orientThisNode == A_B){
        endOffset = thisNode->offsetBEnd;
      }else{
        endOffset = thisNode->offsetAEnd;
      }
      bpLength.mean = complement * (endOffset.mean - beginOffset.mean);
      bpLength.variance = complement * (endOffset.variance -
                                        beginOffset.variance);
      mergeTargetPtrA->found = FALSE;
      mergeTargetPtrA->id = thisNode->id;
      mergeTargetPtrA->lo = (double)(- MAX_SCAFFOLD_GAP_OVERLAP) +
        bpLength.mean - (3.0 * sqrt(bpLength.variance));
      mergeTargetPtrA->hi = sortGapsPtr->maxGapSize;
      if(orientEndNodeAnchor == A_B){
        if(orientThisNode == A_B){
          mergeTargetPtrA->orient = AB_AB;
          anchorNodeEnd = B_END;
          *simCoordsPtr = abs(thisNode->bEndCoord -
                              (sortGapsPtr->nodeA)->bEndCoord);
        }else{
          mergeTargetPtrA->orient = AB_BA;
          anchorNodeEnd = B_END;
          *simCoordsPtr = abs(thisNode->aEndCoord -
                              (sortGapsPtr->nodeA)->bEndCoord);
        }
      }else{
        if(orientThisNode == A_B){
          mergeTargetPtrA->orient = BA_AB;
          anchorNodeEnd = A_END;
          *simCoordsPtr = abs(thisNode->bEndCoord -
                              (sortGapsPtr->nodeA)->aEndCoord);
        }else{
          mergeTargetPtrA->orient = BA_BA;
          anchorNodeEnd = A_END;
          *simCoordsPtr = abs(thisNode->aEndCoord -
                              (sortGapsPtr->nodeA)->aEndCoord);
        }
      }
      if(verbose){
        fprintf(GlobalData->stderrc, "mergeTargetA " F_CID " hi %f lo %f orient %c\n",
                mergeTargetPtrA->id, mergeTargetPtrA->hi, mergeTargetPtrA->lo,
                (char)mergeTargetPtrA->orient);
      }
    }
    if(verbose){
      fprintf(GlobalData->stderrc, "Find_Olap_Path nodeA " F_CID " anchorNodeEnd %d numTargets %d gapBound %f\n",
              (sortGapsPtr->nodeA)->id, anchorNodeEnd,
              sortGapsPtr->numToMerge, sortGapsPtr->maxSpace);
    }
    if(Find_Olap_Path(sortGapsPtr->nodeA, anchorNodeEnd,
                      (ChunkInstanceT *)NULL, sortGapsPtr->numToMerge,
                      mergeTargetsA,
                      sortGapsPtr->maxSpace, &firstTargetA,
                      &numTargetsBestPathA, &firstTargetBestPathA,
                      &targetPosition, USE_TANDEM_OLAPS)){
      int targetIndex;
      int expectedOrder;
      
      pathFromA = TRUE;
      for(expectedOrder = targetIndex = firstTargetBestPathA;
          targetIndex != NULLINDEX;
          targetIndex = mergeTargetsA[targetIndex].next,
            expectedOrder++){
        fprintf(GlobalData->stderrc, "ScafAnodeA:" F_CID "," F_CID ",%c ScafM:" F_CID "," F_CID ",%c %f(%d) %d\n",
                scaffoldAnchor->id, (sortGapsPtr->nodeA)->id,
                orientEndNodeAnchor, scaffoldMerge->id,
                mergeTargetsA[targetIndex].id,
                mergeTargetsA[targetIndex].orient,
                mergeTargetsA[targetIndex].where,
                simCoordsWhere[targetIndex],
                mergeTargetsA[targetIndex].found);
        mergeTargetsA[targetIndex].found = TRUE;
        if(expectedOrder != targetIndex){
          fprintf(GlobalData->stderrc, "WARNING: mergeTargets in path found are not in same order as scaffold!\n");
        }
      }
      if(!mergeTargetsA[0].found){
        if(verbose){
          fprintf(GlobalData->stderrc, "Found a path from node A that did not include the end node.\n");
        }
        // Found a path that did not include the end node - probably
        // need the next gap over.
        continue;
      }
      FoundPathToEndNodeMerge = TRUE;
      for(targetIndex = sortGapsPtr->numToMerge - 1; targetIndex > 0;
          targetIndex--){
        if(mergeTargetsA[targetIndex].found){
          break;
        }
      }
      numNodesInGap = targetIndex + 1;
      endNodeMergeGapOffsetA = mergeTargetsA[0].where -
        endNodeMerge->bpLength.mean;
      if(numNodesInGap < sortGapsPtr->numToMerge){
        double minAnchorGapSize;
        double gapVariance;
	
        gapVariance = sortGapsPtr->gapSize.variance +
          (sortGapsPtr->nodeB)->bpLength.variance;
        minAnchorGapSize = sortGapsPtr->gapSize.mean +
          (sortGapsPtr->nodeB)->bpLength.mean - endNodeMergeGapOffsetA;
        minAnchorGapSize -= (3.0 * sqrt(gapVariance));
        InitCIScaffoldTIteratorFromCI(graph, scaffoldMerge,
                                      mergeTargetsA[numNodesInGap].id,
                                      (orientMerge == A_B), FALSE,
                                      &mergeNodes);
        for(numToMerge = sortGapsPtr->numToMerge - numNodesInGap;
            numToMerge > 0; numToMerge--){
          ChunkOrient orientThisNode;
          LengthT bpLength;
          
          if((thisNode = NextCIScaffoldTIterator(&mergeNodes)) == NULL)
            assert(0);
          
          orientThisNode = GetNodeOrient(thisNode);
          if(orientMerge == B_A){
            orientThisNode = FlipNodeOrient(orientThisNode);
          }
          if(orientThisNode == A_B){
            endOffset = thisNode->offsetBEnd;
          }else{
            endOffset = thisNode->offsetAEnd;
          }
          bpLength.mean = complement * (endOffset.mean - beginOffset.mean);
          bpLength.variance = (complement * (endOffset.variance -
                                             beginOffset.variance)) +
            gapVariance;
          if((bpLength.mean + (3.0 * sqrt(bpLength.variance))) >
             minAnchorGapSize){
            break;
          }
          numNodesInGap++;
        }
      }
    }
    if(sortGapsPtr->nodeB != (NodeCGW_T *)NULL){
      /* Try to find a path from the "B" side of the anchor scaffold
         gap to the merge scaffold end node. */
      orientEndNodeAnchor = GetNodeOrient(sortGapsPtr->nodeB);
      if(orientAnchor == A_B){
        orientEndNodeAnchor = FlipNodeOrient(orientEndNodeAnchor);
      }
      InitCIScaffoldTIteratorFromCI(graph, scaffoldMerge,
                                    endNodeMerge->id,
                                    (orientMerge == A_B),
                                    FALSE, &mergeNodes);
      for(numToMerge = sortGapsPtr->numToMerge,
            mergeTargetPtrB = mergeTargetsB,
            simCoordsPtr = simCoordsWhere;
          numToMerge > 0; numToMerge--, mergeTargetPtrB++,
            simCoordsPtr++){
        ChunkOrient orientThisNode;
        LengthT bpLength;
        if((thisNode = NextCIScaffoldTIterator(&mergeNodes)) == NULL)
          assert(0);
        assert(mergeTargetPtrB < mergeTargetEndB);
        orientThisNode = GetNodeOrient(thisNode);
        if(orientMerge == A_B){
          orientThisNode = FlipNodeOrient(orientThisNode);
        }
        if(orientThisNode == A_B){
          endOffset = thisNode->offsetBEnd;
        }else{
          endOffset = thisNode->offsetAEnd;
        }
        bpLength.mean = complement * (endOffset.mean - beginOffset.mean);
        bpLength.variance = complement * (endOffset.variance -
                                          beginOffset.variance);
        mergeTargetPtrB->found = FALSE;
        mergeTargetPtrB->id = thisNode->id;
        mergeTargetPtrB->lo = (double)(- MAX_SCAFFOLD_GAP_OVERLAP) +
          thisNode->bpLength.mean -
          (3.0 * sqrt(thisNode->bpLength.variance));
        mergeTargetPtrB->hi = sortGapsPtr->maxGapSize - (bpLength.mean -
                                                         (3.0 * sqrt(bpLength.variance)));
        if(orientEndNodeAnchor == A_B){
          if(orientThisNode == A_B){
            mergeTargetPtrB->orient = AB_AB;
            anchorNodeEnd = B_END;
            *simCoordsPtr = abs(thisNode->bEndCoord -
                                (sortGapsPtr->nodeB)->bEndCoord);
          }else{
            mergeTargetPtrB->orient = AB_BA;
            anchorNodeEnd = B_END;
            *simCoordsPtr = abs(thisNode->aEndCoord -
                                (sortGapsPtr->nodeB)->bEndCoord);
          }
        }else{
          if(orientThisNode == A_B){
            mergeTargetPtrB->orient = BA_AB;
            anchorNodeEnd = A_END;
            *simCoordsPtr = abs(thisNode->bEndCoord -
                                (sortGapsPtr->nodeB)->aEndCoord);
          }else{
            mergeTargetPtrB->orient = BA_BA;
            anchorNodeEnd = A_END;
            *simCoordsPtr = abs(thisNode->aEndCoord -
                                (sortGapsPtr->nodeB)->aEndCoord);
          }
        }
        if(verbose){
          fprintf(GlobalData->stderrc, "mergeTargetB " F_CID " hi %f lo %f orient %c\n",
                  mergeTargetPtrB->id, mergeTargetPtrB->hi,
                  mergeTargetPtrB->lo, (char)mergeTargetPtrB->orient);
        }
      }
      if(verbose){
        fprintf(GlobalData->stderrc, "Find_Olap_Path nodeA " F_CID " anchorNodeEnd %d numTargets %d gapBound %f\n",
                (sortGapsPtr->nodeB)->id, anchorNodeEnd,
                sortGapsPtr->numToMerge, sortGapsPtr->maxSpace);
      }
      if(Find_Olap_Path(sortGapsPtr->nodeB, anchorNodeEnd,
                        (ChunkInstanceT *)NULL, sortGapsPtr->numToMerge,
                        mergeTargetsB,
                        sortGapsPtr->maxSpace, &firstTargetB,
                        &numTargetsBestPathB, &firstTargetBestPathB,
                        &targetPosition, USE_TANDEM_OLAPS)){
        int targetIndex;
        int expectedOrder;
        NodeCGW_T *lastNode;
        
        pathFromB = TRUE;
        for(expectedOrder = targetIndex = firstTargetBestPathB;
            targetIndex != NULLINDEX;
            targetIndex = mergeTargetsB[targetIndex].next,
              expectedOrder--){
          fprintf(GlobalData->stderrc, "ScafAnodeB:" F_CID "," F_CID ",%c ScafM:" F_CID "," F_CID ",%c %f(%d) %d\n",
                  scaffoldAnchor->id, (sortGapsPtr->nodeB)->id,
                  orientEndNodeAnchor, scaffoldMerge->id,
                  mergeTargetsB[targetIndex].id,
                  mergeTargetsB[targetIndex].orient,
                  mergeTargetsB[targetIndex].where,
                  simCoordsWhere[targetIndex],
                  mergeTargetsB[targetIndex].found);
          mergeTargetsB[targetIndex].found = TRUE;
          if(expectedOrder != targetIndex){
            fprintf(GlobalData->stderrc, "WARNING: mergeTargets in path found are not in same order as scaffold!\n");
          }
        }
        for(targetIndex = sortGapsPtr->numToMerge - 1; targetIndex > 0;
            targetIndex--){
          if(mergeTargetsB[targetIndex].found){
            break;
          }
        }
        numNodesInGap = targetIndex + 1;
        lastNode = GetGraphNode(graph->RezGraph,
                                mergeTargetsB[targetIndex].id);
        endNodeMergeGapOffsetB = mergeTargetsB[targetIndex].where -
          lastNode->bpLength.mean;
        if(!FoundPathToEndNodeMerge && !mergeTargetsB[0].found){
          double minAnchorGapSize;
          double gapVariance;
          ChunkOrient orientLastNode;
          LengthT bpLength;
          if(verbose){
            fprintf(GlobalData->stderrc, "Found a path from node B that did not include the end node.\n");
          }
          // Found a path that did not include the end node - probably
          // need the next gap over but check if endNodeMerge must go here.
          gapVariance = sortGapsPtr->gapSize.variance +
            (sortGapsPtr->nodeA)->bpLength.variance;
          minAnchorGapSize = sortGapsPtr->gapSize.mean +
            (sortGapsPtr->nodeA)->bpLength.mean - endNodeMergeGapOffsetB;
          minAnchorGapSize -= (3.0 * sqrt(gapVariance));
          orientLastNode = GetNodeOrient(lastNode);
          if(orientMerge == B_A){
            orientLastNode = FlipNodeOrient(orientLastNode);
          }
          if(orientLastNode == A_B){
            endOffset = lastNode->offsetBEnd;
          }else{
            endOffset = lastNode->offsetAEnd;
          }
          bpLength.mean = complement * (endOffset.mean - beginOffset.mean);
          bpLength.variance = (complement * (endOffset.variance -
                                             beginOffset.variance)) +
            gapVariance;
          if((bpLength.mean + (3.0 * sqrt(bpLength.variance))) >
             minAnchorGapSize){
            continue;//Does not have to fit here.
          }
        }
        FoundPathToEndNodeMerge = TRUE;
      }
    }
    
    safe_free(mergeTargetsA);
    safe_free(mergeTargetsB);
    safe_free(simCoordsWhere);
    if(!FoundPathToEndNodeMerge){
      if(verbose){
        fprintf(GlobalData->stderrc, "Did not find a path to the end node.\n");
      }
      continue;
    }
    if(sortGapsPtr->nodeB == (NodeCGW_T *)NULL){
      const LengthT nullLength = {0.0, 0.0};
      SetOrderInfoForGap(orderAnchor, (sortGapsPtr->nodeA)->indexInScaffold,
                         firstTime, mergeGapVariance, sortGapsPtr->gapSize,
                         0.0, endNodeMergeGapOffsetA, pathFromA, nullLength,
                         endNodeMergeGapOffsetB, pathFromB,
                         numNodesInGap, TRUE);
      if(verbose){
        fprintf(GlobalData->stderrc, "Found a path to the end node at end of anchor scaffold.\n");
      }
      return(TRUE);
    }
    {
      double gapVariance;
      double maxMergeGapSize;
      double newMergeGapVariance;
      double endNodeMergeGapOffset;
      LengthT mergeInsertSize;
      LengthT anchorOffset, lastAnchorOffset;
      LengthT mergeGapSize;
      NodeCGW_T *lastNodeMerge = NULL;
      LengthT lastNodeMergeOffset;
      ChunkOrient orientLastNodeMerge;
      
      gapVariance = sortGapsPtr->gapSize.variance +
        (sortGapsPtr->nodeB)->bpLength.variance;
      InitCIScaffoldTIteratorFromCI(graph, scaffoldMerge, endNodeMerge->id,
                                    (orientMerge == A_B), FALSE, &mergeNodes);
      for(numToMerge = numNodesInGap;
          numToMerge > 0; numToMerge--){
        lastNodeMerge = NextCIScaffoldTIterator(&mergeNodes);
      }
      newEndNodeMerge = NextCIScaffoldTIterator(&mergeNodes);
      orientLastNodeMerge = GetNodeOrient(lastNodeMerge);
      if(orientMerge == B_A){
        orientLastNodeMerge = FlipNodeOrient(orientLastNodeMerge);
      }
      if(orientLastNodeMerge == A_B){
        lastNodeMergeOffset = lastNodeMerge->offsetBEnd;
      }else{
        lastNodeMergeOffset = lastNodeMerge->offsetAEnd;
      }
      mergeInsertSize.mean = complement * (lastNodeMergeOffset.mean -
                                           beginOffset.mean);
      mergeInsertSize.variance = complement * (lastNodeMergeOffset.variance -
                                               beginOffset.variance);
      if(newEndNodeMerge == (NodeCGW_T *)NULL){
        SetOrderInfoForGap(orderAnchor,
                           (sortGapsPtr->nodeA)->indexInScaffold, firstTime,
                           mergeGapVariance, sortGapsPtr->gapSize, 0.0,
                           endNodeMergeGapOffsetA, pathFromA,
                           mergeInsertSize, endNodeMergeGapOffsetB,
                           pathFromB, numNodesInGap, TRUE);
        if(verbose){
          fprintf(GlobalData->stderrc, "Found a path to the end node at end of merge scaffold.\n");
        }
        return(TRUE);
      }
      orientNewEndNodeMerge = GetNodeOrient(newEndNodeMerge);
      if(orientMerge == B_A){
        orientNewEndNodeMerge = FlipNodeOrient(orientNewEndNodeMerge);
      }
      if(orientNewEndNodeMerge == A_B){
        endOffset = newEndNodeMerge->offsetAEnd;
      }else{
        endOffset = newEndNodeMerge->offsetBEnd;
      }
      mergeGapSize.mean = complement * (endOffset.mean - beginOffset.mean);
      mergeGapSize.variance = complement * (endOffset.variance -
                                            beginOffset.variance);
      newMergeGapVariance = complement * (endOffset.variance -
                                          lastNodeMergeOffset.variance);
      maxMergeGapSize = (complement * (endOffset.mean -
                                       lastNodeMergeOffset.mean)) +
        (3.0 * sqrt(newMergeGapVariance)) + MAX_SCAFFOLD_GAP_OVERLAP;
      orientEndNodeAnchor = GetNodeOrient(sortGapsPtr->nodeB);
      if(orientAnchor == B_A){
        orientEndNodeAnchor = FlipNodeOrient(orientEndNodeAnchor);
      }
      if(orientEndNodeAnchor == A_B){
        lastAnchorOffset = (sortGapsPtr->nodeB)->offsetAEnd;
      }else{
        lastAnchorOffset = (sortGapsPtr->nodeB)->offsetBEnd;
      }
      orientEndNodeAnchor = GetNodeOrient(sortGapsPtr->nodeA);
      if(orientAnchor == B_A){
        orientEndNodeAnchor = FlipNodeOrient(orientEndNodeAnchor);
      }
      if(orientEndNodeAnchor == A_B){
        anchorOffset = (sortGapsPtr->nodeA)->offsetBEnd;
      }else{
        anchorOffset = (sortGapsPtr->nodeA)->offsetAEnd;
      }
      if(pathFromA && pathFromB){
        endNodeMergeGapOffset = (endNodeMergeGapOffsetA +
                                 (sortGapsPtr->gapSize.mean -
                                  (endNodeMergeGapOffsetB +
                                   mergeInsertSize.mean))) / 2.0;
      }else if(pathFromA){
        endNodeMergeGapOffset = endNodeMergeGapOffsetA;
      }else if(pathFromB){
        endNodeMergeGapOffset = sortGapsPtr->gapSize.mean -
          (endNodeMergeGapOffsetB + mergeInsertSize.mean);
      }else{
        assert(FALSE);// must have pathFromAorB
      }
      if(orientAnchor == B_A){// traversing from A end
        intervalGapBegin = anchorOffset.mean -
          (endNodeMergeGapOffset + mergeGapSize.mean +
           newEndNodeMerge->bpLength.mean);
        intervalBegin = intervalGapBegin -
          (3.0 * sqrt(gapVariance + mergeGapSize.variance +
                      newEndNodeMerge->bpLength.variance));
        intervalGapEnd = anchorOffset.mean - (endNodeMergeGapOffset +
                                              mergeGapSize.mean);
        intervalEnd = intervalGapEnd +
          (3.0 * sqrt(gapVariance + mergeGapSize.variance));
        if(intervalBegin > lastAnchorOffset.mean){
          if(verbose){
            fprintf(GlobalData->stderrc, "newEndNodeMerge did not fit in any of scaffoldAnchor's gaps.\n");
          }
          return(FALSE);
        }
        if(intervalEnd > lastAnchorOffset.mean){
          intervalEnd = lastAnchorOffset.mean - 1.0;
        }
      }else{// traversing from B end
        intervalGapBegin = anchorOffset.mean +
          (endNodeMergeGapOffset + mergeGapSize.mean +
           newEndNodeMerge->bpLength.mean);
        intervalBegin = intervalGapBegin +
          (3.0 * sqrt(gapVariance + mergeGapSize.variance +
                      newEndNodeMerge->bpLength.variance));
        intervalGapEnd = anchorOffset.mean + (endNodeMergeGapOffset +
                                              mergeGapSize.mean);
        intervalEnd = intervalGapEnd -
          (3.0 * sqrt(gapVariance + mergeGapSize.variance));
        if(intervalBegin < lastAnchorOffset.mean){
          if(verbose){
            fprintf(GlobalData->stderrc, "newEndNodeMerge did not fit in any of scaffoldAnchor's gaps.\n");
          }
          return(FALSE);
        }
        if(intervalEnd < lastAnchorOffset.mean){
          intervalEnd = lastAnchorOffset.mean + 1.0;
        }
      }
      if(FindScaffoldMerge(graph, scaffoldAnchor, scaffoldMerge,
                           orientAnchor, orientMerge, newEndNodeMerge,
                           orientNewEndNodeMerge, intervalBegin, intervalEnd,
                           intervalGapBegin, intervalGapEnd, FALSE,
                           maxMergeGapSize, newMergeGapVariance,
                           lastAnchorOffset, orderAnchor, verbose)){
        safe_free(sortGaps);
        SetOrderInfoForGap(orderAnchor, (sortGapsPtr->nodeA)->indexInScaffold,
                           firstTime, mergeGapVariance, sortGapsPtr->gapSize,
                           newMergeGapVariance, endNodeMergeGapOffsetA,
                           pathFromA, mergeInsertSize, endNodeMergeGapOffsetB,
                           pathFromB, numNodesInGap, FALSE);
        return(TRUE);
      }
    }
  }
  safe_free(sortGaps);
  return(FALSE);
}

#define EMIT_STATS 1

static  FILE *fMergeWeight, *fInterleavedMergeWeight;
static  FILE *fMergeDistance, *fInterleavedMergeDistance;



/****************************************************************************/
void OpenStatsFiles(void ) {
  char buffer[2048];

  AS_UTL_mkdir("stat");

  sprintf(buffer,"stat/ScaffoldMerge.mergeWeight.cgm");
  fMergeWeight = fopen(buffer,"w");
  AssertPtr(fMergeWeight);
  fprintf(fMergeWeight,"Non-Interleaved Scaffold Weight edge weights\n");
  
  sprintf(buffer,"stat/ScaffoldMerge.interleavedMergeWeight.cgm");
  fInterleavedMergeWeight = fopen(buffer,"w");
  AssertPtr(fInterleavedMergeWeight);
  fprintf(fInterleavedMergeWeight,"Interleaved Scaffold Weight edge weights\n");
  
  sprintf(buffer,"stat/ScaffoldMerge.mergeDistance.cgm");
  fMergeDistance = fopen(buffer,"w");
  AssertPtr(fMergeDistance);
  fprintf(fMergeDistance,"Non-Interleaved Scaffold Weight edge distances\n");
  
  sprintf(buffer,"stat/ScaffoldMerge.interleavedMergeDistance.cgm");
  fInterleavedMergeDistance = fopen(buffer,"w");
  AssertPtr(fInterleavedMergeDistance);
  fprintf(fInterleavedMergeDistance,"Interleaved Scaffold Weight edge distances\n");
  
}

/****************************************************************************/
void CloseStatsFiles(void){
  fclose(fMergeWeight);
  fclose(fMergeDistance);
  fclose(fInterleavedMergeWeight);
  fclose(fInterleavedMergeDistance);
}


/****************************************************************************/
#define CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD 2

int TouchesMarkedScaffolds(SEdgeT *curEdge){
  CIScaffoldT *scaffoldA = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idA);
  CIScaffoldT *scaffoldB = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idB);
  if(curEdge->orient == AB_AB){
    if((scaffoldA->numEssentialB == 0) &&
       (scaffoldB->numEssentialA == 0)){
      assert(scaffoldA->essentialEdgeB == NULLINDEX);
      assert(scaffoldB->essentialEdgeA == NULLINDEX);
      return FALSE;
    }else{
      return TRUE;
    }
  }else if(curEdge->orient == AB_BA){
    if((scaffoldA->numEssentialB == 0) &&
       (scaffoldB->numEssentialB == 0)){
      assert(scaffoldA->essentialEdgeB == NULLINDEX);
      assert(scaffoldB->essentialEdgeB == NULLINDEX);
      return FALSE;
    }else{
      return TRUE;
    }
  }else if(curEdge->orient == BA_AB){
    if((scaffoldA->numEssentialA == 0) &&
       (scaffoldB->numEssentialA == 0)){
      assert(scaffoldA->essentialEdgeA == NULLINDEX);
      assert(scaffoldB->essentialEdgeA == NULLINDEX);
      return FALSE;
    }else{
      return TRUE;
    }
  }else{//curEdge->orient == BA_BA
    if((scaffoldA->numEssentialA == 0) &&
       (scaffoldB->numEssentialB == 0)){
      assert(scaffoldA->essentialEdgeA == NULLINDEX);
      assert(scaffoldB->essentialEdgeB == NULLINDEX);
      return FALSE;
    }else{
      return TRUE;
    }
  }
}


#define EDGE_QUANTA 5
#define OVERLAP_QUANTA -10000.

int isLargeOverlapSize(LengthT *overlap, int32 numEdges){
  double maxGap, minGap;
  double numEdgeQuanta = numEdges/EDGE_QUANTA;
  double numOverlapQuanta;
  
  maxGap = overlap->mean + 3.0 * sqrt(overlap->variance);
  minGap = overlap->mean - 3.0 * sqrt(overlap->variance);
  
  numOverlapQuanta = minGap/OVERLAP_QUANTA;
  
  if(maxGap > OVERLAP_QUANTA)
    return FALSE;
  
  return(numEdgeQuanta < numOverlapQuanta);
}

int isLargeOverlap(CIEdgeT *curEdge){
  int ret =  isLargeOverlapSize(&(curEdge->distance), curEdge->edgesContributing);
  
  return ret;
  
}

#define EDGE_STRENGTH_FACTOR  2

/*
  it would be better to do this during scaffold edge creation
*/
int ThereIsAStrongerEdgeToSameScaffold(CDS_CID_t scfIID, SEdgeT * curSEdge)
{
  SEdgeTIterator SEdges;
  SEdgeT * sEdge;
  int32 orientValue;
  CDS_CID_t otherScaffoldID;
  int32 retVal = 0;
  
  otherScaffoldID = (scfIID == curSEdge->idA) ? curSEdge->idB : curSEdge->idA;
  orientValue =
    (curSEdge->orient == AB_AB || curSEdge->orient == BA_BA) ? 1 : -1;
  
  /*
    iterate over otherEnd's merged edges to see if there is a
    stronger one to the other scaffold
  */
  InitSEdgeTIterator(ScaffoldGraph, scfIID,
                     FALSE, FALSE, ALL_END, FALSE, &SEdges);
  while((sEdge = NextSEdgeTIterator(&SEdges)) != NULL)
    {
    
      if(sEdge->idA == otherScaffoldID ||
         (sEdge->idB == otherScaffoldID && sEdge != curSEdge))
        {
          ChunkOrientationType newOrientValue =
            (sEdge->orient == AB_AB || sEdge->orient == BA_BA) ? 1 : -1;
      
          /*
            Does the pair of edges agree (shift) or disagree (reversal)?
            If reversal, assume greater weight one is correct
            If shift,
            if both are positive, take one with greater weight
            if one is negative, take the postive one
            if both are negative, take the shortest one (smallest overlap)
          */
          if(orientValue != newOrientValue)
            {
              // reversal
              if(sEdge->edgesContributing > curSEdge->edgesContributing)
                {
                  if(retVal == 0)
                    {
                      fprintf(GlobalData->stderrc,
                              "SCF MERGE CONFLICT: " F_CID "," F_CID "  %s  %dbp  %dvar  %dec\n",
                              curSEdge->idA, curSEdge->idB,
                              ((curSEdge->orient == AB_AB) ? "AB_AB" :
                               ((curSEdge->orient == AB_BA) ? "AB_BA" :
                                ((curSEdge->orient == BA_AB) ? "BA_AB" : "BA_BA"))),
                              (int) curSEdge->distance.mean,
                              (int) curSEdge->distance.variance,
                              curSEdge->edgesContributing);
                    }
                  fprintf(GlobalData->stderrc, "\t" F_CID "," F_CID ", %s, %dbp  %dvar  %dec\n",
                          sEdge->idA, sEdge->idB,
                          ((sEdge->orient == AB_AB) ? "AB_AB" :
                           ((sEdge->orient == AB_BA) ? "AB_BA" :
                            ((sEdge->orient == BA_AB) ? "BA_AB" : "BA_BA"))),
                          (int) sEdge->distance.mean,
                          (int) curSEdge->distance.variance,
                          sEdge->edgesContributing);
                  retVal++;
                }
            }
          else
            {
              // shift
              if(curSEdge->distance.mean > 0)
                {
                  if(sEdge->distance.mean > 0)
                    {
                      // both are positive, prefer stronger
                      if(sEdge->edgesContributing > curSEdge->edgesContributing)
                        {
                          if(retVal == 0)
                            {
                              fprintf(GlobalData->stderrc,
                                      "SCF MERGE CONFLICT: " F_CID "," F_CID "  %s  %dbp  %dvar  %dec\n",
                                      curSEdge->idA, curSEdge->idB,
                                      ((curSEdge->orient == AB_AB) ? "AB_AB" :
                                       ((curSEdge->orient == AB_BA) ? "AB_BA" :
                                        ((curSEdge->orient == BA_AB) ? "BA_AB" : "BA_BA"))),
                                      (int) curSEdge->distance.mean,
                                      (int) curSEdge->distance.variance,
                                      curSEdge->edgesContributing);
                            }
                          fprintf(GlobalData->stderrc, "\t" F_CID "," F_CID ", %s, %dbp  %dvar  %dec\n",
                                  sEdge->idA, sEdge->idB,
                                  ((sEdge->orient == AB_AB) ? "AB_AB" :
                                   ((sEdge->orient == AB_BA) ? "AB_BA" :
                                    ((sEdge->orient == BA_AB) ? "BA_AB" : "BA_BA"))),
                                  (int) sEdge->distance.mean,
                                  (int) curSEdge->distance.variance,
                                  sEdge->edgesContributing);
                          retVal++;
                        }
                    }
                  else
                    {
                      /*
                        curSEdge is positive, sEdge is negative
                        prefer curSEdge unless sEdge is much stronger
                      */
                      if(sEdge->edgesContributing >
                         EDGE_STRENGTH_FACTOR * curSEdge->edgesContributing)
                        {
                          if(retVal == 0)
                            {
                              fprintf(GlobalData->stderrc,
                                      "SCF MERGE CONFLICT: " F_CID "," F_CID "  %s  %dbp  %dvar  %dec\n",
                                      curSEdge->idA, curSEdge->idB,
                                      ((curSEdge->orient == AB_AB) ? "AB_AB" :
                                       ((curSEdge->orient == AB_BA) ? "AB_BA" :
                                        ((curSEdge->orient == BA_AB) ? "BA_AB" : "BA_BA"))),
                                      (int) curSEdge->distance.mean,
                                      (int) curSEdge->distance.variance,
                                      curSEdge->edgesContributing);
                            }
                          fprintf(GlobalData->stderrc, "\t" F_CID "," F_CID ", %s, %dbp  %dvar  %dec\n",
                                  sEdge->idA, sEdge->idB,
                                  ((sEdge->orient == AB_AB) ? "AB_AB" :
                                   ((sEdge->orient == AB_BA) ? "AB_BA" :
                                    ((sEdge->orient == BA_AB) ? "BA_AB" : "BA_BA"))),
                                  (int) sEdge->distance.mean,
                                  (int) curSEdge->distance.variance,
                                  sEdge->edgesContributing);
                          retVal++;
                        }
                    }
                }
              else if(sEdge->distance.mean > 0)
                {
                  /*
                    sEdge is positive, curSEdge is negative
                    prefer sEdge unless curSEdge is much stronger
                  */
                  if(curSEdge->edgesContributing <
                     EDGE_STRENGTH_FACTOR * sEdge->edgesContributing)
                    {
                      if(retVal == 0)
                        {
                          fprintf(GlobalData->stderrc,
                                  "SCF MERGE CONFLICT: " F_CID "," F_CID "  %s  %dbp  %dvar  %dec\n",
                                  curSEdge->idA, curSEdge->idB,
                                  ((curSEdge->orient == AB_AB) ? "AB_AB" :
                                   ((curSEdge->orient == AB_BA) ? "AB_BA" :
                                    ((curSEdge->orient == BA_AB) ? "BA_AB" : "BA_BA"))),
                                  (int) curSEdge->distance.mean,
                                  (int) curSEdge->distance.variance,
                                  curSEdge->edgesContributing);
                        }
                      fprintf(GlobalData->stderrc, "\t" F_CID "," F_CID ", %s, %dbp  %dvar  %dec\n",
                              sEdge->idA, sEdge->idB,
                              ((sEdge->orient == AB_AB) ? "AB_AB" :
                               ((sEdge->orient == AB_BA) ? "AB_BA" :
                                ((sEdge->orient == BA_AB) ? "BA_AB" : "BA_BA"))),
                              (int) sEdge->distance.mean,
                              (int) curSEdge->distance.variance,
                              sEdge->edgesContributing);
                      retVal++;
                    }
                }
              else
                {
                  /*
                    both negative
                    prefer shorter overlap unless longer is much stronger
                    prefer sEdge if much stronger, or shorter & strong enough
                  */
                  if(sEdge->edgesContributing >
                     EDGE_STRENGTH_FACTOR * curSEdge->edgesContributing ||
                     (sEdge->distance.mean > curSEdge->distance.mean &&
                      EDGE_STRENGTH_FACTOR * sEdge->edgesContributing >
                      curSEdge->edgesContributing))
                    {
                      if(retVal == 0)
                        {
                          fprintf(GlobalData->stderrc,
                                  "SCF MERGE CONFLICT: " F_CID "," F_CID "  %s  %dbp  %dvar  %dec\n",
                                  curSEdge->idA, curSEdge->idB,
                                  ((curSEdge->orient == AB_AB) ? "AB_AB" :
                                   ((curSEdge->orient == AB_BA) ? "AB_BA" :
                                    ((curSEdge->orient == BA_AB) ? "BA_AB" : "BA_BA"))),
                                  (int) curSEdge->distance.mean,
                                  (int) curSEdge->distance.variance,
                                  curSEdge->edgesContributing);
                        }
                      fprintf(GlobalData->stderrc, "\t" F_CID "," F_CID ", %s, %dbp  %dvar  %dec\n",
                              sEdge->idA, sEdge->idB,
                              ((sEdge->orient == AB_AB) ? "AB_AB" :
                               ((sEdge->orient == AB_BA) ? "AB_BA" :
                                ((sEdge->orient == BA_AB) ? "BA_AB" : "BA_BA"))),
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


int OtherEndHasStrongerEdgeToSameScaffold(CDS_CID_t scfIID, SEdgeT * curSEdge)
{
  SEdgeTIterator SEdges;
  SEdgeT * sEdge;
  int otherEnd;
  CDS_CID_t otherScaffoldID;
  
  if(scfIID == curSEdge->idA)
    {
      otherEnd =
        (curSEdge->orient == AB_AB || curSEdge->orient == AB_BA) ? A_END : B_END;
      otherScaffoldID = curSEdge->idB;
    }
  else
    {
      otherEnd =
        (curSEdge->orient == AB_BA || curSEdge->orient == BA_BA) ? A_END : B_END;
      otherScaffoldID = curSEdge->idA;
    }
  
  /*
    iterate over otherEnd's merged edges to see if there is a
    stronger one to curSEdge->idB
  */
  InitSEdgeTIterator(ScaffoldGraph, scfIID,
                     FALSE, FALSE, otherEnd, FALSE, &SEdges);
  while((sEdge = NextSEdgeTIterator(&SEdges)) != NULL)
    {
      if((sEdge->idA == otherScaffoldID || sEdge->idB == otherScaffoldID) &&
         sEdge->edgesContributing > curSEdge->edgesContributing)
        return 1;
    }
  return 0;
}


// Find all merge candidates incident on scaffoldA
// Returns TRUE if marked edges are encountered
//
int FindAllMergeCandidates(VA_TYPE(PtrT) *sEdges,
                           VA_TYPE(PtrT) *overlapSEdges,
                           CIScaffoldT *fromScaffold,
                           int fromEnd,
                           CIScaffoldT *ignoreToScaffold,
                           int canonicalOnly, 
                           int minWeight, int doInterleaving,
                           int verbose){
  SEdgeTIterator SEdges;
  SEdgeT *curSEdge;
  CIScaffoldT *otherScaffold;
  CDS_CID_t otherScaffoldID;
  
  if(verbose){
    fprintf(GlobalData->stderrc,"*FindAllMergeCandidates from scaffold " F_CID " end:%d ignore:" F_CID " canonicalOnly:%d minWeight:%d\n",
            fromScaffold->id, fromEnd, (ignoreToScaffold?ignoreToScaffold->id:NULLINDEX), canonicalOnly, minWeight);
  }
  InitSEdgeTIterator(ScaffoldGraph, fromScaffold->id, FALSE, FALSE, fromEnd, FALSE,
                     &SEdges);
  while((curSEdge = NextSEdgeTIterator(&SEdges)) != NULL){
    
    //PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph, "curSEdge", curSEdge, curSEdge->idA);
    
    if(curSEdge->flags.bits.isBogus){ // This edge has already been visited by the recursion
      if(verbose)fprintf(GlobalData->stderrc,"* Hit isBogus edge...continuing\n");
      continue;
    }
    
    if( curSEdge->idA != fromScaffold->id){
      if(canonicalOnly){
        continue;
      }
      otherScaffoldID = curSEdge->idA;
    }else{
      otherScaffoldID = curSEdge->idB;
    }
    otherScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                 otherScaffoldID);
    if(otherScaffold->flags.bits.smoothSeenAlready)
      continue;
    
    if(otherScaffold == ignoreToScaffold){
      if(verbose)fprintf(GlobalData->stderrc,"* otherScaffold " F_CID " == ignoreScaffold...continue\n", otherScaffold->id);
      continue;
    }
    
    if(TouchesMarkedScaffolds(curSEdge)){
      if(verbose)fprintf(GlobalData->stderrc,"* Edge (" F_CID "," F_CID ",%c) touches marked scaffolds\n",
                         curSEdge->idA, curSEdge->idB, curSEdge->orient);
      return TRUE;
      continue;
    }
    
#ifdef OTHER_END_CHECK
    if(OtherEndHasStrongerEdgeToSameScaffold(fromScaffold->id, curSEdge))
      continue;
#endif
    
#ifdef GENERAL_STRONGER_CHECK
    if(ThereIsAStrongerEdgeToSameScaffold(fromScaffold->id, curSEdge))
      continue;
#endif
    
    if(curSEdge->flags.bits.isDeleted ||
       isDeadCIScaffoldT(otherScaffold) ||
       otherScaffold->type != REAL_SCAFFOLD)
      continue;
    
    assert((curSEdge->idA != NULLINDEX) && (curSEdge->idB != NULLINDEX));
    if((curSEdge->edgesContributing - (isOverlapEdge(curSEdge) ? 1 : 0)) < minWeight){
      //      if(verbose)fprintf(GlobalData->stderrc,"* Edge too weak, less than weight %d\n", minWeight);
      continue;
    }
    
    if(!doInterleaving && isLargeOverlap(curSEdge)){
      if(overlapSEdges &&
         curSEdge->edgesContributing >= CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD)
        
        AppendPtrT(overlapSEdges, (const void *)&curSEdge);
      continue;
    }
    
    AppendPtrT(sEdges, (const void *)&curSEdge);
    if(verbose){
      PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph, "\t", curSEdge, fromScaffold->id);
      fflush(GlobalData->stderrc);
    }
  }
  
  return FALSE;
}



// Initialize the numEssentialX and essentialEdgeX fields in all scaffolds
// Collect a list of all "confirmed" scaffold edges thata re candidates for 
// merging.
void SortSEdges(VA_TYPE(PtrT) * sEdges, int32 verbose)
{
  SEdgeT ** sEdge = (SEdgeT **)GetPtrT(sEdges,0);
  
  if(verbose)
    {
      fprintf(GlobalData->stderrc,
              "* Sorting %d sEdges of size " F_SIZE_T " sEdge = %p\n",
              (int) GetNumPtrTs(sEdges), sizeof(*sEdge), sEdge);
      fflush(NULL);
    }
  
#ifdef SORT_BY_EDGE_WEIGHTS
  qsort((void *)sEdge,
        GetNumPtrTs(sEdges),
        sizeof(*sEdge),
        CompareSEdgesContributing);
#else
  qsort((void *)sEdge,
        GetNumPtrTs(sEdges),
        sizeof(*sEdge),
        CompareSEdgeGaps);
#endif
  
  if(verbose)
    {
      int i;
      fprintf(GlobalData->stderrc,">>>>>>>> Usable Edges >>>>>>>\n");
      for(i = 0; i < GetNumPtrTs(sEdges); i++)
        {
          SEdgeT * edge = (SEdgeT *) *GetPtrT(sEdges, i);
          PrintGraphEdge(GlobalData->stderrc,
                         ScaffoldGraph->ScaffoldGraph,
                         "usable> ",
                         edge, edge->idA);	
        }
      fprintf(GlobalData->stderrc,">>>>>>>>>>>>>>>>>>>>>\n");
    }
}


void BuildUsableSEdges(VA_TYPE(PtrT) *sEdges,
                       VA_TYPE(PtrT) *overlapSEdges,
                       int32 verbose)
{
  int32 numRealScaffolds = 0;
  int32 numScaffolds = 0;
  int32 numSEdges;
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;
  
  AssertPtr(sEdges);
  
  numScaffolds = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
  numSEdges = GetNumGraphEdges(ScaffoldGraph->ScaffoldGraph);
  
  InitGraphNodeIterator(&scaffolds,
                        ScaffoldGraph->ScaffoldGraph,
                        GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL)
    {
      if(isDeadCIScaffoldT(scaffold) || (scaffold->type != REAL_SCAFFOLD))
        continue;
    
      scaffold->numEssentialA = 0;
      scaffold->numEssentialB = 0;
      scaffold->essentialEdgeA = NULLINDEX;
      scaffold->essentialEdgeB = NULLINDEX;
      scaffold->setID = NULLINDEX;
      scaffold->flags.bits.smoothSeenAlready = 0;
      scaffold->flags.bits.walkedAlready = 0;
    }
  
  InitGraphNodeIterator(&scaffolds,
                        ScaffoldGraph->ScaffoldGraph,
                        GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL)
    {
      if(isDeadCIScaffoldT(scaffold) || (scaffold->type != REAL_SCAFFOLD))
        continue;
      numRealScaffolds++;
    
      FindAllMergeCandidates(sEdges, overlapSEdges,
                             scaffold,
                             ALL_END, NULL,
                             TRUE,
                             CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD,
                             GlobalData->doInterleavedScaffoldMerging,
                             verbose);
    }
  
  if(verbose){
    fprintf(GlobalData->stderrc, "numSEdges %d numScaffolds %d\nnumRealScaffolds %d numRealSEdges %d\n",
            numSEdges, numScaffolds, numRealScaffolds, (int) GetNumPtrTs(sEdges));
  }
  if(verbose)
    {
      int i;
      fprintf(GlobalData->stderrc,">>>>>>>> Overlap Edges >>>>>>>\n");
      for(i = 0; i < GetNumPtrTs(overlapSEdges); i++){
        SEdgeT *edge = (SEdgeT *) *GetPtrT(overlapSEdges, i);
        PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph, "overlap> ", edge, edge->idA);	
      }
      fprintf(GlobalData->stderrc,">>>>>>>>>>>>>>>>>>>>>\n");
    }
  if(GetNumPtrTs(sEdges) == 0){
    if(verbose){
      fprintf(GlobalData->stderrc, ">>>>No Usable SEdges!!! (there are %d overlap SEdges)\n", (int) GetNumPtrTs(overlapSEdges));
    }
    return;
  }
  
  SortSEdges(sEdges, verbose);
}


ChunkInstanceT * GetTrueBEndCI(ScaffoldGraphT * graph, CIScaffoldT * scaffold)
{
  ChunkInstanceT * returnContig =
    GetGraphNode(graph->ContigGraph, scaffold->info.Scaffold.BEndCI);
  CIScaffoldTIterator contigIterator;
  ChunkInstanceT * myContig;
  
  InitCIScaffoldTIterator(graph, scaffold, FALSE, FALSE, &contigIterator);
  while((myContig = NextCIScaffoldTIterator(&contigIterator)) != NULL)
    {
      if(MAX(returnContig->offsetAEnd.mean, returnContig->offsetBEnd.mean) <
         MAX(myContig->offsetAEnd.mean, myContig->offsetBEnd.mean))
        returnContig = myContig;
      else if(MIN(returnContig->offsetAEnd.mean, returnContig->offsetBEnd.mean) >
              MAX(myContig->offsetAEnd.mean, myContig->offsetBEnd.mean))
        break;
    }
  return returnContig;
}


/****************************************************************************/
// TranslateScaffoldOverlapToContigOverlap
//       Prepare the groundwork for computing an overlap between the appropriate extremal
//       contigs in a pair of scaffolds.
//
void TranslateScaffoldOverlapToContigOverlap(CIScaffoldT *scaffoldA, CIScaffoldT *scaffoldB, 
                                             ChunkOrientationType scaffoldEdgeOrient,
                                             NodeCGW_T **endNodeA, NodeCGW_T **endNodeB, 
                                             ChunkOrient *orientEndNodeA,  ChunkOrient *orientEndNodeB,
                                             ChunkOrientationType *edgeEndsOrient,
                                             double *extremalGapSizeA, double *extremalGapSizeB){
  NodeCGW_T *nextNodeA, *nextNodeB;
  int AGapTowardAEnd, BGapTowardAEnd;
#if 0
  fprintf(GlobalData->stderrc,"* Translate (" F_CID "," F_CID ",%c)   " F_CID " has (" F_CID "," F_CID ")  " F_CID " has (" F_CID "," F_CID ")\n",
          scaffoldA->id, scaffoldB->id, scaffoldEdgeOrient,  
          scaffoldA->id, scaffoldA->info.Scaffold.AEndCI, scaffoldA->info.Scaffold.BEndCI,
          scaffoldB->id, scaffoldB->info.Scaffold.AEndCI, scaffoldB->info.Scaffold.BEndCI);
#endif
  if(scaffoldEdgeOrient == AB_AB){
    // BEndCI may be contained
    *endNodeA = GetTrueBEndCI(ScaffoldGraph, scaffoldA);
    *endNodeB = GetGraphNode(ScaffoldGraph->ContigGraph,
                             scaffoldB->info.Scaffold.AEndCI);
    *orientEndNodeA = GetNodeOrient(*endNodeA);
    *orientEndNodeB = GetNodeOrient(*endNodeB);
    
    AGapTowardAEnd = TRUE;
    BGapTowardAEnd = FALSE;
  }else if(scaffoldEdgeOrient == AB_BA){
    // BendCI may be contained
    *endNodeA = GetTrueBEndCI(ScaffoldGraph, scaffoldA);
    // BendCI may be contained
    *endNodeB = GetTrueBEndCI(ScaffoldGraph, scaffoldB);
    *orientEndNodeA = GetNodeOrient(*endNodeA);
    *orientEndNodeB = FlipNodeOrient(GetNodeOrient(*endNodeB));
    AGapTowardAEnd = TRUE;
    BGapTowardAEnd = TRUE;
  }else if(scaffoldEdgeOrient == BA_AB){
    *endNodeA = GetGraphNode(ScaffoldGraph->ContigGraph,
                             scaffoldA->info.Scaffold.AEndCI);
    *endNodeB = GetGraphNode(ScaffoldGraph->ContigGraph,
                             scaffoldB->info.Scaffold.AEndCI);
    *orientEndNodeA = FlipNodeOrient(GetNodeOrient(*endNodeA));
    *orientEndNodeB = GetNodeOrient(*endNodeB);
    AGapTowardAEnd = FALSE;
    BGapTowardAEnd = FALSE;
  }else{//curEdge->orient == BA_BA
    *endNodeA = GetGraphNode(ScaffoldGraph->ContigGraph,
                             scaffoldA->info.Scaffold.AEndCI);
    *endNodeB = GetTrueBEndCI(ScaffoldGraph, scaffoldB);
    *orientEndNodeA = FlipNodeOrient(GetNodeOrient(*endNodeA));
    *orientEndNodeB = FlipNodeOrient(GetNodeOrient(*endNodeB));
    AGapTowardAEnd = FALSE;
    BGapTowardAEnd = TRUE;
  }
  if(*orientEndNodeA == A_B){
    if(*orientEndNodeB == A_B){
      *edgeEndsOrient = AB_AB;
    }else{//(orientEndNodeB == B_A
      *edgeEndsOrient = AB_BA;
    }
  }else{//(orientEndNodeA == B_A
    if(*orientEndNodeB == A_B){
      *edgeEndsOrient = BA_AB;
    }else{//(orientEndNodeB == B_A
      *edgeEndsOrient = BA_BA;
    }
  }
  nextNodeA = GetGraphNode(ScaffoldGraph->ContigGraph, (AGapTowardAEnd?(*endNodeA)->AEndNext:(*endNodeA)->BEndNext));
  nextNodeB = GetGraphNode(ScaffoldGraph->ContigGraph, (BGapTowardAEnd?(*endNodeB)->AEndNext:(*endNodeB)->BEndNext));
  if(nextNodeA){
    if(AGapTowardAEnd){
      *extremalGapSizeA = MIN((*endNodeA)->offsetAEnd.mean, (*endNodeA)->offsetBEnd.mean) -
        MAX(nextNodeA->offsetAEnd.mean, nextNodeA->offsetBEnd.mean);
    }else{
      *extremalGapSizeA = 
        MIN(nextNodeA->offsetAEnd.mean, nextNodeA->offsetBEnd.mean)-
        MAX((*endNodeA)->offsetAEnd.mean, (*endNodeA)->offsetBEnd.mean);
    }
  }else{
    *extremalGapSizeA = 1000000.0; // large number, anything will fit
  }
  if(nextNodeB){
    if(BGapTowardAEnd){
      *extremalGapSizeB = MIN((*endNodeB)->offsetAEnd.mean, (*endNodeB)->offsetBEnd.mean) -
        MAX(nextNodeB->offsetAEnd.mean, nextNodeB->offsetBEnd.mean);
    }else{
      *extremalGapSizeB = 
        MIN(nextNodeB->offsetAEnd.mean, nextNodeB->offsetBEnd.mean)-
        MAX((*endNodeB)->offsetAEnd.mean, (*endNodeB)->offsetBEnd.mean);
    }
  }else{
    *extremalGapSizeB = 1000000.0; // large number, anything will fit
  }
  
}



void MarkScaffoldsForMerging(SEdgeT *curEdge, int markForMerging){
  int edgeSelected = FALSE;
  CDS_CID_t edgeIndex = NULLINDEX;
  CIScaffoldT *scaffoldA = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idA);
  CIScaffoldT *scaffoldB = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idB);
  
  
  if(markForMerging){
    // ignore these two scaffolds for the rest of this iteration
    scaffoldB->flags.bits.walkedAlready  = 1;
    scaffoldA->flags.bits.walkedAlready  = 1;
    edgeIndex = (CDS_CID_t)GetVAIndex_SEdgeT(ScaffoldGraph->ScaffoldGraph->edges, curEdge);
  }
  
#if 1
  fprintf(GlobalData->stderrc,"* MarkScaffolds " F_CID " " F_CID " %c (%g,%g) for merging %d\n",
          scaffoldA->id, scaffoldB->id, curEdge->orient, curEdge->distance.mean, curEdge->distance.variance,markForMerging);
#endif
  if(curEdge->orient == AB_AB){
    if(markForMerging){
      assert(scaffoldA->essentialEdgeB == NULLINDEX);
      assert(scaffoldB->essentialEdgeA == NULLINDEX);
      scaffoldA->essentialEdgeB = edgeIndex;
      scaffoldB->essentialEdgeA = edgeIndex;
      edgeSelected = TRUE;
    }
#if 0
    fprintf(GlobalData->stderrc,"$$$$ Marking scaffold " F_CID "B as Not For Merging\n", scaffoldA->id);
    fprintf(GlobalData->stderrc,"$$$$ Marking scaffold " F_CID "A as Not For Merging\n", scaffoldB->id);
#endif
    scaffoldA->numEssentialB = 1;
    scaffoldB->numEssentialA = 1;
  }else if(curEdge->orient == AB_BA){
    if(markForMerging){
      assert(scaffoldA->essentialEdgeB == NULLINDEX);
      assert(scaffoldB->essentialEdgeB == NULLINDEX);
      scaffoldA->essentialEdgeB = edgeIndex;
      scaffoldB->essentialEdgeB = edgeIndex;
      edgeSelected = TRUE;
    }
#if 0
    fprintf(GlobalData->stderrc,"$$$$ Marking scaffold " F_CID "B as Not For Merging\n", scaffoldA->id);
    fprintf(GlobalData->stderrc,"$$$$ Marking scaffold " F_CID "B as Not For Merging\n", scaffoldB->id);
#endif
    scaffoldA->numEssentialB = 1;
    scaffoldB->numEssentialB = 1;
  }else if(curEdge->orient == BA_AB){
    if(markForMerging){
      assert(scaffoldA->essentialEdgeA == NULLINDEX);
      assert(scaffoldB->essentialEdgeA == NULLINDEX);
      scaffoldA->essentialEdgeA = edgeIndex;
      scaffoldB->essentialEdgeA = edgeIndex;
      edgeSelected = TRUE;
    }
    scaffoldA->numEssentialA = 1;
    scaffoldB->numEssentialA = 1;
#if 0
    fprintf(GlobalData->stderrc,"$$$$ Marking scaffold " F_CID "A as Not For Merging\n", scaffoldA->id);
    fprintf(GlobalData->stderrc,"$$$$ Marking scaffold " F_CID "A as Not For Merging\n", scaffoldB->id);
#endif
  }else{//curEdge->orient == BA_BA
    if(markForMerging){
      assert(scaffoldA->essentialEdgeA == NULLINDEX);
      assert(scaffoldB->essentialEdgeB == NULLINDEX);
      scaffoldA->essentialEdgeA = edgeIndex;
      scaffoldB->essentialEdgeB = edgeIndex;
      edgeSelected = TRUE;
    }
#if 0
    fprintf(GlobalData->stderrc,"$$$$ Marking scaffold " F_CID "A as Not For Merging\n", scaffoldA->id);
    fprintf(GlobalData->stderrc,"$$$$ Marking scaffold " F_CID "B as Not For Merging\n", scaffoldB->id);
#endif
    scaffoldA->numEssentialA = 1;
    scaffoldB->numEssentialB = 1;
  }
  if(edgeSelected){
#if EMIT_STATS
    //	    PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->ScaffoldGraph, " Non-Interleaved Merging based on ", curEdge, curEdge->idA);
    fprintf(fMergeWeight,"%d\n", curEdge->edgesContributing);
    fprintf(fMergeDistance,"%d\n", (int)curEdge->distance.mean);
#endif
  }
}


void MarkScaffoldForNotMerging(CIScaffoldT *scaffoldA){
  
#if 0
  fprintf(GlobalData->stderrc,"@@@@ Marking scaffold " F_CID " as Not For Merging\n", scaffoldA->id);
#endif
  scaffoldA->numEssentialB = 1;
  scaffoldA->numEssentialA = 1;
}


/*****************************************************************************/

int DoesScaffoldCFit(CIScaffoldT *scaffoldA,
                     CIScaffoldT *scaffoldB,
                     CIScaffoldT *scaffoldC,
                     SEdgeT *edgeAB,
                     SEdgeT *edgeToC,
                     VA_TYPE(PtrT) *edgesFromA,
                     VA_TYPE(PtrT) *edgesFromB,
                     int32 *overlapFound,
                     int checkForTinyScaffolds,
                     int32 verbose){
  int isEdgeFromA = (edgeToC->idA == scaffoldA->id) || (edgeToC->idB == scaffoldA->id);
  float chiSquaredValue;
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
  
  
  if(verbose){
    fprintf(GlobalData->stderrc,"* DoesScaffold " F_CID " fit of length %f fit between scaffolds (" F_CID "," F_CID ") in a gap of size %f +/- %f?\n",
            scaffoldC->id, scaffoldC->bpLength.mean, scaffoldA->id, scaffoldB->id, edgeAB->distance.mean, sqrt(edgeAB->distance.variance));
    DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffoldC, FALSE);
  }
  
  
  if(isEdgeFromA){
    otherEdgeList = edgesFromB;
    otherScaffold = scaffoldA;
  }else{
    otherEdgeList = edgesFromA;
    otherScaffold = scaffoldB;
  }   
  
  
  // We don't want to stick a teeny tiny element in the middle of a gap between two
  // giants.  Often this will later prevent the giants from merging.  So, we
  // prevent it.  We only place stuff such that the size of the element being placed
  // is at least 20% of the size of the implied gap.  So, for a 50k gap, we need a 10k
  // scaffold.  For a 10k gap, we need a 2k scaffold, etc.
  lengthC_to_dist = 1.0;
  if(edgeToC->distance.mean > 0){
    lengthC_to_dist = (scaffoldC->bpLength.mean / edgeToC->distance.mean);
  }
  
  if(checkForTinyScaffolds &&
     scaffoldC->bpLength.mean < 5000 &&
     lengthC_to_dist < 0.20){
    if(verbose){
      fprintf(GlobalData->stderrc,"Scaffold C is too short (%g) relative to edge length (%g)\n",
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
  
  if(verbose){
    fprintf(GlobalData->stderrc,"Gap = (%f +/- %f)  Alternate = (%f +/- %f)\n",
            gap.mean, sqrt(gap.variance), alternateGap.mean, sqrt(alternateGap.variance));
    
    fprintf(GlobalData->stderrc,"*DoesScaffoldCFit  CLength = %g  gap estimate = %g +/- %g\n",
            scaffoldC->bpLength.mean, gap.mean, gap.variance);
  }
  
  if(gap.mean - 3.0 * sqrt(gap.variance)< -5000){
    if(verbose){
      fprintf(GlobalData->stderrc,"* DoesScaffoldCFit fails on gap.mean < -5000  gap = %g\n", gap.mean);
    }
    return FALSE;
  }
  // Look through all of the edges in the list, and try to find one that fits
  //
  if(!edgeToC->flags.bits.isProbablyBogus){ //  we should find one if this flag isn't set
    for(i = 0; i < GetNumPtrTs(otherEdgeList); i++){
      otherEdgeToC = (SEdgeT *) *GetPtrT(otherEdgeList,i);
      if(otherEdgeToC->idA == scaffoldC->id ||
         otherEdgeToC->idB == scaffoldC->id ){
        break;
      }else{
        otherEdgeToC = NULL;
      }
    }
    assert(otherEdgeToC != NULL); // we should find one!
  }
  
  // If we found an edge from A (B) to C, then test its length
  if(otherEdgeToC){
    if(!PairwiseChiSquare((float)otherEdgeToC->distance.mean,    
                          (float)otherEdgeToC->distance.variance,
                          (float)gap.mean,
                          (float)gap.variance, (LengthT *)NULL,
                          &chiSquaredValue, SCAFFOLD_MERGE_CHI2_THRESHHOLD)){ // fails Chi-square
      if(verbose){
        fprintf(GlobalData->stderrc,"* 2 DoesScaffoldCFit fails pairwise chi square %f\n", chiSquaredValue);
      }
      return FALSE;
    }
  }
  
  
  // If we got here, either there is no edge between scaffoldC and scaffoldA or scaffoldB, or
  // the edge length is chi-square compatible with the gap estimate.
  // Now we refine further to handle the other cases.
  
  // If C can't overlap with A or B, we are done
  if((gap.mean - 3.0 * sqrt(gap.variance)) > - CGW_MISSED_OVERLAP){
    if(verbose){
      fprintf(GlobalData->stderrc,"* gap %g +/- %g ...can't overlap... returning TRUE\n", gap.mean, gap.variance);
    }
    return TRUE;
  }else{
    if(verbose)fprintf(GlobalData->stderrc,"* Checking overlap\n");
  }
  
  // C need not overlap with B
  if(PairwiseChiSquare((float)gap.mean,    
                       (float)gap.variance,
                       (float)-CGW_MISSED_OVERLAP,
                       (float)1.0, (LengthT *)NULL,
                       &chiSquaredValue, SCAFFOLD_MERGE_CHI2_THRESHHOLD)){ // passes Chi-square
    if(verbose)
      fprintf(GlobalData->stderrc,"* gap.mean is chi-squ compatible with  %g +/- %g ... returning TRUE\n", (float)-CGW_MISSED_OVERLAP, 1.0);
    needNotOverlap = TRUE;
    // If this is a one sided edge, see if the alternate placement is also compatible with
    // a no overlap scenario.  If so, bail.
    if(edgeToC->flags.bits.isProbablyBogus){
      if(PairwiseChiSquare((float)alternateGap.mean,    
                           (float)alternateGap.variance,
                           (float)-CGW_MISSED_OVERLAP,
                           (float)1.0, (LengthT *)NULL,
                           &chiSquaredValue, (float)SCAFFOLD_MERGE_CHI2_THRESHHOLD)){ // passes Chi-square
        if(verbose)
          fprintf(GlobalData->stderrc,"* One sided edge failed alternate chi square test\n");
        return FALSE;
      }
    }
    
  }else{
    needNotOverlap = FALSE;
    if(verbose)
      fprintf(GlobalData->stderrc,"* gap.mean is NOT chi-squ compatible with  %g +/- %g ... \n", (float)-CGW_MISSED_OVERLAP, 1.0);
  }
  
  
  // C MUST overlap
  // look for overlap
  {
    ChunkOrient orientEndNodeA, orientEndNodeB;
    ChunkOrientationType edgeEndsOrient;
    NodeCGW_T *endNodeA, *endNodeB;
    double aGapSize, bGapSize;
    int alternate;
    
    TranslateScaffoldOverlapToContigOverlap(otherScaffold, scaffoldC, edgeToC->orient, &endNodeA, &endNodeB, &orientEndNodeA, &orientEndNodeB, &edgeEndsOrient, 
                                            &aGapSize, &bGapSize);
    
    if(verbose)
      fprintf(GlobalData->stderrc,"* Looking for overlap nodeA:" F_CID " nodeB: " F_CID ", endAOrient:%c endBOrient:%c orient:%c distance:%g\n",
              endNodeA->id, endNodeB->id, orientEndNodeA, orientEndNodeB, edgeEndsOrient, edgeToC->distance.mean);
    
    overlapEdge = FindOverlapEdgeChiSquare(ScaffoldGraph, endNodeA, endNodeB->id,
                                           edgeEndsOrient, edgeToC->distance.mean,
                                           edgeToC->distance.variance, &chiSquaredValue,
                                           (float)SCAFFOLD_MERGE_CHI2_THRESHHOLD, &alternate, verbose);
    
    // assert(!alternate); // shouldn't get the wrong orientation, should we? What about interleaving?
    if(alternate)
      fprintf( GlobalData->stderrc, "Warning: got an alternate edge orientation in DoesScaffoldCFit!!!!\n");
  }
  if(!overlapEdge){
    if(verbose)
      fprintf(GlobalData->stderrc,"* 3 DoesScaffoldCFit fails %s overlap and doesn't\n",(needNotOverlap?"need not ":" must "));
    if(!needNotOverlap)
      {
        SaveEdgeMeanForLater(edgeToC);
        edgeToC->distance.mean = - CGW_MISSED_OVERLAP;
      }
    return needNotOverlap;//  C Doesn't Fit
  }
  *overlapFound = TRUE;
  if(verbose)
    PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->ContigGraph, " Overlap Found! ", overlapEdge, overlapEdge->idA);
  
  // If overlap is X-square with gap ==> it fits
  if(PairwiseChiSquare(gap.mean,    
                       gap.variance,
                       overlapEdge->distance.mean,
                       overlapEdge->distance.variance,
                       (LengthT *)NULL,
                       &chiSquaredValue, (float)SCAFFOLD_MERGE_CHI2_THRESHHOLD)){ // passes Chi-square
    
    SaveEdgeMeanForLater(edgeToC);
    edgeToC->distance.mean = overlapEdge->distance.mean;
    
    if(verbose)
      fprintf(GlobalData->stderrc,"* Overlap is chi-sq OK...returning TRUE\n");
    return TRUE;
  }else{
    if(verbose)fprintf(GlobalData->stderrc,"* Overlap is NOT chi-sq OK.\n");
    return needNotOverlap;
  }
  
  return FALSE;
}

/********************************************************************/
static CDS_CID_t thisID = NULLINDEX;
static int CompareSEdgesByOtherID(const void *c1, const void *c2){
  SEdgeT *s1 = *(SEdgeT **)c1;
  SEdgeT *s2 = *(SEdgeT **)c2;
  CDS_CID_t other1 = (s1->idA == thisID?s1->idB:s1->idA);
  CDS_CID_t other2 = (s2->idA == thisID?s2->idB:s2->idA);
  int32 diff = other2 - other1;
  
  if(diff > 0)
    return (int)1;
  else if (diff < 0)
    return -1;
  
  diff = s1->edgesContributing - s2->edgesContributing;
  
  if(diff > 0)
    return (int)1;
  else if (diff < 0)
    return -1;
  
  return 0;
  
  
}

/********************************************************************/
void SortSEdgesByOtherID(VA_TYPE(PtrT) *edges, CDS_CID_t ThisID){
  thisID = ThisID;
  qsort((void *)GetPtrT(edges,0),GetNumPtrTs(edges),
        sizeof(PtrT), CompareSEdgesByOtherID);
}


/********************************************************************/
void ZeroEdgeWeights(VA_TYPE(PtrT) *edges){
  SEdgeT *edge;
  int i;
  for(i = 0; i < GetNumPtrTs(edges); i++){
    edge = (SEdgeT *) *GetPtrT(edges,i);
    edge->quality = 0.0;
    edge->flags.bits.isBogus = FALSE;
    edge->flags.bits.isProbablyBogus = FALSE;
  }
}
/********************************************************************/
int  AssignEdgeWeights(VA_TYPE(PtrT) *mergedEdges,
                       SEdgeT *curEdge,
                       CIScaffoldT *scaffoldA,
                       VA_TYPE(PtrT) *edgesFromA,
                       CIScaffoldT *scaffoldB,
                       VA_TYPE(PtrT) *edgesFromB,
                       int doInterleaving,
                       int verbose){
  SEdgeT *edge;  // Edge from A->C
  int i;
  
  if(verbose)
    fprintf(GlobalData->stderrc,"* Assign Edge Weights\n");
  
  // Iterate through edgesFromA. (A,C)  These are sorted by increasing otherID and then by decreasing weight
  if(verbose)
    fprintf(GlobalData->stderrc,"* Looking through %d edgesFromA\n",
            (int) GetNumPtrTs(edgesFromA));
  
  for(i = 0; i < GetNumPtrTs(edgesFromA); i++){
    int j;
    CDS_CID_t otherID;
    int foundBC = FALSE;
    edge = (SEdgeT *) *GetPtrT(edgesFromA,i);
    if(edge->idA == scaffoldA->id){
      otherID = edge->idB;
    }else{
      otherID = edge->idA;
    }
    if(!doInterleaving && isLargeOverlap(edge))
      continue;
    
    if(verbose){
      fprintf(GlobalData->stderrc,"* Found an edge AC (" F_CID "," F_CID ",%c)\n",
              scaffoldA->id, otherID, GetEdgeOrientationWRT(edge, scaffoldA->id));
      
      
      fprintf(GlobalData->stderrc,"* Looking through %d edgesFromB\n",
              (int) GetNumPtrTs(edgesFromB));
    }
    // For each edge, look through edgesFromB for an edge to the same Contig from (B,C)
    // If present, weight = weight(A,C), add the heavier of (A,C) (B,C) to mergedEdges.
    // If absent, add (A,C) to mergedEdges
    for(j = 0; j < GetNumPtrTs(edgesFromB); j++){
      CDS_CID_t otherBID;
      SEdgeT *edgeB = (SEdgeT *) *GetPtrT(edgesFromB,j);
      
      if(verbose)PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph, "EdgeFromBC? ", edgeB, edgeB->idA);
      
      if(edgeB->quality > 0.0){
        if(verbose)fprintf(GlobalData->stderrc,"*Skipping (quality = %g)\n", edgeB->quality);
        continue;   // This edge is already spoken for
      }
      if(edgeB->idA == scaffoldB->id)
        otherBID = edgeB->idB;
      else
        otherBID = edgeB->idA;
      
      // Make sure we have orientation correct
      if(otherBID == otherID){// We have (A,B), (A,C) AND (C,B)
        switch(GetEdgeOrientationWRT(curEdge, scaffoldA->id)){
          case AB_AB:
            switch(GetEdgeOrientationWRT(edge,scaffoldA->id)){ // orientation of AC
              case AB_BA:
                // ------------>   <----------------    ------------>
                //      A                  C                  B
                if(GetEdgeOrientationWRT(edgeB, scaffoldB->id) != BA_AB){ // orientation of BC
                  if(verbose)fprintf(GlobalData->stderrc,"* 1 Orientation is %c should be BA_AB\n",
                                     GetEdgeOrientationWRT(edgeB,scaffoldB->id));
                  continue;
                }
                break;
              case AB_AB:
                // ------------>   ---------------->    ------------>
                //      A                  C                  B
                if(GetEdgeOrientationWRT(edgeB, scaffoldB->id) != BA_BA){
                  if(verbose)fprintf(GlobalData->stderrc,"* 1 Orientation is %c should be BA_BA\n",
                                     GetEdgeOrientationWRT(edgeB,scaffoldB->id));
                  continue;
                }
                break;
              case BA_BA:
                if(verbose)
                  fprintf(GlobalData->stderrc,"* 1 Wrong orientation BA_BA\n");
                continue;  // wrong orientation
                
              case BA_AB:
                if(verbose)
                  fprintf(GlobalData->stderrc,"* 1 Wrong orientation BA_AB\n");
                continue;  // wrong orientation
              default:
                assert(0);
                break;
            }
            break;
          case AB_BA:
            switch(GetEdgeOrientationWRT(edge,scaffoldA->id)){
              case AB_BA:
                // ------------>   <----------------    <------------
                //      A                  C                  B
                if(GetEdgeOrientationWRT(edgeB, scaffoldB->id) != AB_AB){
                  if(verbose)
                    fprintf(GlobalData->stderrc,"* 2 Orientation is %c should be AB_AB\n",
                            GetEdgeOrientationWRT(edgeB,scaffoldB->id));
                  continue;
                }
                break;
              case AB_AB:
                // ------------>   ---------------->    <------------
                //      A                  C                  B
                if(GetEdgeOrientationWRT(edgeB, scaffoldB->id) != AB_BA){
                  if(verbose)
                    fprintf(GlobalData->stderrc,"* 2 Orientation is %c should be AB_BA\n",
                            GetEdgeOrientationWRT(edgeB,scaffoldB->id));
                  continue;
                }
                break;
              case BA_BA:
                if(verbose)
                  fprintf(GlobalData->stderrc,"* 2 Wrong orientation BA_BA\n");
                continue;  // wrong orientation
              case BA_AB:
                if(verbose)
                  fprintf(GlobalData->stderrc,"* 2 Wrong orientation BA_AB\n");
                continue;  // wrong orientation
              default:
                assert(0);
                break;
            }
            break;
          case BA_AB:
            switch(GetEdgeOrientationWRT(edge,scaffoldA->id)){
              case BA_BA:
                // <------------   <----------------    ------------>
                //      A                  C                  B
                if(GetEdgeOrientationWRT(edgeB, scaffoldB->id) != BA_AB){
                  if(verbose)
                    fprintf(GlobalData->stderrc,"* 3 Orientation is %c should be BA_AB\n",
                            GetEdgeOrientationWRT(edgeB,scaffoldB->id));
                  continue;
                }
                break;
              case BA_AB:
                // <------------   ---------------->    ------------>
                //      A                  C                  B
                if(GetEdgeOrientationWRT(edgeB, scaffoldB->id) != BA_BA){
                  if(verbose)
                    fprintf(GlobalData->stderrc,"* 3 Orientation is %c should be BA_BA\n",
                            GetEdgeOrientationWRT(edgeB,scaffoldB->id));
                  continue;
                }
                break;
              case AB_BA:
                if(verbose)
                  fprintf(GlobalData->stderrc,"* 3 Wrong orientation AB_BA\n");
                continue;
              case AB_AB:
                if(verbose)
                  fprintf(GlobalData->stderrc,"* 3 Wrong orientation AB_AB\n");
                continue;  // wrong orientation
              default:
                assert(0);
                break;
            }
            break;
          case BA_BA:
            switch(GetEdgeOrientationWRT(edge,scaffoldA->id)){ // orientation of AC
              case BA_BA:
                // <------------   <----------------    <-----------
                //      A                  C                  B
                if(GetEdgeOrientationWRT(edgeB, scaffoldB->id) != AB_AB){ // BC
                  if(verbose)
                    fprintf(GlobalData->stderrc,"* rOrientation is %c should be AB_AB\n",
                            GetEdgeOrientationWRT(edgeB,scaffoldB->id));
                  continue;
                }
                break;
              case BA_AB:
                // <------------   ---------------->     <-----------
                //      A                  C                  B
                if(GetEdgeOrientationWRT(edgeB, scaffoldB->id) != AB_AB){ // BC
                  if(verbose)
                    fprintf(GlobalData->stderrc,"* 4 Orientation is %c should be AB_AB\n",
                            GetEdgeOrientationWRT(edgeB,scaffoldB->id));
                  continue;
                }
                break;
              case AB_AB:
                if(verbose)
                  fprintf(GlobalData->stderrc,"* 4 Wrong orientation AB_AB\n");
                continue;
              case AB_BA:
                if(verbose)
                  fprintf(GlobalData->stderrc,"* 4 Wrong orientation AB_BA\n");
                continue;  // wrong orientation
              default:
                assert(0);
                break;
            }
            break;
          default:
            assert(0);
            break;
        }
        if(verbose)
          fprintf(GlobalData->stderrc,"* foundBC = TRUE!\n");
        foundBC = TRUE;
        if(verbose){
          fprintf(GlobalData->stderrc,"* Found all three edges (" F_CID "," F_CID ",%c)%d, (" F_CID "," F_CID ",%c)%d AND (" F_CID "," F_CID ",%c)%d\n",
                  scaffoldA->id, scaffoldB->id, GetEdgeOrientationWRT(curEdge, scaffoldA->id), curEdge->edgesContributing,
                  scaffoldA->id, otherID, GetEdgeOrientationWRT(edge, scaffoldA->id), edge->edgesContributing,
                  otherBID,scaffoldB->id, GetEdgeOrientationWRT(edgeB, otherBID), edgeB->edgesContributing);
        }
        if(edge->edgesContributing > edgeB->edgesContributing){
          edge->quality = edge->edgesContributing + edgeB->edgesContributing;
          edgeB->quality = -1.0;
          // If we have a singleton from each side, we accept this
          if(edge->quality +0.1 >= (float)CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD){
            
	    //AppendPtrT(mergedEdges, (void *)&edge);
            AppendPtrT(mergedEdges, (const void *) &edge);
          }
          
        }else{
          edgeB->quality = edge->edgesContributing + edgeB->edgesContributing;
          edge->quality = -1.0;
          if(edge->quality + 0.1 >= (float)CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD){
            
	    //AppendPtrT(mergedEdges, (void *)&edge);
            AppendPtrT(mergedEdges, (const void *) &edge);
          }
          
	  //AppendPtrT(mergedEdges, (void *)&edgeB);
          AppendPtrT(mergedEdges, (const void *)&edgeB);
        }
      }else{
        if(verbose)
          fprintf(GlobalData->stderrc,"* otherBID(" F_CID ") != otherID(" F_CID ")\n", otherBID, otherID);
      }
      
    }
    if(foundBC == FALSE){
      if(verbose)
        fprintf(GlobalData->stderrc,"* Found %sconfirmed AC edge only\n", (edge->edgesContributing < CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD?"UN":""));
      if(edge->edgesContributing < CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD)
        continue;
      
      //AppendPtrT(mergedEdges, (void *)&edge);
      AppendPtrT(mergedEdges, (const void *) &edge);
      
      edge->quality = edge->edgesContributing;
      edge->flags.bits.isProbablyBogus = TRUE; // overloaded for one-sided
      
    }
  }
  
  // Repeat, iterating over edgesFromB
  // If an edge has already been weighted (weight != 0), skip it
  
  if(verbose)
    fprintf(GlobalData->stderrc,"* Iterating over %d edges from B\n",
            (int) GetNumPtrTs(edgesFromB));
  for(i = 0; i < GetNumPtrTs(edgesFromB); i++){
    CDS_CID_t otherID;
    edge = (SEdgeT *) *GetPtrT(edgesFromB,i);
    // We already dealt with this from the other direction
    if(edge->quality != 0.0)
      continue;
    if(!doInterleaving && isLargeOverlap(edge))
      continue;
    if(edge->idA == scaffoldB->id){
      otherID = edge->idB;
    }else{
      otherID = edge->idA;
    }
    if(edge->edgesContributing < CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD)
      continue;
    if(verbose)
      fprintf(GlobalData->stderrc,"* Found a BC edge (" F_CID "," F_CID ",%c)\n",
              scaffoldB->id, otherID, GetEdgeOrientationWRT(edge, scaffoldB->id));
    
    
    //AppendPtrT(mergedEdges, (void *)&edge);
    AppendPtrT(mergedEdges, (const void *) &edge);
    
    edge->quality = edge->edgesContributing;
    edge->flags.bits.isProbablyBogus = TRUE; // overloaded for one-sided
  }
  
  
  return GetNumPtrTs(mergedEdges);
}


/********************************************************************/
SEdgeT *FindMoreAttractiveMergeEdge(SEdgeT *curEdge,
                                    CIScaffoldT *scaffoldA,
                                    CIScaffoldT *scaffoldB,
                                    int32 recDepth,
                                    int doInterleaving,
                                    int checkForTinyScaffolds,
                                    int32 verbose){
  VA_TYPE(PtrT) *edgesFromA = CreateVA_PtrT(10);
  VA_TYPE(PtrT) *edgesFromB = CreateVA_PtrT(10);
  VA_TYPE(PtrT) *mergedEdges = CreateVA_PtrT(10);
  int32 endA, endB;
  int32 markedA, markedB;
  
  assert(scaffoldA->id == curEdge->idA || scaffoldA->id == curEdge->idB);
  assert(scaffoldB->id == curEdge->idA || scaffoldB->id == curEdge->idB);
  
  if(verbose)
    fprintf(GlobalData->stderrc,"* FMA (" F_CID "," F_CID ") depth:%d\n", scaffoldA->id, scaffoldB->id, recDepth);
  
  if(scaffoldA->flags.bits.smoothSeenAlready &&
     scaffoldB->flags.bits.smoothSeenAlready){
    
    if(verbose)
      fprintf(GlobalData->stderrc,"* We've already seen both scaffolds...returning\n");
    return NULL;
  }
  // Cut off recursion back to these nodes
  scaffoldA->flags.bits.smoothSeenAlready = 1;
  scaffoldB->flags.bits.smoothSeenAlready = 1;
  curEdge->flags.bits.isBogus = TRUE; // We mark all edges visited by the recursion so we don't retrace our steps
  
  {
    ChunkOrient orientEndNodeA, orientEndNodeB;
    ChunkOrientationType edgeEndsOrient;
    NodeCGW_T *endNodeA, *endNodeB;
    double aGapSize, bGapSize;
    
    if(verbose)
      fprintf(GlobalData->stderrc,"* FindMoreAttractiveMergeEdge %d\n", recDepth);
    TranslateScaffoldOverlapToContigOverlap(scaffoldA, scaffoldB, curEdge->orient, &endNodeA, &endNodeB, &orientEndNodeA, &orientEndNodeB, &edgeEndsOrient,
                                            &aGapSize, &bGapSize);
  }
  
  
  switch(curEdge->orient){
    case AB_AB:
      endA = B_END;
      endB = A_END;
      break;
    case AB_BA:
      endA = B_END;
      endB = B_END;
      break;
    case BA_AB:
      endA = A_END;
      endB = A_END;
      break;
    case BA_BA:
      endA = A_END;
      endB = B_END;
      break;
    default:
      assert(0);
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
  
  if(markedA || markedB){
    if(verbose)
      fprintf(GlobalData->stderrc,"* Marked Edges Encountered...returning\n");
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
  if(GetNumPtrTs(edgesFromA) == 0 &&
     GetNumPtrTs(edgesFromB) == 0){
    if(verbose){
      fprintf(GlobalData->stderrc,"* No edges from from A OR B!\n");
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
  
  
  
  
  if(AssignEdgeWeights(mergedEdges, curEdge, scaffoldA,edgesFromA, scaffoldB, edgesFromB, doInterleaving, verbose))
    {
      int i;
      //for all edges curEdge in list
      for(i = 0; i < GetNumPtrTs(mergedEdges); i++){
        // fEdge may be an edge from EITHER scaffoldA OR ScaffoldB to
        // a scaffoldC that is in between A and B
        SEdgeT *fEdge = (SEdgeT *)*GetPtrT(mergedEdges,i);
        CDS_CID_t sourceID, targetID;
        int32 overlapFound;
        int32 isFromA;
        CIScaffoldT *targetScaffold, *sourceScaffold;
      
        if(verbose){
          PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph, "AEW", fEdge, fEdge->idA);
          fprintf(GlobalData->stderrc," Edge has quality %g\n", fEdge->quality);
        }
        if(fEdge->idA == scaffoldA->id || fEdge->idA == scaffoldB->id){
          targetID = fEdge->idB;
          sourceID = fEdge->idA;
        }else{
          targetID = fEdge->idA;
          sourceID = fEdge->idB;
        }
      
        if(sourceID == scaffoldA->id){
          isFromA = TRUE;
        }else{
          isFromA = FALSE;
        }
      
        targetScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,targetID);
        sourceScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,sourceID);
      
        if(verbose)
          PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph, "fedge ", fEdge, fEdge->idA);
      
        if(DoesScaffoldCFit(scaffoldA, scaffoldB, targetScaffold,curEdge, fEdge, edgesFromA, edgesFromB, &overlapFound, checkForTinyScaffolds, verbose)){
          SEdgeT *gEdge;
        
          if(verbose){
            fprintf(GlobalData->stderrc,"*### Scaffold " F_CID " fits between (" F_CID "," F_CID ")  curEdge = (" F_CID "," F_CID ",%c)   fEdge = (" F_CID "," F_CID ",%c)\n",
                    targetID, scaffoldA->id, scaffoldB->id, 
                    curEdge->idA, curEdge->idB, curEdge->orient,
                    fEdge->idA, fEdge->idB, fEdge->orient);
          
            fprintf(GlobalData->stderrc,"* Calling FindMoreAttractive recursively on " F_CID " " F_CID "\n", sourceScaffold->id, targetScaffold->id);
          }
        
          gEdge = FindMoreAttractiveMergeEdge(fEdge,sourceScaffold, targetScaffold, recDepth + 1, doInterleaving, checkForTinyScaffolds, verbose);
        
        
          if(gEdge != NULL){
            if(verbose){
              PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph, "FMA returned gedge ", gEdge,gEdge->idA);
            }
            MarkScaffoldsForMerging(fEdge, FALSE);
            MarkScaffoldForNotMerging(targetScaffold);
          }else{
            if(verbose)
              fprintf(GlobalData->stderrc,"*gedge = NULL returning fEdge\n");
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
        }else{
          if(verbose)
            fprintf(GlobalData->stderrc,"*XXX Scaffold " F_CID " DOES NOT fit between (" F_CID "," F_CID ")  curEdge = (" F_CID "," F_CID ",%c)   fEdge = (" F_CID "," F_CID ",%c)\n",
                    targetID, scaffoldA->id, scaffoldB->id, 
                    curEdge->idA, curEdge->idB, curEdge->orient,
                    fEdge->idA, fEdge->idB, fEdge->orient);
        
          //	MarkScaffoldForNotMerging(targetScaffold);
        }
      
      }
    
    }else{
      if(verbose){
        fprintf(GlobalData->stderrc,"No edges found\n");
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


#define MAX_SLOP_IN_STD 3.5

int isQualityScaffoldMergingEdge(SEdgeT * curEdge,
                                 CIScaffoldT * scaffoldA,
                                 CIScaffoldT * scaffoldB,
                                 ScaffoldInstrumenter * si,
                                 VA_TYPE(PtrT) * MIs,
                                 float minSatisfied,
                                 float maxDelta)
{
  fprintf(stderr,"***** MERGE CANDIDATE: Trying to merge scaffold " F_CID " (%.1f bp) with scaffold " F_CID " (%.1f bp) based on gap of %.1f from %d weight edge\n",
	  scaffoldA->id,
	  scaffoldA->bpLength.mean,
	  scaffoldB->id,
	  scaffoldB->bpLength.mean,
	  curEdge->distance.mean,
	  curEdge->edgesContributing);
	  
  if(si != NULL && (minSatisfied > 0.0 || maxDelta > 0.0))
    {
      MateInstrumenter matesAfter;
      MateInstrumenter * sABefore;
      MateInstrumenter * sBBefore;
      float fractMatesHappyAfter;
      MateInstrumenter matesBefore;
      float fractMatesHappyBefore;
    
      if(GetNumVA_PtrT(MIs) > scaffoldA->id &&
         *GetVA_PtrT(MIs, scaffoldA->id) != NULL)
        {
          sABefore = (MateInstrumenter *) *GetVA_PtrT(MIs, scaffoldA->id);
        }
      else
        {
          sABefore = (MateInstrumenter *) safe_calloc(1, sizeof(MateInstrumenter));
          InstrumentScaffold(ScaffoldGraph,
                             scaffoldA,
                             si,
                             InstrumenterVerbose2,
                             GlobalData->stderrc);
          GetMateInstrumenterFromScaffoldInstrumenter(sABefore, si);
          SetVA_PtrT(MIs, scaffoldA->id, (void *) &sABefore);
        }
      ResetMateInstrumenterCounts(&matesBefore);
      AddMateInstrumenterCounts(&matesBefore, sABefore);
    
    
      if(GetNumVA_PtrT(MIs) > scaffoldB->id &&
         *GetVA_PtrT(MIs, scaffoldB->id) != NULL)
        {
          sBBefore = (MateInstrumenter *) *GetVA_PtrT(MIs, scaffoldB->id);
        }
      else
        {
          sBBefore = (MateInstrumenter *) safe_calloc(1, sizeof(MateInstrumenter));
          InstrumentScaffold(ScaffoldGraph,
                             scaffoldB,
                             si,
                             InstrumenterVerbose2,
                             GlobalData->stderrc);
          GetMateInstrumenterFromScaffoldInstrumenter(sBBefore, si);
          SetVA_PtrT(MIs, scaffoldB->id, (void *) &sBBefore);
        }
      AddMateInstrumenterCounts(&matesBefore, sBBefore);
    
      InstrumentScaffoldPair(ScaffoldGraph,
                             curEdge,
                             si,
                             InstrumenterVerbose2,
                             GlobalData->stderrc);
      GetMateInstrumenterFromScaffoldInstrumenter(&matesAfter,
                                                  si);
    
      fractMatesHappyAfter =
        ((float) (GetMateStatsHappy(&(matesAfter.intra)) +
                  GetMateStatsHappy(&(matesAfter.inter)))) /
        (GetMateStatsBad(&(matesAfter.intra)) +
         GetMateStatsBad(&(matesAfter.inter)) +
         GetMateStatsHappy(&(matesAfter.intra)) +
         GetMateStatsHappy(&(matesAfter.inter)));
    

      if(GetMateStatsBad(&(matesBefore.intra)) +
         GetMateStatsBad(&(matesBefore.inter)) +
         GetMateStatsHappy(&(matesBefore.intra)) +
         GetMateStatsHappy(&(matesBefore.inter)) > 0)
        {
          fractMatesHappyBefore =
            ((float) (GetMateStatsHappy(&(matesBefore.intra)) +
                      GetMateStatsHappy(&(matesBefore.inter)))) /
            (GetMateStatsBad(&(matesBefore.intra)) +
             GetMateStatsBad(&(matesBefore.inter)) +
             GetMateStatsHappy(&(matesBefore.intra)) +
             GetMateStatsHappy(&(matesBefore.inter)));
          fprintf(GlobalData->stderrc,
                  "* %.3f mates satisfied if scaffolds " F_CID "," F_CID " instrumented separately\n",
                  fractMatesHappyBefore, curEdge->idA, curEdge->idB);
        }
      else
        {
          fractMatesHappyBefore = 1.0;
          fprintf(GlobalData->stderrc,
                  "* %.3f mates satisfied if scaffolds " F_CID "," F_CID " instrumented separately (no mates)\n",
                  fractMatesHappyBefore, curEdge->idA, curEdge->idB);
        }
    
      fprintf(GlobalData->stderrc,
              "* %.3f mates satisfied if scaffolds " F_CID "," F_CID " merged via (%.2f,%.2f,%s)\n",
              fractMatesHappyAfter, curEdge->idA, curEdge->idB,
              curEdge->distance.mean, curEdge->distance.variance,
              ((curEdge->orient == AB_AB) ? "AB_AB" :
               ((curEdge->orient == AB_BA) ? "AB_BA" :
                ((curEdge->orient == BA_AB) ? "BA_AB" : "BA_BA"))));

#if 0    
      fprintf(GlobalData->stderrc,
	      "******** after bad mates: %d\n",
	      GetMateStatsBad(&(matesAfter.inter))+GetMateStatsBad(&(matesAfter.intra)));
      fprintf(GlobalData->stderrc,
	      "******** after good mates: %d\n",
	      GetMateStatsHappy(&(matesAfter.inter))+GetMateStatsHappy(&(matesAfter.intra)));
      fprintf(GlobalData->stderrc,
	      "******** before bad mates: %d\n",
	      GetMateStatsBad(&(matesBefore.inter))+GetMateStatsBad(&(matesBefore.intra)));
      fprintf(GlobalData->stderrc,
	      "******** before good mates: %d\n",
	      GetMateStatsHappy(&(matesBefore.inter))+GetMateStatsHappy(&(matesBefore.intra)));
#endif

      if(fractMatesHappyAfter < minSatisfied &&
         fractMatesHappyAfter < fractMatesHappyBefore
#ifdef REQUIRE_MORE_THAN_ONE_BAD_TO_REJECT
         &&  GetMateStatsBad(&(matesAfter.inter)) +  GetMateStatsBad(&(matesAfter.intra))  - 
         (GetMateStatsBad(&(matesBefore.inter)) +  GetMateStatsBad(&(matesBefore.intra)))
         >1
#endif
#ifdef REQUIRE_BAD_APPROACHING_HAPPY_TO_REJECT
         &&( (          (float) (GetMateStatsHappy(&(matesAfter.inter)) +  GetMateStatsHappy(&(matesAfter.intra))  - 
                                 GetMateStatsHappy(&(matesBefore.inter)) +  GetMateStatsHappy(&(matesBefore.intra)))
                        <= 0.) ? TRUE : 
             ((float) (GetMateStatsBad(&(matesAfter.inter)) +  GetMateStatsBad(&(matesAfter.intra))  - 
                       GetMateStatsBad(&(matesBefore.inter)) +  GetMateStatsBad(&(matesBefore.intra))) 
              /
              (float) (GetMateStatsHappy(&(matesAfter.inter)) +  GetMateStatsHappy(&(matesAfter.intra))  - 
                       GetMateStatsHappy(&(matesBefore.inter)) +  GetMateStatsHappy(&(matesBefore.intra))) 
              > MAX_FRAC_BAD_TO_GOOD)
             )
#endif
         )
        {
          fprintf(GlobalData->stderrc,
                  "***** Merging would result in too low a satisfied mate fraction (%.3f < %.3f) - shouldn't merge\n",
                  fractMatesHappyAfter, minSatisfied);
          fprintf(GlobalData->stderrc,
                  "******** inter bad mates: %d\n",GetMateStatsBad(&(matesAfter.inter)));
          fprintf(GlobalData->stderrc,
                  "******** inter good mates: %d\n",GetMateStatsHappy(&(matesAfter.inter)));
#ifdef DEBUG_BAD_MATE_RATIO
          DumpACIScaffoldNew(stderr,ScaffoldGraph,scaffoldA,FALSE);
          DumpACIScaffoldNew(stderr,ScaffoldGraph,scaffoldB,FALSE);
#endif
#ifdef DRAW_BAD_MATE_CAMS
          if(curEdge->edgesContributing>=DRAW_BAD_MATE_CAMS)
            {
              int64 scaffoldAEndCoord = 0, scaffoldBEndCoord = 0, endcoord=0;
              char camname[1000];
              static int lowmate_count=0;
              FILE *camfile=NULL;

	
              if(++lowmate_count < 500){
                AS_UTL_mkdir("MergeCams");
                sprintf(camname,"MergeCams/lowmate_failure_%d_%d.cam",scaffoldA->id,scaffoldB->id);
                camfile = fopen(camname,"w");
                assert(camfile!=NULL);
                DumpCelamyColors(camfile);
                DumpCelamyMateColors(camfile);

                do_draw_frags_in_CelamyScaffold=1;
                DumpCelamyFragColors(camfile);
                //	if(endcoord!=0)endcoord+=1000000;
                if(scaffoldA->bpLength.mean < -curEdge->distance.mean){
                  endcoord+=-curEdge->distance.mean-scaffoldA->bpLength.mean;
                }
                if(curEdge->orient == AB_AB || curEdge->orient == AB_BA){
                  scaffoldAEndCoord = endcoord;
                  scaffoldBEndCoord = endcoord + scaffoldA->bpLength.mean;
                } else {
                  scaffoldBEndCoord = endcoord;
                  scaffoldAEndCoord = endcoord + scaffoldA->bpLength.mean;
                }
                CelamyScaffold(camfile,scaffoldA,scaffoldAEndCoord,scaffoldBEndCoord);
                endcoord += scaffoldA->bpLength.mean;
                endcoord += curEdge->distance.mean;
                if(curEdge->orient == AB_AB || curEdge->orient == BA_AB){
                  scaffoldAEndCoord = endcoord;
                  scaffoldBEndCoord = endcoord + scaffoldB->bpLength.mean;
                } else {
                  scaffoldBEndCoord = endcoord;
                  scaffoldAEndCoord = endcoord + scaffoldB->bpLength.mean;
                }
                CelamyScaffold(camfile,scaffoldB,scaffoldAEndCoord,scaffoldBEndCoord);
                do_draw_frags_in_CelamyScaffold=0;
                //	endcoord += scaffoldB->bpLength.mean;
                PrintScaffoldInstrumenterMateDetails(si,camfile,PRINTCELAMY);
                PrintExternalMateDetailsAndDists(ScaffoldGraph,si->bookkeeping.wExtMates,"\t",camfile,PRINTCELAMY);
                PrintUnmatedDetails(si,camfile,PRINTCELAMY);


                fclose(camfile);
              }
            }
#endif
          return FALSE;
        }
    
      if(maxDelta > 0.0 &&
         fractMatesHappyBefore - fractMatesHappyAfter > maxDelta)
        /*
          if(CompareMateInstrumenters(&matesBefore, &matesAfter,
          InstrumenterVerbose2,
          GlobalData->stderrc) == InstrumenterWorse)
        */
        {
          fprintf(GlobalData->stderrc,
                  "***** Merging would decrease satisfied mate fraction by too much (%.3f > %.3f) - shouldn't merge\n",
                  fractMatesHappyBefore - fractMatesHappyAfter, maxDelta);
          return FALSE;
        }
    }
  return TRUE;
}

void SaveBadScaffoldMergeEdge(SEdgeT * edge,
                              ChunkOverlapperT * overlapper)
{
  ChunkOverlapCheckT * lookup;
  ChunkOverlapSpecT spec;
  double delta = sqrt(edge->distance.variance) * 3.;
  
  InitCanonicalOverlapSpec(edge->idA, edge->idB, edge->orient, &spec);
  if((lookup = LookupCanonicalOverlap(overlapper, &spec)) == NULL)
    {
      ChunkOverlapCheckT olap = {0};
      // create
      FillChunkOverlapWithEdge(edge, &olap);
      olap.minOverlap = -edge->distance.mean - delta;
      olap.maxOverlap = -edge->distance.mean + delta;
      InsertChunkOverlap(overlapper, &olap);
    }
  else
    {
      // update
      lookup->minOverlap = MIN(lookup->minOverlap, -edge->distance.mean - delta);
      lookup->maxOverlap = MAX(lookup->maxOverlap, -edge->distance.mean + delta);
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

#ifdef  AGGRESSIVE_ABUTTING
#define MAX_OVERLAP_TO_ABUT                          2000
#define MAX_PERC_SCAFFOLD_LEN                        0.5
#endif

int LooseAbuttingCheck(SEdgeT * curEdge,
                       CIScaffoldT * scaffoldA,
                       CIScaffoldT * scaffoldB)
{
  double stddevsPerWeight;
  double edgeMinScaffoldLengthRatio;

  stddevsPerWeight = -(curEdge->distance.mean + CGW_MISSED_OVERLAP) /
    (curEdge->edgesContributing * sqrt(MAX(0.01,curEdge->distance.variance)));

  edgeMinScaffoldLengthRatio = 
    -curEdge->distance.mean /
    MIN(scaffoldA->bpLength.mean, scaffoldB->bpLength.mean);
  
#ifdef DEBUG1
  fprintf(GlobalData->stderrc, "Performing loose abutting check:\n");
  PrintSEdgeT(GlobalData->stderrc, ScaffoldGraph,
              "sEdge", curEdge, curEdge->idA);
  fprintf(GlobalData->stderrc,
          "Scaffold A: " F_CID ", length: %f\tScaffold B: " F_CID ", length: %f\n",
          scaffoldA->id, scaffoldA->bpLength.mean,
          scaffoldB->id, scaffoldB->bpLength.mean);
#endif

  if ((stddevsPerWeight < STDDEVS_PER_WEIGHT_THRESHOLD &&
       edgeMinScaffoldLengthRatio < EDGE_PER_MIN_SCAFFOLD_LENGTH_THRESHOLD)
#ifdef  AGGRESSIVE_ABUTTING
      ||
      (edgeMinScaffoldLengthRatio < MAX_PERC_SCAFFOLD_LEN &&
       -curEdge->distance.mean < MAX_OVERLAP_TO_ABUT)
#endif
      )
    {
    
#ifdef DEBUG1
      fprintf(GlobalData->stderrc, "Loose abutting check passed.\n");
#endif
    
      return TRUE;
    }
  
#ifdef DEBUG1
  fprintf(GlobalData->stderrc, "Loose abutting check failed.\n");
#endif


#ifdef DRAW_BAD_MATE_CAMS
  if(curEdge->edgesContributing>=DRAW_BAD_MATE_CAMS
     &&curEdge->distance.mean > -5000
     &&scaffoldA->bpLength.mean > curEdge->distance.mean + 10000
     &&scaffoldB->bpLength.mean > curEdge->distance.mean + 10000
     )
    {

      int64 scaffoldAEndCoord = 0, scaffoldBEndCoord = 0, endcoord=0;
      char camname[1000];
      static int looseabut_count=0;
      FILE *camfile=NULL;

      if(++looseabut_count < 500){
        AS_UTL_mkdir("MergeCams");
        sprintf(camname,"MergeCams/looseabut_failure_%d_%d.cam",scaffoldA->id,scaffoldB->id);
        camfile = fopen(camname,"w");
        assert(camfile!=NULL);
        DumpCelamyColors(camfile);
        DumpCelamyMateColors(camfile);

        do_draw_frags_in_CelamyScaffold=1;
        DumpCelamyFragColors(camfile);
        //	if(endcoord!=0)endcoord+=1000000;
        if(scaffoldA->bpLength.mean < -curEdge->distance.mean){
          endcoord+=-curEdge->distance.mean-scaffoldA->bpLength.mean;
        }
        if(curEdge->orient == AB_AB || curEdge->orient == AB_BA){
          scaffoldAEndCoord = endcoord;
          scaffoldBEndCoord = endcoord + scaffoldA->bpLength.mean;
        } else {
          scaffoldBEndCoord = endcoord;
          scaffoldAEndCoord = endcoord + scaffoldA->bpLength.mean;
        }
        CelamyScaffold(camfile,scaffoldA,scaffoldAEndCoord,scaffoldBEndCoord);
        endcoord += scaffoldA->bpLength.mean;
        endcoord += curEdge->distance.mean;
        if(curEdge->orient == AB_AB || curEdge->orient == BA_AB){
          scaffoldAEndCoord = endcoord;
          scaffoldBEndCoord = endcoord + scaffoldB->bpLength.mean;
        } else {
          scaffoldBEndCoord = endcoord;
          scaffoldAEndCoord = endcoord + scaffoldB->bpLength.mean;
        }
        CelamyScaffold(camfile,scaffoldB,scaffoldAEndCoord,scaffoldBEndCoord);
        do_draw_frags_in_CelamyScaffold=0;
        //	endcoord += scaffoldB->bpLength.mean;

        {
          static ScaffoldInstrumenter *si=NULL;
          if(si==NULL){
            si = CreateScaffoldInstrumenter(ScaffoldGraph,INST_OPT_ALL);
            assert(si!=NULL);
          }
          InstrumentScaffoldPair(ScaffoldGraph,
                                 curEdge,
                                 si,
                                 InstrumenterVerbose2,
                                 GlobalData->stderrc);

          PrintScaffoldInstrumenterMateDetails(si,camfile,PRINTCELAMY);
          PrintExternalMateDetailsAndDists(ScaffoldGraph,si->bookkeeping.wExtMates,"\t",camfile,PRINTCELAMY);
          PrintUnmatedDetails(si,camfile,PRINTCELAMY);

        }

        fclose(camfile);
      }
    }
#endif

  
  return FALSE; 
}

/*
  Prior to marking scaffolds for merging, mean may be changed but not variance
  So, see if a change of mean to -CGW_MISSED_OVERLAP will allow abutting
  after which edge can be marked as trusted for recomputing offsets
*/
int AbuttingWillWork(SEdgeT * curEdge,
                     CIScaffoldT * scaffoldA,
                     CIScaffoldT * scaffoldB,
                     InterleavingSpec * iSpec)
{
  int isAbuttable = FALSE;
  float chiSquaredValue;
  
  if(curEdge->distance.mean >= -CGW_MISSED_OVERLAP)
    return TRUE;

  // the edge must be compatible with a -20 edge
  // or be strong & short relative to the scaffold lengths
  if(PairwiseChiSquare(-CGW_MISSED_OVERLAP,
                       curEdge->distance.variance,
                       curEdge->distance.mean,
                       curEdge->distance.variance,
                       (LengthT *) NULL,
                       &chiSquaredValue,
                       PAIRWISECHI2THRESHOLD_CGW) ||
     LooseAbuttingCheck(curEdge, scaffoldA, scaffoldB))
    {
      double originalMean = curEdge->distance.mean;
    
      curEdge->distance.mean = -CGW_MISSED_OVERLAP;

      // a -20 edge must have high satisfied mate %
      if(isQualityScaffoldMergingEdge(curEdge,
                                      scaffoldA,
                                      scaffoldB,
                                      iSpec->sai->scaffInst,
                                      iSpec->MIs,
                                      iSpec->minSatisfied,
                                      iSpec->maxDelta))
        {
          isAbuttable = TRUE;
        }
      curEdge->distance.mean = originalMean;
    }

  if(!isAbuttable)
    SaveBadScaffoldMergeEdge(curEdge, iSpec->badSEdges);
  
  return isAbuttable;
}


int isBadScaffoldMergeEdge(SEdgeT * edge, ChunkOverlapperT * overlapper)
{
  ChunkOverlapCheckT * lookup;
  ChunkOverlapSpecT spec;
  double delta = sqrt(edge->distance.variance) * 3.;
  CDS_COORD_t minOverlap = (CDS_COORD_t) -edge->distance.mean - delta;
  CDS_COORD_t maxOverlap = (CDS_COORD_t) -edge->distance.mean + delta;
  
  InitCanonicalOverlapSpec(edge->idA, edge->idB, edge->orient, &spec);
  if((lookup = LookupCanonicalOverlap(overlapper, &spec)) != NULL)
    {
      if(minOverlap >= lookup->minOverlap - 10 &&
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



void ExamineSEdgeForUsability(VA_TYPE(PtrT) * sEdges,
                              SEdgeT * curEdge, InterleavingSpec * iSpec,
                              int verbose)
{
  SEdgeT *mergeEdge = NULL;
  CIScaffoldT *scaffoldA, *scaffoldB;
  double mergeDistance,maxMergeDistance, minMergeDistance;
  int mayOverlap, mustOverlap;
  double length_to_dist;
  double min_scaffold_length;
  
  if(CONFIRMED_SCAFFOLD_EDGE_THRESHHOLD > curEdge->edgesContributing)
    return;
  
  if(GlobalData->doInterleavedScaffoldMerging == FALSE && isLargeOverlap(curEdge)){   /// too big an overlap
    MarkScaffoldsForMerging(curEdge, FALSE  /* mark for simply exclusion from other merges */);
    return;
  }
  
  scaffoldA = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idA);
  scaffoldB = GetGraphNode(ScaffoldGraph->ScaffoldGraph, curEdge->idB);
  if(verbose)
    PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph,
                   "S ", curEdge, curEdge->idA);
  
  // We don't want to stick a teeny tiny element in the middle of a gap between two
  // giants.  Often this will later prevent the giants from merging.  So, we
  // prevent it.  We only place stuff such that the size of the element being placed
  // is at least 20% of the size of the implied gap.  So, for a 50k gap, we need a 10k
  // scaffold.  For a 10k gap, we need a 2k scaffold, etc.
  length_to_dist = 1.0;
  min_scaffold_length = MIN(scaffoldA->bpLength.mean,
                            scaffoldB->bpLength.mean);
  if(curEdge->distance.mean > 0){
    length_to_dist = ( min_scaffold_length / curEdge->distance.mean);
  }
  
  if(iSpec->checkForTinyScaffolds &&
     min_scaffold_length < 5000 &&
     length_to_dist < 0.20){
    if(verbose)
      fprintf(GlobalData->stderrc,
              "Scaffolds are too short (%g,%g) relative to edge length (%g)\n",
              scaffoldA->bpLength.mean, 
              scaffoldB->bpLength.mean, 
              curEdge->distance.mean);
    return;
  }
  
  if(scaffoldA->flags.bits.walkedAlready ||
     scaffoldB->flags.bits.walkedAlready ){
    //  We want to make sure that a week link doesn't preempt a strong link
    scaffoldA->flags.bits.walkedAlready =
      scaffoldB->flags.bits.walkedAlready = 1;
    //      fprintf(GlobalData->stderrc,"* Skipping edge due to marking!\n");
    return;
  }
  
  /******	******	******	******	******	******	******	******	
   ******	New Recursive Code      ******	******	******	******	
   ******	******	******	******	******	******	******	******	*/
  if(TouchesMarkedScaffolds(curEdge)){
    if(verbose)
      fprintf(GlobalData->stderrc,
              "* Edge (" F_CID "," F_CID ",%c) touches marked scaffolds\n",
              curEdge->idA, curEdge->idB, curEdge->orient);
    return;
  }
  
  
  // See if we already want to merge these two scaffolds, but in an opposite orientation
  // For instance, sometimes there will be evidence suggesting a merge AB_BA and a merge
  // BA_BA between a given pair of scaffolds.
  if((scaffoldA->essentialEdgeB != NULLINDEX &&
      (scaffoldA->essentialEdgeB == scaffoldB->essentialEdgeA ||
       scaffoldA->essentialEdgeB == scaffoldB->essentialEdgeB )) ||
     (scaffoldA->essentialEdgeA != NULLINDEX &&
      (scaffoldA->essentialEdgeA == scaffoldB->essentialEdgeA ||
       scaffoldA->essentialEdgeA == scaffoldB->essentialEdgeB ))){
    
    fprintf(GlobalData->stderrc,
            "* We're already trying to merge scaffold (" F_CID "," F_CID ") ...back off!\n",
            scaffoldA->id, scaffoldB->id);
    return;
  }
  
  
  if(verbose){
    fprintf(GlobalData->stderrc,
            "* Top Level call to FindMoreAttractiveMergeEdge (" F_CID "," F_CID ",%c) (" F_CID "," F_CID ") gap = %g\n",
            curEdge->idA, curEdge->idB, curEdge->orient,
            scaffoldA->id, scaffoldB->id, curEdge->distance.mean);
  }
  
  mergeDistance = curEdge->distance.mean;
  minMergeDistance =
    mergeDistance - MAX_SLOP_IN_STD * sqrt(curEdge->distance.variance);
  maxMergeDistance =
    mergeDistance + MAX_SLOP_IN_STD * sqrt(curEdge->distance.variance);
  
  if(verbose){
    fprintf(GlobalData->stderrc,
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
  if(mustOverlap && curEdge->edgesContributing > EDGE_QUANTA){
    mustOverlap = 0;
    mayOverlap = 1;
  }
  
  
  if (GlobalData->doInterleavedScaffoldMerging == FALSE)
    {
      if(mustOverlap){
        // All edges that got this far passed the !isLargeOverlap test, so we will
        // abut them if we can't find the overlap
        mergeEdge = curEdge;
      
      }else{
      
        mergeEdge = FindMoreAttractiveMergeEdge(curEdge, scaffoldA, scaffoldB, 0, GlobalData->doInterleavedScaffoldMerging, iSpec->checkForTinyScaffolds, verbose);
        if(verbose){
          if(!mergeEdge){
            fprintf(GlobalData->stderrc,"*(NONE) Return from Top Level call to FindMoreAttractiveMergeEdge (" F_CID "," F_CID ",%c) (" F_CID "," F_CID ") ==> NONE!\n",
                    curEdge->idA, curEdge->idB, curEdge->orient, scaffoldA->id, scaffoldB->id);
          }else{
            fprintf(GlobalData->stderrc,"*(%s) Return from Top Level call to FindMoreAttractiveMergeEdge (" F_CID "," F_CID ",%c) (" F_CID "," F_CID ") ==> merge (" F_CID "," F_CID ",%c)\n",
                    (curEdge == mergeEdge?"N":"Y"),
                    curEdge->idA, curEdge->idB, curEdge->orient, scaffoldA->id, scaffoldB->id,
                    mergeEdge->idA, mergeEdge->idB, mergeEdge->orient);
          }
        }
      
        // If we backed off on the merge, continue
        if(mergeEdge == NULL){
          // Should we propagate the marks to prevent anything from happening in this neighborhood?
          MarkScaffoldsForMerging(curEdge,  FALSE  /* mark for simply exclusion from other merges */);
          return;
        }
      
        mergeDistance = mergeEdge->distance.mean;
        minMergeDistance = mergeDistance - MAX_SLOP_IN_STD * sqrt(mergeEdge->distance.variance);
        maxMergeDistance = mergeDistance + MAX_SLOP_IN_STD * sqrt(mergeEdge->distance.variance);
      
        if(verbose){
          fprintf(GlobalData->stderrc,"* mergeDistance = (%g,%g) min:%g max:%g\n",
                  mergeDistance, mergeEdge->distance.variance,
                  minMergeDistance, maxMergeDistance);
        }
        // Look for an overlap
      
        mayOverlap = (minMergeDistance < CGW_MISSED_OVERLAP &&    maxMergeDistance > CGW_MISSED_OVERLAP);
        mustOverlap =(minMergeDistance < CGW_MISSED_OVERLAP &&    maxMergeDistance < CGW_MISSED_OVERLAP);
        //      fprintf(GlobalData->stderrc,"* mayOverlap = %d mustOverlap %d\n", mayOverlap, mustOverlap);
      
      
        if(mustOverlap && mergeEdge->edgesContributing > EDGE_QUANTA){
          mustOverlap = 0;
          mayOverlap = 1;
        }
      
        // Reset scaffoldA or scaffoldB so scaffoldA is to the left of scaffoldB
        if(scaffoldA->id == mergeEdge->idA || scaffoldA->id == mergeEdge->idB)
          {
            scaffoldB = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                     (scaffoldA->id == mergeEdge->idA) ?
                                     mergeEdge->idB : mergeEdge->idA);
          }
        else
          {
            scaffoldA = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                     (scaffoldB->id == mergeEdge->idA) ?
                                     mergeEdge->idB : mergeEdge->idA);
          }
      }
    
      if(isBadScaffoldMergeEdge(curEdge, iSpec->badSEdges))
        return;
    
      if(!isQualityScaffoldMergingEdge(curEdge,
                                       scaffoldA,
                                       scaffoldB,
                                       iSpec->sai->scaffInst,
                                       iSpec->MIs,
                                       iSpec->minSatisfied,
                                       iSpec->maxDelta))
        {
          SaveBadScaffoldMergeEdge(curEdge, iSpec->badSEdges);
          return;
        }
    
      if(mayOverlap || mustOverlap){ // may overlap
        ChunkOrient orientEndNodeA, orientEndNodeB;
        ChunkOrientationType edgeEndsOrient;
        NodeCGW_T *endNodeA, *endNodeB;
        double aGapSize, bGapSize;
        EdgeCGW_T *overlapEdge;
        float chiSquaredValue;
        int alternate;
      
        TranslateScaffoldOverlapToContigOverlap(scaffoldA, scaffoldB, mergeEdge->orient, &endNodeA, &endNodeB, &orientEndNodeA, &orientEndNodeB, &edgeEndsOrient,
                                                &aGapSize, &bGapSize);
      
        overlapEdge = FindOverlapEdgeChiSquare(ScaffoldGraph, endNodeA, endNodeB->id,
                                               edgeEndsOrient, mergeEdge->distance.mean,
                                               mergeEdge->distance.variance, &chiSquaredValue,
                                               (float)SCAFFOLD_MERGE_CHI2_THRESHHOLD, &alternate, verbose);
        if(overlapEdge == (EdgeCGW_T *)NULL){
          double minSizeScaffoldA = scaffoldA->bpLength.mean -
            (MAX_SLOP_IN_STD * sqrt(scaffoldA->bpLength.variance));
          double minSizeScaffoldB = scaffoldB->bpLength.mean -
            (MAX_SLOP_IN_STD * sqrt(scaffoldB->bpLength.variance));
          double maxGapScaffoldA = FindMaxGapInScaffold(ScaffoldGraph, scaffoldA);
          double maxGapScaffoldB = FindMaxGapInScaffold(ScaffoldGraph, scaffoldB);
        
          if(iSpec->checkAbutting &&
             !AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec))
            return;
        
          if(verbose){
            fprintf(GlobalData->stderrc, "Could not find overlap " F_CID "(" F_CID ") %c " F_CID "(" F_CID ")\n",
                    scaffoldA->id, endNodeA->id, (char)edgeEndsOrient,
                    scaffoldB->id, endNodeB->id);
          }
          if(mustOverlap ){
          
            if((minSizeScaffoldA > maxGapScaffoldB) &&
               (minSizeScaffoldB > maxGapScaffoldA)){
              SaveEdgeMeanForLater(mergeEdge);
              mergeEdge->distance.mean = MAX(- CGW_MISSED_OVERLAP, mergeEdge->distance.mean); 
              fprintf(GlobalData->stderrc,"* Must overlap edge -- can't confirm --abutting..going ahead distance is %g (%g,%g)\n", mergeEdge->distance.mean, minMergeDistance, maxMergeDistance);
              PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph, "  MustOverlap ", mergeEdge, mergeEdge->idA);
            }else{
              SaveEdgeMeanForLater(mergeEdge);
              mergeEdge->distance.mean = MAX(- CGW_MISSED_OVERLAP, mergeEdge->distance.mean); 
              fprintf(GlobalData->stderrc,"* Must overlap edge -- can't confirm, may be containment--abutting..going ahead distance is %g (%g,%g)\n", mergeEdge->distance.mean, minMergeDistance, maxMergeDistance);
              PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph, "  MustOverlap ", mergeEdge, mergeEdge->idA);
            }
          
          }else{
            // we don't know what to do, so 
            if(mergeEdge->edgesContributing > EDGE_QUANTA ||
               ((minSizeScaffoldA > maxGapScaffoldB) &&
                (minSizeScaffoldB > maxGapScaffoldA)) ){
              SaveEdgeMeanForLater(mergeEdge);
              mergeEdge->distance.mean = MAX(- CGW_MISSED_OVERLAP, mergeEdge->distance.mean); 
              if(verbose){
                fprintf(GlobalData->stderrc,"* May overlap edge -- can't confirm ...going ahead distance is %g (%g,%g)\n", mergeEdge->distance.mean, minMergeDistance, maxMergeDistance);
                PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph, "  MayOverlap ", mergeEdge, curEdge->idA);
              }
            }else{
              SaveEdgeMeanForLater(mergeEdge);
              mergeEdge->distance.mean = MAX(- CGW_MISSED_OVERLAP, mergeEdge->distance.mean); 
              if(verbose){
                fprintf(GlobalData->stderrc,"* May overlap edge -- can't confirm --abutting..going ahead distance is %g (%g,%g)\n", mergeEdge->distance.mean, minMergeDistance, maxMergeDistance);
                PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->ScaffoldGraph, "  MustOverlap ", mergeEdge, mergeEdge->idA);
              }
            }
          }
        }else if(alternate){
          // This means we were expecting the extremal contigs to overlap as follows:
          //    Expected                                      Found
          // -------A---------> 	                            -------A--------->
          //           -------B--------->           -------B--------->
          //           |-------|                    |--------------------------|
        
        
          fprintf(GlobalData->stderrc,"*** %s:Alternate found: NEED TO FIX THIS OVERLAP (" F_CID "," F_CID ",%c) is really (" F_CID "," F_CID ",%c)...setting to CGW_MISSED_OVERLAP\n",
                  __FILE__, endNodeA->id, endNodeB->id, edgeEndsOrient, overlapEdge->idA, overlapEdge->idB, overlapEdge->orient);
        
          // To deal with this case we need to handle all of the ramifications.  once A and B have slid by each other, B may overlap with
          // the contig to the left of A, and A may overlap with the contig to the right of B.  If either of these overlaps fails to exist,
          // we need to update the scaffold positions of these neighbors so that they do not cause a merge multialigns failure.
        
          // As a first cut, see if the extremal contig will fit in the gap of the other scaffold, so no positions need to be updated
          if(endNodeB->bpLength.mean + overlapEdge->distance.mean < aGapSize &&
             endNodeA->bpLength.mean + overlapEdge->distance.mean < bGapSize){
            fprintf(GlobalData->stderrc,"* We can safely interleave them endalength %g  endblength %g  aGap %g bGap %g\n",
                    endNodeA->bpLength.mean,
                    endNodeB->bpLength.mean,
                    aGapSize, bGapSize);
            SaveEdgeMeanForLater(mergeEdge);
            mergeEdge->distance.mean = -(endNodeA->bpLength.mean + endNodeB->bpLength.mean + overlapEdge->distance.mean);
          
            fprintf(GlobalData->stderrc,"*** Fixed so that overlap is length " F_CID " (%g) + length " F_CID " (%g) + overlap (%g) = %g\n",
                    endNodeA->id,endNodeA->bpLength.mean, endNodeB->id,endNodeB->bpLength.mean,overlapEdge->distance.mean, mergeEdge->distance.mean);
          }else{
            fprintf(GlobalData->stderrc,"*** Couldn't easily fix...gaps too small\n");
            if(iSpec->checkAbutting &&
               !AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec))
              return;
            SaveEdgeMeanForLater(mergeEdge);
            mergeEdge->distance.mean = -CGW_MISSED_OVERLAP;
            fprintf(GlobalData->stderrc,"* We CAN'T safely interleave them, gaps too small:  endalength %g  endblength %g  aGap %g bGap %g\n",
                    endNodeA->bpLength.mean,
                    endNodeB->bpLength.mean,
                    aGapSize, bGapSize);
          }
        }else{
          fprintf(GlobalData->stderrc,"* Confirmed overlap (" F_CID "," F_CID ",%c) overlapEdge:%g   mergeEdge:%g\n",
                  overlapEdge->idA, overlapEdge->idB, overlapEdge->orient, overlapEdge->distance.mean, mergeEdge->distance.mean);
          if(overlapEdge->orient != edgeEndsOrient && overlapEdge->idA == endNodeA->id)
            {
              assert(overlapEdge->orient == InvertEdgeOrient(edgeEndsOrient));
              overlapEdge->distance.mean =
                endNodeA->bpLength.mean + endNodeB->bpLength.mean +
                overlapEdge->distance.mean;
            }
          else
            mergeEdge->distance.mean = overlapEdge->distance.mean;
          if(abs(overlapEdge->distance.mean - mergeEdge->distance.mean) >
             MAX_OVERLAP_SLOP_CGW){
            SaveEdgeMeanForLater(mergeEdge);
            mergeEdge->distance.mean = overlapEdge->distance.mean;
          }
        }
      }
      MarkScaffoldsForMerging(mergeEdge,  TRUE /* mark for merging, not simply exclusion from other merges */);
    }
  else
    {
      ScaffoldAlignmentInterface * sai = iSpec->sai;
    
      if(minMergeDistance < -1000000. &&
         scaffoldA->bpLength.mean > 1000000. &&
         scaffoldB->bpLength.mean > 1000000.)
        {
          // edge is too negative given scaffold lengths
          if(verbose)
            fprintf(GlobalData->stderrc,
                    "Edge is too negative for scaffold lengths\n");
          return;
        }
    
#ifdef DEBUG1
      verbose = 1;
#endif
    
      if(isBadScaffoldMergeEdge(curEdge, iSpec->badSEdges))
        {
          if(verbose)
            fprintf(GlobalData->stderrc,
                    "Edge previously marked as bad for merging.\n");
          return;
        }
    
      if(!isQualityScaffoldMergingEdge(curEdge,
                                       scaffoldA,
                                       scaffoldB,
                                       iSpec->sai->scaffInst,
                                       iSpec->MIs,
                                       iSpec->minSatisfied,
                                       iSpec->maxDelta))
        {
          SaveBadScaffoldMergeEdge(curEdge, iSpec->badSEdges);
          return;
        }
    
      mergeEdge = curEdge;
      if(mayOverlap || mustOverlap)
        {
          // may overlap
          ChunkOrient orientEndNodeA, orientEndNodeB;
          ChunkOrientationType edgeEndsOrient;
          NodeCGW_T *endNodeA, *endNodeB;
          double aGapSize, bGapSize;
          EdgeCGW_T *overlapEdge = NULL;
          float chiSquaredValue;
          int alternate = FALSE;
      
          TranslateScaffoldOverlapToContigOverlap(scaffoldA, scaffoldB,
                                                  mergeEdge->orient,
                                                  &endNodeA, &endNodeB,
                                                  &orientEndNodeA,
                                                  &orientEndNodeB,
                                                  &edgeEndsOrient,
                                                  &aGapSize, &bGapSize);
      
          if(verbose)
            {
              fprintf(GlobalData->stderrc,
                      "scaffoldA length: %d, elements: %d, end node + end gap: %d\n",
                      (int) scaffoldA->bpLength.mean,
                      scaffoldA->info.Scaffold.numElements,
                      (int) (endNodeA->bpLength.mean + aGapSize));
              fprintf(GlobalData->stderrc,
                      "scaffoldB length: %d, elements: %d, end node + end gap: %d\n",
                      (int) scaffoldB->bpLength.mean,
                      scaffoldB->info.Scaffold.numElements,
                      (int) (endNodeB->bpLength.mean + bGapSize));
              fprintf(GlobalData->stderrc,
                      "edge length: %.f, variance: %.f, weight: %d\n",
                      mergeEdge->distance.mean, mergeEdge->distance.variance,
                      mergeEdge->edgesContributing);
            }
      
          // take the easy road if the edge is shorter than end nodes/gaps
          if(endNodeA->bpLength.mean + aGapSize > -minMergeDistance &&
             endNodeB->bpLength.mean + bGapSize > -minMergeDistance)
            {
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
              if(verbose)
                {
                  fprintf(GlobalData->stderrc,
                          "Small overlap - work with end nodes\n");
                  if(overlapEdge)
                    fprintf(GlobalData->stderrc, "End node overlap found\n");
                  else
                    fprintf(GlobalData->stderrc, "End node overlap not found\n");
                }
        
              if(overlapEdge != NULL)
                {
                  // overlap edge found. use it
                  SaveEdgeMeanForLater(mergeEdge);

#ifdef DEBUG_MERGE_EDGE_INVERT	  
                  fprintf(stderr,"MERGE EDGE DEBUG: overlapEdge->orient %d edgeEndsOrient %d overlapEdge->idA " F_CID " endNodeA->id " F_CID " ctg lengths %f %f\n",
                          overlapEdge->orient,
                          edgeEndsOrient,
                          overlapEdge->idA,
                          endNodeA->id,
                          (endNodeA->bpLength).mean,
                          (endNodeB->bpLength).mean );
	  
                  PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph,
                                 "PRE OVERLAPEDGE",overlapEdge,overlapEdge->idA);


                  if(edgeEndsOrient!=GetEdgeOrientationWRT(overlapEdge,endNodeA->id)){
                    fprintf(stderr,"SAULS PROPOSAL: need to invert; result would be %f\n",
                            -(endNodeA->bpLength.mean + endNodeB->bpLength.mean + overlapEdge->distance.mean));
                  } else {
                    fprintf(stderr,"SAULS PROPOSAL: no need to invert; result would be %f\n",
                            overlapEdge->distance.mean);
                  }
#endif	  

#ifdef ORIG_MERGE_EDGE_INVERT
                  if(overlapEdge->orient != edgeEndsOrient &&
                     overlapEdge->idA == endNodeA->id)
                    {
                      assert(overlapEdge->orient == InvertEdgeOrient(edgeEndsOrient));
                      overlapEdge->distance.mean =
                        endNodeA->bpLength.mean + endNodeB->bpLength.mean +
                        overlapEdge->distance.mean;
                    }
                  else
                    mergeEdge->distance.mean = overlapEdge->distance.mean;
#else // new merge edge invert, courtesy of TCAG
                  if(edgeEndsOrient!=GetEdgeOrientationWRT(overlapEdge,endNodeA->id))
                    {
                      overlapEdge->distance.mean = 
                        -(endNodeA->bpLength.mean + endNodeB->bpLength.mean +
                          overlapEdge->distance.mean);
                    }
                  mergeEdge->distance.mean = overlapEdge->distance.mean;
#endif
          
#ifdef DEBUG_MERGE_EDGE_INVERT	  
                  PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph,
                                 "POST OVERLAPEDGE",overlapEdge,overlapEdge->idA);
                  PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph,
                                 "MERGEEDGE",mergeEdge,mergeEdge->idA);
#endif          
                  MarkScaffoldsForMerging(mergeEdge, TRUE);
                }
              else if(!iSpec->checkAbutting ||
                      AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec))
                {
                  // no overlap edge, but abutting will work
                  SaveEdgeMeanForLater(mergeEdge);
                  mergeEdge->distance.mean = MAX(-CGW_MISSED_OVERLAP,
                                                 mergeEdge->distance.mean);
                  MarkScaffoldsForMerging(mergeEdge, TRUE);
                }
              else
                {
                  // else don't prevent scaffolds from merging via other edges
                  SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
                }
            }
          else
            {
              // large negative edge - may involve interleaving
              if(verbose)
                fprintf(GlobalData->stderrc,
                        "Large overlap - attempt interleaved merging\n");
        
              //more contigs may be involved in one or both scaffolds
              if(PopulateScaffoldAlignmentInterface(scaffoldA, scaffoldB,
                                                    mergeEdge, sai) == 0)
                {
          
          
                  if(verbose)
                    {
                      fprintf(GlobalData->stderrc, "Populated ScaffoldAlignmentInterface\n");
                      fprintf(GlobalData->stderrc, "%d segments in possible scaffold overlap.\n", (int) GetNumSegmentsInList(sai->segmentList));
                    }
          
#if 0
                  // punt if the edge distance is too negative & there are no overlaps
                  if( (int) GetNumSegmentsInList(sai->segmentList) == 0 &&
                      endNodeA->bpLength.mean + aGapSize < -minMergeDistance &&
                      endNodeB->bpLength.mean + bGapSize < -minMergeDistance)
                    {
                      if(verbose)
                        {
                          fprintf(GlobalData->stderrc,
                                  "Large scaffold overlap with no contig overlaps.\n");
                          fprintf(GlobalData->stderrc, "Edge not usable for merging.\n");
                        }
                      SaveBadScaffoldMergeEdge(mergeEdge,
                                               iSpec->badSEdges);
                      return;
                    }
#endif
          
                  if(verbose)
                    fprintf(stderr,
                            "Trying to find overlap with band %d,%d\n",
                            sai->scaffoldA->bandBeg,sai->scaffoldA->bandEnd);

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

                  if(sai->segmentList != NULL ||             sai->best >= 0)
                    {
                      if(verbose)
                        {
                          fprintf(GlobalData->stderrc, "Align_Scaffold returned best = %d%s\n",
                                  sai->best,
                                  (sai->best == 0) ? " (this should mean pure interleaving)" : "");
                          fprintf(GlobalData->stderrc,
                                  "Adjusting scaffold contig positions\n");
                        }
            
                      if(sai->segmentList != NULL ||
                         (iSpec->checkAbutting &&
                          !AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec)))
                        {
                          if(verbose)
                            {
                              fprintf(GlobalData->stderrc,
                                      "%d segments in scaffold overlap.\n",
                                      GetNumSegmentsInList(sai->segmentList));
                            }
                          // if there are overlaps or abutting isn't an option
                          overlapEdge = MakeScaffoldAlignmentAdjustments(scaffoldA,
                                                                         scaffoldB,
                                                                         mergeEdge, sai);

                          if(overlapEdge != NULL)
                            {
                              SaveEdgeMeanForLater(mergeEdge);
                              mergeEdge->distance.mean = overlapEdge->distance.mean;
                              if(verbose)
                                fprintf(GlobalData->stderrc,
                                        "Made scaffold adjustments for merging.\n");
                              MarkScaffoldsForMerging(mergeEdge, TRUE);
                            }
                          else
                            {
                              // else don't prevent scaffolds from merging via other edges
                              if(verbose)
                                fprintf(GlobalData->stderrc,
                                        "Failed to make scaffold adjustments for merging.\n");
                              SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
                            }
                        }
                      else
                        {
                          // no overlaps && abutting will work
                          SaveEdgeMeanForLater(mergeEdge);
                          mergeEdge->distance.mean = MAX(-CGW_MISSED_OVERLAP, mergeEdge->distance.mean);
                          MarkScaffoldsForMerging(mergeEdge, TRUE);
                        }
                    }
                  else
                    {
                      // if here, no overlap or interleaving possible
                      if(verbose)
                        fprintf(GlobalData->stderrc, "Align_Scaffold returned best = %d\n",
                                sai->best);
            
                      /* if edge is not highly negative, abutt. Otherwise abort.
                         criterion is if -20 abutting still leaves edge as trusted
                      */
                      if(!iSpec->checkAbutting ||
                         AbuttingWillWork(mergeEdge, scaffoldA, scaffoldB, iSpec))
                        {
                          SaveEdgeMeanForLater(mergeEdge);
                          mergeEdge->distance.mean = MAX(-CGW_MISSED_OVERLAP,
                                                         mergeEdge->distance.mean);
                          if(verbose)
                            fprintf(GlobalData->stderrc, "Abutting will work.\n");
                          MarkScaffoldsForMerging(mergeEdge, TRUE);
                        }
                      else
                        {
                          // else don't prevent scaffolds from merging via other edges
                          // record that this scaffold overlap is bad
                          if(verbose)
                            fprintf(GlobalData->stderrc, "Edge not usable for merging.\n");
                          SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
                        }
                    }
                }
              else
                {
                  fprintf(GlobalData->stderrc,
                          "Failed to populate scaffold alignment interface!\n");
                  SaveBadScaffoldMergeEdge(mergeEdge, iSpec->badSEdges);
                }
            }
        }
      else
        {
          // if here, edge was non-negative
          MarkScaffoldsForMerging(mergeEdge, TRUE);
        }
    }
}


void ExamineUsableSEdgeSet(VA_TYPE(PtrT) *sEdges,
                           int minWeightThreshold,
                           InterleavingSpec * iSpec,
                           int verbose)
{
  static int32 count = 0;
  int i;
  
  for(i = 0; i < GetNumPtrTs(sEdges); i++, count++){
    SEdgeT *curEdge = *(SEdgeT **)GetPtrT(sEdges,i);
    
    if(count != 0 && count % 10000 == 0){
      CDS_CID_t edgeIndex =
        (CDS_CID_t)GetVAIndex_SEdgeT(ScaffoldGraph->ScaffoldGraph->edges, curEdge);
      fprintf(GlobalData->stderrc,
              "* MergeScaffoldsAggressive examining edge %d  " F_CID "\n",
              count, edgeIndex);
    }
    
    if(curEdge->edgesContributing < minWeightThreshold)
      continue;
    
    ExamineSEdgeForUsability(sEdges, curEdge, iSpec, verbose);
#ifdef CHECKPOINT_ON_DEMAND
    CheckpointOnDemand(EXIT_AFTER_CHECKPOINTING);
#endif
  }
}

#define EDGE_WEIGHT_FACTOR  2

void ExamineUsableSEdges(VA_TYPE(PtrT) *sEdges,
                         int minWeightThreshold,
                         InterleavingSpec * iSpec,
                         int verbose)
{
  if (GetNumPtrTs(sEdges) > 0)
    {
      SEdgeT ** sEdge = (SEdgeT **)GetPtrT(sEdges,0);
      int32 maxWeightEdge = 0;
      int i;
    
      if(! GlobalData->doInterleavedScaffoldMerging)
        {
          // find the max weight non-negative distance edge
          for(i = 0; i < GetNumVA_PtrT(sEdges); i++)
            {
              if(sEdge[i]->distance.mean > 0)
                {
                  if(isBadScaffoldMergeEdge(sEdge[i], iSpec->badSEdges))continue;

                  maxWeightEdge = sEdge[i]->edgesContributing;
                  break;
                }
            }
        }
      else
        {

          for(i = 0; i < GetNumVA_PtrT(sEdges); i++)
            {
              if(isBadScaffoldMergeEdge(sEdge[i], iSpec->badSEdges))continue;
              maxWeightEdge = sEdge[i]->edgesContributing;
              break;
            }
      
        }
    
      if ( maxWeightEdge/EDGE_WEIGHT_FACTOR < minWeightThreshold )
        {
          minWeightThreshold = maxWeightEdge / EDGE_WEIGHT_FACTOR;
        }
      minWeightThreshold = MAX( minWeightThreshold, EDGE_WEIGHT_FACTOR);
      fprintf(GlobalData->stderrc,
              "* Considering edges with weight >= %d (maxWeightEdge/%d: %d)\n",
              minWeightThreshold, (int) EDGE_WEIGHT_FACTOR,
              (int) (maxWeightEdge/EDGE_WEIGHT_FACTOR));
    }
  
  // examine the edges
  ExamineUsableSEdgeSet(sEdges, minWeightThreshold, iSpec, verbose);
}


void BuildMergedScaffoldEdges(ScaffoldGraphT * graph)
{
  // create a scaffold edge for every inter-scaffold contig edge
  fprintf(GlobalData->stderrc, "* Building scaffold edges from scratch\n");
  BuildSEdges(graph, GlobalData->doInterleavedScaffoldMerging);
  
  // merge all inter-scaffold edges that are compatible with each other
  fprintf(GlobalData->stderrc, "* Merging scaffold edges from scratch\n");
  MergeAllGraphEdges(graph->ScaffoldGraph, TRUE);// Merge 'em
}


int BuildSEdgesForMerging(ScaffoldGraphT * graph,
                          VA_TYPE(PtrT) ** sEdges,
                          VA_TYPE(PtrT) ** overlapSEdges,
                          float * minWeightThreshold,
                          int adjustThreshold,
                          InterleavingSpec * iSpec,
                          int32 verbose)
{
  int firstTime;
  
  fprintf(GlobalData->stderrc, "* Building SEdges for merging\n");
  fflush(GlobalData->stderrc);
  
  if(GetNumGraphEdges(graph->ScaffoldGraph) == 0)
    {
      if(verbose){
        fprintf(GlobalData->stderrc, "No SEdges!!!\n");
      }
      return 1;
    }
  
  if(*sEdges == NULL){
    *sEdges = CreateVA_PtrT(100);
    *overlapSEdges = CreateVA_PtrT(100);
    firstTime = TRUE;
  }else{
    ResetVA_PtrT(*sEdges);
    ResetVA_PtrT(*overlapSEdges);
    firstTime = FALSE;
  }
  
  // loop over scaffolds to find all merge candidates and
  // create a list of pointers to scaffold edges that
  // may be used for merging
  fprintf(GlobalData->stderrc, "* Building usable SEdges\n");
  fflush(GlobalData->stderrc);
  BuildUsableSEdges(*sEdges, *overlapSEdges, verbose);
  
  if(adjustThreshold)
    {
      if (firstTime == TRUE && GetPtrT(*sEdges, 0) != NULL)
        {
          SEdgeT ** sEdge = (SEdgeT **)GetPtrT(*sEdges,0);
          if (GlobalData->doInterleavedScaffoldMerging)
            {
              int j;
              for(j = 0; j < GetNumVA_PtrT(*sEdges); j++)
                {
                  if(sEdge[j]->distance.mean > 0)
                    {
                      *minWeightThreshold =
                        sEdge[j]->edgesContributing / EDGE_WEIGHT_FACTOR;
                      break;
                    }
                }
            }
          else
            {
              *minWeightThreshold = sEdge[0]->edgesContributing / EDGE_WEIGHT_FACTOR;
            }
          fprintf( GlobalData->stderrc,
                   "initially setting minWeightThreshold to %f\n",
                   *minWeightThreshold);
        }
      else
        *minWeightThreshold -= 0.2;

      // don't let it get to small and go negative
      *minWeightThreshold = MAX( *minWeightThreshold, EDGE_WEIGHT_FACTOR);
    }

  // loop over sEdges
  // mark scaffolds for merging & associate scaffold ends
  fprintf(GlobalData->stderrc, "* Examining usable SEdges\n");
  fflush(GlobalData->stderrc);
  ExamineUsableSEdges(*sEdges, (int) *minWeightThreshold, iSpec, verbose);
  
  return 0;
}


int SetUpSEdges(ScaffoldGraphT * graph,
                VA_TYPE(PtrT) ** sEdges,
                VA_TYPE(PtrT) ** overlapSEdges,
                float * minWeightThreshold,
                int adjustThreshold,
                InterleavingSpec * iSpec,
                int32 verbose)
  
{
  int retVal;
  
  BuildMergedScaffoldEdges(graph);
  
  retVal = BuildSEdgesForMerging(graph,
                                 sEdges, overlapSEdges,
                                 minWeightThreshold, adjustThreshold,
                                 iSpec, verbose);
  
  return retVal;
}


void RemoveDeadRefsFromSEdge(ScaffoldGraphT * graph, SEdgeT * sEdge)
{
  if(!sEdge->flags.bits.isRaw)
    {
      SEdgeT * raw = sEdge;
      SEdgeT * lastRaw = sEdge;
      while((raw = GetGraphEdge(graph->ScaffoldGraph, raw->nextRawEdge)) != NULL)
        {
          // referenceEdge references inducing contig edge
          // topLevelEdge references top contig edge
          CIEdgeT * ciEdge = GetGraphEdge(graph->RezGraph, raw->referenceEdge);
          if(ciEdge->idA == NULLINDEX || ciEdge->idB == NULLINDEX)
            {
              lastRaw->nextRawEdge = raw->nextRawEdge;
              raw = lastRaw;
            }
        }
    }
  else
    {
      // do what?
    }
}


/**********************************************************/
int MergeScaffolds(VA_TYPE(CDS_CID_t) * deadScaffoldIDs,
                   InterleavingSpec * iSpec,
                   int32 verbose)
{
  int mergedSomething = FALSE;
  GraphNodeIterator scaffolds;
  CIScaffoldT *thisScaffold;
  CDS_CID_t thisScaffoldID; /* The index of the thisScaffold. We have to be careful about thisScaffold since
                               it is a pointer to a an element of the scaffolds array, which may be reallocated
                               during the loop. */
  CDS_CID_t currentSetID = 0;
  int cntScaffold = 1;
  int cnt;
#ifdef INSTRUMENT_CGW
  float fractMatesHappyBefore;
  float fractMatesHappyAfter;
  MateInstrumenter matesBefore;
  MateInstrumenter matesAfter;
  ScaffoldInstrumenter * scaff_inst;
#endif
#ifdef CHECK_CONTIG_ORDERS
  ContigOrientChecker * coc;
#endif
  
#ifdef INSTRUMENT_CGW
  scaff_inst = CreateScaffoldInstrumenter(ScaffoldGraph, INST_OPT_ALL);
  assert(scaff_inst != NULL);
#endif
  
#ifdef CHECK_CONTIG_ORDERS
  coc = CreateContigOrientChecker();
  assert(coc != NULL);
#endif
  
  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph,
                        GRAPH_NODE_DEFAULT);
  while((thisScaffold = NextGraphNodeIterator(&scaffolds)) != NULL)
    {
      CIScaffoldT *AendScaffold, *BendScaffold;
      SEdgeT *edgeA, *edgeB, *edge;
      CIScaffoldT *neighbor;
      CDS_CID_t neighborID;
      const LengthT nullLength = {0.0, 0.0};
      LengthT currentOffset = nullLength;
      CIScaffoldT CIScaffold;
      ChunkOrient orientCI;
      CDS_CID_t newScaffoldID;
      int numMerged;
    
      thisScaffoldID = thisScaffold->id;
    
      if(((cntScaffold++) % 10000) == 0){
        fprintf(GlobalData->stderrc,"* Examining scaffold %d " F_CID "\n",
                cntScaffold, thisScaffoldID);
        fflush(GlobalData->stderrc);
      
      }
      if(thisScaffold->type != REAL_SCAFFOLD){
        continue;
      }
      if(thisScaffold->setID != NULLINDEX){
        continue;// This Scaffold has already been placed in a Scaffold.
      }
      if(thisScaffold->essentialEdgeA != NULLINDEX){
        edgeA = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeA);
        AendScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                    (edgeA->idA == thisScaffold->id) ? edgeA->idB : edgeA->idA);
      }else{
        edgeA = (SEdgeT *)NULL;
        AendScaffold = (CIScaffoldT *)NULL;
      }
      if(thisScaffold->essentialEdgeB != NULLINDEX){
        edgeB = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeB);
        BendScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                    (edgeB->idA == thisScaffold->id) ? edgeB->idB : edgeB->idA);
      }else{
        edgeB = (SEdgeT *)NULL;
        BendScaffold = (CIScaffoldT *)NULL;
      }
      if((AendScaffold != (CIScaffoldT *)NULL) &&
         (BendScaffold != (CIScaffoldT *)NULL)){
        continue;// This CI is not a starting point for a Scaffold.
      }
      if(BendScaffold != (CIScaffoldT *)NULL){
        orientCI = A_B;
        edge = edgeB;
        neighbor = BendScaffold;
        neighborID = neighbor->id;
      }else if(AendScaffold != (CIScaffoldT *)NULL){
        orientCI = B_A;
        edge = edgeA;
        neighbor = AendScaffold;
        neighborID = neighbor->id;
      }else{// Singleton Scaffold
        continue;
      }

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
      CIScaffold.aEndCoord = CIScaffold.bEndCoord = -1;
      CIScaffold.numEssentialA = CIScaffold.numEssentialB = 0;
      CIScaffold.essentialEdgeB = CIScaffold.essentialEdgeA = NULLINDEX;
      CIScaffold.microhetScore = 0.0;
      CIScaffold.setID = NULLINDEX;
      thisScaffold->flags.bits.isDead = TRUE;  // Mark the old scaffold dead

      AppendGraphNode(ScaffoldGraph->ScaffoldGraph, &CIScaffold);  /* Potential realloc of ScaffoldGraph->ScaffoldGraph->nodes */

      thisScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, thisScaffoldID);
      neighbor     = GetGraphNode(ScaffoldGraph->ScaffoldGraph, neighborID);
    
      if(verbose){
        fprintf(GlobalData->stderrc,"* START: Inserting scaffold " F_CID " into scaffold " F_CID "\n", thisScaffold->id, newScaffoldID);
      }
      assert(thisScaffold->bpLength.variance >= 0);
      assert(neighbor->bpLength.variance >= 0);
    
      cnt = 0;
    
#ifdef INSTRUMENT_CGW
      {
        fprintf(GlobalData->stderrc,
                "* BEFORE Inserting scaffolds " F_CID " and " F_CID " into scaffold " F_CID "\n",
                thisScaffold->id, neighbor->id, newScaffoldID);
        InstrumentScaffold(ScaffoldGraph, thisScaffold,
                           scaff_inst,
                           InstrumenterVerbose2, GlobalData->stderrc);
        GetMateInstrumenterFromScaffoldInstrumenter(&matesBefore,
                                                    scaff_inst);
      
        InstrumentScaffold(ScaffoldGraph, neighbor,
                           scaff_inst,
                           InstrumenterVerbose2, GlobalData->stderrc);
        AddMateInstrumenterCounts(&matesBefore, &(scaff_inst->mates));
      
        if(GetMateStatsBad(&(matesBefore.intra)) +
           GetMateStatsBad(&(matesBefore.inter)) +
           GetMateStatsHappy(&(matesBefore.intra)) +
           GetMateStatsHappy(&(matesBefore.inter)) > 0)
          {
            fractMatesHappyBefore =
              ((float) (GetMateStatsHappy(&(matesBefore.intra)) +
                        GetMateStatsHappy(&(matesBefore.inter)))) /
              (GetMateStatsBad(&(matesBefore.intra)) +
               GetMateStatsBad(&(matesBefore.inter)) +
               GetMateStatsHappy(&(matesBefore.intra)) +
               GetMateStatsHappy(&(matesBefore.inter)));
          }
        else
          {
            fractMatesHappyBefore = 1.0;
          }
      }
#endif
    
#ifdef CHECK_CONTIG_ORDERS
      ResetContigOrientChecker(coc);
      AddScaffoldToContigOrientChecker(ScaffoldGraph, thisScaffold, coc);
#endif
    
      AppendVA_CDS_CID_t(deadScaffoldIDs, &(thisScaffold->id));
    
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
      mergedSomething = TRUE;
      while(neighbor != (CIScaffoldT *)NULL){
        ChunkOrientationType edgeOrient = GetEdgeOrientationWRT(edge,
                                                                thisScaffold->id);
      
        cnt++;
        assert(orientCI == A_B || orientCI == B_A);
        assert(edge->distance.variance >= 0);
        assert(currentOffset.mean >= 0);
        assert(currentOffset.variance >= 0);
        assert(thisScaffold->bpLength.variance >= 0);
        assert(neighbor->bpLength.variance >= 0);
      
        neighborID = neighbor->id;
      
        if(orientCI == A_B){
          if(edgeOrient == AB_AB){
            orientCI = A_B;
          }else if(edgeOrient == AB_BA){
            orientCI = B_A;
          }else{
            assert(0);
          }
        }else{// orientCI == B_A
          assert(orientCI == B_A);
          if(edgeOrient == BA_AB){
            orientCI = A_B;
          }else if(edgeOrient == BA_BA){
            orientCI = B_A;
          }else{
            assert(0);
          }
        }
        thisScaffold = neighbor;
        thisScaffoldID = thisScaffold->id;
      
        assert(thisScaffold->bpLength.variance >= 0);
        assert(edge->distance.variance >= 0);
        assert(currentOffset.mean >= 0);
        assert(currentOffset.variance >= 0);
        assert(orientCI == A_B || orientCI == B_A);
      
        currentOffset.mean += edge->distance.mean;
        currentOffset.variance += edge->distance.variance;
      
        if(edge->distance.mean < -CGW_MISSED_OVERLAP)
          {
            CIScaffoldT * newScaffold =
              GetCIScaffoldT(ScaffoldGraph->CIScaffolds, newScaffoldID);
            if(currentOffset.mean < 0.0)
              {
                /*
                  0
                  newScaffold:                --------------->
                  thisScaffold:      ------------------
                */
                CIScaffoldTIterator CIs;
                ChunkInstanceT * CI;
                double variance = GetVarianceOffset(thisScaffold,
                                                     -currentOffset.mean,
                                                     orientCI == A_B);
          
                fprintf(stderr, "Adjusting newScaffold offsets by mean - %f variance + %f\n", currentOffset.mean, variance);
          
                InitCIScaffoldTIterator(ScaffoldGraph, newScaffold,
                                        TRUE, FALSE, &CIs);
                while((CI = NextCIScaffoldTIterator(&CIs)) != NULL)
                  {
                    CI->offsetAEnd.mean -= currentOffset.mean;
                    CI->offsetAEnd.variance += variance;
                    CI->offsetBEnd.mean -= currentOffset.mean;
                    CI->offsetBEnd.variance += variance;
                  }
                currentOffset.mean = currentOffset.variance = 0.0;
              }
            else
              {
                /*
                  0
                  newScaffold:     --------------->
                  thisScaffold:              ------------------
                */
                currentOffset.variance = GetVarianceOffset(newScaffold,
                                                           currentOffset.mean,
                                                           TRUE);
              }
          }
        assert(currentOffset.variance >= 0);
      
        if(verbose){
          fprintf(GlobalData->stderrc,"* Adding to scaffold " F_CID " scaffold " F_CID " at orient %c offset %g\n",
                  newScaffoldID, thisScaffoldID, orientCI, currentOffset.mean);
        }
#ifdef CHECK_CONTIG_ORDERS
        AddScaffoldToContigOrientChecker(ScaffoldGraph, thisScaffold, coc);
#endif
      
        AppendVA_CDS_CID_t(deadScaffoldIDs, &(thisScaffold->id));
      
        if (0 == InsertScaffoldContentsIntoScaffold(ScaffoldGraph,
                                                    newScaffoldID, thisScaffold->id,
                                                    orientCI, &currentOffset,
                                                    iSpec->contigNow)) {

          //  Shucks, we failed to insert the scaffold.  Remove all the
          //  work we did, and return that nothing got merged.
          //
          //  OK, it might just be easier to reconstruct this function,
          //  but only have it return success/fail, and do no other
          //  modifications.

        }

        if(iSpec->contigNow != NO_CONTIGGING)
          {
            // remove references to edges of dead contigs from scaffold edge
            RemoveDeadRefsFromSEdge(ScaffoldGraph, edge);
          }
        else
          {
            // need to set edge statuses so scaffold is clone-connected
            MarkUnderlyingCIEdgesTrusted(ScaffoldGraph, edge);
          }
      
        neighbor = GetGraphNode(ScaffoldGraph->ScaffoldGraph, neighborID);
      
        assert(orientCI == A_B || orientCI == B_A);
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
        if(orientCI == A_B){
          if(thisScaffold->essentialEdgeB != NULLINDEX){
            edge = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeB);
            neighbor = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                    (edge->idA == thisScaffold->id) ? edge->idB : edge->idA);
            neighborID = neighbor->id;
          }else{// End of Scaffold
            edge = (SEdgeT *)NULL;
            neighbor = (CIScaffoldT *)NULL;
            neighborID = NULLINDEX;
          }
        }else{// orientCI == B_A
          if(thisScaffold->essentialEdgeA != NULLINDEX){
            edge = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeA);
            neighbor = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                    (edge->idA == thisScaffold->id) ? edge->idB : edge->idA);
            neighborID = neighbor->id;
          }else{// End of Scaffold
            edge = (SEdgeT *)NULL;
            neighbor = (CIScaffoldT *)NULL;
            neighborID = NULLINDEX;
          }
        }
      }  //  while(neighbor != (CIScaffoldT *)NULL)

      // New scaffold fully popuated now.
    
      if(iSpec->contigNow != NO_CONTIGGING &&
         GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                      newScaffoldID)->info.Scaffold.numElements > 1)
        {
          int status = RECOMPUTE_SINGULAR;
          int recomputeIteration = 0;
          while(recomputeIteration < 3 &&
                (status == RECOMPUTE_SINGULAR ||
                 status == RECOMPUTE_CONTIGGED_CONTAINMENTS))
            {
              // need to make sure scaffold is connected with trusted raw edges

#ifdef BPWDEBUG
              CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                                   newScaffoldID);
              dumpTrustedEdges(ScaffoldGraph, scaffold, ALL_TRUSTED_EDGES);
#endif

              MarkInternalEdgeStatus(ScaffoldGraph,
                                     GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                                  newScaffoldID),
                                     PAIRWISECHI2THRESHOLD_CGW,
                                     1000.0 * SLOPPY_EDGE_VARIANCE_THRESHHOLD,
                                     TRUE, TRUE, 0, TRUE);

#ifdef BPWDEBUG
              dumpTrustedEdges(ScaffoldGraph, scaffold, ALL_TRUSTED_EDGES);
#endif

#ifdef CHECKCONNECTED
              assert(IsScaffoldInternallyConnected(ScaffoldGraph,
                                                   GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                                                newScaffoldID),
                                                   ALL_EDGES));
#endif

              status =
                RecomputeOffsetsInScaffold(ScaffoldGraph,
                                           GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                                        newScaffoldID),
                                           TRUE, TRUE, FALSE);
              recomputeIteration++;
            }
          if(status != RECOMPUTE_OK)
            {
              fprintf(GlobalData->stderrc,
                      "ReomputeOffsetsInScaffold failed (%d) "
                      "for scaffold " F_CID " in MergeScaffolds\n",
                      status, newScaffoldID);
            }
        }  //  if (iSpec->contigNow != NO_CONTIGGING && .....)
    
#if defined(INSTRUMENT_CGW) || defined(CHECK_CONTIG_ORDERS)
      {
        CIScaffoldT *newScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                                newScaffoldID);
      
#ifdef INSTRUMENT_CGW
        {
          fprintf(GlobalData->stderrc,"* AFTER inserting scaffold\n");
        
          // instrument the resulting scaffold
          InstrumentScaffold(ScaffoldGraph, newScaffold,
                             scaff_inst,
                             InstrumenterVerbose2, GlobalData->stderrc);
          GetMateInstrumenterFromScaffoldInstrumenter(&matesAfter,
                                                      scaff_inst);
        
          fractMatesHappyAfter =
            ((float) (GetMateStatsHappy(&(matesAfter.intra)) +
                      GetMateStatsHappy(&(matesAfter.inter)))) /
            (GetMateStatsBad(&(matesAfter.intra)) +
             GetMateStatsBad(&(matesAfter.inter)) +
             GetMateStatsHappy(&(matesAfter.intra)) +
             GetMateStatsHappy(&(matesAfter.inter)));
        
          fprintf(GlobalData->stderrc,
                  "%.2f%%, %.2f%% mates satisfied before, after\n",
                  100 * fractMatesHappyBefore,
                  100 * fractMatesHappyAfter);
        
          if(iSpec->maxDelta > 0.0 &&
             (fractMatesHappyAfter + iSpec->maxDelta < fractMatesHappyBefore) ||
             (fractMatesHappyAfter < fractMatesHappyBefore &&
              fractMatesHappyAfter < iSpec->minSatisfied
#ifdef REQUIRE_MORE_THAN_ONE_BAD_TO_REJECT
              &&  GetMateStatsBad(&(matesAfter.inter)) +  GetMateStatsBad(&(matesAfter.intra))  - 
              (GetMateStatsBad(&(matesBefore.inter)) +  GetMateStatsBad(&(matesBefore.intra)))
              >1
#endif
#ifdef REQUIRE_BAD_APPROACHING_HAPPY_TO_REJECT
              &&( (          (float) (GetMateStatsHappy(&(matesAfter.inter)) +  GetMateStatsHappy(&(matesAfter.intra))  - 
                                      GetMateStatsHappy(&(matesBefore.inter)) +  GetMateStatsHappy(&(matesBefore.intra)))
                             <= 0.) ? TRUE : 
                  ((float) (GetMateStatsBad(&(matesAfter.inter)) +  GetMateStatsBad(&(matesAfter.intra))  - 
                            GetMateStatsBad(&(matesBefore.inter)) +  GetMateStatsBad(&(matesBefore.intra))) 
                   /
                   (float) (GetMateStatsHappy(&(matesAfter.inter)) +  GetMateStatsHappy(&(matesAfter.intra))  - 
                            GetMateStatsHappy(&(matesBefore.inter)) +  GetMateStatsHappy(&(matesBefore.intra))) 
                   > MAX_FRAC_BAD_TO_GOOD)
                  )
#endif  // REQUIRE_BAD
              )
             )
            {
              fprintf(GlobalData->stderrc,
                      "Scaffold merging results did not meet expectations!\n");
            }
        }
#endif // INSTRUMENT_CGW
      
#ifdef CHECK_CONTIG_ORDERS
        CompareNewOrientations(ScaffoldGraph, newScaffold, coc);
#endif
      }
#endif // INSTRUMENT_CGW || CHECK_CONTIG_ORDERS
    
      ScaffoldGraph->numLiveScaffolds += (1 - numMerged);
      currentSetID++;
    }
  
#ifdef INSTRUMENT_CGW
  DestroyScaffoldInstrumenter(scaff_inst);
#endif
#ifdef CHECK_CONTIG_ORDERS
  DestroyContigOrientChecker(coc);
#endif
  
  return mergedSomething;
}


void DeleteScaffoldEdgesForScaffold(ScaffoldGraphT * graph,
                                    CIScaffoldT * scaffold)
{
  while(scaffold->edgeHead != NULLINDEX)
    {
      DeleteGraphEdge(graph->ScaffoldGraph,
                      GetGraphEdge(graph->ScaffoldGraph,
                                   scaffold->edgeHead));
    }
}


void DeleteDeadScaffoldEdges(ScaffoldGraphT * graph,
                             VA_TYPE(CDS_CID_t) * deadScaffoldIDs)
{
  int i;
  
  for(i = 0; i < GetNumVA_int32(deadScaffoldIDs); i++)
    {
      CDS_CID_t * scaffoldID = GetVA_CDS_CID_t(deadScaffoldIDs, i);
      CIScaffoldT * scaffold = GetCIScaffoldT(graph->CIScaffolds, *scaffoldID);
    
      assert(scaffold->flags.bits.isDead);
    
      DeleteScaffoldEdgesForScaffold(graph, scaffold);
    }
  RecycleDeletedGraphElements(graph->ScaffoldGraph);
}


void ResetScaffoldAndEdgeFlags(ScaffoldGraphT * graph)
{
  // ZeroEdgeWeights is called in FindMoreAttractiveMergeEdge
  int i;
  
  for(i = 0; i < GetNumGraphEdges(graph->ScaffoldGraph); i++)
    {
      EdgeCGW_T * edge = GetGraphEdge(graph->ScaffoldGraph, i);
      edge->quality = 1.f;
      edge->flags.bits.isBogus = 0;
      edge->flags.bits.isProbablyBogus = 0;
      if(edge->flags.bits.MeanChangedByWalking == TRUE)
        {
          edge->flags.bits.MeanChangedByWalking = FALSE;
          edge->distance.mean = edge->minDistance;
        }
    }
  
  // Reset merge marks
  for(i = 0; i < GetNumGraphNodes(graph->ScaffoldGraph); i++)
    {
      CIScaffoldT * scaffold = GetCIScaffoldT(graph->CIScaffolds, i);
      scaffold->essentialEdgeA = scaffold->essentialEdgeB = NULLINDEX;
      scaffold->numEssentialA = scaffold->numEssentialB = 0;
      scaffold->flags.bits.smoothSeenAlready = FALSE;
      scaffold->flags.bits.walkedAlready = FALSE;
      scaffold->setID = NULLINDEX;
    }
}


void BuildNewScaffoldEdges(ScaffoldGraphT * graph,
                           CDS_CID_t firstScaffoldID)
{
  CDS_CID_t i;
  SEdgeBuildStats buildStats;
  
  // Reset scaffold edge heads
  for(i = firstScaffoldID; i < GetNumGraphNodes(graph->ScaffoldGraph); i++)
    {
      CIScaffoldT * newScaffold = GetCIScaffoldT(graph->CIScaffolds, i);
      newScaffold->edgeHead = NULLINDEX;
    }
  
  // build raw edges
  memset(&buildStats, 0, sizeof(buildStats));
  for(i = firstScaffoldID; i < GetNumGraphNodes(graph->ScaffoldGraph); i++)
    {
      CIScaffoldT * newScaffold = GetCIScaffoldT(graph->CIScaffolds, i);
    
      assert(newScaffold != NULL);

      if ((newScaffold->flags.bits.isDead) ||
          (newScaffold->type != REAL_SCAFFOLD))
        continue;
    
      BuildSEdgesForScaffold(graph, newScaffold, FALSE,
                             GlobalData->doInterleavedScaffoldMerging,
                             &buildStats);
      MergeNodeGraphEdges(graph->ScaffoldGraph, newScaffold, TRUE, TRUE, FALSE);
    }
  fprintf(stderr, "Added %d raw scaffold edges\n", buildStats.edgesSucceeded);
}


int MergeScaffoldsExhaustively(ScaffoldGraphT * graph,
                               InterleavingSpec * iSpec,
                               int logicalcheckpointnumber,
                               int verbose)
{
  static VA_TYPE(PtrT) *sEdges = NULL;
  static VA_TYPE(PtrT) *overlapSEdges = NULL;
  static VA_TYPE(CDS_CID_t) * deadScaffoldIDs = NULL;
  int mergedSomething = FALSE;
  int32 iterations = 0;
  time_t t;
  float minWeightThreshold = 0.0;
  CDS_CID_t prevFirstNewScaffoldID = NULLINDEX;
  int32 totalMerged = 0;

  int    buildEdgeCounter = -1;
  time_t lastCkpTime      = time(0) - 90 * 60;

  if(deadScaffoldIDs == NULL)
    deadScaffoldIDs = CreateVA_CDS_CID_t(10000);

  // loop until nothing gets merged
  do{
    CDS_CID_t currFirstNewScaffoldID;
    
    //  Checkpoint periodically - every two hours seems nice!  The
    //  first checkpoint is done after 30 minutes of work here,
    //  though.
    //
    if (time(0) - lastCkpTime > 120 * 60) {
      CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);
      fprintf(GlobalData->timefp,"Checkpoint %d written written during MergeScaffoldsAggressive at iteration %d\n",
              ScaffoldGraph->checkPointIteration, iterations);
      CheckpointScaffoldGraph(ScaffoldGraph, logicalcheckpointnumber);
      lastCkpTime = time(0);
    }
    currFirstNewScaffoldID = GetNumGraphNodes(graph->ScaffoldGraph);
    
    t = time(0);
    fprintf(GlobalData->stderrc, "* MergeScaffoldsAggressive iteration %d at %s\n", 
            iterations, ctime(&t));
    fprintf(GlobalData->stderrc, "* " F_CID " scaffolds and %d scaffold edges in the graph\n",
            currFirstNewScaffoldID, (int) GetNumGraphEdges(graph->ScaffoldGraph));
    fflush(GlobalData->stderrc);

    iterations++;

    if(verbose){
      char buffer[2048];
      sprintf(buffer,"MSA_%d", iterations);
      CelamyAssembly(buffer);
      CelamyCIScaffolds(buffer, ScaffoldGraph);
    }
    
    //  Initially & periodically build scaffold edges from scratch --
    //  assumes that buildEdgeCounter is initialized to -1.  When
    //  switching the checkpoint above from a fixed iteration to time
    //  based, we needed to change the logic here.  Versions previous
    //  to 2005-09-22 (the date of the checkpoint change) would
    //  SetUpSEdges() on the 1st, 8th, 16th, etc iteration, which we
    //  preserve with the following unnatural test.
    //
    buildEdgeCounter++;
    if ((buildEdgeCounter == 0) || (buildEdgeCounter == 8))
      {
        buildEdgeCounter = 0;
        if(SetUpSEdges(graph, &sEdges, &overlapSEdges, &minWeightThreshold,
                       TRUE, iSpec, verbose))
          {
            fprintf(GlobalData->stderrc,
                    "No additional scaffold merging is possible.\n");
            break;
          }
      }
    else 
      {
        // delete all edges to dead (merged) scaffolds
        DeleteDeadScaffoldEdges(graph, deadScaffoldIDs);

        // raw & merged edges from new scaffolds to all other scaffolds
        BuildNewScaffoldEdges(graph, prevFirstNewScaffoldID);

        // flags used in marking scaffolds for merging
        ResetScaffoldAndEdgeFlags(graph);
      
#ifdef CHECKPOINT_ON_DEMAND
        CheckpointOnDemand(RETURN_AFTER_CHECKPOINTING);
#endif
  
        if(BuildSEdgesForMerging(graph,
                                 &sEdges, &overlapSEdges,
                                 &minWeightThreshold, TRUE,
                                 iSpec,
                                 verbose))
          {
            fprintf(GlobalData->stderrc,
                    "No additional scaffold merging is possible.\n");
            break;
          }
      }

    ResetVA_CDS_CID_t(deadScaffoldIDs);

    mergedSomething = MergeScaffolds(deadScaffoldIDs, iSpec, verbose);
    totalMerged += mergedSomething;
    prevFirstNewScaffoldID = currFirstNewScaffoldID;
    if(mergedSomething == 0 && minWeightThreshold > 2.0)
      {
        mergedSomething = 1;
        minWeightThreshold -= 1.0;
      }
  }while(mergedSomething);
  return totalMerged;
}



void DeleteContigOverlapEdges(void)
{
  int numNulls = 0;
  int32 i;
  int32 numEdges = GetNumGraphEdges(ScaffoldGraph->RezGraph);
  for(i = 0; i < numEdges; i++)
    {
      EdgeCGW_T * edge = GetGraphEdge(ScaffoldGraph->RezGraph, i);
      if(!edge->flags.bits.isDeleted)
        {
          if(edge->fragA == NULLINDEX && edge->fragB == NULLINDEX)
            {
              numNulls++;
              DeleteGraphEdge(ScaffoldGraph->RezGraph, edge);
            }
        }
    }
  
  fprintf(stderr, "Number of null contig edges: %d\n", numNulls);
}

void MergeScaffoldsAggressive(ScaffoldGraphT *graph, int logicalcheckpointnumber, int verbose)
{
  time_t t;
  InterleavingSpec iSpec;
  
  if(verbose)
    fprintf(GlobalData->stderrc, "Starting MergeScaffoldsAggressive\n");

  CheckCIScaffoldTs(ScaffoldGraph);
  CheckCIScaffoldTLengths(ScaffoldGraph);
  fprintf(GlobalData->stderrc, "* Successfully passed checks at beginning of scaffold merging\n");
  fflush(GlobalData->stderrc);
  
#if EMIT_STATS
  OpenStatsFiles();
#endif
  
  iSpec.sai = CreateScaffoldAlignmentInterface();
  iSpec.MIs = CreateVA_PtrT(GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));
  iSpec.maxDelta = -1; //0.005;

  iSpec.minSatisfied = MINSATISFIED_CUTOFF;

  iSpec.badSEdges = CreateChunkOverlapper();
  iSpec.checkAbutting = TRUE;
  iSpec.contigNow = ALL_CONTIGGING;

  if (GlobalData->doInterleavedScaffoldMerging) {

    // merge scaffolds with interleaving

    fprintf(GlobalData->stderrc, "** Merging scaffolds with interleaving.\n");
    fprintf(GlobalData->stderrc, "** MinSatisfied: %f, MaxDelta: %f\n",
            iSpec.minSatisfied, iSpec.maxDelta);

    //GlobalData->aligner = Local_Overlap_AS_forCNS;
    iSpec.checkForTinyScaffolds = FALSE;
    //LeastSquaresGapEstimates(graph, TRUE, FALSE, TRUE, TRUE, FALSE);
    MergeScaffoldsExhaustively(graph, &iSpec, logicalcheckpointnumber, verbose);
    //GlobalData->aligner = DP_Compare;

  } else {

    // merge scaffolds without interleaving

    fprintf(GlobalData->stderrc, "** Merging scaffolds without interleaving.\n");
    fprintf(GlobalData->stderrc, "** MinSatisfied: %f, MaxDelta: %f\n",
            iSpec.minSatisfied, iSpec.maxDelta);

    iSpec.checkForTinyScaffolds = TRUE;
    MergeScaffoldsExhaustively(graph, &iSpec, logicalcheckpointnumber, verbose);
  }
  
 
  DeleteScaffoldAlignmentInterface(iSpec.sai);
  {
    int32 i;
    for(i = 0; i < GetNumVA_PtrT(iSpec.MIs); i++) {
      MateInstrumenter *p = *GetVA_PtrT(iSpec.MIs, i);
      safe_free(p);
    }
    DeleteVA_PtrT(iSpec.MIs);
  }
  DestroyChunkOverlapper(iSpec.badSEdges);
  
  CheckCIScaffoldTs(ScaffoldGraph);
  t = time(0);
  fprintf(GlobalData->stderrc,"* Exiting MSA at %s *\n", ctime(&t));
  fflush(GlobalData->stderrc);
  
#if EMIT_STATS
  CloseStatsFiles();
#endif
  return;
}

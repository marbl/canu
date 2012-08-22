
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
static char *rcsid = "$Id: MergeEdges_CGW.c,v 1.35 2012-08-22 02:00:57 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_CGW_dataTypes.h"
#include "AS_UTL_interval.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ChiSquareTest_CGW.h"

#include <vector>

#define C_SORT

#ifndef C_SORT
#include <algorithm>
#endif

using namespace std;

/* Check that another fragment besides the terminal fragment extends
   by at least the minimal detectable amount into the overlap region. */
static int ConfirmAnotherFragmentOverlap(GraphCGW_T *graph,
                                         CDS_CID_t chunkID,
                                         int endB,
                                         int32 overlap){
  ChunkInstanceT *CI;
  int32 minOffset, maxOffset;

  CI = GetGraphNode(graph, chunkID);

  // If the chunk is a singleton, there are no other fragments!
  if(ScaffoldGraph->tigStore->getNumFrags(CI->id, (graph->type == CI_GRAPH)) == 1)
    return FALSE;

  if(endB){
    maxOffset = (int32)CI->bpLength.mean;
    minOffset = maxOffset + overlap;
  }else{
    minOffset = 0;
    maxOffset = -overlap;
  }
#ifdef DEBUG_CONFIRM
  fprintf(stderr,"* chunk = "F_CID" endB = %d minOffset = "F_S32" maxoffset = "F_S32" overlap = "F_S32"\n",
          chunkID, endB, minOffset, maxOffset,overlap);
#endif

  /* Now we have to find if there are intraChunk or terminal fragments, other
     than the fragment at endB, that overlap the (minOffset,maxOffset)
     interval */
  {
    CIFragT *frag;
    MultiAlignT *ma = ScaffoldGraph->tigStore->loadMultiAlign(CI->id, (graph->type == CI_GRAPH));

    int32 overlap;
    for(uint32 i = 0; i < GetNumIntMultiPoss(ma->f_list); i++){
      IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
      frag = GetCIFragT(ScaffoldGraph->CIFrags, mp->ident);

      /* Skip the extremal read that we're trying to check */
      if((!endB && (frag->flags.bits.chunkLabel == AS_INTERCHUNK_A)) ||
         (endB && (frag->flags.bits.chunkLabel == AS_INTERCHUNK_B)) ||
         frag->flags.bits.chunkLabel == AS_SINGLETON){
        continue;
      }
      if((overlap = IntervalsOverlap((int32)frag->offset3p.mean,
                                     (int32)frag->offset5p.mean,
                                     minOffset,
                                     maxOffset,
                                     CGW_DP_MINLEN)) != 0)
        return TRUE;
    }
  }
  return FALSE;
}


/************************************************************/
/* Check that a weight 2 edge is correctly confirmed        */
static int ConfirmOverlap(GraphCGW_T *graph,
                          CIEdgeT *overlapEdge,
                          CIEdgeT *mateEdge){
  int32 overlap;
  int endB;
  int goodOverlapA = FALSE, goodOverlapB = FALSE;

  if((!mateEdge->flags.bits.hasExtremalAFrag &&
      !mateEdge->flags.bits.hasExtremalBFrag)){
    //fprintf(stderr,"*1*  Confirming edge between chunks ("F_CID","F_CID") -- no extremal fragment\n",
    //    mateEdge->idA, mateEdge->idB);
    return TRUE;
  }

  overlap = (int32)overlapEdge->distance.mean;

  if(mateEdge->flags.bits.hasExtremalAFrag){
    if(mateEdge->orient.isAB_BA() || mateEdge->orient.isAB_AB()){
      endB = TRUE; /* The fragment is on the B end of Chunk A */
    }else{
      endB = FALSE; /* On the A end of Chunk A */
    }
    goodOverlapA = ConfirmAnotherFragmentOverlap(graph, mateEdge->idA, endB, overlap);
  }
  if(mateEdge->flags.bits.hasExtremalBFrag){
    if(mateEdge->orient.isAB_BA() || mateEdge->orient.isBA_BA()){
      endB = TRUE;  /* B end of Chunk B */
    }else{
      endB = FALSE;  /* A end of Chunk B */
    }
    goodOverlapB = ConfirmAnotherFragmentOverlap(graph, mateEdge->idB, endB, overlap);
  }
  return(goodOverlapA && goodOverlapB);
}


//  From numerical recipes in C
//  Produces p-values for Chi-square distributions

#define ITMAX 1000
#define EPS   3.0e-7
#define FPMIN 1.0e-30

static
double
gser(double a, double x) {
  double gln = lgammaf(a);
  double ap  = a;
  double del = 1.0 / a;
  double sum = del;
  int32 n   = 1;

  for (n=1; n<=ITMAX; n++) {
    ++ap;
    del *= x/ap;
    sum += del;
    if (fabs(del) < fabs(sum)*EPS)
      break;
  }

  //fprintf(stderr, "gser(%f, %f)-- OK   sum*exp(-x+a*log(x)-gln) = %f * exp(%f + %f * %f - %f) = %f\n",
  //        a, x, sum, -x, a, log(x), gln, sum * exp( -x + a * log(x) - gln));
  if (n <= ITMAX)
    return sum * exp(-x + a * log(x) - gln);

  fprintf(stderr, "gser(%f, %f)-- WARN sum*exp(-x+a*log(x)-gln) = %f * exp(%f + %f * %f - %f) = %f\n",
          a, x, sum, -x, a, log(x), gln, sum * exp( -x + a * log(x) - gln));
  return(1.0);
  assert(0);
}

static
double
gcf(double a, double x) {
  double gln = lgammaf(a);
  double b   = x + 1.0 - a;
  double c   = 1.0 / FPMIN;
  double d   = 1.0 / b;
  double h   = d;
  double an  = 0;
  double del = 0;
  int32 i   = 1;

  for (i=1; i<=ITMAX; i++) {
    an = -i * (i-a);
    b += 2.0;
    d = an * d + b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c = b + an / c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs(del - 1.0) < EPS)
      break;
  }

  //fprintf(stderr, "gcf(%f, %f)-- OK   exp(-x+a*log(x)-gln)*h = %f * exp(%f + %f * %f - %f) = %f\n",
  //        a, x, -x, a, log(x), h, gln, exp( -x + a * log(x) - gln) * h);
  if (i <= ITMAX)
    return exp(-x + a * log(x) - gln) * h;

  fprintf(stderr, "gcf(%f, %f)-- WARN exp(-x+a*log(x)-gln)*h = %f * exp(%f + %f * %f - %f) = %f\n",
          a, x, -x, a, log(x), h, gln, exp( -x + a * log(x) - gln) * h);
  return(0.0);
  assert(0);
}

static
double
gammq(double a, double x) {

  assert(a >  0.0);
  assert(x >= 0.0);

  if (x == 0.0)
    return(0.0);

  if (x < (a+1.0))
    return(1.0 - gser(a,x));
  else
    return(gcf(a,x));
}


//  Computes a chi squared test on the distribution of the edge distances around the computed mean.
//  If skip is >= 0 and < numEdges then the edge corresponding to that index is ignored. Because we
//  are computing the mean from the edges the number of degrees of freedom is (numEdges - 1) unless
//  we are skipping an edge then it is (numEdges - 2).
//
//  Returns 1 if test SUCCEEDS
//
static
int
ComputeChiSquared(Chi2ComputeT *edges,
                  int           numEdges,
                  CDS_CID_t     skip,
                  LengthT      *distance,
                  double        *score) {
  double  cumScore           = 0.0;
  double  cumWeightedMean    = 0.0;
  double  cumInverseVariance = 0.0;
  double  leastSquaresMean   = 0.0;
  int     numEdgesCorrection = 1;

  for (int32 edgeIndex = 0; edgeIndex < numEdges; edgeIndex++) {
    if (edgeIndex == skip) {
      numEdgesCorrection++;
    } else {
      cumWeightedMean    += edges[edgeIndex].weightedMean;
      cumInverseVariance += edges[edgeIndex].inverseVariance;
    }
  }

  leastSquaresMean = cumWeightedMean / cumInverseVariance;

  distance->mean     = leastSquaresMean;
  distance->variance = 1.0 / cumInverseVariance;

  for (int32 edgeIndex = 0; edgeIndex < numEdges; edgeIndex++) {
    if (edgeIndex != skip)
      cumScore += (((edges[edgeIndex].distance.mean - leastSquaresMean) * 
                    (edges[edgeIndex].distance.mean - leastSquaresMean)) /
                   edges[edgeIndex].distance.variance);
  }

  *score = cumScore;

  //  It'd be a bad day if this fails.
  assert(numEdges > numEdgesCorrection);

  //  This one, on the other hand, fails occasionally.
  if (cumScore < 0.0) {
    fprintf(stderr, "ComputeChiSquared()-- WARNING:  negative variance fails test.\n");
    for (int32 edgeIndex = 0; edgeIndex < numEdges; edgeIndex++) {
      fprintf(stderr, "ComputeChiSquared()--    (cumWeightedMean += %f;  cumInverseVariance += %f\n",
              edges[edgeIndex].weightedMean, edges[edgeIndex].inverseVariance);
    }
    fprintf(stderr, "ComputeChiSquared()--    distance->mean = leastSquaresMean = cumWeightedMean / cumInverseVariance = %f;  distance->variance = %f\n",
            leastSquaresMean, 1.0 / cumInverseVariance);
    for (int32 edgeIndex = 0; edgeIndex < numEdges; edgeIndex++) {
      fprintf(stderr, "ComputeChiSquared()--    (mean - leastSquaresMean) ^ 2 / variance = (%f - %f) ^ 2 / %f = %f%s\n",
              edges[edgeIndex].distance.mean, leastSquaresMean, edges[edgeIndex].distance.variance,
              (((edges[edgeIndex].distance.mean - leastSquaresMean) * 
                (edges[edgeIndex].distance.mean - leastSquaresMean)) /
               edges[edgeIndex].distance.variance),
              (edgeIndex != skip) ? "" : " (SKIPPED)");
    }

    return(FALSE);
  }

  return(gammq(numEdges - numEdgesCorrection, cumScore) > 0.005);
}


/************************************************************************/
static
void
InitializeMergedEdge(CIEdgeT            *newEdge,
                     CIEdgeT            *overlapEdgeAll,
                     LengthT             distance,
                     GraphCGW_T         *graph,
                     vector<CDS_CID_t>  &inputEdges) {

  CIEdgeT *overlapEdge    = NULL;
  CIEdgeT *nonOverlapEdge = NULL;

  //  Fill out the new merged edge

  *newEdge = *GetGraphEdge(graph, inputEdges[0]);
  newEdge->topLevelEdge     = GetVAIndex_EdgeCGW_T(graph->edges, newEdge);

  CIEdgeT  *prevEdge = newEdge;

  //  Move all the edges in inputEdges to a list of raw edges.

  for (uint32 ei=0; ei<inputEdges.size(); ei++){
    prevEdge->nextRawEdge = inputEdges[ei];

    CIEdgeT   *edge   = GetGraphEdge(graph, inputEdges[ei]);

    edge->prevALink = NULLINDEX;
    edge->nextALink = NULLINDEX;
    edge->prevBLink = NULLINDEX;
    edge->nextBLink = NULLINDEX;


    if (isOverlapEdge(edge))
      overlapEdge    = edge;
    else
      nonOverlapEdge = edge;

    newEdge->flags.all |= edge->flags.all;  //  New merged edge has all flags the raw edges have.

    edge->topLevelEdge = newEdge->topLevelEdge;
    edge->nextRawEdge  = NULLINDEX;

    prevEdge = edge;
  }

  newEdge->flags.bits.isRaw = FALSE;

  //  Massage the status of the new merged edge.  If it is a pure status, leave it as is.
  //
  EdgeStatus status = GetEdgeStatus(newEdge);

  switch (status) {
    case UNKNOWN_EDGE_STATUS:
    case UNTRUSTED_EDGE_STATUS:
    case TENTATIVE_UNTRUSTED_EDGE_STATUS:
    case TENTATIVE_TRUSTED_EDGE_STATUS:
    case TRUSTED_EDGE_STATUS:
    case LARGE_VARIANCE_EDGE_STATUS:
    case INTER_SCAFFOLD_EDGE_STATUS:
      break;

    default:
      if      (status & UNTRUSTED_EDGE_STATUS)
        SetGraphEdgeStatus(graph, newEdge, UNTRUSTED_EDGE_STATUS);
      else if (status & TRUSTED_EDGE_STATUS)
        SetGraphEdgeStatus(graph, newEdge, TRUSTED_EDGE_STATUS);
      else
        SetGraphEdgeStatus(graph, newEdge, UNKNOWN_EDGE_STATUS);
      break;
  }

  //  We can have at most one overlap edge
  //  NOT AT ALL WHAT THIS TESTS.
  //assert((overlapEdge == NULL) || (overlapEdge && nonOverlapEdge));

  if ((overlapEdge != NULL) &&
      (inputEdges.size() == 2) &&
      (ConfirmOverlap(graph, overlapEdge, nonOverlapEdge) == FALSE))
    newEdge->flags.bits.isPossibleChimera = TRUE;

  newEdge->flags.bits.mustOverlap = TRUE;

  newEdge->edgesContributing = inputEdges.size();

  newEdge->minDistance = distance.mean - 3.0 * sqrt(distance.variance);
  newEdge->distance    = distance;

  if ((newEdge->flags.bits.isProbablyBogus == TRUE) && (isOverlapEdge(newEdge) == TRUE))
    //  Not bogus if an overlap is present.
    newEdge->flags.bits.isProbablyBogus = FALSE;
}

/* CompareChi2Scores sorts Chi squared scores in ascending order
   but also places active and passed scores before others */

#ifdef C_SORT

static int CompareChi2Scores (const void *c1, const void *c2){
  ClusterScoreChi2T *s1 = *(ClusterScoreChi2T **)c1;
  ClusterScoreChi2T *s2 = *(ClusterScoreChi2T **)c2;

  if(s1->active != s2->active)
    return((s1->active != 0) ? -1 : 1);
  if(s1->passed != s2->passed)
    return((s1->passed != 0) ? -1 : 1);
  return((s1->score < s2->score) ? -1 : 1);
}

#else

bool
operator<(ClusterScoreChi2T &s1, ClusterScoreChi2T &s2) {
  if(s1.active != s2.active)
    return(s1.active != 0);
  if(s1.passed != s2.passed)
    return(s1.passed != 0);
  return (s1.score < s2.score);
}

#endif


/* MergeGraphEdges
   Input:
   VA_TYPE(CDS_CID_t) *inputEdges is the array of references to CIEdges
   that link the same two elements (unitigs/contigs/etc) with the
   same orientation.
   GraphT *graph is the graph containing the edges and is
   needed to detect possible chimeric edges where an overlap edge
   and a mate edge are the only edges and they both depend on the
   same fragment which may be chimeric and so the edge is not
   being independently supported by two kinds of evidence.

   The edges have associated distance records which include mean and variance
   statistics. This routine merges edges which when merged pass a Chi Squared
   test indicating that they are all consistent with a new combined estimate
   of the mean and variance. First all edges are merged and Chi Squared tested.
   If this fails each edge in turn is left out of the merged set and the set is
   Chi Squared tested with the lowest passing set chosen as the merged set and
   the leftover edge remaining as a singleton. If this fails an agglomerative
   (bottom-up) nearest neighbor clustering is performed using a pairwise Chi
   Squared criteria. Each candidate cluster although chosen based on a pairwise
   Chi Squared criterion must pass the full Chi Squared test to be accepted.
   The resulting clusters are output as the merged edges. To compute the best
   least squares estimate of the mean and variance of a merged edge the
   formulation from Chapter 14 of Numerical Recipes in C p523 is used. See
   also pages 66-71 of Data Reduction and Error Analysis for the Physical
   Sciences by Bevington.

   The merged CIEdgeTs reference the raw CIEdges they incorporate by a singly
   linked list via the nextRawEdge field.  The merged edges are marked as not
   raw.  The merged edges are APPENDED to the edges array. The client of this
   routine also needs to know which raw edges were not merged with any other
   edge and therefore need to be retained. MergeEdges will reset the
   VA_TYPE(CDS_CID_t) *inputEdges variable array in order to return the indices
   of the merged edges in addition to the indices of these unmerged edges in
   this variable array.

   The return value is the number of edges that this routine appended to the
   edges array - returns 0 when no merging was done and -1 on failure.

   The client of this routine must INSERT the resulting edges into the
   appropriate graph.
*/

int
MergeGraphEdges(GraphCGW_T        *graph,
                vector<CDS_CID_t> &inputEdges){

  CIEdgeT *newEdge, *edge;
  int numEdges = inputEdges.size();
  int numEdgesAdded = 0;
  int numMergedEdges = 0;
  CDS_CID_t edgeIndex;
  int confirmable;
  CIEdgeT *overlapEdge;
  LengthT distance;
  double chiSquareScore;
  Chi2ComputeT *edgeChi2ComputePtr;

  Chi2ComputeT *edgeChi2Compute = (Chi2ComputeT *)safe_malloc(numEdges * sizeof(*edgeChi2Compute));

  /* Create an array of Chi2ComputeT records one per edge which contain
     a distance record, the inverseVariance = 1 / distance.variance and
     the weightedMean = distance.mean / distance.variance which are used
     in determining the least squares estimator for the mean and variance
     of merged edges and for Chi Squared tests */
  for (overlapEdge = NULL, edgeChi2ComputePtr = edgeChi2Compute, edgeIndex = 0;
       edgeIndex < numEdges; edgeIndex++, edgeChi2ComputePtr++){
    CDS_CID_t eid = inputEdges[edgeIndex];
    edge = GetGraphEdge(graph, eid);
    edge->topLevelEdge = eid; // TEST
    edgeChi2ComputePtr->distance = edge->distance;
    edgeChi2ComputePtr->inverseVariance =
      1./edgeChi2ComputePtr->distance.variance;
    edgeChi2ComputePtr->weightedMean = edgeChi2ComputePtr->distance.mean *
      edgeChi2ComputePtr->inverseVariance;
    /* Check if edge is an overlap edge */
    if(isOverlapEdge(edge)){
      if(overlapEdge != NULL){
        /* It is a violation to have more than one overlap edge in the
           inputEdges so we abort when this occurs */
        fprintf(stderr,"* Found 2 overlap edges out of %d\n",
                numEdges);
        PrintGraphEdge(stderr,graph," e  ", edge, edge->idA);
        PrintGraphEdge(stderr,graph," o  ", overlapEdge, overlapEdge->idA);
        fprintf(stderr,"*Aborting edge mate merge \n");
        safe_free(edgeChi2Compute);
        return(-1);
      }else{
        overlapEdge = edge;
      }
    }
  }
#if 0
  /* This code is used if we don't want to confirm a mate-link edge by
     merging it with a repeat or tandem Overlap edge */
  if(overlapEdge){
    confirmable = (numEdges > 2 ||
                   ((numEdges == 2) && (overlapEdge == NULL)) ||
                   ((numEdges == 2) && overlapEdge->flags.bits.hasContributingOverlap));
  }else{
    confirmable = numEdges >= 2;
  }
#else
  confirmable = numEdges >= 2;
#endif

  /* Check to see if the entire set of inputEdges passes the Chi Squared
     test. */
  if(confirmable && ComputeChiSquared(edgeChi2Compute, numEdges, (CDS_CID_t)(numEdges + 1),
                                      &distance, &chiSquareScore)){
    /* We passed */
    newEdge = GetFreeGraphEdge(graph);
    InitializeMergedEdge(newEdge, overlapEdge, distance, graph, inputEdges);
    inputEdges.clear();

    //    edgeIndex = GetNumCIEdgeTs(edges);  // The index of the last edge is one less than the number of edges
    //    AppendGraphEdge(graph, &newEdge);

    edgeIndex = GetVAIndex_EdgeCGW_T(graph->edges, newEdge);
    inputEdges.push_back(edgeIndex);
    numMergedEdges = 1;
    numEdgesAdded = 1;
    safe_free(edgeChi2Compute);
    return(numMergedEdges);
  }
  if(numEdges > 2){
    /* Check to see if any set of the inputEdges of size numEdges - 1 passes
       the Chi Squared Test and if so merge the set with the best score. */
    double minScore;
    CDS_CID_t skipEdgeIndex;
    Chi2ResultT  *edgeChi2ResultPtr;
    vector<CDS_CID_t> clusterEdges;
    Chi2ResultT * edgeChi2Result = (Chi2ResultT *)safe_malloc(numEdges * sizeof(*edgeChi2Result));
    for (skipEdgeIndex = -1, minScore = FLT_MAX,
           edgeChi2ResultPtr = edgeChi2Result, edgeIndex = 0;
         edgeIndex < numEdges; edgeIndex++, edgeChi2ResultPtr++){
      if((edgeChi2ResultPtr->passed =
          ComputeChiSquared(edgeChi2Compute, numEdges, edgeIndex,
                            &(edgeChi2ResultPtr->distance),
                            &(edgeChi2ResultPtr->score))) != FALSE){
        if(edgeChi2ResultPtr->score < minScore){
          minScore = edgeChi2ResultPtr->score;
          skipEdgeIndex = edgeIndex;
        }
      }
    }
    if(skipEdgeIndex >= 0){
      //  One of the edge sets passed and skipped is set to the edge to be left out of the merged
      //  set.

      CDS_CID_t skipped;

      edgeChi2ResultPtr = edgeChi2Result + skipEdgeIndex;

      for (edgeIndex = 0; edgeIndex < numEdges; edgeIndex++){
        if(edgeIndex == skipEdgeIndex){
          skipped = inputEdges[edgeIndex];
        }else{
          clusterEdges.push_back(inputEdges[edgeIndex]);
        }
      }
      newEdge = GetFreeGraphEdge(graph);
      InitializeMergedEdge(newEdge, overlapEdge, edgeChi2ResultPtr->distance, graph, clusterEdges);

      inputEdges.clear();
      inputEdges.push_back(GetVAIndex_EdgeCGW_T(graph->edges, newEdge));
      inputEdges.push_back(skipped);

      numMergedEdges = 1;
      numEdgesAdded = 2;

      safe_free(edgeChi2Result);
      safe_free(edgeChi2Compute);

      return(numMergedEdges);
    }

    safe_free(edgeChi2Result);
  }

  if(numEdges > 3){
    /* Check to see if any clusters can be formed which pass the Chi Squared
       Test by combining pairs of edges/clusters which have the best pairwise
       Chi Squared Test value and also pass the full Chi Squared Test. This
       is done in a bottom up fashion until all merges which pass have been
       performed. */
    vector<CDS_CID_t> tmpInputEdges;
    vector<CDS_CID_t> clusterEdges;

    int rowIndex, colIndex;
    int numPairs = ((numEdges - 1) * numEdges) / 2;
    int numInCluster;
    int numClusters;
    /* edgeClusterChi2 is an array of numEdges elements which represent
       each initial edge (cluster of size 1) or merged cluster which has
       this edge as the smallest index based on the 0 - (numEdges - 1)
       indexing of the inputEdges array. Each of these elements points
       to an array of ((numEdges - 1) - index) of pairwise Chi Squared
       scores which taken together are the upper right triangle of the
       numEdges x numEdges matrix of cluster versus cluster. */
    ClusterChi2T  *edgeClusterChi2Ptr = NULL;
    /* pairClusterScoreChi2 is an array of numPairs = (numEdges * (numEdges -
       1)) / 2 elements of pairwise Chi Squared scores which taken together
       are the upper right triangle of the numEdges x numEdges matrix of
       cluster versus cluster. */
    ClusterScoreChi2T *pairClusterScoreChi2Ptr = NULL;
    /* sortClusterScoreChi2 is an array for sorting the pairwise Chi Squared
       scores. */
    ClusterScoreChi2T **sortClusterScoreChi2Ptr = NULL;
    /* clusterChi2Compute is an array which must be filled in for each cluster
       in order to perform the full Chi Squared Test for that cluster. */
    Chi2ComputeT *clusterChi2ComputePtr;
    ClusterScoreChi2T *pairClusterScoreChi2 = (ClusterScoreChi2T *)safe_malloc(numPairs * sizeof(*pairClusterScoreChi2));
    ClusterScoreChi2T **sortClusterScoreChi2 = (ClusterScoreChi2T **)safe_malloc(numPairs * sizeof(*sortClusterScoreChi2));
    ClusterChi2T *edgeClusterChi2 = (ClusterChi2T *)safe_malloc(numEdges * sizeof(*edgeClusterChi2));
    Chi2ComputeT *clusterChi2Compute = (Chi2ComputeT *)safe_malloc(numEdges * sizeof(*clusterChi2Compute));

    /* Initialize the arrays. */
    for(edgeChi2ComputePtr = edgeChi2Compute,
          pairClusterScoreChi2Ptr = pairClusterScoreChi2,
          sortClusterScoreChi2Ptr = sortClusterScoreChi2,
          edgeClusterChi2Ptr = edgeClusterChi2, rowIndex = 0;
        rowIndex < numEdges;
        rowIndex++, edgeClusterChi2Ptr++, edgeChi2ComputePtr++){
      /* Set the pointer to the pairwise Chi Squared scores for this row. */
      edgeClusterChi2Ptr->pairwiseScores = pairClusterScoreChi2Ptr;
      /* replacedBy points to the cluster which contains this edge or has in
         turn been replaced by another cluster in which case the pointer
         must be followed until it is NULL. Currently only used to indicate
         that this edge is in another cluster. */
      edgeClusterChi2Ptr->replacedBy = NULL;
      /* replaced is a linked list of additional edges which are in this
         cluster. */
      edgeClusterChi2Ptr->replaced = NULL;
      edgeClusterChi2Ptr->numInCluster = 1;
      edgeClusterChi2Ptr->distance = edgeChi2ComputePtr->distance;
      for(colIndex = rowIndex + 1; colIndex < numEdges;
          colIndex++, pairClusterScoreChi2Ptr++, sortClusterScoreChi2Ptr++){
        *sortClusterScoreChi2Ptr = pairClusterScoreChi2Ptr;
        /* rowIndex and colIndex are the cordinates in the pairwise Chi
           Squared matrix. rowIndex < colIndex because we only store the
           upper right triangle of the matrix. */
        pairClusterScoreChi2Ptr->rowIndex = rowIndex;
        pairClusterScoreChi2Ptr->colIndex = colIndex;
        /* active indicates that the clusters for this row and column have
           not been replaced by previous merge steps. When clusters are
           merged only the lower index cluster remains active by absorbing
           the higher index cluster. */
        pairClusterScoreChi2Ptr->active = TRUE;
        /* passed indicates that this pair of clusters exceeded the threshold
           for passing the pairwise Chi Squared Test. */
        pairClusterScoreChi2Ptr->passed =
          PairwiseChiSquare((double)edgeClusterChi2Ptr->distance.mean,
                            edgeClusterChi2Ptr->distance.variance,
                            (double)(edgeChi2Compute + colIndex)->distance.mean,
                            (edgeChi2Compute + colIndex)->distance.variance,
                            &(pairClusterScoreChi2Ptr->distance),
                            &(pairClusterScoreChi2Ptr->score),
                            (double)PAIRWISECHI2THRESHOLD_CGW);
        /* We want all of the clusters to try to pass the full chi squared test so we set this to TRUE */
        pairClusterScoreChi2Ptr->passed = TRUE;
      }
    }
    /* Sort potential merge candidates by their pairwise Chi Squared scores. */

#ifdef C_SORT
    qsort((void *)sortClusterScoreChi2, numPairs, sizeof(*sortClusterScoreChi2), CompareChi2Scores);
#else
    std::sort(sortClusterScoreChi2, sortClusterScoreChi2 + numPairs);
#endif

    sortClusterScoreChi2Ptr = sortClusterScoreChi2;
    pairClusterScoreChi2Ptr = *sortClusterScoreChi2Ptr;
    /* While there are still potential merge candidates try the next merge. */
    while(pairClusterScoreChi2Ptr->passed &&
          pairClusterScoreChi2Ptr->active){
      ClusterChi2T *lastRowClusterChi2Ptr = NULL;
      /* Fill in the array for the full Chi Squared Test for the row
         cluster. */
      for(edgeClusterChi2Ptr = edgeClusterChi2 +
            pairClusterScoreChi2Ptr->rowIndex,
            numInCluster = edgeClusterChi2Ptr->numInCluster,
            clusterChi2ComputePtr = clusterChi2Compute;
          edgeClusterChi2Ptr !=  NULL;
          edgeClusterChi2Ptr = edgeClusterChi2Ptr->replaced,
            clusterChi2ComputePtr++ ){
        /* Copy this structure from already computed structure at beginning
           of this routine for each edge in the row cluster. */
        *clusterChi2ComputePtr =
          edgeChi2Compute[edgeClusterChi2Ptr - edgeClusterChi2];
        /* Save this pointer at end of linked list so that if we merge these
           two clusters we can join the two linked lists. */
        lastRowClusterChi2Ptr = edgeClusterChi2Ptr;
      }
      /* Fill in the array for the full Chi Squared Test for the column
         cluster. */
      for(edgeClusterChi2Ptr = edgeClusterChi2 +
            pairClusterScoreChi2Ptr->colIndex,
            numInCluster += edgeClusterChi2Ptr->numInCluster;
          edgeClusterChi2Ptr !=  NULL;
          edgeClusterChi2Ptr = edgeClusterChi2Ptr->replaced,
            clusterChi2ComputePtr++ ){
        /* Copy this structure from already computed structure at beginning
           of this routine for each edge in the column cluster. */
        *clusterChi2ComputePtr =
          edgeChi2Compute[edgeClusterChi2Ptr - edgeClusterChi2];
      }
      /* Run the full Chi Squared Test. */
      if(ComputeChiSquared(clusterChi2Compute, numInCluster,
                           (CDS_CID_t)(numInCluster + 1), &distance, &chiSquareScore)){
        /* We passed so we need to update the clusters by merging the
           row and column clusters and computing new pairwise Chi Squared
           scores for this new cluster. */
        ClusterChi2T *rowClusterChi2Ptr = edgeClusterChi2 +
          pairClusterScoreChi2Ptr->rowIndex;
        ClusterChi2T *colClusterChi2Ptr = edgeClusterChi2 +
          pairClusterScoreChi2Ptr->colIndex;
        ClusterScoreChi2T *rowClusterScoreChi2Ptr;
        /* Merge the column cluster into the row cluster. */
        rowClusterChi2Ptr->numInCluster += colClusterChi2Ptr->numInCluster;
        rowClusterChi2Ptr->distance = distance;
        lastRowClusterChi2Ptr->replaced = colClusterChi2Ptr;
        colClusterChi2Ptr->replacedBy = rowClusterChi2Ptr;
        /* Deactivate pairwise Chi Squared scores for the column cluster.
           First the column of scores and ... */
        for(edgeClusterChi2Ptr = edgeClusterChi2,
              colIndex = pairClusterScoreChi2Ptr->colIndex - 1;
            colIndex >= 0 ; edgeClusterChi2Ptr++, colIndex--){
          (edgeClusterChi2Ptr->pairwiseScores + colIndex)->active = FALSE;
        }
        /* then the row of scores. */
        for(rowClusterScoreChi2Ptr = colClusterChi2Ptr->pairwiseScores,
              colIndex = pairClusterScoreChi2Ptr->colIndex + 1;
            colIndex < numEdges; colIndex++, rowClusterScoreChi2Ptr++){
          rowClusterScoreChi2Ptr->active = FALSE;
        }
        /* Now compute the new pairwise Chi Squared scores the merged
           cluster. First the column of scores and ... */
        for(edgeClusterChi2Ptr = edgeClusterChi2,
              colIndex = pairClusterScoreChi2Ptr->rowIndex - 1;
            colIndex >= 0; edgeClusterChi2Ptr++, colIndex--){
          rowClusterScoreChi2Ptr = edgeClusterChi2Ptr->pairwiseScores +
            colIndex;
          if(rowClusterScoreChi2Ptr->active){
            rowClusterScoreChi2Ptr->passed =
              PairwiseChiSquare((double)edgeClusterChi2Ptr->distance.mean,
                                edgeClusterChi2Ptr->distance.variance,
                                (double)rowClusterChi2Ptr->distance.mean,
                                rowClusterChi2Ptr->distance.variance,
                                &(rowClusterScoreChi2Ptr->distance),
                                &(rowClusterScoreChi2Ptr->score),
                                (double)PAIRWISECHI2THRESHOLD_CGW);
            /* We want all of the clusters to try to pass the full chi squared test so we set this to TRUE */
            rowClusterScoreChi2Ptr->passed = TRUE;
          }
        }
        /* then the row of scores. */
        for(rowClusterScoreChi2Ptr = rowClusterChi2Ptr->pairwiseScores,
              colIndex = pairClusterScoreChi2Ptr->rowIndex + 1;
            colIndex < numEdges; colIndex++, rowClusterScoreChi2Ptr++){
          edgeClusterChi2Ptr = edgeClusterChi2 + colIndex;
          if(rowClusterScoreChi2Ptr->active){
            rowClusterScoreChi2Ptr->passed =
              PairwiseChiSquare((double)edgeClusterChi2Ptr->distance.mean,
                                edgeClusterChi2Ptr->distance.variance,
                                (double)rowClusterChi2Ptr->distance.mean,
                                rowClusterChi2Ptr->distance.variance,
                                &(rowClusterScoreChi2Ptr->distance),
                                &(rowClusterScoreChi2Ptr->score),
                                (double)PAIRWISECHI2THRESHOLD_CGW);
            /* We want all of the clusters to try to pass the full chi squared test so we set this to TRUE */
            rowClusterScoreChi2Ptr->passed = TRUE;
          }
        }
        /* Resort the pairwise Chi Squared scores. */

#ifdef C_SORT
        qsort((void *)sortClusterScoreChi2, numPairs, sizeof(*sortClusterScoreChi2), CompareChi2Scores);
#else
        std::sort(sortClusterScoreChi2, sortClusterScoreChi2 + numPairs);
#endif

        sortClusterScoreChi2Ptr = sortClusterScoreChi2;
      }else{
        /* Failed the full Chi Squared Test so see if there are any more
           potential merges and try again. */
        pairClusterScoreChi2Ptr->passed = FALSE;
        if(sortClusterScoreChi2Ptr < (sortClusterScoreChi2 + (numPairs - 1))){
          sortClusterScoreChi2Ptr++;
        }else{
          break;
        }
      }
      pairClusterScoreChi2Ptr = *sortClusterScoreChi2Ptr;
    }
    /* For each cluster either compute the merged edge or if there is only
       one edge in the cluster return it as an unmerged edge. */
    for(numClusters = 0, edgeClusterChi2Ptr = edgeClusterChi2, rowIndex = 0;
        rowIndex < numEdges; rowIndex++, edgeClusterChi2Ptr++){
      if(edgeClusterChi2Ptr->replacedBy == NULL){
        ClusterChi2T *groupClusterChi2Ptr;
        /* Put the edge indices for a cluster into clusterEdges. */
        clusterEdges.clear();
        for(groupClusterChi2Ptr = edgeClusterChi2Ptr, numInCluster = 0;
            groupClusterChi2Ptr != NULL;
            groupClusterChi2Ptr = groupClusterChi2Ptr->replaced,
              numInCluster++){
          clusterEdges.push_back(inputEdges[groupClusterChi2Ptr - edgeClusterChi2]);
        }
        assert(numInCluster == edgeClusterChi2Ptr->numInCluster);
        if(numInCluster > 1){
          newEdge = GetFreeGraphEdge(graph);
          InitializeMergedEdge(newEdge, overlapEdge, edgeClusterChi2Ptr->distance, graph, clusterEdges);

          edgeIndex = GetVAIndex_EdgeCGW_T(graph->edges, newEdge);

          tmpInputEdges.push_back(edgeIndex);

          numMergedEdges++;
          numEdgesAdded++;
        }else{
          tmpInputEdges.push_back(clusterEdges[0]);
          numEdgesAdded++;
        }
        numClusters++;
      }
    }

    inputEdges.clear();

    for(uint32 ei = 0; ei < tmpInputEdges.size(); ei++){
      inputEdges.push_back(tmpInputEdges[ei]);
    }
    edge = GetGraphEdge(graph, inputEdges[0]);

    if (GlobalData->verbose) {
      fprintf(stderr,"**** Couldn't merge these \n");
      fprintf(stderr,"**** MORE THAN ONE (%d)  EDGE BETWEEN "F_CID" and "F_CID"\n",
              numClusters, edge->idA, edge->idB);
      switch(graph->type){
        case CI_GRAPH:
          DumpChunkInstance(stderr,ScaffoldGraph,GetGraphNode(graph, edge->idA), FALSE, FALSE, FALSE, FALSE);
          DumpChunkInstance(stderr,ScaffoldGraph,GetGraphNode(graph, edge->idB), FALSE, FALSE, FALSE, FALSE);
          break;
        case CONTIG_GRAPH:
          DumpContig(stderr,ScaffoldGraph,GetGraphNode(graph, edge->idA), FALSE);
          DumpContig(stderr,ScaffoldGraph,GetGraphNode(graph, edge->idB), FALSE);
          break;
        case SCAFFOLD_GRAPH:
          DumpCIScaffold(stderr,ScaffoldGraph,GetGraphNode(graph, edge->idA), FALSE);
          DumpCIScaffold(stderr,ScaffoldGraph,GetGraphNode(graph, edge->idB), FALSE);
          break;
        default:
          break;
      }

      for (uint32 ei=0; ei<inputEdges.size(); ei++) {
        edge = GetGraphEdge(graph, inputEdges[ei]);
        PrintGraphEdge(stderr, graph," *  ", edge, edge->idA);
      }
    }

    safe_free(pairClusterScoreChi2);
    safe_free(sortClusterScoreChi2);
    safe_free(edgeClusterChi2);
    safe_free(clusterChi2Compute);
    safe_free(edgeChi2Compute);

    return(numMergedEdges);
  }

  /* If we reached here none of the inputEdges was mergable so just return them. */

  if (GlobalData->verbose) {
    edge = GetGraphEdge(graph, inputEdges[0]);

    fprintf(stderr,"**** Couldn't merge these \n");
    fprintf(stderr,"**** MORE THAN ONE (%d)  EDGE BETWEEN "F_CID" and "F_CID"\n",
            numEdges, edge->idA, edge->idB);
  }

  numEdgesAdded = numEdges;

  if (GlobalData->verbose) {
    for(int32 ei = 0; ei < numEdges; ei++) {
      edge = GetGraphEdge(graph, inputEdges[ei]);
      PrintGraphEdge(stderr, graph," *  ", edge, edge->idA);
    }
  }

  safe_free(edgeChi2Compute);

  return(numMergedEdges);
}

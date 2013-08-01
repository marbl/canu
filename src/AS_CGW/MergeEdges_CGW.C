
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
 * MERCHANTABILITY or FITNESS FOR A PATICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/
static char *rcsid = "$Id: MergeEdges_CGW.c,v 1.40 2012-11-13 19:08:33 brianwalenz Exp $";

#include "AS_global.H"
#include "AS_CGW_dataTypes.H"
#include "AS_UTL_interval.H"
#include "AS_UTL_Var.H"
#include "UtilsREZ.H"
#include "Globals_CGW.H"
#include "ScaffoldGraph_CGW.H"
#include "ChiSquareTest_CGW.H"

#include <vector>
#include <algorithm>

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

  //fprintf(stderr,"* chunk = "F_CID" endB = %d minOffset = "F_S32" maxoffset = "F_S32" overlap = "F_S32"\n",
  //        chunkID, endB, minOffset, maxOffset,overlap);

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
ComputeChiSquared(vector<Chi2ComputeT> &edges,
                  int                   numEdges,
                  CDS_CID_t             skip,
                  LengthT              &distance,
                  double               &score) {
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

  distance.mean     = leastSquaresMean;
  distance.variance = 1.0 / cumInverseVariance;

  for (int32 edgeIndex = 0; edgeIndex < numEdges; edgeIndex++) {
    if (edgeIndex != skip)
      cumScore += (((edges[edgeIndex].distance.mean - leastSquaresMean) * 
                    (edges[edgeIndex].distance.mean - leastSquaresMean)) /
                   edges[edgeIndex].distance.variance);
  }

  score = cumScore;

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

bool
sortScoreCompare(ClusterScoreChi2T *s1, ClusterScoreChi2T *s2) {
  if (s1->active != s2->active)  return(s1->active == true);
  if (s1->passed != s2->passed)  return(s1->passed == true);

  return(s1->score < s2->score);
}


//  MergeGraphEdges

//  VA_TYPE(CDS_CID_t) *inputEdges is the array of references to CIEdges that link the same two
//  elements (unitigs/contigs/etc) with the same orientation.
//
//  GraphT *graph is the graph containing the edges and is needed to detect possible chimeric edges
//  where an overlap edge and a mate edge are the only edges and they both depend on the same
//  fragment which may be chimeric and so the edge is not being independently supported by two kinds
//  of evidence.
//
//  The edges have associated distance records which include mean and variance statistics. This
//  routine merges edges which when merged pass a Chi Squared test indicating that they are all
//  consistent with a new combined estimate of the mean and variance.
//
//  All edges are merged and Chi Squared tested.
//
//  If this fails each edge in turn is left out of the merged set and the set is Chi Squared tested
//  with the lowest passing set chosen as the merged set and the leftover edge remaining as a
//  singleton.
//
//  If this fails an agglomerative (bottom-up) nearest neighbor clustering is performed using a
//  pairwise Chi Squared criteria. Each candidate cluster although chosen based on a pairwise Chi
//  Squared criterion must pass the full Chi Squared test to be accepted.  The resulting clusters
//  are output as the merged edges.
//
//  To compute the best least squares estimate of the mean and variance of a merged edge the
//  formulation from Chapter 14 of Numerical Recipes in C p523 is used. See also pages 66-71 of Data
//  Reduction and Error Analysis for the Physical Sciences by Bevington.
//
//  The return value is the number of edges that this routine appended to the edges array - returns
//  0 when no merging was done and -1 on failure.
//
//  The client of this routine must INSERT the resulting edges into the appropriate graph.

bool
MergeGraphEdges_All(GraphCGW_T            *graph,
                    vector<CDS_CID_t>     &inputEdges,
                    vector<Chi2ComputeT>  &edgeCompute,
                    CIEdgeT               *overlapEdge) {

  assert(inputEdges.size() >= 2);

  for (int32 ei=0; ei<inputEdges.size(); ei++) {
    CIEdgeT *edge                   = GetGraphEdge(graph, inputEdges[ei]);
    edge->topLevelEdge              = inputEdges[ei];
    edgeCompute[ei].distance        = edge->distance;
    edgeCompute[ei].inverseVariance = 1.0 / edgeCompute[ei].distance.variance;
    edgeCompute[ei].weightedMean    = edgeCompute[ei].distance.mean * edgeCompute[ei].inverseVariance;
  }

  bool     confirmable = (inputEdges.size() >= 2);
  LengthT  distance;
  double   chiSquareScore;

#if 0
  //  If we don't want to confirm a mate-link edge by
  //  merging it with a repeat or tandem Overlap edge
  if (overlapEdge)
    confirmable = ((inputEdges.size() > 2) ||
                   ((inputEdges.size() == 2) && (overlapEdge == NULL)) ||
                   ((inputEdges.size() == 2) && (overlapEdge->flags.bits.hasContributingOverlap == true)));
#endif

  //  Check to see if the entire set of inputEdges passes the Chi Squared test.

  if ((confirmable == true) &&
      (ComputeChiSquared(edgeCompute, inputEdges.size(), inputEdges.size() + 1, distance, chiSquareScore) == true)) {
    CIEdgeT   *ne  = GetFreeGraphEdge(graph);

    InitializeMergedEdge(ne, overlapEdge, distance, graph, inputEdges);

    inputEdges.clear();
    inputEdges.push_back(GetVAIndex_EdgeCGW_T(graph->edges, ne));

    return(true);
  }

  return(false);
}



bool
MergeGraphEdges_AllButOne(GraphCGW_T            *graph,
                          vector<CDS_CID_t>     &inputEdges,
                          vector<Chi2ComputeT>  &edgeCompute,
                          CIEdgeT               *overlapEdge) {

  assert(inputEdges.size() >= 3);

  vector<Chi2ResultT>   edgeChi2Result(inputEdges.size());

  CDS_CID_t minIndex = -1;
  double    minScore = DBL_MAX;

  for (int32 ei=0; ei<inputEdges.size(); ei++) {
    edgeChi2Result[ei].passed = ComputeChiSquared(edgeCompute, inputEdges.size(), ei, edgeChi2Result[ei].distance, edgeChi2Result[ei].score); 

    if ((edgeChi2Result[ei].passed == true) &&
        (edgeChi2Result[ei].score < minScore)) {
      minScore = edgeChi2Result[ei].score;
      minIndex = ei;
    }
  }

  //  If we found a result, minIndex is > -1 and is the edge we left out of the merged set.
  //  Initialize the new edge with all but the skipped edge.  Result[ei] is the distance estimate
  //  with that edge skipped.

  if (minIndex == -1)
    return(false);

  vector<CDS_CID_t> clusterEdges;

  for (int32 ei=0; ei<inputEdges.size(); ei++)
    if (ei != minIndex)
      clusterEdges.push_back(inputEdges[ei]);

  CIEdgeT    *ne = GetFreeGraphEdge(graph);
  CDS_CID_t   ni = GetVAIndex_EdgeCGW_T(graph->edges, ne);
  CDS_CID_t   ex = inputEdges[minIndex];

  InitializeMergedEdge(ne, overlapEdge, edgeChi2Result[minIndex].distance, graph, clusterEdges);

  inputEdges.clear();
  inputEdges.push_back(ni);
  inputEdges.push_back(ex);

  return(true);
}




//  Starting with the best pairwise set, add one edge as long as it is consistent.
//  If no edges are consistent, seed another best pairwise set.
//


//  A pair of edges, with ChiSquared score and distance.
class mgePair {
public:
  double          score;     //  ChiSquared score
  LengthT         distance;  //  estimated distance
  uint32          aIdx;      //  index into inputEdges vector for one of the edges
  uint32          bIdx;      //  index into inputEdges vector for the other edge

  bool operator<(const mgePair &that) const {
    return(score < that.score);
  };
};

//  A set of edges, with ChiSquared score and distance.
class mgeSet {
public:
  mgeSet(mgePair &pair_) {
    score    = pair_.score;
    distance = pair_.distance;

    edges.push_back(pair_.aIdx);
    edges.push_back(pair_.bIdx);
  };

  mgeSet(CDS_CID_t edge_) {
    score             = 0.0;
    distance.mean     = 0.0;
    distance.variance = 0.0;

    edges.push_back(edge_);
  };

  ~mgeSet() {
  };

  void   insertPair(const mgePair &pair_, double score_, double mean_, double variance_) {
    score              = score_;
    distance.mean      = mean_;
    distance.variance  = variance_;

    edges.push_back(pair_.aIdx);
    edges.push_back(pair_.bIdx);
  }

  void   insertSingle(CDS_CID_t edgeid_, double score_, double mean_, double variance_) {
    score              = score_;
    distance.mean      = mean_;
    distance.variance  = variance_;

    edges.push_back(edgeid_);
  }

  double             score;     //  the ChiSquared score
  LengthT            distance;  //  estimate of true size
  vector<CDS_CID_t>  edges;     //  edges involved
};




bool
MergeGraphEdges_Greedy(GraphCGW_T            *graph,
                       vector<CDS_CID_t>     &inputEdges,
                       vector<Chi2ComputeT>  &edgeCompute,
                       CIEdgeT               *overlapEdge) {

  assert(inputEdges.size() >= 4);

  set<CDS_CID_t>           edgesMerged;  //  These edges have been processed (indices into inputEdges)
  vector<mgeSet>           mergedEdges;  //  A list of the merged edges we are building

  //  Compute ChiSquared for all pairs of edges.

  uint32           pairsLen = 0;
  uint32           pairsMax = inputEdges.size() * (inputEdges.size() - 1) / 2;
  vector<mgePair>  pairs(pairsMax);

  for (uint32 a=0; a<inputEdges.size(); a++) {
    for (uint32 b=a+1; b<inputEdges.size(); b++) {
      pairs[pairsLen].score             = 0.0;
      pairs[pairsLen].distance.mean     = 0.0;
      pairs[pairsLen].distance.variance = 0.0;
      pairs[pairsLen].aIdx              = a;
      pairs[pairsLen].bIdx              = b;

      if (PairwiseChiSquare(edgeCompute[a].distance.mean, edgeCompute[a].distance.variance,
                            edgeCompute[b].distance.mean, edgeCompute[b].distance.variance,
                            &pairs[pairsLen].distance,
                            &pairs[pairsLen].score,
                            PAIRWISECHI2THRESHOLD_500) == false)
        pairs[pairsLen].score = DBL_MAX;

      pairsLen++;
    }
  }

  assert(pairsLen == pairsMax);

  sort(pairs.begin(), pairs.end());

  //  Seed the first cluster with the lowest scoring pair....unless there is no valid lowest scoring pair.

  pairsLen = 0;

  if (pairs[pairsLen].score >= PAIRWISECHI2THRESHOLD_500)
    return(false);

  mergedEdges.push_back(mgeSet(pairs[pairsLen]));
  edgesMerged.insert(pairs[pairsLen].aIdx);
  edgesMerged.insert(pairs[pairsLen].bIdx);

  //  While there are unmerged edges, add one to its best matching cluster.
  //
  //  Iterate over the list of sorted pairs.  If a pair is consistent with an existing merged edge, add it to that merged edge.
  //  Otherwise, 

  for (pairsLen=0; pairsLen < pairsMax; pairsLen++) {

    //  Skip over pairs that have already been merged, or pairs that are not consistent.
    while ((pairsLen < pairsMax) &&
           ((pairs[pairsLen].score >= PAIRWISECHI2THRESHOLD_500) ||  //  Probably should use a 'passed' flag...
            (edgesMerged.count(pairs[pairsLen].aIdx) == 1) ||
            (edgesMerged.count(pairs[pairsLen].bIdx) == 1)))
      pairsLen++;

    //  Just get out of here if there are no more pairs left.
    if (pairsLen == pairsMax)
      break;

    //  Try to find a place to put paired-edges pairs[pairsLen], either in an existing cluster, or its own new cluster.
    //  We'll test all clusters, and pick the lowest score.

    uint32  minci         = UINT32_MAX;
    double  minScore      = DBL_MAX;
    LengthT minDistance   = { 0.0, 0.0 };

    vector<Chi2ComputeT>  clstCompute(inputEdges.size());

    for (uint32 ci=0; ci<mergedEdges.size(); ci++) {
      uint32  numEdges = mergedEdges[ci].edges.size() + 2;

      //  Copy the compute data for each edge in the new cluster.
      for (uint32 ii=0; ii<mergedEdges[ci].edges.size(); ii++)
        clstCompute[ii] = edgeCompute[mergedEdges[ci].edges[ii]];

      //  Add in the new edges
      clstCompute[numEdges - 2] = edgeCompute[pairs[pairsLen].aIdx];
      clstCompute[numEdges - 1] = edgeCompute[pairs[pairsLen].bIdx];

      double  newScore      = 0.0;
      LengthT newDistance   = { 0.0, 0.0 };

      //  Test if this pair is consistent with this set, and that it has the lowest score
      if (ComputeChiSquared(clstCompute, numEdges, numEdges + 1, newDistance, newScore)) {
        if (newScore < minScore) {
          minci       = ci;
          minScore    = newScore;
          minDistance = newDistance;
        }
      }
    }

    if (minci < mergedEdges.size()) {
      mergedEdges[minci].insertPair(pairs[pairsLen], minScore, minDistance.mean, minDistance.variance);

    } else {
      mergedEdges.push_back(mgeSet(pairs[pairsLen]));
    }

    //  Pair has been processed.
    edgesMerged.insert(pairs[pairsLen].aIdx);
    edgesMerged.insert(pairs[pairsLen].bIdx);
  }

  //  For all the remaining unmerged edges, try to add them one at a time to a cluster.  If they
  //  don't add, add as singletons.
  //
  //  It is unlikely these will be added; they didn't have a valid pairwise ChiSquared result (otherwise
  //  they would have been added as the pair already), and so probably won't match a larger cluster.

  for (uint32 ie=0; ie<inputEdges.size(); ie++) {
    if (edgesMerged.count(ie) == 1)
      continue;

    uint32  minci    = UINT32_MAX;
    double  minScore      = 0.0;
    LengthT minDistance   = { 0.0, 0.0 };

    for (uint32 ci=0; ci<mergedEdges.size(); ci++) {
      //  Test if this pair is consistent with this set, and that it has the lowest score
    }

    if (minci < mergedEdges.size()) {
      mergedEdges[minci].insertSingle(ie, minScore, minDistance.mean, minDistance.variance);

    } else {
      mergedEdges.push_back(mgeSet(ie));
    }
  }

  //  Convert all the indices stored in mergedEdges from indices into inputEdges to
  //  the edge index in inputEdges[].

  for (uint32 ci=0; ci<mergedEdges.size(); ci++)
    for (uint32 ei=0; ei<mergedEdges[ci].edges.size(); ei++)
      mergedEdges[ci].edges[ei] = inputEdges[mergedEdges[ci].edges[ei]];

  //  Create new edges for each cluster, and copy to inputEdges.

  inputEdges.clear();

  for (uint32 ci=0; ci<mergedEdges.size(); ci++) {
    if (mergedEdges[ci].edges.size() > 1) {
      CIEdgeT *ne = GetFreeGraphEdge(graph);

      InitializeMergedEdge(ne, overlapEdge, mergedEdges[ci].distance, graph, mergedEdges[ci].edges);

      inputEdges.push_back(GetVAIndex_EdgeCGW_T(graph->edges, ne));

      //fprintf(stderr, "MergeGraphEdges()--  greedy - edge "F_SIZE_T" with %d raw edges, %d %d, score=%f %.0f +- %.0f\n",
      //        GetVAIndex_EdgeCGW_T(graph->edges, ne),
      //        ne->edgesContributing,
      //        mergedEdges[ci].edges[0],
      //        mergedEdges[ci].edges[1],
      //        mergedEdges[ci].score,
      //        mergedEdges[ci].distance.mean, mergedEdges[ci].distance.variance);

    } else {
      //fprintf(stderr, "MergeGraphEdges()--  greedy - edge %d singleton\n",
      //        mergedEdges[ci].edges[0]);
      inputEdges.push_back(mergedEdges[ci].edges[0]);
    }
  }

  return(true);
}




int
MergeGraphEdges(GraphCGW_T        *graph,
                vector<CDS_CID_t> &inputEdges){

  if (inputEdges.size() == 1)
    return(1);


  //  Find the overlap edge, if it exists.  There can be at most one.  Abort the merge if we find such a case.

  CIEdgeT             *overlapEdge     = NULL;

  for (int32 ei=0; ei<inputEdges.size(); ei++) {
    CIEdgeT *edge = GetGraphEdge(graph, inputEdges[ei]);

    if (isOverlapEdge(edge) == false)
      continue;

    if (overlapEdge != NULL) {
      fprintf(stderr, "MergeGraphEdges()--  Found multiple overlap edges between %d and %d, removing duplicate edge %d.\n",
              edge->idA, edge->idB, inputEdges[ei]);

      //  Remove the edge from our list.
      inputEdges.erase(inputEdges.begin() + ei);

      //  Remove the edge from the graph.
      //DeleteGraphEdge(graph, edge);

      //UnlinkGraphEdge(graph, edge);
      FreeGraphEdge(graph, edge);

      //  Back up one, so we examine the new ei element.
      ei--;

      continue;
    }

    assert(overlapEdge == NULL);

    overlapEdge = edge;
  }

  uint32  numOlapEdges = 0;

  for (int32 ei=0; ei<inputEdges.size(); ei++)
    if (isOverlapEdge(GetGraphEdge(graph, inputEdges[ei])))
      numOlapEdges++;

  if (numOlapEdges >= 2) {
    for (int32 ei=0; ei<inputEdges.size(); ei++) {
      CIEdgeT *edge = GetGraphEdge(graph, inputEdges[ei]);
      fprintf(stderr, "MergeGraphEdges()--    edge %d %.0f +- %.0f overlap %d\n",
              inputEdges[ei],
              edge->distance.mean, edge->distance.variance,
              isOverlapEdge(edge));
    }
  }
  assert(numOlapEdges < 2);


  vector<Chi2ComputeT> edgeCompute(inputEdges.size());

  if (MergeGraphEdges_All(graph, inputEdges, edgeCompute, overlapEdge))
    return(1);

  if (inputEdges.size() == 2)
    return(0);

  if (MergeGraphEdges_AllButOne(graph, inputEdges, edgeCompute, overlapEdge))
    return(1);

  if (inputEdges.size() == 3)
    return(0);

  //fprintf(stderr, "MergeGraphEdges()--  "F_SIZE_T" raw edges between %d and %d - use greedy merging.\n",
  //        inputEdges.size(),
  //        GetGraphEdge(graph, inputEdges[0])->idA,
  //        GetGraphEdge(graph, inputEdges[0])->idB);

  if (MergeGraphEdges_Greedy(graph, inputEdges, edgeCompute, overlapEdge))
    return(inputEdges.size());

  //  If we reached here none of the inputEdges was mergable so just return them.

  return(inputEdges.size());
}

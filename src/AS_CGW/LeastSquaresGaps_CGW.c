
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
static char CM_ID[] = "$Id: LeastSquaresGaps_CGW.c,v 1.25 2008-05-31 06:49:46 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_interval.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "ChiSquareTest_CGW.h"
#include "ChunkOverlap_CGW.h"


#define FIXED_RECOMPUTE_NOT_ENOUGH_CLONES /* long standing bug: is it fixed yet? it seems to be */
#undef  FIXED_RECOMPUTE_NOT_ENOUGH_CLONES /* nope */

#undef NEG_GAP_VARIANCE_PROBLEM_FIXED  /* if undef'ed, allow processing to continue despite a negative gap variance */

#undef TEST_FIXUPMISORDER


#define MAX_ABSOLUTE_SLOP 10000
#define MAX_SIGMA_SLOP 10000 /* essentially infinite ... */
#define MAX_OVERLAP_SLOP_CGW 10


/* declarations for LAPACK/DXML calls to linear algebra routines */
#define FTN_INT   long int
#define F_FTN_INT    "%ld"

extern int dgemv_(char *, FTN_INT *, FTN_INT *, 
                  double *, double *, FTN_INT *, double *, FTN_INT *, 
                  double *, double *, FTN_INT *);
extern int dpbtrf_(char *, FTN_INT *, FTN_INT *, double *,
                   FTN_INT *, FTN_INT *);
extern int dpbtrs_(char *, FTN_INT *, FTN_INT *, FTN_INT *, double *,
                   FTN_INT *, double *, FTN_INT *, FTN_INT *);


//  Except as noted:
//    0 nothing, 1 warnings, 2 lots of stuff
//
//  The names end in LV for "level" to distinguish them from the
//  function that is usually a very similar name.  This would be
//  incredibly useful if someone wants to, say, merge two classes into
//  one.
//
//  If we ever create a file of serious warnings, look for "//assert" for
//  things to put into it.
//
typedef struct {
  int   fixUpMisorderedContigsLV;
  FILE *fixUpMisorderedContigsFP;

  int   checkInternalEdgeStatusLV;

  int   fixedRecomputeSingluarLV;      //  dump scaffolds AND assert if this bug appears again

  int   markInternalEdgeStatusLV;
  FILE *markInternalEdgeStatusFP;

  int   recomputeOffsetsVerboseLV;      //  gory detail about the gap size compute
  int   recomputeOffsetsLV;             //  about the recomputeOffsetsInScaffold process
  int   leastSquaresGapsLV;             //  info about the gap proces, contig containment, etc

} debugflags_t;

static debugflags_t  debug = {0, 0L, 0, 0, 0, 0L, 0, 1, 0};




//   If we find a pair of contigs that are misordered, we rip out
//   thisCI, move it to the right place based on the implied position
//   in the overlapEdge
//
int
FixUpMisorderedContigs(CIScaffoldT           *scaffold,
                       ContigT               *prevCI,
                       ContigT               *thisCI, 
                       ChunkOrientationType   edgeOrient, 
                       double                 inferredMean,
                       double                 inferredVariance,
                       EdgeCGW_T             *overlapEdge){
  ChunkOrientationType newEdgeOrient = GetEdgeOrientationWRT(overlapEdge, prevCI->id);

  LengthT aEndOffset, bEndOffset;

  //  Same orientation?  I guess we want to try to merge these contigs then...
  //
  if (edgeOrient == newEdgeOrient)
    return(0);

  // 1 == prevCI
  // 2 == thisCI
  // edgeOrient is expected orientation of (1,2)
  // newEdgeOrient is correct orientation of (1,2) 

  DumpContig(stderr,ScaffoldGraph, prevCI, FALSE);
  DumpContig(stderr,ScaffoldGraph, thisCI, FALSE);
  PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, " overlapEdge: ", overlapEdge, overlapEdge->idA);

  fprintf(stderr,"* edgeOrient %c   edge->orient = %c  newEdgeOrient = %c  prevCI = " F_CID "   thisCI = " F_CID " mean:%g\n",
          edgeOrient, overlapEdge->orient, newEdgeOrient,
          prevCI->id, thisCI->id, inferredMean);



  aEndOffset.mean = aEndOffset.variance = -1.0;
  bEndOffset.mean = bEndOffset.variance = -1.0;


  switch (edgeOrient) {
    case AB_AB:  //  aka AS_NORMAL
      assert(newEdgeOrient == BA_BA);
      //           expected                                 actual
      //       ---------1-------->                            ---------1-------->
      //                  ---------2-------->     ---------2-------->    
      //                                                      |<=====|
      //                                                      overlap length

      // overlap is negative
      bEndOffset.mean = prevCI->offsetAEnd.mean -  overlapEdge->distance.mean; 
      bEndOffset.variance = prevCI->offsetAEnd.variance +  overlapEdge->distance.variance;
      aEndOffset.mean = bEndOffset.mean - thisCI->bpLength.mean;
      aEndOffset.variance = bEndOffset.variance - thisCI->bpLength.variance;
      break;
    case AB_BA:  //  aka AS_INNIE
      assert(newEdgeOrient == BA_AB);
      //           expected                                 actual
      //       ---------1-------->                            ---------1-------->
      //                  <--------2--------     <--------2--------    
      //                                                      |<=====|
      //                                                      overlap length


      // overlap is negative
      aEndOffset.mean = prevCI->offsetAEnd.mean -  overlapEdge->distance.mean; 
      aEndOffset.variance = prevCI->offsetAEnd.variance + overlapEdge->distance.variance;
      bEndOffset.mean = aEndOffset.mean - thisCI->bpLength.mean;
      bEndOffset.variance = aEndOffset.variance - thisCI->bpLength.variance;
      break;
    case BA_AB:  //  aka AS_OUTTIE
      assert(newEdgeOrient == AB_BA);
      //           expected                                    actual
      //       <---------1--------                            <---------1--------
      //                  --------2-------->     --------2------->    
      //                                                      |<=====|
      //                                                      overlap length

      // overlap is negative!
      bEndOffset.mean = prevCI->offsetBEnd.mean -  overlapEdge->distance.mean; 
      bEndOffset.variance = prevCI->offsetBEnd.variance -  overlapEdge->distance.variance;
      aEndOffset.mean = bEndOffset.mean - thisCI->bpLength.mean;
      aEndOffset.variance = bEndOffset.variance - thisCI->bpLength.variance;
      break;
    case BA_BA:  //  aka AS_ANTI
      assert(newEdgeOrient == AB_AB);
      //           expected                                 actual
      //       <---------1--------                            <---------1--------
      //                  <---------2--------     <---------2--------    
      //                                                      |<=====|
      //                                                      overlap length

      // overlap is negative
      aEndOffset.mean = prevCI->offsetBEnd.mean - overlapEdge->distance.mean; 
      aEndOffset.variance = prevCI->offsetBEnd.variance +  overlapEdge->distance.variance;
      bEndOffset.mean = aEndOffset.mean - thisCI->bpLength.mean;
      bEndOffset.variance = aEndOffset.variance - thisCI->bpLength.variance;
      break;
    default:
      return(0);
      break;
  }

  fprintf(stderr,"* Overlap is (" F_CID "," F_CID ",%c)  moving " F_CID " from (%g,%g) to (%g,%g)\n",
          overlapEdge->idA,
          overlapEdge->idB,
          overlapEdge->orient,
          thisCI->id,
          thisCI->offsetAEnd.mean,
          thisCI->offsetBEnd.mean,
          aEndOffset.mean,
          bEndOffset.mean);

  thisCI->offsetAEnd = aEndOffset;
  thisCI->offsetBEnd = bEndOffset;

  RemoveCIFromScaffold(ScaffoldGraph, scaffold, thisCI, FALSE);
  InsertCIInScaffold(ScaffoldGraph,
                     thisCI->id, scaffold->id, aEndOffset, bEndOffset,
                     TRUE, FALSE);

  return(1);
}



EdgeCGW_T *FindOverlapEdgeChiSquare(ScaffoldGraphT *graph,
                                    NodeCGW_T *sourceCI,
                                    CDS_CID_t targetId,
                                    ChunkOrientationType edgeOrient,
                                    double inferredMean,
                                    double inferredVariance,
                                    float *chiSquaredValue,
                                    float chiSquareThreshold,
                                    int *alternate, int verbose){
  GraphEdgeIterator edges;
  EdgeCGW_T *edge;
  EdgeCGW_T *bestEdge = (EdgeCGW_T *)NULL;
  int end;
  float bestChiSquaredValue = FLT_MAX;

  *alternate = FALSE;
  if((edgeOrient == AB_AB) || (edgeOrient == AB_BA)){
    /* edgeOrient == AB_XX */
    end = B_END;
  }else{
    /* edgeOrient == BA_XX */
    end = A_END;
  }
  InitGraphEdgeIterator(ScaffoldGraph->RezGraph, sourceCI->id, end, ALL_EDGES,
                        GRAPH_EDGE_RAW_ONLY, &edges);// Use raw edges
  while((edge = NextGraphEdgeIterator(&edges))!= NULL){
    CDS_CID_t otherCID = (edge->idA == sourceCI->id) ? edge->idB : edge->idA;
    if((otherCID == targetId) && isOverlapEdge(edge) &&
       !isContainmentEdge(edge) ){// deal with these later
      if(GetEdgeOrientationWRT(edge, sourceCI->id) == edgeOrient){
        if(PairwiseChiSquare((float)inferredMean, (float)inferredVariance,
                             (float)edge->distance.mean,
                             (float)((MAX_OVERLAP_SLOP_CGW * MAX_OVERLAP_SLOP_CGW) / 9),
                             (LengthT *)NULL, chiSquaredValue,
                             (float)chiSquareThreshold)){
          if(bestEdge == (EdgeCGW_T *)NULL ||
             (*chiSquaredValue < bestChiSquaredValue)){
            bestEdge = edge;
            bestChiSquaredValue = *chiSquaredValue;
          }
        }
      }
    }
  }

  if (bestEdge != NULL)
    return(bestEdge);

  {
    CDS_COORD_t minOverlap, maxOverlap;
    minOverlap = MAX(CGW_MISSED_OVERLAP,
                     -(inferredMean + (3.0 * sqrt(inferredVariance))));
    maxOverlap = -(inferredMean - (3.0 * sqrt(inferredVariance)));
    if(maxOverlap >= CGW_MISSED_OVERLAP){
      float effectiveOlap;
      ChunkOverlapCheckT olap;
      assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
      olap = OverlapChunks(graph->RezGraph,
                           sourceCI->id, targetId,    // handles suspicious
                           edgeOrient, minOverlap, maxOverlap,
                           AS_CGW_ERROR_RATE, FALSE);
      effectiveOlap = -olap.overlap;
      if(olap.suspicious){
        fprintf(stderr,"* FOEXS: SUSPICIOUS Overlap found! Looked for (" F_CID "," F_CID ",%c)[" F_COORD "," F_COORD "] found (" F_CID "," F_CID ",%c) " F_COORD "\n",
                sourceCI->id, targetId, edgeOrient,
                minOverlap, maxOverlap,
                olap.spec.cidA, olap.spec.cidB,
                olap.spec.orientation, olap.overlap);
        effectiveOlap = -(GetGraphNode(ScaffoldGraph->ContigGraph, targetId)->bpLength.mean +
                          sourceCI->bpLength.mean - olap.overlap);
      }
      if(olap.overlap){
        CDS_CID_t edgeIndex;
        edgeIndex = InsertComputedOverlapEdge(graph->RezGraph, &olap);
        edge = GetGraphEdge(graph->RezGraph, edgeIndex);

        // Create an appropriate hash table entry
        CreateChunkOverlapFromEdge(ScaffoldGraph->RezGraph, edge, FALSE);

        if(PairwiseChiSquare((float)inferredMean, (float)inferredVariance,
                             effectiveOlap,
                             (float)((MAX_OVERLAP_SLOP_CGW * MAX_OVERLAP_SLOP_CGW) / 9),
                             (LengthT *)NULL, chiSquaredValue,
                             (float)chiSquareThreshold)){
          *alternate = olap.suspicious;
          return(edge);
        }else{
          fprintf(stderr,"* Failed pairwise test between (%g, %g) and (%g,%g) not returning edge (" F_CID "," F_CID ",%c) %g\n",
                  inferredMean, inferredVariance, effectiveOlap, (float) ((MAX_OVERLAP_SLOP_CGW * MAX_OVERLAP_SLOP_CGW) / 9),
                  edge->idA, edge->idB, edge->orient, edge->distance.mean);
        }
      }
    }
    return((EdgeCGW_T *)NULL);
  }
}

void CheckInternalEdgeStatus(ScaffoldGraphT *graph, CIScaffoldT *scaffold, 
                             float pairwiseChiSquaredThreshhold,
                             float maxVariance,
                             int doNotChange, int verbose){
  CIScaffoldTIterator CIs;
  /* Iterate over all of the CIEdges */
  GraphEdgeIterator edges;
  EdgeCGW_T *edge;
  NodeCGW_T *thisCI;
  int32 numCIs;
  int32 indexCIs;

  if (debug.checkInternalEdgeStatusLV > 1)
    fprintf(stderr, "Checking Edges for Scaffold " F_CID "\n", scaffold->id);

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  for(indexCIs = 0;
      (thisCI = NextCIScaffoldTIterator(&CIs)) != NULL;){
    thisCI->indexInScaffold = indexCIs;
    indexCIs++;
  }
  numCIs = indexCIs;
  if(numCIs != scaffold->info.Scaffold.numElements){
    if (debug.checkInternalEdgeStatusLV > 1)
      fprintf(stderr, "NumElements inconsistent %d,%d\n", numCIs, scaffold->info.Scaffold.numElements);
    scaffold->info.Scaffold.numElements = numCIs;
  }

  assert(indexCIs == numCIs);

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL){

    // Use merged edges
    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, ALL_END,
                          ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);

    while((edge = NextGraphEdgeIterator(&edges))!= NULL){
      int isA = (edge->idA == thisCI->id);
      NodeCGW_T *otherCI =
        GetGraphNode(ScaffoldGraph->RezGraph,
                     (isA? edge->idB: edge->idA));
      ChunkOrientationType edgeOrient;
      FragOrient thisCIorient, otherCIorient;
      LengthT gapDistance;
      float chiSquareResult;

      /* We do not want to check edges with certain
         labels as specified in the doNotChange mask. */
      if(edge->flags.bits.edgeStatus & doNotChange){
        continue;
      }
      /* Only edges between CIs in the same scaffold should be trusted. */
      if(otherCI->scaffoldID != thisCI->scaffoldID){
        if (debug.checkInternalEdgeStatusLV > 0){
          EdgeStatus edgeStatus = GetEdgeStatus(edge);
          if ((edgeStatus == TRUSTED_EDGE_STATUS) ||
              (edgeStatus == TENTATIVE_TRUSTED_EDGE_STATUS))
            fprintf(stderr, "[" F_CID "." F_CID "," F_CID "." F_CID "]Trusted edge really interscaffold edge.\n",
                    thisCI->id, thisCI->scaffoldID, otherCI->id,
                    otherCI->scaffoldID);
          else if (edgeStatus != INTER_SCAFFOLD_EDGE_STATUS)
            fprintf(stderr, "[" F_CID "." F_CID "," F_CID "." F_CID "]Interscaffold edge marked as %d.\n",
                    thisCI->id, thisCI->scaffoldID, otherCI->id,
                    otherCI->scaffoldID, edgeStatus);
        }
        continue;
      }
      /* We only want to check an edge once and the
         iterator will visit intra scaffold edges twice so we use
         the condition that the edge is canonical - that is the
         condition thisCI->info.CI.indexInScaffold < otherCI->info.CI.indexInScaffold
         to guarantee labelling the edge only once. */
      if(otherCI->indexInScaffold <= thisCI->indexInScaffold){
        continue;
      }

      /* Check that the edge orientation is consistent with the CI positions
         and orientations within the scaffold. */
      edgeOrient = GetEdgeOrientationWRT(edge, thisCI->id);
      switch(edgeOrient){
        case AB_BA:
          //      thisCI                                        otherCI
          //  A --------------------- B               B --------------------- A
          //    5'----->                                           <------5'
          thisCIorient = A_B;
          otherCIorient = B_A;
          gapDistance.mean = otherCI->offsetBEnd.mean - thisCI->offsetBEnd.mean;
          gapDistance.variance = otherCI->offsetBEnd.variance -
            thisCI->offsetBEnd.variance;
          break;
        case AB_AB:
          //      thisCI                                        otherCI
          //  A --------------------- B               A --------------------- B
          //    5'----->                                           <------5'
          thisCIorient = A_B;
          otherCIorient = A_B;
          gapDistance.mean = otherCI->offsetAEnd.mean - thisCI->offsetBEnd.mean;
          gapDistance.variance = otherCI->offsetAEnd.variance -
            thisCI->offsetBEnd.variance;
          break;
        case BA_BA:
          //      thisCI                                        otherCI
          //  B --------------------- A               B --------------------- A
          //    5'----->                                           <------5'
          thisCIorient = B_A;
          otherCIorient = B_A;
          gapDistance.mean = otherCI->offsetBEnd.mean - thisCI->offsetAEnd.mean;
          gapDistance.variance = otherCI->offsetBEnd.variance -
            thisCI->offsetAEnd.variance;
          break;
        case BA_AB:
          //      thisCI                                        otherCI
          //  B --------------------- A               A --------------------- B
          //    5'----->                                           <------5'
          thisCIorient = B_A;
          otherCIorient = A_B;
          gapDistance.mean = otherCI->offsetAEnd.mean - thisCI->offsetAEnd.mean;
          gapDistance.variance = otherCI->offsetAEnd.variance -
            thisCI->offsetAEnd.variance;
          break;
        default:
          assert(0);
          break;
      }
      if ((GetNodeOrient(thisCI) != thisCIorient) ||
          (GetNodeOrient(otherCI) != otherCIorient)){
        /* Mark as untrusted an edge whose orientation does not agree
           with the orientation of the CIs in the scaffold. */
        if (debug.checkInternalEdgeStatusLV > 0) {
          EdgeStatus edgeStatus = GetEdgeStatus(edge);
          if ((edgeStatus == TRUSTED_EDGE_STATUS) ||
              (edgeStatus == TENTATIVE_TRUSTED_EDGE_STATUS))
            fprintf(stderr, "[" F_CID "." F_CID "," F_CID "." F_CID "]Trusted edge really Bad orientation (%c,%c) (%c,%c).\n",
                    thisCI->id, thisCI->scaffoldID, otherCI->id,
                    otherCI->scaffoldID, GetNodeOrient(thisCI), thisCIorient,
                    GetNodeOrient(otherCI), otherCIorient);
          else if ((edgeStatus != UNTRUSTED_EDGE_STATUS) &&
                   (edgeStatus != TENTATIVE_UNTRUSTED_EDGE_STATUS))
            fprintf(stderr, "[" F_CID "." F_CID "," F_CID "." F_CID "]Bad orientation (%c,%c) (%c,%c) edge marked as %d.\n",
                    thisCI->id, thisCI->scaffoldID, otherCI->id,
                    otherCI->scaffoldID, GetNodeOrient(thisCI), thisCIorient,
                    GetNodeOrient(otherCI), otherCIorient, edgeStatus);
        }
        continue;
      }
      if(gapDistance.variance <= 0.0){
        if (debug.checkInternalEdgeStatusLV > 0) {
          EdgeStatus edgeStatus = GetEdgeStatus(edge);
          fprintf(stderr, "[" F_CID "." F_CID "," F_CID "." F_CID "]Bad Gap Variance (%f,%f) (%f,%f) edge marked as %d.\n",
                  thisCI->id, thisCI->scaffoldID, otherCI->id,
                  otherCI->scaffoldID,
                  gapDistance.mean, gapDistance.variance,
                  edge->distance.mean, edge->distance.variance, edgeStatus);
        }
        if (debug.checkInternalEdgeStatusLV > 1) {
          DumpACIScaffoldNew(stderr,ScaffoldGraph,scaffold,TRUE);
          DumpACIScaffoldNew(stderr,ScaffoldGraph,scaffold,FALSE);
        }
      }else if(!PairwiseChiSquare((float)gapDistance.mean,
                                  gapDistance.variance,
                                  (float)edge->distance.mean,
                                  edge->distance.variance,
                                  (LengthT *)NULL, &chiSquareResult,
                                  pairwiseChiSquaredThreshhold)){
        /* Mark  this edge as untrusted if the distance of the edge is not
           consistent with the estimated gap distance as judged by the
           Chi Squared Test. */
        if (debug.checkInternalEdgeStatusLV > 0) {
          EdgeStatus edgeStatus = GetEdgeStatus(edge);
          if ((edgeStatus == TRUSTED_EDGE_STATUS) ||
              (edgeStatus == TENTATIVE_TRUSTED_EDGE_STATUS))
            fprintf(stderr, "[" F_CID "." F_CID "," F_CID "." F_CID "]Trusted edge really Bad Chi Squared %f (%f,%f) (%f,%f).\n",
                    thisCI->id, thisCI->scaffoldID, otherCI->id,
                    otherCI->scaffoldID,
                    chiSquareResult, gapDistance.mean, gapDistance.variance,
                    edge->distance.mean, edge->distance.variance);
          else if ((edgeStatus != UNTRUSTED_EDGE_STATUS) &&
                   (edgeStatus != TENTATIVE_UNTRUSTED_EDGE_STATUS))
            fprintf(stderr, "[" F_CID "." F_CID "," F_CID "." F_CID "]Bad Chi Squared %f (%f,%f) (%f,%f) edge marked as %d.\n",
                    thisCI->id, thisCI->scaffoldID, otherCI->id,
                    otherCI->scaffoldID,
                    chiSquareResult, gapDistance.mean, gapDistance.variance,
                    edge->distance.mean, edge->distance.variance, edgeStatus);
        }
        continue;
      }
      if(edge->distance.variance > maxVariance){
        if (debug.checkInternalEdgeStatusLV > 0) {
          EdgeStatus edgeStatus = GetEdgeStatus(edge);
          if ((edgeStatus == TRUSTED_EDGE_STATUS) ||
              (edgeStatus == TENTATIVE_TRUSTED_EDGE_STATUS))
            fprintf(stderr, "[" F_CID "." F_CID "," F_CID "." F_CID "]Trusted edge really Variance too large %f.\n",
                    thisCI->id, thisCI->scaffoldID, otherCI->id,
                    otherCI->scaffoldID, edge->distance.variance);
          else if (edgeStatus != LARGE_VARIANCE_EDGE_STATUS)
            fprintf(stderr, "[" F_CID "." F_CID "," F_CID "." F_CID "]Variance too large %f edge marked as %d.\n",
                    thisCI->id, thisCI->scaffoldID, otherCI->id,
                    otherCI->scaffoldID, edge->distance.variance, edgeStatus);
        }
        continue;
      }

      if (debug.checkInternalEdgeStatusLV > 1) {
        EdgeStatus edgeStatus = GetEdgeStatus(edge);
        if ((edgeStatus != TRUSTED_EDGE_STATUS) &&
            (edgeStatus != TENTATIVE_TRUSTED_EDGE_STATUS)){
          fprintf(stderr, "[" F_CID "." F_CID "," F_CID "." F_CID "]Edge marked as %d should be trusted.\n",
                  thisCI->id, thisCI->scaffoldID, otherCI->id,
                  otherCI->scaffoldID, edgeStatus);
          fprintf(stderr, " - Good Chi Squared %f (%f,%f) (%f,%f)\n",
                  chiSquareResult, gapDistance.mean, gapDistance.variance,
                  edge->distance.mean, edge->distance.variance);
        }
      }
    }
  }
  return;
}





//  Dump trusted/raw edges
void
dumpTrustedEdges(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold, int32 edgeTypes) {
  CIEdgeT * edge;
  ChunkInstanceT * chunk;
  GraphEdgeIterator   edges;
  CIScaffoldTIterator CIs;
  int set = 0;

  fprintf(stderr, "TRUSTED_EDGE in scaffold %d\n", scaffold->id);

  InitCIScaffoldTIterator(sgraph, scaffold, TRUE, FALSE, &CIs);

  while ((chunk = NextCIScaffoldTIterator(&CIs)) != NULL) {
    //  We don't have a setID if we're in
    //  BuildScaffoldsFromFirstPriniciples() -> RepeatRez() ->
    //  TidyUpScaffolds() -> LeastSquaresGapEstimates() -> here
    //
    //assert(chunk->setID >= 0);
    InitGraphEdgeIterator(sgraph->RezGraph, chunk->id, 
                          ALL_END, edgeTypes, // ALL_TRUSTED_EDGES, 
                          GRAPH_EDGE_DEFAULT, //GRAPH_EDGE_CONFIRMED_ONLY,
                          &edges);
    while ((edge = NextGraphEdgeIterator(&edges)) != NULL) {

      //
      // get the other end
      //
      ChunkInstanceT * otherChunk = GetGraphNode(sgraph->RezGraph,
                                                 (chunk->id == edge->idA) ?
                                                 edge->idB : edge->idA);

      int32 weight = edge->edgesContributing - (isOverlapEdge(edge));

      assert(otherChunk != NULL);

      // See each edge only once
      if(chunk->id != edge->idA)
        continue;

      if(edge->flags.bits.isBridge){
        fprintf(stderr,"* WARNING: chunk " F_CID " weight = %d bridge edge\n",
                chunk->id, weight);
        PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph,
                       "Bridge ", edge, chunk->id);
#ifdef DEBUG
        EdgeCGW_T *e;
        GraphEdgeIterator Edges;
        InitGraphEdgeIterator(sgraph->ContigGraph,chunk->id,ALL_END,
                              ALL_TRUSTED_EDGES,GRAPH_EDGE_DEFAULT,&Edges);
        fprintf(stderr,"Edges out from " F_CID ":\n",chunk->id); 
        while(NULL!= (e = NextGraphEdgeIterator(&Edges)))
          PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph,
                         "DEBUG Bridge ",e, chunk->id);	
#endif
      }
      if(isSingletonOverlapEdge(edge) ||
         (weight == 1 && edge->flags.bits.isBridge))
        continue;

      //
      // if the other end is not in this scaffold
      // ignore it
      //
      if (chunk->scaffoldID != otherChunk->scaffoldID)
        continue;

      fprintf(stderr, "TRUSTED_EDGE "F_CID"("F_CID") - "F_CID"("F_CID") weight=%d\n",
              chunk->id, chunk->setID,
              otherChunk->id, otherChunk->setID,
              weight);
    }
  }
}





typedef struct {
  int numCIs;
  int numClones;
  int numGaps;
  size_t sizeofLengthCIs;
  size_t sizeofCloneGapStart;
  size_t sizeofCloneGapEnd;
  size_t sizeofGapConstants;
  size_t sizeofGapCoefficients;
  size_t sizeofGapVariance;
  size_t sizeofCloneVariance;
  size_t sizeofCloneMean;
  size_t sizeofSpannedGaps;
  size_t sizeofGapSize;
  size_t sizeofGapSizeVariance;
  size_t sizeofGapsToComputeGaps;
  size_t sizeofComputeGapsToGaps;
  LengthT *lengthCIs;
  CDS_COORD_t *cloneGapStart;
  CDS_COORD_t *cloneGapEnd;
  double *gapConstants;
  double *gapCoefficients;
  double *gapVariance;
  double *cloneVariance;
  double *cloneMean;
  double *spannedGaps;
  double *gapSize;
  double *gapSizeVariance;
  int *gapsToComputeGaps;
  int *computeGapsToGaps;
} RecomputeData;

void freeRecomputeData(RecomputeData *data){
  safe_free(data->lengthCIs);
  safe_free(data->gapConstants);
  safe_free(data->gapCoefficients);
  safe_free(data->gapVariance);
  safe_free(data->cloneGapStart);
  safe_free(data->cloneGapEnd);
  safe_free(data->cloneVariance);
  safe_free(data->cloneMean);
  safe_free(data->spannedGaps);
  safe_free(data->gapSize);
  safe_free(data->gapSizeVariance);
  safe_free(data->gapsToComputeGaps);
  safe_free(data->computeGapsToGaps);
}

void ReportRecomputeData(RecomputeData *data, FILE *stream){
  size_t totalMemorySize = 0;
  totalMemorySize = 
    data->sizeofLengthCIs +
    data->sizeofCloneGapStart +
    data->sizeofCloneGapEnd +
    data->sizeofGapConstants +
    data->sizeofGapCoefficients +
    data->sizeofGapVariance +
    data->sizeofCloneVariance +
    data->sizeofCloneMean +
    data->sizeofSpannedGaps +
    data->sizeofGapSize +
    data->sizeofGapSizeVariance +
    data->sizeofGapsToComputeGaps +
    data->sizeofComputeGapsToGaps;

  if(totalMemorySize > 1<<30) // if > 1GB
    fprintf(stream, "* Recompute Offsets CIs:%d Clones:%d Gaps:%d allocated " F_SIZE_T " bytes\n",
            data->numCIs,
            data->numClones,
            data->numGaps,
            totalMemorySize);
}

RecomputeOffsetsStatus RecomputeOffsetsInScaffold(ScaffoldGraphT *graph,
                                                  CIScaffoldT *scaffold,
                                                  int allowOrderChanges,
                                                  int forceNonOverlaps,
                                                  int verbose){

  RecomputeData data;
  CIScaffoldTIterator CIs;
  /* Iterate over all of the "trusted" CIEdges */
  GraphEdgeIterator edges;
  EdgeCGW_T *edge;
  NodeCGW_T *thisCI, *prevCI;
  int32 numCIs;
  int32 indexCIs;
  int standardEdgeStatusFails =0;

  int numGaps, numComputeGaps;
  LengthT *lengthCIs, *lengthCIsPtr;
  int maxDiagonals = 1;
  int numClones = 0;
  CDS_CID_t indexClones;
  CDS_COORD_t *cloneGapStart, *cloneGapEnd;
  int *gapsToComputeGaps, *computeGapsToGaps;
  double *gapCoefficients, *gapConstants;
  double *gapVariance, *cloneVariance;
  double *cloneMean;
  double *spannedGaps;
  double *gapSize, *gapSizeVariance;
  double squaredError;
  LengthT *maxOffset = NULL;
  int hardConstraintSet;

  data.lengthCIs = NULL;
  data.cloneGapStart = NULL;
  data.cloneGapEnd = NULL;
  data.gapConstants = NULL;
  data.gapCoefficients = NULL;
  data.gapVariance = NULL;
  data.cloneVariance = NULL;
  data.cloneMean = NULL;
  data.spannedGaps = NULL;
  data.gapSize = NULL;
  data.gapSizeVariance = NULL;
  data.gapsToComputeGaps = NULL;
  data.computeGapsToGaps = NULL;

  StartTimerT(&GlobalData->RecomputeOffsetsTimer);   /*  START */

  CheckInternalEdgeStatus(graph, scaffold, PAIRWISECHI2THRESHOLD_CGW, 100000000000.0, 0, FALSE);

  if(IsScaffoldInternallyConnected(ScaffoldGraph,scaffold,ALL_TRUSTED_EDGES)!=1){
    standardEdgeStatusFails = 1;
    if (debug.recomputeOffsetsLV > 0) {
      fprintf(stderr, "RecomputeOffsetsInScaffold()- WARNING: scaffold " F_CID " is not internally connected using\n", scaffold->id);
      fprintf(stderr, "                              ALL_TRUSTED_EDGES will proceed with edge set determined by \n");
      fprintf(stderr, "                              IsInternalEdgeStatusVaguelyOK instead of PairwiseChiSquare test\n");
    }
    //assert(0);
  }

  numCIs = scaffold->info.Scaffold.numElements;
  numGaps = numCIs - 1;
  if(numGaps < 1){
    freeRecomputeData(&data);
    return (RECOMPUTE_NO_GAPS);
  }
  data.numCIs = numCIs;
  data.sizeofLengthCIs = numCIs * sizeof(*lengthCIs);
  data.lengthCIs = lengthCIs = (LengthT *)safe_malloc(data.sizeofLengthCIs);

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  for(indexCIs = 0, lengthCIsPtr = lengthCIs;
      (thisCI = NextCIScaffoldTIterator(&CIs)) != NULL;){
    thisCI->indexInScaffold = indexCIs;
    *lengthCIsPtr = thisCI->bpLength;

    if (debug.recomputeOffsetsVerboseLV > 1)
      fprintf(stderr, "Length of CI %d," F_CID " %f\n",
              indexCIs, thisCI->id, lengthCIsPtr->mean);

    indexCIs++;
    lengthCIsPtr++;
  }
  assert(indexCIs == numCIs);

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL){

    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, 
                          ALL_END,
                          (standardEdgeStatusFails ? ALL_INTERNAL_EDGES : ALL_TRUSTED_EDGES),
                          GRAPH_EDGE_RAW_ONLY,
                          //        GRAPH_EDGE_RAW_ONLY | (standardEdgeStatusFails ? GRAPH_EDGE_VERBOSE : 0),
                          &edges);// ONLY RAW

    while((edge = NextGraphEdgeIterator(&edges))!= NULL){
      int isA = (edge->idA == thisCI->id);
      NodeCGW_T *otherCI = GetGraphNode(ScaffoldGraph->RezGraph, (isA ? edge->idB : edge->idA));

      if (otherCI->scaffoldID == -1) {
        if (debug.recomputeOffsetsLV > 1) {
          fprintf(stderr, "RecomputeOffsetsInScaffold()--  WARNING!  otherCI is in scaffold -1!\n");
          fprintf(stderr, "   edge:  idA:%d idB:%d\n", edge->idA, edge->idB);
          fprintf(stderr, " thisCI:  scaffoldID:%d is idA:%d\n", thisCI->scaffoldID, edge->idA == thisCI->id);
          fprintf(stderr, "otherCI:  scaffoldID:%d is idA:%d\n", otherCI->scaffoldID, edge->idA == otherCI->id);
        }
        continue;
      }

      // RAW EDGES ONLY
      assert(edge->flags.bits.isRaw);

      if(otherCI->indexInScaffold <= thisCI->indexInScaffold){
        continue; // Only interested in looking at an edge once
      }


      // the following is invoked if we are looking at ALL_EDGES
      // rather than ALL_TRUSTED_EDGES; this occurs as a
      // desperate attempt to resurrect some scaffold edges that
      // will reconnect a scaffold that is otherwise
      // disconnected and results in a singularity; however, if
      // this edge is really really bad, we need to throw it
      // away ...
      //
      if(standardEdgeStatusFails && !IsInternalEdgeStatusVaguelyOK(edge,thisCI->id)){
        continue;
      }


      numClones++;
      if (debug.recomputeOffsetsVerboseLV > 1)
        fprintf(stderr,"RecomputeOffsets: adding clone between %d and %d\n",
                thisCI->id,otherCI->id);

      if((otherCI->indexInScaffold - thisCI->indexInScaffold) >
         maxDiagonals){
        maxDiagonals = otherCI->indexInScaffold - thisCI->indexInScaffold;

        if (debug.recomputeOffsetsVerboseLV > 1) {
          fprintf(stderr, "Max Diagonals %d (%d,%d) [" F_CID "." F_CID "," F_CID "." F_CID "]\n",
                  maxDiagonals, thisCI->indexInScaffold,
                  otherCI->indexInScaffold, thisCI->scaffoldID,
                  thisCI->id, otherCI->scaffoldID, otherCI->id);
          fprintf(stderr, "Max Diagonals %d (%d,%d) [" F_CID "." F_CID "," F_CID "." F_CID "]\n",
                  maxDiagonals, thisCI->indexInScaffold,
                  otherCI->indexInScaffold, thisCI->scaffoldID,
                  thisCI->id, otherCI->scaffoldID, otherCI->id);
        }
      }
    }
  }
  if(numClones < numGaps){
    freeRecomputeData(&data);
#ifdef  FIXED_RECOMPUTE_NOT_ENOUGH_CLONES
    assert(0 /* Not enough clones */);
#endif
    return (RECOMPUTE_NOT_ENOUGH_CLONES);
  }
  {
    double *gapEnd, *gapPtr;

    data.sizeofCloneGapStart = numClones * sizeof(*cloneGapStart);
    data.cloneGapStart = cloneGapStart = (int32 *)safe_malloc(data.sizeofCloneGapStart);

    data.sizeofCloneGapEnd = numClones * sizeof(*cloneGapEnd);
    data.cloneGapEnd = cloneGapEnd = (int32 *)safe_malloc(data.sizeofCloneGapEnd);

    data.sizeofGapConstants = numGaps * sizeof(*gapConstants);
    data.gapConstants = gapConstants = (double *)safe_malloc(data.sizeofGapConstants);

    for(gapPtr = gapConstants, gapEnd = gapPtr + numGaps; gapPtr < gapEnd; gapPtr++){
      *gapPtr = 0.0;
    }
    data.sizeofGapCoefficients = (maxDiagonals * numGaps) * sizeof(*gapCoefficients);
    data.gapCoefficients = gapCoefficients = (double *)safe_malloc(data.sizeofGapCoefficients);

    for(gapPtr = gapCoefficients, gapEnd = gapPtr + (maxDiagonals * numGaps);
        gapPtr < gapEnd; gapPtr++){
      *gapPtr = 0.0;
    }
    data.numClones = numClones;
    data.numGaps = numGaps;
    data.sizeofCloneMean = numClones * sizeof(*cloneMean);
    data.cloneMean = cloneMean = (double *)safe_malloc(data.sizeofCloneMean);

    data.sizeofCloneVariance = numClones * sizeof(*cloneVariance);
    data.cloneVariance = cloneVariance = (double *)safe_malloc(data.sizeofCloneVariance);

    data.sizeofGapVariance = numGaps * sizeof(*gapVariance);
    data.gapVariance = gapVariance = (double *)safe_malloc(data.sizeofGapVariance);

    for(gapPtr = gapVariance, gapEnd = gapPtr + numGaps; gapPtr < gapEnd; gapPtr++){
      *gapPtr = 0.0;
    }
    data.sizeofSpannedGaps = numGaps * sizeof(*spannedGaps);
    data.spannedGaps = spannedGaps = (double *)safe_malloc(data.sizeofSpannedGaps);

    data.sizeofGapSize = numGaps * sizeof(*gapSize);
    data.gapSize = gapSize = (double *)safe_malloc(data.sizeofGapSize);

    data.sizeofGapSizeVariance = numGaps * sizeof(*gapSizeVariance);
    data.gapSizeVariance = gapSizeVariance = (double *)safe_malloc(data.sizeofGapSizeVariance);

    data.sizeofGapsToComputeGaps = numGaps * sizeof(*gapsToComputeGaps);
    data.gapsToComputeGaps = gapsToComputeGaps = (int32 *)safe_malloc(data.sizeofGapsToComputeGaps);

    data.sizeofComputeGapsToGaps = numGaps * sizeof(*computeGapsToGaps);
    data.computeGapsToGaps = computeGapsToGaps = (int32 *)safe_malloc(data.sizeofComputeGapsToGaps);

    for(numComputeGaps = 0; numComputeGaps < numGaps; numComputeGaps++){
      gapsToComputeGaps[numComputeGaps] =
        computeGapsToGaps[numComputeGaps] = numComputeGaps;
    }
    ReportRecomputeData(&data, stderr);
  }
  /* The following code solves a set of linear equations in order to
     find a least squares minimal solution for the length of the gaps
     between CIs within this scaffold. The squared error to be
     minimized is defined to be the expected size of the gaps spanned
     by a clone minus the gap sizes we are solving for spanned by the
     clone squared divided by the variance of the expected size summed
     over all clones. The expected size of the gaps spanned by a clone
     is computed as follows: the expected/mean length of the clone is
     provided as an input parameter based on what DNA library the clone
     is from, from this mean length we then subtract the portions of
     the CIs which contain the clone end fragments (this step has already
     been done for us and is encoded in the edge->distance record as
     the mean - the variance is also previously computed based on the
     assumption that the two random variables are independent so that
     the variances are additive), in addition we subtract the lengths of
     CIs entirely spanned by the clone (this depends on knowing the order
     of the CIs which is provided by the scaffold) and again add the
     variances assuming independence. In order to find the least squares
     minimal solution we take the partial derivatives of the squared
     error with respect to the gap sizes we are solving for and setting
     them to zero resulting in numGaps equations with numGaps unknowns.
     We use the LAPACK tools to solve this set of equations. Note that
     each term in the squared error sum contributes to a particular
     partial derivative iff the clone for that term spans the gap for
     that partial derivative. */
  do{
    int maxClone;
    indexClones = 0;
    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

    while((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL){

      InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, 
                            ALL_END,
                            (standardEdgeStatusFails ? ALL_INTERNAL_EDGES : ALL_TRUSTED_EDGES),
                            GRAPH_EDGE_RAW_ONLY,
                            //          GRAPH_EDGE_RAW_ONLY | (standardEdgeStatusFails ? GRAPH_EDGE_VERBOSE : 0),
                            &edges);// ONLY RAW

      while((edge = NextGraphEdgeIterator(&edges))!= NULL){
        int isA = (edge->idA == thisCI->id);
        NodeCGW_T *otherCI =
          GetGraphNode(ScaffoldGraph->RezGraph,
                       (isA? edge->idB: edge->idA));
        double constant, constantVariance, inverseVariance;
        int lengthCIsIndex, gapIndex;
        int colIndex;

        // RAW EDGES ONLY
        assert(edge->flags.bits.isRaw);

        if(otherCI->indexInScaffold <= thisCI->indexInScaffold){
          continue; // Only interested in looking at an edge once
        }

        // the following is paired up with a similar test above--see comment there
        if(standardEdgeStatusFails && !IsInternalEdgeStatusVaguelyOK(edge,thisCI->id)){
          continue;
        }


        if(indexClones>=numClones){

          if (debug.recomputeOffsetsVerboseLV > 1)
            fprintf(stderr,"ROIS: Enlarging clone-dependent arrays -- must have improved layout enough to rescue some more clones\n");

          numClones*=1.2;
          data.numClones=numClones;

          data.sizeofCloneMean = numClones * sizeof(*cloneMean);
          data.cloneMean = cloneMean = (double *)safe_realloc(cloneMean,data.sizeofCloneMean);

          data.sizeofCloneVariance = numClones * sizeof(*cloneVariance);
          data.cloneVariance = cloneVariance = (double *)safe_realloc(cloneVariance,data.sizeofCloneVariance);

          data.sizeofCloneGapStart = numClones * sizeof(*cloneGapStart);
          data.cloneGapStart = cloneGapStart = (int32 *)safe_realloc(cloneGapStart,data.sizeofCloneGapStart);

          data.sizeofCloneGapEnd = numClones * sizeof(*cloneGapEnd);
          data.cloneGapEnd = cloneGapEnd = (int32 *)safe_realloc(cloneGapEnd,data.sizeofCloneGapEnd);
        }


        /* We compute the mean and variance for the estimated total gap size
           for this particular clone. We start with the edge mean and variance
           which already takes into account the clone mean and variance and
           the portions of the CIs containing the clone ends mean and variance.
           Next we subtract the length of all of the CIs spanned by the clone.
           Again we assume the variances are additive based on independence. */
        for(constant = edge->distance.mean,
              constantVariance = edge->distance.variance,
              lengthCIsIndex = thisCI->indexInScaffold + 1;
            lengthCIsIndex < otherCI->indexInScaffold; lengthCIsIndex++){
          constant -= lengthCIs[lengthCIsIndex].mean;
          constantVariance += lengthCIs[lengthCIsIndex].variance;
        }
        /* If we are recomputing gap sizes after setting some of the gaps to
           a fixed size based on the lack of an expected overlap then we need
           to take these fixed gaps and their variances into account -
           otherwise this loop is a no-op. The question is how to adjust the
           existing variance if at all. Adding it in produces huge variances
           which seems wrong but not doing anything seems wrong too. */
        for(gapIndex = thisCI->indexInScaffold;
            gapIndex < otherCI->indexInScaffold; gapIndex++){
          if(gapsToComputeGaps[gapIndex] == NULLINDEX){
            constant -= gapSize[gapIndex];
            //constantVariance += gapSizeVariance[gapIndex];
          }
        }
        /* cloneMean and cloneVariance are the statistics for the estimated
           total size of the gaps spanned by this clone. */
        cloneMean[indexClones] = constant;
        cloneVariance[indexClones] = constantVariance;

        if (debug.recomputeOffsetsVerboseLV > 1)
          fprintf(stderr, "Gap clone %f,%f (%d,%d)\n",
                  constant, sqrt(constantVariance),
                  thisCI->indexInScaffold, otherCI->indexInScaffold);

        constant /= constantVariance;
        inverseVariance = 1.0 / constantVariance;
        /* Store which gaps each clone spans so that we can iterate over
           these gaps when we calculate the gap variances and the
           squared error. */
        cloneGapStart[indexClones] = thisCI->indexInScaffold;
        cloneGapEnd[indexClones] = otherCI->indexInScaffold;
        /* Below we incrementally add to the matrices and vector we need for
           solving our equations. When we take the partial derivatives and
           set them to zero we get numGaps equations which we can represent
           as a vector on one side of equation by moving the constant terms
           to one side and a matrix times our set of gap size variables on
           the other. The vector is called gapConstants and the matrix
           gapCoefficients. As expected gapConstants is stored as a one
           dimensional array. The storage for gapCoefficients is also a
           one dimensional array but it represents a more complicated
           data structure. First due to the local effects of the clones
           on the scaffold the array is usually banded so for efficiency
           we only store the nonzero bands and in addition the matrix is
           symmetric so we only store the lower bands (subdiagonals) plus
           the main diagonal. The LAPACK interface expects the subdiagonals
           to be padded out to the same length as the diagonal and to be in
           column major order with the diagonals stored as rows so we
           comply. */
        for(colIndex = thisCI->indexInScaffold;
            colIndex < otherCI->indexInScaffold; colIndex++){
          int rowIndex;
          int colComputeIndex = gapsToComputeGaps[colIndex];
          /* For each gap that the clone spans it contributes the same
             constant value to the gapConstants vector which is equal
             to the mean total gap size for that clone divided by the
             variance of the total gap size. */
          if(colComputeIndex == NULLINDEX){
            continue;
          }
          gapConstants[colComputeIndex] += constant;
          for(rowIndex = colIndex;
              rowIndex < otherCI->indexInScaffold; rowIndex++){
            int rowComputeIndex = gapsToComputeGaps[rowIndex];
            /* If the number of gaps spanned by the clone is N then this clone
               contributes to NxN terms in the gapCoefficients matrix, but
               because the matrix is symmetric we only store the lower triangle
               so N*(N+1)/2 terms are affected for this clone. Remember that we
               store the (sub)diagonals as rows in column major order because
               the matrix tends to be banded and to use the LAPACK interface.
               The contribution of this clone to each term is the inverse of
               the variance of the total gap size for that clone. */
            if(rowComputeIndex == NULLINDEX){
              continue;
            }
            gapCoefficients[(colComputeIndex * maxDiagonals)
                            + (rowComputeIndex - colComputeIndex)] += inverseVariance;
          }
        }
        indexClones++;
      }
    }

    maxClone=indexClones;
    {
      FTN_INT nrhs = 1;
      FTN_INT bands = maxDiagonals - 1;
      FTN_INT ldab = maxDiagonals;
      FTN_INT rows = numComputeGaps;
      FTN_INT info = 0;

      if (debug.recomputeOffsetsVerboseLV > 1) {
        int i = 0;
        double *gapEnd, *gapPtr;
        for(gapPtr = gapConstants, gapEnd = gapPtr + numComputeGaps;
            gapPtr < gapEnd; gapPtr++){
          fprintf(stderr, "Gap Constants %g\n", *gapPtr);
        }
        fprintf(stderr, "Gap Coefficients\n");
        for(gapPtr = gapCoefficients, gapEnd = gapPtr + (maxDiagonals * numComputeGaps);
            gapPtr < gapEnd; gapPtr++){
          fprintf(stderr, "%g", *gapPtr);
          i++;
          if(i == maxDiagonals){
            fprintf(stderr, "\n");
            i = 0;
          }else{
            fprintf(stderr, "\t\t");
          }
        }

        fprintf(stderr, "rows " F_FTN_INT " bands " F_FTN_INT " ldab " F_FTN_INT " info " F_FTN_INT "\n",
                rows, bands, ldab, info);
      }

      dpbtrf_("L", &rows, &bands, gapCoefficients, &ldab, &info);
      if (debug.recomputeOffsetsVerboseLV > 1)
        fprintf(stderr, "dpbtrf: rows " F_FTN_INT " bands " F_FTN_INT " ldab " F_FTN_INT " info " F_FTN_INT "\n",
                rows, bands, ldab, info);
      if(info < 0){
        freeRecomputeData(&data);
        assert(0 /* RECOMPUTE_LAPACK */);
        return (RECOMPUTE_LAPACK);
      }else if(info > 0){
        freeRecomputeData(&data);

        fprintf(stderr,"SOMEBODY IS SCREWING UP SCAFFOLDING -- RecomputeOffsetsInScaffold has a singularity -- assert skipped!\n");

        // mjf 3/9/2001
        // this assert was causing trouble in the mouse_20010307 run, commented it out
        // and the run proceeded w/o further trouble
        // need to figure out why scaffolds that were apparently connected go singular
        //
        if (debug.fixedRecomputeSingluarLV) {
          DumpACIScaffoldNew(stderr,ScaffoldGraph,scaffold,TRUE);
          DumpACIScaffoldNew(stderr,ScaffoldGraph,scaffold,FALSE);
          assert(0 /* RECOMPUTE_SINGULAR */);
        }

        return (RECOMPUTE_SINGULAR);
      }
      /* Call an LAPACK routine to multiply the inverse of the gapCoefficients
         matrix by the gapConstants vector resulting in the least squares
         minimal solution of the gap sizes being returned in the gapConstants
         vector. */
      dpbtrs_("L", &rows, &bands, &nrhs, gapCoefficients, &ldab,
              gapConstants, &rows, &info);
      if (debug.recomputeOffsetsVerboseLV > 1)
        fprintf(stderr, "dpbtrs (call1): rows " F_FTN_INT " bands " F_FTN_INT " ldab " F_FTN_INT " nrhs " F_FTN_INT " info " F_FTN_INT "\n",
                rows, bands, ldab, nrhs, info);
      if(info < 0){
        freeRecomputeData(&data);
        assert(0 /* RECOMPUTE_LAPACK */);
        return (RECOMPUTE_LAPACK);
      }else if(info > 0){
        freeRecomputeData(&data);
        assert(0 /* RECOMPUTE_SINGULAR */);
        return (RECOMPUTE_SINGULAR);
      }
    }

    squaredError = 0;
    for(indexClones = 0; indexClones < maxClone; indexClones++){
      int gapIndex;
      int contributesToVariance = FALSE;

      /* We compute the squared error and gap size variances incrementally
         by adding the contribution from each clone. */
      for(gapIndex = 0; gapIndex < numComputeGaps; gapIndex++){
        spannedGaps[gapIndex] = 0.0;
      }
      for(gapIndex = cloneGapStart[indexClones];
          gapIndex < cloneGapEnd[indexClones]; gapIndex++){
        /* Compute the expected total gap size for this clone minus the solved
           for gap sizes that this clone spans. */
        if(gapsToComputeGaps[gapIndex] != NULLINDEX){
          cloneMean[indexClones] -= gapConstants[gapsToComputeGaps[gapIndex]];
          /* Finish creating a vector whose components are 0.0 for gaps not
             spanned by this clone and 1.0 for gaps that are. */
          spannedGaps[gapsToComputeGaps[gapIndex]] = 1.0;
          contributesToVariance = TRUE;
        }
      }
      /* To compute the squared error we square the difference between
         the expected total gap size for this clone minus the solved
         for gap sizes that this clone spans and divide by the clone
         variance. */
      squaredError += (cloneMean[indexClones] * cloneMean[indexClones]) /
        cloneVariance[indexClones];
      if(contributesToVariance){
        FTN_INT nrhs = 1;
        FTN_INT bands = maxDiagonals - 1;
        FTN_INT ldab = maxDiagonals;
        FTN_INT rows = numComputeGaps;
        FTN_INT info = 0;
        double *gapEnd, *gapPtr, *gapPtr2;

        /* Multiply the inverse of the gapCoefficients matrix times the vector
           of which gaps were spanned by this clone to produce the derivative
           of the gap sizes with respect to this clone (actually we would need
           to divide by the total gap variance for this clone
           but we correct for this below).
           This is computed in order to get an estimate of the variance for
           the gap sizes we have determined as outlined in equation 5-7 page
           70 of Data Reduction and Error Analysis for the Physical Sciences
           by Philip R. Bevington. */
        dpbtrs_("L", &rows, &bands, &nrhs, gapCoefficients, &ldab,
                spannedGaps, &rows, &info);
        if (debug.recomputeOffsetsVerboseLV > 1)
          fprintf(stderr, "dpbtrs (call2): rows " F_FTN_INT " bands " F_FTN_INT " ldab " F_FTN_INT " nrhs " F_FTN_INT " info " F_FTN_INT "\n",
                  rows, bands, ldab, nrhs, info);
        if(info < 0){
          freeRecomputeData(&data);
          assert(0 /* RECOMPUTE_LAPACK */);
          return (RECOMPUTE_LAPACK);
        }else if(info > 0){
          freeRecomputeData(&data);
          assert(0 /* RECOMPUTE_LAPACK */);
          return (RECOMPUTE_LAPACK);
        }
        for(gapPtr = spannedGaps, gapEnd = gapPtr + numComputeGaps,
              gapPtr2 = gapVariance;
            gapPtr < gapEnd; gapPtr++, gapPtr2++){
          /* According to equation 5-7 we need to square the derivative and
             multiply by the total gap variance for this clone but instead
             we end up dividing by the total gap variance for this clone
             because we neglected to divide by it before squaring and so
             the net result is to need to divide by it. */
          double term;
          term = *gapPtr;
          term *= term;
          term /= cloneVariance[indexClones];
          *gapPtr2 += term;
        }
      }
    }

    {
      int gapIndex, computeGapIndex;
      for(gapIndex = 0; gapIndex < numComputeGaps; gapIndex++){
        gapSize[computeGapsToGaps[gapIndex]] = gapConstants[gapIndex];
        gapSizeVariance[computeGapsToGaps[gapIndex]] = gapVariance[gapIndex];
        if (debug.recomputeOffsetsVerboseLV > 1)
          fprintf(stderr,"GapSize(%d:%d) %f:%f\n", gapIndex,
                  computeGapsToGaps[gapIndex], gapConstants[gapIndex],
                  sqrt(gapVariance[gapIndex]));
      }

      if(forceNonOverlaps){
        int alternate;

        hardConstraintSet = FALSE;
        InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
        for(gapIndex = 0, computeGapIndex = 0,
              prevCI = NextCIScaffoldTIterator(&CIs);
            (thisCI = NextCIScaffoldTIterator(&CIs)) != NULL;
            prevCI = thisCI, gapIndex++){
          ChunkOrientationType edgeOrient;
          float chiSquaredValue;
          EdgeCGW_T *overlapEdge;
          if(gapsToComputeGaps[gapIndex] == NULLINDEX){
            continue;
          }
          if(gapSize[gapIndex] >= - CGW_MISSED_OVERLAP){
            computeGapsToGaps[computeGapIndex] = gapIndex;
            gapsToComputeGaps[gapIndex] = computeGapIndex;
            computeGapIndex++;
            continue;
          }
          if(GetNodeOrient(thisCI) == A_B){
            if(GetNodeOrient(prevCI) == A_B){
              edgeOrient = AB_AB;
            }else{//GetNodeOrient(prevCI) == B_A
              edgeOrient = BA_AB;
            }
          }else{//GetNodeOrient(thisCI) == B_A
            if(GetNodeOrient(prevCI) == A_B){
              edgeOrient = AB_BA;
            }else{//GetNodeOrient(prevCI) == B_A
              edgeOrient = BA_BA;
            }
          }
          overlapEdge = FindOverlapEdgeChiSquare(graph, prevCI, thisCI->id,
                                                 edgeOrient, gapSize[gapIndex],
                                                 gapSizeVariance[gapIndex],
                                                 &chiSquaredValue,
                                                 (float)PAIRWISECHI2THRESHOLD_CGW,
                                                 &alternate, verbose);

          if(overlapEdge && alternate){ // found a node that is out of order!!!!

#ifdef TEST_FIXUPMISORDER
            //  XXXX:  BPW is slightly confused as to why we are not
            //  FixUpMisorderedContigs(), instead we are merging the contig.
            //  This was if'd out before I got here.
            //
            //  OK, so test it!
            //
            //  If it succeeds, return, if not, continue with the original method.
            //
            fprintf(stderr, "FixUpMisorderedContigs\n");
            if (FixUpMisorderedContigs(scaffold, prevCI, thisCI, edgeOrient, gapSize[gapIndex], gapSizeVariance[gapIndex], overlapEdge))
              return(RECOMPUTE_FAILED_REORDER_NEEDED);
#endif

            if (debug.leastSquaresGapsLV > 0) {
              fprintf(stderr,"*** Least Squares found the following alternate edge...contig NOW!\n");
              PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, " alternateEdge: ", overlapEdge, overlapEdge->idA);
            }
            if (debug.leastSquaresGapsLV > 1) {
              fprintf(stderr, "BEFORE ContigContainment, scaffold "F_CID" %s connected\n",
                      scaffold->id,
                      IsScaffoldInternallyConnected(ScaffoldGraph,
                                                    scaffold, ALL_TRUSTED_EDGES) ? "is" : "is NOT");
              dumpTrustedEdges(ScaffoldGraph, scaffold, ALL_TRUSTED_EDGES);
              DumpACIScaffold(stderr, graph, scaffold, FALSE);
              CheckInternalEdgeStatus(graph, scaffold, PAIRWISECHI2THRESHOLD_CGW, 100000000000.0, 0, FALSE);
            }

            // We want to merge the two contigs immediately, since these are problematic,
            // but we know we want to overlap these guys.
            //
            if (ContigContainment(scaffold, prevCI, thisCI, overlapEdge, TRUE) != TRUE) {
              fprintf(stderr, "ContigContainment failed.\n");
              assert(0);
            }

            if (debug.leastSquaresGapsLV > 1) {
              fprintf(stderr, "AFTER ContigContainment, scaffold "F_CID" %s connected\n",
                      scaffold->id,
                      IsScaffoldInternallyConnected(ScaffoldGraph,
                                                    scaffold, ALL_TRUSTED_EDGES) ? "is" : "is NOT");
              dumpTrustedEdges(ScaffoldGraph, scaffold, ALL_TRUSTED_EDGES);
              DumpACIScaffold(stderr, graph, scaffold, FALSE);
              CheckInternalEdgeStatus(graph, scaffold, PAIRWISECHI2THRESHOLD_CGW, 100000000000.0, 0, FALSE);
            }

            return RECOMPUTE_CONTIGGED_CONTAINMENTS;
          }
          if(overlapEdge && isContainmentEdge(overlapEdge)){

            if (debug.leastSquaresGapsLV > 0) {
              fprintf(stderr,"*** Least Squares found the following containment edge...contig NOW!\n");
              fprintf(stderr,"*** " F_CID " (length %g) should be contained within " F_CID " (length %g)\n",
                      thisCI->id, thisCI->bpLength.mean,
                      prevCI->id, prevCI->bpLength.mean);
              PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, " overlapEdge: ", overlapEdge, overlapEdge->idA);
            }
            if (debug.leastSquaresGapsLV > 1) {
              fprintf(stderr, "BEFORE ContigContainment, scaffold "F_CID" %s connected\n",
                      scaffold->id,
                      IsScaffoldInternallyConnected(ScaffoldGraph,
                                                    scaffold, ALL_TRUSTED_EDGES) ? "is" : "is NOT");
              dumpTrustedEdges(ScaffoldGraph, scaffold, ALL_TRUSTED_EDGES);
              DumpACIScaffold(stderr, graph, scaffold, FALSE);
              CheckInternalEdgeStatus(graph, scaffold, PAIRWISECHI2THRESHOLD_CGW, 100000000000.0, 0, FALSE);
            }

            // When we find a containment relationship in a scaffold we want to merge the two contigs
            // immediately, since containments have the potential to induce situations that are confusing
            // for least squares
            if (ContigContainment(scaffold, prevCI, thisCI, overlapEdge, TRUE) != TRUE) {
              fprintf(stderr, "ContigContainment failed.\n");
              assert(0);
            }

            if (debug.leastSquaresGapsLV > 1) {
              fprintf(stderr, "AFTER ContigContainment, scaffold "F_CID" %s connected\n",
                      scaffold->id,
                      IsScaffoldInternallyConnected(ScaffoldGraph,
                                                    scaffold, ALL_TRUSTED_EDGES) ? "is" : "is NOT");
              dumpTrustedEdges(ScaffoldGraph, scaffold, ALL_TRUSTED_EDGES);
              DumpACIScaffold(stderr, graph, scaffold, FALSE);
              CheckInternalEdgeStatus(graph, scaffold, PAIRWISECHI2THRESHOLD_CGW, 100000000000.0, 0, FALSE);
            }

            return RECOMPUTE_CONTIGGED_CONTAINMENTS;
          }


          if(overlapEdge == (EdgeCGW_T *)NULL){
            // gapsize is negative, so is leftedge
            double leftEdge = (gapSize[gapIndex] - 3 * sqrt(gapSizeVariance[gapIndex]));
            double newLeftEdge = leftEdge;
            double newStd = (newLeftEdge - (-CGW_MISSED_OVERLAP))/3.0;
            double newVariance = newStd * newStd;

            if (debug.recomputeOffsetsVerboseLV > 1)
              fprintf(stderr,"GapChange Gap(%d:%d) CIs: " F_CID "," F_CID " new:(%f,%f) old(%f,%f)\n",
                      gapsToComputeGaps[gapIndex], gapIndex,
                      prevCI->id, thisCI->id,
                      (float)(- CGW_MISSED_OVERLAP), newVariance,
                      gapSize[gapIndex], gapSizeVariance[gapIndex]);

            //
            // Adjust the gap variance so that the least squares computed mean position is within
            // 3 sigma of the CGW_MISSED_OVERLAP positioning.  Otherwise, mate-induced edges between
            // scaffold neighbors are not X-squared compatible with the scaffold positioning, and are
            // not included in the set of edges used for subsequence scaffold construction.
            //    ---------------------------
            //               -------------------------------------------
            //                  adjustment
            //      |        |<---------->|
            //      |     original      forced
            //      |      position      position
            //      | 
            //   original left edge of 3-sigma interval
            // goal is to adjust forced position variance, so the original left edge of 3 sigma interval
            // is within an adjusted 3 sigma range of the new mean position.


            gapSize[gapIndex] = - CGW_MISSED_OVERLAP;
            gapSizeVariance[gapIndex] = newVariance;
            gapsToComputeGaps[gapIndex] = NULLINDEX;
            hardConstraintSet = TRUE;
            continue;
          }
          if(abs(overlapEdge->distance.mean - gapSize[gapIndex]) <=
             MAX_OVERLAP_SLOP_CGW){
            computeGapsToGaps[computeGapIndex] = gapIndex;
            gapsToComputeGaps[gapIndex] = computeGapIndex;
            computeGapIndex++;
            continue;
          }

          if (debug.recomputeOffsetsVerboseLV > 1)
            fprintf(stderr,"GapChange(%d:%d) %f:%f\n",
                    gapsToComputeGaps[gapIndex], gapIndex,
                    overlapEdge->distance.mean, gapSize[gapIndex]);

          gapSize[gapIndex] = overlapEdge->distance.mean;
          gapsToComputeGaps[gapIndex] = NULLINDEX;
          hardConstraintSet = TRUE;
        }
        numComputeGaps = computeGapIndex;
        if(hardConstraintSet){
          double *gapEnd, *gapPtr;
          for(gapPtr = gapConstants, gapEnd = gapPtr + numComputeGaps;
              gapPtr < gapEnd; gapPtr++){
            *gapPtr = 0.0;
          }
          for(gapPtr = gapCoefficients,
                gapEnd = gapPtr + (maxDiagonals * numComputeGaps);
              gapPtr < gapEnd; gapPtr++){
            *gapPtr = 0.0;
          }
          for(gapPtr = gapVariance, gapEnd = gapPtr + numComputeGaps;
              gapPtr < gapEnd; gapPtr++){
            *gapPtr = 0.0;
          }
        }
      }
    }
  }while(forceNonOverlaps && hardConstraintSet && (numComputeGaps > 0));

  if (debug.recomputeOffsetsVerboseLV > 1)
    fprintf(stderr,"LSE: %f,%f #clones: %d,%d\n",
            scaffold->info.Scaffold.leastSquareError, squaredError,
            scaffold->info.Scaffold.numLeastSquareClones, numClones);

  scaffold->info.Scaffold.leastSquareError = squaredError;
  scaffold->info.Scaffold.numLeastSquareClones = numClones;
  {// Start

    double *gapPtr, *gapVarPtr;
    LengthT *prevLeftEnd, *prevRightEnd, *thisLeftEnd, *thisRightEnd;

    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

    if (debug.recomputeOffsetsVerboseLV > 1)
      fprintf(stderr,"Reestimate gaps for scaffold\n");
    prevCI = NextCIScaffoldTIterator(&CIs);
    if(GetNodeOrient(prevCI) == A_B){
      prevLeftEnd = &(prevCI->offsetAEnd);
      prevRightEnd = &(prevCI->offsetBEnd);
    }else{
      prevLeftEnd = &(prevCI->offsetBEnd);
      prevRightEnd = &(prevCI->offsetAEnd);
    }
    if (debug.recomputeOffsetsVerboseLV > 1)
      fprintf(stderr, "Old %f,%f ",
              prevLeftEnd->mean, sqrt(prevLeftEnd->variance));
    prevLeftEnd->mean = 0.0;
    prevLeftEnd->variance = 0.0;

    if (debug.recomputeOffsetsVerboseLV > 1)
      fprintf(stderr, "New %f,%f\n",
              prevLeftEnd->mean, sqrt(prevLeftEnd->variance));

    for(gapPtr = gapSize, gapVarPtr = gapSizeVariance;
        (thisCI = NextCIScaffoldTIterator(&CIs)) != NULL;
        prevCI = thisCI, prevLeftEnd = thisLeftEnd,
          prevRightEnd = thisRightEnd, gapPtr++, gapVarPtr++){
      LengthT gapDistance;
      CDS_COORD_t realDistance;

      if(GetNodeOrient(thisCI) == A_B){
        thisLeftEnd = &(thisCI->offsetAEnd);
        thisRightEnd = &(thisCI->offsetBEnd);
      }else{
        thisLeftEnd = &(thisCI->offsetBEnd);
        thisRightEnd = &(thisCI->offsetAEnd);
      }

      // Keep track of the biggest offset we've seen so far
      if(!maxOffset || 
         maxOffset->mean < thisRightEnd->mean){
        maxOffset = thisRightEnd;
      }  


      gapDistance.mean = thisLeftEnd->mean - prevRightEnd->mean;
      gapDistance.variance = thisLeftEnd->variance - prevRightEnd->variance;

      if (debug.recomputeOffsetsVerboseLV > 1)
        fprintf(stderr, "Old %f,%f ",
                prevRightEnd->mean, sqrt(prevRightEnd->variance));

      prevRightEnd->mean = prevLeftEnd->mean + prevCI->bpLength.mean;
      prevRightEnd->variance = prevLeftEnd->variance +
        prevCI->bpLength.variance;

      if (debug.recomputeOffsetsVerboseLV > 1) {
        fprintf(stderr, "New %f,%f\n",
                prevRightEnd->mean, sqrt(prevRightEnd->variance));
        fprintf(stderr, "Old %f,%f ",
                thisLeftEnd->mean, sqrt(thisLeftEnd->variance));
      }

      thisLeftEnd->mean = prevRightEnd->mean + *gapPtr;
      thisLeftEnd->variance = prevRightEnd->variance + *gapVarPtr;
      if (debug.recomputeOffsetsVerboseLV > 1)
        fprintf(stderr, "New %f,%f\n",
                thisLeftEnd->mean, sqrt(thisLeftEnd->variance));
      if(GetNodeOrient(thisCI) == A_B){
        if(GetNodeOrient(prevCI) == A_B){
          realDistance = thisCI->aEndCoord - prevCI->bEndCoord;
        }else{//GetNodeOrient(prevCI) == B_A
          realDistance = thisCI->aEndCoord - prevCI->aEndCoord;
        }
      }else{//GetNodeOrient(thisCI) == B_A
        if(GetNodeOrient(prevCI) == A_B){
          realDistance = thisCI->bEndCoord - prevCI->bEndCoord;
        }else{//GetNodeOrient(prevCI) == B_A
          realDistance = thisCI->bEndCoord - prevCI->aEndCoord;
        }
      }
      if (debug.recomputeOffsetsVerboseLV > 1)
        fprintf(stderr, "Old %f New %f Real " F_COORD " StdDev %f,%f\n",
                gapDistance.mean, *gapPtr, realDistance,
                sqrt(gapDistance.variance), sqrt(*gapVarPtr));

    }
    if (debug.recomputeOffsetsVerboseLV > 1)
      fprintf(stderr, "Old %f,%f ",
              prevRightEnd->mean, sqrt(prevRightEnd->variance));
    prevRightEnd->mean = prevLeftEnd->mean + prevCI->bpLength.mean;
    prevRightEnd->variance = prevLeftEnd->variance +
      prevCI->bpLength.variance;
    if (debug.recomputeOffsetsVerboseLV > 1)
      fprintf(stderr, "New %f,%f\n",
              prevRightEnd->mean, sqrt(prevRightEnd->variance));

    if (debug.recomputeOffsetsVerboseLV > 2) {
      fprintf(stderr,"Scaffold " F_CID " Length %f,%f-%f,%f\n",
              scaffold->id,
              scaffold->bpLength.mean, sqrt(scaffold->bpLength.variance),
              maxOffset->mean, sqrt(maxOffset->variance));
      DumpCIScaffold(stderr, graph, scaffold, FALSE);
    }

    scaffold->bpLength = *maxOffset;

    SetCIScaffoldTLength(ScaffoldGraph, scaffold, TRUE); // recompute scaffold length, just to be sure

  }

  StopTimerT(&GlobalData->RecomputeOffsetsTimer);   /*  STOP */

  freeRecomputeData(&data);
  return (RECOMPUTE_OK);
}

void MarkInternalEdgeStatus(ScaffoldGraphT *graph, CIScaffoldT *scaffold, 
                            float pairwiseChiSquaredThreshhold,
                            float maxVariance,
                            int markTrusted, int markUntrusted,
                            int doNotChange,
                            int operateOnMerged){
  CIScaffoldTIterator CIs;
  /* Iterate over all of the CIEdges */
  GraphEdgeIterator edges;
  EdgeCGW_T *edge;
  NodeCGW_T *thisCI;
  int32 numCIs;
  int32 indexCIs;
  int internalEdges = 0;  /* Number of merged edges (not including UNTRUSTED)
                             that are internal to scaffold */
  int confirmedInternalEdges = 0; /* Number of merged edges confirmed by
                                     current CI positions */

  //  Figure out where to put our debug info
  //
  FILE  *debugfp = debug.markInternalEdgeStatusFP;
  if ((debugfp == NULL) && (debug.markInternalEdgeStatusLV > 0))
    debugfp = stderr;

  if (debug.markInternalEdgeStatusLV > 1)
    fprintf(debugfp, "Marking Edges for Scaffold " F_CID " markTrusted=%d markUntrusted=%d\n",
            scaffold->id, markTrusted, markUntrusted);

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  for(indexCIs = 0;
      (thisCI = NextCIScaffoldTIterator(&CIs)) != NULL;){
    thisCI->indexInScaffold = indexCIs;
    indexCIs++;
  }
  numCIs = indexCIs;
  if(numCIs != scaffold->info.Scaffold.numElements){
    if (debug.markInternalEdgeStatusLV > 0)
      fprintf(debugfp, "NumElements inconsistent.  Reset from %d to %d\n", scaffold->info.Scaffold.numElements, numCIs);
    scaffold->info.Scaffold.numElements = numCIs;
  }

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL){
    int flags = GRAPH_EDGE_DEFAULT;

    if(!operateOnMerged)
      flags |= GRAPH_EDGE_RAW_ONLY;

    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, 
                          ALL_END,
                          ALL_EDGES, 
                          flags,
                          &edges);// Use merged edges

    while((edge = NextGraphEdgeIterator(&edges))!= NULL){
      int isA = (edge->idA == thisCI->id);
      NodeCGW_T *otherCI =
        GetGraphNode(ScaffoldGraph->RezGraph,
                     (isA? edge->idB: edge->idA));
      ChunkOrientationType edgeOrient;
      FragOrient thisCIorient, otherCIorient;
      LengthT gapDistance;
      float chiSquareResult;


      if (debug.markInternalEdgeStatusLV > 1)
        fprintf(debugfp, " examining edge [" F_CID "," F_CID "] (%f,%f) weight %d ori: %c\n",
                thisCI->id, otherCI->id,
                edge->distance.mean, edge->distance.variance,
                edge->edgesContributing,
                GetEdgeOrientationWRT(edge, thisCI->id));

      /* We do not want to change the labels for edges with certain
         labels as specified in the doNotChange mask. */
      if(edge->flags.bits.edgeStatus & doNotChange){
        continue;
      }
      /* We only want to label edges between CIs in the same scaffold. */
      if(otherCI->scaffoldID != thisCI->scaffoldID){
        SetEdgeStatus(graph->RezGraph, edge, INTER_SCAFFOLD_EDGE_STATUS);
        continue;
      }
      /* We only want to label an edge once and the
         iterator will visit intra scaffold edges twice so we use
         the condition that the edge is canonical - that is the
         condition thisCI->indexInScaffold < otherCI->indexInScaffold
         to guarantee labelling the edge only once. */
      if(otherCI->indexInScaffold <= thisCI->indexInScaffold){
        continue;
      }

      /* Check that the edge orientation is consistent with the CI positions
         and orientations within the scaffold. */
      edgeOrient = GetEdgeOrientationWRT(edge, thisCI->id);
      switch(edgeOrient){
        case AB_BA:
          //      thisCI                                        otherCI
          //  A --------------------- B               B --------------------- A
          //    5'----->                                           <------5'
          thisCIorient = A_B;
          otherCIorient = B_A;
          gapDistance.mean = otherCI->offsetBEnd.mean - thisCI->offsetBEnd.mean;
          gapDistance.variance = otherCI->offsetBEnd.variance -
            thisCI->offsetBEnd.variance;
          break;
        case AB_AB:
          //      thisCI                                        otherCI
          //  A --------------------- B               A --------------------- B
          //    5'----->                                           <------5'
          thisCIorient = A_B;
          otherCIorient = A_B;
          gapDistance.mean = otherCI->offsetAEnd.mean - thisCI->offsetBEnd.mean;
          gapDistance.variance = otherCI->offsetAEnd.variance -
            thisCI->offsetBEnd.variance;
          break;
        case BA_BA:
          //      thisCI                                        otherCI
          //  B --------------------- A               B --------------------- A
          //    5'----->                                           <------5'
          thisCIorient = B_A;
          otherCIorient = B_A;
          gapDistance.mean = otherCI->offsetBEnd.mean - thisCI->offsetAEnd.mean;
          gapDistance.variance = otherCI->offsetBEnd.variance -
            thisCI->offsetAEnd.variance;
          break;
        case BA_AB:
          //      thisCI                                        otherCI
          //  B --------------------- A               A --------------------- B
          //    5'----->                                           <------5'
          thisCIorient = B_A;
          otherCIorient = A_B;
          gapDistance.mean = otherCI->offsetAEnd.mean - thisCI->offsetAEnd.mean;
          gapDistance.variance = otherCI->offsetAEnd.variance -
            thisCI->offsetAEnd.variance;
          break;
        default:
          assert(0);
          break;
      }
      if((GetNodeOrient(thisCI) != thisCIorient) ||
         (GetNodeOrient(otherCI) != otherCIorient)){
        /* Mark as untrusted an edge whose orientation does not agree
           with the orientation of the CIs in the scaffold. */

        if (debug.markInternalEdgeStatusLV > 0)
          fprintf(debugfp, "[" F_CID "," F_CID "]Bad orientation (%c,%c) (%c,%c)\n",
                  thisCI->id, otherCI->id,
                  GetNodeOrient(thisCI), thisCIorient, GetNodeOrient(otherCI),
                  otherCIorient);

        SetEdgeStatus(graph->RezGraph, edge, markUntrusted ? UNTRUSTED_EDGE_STATUS :
                      TENTATIVE_UNTRUSTED_EDGE_STATUS);
        continue;
      }
      if(gapDistance.variance <= 0.0){
        /* This condition should not occur but when it does it causes an
           assert in the PairwiseChiSquare subroutine so we will just mark
           the edge as untrusted so the program can keep going but this
           needs to be investigated and fixed!!!!
        */

        //  If you're serious about debugging this, enable these dumps, otherwise,
        //  don't fill up the cgwlog with useless crud.
        if (debug.markInternalEdgeStatusLV > 1) {
          DumpACIScaffoldNew(debugfp,ScaffoldGraph,scaffold,TRUE);
          DumpACIScaffoldNew(debugfp,ScaffoldGraph,scaffold,FALSE);
        }

        if (debug.markInternalEdgeStatusLV > 0)
          fprintf(debugfp, "[" F_CID "." F_CID "," F_CID "." F_CID "]Bad Gap Variance (%f,%f) (%f,%f) DANGER WILL ROBINSON!!!\n",
                  thisCI->id, thisCI->scaffoldID, otherCI->id,
                  otherCI->scaffoldID,
                  gapDistance.mean, gapDistance.variance,
                  edge->distance.mean, edge->distance.variance);

        SetEdgeStatus(graph->RezGraph, edge, markUntrusted ? UNTRUSTED_EDGE_STATUS :
                      TENTATIVE_UNTRUSTED_EDGE_STATUS);
        continue;
      }
      if(!PairwiseChiSquare((float)gapDistance.mean,
                            gapDistance.variance,
                            (float)edge->distance.mean,
                            edge->distance.variance,
                            (LengthT *)NULL, &chiSquareResult,
                            pairwiseChiSquaredThreshhold)){
        /* Mark  this edge as untrusted if the distance of the edge is not
           consistent with the estimated gap distance as judged by the
           Chi Squared Test. */
        if (debug.markInternalEdgeStatusLV > 1)
          fprintf(debugfp, "[" F_CID "," F_CID "]Bad Chi Squared %f (%f,%f) (%f,%f)\n",
                  thisCI->id, otherCI->id,
                  chiSquareResult, gapDistance.mean, gapDistance.variance,
                  edge->distance.mean, edge->distance.variance);

        SetEdgeStatus(graph->RezGraph, edge, markUntrusted ? UNTRUSTED_EDGE_STATUS :
                      TENTATIVE_UNTRUSTED_EDGE_STATUS);
        continue;
      }

      if (debug.markInternalEdgeStatusLV > 1)
        fprintf(debugfp, "[" F_CID "," F_CID "]Good Chi Squared %f (%f,%f) (%f,%f)\n",
                thisCI->id, otherCI->id,
                chiSquareResult, gapDistance.mean, gapDistance.variance,
                edge->distance.mean, edge->distance.variance);

      if(edge->distance.variance > maxVariance){
        if (debug.markInternalEdgeStatusLV > 1)
          fprintf(debugfp, "[" F_CID "," F_CID "]Variance too large %f\n",
                  thisCI->id, otherCI->id,
                  edge->distance.variance);

        SetEdgeStatus(graph->RezGraph,edge, LARGE_VARIANCE_EDGE_STATUS);
        continue;
      }

      if (debug.markInternalEdgeStatusLV > 1)
        fprintf(debugfp, "[" F_CID "," F_CID "]Trusted Edge!  Mean %f Variance %f\n",
                thisCI->id, otherCI->id,
                edge->distance.mean,
                edge->distance.variance);

      SetEdgeStatus(graph->RezGraph,edge, markTrusted ? TRUSTED_EDGE_STATUS :
                    TENTATIVE_TRUSTED_EDGE_STATUS);
    }
  }

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL){

    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, 
                          ALL_END,
                          ALL_EDGES, 
                          GRAPH_EDGE_DEFAULT,
                          &edges);// Use merged edges

    while((edge = NextGraphEdgeIterator(&edges))!= NULL){
      int isA = (edge->idA == thisCI->id);
      NodeCGW_T *otherCI =
        GetGraphNode(ScaffoldGraph->RezGraph,
                     (isA? edge->idB: edge->idA));
      EdgeStatus edgeStatus;

      if(otherCI->indexInScaffold <= thisCI->indexInScaffold){
        continue; // Only interested in looking at an edge once
      }
      edgeStatus = GetEdgeStatus(edge);
      if((edgeStatus == TRUSTED_EDGE_STATUS) ||
         (edgeStatus == TENTATIVE_TRUSTED_EDGE_STATUS)){
        internalEdges++;
        confirmedInternalEdges++;
      }else if(edgeStatus == TENTATIVE_UNTRUSTED_EDGE_STATUS){
        internalEdges++;
      }
    }
  }

  if (debug.markInternalEdgeStatusLV > 1)
    fprintf(debugfp,"#Internal Edges %d,%d confirmed %d,%d\n",
            scaffold->info.Scaffold.internalEdges, internalEdges,
            scaffold->info.Scaffold.confirmedInternalEdges, confirmedInternalEdges);

  scaffold->info.Scaffold.internalEdges = internalEdges;
  scaffold->info.Scaffold.confirmedInternalEdges = confirmedInternalEdges;
  return;
}


static int EdgeAndGapAreVaguelyCompatible(double distDiff,
                                          double diffVar,
                                          double maxDiffSlop,
                                          double maxSigmaSlop){
  if(distDiff>maxDiffSlop)return FALSE;
  if(diffVar<0){
    fprintf(stderr,"WARNING: variance difference is negative -- probably trouble with variances after interleaving\n"); 
  }
  if(abs(distDiff)/sqrt(abs(diffVar))>maxSigmaSlop)return FALSE;
  return TRUE;
}


/// New test code to help handle slightly messier cases!
int IsInternalEdgeStatusVaguelyOK(EdgeCGW_T *edge,CDS_CID_t thisCIid){

  NodeCGW_T *thisCI = GetGraphNode(ScaffoldGraph->RezGraph,thisCIid);
  GraphEdgeIterator edges;

  //  Figure out where to put our debug info
  //
  FILE  *debugfp = debug.markInternalEdgeStatusFP;
  if ((debugfp == NULL) && (debug.markInternalEdgeStatusLV > 0))
    debugfp = stderr;

  // WE ASSUME edge COMES FROM SOMETHING LIKE THE FOLLOWING:
  //
  //
  //  InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, 
  //                        ALL_END,
  //                        ALL_EDGES, 
  //                        GRAPH_EDGES_RAW_ONLY,
  //                        &edges); 
  //  while((edge = NextGraphEdgeIterator(&edges))!= NULL){
  //

  int isA = (edge->idA == thisCI->id);
  NodeCGW_T *otherCI =
    GetGraphNode(ScaffoldGraph->RezGraph,
                 (isA? edge->idB: edge->idA));
  ChunkOrientationType edgeOrient;
  FragOrient thisCIorient, otherCIorient;
  LengthT gapDistance;
  float chiSquareResult;

  /* We only want to label edges between CIs in the same scaffold; this is up to the caller. */
  if(otherCI->scaffoldID != thisCI->scaffoldID){

    //  This is a good candidate for including in a file of serious warnings!

    fprintf(stderr, "WARNING:  IsInternalEdgeStatusVaguelyOK() got an edge between different scaffolds.\n");
    fprintf(stderr, "   edge:  idA:%d idB:%d\n", edge->idA, edge->idB);
    fprintf(stderr, " thisCI:  scaffoldID:%d isA:%d\n", thisCI->scaffoldID, edge->idA == thisCI->id);
    fprintf(stderr, "otherCI:  scaffoldID:%d isA:%d\n", otherCI->scaffoldID, edge->idA == otherCI->id);

    //  With motivation from MarkInternalEdgeStatus(), we simply mark
    //  this edge as inter scaffold and return.
    //
    SetEdgeStatus(ScaffoldGraph->RezGraph, edge, INTER_SCAFFOLD_EDGE_STATUS);
    //assert(0);
    return(0);
  }

  /* We only want to do this to an edge once; the calling routine should
     ensure that the edge is canonical - that is the
     condition thisCI->indexInScaffold < otherCI->indexInScaffold;
     it would probably be harmless to allow this through, but why bother to let the 
     programmer off easy?
  */
  if(otherCI->indexInScaffold <= thisCI->indexInScaffold){
    assert(0);
  }
      
  /* Check that the edge orientation is consistent with the CI positions
     and orientations within the scaffold. */
  edgeOrient = GetEdgeOrientationWRT(edge, thisCI->id);
  switch(edgeOrient){
    case AB_BA:
      //      thisCI                                        otherCI
      //  A --------------------- B               B --------------------- A
      //    5'----->                                           <------5'
      thisCIorient = A_B;
      otherCIorient = B_A;
      gapDistance.mean = otherCI->offsetBEnd.mean - thisCI->offsetBEnd.mean;
      gapDistance.variance = otherCI->offsetBEnd.variance -
        thisCI->offsetBEnd.variance;
      break;
    case AB_AB:
      //      thisCI                                        otherCI
      //  A --------------------- B               A --------------------- B
      //    5'----->                                           <------5'
      thisCIorient = A_B;
      otherCIorient = A_B;
      gapDistance.mean = otherCI->offsetAEnd.mean - thisCI->offsetBEnd.mean;
      gapDistance.variance = otherCI->offsetAEnd.variance -
        thisCI->offsetBEnd.variance;
      break;
    case BA_BA:
      //      thisCI                                        otherCI
      //  B --------------------- A               B --------------------- A
      //    5'----->                                           <------5'
      thisCIorient = B_A;
      otherCIorient = B_A;
      gapDistance.mean = otherCI->offsetBEnd.mean - thisCI->offsetAEnd.mean;
      gapDistance.variance = otherCI->offsetBEnd.variance -
        thisCI->offsetAEnd.variance;
      break;
    case BA_AB:
      //      thisCI                                        otherCI
      //  B --------------------- A               A --------------------- B
      //    5'----->                                           <------5'
      thisCIorient = B_A;
      otherCIorient = A_B;
      gapDistance.mean = otherCI->offsetAEnd.mean - thisCI->offsetAEnd.mean;
      gapDistance.variance = otherCI->offsetAEnd.variance -
        thisCI->offsetAEnd.variance;
      break;
    default:
      assert(0);
      break;
  }
  if((GetNodeOrient(thisCI) != thisCIorient) ||
     (GetNodeOrient(otherCI) != otherCIorient)){
    /* edge orientation does not agree
       with the orientation of the CIs in the scaffold. */

    if (debug.markInternalEdgeStatusLV > 1)
      fprintf(debugfp, "[" F_CID "," F_CID "]Bad orientation (%c,%c) (%c,%c)\n",
              thisCI->id, otherCI->id,
              GetNodeOrient(thisCI), thisCIorient, GetNodeOrient(otherCI),
              otherCIorient);

    return FALSE;
  }
  if(gapDistance.variance <= 0.0){
    /* This condition should not occur, so kill it now!
     */

    if (debug.markInternalEdgeStatusLV > 0)
      fprintf(debugfp, "[" F_CID "." F_CID "," F_CID "." F_CID "]Bad Gap Variance (%f,%f) (%f,%f) DANGER WILL ROBINSON!!!\n",
              thisCI->id, thisCI->scaffoldID, otherCI->id,
              otherCI->scaffoldID,
              gapDistance.mean, gapDistance.variance,
              edge->distance.mean, edge->distance.variance);

#ifdef NEG_GAP_VARIANCE_PROBLEM_FIXED
    DumpACIScaffoldNew(debugfp,ScaffoldGraph,GetGraphNode(ScaffoldGraph->ScaffoldGraph,thisCI->scaffoldID),TRUE);
    DumpACIScaffoldNew(debugfp,ScaffoldGraph,GetGraphNode(ScaffoldGraph->ScaffoldGraph,thisCI->scaffoldID),FALSE);
    assert(0);
#endif
  }
  if(edge->distance.variance <= 0.0){
    /* This condition should not occur, so kill it now!
     */

    if (debug.markInternalEdgeStatusLV > 0)
      fprintf(debugfp, "[" F_CID "." F_CID "," F_CID "." F_CID "]Bad Edge Variance (%f,%f) (%f,%f) DANGER WILL ROBINSON!!!\n",
              thisCI->id, thisCI->scaffoldID, otherCI->id,
              otherCI->scaffoldID,
              gapDistance.mean, gapDistance.variance,
              edge->distance.mean, edge->distance.variance);

    DumpACIScaffoldNew(debugfp,ScaffoldGraph,GetGraphNode(ScaffoldGraph->ScaffoldGraph,thisCI->scaffoldID),TRUE);
    DumpACIScaffoldNew(debugfp,ScaffoldGraph,GetGraphNode(ScaffoldGraph->ScaffoldGraph,thisCI->scaffoldID),FALSE);
    assert(0);
  }
  if(!EdgeAndGapAreVaguelyCompatible((double)gapDistance.mean-edge->distance.mean,
                                     (double)gapDistance.variance+edge->distance.variance,
                                     MAX_ABSOLUTE_SLOP,
                                     MAX_SIGMA_SLOP)){
    /* Mark  this edge as untrusted if the distance of the edge is not
       consistent with the estimated gap distance as judged by the
       Chi Squared Test. */

    if (debug.markInternalEdgeStatusLV > 1)
      fprintf(debugfp, "[" F_CID "," F_CID "]Not vaguely compatible %f (%f,%f) (%f,%f)\n",
              thisCI->id, otherCI->id,
              chiSquareResult, gapDistance.mean, gapDistance.variance,
              edge->distance.mean, edge->distance.variance);

    return FALSE;
  }

  if (debug.markInternalEdgeStatusLV > 1)
    fprintf(debugfp, "[" F_CID "," F_CID "]At least vaguely compatible %f (%f,%f) (%f,%f)\n",
            thisCI->id, otherCI->id,
            chiSquareResult, gapDistance.mean, gapDistance.variance,
            edge->distance.mean, edge->distance.variance);

  return TRUE;
}



void  CheckLSScaffoldWierdnesses(char *string, ScaffoldGraphT *graph, CIScaffoldT *scaffold){
  CIScaffoldTIterator CIs;
  ChunkInstanceT *firstCI, *secondCI;
  LengthT delta, *minOffsetp, *maxOffsetp;
  
  InitCIScaffoldTIterator(graph, scaffold, TRUE,  FALSE, &CIs);
  firstCI = NextCIScaffoldTIterator(&CIs);
  if(firstCI->offsetAEnd.mean < firstCI->offsetBEnd.mean){
    minOffsetp = &firstCI->offsetAEnd;;
    maxOffsetp = &firstCI->offsetBEnd;
  }else{
    minOffsetp = &firstCI->offsetBEnd;;
    maxOffsetp = &firstCI->offsetAEnd;
  }
  
  if(minOffsetp->mean < 0.0){
    delta.mean = -minOffsetp->mean;
    delta.variance = minOffsetp->variance;

    if (debug.leastSquaresGapsLV > 0)
      fprintf(stderr,"CheckLSScaffoldWierdnesses < 0 %s for scaffold " F_CID ", shifting by (%g,%g)...fixing...\n",
              string, scaffold->id,
              delta.mean, delta.variance);
    if (debug.leastSquaresGapsLV > 1)
      DumpCIScaffold(stderr,graph, scaffold, FALSE);

    minOffsetp->mean = minOffsetp->variance = 0.0;
    maxOffsetp->mean += delta.mean;
    maxOffsetp->variance += -delta.variance;
  }else if(minOffsetp->mean > 0.0){
    delta.mean = -minOffsetp->mean;
    delta.variance = -minOffsetp->variance;

    if (debug.leastSquaresGapsLV > 0)
      fprintf(stderr,"CheckLSScaffoldWierdnesses > 0 %s for scaffold " F_CID ", shifting by (%g,%g)...fixing...\n",
              string, scaffold->id,
              delta.mean, delta.variance);
    if (debug.leastSquaresGapsLV > 1)
      DumpCIScaffold(stderr,graph, scaffold, FALSE);

    minOffsetp->mean = minOffsetp->variance = 0.0;
    maxOffsetp->mean += delta.mean;
    maxOffsetp->variance += delta.variance;
  }else{
    return; // we are done
  }
  secondCI = NextCIScaffoldTIterator(&CIs);
  
  if(secondCI == NULL){
    scaffold->bpLength = (firstCI->offsetAEnd.mean < firstCI->offsetBEnd.mean) ?
      firstCI->offsetBEnd : firstCI->offsetAEnd;
    return;
  }
  
  AddDeltaToScaffoldOffsets(graph, scaffold->id,  secondCI->id, TRUE, FALSE, delta);
  
  if (debug.leastSquaresGapsLV > 1) {
    fprintf(stderr, "Done!! Scaffold after is:\n"); 
    DumpCIScaffold(stderr,graph, scaffold, FALSE);
  }
}




void LeastSquaresGapEstimates(ScaffoldGraphT *graph, int markEdges,
                              int useGuides, int forceNonOverlaps,
                              int checkConnectivity, int verbose){

  RecomputeOffsetsStatus status;
  int numScaffolds;
  int redo = FALSE;
  int cnt = 0;
  int sID;

  if (debug.leastSquaresGapsLV > 1)
    fprintf(stderr,"* Start of LeastSquaresGapEstimates() markEdges:%d useGuides:%d*\n",
            markEdges, useGuides);

  // if we are marking edges, we don't need to check now
  if(!markEdges)
    CheckAllTrustedEdges(graph);

  /*
    20050819 IMD
    Replaced iterator with incrementing index
    CheckScaffoldConnectivityAndSplit is called in the loop and
    can end up reallocating the scaffold graph, which can relocate
    the scaffold graph in memory, which will screw up the iterator
  */
  for(sID = 0;
      sID < (numScaffolds = GetNumCIScaffoldTs(graph->CIScaffolds));
      sID += (redo?0:1))
    {
      CIScaffoldT * scaffold = GetCIScaffoldT(graph->CIScaffolds, sID);

      if (debug.leastSquaresGapsLV > 1) {
        fprintf(stderr, "LeastSquaresGapEstimates() begins another iteration on scaffold %d, redo=%d\n", sID, redo);
        DumpCIScaffold(stderr, graph, scaffold, TRUE);
        dumpTrustedEdges(graph, scaffold, ALL_TRUSTED_EDGES);
      }

      if(isDeadCIScaffoldT(scaffold) || scaffold->type != REAL_SCAFFOLD){
        assert(!redo);
        continue;
      }

      if(++cnt % 10000 == 0)
        fprintf(stderr," LeastSquaresGapEstimates %d   scaffold " F_CID "/%d\n", cnt - 1, scaffold->id, numScaffolds);

      redo = FALSE;

      if(markEdges){
        MarkInternalEdgeStatus(graph, scaffold, PAIRWISECHI2THRESHOLD_CGW,
                               (useGuides ? (1000.0 * SLOPPY_EDGE_VARIANCE_THRESHHOLD) : SLOPPY_EDGE_VARIANCE_THRESHHOLD),
                               TRUE, TRUE, 0, TRUE);

        if (debug.leastSquaresGapsLV > 1)
          dumpTrustedEdges(graph, scaffold, ALL_TRUSTED_EDGES);
      }

      // Check that the scaffold is connected by trusted edges - otherwise
      // RecomputeOffsetsInScaffold will fail due to a singularity in the matrix
      // that it needs to invert - if it is not, break it into a set of maximal
      // scaffolds which are connected.
      //
      if(checkConnectivity){
        int numComponents = CheckScaffoldConnectivityAndSplit(graph, sID, ALL_TRUSTED_EDGES, verbose);

        if (debug.leastSquaresGapsLV > 1)
          fprintf(stderr, "* Scaffold "F_CID" has %d components.\n", sID, numComponents);

        if(numComponents > 1){ // we split the scaffold because it wasn't connected
          fprintf(stderr,"* Scaffold not connected: Split scaffold " F_CID " into %d pieces\n",
                  sID, numComponents);
          continue;
        }else{
          if (debug.leastSquaresGapsLV > 1)
            fprintf(stderr,"* BPW Scaffold connected " F_CID " hooray!\n",
                    sID, numComponents);

          if(!IsScaffold2EdgeConnected(ScaffoldGraph, scaffold)){
            if (debug.leastSquaresGapsLV > 0)
              fprintf(stderr,"*###### Scaffold " F_CID " is not 2-edge connected... SPLIT IT!\n",
                      sID);

            numComponents = CheckScaffoldConnectivityAndSplit(ScaffoldGraph, sID, ALL_TRUSTED_EDGES, FALSE);
            if (numComponents > 1) {
              if (debug.leastSquaresGapsLV > 0)
                fprintf(stderr,"* Scaffold not 2 edge-connected: Split scaffold " F_CID " into %d pieces\n",
                        sID, numComponents);
              continue;
            }
          }
        }
      }


      {
        int i;
        for(i = 0; i < 100; i++){
          CheckLSScaffoldWierdnesses("BEFORE", graph, scaffold);
          status =  RecomputeOffsetsInScaffold(graph, scaffold, TRUE, forceNonOverlaps, verbose);

          if(status == RECOMPUTE_CONTIGGED_CONTAINMENTS) {
            // We want to restart from the top of the loop, including edge marking
            if (debug.leastSquaresGapsLV > 1)
              fprintf(stderr, "RecomputeOffsetsInScaffold() returned RECOMPUTE_CONTIGGED_CONTAINMENTS, begin anew!\n");
            redo = TRUE;
            break;
          }

          if(status != RECOMPUTE_FAILED_REORDER_NEEDED) {
            if (debug.leastSquaresGapsLV > 1)
              fprintf(stderr, "RecomputeOffsetsInScaffold() returned OK!\n");
            break;
          }

          // We want to simply try again, since we just changed the scaffold order

          if (debug.leastSquaresGapsLV > 1)
            fprintf(stderr,"* RecomputeOffsetsInScaffold " F_CID " attempt %d failed...iterating\n", sID, i);
        }
      }


      CheckLSScaffoldWierdnesses("AFTER", graph, scaffold);


      if(status != RECOMPUTE_OK){
        if (debug.leastSquaresGapsLV > 0)
          fprintf(stderr, "RecomputeOffsetsInScaffold failed (%d) for scaffold " F_CID "\n", status, scaffold->id);

        //  If we are 'redo'ing stuff, we have just merged
        //  contigs.  The merge process does not rebuild trusted
        //  edges, so the scaffold might be disconnected.
        //  CheckCIScaffoldT() includes a call to
        //  RecomputeOffsetsInScaffold(), which needs a connected
        //  scaffold.  Thus, we skip it in this case.
        //
        if (!redo)
          CheckCIScaffoldT(graph, scaffold);
      }

      if ((sID % 100000) == 0)
        clearCacheSequenceDB(graph->sequenceDB);

    }  //  end of sID loop
  return;
}

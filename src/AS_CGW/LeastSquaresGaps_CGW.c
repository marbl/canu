
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
static char *rcsid = "$Id: LeastSquaresGaps_CGW.c,v 1.71 2012-08-31 20:20:56 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "ChiSquareTest_CGW.h"

#include "CIScaffoldT_Analysis.H"

#define FIXED_RECOMPUTE_NOT_ENOUGH_CLONES /* long standing bug: is it fixed yet? it seems to be */
#undef  FIXED_RECOMPUTE_NOT_ENOUGH_CLONES /* nope */

#undef NEG_GAP_VARIANCE_PROBLEM_FIXED  /* if undef'ed, allow processing to continue despite a negative gap variance */

#undef TEST_FIXUPMISORDER


/* declarations for LAPACK/DXML calls to linear algebra routines */
#define FTN_INT   long int
#define F_FTN_INT    "%ld"

extern "C" {
extern int dgemv_(char *, FTN_INT *, FTN_INT *,
                  double *, double *, FTN_INT *, double *, FTN_INT *,
                  double *, double *, FTN_INT *);

extern int dpbtrf_(char *, FTN_INT *, FTN_INT *, double *,
                   FTN_INT *, FTN_INT *);
extern int dpbtrs_(char *, FTN_INT *, FTN_INT *, FTN_INT *, double *,
                   FTN_INT *, double *, FTN_INT *, FTN_INT *);

extern int dgbtrf_(FTN_INT *m, FTN_INT *n, FTN_INT *kl, FTN_INT *ku, 
                   double *ab, FTN_INT *ldab, FTN_INT *ipiv, FTN_INT *info);
extern int dgbtrs_(char *trans, FTN_INT *n, FTN_INT *kl, FTN_INT *ku,
                   FTN_INT *nrhs, double *ab, FTN_INT *ldab, FTN_INT *ipiv, 
                   double *b, FTN_INT *ldb, FTN_INT *info);
}


typedef enum {
  RECOMPUTE_OK = 0,
  RECOMPUTE_SINGULAR = 1,
  RECOMPUTE_LAPACK = 2,
  RECOMPUTE_NO_GAPS = 3,
  RECOMPUTE_FAILED_REORDER_NEEDED = 4,
  RECOMPUTE_NOT_ENOUGH_CLONES = 5,
  RECOMPUTE_CONTIGGED_CONTAINMENTS = 6,
  RECOMPUTE_FAILED_CONTIG_DELETED  = 7
}RecomputeOffsetsStatus;



//   If we find a pair of contigs that are misordered, we rip out
//   thisCI, move it to the right place based on the implied position
//   in the overlapEdge
//
int
FixUpMisorderedContigs(CIScaffoldT           *scaffold,
                       ContigT               *prevCI,
                       ContigT               *thisCI,
                       PairOrient   edgeOrient,
                       double                 inferredMean,
                       double                 inferredVariance,
                       EdgeCGW_T             *overlapEdge){
  PairOrient newEdgeOrient = GetEdgeOrientationWRT(overlapEdge, prevCI->id);

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

  fprintf(stderr,"* edgeOrient %c   edge->orient = %c  newEdgeOrient = %c  prevCI = "F_CID"   thisCI = "F_CID" mean:%g\n",
          edgeOrient.toLetter(), overlapEdge->orient.toLetter(), newEdgeOrient.toLetter(),
          prevCI->id, thisCI->id, inferredMean);



  aEndOffset.mean = aEndOffset.variance = -1.0;
  bEndOffset.mean = bEndOffset.variance = -1.0;

  assert(edgeOrient.isUnknown() == false);

  if (edgeOrient.isAB_AB()) {
    assert(newEdgeOrient.isBA_BA());
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
  }

  if (edgeOrient.isAB_BA()) {
    assert(newEdgeOrient.isBA_AB());
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
  }

  if (edgeOrient.isBA_AB()) {
    assert(newEdgeOrient.isAB_BA());
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
  }

  if (edgeOrient.isBA_BA()) {
    assert(newEdgeOrient.isAB_AB());
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
  }

  fprintf(stderr,"* Overlap is ("F_CID","F_CID",%c)  moving "F_CID" from (%g,%g) to (%g,%g)\n",
          overlapEdge->idA,
          overlapEdge->idB,
          overlapEdge->orient.toLetter(),
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



EdgeCGW_T *
FindOverlapEdgeChiSquare(ScaffoldGraphT *graph,
                         NodeCGW_T *sourceCI,
                         CDS_CID_t targetId,
                         PairOrient edgeOrient,
                         double inferredMean,
                         double inferredVariance,
                         double *chiSquaredValue,
                         double chiSquareThreshold,
                         int *alternate, int verbose){
  EdgeCGW_T *bestEdge = NULL;
  double     bestChiSquaredValue = FLT_MAX;
  int32      end = (edgeOrient.isAB_AB() || edgeOrient.isAB_BA()) ? B_END : A_END;;

  *alternate = FALSE;

  GraphEdgeIterator edges(ScaffoldGraph->ContigGraph, sourceCI->id, end, ALL_EDGES);
  EdgeCGW_T        *edge;

  while ((edge = edges.nextRaw()) != NULL){
    CDS_CID_t otherCID = (edge->idA == sourceCI->id) ? edge->idB : edge->idA;

    if (otherCID != targetId)
      continue;

    if (isOverlapEdge(edge) == false)
      continue;

    if (isContainmentEdge(edge) == true)
      continue;

    if (GetEdgeOrientationWRT(edge, sourceCI->id) != edgeOrient)
      continue;

    if (PairwiseChiSquare(inferredMean, inferredVariance, edge->distance.mean, ((MAX_OVERLAP_SLOP_CGW * MAX_OVERLAP_SLOP_CGW) / 9),
                          NULL, chiSquaredValue, chiSquareThreshold) == false)
      continue;

    if ((bestEdge == NULL) || (*chiSquaredValue < bestChiSquaredValue)) {
      bestEdge = edge;
      bestChiSquaredValue = *chiSquaredValue;
    }
  }

  if (bestEdge != NULL)
    return(bestEdge);

  int32 minOverlap = MAX(CGW_MISSED_OVERLAP, -(inferredMean + (3.0 * sqrt(inferredVariance))));
  int32 maxOverlap = -(inferredMean - (3.0 * sqrt(inferredVariance)));

  if (maxOverlap < CGW_MISSED_OVERLAP)
    return(NULL);

  assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));

  ChunkOverlapCheckT olap = OverlapChunks(graph->ContigGraph,
                                          sourceCI->id, targetId,    // handles suspicious
                                          edgeOrient, minOverlap, maxOverlap,
                                          AS_CGW_ERROR_RATE, FALSE);
  double effectiveOlap = -olap.overlap;

  if(olap.suspicious){
    fprintf(stderr,"* FOEXS: SUSPICIOUS Overlap found! Looked for ("F_CID","F_CID",%c)["F_S32","F_S32"] found ("F_CID","F_CID",%c) "F_S32"\n",
            sourceCI->id, targetId, edgeOrient.toLetter(),
            minOverlap, maxOverlap,
            olap.spec.cidA, olap.spec.cidB,
            olap.spec.orientation.toLetter(), olap.overlap);
    effectiveOlap = -(GetGraphNode(ScaffoldGraph->ContigGraph, targetId)->bpLength.mean + sourceCI->bpLength.mean - olap.overlap);
  }

  if (olap.overlap == 0)
    return(NULL);

  edge = GetGraphEdge(graph->ContigGraph, MakeComputedOverlapEdge(graph->ContigGraph, &olap, TRUE));

  // Create an appropriate hash table entry
  CreateChunkOverlapFromEdge(ScaffoldGraph->ContigGraph, edge);

  if (PairwiseChiSquare(inferredMean, inferredVariance, effectiveOlap, ((MAX_OVERLAP_SLOP_CGW * MAX_OVERLAP_SLOP_CGW) / 9),
                        NULL, chiSquaredValue, chiSquareThreshold)) {
    *alternate = olap.suspicious;
    return(edge);
  }

  fprintf(stderr,"* Failed pairwise test between (%g, %g) and (%g,%g) not returning edge ("F_CID","F_CID",%c) %g\n",
          inferredMean, inferredVariance, effectiveOlap, (double) ((MAX_OVERLAP_SLOP_CGW * MAX_OVERLAP_SLOP_CGW) / 9),
          edge->idA, edge->idB, edge->orient.toLetter(), edge->distance.mean);

  return(NULL);
}




//  Dump trusted/raw edges
void
dumpTrustedEdges(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold) {
  ChunkInstanceT * chunk;
  CIScaffoldTIterator CIs;
  int set = 0;

  fprintf(stderr, "TRUSTED_EDGE in scaffold %d\n", scaffold->id);

  InitCIScaffoldTIterator(sgraph, scaffold, TRUE, FALSE, &CIs);

  while ((chunk = NextCIScaffoldTIterator(&CIs)) != NULL) {
    GraphEdgeIterator   edges(sgraph->ContigGraph, chunk->id, ALL_END, ALL_TRUSTED_EDGES);
    CIEdgeT            *edge;

    while ((edge = edges.nextRaw()) != NULL) {
      ChunkInstanceT * otherChunk = GetGraphNode(sgraph->ContigGraph,
                                                 (chunk->id == edge->idA) ?
                                                 edge->idB : edge->idA);

      int32 weight = edge->edgesContributing - (isOverlapEdge(edge));

      assert(otherChunk != NULL);

      // See each edge only once
      if(chunk->id != edge->idA)
        continue;

      if(edge->flags.bits.isBridge){
        fprintf(stderr,"* WARNING: chunk "F_CID" weight = %d bridge edge\n",
                chunk->id, weight);
        PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph,
                       "Bridge ", edge, chunk->id);
#ifdef DEBUG
        EdgeCGW_T *e;
        GraphEdgeIterator Edges(sgraph->ContigGraph,chunk->id,ALL_END, ALL_TRUSTED_EDGES);
        fprintf(stderr,"Edges out from "F_CID":\n",chunk->id);
        while(NULL!= (e = Edges.nextMerged()))
          PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph,
                         "DEBUG Bridge ",e, chunk->id);
#endif
      }
      if(isSingletonOverlapEdge(edge) ||
         (weight < MIN_EDGES && edge->flags.bits.isBridge))
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
  size_t sizeofGapCoefficientsAlt;
  size_t sizeofIPIV;
  size_t sizeofGapVariance;
  size_t sizeofCloneVariance;
  size_t sizeofCloneMean;
  size_t sizeofSpannedGaps;
  size_t sizeofGapSize;
  size_t sizeofGapSizeVariance;
  size_t sizeofGapsToComputeGaps;
  size_t sizeofComputeGapsToGaps;
  LengthT *lengthCIs;
  int32 *cloneGapStart;
  int32 *cloneGapEnd;
  double *gapConstants;
  double *gapCoefficients;
  double *gapCoefficientsAlt;
  FTN_INT *IPIV;
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
  safe_free(data->gapCoefficientsAlt);
  safe_free(data->IPIV);
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


static
void
dumpGapCoefficients(double *gapCoefficients,
                    int32   maxDiagonals,
                    int32   numComputeGaps,
                    int32   rows,
                    int32   bands) {
  int i = 0;
  double *gapEnd, *gapPtr;
  //for(gapPtr = gapConstants, gapEnd = gapPtr + numComputeGaps;
  //    gapPtr < gapEnd; gapPtr++){
  //  fprintf(stderr, "Gap Constants %g\n", *gapPtr);
  //}
  fprintf(stderr, "Gap Coefficients\n");
  for(gapPtr = gapCoefficients, gapEnd = gapPtr + (maxDiagonals * numComputeGaps);
      gapPtr < gapEnd; gapPtr++){
    fprintf(stderr, "%11.4g", *gapPtr);
    i++;
    if(i == maxDiagonals){
      fprintf(stderr, "\n");
      i = 0;
    }else{
      fprintf(stderr, "\t\t");
    }
  }
}


static
void
dumpGapCoefficientsAlt(double *gapCoefficientsAlt,
                       int32   maxDiagonals,
                       int32   numComputeGaps,
                       int32   rows,
                       int32   bands) {

  int i = 0;
  double *gapEnd, *gapPtr;
  //for(gapPtr = gapConstants, gapEnd = gapPtr + numComputeGaps;
  //    gapPtr < gapEnd; gapPtr++){
  //  fprintf(stderr, "Gap Constants %g\n", *gapPtr);
  //}
  fprintf(stderr, "Gap Coefficients Alternate\n");
  for(gapPtr = gapCoefficientsAlt, gapEnd = gapPtr + ((bands + bands + 1 + bands) * numComputeGaps); gapPtr < gapEnd; gapPtr++){
    fprintf(stderr, "%11.4g", *gapPtr);
    i++;
    if(i == bands + bands + 1 + bands){
      fprintf(stderr, "\n");
      i = 0;
    }else{
      fprintf(stderr, "\t");
    }
  }
}


void
dumpScaffoldContigPositions(ScaffoldGraphT *graph, CIScaffoldT *scaffold, char *label) {
  CIScaffoldTIterator      CIs;
  ChunkInstanceT          *CI;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  CI = NextCIScaffoldTIterator(&CIs);

  CDS_CID_t  prevCIid     = CI->id;
  LengthT    thisLeftEnd  = { 0, 0 };
  LengthT    thisRightEnd = { 0, 0 };
  LengthT    prevLeftEnd  = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;
  LengthT    prevRightEnd = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetBEnd : CI->offsetAEnd;

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    thisLeftEnd  = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;
    thisRightEnd = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetBEnd : CI->offsetAEnd;

    fprintf(stderr, "%s Contig %8d at %9.0f +- %-11.0f to %9.0f +- %-11.0f  ctg len %9.0f  gap to next %9.0f +- %-11.0f\n",
            label,
            prevCIid,
            prevLeftEnd.mean,  prevLeftEnd.variance,
            prevRightEnd.mean, prevRightEnd.variance,
            prevRightEnd.mean - prevLeftEnd.mean,
            thisLeftEnd.mean - prevRightEnd.mean,
            thisLeftEnd.variance - prevRightEnd.variance);

    prevCIid      = CI->id;
    prevLeftEnd   = thisLeftEnd;
    prevRightEnd  = thisRightEnd;
  }

  fprintf(stderr, "%s Contig %8d at %9.0f +- %-11.0f to %9.0f +- %-11.0f  ctg len %9.0f\n",
          label,
          prevCIid,
          prevLeftEnd.mean,  prevLeftEnd.variance,
          prevRightEnd.mean, prevRightEnd.variance,
          prevRightEnd.mean - prevLeftEnd.mean);
}



static
void
RebuildScaffoldGaps(ScaffoldGraphT  *graph,
                    CIScaffoldT     *scaffold,
                    double          *gapMean,
                    double          *gapVariance,
                    bool             beVerbose) {
  LengthT *prevLeftEnd, *prevRightEnd;
  LengthT *thisLeftEnd, *thisRightEnd;

  if (beVerbose)
    fprintf(stderr, "ROIS()-- Rebuild contig and gap positions for scaffold %d\n", scaffold->id);

  if (((gapMean == NULL) || (gapVariance == NULL)) &&
      (scaffold->info.Scaffold.numElements > 1)) {
    fprintf(stderr, "ROIS()-- NULL array!  mean %p variance %p\n", gapMean, gapVariance);
    dumpScaffoldContigPositions(graph, scaffold, "ROIS()-- NULL array!  ");
    return;
  }

  CIScaffoldTIterator CIs;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  NodeCGW_T *prevCI = NextCIScaffoldTIterator(&CIs);
  NodeCGW_T *thisCI = NULL;

  if (GetNodeOrient(prevCI).isForward()) {
    prevLeftEnd  = &prevCI->offsetAEnd;
    prevRightEnd = &prevCI->offsetBEnd;
  } else {
    prevLeftEnd  = &prevCI->offsetBEnd;
    prevRightEnd = &prevCI->offsetAEnd;
  }

  prevLeftEnd->mean     = 0.0;  //  Scaffold starts at zero, no?
  prevLeftEnd->variance = 0.0;

  scaffold->bpLength.mean     = 0.0;  //  Scaffold has no contigs yet, no?
  scaffold->bpLength.variance = 0.0;

  double   *gapPtr      = gapMean;
  double   *gapVarPtr   = gapVariance;

  while ((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL) { 
    if (GetNodeOrient(thisCI).isForward()) {
      thisLeftEnd  = &thisCI->offsetAEnd;
      thisRightEnd = &thisCI->offsetBEnd;
    } else {
      thisLeftEnd  = &thisCI->offsetBEnd;
      thisRightEnd = &thisCI->offsetAEnd;
    }

    double prevGap = thisLeftEnd->mean     - prevRightEnd->mean;
    double prevVar = thisLeftEnd->variance - prevRightEnd->variance;

    //  Set then end of the last contig.

    prevRightEnd->mean     = prevLeftEnd->mean      + prevCI->bpLength.mean;
    prevRightEnd->variance = prevLeftEnd->variance  + prevCI->bpLength.variance;

    if (scaffold->bpLength.mean     < prevRightEnd->mean)       scaffold->bpLength.mean     = prevRightEnd->mean;
    if (scaffold->bpLength.variance < prevRightEnd->variance)   scaffold->bpLength.variance = prevRightEnd->variance;

    if (beVerbose)
      fprintf(stderr, "ROIS()-- Contig %8d at %9.0f +- %-11.0f to %9.0f +- %-11.0f  ctg len %9.0f  gap to next %9.0f +- %-11.0f (was %9.0f +- %-11.0f)\n",
              prevCI->id,
              prevLeftEnd->mean,  prevLeftEnd->variance,
              prevRightEnd->mean, prevRightEnd->variance,
              prevRightEnd->mean - prevLeftEnd->mean,
              *gapPtr, *gapVarPtr,
              prevGap, prevVar);

    //  Adjust for implausible variances.

    if ((*gapPtr > 0) &&
        (*gapVarPtr <= 0)) {
      fprintf(stderr, "ROIS()-- Negative variance %f for gap size %f; reset to %f\n", *gapVarPtr, *gapPtr, (0.01 * *gapPtr) * (0.01 * *gapPtr));
      *gapVarPtr = (0.01 * *gapPtr) * (0.01 * *gapPtr);
    }

    if (*gapVarPtr > 50000.0 * 50000.0) {
      fprintf(stderr, "ROIS()-- Implausibly large variance %f for gap size %f; reset to %f\n", *gapVarPtr, *gapPtr, 50000.0 * 50000.0);
      *gapVarPtr = 50000.0 * 50000.0;
    }

    //  Now set the begin of the next contig.

    thisLeftEnd->mean      = prevRightEnd->mean     + *gapPtr;
    thisLeftEnd->variance  = prevRightEnd->variance + *gapVarPtr;

    prevCI       = thisCI;
    prevLeftEnd  = thisLeftEnd;
    prevRightEnd = thisRightEnd;

    gapPtr++;
    gapVarPtr++;
  }

  prevRightEnd->mean      = prevLeftEnd->mean     + prevCI->bpLength.mean;
  prevRightEnd->variance  = prevLeftEnd->variance + prevCI->bpLength.variance;

  if (scaffold->bpLength.mean     < prevRightEnd->mean)       scaffold->bpLength.mean     = prevRightEnd->mean;
  if (scaffold->bpLength.variance < prevRightEnd->variance)   scaffold->bpLength.variance = prevRightEnd->variance;

  if (beVerbose)
    fprintf(stderr, "ROIS()-- Contig %8d at %9.0f +- %-11.0f to %9.0f +- %-11.0f  ctg len %9.0f\n",
            prevCI->id,
            prevLeftEnd->mean,  prevLeftEnd->variance,
            prevRightEnd->mean, prevRightEnd->variance,
            prevRightEnd->mean - prevLeftEnd->mean);

  //dumpScaffoldContigPositions(graph, scaffold, "ROIS()--");  -  doesn't show previous gap size
}



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
void
LS_IncrementGapsCoveredByOneClone(
				  int thisCI_index, 
				  int otherCI_index,
				  int* _gapsToComputeGaps,
				  double* _gapConstants,
				  double* _gapCoefficients,
				  int32 _maxDiagonals,
				  double clone_constantMean,
				  double clone_inverseVariance) {
  //fprintf (stderr, "LS_IncrementGapsCoveredByOneClone(%d,%d,...%f,%f)\n",
  //         thisCI_index,otherCI_index,clone_constantMean,clone_inverseVariance);
  int colIndex;
  for(colIndex = thisCI_index;
      colIndex < otherCI_index; 
      colIndex++){
    int rowIndex, arrayIndex;
    int colComputeIndex = _gapsToComputeGaps[colIndex];
    /* For each gap that the clone spans it contributes the same
       constant value to the gapConstants vector which is equal
       to the mean total gap size for that clone divided by the
       variance of the total gap size. */
    if(colComputeIndex == NULLINDEX){
      continue;
    }
    _gapConstants[colComputeIndex] += clone_constantMean;

    for(rowIndex = colIndex;
	rowIndex < otherCI_index; 
	rowIndex++){
      int rowComputeIndex = _gapsToComputeGaps[rowIndex];
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
      arrayIndex = (colComputeIndex * _maxDiagonals) + (rowComputeIndex - colComputeIndex);
      _gapCoefficients [arrayIndex] += clone_inverseVariance;
    }
  }
}



static
RecomputeOffsetsStatus
RecomputeOffsetsInScaffold(ScaffoldGraphT *graph,
                           CDS_CID_t scaffoldID,
                           int allowOrderChanges,
                           int forceNonOverlaps,
                           int verbose){

  RecomputeData data;
  CIScaffoldTIterator CIs;
  /* Iterate over all of the "trusted" CIEdges */
  EdgeCGW_T *edge;
  NodeCGW_T *thisCI, *prevCI;
  int32 numCIs;
  int32 indexCIs;

  int numGaps, numComputeGaps;
  LengthT *lengthCIs, *lengthCIsPtr;
  int maxDiagonals = 1;
  int numClones = 0;
  CDS_CID_t indexClones;
  int32 *cloneGapStart, *cloneGapEnd;
  int *gapsToComputeGaps, *computeGapsToGaps;
  double *gapCoefficients, *gapCoefficientsAlt, *gapConstants;
  FTN_INT *IPIV;
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
  data.gapCoefficientsAlt = NULL;
  data.IPIV = NULL;
  data.gapVariance = NULL;
  data.cloneVariance = NULL;
  data.cloneMean = NULL;
  data.spannedGaps = NULL;
  data.gapSize = NULL;
  data.gapSizeVariance = NULL;
  data.gapsToComputeGaps = NULL;
  data.computeGapsToGaps = NULL;

  CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, scaffoldID);

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

    //fprintf(stderr, "Length of CI "F_CID" is %f indexInScaffold %d\n", thisCI->id, thisCI->bpLength.mean, thisCI->indexInScaffold);

    indexCIs++;
    lengthCIsPtr++;
  }
  assert(indexCIs == numCIs);

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL){
    GraphEdgeIterator  edges(ScaffoldGraph->ContigGraph, thisCI->id, ALL_END, ALL_TRUSTED_EDGES);

    while((edge = edges.nextRaw())!= NULL){
      int        isA = (edge->idA == thisCI->id);
      NodeCGW_T *otherCI = GetGraphNode(ScaffoldGraph->ContigGraph, (isA ? edge->idB : edge->idA));

      assert(thisCI->scaffoldID == otherCI->scaffoldID);
      assert(otherCI->indexInScaffold <= numGaps);
      assert(thisCI->indexInScaffold  <= numGaps);

      if (otherCI->scaffoldID == -1)
        continue;

      if (otherCI->indexInScaffold <= thisCI->indexInScaffold)
        continue; // Only interested in looking at an edge once

      numClones++;

      //fprintf(stderr,"RecomputeOffsets: adding clone between %d and %d\n",
      //        thisCI->id,otherCI->id);

      if((otherCI->indexInScaffold - thisCI->indexInScaffold) >
         maxDiagonals){
        maxDiagonals = otherCI->indexInScaffold - thisCI->indexInScaffold;

#if 0
        fprintf(stderr, "Max Diagonals %d (%d,%d) ["F_CID"."F_CID","F_CID"."F_CID"]\n",
                maxDiagonals, thisCI->indexInScaffold,
                otherCI->indexInScaffold, thisCI->scaffoldID,
                thisCI->id, otherCI->scaffoldID, otherCI->id);
        fprintf(stderr, "Max Diagonals %d (%d,%d) ["F_CID"."F_CID","F_CID"."F_CID"]\n",
                maxDiagonals, thisCI->indexInScaffold,
                otherCI->indexInScaffold, thisCI->scaffoldID,
                thisCI->id, otherCI->scaffoldID, otherCI->id);
#endif
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

    memset(gapConstants, 0, data.sizeofGapConstants);

    data.sizeofGapCoefficients = (maxDiagonals * numGaps) * sizeof(*gapCoefficients);
    data.gapCoefficients = gapCoefficients = (double *)safe_malloc(data.sizeofGapCoefficients);

    memset(gapCoefficients, 0, data.sizeofGapCoefficients);

    data.sizeofGapCoefficientsAlt = ((3 * maxDiagonals + 1) * numGaps) * sizeof(*gapCoefficientsAlt);
    data.gapCoefficientsAlt = gapCoefficientsAlt = (double *)safe_malloc(data.sizeofGapCoefficientsAlt);

    memset(gapCoefficientsAlt, 0, data.sizeofGapCoefficientsAlt);

    data.sizeofIPIV = (numGaps) * sizeof(*IPIV);
    data.IPIV = IPIV = (FTN_INT *)safe_malloc(data.sizeofIPIV);

    data.numClones = numClones;
    data.numGaps = numGaps;
    data.sizeofCloneMean = numClones * sizeof(*cloneMean);
    data.cloneMean = cloneMean = (double *)safe_malloc(data.sizeofCloneMean);

    data.sizeofCloneVariance = numClones * sizeof(*cloneVariance);
    data.cloneVariance = cloneVariance = (double *)safe_malloc(data.sizeofCloneVariance);

    data.sizeofGapVariance = numGaps * sizeof(*gapVariance);
    data.gapVariance = gapVariance = (double *)safe_malloc(data.sizeofGapVariance);

    memset(gapVariance, 0, data.sizeofGapVariance);

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

    for(numComputeGaps = 0; numComputeGaps < numGaps; numComputeGaps++)
      gapsToComputeGaps[numComputeGaps] = computeGapsToGaps[numComputeGaps] = numComputeGaps;
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
      GraphEdgeIterator edges(ScaffoldGraph->ContigGraph, thisCI->id, ALL_END, ALL_TRUSTED_EDGES);

      while((edge = edges.nextRaw())!= NULL){
        int isA = (edge->idA == thisCI->id);
        NodeCGW_T *otherCI = GetGraphNode(ScaffoldGraph->ContigGraph, (isA? edge->idB: edge->idA));

        assert(thisCI->scaffoldID == otherCI->scaffoldID);
        assert(otherCI->indexInScaffold <= numGaps);
        assert(thisCI->indexInScaffold  <= numGaps);

        double constantMean, constantVariance, inverseVariance;
        int lengthCIsIndex, gapIndex;

        if (otherCI->scaffoldID == -1)
          continue;

        if (otherCI->indexInScaffold <= thisCI->indexInScaffold)
          continue; // Only interested in looking at an edge once

        if(indexClones>=numClones){
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
        for(constantMean = edge->distance.mean,
              constantVariance = edge->distance.variance,
              lengthCIsIndex = thisCI->indexInScaffold + 1;
            lengthCIsIndex < otherCI->indexInScaffold; lengthCIsIndex++){
          constantMean -= lengthCIs[lengthCIsIndex].mean;
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
          assert(gapIndex < numGaps);
          if(gapsToComputeGaps[gapIndex] == NULLINDEX){
            constantMean -= gapSize[gapIndex];
            //constantVariance += gapSizeVariance[gapIndex];
          }
        }
        /* cloneMean and cloneVariance are the statistics for the estimated
           total size of the gaps spanned by this clone. */
        cloneMean[indexClones] = constantMean;
        cloneVariance[indexClones] = constantVariance;

        //fprintf(stderr, "Gap clone %f,%f (%d,%d)\n",
        //        constantMean, sqrt(constantVariance), thisCI->indexInScaffold, otherCI->indexInScaffold);

        constantMean /= constantVariance;
        inverseVariance = 1.0 / constantVariance;
        /* Store which gaps each clone spans so that we can iterate over
           these gaps when we calculate the gap variances and the
           squared error. */
        cloneGapStart[indexClones] = thisCI->indexInScaffold;
        cloneGapEnd[indexClones] = otherCI->indexInScaffold;

	LS_IncrementGapsCoveredByOneClone
	  (thisCI->indexInScaffold, 
	   otherCI->indexInScaffold,
	   gapsToComputeGaps,
	   gapConstants,
	   gapCoefficients,
	   maxDiagonals,
	   constantMean,
	   inverseVariance);
        indexClones++;
      }
    }

    maxClone=indexClones;

    bool    isCholesky = true;
    FTN_INT bands      = maxDiagonals - 1;
    FTN_INT rows       = numComputeGaps;

    //dpbtrf - Computes the Cholesky factorization of a symmetric/Hermitian positive definite band matrix
    //dpbtrs - Solves a symmetric/Hermitian positive definite banded system of linear equations AX=B, using the Cholesky factorization computed by SPBTRF/CPBTRF
    //
    //dgbtrf - Computes an LU factorization of a general band matrix, using partial pivoting with row interchanges
    //dgbtrs - Solves a general banded system of linear equations AX=B, A**T X=B or A**H X=B, using the LU factorization computed by SGBTRF/CGBTRF

    if (isCholesky == true) {
      FTN_INT ldab = maxDiagonals;
      FTN_INT info = 0;

      //      fprintf(stderr, "Calling dpbtrf() \n");
      //      dumpGapCoefficients(gapCoefficients, maxDiagonals, numComputeGaps, rows, bands); // debug

      dpbtrf_("L", &rows, &bands, gapCoefficients, &ldab, &info);

      //      dumpGapCoefficients(gapCoefficients, maxDiagonals, numComputeGaps, rows, bands); // debug
      //      fprintf(stderr, "dpbtrf: ldab "F_FTN_INT" info "F_FTN_INT"\n", ldab, info);

      if (info < 0) {
        //  The -info'th argument had an illegal value.
        fprintf(stderr, "dpbtrf failed; arg "F_FTN_INT" is illegal.\n", -info);
        isCholesky = false;
      }
      assert(info >= 0);

      if (info > 0) {
        //  Leading minor of order 'info' is not positive definite; factorization could not be completed.
        fprintf(stderr, "dpbtrf failed with info="F_FTN_INT"; no solution found, giving up.\n", info);
        isCholesky = false;
      }
    }

#define LU_BUSTED  // debug
#ifdef LU_BUSTED
    if (isCholesky == false)
      return(RECOMPUTE_LAPACK);
#endif

    //  Force LU
    //isCholesky = false;

    if (isCholesky == false) {
      FTN_INT ldab = maxDiagonals-1 + maxDiagonals-1 + 1 +maxDiagonals-1;
      FTN_INT info = 0;

      //   fprintf(stderr, "Calling dgbtrf() \n");
      //   dumpGapCoefficientsAlt(gapCoefficientsAlt, maxDiagonals, numComputeGaps, rows, bands);  // debug

      dgbtrf_(&rows, &rows, &bands, &bands, gapCoefficientsAlt, &ldab, IPIV, &info);

      //    dumpGapCoefficientsAlt(gapCoefficientsAlt, maxDiagonals, numComputeGaps, rows, bands); // debug
      //    fprintf(stderr, "dgbtrf: ldab "F_FTN_INT" info "F_FTN_INT"\n", ldab, info);

      if (info < 0) {
        //  The -info'th argument had an illegal value.
        fprintf(stderr, "dgbtrf failed; arg "F_FTN_INT" is illegal.\n", -info);
      }
      assert(info >= 0);

      if (info > 0) {
        fprintf(stderr, "dgbtrf failed with info="F_FTN_INT"; a singularity will result in divide by zero, giving up.\n", info);
        freeRecomputeData(&data);
        return(RECOMPUTE_SINGULAR);
      }
    }

    //  multiply the inverse of the gapCoefficients matrix by the gapConstants vector resulting in
    //  the least squares minimal solution of the gap sizes being returned in the gapConstants
    //  vector.

    if (isCholesky) {
      FTN_INT ldab = maxDiagonals;
      FTN_INT info = 0;
      FTN_INT nrhs = 1;

      dpbtrs_("L", &rows, &bands, &nrhs, gapCoefficients, &ldab, gapConstants, &rows, &info);
      assert(info == 0);
    } else {
      FTN_INT ldab = maxDiagonals-1 + maxDiagonals-1 + 1 +maxDiagonals-1;
      FTN_INT info = 0;
      FTN_INT nrhs = 1;

      dgbtrs_("N", &rows, &bands, &bands, &nrhs, gapCoefficientsAlt, &ldab, IPIV, gapConstants, &rows, &info);
      assert(info == 0);
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

        if (isCholesky) {
          FTN_INT ldab = maxDiagonals;
          FTN_INT info = 0;
          FTN_INT nrhs = 1;

          dpbtrs_("L", &rows, &bands, &nrhs, gapCoefficients, &ldab, spannedGaps, &rows, &info);
          assert(info == 0);
        } else {
          FTN_INT ldab = maxDiagonals-1 + maxDiagonals-1 + 1 +maxDiagonals-1;
          FTN_INT info = 0;
          FTN_INT nrhs = 1;

          dgbtrs_("N", &rows, &bands, &bands, &nrhs, gapCoefficientsAlt, &ldab, IPIV, spannedGaps, &rows, &info);
          assert(info == 0);
        }

        /* According to equation 5-7 we need to square the derivative and
           multiply by the total gap variance for this clone but instead
           we end up dividing by the total gap variance for this clone
           because we neglected to divide by it before squaring and so
           the net result is to need to divide by it. */
        for (double *gapPtr = spannedGaps,
                    *gapEnd = gapPtr + numComputeGaps,
                    *gapPtr2 = gapVariance; gapPtr < gapEnd; gapPtr++, gapPtr2++)
          *gapPtr2 += (*gapPtr) * (*gapPtr) / cloneVariance[indexClones];
      }
    }

    {
      int gapIndex, computeGapIndex;
      for(gapIndex = 0; gapIndex < numComputeGaps; gapIndex++){
        gapSize[computeGapsToGaps[gapIndex]] = gapConstants[gapIndex];
        gapSizeVariance[computeGapsToGaps[gapIndex]] = gapVariance[gapIndex];
        //fprintf(stderr,"GapSize(%d:%d) %f:%f\n", gapIndex,
        //        computeGapsToGaps[gapIndex], gapConstants[gapIndex],
        //        sqrt(gapVariance[gapIndex]));
      }

      if(forceNonOverlaps){
        int alternate;

        hardConstraintSet = FALSE;
        InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
        for(gapIndex = 0, computeGapIndex = 0,
              prevCI = NextCIScaffoldTIterator(&CIs);
            (thisCI = NextCIScaffoldTIterator(&CIs)) != NULL;
            prevCI = thisCI, gapIndex++){
          PairOrient edgeOrient;
          double chiSquaredValue;
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
          if(GetNodeOrient(thisCI).isForward()){
            if(GetNodeOrient(prevCI).isForward()){
              edgeOrient.setIsAB_AB();
            }else{//GetNodeOrient(prevCI) == B_A
              edgeOrient.setIsBA_AB();
            }
          }else{//GetNodeOrient(thisCI) == B_A
            if(GetNodeOrient(prevCI).isForward()){
              edgeOrient.setIsAB_BA();
            }else{//GetNodeOrient(prevCI) == B_A
              edgeOrient.setIsBA_BA();
            }
          }
          overlapEdge = FindOverlapEdgeChiSquare(graph, prevCI, thisCI->id,
                                                 edgeOrient, gapSize[gapIndex],
                                                 gapSizeVariance[gapIndex],
                                                 &chiSquaredValue,
                                                 (double)PAIRWISECHI2THRESHOLD_CGW,
                                                 &alternate, verbose);

         if(GlobalData->removeNonOverlapingContigsFromScaffold && overlapEdge == (EdgeCGW_T *)NULL){
            // try to find next contig to join so we can skip this one
            // can call findoverlapedge on next node
            // until we run out of nodes or we find one
            // if we find one, throw out the old contig and return recompute
            CDS_CID_t    startingCI          = thisCI->id;
            NodeCGW_T   *skippedCI           = thisCI;
            int32        startingGapIndex    = gapIndex;
            double        currGapSize         = gapSize[gapIndex];
            double        currGapSizeVariance = gapSizeVariance[gapIndex];
            
            while (overlapEdge == (EdgeCGW_T *)NULL) {
               currGapSize += thisCI->bpLength.mean + gapSize[gapIndex+1];
               currGapSizeVariance += thisCI->bpLength.variance + gapSizeVariance[gapIndex+1];
               gapIndex++; 
               thisCI = NextCIScaffoldTIterator(&CIs);
               if (thisCI == NULL || currGapSize > - CGW_MISSED_OVERLAP) {
                  break;
               }
               edgeOrient = GetChunkPairOrientation(GetNodeOrient(prevCI), GetNodeOrient(thisCI));
               overlapEdge = FindOverlapEdgeChiSquare(graph, prevCI, thisCI->id,
                                                      edgeOrient, currGapSize,
                                                      currGapSizeVariance,
                                                      &chiSquaredValue,
                                                      (double)PAIRWISECHI2THRESHOLD_CGW,
                                                      &alternate, verbose);
               if (overlapEdge) {
                  edgeOrient = GetChunkPairOrientation(GetNodeOrient(skippedCI), GetNodeOrient(thisCI));
                  EdgeCGW_T *skippedToNextEdge = FindOverlapEdgeChiSquare(graph, skippedCI, thisCI->id,
                                                         edgeOrient, gapSize[startingGapIndex+1],
                                                         gapSizeVariance[startingGapIndex+1],
                                                         &chiSquaredValue,
                                                         (double)PAIRWISECHI2THRESHOLD_CGW,
                                                         &alternate, verbose);
                  if (skippedToNextEdge != NULL)
                     overlapEdge = NULL;
               }
            }
            // restore our pre-search positions
            InitCIScaffoldTIteratorFromCI(graph, scaffold, startingCI, TRUE, FALSE, &CIs);
            gapIndex = startingGapIndex;
         
            // if we did find an overlap to a subsequent neighbor, throw intermediate contigs out and rerun 
            // otherwise advance past the contig analyzed above and continue
            if (overlapEdge) {
               NodeCGW_T  *toDelete        = NULL;
               CDS_CID_t   newScaffoldID   = NULLINDEX;
               uint32      seenFirstOffset = FALSE;
               LengthT     firstOffset     = {0.0, 0.0};
               LengthT     offsetAEnd      = {0.0, 0.0};
               LengthT     offsetBEnd      = {0.0, 0.0};
               
               while((toDelete = NextCIScaffoldTIterator(&CIs)) != NULL && toDelete->id != thisCI->id){
                  if(!seenFirstOffset){
                     CIScaffoldT CIScaffold;
                     InitializeScaffold(&CIScaffold, REAL_SCAFFOLD);
                     CIScaffold.info.Scaffold.AEndCI = NULLINDEX;
                     CIScaffold.info.Scaffold.BEndCI = NULLINDEX;
                     CIScaffold.info.Scaffold.numElements = 0;
                     CIScaffold.bpLength = firstOffset;
                     newScaffoldID = CIScaffold.id = GetNumGraphNodes(graph->ScaffoldGraph);
                     CIScaffold.flags.bits.isDead = FALSE;
                     CIScaffold.numEssentialA = CIScaffold.numEssentialB = 0;
                     CIScaffold.essentialEdgeB = CIScaffold.essentialEdgeA = NULLINDEX;

                     AppendGraphNode(graph->ScaffoldGraph, &CIScaffold);

                     //  Ensure that there are no edges, and that the edgeList is allocated.
                     assert(ScaffoldGraph->ScaffoldGraph->edgeLists[CIScaffold.id].empty() == true);
                               
                     scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, scaffoldID);

                     if(GetNodeOrient(toDelete).isForward()){
                       firstOffset = toDelete->offsetAEnd;
                     }else{
                       firstOffset = toDelete->offsetBEnd;
                     }
                     seenFirstOffset = TRUE;
                   }
                   offsetAEnd.mean     = toDelete->offsetAEnd.mean     - firstOffset.mean;
                   offsetAEnd.variance = toDelete->offsetAEnd.variance - firstOffset.variance;
                   offsetBEnd.mean     = toDelete->offsetBEnd.mean     - firstOffset.mean;
                   offsetBEnd.variance = toDelete->offsetBEnd.variance - firstOffset.variance;
         
                   RemoveCIFromScaffold(graph, scaffold, toDelete, FALSE);
                   InsertCIInScaffold(graph, toDelete->id, newScaffoldID, offsetAEnd, offsetBEnd, TRUE, FALSE);
                   fprintf(stderr, "KickOutNonOverlappingContig: Removing contig %d to scaffold %d because we found no overlaps to it\n", toDelete->id, newScaffoldID);
               }
               freeRecomputeData(&data);
               return(RECOMPUTE_FAILED_CONTIG_DELETED);
            } else {
               thisCI = NextCIScaffoldTIterator(&CIs);
            }
         }
         
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

            //fprintf(stderr,"*** Least Squares found the following alternate edge...contig NOW!\n");
            //PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, " alternateEdge: ", overlapEdge, overlapEdge->idA);
            //fprintf(stderr, "BEFORE ContigContainment, scaffold "F_CID" %s connected\n",
            //        scaffold->id,
            //        IsScaffoldInternallyConnected(ScaffoldGraph, scaffold, ALL_TRUSTED_EDGES) ? "is" : "is NOT");
            //dumpTrustedEdges(ScaffoldGraph, scaffold);
            //DumpACIScaffold(stderr, graph, scaffold, FALSE);

            // We want to merge the two contigs immediately, since these are problematic,
            // but we know we want to overlap these guys.
            //
            if (ContigContainment(scaffold, prevCI, thisCI, overlapEdge, TRUE) != TRUE) {
              fprintf(stderr, "ContigContainment failed.\n");
              assert(0);
            }

            //fprintf(stderr, "AFTER ContigContainment, scaffold "F_CID" %s connected\n",
            //        scaffold->id,
            //        IsScaffoldInternallyConnected(ScaffoldGraph, scaffold, ALL_TRUSTED_EDGES) ? "is" : "is NOT");
            //dumpTrustedEdges(ScaffoldGraph, scaffold);
            //DumpACIScaffold(stderr, graph, scaffold, FALSE);

            freeRecomputeData(&data);
            return RECOMPUTE_CONTIGGED_CONTAINMENTS;
          }
          if(overlapEdge && isContainmentEdge(overlapEdge)){

            //fprintf(stderr,"*** Least Squares found the following containment edge...contig NOW!\n");
            //fprintf(stderr,"*** "F_CID" (length %g) should be contained within "F_CID" (length %g)\n",
            //        thisCI->id, thisCI->bpLength.mean,
            //        prevCI->id, prevCI->bpLength.mean);
            //PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, " overlapEdge: ", overlapEdge, overlapEdge->idA);
            
            //fprintf(stderr, "BEFORE ContigContainment, scaffold "F_CID" %s connected\n",
            //        scaffold->id,
            //        IsScaffoldInternallyConnected(ScaffoldGraph, scaffold, ALL_TRUSTED_EDGES) ? "is" : "is NOT");
            //dumpTrustedEdges(ScaffoldGraph, scaffold);
            //DumpACIScaffold(stderr, graph, scaffold, FALSE);

            // When we find a containment relationship in a scaffold we want to merge the two contigs
            // immediately, since containments have the potential to induce situations that are confusing
            // for least squares
            if (ContigContainment(scaffold, prevCI, thisCI, overlapEdge, TRUE) != TRUE) {
              fprintf(stderr, "ContigContainment failed.\n");
              assert(0);
            }

            //fprintf(stderr, "AFTER ContigContainment, scaffold "F_CID" %s connected\n",
            //        scaffold->id,
            //        IsScaffoldInternallyConnected(ScaffoldGraph, scaffold, ALL_TRUSTED_EDGES) ? "is" : "is NOT");
            //dumpTrustedEdges(ScaffoldGraph, scaffold);
            //DumpACIScaffold(stderr, graph, scaffold, FALSE);

            freeRecomputeData(&data);
            return RECOMPUTE_CONTIGGED_CONTAINMENTS;
          }


          if(overlapEdge == (EdgeCGW_T *)NULL){
            // gapsize is negative, so is leftedge
            double leftEdge = (gapSize[gapIndex] - 3 * sqrt(gapSizeVariance[gapIndex]));
            double newLeftEdge = leftEdge;
            double newStd = (newLeftEdge - (-CGW_MISSED_OVERLAP))/3.0;
            double newVariance = newStd * newStd;

            //fprintf(stderr,"GapChange Gap(%d:%d) CIs: "F_CID","F_CID" new:(%f,%f) old(%f,%f)\n",
            //        gapsToComputeGaps[gapIndex], gapIndex,
            //        prevCI->id, thisCI->id,
            //        (double)(- CGW_MISSED_OVERLAP), newVariance,
            //        gapSize[gapIndex], gapSizeVariance[gapIndex]);

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
          if(fabs(overlapEdge->distance.mean - gapSize[gapIndex]) >
             MAX_OVERLAP_SLOP_CGW){
            computeGapsToGaps[computeGapIndex] = gapIndex;
            gapsToComputeGaps[gapIndex] = computeGapIndex;
            computeGapIndex++;
            continue;
          }

          //fprintf(stderr,"GapChange(%d:%d) %f:%f\n",
          //        gapsToComputeGaps[gapIndex], gapIndex,
          //        overlapEdge->distance.mean, gapSize[gapIndex]);

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

  //fprintf(stderr,"LSE: %f,%f #clones: %d,%d\n",
  //        scaffold->info.Scaffold.leastSquareError, squaredError,
  //        scaffold->info.Scaffold.numLeastSquareClones, numClones);

  scaffold->info.Scaffold.leastSquareError = squaredError;
  scaffold->info.Scaffold.numLeastSquareClones = numClones;

  RebuildScaffoldGaps(graph, scaffold, gapSize, gapSizeVariance, false);

  freeRecomputeData(&data);
  return (RECOMPUTE_OK);
}



bool
AdjustNegativePositions(ScaffoldGraphT *graph, CIScaffoldT *scaffold) {
  CIScaffoldTIterator      CIs;
  ChunkInstanceT          *CI;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  LengthT  minPos = { DBL_MAX, DBL_MAX };

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    LengthT thisLeftEnd  = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;

    if (thisLeftEnd.mean < minPos.mean)
      minPos = thisLeftEnd;
  }

  if (minPos.mean == 0.0)
    return(false);

  //dumpScaffoldContigPositions(graph, scaffold, "AdjPos(pre)--");

  fprintf(stderr, "AdjNegPos()-- adjust scaffold start by %f +- %f\n", minPos.mean, minPos.variance);

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    CI->offsetAEnd.mean     -= minPos.mean;
    CI->offsetAEnd.variance -= minPos.variance;
    CI->offsetBEnd.mean     -= minPos.mean;
    CI->offsetBEnd.variance -= minPos.variance;
  }

  //dumpScaffoldContigPositions(graph, scaffold, "AdjPos(post)--");

  return(true);
}


//  Update the variances.  This isn't even close to being correct, but it's the best we have.
//
//  In BounceOutOfOrderContigs(), the contig variances are set to zero.  This looks for zero
//  variances, and sets them to the average of the surrounding valid positions.
//  
static
bool
BounceOutOfOrderSort(LengthT const *a, LengthT const *b) {
  return(a->mean < b->mean);
}

bool
AdjustNegativeVariances(ScaffoldGraphT *graph, CIScaffoldT *scaffold) {
  CIScaffoldTIterator      CIs;
  ChunkInstanceT          *CI;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  bool                       zeroFound = false;

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    if (((CI->offsetAEnd.variance <= 0.0) && (CI->offsetAEnd.mean != 0.0)) ||
        ((CI->offsetBEnd.variance <= 0.0) && (CI->offsetBEnd.mean != 0.0))) {
      zeroFound = true;
      break;
    }
  }

  if (zeroFound == false)
    return(false);

  //dumpScaffoldContigPositions(graph, scaffold, "AdjVar(pre)--");

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  vector<LengthT *>          positions;

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    positions.push_back(&CI->offsetAEnd);
    positions.push_back(&CI->offsetBEnd);
  }

  sort(positions.begin(), positions.end(), BounceOutOfOrderSort);

  for (uint32 pi=1; pi<positions.size()-1; pi++) {
    if (positions[pi]->variance > 0.0)
      continue;

    uint32  pa = pi - 1;  //  Should always have good variance, even if it is an average
    uint32  pb = pi + 1;  //  Might have zero variance

    while ((pb < positions.size()) && (positions[pb]->variance <= 0.0))
      pb++;

    //  The -20 gaps screw this up.  It is close, and probably better than leaving it alone.

    double  adj = 0;

    if (pb == positions.size())
      adj = 500 * 500;
    else
      adj = (positions[pa]->variance + positions[pb]->variance) / (pb - pa);

    //  Adjust all variances from pi to pb.

    for (; pi < pb; pi++) {
      assert(positions[pi]->variance <= 0);

      positions[pi]->variance = positions[pi-1]->variance + adj;

      if (pb == positions.size())
        fprintf(stderr, "AdjNegVar()-- adjust variance for position %f to %f (actual %f) based on %f +- %f\n",
                positions[pi]->mean, positions[pi]->variance, adj,
                positions[pa]->mean, positions[pa]->variance);
      else
        fprintf(stderr, "AdjNegVar()-- adjust variance for position %f to %f (actual %f) based on %f +- %f and %f +- %f and %d gaps\n",
                positions[pi]->mean, positions[pi]->variance, adj,
                positions[pa]->mean, positions[pa]->variance,
                positions[pb]->mean, positions[pb]->variance,
                pb - pa);
    }
  }

  //  Last pass is to reset the end contig variance to something reasonable.

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    LengthT *thisMin = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? &CI->offsetAEnd : &CI->offsetBEnd;
    LengthT *thisMax = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? &CI->offsetBEnd : &CI->offsetAEnd;

    thisMax->variance = thisMin->variance + ComputeFudgeVariance(thisMax->mean - thisMin->mean);
  }

  //dumpScaffoldContigPositions(graph, scaffold, "AdjVar(post)--");

  return(true);
}



void
AdjustNegativeGapVariances(ScaffoldGraphT *graph, CIScaffoldT *scaffold) {

  CIScaffoldTIterator  CIs;
  ChunkInstanceT      *CI;

  LengthT lastMax     = {0, 0};

  double  varAdj      = 0.0;
  double  varOff      = 0.0;
  bool    fixesNeeded = false;

#undef LOG_NEG_GAP_FIX
#ifdef LOG_NEG_GAP_FIX
  //  If no fixes, don't bother logging.
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    LengthT thisMin = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;
    LengthT thisMax = (CI->offsetAEnd.mean > CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;

    if (lastMax.mean < thisMin.mean)
      if (lastMax.variance >= thisMin.variance)
        fixesNeeded = true;

    lastMax.mean     = MAX(lastMax.mean,     thisMax.mean);
    lastMax.variance = MAX(lastMax.variance, thisMax.variance);
  }

  if (fixesNeeded == false)
    return;
#endif

  lastMax.mean     = 0.0;
  lastMax.variance = 0.0;

#ifdef LOG_NEG_GAP_FIX
  dumpScaffoldContigPositions(graph, scaffold, "AdjGapVar(pre)--");
#endif

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    LengthT *thisMin = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? &CI->offsetAEnd : &CI->offsetBEnd;
    LengthT *thisMax = (CI->offsetAEnd.mean > CI->offsetBEnd.mean) ? &CI->offsetAEnd : &CI->offsetBEnd;

    //  Adjust variances for this contig.
    thisMin->variance += varAdj;
    thisMax->variance += varAdj;

    //  Do we need more adjustment?
    if ((lastMax.mean < thisMin->mean) &&
        (lastMax.variance + 1 >= thisMin->variance)) {

      //  The final gap variance
      double  newVariance = 0.10 * (thisMin->mean - lastMax.mean) * 0.10 * (thisMin->mean - lastMax.mean);

      fprintf(stderr, "AdjGapVar()--  contig %d in scaffold %d -- negative gap variance %f on positive gap size %f -- reset to %f\n",
              CI->id, scaffold->id, thisMin->variance - lastMax.variance, thisMin->mean - lastMax.mean, newVariance);

      //  But we also need to account for the currently negative variance for this gap.
      newVariance += lastMax.variance - thisMin->variance;

      //  Adjust again.
      thisMin->variance += newVariance;
      thisMax->variance += newVariance;

      varAdj += newVariance;
    }

    //  Save the highest we've seen.
    lastMax.mean     = MAX(lastMax.mean,     thisMax->mean);
    lastMax.variance = MAX(lastMax.variance, thisMax->variance);
  }

  //  Reset scaffold length
  scaffold->bpLength.mean     = lastMax.mean;
  scaffold->bpLength.variance = lastMax.variance;

#ifdef LOG_NEG_GAP_FIX
  dumpScaffoldContigPositions(graph, scaffold, "AdjGapVar(post)--");
#endif
}



//  This bit of bizarreness was in CheckCIScaffoldT(), which was dropped in favor of
//  ScaffoldSanity().  It serves to remove apparent out of order contigs, and reinsert them.  The
//  reinsert supposedly orders contigs correctly.
//
//  Unlike the original version, this one removes and reinserts ALL contigs.  The original tried to
//  detect which ones were out of order.  The original would iterate up to 20 times detecting and
//  reinserting contigs.
//
static
bool
BounceOutOfOrderContigs(ScaffoldGraphT *graph, CIScaffoldT *scaffold) {
  double                   lastMin = 0.0;

  CIScaffoldTIterator      CIs;
  ChunkInstanceT          *CI;  //  Current

  //  Test if anything is out of order.

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    double thisMin = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetAEnd.mean : CI->offsetBEnd.mean;

    if (thisMin < lastMin)
      break;

    lastMin = thisMin;
  }

  //  Anything?

  if (CI == NULL)
    return(false);

  //  Populate an array with the contigs in the current order.  We'll search through it to find
  //  those that are out of order.

  dumpScaffoldContigPositions(graph, scaffold, "BOOOC(pre)--");

  vector<ChunkInstanceT *> contigList;
  vector<ChunkInstanceT *> bounceList;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs)))
    contigList.push_back(CI);

  //  Search for out of order.  In order is: prev < this < next
  //
  //  Out of order is thus
  //     prev >= this or
  //     this >= next

  for (uint32 i=0; i<contigList.size(); ) {
    double prevMin = -DBL_MAX;
    double thisMin = (contigList[i]->offsetAEnd.mean < contigList[i]->offsetBEnd.mean) ? contigList[i]->offsetAEnd.mean : contigList[i]->offsetBEnd.mean;
    double nextMin = DBL_MAX;

    if (i > 0)
      prevMin = (contigList[i-1]->offsetAEnd.mean < contigList[i-1]->offsetBEnd.mean) ? contigList[i-1]->offsetAEnd.mean : contigList[i-1]->offsetBEnd.mean;

    if (i < contigList.size() - 1)
      nextMin = (contigList[i+1]->offsetAEnd.mean < contigList[i+1]->offsetBEnd.mean) ? contigList[i+1]->offsetAEnd.mean : contigList[i+1]->offsetBEnd.mean;

    if ((prevMin >= thisMin) ||
        (thisMin >= nextMin)) {
      bounceList.push_back(contigList[i]);
      contigList.erase(contigList.begin() + i);

    } else {
      i++;
    }
  }

  //  Eject and reinsert them.  This will place it in the contig list at the correct location.

  for (uint32 bi=0; bi<bounceList.size(); bi++) {
    CI = bounceList[bi];

    fprintf(stderr, "BounceOutOfOrderContigs()-- Bounce contig %d in scaffold %d\n", CI->id, scaffold->id);

    RemoveCIFromScaffold(graph, scaffold, CI, FALSE);

    CI->offsetAEnd.variance = 0.0;
    CI->offsetBEnd.variance = 0.0;

    InsertCIInScaffold(graph, CI->id, scaffold->id, CI->offsetAEnd, CI->offsetBEnd, TRUE, FALSE);
  }

  //  Report the damage

  dumpScaffoldContigPositions(graph, scaffold, "BOOOC(post)--");

  return(true);
}





bool
LeastSquaresGapEstimates(ScaffoldGraphT *graph,
                         CIScaffoldT    *scaffold,
                         uint32          LSFlags, 
                         uint32          bounceIteration) {
  RecomputeOffsetsStatus   status = RECOMPUTE_SINGULAR;

  if ((isDeadCIScaffoldT(scaffold)) ||
      (scaffold->type != REAL_SCAFFOLD))
    return(false);

  if (scaffold->info.Scaffold.numElements == 1) {
    RebuildScaffoldGaps(graph, scaffold, NULL, NULL, false);
    return(true);
  }

  int32  sID = scaffold->id;

  for (uint32 iter=0; iter<5; iter++) {

    //  Even though we only use raw edges, still mark the merged edges.

    MarkInternalEdgeStatus(graph, scaffold, 0, TRUE,  PAIRWISECHI2THRESHOLD_CGW, SLOPPY_EDGE_VARIANCE_THRESHHOLD);  //  Merged
    MarkInternalEdgeStatus(graph, scaffold, 0, FALSE, PAIRWISECHI2THRESHOLD_CGW, SLOPPY_EDGE_VARIANCE_THRESHHOLD);  //  Raw

    //  Don't check variance (false = use raw, true = use trusted)
#if 1
    if (IsScaffoldInternallyConnected(ScaffoldGraph, scaffold, false, true) != 1) {
      MarkInternalEdgeStatus(graph, scaffold, 1, TRUE,  PAIRWISECHI2THRESHOLD_CGW, SLOPPY_EDGE_VARIANCE_THRESHHOLD);  //  Merged
      MarkInternalEdgeStatus(graph, scaffold, 1, FALSE, PAIRWISECHI2THRESHOLD_CGW, SLOPPY_EDGE_VARIANCE_THRESHHOLD);  //  Raw
    }
#endif

    //  Don't check chiSquared - leads to a massive misjoin in pging test set
#if 0
    if (IsScaffoldInternallyConnected(ScaffoldGraph, scaffold, false, true) != 1) {
      MarkInternalEdgeStatus(graph, scaffold, 2, TRUE,  PAIRWISECHI2THRESHOLD_CGW, SLOPPY_EDGE_VARIANCE_THRESHHOLD);  //  Merged
      MarkInternalEdgeStatus(graph, scaffold, 2, FALSE, PAIRWISECHI2THRESHOLD_CGW, SLOPPY_EDGE_VARIANCE_THRESHHOLD);  //  Raw
    }
#endif


    //  If the scaffold isn't 2-edge connected, break it.  We don't bother estimating gaps for any
    //  of the new scaffolds, it is up to the client.  In other words, the client is assumed to be
    //  iterating over all scaffolds.
    //
    if ((LSFlags & LeastSquares_Split) &&
        (CheckScaffoldConnectivityAndSplit(ScaffoldGraph, sID, ALL_TRUSTED_EDGES, FALSE) > 1))
      return(false);

    //dumpScaffoldContigPositions(graph, scaffold, "LSGE(pre)--");

    status = RecomputeOffsetsInScaffold(graph, sID, TRUE, TRUE, FALSE);

    //dumpScaffoldContigPositions(graph, scaffold, "LSGE(post)--");

    //  Well, it's connected at least.  I guess lapack failed to converge.
    if (status == RECOMPUTE_SINGULAR) {
      fprintf(stderr, "RecomputeOffsetsInScaffold() returned RECOMPUTE_SINGULAR on scaffold "F_CID" -- we're stuck\n", sID);
      break;
    }
    assert(status != RECOMPUTE_SINGULAR);

    //  We merged contigs, so remark edge status and repeat.
    if ((status == RECOMPUTE_CONTIGGED_CONTAINMENTS) ||
        (status == RECOMPUTE_FAILED_CONTIG_DELETED)) {
      fprintf(stderr, "RecomputeOffsetsInScaffold() returned RECOMPUTE_CONTIGGED_CONTAINMENTS on scaffold "F_CID"\n", sID);
      continue;
    }

    //  We changed scaffold order, so remark edge status and repeat.
    if (status == RECOMPUTE_FAILED_REORDER_NEEDED) {
      fprintf(stderr,"RecomputeOffsetsInScaffold() returned RECOMPUTE_FAILED_REORDER_NEEDED on scaffold "F_CID"\n", sID);
      continue;
    }

    if (status != RECOMPUTE_OK)
      fprintf(stderr, "RecomputeOffsetsInScaffold() on scaffold "F_CID" failed with code %d\n", sID, status);

    //  No more attempts will help.
    break;
  }

  //  If Least Squares fails, grab the current gap sizes from the scaffold, and rebuild it.  This
  //  might seem like it will do anything, but it resets the size of each contig to truth.  This
  //  seems to wander when we fail LeastSquares enough.

  if (status != RECOMPUTE_OK) {
    vector<double>  GAPmean;
    vector<double>  GAPvari;

    LengthT  lastEnd = { 0.0, 0.0 };

    CIScaffoldTIterator CIs;
    ChunkInstanceT     *CI;

    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
    while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
      LengthT minOffset = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;
      LengthT maxOffset = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetBEnd : CI->offsetAEnd;

      if (lastEnd.mean > 0) {
        GAPmean.push_back(minOffset.mean     - lastEnd.mean);
        GAPvari.push_back(minOffset.variance - lastEnd.variance);
      }

      lastEnd = maxOffset;
    }

    RebuildScaffoldGaps(graph, scaffold, &GAPmean[0], &GAPvari[0], false);
  }


  //  If Least Squares fails, do a greedy size estimate.  This DOES NOT WORK, especially in heavily interleaved scaffolds.  It tries
  //  to compute gap by gap, using mates that span a specific gap.  If those mates span other gaps, and those gap sizes are incorrect,
  //  then the estimate for this gap size is incorrect.
#if 0
  if (status != RECOMPUTE_OK) {
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

    fprintf(stderr, "LeastSquares failed; greedy gap size estimates used.\n");

    instrumentSCF   I(scaffold);

    I.analyze(libs);
    fprintf(stderr, "LeastSquares()--   scaffold (new) instrumenter happy %.1f gap %.1f misorient close %.1f correct %.1f far %.1f oriented close %.1f far %.1f missing %.1f external %.1f\n",
            I.numHappy, I.numGap, I.numMisClose, I.numMis, I.numMisFar, I.numTooClose, I.numTooFar, I.numMissing, I.numExternal);

    I.estimateGaps(libs, false);

    I.analyze(libs);
    fprintf(stderr, "LeastSquares()--   scaffold (new) instrumenter happy %.1f gap %.1f misorient close %.1f correct %.1f far %.1f oriented close %.1f far %.1f missing %.1f external %.1f\n",
            I.numHappy, I.numGap, I.numMisClose, I.numMis, I.numMisFar, I.numTooClose, I.numTooFar, I.numMissing, I.numExternal);

    RebuildScaffoldGaps(graph, scaffold, &I.GAPmean[0], &I.GAPvari[0], true);
    assert(0);
  }
#endif

  bool  redo = false;

  if ((LSFlags & LeastSquares_NoBounce) == 0)
    redo = BounceOutOfOrderContigs(graph, scaffold);

  if (bounceIteration > 5)
    redo = false;

  AdjustNegativePositions(graph, scaffold);
  AdjustNegativeVariances(graph, scaffold);
  AdjustNegativeGapVariances(graph, scaffold);

  if ((redo == false) && ((LSFlags & LeastSquares_Cleanup) == 0))
    //  No cleanup requested, and it didn't bounce, we're done.
    return(true);

  //   Otherwise, cleanup

  if (CleanupAScaffold(graph, scaffold, FALSE, NULLINDEX, FALSE) > 0)
    redo |= true;

  //  If we changed something, recompute least squares.  Note that even if we've reached the bounce
  //  iteration limit, if contigs were merged, we recompute.
  if (redo)
    return(LeastSquaresGapEstimates(graph, scaffold, LSFlags, bounceIteration+1));

  return(true);
}

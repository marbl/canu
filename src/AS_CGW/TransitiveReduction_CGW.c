
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
static char *rcsid = "$Id: TransitiveReduction_CGW.c,v 1.29 2010-01-07 01:09:25 brianwalenz Exp $";

//#define INSTRUMENT_CGW
//#define INSTRUMENT_SMOOTHED
#define INSTRUMENT_TRANS_REDUCED

#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "AS_global.h"
#include "AS_UTL_heap.h"
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "UnionFind_AS.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "ChiSquareTest_CGW.h"
#include "DataTypesREZ.h"
#include "CommonREZ.h"
#include "Instrument_CGW.h"
#include "UtilsREZ.h"

static VA_TYPE(CIEdgeT)        *graphCIEdges;

static
int
CompareCIEdgeTMeans (const void *c1, const void *c2) {
  CIEdgeT *s1 = GetCIEdgeT(graphCIEdges, *(CDS_CID_t *)c1);
  CIEdgeT *s2 = GetCIEdgeT(graphCIEdges, *(CDS_CID_t *)c2);

  if (s1->distance.mean > s2->distance.mean)
    return((int)1);
  else
    return((int)-1);

  return 0;
}

static
ChunkOrientationType
InferredEdgeOrientation(ChunkOrientationType leftOrient,
                        ChunkOrientationType rightOrient) {

  if ((leftOrient == AB_AB) || (leftOrient == BA_AB))
    return(((rightOrient == AB_AB) || (rightOrient == BA_AB)) ? AB_AB : AB_BA);
  else
    return(((rightOrient == AB_AB) || (rightOrient == BA_AB)) ? BA_AB : BA_BA);
}

static
ChunkOrientationType
TransitiveEdgeOrientation(ChunkOrientationType leftOrient,
                          ChunkOrientationType rightOrient) {

  if ((leftOrient == AB_AB) || (leftOrient == AB_BA))
    return(((rightOrient == AB_AB) || (rightOrient == BA_AB)) ? AB_AB : AB_BA);
  else
    return(((rightOrient == AB_AB) || (rightOrient == BA_AB)) ? BA_AB : BA_BA);
}



#define AS_CGW_MAX_FTEP_RECURSE_DEPTH 15

static
int
FoundTransitiveEdgePath(ScaffoldGraphT *graph,
                        CIEdgeT *pathEdge,
                        CIEdgeT *edge,  /* edge we're trying to transitively remove */
                        ChunkInstanceT *startCI, /* Start of path */
                        ChunkInstanceT *endCI,   /* End of Path */
                        ChunkInstanceT *thisCI, /* Where we are now */
                        int recurseDepth) {

  CIEdgeT localPathEdge;
  CIEdgeT *transEdge;
  GraphEdgeIterator transEdges;
  float chiSquaredValue;
  ChunkOrientationType edgeOrient, pathEdgeOrient;
  int returnVal = FALSE;

  assert(!isSloppyEdge(edge));

  if (++recurseDepth > AS_CGW_MAX_FTEP_RECURSE_DEPTH)
    return(FALSE);

  edgeOrient = GetEdgeOrientationWRT(edge, startCI->id);
  pathEdgeOrient = GetEdgeOrientationWRT(pathEdge, startCI->id);

  if (thisCI == endCI)
    // we got there
    if (edgeOrient == pathEdgeOrient)
      // correct orientation
      return(PairwiseChiSquare((float)pathEdge->distance.mean,
                               (float)pathEdge->distance.variance,
                               (float)edge->distance.mean,
                               (float)edge->distance.variance, (LengthT *)NULL,
                               &chiSquaredValue, (float)PAIRWISECHI2THRESHOLD_CGW));
    else
      // incorrect orientation
      return(FALSE);

  if (pathEdge->distance.mean > edge->distance.mean)
    if (!PairwiseChiSquare((float)pathEdge->distance.mean,
                           (float)pathEdge->distance.variance,
                           (float)edge->distance.mean,
                           (float)edge->distance.variance, (LengthT *)NULL,
                           &chiSquaredValue, (float)PAIRWISECHI2THRESHOLD_CGW))
      return(FALSE);

  if ((pathEdgeOrient == AB_AB) || (pathEdgeOrient == BA_AB))
    //  transEdge must be a B_END edge
    InitGraphEdgeIterator(graph->ContigGraph, thisCI->id, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &transEdges);// Use merged edges
  else
    //  transEdge must be a A_END edge
    InitGraphEdgeIterator(graph->ContigGraph, thisCI->id, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &transEdges);// Use merged edges

  while((transEdge = NextGraphEdgeIterator(&transEdges))!= NULL) {
    ChunkInstanceT *nextCI;

    if (!transEdge->flags.bits.isActive || !transEdge->flags.bits.isUniquetoUnique)
      //  Do not use inactive edges or edges to nonUnique CIs.
      continue;

    if (isSingletonOverlapEdge(transEdge))
      //  Do not use overlap edges in path walks
      continue;

    assert(!isSloppyEdge(transEdge));

    nextCI = GetGraphNode(graph->ContigGraph,
			  ((transEdge->idA == thisCI->id) ?
			   transEdge->idB : transEdge->idA));
    localPathEdge = *pathEdge;
    localPathEdge.distance.mean += thisCI->bpLength.mean + transEdge->distance.mean;
    localPathEdge.distance.variance += thisCI->bpLength.variance + transEdge->distance.variance;
    localPathEdge.idB = (transEdge->idA == thisCI->id) ? transEdge->idB : transEdge->idA;
    localPathEdge.orient = TransitiveEdgeOrientation(pathEdgeOrient, GetEdgeOrientationWRT(transEdge, thisCI->id));

    if (FoundTransitiveEdgePath(graph, &localPathEdge, edge, startCI, endCI, nextCI, recurseDepth)) {
      transEdge->flags.bits.isConfirmed = TRUE;
      returnVal = TRUE;

      // For speed we are going to return when we find the first path and not find all paths so as
      // to possibly confirm some extra edges.
      //
      return(returnVal);
    }
  }
  return(returnVal);
}

static
void
MarkPathRemovedEdgesOneEnd(ScaffoldGraphT *graph, ChunkInstanceT *thisCI,
                           int end) {
  GraphEdgeIterator edges;
  CIEdgeT *edge;
  InitGraphEdgeIterator(graph->ContigGraph, thisCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
  while((edge = NextGraphEdgeIterator(&edges))!= NULL) {
    GraphEdgeIterator transEdges;
    CIEdgeT pathEdge;
    CIEdgeT *transEdge;
    ChunkInstanceT *endCI;

    if (edge->idA != thisCI->id)
      // Only check an edge once in canonical direction.
      continue;

    if (!edge->flags.bits.isActive || !edge->flags.bits.isUniquetoUnique)
      // No need to check inactive edges - only interested in edges between Unique CIs.
      continue;

    assert(!isSloppyEdge(edge));
    endCI = GetGraphNode(graph->ContigGraph,
                         ((edge->idA == thisCI->id) ? edge->idB : edge->idA));
    // Only use edges on the same end of thisCI to transitively remove an edge.
    InitGraphEdgeIterator(graph->ContigGraph, thisCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &transEdges);// Use merged edges
    while((transEdge = NextGraphEdgeIterator(&transEdges))!= NULL) {
      ChunkInstanceT *nextCI;
      CDS_CID_t nextCID;
      int recurseDepth;

      if (edge == transEdge)
	// Do not remove an edge with itself.
	continue;

      if (!transEdge->flags.bits.isActive || !transEdge->flags.bits.isUniquetoUnique)
	// Do not use inactive edges - only interested in edges between Unique CIs.
	continue;

      assert(!isSloppyEdge(transEdge));

      nextCID = (transEdge->idA == thisCI->id) ? transEdge->idB : transEdge->idA;
      nextCI = GetGraphNode(graph->ContigGraph, nextCID);
      pathEdge = *transEdge;
      pathEdge.idA = thisCI->id;
      pathEdge.idB = nextCID;
      pathEdge.orient = GetEdgeOrientationWRT(transEdge, thisCI->id);
      recurseDepth = 0;
      if (FoundTransitiveEdgePath(graph, &pathEdge, edge, thisCI, endCI,
                                  nextCI, recurseDepth)) {
	edge->flags.bits.isTransitivelyRemoved = TRUE;
	edge->flags.bits.isConfirmed = TRUE;
	transEdge->flags.bits.isConfirmed = TRUE;
      }
    }
  }
}

static
void
MarkPathRemovedEdges(ScaffoldGraphT *graph) {
  GraphNodeIterator nodes;
  ChunkInstanceT   *thisCI;

  InitGraphNodeIterator(&nodes, graph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL) {
    MarkPathRemovedEdgesOneEnd(graph, thisCI, A_END);
    MarkPathRemovedEdgesOneEnd(graph, thisCI, B_END);
  }
}

static
void
MarkTwoHopConfirmedEdgesOneEnd(ScaffoldGraphT *graph,
                               ChunkInstanceT *thisCI,
                               int end) {
  GraphEdgeIterator edges;
  CIEdgeT          *edge;
  InitGraphEdgeIterator(graph->ContigGraph, thisCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges

  while((edge = NextGraphEdgeIterator(&edges))!= NULL) {
    GraphEdgeIterator hop2Edges;
    CIEdgeT pathEdge;
    CIEdgeT *hop2Edge;
    CDS_CID_t middleCID;
    ChunkInstanceT *middleCI;

    if (!edge->flags.bits.isActive || !edge->flags.bits.isUniquetoUnique)
      // No need to check inactive edges - only interested in edges between Unique CIs.
      continue;

    if (isSingletonOverlapEdge(edge))
      //  Do not use overlap edges in two hops
      continue;

    middleCID = (edge->idA == thisCI->id) ? edge->idB : edge->idA;
    middleCI = GetGraphNode(graph->ContigGraph, middleCID);

    pathEdge = *edge;
    pathEdge.idA = thisCI->id;
    pathEdge.idB = middleCID;
    pathEdge.orient = GetEdgeOrientationWRT(edge, thisCI->id);
    pathEdge.distance.mean += middleCI->bpLength.mean;
    pathEdge.distance.variance += middleCI->bpLength.variance;

    if ((pathEdge.orient == AB_AB) ||
        (pathEdge.orient == BA_AB))
      //  hop2Edge must be a B_END edge
      InitGraphEdgeIterator(graph->ContigGraph, middleCID, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &hop2Edges);// Use merged edges
    else
      //  hop2Edge must be a A_END edge
      InitGraphEdgeIterator(graph->ContigGraph, middleCID, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &hop2Edges);// Use merged edges

    while((hop2Edge = NextGraphEdgeIterator(&hop2Edges))!= NULL) {
      ChunkInstanceT *endCI;
      ChunkOrientationType combinedOrient;
      double combinedMean, combinedVariance;
      GraphEdgeIterator confirmEdges;
      CIEdgeT *confirmEdge;

      if (!hop2Edge->flags.bits.isActive || !hop2Edge->flags.bits.isUniquetoUnique)
	// Do not use inactive edges - only interested in edges between Unique CIs.
	continue;

      if (isSingletonOverlapEdge(hop2Edge))
        //  Do not use overlap edges in two hops
	continue;

      endCI = GetGraphNode(graph->ContigGraph, ((hop2Edge->idA == middleCID) ? hop2Edge->idB : hop2Edge->idA));
      combinedOrient = TransitiveEdgeOrientation(pathEdge.orient, GetEdgeOrientationWRT(hop2Edge, middleCID));
      combinedMean = pathEdge.distance.mean + hop2Edge->distance.mean;
      combinedVariance = pathEdge.distance.variance + hop2Edge->distance.variance;

      // Only use edges on the same end of thisCI to confirm a two hop edge.
      // Use a copy of the original iterator so that we do not find the
      // same pair of two hop paths twice.

      confirmEdges = edges;
      while((confirmEdge = NextGraphEdgeIterator(&confirmEdges))!= NULL) {
	GraphEdgeIterator confirmHop2Edges;
	CIEdgeT confirmPathEdge;
	CIEdgeT *confirmHop2Edge;
	CDS_CID_t confirmMiddleCID;
	ChunkInstanceT *confirmMiddleCI;

	if (!confirmEdge->flags.bits.isActive || !confirmEdge->flags.bits.isUniquetoUnique)
	  // Do not use inactive edges - only interested in edges between Unique CIs.
	  continue;

	if (isSingletonOverlapEdge(confirmEdge))
          //  Do not use overlap edges in two hops
	  continue;

	confirmMiddleCID = (confirmEdge->idA == thisCI->id) ? confirmEdge->idB : confirmEdge->idA;
	confirmMiddleCI = GetGraphNode(graph->ContigGraph, confirmMiddleCID);

	confirmPathEdge = *confirmEdge;
	confirmPathEdge.idA = thisCI->id;
	confirmPathEdge.idB = confirmMiddleCID;
	confirmPathEdge.orient = GetEdgeOrientationWRT(confirmEdge, thisCI->id);
	confirmPathEdge.distance.mean += confirmMiddleCI->bpLength.mean;
	confirmPathEdge.distance.variance += confirmMiddleCI->bpLength.variance;

	if ((confirmPathEdge.orient == AB_AB) ||
            (confirmPathEdge.orient == BA_AB))
	  //  hop2Edge must be a B_END edge
	  InitGraphEdgeIterator(graph->ContigGraph, confirmMiddleCID, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &confirmHop2Edges);// Use merged edges
	else
	  //  hop2Edge must be a A_END edge
	  InitGraphEdgeIterator(graph->ContigGraph, confirmMiddleCID, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &confirmHop2Edges);// Use merged edges

	while((confirmHop2Edge = NextGraphEdgeIterator(&confirmHop2Edges))!= NULL) {
	  ChunkInstanceT *confirmEndCI;
	  ChunkOrientationType confirmOrient;
	  double confirmMean, confirmVariance;
	  float chiSquaredValue;

	  if (!confirmHop2Edge->flags.bits.isActive || !confirmHop2Edge->flags.bits.isUniquetoUnique)
	    // Do not use inactive edges - only interested in edges between Unique CIs.
	    continue;

	  if (isSingletonOverlapEdge(confirmHop2Edge))
            //  Do not use overlap edges in two hops
	    continue;

	  confirmEndCI = GetGraphNode(graph->ContigGraph,
                                      ((confirmHop2Edge->idA == confirmMiddleCID) ? confirmHop2Edge->idB : confirmHop2Edge->idA));
	  confirmOrient = TransitiveEdgeOrientation(confirmPathEdge.orient,
                                                    GetEdgeOrientationWRT(confirmHop2Edge, confirmMiddleCID));

	  confirmMean = confirmPathEdge.distance.mean + confirmHop2Edge->distance.mean;
	  confirmVariance = confirmPathEdge.distance.variance + confirmHop2Edge->distance.variance;

	  if ((endCI == confirmEndCI) && (combinedOrient == confirmOrient) &&
              PairwiseChiSquare((float)combinedMean, (float)combinedVariance,
                                (float)confirmMean, (float)confirmVariance,
                                (LengthT *)NULL, &chiSquaredValue,
                                (float)PAIRWISECHI2THRESHOLD_CGW)) {

	    edge->flags.bits.isConfirmed = TRUE;
	    hop2Edge->flags.bits.isConfirmed = TRUE;
	    confirmEdge->flags.bits.isConfirmed = TRUE;
	    confirmHop2Edge->flags.bits.isConfirmed = TRUE;

            //PrintGraphEdge(stderr, graph->ContigGraph, "2Hop\t1 ", edge, thisCI->id);
            //PrintGraphEdge(stderr, graph->ContigGraph, "\t2 ", hop2Edge, middleCID);
            //PrintGraphEdge(stderr, graph->ContigGraph, "\tC ", confirmEdge, thisCI->id);
            //PrintGraphEdge(stderr, graph->ContigGraph, "\t2 ", confirmHop2Edge, confirmMiddleCID);
	  }
	}
      }
    }
  }
}

static
void
MarkTwoHopConfirmedEdges(ScaffoldGraphT *graph) {
  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;

  InitGraphNodeIterator(&nodes, graph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL) {
    MarkTwoHopConfirmedEdgesOneEnd(graph, thisCI, A_END);
    MarkTwoHopConfirmedEdgesOneEnd(graph, thisCI, B_END);
  }
}

static
int
CompareBestCIEdgeT(const void *c1, const void *c2) {
  CIEdgeT *s1 = (CIEdgeT *)c1;
  CIEdgeT *s2 = (CIEdgeT *)c2;

  if (s1->flags.bits.isEssential != s2->flags.bits.isEssential)
    return s1->flags.bits.isEssential ? 1 : -1;

  if (s1->flags.bits.isConfirmed != s2->flags.bits.isConfirmed)
    return s1->flags.bits.isConfirmed ? 1 : -1;

  if (s1->flags.bits.isLeastSquares != s2->flags.bits.isLeastSquares)
    return s1->flags.bits.isLeastSquares ? 1 : -1;

  if (s1->flags.bits.isActive != s2->flags.bits.isActive)
    return s1->flags.bits.isActive ? 1 : -1;

  if (s1->flags.bits.isProbablyBogus != s2->flags.bits.isProbablyBogus)
    return s1->flags.bits.isProbablyBogus ? -1 : 1;

  if (s1->edgesContributing < s2->edgesContributing)
    return -1;

  if (s1->edgesContributing > s2->edgesContributing)
    return 1;

  if (isOverlapEdge(s1) && !isOverlapEdge(s2))
    return -1;

  if (isOverlapEdge(s2) && !isOverlapEdge(s1))
    return 1;

  if (s1->distance.variance > s2->distance.variance)
    return -1;

  if (s1->distance.variance < s2->distance.variance)
    return 1;

  if (s1->distance.mean > s2->distance.mean)
    return -1;

  if (s1->distance.mean < s2->distance.mean)
    return 1;

  if (s1->minDistance > s2->minDistance)
    return -1;

  if (s1->minDistance < s2->minDistance)
    return 1;

  return 0;
}

static
void
MarkRedundantUniqueToUniqueEdges(ScaffoldGraphT *graph) {

  CDS_CID_t cid;
  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;

  InitGraphNodeIterator(&nodes, graph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL) {
    GraphEdgeIterator edges;
    CIEdgeT *edge;

    cid = thisCI->id;

    InitGraphEdgeIterator(graph->ContigGraph, thisCI->id , ALL_END,	ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
    edge = NextGraphEdgeIterator(&edges);
    while(edge != NULL) {
      CIEdgeT *bestEdge;

      if ((edge->idA != thisCI->id) || !edge->flags.bits.isActive || !edge->flags.bits.isUniquetoUnique) {
	// Only check canonical edges so as not to repeat the work.  No need to check inactive edges
	// - only interested in edges between Unique CIs.
	edge = NextGraphEdgeIterator(&edges);
	continue;
      }

      bestEdge = edge;

      // We make the assumption that the edges are sorted such that the edges
      // between a pair of CIs are contiguous.
      while(((edge = NextGraphEdgeIterator(&edges)) != NULL) &&
	    (edge->idB == bestEdge->idB)) {
	if (!edge->flags.bits.isActive)
	  // No need to check inactive edges.
	  continue;
	if (CompareBestCIEdgeT(bestEdge, edge) < 0) {
	  bestEdge->flags.bits.isActive = FALSE;
	  bestEdge->flags.bits.isRedundantRemoved = TRUE;
	  bestEdge = edge;
	}else{
	  edge->flags.bits.isActive = FALSE;
	  edge->flags.bits.isRedundantRemoved = TRUE;
	}
      }
    }
  }
}

static
EdgeCGW_T *
FindEdgeBetweenCIsChiSquare(GraphCGW_T *graph,
                            ChunkInstanceT *sourceCI, CDS_CID_t targetId,
                            ChunkOrientationType edgeOrient,
                            double inferredMean, double inferredVariance,
                            float *returnChiSquaredValue,
                            float chiSquareThreshold, int *isEssential) {
  GraphEdgeIterator edges;
  CIEdgeT *edge;
  CIEdgeT *bestEdge = (CIEdgeT *)NULL;
  int overlapEdgeExists = FALSE;
  int end;
  float chiSquaredValue, bestChiSquaredValue = FLT_MAX;

  *isEssential = FALSE;
  /* Search all orientations in case an essential edge exists between
     the CIs which is not in the correct orientation. */
  InitGraphEdgeIterator(graph, sourceCI->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
  while((edge = NextGraphEdgeIterator(&edges))!= NULL) {
    CDS_CID_t otherCID = (edge->idA == sourceCI->id) ? edge->idB : edge->idA;

    if (isSloppyEdge(edge))
      // We aren't interested in high variance (sloppy) edges, since they will pass the chi square too easily
      continue;

    if (otherCID == targetId) {
      // If there is an overlap edge with the 'right' orientation, remember it
      if (isOverlapEdge(edge) && (GetEdgeOrientationWRT(edge, sourceCI->id) == edgeOrient))
	overlapEdgeExists = TRUE;

      // Is this the best edge irrespective of orientation?
      if ((bestEdge == (CIEdgeT *)NULL) || (CompareBestCIEdgeT((const void *)edge, (const void *)bestEdge) > 0))
	bestEdge = edge;
    }
  }
  // Were any edges found irrespective of orientation between source,target?
  if (bestEdge != (CIEdgeT *)NULL) {
    // This is a return value that, if true, causes the caller NOT to add an inferred
    // edge when this function returns NULL
    *isEssential = getEssentialEdgeStatus(bestEdge);

    // If the orientation was 'wrong'
    if (GetEdgeOrientationWRT(bestEdge, sourceCI->id) != edgeOrient) {
      // and the edge is essential, return NULL
      if (getEssentialEdgeStatus(bestEdge) == TRUE)
	return((CIEdgeT *)NULL);

    }else // If the orientation is 'right', check for Chi Square compatibility
#if 0
      if ( inferredVariance > FLT_MAX ) inferredVariance = FLT_MAX;  // hack instituted on 5/15/01 to help mouse_20010508
#endif
    if (PairwiseChiSquare((float)inferredMean, (float)inferredVariance,
                          (float)bestEdge->distance.mean,
                          (float)bestEdge->distance.variance,
                          (LengthT *)NULL, &chiSquaredValue,
                          (float)chiSquareThreshold)) {
      *returnChiSquaredValue = chiSquaredValue;
      return(bestEdge);
    }else{ // If right orientation, but not chi square compatible
      // If marked essential, return NULL
      if (getEssentialEdgeStatus(bestEdge) == TRUE)
	return((CIEdgeT *)NULL);
    }
  }


  // If there isn't an overlap edge in the right orientation, look for one
  if (!overlapEdgeExists) {
    int32 minOverlap, maxOverlap;
    int32 minVariance;
    int32 delta;

    if (inferredVariance < 1.0)
      minVariance = 1;
    else
      minVariance = (int32) sqrt(inferredVariance);

    delta = 3 * minVariance;
    minOverlap = MAX(CGW_MISSED_OVERLAP,
		     -(inferredMean + delta));
    maxOverlap = -(inferredMean - delta);
    if (maxOverlap >= CGW_MISSED_OVERLAP) {
      ChunkOverlapCheckT olap;

      // The following has side effect of adding an overlap edge to the graph if found
      // which can cause a realloc of the edge array and impact the calling
      // routine. In addition, currently an edge returned by this function
      // will point to the current edge array but this should be assured
      // when any future code changes are made.
      assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
      olap = OverlapChunks(graph,
                           sourceCI->id, targetId, // handles suspicious
			   edgeOrient, minOverlap, maxOverlap,
                           AS_CGW_ERROR_RATE, TRUE);
      // We found an overlap between source and target in the right orientation
      if (olap.suspicious) {
        fprintf(stderr,"* TR: SUSPICIOUS Overlap found! Looked for ("
                F_CID ","F_CID ",%c)["F_S32","F_S32"] found ("
                F_CID ","F_CID ",%c) "F_S32" (would have been "
                F_S32")\n",
                sourceCI->id, targetId, edgeOrient,
                minOverlap, maxOverlap,
                olap.spec.cidA, olap.spec.cidB, olap.spec.orientation,
                olap.overlap,olap.overlap+olap.ahg+olap.bhg);
      }else if (olap.overlap) {
      }else if (bestEdge == (CIEdgeT *)NULL) { // If we didn't find an overlap, and there was no other edge return NULL
	return((CIEdgeT *)NULL);
      }
    }else if (bestEdge == (CIEdgeT *)NULL) { // There was no potential overlap, and there was no other edge return NULL
      return((CIEdgeT *)NULL);
    }
  }

  // If the 'best' edge didn't pass the above tests (orientation and chi square), we now look for an edge with proper orientation
  // that passes chi square
  // It is important to set bestEdge to NULL here so that any edge found and
  // returned by this function points to the current edge array which might
  // have been reallocated when an overlap edge was added in the code above.

  bestEdge = (CIEdgeT *)NULL;
  if ((edgeOrient == AB_AB) || (edgeOrient == AB_BA)) {
    /* edgeOrient == AB_XX */
    end = B_END;
  }else{
    /* edgeOrient == BA_XX */
    end = A_END;
  }
  InitGraphEdgeIterator(graph, sourceCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT,  &edges);// Use merged edges
  while((edge = NextGraphEdgeIterator(&edges))!= NULL) {
    CDS_CID_t otherCID = (edge->idA == sourceCI->id) ? edge->idB : edge->idA;

    if (isSloppyEdge(edge))
      // We aren't interested in high variance (sloppy) edges, since they will pass the chi square too easily
      continue;

    if ((otherCID == targetId) &&
        (GetEdgeOrientationWRT(edge, sourceCI->id) == edgeOrient) &&
        PairwiseChiSquare((float)inferredMean,
                          (float)inferredVariance,
#if 0
                          (float) ( inferredVariance > FLT_MAX ? FLT_MAX : inferredVariance),  // hack instituted on 5/15/01
                          // to help mouse_20010508
#endif
                          (float)edge->distance.mean,
                          (float)edge->distance.variance,
                          (LengthT *)NULL, &chiSquaredValue,
                          (float)chiSquareThreshold)) {
      if (bestEdge == (CIEdgeT *)NULL ||
          (chiSquaredValue < bestChiSquaredValue)) {
	bestEdge = edge;
	bestChiSquaredValue = chiSquaredValue;
	*isEssential = getEssentialEdgeStatus(bestEdge);
      }
    }
  }
  *returnChiSquaredValue = bestChiSquaredValue;
  return(bestEdge);
}


static
int
RecursiveSmoothWithInferredEdges(ScaffoldGraphT *graph,
                                 ChunkInstanceT *thisCI,
                                 int end,
                                 int * numInferredAdded,
                                 CDS_CID_t firstID,
                                 int instFirstEnd,
                                 CDS_CID_t *lastID,
                                 ScaffoldInstrumenter * smoothed_si,
                                 FILE * instSmoothSuccessfp,
                                 FILE * instSmoothFailurefp) {
  int numBranch;
  CDS_CID_t *branchEdges;
  CDS_CID_t *existingEdges;
  CDS_CID_t *addedInferred;
  CDS_CID_t *branchSource;
  CDS_CID_t *branchTarget;
  CDS_CID_t *branchEnd;
  CDS_CID_t *addedEnd;
  CDS_CID_t *addedPtr;
  CDS_CID_t *existingEnd;
  CDS_CID_t *existingPtr;
  CIEdgeT *sourceEdge = NULL;
  ChunkOrientationType sourceEdgeOrient;
  ChunkInstanceT *sourceCI;
  double sourceMean, sourceVariance;
  int failedInfer;
  int numEssentialAdded;
  int numEssentialRemoved;
  GraphEdgeIterator edges;
  CIEdgeT *edge;
  Target_Info_t *branchTargets, *branchTargetsEnd;
  double maxBound;
#if 0
  int firstTarget, numTargetsBestPath, firstTargetBestPath;
  LengthT  targetPosition;
#endif

  // get the number of essential edges off end of thisCI
  if (end == A_END)
    numBranch = thisCI->numEssentialA;
  else
    numBranch = thisCI->numEssentialB;

  /*
    if there is one essential edge off end of thisCI, simple resolution:
    1. get the next CI (sourceCI) and the edge to it (sourceEdge)
    2. set sourceCI flags so it can be smoothed later
    3. set the essentialX of thisCI & sourceCI to the index of sourceEdge
    4. for instrumenting, if lastID isn't set yet, set it to sourceCI->id
  */
  if (numBranch < 2) {
    CDS_CID_t sourceEdgeIndex;
    sourceCI = (ChunkInstanceT *)NULL;
    assert(numBranch == 1);
    InitGraphEdgeIterator(graph->ContigGraph, thisCI->id, end,ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
    while((edge = NextGraphEdgeIterator(&edges))!= NULL) {
      if (!getEssentialEdgeStatus(edge))
	// Skip nonessential edges.
	continue;

      assert(sourceCI == (ChunkInstanceT *)NULL);
      sourceCI = GetGraphNode(graph->ContigGraph, ((edge->idA == thisCI->id) ? edge->idB : edge->idA));
      sourceEdge = edge;
    }
    assert(sourceCI != (ChunkInstanceT *)NULL);

    // instrument the scaffold to see if it's acceptable
    // before modifying anything
    *lastID = (*lastID == NULLINDEX) ? sourceCI->id : *lastID;

#ifdef INSTRUMENT_TRANS_REDUCED
    {
      double badMates;
      double allMates;

      InstrumentContigPath(graph, smoothed_si,
                           firstID, instFirstEnd, *lastID);
      badMates = GetMateStatsBad(&(smoothed_si->mates.inter));
      allMates = badMates + GetMateStatsHappy(&(smoothed_si->mates.inter));

      if ((int) allMates > 0 && badMates / allMates > .05)
        {
          // bad 'scaffold' - return failure
#ifdef INSTRUMENT_SMOOTHED
          PrintScaffoldInstrumenter(graph, smoothed_si,
                                    InstrumenterVerbose2,
                                    "", instSmoothFailurefp);
#endif
          return(FALSE);
        }
      else
        {
          // since we've done the instrumenting, print successes here
#ifdef INSTRUMENT_SMOOTHED
          PrintScaffoldInstrumenter(graph, smoothed_si,
                                    InstrumenterVerbose2,
                                    "", instSmoothSuccessfp);
#endif
        }
    }
#endif

    sourceCI->flags.bits.smoothSeenAlready = FALSE;
    sourceCI->smoothExpectedCID = NULLINDEX;
    sourceEdgeIndex = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->ContigGraph->edges, sourceEdge);
    if (end == B_END) {
      thisCI->essentialEdgeB = sourceEdgeIndex;
    }else{
      thisCI->essentialEdgeA = sourceEdgeIndex;
    }
    sourceEdgeOrient = GetEdgeOrientationWRT(sourceEdge, thisCI->id);
    if ((sourceEdgeOrient == AB_AB) || (sourceEdgeOrient == BA_AB)) {
      sourceCI->essentialEdgeA = sourceEdgeIndex;
    }else{
      sourceCI->essentialEdgeB = sourceEdgeIndex;
    }

    return(TRUE);
  }

  {
    branchEdges = (CDS_CID_t *)safe_malloc(numBranch * sizeof(*branchEdges));
    existingEdges = (CDS_CID_t *)safe_malloc((numBranch - 1) * sizeof(*existingEdges));
    addedInferred = (CDS_CID_t *)safe_malloc((numBranch - 1) * sizeof(*addedInferred));
    branchTargets = (Target_Info_t *)safe_malloc(numBranch * sizeof(*branchTargets));

    InitGraphEdgeIterator(graph->ContigGraph, thisCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
    maxBound = 0.0;
    branchEnd = branchEdges;
    existingEnd = existingEdges;
    addedEnd = addedInferred;
    branchTargetsEnd = branchTargets;

    /*
      Loop over all essential edges off end of thisCI to attempt to smooth:
      1. identify the CI at the other end of each essential edge (targetCI)
      if it can't be smoothed, then fail
      2. otherwise store the edge indices, orientations, and compute acceptable length intervals
      3. sort these edge references by near to far - nearest = sourceCI, farther ones = targetCI
    */
    while((edge = NextGraphEdgeIterator(&edges))!= NULL) {
      ChunkInstanceT *targetCI;
      if (!getEssentialEdgeStatus(edge))
        // Skip nonessential edges.
        continue;

      targetCI = GetGraphNode(graph->ContigGraph,
                              ((edge->idA == thisCI->id) ? edge->idB : edge->idA));
      if (targetCI->flags.bits.smoothSeenAlready ||
          ((targetCI->smoothExpectedCID != NULLINDEX) &&
           (targetCI->smoothExpectedCID != thisCI->id))) {

#ifdef INSTRUMENT_SMOOTHED
        // instrument the failed 'scaffold' before returning failure
        InstrumentContigPath(graph, smoothed_si,
                             firstID, instFirstEnd, thisCI->id);
        PrintScaffoldInstrumenter(graph, smoothed_si,
                                  InstrumenterVerbose2,
                                  "", instSmoothFailurefp);
#endif

        safe_free(branchEdges);
        safe_free(existingEdges);
        safe_free(addedInferred);
        safe_free(branchTargets);

        return(FALSE);
      }
      // We need to save indices instead of pointers because the Variable
      // Array graph->CIEdges can get reallocated when we insert inferred
      // edges.
      *branchEnd++ = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->ContigGraph->edges, edge);
      branchTargetsEnd->id = targetCI->id;
      branchTargetsEnd->lo = (edge->distance.mean -
                              (3 * sqrt(edge->distance.variance))) +
        (targetCI->bpLength.mean - (3 * sqrt(targetCI->bpLength.variance)));
      branchTargetsEnd->hi = (edge->distance.mean +
                              (3 * sqrt(edge->distance.variance))) +
        (targetCI->bpLength.mean + (3 * sqrt(targetCI->bpLength.variance)));
      if (branchTargetsEnd->hi > maxBound) {
        maxBound = branchTargetsEnd->hi;
      }
      branchTargetsEnd->orient = GetEdgeOrientationWRT(edge, thisCI->id);
      branchTargetsEnd++;
    }
    if ((branchEnd - branchEdges) != numBranch) {
      fprintf(stderr, "RecSub numBranch %d,%d\n",
	      (int)(branchEnd - branchEdges), numBranch);
      DumpChunkInstance(stderr, graph, thisCI, FALSE, TRUE, TRUE, FALSE);
      fflush(NULL);
      assert(FALSE);
    }

    graphCIEdges = graph->ContigGraph->edges;
    qsort((void *)branchEdges, numBranch, sizeof(*branchEdges), CompareCIEdgeTMeans);

    /*
      After sorting essential edges from nearest to farthest,
      1. set sourceCI to nearest CI - intent is to keep existing edge between thisCI & sourceCI
      2. iterate over all other (longer) edges as targetEdge's to targetCI's
      attempt to establish inferred or real edge with acceptable chi2 between sourceCI & targetCI
      if that works, depracate the edge between thisCI and targetCI
    */
    branchSource = branchEdges;
    sourceEdge = GetGraphEdge(graph->ContigGraph, *branchSource);
    sourceCI = GetGraphNode(graph->ContigGraph,
                            ((sourceEdge->idA == thisCI->id) ? sourceEdge->idB :
                             sourceEdge->idA));
    sourceMean = sourceCI->bpLength.mean + sourceEdge->distance.mean;
    sourceVariance = sourceCI->bpLength.variance + sourceEdge->distance.variance;
    sourceEdgeOrient = GetEdgeOrientationWRT(sourceEdge, thisCI->id);
    sourceCI->flags.bits.smoothSeenAlready = TRUE;
    for(branchTarget = branchEdges + 1, failedInfer = FALSE,
          numEssentialAdded = 0, numEssentialRemoved = 0;
        branchTarget < branchEnd;
        branchTarget++) {
      CIEdgeT *targetEdge = GetGraphEdge(graph->ContigGraph, *branchTarget);
      ChunkInstanceT *targetCI = GetGraphNode(graph->ContigGraph,
                                              ((targetEdge->idA == thisCI->id) ? targetEdge->idB : targetEdge->idA));
      CIEdgeT *existingEdge;
      CIEdgeT *inferredEdge;
      CDS_CID_t inferredEdgeIndex;
      double inferredMean, inferredVariance;
      LengthT inferredDistance;
      float chiSquaredValue;
      ChunkOrientationType inferredEdgeOrient;
      ChunkOrientationType targetEdgeOrient;
      int isOverlapInferred = FALSE;
      int existingIsEssential;

      /* Attempt to create inferred edge between sourceCI & targetCI:
         1. compute distance stats of potential inferred edge
         2. see if there is an existing edge compatible with inferred distance
      */
      targetEdgeOrient = GetEdgeOrientationWRT(targetEdge, thisCI->id);
      targetCI->smoothExpectedCID = sourceCI->id;
      inferredDistance.mean =
        inferredMean = targetEdge->distance.mean - sourceMean;
      inferredDistance.variance =
        inferredVariance = targetEdge->distance.variance + sourceVariance;
      inferredEdgeOrient = InferredEdgeOrientation(sourceEdgeOrient,
                                                   targetEdgeOrient);

      existingEdge = FindEdgeBetweenCIsChiSquare(graph->ContigGraph, sourceCI,
                                                 (targetEdge->idA == thisCI->id) ? targetEdge->idB :
                                                 targetEdge->idA,
                                                 inferredEdgeOrient, inferredMean, inferredVariance,
                                                 &chiSquaredValue, (float)PAIRWISECHI2THRESHOLD_CGW,
                                                 &existingIsEssential);
      targetEdge = GetGraphEdge(graph->ContigGraph, *branchTarget); // FindEdgeBetweenCIsChiSquare may cause edges to realloc
      if (existingEdge != (CIEdgeT *)NULL) {
        existingEdge->flags.bits.wasEssential = getEssentialEdgeStatus(existingEdge);

        if (!existingEdge->flags.bits.isEssential) {
          setEssentialEdgeStatus(existingEdge, TRUE);
          //	existingEdge->flags.bits.isEssential = TRUE;
          //	assert(!isSloppyEdge(existingEdge));
          numEssentialAdded++;
        }
        targetEdge->flags.bits.isInferredRemoved = TRUE;
        setEssentialEdgeStatus(targetEdge, FALSE);
        //targetEdge->flags.bits.isEssential = FALSE;
        if (targetEdge->flags.bits.isInferred && targetEdge->flags.bits.isTentative) {
          DeleteGraphEdge(graph->ContigGraph, targetEdge);
        }
        numEssentialRemoved++;
        *existingEnd++ = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->ContigGraph->edges, existingEdge);
        continue;
      }else if (existingIsEssential) {
        failedInfer = TRUE;
        break;
      }
      if (inferredMean < - CGW_MISSED_OVERLAP) {
        // adjust mean and variance so things don't get sloppy (dewim 09/12/01)
        double new_stddev = (inferredDistance.mean + CGW_MISSED_OVERLAP + 3. * sqrt(inferredDistance.variance)) / 3.;
        new_stddev = (new_stddev < 1.) ? 1. : new_stddev;
        inferredDistance.variance = new_stddev * new_stddev;
        inferredDistance.mean = - CGW_MISSED_OVERLAP;

#if 0
        if ( inferredVariance > FLT_MAX ) inferredVariance = FLT_MAX;  // hack instituted on 5/15/01 to help mouse_20010508
#endif

        if (!PairwiseChiSquare((float)inferredMean, (float)inferredVariance,
                               (float)(- CGW_MISSED_OVERLAP), (float)1.0,
                               (LengthT *)NULL, &chiSquaredValue,
                               (float)PAIRWISECHI2THRESHOLD_CGW)) {
          failedInfer = TRUE;
          break;
        }
      }
      targetEdge->flags.bits.isInferredRemoved = TRUE;
      setEssentialEdgeStatus(targetEdge, FALSE);
      if (targetEdge->flags.bits.isInferred && targetEdge->flags.bits.isTentative) {
        DeleteGraphEdge(graph->ContigGraph, targetEdge);
      }
      numEssentialRemoved++;
      inferredEdgeIndex = AddGraphEdge(graph->ContigGraph, sourceCI->id, targetCI->id,
                                       NULLINDEX, NULLINDEX, // No fragments
                                       NULLINDEX, // No distance record
                                       inferredDistance, 1, // fudgeDistance???
                                       1.0, // quality
                                       inferredEdgeOrient,
                                       FALSE, // not unknown orientation
                                       isOverlapInferred,
                                       FALSE, // isAContainsB
                                       FALSE, // isBContainsA
                                       FALSE, // isTransChunk
                                       FALSE, // not extremalA
                                       FALSE, // not extremalB
                                       UNKNOWN_EDGE_STATUS,
                                       FALSE, // do not collect overlap
                                       TRUE); // insert into graph
      (*numInferredAdded)++;
      *addedEnd++ = inferredEdgeIndex;
      inferredEdge = GetGraphEdge(graph->ContigGraph, inferredEdgeIndex);
      AssertPtr(inferredEdge);
      numEssentialAdded++;
      inferredEdge->flags.bits.isInferred = TRUE;
      inferredEdge->flags.bits.isTentative = TRUE;
      setEssentialEdgeStatus(inferredEdge, TRUE);

#if 0
      fprintf(stderr,"* Added inferred GraphEdge v2 (%d,%d,%c) isAContainsB:%d isBContainsA:%d fragA %d fragB %d isInferred %d isRaw %d\n",
              inferredEdge->idA, inferredEdge->idB,
              inferredEdge->orient,
              inferredEdge->flags.bits.aContainsB,
              inferredEdge->flags.bits.bContainsA,
              inferredEdge->fragA,inferredEdge->fragB,
              inferredEdge->flags.bits.isInferred,
              inferredEdge->flags.bits.isRaw);
#endif


    }
    if (failedInfer) {

#ifdef INSTRUMENT_SMOOTHED
      // instrument the failed 'scaffold' before changing essential edges back
      InstrumentContigPath(graph, smoothed_si,
                           firstID, instFirstEnd, sourceCI->id);
      PrintScaffoldInstrumenter(graph, smoothed_si,
                                InstrumenterVerbose2,
                                "", instSmoothFailurefp);
#endif

      sourceCI->flags.bits.smoothSeenAlready = FALSE;
      sourceCI->smoothExpectedCID = NULLINDEX;
      for(addedPtr = addedInferred; addedPtr < addedEnd; addedPtr++) {
        CIEdgeT *addedEdge = GetGraphEdge(graph->ContigGraph, *addedPtr);
        if (!addedEdge->flags.bits.isDeleted) {
          DeleteGraphEdge(graph->ContigGraph, addedEdge);
        }
      }
      for(existingPtr = existingEdges; existingPtr < existingEnd; existingPtr++) {
        CIEdgeT *existingEdge = GetGraphEdge(graph->ContigGraph, *existingPtr);
        setEssentialEdgeStatus(existingEdge, existingEdge->flags.bits.wasEssential);
        //      existingEdge->flags.bits.isEssential =	existingEdge->flags.bits.wasEssential;
      }
      for(branchTarget = branchEdges + 1; branchTarget < branchEnd;
          branchTarget++) {
        CIEdgeT *targetEdge = GetGraphEdge(graph->ContigGraph, *branchTarget);
        ChunkInstanceT *targetCI;

        if (targetEdge->flags.bits.isDeleted)
          continue;

        targetCI = GetGraphNode(graph->ContigGraph, ((targetEdge->idA == thisCI->id) ? targetEdge->idB : targetEdge->idA));
        targetCI->flags.bits.smoothSeenAlready = FALSE;
        targetCI->smoothExpectedCID = NULLINDEX;
        targetEdge->flags.bits.isInferredRemoved = FALSE;
        setEssentialEdgeStatus(targetEdge, TRUE);
        //      targetEdge->flags.bits.isEssential = TRUE;
        //      assert(!isSloppyEdge(targetEdge));

      }

      safe_free(branchEdges);
      safe_free(existingEdges);
      safe_free(addedInferred);
      safe_free(branchTargets);

      return(FALSE);
    }else{
      CDS_CID_t sourceEnd;
      if (end == B_END) {
        thisCI->numEssentialB -= numEssentialRemoved;
      }else{
        thisCI->numEssentialA -= numEssentialRemoved;
      }
      if ((sourceEdgeOrient == AB_AB) || (sourceEdgeOrient == BA_AB)) {
        sourceEnd = B_END;
        sourceCI->numEssentialB += numEssentialAdded;
      }else{
        sourceEnd = A_END;
        sourceCI->numEssentialA += numEssentialAdded;
      }
      if (!RecursiveSmoothWithInferredEdges(graph, sourceCI, sourceEnd,
                                            numInferredAdded,
                                            firstID, instFirstEnd, lastID,
                                            smoothed_si,
                                            instSmoothSuccessfp,
                                            instSmoothFailurefp)) {
        if (end == B_END) {
          thisCI->numEssentialB += numEssentialRemoved;
        }else{
          thisCI->numEssentialA += numEssentialRemoved;
        }
        if (sourceEnd == B_END) {
          sourceCI->numEssentialB -= numEssentialAdded;
        }else{
          sourceCI->numEssentialA -= numEssentialAdded;
        }
        sourceCI->flags.bits.smoothSeenAlready = FALSE;
        sourceCI->smoothExpectedCID = NULLINDEX;
        for(addedPtr = addedInferred; addedPtr < addedEnd; addedPtr++) {
          CIEdgeT *addedEdge = GetGraphEdge(graph->ContigGraph, *addedPtr);
          if (!addedEdge->flags.bits.isDeleted) {
            DeleteGraphEdge(graph->ContigGraph, addedEdge);
          }
        }
        for(existingPtr = existingEdges; existingPtr < existingEnd;
            existingPtr++) {
          CIEdgeT *existingEdge = GetGraphEdge(graph->ContigGraph, *existingPtr);
          setEssentialEdgeStatus(existingEdge, existingEdge->flags.bits.wasEssential);
          //	existingEdge->flags.bits.isEssential =
          //	  existingEdge->flags.bits.wasEssential;
        }
        for(branchTarget = branchEdges + 1; branchTarget < branchEnd;
            branchTarget++) {
          CIEdgeT *targetEdge = GetGraphEdge(graph->ContigGraph, *branchTarget);
          ChunkInstanceT *targetCI;

          if (targetEdge->flags.bits.isDeleted)
            continue;

          targetCI = GetGraphNode(graph->ContigGraph, ((targetEdge->idA == thisCI->id) ? targetEdge->idB : targetEdge->idA));
          targetCI->flags.bits.smoothSeenAlready = FALSE;
          targetCI->smoothExpectedCID = NULLINDEX;
          targetEdge->flags.bits.isInferredRemoved = FALSE;
          assert(!isSloppyEdge(targetEdge));
          setEssentialEdgeStatus(targetEdge, TRUE);
          //	targetEdge->flags.bits.isEssential = TRUE;
        }

        safe_free(branchEdges);
        safe_free(existingEdges);
        safe_free(addedInferred);
        safe_free(branchTargets);

        return(FALSE);
      }else{
        sourceCI->flags.bits.smoothSeenAlready = FALSE;
        sourceCI->smoothExpectedCID = NULLINDEX;
        if (end == B_END) {
          thisCI->essentialEdgeB = *branchSource;
        }else{
          thisCI->essentialEdgeA = *branchSource;
        }
        if ((sourceEdgeOrient == AB_AB) || (sourceEdgeOrient == BA_AB)) {
          sourceCI->essentialEdgeA = *branchSource;
        }else{
          sourceCI->essentialEdgeB = *branchSource;
        }
        for(addedPtr = addedInferred; addedPtr < addedEnd; addedPtr++) {
          CIEdgeT *addedEdge = GetGraphEdge(graph->ContigGraph, *addedPtr);
          if (!addedEdge->flags.bits.isDeleted) {
            addedEdge->flags.bits.isTentative = FALSE;
          }
        }
        for(existingPtr = existingEdges; existingPtr < existingEnd;
            existingPtr++) {
          CIEdgeT *existingEdge = GetGraphEdge(graph->ContigGraph, *existingPtr);
          if (existingEdge->flags.bits.wasEssential) {
            ChunkInstanceT *targetCI =
              GetGraphNode(graph->ContigGraph,
                           ((existingEdge->idA == sourceCI->id) ?
                            existingEdge->idB : existingEdge->idA));
            ChunkOrientationType existingEdgeOrient =
              GetEdgeOrientationWRT(existingEdge, sourceCI->id);
            // Fix numEssential
            if ((existingEdgeOrient == AB_AB) || (existingEdgeOrient == BA_AB)) {
              targetCI->numEssentialA--;
            }else{
              targetCI->numEssentialB--;
            }
          }
        }

        safe_free(branchEdges);
        safe_free(existingEdges);
        safe_free(addedInferred);
        safe_free(branchTargets);

        return(TRUE);
      }
    }
  }
}

static
int
CountEssentialEdgesOneEnd(ScaffoldGraphT *graph, ChunkInstanceT *thisCI,
                          int end, CDS_CID_t *essentialEdge) {

  int count = 0;
  GraphEdgeIterator edges;
  CIEdgeT *edge;
  InitGraphEdgeIterator(graph->ContigGraph, thisCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
  *essentialEdge = (CDS_CID_t)NULLINDEX;
  while((edge = NextGraphEdgeIterator(&edges))!= NULL) {
    if (!edge->flags.bits.isActive || !edge->flags.bits.isConfirmed ||
        !edge->flags.bits.isUniquetoUnique ||
        edge->flags.bits.isTransitivelyRemoved ||
        edge->flags.bits.isRedundantRemoved) {
      // Inactive, unconfirmed, and edges between nonUnique CIs are not essential.
      setEssentialEdgeStatus(edge, FALSE);
      continue;
    }
    setEssentialEdgeStatus(edge, TRUE);
    count++;
    *essentialEdge = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->ContigGraph->edges, edge);
  }
  return(count);
}



static
void
SmoothWithInferredEdges(ScaffoldGraphT *graph,
                        int markShakyBifurcations) {

  GraphNodeIterator nodes;
  ContigT *thisCI;
  int cnt = 1;
  int smooth_success;
  int numInferredAdded;

  InitGraphNodeIterator(&nodes, graph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY );
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL) {
    thisCI->smoothExpectedCID = NULLINDEX;
    thisCI->flags.bits.smoothSeenAlready = FALSE;

    // Do the A_END edges.
    thisCI->numEssentialA = CountEssentialEdgesOneEnd(graph, thisCI, A_END, &(thisCI->essentialEdgeA));
    assert((thisCI->numEssentialA == 0) || (thisCI->essentialEdgeA != NULLINDEX));

    // Do the B_END edges.
    thisCI->numEssentialB = CountEssentialEdgesOneEnd(graph, thisCI, B_END, &(thisCI->essentialEdgeB));
    assert((thisCI->numEssentialB == 0) || (thisCI->essentialEdgeB != NULLINDEX));
  }

  if (!markShakyBifurcations) {
    // for instrumenting successfully/failed smoothings - start & end contigs
    CDS_CID_t firstID;
    CDS_CID_t lastID;
    int instFirstEnd;
    // instrumenter to use
#ifdef INSTRUMENT_TRANS_REDUCED
    ScaffoldInstrumenter * smoothed_si =
      CreateScaffoldInstrumenter(graph, INST_OPT_INTER_MATES);
#else
    ScaffoldInstrumenter * smoothed_si = NULL;
#endif
    // where to write them
    FILE * instSmoothSuccessfp = NULL;
    FILE * instSmoothFailurefp = NULL;
#ifdef INSTRUMENT_CGW
    ScaffoldInstrumenter * si_before =
      CreateScaffoldInstrumenter(graph, INST_OPT_INTER_MATES);
    ScaffoldInstrumenter * si_after =
      CreateScaffoldInstrumenter(graph, INST_OPT_INTER_MATES);
    int num_before, num_after;
#endif
#ifdef INSTRUMENT_SMOOTHED
    {
      int file_count;
      char filename[1024];
      for(file_count = 0; file_count < 1000; file_count++)
        {
          sprintf(filename, "smoothSuccess_%05d.itxt", file_count);
          instSmoothSuccessfp = fopen(filename, "r");
          if (instSmoothSuccessfp == NULL)
            {
              instSmoothSuccessfp = fopen(filename, "w");
              sprintf(filename, "smoothFailure_%05d.itxt", file_count);
              instSmoothFailurefp = fopen(filename, "w");
              break;
            }
          else
            fclose(instSmoothSuccessfp);
        }
      assert(instSmoothSuccessfp != NULL && instSmoothFailurefp != NULL);
    }
#endif
#ifdef INSTRUMENT_TRANS_REDUCED
    assert(smoothed_si != NULL);
#endif

    InitGraphNodeIterator(&nodes, graph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY );
    while((thisCI = NextGraphNodeIterator(&nodes)) != NULL) {
      /* Substitute inferred edges where the graph branches in an attempt to
         eliminate as many branches as possible. */
      if (thisCI->numEssentialA > 1) {
        thisCI->flags.bits.smoothSeenAlready = TRUE;

#ifdef INSTRUMENT_CGW
        num_before = InstrumentContigEnd(graph, si_before, thisCI, A_END);
#endif
        firstID = thisCI->id;
        lastID = NULLINDEX;
        instFirstEnd = A_END;

        numInferredAdded = 0;
        smooth_success =
          RecursiveSmoothWithInferredEdges(graph, thisCI, A_END,
                                           &numInferredAdded,
                                           firstID, instFirstEnd, &lastID,
                                           smoothed_si,
                                           instSmoothSuccessfp,
                                           instSmoothFailurefp);
        thisCI->flags.bits.smoothSeenAlready = FALSE;

#ifdef INSTRUMENT_CGW
        if (smooth_success == TRUE)
          {
            num_after = InstrumentContigEnd(graph, si_after, thisCI, A_END);
            fprintf(stderr, "Smoothed CI "F_CID " on A_END\n",
                    thisCI->id);
            fprintf(stderr,
                    "contigs: %d before (may include duplicates), %d after\n",
                    num_before, num_after);
            if (CompareMateInstrumenters(&(si_before->mates),
                                         &(si_after->mates),
                                         InstrumenterVerbose5,
                                         stderr) ==
                InstrumenterWorse)
              {
                int num_partial;
                for(num_partial = num_after - 1; num_partial > 1; num_partial--)
                  {
                    // assumes there is one essential edge off end
                    InstrumentContigEndPartial(graph, si_after, thisCI, A_END,
                                               num_partial);
                    fprintf(stderr, "re-instrumenting %d contigs\n",
                            num_partial);
                    if (CompareMateInstrumenters(&(si_before->mates),
                                                 &(si_after->mates),
                                                 InstrumenterVerbose5,
                                                 stderr) !=
                        InstrumenterWorse)
                      {
                        break;
                      }
                  }
              }
          }
#endif
      } // if (thisCI->numEssentialA > 1) {
      if (thisCI->numEssentialB > 1) {
        thisCI->flags.bits.smoothSeenAlready = TRUE;
#ifdef INSTRUMENT_CGW
        num_before = InstrumentContigEnd(graph, si_before, thisCI, B_END);
#endif
        firstID = thisCI->id;
        lastID = NULLINDEX;
        instFirstEnd = B_END;

        numInferredAdded = 0;
        smooth_success =
          RecursiveSmoothWithInferredEdges(graph, thisCI, B_END,
                                           &numInferredAdded,
                                           firstID, instFirstEnd, &lastID,
                                           smoothed_si,
                                           instSmoothSuccessfp,
                                           instSmoothFailurefp);
        thisCI->flags.bits.smoothSeenAlready = FALSE;

#ifdef INSTRUMENT_CGW
        if (smooth_success == TRUE)
          {
            num_after = InstrumentContigEnd(graph, si_after, thisCI, B_END);
            fprintf(stderr, "Smoothed CI "F_CID " on B_END\n",
                    thisCI->id);
            fprintf(stderr,
                    "contigs: %d before (may include duplicates), %d after\n",
                    num_before, num_after);
            if (CompareMateInstrumenters(&(si_before->mates),
                                         &(si_after->mates),
                                         InstrumenterVerbose5,
                                         stderr) ==
                InstrumenterWorse)
              {
                int num_partial;
                for(num_partial = num_after - 1; num_partial > 1; num_partial--)
                  {
                    // assumes there is one essential edge off end
                    InstrumentContigEndPartial(graph, si_after, thisCI, B_END,
                                               num_partial);
                    fprintf(stderr, "re-instrumenting %d contigs\n",
                            num_partial);
                    if (CompareMateInstrumenters(&(si_before->mates),
                                                 &(si_after->mates),
                                                 InstrumenterVerbose5,
                                                 stderr) !=
                        InstrumenterWorse)
                      {
                        break;
                      }
                  }
              }
          }
#endif
      } //       if (thisCI->numEssentialB > 1) {
    } // while(thisCI = NextGraphNodeIterator(&nodes)) {
#ifdef INSTRUMENT_CGW
    DestroyScaffoldInstrumenter(si_before);
    DestroyScaffoldInstrumenter(si_after);
#endif
#ifdef INSTRUMENT_SMOOTHED
    fclose(instSmoothSuccessfp);
    fclose(instSmoothFailurefp);
#endif
#ifdef INSTRUMENT_TRANS_REDUCED
    DestroyScaffoldInstrumenter(smoothed_si);
#endif
  } // if (!markShakyBifurcations) {
  return;
}


static
void
SymmetricNeighbors(ChunkInstanceT *thisCI, ChunkOrientationType orient, CDS_CID_t edgeID) {

  if ((orient == AB_AB) || (orient == BA_AB)) {
    // A end of this CI
    assert(thisCI->numEssentialA == 1);
    assert(edgeID == thisCI->essentialEdgeA);
  } else {
    // B end of this CI
    assert(thisCI->numEssentialB == 1);
    assert(edgeID == thisCI->essentialEdgeB);
  }
}


static
int
NeighborBranches(ChunkInstanceT *thisCI,
                 ChunkOrientationType orient,
                 CDS_CID_t edgeID) {

  if ((orient == AB_AB) || (orient == BA_AB)) {
    // A end of this CI
    if (thisCI->numEssentialA != 1) {
      return(TRUE);
    } else {
      assert(edgeID == thisCI->essentialEdgeA);
      return(FALSE);
    }
  } else {
    // B end of this CI
    assert(thisCI->numEssentialB >= 1);
    if (thisCI->numEssentialB != 1) {
      return(TRUE);
    } else {
      assert(edgeID == thisCI->essentialEdgeB);
      return(FALSE);
    }
  }
}


static
void
DetectScaffoldCycles(CDS_CID_t *currentScaffoldID) {
  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL) {
    ChunkInstanceT *AendCI, *BendCI, *startingCI;
    CIEdgeT *edgeA, *edgeB, *edge;
    ChunkInstanceT *neighbor;
    ChunkOrient orientCI;

    // the following test identifies and skips CIs that were
    // assigned to a scaffold by LabelCIScaffolds, or were
    // already operated on below.  These are
    // by construction not in cycles, either because of
    // LabelCIScaffolds or because any cycle they might have
    // originally been in have already been busted up (see below).
    // (In LabelCIScaffolds, all bifurcations were
    // removed before calling LabelCIScaffolds; LabelCIScaffolds
    // works by finding a CI that doesn't have edges off both
    // ends and labeling the whole resulting component).

    if (thisCI->scaffoldID != NULLINDEX) {
      // This CI has already been placed in a Scaffold.
      thisCI->scaffoldID = NULLINDEX;
      continue;
    }

    // what this leaves behind is cycles; so any remaining unique CIs
    // must be in a cycle.
    assert(thisCI->essentialEdgeA != NULLINDEX);
    edgeA = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeA);
    AendCI = GetGraphNode(ScaffoldGraph->ContigGraph,
                          (edgeA->idA == thisCI->id) ? edgeA->idB : edgeA->idA);
    assert(thisCI->essentialEdgeB != NULLINDEX);
    edgeB = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeB);
    BendCI = GetGraphNode(ScaffoldGraph->ContigGraph,
                          (edgeB->idA == thisCI->id) ? edgeB->idB : edgeB->idA);
    fprintf(stderr, "*** Found a Cycle in Essential Edge Graph\ncId:"F_CID " Aend:"F_CID " Bend:"F_CID "\n",
            thisCI->id, AendCI->id, BendCI->id);

    orientCI = A_B;
    edge = edgeB;
    neighbor = BendCI;
    startingCI = thisCI;

    // for all nodes in a cycle (terminated by getting back to the
    // starting CI), remove their essential edges ... i.e. completely
    // break the cycle into singleton bits
    do{
      ChunkOrientationType edgeOrient = GetEdgeOrientationWRT(edge,
                                                              thisCI->id);
      // Dump
      DumpContig(stderr, ScaffoldGraph, thisCI,  FALSE);


      assert(neighbor != NULL);
      if (orientCI == A_B) {
        if (edgeOrient == AB_AB) {
          orientCI = A_B;
        }else if (edgeOrient == AB_BA) {
          orientCI = B_A;
        }else{
          assert(0);
        }
      }else{// orientCI == B_A
        if (edgeOrient == BA_AB) {
          orientCI = A_B;
        }else if (edgeOrient == BA_BA) {
          orientCI = B_A;
        }else{
          assert(0);
        }
      }
      thisCI->essentialEdgeA = NULLINDEX;
      thisCI->essentialEdgeB = NULLINDEX;

      // the actual id assigned here is immaterial; what matters is
      // that all contigs in the cycle get some value not NULLINDEX
      // so that as each one is visited it will pass the test for
      // id != NULLINDEX, have the id set to NULLINDEX, and then
      // continue (see top of outer while loop)
      thisCI->scaffoldID = *currentScaffoldID;

      thisCI = neighbor;

      if (orientCI == A_B) {
        if (thisCI->essentialEdgeB != NULLINDEX) {
          edge = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeB);
          neighbor = GetGraphNode(ScaffoldGraph->ContigGraph,
                                  (edge->idA == thisCI->id) ? edge->idB : edge->idA);
        }else{// End of Scaffold
          edge = (CIEdgeT *)NULL;
          neighbor = (ChunkInstanceT *)NULL;
        }
      }else{// orientCI == B_A
        if (thisCI->essentialEdgeA != NULLINDEX) {
          edge = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeA);
          neighbor = GetGraphNode(ScaffoldGraph->ContigGraph,
                                  (edge->idA == thisCI->id) ? edge->idB : edge->idA);
        }else{// End of Scaffold
          edge = (CIEdgeT *)NULL;
          neighbor = (ChunkInstanceT *)NULL;
        }
      }
    }while(thisCI != startingCI);
    // take the original node that gave us access to the cycle and
    // mark it as processed
    startingCI->scaffoldID = NULLINDEX;
    (*currentScaffoldID)++;
  }
}

static
void
LabelCIScaffolds(CDS_CID_t *currentScaffoldID) {
  /* Label which scaffold a CI would go into but do not create the
     scaffold as a way to detect cycles. */
  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL) {
    ChunkInstanceT *AendCI, *BendCI;
    CIEdgeT *edgeA, *edgeB, *edge;
    ChunkInstanceT *neighbor;
    ChunkOrient orientCI;

    if (thisCI->scaffoldID != NULLINDEX)
      // This CI has already been placed in a Scaffold.
      continue;

    if (thisCI->essentialEdgeA != NULLINDEX) {
      edgeA = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeA);
      AendCI = GetGraphNode(ScaffoldGraph->ContigGraph,
                            (edgeA->idA == thisCI->id) ? edgeA->idB : edgeA->idA);
    }else{
      edgeA = (CIEdgeT *)NULL;
      AendCI = (ChunkInstanceT *)NULL;
    }
    if (thisCI->essentialEdgeB != NULLINDEX) {
      edgeB = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeB);
      BendCI = GetGraphNode(ScaffoldGraph->ContigGraph,
                            (edgeB->idA == thisCI->id) ? edgeB->idB : edgeB->idA);
    }else{
      edgeB = (CIEdgeT *)NULL;
      BendCI = (ChunkInstanceT *)NULL;
    }
    if ((AendCI != (ChunkInstanceT *)NULL) &&
        (BendCI != (ChunkInstanceT *)NULL)) {
      // This CI is not a starting point for a Scaffold.
      continue;
    }


    if (BendCI != (ChunkInstanceT *)NULL) {
      orientCI = A_B;
      edge = edgeB;
      neighbor = BendCI;
    }else if (AendCI != (ChunkInstanceT *)NULL) {
      orientCI = B_A;
      edge = edgeA;
      neighbor = AendCI;
    }else{// Singleton Scaffold
      orientCI = A_B;
      edge = (CIEdgeT *)NULL;
      neighbor = (ChunkInstanceT *)NULL;
    }

    thisCI->scaffoldID = *currentScaffoldID;
    while(neighbor != (ChunkInstanceT *)NULL) {
      ChunkOrientationType edgeOrient = GetEdgeOrientationWRT(edge,
							      thisCI->id);
      if (orientCI == A_B) {
	if (edgeOrient == AB_AB) {
	  orientCI = A_B;
	}else if (edgeOrient == AB_BA) {
	  orientCI = B_A;
	}else{
	  assert(0);
	}
      }else{// orientCI == B_A
	if (edgeOrient == BA_AB) {
	  orientCI = A_B;
	}else if (edgeOrient == BA_BA) {
	  orientCI = B_A;
	}else{
	  assert(0);
	}
      }
      thisCI = neighbor;
      thisCI->scaffoldID = *currentScaffoldID;
      if (orientCI == A_B) {
	if (thisCI->essentialEdgeB != NULLINDEX) {
	  edge = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeB);
	  neighbor = GetGraphNode(ScaffoldGraph->ContigGraph,
				  (edge->idA == thisCI->id) ? edge->idB : edge->idA);
	}else{// End of Scaffold
	  edge = (CIEdgeT *)NULL;
	  neighbor = (ChunkInstanceT *)NULL;
	}
      }else{// orientCI == B_A
	if (thisCI->essentialEdgeA != NULLINDEX) {
	  edge = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeA);
	  neighbor = GetGraphNode(ScaffoldGraph->ContigGraph,
				  (edge->idA == thisCI->id) ? edge->idB : edge->idA);
	}else{// End of Scaffold
	  edge = (CIEdgeT *)NULL;
	  neighbor = (ChunkInstanceT *)NULL;
	}
      }
    }
    (*currentScaffoldID)++;
  }
}

static
void
ActuallyInsertCIsIntoScaffolds(CDS_CID_t *currentScaffoldID) {
  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;
  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY);

  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL) {
    ChunkInstanceT *AendCI, *BendCI;
    CIEdgeT *edgeA, *edgeB, *edge;
    ChunkInstanceT *neighbor;
    LengthT NullLength = {0.0, 0.0};
    LengthT aEndOffset, bEndOffset, currentOffset;
    CIScaffoldT CIScaffold;
    ChunkOrient orientCI;

    memset(&CIScaffold, 0, sizeof(CIScaffoldT));

    if (thisCI->scaffoldID != NULLINDEX)
      // This CI has already been placed in a Scaffold.
      continue;

    if (thisCI->essentialEdgeA != NULLINDEX) {
      edgeA = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeA);
      assert(!isSloppyEdge(edgeA));
      AendCI = GetGraphNode(ScaffoldGraph->ContigGraph, (edgeA->idA == thisCI->id) ? edgeA->idB : edgeA->idA);
      SymmetricNeighbors(AendCI, GetEdgeOrientationWRT(edgeA, thisCI->id),thisCI->essentialEdgeA );
    }else{
      edgeA = (CIEdgeT *)NULL;
      AendCI = (ChunkInstanceT *)NULL;
    }
    if (thisCI->essentialEdgeB != NULLINDEX) {
      edgeB = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeB);
      assert(!isSloppyEdge(edgeB));
      assert(edgeB->idA == thisCI->id || edgeB->idB == thisCI->id);
      BendCI = GetGraphNode(ScaffoldGraph->ContigGraph, (edgeB->idA == thisCI->id) ? edgeB->idB : edgeB->idA);
      SymmetricNeighbors(BendCI, GetEdgeOrientationWRT(edgeB, thisCI->id),thisCI->essentialEdgeB );
    }else{
      edgeB = (CIEdgeT *)NULL;
      BendCI = (ChunkInstanceT *)NULL;
    }

    if ((AendCI != (ChunkInstanceT *)NULL) &&
        (BendCI != (ChunkInstanceT *)NULL))
      // This CI is not a starting point for a Scaffold.
      continue;

    if (BendCI != (ChunkInstanceT *)NULL) {
      orientCI = A_B;
      edge = edgeB;
      neighbor = BendCI;
    }else if (AendCI != (ChunkInstanceT *)NULL) {
      orientCI = B_A;
      edge = edgeA;
      neighbor = AendCI;
    }else{// Singleton Scaffold
      orientCI = A_B;
      edge = (CIEdgeT *)NULL;
      neighbor = (ChunkInstanceT *)NULL;
    }
    if (*currentScaffoldID == 404) {
      fprintf(stderr,"* Scaffold 404\n");
    }
    InitializeScaffold(&CIScaffold, REAL_SCAFFOLD);
    CIScaffold.info.Scaffold.AEndCI = NULLINDEX;
    CIScaffold.info.Scaffold.BEndCI = NULLINDEX;
    CIScaffold.info.Scaffold.numElements = 0;
    CIScaffold.info.Scaffold.leastSquareError = 0.0;
    CIScaffold.edgeHead = NULLINDEX;
    CIScaffold.bpLength = NullLength;
    thisCI->scaffoldID = CIScaffold.id = *currentScaffoldID;
    CIScaffold.flags.bits.isDead = FALSE;
    CIScaffold.numEssentialA = CIScaffold.numEssentialB = 0;
    CIScaffold.essentialEdgeB = CIScaffold.essentialEdgeA = NULLINDEX;
    AppendCIScaffoldT(ScaffoldGraph->CIScaffolds, &CIScaffold);
    assert(*currentScaffoldID == (GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds) - 1));

    aEndOffset = bEndOffset = currentOffset = NullLength;

    if (orientCI == A_B) {
      aEndOffset = NullLength;
      bEndOffset = thisCI->bpLength;
    }else{// orientCI == B_A
      bEndOffset = NullLength;
      aEndOffset = thisCI->bpLength;
    }
    /* DON'T ContigNOW!!! */
    InsertCIInScaffold(ScaffoldGraph, thisCI->id, *currentScaffoldID,
                       aEndOffset, bEndOffset,  TRUE /* Should be FALSE */, FALSE);
    currentOffset = thisCI->bpLength;

    while(neighbor != (ChunkInstanceT *)NULL) {
      ChunkOrientationType edgeOrient = GetEdgeOrientationWRT(edge,
                                                              thisCI->id);
      if (orientCI == A_B) {
        if (edgeOrient == AB_AB) {
          orientCI = A_B;
        }else if (edgeOrient == AB_BA) {
          orientCI = B_A;
        }else{
          assert(0);
        }
      }else{// orientCI == B_A
        if (edgeOrient == BA_AB) {
          orientCI = A_B;
        }else if (edgeOrient == BA_BA) {
          orientCI = B_A;
        }else{
          assert(0);
        }
      }
      thisCI = neighbor;
      thisCI->scaffoldID = *currentScaffoldID;
      currentOffset.mean += edge->distance.mean;
      currentOffset.variance += edge->distance.variance;
      if (orientCI == A_B) {
        aEndOffset = currentOffset;
      }else{// orientCI == B_A
        bEndOffset = currentOffset;
      }
      currentOffset.mean += thisCI->bpLength.mean;
      currentOffset.variance += thisCI->bpLength.variance;
      if (orientCI == A_B) {
        bEndOffset = currentOffset;
      }else{// orientCI == B_A
        aEndOffset = currentOffset;
      }
      /* Don't Contig NOW!!! */
      InsertCIInScaffold(ScaffoldGraph, thisCI->id, *currentScaffoldID,
                         aEndOffset, bEndOffset,  TRUE /* Should be FALSE */, FALSE);
      if (orientCI == A_B) {
        if (thisCI->essentialEdgeB != NULLINDEX) {
          edge = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeB);
          assert(!isSloppyEdge(edge));
          neighbor = GetGraphNode(ScaffoldGraph->ContigGraph,
                                  (edge->idA == thisCI->id) ? edge->idB : edge->idA);
        }else{// End of Scaffold
          edge = (CIEdgeT *)NULL;
          neighbor = (ChunkInstanceT *)NULL;
        }
      }else{// orientCI == B_A
        if (thisCI->essentialEdgeA != NULLINDEX) {
          edge = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeA);
          assert(!isSloppyEdge(edge));
          neighbor = GetGraphNode(ScaffoldGraph->ContigGraph,
                                  (edge->idA == thisCI->id) ? edge->idB : edge->idA);
        }else{// End of Scaffold
          edge = (CIEdgeT *)NULL;
          neighbor = (ChunkInstanceT *)NULL;
        }
      }
    }

    {/***** Check that scaffold is connected ****/
      CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, *currentScaffoldID);
      if (!IsScaffoldInternallyConnected(ScaffoldGraph,  scaffold, ALL_TRUSTED_EDGES)) {
        fprintf(stderr,"* Scaffold "F_CID
                " is DISCONNECTED IMMEDIATELY AFTER INITIAL CONSTRUCTION!!!!\n",
                *currentScaffoldID);
        DumpACIScaffold(stderr,ScaffoldGraph, scaffold, FALSE);
        assert(0);
      }
    }
    (*currentScaffoldID)++;
  }
}


static
void
MarkBifurcations(ScaffoldGraphT *graph) {
  GraphNodeIterator nodes;
  ContigT *thisCI;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL) {

    thisCI->scaffoldID = NULLINDEX;

    switch(thisCI->numEssentialA) {
      case 0:
        break;
      case 1:
        if (thisCI->essentialEdgeA != NULLINDEX) {
          CIEdgeT *edge = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeA);
          if (NeighborBranches(GetGraphNode(ScaffoldGraph->ContigGraph, ((edge->idA == thisCI->id) ? edge->idB : edge->idA)),
                               GetEdgeOrientationWRT(edge, thisCI->id),
                               thisCI->essentialEdgeA)) {
            thisCI->essentialEdgeA = NULLINDEX;
          }
        }
        break;
      default:
        thisCI->essentialEdgeA = NULLINDEX;
        break;
    }

    switch(thisCI->numEssentialB) {
      case 0:
        break;
      case 1:
        if (thisCI->essentialEdgeB != NULLINDEX) {
          CIEdgeT *edge = GetGraphEdge(ScaffoldGraph->ContigGraph, thisCI->essentialEdgeB);
          if (NeighborBranches(GetGraphNode(ScaffoldGraph->ContigGraph,
                                            ((edge->idA == thisCI->id) ?
                                             edge->idB : edge->idA)),
                               GetEdgeOrientationWRT(edge, thisCI->id),
                               thisCI->essentialEdgeB)) {
            thisCI->essentialEdgeB = NULLINDEX;
          }
        }
	break;
      default:
        thisCI->essentialEdgeB = NULLINDEX;
        break;
    }
  }
}

static
CDS_CID_t
MarkShakyBifurcations(ScaffoldGraphT *graph) {
  int numShaky = 0;
  GraphNodeIterator nodes;
  ContigT *contig;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((contig = NextGraphNodeIterator(&nodes)) != NULL) {
    if (IsShakyContigAtScaffoldEnd(contig)) {
      SetNodeType(GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI),
                  UNRESOLVEDCHUNK_CGW);

      contig->flags.bits.isUnique = 0;
      contig->scaffoldID = NULLINDEX;
      contig->AEndNext = contig->BEndNext = NULLINDEX;

      numShaky++;
    }
  }

  return numShaky;
}



static
void
CreateScaffolds(ScaffoldGraphT *graph) {

  CDS_CID_t currentScaffoldID;

  /* Recycle the CIScaffolds VA */
  ResetNodeCGW_T(graph->ScaffoldGraph->nodes);

  MarkBifurcations(graph);

  currentScaffoldID = 0;
  LabelCIScaffolds(&currentScaffoldID);

  /* We're looking for cycles where some of the Unique CIs have not been
     labeled as belonging to a scaffold. */

  DetectScaffoldCycles(&currentScaffoldID);

  // Adjust scaffold labeling, based on instrumenting
  // AdjustCIScaffoldLabels(graph, &currentScaffoldID);

  /* Actually insert CIs into scaffolds. */
  currentScaffoldID = 0;
  ActuallyInsertCIsIntoScaffolds(&currentScaffoldID);

  graph->numLiveScaffolds = GetNumCIScaffoldTs(graph->CIScaffolds);
  assert(graph->numLiveScaffolds == currentScaffoldID);
}

static
void
DeleteInferredEdges(ScaffoldGraphT *graph) {

  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;

  //  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY);
  // An inferred edge may exist between two shaky nodes that have
  // been demoted
  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL) {
    GraphEdgeIterator edges;
    CIEdgeT *edge;

    InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, thisCI->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
    while((edge = NextGraphEdgeIterator(&edges)) != NULL) {
      if (edge->flags.bits.isInferred) {
	DeleteGraphEdge(ScaffoldGraph->ContigGraph, edge);
      }
    }
  }
}



//  This computes edges between the CIs in a scaffold, inferred by the relative positions of the CIs
//  in the scaffold.
//
static
void
AddScaffoldInferredEdges(ScaffoldGraphT *graph) {
  GraphNodeIterator scaffolds;
  CDS_CID_t sid;
  CIScaffoldT *scaffold;
  fprintf(stderr,"* AddScaffoldInferredEdges   scaffolds = %d\n",
	  (int) GetNumGraphNodes(graph->ScaffoldGraph));

  InitGraphNodeIterator(&scaffolds, graph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL) {
    CIScaffoldTIterator CIs;
    ChunkInstanceT *thisCI, *prevCI;
    sid = scaffold->id;

    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
    prevCI = NextCIScaffoldTIterator(&CIs);
    while((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL) {
      LengthT inferredDistance;
      ChunkOrientationType inferredEdgeOrient;
      CDS_CID_t inferredEdgeIndex;
      CIEdgeT *inferredEdge;
      int isOverlapInferred = FALSE;

      assert(!thisCI->flags.bits.isDead &&
             !prevCI->flags.bits.isDead);

      if (GetNodeOrient(prevCI) == A_B) {
	if (GetNodeOrient(thisCI) == A_B) {
	  inferredEdgeOrient = AB_AB;
	  inferredDistance.mean = thisCI->offsetAEnd.mean - prevCI->offsetBEnd.mean;
	  inferredDistance.variance = thisCI->offsetAEnd.variance - prevCI->offsetBEnd.variance;
	}else{// GetNodeOrient(thisCI) == B_A
	  inferredEdgeOrient = AB_BA;
	  inferredDistance.mean = thisCI->offsetBEnd.mean - prevCI->offsetBEnd.mean;
	  inferredDistance.variance = thisCI->offsetBEnd.variance - prevCI->offsetBEnd.variance;
	}
      }else{// GetNodeOrient(prevCI) == B_A
	if (GetNodeOrient(thisCI) == A_B) {
	  inferredEdgeOrient = BA_AB;
	  inferredDistance.mean = thisCI->offsetAEnd.mean - prevCI->offsetAEnd.mean;
	  inferredDistance.variance = thisCI->offsetAEnd.variance - prevCI->offsetAEnd.variance;
	}else{// GetNodeOrient(thisCI) == B_A
	  inferredEdgeOrient = BA_BA;
	  inferredDistance.mean = thisCI->offsetBEnd.mean - prevCI->offsetAEnd.mean;
	  inferredDistance.variance = thisCI->offsetBEnd.variance - prevCI->offsetAEnd.variance;
	}
      }

      if ( inferredDistance.variance > 0.0 ) {// SAK HACK!
        inferredEdgeIndex = AddGraphEdge(ScaffoldGraph->ContigGraph, prevCI->id, thisCI->id,
                                         NULLINDEX, NULLINDEX, // No fragments
                                         NULLINDEX, // No distance record
                                         inferredDistance, 1, // fudgeDistance???
                                         1.0, // quality
                                         inferredEdgeOrient,
                                         FALSE, // not unknown orientation
                                         isOverlapInferred,
                                         FALSE, // isAContainsB
                                         FALSE, // isBContainsA
                                         FALSE, // isTransChunk
                                         FALSE, // not extremalA
                                         FALSE, // not extremalB
                                         UNKNOWN_EDGE_STATUS,
                                         FALSE, // do not collect overlap
                                         TRUE); // insert into graph
        inferredEdge = GetGraphEdge(ScaffoldGraph->ContigGraph, inferredEdgeIndex);
        inferredEdge->flags.bits.isInferred = TRUE;
        inferredEdge->flags.bits.isTentative = FALSE;
        inferredEdge->flags.bits.isActive = TRUE;
        inferredEdge->flags.bits.isConfirmed = TRUE;
        inferredEdge->flags.bits.isLeastSquares = TRUE;

#if 0
        fprintf(stderr,"* Added inferred GraphEdge v1 (%d,%d,%c) isAContainsB:%d isBContainsA:%d fragA %d fragB %d isInferred %d isRaw %d\n",
                inferredEdge->idA, inferredEdge->idB,
                inferredEdge->orient,
                inferredEdge->flags.bits.aContainsB,
                inferredEdge->flags.bits.bContainsA,
                inferredEdge->fragA,inferredEdge->fragB,
                inferredEdge->flags.bits.isInferred,
                inferredEdge->flags.bits.isRaw);
#endif

      }else{
	// This is a containment edge
	fprintf(stderr,"* Did NOT add inferred edge ("F_CID ","F_CID ",%c) offsets [%g,%g] [%g,%g]... a containment?\n",
		prevCI->id, thisCI->id, inferredEdgeOrient,
		thisCI->offsetAEnd.mean, thisCI->offsetBEnd.mean,
		prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean);
      }
      prevCI = thisCI;
    }
  }
  return;
}


static
void
ResetEdgeStatus(ScaffoldGraphT *graph) {

  CDS_CID_t cid;
  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL) {
    GraphEdgeIterator edges;
    CIEdgeT *edge;

    cid = thisCI->id;

    InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, thisCI->id, ALL_END, ALL_EDGES,
			  GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
    while((edge = NextGraphEdgeIterator(&edges)) != NULL) {
      ChunkInstanceT *otherCI = GetGraphNode(ScaffoldGraph->ContigGraph,
					     (edge->idA == thisCI->id) ?
					     edge->idB : edge->idA);

      if ((edge->idA != thisCI->id) && otherCI->flags.bits.isUnique)
	// Only reset canonical edges so as not to repeat the work.
	continue;

      // Initialize the status flags for this edge
      InitGraphEdgeFlags(ScaffoldGraph->ContigGraph, edge);
    }
  }
  return;
}


void MarkEssentialEdges(ScaffoldGraphT *graph,
                        int markShakyBifurcations) {
  ResetEdgeStatus(graph);
  AddScaffoldInferredEdges(graph);
  MarkRedundantUniqueToUniqueEdges(graph);
  MarkTwoHopConfirmedEdges(graph);
  MarkPathRemovedEdges(graph);
  SmoothWithInferredEdges(graph, markShakyBifurcations);
}






void
BuildUniqueCIScaffolds(ScaffoldGraphT *graph,
                       int markShakyBifurcations,
                       int verbose) {
  int iteration = 0;
  int todo = TRUE;

  // on_on
  // this setting used for mouse_20010730 cgw run
  // markShakyBifurcations = FALSE;

  // off_on
  markShakyBifurcations = TRUE;  // ALH, 6/20/2004: why turn this on?

  // Mark essential edges, and look for shaky bifurcations
  // Iterate, until there are no shaky bifucations left
  //
  if (markShakyBifurcations) {
    for(iteration = 0; todo == TRUE ; iteration++) {
      int numShaky = 0;
      MarkEssentialEdges(graph,markShakyBifurcations);

      numShaky = MarkShakyBifurcations(graph);
      todo = (numShaky > 0);
      if (todo)
	DeleteInferredEdges(graph);
    }
  }
  DeleteInferredEdges(graph);
  MarkEssentialEdges(graph,FALSE);

  CreateScaffolds(graph);

  DeleteInferredEdges(graph);
}

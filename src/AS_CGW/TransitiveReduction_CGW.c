
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
static char *rcsid = "$Id: TransitiveReduction_CGW.c,v 1.22 2008-11-02 06:27:13 brianwalenz Exp $";

// This file contains the code for computing the candidate
// chunks of scaffolds.


//#define INSTRUMENT_CGW
//#define MEASURE_GRAPH_COMPLEXITY
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

#if 1
static int CompareCIEdgeTMeans (const void *c1, const void *c2){
  CIEdgeT *s1 = GetCIEdgeT(graphCIEdges, *(CDS_CID_t *)c1);
  CIEdgeT *s2 = GetCIEdgeT(graphCIEdges, *(CDS_CID_t *)c2);

  if(s1->distance.mean > s2->distance.mean)
    {
      return((int)1);
    }
  else
    {
      return((int)-1);
    }
  return 0;
}

#else

static int CompareCIEdgeTMeans (const void *c1, const void *c2)
{
  CIEdgeT *s1 = GetCIEdgeT(graphCIEdges, *(CDS_CID_t *)c1);
  CIEdgeT *s2 = GetCIEdgeT(graphCIEdges, *(CDS_CID_t *)c2);
  CDS_CID_t s1_contigA, s1_contigB, s2_contigA, s2_contigB;
  CDS_CID_t s1_unsharedContig, s2_unsharedContig;
  double s1_bpLength, s2_bpLength;

  s1_contigA = s1->idA;
  s1_contigB = s1->idB;
  s2_contigA = s2->idA;
  s2_contigB = s2->idB;
  if ( s1_contigA == s2_contigA )
    {
      s1_unsharedContig = s1_contigB;
      s2_unsharedContig = s2_contigB;
    }
  else if ( s1_contigA == s2_contigB )
    {
      s1_unsharedContig = s1_contigB;
      s2_unsharedContig = s2_contigA;
    }
  else if ( s1_contigB == s2_contigA )
    {
      s1_unsharedContig = s1_contigA;
      s2_unsharedContig = s2_contigB;
    }
  else if ( s1_contigB == s2_contigB )
    {
      s1_unsharedContig = s1_contigA;
      s2_unsharedContig = s2_contigA;
    }
  else
    assert(0);

  s1_bpLength = GetGraphNode(ScaffoldGraph->ContigGraph,
                             s1_unsharedContig)->bpLength.mean;
  s2_bpLength = GetGraphNode(ScaffoldGraph->ContigGraph,
                             s2_unsharedContig)->bpLength.mean;

  if (s1->distance.mean + sqrt(s1->distance.variance) + s1_bpLength >
      s2->distance.mean + sqrt(s2->distance.variance) + s2_bpLength )
    return((int)1);
  else
    return((int)-1);
  return 0;
}
#endif

void PrintLocalEdges(ScaffoldGraphT *graph, ChunkInstanceT *thisCI, int degree)
{
  CIEdgeT * edge;
  GraphEdgeIterator edges;

  fprintf(stdout, "Local Edges around CI " F_CID "\n", thisCI->id);
  fprintf(stdout, "\tA end:\n");
  InitGraphEdgeIterator(graph->RezGraph,
                        thisCI->id,
                        A_END,
                        ALL_EDGES,
                        GRAPH_EDGE_DEFAULT,
                        &edges);
  while((edge = NextGraphEdgeIterator(&edges))!= NULL)
    {
      if(!getEssentialEdgeStatus(edge))
        continue;
      fprintf(stdout, "\t" F_CID "\n",
              (thisCI->id == edge->idB) ? edge->idA : edge->idB);
    }

  fprintf(stdout, "\tB end:\n");
  InitGraphEdgeIterator(graph->RezGraph,
                        thisCI->id,
                        B_END,
                        ALL_EDGES,
                        GRAPH_EDGE_DEFAULT,
                        &edges);
  while((edge = NextGraphEdgeIterator(&edges))!= NULL)
    {
      if(!getEssentialEdgeStatus(edge))
        continue;
      fprintf(stdout, "\t" F_CID "\n",
              (thisCI->id == edge->idB) ? edge->idA : edge->idB);
    }

  if(degree < 2)
    {
      InitGraphEdgeIterator(graph->RezGraph,
                            thisCI->id,
                            A_END,
                            ALL_EDGES,
                            GRAPH_EDGE_DEFAULT,
                            &edges);
      while((edge = NextGraphEdgeIterator(&edges))!= NULL)
        {
          if(!getEssentialEdgeStatus(edge))
            continue;
          PrintLocalEdges(graph,
                          GetGraphNode(graph->RezGraph,
                                       (thisCI->id == edge->idA) ?
                                       edge->idB : edge->idA),
                          degree+1);
        }

      InitGraphEdgeIterator(graph->RezGraph,
                            thisCI->id,
                            B_END,
                            ALL_EDGES,
                            GRAPH_EDGE_DEFAULT,
                            &edges);
      while((edge = NextGraphEdgeIterator(&edges))!= NULL)
        {
          if(!getEssentialEdgeStatus(edge))
            continue;
          PrintLocalEdges(graph,
                          GetGraphNode(graph->RezGraph,
                                       (thisCI->id == edge->idA) ?
                                       edge->idB : edge->idA),
                          degree+1);
        }
    }
}


/****************************************************************************
 *  InferredEdgeOrientation
 *
 ****************************************************************************/
ChunkOrientationType InferredEdgeOrientation(ChunkOrientationType leftOrient,
                                             ChunkOrientationType rightOrient){

  if((leftOrient == AB_AB) || (leftOrient == BA_AB)){
    /* leftOrient == XX_AB */
    if((rightOrient == AB_AB) || (rightOrient == BA_AB)){
      /* rightOrient == XX_AB */
      return(AB_AB);
    }else{
      /* rightOrient == XX_BA */
      return(AB_BA);
    }
  }else{
    /* leftOrient == XX_BA */
    if((rightOrient == AB_AB) || (rightOrient == BA_AB)){
      /* rightOrient == XX_AB */
      return(BA_AB);
    }else{
      /* rightOrient == XX_BA */
      return(BA_BA);
    }
  }
}

/****************************************************************************
 *  TransitiveEdgeOrientation
 *
 ****************************************************************************/
ChunkOrientationType TransitiveEdgeOrientation(ChunkOrientationType leftOrient,
                                               ChunkOrientationType rightOrient){

  if((leftOrient == AB_AB) || (leftOrient == AB_BA)){
    /* leftOrient == AB_XX */
    if((rightOrient == AB_AB) || (rightOrient == BA_AB)){
      /* rightOrient == XX_AB */
      return(AB_AB);
    }else{
      /* rightOrient == XX_BA */
      return(AB_BA);
    }
  }else{
    /* leftOrient == BA_XX */
    if((rightOrient == AB_AB) || (rightOrient == BA_AB)){
      /* rightOrient == XX_AB */
      return(BA_AB);
    }else{
      /* rightOrient == XX_BA */
      return(BA_BA);
    }
  }
}

/****************************************************************************
 *  FoundTransitiveEdgePath
 *
 ****************************************************************************/
#define AS_CGW_MAX_FTEP_RECURSE_DEPTH 15

int FoundTransitiveEdgePath(ScaffoldGraphT *graph,
			    CIEdgeT *pathEdge,
			    CIEdgeT *edge,  /* edge we're trying to transitively remove */
			    ChunkInstanceT *startCI, /* Start of path */
			    ChunkInstanceT *endCI,   /* End of Path */
			    ChunkInstanceT *thisCI, /* Where we are now */
			    int verbose, int recurseDepth){

  CIEdgeT localPathEdge;
  CIEdgeT *transEdge;
  GraphEdgeIterator transEdges;
  float chiSquaredValue;
  ChunkOrientationType edgeOrient, pathEdgeOrient;
  int returnVal = FALSE;


  assert(!isSloppyEdge(edge));
  if(++recurseDepth > AS_CGW_MAX_FTEP_RECURSE_DEPTH){
    if(verbose)
      fprintf(stderr,"* Max recursion Returning false\n");
    return(FALSE);
  }
  if(verbose)
    fprintf(stderr,"* FoundTransEdgepath path (" F_CID "," F_CID ") edge (" F_CID "," F_CID ") start:" F_CID " end:" F_CID " this:" F_CID " depth %d\n",
	    pathEdge->idA, pathEdge->idB,
	    edge->idA, edge->idB,
	    startCI->id, endCI->id, thisCI->id,recurseDepth);

  edgeOrient = GetEdgeOrientationWRT(edge, startCI->id);
  pathEdgeOrient = GetEdgeOrientationWRT(pathEdge, startCI->id);
  if(thisCI == endCI){  // we got there
    if(edgeOrient == pathEdgeOrient){ // correct orientation
      if(PairwiseChiSquare((float)pathEdge->distance.mean,
			   (float)pathEdge->distance.variance,
			   (float)edge->distance.mean,
			   (float)edge->distance.variance, (LengthT *)NULL,
			   &chiSquaredValue, (float)PAIRWISECHI2THRESHOLD_CGW)){ // passes Chi-square
	if(verbose)
	  fprintf(stderr,"* Returning true\n");
	return(TRUE);
      }else{
	if(verbose){
	  fprintf(GlobalData->stderrc, "Path Failed PairwiseChiSquare cid:" F_CID "-" F_CID "\nReal(%d %f) Implied(%d %f) ChiSquared %f\n",
		  startCI->id, endCI->id,
		  (int)edge->distance.mean, sqrt(edge->distance.variance),
		  (int)pathEdge->distance.mean,
		  sqrt(pathEdge->distance.variance),
		  chiSquaredValue);
	}
	if(verbose)
	  fprintf(stderr,"* Returning false1\n");
	return(FALSE);
      }
    }else{
      if(verbose){
	fprintf(GlobalData->stderrc, "Path Failed Orientation cid:" F_CID "-" F_CID " Orientation:%c-%c\n",
		startCI->id, endCI->id, (char)edgeOrient, (char)pathEdgeOrient);
	if(verbose)
	  fprintf(stderr,"* Returning FALSE2\n");
	return(FALSE);
      }
    }
  }
  if(pathEdge->distance.mean > edge->distance.mean){
    if(!PairwiseChiSquare((float)pathEdge->distance.mean,
			  (float)pathEdge->distance.variance,
			  (float)edge->distance.mean,
			  (float)edge->distance.variance, (LengthT *)NULL,
			  &chiSquaredValue, (float)PAIRWISECHI2THRESHOLD_CGW)){
      if(verbose)
	fprintf(stderr,"* Returning FALSE3\n");
      return(FALSE);
    }
  }
  if((pathEdgeOrient == AB_AB) ||
     (pathEdgeOrient == BA_AB)){
    /* transEdge must be a B_END edge */
    InitGraphEdgeIterator(graph->RezGraph, thisCI->id, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &transEdges);// Use merged edges
  }else{
    /* transEdge must be a A_END edge */
    InitGraphEdgeIterator(graph->RezGraph, thisCI->id, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &transEdges);// Use merged edges
  }
  while((transEdge = NextGraphEdgeIterator(&transEdges))!= NULL){
    ChunkInstanceT *nextCI;
    if(!transEdge->flags.bits.isActive ||
       !transEdge->flags.bits.isUniquetoUnique){
      // Do not use inactive edges or edges to nonUnique CIs.
      continue;
    }
    // SAK HACK
    /* Do not use overlap edges in path walks */
    if(isSingletonOverlapEdge(transEdge))
      continue;

    assert(!isSloppyEdge(transEdge));

    nextCI = GetGraphNode(graph->RezGraph,
			  ((transEdge->idA == thisCI->id) ?
			   transEdge->idB : transEdge->idA));
    localPathEdge = *pathEdge;
    localPathEdge.distance.mean += thisCI->bpLength.mean +
      transEdge->distance.mean;
    localPathEdge.distance.variance += thisCI->bpLength.variance +
      transEdge->distance.variance;
    localPathEdge.idB = (transEdge->idA == thisCI->id) ?
      transEdge->idB : transEdge->idA;
    localPathEdge.orient = TransitiveEdgeOrientation(pathEdgeOrient,
						     GetEdgeOrientationWRT(transEdge, thisCI->id));
    if(verbose)
      fprintf(stderr,"* Recurse\n");
    if(FoundTransitiveEdgePath(graph, &localPathEdge, edge, startCI, endCI,
			       nextCI, verbose, recurseDepth)){
      if(verbose)
        fprintf(stderr,"* Set isConfirmed and returnVal to TRUE\n");
      transEdge->flags.bits.isConfirmed = TRUE;
      if(verbose)
	PrintGraphEdge(GlobalData->stderrc, graph->RezGraph, "\t", transEdge,
                       thisCI->id);
      returnVal = TRUE;
      /* For speed we are going to return when we find the first path
	 and not find all paths so as to possibly confirm some extra
	 edges.
      */
      return(returnVal);
    }
  }
  if(verbose)
    fprintf(stderr,"* Returning %s\n", (returnVal?"TRUEF":"FALSEF"));
  return(returnVal);
}

/****************************************************************************
 *  MarkPathRemovedEdgesOneEnd
 *
 ****************************************************************************/
void MarkPathRemovedEdgesOneEnd(ScaffoldGraphT *graph, ChunkInstanceT *thisCI,
				int end, int verbose){
  GraphEdgeIterator edges;
  CIEdgeT *edge;
  InitGraphEdgeIterator(graph->RezGraph, thisCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
  while((edge = NextGraphEdgeIterator(&edges))!= NULL){
    GraphEdgeIterator transEdges;
    CIEdgeT pathEdge;
    CIEdgeT *transEdge;
    ChunkInstanceT *endCI;
    if(edge->idA != thisCI->id){
      // Only check an edge once in canonical direction.
      continue;
    }
    if(!edge->flags.bits.isActive || !edge->flags.bits.isUniquetoUnique){
      // No need to check inactive edges - only interested in edges between
      // Unique CIs.
      continue;
    }
    assert(!isSloppyEdge(edge));
    endCI = GetGraphNode(graph->RezGraph,
                         ((edge->idA == thisCI->id) ? edge->idB : edge->idA));
    // Only use edges on the same end of thisCI to transitively remove an edge.
    InitGraphEdgeIterator(graph->RezGraph, thisCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &transEdges);// Use merged edges
    while((transEdge = NextGraphEdgeIterator(&transEdges))!= NULL){
      ChunkInstanceT *nextCI;
      CDS_CID_t nextCID;
      int recurseDepth;
      if(edge == transEdge){
	// Do not remove an edge with itself.
	continue;
      }
      if(!transEdge->flags.bits.isActive ||
	 !transEdge->flags.bits.isUniquetoUnique){
	// Do not use inactive edges - only interested in edges between
	// Unique CIs.
	continue;
      }
      assert(!isSloppyEdge(transEdge));
      if(verbose)
        fprintf(stderr,"* Trying to transitively remove (" F_CID "," F_CID ",%c) with (" F_CID "," F_CID ",%c)\n",
                edge->idA, edge->idB, GetEdgeOrientationWRT(edge, thisCI->id),
                transEdge->idA, transEdge->idB, GetEdgeOrientationWRT(transEdge, thisCI->id));

      nextCID = (transEdge->idA == thisCI->id) ? transEdge->idB :
        transEdge->idA;
      nextCI = GetGraphNode(graph->RezGraph, nextCID);
      pathEdge = *transEdge;
      pathEdge.idA = thisCI->id;
      pathEdge.idB = nextCID;
      pathEdge.orient = GetEdgeOrientationWRT(transEdge, thisCI->id);
      recurseDepth = 0;
      if(FoundTransitiveEdgePath(graph, &pathEdge, edge, thisCI, endCI,
				 nextCI, verbose, recurseDepth)){
	edge->flags.bits.isTransitivelyRemoved = TRUE;
	edge->flags.bits.isConfirmed = TRUE;
	transEdge->flags.bits.isConfirmed = TRUE;
	if(verbose){
	  PrintGraphEdge(GlobalData->stderrc, graph->RezGraph, "\t", transEdge,
                         thisCI->id);
	  PrintGraphEdge(GlobalData->stderrc, graph->RezGraph, "TR\t", edge,
                         thisCI->id);
	}
      }
    }
  }
  return;
}

/*****************************************************************************
 *  MarkPathRemovedEdges
 *
 ****************************************************************************/
void MarkPathRemovedEdges(ScaffoldGraphT *graph, int verbose){
  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;
  int cnt = 1;
  InitGraphNodeIterator(&nodes, graph->RezGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL){

    if(verbose){
      if((cnt++ % 1000) == 0){
        fprintf(stderr,"* %d MarkPathRemoved for node " F_CID "\n",
                cnt, thisCI->id);
        fflush(stderr);
      }
    }
    // Do the A_END edges.
    MarkPathRemovedEdgesOneEnd(graph, thisCI, A_END, verbose);
    // Do the B_END edges.
    MarkPathRemovedEdgesOneEnd(graph, thisCI, B_END, verbose);
  }
  return;
}

/****************************************************************************
 *  MarkTwoHopConfirmedEdgesOneEnd
 *
 ****************************************************************************/
void MarkTwoHopConfirmedEdgesOneEnd(ScaffoldGraphT *graph,
				    ChunkInstanceT *thisCI,
				    int end, int verbose){
  GraphEdgeIterator edges;
  CIEdgeT *edge;
  InitGraphEdgeIterator(graph->RezGraph, thisCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges

  while((edge = NextGraphEdgeIterator(&edges))!= NULL){
    GraphEdgeIterator hop2Edges;
    CIEdgeT pathEdge;
    CIEdgeT *hop2Edge;
    CDS_CID_t middleCID;
    ChunkInstanceT *middleCI;
    if(!edge->flags.bits.isActive || !edge->flags.bits.isUniquetoUnique){
      // No need to check inactive edges - only interested in edges between
      // Unique CIs.
      continue;
    }
    /* Do not use overlap edges in two hops */
    if(isSingletonOverlapEdge(edge))
      continue;

    middleCID = (edge->idA == thisCI->id) ? edge->idB : edge->idA;
    middleCI = GetGraphNode(graph->RezGraph, middleCID);
    pathEdge = *edge;
    pathEdge.idA = thisCI->id;
    pathEdge.idB = middleCID;
    pathEdge.orient = GetEdgeOrientationWRT(edge, thisCI->id);
    pathEdge.distance.mean += middleCI->bpLength.mean;
    pathEdge.distance.variance += middleCI->bpLength.variance;
    if((pathEdge.orient == AB_AB) ||
       (pathEdge.orient == BA_AB)){
      /* hop2Edge must be a B_END edge */
      InitGraphEdgeIterator(graph->RezGraph, middleCID, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &hop2Edges);// Use merged edges
    }else{
      /* hop2Edge must be a A_END edge */
      InitGraphEdgeIterator(graph->RezGraph, middleCID, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &hop2Edges);// Use merged edges
    }
    while((hop2Edge = NextGraphEdgeIterator(&hop2Edges))!= NULL){
      ChunkInstanceT *endCI;
      ChunkOrientationType combinedOrient;
      double combinedMean, combinedVariance;
      GraphEdgeIterator confirmEdges;
      CIEdgeT *confirmEdge;
      if(!hop2Edge->flags.bits.isActive ||
	 !hop2Edge->flags.bits.isUniquetoUnique){
	// Do not use inactive edges - only interested in edges between
	// Unique CIs.
	continue;
      }
      /* Do not use overlap edges in two hops */
      if(isSingletonOverlapEdge(hop2Edge))
	continue;

      endCI = GetGraphNode(graph->RezGraph,
                           ((hop2Edge->idA == middleCID) ?
                            hop2Edge->idB : hop2Edge->idA));
      combinedOrient = TransitiveEdgeOrientation(pathEdge.orient,
                                                 GetEdgeOrientationWRT(hop2Edge, middleCID));
      combinedMean = pathEdge.distance.mean + hop2Edge->distance.mean;
      combinedVariance = pathEdge.distance.variance +
        hop2Edge->distance.variance;
      // Only use edges on the same end of thisCI to confirm a two hop edge.
      // Use a copy of the original iterator so that we do not find the
      // same pair of two hop paths twice.
      confirmEdges = edges;
      while((confirmEdge = NextGraphEdgeIterator(&confirmEdges))!= NULL){
	GraphEdgeIterator confirmHop2Edges;
	CIEdgeT confirmPathEdge;
	CIEdgeT *confirmHop2Edge;
	CDS_CID_t confirmMiddleCID;
	ChunkInstanceT *confirmMiddleCI;
	if(!confirmEdge->flags.bits.isActive ||
	   !confirmEdge->flags.bits.isUniquetoUnique){
	  // Do not use inactive edges - only interested in edges between
	  // Unique CIs.
	  continue;
	}
	/* Do not use overlap edges in two hops */
	if(isSingletonOverlapEdge(confirmEdge))
	  continue;
	confirmMiddleCID = (confirmEdge->idA == thisCI->id) ?
          confirmEdge->idB : confirmEdge->idA;
	confirmMiddleCI = GetGraphNode(graph->RezGraph,
                                       confirmMiddleCID);
	confirmPathEdge = *confirmEdge;
	confirmPathEdge.idA = thisCI->id;
	confirmPathEdge.idB = confirmMiddleCID;
	confirmPathEdge.orient = GetEdgeOrientationWRT(confirmEdge,
                                                       thisCI->id);
	confirmPathEdge.distance.mean += confirmMiddleCI->bpLength.mean;
	confirmPathEdge.distance.variance += confirmMiddleCI->bpLength.variance;
	if((confirmPathEdge.orient == AB_AB) ||
	   (confirmPathEdge.orient == BA_AB)){
	  /* hop2Edge must be a B_END edge */
	  InitGraphEdgeIterator(graph->RezGraph, confirmMiddleCID, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &confirmHop2Edges);// Use merged edges
	}else{
	  /* hop2Edge must be a A_END edge */
	  InitGraphEdgeIterator(graph->RezGraph, confirmMiddleCID, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &confirmHop2Edges);// Use merged edges
	}
	while((confirmHop2Edge = NextGraphEdgeIterator(&confirmHop2Edges))!= NULL){
	  ChunkInstanceT *confirmEndCI;
	  ChunkOrientationType confirmOrient;
	  double confirmMean, confirmVariance;
	  float chiSquaredValue;
	  if(!confirmHop2Edge->flags.bits.isActive ||
	     !confirmHop2Edge->flags.bits.isUniquetoUnique){
	    // Do not use inactive edges - only interested in edges between
	    // Unique CIs.
	    continue;
	  }
	  /* Do not use overlap edges in two hops */
	  if(isSingletonOverlapEdge(confirmHop2Edge))
	    continue;

	  confirmEndCI = GetGraphNode(graph->RezGraph,
                                      ((confirmHop2Edge->idA == confirmMiddleCID) ?
                                       confirmHop2Edge->idB : confirmHop2Edge->idA));
	  confirmOrient = TransitiveEdgeOrientation(confirmPathEdge.orient,
                                                    GetEdgeOrientationWRT(confirmHop2Edge,
                                                                          confirmMiddleCID));
	  confirmMean = confirmPathEdge.distance.mean +
            confirmHop2Edge->distance.mean;
	  confirmVariance = confirmPathEdge.distance.variance +
            confirmHop2Edge->distance.variance;
	  if((endCI == confirmEndCI) && (combinedOrient == confirmOrient) &&
	     PairwiseChiSquare((float)combinedMean, (float)combinedVariance,
			       (float)confirmMean, (float)confirmVariance,
			       (LengthT *)NULL, &chiSquaredValue,
			       (float)PAIRWISECHI2THRESHOLD_CGW)){
	    edge->flags.bits.isConfirmed = TRUE;
	    hop2Edge->flags.bits.isConfirmed = TRUE;
	    confirmEdge->flags.bits.isConfirmed = TRUE;
	    confirmHop2Edge->flags.bits.isConfirmed = TRUE;
	    if(verbose){
	      PrintGraphEdge(GlobalData->stderrc, graph->RezGraph, "2Hop\t1 ", edge,
			     thisCI->id);
	      PrintGraphEdge(GlobalData->stderrc, graph->RezGraph, "\t2 ", hop2Edge,
			     middleCID);
	      PrintGraphEdge(GlobalData->stderrc, graph->RezGraph, "\tC ", confirmEdge,
			     thisCI->id);
	      PrintGraphEdge(GlobalData->stderrc, graph->RezGraph, "\t2 ", confirmHop2Edge,
			     confirmMiddleCID);
	    }
	  }
	}
      }
    }
  }
  return;
}

/****************************************************************************
 *  MarkTwoHopConfirmedEdges
 *
 ****************************************************************************/
void MarkTwoHopConfirmedEdges(ScaffoldGraphT *graph, int verbose){

  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;

  InitGraphNodeIterator(&nodes, graph->RezGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL){

    if(verbose){
      if((thisCI->id % 10000) == 0){
        fprintf(stderr,"* MarkTwoHopConfirmedEdges for node " F_CID "\n",
                thisCI->id);
        fflush(stderr);
      }
    }

    // Do the A_END edges.
    MarkTwoHopConfirmedEdgesOneEnd(graph, thisCI, A_END, verbose);
    // Do the B_END edges.
    MarkTwoHopConfirmedEdgesOneEnd(graph, thisCI, B_END, verbose);
  }
  return;
}

/****************************************************************************
 *  CompareBestCIEdgeT
 *
 ****************************************************************************/
static int CompareBestCIEdgeT(const void *c1, const void *c2){
  CIEdgeT *s1 = (CIEdgeT *)c1;
  CIEdgeT *s2 = (CIEdgeT *)c2;

  if(s1->flags.bits.isEssential != s2->flags.bits.isEssential){
    return s1->flags.bits.isEssential ? 1 : -1;
  }
  if(s1->flags.bits.isConfirmed != s2->flags.bits.isConfirmed){
    return s1->flags.bits.isConfirmed ? 1 : -1;
  }
  if(s1->flags.bits.isLeastSquares != s2->flags.bits.isLeastSquares){
    return s1->flags.bits.isLeastSquares ? 1 : -1;
  }
  if(s1->flags.bits.isActive != s2->flags.bits.isActive){
    return s1->flags.bits.isActive ? 1 : -1;
  }
  if(s1->flags.bits.isProbablyBogus != s2->flags.bits.isProbablyBogus){
    return s1->flags.bits.isProbablyBogus ? -1 : 1;
  }
  if(s1->edgesContributing < s2->edgesContributing){
    return -1;
  }
  if(s1->edgesContributing > s2->edgesContributing){
    return 1;
  }
  if(isOverlapEdge(s1) && !isOverlapEdge(s2)){
    return -1;
  }
  if(isOverlapEdge(s2) && !isOverlapEdge(s1)){
    return 1;
  }
  if(s1->distance.variance > s2->distance.variance){
    return -1;
  }
  if(s1->distance.variance < s2->distance.variance){
    return 1;
  }
  if(s1->distance.mean > s2->distance.mean){
    return -1;
  }
  if(s1->distance.mean < s2->distance.mean){
    return 1;
  }
  if(s1->minDistance > s2->minDistance){
    return -1;
  }
  if(s1->minDistance < s2->minDistance){
    return 1;
  }
#if 0
  if(s1->maxDistance > s2->maxDistance){
    return -1;
  }
  if(s1->maxDistance < s2->maxDistance){
    return 1;
  }
#endif
  return 0;
}

/****************************************************************************
 *  MarkRedundantUniqueToUniqueEdges
 *
 ****************************************************************************/
void MarkRedundantUniqueToUniqueEdges(ScaffoldGraphT *graph, int verbose){

  CDS_CID_t cid;
  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;

  InitGraphNodeIterator(&nodes, graph->RezGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL){
    GraphEdgeIterator edges;
    CIEdgeT *edge;

    if(verbose)
      fprintf(stderr,"*MarkRedundant on node " F_CID "\n",   thisCI->id);

    cid = thisCI->id;

    InitGraphEdgeIterator(graph->RezGraph, thisCI->id , ALL_END,	ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
    edge = NextGraphEdgeIterator(&edges);
    while(edge != NULL){
      CIEdgeT *bestEdge;

      if((edge->idA != thisCI->id) ||
	 !edge->flags.bits.isActive || !edge->flags.bits.isUniquetoUnique){
	// Only check canonical edges so as not to repeat the work.
	// No need to check inactive edges - only interested in edges between
	// Unique CIs.
	if(verbose && ((edge->idA == thisCI->id) ||
		       !edge->flags.bits.isUniquetoUnique)){
	  PrintGraphEdge(GlobalData->stderrc, graph->RezGraph,
                         edge->flags.bits.isUniquetoUnique ? "\tIA " : "\tNU ",
                         edge, thisCI->id);
	}
	edge = NextGraphEdgeIterator(&edges);
	continue;
      }
      bestEdge = edge;
      // We make the assumption that the edges are sorted such that the edges
      // between a pair of CIs are contiguous.
      while(((edge = NextGraphEdgeIterator(&edges)) != NULL) &&
	    (edge->idB == bestEdge->idB)){
	if(!edge->flags.bits.isActive){
	  if(verbose){
	    PrintGraphEdge(GlobalData->stderrc, graph->RezGraph,"\tIA ", edge, thisCI->id);
	  }
	  // No need to check inactive edges.
	  continue;
	}
	if(CompareBestCIEdgeT(bestEdge, edge) < 0){
	  bestEdge->flags.bits.isActive = FALSE;
	  bestEdge->flags.bits.isRedundantRemoved = TRUE;
	  if(verbose){
	    PrintGraphEdge(GlobalData->stderrc, graph->RezGraph,"\tR ", bestEdge, thisCI->id);
	  }
	  bestEdge = edge;
	}else{
	  edge->flags.bits.isActive = FALSE;
	  edge->flags.bits.isRedundantRemoved = TRUE;
	  if(verbose){
	    PrintGraphEdge(GlobalData->stderrc, graph->RezGraph, "\tR ", edge, thisCI->id);
	  }
	}
      }
      if(verbose){
	PrintGraphEdge(GlobalData->stderrc, graph->RezGraph, "A ", bestEdge, thisCI->id);
      }
    }
  }
  return;
}

EdgeCGW_T *FindEdgeBetweenCIsChiSquare(GraphCGW_T *graph,
                                       ChunkInstanceT *sourceCI, CDS_CID_t targetId,
                                       ChunkOrientationType edgeOrient,
                                       double inferredMean, double inferredVariance,
                                       float *returnChiSquaredValue,
                                       float chiSquareThreshold, int *isEssential,
                                       int verbose){
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
  while((edge = NextGraphEdgeIterator(&edges))!= NULL){
    CDS_CID_t otherCID = (edge->idA == sourceCI->id) ? edge->idB : edge->idA;
    // We aren't interested in high variance (sloppy) edges, since they will pass the chi square too easily
    if(isSloppyEdge(edge))
      continue;
    if(otherCID == targetId){
      // If there is an overlap edge with the 'right' orientation, remember it
      if(isOverlapEdge(edge) &&
	 (GetEdgeOrientationWRT(edge, sourceCI->id) == edgeOrient)){
	overlapEdgeExists = TRUE;
      }
      // Is this the best edge irrespective of orientation?
      if((bestEdge == (CIEdgeT *)NULL) ||
	 (CompareBestCIEdgeT((const void *)edge, (const void *)bestEdge) > 0)){
	bestEdge = edge;
      }
    }
  }
  // Were any edges found irrespective of orientation between source,target?
  if(bestEdge != (CIEdgeT *)NULL){
    if(verbose){
      PrintGraphEdge(GlobalData->stderrc, graph, "FindChi ", bestEdge,
		     sourceCI->id);
    }
    // This is a return value that, if true, causes the caller NOT to add an inferred
    // edge when this function returns NULL
    *isEssential = getEssentialEdgeStatus(bestEdge);

    // If the orientation was 'wrong'
    if(GetEdgeOrientationWRT(bestEdge, sourceCI->id) != edgeOrient){
      if(verbose){
	fprintf(GlobalData->stderrc, "RecSub(%d) Orientation:%c,%c SoId:" F_CID " TaId:" F_CID "\n",
		bestEdge->flags.bits.isEssential,
		(char)GetEdgeOrientationWRT(bestEdge, sourceCI->id),
		(char)edgeOrient, sourceCI->id, targetId);
      }
      // and the edge is essential, return NULL
      if(getEssentialEdgeStatus(bestEdge) == TRUE){
	return((CIEdgeT *)NULL);
      }
    }else // If the orientation is 'right', check for Chi Square compatibility
#if 0
      if ( inferredVariance > FLT_MAX ) inferredVariance = FLT_MAX;  // hack instituted on 5/15/01 to help mouse_20010508
#endif
    if(PairwiseChiSquare((float)inferredMean, (float)inferredVariance,
                         (float)bestEdge->distance.mean,
                         (float)bestEdge->distance.variance,
                         (LengthT *)NULL, &chiSquaredValue,
                         (float)chiSquareThreshold)){
      if(verbose){
	fprintf(GlobalData->stderrc, "RecSub(%d) PassedChi2:%f SoId:" F_CID " TaId:" F_CID "\n",
		bestEdge->flags.bits.isEssential,
		chiSquaredValue, sourceCI->id, targetId);
      }
      *returnChiSquaredValue = chiSquaredValue;
      return(bestEdge);
    }else{ // If right orientation, but not chi square compatible
      if(verbose){
	fprintf(GlobalData->stderrc, "RecSub(%d) ChiSquare:%f,%f,%f,%f,%f SoId:" F_CID " TaId:" F_CID "\n",
		bestEdge->flags.bits.isEssential, chiSquaredValue,
		inferredMean, sqrt(inferredVariance),
		bestEdge->distance.mean,
		sqrt(bestEdge->distance.variance),
		sourceCI->id, targetId);
      }
      // If marked essential, return NULL
      if(getEssentialEdgeStatus(bestEdge) == TRUE){
	return((CIEdgeT *)NULL);
      }
    }
  }


  // If there isn't an overlap edge in the right orientation, look for one
  if(!overlapEdgeExists){
    CDS_COORD_t minOverlap, maxOverlap;
    CDS_COORD_t minVariance;
    CDS_COORD_t delta;

    if(inferredVariance < 1.0)
      minVariance = 1;
    else
      minVariance = (CDS_COORD_t) sqrt(inferredVariance);

    delta = 3 * minVariance;
    minOverlap = MAX(CGW_MISSED_OVERLAP,
		     -(inferredMean + delta));
    maxOverlap = -(inferredMean - delta);
    if(maxOverlap >= CGW_MISSED_OVERLAP){
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
      if(olap.suspicious){
        fprintf(stderr,"* TR: SUSPICIOUS Overlap found! Looked for ("
                F_CID "," F_CID ",%c)[" F_COORD "," F_COORD "] found ("
                F_CID "," F_CID ",%c) " F_COORD " (would have been "
                F_COORD ")\n",
                sourceCI->id, targetId, edgeOrient,
                minOverlap, maxOverlap,
                olap.spec.cidA, olap.spec.cidB, olap.spec.orientation,
                olap.overlap,olap.overlap+olap.ahg+olap.bhg);
      }else if(olap.overlap){
	if(verbose){
	  fprintf(GlobalData->stderrc, "RecSub Transitive reduction inserting edge (" F_CID "," F_CID ",%c) " F_COORD " (searched on range " F_COORD "," F_COORD ", inf:%g del:" F_COORD ")\n",
		  olap.spec.cidA, olap.spec.cidB, olap.spec.orientation,
		  olap.overlap, minOverlap, maxOverlap, inferredMean, delta);
	}
      }else if(bestEdge == (CIEdgeT *)NULL){ // If we didn't find an overlap, and there was no other edge return NULL
	return((CIEdgeT *)NULL);
      }
    }else if(bestEdge == (CIEdgeT *)NULL){ // There was no potential overlap, and there was no other edge return NULL
      return((CIEdgeT *)NULL);
    }
  }

  // If the 'best' edge didn't pass the above tests (orientation and chi square), we now look for an edge with proper orientation
  // that passes chi square
  // It is important to set bestEdge to NULL here so that any edge found and
  // returned by this function points to the current edge array which might
  // have been reallocated when an overlap edge was added in the code above.

  bestEdge = (CIEdgeT *)NULL;
  if((edgeOrient == AB_AB) || (edgeOrient == AB_BA)){
    /* edgeOrient == AB_XX */
    end = B_END;
  }else{
    /* edgeOrient == BA_XX */
    end = A_END;
  }
  InitGraphEdgeIterator(graph, sourceCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT,  &edges);// Use merged edges
  while((edge = NextGraphEdgeIterator(&edges))!= NULL){
    CDS_CID_t otherCID = (edge->idA == sourceCI->id) ? edge->idB : edge->idA;

    // We aren't interested in high variance (sloppy) edges, since they will pass the chi square too easily
    if(isSloppyEdge(edge))
      continue;

    if((otherCID == targetId) &&
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
			 (float)chiSquareThreshold)){
      if(bestEdge == (CIEdgeT *)NULL ||
	 (chiSquaredValue < bestChiSquaredValue)){
	bestEdge = edge;
	bestChiSquaredValue = chiSquaredValue;
	*isEssential = getEssentialEdgeStatus(bestEdge);
	if(verbose){
	  fprintf(GlobalData->stderrc, "RecSub(%d) PassedChi2Loop:%f SoId:" F_CID " TaId:" F_CID "\n",
		  bestEdge->flags.bits.isEssential,
		  chiSquaredValue, sourceCI->id, targetId);
	}
      }
    }
  }
  *returnChiSquaredValue = bestChiSquaredValue;
  return(bestEdge);
}

int RecursiveSmoothWithInferredEdges(ScaffoldGraphT *graph,
				     ChunkInstanceT *thisCI,
				     int end, int verbose,
                                     int * numInferredAdded,
                                     CDS_CID_t firstID,
                                     int instFirstEnd,
                                     CDS_CID_t *lastID,
                                     ScaffoldInstrumenter * smoothed_si,
                                     FILE * instSmoothSuccessfp,
                                     FILE * instSmoothFailurefp){
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
  if(end == A_END){
    numBranch = thisCI->numEssentialA;
  }else{
    numBranch = thisCI->numEssentialB;
  }
  if(verbose){
    fprintf(GlobalData->stderrc, "RecSub(%d," F_CID ") numBranch:%d cId:" F_CID "(%d{" F_CID "},%d{" F_CID "})\n",
	    thisCI->flags.bits.smoothSeenAlready,
	    thisCI->smoothExpectedCID, numBranch, thisCI->id,
	    thisCI->numEssentialA, thisCI->essentialEdgeA,
	    thisCI->numEssentialB, thisCI->essentialEdgeB);
  }

  /*
    if there is one essential edge off end of thisCI, simple resolution:
    1. get the next CI (sourceCI) and the edge to it (sourceEdge)
    2. set sourceCI flags so it can be smoothed later
    3. set the essentialX of thisCI & sourceCI to the index of sourceEdge
    4. for instrumenting, if lastID isn't set yet, set it to sourceCI->id
  */
  if(numBranch < 2){
    CDS_CID_t sourceEdgeIndex;
    sourceCI = (ChunkInstanceT *)NULL;
    if(verbose){
      fprintf(GlobalData->stderrc, "RecSub numBranch<2:%d cId:" F_CID "\n",
	      numBranch, thisCI->id);
    }
    assert(numBranch == 1);
    InitGraphEdgeIterator(graph->RezGraph, thisCI->id, end,ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
    while((edge = NextGraphEdgeIterator(&edges))!= NULL){
      if(!getEssentialEdgeStatus(edge)){
	// Skip nonessential edges.
	continue;
      }
      assert(sourceCI == (ChunkInstanceT *)NULL);
      sourceCI = GetGraphNode(graph->RezGraph,
                              ((edge->idA == thisCI->id) ? edge->idB : edge->idA));
      sourceEdge = edge;
    }
    assert(sourceCI != (ChunkInstanceT *)NULL);

    // instrument the scaffold to see if it's acceptable
    // before modifying anything
    *lastID = (*lastID == NULLINDEX) ? sourceCI->id : *lastID;
    if(verbose)
      fprintf(stderr, F_CID " ", sourceCI->id);

#ifdef INSTRUMENT_TRANS_REDUCED
    {
      double badMates;
      double allMates;

      InstrumentContigPath(graph, smoothed_si,
                           firstID, instFirstEnd, *lastID);
      badMates = GetMateStatsBad(&(smoothed_si->mates.inter));
      allMates = badMates + GetMateStatsHappy(&(smoothed_si->mates.inter));

      if((int) allMates > 0 && badMates / allMates > .05)
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
    sourceEdgeIndex = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->RezGraph->edges, sourceEdge);
    if(end == B_END){
      thisCI->essentialEdgeB = sourceEdgeIndex;
    }else{
      thisCI->essentialEdgeA = sourceEdgeIndex;
    }
    sourceEdgeOrient = GetEdgeOrientationWRT(sourceEdge, thisCI->id);
    if((sourceEdgeOrient == AB_AB) || (sourceEdgeOrient == BA_AB)){
      sourceCI->essentialEdgeA = sourceEdgeIndex;
    }else{
      sourceCI->essentialEdgeB = sourceEdgeIndex;
    }

    return(TRUE);
  }

  { // Allocate memory
#ifdef ALLOC_MEMORY_ON_STACK
    CDS_CID_t branchEdges[numBranch];           // indices of essential edges off end of thisCI
    CDS_CID_t existingEdges[numBranch -1];      // indices of best chi2 edges between thisCI and each CI off end
    CDS_CID_t addedInferred[numBranch -1];      // indices of inferred edges added
    Target_Info_t branchTargets[numBranch]; //

#else
    branchEdges = (CDS_CID_t *)safe_malloc(numBranch * sizeof(*branchEdges));
    AssertPtr(branchEdges);
    existingEdges = (CDS_CID_t *)safe_malloc((numBranch - 1) * sizeof(*existingEdges));
    AssertPtr(existingEdges);
    addedInferred = (CDS_CID_t *)safe_malloc((numBranch - 1) * sizeof(*addedInferred));
    AssertPtr(addedInferred);
    branchTargets = (Target_Info_t *)safe_malloc(numBranch * sizeof(*branchTargets));
    AssertPtr(branchTargets);
#endif
    InitGraphEdgeIterator(graph->RezGraph, thisCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
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
    while((edge = NextGraphEdgeIterator(&edges))!= NULL){
      ChunkInstanceT *targetCI;
      if(verbose){
        PrintGraphEdge(GlobalData->stderrc, graph->RezGraph, "Tgt ", edge, thisCI->id);
      }
      if(!getEssentialEdgeStatus(edge)){
        // Skip nonessential edges.
        continue;
      }
      targetCI = GetGraphNode(graph->RezGraph,
                              ((edge->idA == thisCI->id) ? edge->idB : edge->idA));
      if(targetCI->flags.bits.smoothSeenAlready ||
         ((targetCI->smoothExpectedCID != NULLINDEX) &&
          (targetCI->smoothExpectedCID != thisCI->id))){

#ifdef INSTRUMENT_SMOOTHED
        // instrument the failed 'scaffold' before returning failure
        InstrumentContigPath(graph, smoothed_si,
                             firstID, instFirstEnd, thisCI->id);
        PrintScaffoldInstrumenter(graph, smoothed_si,
                                  InstrumenterVerbose2,
                                  "", instSmoothFailurefp);
#endif

        if(verbose){
          fprintf(GlobalData->stderrc, "RecSub Loop(%d," F_CID ") cId:" F_CID " TcId:" F_CID "\n",
                  targetCI->flags.bits.smoothSeenAlready,
                  targetCI->smoothExpectedCID, thisCI->id, targetCI->id);
        }
#ifndef ALLOC_MEMORY_ON_STACK
        safe_free(branchEdges);
        safe_free(existingEdges);
        safe_free(addedInferred);
        safe_free(branchTargets);
#endif
        return(FALSE);
      }
      // We need to save indices instead of pointers because the Variable
      // Array graph->CIEdges can get reallocated when we insert inferred
      // edges.
      *branchEnd++ = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->RezGraph->edges, edge);
      branchTargetsEnd->id = targetCI->id;
      branchTargetsEnd->lo = (edge->distance.mean -
                              (3 * sqrt(edge->distance.variance))) +
        (targetCI->bpLength.mean - (3 * sqrt(targetCI->bpLength.variance)));
      branchTargetsEnd->hi = (edge->distance.mean +
                              (3 * sqrt(edge->distance.variance))) +
        (targetCI->bpLength.mean + (3 * sqrt(targetCI->bpLength.variance)));
      if(branchTargetsEnd->hi > maxBound){
        maxBound = branchTargetsEnd->hi;
      }
      branchTargetsEnd->orient = GetEdgeOrientationWRT(edge, thisCI->id);
      branchTargetsEnd++;
    }
    if((branchEnd - branchEdges) != numBranch){
      fprintf(GlobalData->stderrc, "RecSub numBranch %d,%d\n",
	      (int)(branchEnd - branchEdges), numBranch);
      DumpChunkInstance(GlobalData->stderrc, graph, thisCI, FALSE, TRUE, TRUE, FALSE);
      fflush(NULL);
      assert(FALSE);
    }
    graphCIEdges = graph->RezGraph->edges;
    qsort((void *)branchEdges, numBranch,
          sizeof(*branchEdges), CompareCIEdgeTMeans);
#if 0
    if(Find_Olap_Path(thisCI, end, (ChunkInstanceT *)NULL, numBranch,
                      branchTargets, maxBound, &firstTarget, &numTargetsBestPath,
                      &firstTargetBestPath, &targetPosition, USE_TANDEM_OLAPS)){
      CDS_CID_t targetIndex;
      int outOfOrder = FALSE;
      assert(numTargetsBestPath > 0);
      assert(numTargetsBestPath <= numBranch);
      assert(firstTargetBestPath >= 0);
      assert(firstTargetBestPath < numBranch);
      if(numTargetsBestPath < numBranch){
        fprintf(GlobalData->stderrc, "RecSub numTargetsBestPath(%d) < numBranch(%d)\n",
                numTargetsBestPath, numBranch);
      }else{
        fprintf(GlobalData->stderrc, "RecSub numTargetsBestPath(%d) = numBranch(%d)\n",
                numTargetsBestPath, numBranch);
      }
      for(targetIndex = firstTargetBestPath, branchTarget = branchEdges;
          targetIndex != NULLINDEX;
          targetIndex = branchTargets[targetIndex].next, branchTarget++){
        CIEdgeT *targetEdge = GetGraphEdge(graph->RezGraph, *branchTarget);
        ChunkInstanceT *sortTargetCI = GetGraphNode(graph->RezGraph,
                                                    ((targetEdge->idA == thisCI->id) ? targetEdge->idB : targetEdge->idA));

        if(sortTargetCI->id != branchTargets[targetIndex].id){
          outOfOrder = TRUE;
        }
        fprintf(GlobalData->stderrc, "Path: %d %f\n",
                branchTargets[targetIndex].id, branchTargets[targetIndex].where);
        fprintf(GlobalData->stderrc, "Sort: " F_CID " %f %f\n",
                sortTargetCI->id,
                sortTargetCI->bpLength.mean + targetEdge->distance.mean,
                sqrt(sortTargetCI->bpLength.variance +
                     targetEdge->distance.variance));
      }
      for(; branchTarget < branchEnd; branchTarget++){
        CIEdgeT *targetEdge = GetGraphEdge(graph->RezGraph, *branchTarget);
        ChunkInstanceT *sortTargetCI = GetGraphNode(graph->RezGraph,
                                                    ((targetEdge->idA == thisCI->id) ? targetEdge->idB : targetEdge->idA));

        fprintf(GlobalData->stderrc, "Sort: " F_CID " %f %f\n",
                sortTargetCI->id,
                sortTargetCI->bpLength.mean + targetEdge->distance.mean,
                sqrt(sortTargetCI->bpLength.variance +
                     targetEdge->distance.variance));
      }
      if(outOfOrder){
        fprintf(GlobalData->stderrc, "RecSub Targets on Best Path differ from sort.\n");
      }
    }
#endif

    /*
      After sorting essential edges from nearest to farthest,
      1. set sourceCI to nearest CI - intent is to keep existing edge between thisCI & sourceCI
      2. iterate over all other (longer) edges as targetEdge's to targetCI's
      attempt to establish inferred or real edge with acceptable chi2 between sourceCI & targetCI
      if that works, depracate the edge between thisCI and targetCI
    */
    branchSource = branchEdges;
    sourceEdge = GetGraphEdge(graph->RezGraph, *branchSource);
    sourceCI = GetGraphNode(graph->RezGraph,
                            ((sourceEdge->idA == thisCI->id) ? sourceEdge->idB :
                             sourceEdge->idA));
    sourceMean = sourceCI->bpLength.mean + sourceEdge->distance.mean;
    sourceVariance = sourceCI->bpLength.variance + sourceEdge->distance.variance;
    sourceEdgeOrient = GetEdgeOrientationWRT(sourceEdge, thisCI->id);
    if(verbose){
      fprintf(GlobalData->stderrc, "RecSub source SoId:" F_CID "(%d,%d) cId:" F_CID " orient:%c\n",
              sourceCI->id, sourceCI->numEssentialA, sourceCI->numEssentialB,
              thisCI->id, (char)sourceEdgeOrient);
    }
    sourceCI->flags.bits.smoothSeenAlready = TRUE;
    for(branchTarget = branchEdges + 1, failedInfer = FALSE,
          numEssentialAdded = 0, numEssentialRemoved = 0;
        branchTarget < branchEnd;
        branchTarget++){
      CIEdgeT *targetEdge = GetGraphEdge(graph->RezGraph, *branchTarget);
      ChunkInstanceT *targetCI = GetGraphNode(graph->RezGraph,
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

      existingEdge = FindEdgeBetweenCIsChiSquare(graph->RezGraph, sourceCI,
                                                 (targetEdge->idA == thisCI->id) ? targetEdge->idB :
                                                 targetEdge->idA,
                                                 inferredEdgeOrient, inferredMean, inferredVariance,
                                                 &chiSquaredValue, (float)PAIRWISECHI2THRESHOLD_CGW,
                                                 &existingIsEssential, verbose);
      targetEdge = GetGraphEdge(graph->RezGraph, *branchTarget); // FindEdgeBetweenCIsChiSquare may cause edges to realloc
      if(existingEdge != (CIEdgeT *)NULL){
        if(verbose){
          fprintf(GlobalData->stderrc, "RecSub targetE TaId:" F_CID "(%d,%d) cId:" F_CID " orient:%c\n",
                  targetCI->id, targetCI->numEssentialA, targetCI->numEssentialB,
                  thisCI->id, (char)targetEdgeOrient);
        }
        existingEdge->flags.bits.wasEssential =
          getEssentialEdgeStatus(existingEdge);
        if(!existingEdge->flags.bits.isEssential){
          setEssentialEdgeStatus(existingEdge, TRUE);
          //	existingEdge->flags.bits.isEssential = TRUE;
          //	assert(!isSloppyEdge(existingEdge));
          numEssentialAdded++;
        }
        if(verbose){
          PrintGraphEdge(GlobalData->stderrc,graph->RezGraph,"tE Before", targetEdge, thisCI->id);
        }
        targetEdge->flags.bits.isInferredRemoved = TRUE;
        setEssentialEdgeStatus(targetEdge, FALSE);
        //targetEdge->flags.bits.isEssential = FALSE;
        if(verbose){
          PrintGraphEdge(GlobalData->stderrc,graph->RezGraph,"tE After", targetEdge, thisCI->id);
        }
        if(targetEdge->flags.bits.isInferred && targetEdge->flags.bits.isTentative){
          DeleteGraphEdge(graph->RezGraph, targetEdge);
        }
        numEssentialRemoved++;
        *existingEnd++ = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->RezGraph->edges, existingEdge);
        continue;
      }else if(existingIsEssential){
        failedInfer = TRUE;
        break;
      }
      if(inferredMean < - CGW_MISSED_OVERLAP){
        // adjust mean and variance so things don't get sloppy (dewim 09/12/01)
        double new_stddev = (inferredDistance.mean + CGW_MISSED_OVERLAP + 3. * sqrt(inferredDistance.variance)) / 3.;
        new_stddev = (new_stddev < 1.) ? 1. : new_stddev;
        inferredDistance.variance = new_stddev * new_stddev;
        inferredDistance.mean = - CGW_MISSED_OVERLAP;

#if 0
        if ( inferredVariance > FLT_MAX ) inferredVariance = FLT_MAX;  // hack instituted on 5/15/01 to help mouse_20010508
#endif

        if(!PairwiseChiSquare((float)inferredMean, (float)inferredVariance,
                              (float)(- CGW_MISSED_OVERLAP), (float)1.0,
                              (LengthT *)NULL, &chiSquaredValue,
                              (float)PAIRWISECHI2THRESHOLD_CGW)){
          failedInfer = TRUE;
          if(verbose){
            fprintf(GlobalData->stderrc, "RecSub ChiSquareN:%f,%f,%f,%f,%f SoId:" F_CID " TaId:" F_CID "\n",
                    chiSquaredValue, inferredMean,
                    (float)sqrt((double)inferredVariance),
                    (float)(- CGW_MISSED_OVERLAP), (float)1.0,
                    sourceCI->id, targetCI->id);
          }
          break;
        }
      }
      if(verbose){
        fprintf(GlobalData->stderrc, "RecSub target TaId:" F_CID "(%d,%d) cId:" F_CID " orient:%c\n",
                targetCI->id, targetCI->numEssentialA, targetCI->numEssentialB,
                thisCI->id, (char)targetEdgeOrient);
      }
      targetEdge->flags.bits.isInferredRemoved = TRUE;
      setEssentialEdgeStatus(targetEdge, FALSE);
      if(targetEdge->flags.bits.isInferred && targetEdge->flags.bits.isTentative){
        DeleteGraphEdge(graph->RezGraph, targetEdge);
      }
      numEssentialRemoved++;
      inferredEdgeIndex = AddGraphEdge(graph->RezGraph, sourceCI->id, targetCI->id,
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
      inferredEdge = GetGraphEdge(graph->RezGraph, inferredEdgeIndex);
      AssertPtr(inferredEdge);
      numEssentialAdded++;
      inferredEdge->flags.bits.isInferred = TRUE;
      inferredEdge->flags.bits.isTentative = TRUE;
      setEssentialEdgeStatus(inferredEdge, TRUE);

#if 0
      fprintf(GlobalData->stderrc,"* Added inferred GraphEdge v2 (%d,%d,%c) isAContainsB:%d isBContainsA:%d fragA %d fragB %d isInferred %d isRaw %d\n",
              inferredEdge->idA, inferredEdge->idB,
              inferredEdge->orient,
              inferredEdge->flags.bits.aContainsB,
              inferredEdge->flags.bits.bContainsA,
              inferredEdge->fragA,inferredEdge->fragB,
              inferredEdge->flags.bits.isInferred,
              inferredEdge->flags.bits.isRaw);
#endif


    }
    if(failedInfer){

#ifdef INSTRUMENT_SMOOTHED
      // instrument the failed 'scaffold' before changing essential edges back
      InstrumentContigPath(graph, smoothed_si,
                           firstID, instFirstEnd, sourceCI->id);
      PrintScaffoldInstrumenter(graph, smoothed_si,
                                InstrumenterVerbose2,
                                "", instSmoothFailurefp);
#endif

      if(verbose){
        fprintf(GlobalData->stderrc,"Failed RecSub on node " F_CID "\n", thisCI->id);
      }
      sourceCI->flags.bits.smoothSeenAlready = FALSE;
      sourceCI->smoothExpectedCID = NULLINDEX;
      for(addedPtr = addedInferred; addedPtr < addedEnd; addedPtr++){
        CIEdgeT *addedEdge = GetGraphEdge(graph->RezGraph, *addedPtr);
        if(!addedEdge->flags.bits.isDeleted){
          DeleteGraphEdge(graph->RezGraph, addedEdge);
        }
      }
      for(existingPtr = existingEdges; existingPtr < existingEnd; existingPtr++){
        CIEdgeT *existingEdge = GetGraphEdge(graph->RezGraph, *existingPtr);
        setEssentialEdgeStatus(existingEdge, existingEdge->flags.bits.wasEssential);
        //      existingEdge->flags.bits.isEssential =	existingEdge->flags.bits.wasEssential;
      }
      for(branchTarget = branchEdges + 1; branchTarget < branchEnd;
          branchTarget++){
        CIEdgeT *targetEdge = GetGraphEdge(graph->RezGraph, *branchTarget);
        ChunkInstanceT *targetCI;

        if(targetEdge->flags.bits.isDeleted){ // SAK
          //	fprintf(stderr,"* Encountered Deleted Edge at failedInfer in TransitiveReduction\n");
          continue;
        }
        targetCI = GetGraphNode(graph->RezGraph,
                                ((targetEdge->idA == thisCI->id) ?
                                 targetEdge->idB : targetEdge->idA));
        targetCI->flags.bits.smoothSeenAlready = FALSE;
        targetCI->smoothExpectedCID = NULLINDEX;
        targetEdge->flags.bits.isInferredRemoved = FALSE;
        setEssentialEdgeStatus(targetEdge, TRUE);
        //      targetEdge->flags.bits.isEssential = TRUE;
        //      assert(!isSloppyEdge(targetEdge));

      }
#ifndef ALLOC_MEMORY_ON_STACK
      safe_free(branchEdges);
      safe_free(existingEdges);
      safe_free(addedInferred);
      safe_free(branchTargets);
#endif
      if(verbose){
        fprintf(GlobalData->stderrc, "RecSub FailedInfer scId:" F_CID "(%d,%d) cId:" F_CID "(%d,%d)\n",
                sourceCI->id, sourceCI->numEssentialA, sourceCI->numEssentialB,
                thisCI->id, thisCI->numEssentialA, thisCI->numEssentialB);
      }
      return(FALSE);
    }else{
      CDS_CID_t sourceEnd;
      if(end == B_END){
        thisCI->numEssentialB -= numEssentialRemoved;
      }else{
        thisCI->numEssentialA -= numEssentialRemoved;
      }
      if((sourceEdgeOrient == AB_AB) || (sourceEdgeOrient == BA_AB)){
        sourceEnd = B_END;
        sourceCI->numEssentialB += numEssentialAdded;
      }else{
        sourceEnd = A_END;
        sourceCI->numEssentialA += numEssentialAdded;
      }
      if(!RecursiveSmoothWithInferredEdges(graph, sourceCI, sourceEnd, verbose,
                                           numInferredAdded,
                                           firstID, instFirstEnd, lastID,
                                           smoothed_si,
                                           instSmoothSuccessfp,
                                           instSmoothFailurefp)){
        if(end == B_END){
          thisCI->numEssentialB += numEssentialRemoved;
        }else{
          thisCI->numEssentialA += numEssentialRemoved;
        }
        if(sourceEnd == B_END){
          sourceCI->numEssentialB -= numEssentialAdded;
        }else{
          sourceCI->numEssentialA -= numEssentialAdded;
        }
        sourceCI->flags.bits.smoothSeenAlready = FALSE;
        sourceCI->smoothExpectedCID = NULLINDEX;
        for(addedPtr = addedInferred; addedPtr < addedEnd; addedPtr++){
          CIEdgeT *addedEdge = GetGraphEdge(graph->RezGraph, *addedPtr);
          if(!addedEdge->flags.bits.isDeleted){
            DeleteGraphEdge(graph->RezGraph, addedEdge);
          }
        }
        for(existingPtr = existingEdges; existingPtr < existingEnd;
            existingPtr++){
          CIEdgeT *existingEdge = GetGraphEdge(graph->RezGraph, *existingPtr);
          setEssentialEdgeStatus(existingEdge, existingEdge->flags.bits.wasEssential);
          //	existingEdge->flags.bits.isEssential =
          //	  existingEdge->flags.bits.wasEssential;
        }
        for(branchTarget = branchEdges + 1; branchTarget < branchEnd;
            branchTarget++){
          CIEdgeT *targetEdge = GetGraphEdge(graph->RezGraph, *branchTarget);
          ChunkInstanceT *targetCI;

          if(targetEdge->flags.bits.isDeleted){ // SAK
            //	  fprintf(stderr,"* Encountered Deleted Edge at cleanup in TransitiveReduction\n");
            continue;
          }
          targetCI = GetGraphNode(graph->RezGraph,
                                  ((targetEdge->idA == thisCI->id) ?
                                   targetEdge->idB : targetEdge->idA));
          targetCI->flags.bits.smoothSeenAlready = FALSE;
          targetCI->smoothExpectedCID = NULLINDEX;
          targetEdge->flags.bits.isInferredRemoved = FALSE;
          assert(!isSloppyEdge(targetEdge));
          setEssentialEdgeStatus(targetEdge, TRUE);
          //	targetEdge->flags.bits.isEssential = TRUE;
        }
#ifndef ALLOC_MEMORY_ON_STACK
        safe_free(branchEdges);
        safe_free(existingEdges);
        safe_free(addedInferred);
        safe_free(branchTargets);
#endif
        if(verbose){
          fprintf(GlobalData->stderrc, "RecSub Failed scId:" F_CID "(%d,%d) cId:" F_CID "(%d,%d)\n",
                  sourceCI->id, sourceCI->numEssentialA, sourceCI->numEssentialB,
                  thisCI->id, thisCI->numEssentialA, thisCI->numEssentialB);
        }
        return(FALSE);
      }else{
        sourceCI->flags.bits.smoothSeenAlready = FALSE;
        sourceCI->smoothExpectedCID = NULLINDEX;
        if(end == B_END){
          thisCI->essentialEdgeB = *branchSource;
        }else{
          thisCI->essentialEdgeA = *branchSource;
        }
        if((sourceEdgeOrient == AB_AB) || (sourceEdgeOrient == BA_AB)){
          sourceCI->essentialEdgeA = *branchSource;
        }else{
          sourceCI->essentialEdgeB = *branchSource;
        }
        for(addedPtr = addedInferred; addedPtr < addedEnd; addedPtr++){
          CIEdgeT *addedEdge = GetGraphEdge(graph->RezGraph, *addedPtr);
          if(!addedEdge->flags.bits.isDeleted){
            addedEdge->flags.bits.isTentative = FALSE;
          }
        }
        for(existingPtr = existingEdges; existingPtr < existingEnd;
            existingPtr++){
          CIEdgeT *existingEdge = GetGraphEdge(graph->RezGraph, *existingPtr);
          if(existingEdge->flags.bits.wasEssential){
            ChunkInstanceT *targetCI =
              GetGraphNode(graph->RezGraph,
                           ((existingEdge->idA == sourceCI->id) ?
                            existingEdge->idB : existingEdge->idA));
            ChunkOrientationType existingEdgeOrient =
              GetEdgeOrientationWRT(existingEdge, sourceCI->id);
            // Fix numEssential
            if((existingEdgeOrient == AB_AB) || (existingEdgeOrient == BA_AB)){
              targetCI->numEssentialA--;
            }else{
              targetCI->numEssentialB--;
            }
            if(verbose){
              fprintf(GlobalData->stderrc, "RecSub numEssentialFixed-- taId:" F_CID "(%d,%d)\n",
                      targetCI->id, targetCI->numEssentialA,
                      targetCI->numEssentialB);
            }
          }
        }
        if(verbose){
          fprintf(GlobalData->stderrc, "RecSub numEssentialAdded:%d scId:" F_CID "(%d,%d) cId:" F_CID "(%d,%d)\n",
                  numEssentialAdded, sourceCI->id, sourceCI->numEssentialA,
                  sourceCI->numEssentialB, thisCI->id, thisCI->numEssentialA,
                  thisCI->numEssentialB);
        }
#ifndef ALLOC_MEMORY_ON_STACK
        safe_free(branchEdges);
        safe_free(existingEdges);
        safe_free(addedInferred);
        safe_free(branchTargets);
#endif
        return(TRUE);
      }
    }
  }
}

int CountEssentialEdgesOneEnd(ScaffoldGraphT *graph, ChunkInstanceT *thisCI,
                              int end, CDS_CID_t *essentialEdge){

  int count = 0;
  GraphEdgeIterator edges;
  CIEdgeT *edge;
  InitGraphEdgeIterator(graph->RezGraph, thisCI->id, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
  *essentialEdge = (CDS_CID_t)NULLINDEX;
  while((edge = NextGraphEdgeIterator(&edges))!= NULL){
    if(!edge->flags.bits.isActive || !edge->flags.bits.isConfirmed ||
       !edge->flags.bits.isUniquetoUnique ||
       edge->flags.bits.isTransitivelyRemoved ||
       edge->flags.bits.isRedundantRemoved){
      // Inactive, unconfirmed, and edges between nonUnique CIs are not
      // essential.
      setEssentialEdgeStatus(edge, FALSE);
      continue;
    }
    setEssentialEdgeStatus(edge, TRUE);
    count++;
    *essentialEdge = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->RezGraph->edges, edge);
  }
  return(count);
}


void WriteEssentialEdgeCGM(ScaffoldGraphT * graph,
                           char * suffix)
{
  GraphNodeIterator nodes;
  ContigT *thisCI;
  char filename[1024];
  FILE * essentialfp;

  sprintf(filename, "essentialEdges%s.cgm", suffix);
  essentialfp = fopen(filename, "w");
  assert(essentialfp != NULL);

  fprintf(essentialfp,
          "#Essential edges for A & B ends of contigs\n");

  InitGraphNodeIterator(&nodes, graph->RezGraph,
                        GRAPH_NODE_UNIQUE_ONLY );
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL)
    {
      fprintf(essentialfp,  F_CID "\t%d\t%d\n",
              thisCI->id,
              thisCI->numEssentialA,
              thisCI->numEssentialB);
      /*
        if(thisCI->numEssentialA > 1)
        fprintf(essentialfp, "%d\n", thisCI->numEssentialA);
        if(thisCI->numEssentialB > 1)
        fprintf(essentialfp, "%d\n", thisCI->numEssentialB);
      */
    }
  /*
    ChunkInstanceT * contig;
    GraphNodeIterator contigIterator;
    ChunkInstanceT * unitig;
    ContigTIterator unitigIterator;
    int numContigs = 0;
    int numRepetitiveContigs = 0;
    VA_TYPE(CDS_CID_t) * aEndIDs;
    VA_TYPE(CDS_CID_t) * bEndIDs;
    int numOffAEnd;
    int numOffBEnd;
    int maxOffEnd;

    aEndIDs = CreateVA_CDS_CID_t(10000);
    bEndIDs = CreateVA_CDS_CID_t(10000);

    InitGraphNodeIterator(&contigIterator,
    graph->ContigGraph,
    GRAPH_NODE_UNIQUE_ONLY);
    while((contig = NextGraphNodeIterator(&contigIterator)) != NULL)
    {
    int coverageStat;
    GraphEdgeIterator edges;
    CIEdgeT * edge;

    numContigs++;
    ResetVA_CDS_CID_t(aEndIDs);
    ResetVA_CDS_CID_t(bEndIDs);

    // loop over contig edges
    InitGraphEdgeIterator(graph->ContigGraph,
    contig->id,
    ALL_END,
    ALL_EDGES,
    GRAPH_EDGE_ESSENTIAL_ONLY,
    &edges);
    while((edge = NextGraphEdgeIterator(&edges)) != NULL)
    {
    int isA = (edge->idA == contig->id);
    CDS_CID_t othercid = (isA? edge->idB: edge->idA);

    if((isA && (edge->orient == AB_AB || edge->orient == AB_BA)) ||
    (!isA && (edge->orient == AB_BA || edge->orient == BA_BA)))
    AppendVA_CDS_CID_t(bEndIDs, &othercid);
    else
    AppendVA_CDS_CID_t(aEndIDs, &othercid);
    }

    // now sort & examine edges
    if(GetNumVA_CDS_CID_t(aEndIDs) > 1)
    numOffAEnd = CountUniqueIDs(aEndIDs);
    else
    numOffAEnd = GetNumVA_CDS_CID_t(aEndIDs);

    if(GetNumVA_CDS_CID_t(bEndIDs) > 1)
    numOffBEnd = CountUniqueIDs(bEndIDs);
    else
    numOffBEnd = GetNumVA_CDS_CID_t(bEndIDs);

    fprintf(essentialfp,  F_CID "\t%d\t%d\n",
    contig->id, numOffAEnd, numOffBEnd);
    }

    DeleteVA_CDS_CID_t(aEndIDs);
    DeleteVA_CDS_CID_t(bEndIDs);
  */
  fclose(essentialfp);
}


void SmoothWithInferredEdges(ScaffoldGraphT *graph,
                             int markShakyBifurcations,
                             int verbose){

  GraphNodeIterator nodes;
  ContigT *thisCI;
  int cnt = 1;
  int smooth_success;
  int numInferredAdded;

  InitGraphNodeIterator(&nodes, graph->RezGraph, GRAPH_NODE_UNIQUE_ONLY );
  if(verbose){
    fprintf(GlobalData->stderrc, "Smoothing: Unique CIs:\n");
  }
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL){

    if (verbose) {
      if((cnt++ % 10000) == 0){
        fprintf(stderr,"* %d SmoothWithInferredEdges " F_CID "\n",
                cnt, thisCI->id);
        fflush(stderr);
      }
    }
    thisCI->smoothExpectedCID = NULLINDEX;
    thisCI->flags.bits.smoothSeenAlready = FALSE;
    // Do the A_END edges.
    thisCI->numEssentialA = CountEssentialEdgesOneEnd(graph, thisCI, A_END,
                                                      &(thisCI->essentialEdgeA));
    assert((thisCI->numEssentialA == 0) ||
	   (thisCI->essentialEdgeA != NULLINDEX));
    // Do the B_END edges.
    thisCI->numEssentialB = CountEssentialEdgesOneEnd(graph, thisCI, B_END,
                                                      &(thisCI->essentialEdgeB));
    assert((thisCI->numEssentialB == 0) ||
	   (thisCI->essentialEdgeB != NULLINDEX));
    if(verbose){
      fprintf(GlobalData->stderrc, F_CID "(%d{" F_CID "},%d{" F_CID "})\n",
	      thisCI->id, thisCI->numEssentialA, thisCI->essentialEdgeA,
	      thisCI->numEssentialB, thisCI->essentialEdgeB);
    }
  }
  if(!markShakyBifurcations)
    {
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
#ifdef MEASURE_GRAPH_COMPLEXITY
      FILE * successfp;
      FILE * failurefp;
      int file_count;
      char filename[1024];

      for(file_count = 0; file_count < 1000; file_count++)
        {
          // create cgms of inferred edges added in successes AND failures
          sprintf(filename, "inferred_edges_added_%05d.cgm", file_count);
          if((successfp = fopen(filename, "r")) == NULL)
            {
              successfp = fopen(filename, "w");
              sprintf(filename, "inferred_edges_not_added_%05d.cgm", file_count);
              failurefp = fopen(filename, "w");
              break;
            }
          else
            fclose(successfp);
        }
      assert(successfp != NULL && failurefp != NULL);
      fprintf(successfp, "Inferred Edges Added in Smoothing\n");
      fprintf(failurefp, "Inferred Edges In Failed Smoothing Attempt\n");
      /* also, create cgm of num essential edges off each end of each contig
         before smoothing, then after
      */
      sprintf(filename, "%05d_before", file_count);
      WriteEssentialEdgeCGM(graph, filename);
#endif
#ifdef INSTRUMENT_SMOOTHED
      {
        int file_count;
        char filename[1024];
        for(file_count = 0; file_count < 1000; file_count++)
          {
            sprintf(filename, "smoothSuccess_%05d.itxt", file_count);
            instSmoothSuccessfp = fopen(filename, "r");
            if(instSmoothSuccessfp == NULL)
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

      InitGraphNodeIterator(&nodes, graph->RezGraph, GRAPH_NODE_UNIQUE_ONLY );
      while((thisCI = NextGraphNodeIterator(&nodes)) != NULL){
        /* Substitute inferred edges where the graph branches in an attempt to
           eliminate as many branches as possible. */
        if(thisCI->numEssentialA > 1){
          thisCI->flags.bits.smoothSeenAlready = TRUE;
          if(verbose){
            fprintf(GlobalData->stderrc,"Calling RecSmooth on node " F_CID " A_END\n",
                    thisCI->id);
          }

#ifdef INSTRUMENT_CGW
          num_before = InstrumentContigEnd(graph, si_before, thisCI, A_END);
#endif
          firstID = thisCI->id;
          lastID = NULLINDEX;
          instFirstEnd = A_END;
          if(verbose)
            fprintf(stderr,  F_CID ": ", thisCI->id);


          numInferredAdded = 0;
          smooth_success =
            RecursiveSmoothWithInferredEdges(graph, thisCI, A_END, verbose,
                                             &numInferredAdded,
                                             firstID, instFirstEnd, &lastID,
                                             smoothed_si,
                                             instSmoothSuccessfp,
                                             instSmoothFailurefp);
          thisCI->flags.bits.smoothSeenAlready = FALSE;


#ifdef MEASURE_GRAPH_COMPLEXITY
          if(smooth_success)
            fprintf(successfp, F_CID "\n", numInferredAdded);
          else
            fprintf(failurefp, F_CID "\n", numInferredAdded);
#endif

#ifdef INSTRUMENT_CGW
          if(smooth_success == TRUE)
            {
              num_after = InstrumentContigEnd(graph, si_after, thisCI, A_END);
              fprintf(GlobalData->stderrc, "Smoothed CI " F_CID " on A_END\n",
                      thisCI->id);
              fprintf(GlobalData->stderrc,
                      "contigs: %d before (may include duplicates), %d after\n",
                      num_before, num_after);
              if(CompareMateInstrumenters(&(si_before->mates),
                                          &(si_after->mates),
                                          InstrumenterVerbose5,
                                          GlobalData->stderrc) ==
                 InstrumenterWorse)
                {
                  int num_partial;
                  for(num_partial = num_after - 1; num_partial > 1; num_partial--)
                    {
                      // assumes there is one essential edge off end
                      InstrumentContigEndPartial(graph, si_after, thisCI, A_END,
                                                 num_partial);
                      fprintf(GlobalData->stderrc, "re-instrumenting %d contigs\n",
                              num_partial);
                      if(CompareMateInstrumenters(&(si_before->mates),
                                                  &(si_after->mates),
                                                  InstrumenterVerbose5,
                                                  GlobalData->stderrc) !=
                         InstrumenterWorse)
                        {
                          break;
                        }
                    }
                }
            }
#endif
        } // if(thisCI->numEssentialA > 1){
        if(thisCI->numEssentialB > 1){
          thisCI->flags.bits.smoothSeenAlready = TRUE;
          if(verbose){
            fprintf(GlobalData->stderrc,"Calling RecSmooth on node " F_CID " B_END\n",
                    thisCI->id);
          }
#ifdef INSTRUMENT_CGW
          num_before = InstrumentContigEnd(graph, si_before, thisCI, B_END);
#endif
          firstID = thisCI->id;
          lastID = NULLINDEX;
          instFirstEnd = B_END;
          if(verbose)
            fprintf(stderr, F_CID ": ", thisCI->id);

          numInferredAdded = 0;
          smooth_success =
            RecursiveSmoothWithInferredEdges(graph, thisCI, B_END, verbose,
                                             &numInferredAdded,
                                             firstID, instFirstEnd, &lastID,
                                             smoothed_si,
                                             instSmoothSuccessfp,
                                             instSmoothFailurefp);
          thisCI->flags.bits.smoothSeenAlready = FALSE;


#ifdef MEASURE_GRAPH_COMPLEXITY
          if(smooth_success)
            fprintf(successfp, "%d\n", numInferredAdded);
          else
            fprintf(failurefp,  "%\n", numInferredAdded);
#endif

#ifdef INSTRUMENT_CGW
          if(smooth_success == TRUE)
            {
              num_after = InstrumentContigEnd(graph, si_after, thisCI, B_END);
              fprintf(GlobalData->stderrc, "Smoothed CI " F_CID " on B_END\n",
                      thisCI->id);
              fprintf(GlobalData->stderrc,
                      "contigs: %d before (may include duplicates), %d after\n",
                      num_before, num_after);
              if(CompareMateInstrumenters(&(si_before->mates),
                                          &(si_after->mates),
                                          InstrumenterVerbose5,
                                          GlobalData->stderrc) ==
                 InstrumenterWorse)
                {
                  int num_partial;
                  for(num_partial = num_after - 1; num_partial > 1; num_partial--)
                    {
                      // assumes there is one essential edge off end
                      InstrumentContigEndPartial(graph, si_after, thisCI, B_END,
                                                 num_partial);
                      fprintf(GlobalData->stderrc, "re-instrumenting %d contigs\n",
                              num_partial);
                      if(CompareMateInstrumenters(&(si_before->mates),
                                                  &(si_after->mates),
                                                  InstrumenterVerbose5,
                                                  GlobalData->stderrc) !=
                         InstrumenterWorse)
                        {
                          break;
                        }
                    }
                }
            }
#endif
        } //       if(thisCI->numEssentialB > 1){
      } // while(thisCI = NextGraphNodeIterator(&nodes)){
#ifdef INSTRUMENT_CGW
      DestroyScaffoldInstrumenter(si_before);
      DestroyScaffoldInstrumenter(si_after);
#endif
#ifdef MEASURE_GRAPH_COMPLEXITY
      fclose(successfp);
      fclose(failurefp);
      sprintf(filename, "%05d_before", file_count);
      sprintf(filename, "%05d_after", file_count);
      WriteEssentialEdgeCGM(graph, filename);
#endif
#ifdef INSTRUMENT_SMOOTHED
      fclose(instSmoothSuccessfp);
      fclose(instSmoothFailurefp);
#endif
#ifdef INSTRUMENT_TRANS_REDUCED
      DestroyScaffoldInstrumenter(smoothed_si);
#endif
    } // if(!markShakyBifurcations){
  return;
}

void DeleteInferredEdgesForThisCI(ScaffoldGraphT *graph,
                                  CDS_CID_t cid,
				  int verbose){

  GraphEdgeIterator edges;
  EdgeCGW_T *edge;

  // Use merged edges
  InitGraphEdgeIterator(graph->RezGraph, cid, ALL_END,
                        ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
  while((edge = NextGraphEdgeIterator(&edges)) != NULL){
    if(edge->flags.bits.isInferred){
      DeleteGraphEdge(graph->RezGraph, edge);
    }
  }
  return;
}

void SymmetricNeighbors(ChunkInstanceT *thisCI, ChunkOrientationType orient, CDS_CID_t edgeID){

  if((orient == AB_AB) || (orient == BA_AB)){// A end of this CI
    assert(thisCI->numEssentialA == 1);
    assert(edgeID == thisCI->essentialEdgeA);
  }else{// B end of this CI
    assert(thisCI->numEssentialB == 1);
    assert(edgeID == thisCI->essentialEdgeB);
  }
}


int NeighborBranches(ChunkInstanceT *thisCI,
                     ChunkOrientationType orient,
                     CDS_CID_t edgeID){

  if((orient == AB_AB) || (orient == BA_AB)){// A end of this CI
    if(thisCI->numEssentialA != 1){
      return(TRUE);
    }else{
      assert(edgeID == thisCI->essentialEdgeA);
      return(FALSE);
    }
  }else{// B end of this CI
    assert(thisCI->numEssentialB >= 1);
    if(thisCI->numEssentialB != 1){
      return(TRUE);
    }else{
      assert(edgeID == thisCI->essentialEdgeB);
      return(FALSE);
    }
  }
}


/****************************************************************************/
void DetectScaffoldCycles(CDS_CID_t *currentScaffoldID){
  {/* Look for cycles */
    ChunkInstanceT *thisCI;
    GraphNodeIterator nodes;

    InitGraphNodeIterator(&nodes, ScaffoldGraph->RezGraph, GRAPH_NODE_UNIQUE_ONLY);
    while((thisCI = NextGraphNodeIterator(&nodes)) != NULL){
      ChunkInstanceT *AendCI, *BendCI, *startingCI;
      CIEdgeT *edgeA, *edgeB, *edge;
      ChunkInstanceT *neighbor;
      ChunkOrient orientCI;
      CDS_CID_t cid;

      cid = thisCI->id;

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

      if(thisCI->scaffoldID != NULLINDEX){
	thisCI->scaffoldID = NULLINDEX;
	continue;// This CI has already been placed in a Scaffold.
      }

      // what this leaves behind is cycles; so any remaining unique CIs
      // must be in a cycle.
      assert(thisCI->essentialEdgeA != NULLINDEX);
      edgeA = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeA);
      AendCI = GetGraphNode(ScaffoldGraph->RezGraph,
			    (edgeA->idA == thisCI->id) ? edgeA->idB : edgeA->idA);
      assert(thisCI->essentialEdgeB != NULLINDEX);
      edgeB = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeB);
      BendCI = GetGraphNode(ScaffoldGraph->RezGraph,
			    (edgeB->idA == thisCI->id) ? edgeB->idB : edgeB->idA);
      fprintf(stderr, "*** Found a Cycle in Essential Edge Graph\ncId:" F_CID " Aend:" F_CID " Bend:" F_CID "\n",
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
	if(orientCI == A_B){
	  if(edgeOrient == AB_AB){
	    orientCI = A_B;
	  }else if(edgeOrient == AB_BA){
	    orientCI = B_A;
	  }else{
	    assert(0);
	  }
	}else{// orientCI == B_A
	  if(edgeOrient == BA_AB){
	    orientCI = A_B;
	  }else if(edgeOrient == BA_BA){
	    orientCI = B_A;
	  }else{
	    assert(0);
	  }
	}
	thisCI->essentialEdgeA = NULLINDEX;
	thisCI->essentialEdgeB = NULLINDEX;
	//	if(GetNodeType(thisCI) == UNIQUECHUNK_CGW){
	//	  thisCI->scaffoldID = NULLINDEX;
	//	  ResetNodeType(thisCI);
	//	  DeleteInferredEdgesForThisCI(graph, thisCI->id, verbose);
	//	}else{
	//	  assert(GetNodeType(thisCI) == DISCRIMINATORUNIQUECHUNK_CGW);

	// the actual id assigned here is immaterial; what matters is
	// that all contigs in the cycle get some value not NULLINDEX
	// so that as each one is visited it will pass the test for
	// id != NULLINDEX, have the id set to NULLINDEX, and then
	// continue (see top of outer while loop)
        thisCI->scaffoldID = *currentScaffoldID;
        //	}
	thisCI = neighbor;


	if(orientCI == A_B){
	  if(thisCI->essentialEdgeB != NULLINDEX){
	    edge = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeB);
	    neighbor = GetGraphNode(ScaffoldGraph->RezGraph,
				    (edge->idA == thisCI->id) ? edge->idB : edge->idA);
	  }else{// End of Scaffold
	    edge = (CIEdgeT *)NULL;
	    neighbor = (ChunkInstanceT *)NULL;
	  }
	}else{// orientCI == B_A
	  if(thisCI->essentialEdgeA != NULLINDEX){
	    edge = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeA);
	    neighbor = GetGraphNode(ScaffoldGraph->RezGraph,
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
}

/****************************************************************************/
void LabelCIScaffolds(CDS_CID_t *currentScaffoldID, int verbose){
  /* Label which scaffold a CI would go into but do not create the
     scaffold as a way to detect cycles. */
  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->RezGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL){
    ChunkInstanceT *AendCI, *BendCI;
    CIEdgeT *edgeA, *edgeB, *edge;
    ChunkInstanceT *neighbor;
    ChunkOrient orientCI;
    CDS_CID_t cid;

    cid = thisCI->id;
    if(thisCI->scaffoldID != NULLINDEX){
      continue;// This CI has already been placed in a Scaffold.
    }
    if(thisCI->essentialEdgeA != NULLINDEX){
      edgeA = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeA);
      AendCI = GetGraphNode(ScaffoldGraph->RezGraph,
			    (edgeA->idA == thisCI->id) ? edgeA->idB : edgeA->idA);
    }else{
      edgeA = (CIEdgeT *)NULL;
      AendCI = (ChunkInstanceT *)NULL;
    }
    if(thisCI->essentialEdgeB != NULLINDEX){
      edgeB = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeB);
      BendCI = GetGraphNode(ScaffoldGraph->RezGraph,
			    (edgeB->idA == thisCI->id) ? edgeB->idB : edgeB->idA);
    }else{
      edgeB = (CIEdgeT *)NULL;
      BendCI = (ChunkInstanceT *)NULL;
    }
    if(verbose){
      fprintf(GlobalData->stderrc, "ScId:" F_CID " cId:" F_CID " Aend:" F_CID " Bend:" F_CID "\n",
	      thisCI->scaffoldID, thisCI->id,
	      AendCI != (ChunkInstanceT *)NULL ? AendCI->id : NULLINDEX,
	      BendCI != (ChunkInstanceT *)NULL ? BendCI->id : NULLINDEX);
    }
    if((AendCI != (ChunkInstanceT *)NULL) &&
       (BendCI != (ChunkInstanceT *)NULL)){
      continue;// This CI is not a starting point for a Scaffold.
    }


    if(BendCI != (ChunkInstanceT *)NULL){
      orientCI = A_B;
      edge = edgeB;
      neighbor = BendCI;
    }else if(AendCI != (ChunkInstanceT *)NULL){
      orientCI = B_A;
      edge = edgeA;
      neighbor = AendCI;
    }else{// Singleton Scaffold
      orientCI = A_B;
      edge = (CIEdgeT *)NULL;
      neighbor = (ChunkInstanceT *)NULL;
    }

    thisCI->scaffoldID = *currentScaffoldID;
    while(neighbor != (ChunkInstanceT *)NULL){
      ChunkOrientationType edgeOrient = GetEdgeOrientationWRT(edge,
							      thisCI->id);
      if(orientCI == A_B){
	if(edgeOrient == AB_AB){
	  orientCI = A_B;
	}else if(edgeOrient == AB_BA){
	  orientCI = B_A;
	}else{
	  assert(0);
	}
      }else{// orientCI == B_A
	if(edgeOrient == BA_AB){
	  orientCI = A_B;
	}else if(edgeOrient == BA_BA){
	  orientCI = B_A;
	}else{
	  assert(0);
	}
      }
      thisCI = neighbor;
      thisCI->scaffoldID = *currentScaffoldID;
      if(orientCI == A_B){
	if(thisCI->essentialEdgeB != NULLINDEX){
	  edge = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeB);
	  neighbor = GetGraphNode(ScaffoldGraph->RezGraph,
				  (edge->idA == thisCI->id) ? edge->idB : edge->idA);
	}else{// End of Scaffold
	  edge = (CIEdgeT *)NULL;
	  neighbor = (ChunkInstanceT *)NULL;
	}
      }else{// orientCI == B_A
	if(thisCI->essentialEdgeA != NULLINDEX){
	  edge = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeA);
	  neighbor = GetGraphNode(ScaffoldGraph->RezGraph,
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

/********************************************************************/
void ActuallyInsertCIsIntoScaffolds(CDS_CID_t *currentScaffoldID,  int verbose)
{
  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;
  InitGraphNodeIterator(&nodes, ScaffoldGraph->RezGraph, GRAPH_NODE_UNIQUE_ONLY);

  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL){
    ChunkInstanceT *AendCI, *BendCI;
    CIEdgeT *edgeA, *edgeB, *edge;
    ChunkInstanceT *neighbor;
    LengthT NullLength = {0.0, 0.0};
    LengthT aEndOffset, bEndOffset, currentOffset;
    CIScaffoldT CIScaffold = {0};
    ChunkOrient orientCI;
    CDS_CID_t cid;
    cid = thisCI->id;

    if(thisCI->scaffoldID != NULLINDEX){
      continue;// This CI has already been placed in a Scaffold.
    }
    if(thisCI->essentialEdgeA != NULLINDEX){
      edgeA = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeA);
      assert(!isSloppyEdge(edgeA));
      AendCI = GetGraphNode(ScaffoldGraph->RezGraph,
                            (edgeA->idA == thisCI->id) ? edgeA->idB : edgeA->idA);
      SymmetricNeighbors(AendCI, GetEdgeOrientationWRT(edgeA, thisCI->id),thisCI->essentialEdgeA );

    }else{
      edgeA = (CIEdgeT *)NULL;
      AendCI = (ChunkInstanceT *)NULL;
    }
    if(thisCI->essentialEdgeB != NULLINDEX){
      edgeB = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeB);
      assert(!isSloppyEdge(edgeB));
      assert(edgeB->idA == thisCI->id || edgeB->idB == thisCI->id);
      BendCI = GetGraphNode(ScaffoldGraph->RezGraph,
                            (edgeB->idA == thisCI->id) ? edgeB->idB : edgeB->idA);
      SymmetricNeighbors(BendCI, GetEdgeOrientationWRT(edgeB, thisCI->id),thisCI->essentialEdgeB );
    }else{
      edgeB = (CIEdgeT *)NULL;
      BendCI = (ChunkInstanceT *)NULL;
    }
    if(verbose){
      fprintf(GlobalData->stderrc, "ScId:" F_CID " cId:" F_CID " Aend:" F_CID " Bend:" F_CID "\n",
              thisCI->scaffoldID, thisCI->id,
              AendCI != (ChunkInstanceT *)NULL ? AendCI->id : NULLINDEX,
              BendCI != (ChunkInstanceT *)NULL ? BendCI->id : NULLINDEX);
    }
    if((AendCI != (ChunkInstanceT *)NULL) &&
       (BendCI != (ChunkInstanceT *)NULL)){
      continue;// This CI is not a starting point for a Scaffold.
    }

    if(BendCI != (ChunkInstanceT *)NULL){
      orientCI = A_B;
      edge = edgeB;
      neighbor = BendCI;
    }else if(AendCI != (ChunkInstanceT *)NULL){
      orientCI = B_A;
      edge = edgeA;
      neighbor = AendCI;
    }else{// Singleton Scaffold
      orientCI = A_B;
      edge = (CIEdgeT *)NULL;
      neighbor = (ChunkInstanceT *)NULL;
    }
    if(*currentScaffoldID == 404){
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

    if(orientCI == A_B){
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

    while(neighbor != (ChunkInstanceT *)NULL){
      ChunkOrientationType edgeOrient = GetEdgeOrientationWRT(edge,
                                                              thisCI->id);
      if(orientCI == A_B){
        if(edgeOrient == AB_AB){
          orientCI = A_B;
        }else if(edgeOrient == AB_BA){
          orientCI = B_A;
        }else{
          assert(0);
        }
      }else{// orientCI == B_A
        if(edgeOrient == BA_AB){
          orientCI = A_B;
        }else if(edgeOrient == BA_BA){
          orientCI = B_A;
        }else{
          assert(0);
        }
      }
      thisCI = neighbor;
      thisCI->scaffoldID = *currentScaffoldID;
      currentOffset.mean += edge->distance.mean;
      currentOffset.variance += edge->distance.variance;
      if(orientCI == A_B){
        aEndOffset = currentOffset;
      }else{// orientCI == B_A
        bEndOffset = currentOffset;
      }
      currentOffset.mean += thisCI->bpLength.mean;
      currentOffset.variance += thisCI->bpLength.variance;
      if(orientCI == A_B){
        bEndOffset = currentOffset;
      }else{// orientCI == B_A
        aEndOffset = currentOffset;
      }
      /* Don't Contig NOW!!! */
      InsertCIInScaffold(ScaffoldGraph, thisCI->id, *currentScaffoldID,
                         aEndOffset, bEndOffset,  TRUE /* Should be FALSE */, FALSE);
      if(orientCI == A_B){
        if(thisCI->essentialEdgeB != NULLINDEX){
          edge = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeB);
          assert(!isSloppyEdge(edge));
          neighbor = GetGraphNode(ScaffoldGraph->RezGraph,
                                  (edge->idA == thisCI->id) ? edge->idB : edge->idA);
        }else{// End of Scaffold
          edge = (CIEdgeT *)NULL;
          neighbor = (ChunkInstanceT *)NULL;
        }
      }else{// orientCI == B_A
        if(thisCI->essentialEdgeA != NULLINDEX){
          edge = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeA);
          assert(!isSloppyEdge(edge));
          neighbor = GetGraphNode(ScaffoldGraph->RezGraph,
                                  (edge->idA == thisCI->id) ? edge->idB : edge->idA);
        }else{// End of Scaffold
          edge = (CIEdgeT *)NULL;
          neighbor = (ChunkInstanceT *)NULL;
        }
      }
    }

    {/***** Check that scaffold is connected ****/
      CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, *currentScaffoldID);
      if (!IsScaffoldInternallyConnected(ScaffoldGraph,  scaffold, ALL_TRUSTED_EDGES)){
        fprintf(stderr,"* Scaffold " F_CID
                " is DISCONNECTED IMMEDIATELY AFTER INITIAL CONSTRUCTION!!!!\n",
                *currentScaffoldID);
        DumpACIScaffold(stderr,ScaffoldGraph, scaffold, FALSE);
        assert(0);
      }
    }
    (*currentScaffoldID)++;
  }
}


/**********************************************************************/
void MarkBifurcations(ScaffoldGraphT *graph,  int verbose){
  GraphNodeIterator nodes;
  ContigT *thisCI;


  InitGraphNodeIterator(&nodes, ScaffoldGraph->RezGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL){

    thisCI->scaffoldID = NULLINDEX;

    switch(thisCI->numEssentialA){
      case 0:
        break;
      case 1:
        if(thisCI->essentialEdgeA != NULLINDEX){
          CIEdgeT *edge = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeA);
          if(NeighborBranches(GetGraphNode(ScaffoldGraph->RezGraph, ((edge->idA == thisCI->id) ? edge->idB : edge->idA)),
                              GetEdgeOrientationWRT(edge, thisCI->id),
                              thisCI->essentialEdgeA)){
            thisCI->essentialEdgeA = NULLINDEX;
          }
        }
        break;
      default:
        thisCI->essentialEdgeA = NULLINDEX;
        break;
    }

    switch(thisCI->numEssentialB){
      case 0:
        break;
      case 1:
        if(thisCI->essentialEdgeB != NULLINDEX){
          CIEdgeT *edge = GetGraphEdge(ScaffoldGraph->RezGraph, thisCI->essentialEdgeB);
          if(NeighborBranches(GetGraphNode(ScaffoldGraph->RezGraph,
                                           ((edge->idA == thisCI->id) ?
					    edge->idB : edge->idA)),
                              GetEdgeOrientationWRT(edge, thisCI->id),
                              thisCI->essentialEdgeB)){
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

/**********************************************************************/
CDS_CID_t MarkShakyBifurcations(ScaffoldGraphT *graph,  int verbose){
  int numShaky = 0;
  GraphNodeIterator nodes;
  ContigT *contig;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((contig = NextGraphNodeIterator(&nodes)) != NULL){
    if (IsShakyContigAtScaffoldEnd(contig)) {
      SetNodeType(GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI),
                  UNRESOLVEDCHUNK_CGW);

      contig->flags.bits.isUnique = 0;
      contig->scaffoldID = NULLINDEX;
      contig->AEndNext = contig->BEndNext = NULLINDEX;

      numShaky++;

      if (verbose) {
        fprintf(stderr,"* Didn't include the following shaky contig  A end\n");
        DumpContig(stderr,ScaffoldGraph, contig, FALSE);
      }
    }
  }

  return numShaky;
}






/****************************************************************************/
void CreateScaffolds(ScaffoldGraphT *graph, int verbose){

  CDS_CID_t currentScaffoldID;

  /* Recycle the CIScaffolds VA */
  ResetNodeCGW_T(graph->ScaffoldGraph->nodes);

  MarkBifurcations(graph,  verbose);

  currentScaffoldID = 0;
  LabelCIScaffolds(&currentScaffoldID, verbose);

  /* We're looking for cycles where some of the Unique CIs have not been
     labeled as belonging to a scaffold. */

  DetectScaffoldCycles(&currentScaffoldID);

  // Adjust scaffold labeling, based on instrumenting
  // AdjustCIScaffoldLabels(graph, &currentScaffoldID);

  /* Actually insert CIs into scaffolds. */
  currentScaffoldID = 0;
  ActuallyInsertCIsIntoScaffolds(&currentScaffoldID,  verbose);

  graph->numLiveScaffolds = GetNumCIScaffoldTs(graph->CIScaffolds);
  assert(graph->numLiveScaffolds == currentScaffoldID);
}

void DeleteInferredEdges(ScaffoldGraphT *graph, int verbose){

  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;

  //  InitGraphNodeIterator(&nodes, ScaffoldGraph->RezGraph, GRAPH_NODE_UNIQUE_ONLY);
  // An inferred edge may exist between two shaky nodes that have
  // been demoted
  InitGraphNodeIterator(&nodes, ScaffoldGraph->RezGraph, GRAPH_NODE_DEFAULT);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL){
    GraphEdgeIterator edges;
    CIEdgeT *edge;

    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
    while((edge = NextGraphEdgeIterator(&edges)) != NULL){
      if(edge->flags.bits.isInferred){
	DeleteGraphEdge(ScaffoldGraph->RezGraph, edge);
      }
    }
  }
  return;
}

/************************************************************
 * AddScaffoldInferredEdges
 *
 *  This computes edges between the CIs in a scaffold, inferred
 *  by the relative positions of the CIs in the scaffold.
 *
 ***********************************************************/

void AddScaffoldInferredEdges(ScaffoldGraphT *graph,  int verbose){
  GraphNodeIterator scaffolds;
  CDS_CID_t sid;
  CIScaffoldT *scaffold;
  fprintf(stderr,"* AddScaffoldInferredEdges   scaffolds = %d\n",
	  (int) GetNumGraphNodes(graph->ScaffoldGraph));

  InitGraphNodeIterator(&scaffolds, graph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    CIScaffoldTIterator CIs;
    ChunkInstanceT *thisCI, *prevCI;
    sid = scaffold->id;

    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
    prevCI = NextCIScaffoldTIterator(&CIs);
    while((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL){
      LengthT inferredDistance;
      ChunkOrientationType inferredEdgeOrient;
      CDS_CID_t inferredEdgeIndex;
      CIEdgeT *inferredEdge;
      int isOverlapInferred = FALSE;

      assert(!thisCI->flags.bits.isDead &&
             !prevCI->flags.bits.isDead);

      if(GetNodeOrient(prevCI) == A_B){
	if(GetNodeOrient(thisCI) == A_B){
	  inferredEdgeOrient = AB_AB;
	  inferredDistance.mean = thisCI->offsetAEnd.mean -
	    prevCI->offsetBEnd.mean;
	  inferredDistance.variance = thisCI->offsetAEnd.variance -
	    prevCI->offsetBEnd.variance;
	}else{// GetNodeOrient(thisCI) == B_A
	  inferredEdgeOrient = AB_BA;
	  inferredDistance.mean = thisCI->offsetBEnd.mean -
	    prevCI->offsetBEnd.mean;
	  inferredDistance.variance = thisCI->offsetBEnd.variance -
	    prevCI->offsetBEnd.variance;
	}
      }else{// GetNodeOrient(prevCI) == B_A
	if(GetNodeOrient(thisCI) == A_B){
	  inferredEdgeOrient = BA_AB;
	  inferredDistance.mean = thisCI->offsetAEnd.mean -
	    prevCI->offsetAEnd.mean;
	  inferredDistance.variance = thisCI->offsetAEnd.variance -
	    prevCI->offsetAEnd.variance;
	}else{// GetNodeOrient(thisCI) == B_A
	  inferredEdgeOrient = BA_BA;
	  inferredDistance.mean = thisCI->offsetBEnd.mean -
	    prevCI->offsetAEnd.mean;
	  inferredDistance.variance = thisCI->offsetBEnd.variance -
	    prevCI->offsetAEnd.variance;
	}
      }

      if( inferredDistance.variance > 0.0 ){// SAK HACK!
	if(verbose)
	  fprintf(stderr,"* Adding inferred graph edge (" F_CID "," F_CID ",%c)\n",
                  prevCI->id, thisCI->id, inferredEdgeOrient);

        inferredEdgeIndex = AddGraphEdge(ScaffoldGraph->RezGraph, prevCI->id, thisCI->id,
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
        inferredEdge = GetGraphEdge(ScaffoldGraph->RezGraph, inferredEdgeIndex);
        inferredEdge->flags.bits.isInferred = TRUE;
        inferredEdge->flags.bits.isTentative = FALSE;
        inferredEdge->flags.bits.isActive = TRUE;
        inferredEdge->flags.bits.isConfirmed = TRUE;
        inferredEdge->flags.bits.isLeastSquares = TRUE;

#if 0
        fprintf(GlobalData->stderrc,"* Added inferred GraphEdge v1 (%d,%d,%c) isAContainsB:%d isBContainsA:%d fragA %d fragB %d isInferred %d isRaw %d\n",
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
	fprintf(stderr,"* Did NOT add inferred edge (" F_CID "," F_CID ",%c) offsets [%g,%g] [%g,%g]... a containment?\n",
		prevCI->id, thisCI->id, inferredEdgeOrient,
		thisCI->offsetAEnd.mean, thisCI->offsetBEnd.mean,
		prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean);
      }
      prevCI = thisCI;
    }
  }
  return;
}


/*****************************************************************************
 *  ResetEdgeStatus
 *
 ****************************************************************************/
void ResetEdgeStatus(ScaffoldGraphT *graph, int verbose){

  CDS_CID_t cid;
  ChunkInstanceT *thisCI;
  GraphNodeIterator nodes;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->RezGraph, GRAPH_NODE_UNIQUE_ONLY);
  while((thisCI = NextGraphNodeIterator(&nodes)) != NULL){
    GraphEdgeIterator edges;
    CIEdgeT *edge;

    if(verbose)
      fprintf(stderr,"*ResetEdgeStatus on node " F_CID "\n",   thisCI->id);

    cid = thisCI->id;

    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, ALL_END, ALL_EDGES,
			  GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
    while((edge = NextGraphEdgeIterator(&edges)) != NULL){
      ChunkInstanceT *otherCI = GetGraphNode(ScaffoldGraph->RezGraph,
					     (edge->idA == thisCI->id) ?
					     edge->idB : edge->idA);
      if((edge->idA != thisCI->id) && otherCI->flags.bits.isUnique){
	// Only reset canonical edges so as not to repeat the work.
	continue;
      }
      // Initialize the status flags for this edge
      InitGraphEdgeFlags(ScaffoldGraph->RezGraph, edge);
      if(verbose){
	PrintGraphEdge(GlobalData->stderrc, ScaffoldGraph->RezGraph, "r ", edge, thisCI->id);
      }
    }
  }
  return;
}


void MarkEssentialEdges(ScaffoldGraphT *graph,
                        int markShakyBifurcations,
                        int verbose){
  ResetEdgeStatus(graph,  verbose);
  AddScaffoldInferredEdges(graph,  verbose);
  MarkRedundantUniqueToUniqueEdges(graph,  verbose);
  MarkTwoHopConfirmedEdges(graph, verbose);
  MarkPathRemovedEdges(graph, verbose);
  SmoothWithInferredEdges(graph, markShakyBifurcations, verbose);
}



/*****************************************************************************
 * BuildUniqueCIScaffolds
 *****************************************************************************/
void BuildUniqueCIScaffolds(ScaffoldGraphT *graph,
                            int markShakyBifurcations,
                            int verbose){
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
  if(markShakyBifurcations){
    for(iteration = 0; todo == TRUE ; iteration++){
      int numShaky = 0;
      MarkEssentialEdges(graph,markShakyBifurcations,verbose);

      numShaky = MarkShakyBifurcations(graph, verbose);
      todo = (numShaky > 0);
      if(todo)
	DeleteInferredEdges(graph, verbose);
    }
  }
  DeleteInferredEdges(graph, verbose);
  MarkEssentialEdges(graph,FALSE,verbose);

  CreateScaffolds(graph,  verbose);

  DeleteInferredEdges(graph, verbose);

  //clearCacheSequenceDB(ScaffoldGraph->sequenceDB);
}


//#include "obsolete/saul_new_transitive_reduction"


/*****************************************************************************
 *  CompareReviveEdgesMean
 *
 ****************************************************************************/
static int CompareReviveEdgesMean(const void *c1, const void *c2){
  RevivedEdgeT *s1 = (RevivedEdgeT *)c1;
  RevivedEdgeT *s2 = (RevivedEdgeT *)c2;
  CDS_COORD_t mean1, mean2;

  mean1 = s1->minGap + s1->maxGap;
  mean2 = s2->minGap + s2->maxGap;
  if(mean1 > mean2){
    return((int)1);
  }else{
    return((int)-1);
  }
}
/*****************************************************************************
 *  CompareReviveEdgesID
 *
 ****************************************************************************/
static int CompareReviveEdgesID(const void *c1, const void *c2){
  RevivedEdgeT *s1 = (RevivedEdgeT *)c1;
  RevivedEdgeT *s2 = (RevivedEdgeT *)c2;

  if(s1->id > s2->id){
    return((int)1);
  }
  if(s1->id < s2->id){
    return((int)-1);
  }
  if((int)s1->orient > (int)s2->orient){
    return((int)1);
  }else{
    return((int)-1);
  }
}

/******************************************************************************
 *  CleanReviveEdges
 *
 *****************************************************************************/
void CleanReviveEdges(RevivedEdgeT *reviveEdgeStore,
		      RevivedEdgeT *endReviveEdge,
		      RevivedEdgeT **curReviveEdge){
  RevivedEdgeT *edgePtr, *packPtr;
  int numEdges = endReviveEdge - reviveEdgeStore;
  int numPackedEdges;
  CDS_COORD_t minGap, maxGap;

  qsort((void *)reviveEdgeStore, numEdges,
	sizeof(*reviveEdgeStore), CompareReviveEdgesID);
  for(edgePtr = packPtr = reviveEdgeStore, minGap = edgePtr->minGap,
	maxGap = edgePtr->maxGap; edgePtr < endReviveEdge; edgePtr++){
    if((edgePtr->id != packPtr->id) || (edgePtr->orient != packPtr->orient)){
      packPtr->minGap = minGap;
      packPtr->maxGap = maxGap;
      packPtr++;
      packPtr->id = edgePtr->id;
      packPtr->orient = edgePtr->orient;
      minGap = edgePtr->minGap;
      maxGap = edgePtr->maxGap;
    }else{
      if(edgePtr->minGap < minGap){
	minGap = edgePtr->minGap;
      }
      if(edgePtr->maxGap > maxGap){
	maxGap = edgePtr->maxGap;
      }
    }
  }
  packPtr->minGap = minGap;
  packPtr->maxGap = maxGap;
  packPtr++;
  numPackedEdges = packPtr - reviveEdgeStore;
  if(numPackedEdges < (numEdges / 2)){
    *curReviveEdge = reviveEdgeStore + numPackedEdges;
    return;
  }
  qsort((void *)reviveEdgeStore, numPackedEdges,
	sizeof(*reviveEdgeStore), CompareReviveEdgesMean);
  *curReviveEdge = reviveEdgeStore + (numEdges / 2);
  return;
}

/****************************************************************************
 *  ReviveTransitiveEdges
 *
 ****************************************************************************/
void ReviveTransitiveEdges(ScaffoldGraphT *graph,
			   LengthT *pathDistance,
			   ChunkOrientationType pathEdgeOrient,
			   ChunkInstanceT *thisCI,
			   RevivedEdgeT *reviveEdgeStore,
			   RevivedEdgeT *endReviveEdge,
			   RevivedEdgeT **curReviveEdge,
			   int end,
			   int verbose){

  CIEdgeT *transEdge;
  GraphEdgeIterator transEdges;

  if(verbose){
    fprintf(stderr,"*ReviveT thisCI " F_CID " distance (%f,%f) orient %c\n",
	    thisCI->id, pathDistance->mean, sqrt(pathDistance->variance),
	    (char)pathEdgeOrient);
  }

  if((pathEdgeOrient == AB_AB) ||
     (pathEdgeOrient == BA_AB)){
    /* transEdge must be a B_END edge */
    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &transEdges);// Use merged edges
  }else{
    /* transEdge must be a A_END edge */
    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &transEdges);// Use merged edges
  }
  while((transEdge = NextGraphEdgeIterator(&transEdges))!= NULL){
    ChunkInstanceT *nextCI;
    nextCI = GetGraphNode(ScaffoldGraph->RezGraph,
			  ((transEdge->idA == thisCI->id) ?
			   transEdge->idB : transEdge->idA));
    if(!isOverlapEdge(transEdge) || edgeContainsCI(transEdge, thisCI->id)){
      // Only interested in overlap edges and only if they are dovetail
      continue;
    }
    pathDistance->mean += thisCI->bpLength.mean +
      transEdge->distance.mean;
    pathDistance->variance += thisCI->bpLength.variance +
      ComputeFudgeVariance(transEdge->distance.mean);
    pathEdgeOrient = GetEdgeOrientationWRT(transEdge, thisCI->id);
    if(pathDistance->mean > - CGW_MISSED_OVERLAP){
      // Only revive edges longer than some minimum length
      continue;
    }
    if(*curReviveEdge == endReviveEdge){
      CleanReviveEdges(reviveEdgeStore, endReviveEdge, curReviveEdge);
    }
    //revive the edge
    (*curReviveEdge)->id = nextCI->id;
    if(end == A_END){
      if((pathEdgeOrient == AB_AB) ||
	 (pathEdgeOrient == BA_AB)){
	(*curReviveEdge)->orient = BA_AB;
      }else{
	(*curReviveEdge)->orient = BA_BA;
      }
    }else{
      if((pathEdgeOrient == AB_AB) ||
	 (pathEdgeOrient == BA_AB)){
	(*curReviveEdge)->orient = AB_AB;
      }else{
	(*curReviveEdge)->orient = AB_BA;
      }
    }
    (*curReviveEdge)->minGap =  (CDS_COORD_t)(pathDistance->mean -
                                              (3.0 * sqrt(pathDistance->variance)));
    (*curReviveEdge)->maxGap =  (CDS_COORD_t)(pathDistance->mean +
                                              (3.0 * sqrt(pathDistance->variance)));
    if(edgeContainsCI(transEdge, nextCI->id)){
      // Only interested in overlap edges and only if they are dovetail
      continue;
    }
    ReviveTransitiveEdges(graph, pathDistance, pathEdgeOrient,
			  nextCI, reviveEdgeStore, endReviveEdge,
			  curReviveEdge, end, verbose);
  }
  return;
}

/****************************************************************************
 *  ReviveEdges
 *
 ***************************************************************************/
int ReviveEdges(ScaffoldGraphT *graph, ChunkInstanceT *thisCI,
                int end, RevivedEdgeT *reviveEdgeStore, int maxEdges,
                int verbose){
  GraphEdgeIterator edges;
  CIEdgeT *edge;
  RevivedEdgeT *curReviveEdge, *endReviveEdge;
  curReviveEdge = reviveEdgeStore;
  endReviveEdge = reviveEdgeStore + maxEdges;
  InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, end, ALL_EDGES,
			GRAPH_EDGE_RAW_ONLY, &edges);// Use raw edges
  while((edge = NextGraphEdgeIterator(&edges))!= NULL){
    LengthT pathDistance;
    ChunkOrientationType pathEdgeOrient;
    ChunkInstanceT *nextCI;
    nextCI = GetGraphNode(ScaffoldGraph->RezGraph, (edge->idA == thisCI->id) ?
			  edge->idB : edge->idA);
    if(!isOverlapEdge(edge) || edgeContainsCI(edge, nextCI->id)){
      // Only interested in overlap edges and only if they are dovetail
      // or contain thisCI (NOT!!! contained by thisCI)
      continue;
    }
    pathEdgeOrient = GetEdgeOrientationWRT(edge, thisCI->id);
    pathDistance.mean = edge->distance.mean;
    pathDistance.variance = ComputeFudgeVariance(edge->distance.mean);
    ReviveTransitiveEdges(graph, &pathDistance, pathEdgeOrient,
			  nextCI, reviveEdgeStore, endReviveEdge,
			  &curReviveEdge, end, verbose);
  }
  CleanReviveEdges(reviveEdgeStore, endReviveEdge, &curReviveEdge);
  return(curReviveEdge - reviveEdgeStore);
}

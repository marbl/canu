
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
static char CM_ID[] = "$Id: CIScaffoldT_Cleanup_CGW.c,v 1.2 2004-09-23 20:25:19 mcschatz Exp $";

#define DEBUG 0
#undef DEBUG_DETAILED


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "UnionFind_AS.h"
#include "UtilsREZ.h"
#include "ChiSquareTest_CGW.h"
#include "MultiAlignment_CNS.h"
#include "DataTypesREZ.h"
#include "CommonREZ.h"
#include "Stats_CGW.h"   // for collecting scaffold merging stats
#include "FbacREZ.h"
#include <time.h>



typedef struct{
  ChunkInstanceT *firstCI;
  ChunkInstanceT *lastCI;
  CDS_CID_t firstCID;
  CDS_CID_t lastCID;
  LengthT minOffset;
  LengthT maxOffset;
  int32 count;
}ContigEndsT;

VA_DEF(ContigEndsT)

VA_TYPE(ContigEndsT) *ContigEnds = NULL;



  // Propagate Containment Overlaps
void  PropagateContainmentOverlapsToNewContig(ContigT *newContig,
                                              VA_TYPE(IntElementPos) *ContigPositions,
                                              CDS_COORD_t contigBase,
                                              int32 verbose){
  int32 i;

  for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions,i);
    ContigT *contig = GetGraphNode(ScaffoldGraph->ContigGraph,pos->ident);
    GraphEdgeIterator edges;
    EdgeCGW_T *edge;
    int flip =( pos->position.bgn > pos->position.end); // do we need to flip edge orientation
    CDS_COORD_t posBgn = pos->position.bgn - contigBase;
    CDS_COORD_t posEnd = pos->position.end - contigBase;

    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, contig->id, ALL_END, ALL_EDGES, 
			  GRAPH_EDGE_RAW_ONLY, 
			  &edges);
    if(verbose) {
      fprintf(GlobalData->stderrc,"* Scanning edges incident on contig " F_CID "\n", pos->ident);
    }

    while((edge = NextGraphEdgeIterator(&edges)) != NULL){
      ChunkOrientationType edgeOrient = GetEdgeOrientationWRT(edge,contig->id);
      ChunkOrientationType orient;
      LengthT distance;
      ContigT *otherContig = GetGraphNode(ScaffoldGraph->RezGraph, (edge->idA == contig->id? edge->idB:edge->idA));
      CDS_COORD_t overlap;


      if( !(edge->flags.bits.aContainsB ||         // only handling containment edges
	    edge->flags.bits.bContainsA) ||
	  otherContig->flags.bits.beingContigged ) // edges to other contigs that are being merged are not interesting
	continue;
          
      
      // We're only interested in edges where the contigs that are being merged CONTAIN
      // the other contig.
      if((otherContig->id == edge->idA && edge->flags.bits.aContainsB) ||
	 (otherContig->id == edge->idB && edge->flags.bits.bContainsA))
	continue;

      distance.variance = max(edge->distance.variance, 10.0);

      if(verbose){
	fprintf(GlobalData->stderrc,"* Contain:  pos " F_CID "  flip:%d edgeOrient:%c beg:" F_COORD " end:" F_COORD " newContigLength:%g\n",
		pos->ident, flip, edgeOrient, posBgn, posEnd, newContig->bpLength.mean);
	PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,"propogating ", edge, edge->idA);
      }

      //  Compute the two edges for the containment
      overlap = (CDS_COORD_t) -edge->distance.mean;

      if(!flip){
	switch(edgeOrient){
	case AB_AB:

	  //  A-------------------------newcontig-----------------------------------------B
	  //        A ======================contig============================B
	  //  -------------------------A ===========other==========B
          //                           |---------------overlap ---------------|
	  //                           |-----AB_AB overlap length-------------------------|  
	  //    
	  orient = AB_AB;
	  distance.mean = -(newContig->bpLength.mean - posEnd + overlap );
	  break;
	case BA_BA:
	  //  A-------------------------newcontig-----------------------------------------B
	  //        A ======================contig============================B
	  //  -------------------------A ===========other==========B
          //        |---------------overlap -----------------------|
	  //                             |-----AB_AB overlap length-----------------------|  
	  //    

	  orient = AB_AB;
	  distance.mean = -(newContig->bpLength.mean - posBgn - overlap + otherContig->bpLength.mean);
	  break;
	case AB_BA:
	  //  A-------------------------newcontig-----------------------------------------B
	  //        A ======================contig============================B
	  //  -------------------------B ===========other==========A
          //                           |---------------overlap ---------------|
	  //                           |-----AB_BA overlap length----------------------------|  
	  //    

	  orient = AB_BA;
	  distance.mean = -(newContig->bpLength.mean - posEnd +  overlap);
	  break;
	case BA_AB:
	  //  A-------------------------newcontig-----------------------------------------B
	  //        A ======================contig============================B
	  //  -------------------------B ===========other==========A
          //        |-------------------overlap -------------------|
	  //                           |---AB_BA overlap length----------------------------|  
	  //    

	  orient = AB_BA;
	  distance.mean = -(newContig->bpLength.mean - posBgn - overlap + otherContig->bpLength.mean);
	  break;
	default:
	  assert(0);
	}
      }else{
	switch(edgeOrient){
	case BA_AB:
	  //  A-------------------------newcontig-----------------------------------------B
	  //        B ======================contig============================A
	  //  -------------------------A ===========other==========B
          //                           |-------------------overlap -----------|
	  //                           |-----AB_AB overlap length-------------------------|  
	  //    
	  //    
	  orient = AB_AB;
	  distance.mean = -(newContig->bpLength.mean - posBgn + overlap);
	  break;
	case AB_BA:
	  //  A-------------------------newcontig-----------------------------------------B
	  //        B ======================contig============================A
	  //  -------------------------A ===========other==========B
          //        |-------------------overlap -------------------|
	  //                           |-----AB_AB overlap length----------------------------|  
	  //    
	  //    
	  orient = AB_AB;
	  distance.mean = -(newContig->bpLength.mean - posEnd - overlap + otherContig->bpLength.mean);
	  break;
	case BA_BA:
	  //  A-------------------------newcontig-----------------------------------------B
	  //        B ======================contig============================A
	  //  -------------------------B ===========other==========A
          //                           |-------------------overlap -----------|
	  //                           |-----AB_BA overlap length-------------------------|  
	  //    
	  //    
	  orient = AB_BA;
	  distance.mean = -(newContig->bpLength.mean - posBgn  +  overlap);
	  break;
	case AB_AB:
	  //  A-------------------------newcontig-----------------------------------------B
	  //        B ======================contig============================A
	  //  -------------------------B ===========other==========A
          //         |-------------------overlap ------------------|
	  //                           |-----AB_BA overlap length-------------------------|  
	  //    
	  //    
	  orient = AB_BA;
	  distance.mean = -(newContig->bpLength.mean - posEnd - overlap + otherContig->bpLength.mean);
	  break;
	default:
	  assert(0);
	}
      }
      // Add the edge to the graph and insert them 
      {
	EdgeCGW_T *newEdge;
	CDS_CID_t eid;
	int flip = (newContig->id > otherContig->id);
	newEdge = GetFreeGraphEdge(ScaffoldGraph->RezGraph);
	eid = GetVAIndex_EdgeCGW_T(ScaffoldGraph->RezGraph->edges, newEdge);
	
	SetGraphEdgeStatus(ScaffoldGraph->RezGraph, newEdge, UNKNOWN_EDGE_STATUS);
	if(!flip){
	  newEdge->idA = newContig->id;
	  newEdge->idB = otherContig->id;
	  newEdge->orient = orient;
	  newEdge->flags.bits.aContainsB = TRUE;	
	  newEdge->flags.bits.hasContainmentOverlap = TRUE;
	}else{
	  newEdge->idA = otherContig->id;
	  newEdge->idB = newContig->id;
	  newEdge->orient = FlipEdgeOrient(orient);
	  newEdge->flags.bits.bContainsA = TRUE;	
	  newEdge->flags.bits.hasContainmentOverlap = TRUE;
	}
	newEdge->distance = distance;
	newEdge->edgesContributing = 1;

	{
	  ChunkOverlapCheckT olap;
	  if(LookupOverlap(ScaffoldGraph->RezGraph, newEdge->idA, newEdge->idB, newEdge->orient, &olap)){
	    FreeGraphEdgeByEID(ScaffoldGraph->RezGraph, eid);
	    continue;
	  }
	}

	{
	  CDS_COORD_t delta =
            (CDS_COORD_t) min((3 * sqrt(newEdge->distance.variance)), 50);

	  ChunkOverlapCheckT olap =
	    OverlapChunks(ScaffoldGraph->RezGraph,
                          newEdge->idA, newEdge->idB,
                          newEdge->orient, // handle suspicious
			  -newEdge->distance.mean - delta, 
			  -newEdge->distance.mean + delta, 
			  CGW_DP_ERATE, FALSE);
          
	  if(olap.overlap == 0 || olap.suspicious){
	    if(olap.suspicious){
	      fprintf(GlobalData->stderrc,"* CIS_C0 SUSPICIOUS Overlap found! Looked for (" F_CID "," F_CID ",%c) [%d,%d] found (" F_CID "," F_CID ",%c) " F_COORD "\n",
		      newEdge->idA, newEdge->idB, newEdge->orient, 
		      (int)(-newEdge->distance.mean - delta),
		      (int) (-newEdge->distance.mean + delta),
		      olap.spec.cidA, olap.spec.cidB, olap.spec.orientation,
		      olap.overlap);
	    }
	    if(verbose) {
	      fprintf(GlobalData->stderrc,"* NO Overlap found!\n");
              PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,"*F* ",
                             newEdge, newEdge->idA);
            }
	    FreeGraphEdgeByEID(ScaffoldGraph->RezGraph,eid);
	    continue;
	  }
	}
	if(verbose) {
	  PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,"*C* ",
                         newEdge, newEdge->idA);
        }
	InsertGraphEdge(ScaffoldGraph->RezGraph, eid, FALSE);

	// Create an appropriate hash table entry
	CreateChunkOverlapFromEdge(ScaffoldGraph->RezGraph, newEdge, FALSE);

      }
    }
  }
}


static int CompareIntElementPosByBgnPos(const void * c1, const void * c2)
{
  IntElementPos * i1 = (IntElementPos *) c1;
  IntElementPos * i2 = (IntElementPos *) c2;

  return(min(i1->position.bgn, i1->position.end) -
         min(i2->position.bgn, i2->position.end));
}

// Sort routine to order edges that are not yet inserted in graph
//
static int CompareEdges(const void *c1, const void *c2){
  CDS_CID_t eid1 = *(CDS_CID_t *)c1;
  CDS_CID_t eid2 = *(CDS_CID_t *)c2;
  EdgeCGW_T *edge1 = GetGraphEdge(ScaffoldGraph->RezGraph, eid1);
  EdgeCGW_T *edge2 = GetGraphEdge(ScaffoldGraph->RezGraph, eid2);
  int diff;
  diff = edge1->idA - edge2->idA;
  if(diff)
    return diff;
  diff = edge1->idB - edge2->idB;
  if(diff)
    return diff;

  diff = edge1->orient - edge2->orient;
  if(diff)
    return diff;

  return 0;
}

// Collection of indices of raw edges reflecting potential overlaps to
// propagate new contig. 
static VA_TYPE(CDS_CID_t) *CollectedEdges = NULL;

// Collection of indices of raw edges that we've checked
// and discarded
static VA_TYPE(CDS_CID_t) *DiscardedEdges = NULL;

// This VA is used to hold indices of tentative edges between the
// new contig and a particular potential overlap neighbor.  Used in
// the interface to MergeEdges.
static VA_TYPE(CDS_CID_t) *TentativeEdges = NULL;

// This VA is used to hold indices of confirmed edges that we
// need to add to the graph.  We actually discard these merged
// edges and insert corresponding raw edge into the graph
static VA_TYPE(CDS_CID_t) *MergedEdges = NULL;


// Propagate overlap edges (other than extremal or containment) incident on extremal contigs to new contig
// These overlaps become implied containments.  Since containments can be viewed as an edge on either the A or
// B end of the containing/contained elements, we generate a pair of edges for each containment.
//
void PropagateInternalOverlapsToNewContig(ContigT *newContig,
                                          VA_TYPE(IntElementPos) *ContigPositions,
                                          CDS_CID_t scaffoldID,
                                          CDS_COORD_t contigBase,
                                          int32 verbose){
  int32 i;
  int32 cnt = 0;

  if(TentativeEdges == NULL){
    TentativeEdges = CreateVA_CDS_CID_t(100);
    CollectedEdges = CreateVA_CDS_CID_t(100);
    DiscardedEdges = CreateVA_CDS_CID_t(100);
    MergedEdges = CreateVA_CDS_CID_t(100);
  }

  ResetVA_CDS_CID_t(TentativeEdges);
  ResetVA_CDS_CID_t(CollectedEdges);
  ResetVA_CDS_CID_t(DiscardedEdges);
  ResetVA_CDS_CID_t(MergedEdges);

  // First,we map all of the non-containment overlaps of the contigs that were merged onto the
  // new contig.  We create an appropriate set of raw edges

  if(verbose) {
    fprintf(GlobalData->stderrc,"* PropagateInternalOverlaps \n");
  }

  for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions,i);
    ContigT *contig = GetGraphNode(ScaffoldGraph->ContigGraph,pos->ident);
    GraphEdgeIterator edges;
    EdgeCGW_T *edge;
    int flip =( pos->position.bgn > pos->position.end); // do we need to flip edge orientation
    CDS_COORD_t posBgn = pos->position.bgn - contigBase;
    CDS_COORD_t posEnd = pos->position.end - contigBase;

    /* We exclude containment since we handle it elsewhere */
    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, contig->id, ALL_END, ALL_EDGES, 
			  GRAPH_EDGE_RAW_ONLY|GRAPH_EDGE_EXCLUDE_CONTAINMENT, 
			  &edges);
    if(verbose) {
      fprintf(GlobalData->stderrc,"* Scanning edges incident on contig " F_CID "\n", pos->ident);
    }

    while((edge = NextGraphEdgeIterator(&edges)) != NULL){
      ChunkOrientationType edgeOrient = GetEdgeOrientationWRT(edge,contig->id);
      ChunkOrientationType orient1, orient2;
      LengthT distance1, distance2;
      ContigT *otherContig = GetGraphNode(ScaffoldGraph->RezGraph, (edge->idA == contig->id? edge->idB:edge->idA));
      CDS_COORD_t overlap;
      double maxDistance;

      distance1.variance = max(edge->distance.variance, 10.0);
      distance2.variance = distance1.variance;
      if(!isOverlapEdge(edge) ||
	 edge->flags.bits.aContainsB ||
	 edge->flags.bits.bContainsA ||
	 otherContig->flags.bits.beingContigged )
	continue;
          
      
      if(verbose){
	fprintf
          (GlobalData->stderrc,
           "* pos " F_CID "  flip:%d edgeOrient:%c beg:" F_COORD " end:" F_COORD " newContigLength:%g\n",
           pos->ident, flip, edgeOrient, posBgn, posEnd,
           newContig->bpLength.mean);
	PrintGraphEdge(GlobalData->stderrc,
                       ScaffoldGraph->RezGraph,"propogating*I* ",
                       edge, edge->idA);
      }

      //  Compute the two edges for the containment
      // NOTE -- not all of these will be containments!!!!

      overlap = -edge->distance.mean;

      if(!flip){
	switch(edgeOrient){
	case AB_AB:

	  //  A-------------------------newcontig-----------------------------------------B
	  //        A ===========contig============B
	  //  -------------------------A ===========other==========B
	  //                           |-----AB_AB overlap length------------------------|  
	  //  |----------------BA_BA overlap length-----------|  
	  //    
	  orient1 = AB_AB;
	  maxDistance = (newContig->bpLength.mean - posEnd + overlap);
	  distance1.mean = -maxDistance;
	  orient2 = BA_BA;
	  distance2.mean = -(posEnd - overlap + min(otherContig->bpLength.mean, maxDistance));
	  break;
	case BA_BA:
	  //  A-------------------------newcontig-----------------------------------------B
	  //                                  A ===========contig============B
	  //  --------------------A ===========other==========B
	  //  |----------------BA_BA overlap length-----------|  
	  //                             |-----------------AB_AB overlap length-----------|  
	  //    
	  orient1 = BA_BA;
	  maxDistance = posBgn + overlap;
	  distance1.mean = -maxDistance;
	  orient2 = AB_AB;
	  distance2.mean = -(newContig->bpLength.mean - posBgn + min(otherContig->bpLength.mean,maxDistance) - overlap);
	  break;
	case AB_BA:
	  //  A-------------------------newcontig-----------------------------------------B
	  //        A ===========contig============B
	  //  -------------------------B ===========other==========A
	  //  |---------------------BA_AB overlap lengtnh-----------|  
	  //                             |-----------------AB_BA overlap length-----------|  
	  //    
	  orient1 = AB_BA;
	  maxDistance = (newContig->bpLength.mean - posEnd + overlap);
	  distance1.mean = -maxDistance;
	  orient2 = BA_AB;
	  distance2.mean = -(posEnd + min(maxDistance,otherContig->bpLength.mean) - overlap);
	  break;
	case BA_AB:
	  //  A-------------------------newcontig-----------------------------------------B
	  //                                    A ===========contig============B
	  //  ----------------B ===========other==========A
	  //  |------------BA_AB overlap length-----------|  
	  //                  |----------------AB_BA overlap length----------------------|  
	  //    
	  orient1 = BA_AB;
	  maxDistance = (posBgn + overlap);
	  distance1.mean = -maxDistance;
	  orient2 = AB_BA;
	  distance2.mean = -(newContig->bpLength.mean - posBgn + min(maxDistance,otherContig->bpLength.mean) - overlap);
	  break;
	default:
	  assert(0);
	}
      }else{
	switch(edgeOrient){
	case BA_AB:

	  //  A-------------------------newcontig-----------------------------------------B
	  //        B ===========contig============A
	  //  -------------------------A ===========other==========B
	  //  |-----BA_BA overlap length---------------------------|  
	  //                             |----------------AB_AB overlap length-----------|  
	  //    
	  orient1 = BA_BA;
	  maxDistance = (posBgn + otherContig->bpLength.mean - overlap);
	  distance1.mean = -maxDistance;
	  orient2 = AB_AB;
	  distance2.mean = -(newContig->bpLength.mean - posBgn + overlap);
	  break;
	case AB_BA:
	  //  A-------------------------newcontig-----------------------------------------B
	  //                                  B ===========contig============A
	  //  --------------------A ===========other==========B
	  //  |----------------BA_BA overlap length-----------|  
	  //                             |-----------------AB_AB overlap length-----------|  
	  //    
	  orient1 = AB_AB;
	  maxDistance = (posEnd + overlap);
	  distance1.mean = -(newContig->bpLength.mean - posEnd + min(maxDistance,otherContig->bpLength.mean) - overlap);
	  orient2 = BA_BA;
	  distance2.mean = -maxDistance;
	  break;
	case BA_BA:
	  //  A-------------------------newcontig-----------------------------------------B
	  //        B ===========contig============A
	  //  -------------------------B ===========other==========A
	  //  |-----BA_AB overlap length---------------------------|  
	  //                             |-----------------AB_BA overlap length-----------|  
	  //    
	  orient1 = AB_BA;
	  maxDistance = (newContig->bpLength.mean - posBgn + overlap);
	  distance1.mean = -maxDistance;
	  orient2 = BA_AB;
	  distance2.mean = -(posBgn + min(maxDistance,otherContig->bpLength.mean) - overlap);
	  break;
	case AB_AB:
	  //  A-------------------------newcontig-----------------------------------------B
	  //                                    B ===========contig============A
	  //  ----------------B ===========other==========A
	  //  |------------BA_AB overlap length-----------|  
	  //                    |--------------------------------AB_BA overlap length-----|  
	  //    
	  maxDistance = (posEnd + overlap);
	  orient1 = BA_AB;
	  distance1.mean = -maxDistance;
	  orient2 = AB_BA;
     	  distance2.mean = -(newContig->bpLength.mean - posEnd + min(maxDistance,otherContig->bpLength.mean) - overlap);
	  break;
	default:
	  assert(0);
	}
      }
      // Add them to the graph, do NOT insert them.  Keep a list of such edges in collectedEdges
      {
	ChunkOverlapCheckT olap;
	// If we don't find either of these edges already in the graph, proceed
	if(!LookupOverlap(ScaffoldGraph->RezGraph, newContig->id, (edge->idA == contig->id? edge->idB: edge->idA), orient1, &olap) &&
	   !LookupOverlap(ScaffoldGraph->RezGraph, newContig->id, (edge->idA == contig->id? edge->idB: edge->idA), orient2, &olap))
	  {
	    EdgeCGW_T *newEdge;
	    CDS_CID_t eid;
	    int32  isContainment = (min(-distance1.mean,-distance2.mean) > otherContig->bpLength.mean);

            CDS_CID_t edgeID = 
                   GetVAIndex_EdgeCGW_T(ScaffoldGraph->RezGraph->edges, edge);
	    newEdge = GetFreeGraphEdge(ScaffoldGraph->RezGraph);
	    eid = GetVAIndex_EdgeCGW_T(ScaffoldGraph->RezGraph->edges, newEdge);
	    edge = GetGraphEdge(ScaffoldGraph->RezGraph, edgeID);
	    AppendCDS_CID_t(CollectedEdges, &eid);

	    newEdge->idA = newContig->id;
	    newEdge->idB = (edge->idA == contig->id? edge->idB: edge->idA);
	    newEdge->orient = orient1;
	    newEdge->distance = distance1;
	    newEdge->edgesContributing = cnt++;
	    
	    newEdge->flags.bits.bContainsA = FALSE;
	    newEdge->flags.bits.aContainsB = isContainment;
	    newEdge->flags.bits.hasContainmentOverlap = isContainment;
	    newEdge->flags.bits.hasContributingOverlap = !isContainment;
	
	    if(verbose) {
	      PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,
                             "    collecting*I* ", newEdge, newEdge->idA);
            }

	    if(isContainment){
	      newEdge = GetFreeGraphEdge(ScaffoldGraph->RezGraph);
	      eid = GetVAIndex_EdgeCGW_T(ScaffoldGraph->RezGraph->edges, newEdge);
	      AppendCDS_CID_t(CollectedEdges, &eid);
	      newEdge->idA = newContig->id;
	      newEdge->idB = (edge->idA == contig->id? edge->idB: edge->idA);
	      newEdge->edgesContributing = cnt++;

	      newEdge->orient = orient2;
	      newEdge->distance = distance2;
	      newEdge->flags.bits.bContainsA = FALSE;
	      newEdge->flags.bits.aContainsB = TRUE;
	      newEdge->flags.bits.hasContainmentOverlap = TRUE;

	      if(verbose) {
		PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,
                               "    collecting*I* ", newEdge, newEdge->idA);
              }
		assert(!newEdge->flags.bits.isDeleted);
	    }else{
	      if(verbose){
		fprintf(GlobalData->stderrc,"* ***NOT A CONTAINMENT!!!!\n");
	      }
	    }
	  }
      }
    }
  }
  {
    CDS_CID_t *collectedEdges = GetCDS_CID_t(CollectedEdges,0);
    EdgeCGW_T *e;
    int32 i;
  // Sort the edge indices, so they group together multiple edges between the new contig and other contigs
      qsort((void *)collectedEdges, GetNumCDS_CID_ts(CollectedEdges), sizeof(CDS_CID_t), CompareEdges);
      if(verbose) {
	fprintf(GlobalData->stderrc,"* Collected the following edges *\n");
      }
      for(i = 0; i < GetNumCDS_CID_ts(CollectedEdges); i++){
	e = GetGraphEdge(ScaffoldGraph->RezGraph, *GetCDS_CID_t(CollectedEdges,i));
	if(verbose) {
	  PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph," ", e, e->idA);
        }
      }
  }
  // Then we try to merge the edges incident on newContig.  Any that successfully merge are CANDIDATES
  // for creating a raw overlap edge.  We screen with dpalign first.


  {
    int32 i;
    int32 numEdges = GetNumCDS_CID_ts(CollectedEdges);
    CDS_CID_t currID = NULLINDEX;
    ChunkOrientationType currOrient = XX_XX;
    ResetVA_CDS_CID_t(TentativeEdges);

    for(i = 0; i < numEdges; ){
      CDS_CID_t eid = *GetCDS_CID_t(CollectedEdges,i);
      EdgeCGW_T *e = GetGraphEdge(ScaffoldGraph->RezGraph, eid);
      int merge = FALSE;
      CDS_CID_t last = i == numEdges -1;
      if(currID == NULLINDEX){
	currID = e->idB;
	currOrient = e->orient;	
	e->flags.all = 0;
	AppendCDS_CID_t(TentativeEdges,&eid);
	i++;
      }else  if(e->idB == currID &&
	 e->orient == currOrient){
	e->flags.all = 0;
	AppendCDS_CID_t(TentativeEdges,&eid);
	i++;
      }else{
	merge = TRUE;
      }

      if( last && GetNumCDS_CID_ts(TentativeEdges) > 1)
				merge = TRUE;
	
      if(merge){
	int addedEdges=-2;
        //The return value of MergeGraphEdges() is the number of edges
        //that this routine appended to the edges array - returns 0
        //when no merging was done and -1 on failure. We will use -2
        //to denote that the routine was never called.

	EdgeCGW_T *ee;
	int i;
	if(GetNumCDS_CID_ts(TentativeEdges) == 1){
	  // keep track of which edges to discard
#ifdef DEBUG_DETAILED
          fprintf(GlobalData->stderrc,"* Edge " F_CID " to Discarded edges\n",
                  *GetCDS_CID_t(TentativeEdges,0));
#endif          
          AppendCDS_CID_t(DiscardedEdges, GetCDS_CID_t(TentativeEdges,0));
	}else{
	  if(GlobalData->debugLevel > 0)
	    fprintf(GlobalData->stderrc,"* Merging " F_CID "," F_CID ",%c   %d edges\n",
		  e->idA, currID, currOrient,
                    (int) GetNumCDS_CID_ts(TentativeEdges));
	  addedEdges = MergeGraphEdges(ScaffoldGraph->RezGraph, TentativeEdges);
	  if(addedEdges> 0){
	    if(GlobalData->debugLevel > 0)
	      fprintf(GlobalData->stderrc,"* Merging (" F_CID "," F_CID ",%c): Added %d edges!!!\n",
		    e->idA, currID, currOrient, addedEdges);

	    ee = GetGraphEdge(ScaffoldGraph->RezGraph, *GetCDS_CID_t(TentativeEdges,0));
	    ee->flags.bits.bContainsA = FALSE;
	    ee->flags.bits.aContainsB = TRUE;
	    ee->flags.bits.hasContainmentOverlap = TRUE;

	    for(i = 0; i < addedEdges; i++){
#ifdef DEBUG_DETAILED
	      //	      if(GlobalData->debugLevel > 0)
		PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,
                               "* Successful merge: ", ee, ee->idA);
#endif                
	      AppendCDS_CID_t(MergedEdges, GetCDS_CID_t(TentativeEdges,i));
	    }
	  }
#ifdef DEBUG_DETAILED
          fprintf(GlobalData->stderrc,"* addedEdges %d total:%d\n",
                  addedEdges, (int) GetNumCDS_CID_ts(TentativeEdges));
#endif          
	}
		
	if(addedEdges > 0)
	for(i = addedEdges; i < GetNumCDS_CID_ts(TentativeEdges) ; i++){
	  // keep track of which edges to discard
	    ee = GetGraphEdge(ScaffoldGraph->RezGraph, *GetCDS_CID_t(TentativeEdges,i));
	    ee->flags.bits.bContainsA = FALSE;
	    ee->flags.bits.aContainsB = TRUE;
	    ee->flags.bits.hasContainmentOverlap = TRUE;

	    PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,"To Discarded Edges", ee, ee->idA);
	  AppendCDS_CID_t(DiscardedEdges, GetCDS_CID_t(TentativeEdges,i));
	}
	currID = NULLINDEX;
	currOrient = XX_XX;
	ResetVA_CDS_CID_t(TentativeEdges);
	if(last)
	  break;
      }

    }


    // Now Delete All of these edges, and insert raw overlap edges for the survivors
    // that also pass dp_align
    {
      int32 i;
      EdgeCGW_T *newEdge;
      if(verbose) {
	fprintf(GlobalData->stderrc,"* Discarded the following DISCARDED edges:\n");
      }
      for(i = 0; i < GetNumCDS_CID_ts(DiscardedEdges); i++){
	CDS_CID_t eid;
	EdgeCGW_T *e;

	eid = *GetCDS_CID_t(DiscardedEdges,i);
	e = GetGraphEdge(ScaffoldGraph->RezGraph, eid);

#ifdef DEBUG
	fprintf(GlobalData->stderrc,"* i = %d eid = " F_CID "\n", i,eid);
#endif        

	if(e->flags.bits.isDeleted){
	  fprintf(GlobalData->stderrc,"* i = %d Edge " F_CID " is ALREADY deleted...WIERD!\n", i,eid);
	}else{
	  if(verbose) {
	    PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph," ", e, e->idA);
            fprintf(GlobalData->stderrc,"* i = %d Edge " F_CID "!\n", i,eid);
          }
	  FreeGraphEdgeByEID(ScaffoldGraph->RezGraph,eid);
	}
      }


      if(verbose) {
	fprintf(GlobalData->stderrc,"* Inserting the following MERGED edges:\n");
      }

      for(i = 0; i < GetNumCDS_CID_ts(MergedEdges); i++){
	CDS_CID_t eid= *GetCDS_CID_t(MergedEdges,i);
	CDS_CID_t neid;
	EdgeCGW_T *e = GetGraphEdge(ScaffoldGraph->RezGraph, eid);
	ChunkOverlapCheckT olap;
	CDS_COORD_t delta = min((3 * sqrt(e->distance.variance)), 50);
	CDS_CID_t idA, idB;
	ChunkOrientationType orient;
	if(verbose) {
	  PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,"dp_align check ",
                         e, e->idA);
        }

	idA = e->idA;
	idB = e->idB;
	orient = e->orient;
	if(e->idA > e->idB){
	  idA = e->idB;
	  idB = e->idA;
	  orient = FlipEdgeOrient(e->orient);
	}

	
	{
	  ChunkOverlapCheckT olap;
	  if(LookupOverlap(ScaffoldGraph->RezGraph, e->idA, e->idB, e->orient, &olap)){
	    if(verbose) {
	      fprintf(GlobalData->stderrc,"* Overlap found already...skipping\n");
            }
	    FreeGraphEdgeByEID(ScaffoldGraph->RezGraph,eid);
	    continue;
	  }
	}


	olap = OverlapChunks(ScaffoldGraph->RezGraph,
                             idA, idB, orient, // handles suspicious
			     -e->distance.mean - delta, 
			     -e->distance.mean + delta, 
			     CGW_DP_ERATE, FALSE);
        


	if(olap.overlap == 0 || olap.suspicious){
	    if(olap.suspicious){
	      fprintf(GlobalData->stderrc,"* CIS_C1 SUSPICIOUS Overlap found! Looked for (" F_CID "," F_CID ",%c) [%d,%d] found (" F_CID "," F_CID ",%c) " F_COORD " \n",
		      idA, idB, orient, 
		      (int)(-e->distance.mean - delta),
		      (int)(-e->distance.mean + delta),
		      olap.spec.cidA, olap.spec.cidB, olap.spec.orientation, olap.overlap);
	    }
	  
	  if(verbose) {
	    fprintf(GlobalData->stderrc,"* NO Overlap found!\n");
          }
	  FreeGraphEdgeByEID(ScaffoldGraph->RezGraph,eid);
	  continue;
	}

	if(verbose) {
	  fprintf(GlobalData->stderrc,"* Looked for [%g,%g] found " F_COORD " overlap\n",
			     e->distance.mean - delta, 
			     e->distance.mean + delta, 
		             olap.overlap);
        }
        
	newEdge = GetFreeGraphEdge(ScaffoldGraph->RezGraph);
	neid = GetVAIndex_EdgeCGW_T(ScaffoldGraph->RezGraph->edges, newEdge);


	newEdge->flags.bits.aContainsB = olap.BContainsA || olap.AContainsB;
	newEdge->flags.bits.bContainsA = FALSE;
	newEdge->flags.bits.hasContainmentOverlap = olap.BContainsA || olap.AContainsB;
	newEdge->flags.bits.hasContributingOverlap = !(olap.BContainsA || olap.AContainsB);

	newEdge->idA = e->idA;
	newEdge->idB = e->idB;
	newEdge->orient = e->orient;
	newEdge->distance.mean = e->distance.mean;
	newEdge->distance.variance = ComputeFudgeVariance(newEdge->distance.mean);

	if(newEdge->idA > newEdge->idB){
	  CDS_CID_t tmp = newEdge->idA;
	  newEdge->flags.bits.bContainsA = newEdge->flags.bits.aContainsB;
	  newEdge->flags.bits.aContainsB = FALSE;
	  newEdge->flags.bits.hasContainmentOverlap = olap.BContainsA || olap.AContainsB;
	  newEdge->flags.bits.hasContributingOverlap = !(olap.BContainsA || olap.AContainsB);
	  newEdge->idA = newEdge->idB;
	  newEdge->idB = tmp;
	  newEdge->orient =  FlipEdgeOrient(newEdge->orient);
	}
	newEdge->edgesContributing = 1;
	
	if(verbose) {
	  PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,
                         "Inserting #I# ", newEdge, newEdge->idA);
        }

	InsertGraphEdge(ScaffoldGraph->RezGraph, neid, FALSE);

	// Create an appropriate hash table entry
	CreateChunkOverlapFromEdge(ScaffoldGraph->RezGraph, newEdge, FALSE);

	
	FreeGraphEdgeByEID(ScaffoldGraph->RezGraph,eid);
      }
    }
  }
  RecycleDeletedGraphElements(ScaffoldGraph->RezGraph);
}


// Propagate overlap edges incident on extremal contigs to new contig
void PropagateExtremalOverlapsToNewContig(CDS_CID_t contigID, int contigEnd, ContigT *newContig, int newContigEnd, CDS_COORD_t contigBase, int32 verbose){
  GraphEdgeIterator edges;
  EdgeCGW_T *edge;
  int first = 0;
  CDS_CID_t rawEdgeID;

  if(verbose){
    fprintf(GlobalData->stderrc,"* PropagateExtremalOverlaps from (" F_CID ",%c) to (" F_CID ",%c)\n",
	  contigID, (contigEnd == A_END?'A':'B'),
	  newContig->id, (newContigEnd == A_END?'A':'B'));

    DumpContig(GlobalData->stderrc,ScaffoldGraph,
               GetGraphNode(ScaffoldGraph->RezGraph, contigID),FALSE);
  }


  InitGraphEdgeIterator(ScaffoldGraph->RezGraph, contigID, contigEnd, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
  while((edge = NextGraphEdgeIterator(&edges)) != NULL){
    EdgeCGW_T *rawEdge, *newEdge;
    ChunkOrientationType newOrient = edge->orient;
    int flip;
    CDS_CID_t otherContig;
    NodeCGW_T *otherCI = GetGraphNode(ScaffoldGraph->RezGraph, (edge->idA == contigID? edge->idB:edge->idA));
    
    if(!isOverlapEdge(edge))
      continue;

    // Handle containment edges elsewhere
    if(edge->flags.bits.aContainsB ||
       edge->flags.bits.bContainsA ||
       otherCI->flags.bits.beingContigged)
      continue;

    rawEdge = edge;
    if(!rawEdge->flags.bits.isRaw){
      rawEdge = GetGraphEdge(ScaffoldGraph->RezGraph, rawEdge->nextRawEdge);
      AssertPtr(rawEdge);
      while(!isOverlapEdge(rawEdge) && rawEdge->nextRawEdge != NULLINDEX)
	rawEdge = GetGraphEdge(ScaffoldGraph->RezGraph, rawEdge->nextRawEdge);
	   
    }
    //    PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph, "Propogating: ", rawEdge,contigID);
    //    fflush(GlobalData->stderrc);
    AssertPtr(rawEdge);

    /*****  rawEdge should be propagated to be incident on the new contig!!! ****/
    if(verbose){
      if(first++ == 0)
	fprintf(GlobalData->stderrc,"* Need to propagate the following edges from (" F_CID ",%c) to the new contig's %c End\n",
	      contigID, 
	      (contigEnd == A_END?'A':'B'),
	      (newContigEnd == A_END?'A':'B')  );
      PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph, "X Prop *R* ", rawEdge, rawEdge->idA);
    }else{
      first++;
    }
    //UnlinkGraphEdge(ScaffoldGraph->RezGraph, rawEdge); // Unlink it from the graph

    rawEdgeID = GetVAIndex_EdgeCGW_T(ScaffoldGraph->RezGraph->edges, rawEdge);
    newEdge = GetFreeGraphEdge(ScaffoldGraph->RezGraph);

    rawEdge = GetGraphEdge(ScaffoldGraph->RezGraph, rawEdgeID);
    newEdge->flags = rawEdge->flags;
    newEdge->idA = rawEdge->idA;
    newEdge->idB = rawEdge->idB;
    newEdge->orient = rawEdge->orient;
    newEdge->distance = rawEdge->distance;

    {
      // Compute the new edge orientation
      // First get edge orientation into canonical form
      flip = !(rawEdge->idA == contigID);
      otherContig = otherCI->id;

#if 0
      // Get the orientation into canonical form, relative to the contig we're replacing
      if(flip){
	newOrient = FlipEdgeOrient(newOrient);
      }
#endif
      newOrient = GetEdgeOrientationWRT(newEdge, contigID);

      // If the overlap is moving from the A End of the old contig to the B End of the
      // new contig, we need to adjust the orientation of the edge
      if(contigEnd != newContigEnd){
	switch( newOrient){
	case AB_BA:
	  assert(contigEnd == B_END);
	  newOrient = BA_BA;
	  break;
	case AB_AB:
	  assert(contigEnd == B_END);
	  newOrient = BA_AB;
	  break;
	case BA_AB:
	  assert(contigEnd == A_END);
	  newOrient = AB_AB;
	  break;
	case BA_BA:
	  assert(contigEnd == A_END);
	  newOrient = AB_BA;
	  break;
	default:
	  assert(0);
      }

      }

#if 0
      // If we flipped to make canonical, flip it back
      if(flip){
	newOrient = FlipEdgeOrient(newOrient);
      }
#endif
    }

    // Modify the ida and idb to reflect the new contig
    {
      CDS_CID_t newEdgeID = GetVAIndex_EdgeCGW_T(ScaffoldGraph->RezGraph->edges, newEdge);

      assert(rawEdge->idA == contigID || rawEdge->idB == contigID);
      if(otherContig < newContig->id){
	newEdge->orient = FlipEdgeOrient(newOrient);
	newEdge->idA = otherContig;
	newEdge->idB = newContig->id;
      }else{
	newEdge->orient = newOrient;
	newEdge->idB = otherContig;
	newEdge->idA = newContig->id;
      }

	{
	  ChunkOverlapCheckT olap;
	  if(LookupOverlap(ScaffoldGraph->RezGraph, newEdge->idA, newEdge->idB, newEdge->orient, &olap)){
	    fprintf(GlobalData->stderrc,"\t* Overlap already exists (" F_CID "," F_CID ",%c)...skip\n",
		    newEdge->idA, newEdge->idB, newEdge->orient);
	    FreeGraphEdge(ScaffoldGraph->RezGraph, newEdge);
	    continue;
	  }
	}

	{
	  CDS_COORD_t delta = min((3 * sqrt(newEdge->distance.variance)), 50);

	  ChunkOverlapCheckT olap =
	    OverlapChunks(ScaffoldGraph->RezGraph,
                          newEdge->idA, newEdge->idB,
                          newEdge->orient, // handles suspicious
			  -newEdge->distance.mean - delta, 
			  -newEdge->distance.mean + delta, 
			  CGW_DP_ERATE, FALSE);
          
	  if(olap.overlap == 0 || olap.suspicious){
	    if(olap.suspicious){
	      fprintf(GlobalData->stderrc,"* CIS_S2 SUSPICIOUS Overlap found! Looked for (" F_CID "," F_CID ",%c) found (" F_CID "," F_CID ",%c)\n",
		      newEdge->idA, newEdge->idB, newEdge->orient,
		      olap.spec.cidA, olap.spec.cidB, olap.spec.orientation);
	    }
	    if(verbose) {
	      PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,
                             "X No Overlap!  ", newEdge, newEdge->idA);
            }
	    olap = OverlapChunks(ScaffoldGraph->RezGraph,
                                 newEdge->idA, newEdge->idB,
                                 FlipEdgeOrient(newEdge->orient), // handles suspicious
                                 -newEdge->distance.mean - delta, 
                                 -newEdge->distance.mean + delta, 
                                 CGW_DP_ERATE, FALSE);
            
	    if(olap.overlap > 0 && !olap.suspicious)
	      fprintf(GlobalData->stderrc,"*** DUMMY, you got the orientation backwards!!!! ****\n");
	    FreeGraphEdge(ScaffoldGraph->RezGraph,newEdge);
	    continue;
	  }
	}

	
	if(verbose){
	  fprintf(GlobalData->stderrc,"*$#$ Inserting the following edge in the graph\n");
	  PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph, "  X ", newEdge, newEdge->idA);
       }
      // Link the edge back into the graph
      InsertGraphEdge(ScaffoldGraph->RezGraph, newEdgeID, FALSE);

      // Create an appropriate hash table entry
      CreateChunkOverlapFromEdge(ScaffoldGraph->RezGraph, newEdge, FALSE);

    }
  }
  if(verbose){
    DumpContig(GlobalData->stderrc,ScaffoldGraph, newContig,FALSE);
  }
  RecycleDeletedGraphElements(ScaffoldGraph->RezGraph);

}

void PropagateOverlapsToNewContig(ContigT *contig,
                                  VA_TYPE(IntElementPos) *ContigPositions,
                                  CDS_CID_t aEndID, int aEndEnd, 
                                  CDS_CID_t bEndID, int bEndEnd,
                                  CDS_CID_t scaffoldID, int32 verbose){
  ContigT *aContig, *bContig;
  CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, scaffoldID);
  CDS_COORD_t contigBase = CDS_COORD_MAX;
  int32 i;


  // Get a baseline for contig offsets
  for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions,i);
    CDS_COORD_t minPos = min(pos->position.bgn, pos->position.end);
    if(contigBase > minPos)
      contigBase = minPos;
  }

  AssertPtr(scaffold);
#ifdef DEBUG
  fprintf(GlobalData->stderrc,"* PropagateOverlaps...initial\n");
  DumpContig(GlobalData->stderrc,ScaffoldGraph, GetGraphNode(ScaffoldGraph->RezGraph, contig->id),FALSE);
  fprintf(GlobalData->stderrc,"* Calling internal *\n");
#endif

  // First, propagate the tandem overlap marks to the new contig
  //
  aContig = GetGraphNode(ScaffoldGraph->ContigGraph, aEndID);
  bContig = GetGraphNode(ScaffoldGraph->ContigGraph, bEndID);


  // First propagate tandem marks and branch opints to the new contig
  if(aEndEnd == A_END){
    if(aContig->flags.bits.tandemOverlaps & AEND_TANDEM_OVERLAP){
      contig->flags.bits.tandemOverlaps |= AEND_TANDEM_OVERLAP;
    }
    contig->info.Contig.branchPointA = aContig->info.Contig.branchPointA;
  }else if(aEndEnd == B_END){
    if(aContig->flags.bits.tandemOverlaps & BEND_TANDEM_OVERLAP){
      contig->flags.bits.tandemOverlaps |= AEND_TANDEM_OVERLAP;
    }
    contig->info.Contig.branchPointA = aContig->info.Contig.branchPointB;
  }else assert(0 /* Should be A_END or B_END */);


  if(bEndEnd == A_END){
    if(bContig->flags.bits.tandemOverlaps & AEND_TANDEM_OVERLAP){
      contig->flags.bits.tandemOverlaps |= BEND_TANDEM_OVERLAP;
    }
    contig->info.Contig.branchPointB = bContig->info.Contig.branchPointA;
  }else if(bEndEnd == B_END){
    if(bContig->flags.bits.tandemOverlaps & BEND_TANDEM_OVERLAP){
      contig->flags.bits.tandemOverlaps |= BEND_TANDEM_OVERLAP;
    }
    contig->info.Contig.branchPointB = bContig->info.Contig.branchPointB;
  }else assert(0 /* Should be A_END or B_END */);



#ifdef DEBUG
  DumpContig(GlobalData->stderrc,ScaffoldGraph, GetGraphNode(ScaffoldGraph->RezGraph, contig->id),FALSE);

  fprintf(GlobalData->stderrc,"* Calling extermal on " F_CID " *\n", aEndID);
#endif

// Propagate overlaps from the contig at the a end of the new contig
  PropagateExtremalOverlapsToNewContig(aEndID, aEndEnd, contig, A_END, contigBase, verbose);

#ifdef DEBUG
  fprintf(GlobalData->stderrc,"* Calling extermal on " F_CID " *\n", bEndID);
  DumpContig(GlobalData->stderrc,ScaffoldGraph, GetGraphNode(ScaffoldGraph->RezGraph, contig->id),FALSE);
#endif
  // Propagate overlaps from the contig at the b end of the enw contig
  PropagateExtremalOverlapsToNewContig(bEndID, bEndEnd, contig, B_END, contigBase, verbose);

  // Propagate Internal, non-containment overlaps
  PropagateInternalOverlapsToNewContig(contig, ContigPositions, scaffold->id, contigBase, verbose);

#ifdef DEBUG
  DumpContig(GlobalData->stderrc,ScaffoldGraph, GetGraphNode(ScaffoldGraph->RezGraph, contig->id),FALSE);
  fprintf(GlobalData->stderrc,"* Calling containment *\n");
#endif
  // Propagate Containment Overlaps
  PropagateContainmentOverlapsToNewContig(contig, ContigPositions, contigBase, verbose);

#ifdef DEBUG
  DumpContig(GlobalData->stderrc,ScaffoldGraph, GetGraphNode(ScaffoldGraph->RezGraph, contig->id),FALSE);
#endif
}





// Insert the new contig into the scaffold, in lieu of the old contigs
void  ReplaceContigsInScaffolds(CIScaffoldT *scaffold, ContigT *newContig, VA_TYPE(IntElementPos) *contigPositions, 
				LengthT offsetAEnd, LengthT offsetBEnd,   // The offsets of the two contig Ends
				LengthT deltaOffsetBEnd                   // The amount to adjust the rest of the scaffold
				){
  int32 i;
  int32 numElements = GetNumIntElementPoss(contigPositions);
  int orientAB = FALSE;
  LengthT *maxOffset, *minOffset;
  assert(numElements> 0);


  if(offsetAEnd.mean < offsetBEnd.mean){
    orientAB = TRUE;
    maxOffset = &offsetBEnd;
    minOffset = &offsetAEnd;
  }else{
    orientAB = FALSE;
    maxOffset = &offsetAEnd;
    minOffset = &offsetBEnd;
  }

#ifdef DEBUG_DETAILED
  fprintf(GlobalData->stderrc,"* Inserting new contig " F_CID " at (%g,%g) "
          "with delta (%g,%g) in place of:\n", 
	  newContig->id, offsetAEnd.mean, offsetBEnd.mean,
          deltaOffsetBEnd.mean, deltaOffsetBEnd.variance);
#endif  
  for(i = 0; i < numElements; i++){
    IntElementPos *pos = GetIntElementPos(contigPositions,i);
    ContigT *contig = GetGraphNode(ScaffoldGraph->ContigGraph,  pos->ident);
#ifdef DEBUG_DETAILED
    fprintf(GlobalData->stderrc,"* Removing contig " F_CID " "
            "from scaffold " F_CID " prev:" F_CID " next:" F_CID " (%g,%g) (%g,%g) ratio:%g\n",
	    contig->id, scaffold->id,  contig->AEndNext,contig->BEndNext,
	    contig->offsetAEnd.mean, contig->offsetAEnd.variance,
	    contig->offsetBEnd.mean, contig->offsetBEnd.variance,
	    abs(contig->offsetAEnd.variance - contig->offsetBEnd.variance)
            /contig->bpLength.variance);
#endif    
    // Remove all of the contigs that were previously scaffolded
    if(contig->scaffoldID == scaffold->id)
      RemoveCIFromScaffold(ScaffoldGraph, scaffold, contig, FALSE);
    contig->info.Contig.AEndCI = contig->info.Contig.BEndCI = NULLINDEX;
    DeleteGraphNode(ScaffoldGraph->ContigGraph, contig);
  }

  // Make sure scaffold length is OK.  This is probably unnecessary
  SetCIScaffoldTLength(ScaffoldGraph, scaffold, TRUE);


  InsertCIInScaffold(ScaffoldGraph, newContig->id, scaffold->id, offsetAEnd, offsetBEnd, TRUE, NO_CONTIGGING);


  if(newContig->AEndNext == NULLINDEX){ // Make sure we atart at the end
    LengthT offset;

    if(orientAB){
      offset.mean = -newContig->offsetAEnd.mean;
      offset.variance = -newContig->offsetAEnd.variance;
    }else{
      offset.mean = -newContig->offsetBEnd.mean;
      offset.variance = -newContig->offsetBEnd.variance;
    }

    fprintf(stderr,"* newContig " F_CID " at AEND offset [%g,%g] \n", newContig->id, newContig->offsetAEnd.mean, newContig->offsetBEnd.mean);
    AddDeltaToScaffoldOffsets(ScaffoldGraph,
			      scaffold->id,
			      newContig->id,
			      TRUE,
			      FALSE,
			      offset);
  }
  

  if( newContig->BEndNext != NULLINDEX ){ // contig is in the beginning or middle
    AddDeltaToScaffoldOffsets(ScaffoldGraph,
			      scaffold->id,
			      newContig->BEndNext,
			      TRUE,
			      FALSE,
			      deltaOffsetBEnd);
  }
  /* Shouldn't need this */
  if(scaffold->bpLength.mean < (maxOffset->mean - minOffset->mean))
    scaffold->bpLength.mean = maxOffset->mean - minOffset->mean;
  if(scaffold->bpLength.variance < maxOffset->variance - minOffset->variance)
    scaffold->bpLength.variance = maxOffset->variance - minOffset->variance;

}

void ReScaffoldPseudoDegenerates(void)
{
  GraphNodeIterator     nodes;
  ContigT		*ctg;
  int numMatches = 0;
  LengthT NullLength = {0.0, 0.0};
  
  InitGraphNodeIterator(&nodes, ScaffoldGraph->RezGraph, GRAPH_NODE_DEFAULT);
  while((ctg = NextGraphNodeIterator(&nodes)) != NULL)
  {

    if(ctg->flags.bits.isDead) continue;
    if(ctg->flags.bits.isChaff) continue;
    if(ctg->info.Contig.numCI > 1) continue;
    
    if(ctg->scaffoldID == NULLINDEX)
    {
      NodeCGW_T * utg = GetGraphNode(ScaffoldGraph->CIGraph,
                                     ctg->info.Contig.AEndCI);
      if(utg->type == UNRESOLVEDCHUNK_CGW &&
         utg->scaffoldID != NULLINDEX &&
         utg->info.CI.numInstances == 0)
      {
        // if here, we've got one
        CIScaffoldT CIScaffold;

        InitializeScaffold(&CIScaffold, REAL_SCAFFOLD);
        CIScaffold.info.Scaffold.AEndCI = NULLINDEX;
        CIScaffold.info.Scaffold.BEndCI = NULLINDEX;
        CIScaffold.info.Scaffold.numElements = 0;
        CIScaffold.edgeHead = NULLINDEX;
        CIScaffold.bpLength = NullLength;
        CIScaffold.id = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
        CIScaffold.flags.bits.isDead = FALSE;
        CIScaffold.aEndCoord = CIScaffold.bEndCoord = -1;
        CIScaffold.numEssentialA = CIScaffold.numEssentialB = 0;
        CIScaffold.essentialEdgeB = CIScaffold.essentialEdgeA = NULLINDEX;
        AppendGraphNode(ScaffoldGraph->ScaffoldGraph, &CIScaffold);
        InsertCIInScaffold(ScaffoldGraph, ctg->id, CIScaffold.id,
                           NullLength, ctg->bpLength, FALSE, FALSE);
        numMatches++;
      }
    }
  }
}
/****************************************************************************/
int CleanupScaffolds(ScaffoldGraphT *sgraph, int lookForSmallOverlaps,
                     int maxContigsInMerge,
                     int deleteUnMergedSurrogates){
  CIScaffoldT *scaffold;
  GraphNodeIterator scaffolds;
  int didSomething = FALSE;

    if(!GlobalData->performCleanupScaffolds){
      return 0;
    }

#ifdef DEBUG_DETAILED
    fprintf(GlobalData->stderrc,"* CleanupAScaffolds\n");
#endif
    
  InitGraphNodeIterator(&scaffolds, sgraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){

    if(scaffold->type == REAL_SCAFFOLD)
      didSomething |= CleanupAScaffold(sgraph,scaffold, lookForSmallOverlaps, maxContigsInMerge, deleteUnMergedSurrogates);
    
    if((scaffold->id % 10000) == 0){
      fprintf(GlobalData->stderrc,"* CleanupScaffolds through scaffold " F_CID "\n", scaffold->id);
      fflush(GlobalData->stderrc);
    }

#if 0
    if(GetSizeOfCurrentSequenceDB(ScaffoldGraph->sequenceDB) > GlobalData->maxSequencedbSize){
      fprintf(GlobalData->timefp,"\n\nCheckpoint %d written during CleanupScaffolds after scaffold " F_CID "\n",ScaffoldGraph->checkPointIteration, scaffold->id);
      CheckpointScaffoldGraph(ScaffoldGraph);
    }
#endif
    CheckScaffoldGraphCache(ScaffoldGraph);
  }

  RecycleDeletedGraphElements(sgraph->RezGraph);
  return didSomething;
}

/****************************************************************************/
int CleanupFailedMergesInScaffolds(ScaffoldGraphT *sgraph){
  CIScaffoldT *scaffold;
  GraphNodeIterator scaffolds;
  int madeChanges = FALSE;
  fprintf(GlobalData->stderrc,"* CleanupFailedMergesInScaffolds\n");

  /*
    yanked rocks/stones may not be in scaffolds
    without this step, some will inappropriately be degenerates
  */
  ReScaffoldPseudoDegenerates();
  
    if(!GlobalData->performCleanupScaffolds){
      return 0;
    }

  InitGraphNodeIterator(&scaffolds, sgraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    int didSomething = TRUE;
    int iteration = 0;
    while(didSomething > 0){
      if(iteration == 0){
	didSomething = CleanupAScaffold(sgraph,scaffold, FALSE, 16 , FALSE);
	if(didSomething == 0)
	  break;
	if(didSomething > 0)
	  madeChanges = TRUE;
	fprintf(GlobalData->stderrc,"* Merge failure(16) in scaffold " F_CID " didSomething = %d \n", scaffold->id, didSomething);
	didSomething = CleanupAScaffold(sgraph,scaffold, FALSE, 4 , FALSE);
	if(didSomething == 0)
	  break;
	if(didSomething > 0)
	  madeChanges = TRUE;
	fprintf(GlobalData->stderrc,"* Merge failure(4) in scaffold " F_CID " didSomething = %d \n", scaffold->id, didSomething);

      }
      didSomething = CleanupAScaffold(sgraph,scaffold, FALSE, 3 , FALSE);
      if(didSomething == 0)
	break;
      if(didSomething > 0)
	  madeChanges = TRUE;
      fprintf(GlobalData->stderrc,"* Merge failure(3) in scaffold " F_CID " didSomething = %d iteration %d\n", scaffold->id, didSomething, iteration);
      didSomething = CleanupAScaffold(sgraph,scaffold, FALSE, 2 , FALSE);
      if(didSomething > 0)
	madeChanges = TRUE;
      fprintf(GlobalData->stderrc,"* After iteration %d  didSomething = %d\n", iteration, didSomething);
      iteration++;
    }
  }

  RecycleDeletedGraphElements(sgraph->RezGraph);
  return madeChanges;
}

	  
/****************************************************************************/
int  DeleteAllSurrogateContigsFromFailedMerges(CIScaffoldT *scaffold,
                                               NodeCGW_T *contig){
  int didSomething = FALSE;
  int32 i;
  MultiAlignT *ma;    
  int numSurrogates;

  fprintf(GlobalData->stderrc,"* DeleteAllSurrogateContigsFromFailedMerges scaffold " F_CID " contig " F_CID "\n",
	  scaffold->id, contig->id);


  if(contig->bpLength.mean > 2000)
     return FALSE;

    numSurrogates = 0;
    ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);
    //    ma = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id);
    for(i = 0; i < GetNumIntUnitigPoss(ma->u_list); i++){
      IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);
      NodeCGW_T *node = GetGraphNode(ScaffoldGraph->CIGraph, pos->ident);
      if(node->flags.bits.isSurrogate)
	numSurrogates++;
    }

    fprintf(GlobalData->stderrc,"* numSurrogates:%d  numCI :%d numElementPos:%d\n",
	    numSurrogates, contig->info.Contig.numCI,
            (int) GetNumIntUnitigPoss(ma->u_list));

      if(numSurrogates == contig->info.Contig.numCI){
	NodeCGW_T *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, contig->scaffoldID);

	didSomething = TRUE;
	fprintf(GlobalData->stderrc,"*** Deleting contig " F_CID " from scaffold " F_CID " since it is a short, all surrogate contigging failue\n",
		contig->id, contig->scaffoldID);
	DumpContig(GlobalData->stderrc,ScaffoldGraph, contig, FALSE);
        
	/* Remove the Contig from the scaffold */
	RemoveCIFromScaffold(ScaffoldGraph, scaffold, contig, FALSE);

	/* Delete the contig */
	DeleteGraphNode(ScaffoldGraph->ContigGraph, contig);

	/* Delete all of the surrogate CIs */
	for(i = 0; i < GetNumIntUnitigPoss(ma->u_list); i++){
	  IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);
	  NodeCGW_T *node = GetGraphNode(ScaffoldGraph->CIGraph, pos->ident);
	  DeleteGraphNode(ScaffoldGraph->CIGraph, node);
	}
	
      }
    
      return didSomething;
}



static VA_TYPE(IntElementPos) *ContigPositions = NULL;

/****************************************************************************/
// Cleanup contigs for a single scaffold
//n
// If after resolving all repeats, we have placed Contigs that overlap by > threshhold,
// we merge them together into a single contig
//
int CleanupAScaffold(ScaffoldGraphT *graph, CIScaffoldT *scaffold,
                     int lookForSmallOverlaps, 
                     int maxContigsInMerge,
                     int deleteUnmergedSurrogates){

  CIScaffoldTIterator CIs;
  ChunkInstanceT *CI;
  ChunkInstanceT *currCI = NULL;
  ChunkInstanceT *prevCI = NULL;
  ContigEndsT contig;
  int32 i;
  int mergesAttempted = 0;
  int allMergesSucceeded = FALSE;
#if 0
  static time_t lastCkpTime = (time_t) -1;
#endif

  if(!GlobalData->performCleanupScaffolds){
    return 0;
  }
#ifdef DEBUG_DETAILED
  fprintf(GlobalData->stderrc,"* CleanupAScaffold " F_CID " (max:%d)\n", scaffold->id, maxContigsInMerge);
#endif
    
#ifdef DEBUG_CONNECTEDNESS
  // THE FOLLOWING IS DEBUG CODE
  // MAKE SURE WE DIDN'T DISCONNECT THE SCAFFOLD
  {
    float maxVariance = 1000000.0;// useGuides ? 100000000000.0 : 1000000.0;
    int numComponents;


    MarkInternalEdgeStatus(graph, scaffold, PAIRWISECHI2THRESHOLD_CGW,
			   maxVariance, TRUE, TRUE, 0);

    numComponents = IsScaffoldInternallyConnected(graph,scaffold);
    assert(numComponents == 1);
  }
#endif


  if(ContigEnds == NULL){
    ContigEnds = CreateVA_ContigEndsT(10);
  }else{
    ResetContigEndsT(ContigEnds);
  }

  contig.firstCID = contig.lastCID = NULLINDEX;
  contig.firstCI = contig.lastCI = NULL;

  /* First, figure out where the contigs are.  We can't contig them right away, since that
       would be mucking with the data structure we're iterating over */
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
  /* Initialize CI . . . ELA */
  contig.firstCI = currCI = CI = NextCIScaffoldTIterator(&CIs);
  contig.firstCID = contig.firstCI->id; //GetVAIndex_ChunkInstanceT(graph->ChunkInstances, contig.firstCI);
  AssertPtr(CI);		// There should be at least one CI . . ELA
  if(currCI->offsetAEnd.mean < currCI->offsetBEnd.mean){
    contig.minOffset = currCI->offsetAEnd;
    contig.maxOffset = currCI->offsetBEnd;
  }else{
    contig.maxOffset = currCI->offsetAEnd;
    contig.minOffset = currCI->offsetBEnd;
  }
  contig.count = 0;
  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
    prevCI = currCI;
    currCI = CI;

    contig.count++;
    contig.lastCI = prevCI;
    contig.lastCID = contig.lastCI->id; //GetVAIndex_ChunkInstanceT(graph->ChunkInstances, contig.lastCI);
    // The first CI in the contig must have the minimum offset, by definition
    // subsequent CIs are sorted by LEFTmost coordinate, so we need to track
    // the MAX span (rightmost coordinate)

    {
      CDS_COORD_t actual;
      // ChunkOrientationType pairwiseOrient;

#undef DEBUG_CONTIG

#ifdef DEBUG_CONTIG
      fprintf(GlobalData->stderrc,"* contig.min = %g contig.max = %g\n",
	      contig.minOffset.mean, contig.maxOffset.mean);
#endif
      actual = IntervalsOverlap(contig.minOffset.mean, contig.maxOffset.mean,
				currCI->offsetAEnd.mean, currCI->offsetBEnd.mean, -15000);

      if(((maxContigsInMerge == NULLINDEX) || (contig.count < maxContigsInMerge)) &&  // artificially break longs merges if != NULLINDEX
	 actual> CGW_DP_MINLEN){ /*  || (lookForSmallOverlaps && SmallOverlapExists(ScaffoldGraph->RezGraph, prevCI, currCI, pairwiseOrient))){ */
#ifdef DEBUG_CONTIG
	fprintf(GlobalData->stderrc,
		"* CI " F_CID " and " F_CID " mean positions overlap by " F_COORD " (%g,%g) (%g,%g)\n",
		prevCI->id,currCI->id, actual, contig.minOffset.mean, 
		contig.maxOffset.mean, currCI->offsetAEnd.mean, 
		currCI->offsetBEnd.mean);
#endif

	if(currCI->offsetAEnd.mean > contig.maxOffset.mean){
	  contig.maxOffset = currCI->offsetAEnd;
	}
	if(currCI->offsetBEnd.mean > contig.maxOffset.mean){
	  contig.maxOffset = currCI->offsetBEnd;
	}

      }else{
#ifdef DEBUG_CONTIG
	fprintf(GlobalData->stderrc,"* CI " F_CID " and " F_CID " mean positions GAP by " F_COORD " (%g,%g) (%g,%g) first " F_CID " last " F_CID " \n",
		prevCI->id,currCI->id, actual, 
		contig.minOffset.mean,  contig.maxOffset.mean, 
		currCI->offsetAEnd.mean, currCI->offsetBEnd.mean, 
		contig.firstCI->id, contig.lastCI->id);
#endif
	if((maxContigsInMerge != NULLINDEX) &&
	   (contig.count == maxContigsInMerge)){
#ifdef DEBUG_DETAILED
	  fprintf(GlobalData->stderrc,"*^^^^^ Merge terminated at %d elements\n", maxContigsInMerge);
#endif
	}
	AppendContigEndsT(ContigEnds, &contig);
	contig.count = 0;
	contig.firstCI = currCI;
	contig.firstCID = contig.firstCI->id; //GetVAIndex_ChunkInstanceT(graph->ChunkInstances, contig.firstCI);

	if(currCI->offsetAEnd.mean < currCI->offsetBEnd.mean){
	  contig.minOffset = currCI->offsetAEnd;
	  contig.maxOffset = currCI->offsetBEnd;
	}else{
	  contig.maxOffset = currCI->offsetAEnd;
	  contig.minOffset = currCI->offsetBEnd;
	}
#ifdef DEBUG_CONTIG
	fprintf(GlobalData->stderrc,"* Reseting contig  firstCI: " F_CID " min:%g max:%g\n",
		contig.firstCI->id, contig.minOffset.mean, contig.maxOffset.mean);
#endif
      }
    }
  }
  /* Add last contig to the list */
  ++contig.count;
  contig.lastCI = currCI;
  contig.lastCID = contig.lastCI->id; //GetVAIndex_ChunkInstanceT(graph->ChunkInstances, contig.lastCI);
  AppendContigEndsT(ContigEnds, &contig);


    // Now we have a collection of contigs to build
    // Work through them one at a time
  {
    ContigEndsT *ctg = GetContigEndsT(ContigEnds,0);
    int32 numContigs = GetNumContigEndsTs(ContigEnds);
#if 0
    /******************************DEAD CODE *********************************/
    ContigT newContig;
    CDS_CID_t firstContig;
    int32 numContigs;

    InitializeContig(&newContig, CONTIG_CGW);
    firstContig = GetNumGraphNodes(graph->ContigGraph);
    numContigs = GetNumContigEndsTs(ContigEnds);

    newContig.flags.bits.isUnique = (scaffold->type == REAL_SCAFFOLD?TRUE:FALSE);
    newContig.flags.bits.cgbType = UU_CGBTYPE;
    newContig.aEndCoord = newContig.bEndCoord = -1;
    newContig.scaffoldID = ctg->firstCI->scaffoldID;
    assert(newContig.scaffoldID != NULLINDEX);
    newContig.edgeHead = NULLINDEX;
    newContig.microhetScore = NULLINDEX;
    /******************************END DEAD CODE *********************************/
#endif
    if(ContigPositions == NULL){
      ContigPositions = CreateVA_IntElementPos(10);
    }

    /* Now the ContigEnds VA contains the starts/ends of all the contigs */
    for(i = 0; i < numContigs; i++){
      ContigT *contig;
      IntElementPos contigPos;
      ResetVA_IntElementPos(ContigPositions);

      ctg = GetContigEndsT(ContigEnds,i);

      for(contig = GetGraphNode(graph->RezGraph, ctg->firstCID);
	  contig != NULL;
	  contig = GetGraphNode(graph->RezGraph, contig->BEndNext))
	{
	  contigPos.ident = contig->id;
	  contigPos.type = AS_CONTIG;
	  contigPos.position.bgn = contig->offsetAEnd.mean;
	  contigPos.position.end = contig->offsetBEnd.mean;
	  AppendIntElementPos(ContigPositions, &contigPos);
	      
	  if(contig->id == ctg->lastCID)
	    break;
	}

      // If there is no merging to be done, continue
      if(GetNumIntElementPoss(ContigPositions) <= 1)
	continue;
      mergesAttempted++;

      if(deleteUnmergedSurrogates){  // Just find the all-surrogate, unmerged contigs, and nuke them
	int32 i;
	for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
	  IntElementPos *pos = GetIntElementPos(ContigPositions, i);
	  contig = GetGraphNode(ScaffoldGraph->ContigGraph, pos->ident);
	  allMergesSucceeded |= DeleteAllSurrogateContigsFromFailedMerges(scaffold, contig);
	}
      }else{
	allMergesSucceeded &= CreateAContigInScaffold(scaffold, ContigPositions, ctg->minOffset, ctg->maxOffset);
      }
    }
  }
#ifdef DEBUG_CONNECTEDNESS
  // THE FOLLOWING IS DEBUG CODE
  // MAKE SURE WE DIDN'T DISCONNECT THE SCAFFOLD
  {
    float maxVariance = 1000000.0;//useGuides ? 100000000000.0 : 1000000.0;
    int numComponents;
    MarkInternalEdgeStatus(graph, scaffold, PAIRWISECHI2THRESHOLD_CGW,
			   maxVariance, TRUE, TRUE, 0);

    numComponents = IsScaffoldInternallyConnected(graph,scaffold);
    assert(numComponents == 1);
  }
#endif

#if 0
  // initialize the timer the first time we cleanup a scaffold
  if (lastCkpTime == (time_t) -1)
  	lastCkpTime = time(NULL);
  
  if (GetSizeOfCurrentSequenceDB(ScaffoldGraph->sequenceDB) > GlobalData->maxSequencedbSize
	  || time(NULL) - lastCkpTime > 3600 * 4)  // if we've got enough to write or it's been 4 hours
  {
	if (time(NULL) - lastCkpTime > 3600 * 3)  // even with enough to write we must go three hours to checkpoint
	{
	  fprintf(GlobalData->timefp,"\n\nCheckpoint %d written during CleanupAScaffold after scaffold " F_CID "\n",
			  ScaffoldGraph->checkPointIteration, scaffold->id);
	  // fprintf( stderr, "would write a checkpoint at time: " F_TIME_T "\n", time(NULL));
	  CheckpointScaffoldGraph(ScaffoldGraph);
	  lastCkpTime = time(NULL);
	}
  }
#endif

  if (GetSizeOfCurrentSequenceDB(ScaffoldGraph->sequenceDB) > GlobalData->maxSequencedbSize)
  {
	  fprintf(GlobalData->timefp,"\n\nCheckpoint %d written during CleanupAScaffold after scaffold " F_CID "\n",
			  ScaffoldGraph->checkPointIteration, scaffold->id);
	  CheckpointScaffoldGraph(ScaffoldGraph);
  }
  
  
  if(mergesAttempted == 0)
    return 0;
  if(allMergesSucceeded)
    return TRUE;
  //else
  return NULLINDEX;
}


#if 1
/***************************************************************************/
// CheckForContigs
// Insert chunk instance ci int scaffold sid at offset with orientation orient.
// offsetFromAEnd = offset of the end of the CI that is closest to the A end of the scaffold
// orient
int CheckForContigs(ScaffoldGraphT *sgraph,
                      CDS_CID_t cid,
                      CDS_CID_t sid,
                      LengthT offsetAEnd,
                      LengthT offsetBEnd)
{
  /*
    iterate over contigs in scaffold sid & register those that intersect
    offsetAEnd, offsetBEnd
    Identify the left-most
    Call CreateAContigInScaffold
   */
  CIScaffoldT * scaffold = GetGraphNode(sgraph->ScaffoldGraph, sid);
  CIScaffoldTIterator contigIterator;
  ChunkInstanceT * contig;
  IntElementPos pos;
  CDS_COORD_t minPos;
  CDS_COORD_t maxPos;
  LengthT myOffsetAEnd;
  LengthT myOffsetBEnd;
#if DEBUG > 0
  fprintf(GlobalData->stderrc,"* CheckForContigs scaffold " F_CID "\n", scaffold->id);
#endif
  if(ContigPositions == NULL){
    ContigPositions = CreateVA_IntElementPos(10);
  }
  ResetVA_IntElementPos(ContigPositions);

  // append the interval of this cid
  pos.type = AS_CONTIG;
  pos.ident = cid;
  pos.position.bgn = offsetAEnd.mean;
  pos.position.end = offsetBEnd.mean;
  AppendVA_IntElementPos(ContigPositions, &pos);
  if(offsetAEnd.mean < offsetBEnd.mean)
  {
    minPos = offsetAEnd.mean;
    myOffsetAEnd = offsetAEnd;
    maxPos = offsetBEnd.mean;
    myOffsetBEnd = offsetBEnd;
  }
  else
  {
    minPos = offsetBEnd.mean;
    myOffsetAEnd = offsetBEnd;
    maxPos = offsetAEnd.mean;
    myOffsetBEnd = offsetAEnd;
  }

#undef DEBUG_CHECKFORCTGS

  // iterate over all contigs in this scaffold to append those that overlap
  // with the offsetAEnd-offsetBEnd interval
  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE,
                          FALSE, &contigIterator);
  while((contig = NextCIScaffoldTIterator(&contigIterator)) != NULL)
  {
#ifdef DEBUG_CHECKFORCTGS
    fprintf(stderr,"  Testing contig %d ... ",contig->id);
#endif
    if(IntervalsOverlap(contig->offsetAEnd.mean, contig->offsetBEnd.mean,
                        minPos, maxPos, CGW_DP_MINLEN))
    {
#ifdef DEBUG_CHECKFORCTGS
      fprintf(stderr," overlap!\n");
#endif
      pos.ident = contig->id;
      pos.position.bgn = contig->offsetAEnd.mean;
      pos.position.end = contig->offsetBEnd.mean;
      if(contig->offsetAEnd.mean < contig->offsetBEnd.mean)
      {
        myOffsetAEnd = (myOffsetAEnd.mean > contig->offsetAEnd.mean) ?
          contig->offsetAEnd : myOffsetAEnd;
        myOffsetBEnd = (myOffsetBEnd.mean < contig->offsetBEnd.mean) ?
          contig->offsetBEnd : myOffsetBEnd;
      }
      else
      {
        myOffsetAEnd = (myOffsetAEnd.mean > contig->offsetBEnd.mean) ?
          contig->offsetBEnd : myOffsetAEnd;
        myOffsetBEnd = (myOffsetBEnd.mean < contig->offsetAEnd.mean) ?
          contig->offsetAEnd : myOffsetBEnd;
      }
      AppendVA_IntElementPos(ContigPositions, &pos);
#ifdef DEBUG_CHECKFORCTGS
      fprintf(stderr," adding overlap with %d [%g,%g]\n",
	      pos.ident,pos.position.bgn,pos.position.end);
#endif
    }
    else {
#ifdef DEBUG_CHECKFORCTGS
      fprintf(stderr," no overlap!\n");
#endif
      if(maxPos < min(contig->offsetAEnd.mean, contig->offsetBEnd.mean)){
#ifdef DEBUG_CHECKFORCTGS
	fprintf(stderr," done checking for contigs that overlap!\n");
#endif
	break;
      }
    }
  }

#ifdef DEBUG_CHECKFORCTGS
  fprintf(stderr,"  should insert a new contig based on %d old ones\n",
	  GetNumIntElementPoss(ContigPositions));
#endif

  if(GetNumIntElementPoss(ContigPositions) > 1)
  {
    qsort(GetIntElementPos(ContigPositions,0),
          GetNumIntElementPoss(ContigPositions),
          sizeof(IntElementPos),
          CompareIntElementPosByBgnPos);
    // set will be shifted to 0 offset in this function
    return CreateAContigInScaffold(scaffold, ContigPositions,
                                   myOffsetAEnd, myOffsetBEnd);
  }

  // if here, no overlapping contigs
  return FALSE;
}

#else


int CheckForContigs(ScaffoldGraphT *sgraph,
                    CDS_CID_t cid, CDS_CID_t sid,
                    LengthT offsetAEnd, LengthT offsetBEnd){
  ChunkInstanceT *CI = GetGraphNode(sgraph->RezGraph, cid);
  GraphEdgeIterator   edges;
  CIScaffoldT *scaffold = GetGraphNode(sgraph->ScaffoldGraph, sid);
  CIEdgeT *edge;
  IntElementPos pos, *basePos;
  CDS_CID_t cidOther;
  int numEdges = 0;
  int numContainmentEdges = 0;
  int numAEdges = 0;
  int numBEdges = 0;

  fprintf(GlobalData->stderrc,"* CheckForContigs scaffold " F_CID "\n", scaffold->id);
  CI->scaffoldID = NULLINDEX; // test!!!

  if(ContigPositions == NULL){
    ContigPositions = CreateVA_IntElementPos(10);
  }
  ResetVA_IntElementPos(ContigPositions);

  pos.type = AS_CONTIG;
  pos.ident = cid;
  // Orientation of inserted chunk is always A_B
  pos.position.bgn = 0;
  pos.position.end = CI->bpLength.mean;
  AppendIntElementPos(ContigPositions, &pos);
  basePos = GetIntElementPos(ContigPositions,0);


  /* Look for overlap edges to Contigs already in this scaffold.
     The edges must be 'mustOverlap' AND they must overlap by mean
     contig position...otherwise, don't contig them.
  */
  InitGraphEdgeIterator(sgraph->RezGraph, cid, 
			ALL_END, ALL_EDGES, 
			GRAPH_EDGE_DEFAULT,
			&edges); 
  while (edge = NextGraphEdgeIterator(&edges)) {
    ChunkInstanceT *CIOther;
    int end;
    CDS_COORD_t overlap;
    ChunkOrientationType orient = GetEdgeOrientationWRT(edge, cid);

    assert(edge != NULL);


    if(cid == edge->idA)
      cidOther = edge->idB;
    else
      cidOther = edge->idA;

    CIOther = GetGraphNode(sgraph->RezGraph, cidOther);
    
    if(CIOther->scaffoldID != sid || !isOverlapEdge(edge) ||
       isSingletonOverlapEdge(edge))
      continue;
    
      overlap = IntervalsOverlap(CIOther->offsetAEnd.mean, CIOther->offsetBEnd.mean,
					 CI->offsetAEnd.mean, CI->offsetBEnd.mean,
					 -15000);

      //      assert(0/*Why Does overlap > 0 screw this up???? */);

    fprintf(GlobalData->stderrc,"* other scaffold id = " F_CID " (" F_CID ") mustOverlap:%d mean:%g std:%g overlap:" F_COORD " [%g,%g] [%g,%g]\n",
	    CIOther->scaffoldID, sid,isMustOverlapEdge(edge),   edge->distance.mean, sqrt(edge->distance.variance), overlap,
	    CIOther->offsetAEnd.mean, CIOther->offsetBEnd.mean,
	    CI->offsetAEnd.mean, CI->offsetBEnd.mean);
    if(/*(overlap > 0) && */ 
       (isMustOverlapEdge(edge) ||
       (3 * sqrt(edge->distance.variance) <  -edge->distance.mean))){

      numEdges++;
      switch(orient){
      case AB_BA:
	end = B_END;
	pos.ident = cidOther;
	pos.position.end = basePos->position.end + edge->distance.mean;
	pos.position.bgn = pos.position.end + CIOther->bpLength.mean;
	AppendIntElementPos(ContigPositions, &pos);
	numBEdges++;
	break;

      case AB_AB:
	end = B_END;
	pos.ident = cidOther;
	pos.position.bgn = basePos->position.end + edge->distance.mean;
	pos.position.end = pos.position.bgn + CIOther->bpLength.mean;
	AppendIntElementPos(ContigPositions, &pos);
	numBEdges++;
	break;
      case BA_BA:
	end = A_END;
	pos.ident = cidOther;
	pos.position.end = basePos->position.bgn - edge->distance.mean;
	pos.position.bgn = pos.position.end - CIOther->bpLength.mean;
	AppendIntElementPos(ContigPositions, &pos);
	numAEdges++;
	break;

      case BA_AB:
	end = A_END;
	pos.ident = cidOther;
	pos.position.bgn = basePos->position.bgn - edge->distance.mean;
	pos.position.end = pos.position.bgn - CIOther->bpLength.mean;
	AppendIntElementPos(ContigPositions, &pos);
	numAEdges++;
	break;
      default:
	assert(0);
      }
      fprintf(GlobalData->stderrc,"* (" F_CID "," F_CID ",%c) [" F_COORD "," F_COORD "]\n",
	      cid, cidOther, orient, pos.position.bgn, pos.position.end);

      if(edge->flags.bits.hasContainmentOverlap)
	numContainmentEdges++;

      fprintf(GlobalData->stderrc,"* Inserting CI " F_CID " in scaffold " F_CID " implies contigging with CI " F_CID " by %s on %c end\n",
	      cid, sid, cidOther, 
	      (edge->flags.bits.hasContainmentOverlap?"Containment":"Overlap"),
	      (end == A_END? 'A':'B'));
      PrintGraphEdge(GlobalData->stderrc,sgraph->RezGraph, " ", edge, cid);
    }
  }

  if(numEdges >= 0){
    fprintf(GlobalData->stderrc,"* Inserting cid " F_CID " into scaffold " F_CID " implied %d contigging operations (%d contains)(%d A, %d B)\n",
	    cid, sid, numEdges, numContainmentEdges, numAEdges, numBEdges);
  }


  if(GetNumIntElementPoss(ContigPositions) > 1)
    { 

      return CreateAContigInScaffold(scaffold, ContigPositions, offsetAEnd, offsetBEnd);

    }
/* This path takes care of the case where there is no contigging
   implied by this insertion */

  return FALSE;

}
#endif

/***************************************************************************/
// CheckForContainmentContigs
// Insert chunk instance ci int scaffold sid at offset with orientation orient.
// offsetFromAEnd = offset of the end of the CI that is closest to the A end of the scaffold
// orient

int CheckForContainmentContigs(ScaffoldGraphT *sgraph, CDS_CID_t cid, CDS_CID_t sid, LengthT offsetAEnd, LengthT offsetBEnd){
  ChunkInstanceT *CI = GetGraphNode(sgraph->RezGraph, cid);
  ChunkInstanceT *otherCI;
  CIScaffoldT *scaffold = GetGraphNode(sgraph->ScaffoldGraph, sid);
  IntElementPos pos, *basePos;
  CIScaffoldTIterator CIs;
  UnitigOverlapType overlapType;
  int foundContainment;
  LengthT aEnd, bEnd;


  CI->offsetAEnd = offsetAEnd;
  CI->offsetBEnd = offsetBEnd;
  fprintf(GlobalData->stderrc,"* CheckForContainmentContigs scaffold " F_CID "\n", scaffold->id);
  CI->scaffoldID = NULLINDEX; // test!!!

  if(ContigPositions == NULL){
    ContigPositions = CreateVA_IntElementPos(10);
  }
  ResetVA_IntElementPos(ContigPositions);

  pos.type = AS_CONTIG;
  pos.ident = cid;
  // Orientation of inserted chunk is always A_B
  pos.position.bgn = offsetAEnd.mean;
  pos.position.end = offsetBEnd.mean;
  AppendIntElementPos(ContigPositions, &pos);
  basePos = GetIntElementPos(ContigPositions,0);
  aEnd = offsetAEnd;
  bEnd = offsetBEnd;

  foundContainment = FALSE;

  /* See if the inserted element is contained, or contains an existing
     scaffold element(s).  If so, contig them. */

  InitCIScaffoldTIterator(sgraph, scaffold, TRUE, FALSE, &CIs);
  while((otherCI = NextCIScaffoldTIterator(&CIs)) != NULL){
    overlapType = existsContainmentRelationship(otherCI,CI);
    if(overlapType == AS_1_CONTAINS_2_OVERLAP ||
       overlapType == AS_2_CONTAINS_1_OVERLAP){
      foundContainment = TRUE;
      pos.ident = otherCI->id;
      pos.position.bgn = otherCI->offsetAEnd.mean;
      pos.position.end = otherCI->offsetBEnd.mean;
      if(pos.position.bgn < pos.position.end){
	if(pos.position.bgn < aEnd.mean){
	  aEnd = otherCI->offsetAEnd;
	}
	if(pos.position.end > bEnd.mean){
	  bEnd = otherCI->offsetBEnd;
	}
      }else{
	if(pos.position.end < aEnd.mean){
	  aEnd = otherCI->offsetBEnd;
	}
	if(pos.position.bgn > bEnd.mean){
	  bEnd = otherCI->offsetAEnd;
	}
      }
      AppendIntElementPos(ContigPositions, &pos);

    }
  }

  if(foundContainment == FALSE)
    return FALSE;


  if(GetNumIntElementPoss(ContigPositions) > 1)  {
	GraphEdgeIterator edges;
	EdgeCGW_T *edge;
	int32 i;

	for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
	  IntElementPos *ctg = GetIntElementPos(ContigPositions,i);
	  fprintf(GlobalData->stderrc,"* Inserting cid:" F_CID " into scaffold" F_CID " implied contigging with " F_CID "\n",
		  cid,sid, ctg->ident);

	  // Save the edge status in the frags so we can reconstitute it later
	  // in BuildGraphEdgesFromMultiAlign

	  InitGraphEdgeIterator(sgraph->RezGraph, ctg->ident, 
				ALL_END, ALL_EDGES, 
				GRAPH_EDGE_DEFAULT,
				&edges); 
          while((edge = NextGraphEdgeIterator(&edges)) != NULL)
            PropagateEdgeStatusToFrag(sgraph->RezGraph, edge);

	}
      return CreateAContigInScaffold(scaffold, ContigPositions, 
				     aEnd, bEnd);

      }
/* This path takes care of the case where there is no contigging
   implied by this insertion */

  return FALSE;

}


void DumpMultiAlignT(FILE * fp, ScaffoldGraphT * graph,
                     ContigT * ctg, MultiAlignT * ma)
{
  int32 i;
  CDS_COORD_t minPos;
  CDS_COORD_t maxPos;
  int32 numFrag = GetNumIntMultiPoss(ma->f_list);
  int32 numUnitig = GetNumIntUnitigPoss(ma->u_list);
  IntUnitigPos *up = GetIntUnitigPos(ma->u_list,0);

  CIScaffoldT *scaffold = GetGraphNode(graph->ScaffoldGraph,
                                       ctg->scaffoldID);

  fprintf(fp,"* Contig(%s) " F_CID " with %d unitigs and %d frags\n",
          scaffold && (scaffold->type == REAL_SCAFFOLD) ?
           "AS_PLACED": "AS_UNPLACED",
          ctg->id, numUnitig, numFrag);
  fprintf(fp, " %dbp gapped, %dbp ungapped\n",
          (int) GetMultiAlignLength(ma),
          (int) GetMultiAlignUngappedLength(ma));

  minPos = CDS_COORD_MAX;
  maxPos = CDS_COORD_MIN;
  for(i = 0; i < numUnitig; i++){
    IntUnitigPos *iup = up + i;
    NodeCGW_T *unitig = GetGraphNode(graph->CIGraph, iup->ident);
    fprintf(fp, "unitig " F_CID ": ", iup->ident);
    if(unitig->type == DISCRIMINATORUNIQUECHUNK_CGW){
      fprintf(fp, "AS_UNIQUE_UNITIG ");
    }else{
      if(unitig->scaffoldID != NULLINDEX){
        if(!unitig->flags.bits.isSurrogate){
          fprintf(fp, "AS_ROCK_UNITIG ");
        }else  if(unitig->flags.bits.isStoneSurrogate){
          fprintf(fp, "AS_STONE_UNITIG ");
        }else{
          fprintf(fp, "AS_PEBBLE_UNITIG ");
        }
      }else{
        fprintf(fp, "AS_SINGLE_UNITIG ");
      }
    }
    fprintf(fp, "pos:(" F_COORD "," F_COORD ") delta:%d\n",
            iup->position.bgn, iup->position.end,
            iup->delta_length);
    minPos = min(minPos, min(iup->position.bgn, iup->position.end));
    maxPos = max(maxPos, max(iup->position.bgn, iup->position.end));
  }
  fprintf(fp, "minPos = " F_COORD ", maxPos = " F_COORD "\n", minPos, maxPos);
  
  // Null out the source field
  for(i = 0; i < numFrag; i++){
    IntMultiPos *mp_i = GetIntMultiPos(ma->f_list,i);
    CIFragT *frag = GetCIFragT(graph->CIFrags,(CDS_CID_t)mp_i->source);
    fprintf(fp, "frag " F_CID ": ", frag->iid);
    switch(frag->type)
    {
      case AS_READ:
        fprintf(fp, "AS_READ ");
        break;
      case AS_EXTR:
        fprintf(fp, "AS_EXTR ");
        break;
      case AS_TRNR:
        fprintf(fp, "AS_TRNR ");
        break;
      case AS_EBAC:
        fprintf(fp, "AS_EBAC ");
        break;
      case AS_LBAC:
        fprintf(fp, "AS_LBAC ");
        break;
      case AS_UBAC:
        fprintf(fp, "AS_UBAC ");
        break;
      case AS_FBAC:
        fprintf(fp, "AS_FBAC ");
        break;
      case AS_STS:
        fprintf(fp, "AS_STS ");
        break;
      case AS_UNITIG:
        fprintf(fp, "AS_UNITIG ");
        break;
      case AS_CONTIG:
        fprintf(fp, "AS_CONTIG ");
        break;
      case AS_BACTIG:
        fprintf(fp, "AS_BACTIG ");
        break;
      case AS_FULLBAC:
        fprintf(fp, "AS_FULLBAC ");
        break;
      case AS_B_READ:
        fprintf(fp, "AS_B_READ ");
        break;
    }
    fprintf(fp, " utg:" F_CID ", pos:(" F_COORD "," F_COORD ")\n", frag->cid,
            mp_i->position.bgn, mp_i->position.end);
    minPos = min(minPos, min(mp_i->position.bgn, mp_i->position.end));
    maxPos = max(maxPos, max(mp_i->position.bgn, mp_i->position.end));
  }
  fprintf(fp, "minPos = " F_COORD ", maxPos = " F_COORD "\n", minPos, maxPos);
}

/***************************************************************************/
int  CreateAContigInScaffold(CIScaffoldT *scaffold,
                             VA_TYPE(IntElementPos) *ContigPositions,
                             LengthT offsetAEnd, LengthT offsetBEnd){
      /* This path takes care of the on-the-fly contigging of the inserted contig
	 with the existing contigs in the scaffld */
      ContigT *contig = CreateNewGraphNode(ScaffoldGraph->ContigGraph);
      LengthT newOffsetAEnd, newOffsetBEnd, oldOffsetBEnd, deltaOffsetBEnd;
      MultiAlignT  *newMultiAlign;
      int32 numElements = GetNumIntElementPoss(ContigPositions);
      CDS_CID_t i;
      CDS_CID_t aEndID = NULLINDEX, bEndID = NULLINDEX;
      CDS_CID_t first,last;
      int32 aEndEnd = NO_END, bEndEnd = NO_END;
      CDS_COORD_t minPos = CDS_COORD_MAX;
      CDS_COORD_t maxPos = CDS_COORD_MIN;
      int32 includesFinishedBacFragments = FALSE;
      contig->offsetAEnd = offsetAEnd;
      contig->offsetBEnd = offsetBEnd;

#ifdef DEBUG_DETAILED
      fprintf(GlobalData->stderrc,"* Create a contig in scaffold " F_CID "\n", scaffold->id);
#endif
      
      // Normalize the element pos values so their minimum is 0
      // Also, determine which ends are extremal, so we can propagate
      // their overlaps to the new contig
      first = GetIntElementPos(ContigPositions,0)->ident;
      last = GetIntElementPos(ContigPositions,numElements -1)->ident;

      for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
	IntElementPos *pos = GetIntElementPos(ContigPositions,i);
	ContigT *lcontig = GetGraphNode(ScaffoldGraph->ContigGraph,pos->ident);

	lcontig->flags.bits.failedToContig = FALSE;
	lcontig->flags.bits.beingContigged = TRUE; // flag the contig as being involved in a contigging operation

	includesFinishedBacFragments |= lcontig->flags.bits.includesFinishedBacFragments;

	if(pos->position.bgn < pos->position.end){
	  if(pos->position.bgn < minPos){
	    minPos = pos->position.bgn;
	    aEndID = pos->ident;
	    aEndEnd = A_END;
	  }
	  if(pos->position.end > maxPos){
	    maxPos = pos->position.end;
	    bEndID = pos->ident;
	    bEndEnd = B_END;
	  }
	}else{
	  if(pos->position.end < minPos){
	    minPos = pos->position.end;
	    aEndID = pos->ident;
	    aEndEnd = B_END;
	  }
	  if(pos->position.bgn > maxPos){
	    maxPos = pos->position.bgn;
	    bEndID = pos->ident;
	    bEndEnd = A_END;
	  }
	}
      }


#if 0
      fprintf(GlobalData->stderrc,"$$$$ * minPos = " F_COORD " maxPos = " F_COORD " extremes are (" F_CID ",%c) and (" F_CID ",%c)\n", 
	      minPos,maxPos,
	      aEndID,aEndEnd,
	      bEndID,bEndEnd );
#endif

      /*
      for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
        IntElementPos *pos = GetIntElementPos(ContigPositions,i);
        
        pos->position.bgn -= minPos;
        pos->position.end -= minPos;
      }
      */
      
	StartTimerT(&GlobalData->ConsensusTimer);

#ifdef DEBUG_DETAILED
        fprintf(GlobalData->stderrc,"* Calling MergeMultiAligns with %d elements*\n",
		(int) GetNumIntElementPoss(ContigPositions));
#endif        
	  {
#ifdef DEBUG_DETAILED
	    fprintf(GlobalData->stderrc,"*** Contigging ContigPositions\n");
#endif            
#ifdef DEBUG_CONTIG
	    for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
	      IntElementPos *pos = GetIntElementPos(ContigPositions,i);
	      fprintf(GlobalData->stderrc,"* Contig " F_CID "  bgn:" F_COORD " end:" F_COORD "\n",
		      pos->ident, pos->position.bgn, pos->position.end);
	    }
#endif
	  }
	

#if 1
	  //	newMultiAlign = MergeMultiAlignsFast(ScaffoldGraph->CIGraph->maStore, 
	  //					 ScaffoldGraph->ContigGraph->maStore,

        
	        newMultiAlign = MergeMultiAlignsFast_new(ScaffoldGraph->sequenceDB,
					 ScaffoldGraph->fragStore,
					 ContigPositions, FALSE, TRUE,GlobalData->aligner);
#else
	newMultiAlign = NULL;
#endif
	StopTimerT(&GlobalData->ConsensusTimer);

	if(!newMultiAlign){
#ifdef DEBUG_DETAILED
	  fprintf(GlobalData->stderrc,"* Calling MergeMultiAligns with %d elements*\n",
                  (int) GetNumIntElementPoss(ContigPositions));
#endif
	  {
	    int32 i;
      
#ifdef DEBUG_DETAILED
	    fprintf(GlobalData->stderrc,"*** Contigging ContigPositions\n");
#endif
	    for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
	      IntElementPos *pos = GetIntElementPos(ContigPositions,i);
	      ContigT *contig = GetGraphNode(ScaffoldGraph->ContigGraph, pos->ident);
              
	      contig->flags.bits.failedToContig = TRUE;
#ifdef DEBUG_CONTIG
              {
	        CDS_CID_t originalID = GetOriginalContigID(contig->id);
                fprintf(GlobalData->stderrc,"* Contig " F_CID " (" F_CID ")   bgn:" F_COORD " end:" F_COORD " [" F_COORD "," F_COORD "]\n",
                        pos->ident, originalID,
                        pos->position.bgn, pos->position.end,
                        contig->aEndCoord, contig->bEndCoord);
              }
#endif
	    }
	  }
	  DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffold, FALSE);
	  fprintf(GlobalData->stderrc,"* MergeMultiAligns failed....bye\n");
	  fflush(NULL);
	  if(GlobalData->failOn_NoOverlapFound)
	    assert(0 /* No overlap found error in merge multialignments */);
	  //	  SetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id, newMultiAlign); // DeleteGraphNode expects this
	  // NOTE: we do not keep the new multi-align in the cache, newMultiAlign is NULL anyways
	  InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, contig->id,FALSE, newMultiAlign, FALSE);
	  DeleteGraphNode(ScaffoldGraph->ContigGraph, contig);
	  return FALSE;
	}

	newMultiAlign->refCnt = 0;
	newMultiAlign->id = contig->id;
	contig->bpLength.mean = GetMultiAlignUngappedLength(newMultiAlign);
	contig->bpLength.variance = ComputeFudgeVariance(contig->bpLength.mean);

	contig->info.Contig.numCI =  GetNumIntUnitigPoss(newMultiAlign->u_list);
	contig->flags.bits.includesFinishedBacFragments = includesFinishedBacFragments;

      {
	NodeCGW_T *extremeA = GetGraphNode(ScaffoldGraph->ContigGraph, aEndID);
	NodeCGW_T *extremeB = GetGraphNode(ScaffoldGraph->ContigGraph, bEndID);

	if(aEndEnd == A_END)
	  newOffsetAEnd = extremeA->offsetAEnd;
	else
	  newOffsetAEnd = extremeA->offsetBEnd;

	// Determine the old mean/variance of the end of the contig
	if(bEndEnd == A_END)
	  oldOffsetBEnd = extremeB->offsetAEnd;
	else
	  oldOffsetBEnd = extremeB->offsetBEnd;

	// Determine the NEW mean/variance of the B end of the contig
	newOffsetBEnd.mean = newOffsetAEnd.mean + contig->bpLength.mean;
	newOffsetBEnd.variance = newOffsetAEnd.variance + contig->bpLength.variance;

	deltaOffsetBEnd.mean = newOffsetBEnd.mean - oldOffsetBEnd.mean;
	deltaOffsetBEnd.variance = newOffsetBEnd.variance - oldOffsetBEnd.variance;

      }


	if(GlobalData->debugLevel > 0)
	  fprintf(GlobalData->stderrc,"* Merged Contig " F_CID " has ungapped length %d\n",
		  contig->id, (int)contig->bpLength.mean);


#ifdef DEBUG_MERGE
	{
	  int32 i,j;
	  int32 numElements = GetNumIntUnitigPoss(newMultiAlign->u_list);
	  int32 numFrags = GetNumIntMultiPoss(newMultiAlign->f_list);
	  for(j = 0; j <numFrags ; j++){
	    IntMultiPos *fpos = GetIntMultiPos(newMultiAlign->f_list, j);
	    fprintf(GlobalData->stderrc,"* Fragment " F_CID " [" F_COORD "," F_COORD "]\n",(CDS_CID_t)fpos->source, fpos->position.bgn, fpos->position.end);
	  }

	  for(j = 0; j < numElements ; j++){
	    IntUnitigPos *upos = GetIntUnitigPos(newMultiAlign->u_list,j);
	    MultiAlignT *uma = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, upos->ident);
	    fprintf(GlobalData->stderrc,"\t* u" F_CID " [" F_COORD "," F_COORD "] ungapped length:%d\n",
		    upos->ident, upos->position.bgn, upos->position.end, (int) GetMultiAlignUngappedLength(uma));
	  }
	}
#endif


	assert(newMultiAlign->refCnt == 0);

	// Sort the Unitigs from left to right
	MakeCanonicalMultiAlignT(newMultiAlign);

	assert(newMultiAlign->refCnt == 0);

	//	SetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id, newMultiAlign);
	  // NOTE: we keep the new multi-align in the cache for a bit, but free it at the end of this routine
	InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, contig->id,FALSE, newMultiAlign, TRUE);

	assert(newMultiAlign->refCnt == 1);
	AddReferenceMultiAlignT(newMultiAlign); // this makes sure we survive any cache flushes

	// Propagate Overlaps, Tandem Marks and Branch Points to the new Contig
	PropagateOverlapsToNewContig(contig, ContigPositions, aEndID, aEndEnd, bEndID, bEndEnd, scaffold->id, FALSE);

	// Insert the new contig into the scaffold, in lieu of the old contigs
	ReplaceContigsInScaffolds(scaffold, contig,  ContigPositions, newOffsetAEnd, newOffsetBEnd, deltaOffsetBEnd);

	// Mark all frags as being members of this Contig, and set their offsets
	UpdateNodeFragments(ScaffoldGraph->RezGraph,contig->id, TRUE, FALSE);
	// Mark all of the Unitigs of this CI and set their offsets
	UpdateNodeUnitigs(newMultiAlign,contig);

	// Update simulator coordinates
	UpdateContigSimCoordinates(contig);
	UpdateScaffoldSimCoordinates(scaffold);

	// Create the raw link-based edges
	{ 
	  GraphEdgeStatT stats;
	  // This will get reclaimed when we flush the cache
	  //
	  //	  MultiAlignT *ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);	  
	  BuildGraphEdgesFromMultiAlign(ScaffoldGraph->ContigGraph, contig, newMultiAlign, &stats, TRUE);
	}

	if(GlobalData->debugLevel > 0){
	  fprintf(GlobalData->stderrc,"* Here is the new contig before merging: [%g,%g]\n",
		  newOffsetAEnd.mean, newOffsetBEnd.mean);
	  DumpContig(GlobalData->stderrc,ScaffoldGraph, contig,FALSE);
	}
	// Merge the edges incident on this contig
	MergeNodeGraphEdges(ScaffoldGraph->ContigGraph, contig, FALSE, TRUE, FALSE);

#ifdef DEBUG_DETAILED
	//	if(GlobalData->debugLevel > 0){
	  fprintf(GlobalData->stderrc,"*  >>> NEW CONTIG " F_CID " (multiAlign size = " F_SIZE_T ") <<< *\n", contig->id,
		  GetMemorySize(newMultiAlign));
	  DumpContig(GlobalData->stderrc,ScaffoldGraph, contig,FALSE);
	  //	}
#endif
	  RemoveReferenceMultiAlignT(newMultiAlign); // we added a reference previously
	  if(newMultiAlign->refCnt == 0){// we are the owner
	    DeleteMultiAlignT(newMultiAlign);
	  }else{ // the cache owns the memory
	    // Free up the cache space from the new multiAlignT
	    UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);
	  }
          return TRUE;
}

/****************************************************************************/
/*  This function is called by LeastSquares to merge together two contigs that have a containment relationship between them
 */
#define DOIT
void   ContigContainment_old(CIScaffoldT *scaffold, NodeCGW_T *prevCI, NodeCGW_T *thisCI, EdgeCGW_T *overlapEdge){
  ContigT *containingContig;
  ContigT *containedContig;
  IntElementPos contigPos;
  int mergeStatus = 0;

  if(ContigPositions == NULL){
    ContigPositions = CreateVA_IntElementPos(10);
  }
  ResetVA_IntElementPos(ContigPositions);

  if(overlapEdge->flags.bits.aContainsB){
    fprintf(GlobalData->stderrc,"* ContigContainment A (" F_CID ") contains B (" F_CID ")  prevCI = " F_CID "  thisCI = " F_CID "\n",
			overlapEdge->idA, overlapEdge->idB, prevCI->id, thisCI->id);
    if(prevCI->id == overlapEdge->idA){
#ifdef DOIT
      containingContig = prevCI;
      containedContig = thisCI;
#endif
    }else{
      fprintf(GlobalData->stderrc,"*  prevCI is contained!!!!!\n");
#ifdef DOIT
      containedContig = prevCI;
      containingContig = thisCI;
#endif
    }
  }else if(overlapEdge->flags.bits.bContainsA){
    fprintf(GlobalData->stderrc,"* ContigContainment B (" F_CID ") contains A (" F_CID ")  prevCI = " F_CID "  thisCI = " F_CID "\n",
			overlapEdge->idB, overlapEdge->idA, prevCI->id, thisCI->id);
    if(prevCI->id == overlapEdge->idB){
#ifdef DOIT
      containingContig = prevCI;
      containedContig = thisCI;
#endif
    }else{
      fprintf(GlobalData->stderrc,"*  prevCI is contained!!!!!\n");
#ifdef DOIT
      containedContig = prevCI;
      containingContig = thisCI;
#endif
    }
  }else{
    assert(0);
  }
#ifndef DOIT
  // Given the way this is called by least squares, the left most contig should be the containing contig --- NOT!
  // If least squares things it is  dovetail, and it is a contain, the order can be reversed
  containedContig = thisCI;
  containingContig = prevCI;
#endif
  fprintf(GlobalData->stderrc,"* Containing contig is " F_CID " contained contig is " F_CID "\n",
	  containingContig->id, containedContig->id);

  /* Now the ContigEnds VA contains the starts/ends of all the contigs */
  {


    double aEndPos = 0., bEndPos = 0.;
    double overlap = -overlapEdge->distance.mean;
    // all calculations assume that the containing contig is in the same direction as the scaffold
    // flip corrects for the opposite case
    int flip = (containingContig->offsetBEnd.mean < containingContig->offsetAEnd.mean) ? -1 : 1;

    switch(GetEdgeOrientationWRT(overlapEdge,containingContig->id)){
	  case AB_AB:
		//              ----------------------->    NORMAL  AB_AB
		//                ---------------->
		//                |--------------------|   overlap
		aEndPos = containingContig->offsetBEnd.mean - flip * overlap;
		bEndPos = aEndPos + flip * containedContig->bpLength.mean;
		break;
	  case BA_BA:
		//               ------------------------>    Anti-NORMAL AB_BA
		//                 --------------->
		//               |----------------|           overlap
		bEndPos = containingContig->offsetAEnd.mean + flip * overlap;
		aEndPos = bEndPos - flip * containedContig->bpLength.mean;
		break;
	  case BA_AB:
		//               ----------------------->    Outie
		//                 <---------------
		//               |----------------|          overlap
		aEndPos = containingContig->offsetAEnd.mean +  flip * overlap;
		bEndPos = aEndPos - flip * containedContig->bpLength.mean;
		break;
	  case AB_BA:
		//              ----------------------->    Innie  AB_BA
		//               <----------------
		//               |---------------------|   overlap
		bEndPos = containingContig->offsetBEnd.mean - flip * overlap;
		aEndPos = bEndPos + flip * containedContig->bpLength.mean;
		break;
          default:
            assert(0);
            break;
    }
    contigPos.ident = containingContig->id;
    contigPos.type = AS_CONTIG;
    contigPos.position.bgn = containingContig->offsetAEnd.mean;
    contigPos.position.end = containingContig->offsetBEnd.mean;
    AppendIntElementPos(ContigPositions, &contigPos);
    contigPos.ident = containedContig->id;
    contigPos.type = AS_CONTIG;
	contigPos.position.bgn = aEndPos;
    contigPos.position.end = bEndPos;
	// contigPos.position.bgn = containingContig->offsetAEnd.mean;
	// contigPos.position.end = containingContig->offsetBEnd.mean;
    AppendIntElementPos(ContigPositions, &contigPos);
	
    fprintf(GlobalData->stderrc,"* Positions:\n\t" F_CID " [%g,%g]\n\t" F_CID " [%g,%g]\n",
			containingContig->id, containingContig->offsetAEnd.mean, containingContig->offsetBEnd.mean,
			containedContig->id, aEndPos, bEndPos);
	
    flip = (containingContig->offsetBEnd.mean < containingContig->offsetAEnd.mean);
    if(flip){
      mergeStatus = CreateAContigInScaffold(scaffold, ContigPositions, containingContig->offsetBEnd, containingContig->offsetAEnd);
    }else{
      mergeStatus = CreateAContigInScaffold(scaffold, ContigPositions, containingContig->offsetAEnd, containingContig->offsetBEnd);
	}
    assert(mergeStatus == TRUE);
  }
}

#define AHANGSLOP 30

#if 0
/****************************************************************************/
/*  
	This function is meant to be called by LeastSquares to merge together two contigs that have a 
	containment relationship between them,
	but it can actually handle dovetails.  If there is no overlap, it asserts.
*/
void   notContigContainment(CIScaffoldT *scaffold, NodeCGW_T *prevCI, NodeCGW_T *thisCI, EdgeCGW_T *overlapEdge)
{
  CDS_COORD_t minAhang, maxAhang;
  int flip;
  CDS_CID_t firstContig;
  IntElementPos contigPos;
  int mergeStatus = 0;
  Overlap *contigOverlap;
  ChunkOrientationType overlapOrientation, actualOverlapOrientation;
  NodeCGW_T *leftContig, *rightContig;

  if ( min( prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean) <= 
	   min( thisCI->offsetAEnd.mean, thisCI->offsetBEnd.mean))
  {
	leftContig = prevCI;
	rightContig = thisCI;
  }
  else
  {
	leftContig = thisCI;
	rightContig = prevCI;
  }

  if(ContigPositions == NULL)
  {
    ContigPositions = CreateVA_IntElementPos(10);
  }
  ResetVA_IntElementPos(ContigPositions);

  if ( leftContig->offsetAEnd.mean < leftContig->offsetBEnd.mean)  // leftContig is AB
  {
	if ( rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean) // rightContig is AB
	{
	  overlapOrientation = AB_AB;
	  minAhang = (CDS_COORD_t) (rightContig->offsetAEnd.mean - leftContig->offsetAEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
	else // rightContig is BA
	{
	  overlapOrientation = AB_BA;
	  minAhang = (CDS_COORD_t) (rightContig->offsetBEnd.mean - leftContig->offsetAEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
  }
  else  // leftContig is BA
  {
	if ( rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean) // rightContig is AB
	{
	  overlapOrientation = BA_AB;
	  minAhang = (CDS_COORD_t) (rightContig->offsetAEnd.mean - leftContig->offsetBEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
	else // rightContig is BA
	{
	  overlapOrientation = BA_BA;
	  minAhang = (CDS_COORD_t) (rightContig->offsetBEnd.mean - leftContig->offsetBEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
  }

  contigOverlap = OverlapContigs( leftContig, rightContig, &overlapOrientation, minAhang, maxAhang, TRUE);
  if (contigOverlap == NULL)
  {
	CDS_COORD_t maxLength;
	  
	fprintf( GlobalData->stderrc, "no overlap found between " F_CID " and " F_CID " with standard AHANGSLOP, retrying...\n",
			 leftContig->id, rightContig->id);

	maxLength = max( leftContig->bpLength.mean, rightContig->bpLength.mean);
	  
	fprintf( stderr, "overlapOrientation: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n", 
			 (char) overlapOrientation, -maxLength, maxLength);

	// this is code to special case a particular contig arrangement that was happening in mouse
	// we really need to use the edge being passed in to this routine, see ContigContainment_new below
	if (0) //(leftContig->id == 3152359 && rightContig->id == 3152361)
	{
	  LengthT NullLength = {0.0, 0.0}, leftOffset, rightOffset;
	  CDS_CID_t newScaffoldID;
	  // create & init a new scaffold
	  CIScaffoldT * newScaffold = (CIScaffoldT *) safe_malloc( sizeof (CIScaffoldT));
	  assert (newScaffold != NULL);  
	  InitializeScaffold( newScaffold, REAL_SCAFFOLD);
	  newScaffold->info.Scaffold.AEndCI = NULLINDEX;
	  newScaffold->info.Scaffold.BEndCI = NULLINDEX;
	  newScaffold->info.Scaffold.numElements = 0;
	  newScaffold->edgeHead = NULLINDEX;
	  newScaffold->bpLength = NullLength;
	  newScaffoldID = newScaffold->id = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
	  newScaffold->flags.bits.isDead = FALSE;
	  newScaffold->aEndCoord = newScaffold->bEndCoord = -1;
	  newScaffold->numEssentialA = newScaffold->numEssentialB = 0;
	  newScaffold->essentialEdgeB = newScaffold->essentialEdgeA = NULLINDEX;
	  AppendGraphNode( ScaffoldGraph->ScaffoldGraph, newScaffold);

	  leftOffset.mean = 0.0;
	  leftOffset.variance = 0.0;
	  rightOffset.mean = leftContig->bpLength.mean;
	  rightOffset.variance = ComputeFudgeVariance(leftContig->bpLength.mean);
	  
	  RemoveCIFromScaffold( ScaffoldGraph, scaffold, leftContig, FALSE);
	  InsertCIInScaffold( ScaffoldGraph, leftContig->id, newScaffoldID, 
						  leftOffset, rightOffset, TRUE, FALSE);
	  
	  leftOffset.mean = 1590.0;
	  leftOffset.variance = ComputeFudgeVariance( leftOffset.mean );
	  rightOffset.mean = leftOffset.mean + rightContig->bpLength.mean;
	  rightOffset.variance = leftOffset.variance + ComputeFudgeVariance(rightContig->bpLength.mean);;

	  RemoveCIFromScaffold( ScaffoldGraph, scaffold, rightContig, FALSE);
	  InsertCIInScaffold( ScaffoldGraph, rightContig->id, newScaffoldID, 
						  leftOffset, rightOffset, TRUE, FALSE);

	  return;
	}
	
	contigOverlap = OverlapContigs( leftContig, rightContig, &overlapOrientation, -maxLength, maxLength, FALSE);

	if (contigOverlap == NULL)
	{
	  fprintf( GlobalData->stderrc, "no overlap found between " F_CID " and " F_CID " with maximum slop, aborting...\n",
			   leftContig->id, rightContig->id);
	  dumpContigInfo(leftContig);
	  dumpContigInfo(rightContig);
	  assert(0);
	}	
  }

  if (contigOverlap->begpos < 0) // contigs need to be reversed
  {
	leftContig = thisCI;
	rightContig = prevCI;
	// adjust Overlap fields for later use in positioning
	contigOverlap->begpos = - contigOverlap->begpos;
	contigOverlap->endpos = - contigOverlap->endpos;	

	switch(overlapOrientation){
	case AB_AB:
	case BA_BA:
	  actualOverlapOrientation = overlapOrientation;	  
	  break;
	case AB_BA:
	  actualOverlapOrientation = BA_AB;
	  break;
	case BA_AB:
	  actualOverlapOrientation = AB_BA;
	  break;
	default:
	  assert(0);
	}
	fprintf(GlobalData->stderrc,"* Switched right-left  orientation went from %c to %c\n",
		overlapOrientation, actualOverlapOrientation);
  }else{

    actualOverlapOrientation = overlapOrientation;
  }
  
  fprintf(GlobalData->stderrc,"* Containing contig is " F_CID " contained contig is " F_CID " ahg:%d bhg:%d orient:%c\n",
		  leftContig->id, rightContig->id, contigOverlap->begpos, contigOverlap->endpos, actualOverlapOrientation);

  fprintf(GlobalData->stderrc,"* Initial Positions:\n\t" F_CID " [%g,%g]\n\t" F_CID " [%g,%g]\n",
		  leftContig->id, 
	          leftContig->offsetAEnd.mean, leftContig->offsetBEnd.mean,
		  rightContig->id, 
	          rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);

  // assume we leave the leftContig where it is
  if ( actualOverlapOrientation == AB_AB)
  {
	rightContig->offsetAEnd.mean = leftContig->offsetAEnd.mean + contigOverlap->begpos;
	rightContig->offsetBEnd.mean = rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
  }
  else if ( actualOverlapOrientation == AB_BA)
  {
	rightContig->offsetBEnd.mean = leftContig->offsetAEnd.mean + contigOverlap->begpos;
	rightContig->offsetAEnd.mean = rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
  }
  else if ( actualOverlapOrientation == BA_AB)
  {
	rightContig->offsetAEnd.mean = leftContig->offsetBEnd.mean + contigOverlap->begpos;
	rightContig->offsetBEnd.mean = rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
  }
  if ( actualOverlapOrientation == BA_BA)
  {
	rightContig->offsetBEnd.mean = leftContig->offsetBEnd.mean + contigOverlap->begpos;
	rightContig->offsetAEnd.mean = rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
  }

  contigPos.ident = leftContig->id;
  contigPos.type = AS_CONTIG;
  contigPos.position.bgn = leftContig->offsetAEnd.mean;
  contigPos.position.end = leftContig->offsetBEnd.mean;
  AppendIntElementPos(ContigPositions, &contigPos);
  contigPos.ident = rightContig->id;
  contigPos.type = AS_CONTIG;
  contigPos.position.bgn = rightContig->offsetAEnd.mean;
  contigPos.position.end = rightContig->offsetBEnd.mean;
  AppendIntElementPos(ContigPositions, &contigPos);
  
  fprintf(GlobalData->stderrc,"* Final Positions:\n\t" F_CID " [%g,%g]\n\t" F_CID " [%g,%g]\n",
		  leftContig->id, 
	          leftContig->offsetAEnd.mean, leftContig->offsetBEnd.mean,
		  rightContig->id, 
	          rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);
  
  flip = (leftContig->offsetBEnd.mean < leftContig->offsetAEnd.mean);
  if(flip)
  {
	mergeStatus = CreateAContigInScaffold(scaffold, ContigPositions, leftContig->offsetBEnd, leftContig->offsetAEnd);
  }
  else
  {
	mergeStatus = CreateAContigInScaffold(scaffold, ContigPositions, leftContig->offsetAEnd, leftContig->offsetBEnd);
  }
  assert(mergeStatus == TRUE);
}
#endif

// this version tries several ways of finding the overlap relationsship if tryHarder is set
void ContigContainment(CIScaffoldT *scaffold, NodeCGW_T *prevCI, NodeCGW_T *thisCI, EdgeCGW_T *overlapEdge,
					   int tryHarder)
{
  CDS_COORD_t minAhang, maxAhang;
  int flip;
  IntElementPos contigPos;
  int mergeStatus = 0;
  Overlap *contigOverlap;
  ChunkOrientationType overlapOrientation, actualOverlapOrientation;
  NodeCGW_T *leftContig, *rightContig;
  static VA_TYPE(IntElementPos) *ContigPositions = NULL;

  if(ContigPositions == NULL)
  {
    ContigPositions = CreateVA_IntElementPos(10);
  }
  ResetVA_IntElementPos(ContigPositions);

  if ( min( prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean) <= 
	   min( thisCI->offsetAEnd.mean, thisCI->offsetBEnd.mean))
  {
	leftContig = prevCI;
	rightContig = thisCI;
  }
  else
  {
	leftContig = thisCI;
	rightContig = prevCI;
  }

  if(ContigPositions == NULL)
  {
    ContigPositions = CreateVA_IntElementPos(10);
  }
  ResetVA_IntElementPos(ContigPositions);

  if ( leftContig->offsetAEnd.mean < leftContig->offsetBEnd.mean)  // leftContig is AB
  {
	if ( rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean) // rightContig is AB
	{
	  overlapOrientation = AB_AB;
	  minAhang = (CDS_COORD_t) (rightContig->offsetAEnd.mean - leftContig->offsetAEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
	else // rightContig is BA
	{
	  overlapOrientation = AB_BA;
	  minAhang = (CDS_COORD_t) (rightContig->offsetBEnd.mean - leftContig->offsetAEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
  }
  else  // leftContig is BA
  {
	if ( rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean) // rightContig is AB
	{
	  overlapOrientation = BA_AB;
	  minAhang = (CDS_COORD_t) (rightContig->offsetAEnd.mean - leftContig->offsetBEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
	else // rightContig is BA
	{
	  overlapOrientation = BA_BA;
	  minAhang = (CDS_COORD_t) (rightContig->offsetBEnd.mean - leftContig->offsetBEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
  }

  fprintf( stderr, "calling OverlapContigs with orient: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n",
		   overlapOrientation, minAhang, maxAhang);
  contigOverlap = OverlapContigs( leftContig, rightContig, &overlapOrientation, minAhang, maxAhang, TRUE);

  if ( tryHarder )
  {
	// if no overlap found, try flipping orientation in case DPCompare asymmetry is biting us
	if (contigOverlap == NULL)
	{
	  fprintf( stderr, "no overlap found between " F_CID " and " F_CID ", retrying with flipped orientation\n",
			   leftContig->id, rightContig->id);
	  
	  // try with the reverse orientation
	  overlapOrientation = InvertEdgeOrient( (const ChunkOrientationType) overlapOrientation );
	  
	  contigOverlap = OverlapContigs( leftContig, rightContig, &overlapOrientation, minAhang, maxAhang, TRUE);
	  
	  if ( contigOverlap != NULL)
	  {
		CDS_COORD_t temp;
		
		temp = -contigOverlap->begpos;
		contigOverlap->begpos = -contigOverlap->endpos;
		contigOverlap->endpos = temp;
	  }
	  // restore the orientation
	  overlapOrientation = InvertEdgeOrient( (const ChunkOrientationType) overlapOrientation );
	}
	
	// if still no overlap found, try maxing out hangs
	if (contigOverlap == NULL)
	{
	  CDS_COORD_t maxLength;
	  
	  fprintf( stderr, "no overlap found between " F_CID " and " F_CID ", retrying with max AHANGSLOP\n",
			   leftContig->id, rightContig->id);
	  
	  maxLength = max( leftContig->bpLength.mean, rightContig->bpLength.mean);
	  
	  fprintf( stderr, "overlapOrientation: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n", 
			   (char) overlapOrientation, -maxLength, maxLength);
	  
	  contigOverlap = OverlapContigs( leftContig, rightContig, &overlapOrientation, -maxLength, maxLength, FALSE);	  
	}
  
	// if still no overlap found, try flipping orientation and maxing out hangs
	if (contigOverlap == NULL)
	{
	  CDS_COORD_t maxLength;
	  
	  fprintf( stderr, 
			   "no overlap found between " F_CID " and " F_CID ", retrying with flipped orientation and max AHANGSLOP\n",
			   leftContig->id, rightContig->id);
	  
	  maxLength = max( leftContig->bpLength.mean, rightContig->bpLength.mean);
	  
	  // try with the reverse orientation
	  overlapOrientation = InvertEdgeOrient( (const ChunkOrientationType) overlapOrientation );
	  
	  fprintf( stderr, "overlapOrientation: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n", 
			   (char) overlapOrientation, -maxLength, maxLength);
	  
	  contigOverlap = OverlapContigs( leftContig, rightContig, &overlapOrientation, -maxLength, maxLength, FALSE);
	  
	  if ( contigOverlap != NULL)
	  {
		CDS_COORD_t temp;
		
		temp = -contigOverlap->begpos;
		contigOverlap->begpos = -contigOverlap->endpos;
		contigOverlap->endpos = temp;
	  }
	  
	  // restore the orientation
	  overlapOrientation = InvertEdgeOrient( (const ChunkOrientationType) overlapOrientation );
	}
  }

  if (contigOverlap == NULL)
  {
	fprintf( stderr, "no overlap found between " F_CID " and " F_CID ", aborting...\n",
			 leftContig->id, rightContig->id);
	dumpContigInfo(leftContig);
	dumpContigInfo(rightContig);
	assert(0);
  }	
  
  if (contigOverlap->begpos < 0) // contigs need to be reversed
  {
	leftContig = thisCI;
	rightContig = prevCI;
	// adjust Overlap fields for later use in positioning
	contigOverlap->begpos = - contigOverlap->begpos;
	contigOverlap->endpos = - contigOverlap->endpos;	

	switch(overlapOrientation)
	{
	  case AB_AB:
	  case BA_BA:
		actualOverlapOrientation = overlapOrientation;	  
		break;
	  case AB_BA:
		actualOverlapOrientation = BA_AB;
		break;
	  case BA_AB:
		actualOverlapOrientation = AB_BA;
		break;
	  default:
		assert(0);
	}
	fprintf( stderr, "* Switched right-left  orientation went from %c to %c\n",
			 overlapOrientation, actualOverlapOrientation);
  }
  else
  {
    actualOverlapOrientation = overlapOrientation;
  }
  
  fprintf( stderr, "* Containing contig is " F_CID " contained contig is " F_CID " ahg:%d bhg:%d orient:%c\n",
		   leftContig->id, rightContig->id, contigOverlap->begpos, contigOverlap->endpos, actualOverlapOrientation);

  fprintf( stderr, "* Initial Positions:\n\t" F_CID " [%g,%g]\n\t" F_CID " [%g,%g]\n",
		   leftContig->id, 
		   leftContig->offsetAEnd.mean, leftContig->offsetBEnd.mean,
		   rightContig->id, 
		   rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);

  // assume we leave the leftContig where it is
  if ( actualOverlapOrientation == AB_AB)
  {
	rightContig->offsetAEnd.mean = leftContig->offsetAEnd.mean + contigOverlap->begpos;
	rightContig->offsetBEnd.mean = rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
  }
  else if ( actualOverlapOrientation == AB_BA)
  {
	rightContig->offsetBEnd.mean = leftContig->offsetAEnd.mean + contigOverlap->begpos;
	rightContig->offsetAEnd.mean = rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
  }
  else if ( actualOverlapOrientation == BA_AB)
  {
	rightContig->offsetAEnd.mean = leftContig->offsetBEnd.mean + contigOverlap->begpos;
	rightContig->offsetBEnd.mean = rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
  }
  if ( actualOverlapOrientation == BA_BA)
  {
	rightContig->offsetBEnd.mean = leftContig->offsetBEnd.mean + contigOverlap->begpos;
	rightContig->offsetAEnd.mean = rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
  }

  contigPos.ident = leftContig->id;
  contigPos.type = AS_CONTIG;
  contigPos.position.bgn = leftContig->offsetAEnd.mean;
  contigPos.position.end = leftContig->offsetBEnd.mean;
  AppendIntElementPos(ContigPositions, &contigPos);
  contigPos.ident = rightContig->id;
  contigPos.type = AS_CONTIG;
  contigPos.position.bgn = rightContig->offsetAEnd.mean;
  contigPos.position.end = rightContig->offsetBEnd.mean;
  AppendIntElementPos(ContigPositions, &contigPos);
  
  fprintf( stderr,"* Final Positions:\n\t" F_CID " [%g,%g]\n\t" F_CID " [%g,%g]\n",
		   leftContig->id, 
		   leftContig->offsetAEnd.mean, leftContig->offsetBEnd.mean,
		   rightContig->id, 
		   rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);
  
  flip = (leftContig->offsetBEnd.mean < leftContig->offsetAEnd.mean);
  if(flip)
  {
	mergeStatus = CreateAContigInScaffold(scaffold, ContigPositions, leftContig->offsetBEnd, leftContig->offsetAEnd);
  }
  else
  {
	mergeStatus = CreateAContigInScaffold(scaffold, ContigPositions, leftContig->offsetAEnd, leftContig->offsetBEnd);
  }
  assert(mergeStatus == TRUE);
}


/******************* Code in progress below **********************************/
/*****************************************************************************/
/*  
	This function is meant to be called by LeastSquares to merge together two contigs that have a 
	containment relationship between them,
	but it can actually handle dovetails.  If there is no overlap, it asserts.
	This version uses the overlapEdge information passed in rather than recomputing overlaps.
*/
#ifdef NEVER
void   ContigContainment_new(CIScaffoldT *scaffold,
                             NodeCGW_T *prevCI,
                             NodeCGW_T *thisCI,
                             EdgeCGW_T *overlapEdge)
{
  CDS_COORD_t minAhang, maxAhang;
  int flip;
  CDS_CID_t firstContig;
  IntElementPos contigPos;
  int mergeStatus = 0;
  Overlap *contigOverlap;
  ChunkOrientationType overlapOrientation, actualOverlapOrientation;
  NodeCGW_T *leftContig, *rightContig;

  fprintf( GlobalData->stderrc, "1 prevCI->id: " F_CID " (%f, %f) and thisCI->id: " F_CID " (%f, %f)\n",
		   prevCI->id, prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean,
		   thisCI->id, thisCI->offsetAEnd.mean, thisCI->offsetBEnd.mean);

  if ( overlapEdge->idA == prevCI->id )
  {
	leftContig = prevCI;
	rightContig = thisCI;
  }
  else
  {
	leftContig = thisCI;
	rightContig = prevCI;
  }

  fprintf( GlobalData->stderrc, "1 leftContig->id: " F_CID " and rightContig->id: " F_CID "\n",
		   leftContig->id, rightContig->id);

  if(ContigPositions == NULL)
  {
    ContigPositions = CreateVA_IntElementPos(10);
  }
  ResetVA_IntElementPos(ContigPositions);

  if ( leftContig->offsetAEnd.mean < leftContig->offsetBEnd.mean)  // leftContig is AB
  {
	if ( rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean) // rightContig is AB
	{
	  overlapOrientation = AB_AB;
	  minAhang = (CDS_COORD_t) (rightContig->offsetAEnd.mean - leftContig->offsetAEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
	else // rightContig is BA
	{
	  overlapOrientation = AB_BA;
	  minAhang = (CDS_COORD_t) (rightContig->offsetBEnd.mean - leftContig->offsetAEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
  }
  else  // leftContig is BA
  {
	if ( rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean) // rightContig is AB
	{
	  overlapOrientation = BA_AB;
	  minAhang = (CDS_COORD_t) (rightContig->offsetAEnd.mean - leftContig->offsetBEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
	else // rightContig is BA
	{
	  overlapOrientation = BA_BA;
	  minAhang = (CDS_COORD_t) (rightContig->offsetBEnd.mean - leftContig->offsetBEnd.mean) - AHANGSLOP;
	  maxAhang = minAhang + (2 * AHANGSLOP);
	}
  }

  fprintf( stderr, "leftContig: " F_CID ", rightContig: " F_CID ", overlapOrientation: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n", 
		   leftContig->id, rightContig->id, (char) overlapOrientation, minAhang, maxAhang);
  contigOverlap = OverlapContigs( leftContig, rightContig, &overlapOrientation, minAhang, maxAhang, TRUE);
  if (contigOverlap == NULL)
  {
	CDS_COORD_t maxLength;
	  
	fprintf( GlobalData->stderrc, "no overlap found between " F_CID " and " F_CID " with standard AHANGSLOP, retrying...\n",
			 leftContig->id, rightContig->id);
	// minAhang = - rightContig->bpLength.mean;
	// maxAhang = leftContig->bpLength.mean;
	
	maxLength = max( leftContig->bpLength.mean, rightContig->bpLength.mean);
	  
	fprintf( stderr, "overlapOrientation: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n", 
			 (char) overlapOrientation, -maxLength, maxLength);

	if (leftContig->id == 3152359 && rightContig->id == 3152361)
	{
	  contigOverlap = (Overlap *) safe_malloc( sizeof (Overlap));
	  if (contigOverlap == NULL)
	  {
		fprintf( stderr, "failed to safe_malloc Overlap struct\n");
		assert(0);
	  }
	  contigOverlap->begpos = 1590;
	}
	else
	  contigOverlap = OverlapContigs( leftContig, rightContig, &overlapOrientation, -maxLength, maxLength, FALSE);

	if (contigOverlap == NULL)
	{
	  fprintf( GlobalData->stderrc, "no overlap found between " F_CID " and " F_CID " with maximum slop, aborting...\n",
			   leftContig->id, rightContig->id);
	  dumpContigInfo(leftContig);
	  dumpContigInfo(rightContig);
	  assert(0);
	}	
  }

  fprintf( GlobalData->stderrc, "2 leftContig->id: " F_CID " and rightContig->id: " F_CID "\n",
		   leftContig->id, rightContig->id);

  if (contigOverlap->begpos < 0) // contigs need to be reversed
  {
	// NodeCGW_T *tempCI;
	// tempCI = leftContig;
	// leftContig = rightContig;
	// rightContig = tempCI;
	
	leftContig = thisCI;
	rightContig = prevCI;
	// adjust Overlap fields for later use in positioning
	contigOverlap->begpos = - contigOverlap->begpos;
	contigOverlap->endpos = - contigOverlap->endpos;	

	switch(overlapOrientation){
	case AB_AB:
	case BA_BA:
	  actualOverlapOrientation = overlapOrientation;	  
	  break;
	case AB_BA:
	  actualOverlapOrientation = BA_AB;
	  break;
	case BA_AB:
	  actualOverlapOrientation = AB_BA;
	  break;
	default:
	  assert(0);
	}
	fprintf(GlobalData->stderrc,"* Switched right-left  orientation went from %c to %c\n",
		overlapOrientation, actualOverlapOrientation);
  }else{

    actualOverlapOrientation = overlapOrientation;
  }
  
  fprintf( GlobalData->stderrc, "1 leftContig->id: " F_CID " and rightContig->id: " F_CID "\n",
		   leftContig->id, rightContig->id);

  fprintf(GlobalData->stderrc,"* Containing contig is " F_CID " contained contig is " F_CID " ahg:" F_COORD " bhg:" F_COORD " orient:%c\n",
		  leftContig->id, rightContig->id, contigOverlap->begpos, contigOverlap->endpos, actualOverlapOrientation);

  fprintf(GlobalData->stderrc,"* Initial Positions:\n\t" F_CID " [%g,%g]\n\t" F_CID " [%g,%g]\n",
		  leftContig->id, 
	          leftContig->offsetAEnd.mean, leftContig->offsetBEnd.mean,
		  rightContig->id, 
	          rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);

  // assume we leave the leftContig where it is
  if ( actualOverlapOrientation == AB_AB)
  {
	rightContig->offsetAEnd.mean = leftContig->offsetAEnd.mean + contigOverlap->begpos;
	rightContig->offsetBEnd.mean = rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
  }
  else if ( actualOverlapOrientation == AB_BA)
  {
	rightContig->offsetBEnd.mean = leftContig->offsetAEnd.mean + contigOverlap->begpos;
	rightContig->offsetAEnd.mean = rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
  }
  else if ( actualOverlapOrientation == BA_AB)
  {
	rightContig->offsetAEnd.mean = leftContig->offsetBEnd.mean + contigOverlap->begpos;
	rightContig->offsetBEnd.mean = rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
  }
  if ( actualOverlapOrientation == BA_BA)
  {
	rightContig->offsetBEnd.mean = leftContig->offsetBEnd.mean + contigOverlap->begpos;
	rightContig->offsetAEnd.mean = rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
  }

  contigPos.ident = leftContig->id;
  contigPos.type = AS_CONTIG;
  contigPos.position.bgn = leftContig->offsetAEnd.mean;
  contigPos.position.end = leftContig->offsetBEnd.mean;
  AppendIntElementPos(ContigPositions, &contigPos);
  contigPos.ident = rightContig->id;
  contigPos.type = AS_CONTIG;
  contigPos.position.bgn = rightContig->offsetAEnd.mean;
  contigPos.position.end = rightContig->offsetBEnd.mean;
  AppendIntElementPos(ContigPositions, &contigPos);
  
  fprintf(GlobalData->stderrc,"* Final Positions:\n\t" F_CID " [%g,%g]\n\t" F_CID " [%g,%g]\n",
		  leftContig->id, 
	          leftContig->offsetAEnd.mean, leftContig->offsetBEnd.mean,
		  rightContig->id, 
	          rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);
  
  flip = (leftContig->offsetBEnd.mean < leftContig->offsetAEnd.mean);
  if(flip)
  {
	mergeStatus = CreateAContigInScaffold(scaffold, ContigPositions, leftContig->offsetBEnd, leftContig->offsetAEnd);
  }
  else
  {
	mergeStatus = CreateAContigInScaffold(scaffold, ContigPositions, leftContig->offsetAEnd, leftContig->offsetBEnd);
  }
  assert(mergeStatus == TRUE);
}
#endif



/****************** DEAD CODE BELOW *****************************************/
/****************************************************************************/
#if 0
int BuildContigs(ScaffoldGraphT *graph){
  CDS_CID_t sid;
  /* For each chunk */
  for(sid = 0; sid < GetNumCIScaffoldTs(graph->CIScaffolds); sid++){
    CIScaffoldT *scaffold = GetGraphNode(graph->ScaffoldGraph,sid);
    if(isDeadCIScaffoldT(scaffold))
      continue;
    ContigAScaffold(graph, sid);
  }
  graph->numContigs = GetNumGraphNodes(graph->ContigGraph);
  return TRUE;
}

/***************************************************************************/
// Build contigs for a single scaffold
// ASSUMPTIONS:
//    If the mean positions of two CIs place them overlapping by > Threshhold,
//    then, they indeed overlap and should be contigged
//
int ContigAScaffold(ScaffoldGraphT *graph, CDS_CID_t sid){

    CIScaffoldTIterator CIs;
    CIScaffoldT *scaffold = GetGraphNode(graph->ScaffoldGraph,sid);
    ChunkInstanceT *CI;
    ChunkInstanceT *currCI = NULL;
    ChunkInstanceT *prevCI = NULL;
    ContigEndsT contig;
    int i;
    int32 initialNumElements = scaffold->info.Scaffold.numElements;
    // fprintf(GlobalData->stderrc,"* ContigAScaffold " F_CID "\n", sid);

    if(ContigEnds == NULL){
      ContigEnds = CreateVA_ContigEndsT(10);
    }else{
      ResetContigEndsT(ContigEnds);
    }

    contig.firstCID = contig.lastCID = NULLINDEX;
    contig.firstCI = contig.lastCI = NULL;

    /* First, figure out where the contigs are.  We can't contig them right away, since that
       would be mucking with the data structure we're iterating over */
    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
    /* Initialize CI . . . ELA */
    contig.firstCI = currCI = CI = NextCIScaffoldTIterator(&CIs);
    contig.firstCID = contig.firstCI->id; //GetVAIndex_ChunkInstanceT(graph->ChunkInstances, contig.firstCI);
    AssertPtr(CI);		// There should be at least one CI . . ELA
    if(currCI->offsetAEnd.mean < currCI->offsetBEnd.mean){
      contig.minOffset = currCI->offsetAEnd;
      contig.maxOffset = currCI->offsetBEnd;
    }else{
      contig.maxOffset = currCI->offsetAEnd;
      contig.minOffset = currCI->offsetBEnd;
    }
    contig.count = 0;
    while(CI = NextCIScaffoldTIterator(&CIs)){
      prevCI = currCI;
      currCI = CI;

      contig.count++;
      contig.lastCI = prevCI;
      contig.lastCID = contig.lastCI->id; //GetVAIndex_ChunkInstanceT(graph->ChunkInstances, contig.lastCI);
      // The first CI in the contig must have the minimum offset, by definition
      // subsequent CIs are sorted by LEFTmost coordinate, so we need to track
      // the MAX span (rightmost coordinate)

      {
	CDS_COORD_t actual;

	//#define DEBUG_CONTIG
#ifdef DEBUG_CONTIG
	fprintf(GlobalData->stderrc,"* contig.min = %g contig.max = %g\n",
		contig.minOffset.mean, contig.maxOffset.mean);
#endif
	actual = IntervalsOverlap(contig.minOffset.mean, contig.maxOffset.mean,
                                  currCI->offsetAEnd.mean, currCI->offsetBEnd.mean, -15000);

	if(actual> CGW_DP_MINLEN){
#ifdef DEBUG_CONTIG
	  fprintf(GlobalData->stderrc,
	       "* CI " F_CID " and " F_CID " mean positions overlap by " F_COORD " (%g,%g) (%g,%g)\n",
		prevCI->id,currCI->id, actual, contig.minOffset.mean, 
		contig.maxOffset.mean, currCI->offsetAEnd.mean, 
		currCI->offsetBEnd.mean);
#endif

	  if(currCI->offsetAEnd.mean > contig.maxOffset.mean){
	    contig.maxOffset = currCI->offsetAEnd;
	  }
	  if(currCI->offsetBEnd.mean > contig.maxOffset.mean){
	    contig.maxOffset = currCI->offsetBEnd;
	  }
	}else{
#ifdef DEBUG_CONTIG
	  fprintf(GlobalData->stderrc,"* CI " F_CID " and " F_CID " mean positions GAP by " F_COORD " (%g,%g) (%g,%g) first " F_CID " last " F_CID " \n",
		  prevCI->id,currCI->id, actual, 
		  contig.minOffset.mean,  contig.maxOffset.mean, 
		  currCI->offsetAEnd.mean, currCI->offsetBEnd.mean, 
		  contig.firstCI->id, contig.lastCI->id);
#endif
	  AppendContigEndsT(ContigEnds, &contig);
	  contig.count = 0;
	  contig.firstCI = currCI;
	  contig.firstCID = contig.firstCI->id; //GetVAIndex_ChunkInstanceT(graph->ChunkInstances, contig.firstCI);

	  if(currCI->offsetAEnd.mean < currCI->offsetBEnd.mean){
	    contig.minOffset = currCI->offsetAEnd;
	    contig.maxOffset = currCI->offsetBEnd;
	  }else{
	    contig.maxOffset = currCI->offsetAEnd;
	    contig.minOffset = currCI->offsetBEnd;
	  }
#ifdef DEBUG_CONTIG
	  fprintf(GlobalData->stderrc,"* Reseting contig  firstCI: " F_CID " min:%g max:%g\n",
	    contig.firstCI->id, contig.minOffset.mean, contig.maxOffset.mean);
#endif
	}
      }
    }
    /* Add last contig to the list */
    ++contig.count;
    contig.lastCI = currCI;
    contig.lastCID = contig.lastCI->id; //GetVAIndex_ChunkInstanceT(graph->ChunkInstances, contig.lastCI);
    AppendContigEndsT(ContigEnds, &contig);

    {
      ContigEndsT *ctg = GetContigEndsT(ContigEnds,0);
      ContigT newContig;
      CDS_CID_t firstContig;
      int32 numContigs;

      InitializeContig(&newContig, CONTIG_CGW);
      firstContig = GetNumGraphNodes(graph->ContigGraph);
      numContigs = GetNumContigEndsTs(ContigEnds);

      newContig.flags.bits.isUnique = (scaffold->type == REAL_SCAFFOLD?TRUE:FALSE);
      newContig.flags.bits.cgbType = UU_CGBTYPE;
      newContig.aEndCoord = newContig.bEndCoord = -1;
      newContig.scaffoldID = ctg->firstCI->scaffoldID;
      assert(newContig.scaffoldID != NULLINDEX);
      newContig.edgeHead = NULLINDEX;
      newContig.microhetScore = NULLINDEX;

    /* Now the ContigEnds VA contains the starts/ends of all the contigs */
      for(i = 0; i < numContigs; i++){
	ctg = GetContigEndsT(ContigEnds,i);

	newContig.offsetAEnd = ctg->minOffset;
	newContig.offsetBEnd = ctg->maxOffset;
	newContig.aEndCoord = newContig.bEndCoord = -1;
	newContig.id = firstContig + i;

	// Maybe things have moved by reallocation, refresh this pointer!
	ctg->firstCI = GetGraphNode(graph->RezGraph, ctg->firstCID);
	ctg->lastCI = GetGraphNode(graph->RezGraph, ctg->lastCID);

	if(ctg->firstCI->flags.bits.cgbType == UU_CGBTYPE){
	  if(GetNodeOrient(ctg->firstCI) == A_B){
	    newContig.aEndCoord = ctg->firstCI->aEndCoord;
	  }else{
	    newContig.aEndCoord = ctg->firstCI->bEndCoord;
	  }
	}
	if(ctg->lastCI->flags.bits.cgbType == UU_CGBTYPE){
	  if(GetNodeOrient(ctg->lastCI) == A_B){
	    newContig.bEndCoord = ctg->lastCI->bEndCoord;
	  }else{
	    newContig.bEndCoord = ctg->lastCI->aEndCoord;
	  }
	}

	if(ctg->firstCI->flags.bits.cgbType == UU_CGBTYPE &&
	   ctg->lastCI->flags.bits.cgbType == UU_CGBTYPE){
	  newContig.flags.bits.cgbType = UU_CGBTYPE;
	}else{
	  newContig.flags.bits.cgbType = RR_CGBTYPE;
	}

        ComputeLength(&newContig.bpLength, &(ctg->minOffset), &(ctg->maxOffset));
#ifdef DEBUG_CONTIG
	fprintf(GlobalData->stderrc,"*** newContig.bpLength = %g endLength:%g contig.length = %g (%g,%g)\n",
		newContig.bpLength.mean, 
		newContig.offsetBEnd.mean -newContig.offsetAEnd.mean,
		ctg->maxOffset.mean - ctg->minOffset.mean,
		ctg->minOffset.mean,ctg->maxOffset.mean);
#endif
	// Connect the contig to its consitutent CIs
        newContig.info.Contig.AEndCI = ctg->firstCID;
        newContig.info.Contig.BEndCI = ctg->lastCID;
	newContig.info.Contig.numCI = ctg->count;


	// Disconnect the CIs from their scaffold links
	ctg->firstCI->AEndNext = NULLINDEX;
	ctg->lastCI->BEndNext = NULLINDEX;

	// Connect the contigs within a scaffold with scaffold links 
	if(i != 0)
	  newContig.AEndNext = newContig.id - 1;
	else
	  newContig.AEndNext = NULLINDEX;

	if(i != numContigs - 1)
	  newContig.BEndNext = newContig.id + 1;
	else
	  newContig.BEndNext = NULLINDEX;
	  

#ifdef DEBUG_CONTIG
	fprintf(GlobalData->stderrc,"* Contig " F_CID " at offset %g,%g has length %g +/- %g  nextA:" F_CID "  nextB:" F_CID "  firstCI:" F_CID "  lastCI:" F_CID "\n",
		newContig.id, 
		newContig.offsetAEnd.mean, 
		newContig.offsetBEnd.mean,
		newContig.bpLength.mean, newContig.bpLength.variance,
		newContig.AEndNext, newContig.BEndNext,
		newContig.info.Contig.AEndCI, newContig.info.Contig.BEndCI);
#endif
	AppendGraphNode(graph->ContigGraph, &newContig);

	{
	  ContigTIterator CIs;
	  ChunkInstanceT *CI;
	  int cnt = 0;
	  //#define DEBUG_SCAFFOLD
#ifdef DEBUG_SCAFFOLD
	  fprintf(GlobalData->stderrc,"* Updating CI offsets for Contig " F_CID " relative to (%d,%d)\n",
		  newContig.id, (int) newContig.offsetAEnd.mean, (int) newContig.offsetBEnd.mean);
#endif
	/* Now iterate through the CIs in the Contig, and mark them as belonging to teh contig */
	InitContigTIterator(graph, newContig.id, TRUE, FALSE, &CIs);
	while(CI = NextContigTIterator(&CIs)){

#ifdef DEBUG_SCAFFOLD
	  fprintf(GlobalData->stderrc,"* Marking CI " F_CID " starting at scaffold offsets (%d,%d) contig offset (%d,%d) anext:" F_CID " bnext:" F_CID "\n",
		  CI->id, 
		  (int) CI->offsetAEnd.mean,
		  (int) CI->offsetBEnd.mean,
		  (int) newContig.offsetAEnd.mean,
		  (int) newContig.offsetBEnd.mean,
		  CI->AEndNext, CI->BEndNext);
#endif
	  cnt++;
	  CI->info.CI.contigID = newContig.id;
	  CI->offsetAEnd.mean -= newContig.offsetAEnd.mean;
	  CI->offsetAEnd.variance -= newContig.offsetAEnd.variance;

	  CI->offsetBEnd.mean -= newContig.offsetAEnd.mean;
	  CI->offsetBEnd.variance -= newContig.offsetAEnd.variance;
	  if( CI->offsetAEnd.variance < 0.0){
	    fprintf(GlobalData->stderrc,"* Negative variance on CI (" F_CID ") ->offsetAEnd (%g)\n",
		    CI->id, CI->offsetAEnd.variance);
	    CI->offsetAEnd.variance = 0.0;
	  }
	  fflush(GlobalData->stderrc);
	  if( CI->offsetBEnd.variance < 0.0){
	    fprintf(GlobalData->stderrc,"* Negative variance on CI (" F_CID ") ->offsetBEnd (%g)\n",
		    CI->id, CI->offsetBEnd.variance);
	    CI->offsetBEnd.variance = 0.0;
	  }
	  fflush(GlobalData->stderrc);

#ifdef DEBUG_SCAFFOLD
	  fprintf(GlobalData->stderrc,"* Marking CI " F_CID " as belonging to contig " F_CID " at offsets (%d,%d)\n",
		  CI->id, newContig.id,
		  (int) CI->offsetAEnd.mean,
		  (int) CI->offsetBEnd.mean );
#endif
	}
	  assert(cnt ==  newContig.info.Contig.numCI);

	}

      }
	/* Reset the scaffold, and set it up to hold the collection of contigs */
      scaffold->info.Scaffold.numElements = numContigs;
      scaffold->info.Scaffold.AEndCI = firstContig;
      scaffold->info.Scaffold.BEndCI = firstContig + numContigs - 1;
      scaffold->flags.bits.containsCIs = FALSE;

#ifdef DEBUG
    fprintf(GlobalData->stderrc,"* Scaffold " F_CID " had %d contigs %d CIs aend:" F_CID " bend:" F_CID "\n",
	    sid, (int) GetNumContigEndsTs(ContigEnds), scaffold->numElements,
	    scaffold->info.Scaffold.AEndCI, scaffold->info.Scaffold.BEndCI);
#endif
    }
    return TRUE;
}
#endif

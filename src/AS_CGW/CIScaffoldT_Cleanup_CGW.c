
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
static char CM_ID[] = "$Id: CIScaffoldT_Cleanup_CGW.c,v 1.17 2006-11-14 17:52:14 eliv Exp $";

#undef DEBUG_CHECKFORCTGS
#undef DEBUG_DETAILED
#undef DEBUG_CONNECTEDNESS


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

VA_DEF(ContigEndsT);
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

      distance.variance = MAX(edge->distance.variance, 10.0);

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

      distance1.variance = MAX(edge->distance.variance, 10.0);
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

	//fprintf(GlobalData->stderrc,"* i = %d eid = " F_CID "\n", i,eid);

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

    // Compute the new edge orientation
    // First get edge orientation into canonical form


    newOrient = GetEdgeOrientationWRT(newEdge, contigID);

    //  If the overlap is moving from the A End of the old contig to
    //  the B End of the new contig, we need to adjust the orientation
    //  of the edge
    //
    if (contigEnd != newContigEnd) {
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

    // Modify the ida and idb to reflect the new contig
    {
      CDS_CID_t newEdgeID = GetVAIndex_EdgeCGW_T(ScaffoldGraph->RezGraph->edges, newEdge);

      assert(rawEdge->idA == contigID || rawEdge->idB == contigID);

      if(otherCI->id < newContig->id){
	newEdge->orient = FlipEdgeOrient(newOrient);
	newEdge->idA = otherCI->id;
	newEdge->idB = newContig->id;
      }else{
	newEdge->orient = newOrient;
	newEdge->idB = otherCI->id;
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
#undef DEBUG_PROPAGATE
#ifdef DEBUG_PROPAGATE
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



#ifdef DEBUG_PROPAGATE
  DumpContig(GlobalData->stderrc,ScaffoldGraph, GetGraphNode(ScaffoldGraph->RezGraph, contig->id),FALSE);

  fprintf(GlobalData->stderrc,"* Calling extermal on " F_CID " *\n", aEndID);
#endif

  // Propagate overlaps from the contig at the a end of the new contig
  PropagateExtremalOverlapsToNewContig(aEndID, aEndEnd, contig, A_END, contigBase, verbose);

#ifdef DEBUG_PROPAGATE
  fprintf(GlobalData->stderrc,"* Calling extermal on " F_CID " *\n", bEndID);
  DumpContig(GlobalData->stderrc,ScaffoldGraph, GetGraphNode(ScaffoldGraph->RezGraph, contig->id),FALSE);
#endif
  // Propagate overlaps from the contig at the b end of the enw contig
  PropagateExtremalOverlapsToNewContig(bEndID, bEndEnd, contig, B_END, contigBase, verbose);

  // Propagate Internal, non-containment overlaps
  PropagateInternalOverlapsToNewContig(contig, ContigPositions, scaffold->id, contigBase, verbose);

#ifdef DEBUG_PROPAGATE
  DumpContig(GlobalData->stderrc,ScaffoldGraph, GetGraphNode(ScaffoldGraph->RezGraph, contig->id),FALSE);
  fprintf(GlobalData->stderrc,"* Calling containment *\n");
#endif
  // Propagate Containment Overlaps
  PropagateContainmentOverlapsToNewContig(contig, ContigPositions, contigBase, verbose);

#ifdef DEBUG_PROPAGATE
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

    if(GetSizeOfCurrentSequenceDB(ScaffoldGraph->sequenceDB) > GlobalData->maxSequencedbSize){
      fprintf(GlobalData->timefp, "\n\nCheckpoint %d written during CleanupScaffolds after scaffold " F_CID "\n",
              ScaffoldGraph->checkPointIteration, scaffold->id);
      fprintf(GlobalData->timefp, "Sorry, no way to know which logical checkpoint this really is!\n");
      CheckpointScaffoldGraph(ScaffoldGraph, -1);
    }

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

//  BPW -- appears to delete a contig from a scaffold if that contig
//  is made up entirely of surrogates and short.  See the printf
//  below.  Returns TRUE if it did this.
//
int  DeleteAllSurrogateContigsFromFailedMerges(CIScaffoldT *scaffold,
                                               NodeCGW_T *contig){
  int didSomething = FALSE;
  int32 i;
  MultiAlignT *ma;
  int numSurrogates;

#if 0
  fprintf(GlobalData->stderrc,"* DeleteAllSurrogateContigsFromFailedMerges scaffold " F_CID " contig " F_CID "\n",
	  scaffold->id, contig->id);
#endif

  if(contig->bpLength.mean > 2000)
    return FALSE;

  numSurrogates = 0;
  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);
  //ma = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id);

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
//
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
			   maxVariance, TRUE, TRUE, 0,TRUE);

    numComponents = IsScaffoldInternallyConnected(graph,scaffold,ALL_TRUSTED_EDGES);
    if(numComponents>1){
      //assert(numComponents == 1);
      fprintf(stderr,"WARNING  CUAS1: scaffold %d has %d components\n",scaffold->id,numComponents);
    }
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

	 // if we are not careful elsewhere, DP_Compare can return overlaps shorter than minlen;
	 // this makes the following tempting ... but since failed overlap gaps are set to -CGW_DP_MINLEN
	 // it seems to be a bad idea in the end ... [ALH, 09/04]
#ifdef ALLOW_SHORT_PERFECT
	 actual> (CGW_DP_MINLEN-(int)ceil( CGW_DP_ERATE *(double)CGW_DP_MINLEN )) ){ /*  || (lookForSmallOverlaps && SmallOverlapExists(ScaffoldGraph->RezGraph, prevCI, currCI, pairwiseOrient))){ */
#else
        actual> CGW_DP_MINLEN ){ /*  || (lookForSmallOverlaps && SmallOverlapExists(ScaffoldGraph->RezGraph, prevCI, currCI, pairwiseOrient))){ */
#endif

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

          //  XXX: DeleteAllSurrogate..() returns true if it deletes a
          //  contig, so if it does just one, we set
          //  allMergesSucceeded to true.  Possibly a bug?
          //
          //  Not sure if this ever gets used, though.
          //
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
			   maxVariance, TRUE, TRUE, 0,TRUE);

    numComponents = IsScaffoldInternallyConnected(graph,scaffold,ALL_TRUSTED_EDGES);
    if(numComponents>1){
      //assert(numComponents == 1);
      fprintf(stderr," CUAS: scaffold %d has %d components\n",scaffold->id,numComponents);
    }
  }
#endif


  if(mergesAttempted == 0)
    return 0;
  if(allMergesSucceeded)
    return TRUE;

  return NULLINDEX;
}


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

  if(ContigPositions == NULL)
    ContigPositions = CreateVA_IntElementPos(10);
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

//#include "obsolete/checkforcontigs"

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
    maxPos = MAX(maxPos, MAX(iup->position.bgn, iup->position.end));
  }
  fprintf(fp, "minPos = " F_COORD ", maxPos = " F_COORD "\n", minPos, maxPos);
  
  // Null out the source field
  for(i = 0; i < numFrag; i++){
    IntMultiPos *mp_i = GetIntMultiPos(ma->f_list,i);
    CIFragT *frag = GetCIFragT(graph->CIFrags,(CDS_CID_t)mp_i->sourceInt);
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
    maxPos = MAX(maxPos, MAX(mp_i->position.bgn, mp_i->position.end));
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
  int32 aEndEnd = NO_END, bEndEnd = NO_END;
  CDS_COORD_t minPos = CDS_COORD_MAX;
  CDS_COORD_t maxPos = CDS_COORD_MIN;
  int32 includesFinishedBacFragments = FALSE;

  contig->offsetAEnd = offsetAEnd;
  contig->offsetBEnd = offsetBEnd;

#ifdef DEBUG_CREATEACONTIG
  fprintf(GlobalData->stderrc,"* Create a contig in scaffold " F_CID "\n", scaffold->id);
#endif
      
  //  Determine which ends are extremal, so we can propagate their
  //  overlaps to the new contig.
  //
  for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions, i);
    ContigT       *ctg = GetGraphNode(ScaffoldGraph->ContigGraph, pos->ident);

    ctg->flags.bits.failedToContig = FALSE;
    ctg->flags.bits.beingContigged = TRUE; // flag the contig as being involved in a contigging operation

    includesFinishedBacFragments |= ctg->flags.bits.includesFinishedBacFragments;

    if (pos->position.bgn < pos->position.end) {
      if (pos->position.bgn < minPos) {
        minPos  = pos->position.bgn;
        aEndID  = pos->ident;
        aEndEnd = A_END;
      }
      if (pos->position.end > maxPos) {
        maxPos  = pos->position.end;
        bEndID  = pos->ident;
        bEndEnd = B_END;
      }
    } else {
      if (pos->position.end < minPos) {
        minPos  = pos->position.end;
        aEndID  = pos->ident;
        aEndEnd = B_END;
      }
      if (pos->position.bgn > maxPos) {
        maxPos  = pos->position.bgn;
        bEndID  = pos->ident;
        bEndEnd = A_END;
      }
    }
  }


#ifdef DEBUG_CREATEACONTIG
  fprintf(GlobalData->stderrc, "minPos = " F_COORD " maxPos = " F_COORD " extremes are (" F_CID ",%c) and (" F_CID ",%c)\n", 
          minPos,maxPos,
          aEndID,aEndEnd,
          bEndID,bEndEnd );
#endif

#if 0
  //  Normalize element positions to zero
  for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions,i);
        
    pos->position.bgn -= minPos;
    pos->position.end -= minPos;
  }
#endif
      
#ifdef DEBUG_CREATEACONTIG
  for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions,i);
    fprintf(GlobalData->stderrc,"* Contig " F_CID "  bgn:" F_COORD " end:" F_COORD "\n",
            pos->ident, pos->position.bgn, pos->position.end);
  }
#endif

  StartTimerT(&GlobalData->ConsensusTimer);

  newMultiAlign = MergeMultiAlignsFast_new(ScaffoldGraph->sequenceDB,
                                           ScaffoldGraph->fragStore,
                                           ContigPositions,
                                           FALSE,
                                           TRUE,
                                           GlobalData->aligner,
                                           NULL);

  StopTimerT(&GlobalData->ConsensusTimer);

  if(!newMultiAlign){
    int32 i;

    for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++) {
      IntElementPos *pos     = GetIntElementPos(ContigPositions,i);
      ContigT        *contig = GetGraphNode(ScaffoldGraph->ContigGraph, pos->ident);

      contig->flags.bits.failedToContig = TRUE;

#ifdef DEBUG_CREATEACONTIG
      fprintf(GlobalData->stderrc,"* Contig " F_CID " (" F_CID ")   bgn:" F_COORD " end:" F_COORD " [" F_COORD "," F_COORD "]\n",
              pos->ident, GetOriginalContigID(contig->id),
              pos->position.bgn, pos->position.end,
              contig->aEndCoord, contig->bEndCoord);
#endif
    }

    DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffold, FALSE);

    fprintf(GlobalData->stderrc,"* MergeMultiAligns failed....bye\n");
    fflush(NULL);

    if(GlobalData->failOn_NoOverlapFound)
      assert(0 /* No overlap found error in merge multialignments */);

    //	  SetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id, newMultiAlign); // DeleteGraphNode expects this

    // NOTE: we do not keep the new multi-align in the cache, newMultiAlign is NULL anyways
    InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE, newMultiAlign, FALSE);
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

    {
      // In the course of interleaved scaffold merging, we may create a new contig
      // from pieces of both contributing scaffolds.  In this event, the adjusted
      // variances of the original pieces may not be monotonic, which was leading
      // to non-monotonic variances on the contigs in the resulting scaffold.  To
      // "fix" this, we ensure that the resulting variance is based on the input piece
      // with the lowest variance.
      // Something similar happens when contigs are reordered (fix up a suspicious overlap)

      int i,errflag=0,errID;
      double minvariance = newOffsetAEnd.variance;
      for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
        IntElementPos *pos = GetIntElementPos(ContigPositions,i);
        ContigT *anOrigcontig = GetGraphNode(ScaffoldGraph->ContigGraph, pos->ident);
        if(anOrigcontig->offsetAEnd.variance<minvariance||
           anOrigcontig->offsetBEnd.variance<minvariance){
          errflag=1;
          errID = pos->ident;
          if(anOrigcontig->offsetAEnd.variance <
             anOrigcontig->offsetBEnd.variance){
            minvariance = anOrigcontig->offsetAEnd.variance;
          } else {
            minvariance = anOrigcontig->offsetBEnd.variance;
          }
        }
      }
      if(errflag){
        double vardiff = newOffsetAEnd.variance-minvariance;

        fprintf(stderr,
                "WARNING: NEG. VARIANCE: Looks like in creating contig %d, the variance of the contig may be screwed up!\n"
                "Set to %e based on extremal orig contig, but there is another, non-extremal orig contig with min variance %e, namely %d\n",
                contig->id,
                newOffsetAEnd.variance,
                minvariance,
                errID);
        fprintf(stderr,"This sometimes happens when ...\n"
                "\t(a) scaffolds are interleaved [variance computation here not really coherent :-(]\n"
                "\t(b) contigs are reordered (FOEXS suspicious overlap ...)\n");

#ifdef DEBUG_NEG_VARIANCE
        DumpACIScaffoldNew(stderr,ScaffoldGraph,scaffold,FALSE);
#endif
        newOffsetAEnd.variance -= vardiff;

        fprintf(stderr,
                "Resetting new contig AEnd variance to %e\n",
                newOffsetAEnd.variance);
      }
    }

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
      fprintf(GlobalData->stderrc,"* Fragment " F_CID " [" F_COORD "," F_COORD "]\n",(CDS_CID_t)fpos->sourceInt, fpos->position.bgn, fpos->position.end);
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
  fprintf(GlobalData->stderrc,"*  >>> NEW CONTIG " F_CID " (multiAlign size = " F_SIZE_T ") <<< *\n", contig->id,
          GetMemorySize(newMultiAlign));
  DumpContig(GlobalData->stderrc,ScaffoldGraph, contig,FALSE);
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


//#include "obsolete/contig_containment_old"
//#include "obsolete/not_contig_containment"

#define AHANGSLOP 30


//  This function is called by LeastSquares to merge together two
//  contigs that have a containment relationship between them.
//
//  This version tries several ways of finding the overlap
//  relationship if tryHarder is set.
//
//  BPW claims ownership.
//
void ContigContainment(CIScaffoldT  *scaffold,
                       NodeCGW_T    *prevCI,
                       NodeCGW_T    *thisCI,
                       EdgeCGW_T    *overlapEdge,
                       int           tryHarder) {
  CDS_COORD_t            minAhang;
  CDS_COORD_t            maxAhang;
  IntElementPos          contigPos;
  Overlap               *contigOverlap;
  ChunkOrientationType   overlapOrientation;
  ChunkOrientationType   actualOverlapOrientation;
  NodeCGW_T             *leftContig;
  NodeCGW_T             *rightContig;

  static VA_TYPE(IntElementPos) *ContigPositions = NULL;

  if (ContigPositions == NULL)
    ContigPositions = CreateVA_IntElementPos(10);

  ResetVA_IntElementPos(ContigPositions);

  if (min(prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean) <= 
      min(thisCI->offsetAEnd.mean, thisCI->offsetBEnd.mean)) {
    fprintf(stderr, "leftContig is prevCI\n");
    leftContig = prevCI;
    rightContig = thisCI;
  } else {
    fprintf(stderr, "leftContig is thisCI (swapped order)\n");
    leftContig = thisCI;
    rightContig = prevCI;
  }


  // Ad hominem critique (ALH): the following seems at this late hour an insane way to
  // set the allowed ahang range: the input edge should be used rather than the 
  // positions the contigs currently think they have in the scaffold ... which can be way off!

  if (leftContig->offsetAEnd.mean < leftContig->offsetBEnd.mean) {  // leftContig is AB
    if (rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean) {  // rightContig is AB
      overlapOrientation = AB_AB;
      minAhang = (CDS_COORD_t) (rightContig->offsetAEnd.mean - leftContig->offsetAEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    } else {  // rightContig is BA
      overlapOrientation = AB_BA;
      minAhang = (CDS_COORD_t) (rightContig->offsetBEnd.mean - leftContig->offsetAEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    }
  } else {  // leftContig is BA
    if (rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean) {  // rightContig is AB
      overlapOrientation = BA_AB;
      minAhang = (CDS_COORD_t) (rightContig->offsetAEnd.mean - leftContig->offsetBEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    } else {  // rightContig is BA
      overlapOrientation = BA_BA;
      minAhang = (CDS_COORD_t) (rightContig->offsetBEnd.mean - leftContig->offsetBEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    }
  }

  fprintf(stderr, "calling OverlapContigs with orient: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n",
          overlapOrientation, minAhang, maxAhang);

  contigOverlap = OverlapContigs(leftContig, rightContig, &overlapOrientation, minAhang, maxAhang, TRUE);

  // if no overlap found, try flipping orientation in case DPCompare
  // asymmetry is biting us
  //
  if ((contigOverlap == NULL) && (tryHarder)) {
    fprintf(stderr, "no overlap found between " F_CID " and " F_CID ", retrying with flipped orientation\n",
            leftContig->id, rightContig->id);
	  
    overlapOrientation = InvertEdgeOrient(overlapOrientation);
    contigOverlap = OverlapContigs(leftContig, rightContig, &overlapOrientation, minAhang, maxAhang, TRUE);
    overlapOrientation = InvertEdgeOrient(overlapOrientation);

    //  correct for flipped orientation
    if (contigOverlap != NULL) {
      CDS_COORD_t temp      = -contigOverlap->begpos;
      contigOverlap->begpos = -contigOverlap->endpos;
      contigOverlap->endpos = temp;
    }
  }
	
  // if still no overlap found, try maxing out hangs
  if ((contigOverlap == NULL) && (tryHarder)) {
    CDS_COORD_t maxLength = MAX(leftContig->bpLength.mean, rightContig->bpLength.mean);
	  
    fprintf(stderr, "no overlap found between " F_CID " and " F_CID ", retrying with max AHANGSLOP\n",
            leftContig->id, rightContig->id);
    fprintf(stderr, "overlapOrientation: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n", 
            (char) overlapOrientation, -maxLength, maxLength);

    contigOverlap = OverlapContigs(leftContig, rightContig, &overlapOrientation, -maxLength, maxLength, FALSE);	  
  }
  
  // if still no overlap found, try flipping orientation and maxing out hangs
  if ((contigOverlap == NULL) && (tryHarder)) {
    CDS_COORD_t maxLength = MAX(leftContig->bpLength.mean, rightContig->bpLength.mean);

    fprintf(stderr, "no overlap found between " F_CID " and " F_CID ", retrying with flipped orientation and max AHANGSLOP\n",
            leftContig->id, rightContig->id);
    fprintf(stderr, "overlapOrientation: %c, minAhang: " F_COORD ", maxAhang: " F_COORD "\n", 
            (char) overlapOrientation, -maxLength, maxLength);

    overlapOrientation = InvertEdgeOrient(overlapOrientation);
    contigOverlap = OverlapContigs(leftContig, rightContig, &overlapOrientation, -maxLength, maxLength, FALSE);
    overlapOrientation = InvertEdgeOrient(overlapOrientation);

    //  correct for flipped orientation
    if (contigOverlap != NULL) {
      CDS_COORD_t temp      = -contigOverlap->begpos;
      contigOverlap->begpos = -contigOverlap->endpos;
      contigOverlap->endpos = temp;
    }
  }

  // Now, a series of attempts involving swapping the two contigs ...
  //
  // to swap, need to:
  //  - change orientation to the other contig's perspective
  //  - change ahang range to other contig's perspective -- unless we just use -/+maxLength
  //  - call OverlapContigs with left and right swapped
  //
  if ((contigOverlap == NULL) && (tryHarder)) {
    CDS_COORD_t maxLength = MAX(leftContig->bpLength.mean, rightContig->bpLength.mean);

    fprintf(stderr,"trying to swap contigs ... " F_CID " and " F_CID " with overlapOrientation %c minAhang: " F_COORD ", maxAhang: " F_COORD "\n",
            rightContig->id,leftContig->id,
            overlapOrientation,-maxLength,maxLength);

    overlapOrientation = EdgeOrientSwap(overlapOrientation);
    contigOverlap = OverlapContigs(rightContig, leftContig, &overlapOrientation, -maxLength, maxLength, FALSE);
    overlapOrientation = EdgeOrientSwap(overlapOrientation);
	    
    // after swapping fragments, need to correct the hangs
    if (contigOverlap != NULL) {
      contigOverlap->begpos *= -1;
      contigOverlap->endpos *= -1;
    } 
  }

  if ((contigOverlap == NULL) && (tryHarder)) {
    CDS_COORD_t maxLength = MAX(leftContig->bpLength.mean, rightContig->bpLength.mean);

    fprintf(stderr,"trying swapped and inverted contigs ... " F_CID " and " F_CID " with overlapOrientation %c minAhang: " F_COORD ", maxAhang: " F_COORD "\n",
            rightContig->id,leftContig->id,
            overlapOrientation,-maxLength,maxLength);

    overlapOrientation = EdgeOrientSwap(overlapOrientation);
    overlapOrientation = InvertEdgeOrient(overlapOrientation);
    contigOverlap = OverlapContigs(rightContig, leftContig, &overlapOrientation, -maxLength, maxLength, FALSE);
    overlapOrientation = InvertEdgeOrient(overlapOrientation);
    overlapOrientation = EdgeOrientSwap(overlapOrientation);

    //  correct for flipped orientation
    if (contigOverlap != NULL) {
      CDS_COORD_t temp      = -contigOverlap->begpos;
      contigOverlap->begpos = -contigOverlap->endpos;
      contigOverlap->endpos = temp;
    }
    // after swapping fragments, need to correct the hangs
    if (contigOverlap != NULL) {
      contigOverlap->begpos *= -1;
      contigOverlap->endpos *= -1;
    } 
  }

  // if we STILL do not have an overlap AND we have been relatively
  // conservative finding overlaps...get out the big guns, redo the whole thing.
  //
  if ((contigOverlap == NULL) && (tryHarder) && (GlobalData->aligner == DP_Compare)) {
    GlobalData->aligner = Local_Overlap_AS_forCNS;
    ContigContainment(scaffold, prevCI, thisCI, overlapEdge, tryHarder);
    GlobalData->aligner = DP_Compare;

    // if we found an overlap, then all the necessary side effects
    // were taken care of in the inner call of the function, so just
    // return
    if (contigOverlap != NULL)
      return;
  }


  if (contigOverlap == NULL) {
    fprintf(stderr, "no overlap found between " F_CID " and " F_CID ", aborting...\n",
            leftContig->id, rightContig->id);
    dumpContigInfo(leftContig);
    dumpContigInfo(rightContig);
    assert(0);
  }	
  
  // contigs need to be reversed
  if (contigOverlap->begpos < 0) {
    NodeCGW_T *t = leftContig;
    leftContig   = rightContig;
    rightContig  = t;
    
    // adjust Overlap fields for later use in positioning
    contigOverlap->begpos = - contigOverlap->begpos;
    contigOverlap->endpos = - contigOverlap->endpos;	

    switch (overlapOrientation) {
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
    fprintf(stderr, "* Switched right-left, orientation went from %c to %c\n",
            overlapOrientation, actualOverlapOrientation);
  } else {
    actualOverlapOrientation = overlapOrientation;
  }

  fprintf(stderr, "* Containing contig is " F_CID " contained contig is " F_CID " ahg:%d bhg:%d orient:%c\n",
          leftContig->id, rightContig->id, contigOverlap->begpos, contigOverlap->endpos, actualOverlapOrientation);

  fprintf(stderr, "* Initial Positions:\n\t" F_CID " [%g,%g]\n\t" F_CID " [%g,%g]\n",
          leftContig->id,  leftContig->offsetAEnd.mean,  leftContig->offsetBEnd.mean,
          rightContig->id, rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);

  // assume we leave the leftContig where it is
  switch (actualOverlapOrientation) {
    case AB_AB:
      rightContig->offsetAEnd.mean = leftContig->offsetAEnd.mean  + contigOverlap->begpos;
      rightContig->offsetBEnd.mean = rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
      break;
    case AB_BA:
      rightContig->offsetBEnd.mean = leftContig->offsetAEnd.mean  + contigOverlap->begpos;
      rightContig->offsetAEnd.mean = rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
      break;
    case BA_AB:
      rightContig->offsetAEnd.mean = leftContig->offsetBEnd.mean  + contigOverlap->begpos;
      rightContig->offsetBEnd.mean = rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
      break;
    case BA_BA:
      rightContig->offsetBEnd.mean = leftContig->offsetBEnd.mean  + contigOverlap->begpos;
      rightContig->offsetAEnd.mean = rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
      break;
    default:
      break;
  }

  contigPos.ident        = leftContig->id;
  contigPos.type         = AS_CONTIG;
  contigPos.position.bgn = leftContig->offsetAEnd.mean;
  contigPos.position.end = leftContig->offsetBEnd.mean;
  AppendIntElementPos(ContigPositions, &contigPos);

  contigPos.ident         = rightContig->id;
  contigPos.type          = AS_CONTIG;
  contigPos.position.bgn = rightContig->offsetAEnd.mean;
  contigPos.position.end = rightContig->offsetBEnd.mean;
  AppendIntElementPos(ContigPositions, &contigPos);

  fprintf(stderr,"* Final Positions:\n\t" F_CID " [%g,%g]\n\t" F_CID " [%g,%g]\n",
          leftContig->id,  leftContig->offsetAEnd.mean,  leftContig->offsetBEnd.mean,
          rightContig->id, rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);

  {
    int mergeStatus = 0;
    int flip        = (leftContig->offsetBEnd.mean < leftContig->offsetAEnd.mean);
    mergeStatus = CreateAContigInScaffold(scaffold,
                                          ContigPositions,
                                          flip ? leftContig->offsetBEnd : leftContig->offsetAEnd,
                                          flip ? leftContig->offsetAEnd : leftContig->offsetBEnd);
    assert(mergeStatus == TRUE);
  }
}

//#include "obsolete/contig_containment_new"
//#include "obsolete/scaffold_cleanup_dead"

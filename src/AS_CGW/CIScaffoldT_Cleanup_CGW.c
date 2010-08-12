
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
static char *rcsid = "$Id: CIScaffoldT_Cleanup_CGW.c,v 1.67 2010-08-12 19:19:48 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
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
#include <time.h>


#undef DEBUG_DETAILED
#undef DEBUG_CONNECTEDNESS
#undef DEBUG_PROPAGATE
#undef DEBUG_CONTIG
#undef DEBUG_CREATEACONTIG
#undef DEBUG_CONTIGCONTAINMENT
#undef DEBUG_NEG_VARIANCE

typedef struct{
  ChunkInstanceT *firstCI;
  ChunkInstanceT *lastCI;
  ChunkInstanceT *maxCI;
  CDS_CID_t firstCID;
  CDS_CID_t lastCID;
  CDS_CID_t maxCID;
  LengthT minOffset;
  LengthT maxOffset;
  int32 count;
}ContigEndsT;

VA_DEF(ContigEndsT);
VA_TYPE(ContigEndsT) *ContigEnds = NULL;


// Propagate Containment Overlaps
void  PropagateContainmentOverlapsToNewContig(ContigT *newContig,
                                              VA_TYPE(IntElementPos) *ContigPositions,
                                              int32 contigBase,
                                              int32 verbose){
  int32 i;

  for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions,i);
    ContigT *contig = GetGraphNode(ScaffoldGraph->ContigGraph,pos->ident);
    GraphEdgeIterator edges;
    EdgeCGW_T *edge;
    int flip =( pos->position.bgn > pos->position.end); // do we need to flip edge orientation
    int32 posBgn = pos->position.bgn - contigBase;
    int32 posEnd = pos->position.end - contigBase;

    InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, contig->id, ALL_END, ALL_EDGES,
                          GRAPH_EDGE_RAW_ONLY,
                          &edges);
    if(verbose) {
      fprintf(stderr,"* Scanning edges incident on contig "F_CID "\n", pos->ident);
    }

    while((edge = NextGraphEdgeIterator(&edges)) != NULL){
      PairOrient edgeOrient = GetEdgeOrientationWRT(edge,contig->id);
      PairOrient orient;
      LengthT distance;
      ContigT *otherContig = GetGraphNode(ScaffoldGraph->ContigGraph, (edge->idA == contig->id? edge->idB:edge->idA));
      int32 overlap;


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
        fprintf(stderr,"* Contain:  pos "F_CID "  flip:%d edgeOrient:%c beg:"F_S32" end:"F_S32" newContigLength:%g\n",
                pos->ident, flip, edgeOrient.toLetter(), posBgn, posEnd, newContig->bpLength.mean);
        PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph,"propogating ", edge, edge->idA);
      }

      //  Compute the two edges for the containment
      overlap = (int32) -edge->distance.mean;

      assert(edgeOrient.isUnknown() == false);

      if(!flip){
        if (edgeOrient.isAB_AB()) {

          //  A-------------------------newcontig-----------------------------------------B
          //        A ======================contig============================B
          //  -------------------------A ===========other==========B
          //                           |---------------overlap ---------------|
          //                           |-----AB_AB overlap length-------------------------|
          //
          orient.setIsAB_AB();
          distance.mean = -(newContig->bpLength.mean - posEnd + overlap );
        }
        if (edgeOrient.isBA_BA()) {
          //  A-------------------------newcontig-----------------------------------------B
          //        A ======================contig============================B
          //  -------------------------A ===========other==========B
          //        |---------------overlap -----------------------|
          //                             |-----AB_AB overlap length-----------------------|
          //
          orient.setIsAB_AB();
          distance.mean = -(newContig->bpLength.mean - posBgn - overlap + otherContig->bpLength.mean);
        }
        if (edgeOrient.isAB_BA()) {
          //  A-------------------------newcontig-----------------------------------------B
          //        A ======================contig============================B
          //  -------------------------B ===========other==========A
          //                           |---------------overlap ---------------|
          //                           |-----AB_BA overlap length----------------------------|
          //
          orient.setIsAB_BA();
          distance.mean = -(newContig->bpLength.mean - posEnd +  overlap);
        }
        if (edgeOrient.isBA_AB()) {
          //  A-------------------------newcontig-----------------------------------------B
          //        A ======================contig============================B
          //  -------------------------B ===========other==========A
          //        |-------------------overlap -------------------|
          //                           |---AB_BA overlap length----------------------------|
          //
          orient.setIsAB_BA();
          distance.mean = -(newContig->bpLength.mean - posBgn - overlap + otherContig->bpLength.mean);
        }
      }else{
        if (edgeOrient.isBA_AB()) {
          //  A-------------------------newcontig-----------------------------------------B
          //        B ======================contig============================A
          //  -------------------------A ===========other==========B
          //                           |-------------------overlap -----------|
          //                           |-----AB_AB overlap length-------------------------|
          //
          orient.setIsAB_AB();
          distance.mean = -(newContig->bpLength.mean - posBgn + overlap);
        }
        if (edgeOrient.isAB_BA()) {
          //  A-------------------------newcontig-----------------------------------------B
          //        B ======================contig============================A
          //  -------------------------A ===========other==========B
          //        |-------------------overlap -------------------|
          //                           |-----AB_AB overlap length----------------------------|
          //
          orient.setIsAB_AB();
          distance.mean = -(newContig->bpLength.mean - posEnd - overlap + otherContig->bpLength.mean);
        }
        if (edgeOrient.isBA_BA()) {
          //  A-------------------------newcontig-----------------------------------------B
          //        B ======================contig============================A
          //  -------------------------B ===========other==========A
          //                           |-------------------overlap -----------|
          //                           |-----AB_BA overlap length-------------------------|
          //
          orient.setIsAB_BA();
          distance.mean = -(newContig->bpLength.mean - posBgn  +  overlap);
        }
        if (edgeOrient.isAB_AB()) {
          //  A-------------------------newcontig-----------------------------------------B
          //        B ======================contig============================A
          //  -------------------------B ===========other==========A
          //         |-------------------overlap ------------------|
          //                           |-----AB_BA overlap length-------------------------|
          //
          orient.setIsAB_BA();
          distance.mean = -(newContig->bpLength.mean - posEnd - overlap + otherContig->bpLength.mean);
        }
      }
      // Add the edge to the graph and insert them
      {
        EdgeCGW_T *newEdge;
        CDS_CID_t eid;
        int flip = (newContig->id > otherContig->id);
        newEdge = GetFreeGraphEdge(ScaffoldGraph->ContigGraph);
        eid = GetVAIndex_EdgeCGW_T(ScaffoldGraph->ContigGraph->edges, newEdge);

        SetGraphEdgeStatus(ScaffoldGraph->ContigGraph, newEdge, UNKNOWN_EDGE_STATUS);
        if(!flip){
          newEdge->idA = newContig->id;
          newEdge->idB = otherContig->id;
          newEdge->orient = orient;
          newEdge->flags.bits.aContainsB = TRUE;
          newEdge->flags.bits.hasContainmentOverlap = TRUE;
        }else{
          newEdge->idA = otherContig->id;
          newEdge->idB = newContig->id;
          newEdge->orient = orient;
          newEdge->orient.flip();
          newEdge->flags.bits.bContainsA = TRUE;
          newEdge->flags.bits.hasContainmentOverlap = TRUE;
        }
        newEdge->distance = distance;
        newEdge->edgesContributing = 1;

        {
          ChunkOverlapCheckT olap = {0};
          if(LookupOverlap(ScaffoldGraph->ContigGraph, newEdge->idA, newEdge->idB, newEdge->orient, &olap)){
            FreeGraphEdgeByEID(ScaffoldGraph->ContigGraph, eid);
            continue;
          }
        }

        {
          int32 delta =
            (int32) MIN((3 * sqrt(newEdge->distance.variance)), 50);

          assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));

          ChunkOverlapCheckT olap =
            OverlapChunks(ScaffoldGraph->ContigGraph,
                          newEdge->idA, newEdge->idB,
                          newEdge->orient, // handle suspicious
                          -newEdge->distance.mean - delta,
                          -newEdge->distance.mean + delta,
                          AS_CGW_ERROR_RATE, FALSE);

          if(olap.overlap == 0 || olap.suspicious){
            if(olap.suspicious){
              fprintf(stderr,"* CIS_C0 SUSPICIOUS Overlap found! Looked for ("F_CID ","F_CID ",%c) [%d,%d] found ("F_CID ","F_CID ",%c) "F_S32"\n",
                      newEdge->idA, newEdge->idB, newEdge->orient.toLetter(),
                      (int)(-newEdge->distance.mean - delta),
                      (int) (-newEdge->distance.mean + delta),
                      olap.spec.cidA, olap.spec.cidB, olap.spec.orientation.toLetter(),
                      olap.overlap);
            }
            if(verbose) {
              fprintf(stderr,"* NO Overlap found!\n");
              PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph,"*F* ",
                             newEdge, newEdge->idA);
            }
            FreeGraphEdgeByEID(ScaffoldGraph->ContigGraph,eid);
            continue;
          }
        }
        if(verbose) {
          PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph,"*C* ",
                         newEdge, newEdge->idA);
        }
        InsertGraphEdge(ScaffoldGraph->ContigGraph, eid, FALSE);

        // Create an appropriate hash table entry
        CreateChunkOverlapFromEdge(ScaffoldGraph->ContigGraph, newEdge);

      }
    }
  }
}


static int CompareIntElementPosByBgnPos(const void * c1, const void * c2)
{
  IntElementPos * i1 = (IntElementPos *) c1;
  IntElementPos * i2 = (IntElementPos *) c2;

  return(MIN(i1->position.bgn, i1->position.end) -
         MIN(i2->position.bgn, i2->position.end));
}

// Sort routine to order edges that are not yet inserted in graph
//
static int CompareEdges(const void *c1, const void *c2){
  CDS_CID_t eid1 = *(CDS_CID_t *)c1;
  CDS_CID_t eid2 = *(CDS_CID_t *)c2;
  EdgeCGW_T *edge1 = GetGraphEdge(ScaffoldGraph->ContigGraph, eid1);
  EdgeCGW_T *edge2 = GetGraphEdge(ScaffoldGraph->ContigGraph, eid2);
  int diff;
  diff = edge1->idA - edge2->idA;
  if(diff)
    return diff;
  diff = edge1->idB - edge2->idB;
  if(diff)
    return diff;

  diff = (int)edge1->orient.toLetter() - (int)edge2->orient.toLetter();
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
                                          int32 contigBase,
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
    fprintf(stderr,"* PropagateInternalOverlaps \n");
  }

  for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions,i);
    ContigT *contig = GetGraphNode(ScaffoldGraph->ContigGraph,pos->ident);
    GraphEdgeIterator edges;
    EdgeCGW_T *edge;
    int flip =( pos->position.bgn > pos->position.end); // do we need to flip edge orientation
    int32 posBgn = pos->position.bgn - contigBase;
    int32 posEnd = pos->position.end - contigBase;

    /* We exclude containment since we handle it elsewhere */
    InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, contig->id, ALL_END, ALL_EDGES,
			  GRAPH_EDGE_RAW_ONLY|GRAPH_EDGE_EXCLUDE_CONTAINMENT,
			  &edges);
    if(verbose) {
      fprintf(stderr,"* Scanning edges incident on contig "F_CID "\n", pos->ident);
    }

    while((edge = NextGraphEdgeIterator(&edges)) != NULL){
      PairOrient edgeOrient = GetEdgeOrientationWRT(edge,contig->id);
      PairOrient orient1, orient2;
      LengthT distance1, distance2;
      ContigT *otherContig = GetGraphNode(ScaffoldGraph->ContigGraph, (edge->idA == contig->id? edge->idB:edge->idA));
      int32 overlap;
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
          (stderr,
           "* pos "F_CID "  flip:%d edgeOrient:%c beg:"F_S32" end:"F_S32" newContigLength:%g\n",
           pos->ident, flip, edgeOrient.toLetter(), posBgn, posEnd,
           newContig->bpLength.mean);
	PrintGraphEdge(stderr,
                       ScaffoldGraph->ContigGraph,"propogating*I* ",
                       edge, edge->idA);
      }

      //  Compute the two edges for the containment
      // NOTE -- not all of these will be containments!!!!

      overlap = -edge->distance.mean;

      if(!flip){
        if (edgeOrient.isAB_AB()) {
          //  A-------------------------newcontig-----------------------------------------B
          //        A ===========contig============B
          //  -------------------------A ===========other==========B
          //                           |-----AB_AB overlap length------------------------|
          //  |----------------BA_BA overlap length-----------|
          //
          orient1.setIsAB_AB();
          maxDistance = (newContig->bpLength.mean - posEnd + overlap);
          distance1.mean = -maxDistance;
          orient2.setIsBA_BA();
          distance2.mean = -(posEnd - overlap + MIN(otherContig->bpLength.mean, maxDistance));
        }
        if (edgeOrient.isBA_BA()) {
          //  A-------------------------newcontig-----------------------------------------B
          //                                  A ===========contig============B
          //  --------------------A ===========other==========B
          //  |----------------BA_BA overlap length-----------|
          //                             |-----------------AB_AB overlap length-----------|
          //
          orient1.setIsBA_BA();
          maxDistance = posBgn + overlap;
          distance1.mean = -maxDistance;
          orient2.setIsAB_AB();
          distance2.mean = -(newContig->bpLength.mean - posBgn + MIN(otherContig->bpLength.mean,maxDistance) - overlap);
        }
        if (edgeOrient.isAB_BA()) {
          //  A-------------------------newcontig-----------------------------------------B
          //        A ===========contig============B
          //  -------------------------B ===========other==========A
          //  |---------------------BA_AB overlap lengtnh-----------|
          //                             |-----------------AB_BA overlap length-----------|
          //
          orient1.setIsAB_BA();
          maxDistance = (newContig->bpLength.mean - posEnd + overlap);
          distance1.mean = -maxDistance;
          orient2.setIsBA_AB();
          distance2.mean = -(posEnd + MIN(maxDistance,otherContig->bpLength.mean) - overlap);
        }
        if (edgeOrient.isBA_AB()) {
          //  A-------------------------newcontig-----------------------------------------B
          //                                    A ===========contig============B
          //  ----------------B ===========other==========A
          //  |------------BA_AB overlap length-----------|
          //                  |----------------AB_BA overlap length----------------------|
          //
          orient1.setIsBA_AB();
          maxDistance = (posBgn + overlap);
          distance1.mean = -maxDistance;
          orient2.setIsAB_BA();
          distance2.mean = -(newContig->bpLength.mean - posBgn + MIN(maxDistance,otherContig->bpLength.mean) - overlap);
        }
      }else{
        if (edgeOrient.isBA_AB()) {
          //  A-------------------------newcontig-----------------------------------------B
          //        B ===========contig============A
          //  -------------------------A ===========other==========B
          //  |-----BA_BA overlap length---------------------------|
          //                             |----------------AB_AB overlap length-----------|
          //
          orient1.setIsBA_BA();
          maxDistance = (posBgn + otherContig->bpLength.mean - overlap);
          distance1.mean = -maxDistance;
          orient2.setIsAB_AB();
          distance2.mean = -(newContig->bpLength.mean - posBgn + overlap);
        }
        if (edgeOrient.isAB_BA()) {
          //  A-------------------------newcontig-----------------------------------------B
          //                                  B ===========contig============A
          //  --------------------A ===========other==========B
          //  |----------------BA_BA overlap length-----------|
          //                             |-----------------AB_AB overlap length-----------|
          //
          orient1.setIsAB_AB();
          maxDistance = (posEnd + overlap);
          distance1.mean = -(newContig->bpLength.mean - posEnd + MIN(maxDistance,otherContig->bpLength.mean) - overlap);
          orient2.setIsBA_BA();
          distance2.mean = -maxDistance;
        }
        if (edgeOrient.isBA_BA()) {
          //  A-------------------------newcontig-----------------------------------------B
          //        B ===========contig============A
          //  -------------------------B ===========other==========A
          //  |-----BA_AB overlap length---------------------------|
          //                             |-----------------AB_BA overlap length-----------|
          //
          orient1.setIsAB_BA();
          maxDistance = (newContig->bpLength.mean - posBgn + overlap);
          distance1.mean = -maxDistance;
          orient2.setIsBA_AB();
          distance2.mean = -(posBgn + MIN(maxDistance,otherContig->bpLength.mean) - overlap);
        }
        if (edgeOrient.isAB_AB()) {
          //  A-------------------------newcontig-----------------------------------------B
          //                                    B ===========contig============A
          //  ----------------B ===========other==========A
          //  |------------BA_AB overlap length-----------|
          //                    |--------------------------------AB_BA overlap length-----|
          //
          maxDistance = (posEnd + overlap);
          orient1.setIsBA_AB();
          distance1.mean = -maxDistance;
          orient2.setIsAB_BA();
          distance2.mean = -(newContig->bpLength.mean - posEnd + MIN(maxDistance,otherContig->bpLength.mean) - overlap);
        }
      }
      // Add them to the graph, do NOT insert them.  Keep a list of such edges in collectedEdges
      {
	ChunkOverlapCheckT olap = {0};
	// If we don't find either of these edges already in the graph, proceed
	if(!LookupOverlap(ScaffoldGraph->ContigGraph, newContig->id, (edge->idA == contig->id? edge->idB: edge->idA), orient1, &olap) &&
	   !LookupOverlap(ScaffoldGraph->ContigGraph, newContig->id, (edge->idA == contig->id? edge->idB: edge->idA), orient2, &olap))
	  {
	    EdgeCGW_T *newEdge;
	    CDS_CID_t eid;
	    int32  isContainment = (MIN(-distance1.mean,-distance2.mean) > otherContig->bpLength.mean);

            CDS_CID_t edgeID =
              GetVAIndex_EdgeCGW_T(ScaffoldGraph->ContigGraph->edges, edge);
	    newEdge = GetFreeGraphEdge(ScaffoldGraph->ContigGraph);
	    eid = GetVAIndex_EdgeCGW_T(ScaffoldGraph->ContigGraph->edges, newEdge);
	    edge = GetGraphEdge(ScaffoldGraph->ContigGraph, edgeID);
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
	      PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph,
                             "    collecting*I* ", newEdge, newEdge->idA);
            }

	    if(isContainment){
	      newEdge = GetFreeGraphEdge(ScaffoldGraph->ContigGraph);
	      eid = GetVAIndex_EdgeCGW_T(ScaffoldGraph->ContigGraph->edges, newEdge);
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
		PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph,
                               "    collecting*I* ", newEdge, newEdge->idA);
              }
              assert(!newEdge->flags.bits.isDeleted);
	    }else{
	      if(verbose){
		fprintf(stderr,"* ***NOT A CONTAINMENT!!!!\n");
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
      fprintf(stderr,"* Collected the following edges *\n");
    }
    for(i = 0; i < GetNumCDS_CID_ts(CollectedEdges); i++){
      e = GetGraphEdge(ScaffoldGraph->ContigGraph, *GetCDS_CID_t(CollectedEdges,i));
      if(verbose) {
        PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph," ", e, e->idA);
      }
    }
  }
  // Then we try to merge the edges incident on newContig.  Any that successfully merge are CANDIDATES
  // for creating a raw overlap edge.  We screen with dpalign first.


  {
    int32 i;
    int32 numEdges = GetNumCDS_CID_ts(CollectedEdges);
    CDS_CID_t currID = NULLINDEX;
    PairOrient currOrient;
    ResetVA_CDS_CID_t(TentativeEdges);

    for(i = 0; i < numEdges; ){
      CDS_CID_t eid = *GetCDS_CID_t(CollectedEdges,i);
      EdgeCGW_T *e = GetGraphEdge(ScaffoldGraph->ContigGraph, eid);
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
          fprintf(stderr,"* Edge "F_CID " to Discarded edges\n",
                  *GetCDS_CID_t(TentativeEdges,0));
#endif
          AppendCDS_CID_t(DiscardedEdges, GetCDS_CID_t(TentativeEdges,0));
	}else{
	  if(GlobalData->debugLevel > 0)
	    fprintf(stderr,"* Merging "F_CID ","F_CID ",%c   %d edges\n",
                    e->idA, currID, currOrient.toLetter(),
                    (int) GetNumCDS_CID_ts(TentativeEdges));
	  addedEdges = MergeGraphEdges(ScaffoldGraph->ContigGraph, TentativeEdges);
	  if(addedEdges> 0){
	    if(GlobalData->debugLevel > 0)
	      fprintf(stderr,"* Merging ("F_CID ","F_CID ",%c): Added %d edges!!!\n",
                      e->idA, currID, currOrient.toLetter(), addedEdges);

	    ee = GetGraphEdge(ScaffoldGraph->ContigGraph, *GetCDS_CID_t(TentativeEdges,0));
	    ee->flags.bits.bContainsA = FALSE;
	    ee->flags.bits.aContainsB = TRUE;
	    ee->flags.bits.hasContainmentOverlap = TRUE;

	    for(i = 0; i < addedEdges; i++){
#ifdef DEBUG_DETAILED
	      //	      if(GlobalData->debugLevel > 0)
              PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph,
                             "* Successful merge: ", ee, ee->idA);
#endif
	      AppendCDS_CID_t(MergedEdges, GetCDS_CID_t(TentativeEdges,i));
	    }
	  }
#ifdef DEBUG_DETAILED
          fprintf(stderr,"* addedEdges %d total:%d\n",
                  addedEdges, (int) GetNumCDS_CID_ts(TentativeEdges));
#endif
	}

	if(addedEdges > 0)
          for(i = addedEdges; i < GetNumCDS_CID_ts(TentativeEdges) ; i++){
            // keep track of which edges to discard
	    ee = GetGraphEdge(ScaffoldGraph->ContigGraph, *GetCDS_CID_t(TentativeEdges,i));
	    ee->flags.bits.bContainsA = FALSE;
	    ee->flags.bits.aContainsB = TRUE;
	    ee->flags.bits.hasContainmentOverlap = TRUE;

	    PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph,"To Discarded Edges", ee, ee->idA);
            AppendCDS_CID_t(DiscardedEdges, GetCDS_CID_t(TentativeEdges,i));
          }
	currID = NULLINDEX;
	currOrient.setIsUnknown();
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
	fprintf(stderr,"* Discarded the following DISCARDED edges:\n");
      }
      for(i = 0; i < GetNumCDS_CID_ts(DiscardedEdges); i++){
	CDS_CID_t eid;
	EdgeCGW_T *e;

	eid = *GetCDS_CID_t(DiscardedEdges,i);
	e = GetGraphEdge(ScaffoldGraph->ContigGraph, eid);

	//fprintf(stderr,"* i = %d eid = "F_CID "\n", i,eid);

	if(e->flags.bits.isDeleted){
	  fprintf(stderr,"* i = %d Edge "F_CID " is ALREADY deleted...WIERD!\n", i,eid);
	}else{
	  if(verbose) {
	    PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph," ", e, e->idA);
            fprintf(stderr,"* i = %d Edge "F_CID "!\n", i,eid);
          }
	  FreeGraphEdgeByEID(ScaffoldGraph->ContigGraph,eid);
	}
      }


      if(verbose) {
	fprintf(stderr,"* Inserting the following MERGED edges:\n");
      }

      for(i = 0; i < GetNumCDS_CID_ts(MergedEdges); i++){
	CDS_CID_t eid= *GetCDS_CID_t(MergedEdges,i);
	CDS_CID_t neid;
	EdgeCGW_T *e = GetGraphEdge(ScaffoldGraph->ContigGraph, eid);
	ChunkOverlapCheckT olap = {0};
	int32 delta = MIN((3 * sqrt(e->distance.variance)), 50);
	CDS_CID_t idA, idB;
	PairOrient orient;

	if(verbose) {
	  PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph,"dp_align check ",
                         e, e->idA);
        }

	idA = e->idA;
	idB = e->idB;
	orient = e->orient;
	if(e->idA > e->idB){
	  idA = e->idB;
	  idB = e->idA;
	  orient = e->orient;
          orient.flip();
	}


	{
	  ChunkOverlapCheckT olap = {0};
	  if(LookupOverlap(ScaffoldGraph->ContigGraph, e->idA, e->idB, e->orient, &olap)){
	    if(verbose) {
	      fprintf(stderr,"* Overlap found already...skipping\n");
            }
	    FreeGraphEdgeByEID(ScaffoldGraph->ContigGraph,eid);
	    continue;
	  }
	}

        assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));

	olap = OverlapChunks(ScaffoldGraph->ContigGraph,
                             idA, idB, orient, // handles suspicious
			     -e->distance.mean - delta,
			     -e->distance.mean + delta,
			     AS_CGW_ERROR_RATE, FALSE);



	if(olap.overlap == 0 || olap.suspicious){
          if(olap.suspicious){
            fprintf(stderr,"* CIS_C1 SUSPICIOUS Overlap found! Looked for ("F_CID ","F_CID ",%c) [%d,%d] found ("F_CID ","F_CID ",%c) "F_S32" \n",
                    idA, idB, orient.toLetter(),
                    (int)(-e->distance.mean - delta),
                    (int)(-e->distance.mean + delta),
                    olap.spec.cidA, olap.spec.cidB, olap.spec.orientation.toLetter(), olap.overlap);
          }

	  if(verbose) {
	    fprintf(stderr,"* NO Overlap found!\n");
          }
	  FreeGraphEdgeByEID(ScaffoldGraph->ContigGraph,eid);
	  continue;
	}

	if(verbose) {
	  fprintf(stderr,"* Looked for [%g,%g] found "F_S32" overlap\n",
                  e->distance.mean - delta,
                  e->distance.mean + delta,
                  olap.overlap);
        }

	newEdge = GetFreeGraphEdge(ScaffoldGraph->ContigGraph);
	neid = GetVAIndex_EdgeCGW_T(ScaffoldGraph->ContigGraph->edges, newEdge);

        //  Refresh the e pointer; GetFreeGraphEdge can reallocate.
	e = GetGraphEdge(ScaffoldGraph->ContigGraph, eid);

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
          newEdge->orient.flip();
	}
	newEdge->edgesContributing = 1;

	if(verbose) {
	  PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph,
                         "Inserting #I# ", newEdge, newEdge->idA);
        }

	InsertGraphEdge(ScaffoldGraph->ContigGraph, neid, FALSE);

	// Create an appropriate hash table entry
	CreateChunkOverlapFromEdge(ScaffoldGraph->ContigGraph, newEdge);


	FreeGraphEdgeByEID(ScaffoldGraph->ContigGraph,eid);
      }
    }
  }
  RecycleDeletedGraphElements(ScaffoldGraph->ContigGraph);
}


// Propagate overlap edges incident on extremal contigs to new contig
void PropagateExtremalOverlapsToNewContig(CDS_CID_t contigID, int contigEnd, ContigT *newContig, int newContigEnd, int32 contigBase, int32 verbose){
  GraphEdgeIterator edges;
  EdgeCGW_T *edge;
  int first = 0;
  CDS_CID_t rawEdgeID;

  if(verbose){
    fprintf(stderr,"* PropagateExtremalOverlaps from ("F_CID ",%c) to ("F_CID ",%c)\n",
            contigID, (contigEnd == A_END?'A':'B'),
            newContig->id, (newContigEnd == A_END?'A':'B'));

    DumpContig(stderr,ScaffoldGraph,
               GetGraphNode(ScaffoldGraph->ContigGraph, contigID),FALSE);
  }


  InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, contigID, contigEnd, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
  while((edge = NextGraphEdgeIterator(&edges)) != NULL){
    EdgeCGW_T *rawEdge, *newEdge;
    PairOrient newOrient = edge->orient;
    NodeCGW_T *otherCI = GetGraphNode(ScaffoldGraph->ContigGraph, (edge->idA == contigID? edge->idB:edge->idA));

    if(!isOverlapEdge(edge))
      continue;

    // Handle containment edges elsewhere
    if(edge->flags.bits.aContainsB ||
       edge->flags.bits.bContainsA ||
       otherCI->flags.bits.beingContigged)
      continue;

    rawEdge = edge;
    if(!rawEdge->flags.bits.isRaw){
      rawEdge = GetGraphEdge(ScaffoldGraph->ContigGraph, rawEdge->nextRawEdge);
      AssertPtr(rawEdge);
      while(!isOverlapEdge(rawEdge) && rawEdge->nextRawEdge != NULLINDEX)
	rawEdge = GetGraphEdge(ScaffoldGraph->ContigGraph, rawEdge->nextRawEdge);

    }
    //    PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph, "Propogating: ", rawEdge,contigID);
    //    fflush(stderr);
    AssertPtr(rawEdge);

    /*****  rawEdge should be propagated to be incident on the new contig!!! ****/
    if(verbose){
      if(first++ == 0)
	fprintf(stderr,"* Need to propagate the following edges from ("F_CID ",%c) to the new contig's %c End\n",
                contigID,
                (contigEnd == A_END?'A':'B'),
                (newContigEnd == A_END?'A':'B')  );
      PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph, "X Prop *R* ", rawEdge, rawEdge->idA);
    }else{
      first++;
    }
    //UnlinkGraphEdge(ScaffoldGraph->ContigGraph, rawEdge); // Unlink it from the graph

    rawEdgeID = GetVAIndex_EdgeCGW_T(ScaffoldGraph->ContigGraph->edges, rawEdge);
    newEdge = GetFreeGraphEdge(ScaffoldGraph->ContigGraph);

    rawEdge = GetGraphEdge(ScaffoldGraph->ContigGraph, rawEdgeID);
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
      assert(newOrient.isUnknown() == false);

      if        (newOrient.isAB_BA()) {
        assert(contigEnd == B_END);
        newOrient.setIsBA_BA();
      } else if (newOrient.isAB_AB()) {
        assert(contigEnd == B_END);
        newOrient.setIsBA_AB();
      } else if (newOrient.isBA_AB()) {
        assert(contigEnd == A_END);
        newOrient.setIsAB_AB();
      } else if (newOrient.isBA_BA()) {
        assert(contigEnd == A_END);
        newOrient.setIsAB_BA();
      }
    }

    // Modify the ida and idb to reflect the new contig
    {
      CDS_CID_t newEdgeID = GetVAIndex_EdgeCGW_T(ScaffoldGraph->ContigGraph->edges, newEdge);

      assert(rawEdge->idA == contigID || rawEdge->idB == contigID);

      if(otherCI->id < newContig->id){
        newEdge->orient = newOrient;
	newEdge->orient.flip();
	newEdge->idA = otherCI->id;
	newEdge->idB = newContig->id;
      }else{
	newEdge->orient = newOrient;
	newEdge->idB = otherCI->id;
	newEdge->idA = newContig->id;
      }

      {
        ChunkOverlapCheckT olap = {0};
        if(LookupOverlap(ScaffoldGraph->ContigGraph, newEdge->idA, newEdge->idB, newEdge->orient, &olap)){
          fprintf(stderr,"\t* Overlap already exists ("F_CID ","F_CID ",%c)...skip\n",
                  newEdge->idA, newEdge->idB, newEdge->orient.toLetter());
          FreeGraphEdge(ScaffoldGraph->ContigGraph, newEdge);
          continue;
        }
      }

      {
        int32 delta = MIN((3 * sqrt(newEdge->distance.variance)), 50);

        assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));

        ChunkOverlapCheckT olap =
          OverlapChunks(ScaffoldGraph->ContigGraph,
                        newEdge->idA, newEdge->idB,
                        newEdge->orient, // handles suspicious
                        -newEdge->distance.mean - delta,
                        -newEdge->distance.mean + delta,
                        AS_CGW_ERROR_RATE, FALSE);

        if(olap.overlap == 0 || olap.suspicious){
          if(olap.suspicious){
            fprintf(stderr,"* CIS_S2 SUSPICIOUS Overlap found! Looked for ("F_CID ","F_CID ",%c) found ("F_CID ","F_CID ",%c)\n",
                    newEdge->idA, newEdge->idB, newEdge->orient.toLetter(),
                    olap.spec.cidA, olap.spec.cidB, olap.spec.orientation.toLetter());
          }
          if(verbose) {
            PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph,
                           "X No Overlap!  ", newEdge, newEdge->idA);
          }

          newEdge->orient.flip();  //  Handles suspicious

          olap = OverlapChunks(ScaffoldGraph->ContigGraph,
                               newEdge->idA, newEdge->idB,
                               newEdge->orient,
                               -newEdge->distance.mean - delta,
                               -newEdge->distance.mean + delta,
                               AS_CGW_ERROR_RATE, FALSE);

          newEdge->orient.flip();

          if(olap.overlap > 0 && !olap.suspicious)
            fprintf(stderr,"*** DUMMY, you got the orientation backwards!!!! ****\n");
          FreeGraphEdge(ScaffoldGraph->ContigGraph,newEdge);
          continue;
        }
      }


      if(verbose){
        fprintf(stderr,"*$#$ Inserting the following edge in the graph\n");
        PrintGraphEdge(stderr,ScaffoldGraph->ContigGraph, "  X ", newEdge, newEdge->idA);
      }

      // Link the edge back into the graph
      InsertGraphEdge(ScaffoldGraph->ContigGraph, newEdgeID, FALSE);

      // Create an appropriate hash table entry
      CreateChunkOverlapFromEdge(ScaffoldGraph->ContigGraph, newEdge);

    }
  }
  if(verbose){
    DumpContig(stderr,ScaffoldGraph, newContig,FALSE);
  }
  RecycleDeletedGraphElements(ScaffoldGraph->ContigGraph);

}

void PropagateOverlapsToNewContig(ContigT *contig,
                                  VA_TYPE(IntElementPos) *ContigPositions,
                                  CDS_CID_t aEndID, int aEndEnd,
                                  CDS_CID_t bEndID, int bEndEnd,
                                  CDS_CID_t scaffoldID, int32 verbose){
  ContigT *aContig, *bContig;
  CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, scaffoldID);
  int32 contigBase = INT32_MAX;
  int32 i;


  // Get a baseline for contig offsets
  for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions,i);
    int32 minPos = MIN(pos->position.bgn, pos->position.end);
    if(contigBase > minPos)
      contigBase = minPos;
  }

  AssertPtr(scaffold);

#ifdef DEBUG_PROPAGATE
  fprintf(stderr,"* PropagateOverlaps...initial\n");
  DumpContig(stderr,ScaffoldGraph, GetGraphNode(ScaffoldGraph->ContigGraph, contig->id),FALSE);
  fprintf(stderr,"* Calling internal *\n");
#endif

  // First, propagate the tandem overlap marks to the new contig
  //
  aContig = GetGraphNode(ScaffoldGraph->ContigGraph, aEndID);
  bContig = GetGraphNode(ScaffoldGraph->ContigGraph, bEndID);



#ifdef DEBUG_PROPAGATE
  DumpContig(stderr,ScaffoldGraph, GetGraphNode(ScaffoldGraph->ContigGraph, contig->id),FALSE);

  fprintf(stderr,"* Calling extermal on "F_CID " *\n", aEndID);
#endif

  // Propagate overlaps from the contig at the a end of the new contig
  PropagateExtremalOverlapsToNewContig(aEndID, aEndEnd, contig, A_END, contigBase, verbose);

#ifdef DEBUG_PROPAGATE
  fprintf(stderr,"* Calling extermal on "F_CID " *\n", bEndID);
  DumpContig(stderr,ScaffoldGraph, GetGraphNode(ScaffoldGraph->ContigGraph, contig->id),FALSE);
#endif
  // Propagate overlaps from the contig at the b end of the enw contig
  PropagateExtremalOverlapsToNewContig(bEndID, bEndEnd, contig, B_END, contigBase, verbose);

  // Propagate Internal, non-containment overlaps
  PropagateInternalOverlapsToNewContig(contig, ContigPositions, scaffold->id, contigBase, verbose);

#ifdef DEBUG_PROPAGATE
  DumpContig(stderr,ScaffoldGraph, GetGraphNode(ScaffoldGraph->ContigGraph, contig->id),FALSE);
  fprintf(stderr,"* Calling containment *\n");
#endif
  // Propagate Containment Overlaps
  PropagateContainmentOverlapsToNewContig(contig, ContigPositions, contigBase, verbose);

#ifdef DEBUG_PROPAGATE
  DumpContig(stderr,ScaffoldGraph, GetGraphNode(ScaffoldGraph->ContigGraph, contig->id),FALSE);
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
  fprintf(stderr,"* Inserting new contig "F_CID " at (%g,%g) "
          "with delta (%g,%g) in place of:\n",
	  newContig->id, offsetAEnd.mean, offsetBEnd.mean,
          deltaOffsetBEnd.mean, deltaOffsetBEnd.variance);
#endif
  for(i = 0; i < numElements; i++){
    IntElementPos *pos = GetIntElementPos(contigPositions,i);
    ContigT *contig = GetGraphNode(ScaffoldGraph->ContigGraph,  pos->ident);
#ifdef DEBUG_DETAILED
    fprintf(stderr,"* Removing contig "F_CID " "
            "from scaffold "F_CID " prev:"F_CID " next:"F_CID " (%g,%g) (%g,%g) ratio:%g\n",
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


  InsertCIInScaffold(ScaffoldGraph, newContig->id, scaffold->id, offsetAEnd, offsetBEnd, TRUE, FALSE);


  if(newContig->AEndNext == NULLINDEX){ // Make sure we atart at the end
    LengthT offset;

    if(orientAB){
      offset.mean = -newContig->offsetAEnd.mean;
      offset.variance = -newContig->offsetAEnd.variance;
    }else{
      offset.mean = -newContig->offsetBEnd.mean;
      offset.variance = -newContig->offsetBEnd.variance;
    }

    //fprintf(stderr,"* newContig "F_CID " at AEND offset [%g,%g] \n", newContig->id, newContig->offsetAEnd.mean, newContig->offsetBEnd.mean);
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

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
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
              //CIScaffold.aEndCoord = CIScaffold.bEndCoord = -1;
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
  fprintf(stderr,"* CleanupAScaffolds\n");
#endif

  InitGraphNodeIterator(&scaffolds, sgraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){

    if(scaffold->type == REAL_SCAFFOLD)
      didSomething |= CleanupAScaffold(sgraph,scaffold, lookForSmallOverlaps, maxContigsInMerge, deleteUnMergedSurrogates);

    if((scaffold->id % 10000) == 0){
      fprintf(stderr,"* CleanupScaffolds through scaffold "F_CID "\n", scaffold->id);
      //ScaffoldGraph->tigStore->flushCache();
    }
  }

  //ScaffoldGraph->tigStore->flushCache();

  RecycleDeletedGraphElements(sgraph->ContigGraph);
  return didSomething;
}

/****************************************************************************/
int CleanupFailedMergesInScaffolds(ScaffoldGraphT *sgraph){
  CIScaffoldT *scaffold;
  GraphNodeIterator scaffolds;
  int madeChanges = FALSE;
  fprintf(stderr,"* CleanupFailedMergesInScaffolds\n");

  //  yanked rocks/stones may not be in scaffolds.  without this step,
  //  some will inappropriately be degenerates
  ReScaffoldPseudoDegenerates();

  if (!GlobalData->performCleanupScaffolds)
    return 0;

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
	//fprintf(stderr,"* Merge failure(16) in scaffold "F_CID " didSomething = %d \n", scaffold->id, didSomething);

	didSomething = CleanupAScaffold(sgraph,scaffold, FALSE, 4 , FALSE);
	if(didSomething == 0)
	  break;
	if(didSomething > 0)
	  madeChanges = TRUE;
	//fprintf(stderr,"* Merge failure(4) in scaffold "F_CID " didSomething = %d \n", scaffold->id, didSomething);
      }

      didSomething = CleanupAScaffold(sgraph,scaffold, FALSE, 3 , FALSE);
      if(didSomething == 0)
	break;
      if(didSomething > 0)
        madeChanges = TRUE;
      //fprintf(stderr,"* Merge failure(3) in scaffold "F_CID " didSomething = %d iteration %d\n", scaffold->id, didSomething, iteration);

      didSomething = CleanupAScaffold(sgraph,scaffold, FALSE, 2 , FALSE);
      if(didSomething > 0)
	madeChanges = TRUE;
      //fprintf(stderr,"* After iteration %d  didSomething = %d\n", iteration, didSomething);

      iteration++;
    }
  }

  //  Flushing after every scaffold is very expensive on larger (metagenomic) assemblies.
  //  Flushing is of debatable value anyway.
  //ScaffoldGraph->tigStore->flushCache();

  RecycleDeletedGraphElements(sgraph->ContigGraph);
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
  fprintf(stderr,"* DeleteAllSurrogateContigsFromFailedMerges scaffold "F_CID " contig "F_CID "\n",
	  scaffold->id, contig->id);
#endif

  if(contig->bpLength.mean > 2000)
    return FALSE;

  numSurrogates = 0;
  ma = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);

  for(i = 0; i < GetNumIntUnitigPoss(ma->u_list); i++){
    IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);
    NodeCGW_T *node = GetGraphNode(ScaffoldGraph->CIGraph, pos->ident);
    if(node->flags.bits.isSurrogate)
      numSurrogates++;
  }

  fprintf(stderr,"* numSurrogates:%d  numCI :%d numElementPos:%d\n",
          numSurrogates, contig->info.Contig.numCI,
          (int) GetNumIntUnitigPoss(ma->u_list));

  if(numSurrogates == contig->info.Contig.numCI){
    NodeCGW_T *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, contig->scaffoldID);

    didSomething = TRUE;

    fprintf(stderr,"*** Deleting contig "F_CID " from scaffold "F_CID " since it is a short, all surrogate contigging failue\n",
            contig->id, contig->scaffoldID);
    DumpContig(stderr,ScaffoldGraph, contig, FALSE);

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
  fprintf(stderr,"* CleanupAScaffold "F_CID " (max:%d)\n", scaffold->id, maxContigsInMerge);
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

  contig.firstCI = contig.maxCI = currCI = CI = NextCIScaffoldTIterator(&CIs);
  contig.firstCID = contig.firstCI->id;
  contig.maxCID = contig.maxCI->id;

  if(currCI->offsetAEnd.mean < currCI->offsetBEnd.mean){
    contig.minOffset = currCI->offsetAEnd;
    contig.maxOffset = currCI->offsetBEnd;
  }else{
    contig.maxOffset = currCI->offsetAEnd;
    contig.minOffset = currCI->offsetBEnd;
  }
  contig.count = 0;
  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
    int32 actual;

    prevCI = currCI;
    currCI = CI;

    contig.count++;
    contig.lastCI = prevCI;
    contig.lastCID = contig.lastCI->id;

    // The first CI in the contig must have the minimum offset, by definition
    // subsequent CIs are sorted by LEFTmost coordinate, so we need to track
    // the MAX span (rightmost coordinate)

#ifdef DEBUG_CONTIG
    fprintf(stderr,"* contig.min = %g contig.max = %g\n",
            contig.minOffset.mean, contig.maxOffset.mean);
#endif
    actual = IntervalsOverlap(contig.minOffset.mean, contig.maxOffset.mean,
                              currCI->offsetAEnd.mean, currCI->offsetBEnd.mean, -15000);

    assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));

    if (GlobalData->doUnjiggleWhenMerging == 1 && 
        actual >= 0 && actual <= CGW_DP_MINLEN && 
        currCI->flags.bits.isJiggled == 1 && currCI->offsetDelta.mean > 0) {
      IntElementPos aPos;
      if (((uint32)(contig.maxOffset.mean - trunc(currCI->offsetAEnd.mean)) == actual) ||
          ((uint32)(contig.maxOffset.mean - trunc(currCI->offsetBEnd.mean)) == actual)) {
         currCI->offsetAEnd.mean -= currCI->offsetDelta.mean;
         currCI->offsetAEnd.variance -= currCI->offsetDelta.variance;
         currCI->offsetBEnd.mean -= currCI->offsetDelta.mean;
         currCI->offsetBEnd.variance -= currCI->offsetDelta.variance;

         int failed = TRUE;
         fprintf(stderr, "CleanupAScaffold() Processing total multialign of %g\n",(fabs(contig.maxCI->offsetBEnd.mean - contig.maxCI->offsetAEnd.mean) + fabs(currCI->offsetBEnd.mean - currCI->offsetAEnd.mean))); 

         actual = IntervalsOverlap(contig.minOffset.mean, contig.maxOffset.mean, currCI->offsetAEnd.mean, currCI->offsetBEnd.mean, -15000);
         if (actual >= 0) {
        	 int32 expected_ahang = (GetNodeOrient(currCI).isForward() ? currCI->offsetAEnd.mean : currCI->offsetBEnd.mean);
        	 expected_ahang -= (GetNodeOrient(contig.maxCI).isForward() ? contig.maxCI->offsetAEnd.mean : contig.maxCI->offsetBEnd.mean);
        	 int32 expected_bhang = (GetNodeOrient(currCI).isForward() ? currCI->offsetBEnd.mean : currCI->offsetAEnd.mean);
        	 expected_bhang -= (GetNodeOrient(contig.maxCI).isForward() ? contig.maxCI->offsetBEnd.mean : contig.maxCI->offsetAEnd.mean);

        	 float maxVar = (GetNodeOrient(contig.maxCI).isForward() ? contig.maxCI->offsetBEnd.variance : contig.maxCI->offsetAEnd.variance);
        	 float currVar = (GetNodeOrient(currCI).isForward() ? currCI->offsetAEnd.variance : currCI->offsetBEnd.variance);
        	 float inferredVar = (maxVar >= 0 ? maxVar : 0);
        	 inferredVar += (currVar >= 0 ? currVar : 0);

        	 int32 minOverlap = MAX(CGW_MISSED_OVERLAP, (actual - (3.0 * sqrt(inferredVar))));
        	 int32 maxOverlap = (actual + (3.0 * sqrt(inferredVar)));
        	 int32 minAhang = contig.maxCI->bpLength.mean - maxOverlap;
        	 int32 maxAhang = contig.maxCI->bpLength.mean - minOverlap;
        	 PairOrient orient = GetChunkPairOrientation(GetNodeOrient(contig.maxCI), GetNodeOrient(currCI));

                 ALNoverlap* ovl = OverlapContigs(contig.maxCI, currCI, &orient, minAhang, maxAhang, FALSE, TRUE, TRUE);
        	 if (ovl != NULL && ScoreOverlap(ovl, actual, expected_ahang, expected_bhang, AS_CGW_ERROR_RATE, NULL, NULL, NULL) != 0) {
        		 failed = FALSE;
        		 fprintf(stderr, "CleanupAScaffold() Undoing jiggling for %d worked.\n", currCI->id);
        	 }
         }

         if (failed) {
            fprintf(stderr, "CleanupAScaffold() Error: Undoing jiggle failed multi-align, backing out changes\n", currCI->id);
            currCI->offsetAEnd.mean += currCI->offsetDelta.mean;
            currCI->offsetBEnd.mean += currCI->offsetDelta.mean;
         }
         currCI->flags.bits.isJiggled = 0;
         actual = IntervalsOverlap(contig.minOffset.mean, contig.maxOffset.mean, currCI->offsetAEnd.mean, currCI->offsetBEnd.mean, -15000);
      } else {
         fprintf(stderr, "CleanupAScaffold() Error: don't know how to unjiggle positions for contig id %d, skipping\n", currCI->id);
      }
    }
    //  artificially break longs merges if != NULLINDEX
    //
    // if we are not careful elsewhere, DP_Compare can return overlaps shorter than minlen;
    // this makes the following tempting ... but since failed overlap gaps are set to -CGW_DP_MINLEN
    // it seems to be a bad idea in the end ... [ALH, 09/04]
    //       actual > (CGW_DP_MINLEN-(int)ceil( AS_CGW_ERROR_RATE *(double)CGW_DP_MINLEN
    //
    if (((maxContigsInMerge == NULLINDEX) || (contig.count < maxContigsInMerge)) && (actual> CGW_DP_MINLEN)) {

#ifdef DEBUG_CONTIG
      fprintf(stderr,
              "* CI "F_CID " and "F_CID " mean positions overlap by "F_S32" (%g,%g) (%g,%g)\n",
              prevCI->id,currCI->id, actual, contig.minOffset.mean,
              contig.maxOffset.mean, currCI->offsetAEnd.mean,
              currCI->offsetBEnd.mean);
#endif

      if(currCI->offsetAEnd.mean > contig.maxOffset.mean){
        contig.maxOffset = currCI->offsetAEnd;
        contig.maxCI = currCI;
        contig.maxCID = contig.maxCI->id;
      }
      if(currCI->offsetBEnd.mean > contig.maxOffset.mean){
        contig.maxOffset = currCI->offsetBEnd;
        contig.maxCI = currCI;
        contig.maxCID = contig.maxCI->id;
      }

    }else{
#ifdef DEBUG_CONTIG
      fprintf(stderr,"* CI "F_CID " and "F_CID " mean positions GAP by "F_S32" (%g,%g) (%g,%g) first "F_CID " last "F_CID " \n",
              prevCI->id,currCI->id, actual,
              contig.minOffset.mean,  contig.maxOffset.mean,
              currCI->offsetAEnd.mean, currCI->offsetBEnd.mean,
              contig.firstCI->id, contig.lastCI->id);
#endif
#ifdef DEBUG_DETAILED
      if((maxContigsInMerge != NULLINDEX) &&
         (contig.count == maxContigsInMerge)){
        fprintf(stderr,"*^^^^^ Merge terminated at %d elements\n", maxContigsInMerge);
      }
#endif
      AppendContigEndsT(ContigEnds, &contig);
      contig.count = 0;
      contig.firstCI = currCI;
      contig.firstCID = contig.firstCI->id;
      contig.maxCI = currCI;
      contig.maxCID = contig.maxCI->id;

      if(currCI->offsetAEnd.mean < currCI->offsetBEnd.mean){
        contig.minOffset = currCI->offsetAEnd;
        contig.maxOffset = currCI->offsetBEnd;
      }else{
        contig.maxOffset = currCI->offsetAEnd;
        contig.minOffset = currCI->offsetBEnd;
      }
#ifdef DEBUG_CONTIG
      fprintf(stderr,"* Reseting contig  firstCI: "F_CID " min:%g max:%g\n",
              contig.firstCI->id, contig.minOffset.mean, contig.maxOffset.mean);
#endif
    }
  }

  /* Add last contig to the list */
  ++contig.count;
  contig.lastCI = currCI;
  contig.lastCID = contig.lastCI->id;
  AppendContigEndsT(ContigEnds, &contig);


  // Now we have a collection of contigs to build
  // Work through them one at a time

    ContigEndsT *ctg = GetContigEndsT(ContigEnds,0);
    int32 numContigs = GetNumContigEndsTs(ContigEnds);

    VA_TYPE(IntElementPos) *ContigPositions = CreateVA_IntElementPos(32);

    /* Now the ContigEnds VA contains the starts/ends of all the contigs */
    for(i = 0; i < numContigs; i++){
      ContigT *contig;
      IntElementPos contigPos;
      ResetVA_IntElementPos(ContigPositions);

      ctg = GetContigEndsT(ContigEnds,i);

      for(contig = GetGraphNode(graph->ContigGraph, ctg->firstCID);
          contig != NULL;
          contig = GetGraphNode(graph->ContigGraph, contig->BEndNext))
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

  Delete_VA(ContigPositions);

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
  int32 minPos;
  int32 maxPos;
  LengthT myOffsetAEnd;
  LengthT myOffsetBEnd;
  int retVal = FALSE;

  VA_TYPE(IntElementPos) *ContigPositions = CreateVA_IntElementPos(10);

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
      if(IntervalsOverlap(contig->offsetAEnd.mean, contig->offsetBEnd.mean,
                          minPos, maxPos, CGW_DP_MINLEN))
        {
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
        }
      else {
        if(maxPos < MIN(contig->offsetAEnd.mean, contig->offsetBEnd.mean)){
          break;
        }
      }
    }

  if(GetNumIntElementPoss(ContigPositions) > 1) {
      qsort(GetIntElementPos(ContigPositions,0),
            GetNumIntElementPoss(ContigPositions),
            sizeof(IntElementPos),
            CompareIntElementPosByBgnPos);
      // set will be shifted to 0 offset in this function
      retVal = CreateAContigInScaffold(scaffold, ContigPositions, myOffsetAEnd, myOffsetBEnd);
    }

  Delete_VA(ContigPositions);

  return(retVal);
}



void DumpMultiAlignT(FILE * fp, ScaffoldGraphT * graph,
                     ContigT * ctg, MultiAlignT * ma)
{
  int32 i;
  int32 minPos;
  int32 maxPos;
  int32 numFrag = GetNumIntMultiPoss(ma->f_list);
  int32 numUnitig = GetNumIntUnitigPoss(ma->u_list);
  IntUnitigPos *up = GetIntUnitigPos(ma->u_list,0);

  CIScaffoldT *scaffold = GetGraphNode(graph->ScaffoldGraph,
                                       ctg->scaffoldID);

  fprintf(fp,"* Contig(%s) "F_CID " with %d unitigs and %d frags\n",
          scaffold && (scaffold->type == REAL_SCAFFOLD) ?
          "AS_PLACED": "AS_UNPLACED",
          ctg->id, numUnitig, numFrag);
  fprintf(fp, " %dbp gapped, %dbp ungapped\n",
          (int) GetMultiAlignLength(ma),
          (int) GetMultiAlignUngappedLength(ma));

  minPos = INT32_MAX;
  maxPos = INT32_MIN;
  for(i = 0; i < numUnitig; i++){
    IntUnitigPos *iup = up + i;
    NodeCGW_T *unitig = GetGraphNode(graph->CIGraph, iup->ident);
    fprintf(fp, "unitig "F_CID ": ", iup->ident);
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
    fprintf(fp, "pos:("F_S32","F_S32") delta:%d\n",
            iup->position.bgn, iup->position.end,
            iup->delta_length);
    minPos = MIN(minPos, MIN(iup->position.bgn, iup->position.end));
    maxPos = MAX(maxPos, MAX(iup->position.bgn, iup->position.end));
  }
  fprintf(fp, "minPos = "F_S32", maxPos = "F_S32"\n", minPos, maxPos);

  // Null out the source field
  for(i = 0; i < numFrag; i++){
    IntMultiPos *mp_i = GetIntMultiPos(ma->f_list,i);
    CIFragT *frag = GetCIFragT(graph->CIFrags, mp_i->ident);
    fprintf(fp, "frag "F_CID ":  utg:"F_CID ", pos:("F_S32","F_S32")\n",
            mp_i->ident, ma->maID, mp_i->position.bgn, mp_i->position.end);
    minPos = MIN(minPos, MIN(mp_i->position.bgn, mp_i->position.end));
    maxPos = MAX(maxPos, MAX(mp_i->position.bgn, mp_i->position.end));
  }
  fprintf(fp, "minPos = "F_S32", maxPos = "F_S32"\n", minPos, maxPos);
}

/***************************************************************************/

//  This path takes care of the on-the-fly contigging of the inserted contig with the existing
//  contigs in the scaffld

int  CreateAContigInScaffold(CIScaffoldT *scaffold,
                             VA_TYPE(IntElementPos) *ContigPositions,
                             LengthT offsetAEnd, LengthT offsetBEnd){
  CDS_CID_t   aEndID = NULLINDEX;
  CDS_CID_t   bEndID = NULLINDEX;
  int32       aEndEnd = NO_END;
  int32       bEndEnd = NO_END;
  int32       minPos = INT32_MAX;
  int32       maxPos = INT32_MIN;


  //  Determine which ends are extremal, so we can propagate their
  //  overlaps to the new contig.
  //
  for(int32 i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions, i);
    ContigT       *ctg = GetGraphNode(ScaffoldGraph->ContigGraph, pos->ident);

    ctg->flags.bits.failedToContig = FALSE;
    ctg->flags.bits.beingContigged = TRUE; // flag the contig as being involved in a contigging operation

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
  fprintf(stderr, "minPos = "F_S32" maxPos = "F_S32" extremes are ("F_CID ",%d) and ("F_CID ",%d)\n",
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
  for(int32 i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions,i);
    fprintf(stderr,"* Contig "F_CID " %c bgn:"F_S32" end:"F_S32"\n",
            pos->ident, pos->type, pos->position.bgn, pos->position.end);
  }

  VERBOSE_MULTIALIGN_OUTPUT = 1;
#endif

  MultiAlignT *newMultiAlign = MergeMultiAlignsFast_new(ContigPositions, NULL);

  if (newMultiAlign == NULL) {
    for(int32 i = 0; i < GetNumIntElementPoss(ContigPositions); i++) {
      IntElementPos *pos = GetIntElementPos(ContigPositions,i);
      ContigT       *ctg = GetGraphNode(ScaffoldGraph->ContigGraph, pos->ident);

      ctg->flags.bits.failedToContig = TRUE;

#ifdef DEBUG_CREATEACONTIG
      fprintf(stderr,"* Contig "F_CID" ("F_CID")   bgn:"F_S32" end:"F_S32"\n",
              pos->ident, GetOriginalContigID(ctg->id),
              pos->position.bgn, pos->position.end);
#endif
    }

#ifdef DEBUG_CREATEACONTIG
    DumpCIScaffold(stderr, ScaffoldGraph, scaffold, FALSE);
    fprintf(stderr,"* MergeMultiAligns failed....bye\n");
#endif

    return FALSE;
  }

  //  We have a multialign, create a new contig.

  ContigT *contig = CreateNewGraphNode(ScaffoldGraph->ContigGraph);

#ifdef DEBUG_CREATEACONTIG
  fprintf(stderr,"* Create a contig "F_CID " in scaffold "F_CID "\n", contig->id, scaffold->id);
#endif

  newMultiAlign->maID = contig->id;

  contig->offsetAEnd = offsetAEnd;
  contig->offsetBEnd = offsetBEnd;

  contig->bpLength.mean     = GetMultiAlignUngappedLength(newMultiAlign);
  contig->bpLength.variance = ComputeFudgeVariance(contig->bpLength.mean);

  contig->info.Contig.numCI =  GetNumIntUnitigPoss(newMultiAlign->u_list);

  // SK: for closure, when we merge together a unique closure unitig, propagate the flag up to the new contig
  for(int32 i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions, i);
    ContigT       *ctg = GetGraphNode(ScaffoldGraph->ContigGraph, pos->ident);

    contig->flags.bits.isClosure |= ctg->flags.bits.isClosure;
  }

#ifdef DEBUG_CREATEACONTIG
  for (int32 i=0; i<GetNumIntUnitigPoss(newMultiAlign->u_list); i++) {
    IntUnitigPos  *iup = GetIntUnitigPos(newMultiAlign->u_list, i);
    fprintf(stderr, "utg %d %d-%d\n", iup->ident, iup->position.bgn, iup->position.end);
  }
#endif

  NodeCGW_T *extremeA = GetGraphNode(ScaffoldGraph->ContigGraph, aEndID);
  NodeCGW_T *extremeB = GetGraphNode(ScaffoldGraph->ContigGraph, bEndID);

  LengthT  newOffsetAEnd = (aEndEnd == A_END) ? extremeA->offsetAEnd : extremeA->offsetBEnd;
  LengthT  oldOffsetBEnd = (bEndEnd == A_END) ? extremeB->offsetAEnd : extremeB->offsetBEnd;

  // In the course of interleaved scaffold merging, we may create a new contig from pieces of both
  // contributing scaffolds.  In this event, the adjusted variances of the original pieces may not
  // be monotonic, which was leading to non-monotonic variances on the contigs in the resulting
  // scaffold.  To "fix" this, we ensure that the resulting variance is based on the input piece
  // with the lowest variance.
  //
  // Something similar happens when contigs are reordered (fix up a suspicious overlap)
  //
  // This fails (errID > 0) sometimes when:
  //     (a) scaffolds are interleaved [variance computation here not really coherent :-(
  //     (b) contigs are reordered (FOEXS suspicious overlap
  //
  int32  errID       = -1;
  double minvariance = newOffsetAEnd.variance;

  for(int32 i = 0; i < GetNumIntElementPoss(ContigPositions); i++){
    IntElementPos *pos = GetIntElementPos(ContigPositions,i);
    ContigT       *ctg = GetGraphNode(ScaffoldGraph->ContigGraph, pos->ident);

    if ((ctg->offsetAEnd.variance < minvariance) ||
        (ctg->offsetBEnd.variance < minvariance)) {
      errID       = pos->ident;
      minvariance = (ctg->offsetAEnd.variance < ctg->offsetBEnd.variance) ? ctg->offsetAEnd.variance : ctg->offsetBEnd.variance;
    }
  }

  if (errID >= 0) {
    fprintf(stderr, "WARNING:  NEGATIVE VARIANCE in creating contig %d.  Set to %e based on extremal orig contig, but contig %d (non-extremal) has min variance %e\n",
            contig->id, newOffsetAEnd.variance, errID, minvariance);

#ifdef DEBUG_NEG_VARIANCE
    DumpACIScaffoldNew(stderr,ScaffoldGraph,scaffold,FALSE);
#endif
    newOffsetAEnd.variance -= (newOffsetAEnd.variance - minvariance);
  }

  // Determine the NEW mean/variance of the B end of the contig
  LengthT newOffsetBEnd;
  LengthT deltaOffsetBEnd;

  newOffsetBEnd.mean     = newOffsetAEnd.mean     + contig->bpLength.mean;
  newOffsetBEnd.variance = newOffsetAEnd.variance + contig->bpLength.variance;

  deltaOffsetBEnd.mean     = newOffsetBEnd.mean     - oldOffsetBEnd.mean;
  deltaOffsetBEnd.variance = newOffsetBEnd.variance - oldOffsetBEnd.variance;


  if(GlobalData->debugLevel > 0)
    fprintf(stderr,"* Merged Contig "F_CID " has ungapped length %d\n",
            contig->id, (int)contig->bpLength.mean);

  // Sort the Unitigs from left to right
  MakeCanonicalMultiAlignT(newMultiAlign);

  // INsert the multialign
  ScaffoldGraph->tigStore->insertMultiAlign(newMultiAlign, FALSE, TRUE);

  // Propagate Overlaps, Tandem Marks and Branch Points to the new Contig
  PropagateOverlapsToNewContig(contig, ContigPositions, aEndID, aEndEnd, bEndID, bEndEnd, scaffold->id, FALSE);

  // Insert the new contig into the scaffold, in lieu of the old contigs
  ReplaceContigsInScaffolds(scaffold, contig,  ContigPositions, newOffsetAEnd, newOffsetBEnd, deltaOffsetBEnd);

  // Mark all frags as being members of this Contig, and set their offsets
  UpdateNodeFragments(ScaffoldGraph->ContigGraph,contig->id, TRUE, FALSE);

  // Mark all of the Unitigs of this CI and set their offsets
  UpdateNodeUnitigs(newMultiAlign,contig);

  // Create the raw link-based edges
  {
    GraphEdgeStatT stats;
    BuildGraphEdgesFromMultiAlign(ScaffoldGraph->ContigGraph, contig, newMultiAlign, &stats, TRUE);
  }

  if(GlobalData->debugLevel > 0){
    fprintf(stderr,"* Here is the new contig before merging: [%g,%g]\n",
            newOffsetAEnd.mean, newOffsetBEnd.mean);
    DumpContig(stderr,ScaffoldGraph, contig,FALSE);
  }

  // Merge the edges incident on this contig
  MergeNodeGraphEdges(ScaffoldGraph->ContigGraph, contig, FALSE, TRUE, FALSE);

#ifdef DEBUG_DETAILED
  fprintf(stderr,"*  >>> NEW CONTIG "F_CID" <<< *\n", contig->id);
  DumpContig(stderr,ScaffoldGraph, contig,FALSE);
#endif

  return TRUE;
}


#define AHANGSLOP 30


//  This function is called by LeastSquares to merge together two
//  contigs that have a containment relationship between them.
//
//  This version tries several ways of finding the overlap
//  relationship if tryHarder is set.
//
//  Returns true if it succeeded.  Previous behavior was to return
//  success or crash.
//
//  BPW claims ownership.
//
int
ContigContainment(CIScaffoldT  *scaffold,
                  NodeCGW_T    *prevCI,
                  NodeCGW_T    *thisCI,
                  EdgeCGW_T    *overlapEdge,
                  int           tryHarder) {
  int32            minAhang;
  int32            maxAhang;
  IntElementPos          contigPos;
  ALNoverlap            *contigOverlap;
  PairOrient             overlapOrientation;
  PairOrient             actualOverlapOrientation;
  NodeCGW_T             *leftContig;
  NodeCGW_T             *rightContig;


  if (MIN(prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean) <=
      MIN(thisCI->offsetAEnd.mean, thisCI->offsetBEnd.mean)) {
    leftContig = prevCI;
    rightContig = thisCI;
  } else {
    leftContig = thisCI;
    rightContig = prevCI;
  }


  // Ad hominem critique (ALH): the following seems at this late hour an insane way to
  // set the allowed ahang range: the input edge should be used rather than the
  // positions the contigs currently think they have in the scaffold ... which can be way off!

  if (leftContig->offsetAEnd.mean < leftContig->offsetBEnd.mean) {  // leftContig is AB
    if (rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean) {  // rightContig is AB
      overlapOrientation.setIsAB_AB();
      minAhang = (int32) (rightContig->offsetAEnd.mean - leftContig->offsetAEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    } else {  // rightContig is BA
      overlapOrientation.setIsAB_BA();
      minAhang = (int32) (rightContig->offsetBEnd.mean - leftContig->offsetAEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    }
  } else {  // leftContig is BA
    if (rightContig->offsetAEnd.mean < rightContig->offsetBEnd.mean) {  // rightContig is AB
      overlapOrientation.setIsBA_AB();
      minAhang = (int32) (rightContig->offsetAEnd.mean - leftContig->offsetBEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    } else {  // rightContig is BA
      overlapOrientation.setIsBA_BA();
      minAhang = (int32) (rightContig->offsetBEnd.mean - leftContig->offsetBEnd.mean) - AHANGSLOP;
      maxAhang = minAhang + (2 * AHANGSLOP);
    }
  }

#ifdef DEBUG_CONTIGCONTAINMENT
  fprintf(stderr, "ContigContainment()-- orient=%d ahang min,max=%d,%d\n",
          overlapOrientation, minAhang, maxAhang);
#endif

  contigOverlap = OverlapContigs(leftContig, rightContig, &overlapOrientation, minAhang, maxAhang, TRUE);

  // if no overlap found, try flipping orientation in case DPCompare
  // asymmetry is biting us
  //
  if ((contigOverlap == NULL) && (tryHarder)) {
#ifdef DEBUG_CONTIGCONTAINMENT
    fprintf(stderr, "ContigContainment()--  Try #2\n");
#endif
    overlapOrientation.invert();
    contigOverlap = OverlapContigs(leftContig, rightContig, &overlapOrientation, minAhang, maxAhang, TRUE);
    overlapOrientation.invert();

    if (contigOverlap != NULL) {
      int32 temp      = -contigOverlap->begpos;
      contigOverlap->begpos = -contigOverlap->endpos;
      contigOverlap->endpos = temp;
    }
  }

  // if still no overlap found, try maxing out hangs
  if ((contigOverlap == NULL) && (tryHarder)) {
#ifdef DEBUG_CONTIGCONTAINMENT
    fprintf(stderr, "ContigContainment()--  Try #3\n");
#endif
    int32 maxLength = MAX(leftContig->bpLength.mean, rightContig->bpLength.mean);

    contigOverlap = OverlapContigs(leftContig, rightContig, &overlapOrientation, -maxLength, maxLength, FALSE);
  }

  // if still no overlap found, try flipping orientation and maxing out hangs
  if ((contigOverlap == NULL) && (tryHarder)) {
#ifdef DEBUG_CONTIGCONTAINMENT
    fprintf(stderr, "ContigContainment()--  Try #4\n");
#endif
    int32 maxLength = MAX(leftContig->bpLength.mean, rightContig->bpLength.mean);

    overlapOrientation.invert();
    contigOverlap = OverlapContigs(leftContig, rightContig, &overlapOrientation, -maxLength, maxLength, FALSE);
    overlapOrientation.invert();

    if (contigOverlap != NULL) {
      int32 temp      = -contigOverlap->begpos;
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

  //  swap
  if ((contigOverlap == NULL) && (tryHarder)) {
#ifdef DEBUG_CONTIGCONTAINMENT
    fprintf(stderr, "ContigContainment()--  Try #5\n");
#endif
    int32 maxLength = MAX(leftContig->bpLength.mean, rightContig->bpLength.mean);

    overlapOrientation.swap();
    contigOverlap = OverlapContigs(rightContig, leftContig, &overlapOrientation, -maxLength, maxLength, FALSE);
    overlapOrientation.swap();

    if (contigOverlap != NULL) {
      contigOverlap->begpos *= -1;
      contigOverlap->endpos *= -1;
    }
  }

  //  swap and flip orientation
  if ((contigOverlap == NULL) && (tryHarder)) {
#ifdef DEBUG_CONTIGCONTAINMENT
    fprintf(stderr, "ContigContainment()--  Try #6\n");
#endif
    int32 maxLength = MAX(leftContig->bpLength.mean, rightContig->bpLength.mean);

    overlapOrientation.swap();
    overlapOrientation.invert();

    contigOverlap = OverlapContigs(rightContig, leftContig, &overlapOrientation, -maxLength, maxLength, FALSE);

    overlapOrientation.invert();
    overlapOrientation.swap();

    if (contigOverlap != NULL) {
      int32 temp      = -contigOverlap->begpos;
      contigOverlap->begpos = -contigOverlap->endpos;
      contigOverlap->endpos = temp;

      contigOverlap->begpos *= -1;
      contigOverlap->endpos *= -1;
    }
  }

  //  High error rate assemblies have trouble in eCR.  eCR finds
  //  overlaps at CGW_ERATE + 0.02, builds a new contig, and then dies
  //  when it calls RecomputeOffsetsInScaffold() (which calls us).
  //
  if ((contigOverlap == NULL) && (tryHarder)) {
#ifdef DEBUG_CONTIGCONTAINMENT
    fprintf(stderr, "ContigContainment()--  Try #7\n");
#endif
    int32 maxLength = MAX(leftContig->bpLength.mean, rightContig->bpLength.mean);

    AS_CGW_ERROR_RATE += 0.02;

    contigOverlap = OverlapContigs(leftContig, rightContig, &overlapOrientation, minAhang, maxAhang, TRUE);

    if (contigOverlap == NULL)
      contigOverlap = OverlapContigs(leftContig, rightContig, &overlapOrientation, -maxLength, maxLength, FALSE);

    AS_CGW_ERROR_RATE -= 0.02;
  }

  if (contigOverlap == NULL) {
    fprintf(stderr, "================================================================================\n");
    dumpContigInfo(leftContig);
    dumpContigInfo(rightContig);
    fprintf(stderr, "* No overlap found between "F_CID" and "F_CID".  Fail.\n", leftContig->id, rightContig->id);
    return(FALSE);
  }
  //assert(contigOverlap != NULL);

  // contigs need to be reversed
  if (contigOverlap->begpos < 0) {
    NodeCGW_T *t = leftContig;
    leftContig   = rightContig;
    rightContig  = t;

    // adjust Overlap fields for later use in positioning
    contigOverlap->begpos = - contigOverlap->begpos;
    contigOverlap->endpos = - contigOverlap->endpos;

    //  AB_AB and BA_BA don't change
    actualOverlapOrientation = overlapOrientation;

    if (overlapOrientation.isAB_BA())
      actualOverlapOrientation.setIsBA_AB();
    if (overlapOrientation.isBA_AB())
      actualOverlapOrientation.setIsAB_BA();

    fprintf(stderr, "* Switched right-left, orientation went from %c to %c\n",
            overlapOrientation.toLetter(), actualOverlapOrientation.toLetter());
  } else {
    actualOverlapOrientation = overlapOrientation;
  }

#ifdef DEBUG_CONTIGCONTAINMENT
  fprintf(stderr, "* Containing contig is "F_CID" contained contig is "F_CID" ahg:%d bhg:%d orient:%c\n",
          leftContig->id, rightContig->id, contigOverlap->begpos, contigOverlap->endpos, actualOverlapOrientation);

  fprintf(stderr, "* Initial Positions:  left:"F_CID" [%g,%g]  right:"F_CID" [%g,%g]\n",
          leftContig->id,  leftContig->offsetAEnd.mean,  leftContig->offsetBEnd.mean,
          rightContig->id, rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);
#endif

  // assume we leave the leftContig where it is
  if (actualOverlapOrientation.isAB_AB()) {
    rightContig->offsetAEnd.mean = leftContig->offsetAEnd.mean  + contigOverlap->begpos;
    rightContig->offsetBEnd.mean = rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
  }
  if (actualOverlapOrientation.isAB_BA()) {
    rightContig->offsetBEnd.mean = leftContig->offsetAEnd.mean  + contigOverlap->begpos;
    rightContig->offsetAEnd.mean = rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
  }
  if (actualOverlapOrientation.isBA_AB()) {
    rightContig->offsetAEnd.mean = leftContig->offsetBEnd.mean  + contigOverlap->begpos;
    rightContig->offsetBEnd.mean = rightContig->offsetAEnd.mean + rightContig->bpLength.mean;
  }
  if (actualOverlapOrientation.isBA_BA()) {
    rightContig->offsetBEnd.mean = leftContig->offsetBEnd.mean  + contigOverlap->begpos;
    rightContig->offsetAEnd.mean = rightContig->offsetBEnd.mean + rightContig->bpLength.mean;
  }

  VA_TYPE(IntElementPos) *ContigPositions = CreateVA_IntElementPos(2);

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

#ifdef DEBUG_CONTIGCONTAINMENT
  fprintf(stderr, "* Final   Positions:  left:"F_CID" [%g,%g]  right:"F_CID" [%g,%g]\n",
          leftContig->id,  leftContig->offsetAEnd.mean,  leftContig->offsetBEnd.mean,
          rightContig->id, rightContig->offsetAEnd.mean, rightContig->offsetBEnd.mean);
#endif

  int flip        = (leftContig->offsetBEnd.mean < leftContig->offsetAEnd.mean);
  int mergeStatus = CreateAContigInScaffold(scaffold,
                                            ContigPositions,
                                            flip ? leftContig->offsetBEnd : leftContig->offsetAEnd,
                                            flip ? leftContig->offsetAEnd : leftContig->offsetBEnd);

  Delete_VA(ContigPositions);

  if (mergeStatus == FALSE) {
    fprintf(stderr, "* CreateAContigInScaffold() failed.\n");
    return(FALSE);
  }

  return(TRUE);
}





// Repair the CI neighbors in a contig AEndNext and BEndNext when removing a surrogate copy
void
RepairContigNeighbors(ChunkInstanceT *surr){
  ChunkInstanceT *AEndNeighbor = (surr->AEndNext >= 0) ? GetGraphNode(ScaffoldGraph->CIGraph, surr->AEndNext) : NULL;
  ChunkInstanceT *BEndNeighbor = (surr->BEndNext >= 0) ? GetGraphNode(ScaffoldGraph->CIGraph, surr->BEndNext) : NULL;

  ContigT *contig = GetGraphNode(ScaffoldGraph->ContigGraph, surr->info.CI.contigID);

  assert(contig != NULL);

  assert(surr->type == RESOLVEDREPEATCHUNK_CGW);

  if(AEndNeighbor != NULL){
    assert(AEndNeighbor->flags.bits.isCI);
    if(AEndNeighbor->AEndNext == surr->id){
      AEndNeighbor->AEndNext = surr->BEndNext;
    }else{
      assert(AEndNeighbor->BEndNext == surr->id);
      AEndNeighbor->BEndNext = surr->BEndNext;
    }
  }else{
    assert(contig->info.Contig.AEndCI == surr->id);
    contig->info.Contig.AEndCI = surr->BEndNext;
  }

  if(BEndNeighbor != NULL){
    assert(BEndNeighbor->flags.bits.isCI);
    if(BEndNeighbor->AEndNext == surr->id){
      BEndNeighbor->AEndNext = surr->AEndNext;
    }else{
      assert(BEndNeighbor->BEndNext == surr->id);
      BEndNeighbor->BEndNext = surr->AEndNext;
    }
  }else{
    assert(contig->info.Contig.BEndCI == surr->id);
    contig->info.Contig.BEndCI = surr->AEndNext;
  }

  // For the contig containing the removed duplicate surrogate, update
  // the SeqDB entry to delete the unitig position

  MultiAlignT *ma;
  int32 i, j;
  int32 delete_index = -1;

  //     1. get the old multialign from the seqDB
  ma =  ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);

  //     2. delete surrogate from multialign
  for(i = j = 0; i < GetNumIntUnitigPoss(ma->u_list); i++){
    IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);

    if(surr->id == pos->ident){
      delete_index = i;
    }else{
      if (i != j)
        SetIntUnitigPos(ma->u_list, j, GetIntUnitigPos(ma->u_list,i));
      j++;
    }
  }

  assert(delete_index >= 0);
  assert((j + 1) == i);
  ResetToRange_IntUnitigPos(ma->u_list,j);

  //     3. update the multialign
  assert(ma->maID == contig->id);
  ScaffoldGraph->tigStore->insertMultiAlign(ma, FALSE, TRUE);
}



// Comparison for qsort in RemoveSurrogateDuplicates

static
int
CompareSurrogatePlacements(const void *c1, const void *c2) {
  ChunkInstanceT *s1 = GetGraphNode(ScaffoldGraph->CIGraph, *(CDS_CID_t *)c1);
  ChunkInstanceT *s2 = GetGraphNode(ScaffoldGraph->CIGraph, *(CDS_CID_t *)c2);

  assert(s1 != NULL);
  assert(s2 != NULL);

  if (s1->info.CI.contigID > s2->info.CI.contigID)
    return((int)1);
  if (s1->info.CI.contigID < s2->info.CI.contigID)
    return((int)-1);

  if (s1->offsetAEnd.mean > s2->offsetAEnd.mean)
    return((int)1);
  if (s1->offsetAEnd.mean < s2->offsetAEnd.mean)
    return((int)-1);

  if (s1->offsetBEnd.mean > s2->offsetBEnd.mean)
    return((int)1);
  if (s1->offsetBEnd.mean < s2->offsetBEnd.mean)
    return((int)-1);

  return((int)0);
}


// Remove copies of surrogates which are placed multiple times in the same place in a contig

void
RemoveSurrogateDuplicates(void) {
  GraphNodeIterator chunks;
  ChunkInstanceT   *curChunk;

  fprintf(stderr, "RemoveSurrogateDuplicates()--\n");

  InitGraphNodeIterator(&chunks, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);

  while ((curChunk = NextGraphNodeIterator(&chunks)) != NULL) {
    if (curChunk->type != UNRESOLVEDCHUNK_CGW)
      continue;

    if (curChunk->info.CI.numInstances < 2)
      continue;

    if (curChunk->info.CI.numInstances == 2) {
      ChunkInstanceT *surr1 = GetGraphNode(ScaffoldGraph->CIGraph, curChunk->info.CI.instances.in_line.instance1);
      ChunkInstanceT *surr2 = GetGraphNode(ScaffoldGraph->CIGraph, curChunk->info.CI.instances.in_line.instance2);

      assert((surr1 != NULL));
      assert((surr2 != NULL));

      fprintf(stderr, "RemoveSurrogateDuplicates()-- surrogate CI=%d contig=%d pos=%f,%f\n",
              surr1->id,
              surr1->info.CI.contigID,
              surr1->offsetAEnd.mean,
              surr1->offsetBEnd.mean);
      fprintf(stderr, "RemoveSurrogateDuplicates()-- surrogate CI=%d contig=%d pos=%f,%f\n",
              surr2->id,
              surr2->info.CI.contigID,
              surr2->offsetAEnd.mean,
              surr2->offsetBEnd.mean);

      if ((surr1->info.CI.contigID >= 0) &&
          (surr1->info.CI.contigID == surr2->info.CI.contigID) &&
          (fabs(surr1->offsetAEnd.mean - surr2->offsetAEnd.mean) < 10.0) &&
          (fabs(surr1->offsetBEnd.mean - surr2->offsetBEnd.mean) < 10.0)) {

        curChunk->info.CI.numInstances = 1;
        curChunk->info.CI.instances.in_line.instance2 = -1;

        fprintf(stderr, "RemoveSurrogateDuplicates()--  Remove surrogate CI=%d from contig=%d\n", surr2->id, surr2->info.CI.contigID);

        RepairContigNeighbors(surr2);
        DeleteGraphNode(ScaffoldGraph->CIGraph, surr2);
      }
    }

    if (curChunk->info.CI.numInstances > 2) {
      if (curChunk->info.CI.numInstances != GetNumCDS_CID_ts(curChunk->info.CI.instances.va))
        fprintf(stderr, "RemoveSurrogateDuplicates()-- CI.numInstances=%d != instances.va=%d\n",
                curChunk->info.CI.numInstances,
                GetNumCDS_CID_ts(curChunk->info.CI.instances.va));
      assert(curChunk->info.CI.numInstances == GetNumCDS_CID_ts(curChunk->info.CI.instances.va));

      int numVaInstances = curChunk->info.CI.numInstances;

      qsort(GetCDS_CID_t(curChunk->info.CI.instances.va, 0), numVaInstances, sizeof(CDS_CID_t), CompareSurrogatePlacements);

      ChunkInstanceT *prevSurr = GetGraphNode(ScaffoldGraph->CIGraph, *GetCDS_CID_t(curChunk->info.CI.instances.va, 0));
      ChunkInstanceT *curSurr  = NULL;
      int i, copyto;

      assert(prevSurr != NULL);

      for(i = 0; i < numVaInstances; i++){
        curSurr = GetGraphNode(ScaffoldGraph->CIGraph, *GetCDS_CID_t(curChunk->info.CI.instances.va, i));

        fprintf(stderr, "RemoveSurrogateDuplicates()-- surrogate CI=%d contig=%d pos=%f,%f\n",
                curSurr->id,
                curSurr->info.CI.contigID,
                curSurr->offsetAEnd.mean,
                curSurr->offsetBEnd.mean);
      }

      for(i = 1, copyto = 1; i < numVaInstances; i++){
        curSurr = GetGraphNode(ScaffoldGraph->CIGraph, *GetCDS_CID_t(curChunk->info.CI.instances.va, i));

        assert(curSurr != NULL);

        if ((prevSurr->info.CI.contigID >= 0) &&
            (prevSurr->info.CI.contigID == curSurr->info.CI.contigID) &&
            (fabs(prevSurr->offsetAEnd.mean - curSurr->offsetAEnd.mean) < 10.0) &&
            (fabs(prevSurr->offsetBEnd.mean - curSurr->offsetBEnd.mean) < 10.0)) {

          fprintf(stderr, "RemoveSurrogateDuplicates()--  Remove surrogate CI=%d from contig=%d\n", curSurr->id, curSurr->info.CI.contigID);

          RepairContigNeighbors(curSurr);
          DeleteGraphNode(ScaffoldGraph->CIGraph, curSurr);

        } else {
          if (copyto != i)
            SetCDS_CID_t(curChunk->info.CI.instances.va, copyto, GetCDS_CID_t(curChunk->info.CI.instances.va, i));

          prevSurr = curSurr;
          copyto++;
        }
      }

      curChunk->info.CI.numInstances = copyto;

      assert(curChunk->info.CI.numInstances > 0);

      if (curChunk->info.CI.numInstances < 3){
        CDS_CID_t  a = *GetCDS_CID_t(curChunk->info.CI.instances.va, 0);
        CDS_CID_t  b = *GetCDS_CID_t(curChunk->info.CI.instances.va, 1);

        DeleteVA_CDS_CID_t(curChunk->info.CI.instances.va);

        curChunk->info.CI.instances.in_line.instance1 = a;
        curChunk->info.CI.instances.in_line.instance2 = (curChunk->info.CI.numInstances == 2) ? b : -1;

      } else {
        ResetToRange_CDS_CID_t(curChunk->info.CI.instances.va, curChunk->info.CI.numInstances);
        assert(curChunk->info.CI.numInstances == GetNumCDS_CID_ts(curChunk->info.CI.instances.va));
      }
    }
  }
}

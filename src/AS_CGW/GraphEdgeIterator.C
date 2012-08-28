
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2012, J. Craig Venter Institute.
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

static char *rcsid = "$Id: GraphEdgeIterator.C,v 1.8 2012-08-28 21:09:39 brianwalenz Exp $";

#include "GraphCGW_T.h"
#include "ScaffoldGraph_CGW.h"
//#include "GraphEdgeIterator.H"

#include <vector>
#include <algorithm>

using namespace std;


GraphEdgeIterator::GraphEdgeIterator(GraphCGW_T *graph_,
                                     CDS_CID_t   cid_,
                                     int         end_,
                                     int         edgeStatusSet_) {

  assert(graph_ != NULL);
  assert(end_ & ALL_END);

  graph              = graph_;
  cid                = cid_;
  end                = end_;
  edgeStatusSet      = edgeStatusSet_;

  noContains         = false;

  reset();

  //fprintf(stderr, "Created GraphEdgeIterator for cid %d end %d status %d nextM %d nextR %d\n",
  //        cid_, end_, edgeStatusSet_, nextM, nextR);
}


//  Return the next raw edge in the traversal.  Precondition: nextM is a pointer to the next merged
//  edge; nextR is to the next raw edge.  currM is not used.
//
EdgeCGW_T *
GraphEdgeIterator::nextRaw(void) {
  EdgeCGW_T *r = NULL;

 getRawEdge:

  //  If we're in a list of raw edges, keep marching down the list.
  if (nextR != NULLINDEX) {
    r = GetGraphEdge(graph, nextR);

    nextR = r->nextRawEdge;

    if ((r->idA != cid) && (r->idB != cid))
      fprintf(stderr, "ERROR: edge %16p idA=%d idB=%d not for cid %d\n",
              r, r->idA, r->idB, cid);

    assert(r->flags.bits.isRaw == true);
    assert(r->flags.bits.isDeleted == false);
    assert((r->idA == cid) || (r->idB == cid));

    if ((GetEdgeStatus(r) & edgeStatusSet) == 0)
      //  If the raw edge isn't marked appropraitely, get another.
      goto getRawEdge;

    //fprintf(stderr, "Return raw edge %d-%d\n", r->idA, r->idB);
    return(r);
  }

  //  No more top level edges?

  if (it == graph->edgeLists[cid].end())
    return(NULL);

  r = GetGraphEdge(graph, *it);

  it++;

  //  If this merged edge in fact is raw (it didn't have anything to merge with),
  //  return it.

  if (r->flags.bits.isRaw == true) {

    //  This assert was triggering in a variety of Salmon assemblies.  The top level edge (r) is
    //  marked as raw, and claimed one contributing edge.  The top level edge was an overlap edge,
    //  while the contributing edge was not.  The IDs did not agree.  There indeed was one
    //  raw edge in the nextRawEdge list.
    //
    //  It is possibly caused by deleting a raw edge from a 2-edge merged edge.
    //
    if (r->nextRawEdge != NULLINDEX) {
      fprintf(stderr, "GraphEdgeIterator()-- WARNING: top level raw edge %p (%d,%d contributing=%d) claims to have raw edges (id %d).  Fixing and ignoring.\n",
              r, r->idA, r->idB, r->edgesContributing, r->nextRawEdge);
      r->nextRawEdge = NULLINDEX;
    }
    assert(r->nextRawEdge == NULLINDEX);

    if ((GetEdgeStatus(r) & edgeStatusSet) == 0)
      //  If the raw edge isn't marked appropraitely, get another.
      goto getRawEdge;

    return(r);
  }

  //  Otherwise, it is a top level merged edge, so dive into the raw edges.

  nextR = r->nextRawEdge;

  goto getRawEdge;
}


EdgeCGW_T *
GraphEdgeIterator::nextMerged(void) {
  EdgeCGW_T *r = NULL;

 getMergedEdge:
  if (it == graph->edgeLists[cid].end())
    return(NULL);

  r = GetGraphEdge(graph, *it);

  assert(r->topLevelEdge == *it);
  assert(r->flags.bits.isDeleted == false);
  assert((r->idA == cid) || (r->idB == cid));

  it++;

  PairOrient orient = GetEdgeOrientationWRT(r, cid);

  assert(orient.isUnknown() == false);

  //  Check for correct end and confirmed status

  if ((GetEdgeStatus(r) & edgeStatusSet) == 0) {
    //fprintf(stderr, "Merged edge eid %d %d-%d status %d not %d\n",
    //        GetVAIndex_EdgeCGW_T(graph->edges, r), r->idA, r->idB, GetEdgeStatus(r), edgeStatusSet);
    goto getMergedEdge;
  }

  if ((end == A_END) && (orient.isAB_BA() || orient.isAB_AB())) {
    //fprintf(stderr, "Merged edge eid %d %d-%d orient %c not A_END\n",
    //        GetVAIndex_EdgeCGW_T(graph->edges, r), r->idA, r->idB, orient.toLetter());
    goto getMergedEdge;
  }

  if ((end == B_END) && (orient.isBA_BA() || orient.isBA_AB())) {
    //fprintf(stderr, "Merged edge eid %d %d-%d orient %c not B_END\n",
    //        GetVAIndex_EdgeCGW_T(graph->edges, r), r->idA, r->idB, orient.toLetter());
    goto getMergedEdge;
  }

  if ((noContains == true) && (isContainmentEdge(r))) {
    //fprintf(stderr, "Merged edge eid %d %d-%d is contained, skip\n",
    //        GetVAIndex_EdgeCGW_T(graph->edges, r), r->idA, r->idB);
    goto getMergedEdge;
  }

  //fprintf(stderr, "Return merged edge eid %d %d-%d\n",
  //        GetVAIndex_EdgeCGW_T(graph->edges, r), r->idA, r->idB);
  return(r);
}




bool
edgeCompareForMerging(EdgeCGW_T const *edge1, EdgeCGW_T const *edge2) {
  int diff;

  assert(edge1->idA <= edge1->idB);
  assert(edge2->idA <= edge2->idB);
  assert(edge1->flags.bits.isDeleted == false);
  assert(edge2->flags.bits.isDeleted == false);

  if (edge1->idA < edge2->idA)
    return(true);
  if (edge1->idA > edge2->idA)
    return(false);

  if (edge1->idB < edge2->idB)
    return(true);
  if (edge1->idB > edge2->idB)
    return(false);

  if (edge1->orient.toLetter() < edge2->orient.toLetter())
    return(true);
  if (edge1->orient.toLetter() > edge2->orient.toLetter())
    return(false);

  // We want guide edges AFTER all other edges
  diff = isSloppyEdge(edge1) - isSloppyEdge(edge2);
  if (diff < 0)
    return(true);
  if (diff > 0)
    return(false);

  // We want overlap edges AFTER non-overlap edges...
  diff = isOverlapEdge(edge1) - isOverlapEdge(edge2);
  if (diff < 0)
    return(true);
  if (diff > 0)
    return(false);

  if (edge1->distance.mean < edge2->distance.mean)
    return(true);
  if (edge1->distance.mean > edge2->distance.mean)
    return(false);

  //  This *is* well behaved.  We should only be comparing edges allocated from the same array,
  //  and this is essentially just comparing the edge IDs.  It is needed so that we can add two
  //  (raw) edges of the same stats.
  return(edge1 < edge2);

  return(false);
}



bool
edgeCompareForStoring(EdgeCGW_T const *edge1, EdgeCGW_T const *edge2) {
  int diff;

  assert(edge1->idA <= edge1->idB);
  assert(edge2->idA <= edge2->idB);
  assert(edge1->flags.bits.isDeleted == false);
  assert(edge2->flags.bits.isDeleted == false);

  if (edge1->idA < edge2->idA)
    return(true);
  if (edge1->idA > edge2->idA)
    return(false);

  if (edge1->idB < edge2->idB)
    return(true);
  if (edge1->idB > edge2->idB)
    return(false);

  if (edge1->orient.toLetter() < edge2->orient.toLetter())
    return(true);
  if (edge1->orient.toLetter() > edge2->orient.toLetter())
    return(false);

  //  This *is* well behaved.  We should only be comparing edges allocated from the same array,
  //  and this is essentially just comparing the edge IDs.  It is needed so that we can add two
  //  (raw) edges of the same stats.
  return(edge1 < edge2);
}

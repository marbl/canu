
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

static char *rcsid = "$Id: GraphEdgeIterator.C,v 1.6 2012-08-19 04:13:52 brianwalenz Exp $";

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

  if (GetGraphNode(graph, cid)->flags.bits.edgesModified)
    sortEdges();

  prevM              = NULLINDEX;
  currM              = NULLINDEX;
  nextM              = GetGraphNode(graph, cid)->edgeHead;

  nextR              = NULLINDEX;

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

  //  If there is no nextM, we're done with the iteration.
  if (nextM == NULLINDEX)
    return(NULL);

  //  Otherwise, grab the nextM and do soemthing with it.  Whatever happens,
  //  we want to also advance to the next nextM.

  r = GetGraphEdge(graph, nextM);

  nextM = (r->idA == cid) ? r->nextALink : r->nextBLink;
  nextR = NULLINDEX;

  //  Previous versions would check that the merged edge was marked appropriately, and if not, skip
  //  ALL raw edges below it.  Starting in v1.2 we test the raw edges separately.

  //  If this merged edge in fact is raw (it didn't have anything to merge with),
  //  return it.

  if (r->flags.bits.isRaw == true) {

    //  This assert was triggering in a variety of Salmon assemblies.  The top level edge (r) is
    //  marked as raw, and claimed one contributing edge.  The top level edge was an overlap edge,
    //  while the contributing edge was not.  The IDs did not agree.  There indeed was one
    //  raw edge in the nextRawEdge list.
    //
    //  Instead of tracking down what caused this, we just 'fix' the top level edge to have
    //  no underlying edges.
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

  //  Like, if it is a merged edge, start iterating down the raw edges.

  nextR = r->nextRawEdge;

  goto getRawEdge;
}


EdgeCGW_T *
GraphEdgeIterator::nextMerged(void) {
  EdgeCGW_T *r = NULL;

 getMergedEdge:
  prevM = currM;
  currM = nextM;
  nextM = NULLINDEX;

  if (currM == NULLINDEX)
    return(NULL);

  r = GetGraphEdge(graph, currM);

  nextM = (r->idA == cid) ? r->nextALink : r->nextBLink;

  if ((r->idA != cid) && (r->idB != cid))
    fprintf(stderr, "ERROR: edge %16p idA=%d idB=%d not for cid %d\n",
            r, r->idA, r->idB, cid);

  //  This is not true - we will return raw edges, as merged edges with only one raw edge.
  //assert(r->flags.bits.isRaw == false);
  assert(r->flags.bits.isDeleted == false);
  assert((r->idA == cid) || (r->idB == cid));

  // Flip orientation to accomodate canonical graph form

  PairOrient orient = GetEdgeOrientationWRT(r, cid);

  assert(orient.isUnknown() == false);

  //  Check for correct end and confirmed status

  if ((GetEdgeStatus(r) & edgeStatusSet) == 0) {
    //fprintf(stderr, "Merged edge %d-%d status %d not %d\n",
    //        r->idA, r->idB, GetEdgeStatus(r), edgeStatusSet);
    goto getMergedEdge;
  }

  if ((end == A_END) && (orient.isAB_BA() || orient.isAB_AB())) {
    //fprintf(stderr, "Merged edge %d-%d orient %c not A_END\n",
    //        r->idA, r->idB, orient.toLetter());
    goto getMergedEdge;
  }

  if ((end == B_END) && (orient.isBA_BA() || orient.isBA_AB())) {
    //fprintf(stderr, "Merged edge %d-%d orient %c not B_END\n",
    //        r->idA, r->idB, orient.toLetter());
    goto getMergedEdge;
  }

  if ((noContains == true) && (isContainmentEdge(r))) {
    //fprintf(stderr, "Merged edge %d-%d is contained, skip\n",
    //        r->idA, r->idB);
    goto getMergedEdge;
  }

  //fprintf(stderr, "Return merged edge %d-%d p=%d c=%d n=%d\n",
  //        r->idA, r->idB, prevM, currM, nextM);
  return(r);
}




bool
edgeCompare(EdgeCGW_T const *edge1, EdgeCGW_T const *edge2) {
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

  diff = edge1->distance.mean - edge2->distance.mean;
  if (diff < 0)
    return(true);

  return(false);
}


void
GraphEdgeIterator::sortEdges(void) {
  NodeCGW_T *node = GetGraphNode(graph, cid);
  EdgeCGW_T *edge = NULL;

  if (node->flags.bits.edgesModified == 0)
    return;

#if 0
  char const *type = "UNKNOWN";
  type = (graph == ScaffoldGraph->CIGraph)       ? "unitig" : type;
  type = (graph == ScaffoldGraph->ContigGraph)   ? "contig" : type;
  type = (graph == ScaffoldGraph->ScaffoldGraph) ? "scaffold" : type;

  fprintf(stderr, "GraphEdgeIterator::sortEdges()--  sorting edges for %s %d\n", type, node->id);
#endif

  node->flags.bits.edgesModified = 0;

  if (node->edgeHead == NULLINDEX)
    return;

  //  Move pointers from the unsorted list to our vector for sorting

  vector<EdgeCGW_T *>  edges;

  for (CDS_CID_t ie=node->edgeHead; ie != NULLINDEX; ) {
    edge = GetGraphEdge(graph, ie);

    edges.push_back(edge);

    if ((edge->idA != cid) && (edge->idB != cid))
      fprintf(stderr, "ERROR: edge %16p idA=%d idB=%d not for cid %d\n",
              edge, edge->idA, edge->idB, cid);
    assert((edge->idA == cid) || (edge->idB == cid));

    ie = (edge->idA == cid) ? edge->nextALink : edge->nextBLink;

    //fprintf(stderr, "GraphEdgeIterator::sortEdges()--   edge %16p idx %d next %d\n",
    //        edge, GetVAIndex_EdgeCGW_T(graph->edges, edge), ie);
  }

  //  Sort

#ifdef _GLIBCXX_PARALLEL
  //  If we have the parallel STL, don't use it!  We don't want the expense of firing up threads here.
  __gnu_sequential::sort(edges.begin(), edges.end(), edgeCompare);
#else
  sort(edges.begin(), edges.end(), edgeCompare);
#endif

  //  Update pointers.

  //for (uint32 ei=0; ei<edges.size(); ei++) {
  //  edges[ei]->prevALink = NULLINDEX;
  //  edges[ei]->nextALink = NULLINDEX;
  //  edges[ei]->prevBLink = NULLINDEX;
  //  edges[ei]->nextBLink = NULLINDEX;
  //}

  CDS_CID_t  last = NULLINDEX;

  //  Move forward, to set backward pointers

  for (uint32 ei=0; ei<edges.size(); ei++) {
    EdgeCGW_T *edge = edges[ei];

    if (edge->idA == cid)
      edge->prevALink = last;
    else
      edge->prevBLink = last;

    last = GetVAIndex_EdgeCGW_T(graph->edges, edge);
  }

  //  Move backward, to set forward pointers

  last = NULLINDEX;

  for (uint32 ei=edges.size(); ei-->0; ) {
    EdgeCGW_T *edge = edges[ei];

    if (edge->idA == cid)
      edge->nextALink = last;
    else
      edge->nextBLink = last;

    last = GetVAIndex_EdgeCGW_T(graph->edges, edge);
  }

  //  Set the list head to the node.

  node->edgeHead = GetVAIndex_EdgeCGW_T(graph->edges, edges[0]);
}

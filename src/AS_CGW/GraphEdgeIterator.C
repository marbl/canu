
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

static char *rcsid = "$Id: GraphEdgeIterator.C,v 1.2 2012-08-03 21:14:14 brianwalenz Exp $";

#include "GraphCGW_T.h"
//#include "GraphEdgeIterator.H"


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

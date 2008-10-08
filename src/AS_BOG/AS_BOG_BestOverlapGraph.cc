
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id: AS_BOG_BestOverlapGraph.cc,v 1.64 2008-10-08 22:02:54 brianwalenz Exp $";

#include<iostream>
#include<vector>
#include<limits>
#include<cmath>

#include "AS_BOG_BestOverlapGraph.hh"

#undef max

//  The overlap, pi, exists between A and B:
//
//  A -------------->
//         |||||||||
//  B      ---------------->
//
//  AEnd(pi) is 3'
//  BEnd(pi) is 5'
//

fragment_end_type BestOverlapGraph::AEnd(const OVSoverlap& olap) {
  if (olap.dat.ovl.a_hang < 0 && olap.dat.ovl.b_hang < 0)
    return FIVE_PRIME;
  if (olap.dat.ovl.a_hang > 0 && olap.dat.ovl.b_hang > 0)
    return THREE_PRIME;

  assert(0); // no contained
}

fragment_end_type BestOverlapGraph::BEnd(const OVSoverlap& olap) {
  if (olap.dat.ovl.a_hang < 0 && olap.dat.ovl.b_hang < 0)
    if (olap.dat.ovl.flipped)
      return FIVE_PRIME;
    else
      return THREE_PRIME;

  if (olap.dat.ovl.a_hang > 0 && olap.dat.ovl.b_hang > 0)
    if (olap.dat.ovl.flipped)
      return THREE_PRIME;
    else
      return FIVE_PRIME;

  assert(0); // no contained
}


// Create BestOverlapGraph as an array of size max fragments.
//     Assuming that our iuids start at index value of 1.
//
// AS_UTG_ERROR_RATE is fraction error, same as AS_CNS_ERROR_RATE.
//
BestOverlapGraph::BestOverlapGraph(FragmentInfo        *fi,
                                   OverlapStore        *ovlStore,
                                   double               AS_UTG_ERROR_RATE) {
  _fi = fi;

  _best_overlaps = new BestFragmentOverlap [fi->numFragments() + 1];
  _best_contains = new BestContainment     [fi->numFragments() + 1];

  memset(_best_overlaps, 0, sizeof(BestFragmentOverlap) * (fi->numFragments() + 1));
  memset(_best_contains, 0, sizeof(BestContainment)     * (fi->numFragments() + 1));

  assert(AS_UTG_ERROR_RATE >= 0.0);
  assert(AS_UTG_ERROR_RATE <= AS_MAX_ERROR_RATE);

  assert(AS_CNS_ERROR_RATE >= 0.0);
  assert(AS_CNS_ERROR_RATE <= AS_MAX_ERROR_RATE);

  fprintf(stderr, "BestOverlapGraph()-- UTG erate %.4f%%, CNS erate %.4f%%\n",
          100.0 * AS_UTG_ERROR_RATE, 100.0 * AS_CNS_ERROR_RATE);

  mismatchCutoff  = AS_OVS_encodeQuality(AS_UTG_ERROR_RATE);
  consensusCutoff = AS_OVS_encodeQuality(AS_CNS_ERROR_RATE);

  // Go through the overlap stream in two passes:
  // first pass finds the containments
  // second pass builds the overlap graph, excluding contained frags

  OVSoverlap olap;

  AS_OVS_resetRangeOverlapStore(ovlStore);
  while  (AS_OVS_readOverlapFromStore(ovlStore, &olap, AS_OVS_TYPE_OVL))
    scoreContainment(olap);

  AS_OVS_resetRangeOverlapStore(ovlStore);
  while  (AS_OVS_readOverlapFromStore(ovlStore, &olap, AS_OVS_TYPE_OVL))
    scoreEdge(olap);

  //  Diagnostic.  Dump the best edges, count the number of contained
  //  reads, etc.
  {
    FILE *BC = fopen("best.contains", "w");
    FILE *BE = fopen("best.edges", "w");

    if ((BC) && (BE)) {
      fprintf(BC, "#fragId\tlibId\tmated\tbestCont\n");
      fprintf(BE, "#fragId\tlibId\tmated\tbest5\tbest3\n");

      for (int id=1; id<_fi->numFragments() + 1; id++) {
        BestContainment *bestcont  = getBestContainer(id);
        BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, FIVE_PRIME);
        BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, THREE_PRIME);

        if (bestcont)
          fprintf(BC, "%d\t%d\t%c\t%d\n", id, _fi->libraryIID(id), (_fi->mateIID(id) > 0) ? 'm' : 'f', bestcont->container);
        else if ((bestedge5->frag_b_id > 0) || (bestedge3->frag_b_id > 0))
          fprintf(BE, "%d\t%d\t%d\t%d\n", id, _fi->libraryIID(id), bestedge5->frag_b_id, bestedge3->frag_b_id);
      }

      fclose(BC);
      fclose(BE);
    }
  }
}

BestOverlapGraph::~BestOverlapGraph(){
  delete[] _best_overlaps;
  delete[] _best_contains;
}


//  Given a fragment IUID and which end, returns pointer to
//  BestOverlap node.
BestEdgeOverlap *BestOverlapGraph::getBestEdgeOverlap(iuid frag_id, fragment_end_type which_end){
  if(which_end == FIVE_PRIME)
    return(&_best_overlaps[frag_id].five_prime);
  if(which_end == THREE_PRIME)
    return(&_best_overlaps[frag_id].three_prime);
  return(NULL);
}

BestEdgeOverlap *BestOverlapGraph::getBestEdgeOverlap(FragmentEnd* end) {
  return getBestEdgeOverlap(end->fragId(),end->fragEnd());
}

void BestOverlapGraph::followOverlap(FragmentEnd* end) {
  BestEdgeOverlap* edge = getBestEdgeOverlap(end);
  *end = FragmentEnd(edge->frag_b_id, (edge->bend == FIVE_PRIME) ? THREE_PRIME : FIVE_PRIME);
}

bool BestOverlapGraph::containHaveEdgeTo(iuid contain, iuid otherRead) {
  BestContainment  *c = &_best_contains[contain];
  bool              r = false;

  if (c->isContained) {
    if (c->olapsLen < 16) {
      for (int i=0; i<c->olapsLen; i++)
        if (c->olaps[i] == otherRead) {
          r = true;
          break;
        }
    } else {
      if (c->olapsSorted == false) {
        std::sort(c->olaps, c->olaps + c->olapsLen);
        c->olapsSorted = true;
      }
      r = std::binary_search(c->olaps, c->olaps + c->olapsLen, otherRead);
    }
  }

  return(r);
}


void BestOverlapGraph::scoreContainment(const OVSoverlap& olap) {

  //  Count the number of overlaps we have -- used by scoreEdge to
  //  keep a list of dovetail overlaps to contained fragments.
  //
  //  GOOFY!  We're counting on A, but the rest of the routine uses B.
  //
  if ((olap.dat.ovl.a_hang < 0 && olap.dat.ovl.b_hang < 10) ||
      (olap.dat.ovl.a_hang > 0 && olap.dat.ovl.b_hang > -10))
    _best_contains[olap.a_iid].olapsMax++;

  //  In the case of no hang, make the lower frag the container
  //
  if ((olap.dat.ovl.a_hang == 0) &&
      (olap.dat.ovl.b_hang == 0) &&
      (olap.a_iid > olap.b_iid))
    return;

  //  We only care if A contains B.

  if ((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0)) {
    float             newScr = scoreOverlap(olap);
    BestContainment  *c = &_best_contains[olap.b_iid];

    if ((newScr > 0) &&
        (newScr > c->contain_score)) {
      //  NOTE!  This is already initialized.  We do not need to, and
      //  it is an error to, initialize olaps to zero!  (We're
      //  counting olapsMax above, see?  This stupid bug took me about
      //  an hour to find, grrr.)
      c->container         = olap.a_iid;
      c->contain_score     = newScr;
      c->a_hang            = olap.dat.ovl.a_hang;
      c->b_hang            = olap.dat.ovl.b_hang;
      c->sameOrientation   = olap.dat.ovl.flipped ? false : true;
      c->isContained       = true;
      c->isPlaced          = false;
      c->olapsSorted       = false;
    }
  }
}

void BestOverlapGraph::scoreEdge(const OVSoverlap& olap) {

  //  Store real edges from contained frags to help with unhappy
  //  mate splitting
  //
  //  If A was previously found to be a containee, and
  //
  //  a < 0           b < 10
  //        ----------..
  //  ......------------
  //
  //  a > 0           b > -10
  //  ......-----------
  //        ---------.......
  //
  //  That looks wrong.
  //
  //  From Eli: These are contained, but close either way.  We're
  //  storing the non-containment edges for this fragment, plus a
  //  few containment edges that are "close" to being dovetails.  "I
  //  think there are cases when a change in the alignemtn
  //  (consensus) will change which one is contained and screw up
  //  the order, so having this 10 base fudge factor helps things
  //  work out."
  //
  if (isContained(olap.a_iid)) {
    if ((olap.dat.ovl.a_hang < 0 && olap.dat.ovl.b_hang < 10) ||
        (olap.dat.ovl.a_hang > 0 && olap.dat.ovl.b_hang > -10)) {
      BestContainment *c = &_best_contains[olap.a_iid];
      if (c->olaps == NULL)
        c->olaps = new iuid [c->olapsMax];
      c->olaps[c->olapsLen++] = olap.b_iid;
    }
    return;
  }

  //  Skip contained fragments.
  if (isContained(olap.a_iid) || isContained(olap.b_iid))
    return;

  //  Skip containment overlaps.  Can this happen?  Yup.  How?
  //  The overlap could be above our allowed error.
  //
  if (((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0)) ||
      ((olap.dat.ovl.a_hang <= 0) && (olap.dat.ovl.b_hang >= 0)))
    return;

  float newScr = scoreOverlap(olap);

  //  If the score is 0, the overlap doesn't pass the scoring
  //  criteria at all so don't store the overlap whether or not
  //  it's dovetailing or containment.
  if (newScr <= 0)
    return;

  //  Dove tailing overlap

  BestEdgeOverlap *best = getBestEdgeOverlap(olap.a_iid, AEnd(olap));
  short            olapLen = olapLength(olap);

  // Store the overlap if:
  //   1.)  The score is better than what is already in the graph
  //   2.)  If the scores are identical, the one with the longer length
  //
  // Since the order of how the overlaps are read in from the overlap
  // store are by A's increasing iuid, by default, if the score and
  // length are the same, the iuid of the lower value will be kept.

  if ((newScr > best->olap_score) ||
      ((newScr == best->olap_score) && (olapLen > best->olap_length))) {
    best->frag_b_id    = olap.b_iid;
    best->olap_score   = newScr;
    best->olap_length  = olapLen;
    best->bend         = BEnd(olap);
    best->ahang        = olap.dat.ovl.a_hang;
    best->bhang        = olap.dat.ovl.b_hang;
  }
}

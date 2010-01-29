
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

static const char *rcsid = "$Id: AS_BOG_BestOverlapGraph.cc,v 1.72 2010-01-29 07:36:36 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
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

uint32 BestOverlapGraph::AEnd(const OVSoverlap& olap) {
  if (olap.dat.ovl.a_hang < 0 && olap.dat.ovl.b_hang < 0)
    return FIVE_PRIME;
  if (olap.dat.ovl.a_hang > 0 && olap.dat.ovl.b_hang > 0)
    return THREE_PRIME;

  assert(0); // no contained
  return(0);
}

uint32 BestOverlapGraph::BEnd(const OVSoverlap& olap) {
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
  return(0);
}


// Create BestOverlapGraph as an array of size max fragments.
//     Assuming that our uint32s start at index value of 1.
//
// AS_UTG_ERROR_RATE is fraction error, same as AS_CNS_ERROR_RATE.
//
BestOverlapGraph::BestOverlapGraph(FragmentInfo        *fi,
                                   OverlapStore        *ovlStoreUniq,
                                   OverlapStore        *ovlStoreRept,
                                   double               AS_UTG_ERROR_RATE,
                                   double               AS_UTG_ERROR_LIMIT) {
  OVSoverlap olap;

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

  mismatchLimit   = AS_UTG_ERROR_LIMIT;

  //  Pass 1 through overlaps -- find the contained fragments.

  //setLogFile("unitigger", "bestoverlapgraph-containments");

  _best_contains_score    = new uint64 [fi->numFragments() + 1];
  memset(_best_contains_score,    0, sizeof(uint64) * (fi->numFragments() + 1));

  AS_OVS_resetRangeOverlapStore(ovlStoreUniq);
  while  (AS_OVS_readOverlapFromStore(ovlStoreUniq, &olap, AS_OVS_TYPE_OVL))
    scoreContainment(olap);

  if (ovlStoreRept) {
    AS_OVS_resetRangeOverlapStore(ovlStoreRept);
    while  (AS_OVS_readOverlapFromStore(ovlStoreRept, &olap, AS_OVS_TYPE_OVL))
      scoreContainment(olap);
  }

  delete [] _best_contains_score;
  _best_contains_score    = NULL;

  //  Report some statistics on overlaps
  {
    uint64  numContainsToSave = 0;

    for (uint32 i=0; i<fi->numFragments(); i++)
      numContainsToSave += _best_contains[i].olapsLen;

    fprintf(stderr, "Need to save "F_U64" near-containment overlaps\n", numContainsToSave);
  }


  //  Pass 2 through overlaps -- find dovetails, build the overlap graph.  For each
  //  contained fragment, remember some of the almost containment overlaps.

  //setLogFile("unitigger", "bestoverlapgraph-dovetails");

  _best_overlaps_5p_score = new uint64 [fi->numFragments() + 1];
  _best_overlaps_3p_score = new uint64 [fi->numFragments() + 1];

  memset(_best_overlaps_5p_score, 0, sizeof(uint64) * (fi->numFragments() + 1));
  memset(_best_overlaps_3p_score, 0, sizeof(uint64) * (fi->numFragments() + 1));

  AS_OVS_resetRangeOverlapStore(ovlStoreUniq);
  while  (AS_OVS_readOverlapFromStore(ovlStoreUniq, &olap, AS_OVS_TYPE_OVL))
    scoreEdge(olap);

  if (ovlStoreRept) {
    AS_OVS_resetRangeOverlapStore(ovlStoreRept);
    while  (AS_OVS_readOverlapFromStore(ovlStoreRept, &olap, AS_OVS_TYPE_OVL))
      scoreEdge(olap);
  }

  delete [] _best_overlaps_5p_score;
  delete [] _best_overlaps_3p_score;

  _best_overlaps_5p_score = NULL;
  _best_overlaps_3p_score = NULL;

  //setLogFile("unitigger", NULL);

  //  Diagnostic.  Dump the best edges, count the number of contained
  //  reads, etc.
  {
    FILE *BC = fopen("best.contains", "w");
    FILE *BE = fopen("best.edges", "w");

    if ((BC) && (BE)) {
      fprintf(BC, "#fragId\tlibId\tmated\tbestCont\n");
      fprintf(BE, "#fragId\tlibId\tmated\tbest5\tbest3\n");

      for (uint32 id=1; id<_fi->numFragments() + 1; id++) {
        BestContainment *bestcont  = getBestContainer(id);
        BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, FIVE_PRIME);
        BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, THREE_PRIME);

        if (bestcont)
          fprintf(BC, "%u\t%u\t%c\t%u\n", id, _fi->libraryIID(id), (_fi->mateIID(id) > 0) ? 'm' : 'f', bestcont->container);
        else if ((bestedge5->frag_b_id > 0) || (bestedge3->frag_b_id > 0))
          fprintf(BE, "%u\t%u\t%u\t%u\n", id, _fi->libraryIID(id), bestedge5->frag_b_id, bestedge3->frag_b_id);
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



void BestOverlapGraph::scoreContainment(const OVSoverlap& olap) {

  if (isOverlapBadQuality(olap))
    return;

  //  Count the number of good containment and near-containment
  //  overlaps the B fragment has -- used by scoreEdge (see the
  //  comment there) to keep a list of dovetail overlaps to contained
  //  fragments.
  //
  if (((olap.dat.ovl.a_hang >= -10) && (olap.dat.ovl.b_hang <=  0)) ||
      ((olap.dat.ovl.a_hang >=   0) && (olap.dat.ovl.b_hang <= 10)))
    _best_contains[olap.b_iid].olapsLen++;

  //  In the case of no hang, make the lower frag the container
  //
  if ((olap.dat.ovl.a_hang == 0) &&
      (olap.dat.ovl.b_hang == 0) &&
      (olap.a_iid > olap.b_iid))
    return;

  //  We only care if A contains B.

  if ((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0)) {
    uint64           newScr = scoreOverlap(olap);
    BestContainment      *c = &_best_contains[olap.b_iid];

    if (newScr > _best_contains_score[olap.b_iid]) {
      //  NOTE!  This is already initialized.  We do not need to, and
      //  it is an error to, initialize olaps to zero!  (We're
      //  counting olapsLen above, see?  This stupid bug took me about
      //  an hour to find, grrr.)
      c->container         = olap.a_iid;
      c->a_hang            = olap.dat.ovl.a_hang;
      c->b_hang            = olap.dat.ovl.b_hang;
      c->sameOrientation   = olap.dat.ovl.flipped ? false : true;
      c->isContained       = true;
      c->isPlaced          = false;
      c->olapsSorted       = false;

      _best_contains_score[olap.b_iid] = newScr;
    }
  }
}

void BestOverlapGraph::scoreEdge(const OVSoverlap& olap) {

  if (isOverlapBadQuality(olap))
    return;

  //  Store edges from contained frags to help with unhappy mate
  //  splitting.
  //
  //  From Eli: These are contained, but close either way.  We're
  //  storing the non-containment edges for this fragment, plus a few
  //  containment edges that are "close" to being dovetails.  "I think
  //  there are cases when a change in the alignemtn (consensus) will
  //  change which one is contained and screw up the order, so having
  //  this 10 base fudge factor helps things work out."
  //
  if (isContained(olap.b_iid)) {
    if (((olap.dat.ovl.a_hang >= -10) && (olap.dat.ovl.b_hang <=  0)) ||
        ((olap.dat.ovl.a_hang >=   0) && (olap.dat.ovl.b_hang <= 10))) {
      BestContainment *c = &_best_contains[olap.b_iid];
      if (c->olaps == NULL) {
        c->olaps    = new uint32 [c->olapsLen];
        c->olapsLen = 0;
      }
      c->olaps[c->olapsLen++] = olap.a_iid;
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

  uint64 newScr = scoreOverlap(olap);

  //  If the score is 0, the overlap doesn't pass the scoring
  //  criteria at all so don't store the overlap whether or not
  //  it's dovetailing or containment.

  if (newScr == 0)
    return;

  //  Dove tailing overlap
  uint32           aend    = AEnd(olap);
  BestEdgeOverlap *best    = getBestEdgeOverlap(olap.a_iid, aend);
  uint64           score   = 0;

  // Store the overlap if:
  //   1.)  The score is better than what is already in the graph
  //   2.)  If the scores are identical, the one with the longer length
  //
  // Since the order of how the overlaps are read in from the overlap
  // store are by A's increasing uint32, by default, if the score and
  // length are the same, the uint32 of the lower value will be kept.

  if (aend == THREE_PRIME)
    score = _best_overlaps_3p_score[olap.a_iid];
  else
    score = _best_overlaps_5p_score[olap.a_iid];

  if (newScr > score) {
    best->frag_b_id    = olap.b_iid;
    best->bend         = BEnd(olap);
    best->ahang        = olap.dat.ovl.a_hang;
    best->bhang        = olap.dat.ovl.b_hang;

    if (aend == THREE_PRIME)
      _best_overlaps_3p_score[olap.a_iid] = newScr;
    else
      _best_overlaps_5p_score[olap.a_iid] = newScr;
  }
}


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

static const char *rcsid = "$Id: AS_BAT_Joining.C,v 1.1 2010-11-24 01:03:31 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_UnitigGraph.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "MultiAlignStore.h"


class joinEntry {
public:
  joinEntry() {
    toFragID = 0;
    frFragID = 0;
    toLen    = 0;
    frLen    = 0;
    ttLen    = 0;
  };
  joinEntry(uint32 toFragID_, uint32 frFragID_, uint32 toLen_, uint32 frLen_) {
    toFragID = toFragID_;
    frFragID = frFragID_;
    toLen    = toLen_;
    frLen    = frLen_;
    ttLen    = toLen_ + frLen_;
  };
  ~joinEntry() {
  };

  uint32            toFragID;
  uint32            frFragID;

  uint32            toLen;
  uint32            frLen;
  uint32            ttLen;

  bool operator<(const joinEntry &that) const { return(ttLen < that.ttLen); };
  bool operator>(const joinEntry &that) const { return(ttLen > that.ttLen); };
};




//  Examine the first (few?) fragments of a unitig, evaluate if they indicate a join should be made.
void
UnitigGraph::joinUnitigs_examineFirst(Unitig *fr,
                                      vector<joinEntry> &joins) {
  ufNode          *frg         = &fr->ufpath[0];
  bool             fragForward = (frg->position.bgn < frg->position.end);
  BestEdgeOverlap *bestEdge    = OG->getBestEdgeOverlap(frg->ident, (fragForward == false));

  if (bestEdge->fragId() == 0)
    //  No best edge?  Skip it.
    return;

  uint32  toID = fr->fragIn(bestEdge->fragId());
  Unitig *to   = unitigs[toID];

  if (to->ufpath.size() == 1)
    //  Joining to something teeny?  Skip it.
    return;

  if (to->id() == fr->id())
    //  Join to myself?  Nope.
    return;

  //  Figure out which fragment we have an edge to, and which end of it is sticking out from the
  //  end of the unitig.

  ufNode *tgt      = NULL;
  bool    tgt3p    = false;

  if (bestEdge->fragId() == to->ufpath[0].ident) {
    tgt      = &to->ufpath[0];
    tgt3p    = (tgt->position.bgn > tgt->position.end);
  }

  if (bestEdge->fragId() == to->ufpath[to->ufpath.size()-1].ident) {
    tgt      = &to->ufpath[to->ufpath.size()-1];
    tgt3p    = (tgt->position.bgn < tgt->position.end);
  }

  if (tgt == NULL)
    //  Attempt to join to a fragment that isn't first or last.  Skip it.
    return;

  if (bestEdge->frag3p() != tgt3p)
    //  Joining to the wrong end of the first/last fragment?  Skip it.
    return;

  //  OK, join.

  joins.push_back(joinEntry(tgt->ident, frg->ident, to->getLength(), fr->getLength()));
}




//  Examine the last (few?) fragments of a unitig, evaluate if they indicate a join should be made.
//  This is complicated by wanting the last non-contained fragment.
void
UnitigGraph::joinUnitigs_examineLast(Unitig *fr,
                                     vector<joinEntry> &joins) {
  ufNode          *frg         = &fr->ufpath[fr->ufpath.size()-1];
  bool             fragForward = (frg->position.bgn < frg->position.end);
  BestEdgeOverlap *bestEdge    = OG->getBestEdgeOverlap(frg->ident, (fragForward == true));

  if (bestEdge->fragId() == 0)
    //  No best edge?  Skip it.
    return;

  uint32  toID = fr->fragIn(bestEdge->fragId());
  Unitig *to   = unitigs[toID];

  if (to->ufpath.size() == 1)
    //  Joining to something teeny?  Skip it.
    return;

  if (to->id() == fr->id())
    //  Join to myself?  Nope.
    return;

  //  Figure out which fragment we have an edge to, and which end of it is sticking out from the
  //  end of the unitig.

  ufNode *tgt      = NULL;
  uint32  tgt3p    = false;

  if (bestEdge->fragId() == to->ufpath[0].ident) {
    tgt      = &to->ufpath[0];
    tgt3p    = (tgt->position.bgn > tgt->position.end);
  }

  if (bestEdge->fragId() == to->ufpath[to->ufpath.size()-1].ident) {
    tgt      = &to->ufpath[to->ufpath.size()-1];
    tgt3p    = (tgt->position.bgn < tgt->position.end);
  }

  if (tgt == NULL)
    //  Attempt to join to a fragment that isn't first or last.  Skip it.
    return;

  if (bestEdge->frag3p() != tgt3p)
    //  Joining to the wrong end of the first/last fragment?  Skip it.
    return;

  //  OK, join.

  joins.push_back(joinEntry(tgt->ident, frg->ident, to->getLength(), fr->getLength()));
}




bool
UnitigGraph::joinUnitigs_confirmJoin(joinEntry *join) {
  uint32    frUnitigID = Unitig::fragIn(join->frFragID);
  uint32    toUnitigID = Unitig::fragIn(join->toFragID);

  if (frUnitigID == toUnitigID)
    //  Circular join?  Nope!  Won't do it.  Should never happen, really.
    return(false);

  Unitig   *frUnitig = unitigs[frUnitigID];
  Unitig   *toUnitig = unitigs[toUnitigID];

  if ((frUnitig == NULL) || (toUnitig == NULL))
    //  Already added one of the unitigs to another one.  Should never happen, really.
    return(false);

  uint32  frIdx = Unitig::pathPosition(join->frFragID);
  uint32  toIdx = Unitig::pathPosition(join->toFragID);

  //  If the fragment is not the first, reverse complement the unitig.  This (should) ensure that
  //  the fragments in the edge we want to join with are both at position 0 in the fragment list.
  //  In particular, we are (supposedly) guaranteeing that there is no contained fragment before it
  //  (like would happen if we left the fragment at the end of a unitig).

  if (frIdx != 0)
    frUnitig->reverseComplement();
  frIdx = Unitig::pathPosition(join->frFragID);

  if (toIdx != 0)
    toUnitig->reverseComplement();
  toIdx = Unitig::pathPosition(join->toFragID);

  //  Now if either Idx is not zero, we must have already joined something onto this unitig.

  if ((frIdx != 0) ||
      (toIdx != 0))
    return(false);

  //  Reverse complement (possibly again) the toUnitig, so we can append new fragments.  Not
  //  strictly needed (we'll usually just end up sorting it again) but nice.

  toUnitig->reverseComplement();
  toIdx = Unitig::pathPosition(join->toFragID);

  //  Now, there had better be an edge between toIdx and frIdx.  If not, we messed up somewhere.  We
  //  can still just abort the merge though.

  ufNode          *frFragment = &frUnitig->ufpath[frIdx];
  bool             frForward = (frFragment->position.bgn < frFragment->position.end);
  BestEdgeOverlap *frBest    = OG->getBestEdgeOverlap(frFragment->ident, (frForward == false));

  ufNode          *toFragment = &toUnitig->ufpath[toIdx];
  bool             toForward = (toFragment->position.bgn < toFragment->position.end);
  BestEdgeOverlap *toBest    = OG->getBestEdgeOverlap(toFragment->ident, (toForward == true));

  //  We forced 'fr' to be the first fragment, and 'to' to be the last fragment.  If 'to' is
  //  forward, then its 3' end is sticking out of the 'to' unitig, and the edge should hit that end.

  if (((frBest->fragId() != toFragment->ident) || (frBest->frag3p() != toForward)) &&
      ((toBest->fragId() != frFragment->ident) || (toBest->frag3p() == frForward))) {
    //  Hmmm.  No best edge between these?
    fprintf(logFile, "to: unitig %d frag %d best: %d/%d'\n", toUnitig->id(), toFragment->ident, toBest->fragId(), toBest->frag3p() ? 3 : 5);
    fprintf(logFile, "fr: unitig %d frag %d best: %d/%d'\n", frUnitig->id(), frFragment->ident, frBest->fragId(), frBest->frag3p() ? 3 : 5);
    return(false);
  }

  //  Compute the offset for appending the frUnitig onto the toUnitig.



  //  We could jump down to overlaps and check all the fragments at the join point.
  //  Instead, we (hope) that we do the same type of check at the end of unitig construction.

  return(true);
}



//  Given a fragment and a unitig, this will return exactly one edge that will
//  associate aFrg to some fragment in tUtg.  If no such edge exists, the function
//  returns false.
//
bool
UnitigGraph::findEdgeToUnitig(ufNode *aFrg, BestEdgeOverlap &a5, BestEdgeOverlap &a3,
                              Unitig *tUtg) {
  return(false);
}


//  Given two fragments that share at least one edge, this will find that edge and construct a new
//  edge to make it mutual.
//
//  For example, if there is a best edge from aFrg 3' to bFrg 5', this will return that edge in a3,
//  and also create the symmetric edge in b5.
//
bool
UnitigGraph::findEdges(ufNode *aFrg, BestEdgeOverlap &a5, BestEdgeOverlap &a3,
                       ufNode *bFrg, BestEdgeOverlap &b5, BestEdgeOverlap &b3) {

  if (OG->isContained(aFrg->ident) ||
      OG->isContained(bFrg->ident))
    return(false);

  //  Grab what edges we have.

  a5 = *OG->getBestEdgeOverlap(aFrg->ident, false);
  a3 = *OG->getBestEdgeOverlap(aFrg->ident, true);
  b5 = *OG->getBestEdgeOverlap(bFrg->ident, false);
  b3 = *OG->getBestEdgeOverlap(bFrg->ident, true);

  //  Erase things that aren't correct

  if (a5.fragId() != bFrg->ident)  a5 = BestEdgeOverlap();
  if (a3.fragId() != bFrg->ident)  a3 = BestEdgeOverlap();
  if (b5.fragId() != aFrg->ident)  b5 = BestEdgeOverlap();
  if (b3.fragId() != aFrg->ident)  b3 = BestEdgeOverlap();

  //  If we have no edges left, there are no edges!

  if ((b5.fragId() != aFrg->ident) && (b3.fragId() != aFrg->ident) &&
      (a5.fragId() != bFrg->ident) && (a3.fragId() != bFrg->ident))
    return(false);

  //  If we found TWO edges for any single fragment....that's madness!  That means the fragment
  //  had best dovetail overlaps to the same other fragment off of both ends.  We'll complain
  //  and return failure.  Ideally, data like this will be cleaned up by OBT, or filtered from
  //  our input.
  //
  if (a5.fragId() == a3.fragId()) {
    fprintf(logFile, "findEdges()-- frag %d has multiple edges to frag %d - a5 %d/%d' a3 %d/%d'\n",
            aFrg->ident, a5.fragId(),
            a5.fragId(), a5.frag3p() ? 3 : 5,
            a5.fragId(), a5.frag3p() ? 3 : 5);
  }

  if (b5.fragId() == b3.fragId()) {
    fprintf(logFile, "findEdges()-- frag %d has multiple edges to frag %d - b5 %d/%d' b3 %d/%d'\n",
            bFrg->ident, b5.fragId(),
            b5.fragId(), b5.frag3p() ? 3 : 5,
            b5.fragId(), b5.frag3p() ? 3 : 5);
  }

  if (((a5.fragId() != 0) && (a5.fragId() == a3.fragId())) ||
      ((b5.fragId() != 0) && (b5.fragId() == b3.fragId()))) {
    a5 = BestEdgeOverlap();
    a3 = BestEdgeOverlap();
    b5 = BestEdgeOverlap();
    b3 = BestEdgeOverlap();
    return(false);
  }

  //  Now, populate the other edges using whatever we have.  Best case is that we have two edges
  //  (because we're done).

  assert(((a5.fragId() == bFrg->ident) +
          (a3.fragId() == bFrg->ident) +
          (b5.fragId() == aFrg->ident) +
          (b3.fragId() == aFrg->ident)) <= 2);

  if (((a5.fragId() == bFrg->ident) || (a3.fragId() == bFrg->ident)) &&
      ((b5.fragId() == aFrg->ident) || (b3.fragId() == aFrg->ident)))
    return(true);

  //  Otherwise, we have exactly one edge, and the other one needs to be created.

  assert(((a5.fragId() == bFrg->ident) +
          (a3.fragId() == bFrg->ident) +
          (b5.fragId() == aFrg->ident) +
          (b3.fragId() == aFrg->ident)) == 1);

  if        (a5.fragId() == bFrg->ident) {
    //assert(a5.fragId() == 0);
    assert(a3.fragId() == 0);
    assert(b5.fragId() == 0);
    assert(b3.fragId() == 0);

    //  Edge off of A's 5' end ('false' below)...
    //  ...to B's 3' end (so ANTI or NORMAL -- negate the hangs)
    //  ...to B's 5' end (so INNIE or OUTTIE -- swap the hangs)
    if (a5.frag3p())
      b3.set(aFrg->ident, false, -a5.ahang(), -a5.bhang());
    else
      b5.set(aFrg->ident, false, a5.bhang(), a5.ahang());

  } else if (a3.fragId() == bFrg->ident) {
    assert(a5.fragId() == 0);
    //assert(a3.fragId() == 0);
    assert(b5.fragId() == 0);
    assert(b3.fragId() == 0);

    //  Edge off of A's 3' end ('true' below)...
    //  ...to B's 3' end (so INNIE or OUTTIE -- swap the hangs)
    //  ...to B's 5' end (so ANTI or NORMAL -- negate the hangs)
    if (a3.frag3p())
      b3.set(aFrg->ident, true, a3.bhang(), a3.ahang());
    else
      b5.set(aFrg->ident, true, -a3.ahang(), -a3.bhang());

  } else if (b5.fragId() == aFrg->ident) {
    assert(a5.fragId() == 0);
    assert(a3.fragId() == 0);
    //assert(b5.fragId() == 0);
    assert(b3.fragId() == 0);

    if (b5.frag3p())
      a3.set(bFrg->ident, false, -b5.ahang(), -b5.bhang());
    else
      a5.set(bFrg->ident, false, b5.bhang(), b5.ahang());


  } else if (b3.fragId() == aFrg->ident) {
    assert(a5.fragId() == 0);
    assert(a3.fragId() == 0);
    assert(b5.fragId() == 0);
    //assert(b3.fragId() == 0);

    if (b3.frag3p())
      a3.set(bFrg->ident, true, b3.bhang(), b3.ahang());
    else
      a5.set(bFrg->ident, true, -b3.ahang(), -b3.bhang());

  } else {
    fprintf(stderr, "UnitigGraph::findEdges()-- Logically impossible!\n");
    assert(0);
  }

  //  And now we should have exactly two edges.

  assert(((a5.fragId() == bFrg->ident) +
          (a3.fragId() == bFrg->ident) +
          (b5.fragId() == aFrg->ident) +
          (b3.fragId() == aFrg->ident)) == 2);

  return(true);
}



//  Join two unitigs by using edges to move fragments from the second unitig to the first.
//  Assumes that the second unitig is appended to the first.
//
//  This method DOES NOT WORK due to ties, imprecise lengths of alignments, sorting and reverse
//  complementing.
//
void
UnitigGraph::joinUnitigs_merge(joinEntry *join) {
  uint32    frUnitigID = Unitig::fragIn(join->frFragID);
  uint32    toUnitigID = Unitig::fragIn(join->toFragID);

  Unitig   *frUnitig = unitigs[frUnitigID];
  Unitig   *toUnitig = unitigs[toUnitigID];

  uint32    frIdx = Unitig::pathPosition(join->frFragID);
  uint32    toIdx = Unitig::pathPosition(join->toFragID);

  if (logFileFlagSet(LOG_INTERSECTION_JOINING))
    fprintf(logFile, "Join from unitig "F_U32" (length %d idx %d/%d) to unitig "F_U32" (length %d idx %d/%d) for a total length of %d\n",
            frUnitigID, join->frLen, frIdx, (int)frUnitig->ufpath.size(),
            toUnitigID, join->toLen, toIdx, (int)toUnitig->ufpath.size(),
            join->ttLen);

  //  Over all fragments in the frUnitig, add them to the toUnitig.

  while (frIdx < frUnitig->ufpath.size()) {
    ufNode  *frFragment = &frUnitig->ufpath[frIdx];
    ufNode  *toFragment = &toUnitig->ufpath[toIdx];

    if (OG->isContained(frFragment->ident)) {
      Unitig::removeFrag(frFragment->ident);
      frIdx++;
      continue;
    }

    BestEdgeOverlap  fr5, fr3;
    BestEdgeOverlap  to5, to3;

  tryAgain:
    findEdges(frFragment, fr5, fr3,
              toFragment, to5, to3);

    //  If we found NO edges, find the edge of to frFragment that will place it in the toUnitig.
    //
    if ((fr5.fragId() == 0) &&
        (fr3.fragId() == 0)) {
      fr5 = *OG->getBestEdgeOverlap(frFragment->ident, false);
      fr3 = *OG->getBestEdgeOverlap(frFragment->ident, true);

      fprintf(logFile, "Get real best edges: fr5 %d/%d' in unitig %d -- fr3 %d/%d' in unitig %d\n",
              fr5.fragId(), fr5.frag3p() ? 3 : 5, Unitig::fragIn(fr5.fragId()),
              fr3.fragId(), fr3.frag3p() ? 3 : 5, Unitig::fragIn(fr3.fragId()));

      if (Unitig::fragIn(fr5.fragId()) != toUnitig->id())
        fr5 = BestEdgeOverlap();

      if (Unitig::fragIn(fr3.fragId()) != toUnitig->id())
        fr3 = BestEdgeOverlap();
    }

    //  And if we STILL found no edge.....do what?
    if ((fr5.fragId() == 0) &&
        (fr3.fragId() == 0) &&
        (toIdx > 0)) {
      toFragment = &toUnitig->ufpath[--toIdx];
      goto tryAgain;
    }

    //  Both edges from a single fragment cannot be set -- both ends of a fragment overlap with
    //  the same fragment?!  Makes no sense at all.
    //
    assert((fr5.fragId() == 0) || (fr3.fragId() == 0));
    assert((to5.fragId() == 0) || (to3.fragId() == 0));

    if (logFileFlagSet(LOG_INTERSECTION_JOINING_DEBUG)) {
      fprintf(logFile, "Adding unitig %d frag %d to unitig %d using frag %d at %d,%d\n",
              frUnitig->id(), frFragment->ident,
              toUnitig->id(), toFragment->ident,
              toFragment->position.bgn, toFragment->position.end);
      fprintf(logFile, "fr5 = %8d %1d %5d %5d\n", fr5.fragId(), fr5.frag3p(), fr5.ahang(), fr5.bhang());
      fprintf(logFile, "fr3 = %8d %1d %5d %5d\n", fr3.fragId(), fr3.frag3p(), fr3.ahang(), fr3.bhang());
      fprintf(logFile, "to5 = %8d %1d %5d %5d\n", to5.fragId(), to5.frag3p(), to5.ahang(), to5.bhang());
      fprintf(logFile, "to3 = %8d %1d %5d %5d\n", to3.fragId(), to3.frag3p(), to3.ahang(), to3.bhang());
    }

    if ((fr5.fragId() == 0) &&
        (fr3.fragId() == 0)) {
      //  No suitable edge found to add frFragment to the toUnitig!
      fprintf(logFile, "Adding unitig %d frag %d to unitig %d using frag %d at %d,%d\n",
              frUnitig->id(), frFragment->ident,
              toUnitig->id(), toFragment->ident,
              toFragment->position.bgn, toFragment->position.end);
      fprintf(logFile, "fr5 = %d %d %d %d\n", fr5.fragId(), fr5.frag3p(), fr5.ahang(), fr5.bhang());
      fprintf(logFile, "fr3 = %d %d %d %d\n", fr3.fragId(), fr3.frag3p(), fr3.ahang(), fr3.bhang());
      fprintf(logFile, "to5 = %d %d %d %d\n", to5.fragId(), to5.frag3p(), to5.ahang(), to5.bhang());
      fprintf(logFile, "to3 = %d %d %d %d\n", to3.fragId(), to3.frag3p(), to3.ahang(), to3.bhang());
      fprintf(logFile, "ERROR:  No edge found!  Can't place fragment.\n");
      assert(0);
    }

    //  Now, just add the fragment.  Simple!

    if (toUnitig->addAndPlaceFrag(frFragment->ident, &fr5, &fr3, 
                                  logFileFlagSet(LOG_INTERSECTION_JOINING)) == false) {
      fprintf(logFile, "ERROR:  Failed to place frag %d into extended unitig.\n", frFragment->ident);
      assert(0);
    }

    //  Move to the next fragment in the frUnitig, and to the last fragment (the one we just added)
    //  in the toUnitig.
    frIdx++;
    toIdx = toUnitig->ufpath.size() - 1;
  }

  //  Finally, delete the frUnitig.  It's now in the toUnitig.

  delete frUnitig;
  unitigs[frUnitigID] = NULL;

  //  And make sure the toUnitig is consistent.

  toUnitig->sort();
}





void
UnitigGraph::joinUnitigs_append(joinEntry *join) {
  uint32    frUnitigID = Unitig::fragIn(join->frFragID);
  uint32    toUnitigID = Unitig::fragIn(join->toFragID);

  Unitig   *frUnitig = unitigs[frUnitigID];
  Unitig   *toUnitig = unitigs[toUnitigID];

  uint32    frIdx = Unitig::pathPosition(join->frFragID);
  uint32    toIdx = Unitig::pathPosition(join->toFragID);

  if (logFileFlagSet(LOG_INTERSECTION_JOINING))
    fprintf(logFile, "Join from unitig "F_U32" (length %d idx %d/%d) to unitig "F_U32" (length %d idx %d/%d) for a total length of %d\n",
            frUnitigID, join->frLen, frIdx, (int)frUnitig->ufpath.size(),
            toUnitigID, join->toLen, toIdx, (int)toUnitig->ufpath.size(),
            join->ttLen);

  //  Compute the offset for our append.  We just need to compute where the join fragment would
  //  appear in the unitig.  The join fragment MUST be the first thing in the frUnitig.

  assert(frIdx == 0);

  ufNode           frF,      toF;
  BestEdgeOverlap  fr5, fr3, to5, to3;
  int32            id5,      id3;

  frF.ident = join->frFragID;
  toF.ident = join->toFragID;

  findEdges(&frF, fr5, fr3,
            &toF, to5, to3);

  toUnitig->placeFrag(frF, id5, &fr5,
                      frF, id3, &fr3);

  //  We should get exactly one placement for this.

  assert((id5 == -1) || (id3 == -1));
  assert((id5 != -1) || (id3 != -1));

  int32 offset = MIN(frF.position.bgn, frF.position.end);

  assert(offset > 0);

  //  Over all fragments in the frUnitig, add them to the toUnitig.

  for (; frIdx < frUnitig->ufpath.size(); frIdx++) {
    ufNode  frFragment = frUnitig->ufpath[frIdx];

    toUnitig->addFrag(frFragment, offset, logFileFlagSet(LOG_INTERSECTION_JOINING));
  }

  //  Finally, delete the frUnitig.  It's now in the toUnitig.

  delete frUnitig;
  unitigs[frUnitigID] = NULL;

  //  And make sure the toUnitig is consistent.

  toUnitig->sort();
}





void
UnitigGraph::joinUnitigs(bool enableJoining) {

  if (enableJoining == false)
    return;

  fprintf(logFile, "==> JOINING SPLIT UNITIGS\n");

  //  Sort unitigs by joined size.  Sort.  Join the largest first.

  vector<joinEntry>  joins;

  //  Over all unitigs, evaluate if a unitig is a candidate for merging onto something.

  for (uint32 frID=0; frID<unitigs.size(); frID++) {
    Unitig        *fr   = unitigs[frID];

    if (fr == NULL)
      //  Ain't no unitig here, mister!
      continue;

    if (fr->ufpath.size() == 1)
      //  Ain't no real unitig here, mister!
      continue;

    joinUnitigs_examineFirst(fr, joins);
    joinUnitigs_examineLast(fr, joins);
  }  //  Over all unitigs.


  fprintf(logFile, "Found %d pairs of unitigs to join.\n", (int)joins.size());
  std::sort(joins.begin(), joins.end(), greater<joinEntry>());


  for (uint32 j=0; j<joins.size(); j++) {
    joinEntry  *join = &joins[j];

    if (joinUnitigs_confirmJoin(join) == false)
      //  Join didn't pass muster.  Maybe one unitig already joined to something, maybe
      //  overlaps didn't agree.
      continue;

    joinUnitigs_append(join);
    //joinUnitigs_merge(join);
  }
}



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

static const char *rcsid = "$Id: AS_BOG_Joining.cc,v 1.4 2010-09-30 05:50:17 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"



class joinEntry {
public:
  int32             toFragID;
  int32             frFragID;

  uint32            toLen;
  uint32            frLen;
  uint32            combinedLength;

  bool operator<(const joinEntry &that) const { return(combinedLength < that.combinedLength); };
  bool operator>(const joinEntry &that) const { return(combinedLength > that.combinedLength); };
};


void
UnitigGraph::joinUnitigs(bool enableJoining) {

  fprintf(logFile, "==> JOINING SPLIT UNITIGS\n");

  //  Sort unitigs by joined size.  Sort.  Join the largest first.

  joinEntry  *joins    = new joinEntry [ 2*unitigs.size() ];
  uint32      joinsLen = 0;

  for (uint32 frID=0; frID<unitigs.size(); frID++) {
    Unitig        *fr   = unitigs[frID];

    if (fr == NULL)
      continue;

    if (fr->dovetail_path_ptr->size() == 1)
      continue;

    //  Examine the first fragment.
    {
      DoveTailNode *frg      = &(*fr->dovetail_path_ptr)[0];
      DoveTailNode *tgt      = NULL;
      uint32        tgtEnd   = 0;

      bool  fragForward = (frg->position.bgn < frg->position.end);

      uint32           bestEnd  = (fragForward) ? FIVE_PRIME : THREE_PRIME;
      BestEdgeOverlap *bestEdge = OG->getBestEdgeOverlap(frg->ident, bestEnd);

      int32            toID = 0;
      Unitig          *to   = NULL;

      //  No best edge?  Skip it.
      if (bestEdge->frag_b_id == 0)
        goto skipFirst;

      toID = fr->fragIn(bestEdge->frag_b_id);
      to   = unitigs[toID];

      //  Joining to something teeny?  Skip it.
      if (to->dovetail_path_ptr->size() == 1)
        goto skipFirst;

      //  Figure out which fragment we have an edge to, and which end of it is sticking out from the
      //  end of the unitig.

      if (bestEdge->frag_b_id == (*to->dovetail_path_ptr)[0].ident) {
        tgt      = &(*to->dovetail_path_ptr)[0];
        tgtEnd   = (tgt->position.bgn < tgt->position.end) ? FIVE_PRIME : THREE_PRIME;
      }

      if (bestEdge->frag_b_id == (*to->dovetail_path_ptr)[to->dovetail_path_ptr->size()-1].ident) {
        tgt      = &(*to->dovetail_path_ptr)[to->dovetail_path_ptr->size()-1];
        tgtEnd   = (tgt->position.bgn < tgt->position.end) ? THREE_PRIME : FIVE_PRIME;
      }

      //  Attempt to join to a fragment that isn't first or last.  Skip it.
      if (tgt == NULL)
        goto skipFirst;

      //  Joining to the wrong end of the first/last fragment?  Skip it.
      if (bestEdge->bend != tgtEnd)
        goto skipFirst;

      //  OK, join.

      joins[joinsLen].toFragID       = tgt->ident;
      joins[joinsLen].frFragID       = frg->ident;

      joins[joinsLen].toLen          = to->getLength();
      joins[joinsLen].frLen          = fr->getLength();
      joins[joinsLen].combinedLength = to->getLength() + fr->getLength();

      joinsLen++;

    skipFirst:
      ;
    }

    //  Examine the last fragment.
    {
      DoveTailNode *frg      = &(*fr->dovetail_path_ptr)[0];
      DoveTailNode *tgt      = NULL;
      uint32        tgtEnd   = 0;

      bool  fragForward = (frg->position.bgn < frg->position.end);

      uint32           bestEnd  = (fragForward) ? THREE_PRIME : FIVE_PRIME;
      BestEdgeOverlap *bestEdge = OG->getBestEdgeOverlap(frg->ident, bestEnd);

      int32            toID = 0;
      Unitig          *to   = NULL;

      //  No best edge?  Skip it.
      if (bestEdge->frag_b_id == 0)
        goto skipLast;

      toID = fr->fragIn(bestEdge->frag_b_id);
      to   = unitigs[toID];

      //  Joining to something teeny?  Skip it.
      if (to->dovetail_path_ptr->size() == 1)
        goto skipLast;

      //  Figure out which fragment we have an edge to, and which end of it is sticking out from the
      //  end of the unitig.

      if (bestEdge->frag_b_id == (*to->dovetail_path_ptr)[0].ident) {
        tgt      = &(*to->dovetail_path_ptr)[0];
        tgtEnd   = (tgt->position.bgn < tgt->position.end) ? FIVE_PRIME : THREE_PRIME;
      }

      if (bestEdge->frag_b_id == (*to->dovetail_path_ptr)[to->dovetail_path_ptr->size()-1].ident) {
        tgt      = &(*to->dovetail_path_ptr)[to->dovetail_path_ptr->size()-1];
        tgtEnd   = (tgt->position.bgn < tgt->position.end) ? THREE_PRIME : FIVE_PRIME;
      }

      //  Attempt to join to a fragment that isn't first or last.  Skip it.
      if (tgt == NULL)
        goto skipLast;

      //  Joining to the wrong end of the first/last fragment?  Skip it.
      if (bestEdge->bend != tgtEnd)
        goto skipLast;

      //  OK, join.

      joins[joinsLen].toFragID       = tgt->ident;
      joins[joinsLen].frFragID       = frg->ident;

      joins[joinsLen].toLen          = to->getLength();
      joins[joinsLen].frLen          = fr->getLength();
      joins[joinsLen].combinedLength = to->getLength() + fr->getLength();

      joinsLen++;

    skipLast:
      ;
    }
  }  //  Over all unitigs.

  fprintf(logFile, "Found %d pairs of unitigs to join.\n", joinsLen);

  std::sort(joins, joins + joinsLen, greater<joinEntry>());

  for (uint32 j=0; j<joinsLen; j++) {
    joinEntry  *join = joins + j;

    int32    frUnitigID = Unitig::fragIn(join->frFragID);
    int32    toUnitigID = Unitig::fragIn(join->toFragID);

    Unitig  *frUnitig = unitigs[frUnitigID];
    Unitig  *toUnitig = unitigs[toUnitigID];

    if ((frUnitig == NULL) || (toUnitig == NULL))
      //  Already added one of the unitigs to another one.  Should never happen, really.
      continue;

    int32  frIdx = Unitig::pathPosition(join->frFragID);
    int32  toIdx = Unitig::pathPosition(join->toFragID);

    //  If either fragment is not the first or last, bail.  We already joined something to
    //  one of these unitigs.

    if ((frIdx != 0) && (frIdx != frUnitig->dovetail_path_ptr->size() - 1))
      continue;

    if ((toIdx != 0) && (toIdx != toUnitig->dovetail_path_ptr->size() - 1))
      continue;

    if (logFileFlagSet(LOG_INTERSECTION_JOINING))
      fprintf(logFile, "Join from unitig %d (length %d idx %d/%d) to unitig %d (length %d idx %d/%d) for a total length of %d\n",
              frUnitigID, join->frLen, frIdx, frUnitig->dovetail_path_ptr->size(),
              toUnitigID, join->toLen, toIdx, toUnitig->dovetail_path_ptr->size(),
              join->combinedLength);

    //  Reverse complement to ensure that we append to the tail of the toUnitig, and grab fragments
    //  from the start of the frUnitig.

    if (frIdx != 0)
      frUnitig->reverseComplement();

    if (toIdx == 0)
      toUnitig->reverseComplement();

    frIdx = 0;
    toIdx = toUnitig->dovetail_path_ptr->size() - 1;

    //  Over all fragments in the frUnitig, add them to the toUnitig.

    while (frIdx < frUnitig->dovetail_path_ptr->size()) {
      DoveTailNode  *frFragment = &(*frUnitig->dovetail_path_ptr)[frIdx];
      DoveTailNode  *toFragment = &(*toUnitig->dovetail_path_ptr)[toIdx];

      //  Construct an edge that will place the frFragment onto the toUnitig.  This edge
      //  can come from either fragment.

      BestEdgeOverlap  *best5fr = OG->getBestEdgeOverlap(frFragment->ident, FIVE_PRIME);
      BestEdgeOverlap  *best3fr = OG->getBestEdgeOverlap(frFragment->ident, THREE_PRIME);
      BestEdgeOverlap  *best5to = OG->getBestEdgeOverlap(toFragment->ident, FIVE_PRIME);
      BestEdgeOverlap  *best3to = OG->getBestEdgeOverlap(toFragment->ident, THREE_PRIME);

      BestEdgeOverlap   joinEdge5;
      BestEdgeOverlap   joinEdge3;

      joinEdge5.frag_b_id = 0;
      joinEdge5.bend      = 0;
      joinEdge5.ahang     = 0;
      joinEdge5.bhang     = 0;

      joinEdge3.frag_b_id = 0;
      joinEdge3.bend      = 0;
      joinEdge3.ahang     = 0;
      joinEdge3.bhang     = 0;

      //  The first two cases are trivial, the edge is already correctly oriented, so just copy it over.

      if ((best5fr->frag_b_id == toFragment->ident) ||
          (Unitig::fragIn(best5fr->frag_b_id) == toUnitigID))
        joinEdge5 = *best5fr;

      if ((best3fr->frag_b_id == toFragment->ident) ||
          (Unitig::fragIn(best3fr->frag_b_id) == toUnitigID))
        joinEdge3 = *best3fr;


      //  The last two cases are trouble.  The edge is from the toFrag to the frFrag, and we need to
      //  reverse it.  In the first case, the edge is off of the 5' end of the toFrag, but which
      //  edge we create depends on what end of the frFrag it hits.
      //
      //  If the fragments are innie/outtie, we need to reverse the overlap to maintain that the
      //  A fragment is forward.
      //
      if (best5to->frag_b_id == frFragment->ident) {
        if (best5to->bend == FIVE_PRIME) {
          joinEdge5.frag_b_id = toFragment->ident;
          joinEdge5.bend      = FIVE_PRIME;
          joinEdge5.ahang     = best5to->bhang;
          joinEdge5.bhang     = best5to->ahang;
        } else {
          joinEdge3.frag_b_id = toFragment->ident;
          joinEdge3.bend      = FIVE_PRIME;
          joinEdge3.ahang     = -best5to->ahang;
          joinEdge3.bhang     = -best5to->bhang;
        }
      }

      if (best3to->frag_b_id == frFragment->ident) {
        if (best3to->bend == FIVE_PRIME) {
          joinEdge5.frag_b_id = toFragment->ident;
          joinEdge5.bend      = THREE_PRIME;
          joinEdge5.ahang     = -best3to->ahang;
          joinEdge5.bhang     = -best3to->bhang;
        } else {
          joinEdge3.frag_b_id = toFragment->ident;
          joinEdge3.bend      = THREE_PRIME;
          joinEdge3.ahang     = best3to->bhang;
          joinEdge3.bhang     = best3to->ahang;
        }
      }

      if (logFileFlagSet(LOG_INTERSECTION_JOINING_DEBUG)) {
        fprintf(logFile, "Adding frag %d using frag %d (%d,%d)\n", frFragment->ident, toFragment->ident, toFragment->position.bgn, toFragment->position.end);
        fprintf(logFile, "best5fr = %8d %1d %5d %5d\n", best5fr->frag_b_id, best5fr->bend, best5fr->ahang, best5fr->bhang);
        fprintf(logFile, "best3fr = %8d %1d %5d %5d\n", best3fr->frag_b_id, best3fr->bend, best3fr->ahang, best3fr->bhang);
        fprintf(logFile, "best5to = %8d %1d %5d %5d\n", best5to->frag_b_id, best5to->bend, best5to->ahang, best5to->bhang);
        fprintf(logFile, "best3to = %8d %1d %5d %5d\n", best3to->frag_b_id, best3to->bend, best3to->ahang, best3to->bhang);
        fprintf(logFile, "join3   = %8d %1d %5d %5d\n", best5to->frag_b_id, best5to->bend, best5to->ahang, best5to->bhang);
        fprintf(logFile, "join5   = %8d %1d %5d %5d\n", best3to->frag_b_id, best3to->bend, best3to->ahang, best3to->bhang);
      }

      if ((joinEdge5.frag_b_id == 0) &&
          (joinEdge3.frag_b_id == 0)) {
        //  No suitable edge found to add frFragment to the toUnitig!
        fprintf(logFile, "ERROR:  No edge found!  Can't place fragment.\n");
        fprintf(logFile, "best5fr = %d %d %d %d\n", best5fr->frag_b_id, best5fr->bend, best5fr->ahang, best5fr->bhang);
        fprintf(logFile, "best3fr = %d %d %d %d\n", best3fr->frag_b_id, best3fr->bend, best3fr->ahang, best3fr->bhang);
        fprintf(logFile, "best5to = %d %d %d %d\n", best5to->frag_b_id, best5to->bend, best5to->ahang, best5to->bhang);
        fprintf(logFile, "best3to = %d %d %d %d\n", best3to->frag_b_id, best3to->bend, best3to->ahang, best3to->bhang);
        assert(0);
      }

      //  Now, just add the fragment.  Simple!

      if (toUnitig->addAndPlaceFrag(frFragment->ident,
                                    (joinEdge5.frag_b_id == 0) ? NULL : &joinEdge5,
                                    (joinEdge3.frag_b_id == 0) ? NULL : &joinEdge3,
                                    logFileFlagSet(LOG_INTERSECTION_JOINING)) == false) {
        fprintf(logFile, "ERROR:  Failed to place frag %d into extended unitig.\n", frFragment->ident);
        assert(0);
      }

      frIdx++;
      toIdx++;
    }

    //  Finally, delete the frUnitig.  It's now in the toUnitig.

    delete frUnitig;
    unitigs[frUnitigID] = NULL;
  }  //  Over all joins.
}



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

static const char *rcsid = "$Id: AS_BOG_Joining.cc,v 1.7 2011-09-06 02:15:18 mkotelbajcvi Exp $";

#include "AS_BOG_Datatypes.H"
#include "AS_BOG_UnitigGraph.H"
#include "AS_BOG_BestOverlapGraph.H"

#include "MultiAlignStore.H"



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

    if (fr->ufpath.size() == 1)
      continue;

    //  Examine the first fragment.
    {
      ufNode *frg      = &fr->ufpath[0];
      ufNode *tgt      = NULL;
      bool    tgt3p    = false;

      bool             fragForward = (frg->position.bgn < frg->position.end);
      BestEdgeOverlap *bestEdge = OG->getBestEdgeOverlap(frg->ident, (fragForward == false));

      int32            toID = 0;
      Unitig          *to   = NULL;

      //  No best edge?  Skip it.
      if (bestEdge->fragId() == 0)
        goto skipFirst;

      toID = fr->fragIn(bestEdge->fragId());
      to   = unitigs[toID];

      //  Joining to something teeny?  Skip it.
      if (to->ufpath.size() == 1)
        goto skipFirst;

      //  Figure out which fragment we have an edge to, and which end of it is sticking out from the
      //  end of the unitig.

      if (bestEdge->fragId() == to->ufpath[0].ident) {
        tgt      = &to->ufpath[0];
        tgt3p    = (tgt->position.bgn > tgt->position.end);
      }

      if (bestEdge->fragId() == to->ufpath[to->ufpath.size()-1].ident) {
        tgt      = &to->ufpath[to->ufpath.size()-1];
        tgt3p    = (tgt->position.bgn < tgt->position.end);
      }

      //  Attempt to join to a fragment that isn't first or last.  Skip it.
      if (tgt == NULL)
        goto skipFirst;

      //  Joining to the wrong end of the first/last fragment?  Skip it.
      if (bestEdge->frag3p() != tgt3p)
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
      ufNode *frg      = &fr->ufpath[0];
      ufNode *tgt      = NULL;
      uint32  tgt3p    = false;

      bool             fragForward = (frg->position.bgn < frg->position.end);
      BestEdgeOverlap *bestEdge = OG->getBestEdgeOverlap(frg->ident, fragForward);

      int32            toID = 0;
      Unitig          *to   = NULL;

      //  No best edge?  Skip it.
      if (bestEdge->fragId() == 0)
        goto skipLast;

      toID = fr->fragIn(bestEdge->fragId());
      to   = unitigs[toID];

      //  Joining to something teeny?  Skip it.
      if (to->ufpath.size() == 1)
        goto skipLast;

      //  Figure out which fragment we have an edge to, and which end of it is sticking out from the
      //  end of the unitig.

      if (bestEdge->fragId() == to->ufpath[0].ident) {
        tgt      = &to->ufpath[0];
        tgt3p    = (tgt->position.bgn > tgt->position.end);
      }

      if (bestEdge->fragId() == to->ufpath[to->ufpath.size()-1].ident) {
        tgt      = &to->ufpath[to->ufpath.size()-1];
        tgt3p    = (tgt->position.bgn < tgt->position.end);
      }

      //  Attempt to join to a fragment that isn't first or last.  Skip it.
      if (tgt == NULL)
        goto skipLast;

      //  Joining to the wrong end of the first/last fragment?  Skip it.
      if (bestEdge->frag3p() != tgt3p)
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

    uint32    frUnitigID = Unitig::fragIn(join->frFragID);
    uint32    toUnitigID = Unitig::fragIn(join->toFragID);

    Unitig  *frUnitig = unitigs[frUnitigID];
    Unitig  *toUnitig = unitigs[toUnitigID];

    if ((frUnitig == NULL) || (toUnitig == NULL))
      //  Already added one of the unitigs to another one.  Should never happen, really.
      continue;

    uint32  frIdx = Unitig::pathPosition(join->frFragID);
    uint32  toIdx = Unitig::pathPosition(join->toFragID);

    //  If either fragment is not the first or last, bail.  We already joined something to
    //  one of these unitigs.

    if ((frIdx != 0) && (frIdx != frUnitig->ufpath.size() - 1))
      continue;

    if ((toIdx != 0) && (toIdx != toUnitig->ufpath.size() - 1))
      continue;

    if (logFileFlagSet(LOG_INTERSECTION_JOINING))
      fprintf(logFile, "Join from unitig "F_U32" (length %d idx %d/"F_SIZE_T") to unitig "F_U32" (length %d idx %d/"F_SIZE_T") for a total length of "F_U32"\n",
              frUnitigID, join->frLen, frIdx, frUnitig->ufpath.size(),
              toUnitigID, join->toLen, toIdx, toUnitig->ufpath.size(),
              join->combinedLength);

    //  Reverse complement to ensure that we append to the tail of the toUnitig, and grab fragments
    //  from the start of the frUnitig.

    if (frIdx != 0)
      frUnitig->reverseComplement();

    if (toIdx == 0)
      toUnitig->reverseComplement();

    frIdx = 0;
    toIdx = toUnitig->ufpath.size() - 1;

    //  Over all fragments in the frUnitig, add them to the toUnitig.

    while (frIdx < frUnitig->ufpath.size()) {
      ufNode  *frFragment = &frUnitig->ufpath[frIdx];
      ufNode  *toFragment = &toUnitig->ufpath[toIdx];

      //  Construct an edge that will place the frFragment onto the toUnitig.  This edge
      //  can come from either fragment.

      BestEdgeOverlap  *best5fr = OG->getBestEdgeOverlap(frFragment->ident, false);
      BestEdgeOverlap  *best3fr = OG->getBestEdgeOverlap(frFragment->ident, true);
      BestEdgeOverlap  *best5to = OG->getBestEdgeOverlap(toFragment->ident, false);
      BestEdgeOverlap  *best3to = OG->getBestEdgeOverlap(toFragment->ident, true);

      BestEdgeOverlap   joinEdge5;
      BestEdgeOverlap   joinEdge3;

      //  The first two cases are trivial, the edge is already correctly oriented, so just copy it over.

      if ((best5fr->fragId() == toFragment->ident) ||
          (Unitig::fragIn(best5fr->fragId()) == toUnitigID))
        joinEdge5 = *best5fr;

      if ((best3fr->fragId() == toFragment->ident) ||
          (Unitig::fragIn(best3fr->fragId()) == toUnitigID))
        joinEdge3 = *best3fr;


      //  The last two cases are trouble.  The edge is from the toFrag to the frFrag, and we need to
      //  reverse it.  In the first case, the edge is off of the 5' end of the toFrag, but which
      //  edge we create depends on what end of the frFrag it hits.
      //
      //  If the fragments are innie/outtie, we need to reverse the overlap to maintain that the
      //  A fragment is forward.
      //
      if (best5to->fragId() == frFragment->ident) {
        if (best5to->frag3p() == false)
          joinEdge5.set(toFragment->ident, false, best5to->bhang(), best5to->ahang());
        else
          joinEdge3.set(toFragment->ident, false, -best5to->ahang(), -best5to->bhang());
      }

      if (best3to->fragId() == frFragment->ident) {
        if (best3to->frag3p() == false)
          joinEdge5.set(toFragment->ident, true, -best3to->ahang(), -best3to->bhang());
        else
          joinEdge3.set(toFragment->ident, true, best3to->bhang(), best3to->ahang());
      }

      if (logFileFlagSet(LOG_INTERSECTION_JOINING_DEBUG)) {
        fprintf(logFile, "Adding frag %d using frag %d (%d,%d)\n", frFragment->ident, toFragment->ident, toFragment->position.bgn, toFragment->position.end);
        fprintf(logFile, "best5fr = %8d %1d %5d %5d\n", best5fr->fragId(), best5fr->frag3p(), best5fr->ahang(), best5fr->bhang());
        fprintf(logFile, "best3fr = %8d %1d %5d %5d\n", best3fr->fragId(), best3fr->frag3p(), best3fr->ahang(), best3fr->bhang());
        fprintf(logFile, "best5to = %8d %1d %5d %5d\n", best5to->fragId(), best5to->frag3p(), best5to->ahang(), best5to->bhang());
        fprintf(logFile, "best3to = %8d %1d %5d %5d\n", best3to->fragId(), best3to->frag3p(), best3to->ahang(), best3to->bhang());
        fprintf(logFile, "join3   = %8d %1d %5d %5d\n", best5to->fragId(), best5to->frag3p(), best5to->ahang(), best5to->bhang());
        fprintf(logFile, "join5   = %8d %1d %5d %5d\n", best3to->fragId(), best3to->frag3p(), best3to->ahang(), best3to->bhang());
      }

      if ((joinEdge5.fragId() == 0) &&
          (joinEdge3.fragId() == 0)) {
        //  No suitable edge found to add frFragment to the toUnitig!
        fprintf(logFile, "ERROR:  No edge found!  Can't place fragment.\n");
        fprintf(logFile, "best5fr = %d %d %d %d\n", best5fr->fragId(), best5fr->frag3p(), best5fr->ahang(), best5fr->bhang());
        fprintf(logFile, "best3fr = %d %d %d %d\n", best3fr->fragId(), best3fr->frag3p(), best3fr->ahang(), best3fr->bhang());
        fprintf(logFile, "best5to = %d %d %d %d\n", best5to->fragId(), best5to->frag3p(), best5to->ahang(), best5to->bhang());
        fprintf(logFile, "best3to = %d %d %d %d\n", best3to->fragId(), best3to->frag3p(), best3to->ahang(), best3to->bhang());
        assert(0);
      }

      //  Now, just add the fragment.  Simple!

      if (toUnitig->addAndPlaceFrag(frFragment->ident,
                                    (joinEdge5.fragId() == 0) ? NULL : &joinEdge5,
                                    (joinEdge3.fragId() == 0) ? NULL : &joinEdge3,
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


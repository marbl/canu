
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

static const char *rcsid = "$Id: AS_BAT_PlaceContains.C,v 1.5 2012-01-05 16:29:26 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_PlaceContains.H"


void
placeContainsUsingBestOverlaps(UnitigVector &unitigs) {
  uint32   fragsPlaced  = 1;
  uint32   fragsPending = 0;

  logFileFlags &= ~LOG_PLACE_FRAG;

  while (fragsPlaced > 0) {
    fragsPlaced  = 0;
    fragsPending = 0;

    fprintf(logFile, "==> PLACING CONTAINED FRAGMENTS\n");

    for (uint32 fid=1; fid<FI->numFragments()+1; fid++) {
      BestContainment *bestcont = OG->getBestContainer(fid);
      Unitig          *utg;

      if (bestcont->isContained == false)
        //  Not a contained fragment.
        continue;

      if (Unitig::fragIn(fid) != 0)
        //  Containee already placed.
        continue;

      if (Unitig::fragIn(bestcont->container) == 0) {
        //  Container not placed (yet).
        fragsPending++;
        continue;
      }

      utg = unitigs[Unitig::fragIn(bestcont->container)];
      utg->addContainedFrag(fid, bestcont, logFileFlagSet(LOG_INITIAL_CONTAINED_PLACEMENT));

      if (utg->id() != Unitig::fragIn(fid))
        fprintf(logFile, "placeContainsUsingBestOverlaps()-- FAILED to add frag %d to unitig %d.\n", fid, bestcont->container);
      assert(utg->id() == Unitig::fragIn(fid));


      fragsPlaced++;
    }

    fprintf(logFile, "==> PLACING CONTAINED FRAGMENTS - placed %d fragments; still need to place %d\n",
            fragsPlaced, fragsPending);

    if ((fragsPlaced == 0) && (fragsPending > 0)) {
      fprintf(logFile, "Stopping contained fragment placement due to zombies.\n");
      fragsPlaced  = 0;
      fragsPending = 0;
    }
  }

  for (uint32 ti=1; ti<unitigs.size(); ti++) {
    Unitig *utg = unitigs[ti];

    if (utg)
      utg->sort();
  }
}



void
placeContainsUsingBestOverlaps(Unitig *target, set<AS_IID> *fragments) {
  uint32   fragsPlaced  = 1;

  logFileFlags &= ~LOG_PLACE_FRAG;

  while (fragsPlaced > 0) {
    fragsPlaced  = 0;

    for (set<AS_IID>::iterator it=fragments->begin(); it != fragments->end(); it++) {
      AS_IID           fid      = *it;
      BestContainment *bestcont = OG->getBestContainer(fid);

      if ((bestcont->isContained == false) ||
          (Unitig::fragIn(fid) != 0) ||
          (Unitig::fragIn(bestcont->container) == 0) ||
          (Unitig::fragIn(bestcont->container) != target->id()))
        //  Not a contained fragment OR
        //  Containee already placed OR
        //  Container not placed (yet) OR
        //  Containee not in the target unitig
        continue;

      target->addContainedFrag(fid, bestcont, false);

      if (target->id() != Unitig::fragIn(fid))
        fprintf(logFile, "placeContainsUsingBestOverlaps()-- FAILED to add frag %d to unitig %d.\n", fid, bestcont->container);
      assert(target->id() == Unitig::fragIn(fid));

      fragsPlaced++;
    }
  }

  target->sort();
}


void
placeContainsUsingAllOverlaps(UnitigVector &unitigs,
                              bool   withMatesToNonContained,
                              bool   withMatesToUnambiguousContain) {

#if 0
  for (uint32 fid=1; fid<FI->numFragments()+1; fid++) {
    ufNode frag;

    if (Unitig::fragIn(fid) > 0)
      //  Fragment placed already.
      continue;

    frag.ident = fid;

    placeFragUsingOverlaps(frag, ovlStoreUniq, ovlStoreRept);
  }
#endif
}

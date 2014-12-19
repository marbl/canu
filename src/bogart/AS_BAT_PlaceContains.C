
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

static const char *rcsid = "$Id$";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_PlaceContains.H"

#include "AS_BAT_PlaceFragUsingOverlaps.H"


void
placeContainsUsingBestOverlaps(UnitigVector &unitigs) {
  uint32   fragsPlaced  = 1;
  uint32   fragsPending = 0;

  uint32  *nReadsPer = new uint32 [unitigs.size()];

  uint32   totalPlaced            = 0;
  uint32   totalPlacedInSingleton = 0;

  logFileFlags &= ~LOG_PLACE_FRAG;

  for (uint32 ii=0; ii<unitigs.size(); ii++)
    nReadsPer[ii] = (unitigs[ii] == NULL) ? 0 : unitigs[ii]->getNumFrags();

  while (fragsPlaced > 0) {
    fragsPlaced  = 0;
    fragsPending = 0;

    writeLog("==> PLACING CONTAINED FRAGMENTS\n");

    for (uint32 fid=1; fid<FI->numFragments()+1; fid++) {
      BestContainment *bestcont = OG->getBestContainer(fid);

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

      uint32  utgid = Unitig::fragIn(bestcont->container);
      Unitig *utg   = unitigs[utgid];

      totalPlaced++;

      if (nReadsPer[utgid] == 1)
        totalPlacedInSingleton++;

      utg->addContainedFrag(fid, bestcont, logFileFlagSet(LOG_INITIAL_CONTAINED_PLACEMENT));

      if (utg->id() != Unitig::fragIn(fid))
        writeLog("placeContainsUsingBestOverlaps()-- FAILED to add frag %d to unitig %d.\n", fid, bestcont->container);
      assert(utg->id() == Unitig::fragIn(fid));


      fragsPlaced++;
    }

    writeLog("placeContainsUsingBestOverlaps()-- Placed %d fragments; still need to place %d\n",
            fragsPlaced, fragsPending);

    if ((fragsPlaced == 0) && (fragsPending > 0))
      writeLog("placeContainsUsingBestOverlaps()-- Stopping contained fragment placement due to zombies.\n");
  }

  writeLog("placeContainsUsingBestOverlaps()-- %u frags placed in unitigs (including singleton unitigs)\n", totalPlaced);
  writeLog("placeContainsUsingBestOverlaps()-- %u frags placed in singleton unitigs\n", totalPlacedInSingleton);
  writeLog("placeContainsUsingBestOverlaps()-- %u frags unplaced\n", fragsPending);

  delete [] nReadsPer;

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
        writeLog("placeContainsUsingBestOverlaps()-- FAILED to add frag %d to unitig %d.\n", fid, bestcont->container);
      assert(target->id() == Unitig::fragIn(fid));

      fragsPlaced++;
    }
  }

  target->sort();
}


void
placeContainsUsingAllOverlaps(UnitigVector &unitigs,
                              bool          withMatesToNonContained,
                              bool          withMatesToUnambiguousContain) {

  //  UNFINISHED.  This results in crashes later in the process.

  for (uint32 fid=1; fid<FI->numFragments()+1; fid++) {
    ufNode frg;
    //ufNode mat;

    if (Unitig::fragIn(fid) > 0)
      //  Fragment placed already.
      continue;

    frg.ident = fid;
    //mat.ident = 0;  //mid;

    overlapPlacement    frgPlacement;
    //overlapPlacement    matPlacement;

    frgPlacement.errors  = 4.0e9;
    frgPlacement.aligned = 1;

    //matPlacement.errors  = 4.0e9;
    //matPlacement.aligned = 1;

    vector<overlapPlacement>   placements;

    //  Place the read.

    placeFragUsingOverlaps(unitigs, NULL, frg.ident, placements);

    //  Search the placements for the highest expect identity placement using all overlaps in the unitig.

    Unitig  *frgTig = NULL;
    Unitig  *matTig = NULL;

    for (uint32 i=0; i<placements.size(); i++) {
      if (placements[i].fCoverage < 0.99)
        continue;

      if (placements[i].errors / placements[i].aligned < frgPlacement.errors / frgPlacement.aligned) {
        frgPlacement = placements[i];
        frgTig       = unitigs[placements[i].tigID];
      }
    }

    frg.ident             = frgPlacement.frgID;
    frg.contained         = 0;
    frg.parent            = 0;
    frg.ahang             = 0;
    frg.bhang             = 0;
    frg.position          = frgPlacement.position;
    frg.containment_depth = 0;

    if ((frg.position.bgn == 0) &&
        (frg.position.end == 0))
      //  Failed to place the contained read anywhere.  We should probably just make a new unitig
      //  for it right here.
      continue;

    //  Place the mate

    //  Add the placed read to the unitig.

    writeLog("placeContainsUsingAllOverlaps()-- frag %u placed in tig %u at %u-%u.\n",
             frg.ident, frgTig->id(), frg.position.bgn, frg.position.end);

    frgTig->addFrag(frg, 0, false);
  }
}

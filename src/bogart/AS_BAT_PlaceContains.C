
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_PlaceContains.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2014-JAN-29
 *      are Copyright 2010-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2015-JUN-03
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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

      utg->addContainedFrag(fid, bestcont, true);  //logFileFlagSet(LOG_INITIAL_CONTAINED_PLACEMENT));

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
placeContainsUsingBestOverlaps(Unitig *target, set<uint32> *fragments) {
  uint32   fragsPlaced  = 1;

  logFileFlags &= ~LOG_PLACE_FRAG;

  while (fragsPlaced > 0) {
    fragsPlaced  = 0;

    for (set<uint32>::iterator it=fragments->begin(); it != fragments->end(); it++) {
      uint32           fid      = *it;
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
                              double        erate) {

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

    placeFragUsingOverlaps(unitigs, erate, NULL, frg.ident, placements);

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

    //  Add the placed read to the unitig.

    writeLog("placeContainsUsingAllOverlaps()-- frag %u placed in tig %u at %u-%u.\n",
             frg.ident, frgTig->id(), frg.position.bgn, frg.position.end);

    frgTig->addFrag(frg, 0, false);
  }
}

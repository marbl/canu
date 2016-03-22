
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
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
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
placeContainsUsingAllOverlaps(UnitigVector &unitigs,
                              double        erate) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  uint32       *placedTig = new uint32      [FI->numFragments() + 1];
  SeqInterval  *placedPos = new SeqInterval [FI->numFragments() + 1];

  memset(placedTig, 0, sizeof(uint32)      * (FI->numFragments() + 1));
  memset(placedPos, 0, sizeof(SeqInterval) * (FI->numFragments() + 1));

  //  Just some logging.

  uint32   nPlaced  = 0;
  uint32   nToPlace = 0;
  uint32   nFailed  = 0;

  for (uint32 fid=1; fid<FI->numFragments()+1; fid++)
    if (Unitig::fragIn(fid) == 0)
      nToPlace++;

  fprintf(stderr, "placeContains()-- placing %u contained reads, with %d threads.\n", nToPlace, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fid=1; fid<FI->numFragments()+1; fid++) {
    if (Unitig::fragIn(fid) > 0)
      continue;

    //  Place the read.

    vector<overlapPlacement>   placements;

    placeFragUsingOverlaps(unitigs, erate, NULL, fid, placements);

    //  Search the placements for the highest expected identity placement using all overlaps in the unitig.

    uint32   b = UINT32_MAX;

    for (uint32 i=0; i<placements.size(); i++) {
#if 0
      writeLog("frag %u placed in tig %u (%u reads) with coverage %f and ident %f\n",
               fid, placements[i].tigID, unitigs[placements[i].tigID]->ufpath.size(),
               placements[i].fCoverage,
               placements[i].errors / placements[i].aligned);
#endif

      if (placements[i].fCoverage < 0.99)                   //  Ignore partially placed reads.
        continue;

      if (unitigs[placements[i].tigID]->ufpath.size() == 1)  //  Ignore placements in singletons.
        continue;

      if ((b == UINT32_MAX) ||
          (placements[i].errors / placements[i].aligned < placements[b].errors / placements[b].aligned))
        b = i;
    }

    //  If we didn't find a best, b will be invalid; set positions for adding to a new tig.
    //  If we did, save both the position it was placed at, and the tigID it was placed in.

    if (b == UINT32_MAX) {
      placedPos[fid].bgn = 0;
      placedPos[fid].end = FI->fragmentLength(fid);
    } else {
      placedTig[fid] = placements[b].tigID;
      placedPos[fid] = placements[b].position;
    }
  }

  //  All reads placed, now just dump them in their correct tigs.

  for (uint32 fid=1; fid<FI->numFragments()+1; fid++) {
    Unitig  *tig = NULL;
    ufNode   frg;

    if (Unitig::fragIn(fid) > 0)
      continue;

    //  If not placed, dump it in a new unitig.  Otherwise, it was placed somewhere, grab the tig.

    if (placedTig[fid] == 0) {
      nFailed++;
      tig = unitigs.newUnitig(false);
    } else {
      nPlaced++;
      tig = unitigs[placedTig[fid]];
    }

    //  Regardless, add it to the tig.

    frg.ident             = fid;
    frg.contained         = 0;
    frg.parent            = 0;
    frg.ahang             = 0;
    frg.bhang             = 0;
    frg.position          = placedPos[fid];
    frg.containment_depth = 0;

    writeLog("placeContainsUsingAllOverlaps()-- frag %u placed in tig %u at %u-%u.\n",
             fid, tig->id(), frg.position.bgn, frg.position.end);

    tig->addFrag(frg, 0, false);
  }

  //  Cleanup.

  delete [] placedPos;
  delete [] placedTig;

  fprintf(stderr, "placeContains()-- placed %u contained reads.  failed to place %u contained reads.\n", nPlaced, nFailed);

  //  But wait!  All the tigs need to be sorted.  Well, not really _all_, but the hard ones to sort
  //  are big, and those quite likely had reads added to them, so it's really not worth the effort
  //  of tracking which ones need sorting, since the ones that don't need it are trivial to sort.

  for (uint32 ti=1; ti<unitigs.size(); ti++) {
    Unitig *utg = unitigs[ti];

    if (utg)
      utg->sort();
  }
}

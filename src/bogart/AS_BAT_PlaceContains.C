
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

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_PlaceContains.H"
#include "AS_BAT_PlaceReadUsingOverlaps.H"

#undef  SHOW_PLACEMENT_DETAIL    //  Reports evidence (too much) for placing reads.
#undef  SHOW_PLACEMENT           //  Reports where the read was placed.


void
breakSingletonTigs(TigVector &tigs) {

  //  For any singleton unitig, eject the read and delete the unitig.  Eventually,
  //  we will stop making singleton tigs.

  uint32   removed = 0;

  for (uint32 ti=1; ti<tigs.size(); ti++) {
    Unitig *utg = tigs[ti];

    if (utg == NULL)
      continue;

    if (utg->ufpath.size() > 1)
      continue;

    if (OG->isZombie(utg->ufpath[0].ident) == true) {
      writeLog("Not breaking sinleton zombie tig %u with read %u\n", ti, utg->ufpath[0].ident);
      continue;
    }

    tigs[ti] = NULL;                           //  Remove the tig from the list
    tigs.registerRead(utg->ufpath[0].ident);   //  Eject the read
    delete utg;                                //  Reclaim space
    removed++;                                 //  Count
  }

  writeStatus("breakSingletonTigs()-- Removed %u singleton tig%s; reads are now unplaced.\n",
              removed, (removed == 1) ? "" : "s");
}



void
placeUnplacedUsingAllOverlaps(TigVector           &tigs,
                              const char   *UNUSED(prefix)) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  uint32       *placedTig = new uint32      [RI->numReads() + 1];
  SeqInterval  *placedPos = new SeqInterval [RI->numReads() + 1];

  memset(placedTig, 0, sizeof(uint32)      * (RI->numReads() + 1));
  memset(placedPos, 0, sizeof(SeqInterval) * (RI->numReads() + 1));

  //  Just some logging.  Count the number of reads we try to place.

  uint32   nToPlaceContained = 0;
  uint32   nToPlace          = 0;
  uint32   nPlacedContained  = 0;
  uint32   nPlaced           = 0;
  uint32   nFailedContained  = 0;
  uint32   nFailed           = 0;

  for (uint32 fid=1; fid<RI->numReads()+1; fid++)
    if (tigs.inUnitig(fid) == 0)    //  I'm NOT ambiguous!
      if (OG->isContained(fid))
        nToPlaceContained++;
      else
        nToPlace++;

  writeStatus("\n");
  writeStatus("placeContains()-- placing %u contained and %u unplaced reads, with %d thread%s.\n",
              nToPlaceContained, nToPlace, numThreads, (numThreads == 1) ? "" : "s");

  //  Do the placing!

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fid=1; fid<RI->numReads()+1; fid++) {
    bool  enableLog = true;

    if (tigs.inUnitig(fid) > 0)
      continue;

    //  Place the read.

    vector<overlapPlacement>   placements;

    placeReadUsingOverlaps(tigs, NULL, fid, placements, placeRead_fullMatch);

    //  If all placements are in singletons, allow them.  If any placement is to a 'real' tig,
    //  ignore singleton placements.

    bool  ignoreSingleton = false;

    for (uint32 i=0; i<placements.size(); i++)
      if ((placements[i].fCoverage >= 0.99) &&
          (tigs[placements[i].tigID]->ufpath.size() > 1))
        ignoreSingleton = true;

    //  Search the placements for the highest expected identity placement using all overlaps in the unitig.

    uint32   b = UINT32_MAX;

    for (uint32 i=0; i<placements.size(); i++) {
      Unitig *tig = tigs[placements[i].tigID];

      if (placements[i].fCoverage < 0.99)     //  Ignore partially placed reads.
        continue;

      if ((ignoreSingleton == true) &&
          (tig->ufpath.size() == 1))          //  Ignore placements in singletons.
        continue;

      uint32  bgn   = placements[i].position.min();
      uint32  end   = placements[i].position.max();

      double  erate = placements[i].errors / placements[i].aligned;

      if (tig->overlapConsistentWithTig(5.0, bgn, end, erate) < 0.5) {
        if ((enableLog == true) && (logFileFlagSet(LOG_PLACE_UNPLACED)))
          writeLog("read %8u tested tig %6u (%6u reads) at %8u-%8u (cov %7.5f erate %6.4f) - HIGH ERROR\n",
                   fid, placements[i].tigID, tig->ufpath.size(), placements[i].position.bgn, placements[i].position.end, placements[i].fCoverage, erate);
        continue;
      }

      if ((enableLog == true) && (logFileFlagSet(LOG_PLACE_UNPLACED)))
        writeLog("read %8u tested tig %6u (%6u reads) at %8u-%8u (cov %7.5f erate %6.4f)\n",
                 fid, placements[i].tigID, tig->ufpath.size(), placements[i].position.bgn, placements[i].position.end, placements[i].fCoverage, erate);

      if ((b == UINT32_MAX) ||
          (placements[i].errors / placements[i].aligned < placements[b].errors / placements[b].aligned))
        b = i;
    }

    //  If we didn't find a best, b will be invalid; set positions for adding to a new tig.
    //  If we did, save both the position it was placed at, and the tigID it was placed in.

    if (b == UINT32_MAX) {
      if ((enableLog == true) && (logFileFlagSet(LOG_PLACE_UNPLACED)))
        writeLog("read %8u remains unplaced\n", fid);
      placedPos[fid].bgn = 0;
      placedPos[fid].end = RI->readLength(fid);
    }

    else {
      if ((enableLog == true) && (logFileFlagSet(LOG_PLACE_UNPLACED)))
        writeLog("read %8u placed tig %6u (%6u reads) at %8u-%8u (cov %7.5f erate %6.4f)\n",
                 fid, placements[b].tigID, tigs[placements[b].tigID]->ufpath.size(),
                 placements[b].position.bgn, placements[b].position.end,
                 placements[b].fCoverage,
                 placements[b].errors / placements[b].aligned);
      placedTig[fid] = placements[b].tigID;
      placedPos[fid] = placements[b].position;
    }
  }

  //  All reads placed, now just dump them in their correct tigs.

  for (uint32 fid=1; fid<RI->numReads()+1; fid++) {
    Unitig  *tig = NULL;
    ufNode   frg;

    if (tigs.inUnitig(fid) > 0)  //  Already placed, just skip it.
      continue;

    //  If not placed, it's garbage.  These reads were not placed in any tig initially, were not
    //  allowed to seed a tig, and now, could find no place to go.  They're garbage.

    if (placedTig[fid] == 0) {
      if (OG->isContained(fid))
        nFailedContained++;
      else
        nFailed++;
    }

    //  Otherwise, it was placed somewhere, grab the tig.

    else {
      if (OG->isContained(fid))
        nPlacedContained++;
      else
        nPlaced++;

      tig = tigs[placedTig[fid]];
    }

    //  Regardless, add it to the tig.  Logging for this is above.

    if (tig) {
      frg.ident             = fid;
      frg.contained         = 0;
      frg.parent            = 0;
      frg.ahang             = 0;
      frg.bhang             = 0;
      frg.position          = placedPos[fid];

      tig->addRead(frg, 0, false);
    }

    //  Update status.

    if (tig)
      RI->setUnplaced(fid);
    else
      RI->setLeftover(fid);
  }

  //  Cleanup.

  delete [] placedPos;
  delete [] placedTig;

  writeStatus("placeContains()-- Placed %u contained reads and %u unplaced reads.\n", nPlacedContained, nPlaced);
  writeStatus("placeContains()-- Failed to place %u contained reads (too high error suspected) and %u unplaced reads (lack of overlaps suspected).\n", nFailedContained, nFailed);

  //  But wait!  All the tigs need to be sorted.  Well, not really _all_, but the hard ones to sort
  //  are big, and those quite likely had reads added to them, so it's really not worth the effort
  //  of tracking which ones need sorting, since the ones that don't need it are trivial to sort.

  for (uint32 ti=1; ti<tigs.size(); ti++) {
    Unitig *utg = tigs[ti];

    if (utg)
      utg->sort();
  }
}

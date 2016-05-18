
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2016-MAY-17
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "intervalList.H"
#include "stddev.H"

#include <vector>

using namespace std;




void
mergeUnitigs_findPlacements(UnitigVector              &unitigs,
                            ufNode                    *rd,
                            double                     deviation,
                            vector<overlapPlacement>  &validPlacements) {
  vector<overlapPlacement>   placements;

  placeFragUsingOverlaps(unitigs, AS_MAX_ERATE, NULL, rd->ident, placements);

  for (uint32 pi=0; pi<placements.size(); pi++) {
    Unitig *tig   = unitigs[placements[pi].tigID];

    uint32  bgn   = placements[pi].position.min();
    uint32  end   = placements[pi].position.max();

    double  erate = placements[pi].errors / placements[pi].aligned;

    if ((rd->position.min() < end) && (bgn < rd->position.max()))  //  Ignore placements to the same place
      continue;

    if ((placements[pi].fCoverage < 0.99) ||  //  Ignore partially placed reads.
        (tig->ufpath.size() == 1)) {          //  Ignore placements in singletons.
      //writeLog("read %8u tig %6u (%8u-%8u) placed -- tig %6u (%6u reads) at %8u-%8u (cov %7.5f erate %6.4f) - LOW COV or SINGLETON\n",
      //         rd->ident, Unitig::fragIn(rd->ident), rd->position.bgn, rd->position.end,
      //         placements[pi].tigID, tig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
      continue;
    }

    if (tig->overlapConsistentWithTig(deviation, bgn, end, erate) < 0.5) {
      //if ((enableLog == true) && (logFileFlagSet(LOG_PLACE_UNPLACED)))
      //  writeLog("read %8u tig %6u (%8u-%8u) placed -- tig %6u (%6u reads) at %8u-%8u (cov %7.5f erate %6.4f) - HIGH ERROR\n",
      //           rd->ident, Unitig::fragIn(rd->ident), rd->position.bgn, rd->position.end,
      //           placements[pi].tigID, tig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
      continue;
    }

    writeLog("read %8u tig %6u (%8u-%8u) placed -- tig %6u (%6u reads) at %8u-%8u (cov %7.5f erate %6.4f)\n",
             rd->ident, Unitig::fragIn(rd->ident), rd->position.bgn, rd->position.end,
             placements[pi].tigID, tig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);

    validPlacements.push_back(placements[pi]);
  }
}





void
mergeUnitigs(UnitigVector &unitigs,
             double        deviation,
             bool          findCircularTigs) {


  //  For every tig, decide if it can merge, end-to-end, with some other tig.  This operation
  //  should occur before bubbles are popped (so that whatever we chop off can be popped as a
  //  bubble) and repeats are split (so that whatever we join can be split if it's not supported).
  //
  //  The basic idea is that the end read on each tig should align to the middle of the other tig.
  //  If the reads between those also align, we can merge.  If they do not align, we should split
  //  off one end and join.  The split off end is either a bubble, or we made bad joins and will
  //  end up with four pieces after repeat breaking.
  //
  //    -----------------------------------
  //                        ^^^^       ----
  //                        ||||       ||||
  //                        ----       vvvv
  //                        ----------------------------------
  //
  //  This is the same basic operation as for finding circular tigs, and those
  //  are found too.  However, this should occur after bubbles and repeats.




  //  Step 1: For every end read, place it.  Save only placements that are full-length and
  //  compatible with the destination tig.

  vector<overlapPlacement>   validPlacements;

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) ||
        (tig->getNumFrags() < 2) ||
        (tig->_isUnassembled == true))
      continue;

    ufNode  *f = tig->firstRead();
    ufNode  *l = tig->lastRead();

    if (f == l)
      continue;

    mergeUnitigs_findPlacements(unitigs, f, deviation, validPlacements);
    mergeUnitigs_findPlacements(unitigs, l, deviation, validPlacements);
  }

  writeLog("Found "F_SIZE_T" valid placements of end reads.\n", validPlacements.size());

  //  Step 2: Find pairs of placements between two tigs.

  vector<pair<uint32, uint32> >   potentialCircles;
  vector<pair<uint32, uint32> >   potentialMerges;

  for (uint32 pa=0; pa<validPlacements.size(); pa++) {
    uint32  paSrcTigID = Unitig::fragIn(validPlacements[pa].frgID);
    uint32  paDstTigID = validPlacements[pa].tigID;

    for (uint32 pb=pa+1; pb<validPlacements.size(); pb++) {
      uint32  pbSrcTigID = Unitig::fragIn(validPlacements[pb].frgID);
      uint32  pbDstTigID = validPlacements[pb].tigID;

      if (validPlacements[pa].frgID == validPlacements[pb].frgID)  //  Whatever we're trying, we can't use the same read twice.
        continue;

      if ((paDstTigID == paSrcTigID) &&   //  pa placed in same tig as it came from
          (pbDstTigID == pbSrcTigID) &&   //  pb placed in same tig as it came from
          (paSrcTigID == pbSrcTigID)) {   //  and both placed in same tig
        potentialCircles.push_back(pair<uint32,uint32>(pa, pb));
        continue;
      }

      if ((paDstTigID == pbSrcTigID) &&   //  pa placed in same tig as pb came from
          (pbDstTigID == paSrcTigID)) {   //  pb placed in same tig as pa came from
        potentialMerges.push_back(pair<uint32,uint32>(pa, pb));
        continue;
      }
    }
  }

  writeLog("Found "F_SIZE_T" potential circular tigs.\n", potentialCircles.size());
  writeLog("Found "F_SIZE_T" potential joins.\n",         potentialMerges.size());


  //  Step 3: For the potential circles, each read needs to be placed with the same orientation as
  //  its source, and the distance between (paSrc,pbDst) and (pbSrc,paDst) needs to be
  //  (approximately) the same.  Then, we should really check the reads between those two points.

  for (uint32 pc=0; pc<potentialCircles.size(); pc++) {
    uint32  pa         = potentialCircles[pc].first;
    uint32  pb         = potentialCircles[pc].second;

    uint32  paReadID   = validPlacements[pa].frgID;
    uint32  pbReadID   = validPlacements[pb].frgID;

    uint32  tigID      = validPlacements[pa].tigID;  //  All reads placed in the same tig, see above.
    Unitig *tig        = unitigs[tigID];

    ufNode *paRead     = &tig->ufpath[Unitig::pathPosition(paReadID)];
    bool    paSrcFwd   = paRead->position.isForward();
    bool    paDstFwd   = validPlacements[pa].position.isForward();

    ufNode *pbRead     = &tig->ufpath[Unitig::pathPosition(pbReadID)];
    bool    pbSrcFwd   = pbRead->position.isForward();
    bool    pbDstFwd   = validPlacements[pb].position.isForward();

    writeLog("TEST CIRCULAR - tig %u - pa=%u pb=%u - reads %u @ %u-%u -> %u-%u and %u @ %u-%u -> %u-%u\n",
             tigID, pa, pb,
             paReadID, paRead->position.bgn, paRead->position.end, validPlacements[pa].position.bgn, validPlacements[pa].position.end,
             pbReadID, pbRead->position.bgn, pbRead->position.end, validPlacements[pb].position.bgn, validPlacements[pb].position.end);

    if ((paSrcFwd != paDstFwd) ||
        (pbSrcFwd != pbDstFwd)) {
      writeLog("not circular - orient mismatch for tig %u pa %u pb %u reads %u and %u\n",
               tigID, pa, pb, paReadID, pbReadID);
      continue;
    }
  }



  for (uint32 pc=0; pc<potentialMerges.size(); pc++) {
    uint32  pa         = potentialMerges[pc].first;
    uint32  pb         = potentialMerges[pc].second;

    uint32  paReadID   = validPlacements[pa].frgID;
    uint32  pbReadID   = validPlacements[pb].frgID;

    uint32  paTigID    = validPlacements[pa].tigID;  //  All reads placed in the same tig, see above.
    uint32  pbTigID    = validPlacements[pb].tigID;  //  All reads placed in the same tig, see above.

    Unitig *paTig      = unitigs[paTigID];
    Unitig *pbTig      = unitigs[pbTigID];

    ufNode *paRead     = &paTig->ufpath[Unitig::pathPosition(paReadID)];
    bool    paSrcFwd   = paRead->position.isForward();
    bool    paDstFwd   = validPlacements[pa].position.isForward();

    ufNode *pbRead     = &pbTig->ufpath[Unitig::pathPosition(pbReadID)];
    bool    pbSrcFwd   = pbRead->position.isForward();
    bool    pbDstFwd   = validPlacements[pb].position.isForward();

    writeLog("TEST JOIN - pa tig %u read %u @ %u-%u -> %u-%u -- pb tig %u read %u @ %u-%u -> %u-%u\n",
             paTigID, paReadID, paRead->position.bgn, paRead->position.end, validPlacements[pa].position.bgn, validPlacements[pa].position.end,
             pbTigID, pbReadID, pbRead->position.bgn, pbRead->position.end, validPlacements[pb].position.bgn, validPlacements[pb].position.end);
  }




  exit(0);
}

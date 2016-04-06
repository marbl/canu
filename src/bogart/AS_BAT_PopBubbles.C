
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
 *    src/AS_BAT/AS_BAT_MergeSplitJoin.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2011-FEB-15 to 2014-MAY-03
 *      are Copyright 2011-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2015-AUG-05
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_IntersectSplit.H"

#include "AS_BAT_Instrumentation.H"

#include "intervalList.H"

#include <vector>
#include <set>
#include <map>

using namespace std;

#define BUBBLE_CONSISTENT  3.0

#undef  LOG_BUBBLE_TESTS
#undef  LOG_BUBBLE_FAILURE
#define LOG_BUBBLE_SUCCESS

#define BUBBLINESS_DETAILS


class candidatePop {
public:
  candidatePop() {
  };
  candidatePop(Unitig *bubble_, Unitig *target_, uint32 bgn_, uint32 end_) {
    bubble = bubble_;
    target = target_;
    bgn    = bgn_;
    end    = end_;
  };
  ~candidatePop() {
  };

  Unitig  *bubble;
  Unitig  *target;

  uint32   bgn;
  uint32   end;

  vector<overlapPlacement>  placed;
};


//  A list of the target unitigs that a bubble could be popped into.
typedef  map<uint32, vector<uint32> >  BubTargetList;



  //  Decide which unitigs can be bubbles.  The first pass finds unitigs that can be potential
  //  bubbles.  Any unitig where every dovetail read has an overlap to some other unitig is a
  //  candidate for bubble popping.
void
findPotentialBubbles(UnitigVector    &unitigs,
                     double           erateBubble,
                     BubTargetList   &potentialBubbles) {
  uint32  tiLimit      = unitigs.size();
  uint32  tiNumThreads = omp_get_max_threads();
  uint32  tiBlockSize  = (tiLimit < 100000 * tiNumThreads) ? tiNumThreads : tiLimit / 99999;

  writeLog("bubbleDetect()-- working on "F_U32" unitigs, with "F_U32" threads.\n", tiLimit, tiNumThreads);

  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) ||               //  Not a tig, ignore it.
        (tig->ufpath.size() == 1))     //  Singleton, handled elsewhere
      continue;

    uint32  nonContainedReads = 0;
    bool    validBubble       = true;

    map<uint32,uint32>  tigOlapsTo;

    uint32  fiLimit      = tig->ufpath.size();
    uint32  fiNumThreads = omp_get_max_threads();
    uint32  fiBlockSize  = (fiLimit < 100 * fiNumThreads) ? fiNumThreads : fiLimit / 99;

    for (uint32 fi=0; (validBubble == true) && (fi<fiLimit); fi++) {
      uint32      rid      = tig->ufpath[fi].ident;

      if (OG->isContained(rid) == true)
        continue;

      nonContainedReads++;

      uint32      ovlLen   = 0;
      BAToverlap *ovl      = OC->getOverlaps(rid, AS_MAX_ERATE, ovlLen);

      set<uint32>  readOlapsTo;

      for (uint32 oi=0; oi<ovlLen; oi++) {
        uint32  ovlTigID = Unitig::fragIn(ovl[oi].b_iid);
        Unitig *ovlTig   = unitigs[ovlTigID];

        //  Skip this overlap if it is to an unplaced read, to a singleton tig, to ourself,
        //  or to a unitig that is shorter than us.  We can not pop this tig as a bubble
        //  in any of those cases.

        if ((ovlTigID == 0) ||
            (ovlTig == NULL) ||
            (ovlTig->ufpath.size() == 1) ||
            (ovlTig->id() == tig->id()) ||
            (ovlTig->getLength() < tig->getLength()))
          continue;

        //  Otherwise, remember that we had an overlap to ovlTig.

        readOlapsTo.insert(ovlTigID);
      }

      //  Transfer the per-read counts to the per-unitig counts.  Decide if we're a valid potential bubble.

      validBubble = false;

      for (set<uint32>::iterator it=readOlapsTo.begin(); it != readOlapsTo.end(); ++it) {
        tigOlapsTo[*it]++;

        if (tigOlapsTo[*it] == nonContainedReads)
          validBubble = true;
      }
    }

    //  If validBubble, then there is a tig that every dovetail read has at least one overlap to.
    //  Save those tigs in potentialBubbles.

    if (validBubble) {
      uint32  nTigs = 0;

      for (map<uint32,uint32>::iterator it=tigOlapsTo.begin(); it != tigOlapsTo.end(); ++it)
        if (it->second == nonContainedReads)
          nTigs++;

      writeLog("potential bubble tig %8u length %9u nReads %7u to %3u tigs:\n",
               tig->id(), tig->getLength(), tig->ufpath.size(), nTigs);

      for (map<uint32,uint32>::iterator it=tigOlapsTo.begin(); it != tigOlapsTo.end(); ++it) {
        Unitig  *dest = unitigs[it->first];

        if (it->second != nonContainedReads)
          continue;

        writeLog("                 tig %8u length %9u nReads %7u\n", dest->id(), dest->getLength(), dest->ufpath.size());

        potentialBubbles[ti].push_back(dest->id());
      }
    }
  }

  flushLog();
}
                     



//  Find filtered placements for all the reads in the potential bubble tigs.

vector<overlapPlacement>  *
findBubbleReadPlacements(UnitigVector    &unitigs,
                         BubTargetList   &potentialBubbles) {
  uint32  fiLimit      = FI->numFragments();
  uint32  fiNumThreads = omp_get_max_threads();
  uint32  fiBlockSize  = (fiLimit < 1000 * fiNumThreads) ? fiNumThreads : fiLimit / 999;

  vector<overlapPlacement>   *placed = new vector<overlapPlacement> [fiLimit + 1];

#pragma omp parallel for schedule(dynamic, fiBlockSize)
  for (uint32 fi=0; fi<fiLimit; fi++) {
    uint32     rdAtigID = Unitig::fragIn(fi);

    if ((rdAtigID == 0) ||                           //  Read not placed in a tig, ignore it.
        (OG->isContained(fi)) ||                     //  Read is contained, ignore it.
        (potentialBubbles.count(rdAtigID) == 0))     //  Read isn't in a potential bubble, ignore it.
      continue;

    Unitig     *rdAtig   = unitigs[rdAtigID];
    ufNode     *rdA      = &rdAtig->ufpath[ Unitig::pathPosition(fi) ];
    bool        rdAfwd   = (rdA->position.bgn < rdA->position.end);
    int32       rdAlo    = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
    int32       rdAhi    = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

    uint32      ovlLen   = 0;
    BAToverlap *ovl      = OC->getOverlaps(rdA->ident, AS_MAX_ERATE, ovlLen);

    set<uint32> intersections;

    if ((fi % 100) == 0)
      fprintf(stderr, "bubbliness()-- read %8u with %6u overlaps - %6.2f%% finished.\r",
              rdA->ident, ovlLen, 100.0 * fi / fiLimit);

    //  Compute all placements for this read.

    vector<overlapPlacement>   placements;

    placeFragUsingOverlaps(unitigs, AS_MAX_ERATE, NULL, rdA->ident, placements);

    //  Weed out placements that aren't for bubbles, or that are for bubbles but are poor quality.

    for (uint32 pi=0; pi<placements.size(); pi++) {
      uint32    rdBtigID = placements[pi].tigID;
      Unitig   *rdBtig   = unitigs[rdBtigID];

      //  Ignore the placement if it is to a non-tig / singleton read, or if it didn't place the
      //  read fully.

      if ((rdBtigID == 0) ||
          (rdBtig   == NULL) ||
          (rdBtig->ufpath.size() == 1) ||
          (placements[pi].fCoverage < 0.99))
        continue;

      //  Ignore the placement if it isn't to one of our bubble-popping candidate unitigs.

      bool             dontcare = true;
      vector<uint32>  &pbubbles = potentialBubbles[rdAtigID];

      for (uint32 pb=0; pb<pbubbles.size(); pb++) {
        if (pbubbles[pb] == rdBtigID)
          dontcare = false;
      }

      if (dontcare)
        continue;

      //  Ignore the placement if it is too diverged from the destination tig.

      uint32  lo = (placements[pi].position.bgn < placements[pi].position.end) ? placements[pi].position.bgn : placements[pi].position.end;
      uint32  hi = (placements[pi].position.bgn < placements[pi].position.end) ? placements[pi].position.end : placements[pi].position.bgn;

      double  erate = placements[pi].errors / placements[pi].aligned;

      if (rdBtig->overlapConsistentWithTig(BUBBLE_CONSISTENT, lo, hi, erate) < 0.5) {
#ifdef SHOW_PLACEMENT_DETAIL
        writeLog("frag %8u placed in tig %6u (%6u reads) at %8u-%8u with coverage %7.5f and erate %6.4f - HIGH ERROR\n",
                 fid, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
#endif
        continue;
      }

#ifdef SHOW_PLACEMENT_DETAIL
      writeLog("frag %8u placed in tig %6u (%6u reads) at %8u-%8u with coverage %7.5f and erate %6.4f\n",
               fid, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
#endif

      //  Good placement!

      placed[fi].push_back(placements[pi]);
    }
  }

  return(placed);
}






//  Bubble popping cannot be done in parallel -- there is a race condition when both unitigs
//  A and B are considering merging in unitig C.

void
popBubbles(UnitigVector &unitigs,
           double UNUSED(erateGraph), double erateBubble, double UNUSED(erateMerge), double UNUSED(erateRepeat),
           const char *UNUSED(prefix),
           uint32 UNUSED(minOverlap),
           uint64 UNUSED(genomeSize)) {

  BubTargetList   potentialBubbles;

  findPotentialBubbles(unitigs, erateBubble, potentialBubbles);

  vector<overlapPlacement>   *placed = findBubbleReadPlacements(unitigs, potentialBubbles);

  //  We now have, in 'placed', a list of all the places that each read could be placed.  Decide if there is a _single_
  //  place for each bubble to be popped.

  uint32  tiLimit      = unitigs.size();
  //uint32  tiNumThreads = omp_get_max_threads();
  //uint32  tiBlockSize  = (tiLimit < 100000 * tiNumThreads) ? tiNumThreads : tiLimit / 99999;

  for (uint32 ti=0; ti<tiLimit; ti++) {
    if (potentialBubbles.count(ti) == 0)   //  Not a potential bubble
      continue;

    Unitig  *bubble = unitigs[ti];

    //  Split the placements into piles for each target and build an interval list for each target.

    map<uint32, intervalList<uint32> *>  targetIntervals;

    //  For each read in the tig, convert the vector of placements into interval lists, one list per target tig.

    for (uint32 fi=0; fi<bubble->ufpath.size(); fi++) {
      uint32  readID  = bubble->ufpath[fi].ident;

      for (uint32 pp=0; pp<placed[readID].size(); pp++) {
        uint32  tid = placed[readID][pp].tigID;

        uint32  bgn = (placed[readID][pp].position.bgn < placed[readID][pp].position.end) ? placed[readID][pp].position.bgn : placed[readID][pp].position.end;
        uint32  end = (placed[readID][pp].position.bgn < placed[readID][pp].position.end) ? placed[readID][pp].position.end : placed[readID][pp].position.bgn;

        if (targetIntervals[tid] == NULL)
          targetIntervals[tid] = new intervalList<uint32>;

        targetIntervals[tid]->add(bgn, end-bgn);
      }
    }

    vector<candidatePop *>    targets;

    //  Squish the intervals.  Create new candidatePops for each interval that isn't too big or
    //  small.  Assign each overlapPlacements to the correct candidatePop.

    for (map<uint32, intervalList<uint32> *>::iterator it=targetIntervals.begin(); it != targetIntervals.end(); ++it) {
      uint32                 targetID = it->first;
      intervalList<uint32>  *IL       = it->second;

      IL->merge();

      //  Discard intervals that are significantly too small or large.  Save the ones that are nicely sized.

      for (uint32 ii=0; ii<IL->numberOfIntervals(); ii++) {
        writeLog("tig %u length %u -> target %u piece %u position %u-%u length %u\n",
                 bubble->id(), bubble->getLength(),
                 targetID, ii, IL->lo(ii), IL->hi(ii), IL->hi(ii) - IL->lo(ii));

        if (IL->hi(ii) - IL->lo(ii) < 0.75 * bubble->getLength())   //  Too small!
          continue;

        if (1.25 * bubble->getLength() < IL->hi(ii) - IL->lo(ii))     //  Too big!
          continue;

        targets.push_back(new candidatePop(bubble, unitigs[targetID], IL->lo(ii), IL->hi(ii)));
      }

      delete IL;
    }

    targetIntervals.clear();

    //  Run through the placements again, and assign them to the correct target.
    //
    //  For each read:
    //  For each acceptable placement:
    //  For each target location:
    //  If the placement is for this target, save it.

    for (uint32 fi=0; fi<bubble->ufpath.size(); fi++) {
      uint32  readID  = bubble->ufpath[fi].ident;

      for (uint32 pp=0; pp<placed[readID].size(); pp++) {
        uint32  tid = placed[readID][pp].tigID;

        uint32  bgn = (placed[readID][pp].position.bgn < placed[readID][pp].position.end) ? placed[readID][pp].position.bgn : placed[readID][pp].position.end;
        uint32  end = (placed[readID][pp].position.bgn < placed[readID][pp].position.end) ? placed[readID][pp].position.end : placed[readID][pp].position.bgn;

        for (uint32 tt=0; tt<targets.size(); tt++)
          if ((targets[tt]->target->id() == tid) &&
              (targets[tt]->bgn < end) && (bgn < targets[tt]->end))
            targets[tt]->placed.push_back(placed[readID][pp]);
      }
    }

    //  Count the number of targets that have all the reads (in the correct order, etc, etc).  Remove those
    //  that don't.

    uint32  nTargets = 0;

    set<uint32>  tigReads;
    set<uint32>  tgtReads;

    for (uint32 fi=0; fi<bubble->ufpath.size(); fi++)
      tigReads.insert(bubble->ufpath[fi].ident);


    for (uint32 tt=0; tt<targets.size(); tt++) {
      tgtReads.clear();

      for (uint32 op=0; op<targets[tt]->placed.size(); op++)
        tgtReads.insert(targets[tt]->placed[op].frgID);

      uint32  n5 = 0;
      uint32  n3 = 0;

      for (uint32 fi=0; fi<bubble->ufpath.size(); fi++)
        if (tgtReads.count(bubble->ufpath[fi].ident) > 0)
          n5++;
        else
          break;

      for (uint32 fi=bubble->ufpath.size(); fi-->0; )
        if (tgtReads.count(bubble->ufpath[fi].ident) > 0)
          n3++;
        else
          break;

      writeLog("target %u piece %u position %u-%u length %u - expected "F_SIZE_T" reads, had "F_SIZE_T" reads.  n5=%u n3=%u\n",
               targets[tt]->target->id(), tt, targets[tt]->bgn, targets[tt]->end, targets[tt]->end - targets[tt]->bgn,
               tigReads.size(),
               tgtReads.size(), n5, n3);

      if (tigReads == tgtReads) {
        nTargets++;
      } else {
        targets[tt]->bubble = NULL;
        targets[tt]->target = NULL;
      }
    }


    if (nTargets == 0) {
      writeLog("no targets\n");
    }

    else if (nTargets == 1) {
      writeLog("one target\n");
    }

    else {
      writeLog("many targets\n");
    }
  }

  delete [] placed;
}

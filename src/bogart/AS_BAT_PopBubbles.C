
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
 *    Brian P. Walenz beginning on 2016-MAR-11
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

#include "AS_BAT_Instrumentation.H"

#include "intervalList.H"

#include <vector>
#include <set>
#include <map>

using namespace std;

#define BUBBLE_READ_FRACTION  0.5

#undef  SHOW_MULTIPLE_PLACEMENTS  //  Reports reads that are placed multiple times in a single target region

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
                     BubTargetList   &potentialBubbles) {
  uint32  tiLimit      = unitigs.size();
  uint32  tiNumThreads = omp_get_max_threads();
  uint32  tiBlockSize  = (tiLimit < 100000 * tiNumThreads) ? tiNumThreads : tiLimit / 99999;

  writeLog("bubbleDetect()-- working on "F_U32" unitigs, with "F_U32" threads.\n", tiLimit, tiNumThreads);

  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) ||               //  Not a tig, ignore it.
        (tig->ufpath.size() == 1))     //  Singleton, handled elsewhere.
      continue;

    uint32  nonContainedReads = 0;
    bool    validBubble       = true;

    map<uint32,uint32>  tigOlapsTo;

    uint32  fiLimit      = tig->ufpath.size();
    uint32  fiNumThreads = omp_get_max_threads();
    uint32  fiBlockSize  = (fiLimit < 100 * fiNumThreads) ? fiNumThreads : fiLimit / 99;

    for (uint32 fi=0; (validBubble == true) && (fi<fiLimit); fi++) {
      uint32      rid      = tig->ufpath[fi].ident;

      if (OG->isContained(rid) == true)  //  Don't need to check contained reads.  If their container
        continue;                        //  passes the tests below, the contained read will too.

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

        //writeLog("tig %u read %u overlap to tig %u read %u\n",
        //         tig->id(), rid, ovlTigID, ovl[oi].b_iid);

        readOlapsTo.insert(ovlTigID);
      }

      //writeLog("tig %8u read %8u has %u olaps\n", tig->id(), rid, readOlapsTo.size());

      //  Transfer the per-read counts to the per-unitig counts:  add one to the counter for each tig
      //  that we have overlaps to.

      for (set<uint32>::iterator it=readOlapsTo.begin(); it != readOlapsTo.end(); ++it)
        tigOlapsTo[*it]++;

      //  Decide if we're a valid potential bubble.  If tig id (in it->first) has overlaps to every
      //  read we've seen so far (nonContainedReads), we're still a valid bubble.
      //
      //  To _attempt_ to have differences in the bubble, we'll accept it if 3/4 of the reads
      //  have overlaps.

      validBubble = false;

      for (map<uint32,uint32>::iterator it=tigOlapsTo.begin(); it != tigOlapsTo.end(); ++it)
        if (it->second >= BUBBLE_READ_FRACTION * nonContainedReads)
          validBubble = true;

      //  If we've not seen that many reads, pretend it's a valid bubble.  It'll get screened out later.

      if (nonContainedReads < 16)
        validBubble = true;
    }

    //  If not validBubble, report.

#if 0
    if (validBubble == false) {
      writeLog("notValidBubble tig %8d expects %6u reads\n", tig->id(), nonContainedReads);

      for (map<uint32,uint32>::iterator it=tigOlapsTo.begin(); it != tigOlapsTo.end(); ++it)
        writeLog("  to tig %8u overlaps %6u\n", it->first, it->second);
    }
#endif

    //  If validBubble, then there is a tig that every dovetail read has at least one overlap to.
    //  Save those tigs in potentialBubbles.

    uint32  nTigs = 0;

    if (validBubble) {
      for (map<uint32,uint32>::iterator it=tigOlapsTo.begin(); it != tigOlapsTo.end(); ++it)
        if (it->second >= BUBBLE_READ_FRACTION * nonContainedReads)
          nTigs++;
    }

    //  ALWAYS log potential bubbles.

    if (nTigs > 0) {
      writeLog("\n");
      writeLog("potential bubble tig %8u length %9u nReads %7u to %3u tigs:\n",
               tig->id(), tig->getLength(), tig->ufpath.size(), nTigs);

      for (map<uint32,uint32>::iterator it=tigOlapsTo.begin(); it != tigOlapsTo.end(); ++it) {
        if (it->second >= BUBBLE_READ_FRACTION * nonContainedReads) {
          Unitig  *dest = unitigs[it->first];

          writeLog("                 tig %8u length %9u nReads %7u\n", dest->id(), dest->getLength(), dest->ufpath.size());

          potentialBubbles[ti].push_back(dest->id());
        }
      }
    }
  }

  flushLog();
}




//  Find filtered placements for all the reads in the potential bubble tigs.

vector<overlapPlacement>  *
findBubbleReadPlacements(UnitigVector    &unitigs,
                         BubTargetList   &potentialBubbles,
                         double           deviationBubble) {
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

    //if ((fi % 100) == 0)
    //  fprintf(stderr, "findBubbleReadPlacements()-- read %8u with %6u overlaps - %6.2f%% finished.\r",
    //          rdA->ident, ovlLen, 100.0 * fi / fiLimit);

    //  Compute all placements for this read.

    vector<overlapPlacement>   placements;

    placeFragUsingOverlaps(unitigs, AS_MAX_ERATE, NULL, rdA->ident, placements);

    //  Weed out placements that aren't for bubbles, or that are for bubbles but are poor quality.  Or are to ourself!

    for (uint32 pi=0; pi<placements.size(); pi++) {
      uint32    rdBtigID = placements[pi].tigID;
      Unitig   *rdBtig   = unitigs[rdBtigID];

      uint32    lo       = (placements[pi].position.bgn < placements[pi].position.end) ? placements[pi].position.bgn : placements[pi].position.end;
      uint32    hi       = (placements[pi].position.bgn < placements[pi].position.end) ? placements[pi].position.end : placements[pi].position.bgn;

      double    erate    = placements[pi].errors / placements[pi].aligned;

      //  Ignore the placement if it is to ourself.

      if (rdAtigID == rdBtigID) {
        //writeLog("tig %6u frag %8u -> tig %6u %6u reads at %8u-%8u (cov %7.5f erate %6.4f) - SAME TIG\n",
        //         rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Ignore the placement if it is to a non-tig / singleton read, or if it didn't place the
      //  read fully.

      if ((rdBtigID == 0) ||
          (rdBtig   == NULL) ||
          (rdBtig->ufpath.size() == 1) ||
          (placements[pi].fCoverage < 0.99)) {
        if (logFileFlagSet(LOG_BUBBLE_DETAIL))
          writeLog("tig %6u frag %8u -> tig %6u %6u reads at %8u-%8u (cov %7.5f erate %6.4f) - PARTIALLY PLACED\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Ignore the placement if it isn't to one of our bubble-popping candidate unitigs.

      bool             dontcare = true;
      vector<uint32>  &pbubbles = potentialBubbles[rdAtigID];

      for (uint32 pb=0; pb<pbubbles.size(); pb++) {
        if (pbubbles[pb] == rdBtigID)
          dontcare = false;
      }

      if (dontcare) {
        if (logFileFlagSet(LOG_BUBBLE_DETAIL))
          writeLog("tig %6u frag %8u -> tig %6u %6u reads at %8u-%8u (cov %7.5f erate %6.4f) - NOT CANDIDATE TIG\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Ignore the placement if it is too diverged from the destination tig.

      if (rdBtig->overlapConsistentWithTig(deviationBubble, lo, hi, erate) < 0.5) {
        if (logFileFlagSet(LOG_BUBBLE_DETAIL))
          writeLog("tig %6u frag %8u -> tig %6u %6u reads at %8u-%8u (cov %7.5f erate %6.4f) - HIGH ERROR\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Good placement!

      if (logFileFlagSet(LOG_BUBBLE_DETAIL))
        writeLog("tig %6u frag %8u -> tig %6u %6u reads at %8u-%8u (cov %7.5f erate %6.4f)\n",
                 rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);

      placed[fi].push_back(placements[pi]);
    }
  }

  return(placed);
}






//  Bubble popping cannot be done in parallel -- there is a race condition when both unitigs
//  A and B are considering merging in unitig C.

void
popBubbles(UnitigVector &unitigs,
           double deviationBubble) {

  BubTargetList   potentialBubbles;

  findPotentialBubbles(unitigs, potentialBubbles);

  writeLog("\n");
  writeLog("Found "F_SIZE_T" potential bubbles.\n", potentialBubbles.size());
  writeLog("\n");

  vector<overlapPlacement>   *placed = findBubbleReadPlacements(unitigs, potentialBubbles, deviationBubble);

  //  We now have, in 'placed', a list of all the places that each read could be placed.  Decide if there is a _single_
  //  place for each bubble to be popped.

  uint32  tiLimit      = unitigs.size();
  //uint32  tiNumThreads = omp_get_max_threads();
  //uint32  tiBlockSize  = (tiLimit < 100000 * tiNumThreads) ? tiNumThreads : tiLimit / 99999;

  //  Clear flags.
  for (uint32 ti=0; ti<tiLimit; ti++) {
    if (unitigs[ti]) {
      unitigs[ti]->_isBubble = false;
      unitigs[ti]->_isRepeat = false;
    }
  }

  //  In parallel, process the placements.

  for (uint32 ti=0; ti<tiLimit; ti++) {
    if (potentialBubbles.count(ti) == 0)   //  Not a potential bubble
      continue;

    //  Scan the bubble, decide if there are _ANY_ read placements.  Log appropriately.

    Unitig  *bubble = unitigs[ti];
    bool     hasPlacements = false;

    for (uint32 fi=0; fi<bubble->ufpath.size(); fi++) {
      uint32  readID  = bubble->ufpath[fi].ident;

      if (placed[readID].size() > 0)
        hasPlacements = true;
    }

    if (hasPlacements == false)
      writeLog("potential bubble %u had no valid placements (all were not contained in target tig)\n", ti);
    else
      writeLog("potential bubble %u\n", ti);

    //  Split the placements into piles for each target and build an interval list for each target.
    //  For each read in the tig, convert the vector of placements into interval lists, one list per target tig.

    map<uint32, intervalList<uint32> *>  targetIntervals;

    for (uint32 fi=0; fi<bubble->ufpath.size(); fi++) {
      uint32  readID  = bubble->ufpath[fi].ident;

      for (uint32 pp=0; pp<placed[readID].size(); pp++) {
        uint32  tid = placed[readID][pp].tigID;

        assert(placed[readID][pp].frgID > 0);

        uint32  bgn = (placed[readID][pp].position.bgn < placed[readID][pp].position.end) ? placed[readID][pp].position.bgn : placed[readID][pp].position.end;
        uint32  end = (placed[readID][pp].position.bgn < placed[readID][pp].position.end) ? placed[readID][pp].position.end : placed[readID][pp].position.bgn;

        if (targetIntervals[tid] == NULL)
          targetIntervals[tid] = new intervalList<uint32>;

        //writeLog("read %u -> tig %u intervals %u-%u\n", readID, tid, bgn, end);

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

      //  Discard intervals that are significantly too small or large.  Save the ones that are
      //  nicely sized.  Logging here isn't terribly useful, it's just repeated (out of order) later
      //  when we try to make sense of the read alignments.

      for (uint32 ii=0; ii<IL->numberOfIntervals(); ii++) {
        if ((IL->hi(ii) - IL->lo(ii) < 0.75 * bubble->getLength()) ||   //  Too small!
            (1.25 * bubble->getLength() < IL->hi(ii) - IL->lo(ii))) {   //  Too big!
          writeLog("tig %8u length %9u -> target %8u piece %2u position %9u-%9u length %8u - size mismatch, discarded\n",
                   bubble->id(), bubble->getLength(),
                   targetID, ii, IL->lo(ii), IL->hi(ii), IL->hi(ii) - IL->lo(ii));
          continue;
        }

        writeLog("tig %8u length %9u -> target %8u piece %2u position %9u-%9u length %8u\n",
                 bubble->id(), bubble->getLength(),
                 targetID, ii, IL->lo(ii), IL->hi(ii), IL->hi(ii) - IL->lo(ii));

        targets.push_back(new candidatePop(bubble, unitigs[targetID], IL->lo(ii), IL->hi(ii)));
      }

      delete IL;
    }

    targetIntervals.clear();

    //  If no targets, nothing to do.

    if (targets.size() == 0)
      continue;

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

    //  Count the number of targets that have all the reads (later: in the correct order, etc, etc).  Remove those
    //  that don't.

    uint32  nTargets = 0;

    set<uint32>  tigReads;  //  Reads in the bubble tig.
    set<uint32>  tgtReads;  //  Reads in the bubble that have a placement in the target.

    //  Remove duplicate placements from each target.

    for (uint32 tt=0; tt<targets.size(); tt++) {
      candidatePop *t = targets[tt];

      //  Detect duplicates, keep the one with lower error.  There are a lot of duplicate
      //  placements, logging isn't terribly useful.

      for (uint32 aa=0; aa<t->placed.size(); aa++) {
        for (uint32 bb=0; bb<t->placed.size(); bb++) {
          if ((aa == bb) ||
              (t->placed[aa].frgID != t->placed[bb].frgID) ||
              (t->placed[aa].frgID == 0) ||
              (t->placed[bb].frgID == 0))
            continue;

          if (t->placed[aa].errors / t->placed[aa].aligned < t->placed[bb].errors / t->placed[bb].aligned) {
#ifdef SHOW_MULTIPLE_PLACEMENTS
            writeLog("duplicate read alignment for tig %u read %u - better %u-%u %.4f - worse %u-%u %.4f\n",
                     t->placed[aa].tigID, t->placed[aa].frgID,
                     t->placed[aa].position.bgn, t->placed[aa].position.end, t->placed[aa].errors / t->placed[aa].aligned,
                     t->placed[bb].position.bgn, t->placed[bb].position.end, t->placed[bb].errors / t->placed[bb].aligned);
#endif
            t->placed[bb] = overlapPlacement();
          } else {
#ifdef SHOW_MULTIPLE_PLACEMENTS
            writeLog("duplicate read alignment for tig %u read %u - better %u-%u %.4f - worse %u-%u %.4f\n",
                     t->placed[aa].tigID, t->placed[aa].frgID,
                     t->placed[bb].position.bgn, t->placed[bb].position.end, t->placed[bb].errors / t->placed[bb].aligned,
                     t->placed[aa].position.bgn, t->placed[aa].position.end, t->placed[aa].errors / t->placed[aa].aligned);
#endif
            t->placed[aa] = overlapPlacement();
          }
        }
      }

      //  Get rid of any now-empty entries.

      for (uint32 aa=t->placed.size(); aa--; ) {
        if (t->placed[aa].frgID == 0) {
          t->placed[aa] = t->placed.back();
          t->placed.pop_back();
        }
      }
    }

    //  Make a set of the reads in the bubble.  We'll compare each target against this to decide if all reads are placed.

    for (uint32 fi=0; fi<bubble->ufpath.size(); fi++)
      tigReads.insert(bubble->ufpath[fi].ident);

    uint32   nOrphan      = 0;   //  Full coverage; bubble can be popped.
    uint32   orphanTarget = 0;

    uint32   nBubble      = 0;   //  Partial coverage, bubble cannot be popped.
    uint32   bubbleTarget = 0;

    for (uint32 tt=0; tt<targets.size(); tt++) {
      tgtReads.clear();

      for (uint32 op=0; op<targets[tt]->placed.size(); op++) {
        if (logFileFlagSet(LOG_BUBBLE_DETAIL))
          writeLog("tig %8u length %9u -> target %8u piece %2u position %9u-%9u length %8u - read %7u at %9u-%9u\n",
                   bubble->id(), bubble->getLength(),
                   targets[tt]->target->id(), tt, targets[tt]->bgn, targets[tt]->end, targets[tt]->end - targets[tt]->bgn,
                   targets[tt]->placed[op].frgID,
                   targets[tt]->placed[op].position.bgn, targets[tt]->placed[op].position.end);

        assert(targets[tt]->placed[op].frgID > 0);
        tgtReads.insert(targets[tt]->placed[op].frgID);
      }

      //  Count the number of consecutive reads from the 5' or 3' end of the bubble that are placed
      //  in the target.
      //
      //  Also, count the number of reads in the bubble that are placed in the target.  Likely the
      //  same as n5 + n3.

      uint32  n5 = 0;
      uint32  n3 = 0;
      uint32  nt = 0;

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


      for (uint32 fi=0; fi<bubble->ufpath.size(); fi++)
        if (tgtReads.count(bubble->ufpath[fi].ident) > 0)
          nt++;


      //  Report now, before we nuke targets[tt] for being not a bubble!

      if ((nt == bubble->ufpath.size()) ||
          ((n5 > 0) && (n3 > 0)))
        writeLog("tig %8u length %9u -> target %8u piece %2u position %9u-%9u length %8u - expected %3"F_SIZE_TP" reads, had %3"F_SIZE_TP" reads.  n5=%3u n3=%3u nt=%3u\n",
                 bubble->id(), bubble->getLength(),
                 targets[tt]->target->id(), tt, targets[tt]->bgn, targets[tt]->end, targets[tt]->end - targets[tt]->bgn,
                 tigReads.size(),
                 tgtReads.size(), n5, n3, nt);

      //  Decide if this is a bubble, orphan from construction, or repeat.

      if (nt == bubble->ufpath.size()) {
        nOrphan++;
        orphanTarget = tt;
      }

      else if ((n5 > 0) && (n3 > 0)) {
        nBubble++;
        bubbleTarget = tt;
      }
    }

    //  If no placements, pbbbt.

    if (nOrphan + nBubble == 0) {
      //writeLog("tig %8u length %8u reads %6u had no bubble or orphan placements.\n", bubble->id(), bubble->getLength(), bubble->ufpath.size());
      continue;
    }

    //  If multiple orphan and/or bubble placements, it's a repeat.

    if (nOrphan + nBubble > 1) {
      writeLog("tig %8u length %8u reads %6u - repeat - %u orphan %u bubble placements.\n",
               bubble->id(), bubble->getLength(), bubble->ufpath.size(),
               nOrphan, nBubble);
      writeLog("\n");
      bubble->_isRepeat = true;
      continue;
    }

    //  If a bubble placement, mark it as a bubble so it can be skipped during repeat detection.

    if (nBubble > 0) {
      writeLog("tig %8u length %8u reads %6u - bubble\n",
               bubble->id(), bubble->getLength(), bubble->ufpath.size());
      writeLog("\n");
      bubble->_isBubble = true;
      continue;
    }

    //  Otherwise, it's an orphan, move the reads to the proper place.

    writeLog("tig %8u length %8u reads %6u - orphan\n", bubble->id(), bubble->getLength(), bubble->ufpath.size());

    for (uint32 op=0, tt=orphanTarget; op<targets[tt]->placed.size(); op++) {
      ufNode  frg;

      frg.ident        = targets[tt]->placed[op].frgID;
      frg.contained    = 0;
      frg.parent       = 0;
      frg.ahang        = 0;
      frg.bhang        = 0;
      frg.position.bgn = targets[tt]->placed[op].position.bgn;
      frg.position.end = targets[tt]->placed[op].position.end;

      writeLog("move read %u from tig %u to tig %u %u-%u\n",
               frg.ident,
               bubble->id(),
               targets[tt]->target->id(), frg.position.bgn, frg.position.end);

      targets[tt]->target->addFrag(frg, 0, false);
    }

    writeLog("\n");

    unitigs[bubble->id()] = NULL;
    delete bubble;
  }  //  Over all bubbles

  writeLog("\n");   //  Needed if no bubbles are popped.

  delete [] placed;

  //  Sort reads in all the tigs.  Overkill, but correct.

  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) ||               //  Not a tig, ignore it.
        (tig->ufpath.size() == 1))     //  Singleton, already sorted.
      continue;

    tig->sort();
  }
}


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
 *    src/bogart/AS_BAT_PopBubbles.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2016-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_PlaceReadUsingOverlaps.H"

#include "AS_BAT_Instrumentation.H"

#include "AS_BAT_MergeOrphans.H"

#include "intervalList.H"

#include <vector>
#include <set>
#include <map>

using namespace std;

#define BUBBLE_READ_FRACTION  0.5

#undef  SHOW_MULTIPLE_PLACEMENTS  //  Reports reads that are placed multiple times in a single target region

class candidatePop {
public:
  candidatePop(Unitig *bubble_, Unitig *target_, uint32 bgn_, uint32 end_) {
    bubble = bubble_;
    target = target_;
    bgn    = bgn_;
    end    = end_;
  };

  Unitig  *bubble;
  Unitig  *target;

  uint32   bgn;
  uint32   end;

  vector<overlapPlacement>  placed;
};


//  A list of the target tigs that a bubble could be popped into.
typedef  map<uint32, vector<uint32> >  BubTargetList;



//  Decide which tigs can be bubbles.  The first pass finds tigs that can be potential
//  bubbles.  Any unitig where every dovetail read has an overlap to some other unitig is a
//  candidate for bubble popping.

void
findPotentialBubbles(TigVector       &tigs,
                     BubTargetList   &potentialBubbles) {
  uint32  tiLimit      = tigs.size();
  uint32  tiNumThreads = omp_get_max_threads();
  uint32  tiBlockSize  = (tiLimit < 100000 * tiNumThreads) ? tiNumThreads : tiLimit / 99999;

  writeStatus("\n");
  writeStatus("bubbleDetect()-- working on " F_U32 " tigs, with " F_U32 " thread%s.\n", tiLimit, tiNumThreads, (tiNumThreads == 1) ? "" : "s");

  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = tigs[ti];

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
      BAToverlap *ovl      = OC->getOverlaps(rid, ovlLen);

      set<uint32>  readOlapsTo;

      for (uint32 oi=0; oi<ovlLen; oi++) {
        uint32  ovlTigID = tigs.inUnitig(ovl[oi].b_iid);
        Unitig *ovlTig   = tigs[ovlTigID];

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
          Unitig  *dest = tigs[it->first];

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
findBubbleReadPlacements(TigVector       &tigs,
                         BubTargetList   &potentialBubbles,
                         double           deviationBubble) {
  uint32  fiLimit      = RI->numReads();
  uint32  fiNumThreads = omp_get_max_threads();
  uint32  fiBlockSize  = (fiLimit < 1000 * fiNumThreads) ? fiNumThreads : fiLimit / 999;

  vector<overlapPlacement>   *placed = new vector<overlapPlacement> [fiLimit + 1];

#pragma omp parallel for schedule(dynamic, fiBlockSize)
  for (uint32 fi=0; fi<fiLimit; fi++) {
    uint32     rdAtigID = tigs.inUnitig(fi);

    if ((rdAtigID == 0) ||                           //  Read not placed in a tig, ignore it.
        (OG->isContained(fi)) ||                     //  Read is contained, ignore it.
        (potentialBubbles.count(rdAtigID) == 0))     //  Read isn't in a potential bubble, ignore it.
      continue;

    Unitig     *rdAtig   = tigs[rdAtigID];
    ufNode     *rdA      = &rdAtig->ufpath[ tigs.ufpathIdx(fi) ];
    bool        rdAfwd   = (rdA->position.bgn < rdA->position.end);
    int32       rdAlo    = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
    int32       rdAhi    = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

    bool        isEnd    = (fi == 0) || (fi == fiLimit-1);

    uint32      ovlLen   = 0;
    BAToverlap *ovl      = OC->getOverlaps(rdA->ident, ovlLen);

    set<uint32> intersections;

    //if ((fi % 100) == 0)
    //  fprintf(stderr, "findBubbleReadPlacements()-- read %8u with %6u overlaps - %6.2f%% finished.\r",
    //          rdA->ident, ovlLen, 100.0 * fi / fiLimit);

    //  Compute all placements for this read.

    vector<overlapPlacement>   placements;

    placeReadUsingOverlaps(tigs, NULL, rdA->ident, placements, placeRead_noExtend);

    //  Weed out placements that aren't for bubbles, or that are for bubbles but are poor quality.  Or are to ourself!

    for (uint32 pi=0; pi<placements.size(); pi++) {
      uint32    rdBtigID = placements[pi].tigID;
      Unitig   *rdBtig   = tigs[rdBtigID];

      uint32    lo       = placements[pi].position.min();
      uint32    hi       = placements[pi].position.max();

      double    erate    = placements[pi].errors / placements[pi].aligned;

      //  Ignore the placement if it is to ourself.

      if (rdAtigID == rdBtigID)
        continue;

      //  Ignore the placement if it is to a non-tig or a singleton read.

      if ((rdBtigID == 0) ||
          (rdBtig   == NULL) ||
          (rdBtig->ufpath.size() == 1))
        continue;

      //  Ignore the placement if it is partial and not a terminal read.

      if ((isEnd == false) &&
          (placements[pi].fCoverage < 0.99)) {
        if (logFileFlagSet(LOG_BUBBLE_DETAIL))
          writeLog("tig %6u read %8u -> tig %6u %6u reads at %8u-%-8u (cov %7.5f erate %6.4f) - PARTIALLY PLACED\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Ignore the placement if it isn't to one of our bubble-popping candidate tigs.

      bool             dontcare = true;
      vector<uint32>  &pbubbles = potentialBubbles[rdAtigID];

      for (uint32 pb=0; pb<pbubbles.size(); pb++) {
        if (pbubbles[pb] == rdBtigID)
          dontcare = false;
      }

      if (dontcare) {
        if (logFileFlagSet(LOG_BUBBLE_DETAIL))
          writeLog("tig %6u read %8u -> tig %6u %6u reads at %8u-%-8u (cov %7.5f erate %6.4f) - NOT CANDIDATE TIG\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Ignore the placement if it is too diverged from the destination tig.

      if (rdBtig->overlapConsistentWithTig(deviationBubble, lo, hi, erate) < 0.5) {
        if (logFileFlagSet(LOG_BUBBLE_DETAIL))
          writeLog("tig %6u read %8u -> tig %6u %6u reads at %8u-%-8u (cov %7.5f erate %6.4f) - HIGH ERROR\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Good placement!

      if (logFileFlagSet(LOG_BUBBLE_DETAIL))
        writeLog("tig %6u read %8u -> tig %6u %6u reads at %8u-%-8u (cov %7.5f erate %6.4f)\n",
                 rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);

      placed[fi].push_back(placements[pi]);
    }
  }

  return(placed);
}






//  Bubble popping cannot be done in parallel -- there is a race condition when both tigs
//  A and B are considering merging in unitig C.

void
mergeOrphans(TigVector &tigs,
             double     deviationBubble) {

  BubTargetList   potentialBubbles;

  findPotentialBubbles(tigs, potentialBubbles);

  writeStatus("mergeOrphans()-- Found " F_SIZE_T " potential bubbles.\n", potentialBubbles.size());

  //if (potentialBubbles.size() == 0)
  //  return;

  writeLog("\n");
  writeLog("Found " F_SIZE_T " potential bubbles.\n", potentialBubbles.size());
  writeLog("\n");

  vector<overlapPlacement>   *placed = findBubbleReadPlacements(tigs, potentialBubbles, deviationBubble);

  //  We now have, in 'placed', a list of all the places that each read could be placed.  Decide if there is a _single_
  //  place for each bubble to be popped.

  uint32  tiLimit      = tigs.size();
  //uint32  tiNumThreads = omp_get_max_threads();
  //uint32  tiBlockSize  = (tiLimit < 100000 * tiNumThreads) ? tiNumThreads : tiLimit / 99999;

  //  Clear flags.
  for (uint32 ti=0; ti<tiLimit; ti++) {
    if (tigs[ti]) {
      tigs[ti]->_isBubble = false;
      tigs[ti]->_isRepeat = false;
    }
  }

  uint32  nUniqOrphan = 0;
  uint32  nReptOrphan = 0;
  uint32  nUniqBubble = 0;
  uint32  nReptBubble = 0;

  //  In parallel, process the placements.

  for (uint32 ti=0; ti<tiLimit; ti++) {
    if (potentialBubbles.count(ti) == 0)   //  Not a potential bubble
      continue;

    writeLog("\n");

    //  Save some interesting bits about our bubble.

    Unitig  *bubble        = tigs[ti];
    uint32   bubbleLen     = bubble->getLength();
    uint32   nReads        = bubble->ufpath.size();

    ufNode  &fRead         = bubble->ufpath.front();
    ufNode  &lRead         = bubble->ufpath.back();

    uint32   fReadID       = fRead.ident;  //  Ident of the first read
    uint32   lReadID       = lRead.ident;

    bool     bubbleInnie   = (fRead.position.isForward() && lRead.position.isReverse());
    bool     bubbleOuttie  = (fRead.position.isReverse() && lRead.position.isForward());
    bool     bubbleFwd     = (fRead.position.isForward() && lRead.position.isForward());
    bool     bubbleRev     = (fRead.position.isReverse() && lRead.position.isReverse());

    //  Scan the bubble, decide if there are _ANY_ read placements.  Log appropriately.

    bool     failedToPlaceAnchor = false;

    {
      char     placedS[128];

      char     placed0 = ((nReads > 0) && (placed[ bubble->ufpath[        0 ].ident ].size() > 0)) ? 't' : '-';
      char     placed1 = ((nReads > 1) && (placed[ bubble->ufpath[        1 ].ident ].size() > 0)) ? 't' : '-';
      char     placedb = ((nReads > 1) && (placed[ bubble->ufpath[ nReads-2 ].ident ].size() > 0)) ? 't' : '-';
      char     placeda = ((nReads > 0) && (placed[ bubble->ufpath[ nReads-1 ].ident ].size() > 0)) ? 't' : '-';

      uint32   placedN = 0;

      if (nReads > 3)
        for (uint32 fi=2; fi<nReads-2; fi++)
          if (placed[bubble->ufpath[fi].ident].size() > 0)
            placedN++;

      switch (nReads) {
        case 0:
          assert(0);
          break;

        case 1:
          snprintf(placedS, 128, "%c", placed0);
          break;

        case 2:
          snprintf(placedS, 128, "%c%c", placed0, placeda);
          break;

        case 3:
          snprintf(placedS, 128, "%c%c%c", placed0, placed1, placeda);
          break;

        case 4:
          snprintf(placedS, 128, "%c%c%c%c", placed0, placed1, placedb, placeda);
          break;

        default:
          snprintf(placedS, 128, "%c%c[%u]%c%c",
                  placed0, placed1, placedN, placedb, placeda);
          break;
      }

      failedToPlaceAnchor = ((placed0 != 't') || (placeda != 't'));

      writeLog("potential bubble tig %8u (reads %5u length %8u) - placed %s%s\n",
               bubble->id(), nReads, bubbleLen, placedS, failedToPlaceAnchor ? " FAILED" : "");
    }

    if (failedToPlaceAnchor)
      continue;


    //  Split the placements into piles for each target and build an interval list for each target.
    //  For each read in the tig, convert the vector of placements into interval lists, one list per target tig.

    map<uint32, intervalList<uint32> *>  targetIntervals;

    //  Add extended intervals for the first read.

    for (uint32 pp=0; pp<placed[fReadID].size(); pp++) {
      uint32  tid = placed[fReadID][pp].tigID;
      uint32  bgn = placed[fReadID][pp].position.min();

      if (targetIntervals[tid] == NULL)
        targetIntervals[tid] = new intervalList<uint32>;

      targetIntervals[tid]->add(bgn, bubbleLen);  //  Don't care if it goes off the high end of the tig.
    }

    //  Add extended intervals for the last read.

    for (uint32 pp=0; pp<placed[lReadID].size(); pp++) {
      uint32  tid = placed[lReadID][pp].tigID;
      uint32  end = placed[lReadID][pp].position.max();

      if (targetIntervals[tid] == NULL)
        targetIntervals[tid] = new intervalList<uint32>;

      if (end < bubbleLen)
        targetIntervals[tid]->add(0, end);  //  Careful!  Negative will underflow!
      else
        targetIntervals[tid]->add(end - bubbleLen, bubbleLen);
    }

    //  For each destination tig:
    //    merge the intervals
    //    for each interval
    //      find which bubble first/last reads map to each interval
    //      ignore if the extent of first/last is too big or small
    //      save otherwise

    vector<candidatePop *>    targets;

    for (map<uint32, intervalList<uint32> *>::iterator it=targetIntervals.begin(); it != targetIntervals.end(); ++it) {
      uint32                 targetID = it->first;
      intervalList<uint32>  *IL       = it->second;

      //  Merge.

      IL->merge();

      //  Figure out if each interval has both the first and last read of some bubble, and if those
      //  are properly sized.

      for (uint32 ii=0; ii<IL->numberOfIntervals(); ii++) {
        bool    noFirst = true;
        bool    noLast  = true;

        uint32  intBgn   = IL->lo(ii);
        uint32  intEnd   = IL->hi(ii);

        SeqInterval    fPos;
        SeqInterval    lPos;

        for (uint32 pp=0; pp<placed[fReadID].size(); pp++) {
          fPos = placed[fReadID][pp].position;

          if ((targetID == placed[fReadID][pp].tigID) &&
              (intBgn <= fPos.min()) && (fPos.max() <= intEnd)) {
            noFirst = false;
            break;
          }
        }

        for (uint32 pp=0; pp<placed[lReadID].size(); pp++) {
          lPos = placed[lReadID][pp].position;

          if ((targetID == placed[lReadID][pp].tigID) &&
              (intBgn <= lPos.min()) && (lPos.max() <= intEnd)) {
            noLast = false;
            break;
          }
        }

        //  Ignore if missing either read.

        if ((noFirst == true) ||
            (noLast  == true)) {
          writeLog("potential bubble tig %8u (length %8u) - target %8u %8u-%-8u (length %8u) - MISSING %s%s%s READ%s\n",
                   bubble->id(), bubble->getLength(),
                   targetID, intBgn, intEnd, intEnd - intBgn,
                   (noFirst) ? "FIRST" : "",
                   (noFirst && noLast) ? " and " : "",
                   (noLast)  ? "LAST"  : "",
                   (noFirst && noLast) ? "S" : "");
          continue;
        }

        writeLog("potential bubble tig %8u (length %8u) - target %8u %8u-%-8u (length %8u) - %8u-%-8u %8u-%-8u\n",
                 bubble->id(), bubble->getLength(),
                 targetID, intBgn, intEnd, intEnd - intBgn,
                 fPos.min(), fPos.max(),
                 lPos.min(), lPos.max());


        //  Ignore if the reads align in inconsistent orientations.

#if 0
        bool  alignFwd   = (fPos.min() < lPos.max()) ? true : false;
        bool  fPosFwd    = fPos.isForward();
        bool  lPosFwd    = lPos.isForward();

        bool  alignInnie  = (alignFwd == true) ? ((fPosFwd == true) && (lPosFwd == false)) : ((fPosFwd == false) && (lPosFwd == true));
        bool  alignOuttie = false;
        bool  alignFwd    = false;
        bool  alignRev    = false;

        bool  alignInnie = (alignFwd  && fPosFwd && !rPosFwd);


        //if ((bubbleInnie  == true) &&
        //if ((bubbleOuttie == true) && ((alignFwd == true) || (fPosFwd == true) || (rPosFwd == false)));
        //if ((bubbleFwd    == true) && ((alignFwd == true) || (fPosFwd == true) || (rPosFwd == false)));
        //if ((bubbleRev    == true) && ((alignFwd == true) || (fPosFwd == true) || (rPosFwd == false)));
#endif

        //  Ignore if the region is too small or too big.

        uint32  regionMin = min(fPos.min(), lPos.min());
        uint32  regionMax = max(fPos.max(), lPos.max());

        if ((regionMax - regionMin < 0.75 * bubbleLen) ||
            (regionMax - regionMin > 1.25 * bubbleLen))
          continue;

        //  Both reads placed, and at about the right size.  We probably should be checking orientation.  Maybe tomorrow.

        targets.push_back(new candidatePop(bubble, tigs[targetID], regionMin, regionMax));
      }  //  Over all intervals for this target
    }  //  Over all targets

    //  Done with the targetIntervals.  Clean up.

    for (map<uint32, intervalList<uint32> *>::iterator it=targetIntervals.begin(); it != targetIntervals.end(); ++it)
      delete it->second;

    targetIntervals.clear();

    //  If no targets, nothing to do.

    if (targets.size() == 0) {
      writeLog("potential bubble tig %8u - generated no targets\n", ti);
      continue;
    }

    //  Run through the placements again, and assign them to the correct target.
    //
    //  For each read:
    //  For each acceptable placement:
    //  For each target location:
    //  If the placement is for this target, save it.

    for (uint32 fi=0; fi<nReads; fi++) {
      uint32  readID  = bubble->ufpath[fi].ident;

      for (uint32 pp=0; pp<placed[readID].size(); pp++) {
        uint32  tid = placed[readID][pp].tigID;

        uint32  bgn = placed[readID][pp].position.min();
        uint32  end = placed[readID][pp].position.max();

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
            writeLog("duplicate read alignment for tig %u read %u - better %u-%-u %.4f - worse %u-%-u %.4f\n",
                     t->placed[aa].tigID, t->placed[aa].frgID,
                     t->placed[aa].position.bgn, t->placed[aa].position.end, t->placed[aa].errors / t->placed[aa].aligned,
                     t->placed[bb].position.bgn, t->placed[bb].position.end, t->placed[bb].errors / t->placed[bb].aligned);
#endif
            t->placed[bb] = overlapPlacement();
          } else {
#ifdef SHOW_MULTIPLE_PLACEMENTS
            writeLog("duplicate read alignment for tig %u read %u - better %u-%-u %.4f - worse %u-%-u %.4f\n",
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

    //  Make a set of the reads in the bubble.

    for (uint32 fi=0; fi<nReads; fi++)
      tigReads.insert(bubble->ufpath[fi].ident);

    //  Compare the bubble against each target.

    uint32   nOrphan      = 0;   //  Full coverage; bubble can be popped.
    uint32   orphanTarget = 0;

    uint32   nBubble      = 0;   //  Partial coverage, bubble cannot be popped.
    uint32   bubbleTarget = 0;

    for (uint32 tt=0; tt<targets.size(); tt++) {
      tgtReads.clear();

      for (uint32 op=0; op<targets[tt]->placed.size(); op++) {
        if (logFileFlagSet(LOG_BUBBLE_DETAIL))
          writeLog("tig %8u length %9u -> target %8u piece %2u position %9u-%-9u length %8u - read %7u at %9u-%-9u\n",
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

      for (uint32 fi=0; fi<nReads; fi++)
        if (tgtReads.count(bubble->ufpath[fi].ident) > 0)
          n5++;
        else
          break;

      for (uint32 fi=nReads; fi-->0; )
        if (tgtReads.count(bubble->ufpath[fi].ident) > 0)
          n3++;
        else
          break;


      for (uint32 fi=0; fi<nReads; fi++)
        if (tgtReads.count(bubble->ufpath[fi].ident) > 0)
          nt++;


      //  Report now, before we nuke targets[tt] for being not a bubble!

      if ((nt == nReads) ||
          ((n5 > 0) && (n3 > 0)))
        writeLog("tig %8u length %9u -> target %8u piece %2u position %9u-%-9u length %8u - expected %3" F_SIZE_TP " reads, had %3" F_SIZE_TP " reads.  n5=%3u n3=%3u nt=%3u\n",
                 bubble->id(), bubble->getLength(),
                 targets[tt]->target->id(), tt, targets[tt]->bgn, targets[tt]->end, targets[tt]->end - targets[tt]->bgn,
                 tigReads.size(),
                 tgtReads.size(), n5, n3, nt);

      //  Decide if this is a bubble, orphan from construction, or repeat.

      if (nt == nReads) {
        nOrphan++;
        orphanTarget = tt;
      }

      else if ((n5 > 0) && (n3 > 0)) {
        nBubble++;
        bubbleTarget = tt;
      }
    }

    //  If no placements, pbbbt, not a whole lot we can do here.  Leave it as is.  It's not even
    //  worth logging (there are many of these).

    if (nOrphan + nBubble == 0) {
    }

    //  If not an orphan, mark it as a bubble.  If multiple bubble placements, mark it as a repeat
    //  so we can use it in repeat detection.
    //
    //  If there are orphan placements also, those placements are superior to the bubble placements,
    //  and we'll place the orphan.

    else if (nOrphan == 0) {
      if (nBubble == 1) {
        nUniqBubble++;
        writeStatus("mergeOrphans()-- tig %8u BUBBLE -> tig %8u\n",
                    bubble->id(),
                    targets[bubbleTarget]->target->id());
      } else {
        nReptBubble++;
        writeStatus("mergeOrphans()-- tig %8u BUBBLE -> repeat\n",
                    bubble->id());
      }

      writeLog("tig %8u length %8u reads %6u - %s.\n",
               bubble->id(), bubble->getLength(), nReads,
               (nBubble == 1) ? "bubble" : "bubble-repeat");
      writeLog("\n");

      bubble->_isRepeat = (nBubble > 1);
      bubble->_isBubble = true;
    }

    //  If a unique orphan placement, place it there.

    else if (nOrphan == 1) {
      nUniqOrphan++;
      writeStatus("mergeOrphans()-- tig %8u ORPHAN -> tig %8u\n",
                  bubble->id(),
                  targets[bubbleTarget]->target->id());

      writeLog("tig %8u length %8u reads %6u - orphan\n", bubble->id(), bubble->getLength(), nReads);

      for (uint32 op=0, tt=orphanTarget; op<targets[tt]->placed.size(); op++) {
        ufNode  frg;

        frg.ident        = targets[tt]->placed[op].frgID;
        frg.contained    = 0;
        frg.parent       = 0;
        frg.ahang        = 0;
        frg.bhang        = 0;
        frg.position.bgn = targets[tt]->placed[op].position.bgn;
        frg.position.end = targets[tt]->placed[op].position.end;

        writeLog("move read %u from tig %u to tig %u %u-%-u\n",
                 frg.ident,
                 bubble->id(),
                 targets[tt]->target->id(), frg.position.bgn, frg.position.end);

        targets[tt]->target->addRead(frg, 0, false);
      }

      writeLog("\n");

      tigs[bubble->id()] = NULL;
      delete bubble;
    }

    //  Otherwise, there are multiple orphan placements.  We can't distinguish between them, and
    //  instead just place reads where they individually decide to go.

    else {
      nReptBubble++;
      writeStatus("mergeOrphans()-- tig %8u ORPHAN -> multiple tigs\n",
                  bubble->id(),
                  targets[bubbleTarget]->target->id());

      writeLog("tig %8u length %8u reads %6u - orphan with multiple placements\n", bubble->id(), bubble->getLength(), nReads);

      for (uint32 fi=0; fi<nReads; fi++) {
        uint32  rr = bubble->ufpath[fi].ident;
        double  er = 1.00;
        uint32  bb = 0;

        for (uint32 pp=0; pp<placed[rr].size(); pp++) {
          double erate = placed[rr][pp].errors / placed[rr][pp].aligned;

          if ((erate < er) &&                             //  Reads placed in 'bubble' by this same method
              (placed[rr][pp].tigID != bubble->id())) {   //  will pick 'bubble' (which is about the be removed)
            er = erate;                                   //  as their best location, so careful to exclude
            bb = pp;                                      //  all self placements!
          }
        }

        ufNode  frg;

        frg.ident        = placed[rr][bb].frgID;
        frg.contained    = 0;
        frg.parent       = 0;
        frg.ahang        = 0;
        frg.bhang        = 0;
        frg.position.bgn = placed[rr][bb].position.bgn;
        frg.position.end = placed[rr][bb].position.end;

        Unitig  *target  = tigs[placed[rr][bb].tigID];

        writeLog("move read %u from tig %u to tig %u %u-%-u\n",
                 frg.ident,
                 bubble->id(),
                 target->id(), frg.position.bgn, frg.position.end);

        assert(target->id() != bubble->id());

        target->addRead(frg, 0, false);
      }

      writeLog("\n");

      tigs[bubble->id()] = NULL;
      delete bubble;
    }

    //  Clean up the targets list.

    for (uint32 tt=0; tt<targets.size(); tt++) {
      delete targets[tt];
      targets[tt] = NULL;
    }

    targets.clear();

  }  //  Over all bubbles

  writeLog("\n");   //  Needed if no bubbles are popped.

  writeStatus("mergeOrphans()-- placed    %5u unique orphan tigs\n", nUniqOrphan);
  writeStatus("mergeOrphans()-- shattered %5u repeat orphan tigs\n", nReptOrphan);
  writeStatus("mergeOrphans()-- marked    %5u unique bubble tigs\n", nUniqBubble);
  writeStatus("mergeOrphans()-- marked    %5u repeat bubble tigs\n", nReptBubble);

  delete [] placed;

  //  Sort reads in all the tigs.  Overkill, but correct.

  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = tigs[ti];

    if ((tig == NULL) ||               //  Not a tig, ignore it.
        (tig->ufpath.size() == 1))     //  Singleton, already sorted.
      continue;

    tig->sort();
  }
}


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
 *    Brian P. Walenz beginning on 2016-MAR-11
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

#undef  SHOW_MULTIPLE_PLACEMENTS  //  Reports reads that are placed multiple times in a single target region

class candidatePop {
public:
  candidatePop(Unitig *orphan_, Unitig *target_, uint32 bgn_, uint32 end_) {
    orphan = orphan_;
    target = target_;
    bgn    = bgn_;
    end    = end_;
  };

  Unitig  *orphan;
  Unitig  *target;

  uint32   bgn;
  uint32   end;

  vector<overlapPlacement>  placed;
};


//  A list of the target tigs that a orphan could be popped into.
typedef  map<uint32, vector<uint32> >  BubTargetList;



//  Decide which tigs can be orphans.  Any unitig where (nearly) every dovetail
//  read has an overlap to some other unitig is a candidate for orphan popping.
//
//  Counts the number of reads that have an overlap to some other tig
//  (tigOlapsTo).  if more than half the reads in the tig have an overlap to
//  some other tig, it is a potential place to pop the bubble.
//
//  Returns BubTargetList, a map of uint32 to vector<uint32>, of the potential
//  places that some tig could be popped into.

void
findPotentialOrphans(TigVector       &tigs,
                     BubTargetList   &potentialOrphans) {

  writeStatus("\n");
  writeStatus("findPotentialOrphans()-- working on " F_U32 " tigs.\n", tigs.size());

  for (uint32 ti=0; ti<tigs.size(); ti++) {
    Unitig  *tig = tigs[ti];

    if ((tig == NULL) ||               //  Not a tig, ignore it.
        (tig->ufpath.size() == 1))     //  Singleton, handled elsewhere.
      continue;

    //  Count the number of reads that have an overlap to some other tig.  tigOlapsTo[otherTig] = count.

    map<uint32,uint32>  tigOlapsTo;
    uint32              nonContainedReads    = 0;
    bool                validOrphan          = true;

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      uint32      rid = tig->ufpath[fi].ident;

      if (OG->isContained(rid) == true)  //  Don't need to check contained reads.  If their container
        continue;                        //  passes the tests below, the contained read will too.

      nonContainedReads++;

      //  Find the list of tigs that we have an overlap to.

      set<uint32>  readOlapsTo;

      uint32      ovlLen   = 0;
      BAToverlap *ovl      = OC->getOverlaps(rid, ovlLen);

      for (uint32 oi=0; oi<ovlLen; oi++) {
        uint32  ovlTigID = tigs.inUnitig(ovl[oi].b_iid);
        Unitig *ovlTig   = tigs[ovlTigID];

        if ((ovlTigID == 0) ||                          //  Skip this overlap if it is to an unplaced read,
            (ovlTig == NULL) ||                         //  to a singleton tig, to ourself, or to a unitig
            (ovlTig->ufpath.size() == 1) ||             //  that is shorter than us.  We can not pop this
            (ovlTig->id() == tig->id()) ||              //  tig as a orphan in any of those cases.
            (ovlTig->getLength() < tig->getLength()))   //
          continue;

        readOlapsTo.insert(ovlTigID);                   //  Otherwise, remember that we had an overlap to ovlTig.
      }

      //  With the list of tigs that this read has an overlap to, add one to each tig in the list of
      //  tigs that this tig has an overlap to.

      for (set<uint32>::iterator it=readOlapsTo.begin(); it != readOlapsTo.end(); ++it)
        tigOlapsTo[*it]++;

      //  Decide if we're a valid potential orphan.  If tig id (in it->first) has overlaps to
      //  (nearly) every read we've seen so far (nonContainedReads), we're still a valid orphan.

      validOrphan = false;

      for (map<uint32,uint32>::iterator it=tigOlapsTo.begin(); it != tigOlapsTo.end(); ++it)
        if (it->second == nonContainedReads)         //  All reads have an overlap to the tig
          validOrphan = true;                        //  at *it, so valid orphan.

      if (validOrphan == false)                      //  If not a valid orphan, bail.  There is no other
        break;                                       //  tig that all of our reads have overlaps to.
    }

    //  If not a valid orphan, just move on to the next tig.

    if (validOrphan == false)
      continue;

    //  Otherwise, a valid orphan!  There is at least one tig that (nearly) every dovetail read has
    //  at least one overlap to.  Save those tigs in potentialOrphans.

    uint32  nTigs = 0;

    for (map<uint32,uint32>::iterator it=tigOlapsTo.begin(); it != tigOlapsTo.end(); ++it)
      if (it->second >= 0.5 * nonContainedReads)
        nTigs++;

    writeLog("findPotentialOrphans()--\n");
    writeLog("findPotentialOrphans()-- potential orphan tig %8u length %9u nReads %7u to %3u tigs:\n",
             tig->id(), tig->getLength(), tig->ufpath.size(), nTigs);

    for (map<uint32,uint32>::iterator it=tigOlapsTo.begin(); it != tigOlapsTo.end(); ++it) {
      if (it->second >= 0.5 * nonContainedReads) {
        Unitig  *dest = tigs[it->first];

        writeLog("findPotentialOrphans()--                  tig %8u length %9u nReads %7u\n", dest->id(), dest->getLength(), dest->ufpath.size());

        potentialOrphans[ti].push_back(dest->id());
      }
    }
  }  //  Over all tigs.

  flushLog();
}



//  Find filtered placements for all the reads in the potential orphan tigs.

vector<overlapPlacement>  *
findOrphanReadPlacements(TigVector       &tigs,
                         BubTargetList   &potentialOrphans,
                         double           deviationOrphan) {
  uint32  fiLimit      = RI->numReads();
  uint32  fiNumThreads = omp_get_max_threads();
  uint32  fiBlockSize  = (fiLimit < 1000 * fiNumThreads) ? fiNumThreads : fiLimit / 999;

  uint64  nReads       = 0;
  uint64  nPlaces      = 0;

  vector<overlapPlacement>   *placed = new vector<overlapPlacement> [fiLimit + 1];

  writeLog("findOrphanReadPlacement()--\n");

#pragma omp parallel for schedule(dynamic, fiBlockSize)
  for (uint32 fi=0; fi<fiLimit; fi++) {
    uint32     rdAtigID = tigs.inUnitig(fi);

    if ((rdAtigID == 0) ||                           //  Read not placed in a tig, ignore it.
        (OG->isContained(fi)) ||                     //  Read is contained, ignore it.
        (potentialOrphans.count(rdAtigID) == 0))     //  Read isn't in a potential orphan, ignore it.
      continue;

#pragma omp atomic
    nReads++;

    Unitig     *rdAtig   = tigs[rdAtigID];
    ufNode     *rdA      = &rdAtig->ufpath[ tigs.ufpathIdx(fi) ];
    bool        rdAfwd   = (rdA->position.bgn < rdA->position.end);
    int32       rdAlo    = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
    int32       rdAhi    = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

    bool        isEnd    = (rdAlo == 0) || (rdAhi == rdAtig->getLength());

    //  Compute all placements for this read.  We ask for only fully placed reads.

    vector<overlapPlacement>   placements;

    placeReadUsingOverlaps(tigs, NULL, rdA->ident, placements, placeRead_fullMatch);

    //  Weed out placements that aren't for orphans, or that are for orphans but are poor quality.  Or are to ourself!

    for (uint32 pi=0; pi<placements.size(); pi++) {
      uint32    rdBtigID = placements[pi].tigID;
      Unitig   *rdBtig   = tigs[rdBtigID];

      uint32    lo       = placements[pi].position.min();
      uint32    hi       = placements[pi].position.max();

      double    erate    = placements[pi].errors / placements[pi].aligned;

      if ((rdAtigID == rdBtigID) ||                    //  To ourself.
          (rdBtigID == 0) ||                           //  To a singleton read.
          (rdBtig   == NULL) ||                        //  To a singleton read.
          (rdBtig->ufpath.size() == 1) ||              //  To a singleton tig.
          (potentialOrphans.count((rdBtigID) > 0)))    //  To a potential orphan tig
        continue;

      //  Ignore the placement if it isn't to one of our orphan-popping candidate tigs.

      bool             dontcare = true;
      vector<uint32>  &porphans = potentialOrphans[rdAtigID];

      for (uint32 pb=0; pb<porphans.size(); pb++)
        if (porphans[pb] == rdBtigID)
          dontcare = false;

      if (dontcare) {
        if (logFileFlagSet(LOG_ORPHAN_DETAIL))
          writeLog("findOrphanReadPlacement()-- tig %6u read %8u -> tig %6u %6u reads at %8u-%-8u (cov %7.5f erate %6.4f) - NOT CANDIDATE TIG\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Ignore the placement if it is to a potential orphan.

      if (potentialOrphans.count(rdBtigID) > 0) {
        if (logFileFlagSet(LOG_ORPHAN_DETAIL))
          writeLog("findOrphanReadPlacement()-- tig %6u read %8u -> tig %6u %6u reads at %8u-%-8u (cov %7.5f erate %6.4f) - INTO POTENTIAL ORPHAN\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Ignore the placement if it is too diverged from the destination tig.

      if (rdBtig->overlapConsistentWithTig(deviationOrphan, lo, hi, erate) < 0.5) {
        if (logFileFlagSet(LOG_ORPHAN_DETAIL))
          writeLog("findOrphanReadPlacement()-- tig %6u read %8u -> tig %6u %6u reads at %8u-%-8u (cov %7.5f erate %6.4f) - HIGH ERROR\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Good placement!

      if (logFileFlagSet(LOG_ORPHAN_DETAIL))
        writeLog("findOrphanReadPlacement()-- tig %6u read %8u -> tig %6u %6u reads at %8u-%-8u (cov %7.5f erate %6.4f)\n",
                 rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);

#pragma omp atomic
      nPlaces++;

      placed[fi].push_back(placements[pi]);
    }
  }

  writeLog("findOrphanReadPlacement()--  placed %u reads into %u locations\n", nReads, nPlaces);

  return(placed);
}




static
bool
failedToPlaceAnchor(Unitig                     *orphan,
                    vector<overlapPlacement>   *placed) {
  uint32   nReads  = orphan->ufpath.size();

  char     placed0 = ((nReads > 0) && (placed[ orphan->ufpath[        0 ].ident ].size() > 0)) ? 't' : '-';
  char     placed1 = ((nReads > 1) && (placed[ orphan->ufpath[        1 ].ident ].size() > 0)) ? 't' : '-';
  char     placedb = ((nReads > 1) && (placed[ orphan->ufpath[ nReads-2 ].ident ].size() > 0)) ? 't' : '-';
  char     placeda = ((nReads > 0) && (placed[ orphan->ufpath[ nReads-1 ].ident ].size() > 0)) ? 't' : '-';

  char     placedS[128];
  uint32   placedN = 0;

  bool     failed  = false;

  if (nReads > 3)
    for (uint32 fi=2; fi<nReads-2; fi++)
      if (placed[orphan->ufpath[fi].ident].size() > 0)
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

  failed = ((placed0 != 't') || (placeda != 't'));

  writeLog("failedToPlaceAnchor()-- potential orphan tig %8u (reads %5u length %8u) - placed %s%s\n",
           orphan->id(), nReads, orphan->getLength(), placedS, failed ? " FAILED" : "");

  return(failed);
}




static
void
addInitialIntervals(Unitig                               *orphan,
                    vector<overlapPlacement>             *placed,
                    uint32                                fReadID,
                    uint32                                lReadID,
                    map<uint32, intervalList<uint32> *>  &targetIntervals) {
  uint32   orphanLen  = orphan->getLength();

  //  Add extended intervals for the first read.
  //
  //    target ---------------------------------------------
  //    read        -------
  //    orphan      -------------------------

  for (uint32 pp=0; pp<placed[fReadID].size(); pp++) {
    uint32  tid = placed[fReadID][pp].tigID;
    uint32  bgn = placed[fReadID][pp].position.min();

    if (targetIntervals[tid] == NULL)
      targetIntervals[tid] = new intervalList<uint32>;

    targetIntervals[tid]->add(bgn, orphanLen);  //  Don't care if it goes off the high end of the tig.
  }

  //  Add extended intervals for the last read.
  //
  //    target ---------------------------------------------
  //    read                          -------
  //    orphan      -------------------------

  for (uint32 pp=0; pp<placed[lReadID].size(); pp++) {
    uint32  tid = placed[lReadID][pp].tigID;
    uint32  end = placed[lReadID][pp].position.max();

    if (targetIntervals[tid] == NULL)
      targetIntervals[tid] = new intervalList<uint32>;

    if (end < orphanLen)
      targetIntervals[tid]->add(0, end);  //  Careful!  Negative will underflow!
    else
      targetIntervals[tid]->add(end - orphanLen, orphanLen);
  }
}



static
void
saveCorrectlySizedInitialIntervals(Unitig                    *orphan,
                                   Unitig                    *target,
                                   intervalList<uint32>      *IL,
                                   uint32                     fReadID,
                                   uint32                     lReadID,
                                   vector<overlapPlacement>  *placed,
                                   vector<candidatePop *>    &targets) {

  IL->merge();  // Merge overlapping initial intervals created above.

  for (uint32 ii=0; ii<IL->numberOfIntervals(); ii++) {
    bool    noFirst = true;
    bool    noLast  = true;

    uint32  intBgn   = IL->lo(ii);
    uint32  intEnd   = IL->hi(ii);

    SeqInterval    fPos;
    SeqInterval    lPos;

    //  Find the read placement in this interval, if it exists.

    for (uint32 pp=0; pp<placed[fReadID].size(); pp++) {
      fPos = placed[fReadID][pp].position;

      if ((target->id() == placed[fReadID][pp].tigID) &&
          (intBgn <= fPos.min()) && (fPos.max() <= intEnd)) {
        noFirst = false;
        break;
      }
    }

    for (uint32 pp=0; pp<placed[lReadID].size(); pp++) {
      lPos = placed[lReadID][pp].position;

      if ((target->id() == placed[lReadID][pp].tigID) &&
          (intBgn <= lPos.min()) && (lPos.max() <= intEnd)) {
        noLast = false;
        break;
      }
    }

    //  Ignore if missing either read.

    if ((noFirst == true) ||
        (noLast  == true)) {
      writeLog("saveCorrectlySizedInitialIntervals()-- potential orphan tig %8u (length %8u) - target %8u %8u-%-8u (length %8u) - MISSING %s%s%s READ%s\n",
               orphan->id(), orphan->getLength(),
               target->id(), intBgn, intEnd, intEnd - intBgn,
               (noFirst) ? "FIRST" : "",
               (noFirst && noLast) ? " and " : "",
               (noLast)  ? "LAST"  : "",
               (noFirst && noLast) ? "S" : "");
      continue;
    }

    writeLog("saveCorrectlySizedInitialIntervals()-- potential orphan tig %8u (length %8u) - target %8u %8u-%-8u (length %8u) - %8u-%-8u %8u-%-8u\n",
             orphan->id(), orphan->getLength(),
             target->id(), intBgn, intEnd, intEnd - intBgn,
             fPos.min(), fPos.max(),
             lPos.min(), lPos.max());

    //  Ignore if the region is too small or too big.

    uint32  regionMin = min(fPos.min(), lPos.min());
    uint32  regionMax = max(fPos.max(), lPos.max());

    if ((regionMax - regionMin < 0.75 * orphan->getLength()) ||
        (regionMax - regionMin > 1.25 * orphan->getLength()))
      continue;

    //  We probably should be checking orientation.  Maybe tomorrow.

    //  Both reads placed, and at about the right size.  Save the candidate position - we can
    //  possibly place 'orphan' in 'tigs[target->id()' at position regionMin-regionMax.

    targets.push_back(new candidatePop(orphan, target, regionMin, regionMax));
  }  //  Over all intervals for this target

  //  We're done with this intervalList, clean it up.  This does leave a dangling pointer in the map<> though.

  delete IL;
}





void
assignReadsToTargets(Unitig                     *orphan,
                     vector<overlapPlacement>   *placed,
                     vector<candidatePop *>     targets) {

  for (uint32 fi=0; fi<orphan->ufpath.size(); fi++) {
    uint32  readID  = orphan->ufpath[fi].ident;

    for (uint32 pp=0; pp<placed[readID].size(); pp++) {
      uint32  tid = placed[readID][pp].tigID;
      uint32  bgn = placed[readID][pp].position.min();
      uint32  end = placed[readID][pp].position.max();

      for (uint32 tt=0; tt<targets.size(); tt++)                           //  For a read placed in tig 'tid' at 'bgn-end',
        if ((targets[tt]->target->id() == tid) &&                          //  if the target is the same tig and the read
            (isContained(bgn, end, targets[tt]->bgn, targets[tt]->end)))   //  is contained in the target position,
          targets[tt]->placed.push_back(placed[readID][pp]);               //    save the position to the target
    }
  }

  //  Remove duplicate placements from each target.
  //
  //  Detect duplicates, keep the one with lower error.
  //  There are a lot of duplicate placements, logging isn't terribly useful.

  uint32  nDup = 0;
  uint32  save;
  uint32  remo;

  for (uint32 tt=0; tt<targets.size(); tt++) {
    candidatePop *t = targets[tt];

    for (uint32 aa=0; aa<t->placed.size(); aa++) {
      for (uint32 bb=0; bb<t->placed.size(); bb++) {
        if ((aa == bb) ||
            (t->placed[aa].frgID != t->placed[bb].frgID) ||
            (t->placed[aa].frgID == 0) ||
            (t->placed[bb].frgID == 0))
          continue;

        nDup++;

        if (t->placed[aa].errors / t->placed[aa].aligned < t->placed[bb].errors / t->placed[bb].aligned) {
          save = aa;
          remo = bb;
        } else {
          save = bb;
          remo = aa;
        }

#ifdef SHOW_MULTIPLE_PLACEMENTS
        writeLog("assignReadsToTargets()-- duplicate read alignment for tig %u read %u - better %u-%-u %.4f - worse %u-%-u %.4f\n",
                 t->placed[save].tigID, t->placed[save].frgID,
                 t->placed[save].position.bgn, t->placed[save].position.end, t->placed[save].errors / t->placed[save].aligned,
                 t->placed[remo].position.bgn, t->placed[remo].position.end, t->placed[remo].errors / t->placed[remo].aligned);
#endif

        t->placed[remo] = overlapPlacement();
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

  writeLog("assignReadsToTargets()-- Removed %u duplicate placements.\n", nDup);
}






void
mergeOrphans(TigVector &tigs,
             double     deviationOrphan) {

  //  Find, for each tig, the list of other tigs that it could potentially be placed into.

  BubTargetList   potentialOrphans;

  findPotentialOrphans(tigs, potentialOrphans);

  writeStatus("mergeOrphans()-- Found " F_SIZE_T " potential orphans.\n", potentialOrphans.size());

  writeLog("\n");
  writeLog("mergeOrphans()-- Found " F_SIZE_T " potential orphans.\n", potentialOrphans.size());
  writeLog("\n");

  //  For any tig that is a potential orphan, find all read placements.

  vector<overlapPlacement>   *placed = findOrphanReadPlacements(tigs, potentialOrphans, deviationOrphan);

  //  We now have, in 'placed', a list of all the places that each read could be placed.  Decide if there is a _single_
  //  place for each orphan to be popped.

  uint32  nUniqOrphan = 0;
  uint32  nReptOrphan = 0;

  for (uint32 ti=0; ti<tigs.size(); ti++) {
    Unitig  *orphan        = tigs[ti];

    if (potentialOrphans.count(ti) == 0)
      continue;

    //  Scan the orphan, decide if there are _ANY_ read placements.  Log appropriately.

    if (failedToPlaceAnchor(orphan, placed) == true)
      continue;

    writeLog("mergeOrphans()-- Processing orphan %u - %u bp %u reads\n", ti, orphan->getLength(), orphan->ufpath.size());

    //  Create intervals for each placed read.
    //
    //    target ---------------------------------------------
    //    read        -------
    //    orphan      -------------------------

    uint32                                fReadID = orphan->ufpath.front().ident;
    uint32                                lReadID = orphan->ufpath.back().ident;
    map<uint32, intervalList<uint32> *>   targetIntervals;

    addInitialIntervals(orphan, placed, fReadID, lReadID, targetIntervals);

    //  Figure out if each interval has both the first and last read of some orphan, and if those
    //  are properly sized.  If so, save a candidatePop.

    vector<candidatePop *>    targets;

    for (map<uint32, intervalList<uint32> *>::iterator it=targetIntervals.begin(); it != targetIntervals.end(); ++it)
      if (tigs[it->first] == NULL)
        writeLog("mergeOrphans()-- orphan %u wants to go into nonexistent tig %u!\n", ti, it->first);
      else
        saveCorrectlySizedInitialIntervals(orphan,
                                           tigs[it->first],     //  The targetID      in targetIntervals
                                           it->second,          //  The interval list in targetIntervals
                                           fReadID,
                                           lReadID,
                                           placed,
                                           targets);

    targetIntervals.clear();   //  intervalList already freed.

    //  If no targets, nothing to do.

    writeLog("mergeOrphans()-- Processing orphan %u - found %u target location%s\n", ti, targets.size(), (targets.size() == 1) ? "" : "s");

    if (targets.size() == 0)
      continue;

    //  Assign read placements to targets.

    assignReadsToTargets(orphan, placed, targets);

    //  Compare the orphan against each target.

    uint32   nOrphan      = 0;   //  Number of targets that have all the reads.
    uint32   orphanTarget = 0;   //  If nOrphan == 1, the target we're popping into.

    for (uint32 tt=0; tt<targets.size(); tt++) {
      uint32  orphanSize = orphan->ufpath.size();
      uint32  targetSize = targets[tt]->placed.size();

      //  Report now, before we nuke targets[tt] for being not a orphan!

      if (logFileFlagSet(LOG_ORPHAN_DETAIL))
        for (uint32 op=0; op<targets[tt]->placed.size(); op++)
          writeLog("mergeOrphans()-- tig %8u length %9u -> target %8u piece %2u position %9u-%-9u length %8u - read %7u at %9u-%-9u\n",
                   orphan->id(), orphan->getLength(),
                   targets[tt]->target->id(), tt, targets[tt]->bgn, targets[tt]->end, targets[tt]->end - targets[tt]->bgn,
                   targets[tt]->placed[op].frgID,
                   targets[tt]->placed[op].position.bgn, targets[tt]->placed[op].position.end);

      writeLog("mergeOrphans()-- tig %8u length %9u -> target %8u piece %2u position %9u-%-9u length %8u - expected %3u reads, had %3u reads.\n",
               orphan->id(), orphan->getLength(),
               targets[tt]->target->id(), tt, targets[tt]->bgn, targets[tt]->end, targets[tt]->end - targets[tt]->bgn,
               orphanSize, targetSize);

      //  If all reads placed, we can merge this orphan into the target.  Preview: if this happens more than once, we just
      //  split the orphan and place reads individually.

      if (orphanSize == targetSize) {
        nOrphan++;
        orphanTarget = tt;
      }
    }

    //  If a unique orphan placement, place it there.

    if (nOrphan == 1) {
      writeLog("mergeOrphans()-- tig %8u length %8u reads %6u - orphan\n", orphan->id(), orphan->getLength(), orphan->ufpath.size());
      nUniqOrphan++;

      for (uint32 op=0, tt=orphanTarget; op<targets[tt]->placed.size(); op++) {
        ufNode  frg;

        frg.ident        = targets[tt]->placed[op].frgID;
        frg.contained    = 0;
        frg.parent       = 0;
        frg.ahang        = 0;
        frg.bhang        = 0;
        frg.position.bgn = targets[tt]->placed[op].position.bgn;
        frg.position.end = targets[tt]->placed[op].position.end;

        writeLog("mergeOrphans()-- move read %u from tig %u to tig %u %u-%-u\n",
                 frg.ident,
                 orphan->id(),
                 targets[tt]->target->id(), frg.position.bgn, frg.position.end);

        targets[tt]->target->addRead(frg, 0, false);
      }

      writeLog("\n");

      tigs[orphan->id()] = NULL;
      delete orphan;
    }

    //  If multiply placed, we can't distinguish between them, and
    //  instead just place reads where they individually decide to go.

    if (nOrphan > 1) {
      writeLog("tig %8u length %8u reads %6u - orphan with multiple placements\n", orphan->id(), orphan->getLength(), orphan->ufpath.size());
      nReptOrphan++;

      for (uint32 fi=0; fi<orphan->ufpath.size(); fi++) {
        uint32  rr = orphan->ufpath[fi].ident;
        double  er = 1.00;
        uint32  bb = 0;

        //  Over all placements for this read, pick the one with lowest error, as long as it isn't
        //  to the orphan.

        for (uint32 pp=0; pp<placed[rr].size(); pp++) {
          double erate = placed[rr][pp].errors / placed[rr][pp].aligned;

          if ((er < erate) ||                           //  Worse placement.
              (placed[rr][pp].tigID == orphan->id()))   //  Self placement.
            continue;

          er = erate;
          bb = pp;
        }

        assert(rr == placed[rr][bb].frgID);

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
                 orphan->id(),
                 placed[rr][bb].tigID, frg.position.bgn, frg.position.end);

        assert(target       != NULL);
        assert(target->id() != orphan->id());

        target->addRead(frg, 0, false);
      }

      writeLog("\n");

      tigs[orphan->id()] = NULL;
      delete orphan;
    }

    //  Clean up the targets list.

    for (uint32 tt=0; tt<targets.size(); tt++) {
      delete targets[tt];
      targets[tt] = NULL;
    }

    targets.clear();

  }  //  Over all orphans

  writeLog("\n");   //  Needed if no orphans are popped.

  writeStatus("mergeOrphans()-- placed    %5u unique orphan tigs\n", nUniqOrphan);
  writeStatus("mergeOrphans()-- shattered %5u repeat orphan tigs\n", nReptOrphan);
  writeStatus("mergeOrphans()--\n");

  delete [] placed;

  //  Sort reads in all the tigs.  Overkill, but correct.

  for (uint32 ti=0; ti<tigs.size(); ti++) {
    Unitig  *tig = tigs[ti];

    if ((tig == NULL) ||               //  Not a tig, ignore it.
        (tig->ufpath.size() == 1))     //  Singleton, already sorted.
      continue;

    tig->sort();
  }
}

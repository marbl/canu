
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
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

//
//
//  Candidate bubble tigs are found by annotating regions on each tig that
//  are covered by read-level overlaps to some other larger tig.  After
//  merging all overlapping regions, if a tig has ten or fewer regions, and
//  at least one end of the tig is covered, the tig will be considered for
//  bubble popping.  Additionally, a list of the tigs with overlapping reads
//  is kept for each candidate bubble tig.
//
//  Every read in a bubble tig is 'placed', using read-level overlaps, in all
//  other tigs.  A read can be 'placed' at a specific location in a tig if
//  the overlaps between it and the reads at that location are of similar
//  quality to the overlaps between just the reads at that location.
//  Additionally, a read can only be placed in a tig previously identified as
//  a potential location for the bubble.
//
//  Each placement of the first and last read in a tig is extended by the
//  length of the candidate bubble tig.  Overlapping placements are merged,
//  and the merged regions are exteded by 25% on each end.  Any region that
//  contains a placement for both the first and last read, correctly oriented
//  and sized, is retained.
//
//  Four outcomes are possible:
//
//  1) A single region is indentified, and every read in the bubble tig has
//     been "placed" in the region.  The candidate tig is merged into the
//     larger tig.
//
//  2) Multiple regions are identified, and every read in the candidate tig
//     are placed in every regoion.  The reads in the candidate tig are
//     individually placed at their best location.
//
//  3) Any number of regions are identified, and both the first and last read
//     in the candidate tig are placed.  The candidate tig is flagged as a
//     "bubble", and excluded from later repeat detection.
//
//  4) None of the above.  The candidate tig remains as is.
//
//


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
                     BubTargetList   &potentialOrphans,
                     bool             isBubble) {

  writeStatus("\n");
  writeStatus("findPotentialOrphans()-- working on " F_U32 " tigs.\n", tigs.size());

  writeLog("\n");
  writeLog("== Finding Potential %s ==\n", (isBubble ? "Bubbles" : "Orphans"));
  writeLog("\n");

  for (uint32 ti=0; ti<tigs.size(); ti++) {
    Unitig               *tig = tigs[ti];

    if ((tig == NULL) ||               //  Not a tig, ignore it.
        (tig->ufpath.size() == 1))     //  Singleton, handled elsewhere.
      continue;

    //  If the first or last read has no best edge, that's it, we're done.

    ufNode   *fRead = tig->firstRead();
    ufNode   *lRead = tig->lastRead();

    bool      fEdgeExists = (fRead->isForward() == true) ? (OG->bestEdgeExists(fRead->ident, false)) : (OG->bestEdgeExists(fRead->ident,  true));
    bool      lEdgeExists = (lRead->isForward() == true) ? (OG->bestEdgeExists(lRead->ident,  true)) : (OG->bestEdgeExists(lRead->ident, false));

#if 0
    if ((fEdgeExists == false) ||
        (lEdgeExists == false))
      continue;
#endif

    //  Count the number of reads that have an overlap to some other tig.  tigOlapsTo[otherTig] = count.

    intervalList<int32> tigCoverage;
    map<uint32,uint32>  tigOlapsTo;
    uint32              nonContainedReads = 0;

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode     *rdA   = &tig->ufpath[fi];
      uint32      rdAid =  tig->ufpath[fi].ident;

      if (OG->isContained(rdAid) == true)  //  Don't need to check contained reads.  If their container
        continue;                          //  passes the tests below, the contained read will too.

      nonContainedReads++;

      //  Find the list of tigs that this read has an overlap to.

      set<uint32>  readOlapsTo;

      uint32      ovlLen   = 0;
      BAToverlap *ovl      = OC->getOverlaps(rdAid, ovlLen);

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

        int32  mincoord = rdA->hangToMinCoord(ovl[oi].a_hang, ovl[oi].b_hang);
        int32  maxcoord = rdA->hangToMaxCoord(ovl[oi].a_hang, ovl[oi].b_hang);

        if (mincoord >= maxcoord)
          fprintf(stderr, "read %u at %u %u olap to read %u hangs %ld %ld -> coords %d %d\n",
                  ovl[oi].a_iid, rdA->position.bgn, rdA->position.end,
                  ovl[oi].b_iid, ovl[oi].a_hang, ovl[oi].b_hang,
                  mincoord, maxcoord);
        assert(mincoord < maxcoord);

        tigCoverage.add(mincoord, maxcoord - mincoord);
      }

      //  With the list of tigs that this read has an overlap to, add one to
      //  each tig in the list of tigs that this tig has an overlap to.

      for (set<uint32>::iterator it=readOlapsTo.begin(); it != readOlapsTo.end(); ++it)
        tigOlapsTo[*it]++;
    }

    //  Squash the tigcoverage down to intervals and decide if enough of this
    //  tig is overlapping anywhere else to consider this an orphan or
    //  bubble.

    tigCoverage.merge();

    uint32   spannedBases = 0;
    uint32   maxUncovered = 0;
    uint32   bgnUncovered = 0;
    uint32   endUncovered = 0;

    for (uint32 ii=0; ii<tigCoverage.numberOfIntervals(); ii++)
      spannedBases += tigCoverage.hi(ii) - tigCoverage.lo(ii);

    for (uint32 ii=1; ii<tigCoverage.numberOfIntervals(); ii++) {
      uint32   uncovered = tigCoverage.lo(ii) - tigCoverage.hi(ii-1);

      if (maxUncovered < uncovered)
        maxUncovered = uncovered;
    }

    if (tigCoverage.numberOfIntervals() > 0) {
      bgnUncovered =                    tigCoverage.lo(0);
      endUncovered = tig->getLength() - tigCoverage.hi( tigCoverage.numberOfIntervals()-1 );
    }

    writeLog("tig %8u length %9u nReads %7u/%7u - %3u regions covering %6.2f uncovered %5u/%6u/%5u edges %d %d -- ",
             tig->id(), tig->getLength(), nonContainedReads, tig->ufpath.size(),
             tigCoverage.numberOfIntervals(),
             100.0 * spannedBases / tig->getLength(),
             bgnUncovered, maxUncovered, endUncovered,
             fEdgeExists, lEdgeExists);

    //  Reject this tig as a potential bubble if
    //    there are more than 10 coverage intervals
    //    both bgn and end uncovered are non-zero
    //    the largest uncovered region ... ??
    //    bubbles don't care about intervals since only the ends must be covered

    if (!isBubble && tigCoverage.numberOfIntervals() > 10) {
      writeLog("too many intervals: is not orphan\n");
      continue;
    }

    if ((bgnUncovered > 0) ||
        (endUncovered > 0)) {
      writeLog("ends uncovered: is not orphan\n");
      continue;
    }

    //  Log the places where this orphan can go, and remember those places.

    writeLog("potential orphan\n");

    for (map<uint32,uint32>::iterator it=tigOlapsTo.begin(); it != tigOlapsTo.end(); ++it) {
      Unitig  *dest = tigs[it->first];

      writeLog("             tig %8u length %9u nReads %7u   %5u reads with overlaps\n",
               dest->id(), dest->getLength(), dest->ufpath.size(), it->second);

      potentialOrphans[ti].push_back(dest->id());
    }
  }  //  Over all tigs.

  writeStatus("findPotentialOrphans()-- found " F_SIZE_T " potential orphans.\n", potentialOrphans.size());

  writeLog("\n");
  writeLog("== Found %u Potential Orphans ==\n", potentialOrphans.size());
  writeLog("\n");

  //  Write non-orphan tigs.

#if 0
  for (uint32 ti=0; ti<tigs.size(); ti++) {
    Unitig               *tig = tigs[ti];

    if ((tig == NULL) ||               //  Not a tig, ignore it.
        (tig->ufpath.size() == 1))     //  Singleton, handled elsewhere.
      continue;

    if (potentialOrphans.count(ti) == 0)
      writeLog("tig %u of length %u with %u reads is NOT an orphan.\n", ti, tig->getLength(), tig->ufpath.size());
  }
#endif

  flushLog();
}



//  Find filtered placements for all the reads in the potential orphan tigs.

vector<overlapPlacement>  *
findOrphanReadPlacements(TigVector       &tigs,
                         BubTargetList   &potentialOrphans,
                         double           deviation,
                         double           similarity,
                         double           coverage) {
  uint32  fiLimit      = RI->numReads();
  uint32  fiNumThreads = omp_get_max_threads();
  uint32  fiBlockSize  = (fiLimit < 1000 * fiNumThreads) ? fiNumThreads : fiLimit / 999;

  vector<overlapPlacement>   *placed = new vector<overlapPlacement> [fiLimit + 1];

  writeLog("\n");
  writeLog("== Finding Read Placements for Potential Orphans ==\n");
  writeLog("\n");

#pragma omp parallel for schedule(dynamic, fiBlockSize)
  for (uint32 fi=1; fi<fiLimit; fi++) {
    uint32     rdAtigID = tigs.inUnitig(fi);

    if ((rdAtigID == 0) ||                           //  Read not placed in a tig, ignore it.
        (OG->isContained(fi)) ||                     //  Read is contained, ignore it.
        (potentialOrphans.count(rdAtigID) == 0))     //  Read isn't in a potential orphan, ignore it.
      continue;

    Unitig     *rdAtig   = tigs[rdAtigID];
    ufNode     *rdA      = &rdAtig->ufpath[ tigs.ufpathIdx(fi) ];
    int32       rdAlo    = rdA->position.min();
    int32       rdAhi    = rdA->position.max();

    bool        isEnd    = (rdAlo == 0) || (rdAhi == rdAtig->getLength());

    //  Compute all placements for this read.  It is critical to search for partial
    //  placements, otherwise we'll generally find no bubbles (only orphans).

    vector<overlapPlacement>   placements;

    placeReadUsingOverlaps(tigs, NULL, rdA->ident, placements, placeRead_noExtend);

    //  Weed out placements that aren't for orphans, or that are for orphans but are poor quality.  Or are to ourself!

    for (uint32 pi=0; pi<placements.size(); pi++) {
      uint32    rdBtigID = placements[pi].tigID;
      Unitig   *rdBtig   = tigs[rdBtigID];

      uint32    lo       = placements[pi].position.min();
      uint32    hi       = placements[pi].position.max();

      double    erate    = placements[pi].erate();

#if 0
      if (rdAtigID == rdBtigID)
        writeLog("tig %6u read %8u -> placed in source tig\n", rdAtigID, placements[pi].frgID);

      if (rdBtigID == 0)
        writeLog("tig %6u read %8u -> placed in singleton read (id == 0)\n", rdAtigID, placements[pi].frgID);

      if (rdBtig   == NULL)
        writeLog("tig %6u read %8u -> placed in singleton read (null ptr)\n", rdAtigID, placements[pi].frgID);

      if (rdBtig->ufpath.size() == 1)
        writeLog("tig %6u read %8u -> placed in singleton tig\n", rdAtigID, placements[pi].frgID);

      if ((potentialOrphans.count(rdBtigID) > 0) && (rdAtigID != rdBtigID))
        writeLog("tig %6u read %8u -> placed in orphan tig %u\n", rdAtigID, placements[pi].frgID, rdBtigID);
#endif

      if ((rdAtigID == rdBtigID) ||                     //  To ourself.
          (rdBtigID == 0) ||                            //  To a singleton read.
          (rdBtig   == NULL) ||                         //  To a singleton read.
          (rdBtig->ufpath.size() == 1) ||               //  To a singleton tig.
          (potentialOrphans.count(rdBtigID) > 0))       //  To a potential orphan tig
        continue;

      //  Ignore the placement if it isn't to one of our orphan-popping candidate tigs.

      assert(potentialOrphans.count(rdAtigID) == 1);

      bool             dontcare = true;
      vector<uint32>  &porphans = potentialOrphans[rdAtigID];

      for (uint32 pb=0; pb<porphans.size(); pb++)
        if (porphans[pb] == rdBtigID)
          dontcare = false;

      if (dontcare) {
        if (logFileFlagSet(LOG_ORPHAN_DETAIL))
          writeLog("tig %6u read %8u -> tig %6u (%6u reads) at %8u-%-8u (cov %7.5f erate %6.4f) - NOT CANDIDATE TIG\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Ignore the placement if it is to a potential orphan.

      if (potentialOrphans.count(rdBtigID) > 0) {
        if (logFileFlagSet(LOG_ORPHAN_DETAIL))
          writeLog("tig %6u read %8u -> tig %6u (%6u reads) at %8u-%-8u (cov %7.5f erate %6.4f) - INTO POTENTIAL ORPHAN\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      if (placements[pi].fCoverage < coverage) {    //  Ignore partially placed reads.
        if (logFileFlagSet(LOG_ORPHAN_DETAIL))
          writeLog("tig %6u read %8u -> tig %6u (%6u reads) at %8u-%-8u (cov %7.5f erate %6.4f) - LOW COVERAGE\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Ignore the placement if it is too diverged from the destination tig.

      double fGood = rdBtig->overlapConsistentWithTig(deviation, lo, hi, erate);

      if ((erate > similarity) && (fGood < 0.5)) {
        if (logFileFlagSet(LOG_ORPHAN_DETAIL))
          writeLog("tig %6u read %8u -> tig %6u (%6u reads) at %8u-%-8u (cov %7.5f erate %6.4f) - HIGH ERROR\n",
                   rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);
        continue;
      }

      //  Good placement!

      if (logFileFlagSet(LOG_ORPHAN_DETAIL))
        writeLog("tig %6u read %8u -> tig %6u (%6u reads) at %8u-%-8u (cov %7.5f erate %6.4f)\n",
                 rdAtigID, placements[pi].frgID, placements[pi].tigID, rdBtig->ufpath.size(), placements[pi].position.bgn, placements[pi].position.end, placements[pi].fCoverage, erate);

      placed[fi].push_back(placements[pi]);
    }
  }

  //  Done with the parallel.  Count things.

  uint32  nZeroTig   = 0;
  uint32  nContain   = 0;
  uint32  nNotOrphan = 0;

  uint32  nReads  = 0;
  uint32  nPlaces = 0;

  for (uint32 fi=1; fi<fiLimit; fi++) {
    uint32     rdAtigID = tigs.inUnitig(fi);

    if (rdAtigID == 0)
      nZeroTig++;

    if (OG->isContained(fi))
      nContain++;

    if (potentialOrphans.count(rdAtigID) == 0)
      nNotOrphan++;

    if ((rdAtigID == 0) ||                           //  Read not placed in a tig, ignore it.
        (OG->isContained(fi)) ||                     //  Read is contained, ignore it.
        (potentialOrphans.count(rdAtigID) == 0))     //  Read isn't in a potential orphan, ignore it.
      continue;

    nReads  += 1;
    nPlaces += placed[fi].size();
  }

  //  And report what we did.

  writeLog("\n");
  writeLog("== Found %u placements for %u reads. ==\n", nPlaces, nReads);
  writeLog("     %8u reads not placed: not in a tig\n",     nZeroTig);
  writeLog("     %8u reads not placed: contained\n",        nContain);
  writeLog("     %8u reads not placed: not in an orhpan\n", nNotOrphan);
  writeLog("\n");

  return(placed);
}




static
bool
placeAnchor(Unitig                     *orphan,
            vector<overlapPlacement>   *placed) {
  uint32   nReads = orphan->ufpath.size();
  ufNode  *fRead  = orphan->firstRead();
  ufNode  *lRead  = orphan->lastRead();

  assert(nReads > 0);

  //  Count the number of reads in the middle that are placed.

  uint32   placedD = 0, totalD = 0;
  uint32   placedC = 0, totalC = 0;

  for (uint32 fi=0; fi<nReads; fi++) {
    ufNode *rd   = &orphan->ufpath[fi];
    uint32  rdId = rd->ident;

    if ((rd == fRead) ||
        (rd == lRead))
      continue;

    if (OG->isContained(rdId)) {
      totalC++;

      if (placed[rdId].size() > 0)
        placedC++;
    }

    else {
      totalD++;

      if (placed[rdId].size() > 0)
        placedD++;
    }
  }

  //  The anchor isn't placed if either terminal read isn't placed.

  uint32   fPlaced = placed[fRead->ident].size();
  uint32   lPlaced = placed[lRead->ident].size();

  writeLog("Find anchors for orphan %u:\n", orphan->id());

  if (fPlaced == 0)   writeLog("  First read %6u - unplaced\n",         fRead->ident);
  else                writeLog("  First read %6u - placed %u time%s\n", fRead->ident, fPlaced, (fPlaced == 1) ? "" : "s");

  if (lPlaced == 0)   writeLog("   Last read %6u - unplaced\n",         lRead->ident);
  else                writeLog("   Last read %6u - placed %u time%s\n", lRead->ident, lPlaced, (lPlaced == 1) ? "" : "s");

  writeLog("     Internal reads - placed %6u/%-6u dovetail reads\n",  placedD, totalD);
  writeLog("                    - placed %6u/%-6u contained reads\n", placedC, totalC);

  return((fPlaced > 0) &&
         (lPlaced > 0));
}




static
void
addInitialIntervals(Unitig                               *orphan,
                    vector<overlapPlacement>             *placed,
                    ufNode                               *fRead,
                    ufNode                               *lRead,
                    map<uint32, intervalList<int32> *>   &targetIntervals) {
  uint32   orphanLen  = orphan->getLength();

  //  Add extended intervals for the first read.
  //
  //    target ---------------------------------------------
  //    read        -------
  //    orphan      -------------------------

  writeLog("\n");
  writeLog("  Intervals (first read):\n");

  for (uint32 pp=0; pp<placed[fRead->ident].size(); pp++) {
    uint32  tid = placed[fRead->ident][pp].tigID;
    uint32  bgn = placed[fRead->ident][pp].position.min();
    uint32  end = placed[fRead->ident][pp].position.max();

    if (targetIntervals[tid] == NULL)
      targetIntervals[tid] = new intervalList<int32>;

    //  If placed in the same orientation as in the tig, the orphan extends
    //  to the right of the min coordinate.
    //
    //  Otherwise, the the orphan extends to the left of the max coordinate.
    if (placed[fRead->ident][pp].position.isForward() == fRead->position.isForward()) {
      writeLog("    tig %8u %9u-%-9u ->\n", tid, bgn, bgn+orphanLen);
      targetIntervals[tid]->add(bgn, orphanLen);
    } else {
      writeLog("    tig %8u %9u-%-9u <-\n", tid, end-orphanLen, end);
      targetIntervals[tid]->add(end - orphanLen, orphanLen);
    }
  }

  //  Add extended intervals for the last read.
  //
  //    target ---------------------------------------------
  //    read                          -------
  //    orphan      -------------------------

  writeLog("\n");
  writeLog("  Intervals (last read):\n");

  for (uint32 pp=0; pp<placed[lRead->ident].size(); pp++) {
    uint32  tid = placed[lRead->ident][pp].tigID;
    uint32  bgn = placed[lRead->ident][pp].position.min();
    uint32  end = placed[lRead->ident][pp].position.max();

    if (targetIntervals[tid] == NULL)
      targetIntervals[tid] = new intervalList<int32>;

    //  Same as above, just backwards.

    if (placed[lRead->ident][pp].position.isForward() == lRead->position.isForward()) {
      writeLog("    tig %8u %9u-%-9u ->\n", tid, end-orphanLen, end);
      targetIntervals[tid]->add(end - orphanLen, orphanLen);
    } else {
      writeLog("    tig %8u %9u-%-9u <-\n", tid, bgn, bgn+orphanLen);
      targetIntervals[tid]->add(bgn, orphanLen);
    }
  }
}



static
void
saveCorrectlySizedInitialIntervals(Unitig                    *orphan,
                                   Unitig                    *target,
                                   intervalList<int32>       *IL,
                                   ufNode                    *fRead,
                                   ufNode                    *lRead,
                                   vector<overlapPlacement>  *placed,
                                   vector<candidatePop *>    &targets) {
  uint32   orphanLen  = orphan->getLength();

  IL->merge();  // Merge overlapping initial intervals.

  writeLog("\n");
  writeLog("Finding intervals for orphan %u placed in tig %u.\n", orphan->id(), target->id());

  //  Search all the intervals we think an orphan/bubble can go, and decide
  //  if both the first and last read in the orphan are:
  //    placed in the interval
  //    oriented
  //    the correct size
  //
  //    ----------------[----------------]--------------   //  initial interval size
  //    -----------[--------------------------]---------   //  interval extended by 50% of the orphan size
  //                   -->                <--              //  reads placed
  //                    --------------------               //  compared to orpan itself
  //
  for (uint32 ii=0; ii<IL->numberOfIntervals(); ii++) {
    int32  intBgn   = IL->lo(ii) - 0.25 * orphanLen;   //  Extend the region by 50% of the
    int32  intEnd   = IL->hi(ii) + 0.25 * orphanLen;   //  orphan length.

    intBgn = max(intBgn, 0);
    intEnd = min(intEnd, target->getLength());

    SeqInterval    fPos;
    SeqInterval    lPos;

    //  Search placements for a valid placement pair.

    vector<overlapPlacement>  &fPlaces = placed[fRead->ident];
    vector<overlapPlacement>  &lPlaces = placed[lRead->ident];

    //  Over all the first read placements...
    for (uint32 fp=0; fp<fPlaces.size(); fp++) {
      if ((target->id() != fPlaces[fp].tigID) ||      //  Placed in wrong tig
          (fPlaces[fp].position.min() < intBgn) ||    //  Read not placed fully
          (intEnd    < fPlaces[fp].position.max()))   //  in the region.
        continue;

      //  First read is in this region.  Decide if the tig should be aligned
      //  forward or reverse based on the alignment of this read.

      bool  fPlaceForward = (fPlaces[fp].position.isForward() == fRead->position.isForward()) ? true : false;

      //  Over all the last read placements...
      for (uint32 lp=0; lp<lPlaces.size(); lp++) {
        if ((target->id() != lPlaces[lp].tigID) ||      //  Placed in wrong tig
            (lPlaces[lp].position.min() < intBgn) ||    //  Read not placed fully
            (intEnd    < lPlaces[lp].position.max()))   //  in the region.
          continue;

        //  Second read is in this region.  Decide if the tig should be
        //  aligned forward or reverse, again based on only this read.

        bool  lPlaceForward = (lPlaces[lp].position.isForward() == lRead->position.isForward()) ? true : false;

        //  If they disagree, this isn't a valid placement.

        bool  misOrient = (fPlaceForward != lPlaceForward) ? true : false;

        //  Decide if their order is correct, and if the length is
        //  appropriate.

        int32  pBgn     = (fPlaceForward) ? (fPlaces[fp].position.min()) : (lPlaces[lp].position.min());
        int32  pEnd     = (fPlaceForward) ? (lPlaces[lp].position.max()) : (fPlaces[fp].position.max());
        int32  length   = pEnd - pBgn;

        bool   misOrder = (length < 0) ? true : false;
        bool   tooSmall = (length < 0.75 * orphan->getLength()) ? true : false;
        bool   tooLarge = (length > 1.25 * orphan->getLength()) ? true : false;

        if (misOrient) {
          writeLog("  %9d-%-9d %7.1f%% of orphan length - first read at %9d-%-9d last read at %9d-%-9d  MIS-ORIENT\n",
                   pBgn, pEnd, 100.0 * length / orphan->getLength(),
                   fPlaces[fp].position.min(), fPlaces[fp].position.max(),
                   lPlaces[lp].position.min(), lPlaces[lp].position.max());
          continue;
        }

        if (misOrder) {
          writeLog("  %9d-%-9d %7.1f%% of orphan length - first read at %9d-%-9d last read at %9d-%-9d  MIS-ORDER\n",
                   pBgn, pEnd, 100.0 * length / orphan->getLength(),
                   fPlaces[fp].position.min(), fPlaces[fp].position.max(),
                   lPlaces[lp].position.min(), lPlaces[lp].position.max());
          continue;
        }

        if (tooSmall) {
          writeLog("  %9d-%-9d %7.1f%% of orphan length - first read at %9d-%-9d last read at %9d-%-9d  TOO SMALL\n",
                   pBgn, pEnd, 100.0 * length / orphan->getLength(),
                   fPlaces[fp].position.min(), fPlaces[fp].position.max(),
                   lPlaces[lp].position.min(), lPlaces[lp].position.max());
          continue;
        }

        if (tooLarge) {
          writeLog("  %9d-%-9d %7.1f%% of orphan length - first read at %9d-%-9d last read at %9d-%-9d  TOO LARGE\n",
                   pBgn, pEnd, 100.0 * length / orphan->getLength(),
                   fPlaces[fp].position.min(), fPlaces[fp].position.max(),
                   lPlaces[lp].position.min(), lPlaces[lp].position.max());
          continue;
        }

        //  A valid placement.

        writeLog("  %9d-%-9d %7.1f%% of orphan length - first read at %9d-%-9d last read at %9d-%-9d  SUCCESS!\n",
                 pBgn, pEnd, 100.0 * length / orphan->getLength(),
                 fPlaces[fp].position.min(), fPlaces[fp].position.max(),
                 lPlaces[lp].position.min(), lPlaces[lp].position.max());

        targets.push_back(new candidatePop(orphan, target, pBgn, pEnd));
      }
    }
  }

  delete IL;
}



void
assignReadsToTargets(Unitig                     *orphan,
                     vector<overlapPlacement>   *placed,
                     vector<candidatePop *>     targets) {

  //  For each read in the orphan,
  //  For each placement of the read,
  //  For each target location
  //    If the target tig is the same as the placement tig
  //    and the placement of the read is contained in the target region
  //      save the placement to a list of placements for this target


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

  writeLog("\n");
  writeLog("Removing duplicate placements.\n");

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

        if (t->placed[aa].erate() < t->placed[bb].erate()) {
          save = aa;
          remo = bb;
        } else {
          save = bb;
          remo = aa;
        }

        writeLog("  duplicate read alignment for tig %u read %u - better %u-%-u %.4f - worse %u-%-u %.4f\n",
                 t->placed[save].tigID, t->placed[save].frgID,
                 t->placed[save].position.bgn, t->placed[save].position.end, t->placed[save].erate(),
                 t->placed[remo].position.bgn, t->placed[remo].position.end, t->placed[remo].erate());

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

  if (nDup > 0) {
    writeLog("\n");
    writeLog("Removed %u duplicate placement%s.\n", nDup, (nDup == 1) ? "" : "s");
  }
}



void
mergeOrphans(TigVector    &tigs,
             double        deviation,
             double        similarity,
             bool          isBubble) {


  //  Find, for each tig, the list of other tigs that it could potentially be placed into.

  BubTargetList   potentialOrphans;

  findPotentialOrphans(tigs, potentialOrphans, isBubble);

  //  If you enable this, all reads with any overlap will get removed from the 'reduced' graph.
#if 0
  for (auto it = potentialOrphans.begin(); it != potentialOrphans.end(); it++) {
    uint32   tid = it->first;
    Unitig  *tig = tigs[tid];

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++)
      OG->setBubble(tig->ufpath[fi].ident);
  }
#endif

  //  For any tig that is a potential orphan, find all read placements.
  //  We don't try to insert orphans if the reads aren't fully contained, for bubbles no minimum threshold as the minimum is the shortest overlap we are willing to consider
  vector<overlapPlacement>   *placed = findOrphanReadPlacements(tigs, potentialOrphans, deviation, similarity, (isBubble ? 0.01 : 0.99));

  //  We now have, in 'placed', a list of all the places that each read could be placed.  Decide if there is a _single_
  //  place for each orphan to be popped.

  uint32        nNeither    = 0, nNeitherReads    = 0;
  uint32        nUniqBubble = 0, nUniqBubbleReads = 0;
  uint32        nUniqOrphan = 0, nUniqOrphanReads = 0;
  uint32        nReptOrphan = 0, nReptOrphanReads = 0;

  for (uint32 ti=0; ti<tigs.size(); ti++) {
    Unitig  *orphan        = tigs[ti];

    if (potentialOrphans.count(ti) == 0)
      continue;

    writeLog("\n");
    writeLog("========================================\n");
    writeLog("Processing potential orphan %u of length %u bp with %u reads\n", ti, orphan->getLength(), orphan->ufpath.size());
    writeLog("\n");

    //  Scan the orphan, decide if there are _ANY_ read placements.  Log appropriately.

    if (placeAnchor(orphan, placed) == false) {
      writeLog("\n");
      writeLog("ANCHOR READS FAILED TO PLACE.\n");
      continue;
    }

    //  Create intervals for each placed read.
    //
    //    target ---------------------------------------------
    //    read        -------
    //    orphan      -------------------------

    ufNode  *fRead = orphan->firstRead();
    ufNode  *lRead = orphan->lastRead();

    map<uint32, intervalList<int32> *>   targetIntervals;

    addInitialIntervals(orphan, placed, fRead, lRead, targetIntervals);

    //  Figure out if each interval has both the first and last read of some orphan, and if those
    //  are properly sized.  If so, save a candidatePop.

    vector<candidatePop *>    targets;

    for (auto it=targetIntervals.begin(); it != targetIntervals.end(); ++it)
      if (tigs[it->first] == NULL)
        writeLog("WARNING: Orphan %u wants to go into nonexistent tig %u!\n", ti, it->first);
      else
        saveCorrectlySizedInitialIntervals(orphan,
                                           tigs[it->first],     //  The targetID      in targetIntervals
                                           it->second,          //  The interval list in targetIntervals
                                           fRead,
                                           lRead,
                                           placed,
                                           targets);

    targetIntervals.clear();   //  intervalList already freed.

    writeLog("\n");
    writeLog("Found %u target location%s\n", targets.size(), (targets.size() == 1) ? "" : "s");

    //  If no targets, nothing to do.

    if (targets.size() == 0)
      continue;

    //  Assign read placements to targets.

    assignReadsToTargets(orphan, placed, targets);

    //  Compare the orphan against each target.

    uint32   nOrphan      = 0;   //  Number of targets that have all the reads.
    uint32   nBubble      = 0;   //  Number of targets that have some reads placed.
    bool     repeatBubble = false;
    uint32   orphanTarget = 0;   //  If nOrphan == 1, the target we're popping into.

    for (uint32 tt=0; tt<targets.size(); tt++) {
      uint32  orphanSize   = orphan->ufpath.size();        //  Size of the orphan tig
      uint32  targetSize   = targets[tt]->placed.size();   //  Number of those reads placed at this target
      uint32  terminalSize = 0;                            //  Number of terminal reads in the orphan placed

      //  Count how many terminal reads are placed.  We don't know where they
      //  are in the list, so need to check every read.

      for (uint32 op=0; op<targets[tt]->placed.size(); op++) {
        if (targets[tt]->placed[op].frgID == fRead->ident)    terminalSize++;
        if (targets[tt]->placed[op].frgID == lRead->ident)    terminalSize++;
      }

      // add interval
      if (targetIntervals[targets[tt]->target->id()] == NULL)
         targetIntervals[targets[tt]->target->id()] = new intervalList<int32>;
      targetIntervals[targets[tt]->target->id()]->add(targets[tt]->bgn, targets[tt]->end-targets[tt]->bgn);

      //  Report now, before we nuke targets[tt] for being not a orphan!

      writeLog("\n");
      writeLog("Placing orphan %u (length %u) into tig %u at position %u-%u (length %u):\n",
               orphan->id(), orphan->getLength(),
               targets[tt]->target->id(), targets[tt]->bgn, targets[tt]->end, targets[tt]->end - targets[tt]->bgn);

      for (uint32 op=0; op<targets[tt]->placed.size(); op++)
        writeLog("    read %7u at %9u-%-9u\n",
                 targets[tt]->placed[op].frgID,
                 targets[tt]->placed[op].position.bgn, targets[tt]->placed[op].position.end);

      writeLog("  out of %u reads, found %u reads placed, %u terminal reads\n", orphanSize, targetSize, terminalSize);

      //  If all reads placed, we can merge this orphan into the target.  But
      //  if this happens more than once, we just split the orphan and merge
      //  reads at their best location.

      if (orphanSize == targetSize && !isBubble) {
        nOrphan++;
        orphanTarget = tt;
      }

      //  If only some of the reads are placed or we're in the bubble rounds, declare this a bubble so we
      //  can ignore reads when finding repeats.

      else if (orphanSize == targetSize || terminalSize == 2) {
        nBubble++;
      }
    }

    // use the merged intervals to check for if this is a repeat
    for (auto it=targetIntervals.begin(); it != targetIntervals.end(); ++it) {
       it->second->merge();
       writeLog("Merged intervals of orphan %u with %d equivalents into tig %u with %u intervals, first is %u-%u (length %u), ratio %f\n", orphan->id(), nBubble, it->first, it->second->numberOfIntervals(), it->second->lo(0), it->second->hi(0), it->second->hi(0)-it->second->lo(0), double(it->second->hi(0) - it->second->lo(0)) / orphan->getLength());
       // if this is being placed in multiple disjoint locations in a tig or in a single but extremely large location, it is a repeat
       if (it->second->numberOfIntervals() > 1 || double((it->second->hi(0) - it->second->lo(0)) / orphan->getLength()) >= 5) {
          repeatBubble = true;
          break;
       }
       delete it->second;
    }
    // lastly even if it met the above criteria, if there weren't too many placements it's not a repat
    // essentially our critera for repeat is lots of possibilities and they're not all to one/a few tig regions
    repeatBubble = repeatBubble && (nBubble >= 50);
    targetIntervals.clear();

    //
    //  If neither, be obnoxious.
    //

    if ((nOrphan == 0) && (nBubble == 0)) {
      writeLog("\n");
      writeLog("Result:\n");
      writeLog("  tig %8u of length %8u with %6u reads - NO GOOD PLACEMENTS\n", orphan->id(), orphan->getLength(), orphan->ufpath.size());

      nNeither      += 1;
      nNeitherReads += orphan->ufpath.size();
    }

    //
    //  If not an orphan, but a bubble, flag the tig and all reads as a bubble.
    //

    if ((nOrphan == 0) && (nBubble > 0)) {
      writeLog("\n");
      writeLog("Result:\n");
      writeLog("  tig %8u of length %8u with %6u reads - BUBBLE%s\n", orphan->id(), orphan->getLength(), orphan->ufpath.size(), repeatBubble ? " TOO MANY PLACEMENTS" : "");

      if (!repeatBubble) {
        nUniqBubble      += 1;
        nUniqBubbleReads += orphan->ufpath.size();

        orphan->_isBubble = true;

        for (uint32 fi=0; fi<orphan->ufpath.size(); fi++)
          OG->setBubble(orphan->ufpath[fi].ident);
      }
    }
    //
    //  If a unique orphan placement, place it there.
    //

    if (nOrphan == 1) {
      writeLog("\n");
      writeLog("Result:\n");
      writeLog("  tig %8u of length %8u with %6u reads - UNIQUELY PLACED ORPHAN\n", orphan->id(), orphan->getLength(), orphan->ufpath.size());

      nUniqOrphan      += 1;
      nUniqOrphanReads += orphan->ufpath.size();

      for (uint32 op=0, tt=orphanTarget; op<targets[tt]->placed.size(); op++)     //  Move all the reads to
        targets[tt]->target->addRead(ufNode(targets[tt]->placed[op].frgID,        //  their placed position.
                                            targets[tt]->placed[op].position));

      for (uint32 fi=0; fi<orphan->ufpath.size(); fi++) {                           //  Flag them as being an orphan.
        OG->setOrphan(orphan->ufpath[fi].ident);
        OG->setBackbone(orphan->ufpath[fi].ident, false);                           // and reset their backbone status, they're not the backbone of the new tig they've been placed into
      }

      tigs[orphan->id()] = NULL;                                                  //  Delete the original tig.
      delete orphan;
    }

    //
    //  If multiply placed, we can't distinguish between them, and
    //  instead just place reads where they individually decide to go.
    //

    if (nOrphan > 1) {
      writeLog("\n");
      writeLog("Result:\n");
      writeLog("  tig %8u of length %8u with %6u reads %6u - MULTIPLY PLACED ORPHAN\n", orphan->id(), orphan->getLength(), orphan->ufpath.size());

      nReptOrphan      += 1;
      nReptOrphanReads += orphan->ufpath.size();

      for (uint32 fi=0; fi<orphan->ufpath.size(); fi++) {
        uint32  rr = orphan->ufpath[fi].ident;
        uint32  bb = 0;

        for (uint32 pp=0; pp<placed[rr].size(); pp++)                //  Over all placements for this read, pick
          if ((placed[rr][pp].erate() < placed[rr][bb].erate()) &&   //  the lowest error one, that isn't in
              (placed[rr][pp].tigID  != orphan->id()))               //  the orphan tig.
            bb = pp;

        assert(rr == placed[rr][bb].frgID);

        Unitig  *target = tigs[placed[rr][bb].tigID];

        assert(target       != NULL);
        assert(target->id() != orphan->id());

        target->addRead(ufNode(placed[rr][bb].frgID,
                               placed[rr][bb].position));

        OG->setOrphan(placed[rr][bb].frgID);
      }

      writeLog("\n");

      tigs[orphan->id()] = NULL;
      delete orphan;
    }

    //
    //  All done with this orphan.  Clean up the targets list.
    //

    for (uint32 tt=0; tt<targets.size(); tt++) {
      delete targets[tt];
      targets[tt] = NULL;
    }

    targets.clear();

  }  //  Over all orphans

  writeLog("\n");   //  Needed if no orphans are popped.

  writeStatus("mergeOrphans()-- flagged   %5u        bubble tigs with %u reads\n", nUniqBubble, nUniqBubbleReads);
  writeStatus("mergeOrphans()-- placed    %5u unique orphan tigs with %u reads\n", nUniqOrphan, nUniqOrphanReads);
  writeStatus("mergeOrphans()-- shattered %5u repeat orphan tigs with %u reads\n", nReptOrphan, nReptOrphanReads);
  writeStatus("mergeOrphans()-- ignored   %5u               tigs with %u reads; failed to place\n", nNeither, nNeitherReads);
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

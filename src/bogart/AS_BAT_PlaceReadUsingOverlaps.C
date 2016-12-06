
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
 *    src/bogart/AS_BAT_PlaceFragUsingOverlaps.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2016-AUG-12
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
#include "AS_BAT_PlaceReadUsingOverlaps.H"

#include "intervalList.H"


#undef TEST_ALT



overlapPlacement *
placeRead_fromOverlaps(TigVector   &tigs,
                       Unitig      *target,
                       uint32       fid,
                       uint32       flags,
                       uint32       ovlLen,
                       BAToverlap   *ovl) {
  overlapPlacement *ovlPlace = new overlapPlacement[ovlLen];

  for (uint32 i=0; i<ovlLen; i++) {
    int32             tigID = tigs.inUnitig(ovl[i].b_iid);
    Unitig           *tig   = tigs[tigID];

    assert(ovl[i].a_iid == fid);

    if (tigID == 0)                           //  Skip if overlapping read isn't in a tig yet - unplaced contained, or garbage read.
      continue;

    if ((target != NULL) && (target != tig))  //  Skip if we requested a specific tig and if this isn't it.
      continue;

    //  Place the read relative to the other read.

    BestEdgeOverlap   edge(ovl[i]);
    ufNode            read;

    if (tig->placeRead(read, fid, ovl[i].AEndIs3prime(), &edge) == false) {
      if (logFileFlagSet(LOG_PLACE_READ))
        writeLog("pRUO()-- WARNING: Failed to place with overlap %u %u hangs %u %u flipped %u\n",
                 ovl[i].a_iid, ovl[i].b_iid, ovl[i].a_hang, ovl[i].b_hang, ovl[i].flipped);
      continue;
    }

    //  Save the placement in our work space.

    uint32  olen = RI->overlapLength(ovl[i].a_iid, ovl[i].b_iid, ovl[i].a_hang, ovl[i].b_hang);
    uint32  flen = RI->readLength(ovl[i].a_iid);

    ovlPlace[i].frgID        = fid;
    ovlPlace[i].refID        = ovl[i].b_iid;
    ovlPlace[i].tigID        = tig->id();
    ovlPlace[i].position     = read.position;
    ovlPlace[i].verified.bgn = INT32_MAX;
    ovlPlace[i].verified.end = INT32_MIN;
    ovlPlace[i].covered.bgn  = (ovl[i].a_hang < 0) ?    0 : ovl[i].a_hang;          //  The portion of the read
    ovlPlace[i].covered.end  = (ovl[i].b_hang > 0) ? flen : ovl[i].b_hang + flen;   //  covered by the overlap.
    ovlPlace[i].clusterID    = 0;
    ovlPlace[i].fCoverage    = 0.0;
    ovlPlace[i].errors       = olen * ovl[i].erate();
    ovlPlace[i].aligned      = ovlPlace[i].covered.end - ovlPlace[i].covered.bgn;
    ovlPlace[i].tigFidx      = UINT32_MAX;
    ovlPlace[i].tigLidx      = 0;

    assert(ovlPlace[i].covered.bgn >= 0);
    assert(ovlPlace[i].covered.end >= 0);
    assert(ovlPlace[i].covered.bgn <= flen);
    assert(ovlPlace[i].covered.end <= flen);
    assert(ovlPlace[i].covered.bgn < ovlPlace[i].covered.end);

    //  Disallow any placements that exceed the boundary of the unitig.  These cannot be confirmed
    //  by overlaps and might be wrong.  Sample cases:
    //    o  sticking a unique/repeat read onto a repeat (leaving the unique uncovered)
    //    o  sticking a chimeric read onto the end of a unitig (leaving the chimeric join uncovered)

    if (((flags & placeRead_fullMatch) ||
         (flags & placeRead_noExtend)) &&
        ((ovlPlace[i].position.min() < 0) ||
         (ovlPlace[i].position.max() > tig->getLength()))) {
      ovlPlace[i] = overlapPlacement();
    }

    //  Report the placement.

    if (logFileFlagSet(LOG_PLACE_READ))
      writeLog("pRUO()-- read %7d (%5d,%5d) in unitig %5d at %8d,%-8d via read %7d at %8d:%-8d hang %6d %6d %s%s\n",
               ovlPlace[i].frgID,
               ovlPlace[i].covered.bgn, ovlPlace[i].covered.end,
               ovlPlace[i].tigID,
               ovlPlace[i].position.bgn, ovlPlace[i].position.end,
               ovl[i].b_iid,
               tig->readFromId(ovl[i].b_iid)->position.bgn,
               tig->readFromId(ovl[i].b_iid)->position.end,
               ovl[i].a_hang, ovl[i].b_hang,
               (ovl[i].flipped == true) ? "<--" : "-->",
               (ovlPlace[i].frgID == 0) ? " DISALLOWED" : "");
  }  //  Over all overlaps.

  return(ovlPlace);
}



void
placeRead_assignEndPointsToCluster(uint32  bgn, uint32  end,
                                   uint32  fid,
                                   overlapPlacement     *ovlPlace,
                                   intervalList<int32>  &bgnPoints,
                                   intervalList<int32>  &endPoints) {
  int32  windowSlop = 0.075 * RI->readLength(fid);

  if (windowSlop < 5)
    windowSlop = 5;

  for (uint32 oo=bgn; oo<end; oo++) {
    bgnPoints.add(ovlPlace[oo].position.bgn - windowSlop, 2 * windowSlop);
    endPoints.add(ovlPlace[oo].position.end - windowSlop, 2 * windowSlop);
  }

  bgnPoints.merge();
  endPoints.merge();

  if (logFileFlagSet(LOG_PLACE_READ)) {
    writeLog("pRUO()-- Using windowSlop %d\n", windowSlop);
    writeLog("pRUO()-- Found %3u bgn interval%s", bgnPoints.numberOfIntervals(), bgnPoints.numberOfIntervals() == 1 ? ":  " : "s: ");
    for (uint32 r=0; r<bgnPoints.numberOfIntervals(); r++)
      writeLog(" %6d:%-6d", bgnPoints.lo(r), bgnPoints.hi(r));
    writeLog("\n");

    writeLog("pRUO()-- Found %3u end interval%s", endPoints.numberOfIntervals(), endPoints.numberOfIntervals() == 1 ? ":  " : "s: ");
    for (uint32 r=0; r<endPoints.numberOfIntervals(); r++)
      writeLog(" %6d:%-6d", endPoints.lo(r), endPoints.hi(r));
    writeLog("\n");
  }
}



void
placeRead_assignPlacementsToCluster(uint32  bgn, uint32  end,
                                    uint32  fid,
                                    overlapPlacement     *ovlPlace,
                                    intervalList<int32>  &bgnPoints,
                                    intervalList<int32>  &endPoints) {
  int32   numBgnPoints = bgnPoints.numberOfIntervals();
  int32   numEndPoints = endPoints.numberOfIntervals();

  for (uint32 oo=bgn; oo<end; oo++) {
    int32   b = ovlPlace[oo].position.bgn;
    int32   e = ovlPlace[oo].position.end;
    int32   c = 0;

    ovlPlace[oo].clusterID = 0;

    for (int32 r=0; r<numBgnPoints; r++)
      if ((bgnPoints.lo(r) <= b) && (b <= bgnPoints.hi(r))) {
        assert(ovlPlace[oo].clusterID == 0);  //  Obvious; we just set it to zero above.
        ovlPlace[oo].clusterID = c = r * numEndPoints + 1;
      }

    for (int32 r=0; r<numEndPoints; r++)
      if ((endPoints.lo(r) <= e) && (e <= endPoints.hi(r))) {
        assert(ovlPlace[oo].clusterID == c);  //  Otherwise, bgn point wasn't placed in a cluster!
        ovlPlace[oo].clusterID += r;
      }
  }
}



void
placeRead_findFirstLastOverlapping(overlapPlacement &op,
                                   Unitig *tig,
                                   uint32 os, uint32 oe,
                                   overlapPlacement *ovlPlace) {
  op.tigFidx = UINT32_MAX;
  op.tigLidx = 0;

  for (uint32 oo=os; oo<oe; oo++) {
    uint32   ord = tig->ufpathIdx(ovlPlace[oo].refID);

    op.tigFidx = min(ord, op.tigFidx);
    op.tigLidx = max(ord, op.tigLidx);

    //if (logFileFlagSet(LOG_PLACE_READ))
    //  writeLog("pRUO()--     find range from os=%u to oe=%u  tig=%u  ord=%u  f=%u l=%u\n",
    //            os, oe, op.tigID, ord, op.tigFidx, op.tigLidx);
  }

  if (op.tigFidx > op.tigLidx)
    writeStatus("pRUO()-- Invalid placement indices: tigFidx %u tigLidx %u\n", op.tigFidx, op.tigLidx);
  assert(op.tigFidx <= op.tigLidx);

  if (logFileFlagSet(LOG_PLACE_READ))
    writeLog("pRUO()--   spans reads #%u (%u) to #%u (%u) in tig %u\n",
             op.tigFidx, tig->ufpath[op.tigFidx].ident,
             op.tigLidx, tig->ufpath[op.tigLidx].ident,
             op.tigID);
}



void
placeRead_computeQualityAndCoverage(overlapPlacement &op,
                                    uint32 os, uint32 oe,
                                    overlapPlacement *ovlPlace) {
  op.errors  = 0;
  op.aligned = 0;

  op.covered.bgn = INT32_MAX;  //  Covered interval is always in
  op.covered.end = INT32_MIN;  //  forward read coordinates

  for (uint32 oo=os; oo<oe; oo++) {
    if ((ovlPlace[oo].position.bgn == 0) &&
        (ovlPlace[oo].position.end == 0)) {
      if (logFileFlagSet(LOG_PLACE_READ))
        writeLog("OLD place=%3d  read %8d ref read %8d - covered %5d:%-5d with %6.1f errors - DELETED\n",
                 op.frgID, ovlPlace[oo].refID, ovlPlace[oo].covered.bgn, ovlPlace[oo].covered.end, ovlPlace[oo].errors);
      continue;
    }

    op.errors      += ovlPlace[oo].errors;
    op.aligned     += ovlPlace[oo].aligned;

    op.covered.bgn  = min(op.covered.bgn, ovlPlace[oo].covered.bgn);
    op.covered.end  = max(op.covered.end, ovlPlace[oo].covered.end);

    //if (logFileFlagSet(LOG_PLACE_READ))
    //  writeLog("OLD place=%3d  read %8d ref read %8d - covered %5d:%-5d with %6.1f errors\n",
    //           oo, op.frgID, ovlPlace[oo].refID, ovlPlace[oo].covered.bgn, ovlPlace[oo].covered.end, ovlPlace[oo].errors);
  }

  op.fCoverage = (op.covered.end - op.covered.bgn) / (double)RI->readLength(op.frgID);
}



void
placeRead_computeQualityAndCoverage(overlapPlacement &op,
                                    BAToverlap       *ovl,
                                    uint32            ovlLen,
                                    set<uint32>      &reads) {
  op.errors  = 0;
  op.aligned = 0;

  op.covered.bgn = INT32_MAX;  //  Covered interval is always in
  op.covered.end = INT32_MIN;  //  forward read coordinates

  //  For reads that have two overlaps to the same other read, we have no way of knowing
  //  which is the correct overlap, just that we have an overlap.
  //
  //  This happens in dros a whole bunch of times, and does change the fCoverave value.

  for (uint32 oo=0; oo<ovlLen; oo++) {
    if (reads.count(ovl[oo].b_iid) == 0)
      continue;

    int32   olen = RI->overlapLength(ovl[oo].a_iid, ovl[oo].b_iid, ovl[oo].a_hang, ovl[oo].b_hang);
    int32   flen = RI->readLength(ovl[oo].a_iid);

    int32   cbgn = (ovl[oo].a_hang < 0) ?    0 : ovl[oo].a_hang;          //  The portion of the read
    int32   cend = (ovl[oo].b_hang > 0) ? flen : ovl[oo].b_hang + flen;   //  covered by the overlap.

    //if (logFileFlagSet(LOG_PLACE_READ))
    //  writeLog("NEW place=%3d  read %8d ref read %8d - covered %5d:%-d with %f errors\n",
    //           op.frgID, ovlPlace[oo].refID, cbgn, cend, olen * ovl[oo].erate());

    op.errors   += olen * ovl[oo].erate();
    op.aligned  += cend - cbgn;

    op.covered.bgn = min(op.covered.bgn, cbgn);
    op.covered.end = max(op.covered.end, cend);
  }

  op.fCoverage = (op.covered.end - op.covered.bgn) / (double)RI->readLength(op.frgID);
}




//  Now that it is placed, estimate the span that is verified by overlaps.
//  Threshold the floating end so it doesn't exceed the placement.
//
//  Annoyingly, the verified placement can, and does, exceed the bounds of the
//  unitig, and we  need to check that threshold too.  Indel in the read and all that.
//
void
placeRead_computeVerified(overlapPlacement &op, uint32 tigLen) {

  //writeLog("computeVer pos %d-%d cov %d-%d\n", op.position.bgn, op.position.end, op.covered.bgn, op.covered.end);

  if (op.position.isForward()) {
    op.verified.bgn = op.position.bgn + op.covered.bgn;
    op.verified.end = op.position.bgn + op.covered.end;

    if (op.verified.end > op.position.end)   //  verified.bgn is always valid if covered.bgn > 0
      op.verified.end = op.position.end;     //

    if (op.verified.bgn < 0)
      op.verified.bgn = 0;
    if (op.verified.end > tigLen)
      op.verified.end = tigLen;

    assert(op.verified.bgn >= op.position.bgn);
    assert(op.verified.end <= op.position.end);
    assert(op.verified.bgn <  op.verified.end);
  }

  else {
    op.verified.bgn = op.position.bgn - op.covered.bgn;  //  High coord
    op.verified.end = op.position.bgn - op.covered.end;  //  Low coord

    if (op.verified.end < op.position.end)  //  verified.bgn is always valid if covered.bgn > 0
      op.verified.end = op.position.end;

    if (op.verified.end < 0)
      op.verified.end = 0;
    if (op.verified.bgn > tigLen)
      op.verified.bgn = tigLen;

    assert(op.verified.end >= op.position.end);
    assert(op.verified.bgn <= op.position.bgn);
    assert(op.verified.end <  op.verified.bgn);
  }

  assert(op.position.isForward() == op.verified.isForward());
}



void
placeRead_computePlacement(overlapPlacement &op,
                           uint32            os,
                           uint32            oe,
                           overlapPlacement *ovlPlace,
                           Unitig           *tig) {
  stdDev<double>   bgnPos;
  stdDev<double>   endPos;

  for (uint32 oo=os; oo<oe; oo++) {
    if ((ovlPlace[oo].position.bgn == 0) &&
        (ovlPlace[oo].position.end == 0))
      continue;

    //writeLog("OLD place %d-%d\n", ovlPlace[oo].position.bgn, ovlPlace[oo].position.end);

    bgnPos.insert(ovlPlace[oo].position.bgn);
    endPos.insert(ovlPlace[oo].position.end);
  }

  bgnPos.finalize();
  endPos.finalize();

  op.position.bgn = bgnPos.mean();
  op.position.end = endPos.mean();

  placeRead_computeVerified(op, tig->getLength());
}




void
placeRead_computePlacement(overlapPlacement &op,
                           BAToverlap       *ovl,
                           uint32            ovlLen,
                           set<uint32>      &reads,
                           uint32            flags,
                           Unitig           *tig) {
  stdDev<double>   bgnPos;
  stdDev<double>   endPos;

  //  For reads that have two overlaps to the same other read, we have no way of knowing
  //  which is the correct overlap, just that we have an overlap.

  for (uint32 oo=0; oo<ovlLen; oo++) {
    if (reads.count(ovl[oo].b_iid) == 0)
      continue;

    BestEdgeOverlap   edge(ovl[oo]);
    ufNode            read;

    if (tig->placeRead(read, op.frgID, ovl[oo].AEndIs3prime(), &edge) == false) {
      if (logFileFlagSet(LOG_PLACE_READ))
        writeLog("pRUO()-- WARNING: Failed to place with overlap %u %u hangs %d %d flipped %u\n",
                 ovl[oo].a_iid, ovl[oo].b_iid, ovl[oo].a_hang, ovl[oo].b_hang, ovl[oo].flipped);
      continue;
    }

    if (((flags & placeRead_fullMatch) ||
         (flags & placeRead_noExtend)) &&
        ((read.position.min() < 0) ||
         (read.position.max() > tig->getLength())))
      continue;

    //writeLog("NEW place %d-%d\n", read.position.bgn, read.position.end);

    bgnPos.insert(read.position.bgn);
    endPos.insert(read.position.end);
  }

  bgnPos.finalize();
  endPos.finalize();

  op.position.bgn = bgnPos.mean();
  op.position.end = endPos.mean();

  placeRead_computeVerified(op, tig->getLength());
}




bool
placeReadUsingOverlaps(TigVector                &tigs,
                       Unitig                   *target,
                       uint32                    fid,
                       vector<overlapPlacement> &placements,
                       uint32                    flags) {

  //if ((fid == 232074) || (fid == 72374) || (fid == 482602))
  //  logFileFlags |= LOG_PLACE_READ;

  if (logFileFlagSet(LOG_PLACE_READ))  //  Nope, not ambiguous.
    if (target)
      writeLog("\npRUO()-- begin for read %d into target tig %d\n", fid, target->id());
    else
      writeLog("\npRUO()-- begin for read %d into all tigs\n", fid);

  assert(fid > 0);
  assert(fid <= RI->numReads());

  //  Grab overlaps we'll use to place this read.

  uint32                ovlLen = 0;
  BAToverlap           *ovl    = OC->getOverlaps(fid, ovlLen);

  //  Grab some work space, and clear the output.

  placements.clear();

  //  Compute placements.  Anything that doesn't get placed is left as 'nowhere', specifically, in
  //  unitig 0 (which doesn't exist).

  overlapPlacement *ovlPlace = placeRead_fromOverlaps(tigs, target, fid, flags, ovlLen, ovl);

  //  We've placed the read in all possible places, or set unitig ID to 0 (an invalid unitig).
  //  Sort all the placements.  Sort order is:
  //    unitig ID (so zero is first)
  //    placed orientation (reverse is first)
  //    position

  sort(ovlPlace, ovlPlace + ovlLen, overlapPlacement_byLocation);


  //  Segregate the overlaps by placement in the unitig.  We want to construct one
  //  overlapPlacement for each distinct placement.  How this is done:
  //
  //  For all overlapping overlaps for a specific unitig and a specific read orientation, the
  //  end points (+- a few bases) are added to an interval list.  The intervals are combined.
  //  Each combined interval now forms the basis of a cluster of overlaps.  A list of the pairs of
  //  clusters hit by each overlap is built.  If there is one clear winner, that is picked.  If
  //  there is no clear winner, the read cannot be placed.
  //
  //  unitig    ==========================================================================
  //  overlaps       --------------------                           ------------------
  //                          -------------------                   -------------------
  //                          -------------------                   -------------------
  //  intervals     x1x      x2x       x3x     x4x                 x5x              x6xx
  //
  //  This unitig has two sets of "overlapping overlaps".  The left set could be from a tandem
  //  repeat.  We'll get two overlaps to the 2,4 pair, and one overlap to the 1,3 pair.  Assuming
  //  that is good enough to be a clear winner, we'll ignore the 1,3 overlap and compute position
  //  based on the other two overlaps.

  uint32         bgn = 0;  //  Range of overlaps with the same unitig/orientation
  uint32         end = 1;

  //  Skip overlaps that didn't generate a placement

  while ((bgn < ovlLen) && (ovlPlace[bgn].tigID == 0))
    bgn++;

  //  Process all placements.

  while (bgn < ovlLen) {

    //  Find the last placement with the same unitig/orientation as the 'bgn' read.

    end = bgn + 1;
    while ((end < ovlLen) &&
           (ovlPlace[bgn].tigID == ovlPlace[end].tigID) &&
           (ovlPlace[bgn].position.isReverse() == ovlPlace[end].position.isReverse()))
      end++;

    if (logFileFlagSet(LOG_PLACE_READ))
      writeLog("\nplaceReadUsingOverlaps()-- Merging placements %u to %u to place the read.\n", bgn, end);

    //  Build interval lists for the begin point and the end point.  Remember, this is all reads
    //  to a single unitig (the whole picture above), not just the overlapping read sets (left
    //  or right blocks).

    intervalList<int32>   bgnPoints;
    intervalList<int32>   endPoints;

    placeRead_assignEndPointsToCluster(bgn, end, fid, ovlPlace, bgnPoints, endPoints);

    //  Now, assign each placement to an end-pair cluster based on the interval ID that the end point falls in.
    //
    //  Count the number of reads that hit each pair of points.  Assign each ovlPlace to an implicit
    //  numbering of each pair of points.

    placeRead_assignPlacementsToCluster(bgn, end, fid, ovlPlace, bgnPoints, endPoints);

    //  Sort these placements by their clusterID.

    sort(ovlPlace + bgn, ovlPlace + end, overlapPlacement_byCluster);

    //  Run through each 'cluster' and compute a final placement for the read.
    //    A cluster extends from placements os to oe.
    //    Each cluster generates one placement.

    for (uint32 os=bgn, oe=bgn+1; os<end; ) {

      //  Find the end ovlPlace, oe, for this cluster, and do a quick check on orientation.

      for (oe=os+1; (oe < end) && (ovlPlace[os].clusterID == ovlPlace[oe].clusterID); oe++) {
        assert(ovlPlace[os].tigID                == ovlPlace[oe].tigID);
        assert(ovlPlace[os].position.isReverse() == ovlPlace[oe].position.isReverse());
      }

      //  Build the set of reads we care about.

#ifdef TEST_ALT
      set<uint32> reads;

      for (uint32 oo=os; oo<oe; oo++)
        reads.insert(ovlPlace[oo].refID);
#endif

      //  Make a new overlapPlacement from the first placement in this cluster.

      if (logFileFlagSet(LOG_PLACE_READ))
        writeLog("pRUO()-- process clusterID %u\n", ovlPlace[os].clusterID);

      overlapPlacement  op(fid, ovlPlace[os]);

      //  Find the first and last read in the unitig that we overlap with.

      placeRead_findFirstLastOverlapping(op, tigs[op.tigID], os, oe, ovlPlace);

      //  Sum the errors and bases aligned for each overlap.
      //  Find the minimum and maximum coordinates covered in the read, use that to compute the
      //  fraction of read coverage.

      placeRead_computeQualityAndCoverage(op, os, oe, ovlPlace);

#ifdef TEST_ALT
      //  Test the alternate qual and cov compute that uses overlaps directly
      {
        double er = op.errors;
        uint32 al = op.aligned;
        double fC = op.fCoverage;
        int32  bg = op.covered.bgn;
        int32  ed = op.covered.end;

        placeRead_computeQualityAndCoverage(op, ovl, ovlLen, reads);

        if ((er - op.errors > 0.0001) ||
            ((int32)al - (int32)op.aligned != 0) ||
            (fC - op.fCoverage > 0.0001) ||
            (bg - op.covered.bgn != 0) ||
            (ed - op.covered.end != 0))
        writeLog("COMPARE er %8.3f %8.3f %8.3f al %7u %7u %7d fC %8.4f %8.4f %8.4f bg %8d %8d %8d ed %8d %8d %8d\n",
                 er, op.errors, er - op.errors,
                 al, op.aligned, (int32)al - (int32)op.aligned,
                 fC, op.fCoverage, fC - op.fCoverage,
                 bg, op.covered.bgn, bg - op.covered.bgn,
                 ed, op.covered.end, ed - op.covered.end);
      }
#endif

      //  Compute placement based on the longest overlap on each end, or the best contain.

      placeRead_computePlacement(op, os, oe, ovlPlace, tigs[op.tigID]);

#ifdef TEST_ALT
      {
        SeqInterval  origpos = op.position;
        SeqInterval  origver = op.verified;

        placeRead_computePlacement(op, ovl, ovlLen, reads, flags, tigs[op.tigID]);

        if ((origpos.bgn - op.position.bgn > 10) ||   //  Placements wobble by a few bases
            (origpos.end - op.position.end > 10) ||
            (origver.bgn - op.verified.bgn > 10) ||
            (origver.end - op.verified.end > 10))
          writeLog("COMPARE pos bgn %d-%d end %d-%d  ver bgn %d-%d end %d-%d\n",
                   origpos.bgn, op.position.bgn,
                   origpos.end, op.position.end,
                   origver.bgn, op.verified.bgn,
                   origver.end, op.verified.end);
      }
#endif

      //  Filter out bogus placements.  There used to be a few more, but they made no sense for long reads.
      //  Reject if either end stddev is high.  It has to be pretty bad before this triggers.

      bool   goodPlacement   = true;

#if 0
      double allowableStdDev = max(2.0, 0.075 * RI->readLength(op.frgID));

      if ((bgnPos.stddev() > allowableStdDev) ||
          (endPos.stddev() > allowableStdDev))
        goodPlacement = false;
#endif

      if ((flags & placeRead_fullMatch) &&
          (op.fCoverage < 0.99))
        goodPlacement = false;

      if ((flags & placeRead_noExtend) &&
          ((op.position.min() < 0) ||
           (op.position.max() > tigs[op.tigID]->getLength())))
        goodPlacement = false;

      if (goodPlacement)
        placements.push_back(op);


      if (logFileFlagSet(LOG_PLACE_READ))
        writeLog("pRUO()--   placements[%u] - PLACE READ %d in tig %d at %d,%d -- verified %d,%d -- covered %d,%d %4.1f%% -- errors %.2f aligned %d novl %d%s\n",
                 placements.size() - 1,
                 op.frgID, op.tigID,
                 op.position.bgn, op.position.end,
                 op.verified.bgn, op.verified.end,
                 op.covered.bgn, op.covered.end,
                 op.fCoverage * 100.0,
                 op.errors, op.aligned, oe - os,
                 (goodPlacement == false) ? " -- INVALID"  : "");

      os = oe;
    }  //  End of segregating overlaps by placement

    //  Move to the next block of overlaps.
    bgn = end;
  }

  delete [] ovlPlace;

  //if ((fid == 232074) || (fid == 72374) || (fid == 482602))
  //  logFileFlags &= ~LOG_PLACE_READ;

  return(true);
}

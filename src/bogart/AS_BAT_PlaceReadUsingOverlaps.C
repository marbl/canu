
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
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_PlaceReadUsingOverlaps.H"

#include "intervalList.H"

#include <vector>
#include <algorithm>



void
placeRead_fromOverlaps(TigVector          &tigs,
                       Unitig             *target,
                       uint32              fid,
                       uint32              flags,
                       double              errorLimit,
                       uint32              ovlLen,
                       BAToverlap         *ovl,
                       uint32             &ovlPlaceLen,
                       overlapPlacement   *ovlPlace) {

  if (logFileFlagSet(LOG_PLACE_READ))
    writeLog("pRUO()-- placements for read %u with %u overlaps\n", fid, ovlLen);

  for (uint32 oo=0; oo<ovlLen; oo++) {
    bool              disallow = false;
    uint32            btID     = tigs.inUnitig(ovl[oo].b_iid);

    assert(ovl[oo].a_iid == fid);

    if ((btID == 0) ||                                  //  Skip if overlapping read isn't in a tig yet - unplaced contained, or garbage read.
        ((target != NULL) && (target->id() != btID)))   //  Skip if we requested a specific tig and if this isn't it.
      continue;

   if (ovl[oo].erate() > errorLimit)                    //  Skip if the error is higher than we were told to use
      continue;

    Unitig           *btig   = tigs[btID];
    ufNode           &bread  = btig->ufpath[ tigs.ufpathIdx(ovl[oo].b_iid) ];

    if (btig->_isUnassembled == true)  //  Skip if overlapping read is in an unassembled contig.
      continue;

    SeqInterval       apos;   //  Position of the read in the btig
    SeqInterval       bver;   //  Bases covered by the overlap in the B tig

    //  Place the read relative to the other read.  The overlap we have is relative to the A read,
    //  so hangs need to be subtracted from the other coordinate.
    //
    //  Pictures all show positive hangs.

    //   A  ------------>     (b)
    //   B  (a)     ------------>
    if ((ovl[oo].flipped == false) && (bread.position.isForward() == true)) {
      apos.bgn = bread.position.min() - ovl[oo].a_hang;
      apos.end = bread.position.max() - ovl[oo].b_hang;

      bver.bgn = bread.position.min() - ((ovl[oo].a_hang > 0) ? 0 : ovl[oo].a_hang);
      bver.end = bread.position.max() - ((ovl[oo].b_hang > 0) ? ovl[oo].b_hang : 0);
    }

    //   A  ------------>     (b)
    //   B  (a)    <-------------
    if ((ovl[oo].flipped == true) && (bread.position.isForward() == false)) {
      apos.bgn = bread.position.min() - ovl[oo].a_hang;
      apos.end = bread.position.max() - ovl[oo].b_hang;

      bver.bgn = bread.position.min() - ((ovl[oo].a_hang > 0) ? 0 : ovl[oo].a_hang);
      bver.end = bread.position.max() - ((ovl[oo].b_hang > 0) ? ovl[oo].b_hang : 0);
    }

    //   A  (b)     <------------
    //   B  ------------>     (a)
    if ((ovl[oo].flipped == true) && (bread.position.isForward() == true)) {
      apos.end = bread.position.min() + ovl[oo].b_hang;
      apos.bgn = bread.position.max() + ovl[oo].a_hang;

      bver.end = bread.position.min() + ((ovl[oo].b_hang > 0) ? ovl[oo].b_hang : 0);
      bver.bgn = bread.position.max() + ((ovl[oo].a_hang > 0) ? 0 : ovl[oo].a_hang);
    }

    //   A  (b)     <------------
    //   B  <------------     (a)
    if ((ovl[oo].flipped == false) && (bread.position.isForward() == false)) {
      apos.end = bread.position.min() + ovl[oo].b_hang;
      apos.bgn = bread.position.max() + ovl[oo].a_hang;

      bver.end = bread.position.min() + ((ovl[oo].b_hang > 0) ? ovl[oo].b_hang : 0);
      bver.bgn = bread.position.max() + ((ovl[oo].a_hang > 0) ? 0 : ovl[oo].a_hang);
    }

    //  HOWEVER, the verified position is all goobered up if the overlapping read
    //  was placed too short.  Imagine a 20k read with a 500bp overlap, so the hangs
    //  are 19.5k.  If we position this read 1k too short, then readLen-hang is negative,
    //  and we end up misorienting the verified coords (not too mention that they're
    //  likely bogus too).  So, if that happens, we just ignore the overlap.

    int32             bposlen = (bread.position.max() - bread.position.min());

    if (ovl[oo].a_hang < 0)
      bposlen += ovl[oo].a_hang;

    if (ovl[oo].b_hang > 0)
      bposlen -= ovl[oo].b_hang;

    if (bposlen < 0) {
      writeLog("WARNING: read %u overlap to read %u in tig %u at %d-%d - hangs %d %d to large for placement, ignoring overlap\n",
               ovl[oo].a_iid,
               ovl[oo].b_iid,
               btID,
               bread.position.bgn, bread.position.end,
               ovl[oo].a_hang, ovl[oo].b_hang);
      disallow = true;
    }

    //  Save the placement in our work space.

    uint32  flen = RI->readLength(ovl[oo].a_iid);

    overlapPlacement  op;

    op.frgID        = fid;
    op.refID        = ovl[oo].b_iid;
    op.tigID        = btig->id();
    op.position     = apos;
    op.verified     = bver;
    op.covered.bgn  = (ovl[oo].a_hang < 0) ?    0 : ovl[oo].a_hang;          //  The portion of the read
    op.covered.end  = (ovl[oo].b_hang > 0) ? flen : ovl[oo].b_hang + flen;   //  covered by the overlap.
    op.clusterID    = 0;
    op.fCoverage    = 0.0;
    op.errors       = RI->overlapLength(ovl[oo].a_iid, ovl[oo].b_iid, ovl[oo].a_hang, ovl[oo].b_hang) * ovl[oo].erate();
    op.aligned      = op.covered.end - op.covered.bgn;
    op.tigFidx      = UINT32_MAX;
    op.tigLidx      = 0;

    //  If we're looking for final placments either contained completely in the tig or covering the
    //  whole read, disallow any placements that exceed the boundary of the unitig.  This is NOT a
    //  filter on these placements; but any placement here that extends past the end of the tig is
    //  guaranteed to not generate a contained/whole-read placement.

    if ((flags & placeRead_noExtend) || (flags & placeRead_fullMatch))
      if ((op.position.min() < 0) ||
          (op.position.max() > btig->getLength()))
        disallow = true;

    if (logFileFlagSet(LOG_PLACE_READ))
      writeLog("pRUO()-- bases %5d-%-5d to tig %5d %8ubp at %8d-%-8d olap %8d-%-8d via read %7d at %8d-%-8d hang %6d %6d erate %.4f %s%s\n",
               op.covered.bgn, op.covered.end,
               btig->id(),
               btig->getLength(),
               op.position.bgn, op.position.end,
               op.verified.bgn, op.verified.end,
               bread.ident,
               bread.position.bgn,
               bread.position.end,
               ovl[oo].a_hang, ovl[oo].b_hang,
               ovl[oo].erate(),
               (ovl[oo].flipped == true) ? "I" : "N",
               (disallow) ? " DISALLOW" : "");

    //  Ensure everything is hunkey dorey and save the overlap.

    if (disallow == false) {
      assert(op.covered.bgn >= 0);
      assert(op.covered.end <= flen);
      assert(op.covered.isForward() == true);

      assert(op.position.isForward() == op.verified.isForward());

      if (op.position.isForward() == true) {
        assert(op.position.bgn <= op.verified.bgn);
        assert(op.verified.end <= op.position.end);
      } else {
        assert(op.position.end <= op.verified.end);
        assert(op.verified.bgn <= op.position.bgn);
      }

      ovlPlace[ovlPlaceLen++] = op;
    }
  }
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
                                   uint32 os, uint32 oe,
                                   overlapPlacement *ovlPlace,
                                   Unitig *tig) {
  op.tigFidx = UINT32_MAX;
  op.tigLidx = 0;

  for (uint32 oo=os; oo<oe; oo++) {
    uint32   ord = tig->ufpathIdx(ovlPlace[oo].refID);

    op.tigFidx = min(ord, op.tigFidx);
    op.tigLidx = max(ord, op.tigLidx);
  }

  if (op.tigFidx > op.tigLidx)
    writeStatus("pRUO()-- Invalid placement indices: tigFidx %u tigLidx %u\n", op.tigFidx, op.tigLidx);
  assert(op.tigFidx <= op.tigLidx);

  if (logFileFlagSet(LOG_PLACE_READ))
    writeLog("pRUO()--  spans reads #%u (%u) to #%u (%u) in tig %u\n",
             op.tigFidx, tig->ufpath[op.tigFidx].ident,
             op.tigLidx, tig->ufpath[op.tigLidx].ident,
             op.tigID);
}



void
placeRead_computePlacement(overlapPlacement &op,
                           uint32            os,
                           uint32            oe,
                           overlapPlacement *ovlPlace,
                           Unitig           *tig) {
  stdDev<double>   bgnPos, endPos;

  bool             isFwd   = ovlPlace[os].position.isForward();
  int32            readLen = RI->readLength(op.frgID);
  int32            tigLen  = tig->getLength();

  op.errors  = 0;
  op.aligned = 0;

  op.verified.bgn = (isFwd) ? INT32_MAX : INT32_MIN;
  op.verified.end = (isFwd) ? INT32_MIN : INT32_MAX;

  int32            bgnVer2 = (isFwd) ? INT32_MAX : INT32_MIN;;
  int32            endVer2 = (isFwd) ? INT32_MIN : INT32_MAX;;

  int32            bgnVer3 = (isFwd) ? INT32_MAX : INT32_MIN;;
  int32            endVer3 = (isFwd) ? INT32_MIN : INT32_MAX;;

  op.covered.bgn = INT32_MAX;  //  Covered interval is always in
  op.covered.end = INT32_MIN;  //  forward read coordinates

  //  Deleted overlaps?  From where?
  for (uint32 oo=os; oo<oe; oo++)
    assert((ovlPlace[oo].position.bgn != 0) ||
           (ovlPlace[oo].position.end != 0));

  if (logFileFlagSet(LOG_PLACE_READ))
    writeLog("pRUO()-- compute placement for os=%u od=%u\n", os, oe);

  //  Over all the placements that support this position:
  //    compute the final position as the mean of the supporting overlaps.
  //    compute the verified position as.....
  //    compute the read bases covered by an overlap as the min/max.
  //
  //  The verified position is a bit annoying.
  //
  //    The first attempt used the mean just as for the position.  But this occasionally
  //    left the verified outside the placed position.  It was thresholded to make it sane.
  //    IT ALSO TOTALLY BREAKS GFA EDGE FINDING.  (I think because the verified overlap position
  //    is too small).
  //
  //    The second attempt set it relative to the position, using the hangs from the
  //    'covered' position on the read.  This failed on am ~8k read placed with a 500bp
  //    overlap.  The overlapping read was placed shorter than expected.  The sum of the overlap
  //    hangs was larger than this placement.  (Largely solved by recomputing positions
  //    after unplaced reads are placed).
  //
  //    The third attempt mirrors what is done for 'covered' -- just take the min/max
  //    of all the overlaps used when placing the read.

  for (uint32 oo=os; oo<oe; oo++) {
    bgnPos.insert(ovlPlace[oo].position.bgn);
    endPos.insert(ovlPlace[oo].position.end);

    //  Third attempt

    if (logFileFlagSet(LOG_PLACE_READ))
      writeLog("pRUO()-- op %3d ovl ver %12d %12d pos %12d %12d\n",
              oo,
              ovlPlace[oo].verified.bgn, ovlPlace[oo].verified.end,
              ovlPlace[oo].position.bgn, ovlPlace[oo].position.end);

#if 1
    if (isFwd) {
      bgnVer3 = min(bgnVer3, ovlPlace[oo].verified.bgn);
      endVer3 = max(endVer3, ovlPlace[oo].verified.end);
    } else {
      bgnVer3 = max(bgnVer3, ovlPlace[oo].verified.bgn);
      endVer3 = min(endVer3, ovlPlace[oo].verified.end);
    }
#endif

    op.errors      += ovlPlace[oo].errors;
    op.aligned     += ovlPlace[oo].aligned;

    op.covered.bgn  = min(op.covered.bgn, ovlPlace[oo].covered.bgn);
    op.covered.end  = max(op.covered.end, ovlPlace[oo].covered.end);
  }

  op.fCoverage = (op.covered.end - op.covered.bgn) / (double)readLen;

  //  Take the mean of the positions as the final position.

  bgnPos.finalize();
  endPos.finalize();

  op.position.bgn = bgnPos.mean();
  op.position.end = endPos.mean();

  //  Second attempt.

#if 1
  if (isFwd) {
    bgnVer2 = op.position.bgn +            op.covered.bgn;
    endVer2 = op.position.end - (readLen - op.covered.end);
  } else {
    bgnVer2 = op.position.bgn -            op.covered.bgn;
    endVer2 = op.position.end + (readLen - op.covered.end);
  }
#endif

  if (logFileFlagSet(LOG_PLACE_READ))
    writeLog("pRUO()-- position %d-%d verified %d-%d %d-%d\n",
             op.position.bgn, op.position.end,
             bgnVer2, endVer2,
             bgnVer3, endVer3);

  //  Results in about 15% fewer contig and 3% fewer unitig edges, compared to v3.
#if 0
  op.verified.bgn = bgnVer2;
  op.verified.end = endVer2;
#endif

  //  On dmel, gives more contig edges than v2, mostly small stuff.
#if 1
  op.verified.bgn = bgnVer3;
  op.verified.end = endVer3;
#endif

  //  Finally, limit verified to be the extent of the tig, or the extent of the placement.

  if (isFwd) {
    if (op.verified.bgn < 0)                  op.verified.bgn = 0;
    if (op.verified.end > tigLen)             op.verified.end = tigLen;
  } else {
    if (op.verified.bgn > tigLen)             op.verified.bgn = tigLen;
    if (op.verified.end < 0)                  op.verified.end = 0;
  }

  if (isFwd) {
    if (op.verified.bgn < op.position.bgn)    op.verified.bgn = op.position.bgn;
    if (op.verified.end > op.position.end)    op.verified.end = op.position.end;
  } else {
    if (op.verified.bgn > op.position.bgn)    op.verified.bgn = op.position.bgn;
    if (op.verified.end < op.position.end)    op.verified.end = op.position.end;
  }

  //  And check that the result is sane.

  assert(op.position.isForward() == isFwd);
  assert(op.position.isForward() == op.verified.isForward());
  assert(op.covered.isForward() == true);

  if (isFwd) {
    assert(op.position.bgn <= op.verified.bgn);
    assert(op.verified.end <= op.position.end);
  } else {
    assert(op.position.end <= op.verified.end);
    assert(op.verified.bgn <= op.position.bgn);
  }
}



//  Compute how much of the read is covered by overlaps to the contig.
//
//  There's a bit of dead code in placeRead_computePlacement() above, that computes this coverage
//  using the extent of overlaps, but ignoring any gaps in coverage in the middle.
//
void
placeRead_computeCoverage(overlapPlacement &op,
                          uint32            os,
                          uint32            oe,
                          overlapPlacement *ovlPlace,
                          Unitig           *tig) {
  intervalList<int32>    readCov;

  //  Recompute op.covered, for no good reason except that the computation above should be removed.

  op.covered.bgn = INT32_MAX;  //  Covered interval is always in
  op.covered.end = INT32_MIN;  //  forward read coordinates

  for (uint32 oo=os; oo<oe; oo++) {
    op.covered.bgn  = min(op.covered.bgn, ovlPlace[oo].covered.bgn);
    op.covered.end  = max(op.covered.end, ovlPlace[oo].covered.end);

    if (logFileFlagSet(LOG_PLACE_READ))
      writeLog("pRUO()-- op %3d covers %8d-%-8d extent %8d-%-8d\n",
               oo,
               ovlPlace[oo].covered.bgn, ovlPlace[oo].covered.end,
               op.covered.bgn, op.covered.end);

    readCov.add(ovlPlace[oo].covered.bgn, ovlPlace[oo].covered.end - ovlPlace[oo].covered.bgn);
  }

  assert(op.covered.isForward() == true);

  readCov.merge();

  int32  readLen = RI->readLength(op.frgID);
  int32  bCov    = 0;

  for (uint32 ii=0; ii<readCov.numberOfIntervals(); ii++)
    bCov += readCov.hi(ii) - readCov.lo(ii);

  double  eCov = (op.covered.end - op.covered.bgn) / (double)readLen;
  double  fCov = (bCov)                            / (double)readLen;

  op.fCoverage = fCov;

  if (logFileFlagSet(LOG_PLACE_READ))
    writeLog("pRUO()-- covered %6.4f extent %6.4f\n",
             fCov, eCov);
}



bool
placeReadUsingOverlaps(TigVector                &tigs,
                       Unitig                   *target,
                       uint32                    fid,
                       vector<overlapPlacement> &placements,
                       uint32                    flags,
                       double                    errorLimit) {

  set<uint32>  verboseEnable;

  //  Enable logging for all reads, or specific tigs, or not at all.
  //verboseEnable.insert(fid);
  //
  //for (uint32 fi=0; (tigs[   2]) && (fi<tigs[   2]->ufpath.size()); fi++)   verboseEnable.insert(tigs[   2]->ufpath[fi].ident);
  //for (uint32 fi=0; (tigs[   3]) && (fi<tigs[   3]->ufpath.size()); fi++)   verboseEnable.insert(tigs[   3]->ufpath[fi].ident);
  //for (uint32 fi=0; (tigs[   4]) && (fi<tigs[   4]->ufpath.size()); fi++)   verboseEnable.insert(tigs[   4]->ufpath[fi].ident);
  //for (uint32 fi=0; (tigs[ 255]) && (fi<tigs[ 255]->ufpath.size()); fi++)   verboseEnable.insert(tigs[ 255]->ufpath[fi].ident);

  if (verboseEnable.count(fid) > 0)
    logFileFlags |= LOG_PLACE_READ;

  if (logFileFlagSet(LOG_PLACE_READ))  //  Nope, not ambiguous.
    if (target)
      writeLog("\npRUO()-- begin for read %u length %u into target tig %d\n", fid, RI->readLength(fid), target->id());
    else
      writeLog("\npRUO()-- begin for read %u length %u into all tigs\n", fid, RI->readLength(fid));

  assert(fid > 0);
  assert(fid <= RI->numReads());

  //  Grab overlaps we'll use to place this read.

  uint32                ovlLen = 0;
  BAToverlap           *ovl    = OC->getOverlaps(fid, ovlLen);

  //  Grab some work space, and clear the output.

  placements.clear();

  //  Compute placements.  Anything that doesn't get placed is left as 'nowhere', specifically, in
  //  unitig 0 (which doesn't exist).

  uint32             ovlPlaceLen = 0;
  overlapPlacement  *ovlPlace    = new overlapPlacement [ovlLen];

  placeRead_fromOverlaps(tigs, target, fid, flags, errorLimit, ovlLen, ovl, ovlPlaceLen, ovlPlace);

  //  Sort all the placements.  Sort order is:
  //    unitig ID
  //    placed orientation (reverse is first)
  //    position

  sort(ovlPlace, ovlPlace + ovlPlaceLen, overlapPlacement_byLocation);

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

  while (bgn < ovlPlaceLen) {

    //  Find the last placement with the same unitig/orientation as the 'bgn' read.

    end = bgn + 1;
    while ((end < ovlPlaceLen) &&
           (ovlPlace[bgn].tigID == ovlPlace[end].tigID) &&
           (ovlPlace[bgn].position.isReverse() == ovlPlace[end].position.isReverse()))
      end++;

    if (logFileFlagSet(LOG_PLACE_READ))
      writeLog("\npRUO()-- Merging placements %u to %u to place the read.\n", bgn, end);

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
      if (logFileFlagSet(LOG_PLACE_READ))
        writeLog("pRUO()-- process clusterID %u\n", ovlPlace[os].clusterID);

      //  Find the end ovlPlace, oe, for this cluster, and do a quick check on orientation.

      for (oe=os+1; (oe < end) && (ovlPlace[os].clusterID == ovlPlace[oe].clusterID); oe++) {
        assert(ovlPlace[os].tigID                == ovlPlace[oe].tigID);
        assert(ovlPlace[os].position.isReverse() == ovlPlace[oe].position.isReverse());
      }

      //  Make a new overlapPlacement from the first placement in this cluster, figure out the first/last tig reads that
      //  have overlaps to it, and figure out final positions.

      overlapPlacement  op(fid, ovlPlace[os]);

      placeRead_findFirstLastOverlapping(op, os, oe, ovlPlace, tigs[op.tigID]);
      placeRead_computePlacement        (op, os, oe, ovlPlace, tigs[op.tigID]);
      placeRead_computeCoverage         (op, os, oe, ovlPlace, tigs[op.tigID]);

      //  Filter out bogus placements.  There used to be a few more, but they made no sense for long reads.
      //  Reject if either end stddev is high.  It has to be pretty bad before this triggers.

      bool   fullMatch = true;
      bool   noExtend  = true;

      if ((flags & placeRead_fullMatch) &&
          (op.fCoverage < 1.0))
        fullMatch = false;

      if ((flags & placeRead_noExtend) &&
          ((op.position.min() < 0) ||
           (op.position.max() > tigs[op.tigID]->getLength())))
        noExtend = false;

      if ((fullMatch == true) &&
          (noExtend  == true)) {
        placements.push_back(op);

        if (logFileFlagSet(LOG_PLACE_READ))
          writeLog("pRUO()--   placements[%u] - PLACE READ %d in tig %d at %d,%d -- verified %d,%d -- covered %d,%d %4.1f%% -- errors %.2f aligned %d novl %d\n",
                   placements.size() - 1,
                   op.frgID, op.tigID,
                   op.position.bgn, op.position.end,
                   op.verified.bgn, op.verified.end,
                   op.covered.bgn, op.covered.end,
                   op.fCoverage * 100.0,
                   op.errors, op.aligned, oe - os);
      } else {
        if (logFileFlagSet(LOG_PLACE_READ))
          writeLog("pRUO()--   placements[%u] - DO NOT PLACE READ %d in tig %d at %d,%d -- verified %d,%d -- covered %d,%d %4.1f%% -- errors %.2f aligned %d novl %d%s%s\n",
                   placements.size() - 1,
                   op.frgID, op.tigID,
                   op.position.bgn, op.position.end,
                   op.verified.bgn, op.verified.end,
                   op.covered.bgn, op.covered.end,
                   op.fCoverage * 100.0,
                   op.errors, op.aligned, oe - os,
                   (fullMatch == false) ? " -- PARTIAL"  : "",
                   (noExtend  == false) ? " -- EXTENDS"  : "");
      }


      os = oe;
    }  //  End of segregating overlaps by placement

    //  Move to the next block of overlaps.
    bgn = end;
  }

  delete [] ovlPlace;

  if (verboseEnable.count(fid) > 0)
    logFileFlags &= ~LOG_PLACE_READ;

  return(true);
}

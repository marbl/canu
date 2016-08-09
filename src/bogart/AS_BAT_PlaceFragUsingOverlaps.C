
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
 *    src/AS_BAT/AS_BAT_PlaceFragUsingOverlaps.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-SEP-08
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2015-AUG-05
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

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "intervalList.H"



bool
placeFragUsingOverlaps(TigVector             &unitigs,
                       Unitig                   *target,
                       uint32                    fid,
                       vector<overlapPlacement> &placements,
                       uint32                    flags) {

  //if (fid == 328)
  //  logFileFlags |= LOG_PLACE_FRAG;

  if (logFileFlagSet(LOG_PLACE_FRAG))  //  Nope, not ambiguous.
    if (target)
      writeLog("\npFUO()-- begin for frag %d into target tig %d\n", fid, target->id());
    else
      writeLog("\npFUO()-- begin for frag %d into all tigs\n", fid);

  assert(fid > 0);
  assert(fid <= FI->numFragments());

  //  Grab overlaps we'll use to place this read.

  uint32                ovlLen = 0;
  BAToverlap           *ovl    = OC->getOverlaps(fid, AS_MAX_ERATE, ovlLen);

  //  Grab some work space, and clear the output.

  overlapPlacement     *ovlPlace = new overlapPlacement[ovlLen];

  placements.clear();

  //  Compute placements.
  //  Anything that doesn't get placed is left as 'nowhere', specifically, in unitig 0 (which doesn't exist).

  for (uint32 i=0; i<ovlLen; i++) {
    int32             tigID = Unitig::fragIn(ovl[i].b_iid);
    Unitig           *tig   = unitigs[tigID];

    assert(ovl[i].a_iid == fid);

    if (tigID == 0)                           //  Skip if overlapping read isn't in a tig yet - unplaced contained, or garbage read.
      continue;

    if ((target != NULL) && (target != tig))  //  Skip if we requested a specific tig and if this isn't it.
      continue;

    //  Place the fragment relative to the other fragment.

    BestEdgeOverlap   edge(ovl[i]);
    ufNode            frag;

    if (tig->placeFrag(frag, fid, ovl[i].AEndIs3prime(), &edge) == false) {
      if (logFileFlagSet(LOG_PLACE_FRAG))
        writeLog("pFUO()-- WARNING: Failed to place with overlap %u %u hangs %u %u flipped %u\n",
                 ovl[i].a_iid, ovl[i].b_iid, ovl[i].a_hang, ovl[i].b_hang, ovl[i].flipped);
      continue;
    }

    //  Save the placement in our work space.

    uint32  olen = FI->overlapLength(ovl[i].a_iid, ovl[i].b_iid, ovl[i].a_hang, ovl[i].b_hang);
    uint32  flen = FI->fragmentLength(ovl[i].a_iid);

    ovlPlace[i].frgID        = fid;
    ovlPlace[i].refID        = ovl[i].b_iid;
    ovlPlace[i].tigID        = tig->id();
    ovlPlace[i].position     = frag.position;
    ovlPlace[i].verified.bgn = INT32_MAX;
    ovlPlace[i].verified.end = INT32_MIN;
    ovlPlace[i].covered.bgn  = (ovl[i].a_hang < 0) ?    0 : ovl[i].a_hang;          //  The portion of the read
    ovlPlace[i].covered.end  = (ovl[i].b_hang > 0) ? flen : ovl[i].b_hang + flen;   //  covered by the overlap.
    ovlPlace[i].bgnStdDev    = 0.0;
    ovlPlace[i].endStdDev    = 0.0;
    ovlPlace[i].clusterID    = 0;
    ovlPlace[i].fCoverage    = 0.0;
    ovlPlace[i].errors       = olen * ovl[i].erate;
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
    //    o  sticking a unique/repeat fragment onto a repeat (leaving the unique uncovered)
    //    o  sticking a chimeric fragment onto the end of a unitig (leaving the chimeric join uncovered)

    if (((flags & placeFrag_fullMatch) ||
         (flags & placeFrag_noExtend)) &&
        ((ovlPlace[i].position.min() < 0) ||
         (ovlPlace[i].position.max() > tig->getLength()))) {
      ovlPlace[i] = overlapPlacement();
    }

    //  Report the placement.


    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("pFUO()-- frag %d in unitig %d at %d,%d (covered %d,%d) from overlap ident %d %d hang %d %d flipped %d%s\n",
               ovlPlace[i].frgID,
               ovlPlace[i].tigID,
               ovlPlace[i].position.bgn, ovlPlace[i].position.end,
               ovlPlace[i].covered.bgn, ovlPlace[i].covered.end,
               ovl[i].a_iid, ovl[i].b_iid, ovl[i].a_hang, ovl[i].b_hang, ovl[i].flipped,
               (ovlPlace[i].frgID == 0) ? " DISALLOWED" : "");
  }  //  Over all overlaps.


  //  We've placed the read in all possible places, or set unitig ID to 0 (an invalid unitig).
  //  Sort all the placements.  Sort order is:
  //    unitig ID (so zero is first)
  //    placed orientation (reverse is first)
  //    position

  sort(ovlPlace, ovlPlace + ovlLen, overlapPlacement_byLocation);



  //  Segregate the overlaps by placement in the unitig.  We want to construct one
  //  overlapPlacement for each distinct placement.  How this is done:
  //
  //  For all overlapping overlaps for a specific unitig and a specific fragment orientation, the
  //  end points (+- a few bases) are added to an interval list.  The intervals are combined.
  //  Each combined interval now forms the basis of a cluster of overlaps.  A list of the pairs of
  //  clusters hit by each overlap is built.  If there is one clear winner, that is picked.  If
  //  there is no clear winner, the fragment cannot be placed.
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

    //  Find the last placement with the same unitig/orientation as the 'bgn' fragment.

    end = bgn + 1;
    while ((end < ovlLen) &&
           (ovlPlace[bgn].tigID == ovlPlace[end].tigID) &&
           (ovlPlace[bgn].position.isReverse() == ovlPlace[end].position.isReverse()))
      end++;

    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("pFUO()-- Merging placements %u to %u to place the read.\n", bgn, end);

    //  Build interval lists for the begin point and the end point.  Remember, this is all fragments
    //  to a single unitig (the whole picture above), not just the overlapping fragment sets (left
    //  or right blocks).

    intervalList<int32>   bgnPoints;
    intervalList<int32>   endPoints;

    int32                 windowSlop = 0.075 * FI->fragmentLength(fid);

    if (windowSlop < 5)
      windowSlop = 5;

    for (uint32 oo=bgn; oo<end; oo++) {
      bgnPoints.add(ovlPlace[oo].position.bgn - windowSlop, 2 * windowSlop);
      endPoints.add(ovlPlace[oo].position.end - windowSlop, 2 * windowSlop);
    }

    bgnPoints.merge();
    endPoints.merge();

    if (logFileFlagSet(LOG_PLACE_FRAG)) {
      writeLog("pFUO()-- Found %u bgn intervals: ", bgnPoints.numberOfIntervals());
      for (uint32 r=0; r<bgnPoints.numberOfIntervals(); r++)
        writeLog(" %u-%u", bgnPoints.lo(r), bgnPoints.hi(r));
      writeLog("\n");

      writeLog("pFUO()-- Found %u end intervals: ", endPoints.numberOfIntervals());
      for (uint32 r=0; r<endPoints.numberOfIntervals(); r++)
        writeLog(" %u-%u", endPoints.lo(r), endPoints.hi(r));
      writeLog("\n");
    }

    //  Now, assign each placement to an end-pair cluster based on the interval ID that the end point falls in.
    //
    //  Count the number of fragments that hit each pair of points.  Assign each ovlPlace to an implicit
    //  numbering of each pair of points.

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

    //  Sort these placements by their clusterID.

    sort(ovlPlace + bgn, ovlPlace + end, overlapPlacement_byCluster);

    //  Run through each 'cluster' and compute a final placement for the read.
    //    A cluster extends from placements os to oe.
    //    Each cluster generates one placement.

    for (uint32 os=bgn, oe=bgn; os<end; ) {
      if (logFileFlagSet(LOG_PLACE_FRAG))
        writeLog("pFUO()-- process clusterID %u\n", ovlPlace[os].clusterID);

      overlapPlacement  op;

      //  Initialize a new op.

      op.frgID          = fid;

      op.refID          = UINT32_MAX;              //  Not valid in the output overlapPlacement.
      op.tigID          = ovlPlace[os].tigID;
      op.position.bgn   = 0;
      op.position.end   = 0;

      op.verified.bgn   = 0;
      op.verified.end   = 0;

      op.covered.bgn    = ovlPlace[os].covered.bgn;
      op.covered.end    = ovlPlace[os].covered.end;

      op.bgnStdDev      = 0.0;
      op.endStdDev      = 0.0;

      op.clusterID      = ovlPlace[os].clusterID;  //  Useless to track forward.

      op.fCoverage      = 0.0;                     //  coverage of the fragment

      op.errors         = 0.0;                     //  sum of the estimated number of errors in all the overlaps
      op.aligned        = 0;                       //  bases aligned?

      op.tigFidx        = UINT32_MAX;
      op.tigLidx        = UINT32_MAX;

      //  Find the end ovlPlace for this cluster.
      //  Do a quick sanity check to make sure all placements are the same tig and the same orientation.

      for (oe=os; (oe < end) && (ovlPlace[os].clusterID == ovlPlace[oe].clusterID); oe++) {
        assert(ovlPlace[os].tigID                == ovlPlace[oe].tigID);
        assert(ovlPlace[os].position.isReverse() == ovlPlace[oe].position.isReverse());
      }

      //  Sum the errors and bases aligned for each overlap.
      //  Find the minimum and maximum coordinates covered in the read, use that to compute the fraction of read coverage.

      for (uint32 oo=os; oo<oe; oo++) {
        if ((ovlPlace[oo].position.bgn == 0) &&
            (ovlPlace[oo].position.end == 0))
          continue;

        op.errors      += ovlPlace[oo].errors;
        op.aligned     += ovlPlace[oo].aligned;

        //  Find min/max covered interval.  These are always forward on the read.

        op.covered.bgn  = min(op.covered.bgn, ovlPlace[oo].covered.bgn);
        op.covered.end  = max(op.covered.end, ovlPlace[oo].covered.end);
      }

      op.fCoverage = (op.covered.end - op.covered.bgn) / (double)FI->fragmentLength(op.frgID);

      //  Find the first and last fragment in the unitig that we overlap with.

      Unitig   *tig    = unitigs[op.tigID];
      uint32    tigLen = tig->getLength();

      op.tigFidx = UINT32_MAX;
      op.tigLidx = 0;

      for (uint32 oo=os; oo<oe; oo++) {
        uint32   ord = tig->pathPosition(ovlPlace[oo].refID);

        op.tigFidx = min(ord, op.tigFidx);
        op.tigLidx = max(ord, op.tigLidx);

        //if (logFileFlagSet(LOG_PLACE_FRAG))
        //  writeLog("pFUO()--     find range from os=%u to oe=%u  tig=%u  ord=%u  f=%u l=%u\n",
        //            os, oe, op.tigID, ord, op.tigFidx, op.tigLidx);
      }

      if (op.tigFidx > op.tigLidx)
        writeStatus("Invalid placement indices: tigFidx %u tigLidx %u\n", op.tigFidx, op.tigLidx);
      assert(op.tigFidx <= op.tigLidx);

      if (logFileFlagSet(LOG_PLACE_FRAG))
        writeLog("pFUO()--   spans reads #%u (%u) to #%u (%u) in tig %u\n",
                 op.tigFidx, tig->ufpath[op.tigFidx].ident,
                 op.tigLidx, tig->ufpath[op.tigLidx].ident,
                 op.tigID);

      //  Compute mean and stddev placement.

      stdDev<double>   bgnPos;
      stdDev<double>   endPos;

      for (uint32 oo=os; oo<oe; oo++) {
        if ((ovlPlace[oo].position.bgn == 0) &&
            (ovlPlace[oo].position.end == 0))
          continue;

        bgnPos.insert(ovlPlace[oo].position.bgn);
        endPos.insert(ovlPlace[oo].position.end);
      }

      bgnPos.finalize();
      endPos.finalize();

      op.position.bgn = bgnPos.mean();
      op.position.end = endPos.mean();

      op.bgnStdDev = bgnPos.stddev();   //  StdDev's are the same because the placement (placeFrag())
      op.endStdDev = endPos.stddev();   //  forces end to be set to bgn+readLen.

      //  Now that it is placed, estimate the span that is verified by overlaps.
      //  Threshold the floating end so it doesn't exceed the placement.
      //
      //  Annoyingly, the verified placement can, and does, exceed the bounds of the
      //  unitig, and we  need to check that threshold too.  Indel in the read and all that.
      
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

      //  Filter out bogus placements.  There used to be a few more, but they made no sense for long reads.
      //  Reject if either end stddev is high.  It has to be pretty bad before this triggers.

      bool   goodPlacement   = true;

      double allowableStdDev = max(2.0, 0.075 * FI->fragmentLength(op.frgID));

      if ((op.bgnStdDev > allowableStdDev) ||
          (op.endStdDev > allowableStdDev))
        goodPlacement = false;

      if ((flags & placeFrag_fullMatch) &&
          (op.fCoverage < 0.99))
        goodPlacement = false;

      if ((flags & placeFrag_noExtend) &&
          ((op.position.min() < 0) ||
           (op.position.max() > tigLen)))
        goodPlacement = false;


      if (goodPlacement)
        placements.push_back(op);


      if (logFileFlagSet(LOG_PLACE_FRAG))
        writeLog("pFUO()--   placements[%u] - PLACE FRAG %d in unitig %d at %d,%d (+- %.2f,%.2f) -- ovl %d,%d -- cov %d,%d %.2f -- errors %.2f aligned %d novl %d%s\n",
                 placements.size() - 1,
                 op.frgID, op.tigID,
                 op.position.bgn, op.position.end, op.bgnStdDev, op.endStdDev,
                 op.verified.bgn, op.verified.end,
                 op.covered.bgn, op.covered.end,
                 op.fCoverage,
                 op.errors, op.aligned, oe - os,
                 (goodPlacement == false) ? " -- INVALID"  : "");

      os = oe;
      oe = oe + 1;
    }  //  End of segregating overlaps by placement

    //  Move to the next block of overlaps.
    bgn = end;
    end = end + 1;
  }

  delete [] ovlPlace;

  //if (fid == 328)
  //  logFileFlags &= ~LOG_PLACE_FRAG;

  return(true);
}

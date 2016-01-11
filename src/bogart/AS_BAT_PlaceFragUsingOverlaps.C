
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "intervalList.H"

//  Report LOTS of details on placement, including evidence.
#undef  VERBOSE_PLACEMENT


//  Given an implicit fragment -- a ufNode with only the 'ident' set -- this will compute the
//  best placement for the fragment in an existing unitig.  ALL overlaps are used, not just
//  the best.
//
//  Ties are broken using overlap identities or arbitrarily.
//
//  Returns true if any placement is found, false otherwise.
//


bool
placeAcontainsB(Unitig *utg, ufNode &frag, BAToverlap &ovl, overlapPlacement &op) {
  BestContainment  best;

  //  The placeFrag() function is expecting the overlap to be from the container to the us fragment,
  //  which is opposite the overlap that we have.  We need to flip the fragments in the overlap --
  //  negate the hangs.

  best.container       = ovl.b_iid;  //  Not really the container...
  best.isContained     = false;      //  ...so mark this as a false BestContainment
  best.a_hang          = ovl.flipped ? ovl.b_hang : -ovl.a_hang;
  best.b_hang          = ovl.flipped ? ovl.a_hang : -ovl.b_hang;
  best.sameOrientation = ovl.flipped ? false : true;

  if (utg->placeFrag(frag, &best) == false)
    return(false);

  uint32   parentOrd = utg->pathPosition(ovl.b_iid);
  ufNode  &parent    = utg->ufpath[parentOrd];

  op.frgID       = frag.ident;
  op.refID       = ovl.b_iid;
  op.tigID       = utg->id();
  op.position    = frag.position;
  op.errors      = FI->fragmentLength(ovl.b_iid) * ovl.erate;
  op.covered.bgn = MIN(parent.position.bgn, parent.position.end);  //  Adjusted by hang later
  op.covered.end = MAX(parent.position.bgn, parent.position.end);
  op.aligned     = op.covered.end - op.covered.bgn;

  assert(op.covered.bgn < op.covered.end);

  //  Compute the portion of the unitig that is actually verified by
  //  the overlap.

  if (op.position.bgn < op.position.end) {
    int32  poslo = op.position.bgn;
    int32  poshi = op.position.end;

    assert(op.position.bgn <= op.covered.bgn);
    assert(op.position.bgn <= op.covered.end);

    op.covered.bgn -= op.position.bgn;
    op.covered.end -= op.position.bgn;

    op.verified.bgn = poslo + op.covered.bgn;
    op.verified.end = poslo + op.covered.end;

    if (op.verified.end > poshi)
      op.verified.end = poshi;

    assert(op.verified.bgn <  op.verified.end);
    assert(poslo           <= op.verified.bgn);
    assert(op.verified.end <= poshi);

  } else {
    int32  poslo = op.position.end;
    int32  poshi = op.position.bgn;

    assert(op.position.end <= op.covered.bgn);
    assert(op.position.end <= op.covered.end);

    op.covered.bgn -= op.position.end;
    op.covered.end -= op.position.end;

    op.verified.bgn = poslo + op.covered.end;
    op.verified.end = poslo + op.covered.bgn;

    if (op.verified.bgn > poshi)
      op.verified.bgn = poshi;

    assert(op.verified.end <  op.verified.bgn);
    assert(poslo           <= op.verified.end);
    assert(op.verified.bgn <= poshi);
  }

  //  Disallow any placements that exceed the boundary of the unitig.  These cannot be confirmed
  //  by overlaps and might be wrong.  Sample cases:
  //    o  sticking a unique/repeat fragment onto a repeat (leaving the unique uncovered)
  //    o  sticking a chimeric fragment onto the end of a unitig (leaving the chimeric join uncovered)

  if ((MIN(op.position.bgn, op.position.end) < 0) ||
      (MAX(op.position.bgn, op.position.end) > utg->getLength())) {
#ifdef VERBOSE_PLACEMENT
    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("placeFragUsingOverlaps()-- (container) - frag %d in unitig %d at %d,%d (verified %d,%d) from overlap ident %d %d hang %d %d flipped %d covered %d,%d DISALLOWED\n",
              frag.ident, utg->id(), op.position.bgn, op.position.end, op.verified.bgn, op.verified.end,
              ovl.a_iid, ovl.b_iid, ovl.a_hang, ovl.b_hang, ovl.flipped,
              op.covered.bgn, op.covered.end);
#endif
    op = overlapPlacement();

  } else {
#ifdef VERBOSE_PLACEMENT
    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("placeFragUsingOverlaps()-- (container) - frag %d in unitig %d at %d,%d (verified %d,%d) from overlap ident %d %d hang %d %d flipped %d covered %d,%d\n",
              frag.ident, utg->id(), op.position.bgn, op.position.end, op.verified.bgn, op.verified.end,
              ovl.a_iid, ovl.b_iid, ovl.a_hang, ovl.b_hang, ovl.flipped,
              op.covered.bgn, op.covered.end);
#endif
  }

  return(true);
}



bool
placeBcontainsA(Unitig *utg, ufNode &frag, BAToverlap &ovl, overlapPlacement &op) {
  BestContainment  best;

  best.container       = ovl.b_iid;
  best.isContained     = true;
  best.a_hang          = ovl.flipped ? ovl.b_hang : -ovl.a_hang;
  best.b_hang          = ovl.flipped ? ovl.a_hang : -ovl.b_hang;
  best.sameOrientation = ovl.flipped ? false : true;

  if (utg->placeFrag(frag, &best) == false)
    return(false);

  op.frgID       = frag.ident;
  op.refID       = ovl.b_iid;
  op.tigID       = utg->id();
  op.position    = frag.position;
  op.errors      = FI->fragmentLength(ovl.a_iid) * ovl.erate;
  op.covered.bgn = 0;
  op.covered.end = FI->fragmentLength(ovl.a_iid);
  op.aligned     = op.covered.end - op.covered.bgn;

  assert(op.covered.bgn < op.covered.end);

  op.verified.bgn = op.position.bgn;
  op.verified.end = op.position.end;

  if ((MIN(op.position.bgn, op.position.end) < 0) ||
      (MAX(op.position.bgn, op.position.end) > utg->getLength())) {
#ifdef VERBOSE_PLACEMENT
    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("placeFragUsingOverlaps()-- (contained) - frag %d in unitig %d at %d,%d (verified %d,%d) from overlap ident %d %d hang %d %d flipped %d covered %d,%d DISALLOWED\n",
              frag.ident, utg->id(), op.position.bgn, op.position.end, op.verified.bgn, op.verified.end,
              ovl.a_iid, ovl.b_iid, ovl.a_hang, ovl.b_hang, ovl.flipped,
              op.covered.bgn, op.covered.end);
#endif
    op = overlapPlacement();

  } else {
#ifdef VERBOSE_PLACEMENT
    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("placeFragUsingOverlaps()-- (contained) - frag %d in unitig %d at %d,%d (verified %d,%d) from overlap ident %d %d hang %d %d flipped %d covered %d,%d\n",
              frag.ident, utg->id(), op.position.bgn, op.position.end, op.verified.bgn, op.verified.end,
              ovl.a_iid, ovl.b_iid, ovl.a_hang, ovl.b_hang, ovl.flipped,
              op.covered.bgn, op.covered.end);
#endif
  }

  return(true);
}



bool
placeDovetail(Unitig *utg, ufNode &frag, BAToverlap &ovl, overlapPlacement &op) {
  BestEdgeOverlap   best;
  int32             plac3, plac5;
  int32             aend3p = AS_BAT_overlapAEndIs3prime(ovl);

  best.set(ovl);

  if (utg->placeFrag(frag, plac5, (aend3p ? NULL  : &best),
                     frag, plac3, (aend3p ? &best : NULL)) == false)
    return(false);

  uint32  olen = FI->overlapLength(ovl.a_iid, ovl.b_iid, ovl.a_hang, ovl.b_hang);
  uint32  flen = FI->fragmentLength(ovl.a_iid);

  op.frgID       = frag.ident;
  op.refID       = ovl.b_iid;
  op.tigID       = utg->id();
  op.position    = frag.position;
  op.errors      = olen * ovl.erate;
  op.covered.bgn = (ovl.a_hang < 0) ?    0 : ovl.a_hang;
  op.covered.end = (ovl.b_hang > 0) ? flen : ovl.b_hang + flen;
  op.aligned     = op.covered.end - op.covered.bgn;

  assert(op.covered.bgn < op.covered.end);

  if (op.position.bgn < op.position.end) {
    op.verified.bgn = op.position.bgn + op.covered.bgn;
    op.verified.end = op.position.bgn + op.covered.end;

    if (op.verified.end > op.position.end)
      op.verified.end = op.position.end;

    assert(op.verified.bgn >= op.position.bgn);
    assert(op.verified.end <= op.position.end);
  } else {
    op.verified.bgn = op.position.bgn - op.covered.bgn;
    op.verified.end = op.position.bgn - op.covered.end;

    if (op.verified.end < op.position.end)
      op.verified.end = op.position.end;

    assert(op.verified.end >= op.position.end);
    assert(op.verified.bgn <= op.position.bgn);
  }

  if ((MIN(op.position.bgn, op.position.end) < 0) ||
      (MAX(op.position.bgn, op.position.end) > utg->getLength())) {
#ifdef VERBOSE_PLACEMENT
    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("placeFragUsingOverlaps()-- (dovetail)  - frag %d in unitig %d at %d,%d (verified %d,%d) from overlap ident %d %d hang %d %d flipped %d covered %d,%d DISALLOWED\n",
              frag.ident, utg->id(), op.position.bgn, op.position.end, op.verified.bgn, op.verified.end,
              ovl.a_iid, ovl.b_iid, ovl.a_hang, ovl.b_hang, ovl.flipped,
              op.covered.bgn, op.covered.end);
#endif
    op = overlapPlacement();

  } else {
#ifdef VERBOSE_PLACEMENT
    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("placeFragUsingOverlaps()-- (dovetail)  - frag %d in unitig %d at %d,%d (verified %d,%d) from overlap ident %d %d hang %d %d flipped %d covered %d,%d\n",
              frag.ident, utg->id(), op.position.bgn, op.position.end, op.verified.bgn, op.verified.end,
              ovl.a_iid, ovl.b_iid, ovl.a_hang, ovl.b_hang, ovl.flipped,
              op.covered.bgn, op.covered.end);
#endif
  }

  return(true);
}







bool
placeFragUsingOverlaps(UnitigVector             &unitigs,
                       double                    erate,
                       Unitig                   *target,
                       uint32                    fid,
                       vector<overlapPlacement> &placements) {

  //logFileFlags |= LOG_PLACE_FRAG;

  if (logFileFlagSet(LOG_PLACE_FRAG))
    writeLog("placeFragUsingOverlaps()-- begin for frag %d into target tig %d\n", fid, target->id());

  assert(fid > 0);
  assert(fid <= FI->numFragments());

  ufNode frag;

  frag.ident             = fid;
  frag.contained         = 0;
  frag.parent            = 0;
  frag.ahang             = 0;
  frag.bhang             = 0;
  frag.position.bgn      = 0;
  frag.position.end      = 0;
  frag.containment_depth = 0;

  placements.clear();

  uint32      ovlLen = 0;
  BAToverlap *ovl    = OC->getOverlaps(frag.ident, erate, ovlLen);

  overlapPlacement   *ovlPlace = new overlapPlacement[ovlLen];
  uint32              nFragmentsNotPlaced = 0;

  //  Initialize placements to nowhere.

  for (uint32 i=0; i<ovlLen; i++)
    ovlPlace[i] = overlapPlacement();

  //  Compute placements.  Anything that doesn't get placed is left as 'nowhere', in particular, unitig == 0.

  for (uint32 i=0; i<ovlLen; i++) {
    int32             utgID = Unitig::fragIn(ovl[i].b_iid);
    Unitig           *utg   = unitigs[utgID];

    assert(ovl[i].a_iid == frag.ident);

    if (utgID == 0)
      //  Fragment not in a unitig yet -- possibly this is a contained fragment that we haven't
      //  placed yet, or have temporarily removed it from a unitig.
      continue;

    if ((target != NULL) && (target != utg))
      //  Requested placement in a specific unitig, and this isn't it.
      continue;

    //  Depending on the type of overlap (containment vs dovetail), place the fragment relative to
    //  the other fragment.

    if        ((ovl[i].a_hang >= 0) && (ovl[i].b_hang <= 0)) {
      //  A (us) contains B (the other fragment)
      if (placeAcontainsB(utg, frag, ovl[i], ovlPlace[i]) == false)
        nFragmentsNotPlaced++;

    } else if ((ovl[i].a_hang <= 0) && (ovl[i].b_hang >= 0)) {
      //  A (us) is contained in B (the other fragment)
      if (placeBcontainsA(utg, frag, ovl[i], ovlPlace[i]) == false)
        nFragmentsNotPlaced++;

    } else {
      //  A dovetail, use the existing placement routine
      if (placeDovetail(utg, frag, ovl[i], ovlPlace[i]) == false)
        nFragmentsNotPlaced++;
    }

    assert((ovlPlace[i].position.bgn < ovlPlace[i].position.end) == (ovlPlace[i].verified.bgn < ovlPlace[i].verified.end));
  }  //  Over all overlaps.


  //  Report if any of the placement routines fail.  This shouldn't happen, but if it does, it is
  //  hardly fatal.

#ifdef VERBOSE_PLACEMENT
  if (nFragmentsNotPlaced > 0)
    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("placeFragUsingOverlaps()-- WARNING: Failed to place %d fragments\n", nFragmentsNotPlaced);
#endif

  //  Sort all the placements.  Any overlap we couldn't place is automatically in Unitig 0, the
  //  invalid unitig.  Sort order is by unitig ID, then by orientation, then by position.
  //
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
    //  Orientation of 'position' and 'verified' is the same, asserted above.

    end = bgn + 1;
    while ((end < ovlLen) &&
           (ovlPlace[bgn].tigID == ovlPlace[end].tigID) &&
           (isReverse(ovlPlace[bgn].verified) == isReverse(ovlPlace[end].verified)))
      end++;

    //  Over all placements with the same unitig/orientation (that'd be from bgn to end), build
    //  interval lists for the begin point and the end point.  Remember, this is all fragments to a
    //  single unitig (the whole picture above), not just the overlapping fragment sets (left or
    //  right blocks).

    intervalList<int32>   bgnPoints;
    intervalList<int32>   endPoints;

    int32                 windowSlop = 0.075 * FI->fragmentLength(frag.ident);

    if (windowSlop < 5)
      windowSlop = 5;

    for (uint32 oo=bgn; oo<end; oo++) {
      assert(ovlPlace[oo].tigID > 0);

      int32   b  = ovlPlace[oo].verified.bgn;
      int32   be = ovlPlace[oo].verified.bgn + windowSlop;
      int32   e  = ovlPlace[oo].verified.end;
      int32   ee = ovlPlace[oo].verified.end + windowSlop;

      b = (b < windowSlop) ? 0 : b - windowSlop;
      e = (e < windowSlop) ? 0 : e - windowSlop;

      bgnPoints.add(b, be - b);
      endPoints.add(e, ee - e);
    }

    bgnPoints.merge();
    endPoints.merge();

    //  Now, assign each placement to a end-pair cluster based on the interval ID that the end point
    //  falls in.
    //
    //  count the number of fragments that hit each pair of points.  We can do this two ways:
    //    1)  With a list of point-pairs that we sort and count           -- O(n) size, O(n log n) time
    //    2)  With an array of all point-pairs that we increment directly -- O(p*p) size, O(n) time
    //  Typically, p is small.

    int32   numBgnPoints = bgnPoints.numberOfIntervals();
    int32   numEndPoints = endPoints.numberOfIntervals();

    for (uint32 oo=bgn; oo<end; oo++) {
      int32   b = ovlPlace[oo].verified.bgn;  //  WAS expected position of read in tig!
      int32   e = ovlPlace[oo].verified.end;
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

    sort(ovlPlace + bgn, ovlPlace + end, overlapPlacement_byCluster);

    //  Run through each 'cluster' and compute the placement.

    for (uint32 os=bgn, oe=bgn; os<end; ) {
      overlapPlacement  op;

      while ((oe < end) && (ovlPlace[os].clusterID == ovlPlace[oe].clusterID))
        oe++;

      //  Overlaps from os to oe are all for a single location.  Examine them to fill out an
      //  overlapPlacement, including scores.
      //
      //  position:   the MAX extent (which is actually exactly what the intervalList computed).  A possibly
      //              better solution is to use the mode.
      //
      //  errors:     sum of the estimated number of errors in all the overlaps
      //
      //  fCoverage:  coverage of the fragment.  Instead of building another interval list, this is approximated
      //              by (max-min) overlap position.

      op.frgID = frag.ident;
      op.refID = ovlPlace[os].refID;
      op.tigID = ovlPlace[os].tigID;

      op.fCoverage   = 0.0;

      op.errors      = 0.0;
      op.aligned     = 0;

      assert((ovlPlace[os].position.bgn < ovlPlace[os].position.end) == (ovlPlace[os].verified.bgn < ovlPlace[os].verified.end));

      //  op.position is not set yet.
      //assert((op.position.bgn < op.position.end)                     == (ovlPlace[os].verified.bgn < ovlPlace[os].verified.end));

      op.verified.bgn = ovlPlace[os].verified.bgn;
      op.verified.end = ovlPlace[os].verified.end;

      op.covered.bgn = ovlPlace[os].covered.bgn;
      op.covered.end = ovlPlace[os].covered.end;

      uint32  nForward = 0;
      uint32  nReverse = 0;

      for (uint32 oo=os; oo<oe; oo++) {
        assert(op.tigID == ovlPlace[oo].tigID);

        op.errors      += ovlPlace[oo].errors;
        op.aligned     += ovlPlace[oo].aligned;

        op.covered.bgn  = MIN(op.covered.bgn, ovlPlace[oo].covered.bgn);
        op.covered.end  = MAX(op.covered.end, ovlPlace[oo].covered.end);

        if (isReverse(ovlPlace[oo].position))
          nReverse++;
        else
          nForward++;
      }

      assert((nReverse == 0) || (nForward == 0));

      op.fCoverage = (op.covered.end - op.covered.bgn) / (double)FI->fragmentLength(op.frgID);

      //  Find the first and last fragment in the unitig that we overlap with.
      //
      //  The first fragment is easy.  We can run through the list of overlaps, ask the unitig for
      //  the ordinal of each fragment, and return the lowest.
      //
      //  The last fragment is similar, but we need to return the highest ordinal that is also the
      //  longest.  (In the first fragment case, we are guaranteed by the construction of the unitig
      //  to have the earliest fragment position).

      Unitig *destTig = unitigs[op.tigID];

      uint32  firstOrdinal  = UINT32_MAX;
      uint32  firstPosition = UINT32_MAX;
      uint32  lastOrdinal   = 0;
      uint32  lastPosition  = 0;

      FragmentEnd  firstEnd;
      FragmentEnd  lastEnd;

      for (uint32 oo=os; oo<oe; oo++) {
        uint32   ordinal = destTig->pathPosition(ovlPlace[oo].refID);
        ufNode  &ovlFrg  = destTig->ufpath[ordinal];
        uint32   minPos  = MIN(ovlFrg.position.bgn, ovlFrg.position.end);
        uint32   maxPos  = MAX(ovlFrg.position.bgn, ovlFrg.position.end);

        //writeLog("placeFragUsingOverlaps()-- PickEnds ordinal %d tigFrg %d pos %d,%d\n",
        //        ordinal, ovlFrg.ident, minPos, maxPos);

        //  For a normal dovetail alignment, this is pretty straight forward.  We pick the
        //  end we align to.
        //
        //  For spur alignments (repeat detection) it is backwards to what we want.  More comments
        //  in repeatJunctionEvidence::repeatJunctionEvidence().
        //
        //         \                        /
        //          -----alignedfragment----
        //          ------....    .....-----
        //

        if (((minPos  <  firstPosition)) ||
            ((minPos  <= firstPosition) && (ordinal < firstOrdinal))) {
          firstOrdinal  = ordinal;
          firstPosition = minPos;
          firstEnd      = FragmentEnd(ovlFrg.ident, (ovlFrg.position.bgn < ovlFrg.position.end));
        }

        if (((maxPos  >  lastPosition)) ||
            ((maxPos  >= lastPosition) && (ordinal > lastOrdinal))) {
          lastOrdinal  = ordinal;
          lastPosition = minPos;
          lastEnd      = FragmentEnd(ovlFrg.ident, (ovlFrg.position.end < ovlFrg.position.bgn));
        }
      }

      if (nForward > 0) {
        op.frag5p = firstEnd;
        op.frag3p = lastEnd;
      } else {
        op.frag5p = lastEnd;
        op.frag3p = firstEnd;
      }

      //  Compute mean and stddev placement.

      uint32  numPlace = 0;
      double  bgnMean = 0;
      double  endMean = 0;

      op.bgnStdDev = 0.0;
      op.endStdDev = 0.0;

      for (uint32 oo=os; oo<oe; oo++) {
        if ((ovlPlace[oo].position.bgn == 0) &&
            (ovlPlace[oo].position.end == 0))
          continue;

        if (ovlPlace[oo].position.bgn < ovlPlace[oo].position.end) {
          if (ovlPlace[oo].verified.bgn >= ovlPlace[oo].verified.end)
            writeLog("placeFragUsingOverlaps()-- frag %d FWD verified placement invalid (bgn,end %d,%d) for position (bgn,end %d,%d)\n",
                    ovlPlace[oo].frgID,
                    ovlPlace[oo].verified.bgn, ovlPlace[oo].verified.end,
                    ovlPlace[oo].position.bgn, ovlPlace[oo].position.end);
          assert(ovlPlace[oo].verified.bgn < ovlPlace[oo].verified.end);

          bgnMean += ovlPlace[oo].position.bgn;
          endMean += ovlPlace[oo].position.end;

          op.verified.bgn  = MIN(op.verified.bgn, ovlPlace[oo].verified.bgn);
          op.verified.end  = MAX(op.verified.end, ovlPlace[oo].verified.end);

        } else {
          if (ovlPlace[oo].verified.bgn < ovlPlace[oo].verified.end)
            writeLog("placeFragUsingOverlaps()-- frag %d REV verified placement invalid (bgn,end %d,%d) for position (bgn,end %d,%d)\n",
                    ovlPlace[oo].frgID,
                    ovlPlace[oo].verified.bgn, ovlPlace[oo].verified.end,
                    ovlPlace[oo].position.bgn, ovlPlace[oo].position.end);
          assert(ovlPlace[oo].verified.bgn >= ovlPlace[oo].verified.end);

          bgnMean += ovlPlace[oo].position.end;
          endMean += ovlPlace[oo].position.bgn;

          op.verified.bgn  = MAX(op.verified.bgn, ovlPlace[oo].verified.bgn);
          op.verified.end  = MIN(op.verified.end, ovlPlace[oo].verified.end);
        }

        numPlace++;
      }

      bgnMean /= numPlace;
      endMean /= numPlace;

      op.position.bgn = (int32)((nReverse == 0) ? bgnMean : endMean);
      op.position.end = (int32)((nReverse == 0) ? endMean : bgnMean);

      for (uint32 oo=os; oo<oe; oo++) {
        if ((ovlPlace[oo].position.bgn == 0) &&
            (ovlPlace[oo].position.end == 0))
          continue;

        if (ovlPlace[oo].position.bgn < ovlPlace[oo].position.end) {
          op.bgnStdDev += (ovlPlace[oo].position.bgn - bgnMean) * (ovlPlace[oo].position.bgn - bgnMean);
          op.endStdDev += (ovlPlace[oo].position.end - endMean) * (ovlPlace[oo].position.end - endMean);
        } else {
          op.bgnStdDev += (ovlPlace[oo].position.end - bgnMean) * (ovlPlace[oo].position.end - bgnMean);
          op.endStdDev += (ovlPlace[oo].position.bgn - endMean) * (ovlPlace[oo].position.bgn - endMean);
        }
      }

      op.bgnStdDev = sqrt(op.bgnStdDev / numPlace);
      op.endStdDev = sqrt(op.endStdDev / numPlace);

      //  Filter out bogus placements.
      //
      //  This placement is invalid if the std.dev is too high on either end.  (Was 3% fragment length before 11 Apr 2013)
      //  This placement is invalid if both nReverse and nForward are more than zero.

      bool   weakStdDev      = false;
      bool   overlappingSpan = false;
      bool   spanBad         = false;
      bool   reject          = false;

      double allowableStdDev = MAX(2.0, 0.075 * FI->fragmentLength(op.frgID));

      if ((op.bgnStdDev > allowableStdDev) ||
          (op.endStdDev > allowableStdDev))
        weakStdDev = true;

      if (((op.position.bgn < op.position.end) && (op.position.bgn + 3 * op.bgnStdDev > op.position.end - 3 * op.endStdDev)) ||
          ((op.position.end < op.position.bgn) && (op.position.end + 3 * op.endStdDev > op.position.bgn - 3 * op.bgnStdDev)))
        overlappingSpan = true;

      int32   poslen = (op.position.end > op.position.bgn) ? (op.position.end - op.position.bgn) : (op.position.bgn - op.position.end);
      int32   trulen = FI->fragmentLength(op.frgID);
      double  scaled = (double)poslen / trulen;

      if ((scaled < 0.3333) ||
          (2.0    < scaled))
        spanBad = true;


      if ((weakStdDev) && (0))
        //  Read is not known to have lots of indel, but the stddev is high.
        reject = true;

      if (overlappingSpan)
        //  Read placements are conflicting and overlapping.
        reject = true;

      if (spanBad)
        //  Bogus placement, more than twice as large as expected, or less than 1/3 expected.
        reject = true;


      if (reject) {
#ifdef VERBOSE_PLACEMENT
        if (logFileFlagSet(LOG_PLACE_FRAG)) {
          writeLog("placeFragUsingOverlaps()-- frag %d in unitig %d at %d,%d (+- %.2f,%.2f) -- cov %.2f (%d,%d) errors %.2f aligned %d novl %d -- INVALID stddev weak %d overlapping %d bad size %d\n",
                   op.frgID, op.tigID, op.position.bgn, op.position.end, op.bgnStdDev, op.endStdDev,
                   op.fCoverage, op.covered.bgn, op.covered.end,
                   op.errors,
                   op.aligned,
                   oe - os,
                   weakStdDev, overlappingSpan, spanBad);
          for (uint32 oo=os; oo<oe; oo++) {
            if ((ovlPlace[oo].position.bgn == 0) &&
                (ovlPlace[oo].position.end == 0))
              continue;

            writeLog("placeFragUsingOverlaps()--   %8u,%8u\n", ovlPlace[oo].position.bgn, ovlPlace[oo].position.end);

          }
        }
#endif

      } else {
        placements.push_back(op);
#ifdef VERBOSE_PLACEMENT
        if (logFileFlagSet(LOG_PLACE_FRAG))
          writeLog("placeFragUsingOverlaps()-- frag %d in unitig %d at %d,%d (+- %.2f,%.2f) -- cov %.2f (%d,%d) errors %.2f aligned %d novl %d\n",
                  op.frgID, op.tigID, op.position.bgn, op.position.end, op.bgnStdDev, op.endStdDev,
                  op.fCoverage, op.covered.bgn, op.covered.end,
                  op.errors,
                  op.aligned,
                  oe - os);
#endif
      }

      os = oe;
      oe = oe + 1;
    }  //  End of segregating overlaps by placement

    //  Move to the next block of overlaps.
    bgn = end;
    end = end + 1;
  }

  delete [] ovlPlace;

  //logFileFlags &= ~LOG_PLACE_FRAG;

  return(true);
}




void
placeFragInBestLocation(UnitigVector   &unitigs,
                        double          erate,
                        uint32          fid) {

  if (Unitig::fragIn(fid) != 0)
    //  Already placed.
    return;

  ufNode                    frg;
  vector<overlapPlacement>  op;

  frg.ident             = fid;
  frg.contained         = 0;
  frg.parent            = 0;
  frg.ahang             = 0;
  frg.bhang             = 0;
  frg.position.bgn      = 0;
  frg.position.end      = FI->fragmentLength(fid);
  frg.containment_depth = 0;

  placeFragUsingOverlaps(unitigs, erate, NULL, fid, op);

  //  Pick the lowest error placement, and of those lowest, the least aligned region.

  double  minError = DBL_MAX;
  uint32  minAlign = UINT32_MAX;
  uint32  bp       = UINT32_MAX;

  for (uint32 pl=0; pl<op.size(); pl++) {
    if (op[pl].fCoverage < 0.99)
      continue;

    double e = op[pl].errors / op[pl].aligned;

    if ((e < minError) ||
        ((e <= minError) && (op[pl].aligned < minAlign))) {
      minError = e;
      minAlign = op[pl].aligned;
      bp       = pl;
    }
  }

  //  No placement?  New unitig!

  if (bp == UINT32_MAX) {
    Unitig  *sing = unitigs.newUnitig(false);
    sing->addFrag(frg, 0, false);
    return;
  }

  //  Place the frag in the unitig at the spot.

  Unitig  *tig = unitigs[op[bp].tigID];

  frg.position.bgn = op[bp].position.bgn;
  frg.position.end = op[bp].position.end;

  tig->addFrag(frg, 0, false);
  tig->bubbleSortLastFrag();
}

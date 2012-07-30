
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static const char *rcsid = "$Id: AS_BAT_PlaceFragUsingOverlaps.C,v 1.17 2012-07-30 01:21:01 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "AS_UTL_intervalList.H"

#define VERBOSE

//  Given an implicit fragment -- a ufNode with only the 'ident' set -- this will compute the
//  best placement for the fragment in an existing unitig.  ALL overlaps are used, not just
//  the best.
//
//  Ties are broken using mate pairs, overlap identities, or arbitrarily.
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

  op.frgID       = frag.ident;
  op.refID       = ovl.b_iid;
  op.tigID       = utg->id();
  op.position    = frag.position;
  op.errors      = FI->fragmentLength(ovl.b_iid) * ovl.error;
  op.covered.bgn = ovl.a_hang;
  op.covered.end = FI->fragmentLength(ovl.a_iid) + ovl.b_hang;
  op.aligned     = op.covered.end - op.covered.bgn;

  assert(op.covered.bgn < op.covered.end);

  //  Compute the portion of the unitig that is actually verified by
  //  the overlap.

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

    if (op.position.bgn < op.covered.bgn)
      op.verified.bgn = 0;

    if (op.verified.end < op.position.end)
      op.verified.end = op.position.end;

    assert(op.verified.end >= op.position.end);
    assert(op.verified.bgn <= op.position.bgn);
  }

  //  Disallow any placements that exceed the boundary of the unitig.  These cannot be confirmed
  //  by overlaps and might be wrong.  Sample cases:
  //    o  sticking a unique/repeat fragment onto a repeat (leaving the unique uncovered)
  //    o  sticking a chimeric fragment onto the end of a unitig (leaving the chimeric join uncovered)

  if ((MIN(op.position.bgn, op.position.end) < 0) ||
      (MAX(op.position.bgn, op.position.end) > utg->getLength())) {
#ifdef VERBOSE
    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("placeFragUsingOverlaps()-- (container) - frag %d in unitig %d at %d,%d (verified %d,%d) from overlap ident %d %d hang %d %d flipped %d covered %d,%d DISALLOWED\n",
              frag.ident, utg->id(), op.position.bgn, op.position.end, op.verified.bgn, op.verified.end,
              ovl.a_iid, ovl.b_iid, ovl.a_hang, ovl.b_hang, ovl.flipped,
              op.covered.bgn, op.covered.end);
#endif
    op = overlapPlacement();

  } else {
#ifdef VERBOSE
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
  op.errors      = FI->fragmentLength(ovl.a_iid) * ovl.error;
  op.covered.bgn = 0;
  op.covered.end = FI->fragmentLength(ovl.a_iid);
  op.aligned     = op.covered.end - op.covered.bgn;

  assert(op.covered.bgn < op.covered.end);

  op.verified.bgn = op.position.bgn;
  op.verified.end = op.position.end;

  if ((MIN(op.position.bgn, op.position.end) < 0) ||
      (MAX(op.position.bgn, op.position.end) > utg->getLength())) {
#ifdef VERBOSE
    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("placeFragUsingOverlaps()-- (contained) - frag %d in unitig %d at %d,%d (verified %d,%d) from overlap ident %d %d hang %d %d flipped %d covered %d,%d DISALLOWED\n",
              frag.ident, utg->id(), op.position.bgn, op.position.end, op.verified.bgn, op.verified.end,
              ovl.a_iid, ovl.b_iid, ovl.a_hang, ovl.b_hang, ovl.flipped,
              op.covered.bgn, op.covered.end);
#endif
    op = overlapPlacement();

  } else {
#ifdef VERBOSE
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
  op.errors      = olen * ovl.error;
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
#ifdef VERBOSE
    if (logFileFlagSet(LOG_PLACE_FRAG))
      writeLog("placeFragUsingOverlaps()-- (dovetail)  - frag %d in unitig %d at %d,%d (verified %d,%d) from overlap ident %d %d hang %d %d flipped %d covered %d,%d DISALLOWED\n",
              frag.ident, utg->id(), op.position.bgn, op.position.end, op.verified.bgn, op.verified.end,
              ovl.a_iid, ovl.b_iid, ovl.a_hang, ovl.b_hang, ovl.flipped,
              op.covered.bgn, op.covered.end);
#endif
    op = overlapPlacement();

  } else {
#ifdef VERBOSE
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
                       Unitig                   *target,
                       AS_IID                    fid,
                       vector<overlapPlacement> &placements) {

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
  BAToverlap *ovl    = OC->getOverlaps(frag.ident, ovlLen);

  overlapPlacement   *ovlPlace = new overlapPlacement[ovlLen];
  uint32              nFragmentsNotPlaced = 0;

  for (uint32 i=0; i<ovlLen; i++)
    ovlPlace[i] = overlapPlacement();

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
  }  //  Over all overlaps.


  //  For whatever reason, the fragment placement routines failed to place a fragment using an overlap.
  //  This shouldn't happen, but if it does, it is hardly fatal.
#if 0
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
  //  This unitig has two sets of "overlapping overlaps".  The left set is probably caused caused
  //  by a tandem repeat.  We'll get two overlaps to the 2,4 pair, and one overlap to the 1,3
  //  pair.  Assuming that is good enough to be a clear winner, we'll ignore the 1,3 overlap and compute
  //  position based on the other two overlaps.



  uint32         bgn = 0;  //  Range of overlaps with the same unitig/orientation
  uint32         end = 1;

  while ((bgn < ovlLen) && (ovlPlace[bgn].tigID == 0))
    //  Skip overlaps that didn't generate a placement
    bgn++;


  //  Process all placements.
  while (bgn < ovlLen) {

    //  Find the last placement with the same unitig/orientation as the 'bgn' fragment.
    end = bgn + 1;
    while ((end < ovlLen) &&
           (ovlPlace[bgn].tigID == ovlPlace[end].tigID) &&
           (isReverse(ovlPlace[bgn].position) == isReverse(ovlPlace[end].position)))
      end++;

#if 0
    writeLog("PROCESS bgn=%d to end=%d\n", bgn, end);
#endif

    //  Over all placements with the same unitig/orientation, build interval lists for the begin
    //  point and the end point.  Remember, this is all fragments to a single unitig (the whole
    //  picture above), not just the overlapping fragment sets (left or right in the above picture
    //  above).
    //
    intervalList   bgnPoints;
    intervalList   endPoints;
   
    for (uint32 oo=bgn; oo<end; oo++) {
      int32   b = ovlPlace[oo].position.bgn;
      int32   e = ovlPlace[oo].position.end;

      assert(ovlPlace[oo].tigID > 0);

      bgnPoints.add(b-5, 10);
      endPoints.add(e-5, 10);
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
    
    //int32   *pointPairCount = new int32 [bgnPoints.numberOfIntervals() * endPoints.numberOfIntervals()];
    //memset(pointPairCount, 0, sizeof(int32) * bgnPoints.numberOfIntervals() * endPoints.numberOfIntervals());

    int32   numBgnPoints = bgnPoints.numberOfIntervals();
    int32   numEndPoints = endPoints.numberOfIntervals();

    for (uint32 oo=bgn; oo<end; oo++) {
      int32   b = ovlPlace[oo].position.bgn;
      int32   e = ovlPlace[oo].position.end;

      int32   c = 0;

      ovlPlace[oo].clusterID = 0;

      for (int32 r=0; r<numBgnPoints; r++)
        if ((bgnPoints.lo(r) <= b) && (b <= bgnPoints.hi(r))) {
          assert(ovlPlace[oo].clusterID == 0);
          ovlPlace[oo].clusterID = c = r * numEndPoints + 1;
        }

      for (int32 r=0; r<numEndPoints; r++)
        if ((endPoints.lo(r) <= e) && (e <= endPoints.hi(r))) {
          assert(ovlPlace[oo].clusterID == c);
          ovlPlace[oo].clusterID += r;
        }
    }

      
    sort(ovlPlace + bgn, ovlPlace + end, overlapPlacement_byCluster);

#if 0
    writeLog("POSTSORT\n");
    for (uint32 oo=bgn; oo<end; oo++)
      writeLog("overlapPlacement [%d] cluster %d frg %d tig %d position %d,%d errors %.2f covered %d,%d aligned %d\n",
              oo,
              ovlPlace[oo].clusterID,
              ovlPlace[oo].frgID,
              ovlPlace[oo].tigID,
              ovlPlace[oo].position.bgn, ovlPlace[oo].position.end,
              ovlPlace[oo].errors,
              ovlPlace[oo].covered.bgn,  ovlPlace[oo].covered.end,
              ovlPlace[oo].aligned);
    fflush(logFile);
#endif

    //  Run through each 'cluster' and compute the placement.

    for (uint32 os=bgn, oe=bgn; os<end; ) {
      overlapPlacement  op;

      while ((oe < end) && (ovlPlace[os].clusterID == ovlPlace[oe].clusterID))
        oe++;

#if 0
      writeLog("PROCESS os=%d to oe=%d in range bgn=%d end=%d\n", os, oe, bgn, end);
#endif

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
          assert(ovlPlace[oo].verified.bgn < ovlPlace[oo].verified.end);

          bgnMean += ovlPlace[oo].position.bgn;
          endMean += ovlPlace[oo].position.end;

          op.verified.bgn  = MIN(op.verified.bgn, ovlPlace[oo].verified.bgn);
          op.verified.end  = MAX(op.verified.end, ovlPlace[oo].verified.end);

        } else {
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
      //  This placement is invalid if the std.dev is too high on either end.
      //  This placement is invalid if both nReverse and nForward are more than zero.

      double allowableStdDev = MAX(2.0, 0.03 * FI->fragmentLength(op.frgID));

      if        ((op.bgnStdDev > allowableStdDev) ||
                 (op.endStdDev > allowableStdDev)) {
#if 0
        if (logFileFlagSet(LOG_PLACE_FRAG))
          writeLog("placeFragUsingOverlaps()-- frag %d in unitig %d at %d,%d (+- %.2f,%.2f) -- cov %.2f (%d,%d) errors %.2f aligned %d novl %d -- INVALID standard deviation too large (min %0.2f)\n",
                  op.frgID, op.tigID, op.position.bgn, op.position.end, op.bgnStdDev, op.endStdDev,
                  op.fCoverage, op.covered.bgn, op.covered.end,
                  op.errors,
                  op.aligned,
                  oe - os,
                  allowableStdDev);
#endif

      } else {
        placements.push_back(op);
#if 0
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

  return(true);
}




void
placeUnmatedFragInBestLocation(UnitigVector   &unitigs,
                               AS_IID          fid) {
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

  placeFragUsingOverlaps(unitigs, NULL, fid, op);

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


void
placeMatedFragInBestLocation(UnitigVector   &unitigs,
                             AS_IID          fid,
                             AS_IID          mid) {
  placeUnmatedFragInBestLocation(unitigs, fid);
}


void
placeFragInBestLocation(UnitigVector   &unitigs,
                        AS_IID          fid) {

  if (Unitig::fragIn(fid) != 0)
    //  Already placed.
    return;

  AS_IID  mid = FI->mateIID(fid);

  if (mid == 0)
    placeUnmatedFragInBestLocation(unitigs, fid);
  else
    placeMatedFragInBestLocation(unitigs, fid, mid);
}

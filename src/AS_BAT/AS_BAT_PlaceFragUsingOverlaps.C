
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

static const char *rcsid = "$Id: AS_BAT_PlaceFragUsingOverlaps.C,v 1.3 2010-12-06 08:03:48 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "AS_UTL_intervalList.H"

#define MAX_OVERLAPS_PER_FRAG   (16 * 1024 * 1024)


//  Given an implicit fragment -- a ufNode with only the 'ident' set -- this will compute the
//  best placement for the fragment in an existing unitig.  ALL overlaps are used, not just
//  the best.
//
//  Ties are broken using mate pairs, overlap identities, or arbitrarily.
//
//  Returns true if any placement is found, false otherwise.
//
bool
placeFragUsingOverlaps(UnitigVector &unitigs,
                       uint32        fid,
                       OverlapStore *ovlStoreUniq,
                       OverlapStore *ovlStoreRept,
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

  ////////////////////////////////////////
  //
  //  Read overlaps from the store

  uint32      ovlMax = MAX_OVERLAPS_PER_FRAG;
  uint32      ovlLen = 0;
  OVSoverlap *ovl    = new OVSoverlap [ovlMax];
  ovlLen  = 0;

  if (ovlStoreUniq) {
    AS_OVS_setRangeOverlapStore(ovlStoreUniq, frag.ident, frag.ident);
    ovlLen += AS_OVS_readOverlapsFromStore(ovlStoreUniq, ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
  }

  if (ovlStoreRept) {
    AS_OVS_setRangeOverlapStore(ovlStoreRept, frag.ident, frag.ident);
    ovlLen += AS_OVS_readOverlapsFromStore(ovlStoreRept, ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
  }

  if (ovlLen == 0) {
    delete [] ovl;
    return(false);
  }

  logFileFlags |= LOG_PLACE_FRAG;

  if (logFileFlagSet(LOG_PLACE_FRAG)) {
    fprintf(logFile, "placeFragUsingOverlaps()-- frag %d has %d overlaps\n", frag.ident, ovlLen);
  }

  ////////////////////////////////////////
  //
  //  Assign overlaps to unitigs, place fragments.

  overlapPlacement   *ovlPlace = new overlapPlacement[ovlLen];

  //for (uint32 i=0; i<ovlLen; i++)
  //  fprintf(logFile, "ovl iid %d %d hang %d %d\n", ovl[i].a_iid, ovl[i].b_iid, ovl[i].dat.ovl.a_hang, ovl[i].dat.ovl.b_hang);

  uint32   nOverlapsNotPlaced          = 0;
  uint32   nOverlapsNotPlacedContained = 0;

  uint32   nFragmentsNotPlaced         = 0;

  for (uint32 i=0; i<ovlLen; i++) {
    int32             utgID = Unitig::fragIn(ovl[i].b_iid);
    Unitig           *utg   = unitigs[utgID];

    assert(ovl[i].a_iid == frag.ident);

    if (utgID == 0) {
      //  Overlapping fragment not placed in a unitig yet.
      nOverlapsNotPlaced++;

      if (OG->isContained(ovl[i].b_iid))
        //  Overlapping fragment is contained, not expected to be placed?
        nOverlapsNotPlacedContained++;

      continue;
    }

    //  Depending on the type of overlap (containment vs dovetail), place the fragment relative to
    //  the other fragment.  The placeFrag() function is expecting the overlap to be from the
    //  container to the us fragment, which is opposite the overlap that we have.  We need to flip the
    //  fragments in the overlap -- negate the hangs.

    if        ((ovl[i].dat.ovl.a_hang >= 0) && (ovl[i].dat.ovl.b_hang <= 0)) {
      //  A (us) contains B (the other fragment)
      BestContainment  best;

      best.container       = ovl[i].b_iid;  //  Not really the container...
      best.isContained     = false;         //  ...so mark this as a false BestContainment
      best.a_hang          = ovl[i].dat.ovl.flipped ? ovl[i].dat.ovl.b_hang : -ovl[i].dat.ovl.a_hang;
      best.b_hang          = ovl[i].dat.ovl.flipped ? ovl[i].dat.ovl.a_hang : -ovl[i].dat.ovl.b_hang;
      best.sameOrientation = ovl[i].dat.ovl.flipped ? false : true;
      best.olapsSorted     = false;
      best.olapsLen        = 0;
      best.olaps           = NULL;

      if (utg->placeFrag(frag, &best)) {
        ovlPlace[i].frgID       = frag.ident;
        ovlPlace[i].tigID       = utgID;
        ovlPlace[i].position    = frag.position;
        ovlPlace[i].errors      = FI->fragmentLength(ovl[i].b_iid) * AS_OVS_decodeQuality(ovl[i].dat.ovl.corr_erate);
        ovlPlace[i].covered.bgn = ovl[i].dat.ovl.a_hang;
        ovlPlace[i].covered.end = FI->fragmentLength(ovl[i].a_iid) + ovl[i].dat.ovl.b_hang;
        ovlPlace[i].aligned     = ovlPlace[i].covered.end - ovlPlace[i].covered.bgn;
        assert(ovlPlace[i].covered.bgn < ovlPlace[i].covered.end);
#ifdef LOG_PLACE
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "placeFragUsingOverlaps()-- (container) - frag %d in unitig %d at %d,%d from overlap ident %d %d hang %d %d flipped %d covered %d,%d\n",
                  frag.ident, utgID, frag.position.bgn, frag.position.end,
                  ovl[i].a_iid, ovl[i].b_iid, ovl[i].dat.ovl.a_hang, ovl[i].dat.ovl.b_hang, ovl[i].dat.ovl.flipped,
                  ovlPlace[i].covered.bgn, ovlPlace[i].covered.end);
#endif
      } else {
        nFragmentsNotPlaced++;          
      }

    } else if ((ovl[i].dat.ovl.a_hang <= 0) && (ovl[i].dat.ovl.b_hang >= 0)) {
      //  A (us) is contained in B (the other fragment)
      BestContainment  best;

      best.container       = ovl[i].b_iid;
      best.isContained     = true;
      best.a_hang          = ovl[i].dat.ovl.flipped ? ovl[i].dat.ovl.b_hang : -ovl[i].dat.ovl.a_hang;
      best.b_hang          = ovl[i].dat.ovl.flipped ? ovl[i].dat.ovl.a_hang : -ovl[i].dat.ovl.b_hang;
      best.sameOrientation = ovl[i].dat.ovl.flipped ? false : true;
      best.olapsSorted     = false;
      best.olapsLen        = 0;
      best.olaps           = NULL;

      if (utg->placeFrag(frag, &best)) {
        ovlPlace[i].frgID       = frag.ident;
        ovlPlace[i].tigID       = utgID;
        ovlPlace[i].position    = frag.position;
        ovlPlace[i].errors      = FI->fragmentLength(ovl[i].a_iid) * AS_OVS_decodeQuality(ovl[i].dat.ovl.corr_erate);
        ovlPlace[i].covered.bgn = 0;
        ovlPlace[i].covered.end = FI->fragmentLength(ovl[i].a_iid);
        ovlPlace[i].aligned     = ovlPlace[i].covered.end - ovlPlace[i].covered.bgn;
        assert(ovlPlace[i].covered.bgn < ovlPlace[i].covered.end);
#ifdef LOG_PLACE
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "placeFragUsingOverlaps()-- (contained) - frag %d in unitig %d at %d,%d from overlap ident %d %d hang %d %d flipped %d covered %d,%d\n",
                  frag.ident, utgID, frag.position.bgn, frag.position.end,
                  ovl[i].a_iid, ovl[i].b_iid, ovl[i].dat.ovl.a_hang, ovl[i].dat.ovl.b_hang, ovl[i].dat.ovl.flipped,
                  ovlPlace[i].covered.bgn, ovlPlace[i].covered.end);
#endif
      } else {
        nFragmentsNotPlaced++;          
      }

    } else {
      //  A dovetail, use the existing placement routine
      BestEdgeOverlap   best;
      int32             plac3, plac5;
      int32             aend3p = AS_OVS_overlapAEndIs3prime(ovl[i]);

      best.set(ovl[i]);

      if (utg->placeFrag(frag, plac5, (aend3p ? NULL  : &best),
                         frag, plac3, (aend3p ? &best : NULL))) {
        uint32  olen = OG->olapLength(ovl[i].a_iid, ovl[i].b_iid, ovl[i].dat.ovl.a_hang, ovl[i].dat.ovl.b_hang);
        uint32  flen = FI->fragmentLength(ovl[i].a_iid);

        ovlPlace[i].frgID       = frag.ident;
        ovlPlace[i].tigID       = utgID;
        ovlPlace[i].position    = frag.position;
        ovlPlace[i].errors      = olen * AS_OVS_decodeQuality(ovl[i].dat.ovl.corr_erate);
        ovlPlace[i].covered.bgn = (ovl[i].dat.ovl.a_hang < 0) ?    0 : ovl[i].dat.ovl.a_hang;
        ovlPlace[i].covered.end = (ovl[i].dat.ovl.b_hang > 0) ? flen : ovl[i].dat.ovl.b_hang + flen;
        ovlPlace[i].aligned     = ovlPlace[i].covered.end - ovlPlace[i].covered.bgn;
        assert(ovlPlace[i].covered.bgn < ovlPlace[i].covered.end);
#ifdef LOG_PLACE
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "placeFragUsingOverlaps()-- (dovetail)  - frag %d in unitig %d at %d,%d from overlap ident %d %d hang %d %d flipped %d covered %d,%d\n",
                  frag.ident, utgID, frag.position.bgn, frag.position.end,
                  ovl[i].a_iid, ovl[i].b_iid, ovl[i].dat.ovl.a_hang, ovl[i].dat.ovl.b_hang, ovl[i].dat.ovl.flipped,
                  ovlPlace[i].covered.bgn, ovlPlace[i].covered.end);
#endif
      } else {
        nFragmentsNotPlaced++;          
      }
    }

    //  Disallow any placements that exceed the boundary of the unitig.  These cannot be confirmed
    //  by overlaps and might be wrong:
    //    o  sticking a unique/repeat fragment onto a repeat
    //    o  sticking a chimeric fragment onto the end of a unitig
    //
    if ((ovlPlace[i].position.bgn < 0) ||
        (ovlPlace[i].position.end < 0) ||
        (ovlPlace[i].position.bgn > utg->getLength()) ||
        (ovlPlace[i].position.end > utg->getLength())) {
#ifdef LOG_PLACE
      if (logFileFlagSet(LOG_PLACE_FRAG))
        fprintf(logFile, "placeFragUsingOverlaps()-- DISALLOW placement; it extends the unitig.\n");
#endif
      ovlPlace[i] = overlapPlacement();
    }
  }

  if (nFragmentsNotPlaced > 0)
    if (logFileFlagSet(LOG_PLACE_FRAG))
      fprintf(logFile, "placeFragUsingOverlaps()-- WARNING: Failed to place %d fragments\n", nFragmentsNotPlaced);

  ////////////////////////////////////////
  //
  //  Sort all the placements.  Any overlap we couldn't place is automatically in Unitig 0, the
  //  invalid unitig.  Sort order is by unitig ID, then by position in the unitig.

  sort(ovlPlace, ovlPlace + ovlLen);

  ////////////////////////////////////////
  //
  //  Examine those placements.  See if there is only one per unitig.
  //

  uint32         bgn = 0;
  uint32         end = 1;

  intervalList   UP;  //  Unitig position
  intervalList   FC;  //  Fragment coverage

  //  Simplify - skip the empty invalid alignents.  This could probably have been done inline
  //  but....they're sorted, and this is so much easier.
  while ((bgn < ovlLen) && (ovlPlace[bgn].tigID == 0))
    bgn++;


  while (bgn < ovlLen) {
    UP.clear();
    FC.clear();

    //  Find the last fragment in this unitig.
    end = bgn + 1;
    while ((end < ovlLen) && (ovlPlace[bgn].tigID == ovlPlace[end].tigID))
      end++;

    //  Populate an interval list with the placements, collapse.
    for (uint32 oo=bgn; oo<end; oo++) {
      int32   b = ovlPlace[oo].position.bgn;
      int32   e = ovlPlace[oo].position.end;

      int32   s = (b < e) ? (b)     : (e);
      int32   l = (b < e) ? (e - b) : (b - e);

      assert(ovlPlace[oo].tigID > 0);

      UP.add(s, l);
    }

    UP.merge();

#if 0
    if (logFileFlagSet(LOG_PLACE_FRAG))
      for (uint32 i=0; i<UP.numberOfIntervals(); i++)
        fprintf(logFile, "placeFragUsingOverlaps()-- frag %d in unitig %d (len %d nfrags %d) at "F_S64","F_S64" from "F_U32" overlaps\n",
                frag.ident,
                ovlPlace[bgn].tigID, unitigs[ovlPlace[bgn].tigID]->getLength(), unitigs[ovlPlace[bgn].tigID]->ufpath.size(),
                UP.lo(i),
                UP.hi(i),
                UP.ct(i));
#endif

    //  Segregate the overlaps by placement in the unitig.  We want to construct one overlapPlacement for each distinct placement.

    for (uint32 il=0, os=bgn, oe=bgn; os<end; il++) {
      //  Find the last overlap that belongs with intervalList il.

      assert(il < UP.numberOfIntervals());

      int32 ilbgn = UP.lo(il);
      int32 ilend = UP.hi(il);
      int32 novl  = UP.ct(il);

      overlapPlacement  op;

      while ((oe < end) && (ovlPlace[oe].position.end <= ilend))
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
      op.tigID = ovlPlace[os].tigID;

      op.fCoverage   = 0.0;

      op.errors      = 0.0;
      op.aligned     = 0;
      op.covered.bgn = ovlPlace[os].covered.bgn;
      op.covered.end = ovlPlace[os].covered.end;

      uint32  nForward = 0;
      uint32  nReverse = 0;

      for (uint32 oo=os; oo<oe; oo++) {
        op.errors      += ovlPlace[oo].errors;
        op.aligned     += ovlPlace[oo].aligned;
        op.covered.bgn  = MIN(op.covered.bgn, ovlPlace[oo].covered.bgn);
        op.covered.end  = MAX(op.covered.end, ovlPlace[oo].covered.end);

        if (isReverse(ovlPlace[oo].position))
          op.nReverse++;
        else
          op.nForward++;
      }

      op.position.bgn = (op.nReverse > 0) ? UP.hi(il) : UP.lo(il);
      op.position.end = (op.nReverse > 0) ? UP.lo(il) : UP.hi(il);

      op.fCoverage = (op.covered.end - op.covered.bgn) / (double)FI->fragmentLength(op.frgID);

      //  This placement is invalid if both nReverse and nForward are more than zero.

      if ((op.nReverse == 0) || (op.nForward == 0)) {
        placements.push_back(op);
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "placeFragUsingOverlaps()-- frag %d in unitig %d at %d,%d -- cov %.2f (%d,%d) errors %.2f aligned %d novl %d\n",
                  op.frgID, op.tigID, op.position.bgn, op.position.end,
                  op.fCoverage, op.covered.bgn, op.covered.end,
                  op.errors,
                  op.aligned,
                  novl);
      } else {
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "placeFragUsingOverlaps()-- frag %d in unitig %d at %d,%d -- cov %.2f (%d,%d) errors %.2f aligned %d novl %d -- INVALID nRev %d nFor %d\n",
                  op.frgID, op.tigID, op.position.bgn, op.position.end,
                  op.fCoverage, op.covered.bgn, op.covered.end,
                  op.errors,
                  op.aligned,
                  novl,
                  op.nReverse, op.nForward);
      }

      os = oe;
      oe = oe + 1;
    }  //  End of segregating overlaps by placement

    //  Move to the next block of overlaps.
    bgn = end;
    end = end + 1;
  }

  logFileFlags &= ~LOG_PLACE_FRAG;

  delete [] ovlPlace;
  delete [] ovl;

  return(true);
}

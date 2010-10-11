
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

static const char *rcsid = "$Id: AS_BOG_PlaceFragUsingOverlaps.cc,v 1.2 2010-10-11 03:43:44 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_Unitig.hh"

#include "AS_UTL_intervalList.H"

#define MAX_OVERLAPS_PER_FRAG   (16 * 1024 * 1024)


class overlapPlacement {
public:
  overlapPlacement() {
    tigID                 = 0;
    frg.ident             = 0;
    frg.contained         = 0;
    frg.parent            = 0;
    frg.ahang             = 0;
    frg.bhang             = 0;
    frg.position.bgn      = 0;
    frg.position.end      = 0;
    frg.containment_depth = 0;
  };
  ~overlapPlacement() {
  };

  bool operator<(const overlapPlacement &that) const {
    return(tigID < that.tigID);
  };

public:
  int32             tigID;
  ufNode            frg;      //  Position computed from the overlap
};



//  Given an implicit fragment -- a ufNode with only the 'ident' set -- this will compute the
//  best placement for the fragment in an existing unitig.  ALL overlaps are used, not just
//  the best.
//
//  Ties are broken using mate pairs, overlap identities, or arbitrarily.
//
//  Returns true if any placement is found, false otherwise.
//
bool
UnitigGraph::placeFragUsingOverlaps(ufNode frag,
                                    OverlapStore *ovlStoreUniq,
                                    OverlapStore *ovlStoreRept) {

  assert(frag.ident > 0);
  assert(frag.ident <= FI->numFragments());

  frag.contained         = 0;
  frag.parent            = 0;
  frag.ahang             = 0;
  frag.bhang             = 0;
  frag.position.bgn      = 0;
  frag.position.end      = 0;
  frag.containment_depth = 0;

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
    fprintf(logFile, "placeFragUsingOverlaps()--\n");
    fprintf(logFile, "placeFragUsingOverlaps()-- frag %d has %d overlaps\n", frag.ident, ovlLen);
  }

  //  Assign overlaps to unitigs, place fragments.

  overlapPlacement   *ovlPlace = new overlapPlacement[ovlLen];

  //for (uint32 i=0; i<ovlLen; i++) {
  //  fprintf(logFile, "ovl iid %d %d hang %d %d\n", ovl[i].a_iid, ovl[i].b_iid, ovl[i].dat.ovl.a_hang, ovl[i].dat.ovl.b_hang);
  //}

  uint32   nOverlapsNotPlaced          = 0;
  uint32   nOverlapsNotPlacedContained = 0;

  for (uint32 i=0; i<ovlLen; i++) {
    int32             utgID = Unitig::fragIn(ovl[i].b_iid);
    Unitig           *utg   = unitigs[utgID];

    assert(ovl[i].a_iid == frag.ident);

    if (utgID == 0) {
      //  Fragment not placed in a unitig yet
      nOverlapsNotPlaced++;

      if (OG->isContained(ovl[i].b_iid))
        //  Fragment is contained, not expected to be placed?
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

      best.container       = ovl[i].b_iid;
      best.isContained     = false;  //  Mark as a false BestContainment
      best.a_hang          = ovl[i].dat.ovl.flipped ? ovl[i].dat.ovl.b_hang : -ovl[i].dat.ovl.a_hang;
      best.b_hang          = ovl[i].dat.ovl.flipped ? ovl[i].dat.ovl.a_hang : -ovl[i].dat.ovl.b_hang;
      best.sameOrientation = ovl[i].dat.ovl.flipped ? false : true;
      best.isPlaced        = false;
      best.olapsSorted     = false;
      best.olapsLen        = 0;
      best.olaps           = NULL;

      frag.contained = ovl[i].b_iid;

      if (utg->placeFrag(frag, &best)) {
        ovlPlace[i].tigID = utgID;
        ovlPlace[i].frg   = frag;
#if 0
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "placeFragUsingOverlaps()-- (1) - frag %d in unitig %d at %d,%d from overlap ident %d %d hang %d %d flipped %d\n",
                  frag.ident, utgID, frag.position.bgn, frag.position.end,
                  ovl[i].a_iid, ovl[i].b_iid, ovl[i].dat.ovl.a_hang, ovl[i].dat.ovl.b_hang, ovl[i].dat.ovl.flipped);
#endif
      }

    } else if ((ovl[i].dat.ovl.a_hang <= 0) && (ovl[i].dat.ovl.b_hang >= 0)) {
      //  A (us) is contained in B (the other fragment)
      BestContainment  best;

      best.container       = ovl[i].b_iid;
      best.isContained     = true;
      best.a_hang          = ovl[i].dat.ovl.flipped ? ovl[i].dat.ovl.b_hang : -ovl[i].dat.ovl.a_hang;
      best.b_hang          = ovl[i].dat.ovl.flipped ? ovl[i].dat.ovl.a_hang : -ovl[i].dat.ovl.b_hang;
      best.sameOrientation = ovl[i].dat.ovl.flipped ? false : true;
      best.isPlaced        = false;
      best.olapsSorted     = false;
      best.olapsLen        = 0;
      best.olaps           = NULL;

      if (utg->placeFrag(frag, &best)) {
        ovlPlace[i].tigID = utgID;
        ovlPlace[i].frg   = frag;
#if 0
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "placeFragUsingOverlaps()-- (2) - frag %d in unitig %d at %d,%d from overlap ident %d %d hang %d %d flipped %d\n",
                  frag.ident, utgID, frag.position.bgn, frag.position.end,
                  ovl[i].a_iid, ovl[i].b_iid, ovl[i].dat.ovl.a_hang, ovl[i].dat.ovl.b_hang, ovl[i].dat.ovl.flipped);
#endif
      }

    } else {
      //  A dovetail, use the existing placement routine
      BestEdgeOverlap   best;
      int32             plac3, plac5;
      int32             aend3p = AS_OVS_overlapAEndIs3prime(ovl[i]);

      best.set(ovl[i]);

      //best.fragId = ovl[i].b_iid;
      //best.frag3p      = Frag3p(ovl[i]);
      //best.ahang     = ovl[i].dat.ovl.a_hang;
      //best.bhang     = ovl[i].dat.ovl.b_hang;

      if (utg->placeFrag(frag, plac5, (aend3p ? NULL  : &best),
                         frag, plac3, (aend3p ? &best : NULL))) {
        ovlPlace[i].tigID = utgID;
        ovlPlace[i].frg   = frag;
#if 0
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "placeFragUsingOverlaps()-- (3) - frag %d in unitig %d at %d,%d from overlap ident %d %d hang %d %d flipped %d\n",
                  frag.ident, utgID, frag.position.bgn, frag.position.end,
                  ovl[i].a_iid, ovl[i].b_iid, ovl[i].dat.ovl.a_hang, ovl[i].dat.ovl.b_hang, ovl[i].dat.ovl.flipped);
#endif
      }
    }
  }

  //  Sort all the placements.  Any overlap we couldn't place is automatically in Unitig 0, the
  //  invalid unitig.

  sort(ovlPlace, ovlPlace + ovlLen);

  //  Examine those placements.  See if there is only one per unitig.

  uint32         bgn = 0;
  uint32         end = 1;

  intervalList   UP;  //  Unitig position
  intervalList   FC;  //  Fragment coverage

  while (bgn < ovlLen) {
    UP.clear();
    FC.clear();

    //  Find the last fragment in this unitig.
    while ((end < ovlLen) && (ovlPlace[bgn].tigID == ovlPlace[end].tigID))
      end++;

    //  Populate an interval list with the placements, collapse.
    for (uint32 i=bgn; i<end; i++) {
      int32   b = ovlPlace[i].frg.position.bgn;
      int32   e = ovlPlace[i].frg.position.end;

      int32   s = (b < e) ? (b)     : (e);
      int32   l = (b < e) ? (e - b) : (b - e);

      if (ovlPlace[i].tigID == 0) {
        assert(b == 0);
        assert(e == 0);
      }

#if 0
      if ((ovlPlace[i].tigID != 0) &&
          (logFileFlagSet(LOG_PLACE_FRAG)))
        fprintf(logFile, "placeFragUsingOverlaps()-- frag %d in unitig %d at %d,%d\n",
                ovlPlace[i].frg.ident, ovlPlace[i].tigID, b, e);
#endif

      UP.add(s, l);
    }

    UP.merge();

    if (logFileFlagSet(LOG_PLACE_FRAG))
      for (uint32 i=0; i<UP.numberOfIntervals(); i++)
        fprintf(logFile, "placeFragUsingOverlaps()-- frag %d in unitig %d (len %d nfrags %d) at "F_S64","F_S64" from "F_U32" overlaps\n",
                frag.ident,
                ovlPlace[bgn].tigID, unitigs[ovlPlace[bgn].tigID]->getLength(), unitigs[ovlPlace[bgn].tigID]->ufpath.size(),
                UP.lo(i),
                UP.hi(i),
                UP.ct(i));

    bgn = end;
    end = end + 1;
  }

  {
    BestContainment *bestcont = OG->getBestContainer(frag.ident);

    //  Fragment has a best container, let's see where it put us.

    if (bestcont != NULL) {
      if (Unitig::fragIn(bestcont->container) > 0) {
        Unitig *utg = unitigs[Unitig::fragIn(bestcont->container)];

        if ((utg->placeFrag(frag, bestcont)) &&
            (logFileFlagSet(LOG_PLACE_FRAG)))
          fprintf(logFile, "placeFragUsingOverlaps()-- frag %d in unitig %d at %d,%d from best container\n",
                  frag.ident, utg->id(), frag.position.bgn, frag.position.end);
      } else {
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "placeFragUsingOverlaps()-- frag %d has unplaced best container\n",
                  frag.ident);
      }
    }
  }

  if (logFileFlagSet(LOG_PLACE_FRAG))
    fprintf(logFile, "placeFragUsingOverlaps()-- frag %d had %u unplacable overlaps, %u because they were to contained fragments\n",
            frag.ident, nOverlapsNotPlaced, nOverlapsNotPlacedContained);

  logFileFlags &= ~LOG_PLACE_FRAG;

  delete [] ovlPlace;
  delete [] ovl;

  return(false);
}

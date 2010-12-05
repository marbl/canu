
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

static const char *rcsid = "$Id: AS_BAT_Unitig_PlaceFragUsingEdges.C,v 1.2 2010-12-05 00:05:46 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#undef DEBUG_PLACE_FRAG

//  Given an implicit fragment -- a ufNode with only the 'ident' set -- and at least one best edge
//  to a fragment in this unitig, compute the position of the fragment in this unitig.  If both
//  edges are given, both will independently compute a placement, which might disagree.  It is up to
//  the client to figure out what to do in this case.
//
//  If a placement is not found for an edge, the corresponding bidx value is set to -1.  Otherwise,
//  it is set to the position in the fragment list of the fragment in this unitig (from above).
//
//  Returns true if any placement is found, false otherwise.
//
bool
Unitig::placeFrag(ufNode &frag5, int32 &bidx5, BestEdgeOverlap *bestedge5,
                  ufNode &frag3, int32 &bidx3, BestEdgeOverlap *bestedge3) {

  bidx5              = -1;
  bidx3              = -1;

  assert(frag5.ident > 0);
  assert(frag5.ident <= FI->numFragments());

  frag5.contained         = 0;
  frag5.parent            = 0;
  frag5.ahang             = 0;
  frag5.bhang             = 0;
  frag5.position.bgn      = 0;
  frag5.position.end      = 0;
  frag5.containment_depth = 0;

  assert(frag3.ident > 0);
  assert(frag3.ident <= FI->numFragments());

  frag3.contained         = 0;
  frag3.parent            = 0;
  frag3.ahang             = 0;
  frag3.bhang             = 0;
  frag3.position.bgn      = 0;
  frag3.position.end      = 0;
  frag3.containment_depth = 0;

  if ((bestedge5) && (bestedge5->fragId() == 0))
    bestedge5 = NULL;

  if ((bestedge3) && (bestedge3->fragId() == 0))
    bestedge3 = NULL;

  //  If we have an incoming edge, AND the fragment for that edge is in this unitig, look up its
  //  index.  Otherwise, discard the edge to prevent placement.
  //
  if ((bestedge5) && (fragIn(bestedge5->fragId()) == id())) {
    bidx5 = pathPosition(bestedge5->fragId());
    assert(bestedge5->fragId() == ufpath[bidx5].ident);
  } else {
    bestedge5 = NULL;
    bidx5     = -1;
  }

  if ((bestedge3) && (fragIn(bestedge3->fragId()) == id())) {
    bidx3 = pathPosition(bestedge3->fragId());;
    assert(bestedge3->fragId() == ufpath[bidx3].ident);
  } else {
    bestedge3 = NULL;
    bidx3     = -1;

  }
  //  Now, just compute the placement based on edges that exist.

  if ((bestedge5) && (bidx5 != -1)) {
    ufNode *parent = &ufpath[bidx5];

    assert(parent->ident == bestedge5->fragId());

    //  Overlap is stored using 'node' as the A frag, and we negate the hangs to make them relative
    //  to the 'parent'.  (This is opposite from how containment edges are saved.)  A special case
    //  exists when we overlap to the 5' end of the other fragment; we need to flip the overlap to
    //  ensure the (new) A frag is forward.

    int ahang = -bestedge5->ahang();
    int bhang = -bestedge5->bhang();

    if (bestedge5->frag3p() == false) {
      ahang = bestedge5->bhang();
      bhang = bestedge5->ahang();
    }

    int  pbgn, pend;
    int  fbgn, fend;

    //  Place the new fragment using the overlap.  We don't worry about the orientation of the new
    //  fragment, only the location.  Orientation of the parent fragment matters (1) to know which
    //  coordinate is the lower, and (2) to decide if the overlap needs to be flipped (again).

    if (parent->position.bgn < parent->position.end) {
      pbgn = parent->position.bgn;
      pend = parent->position.end;
      fbgn = parent->position.bgn + ahang;
      fend = parent->position.end + bhang;
    } else {
      pbgn = parent->position.end;
      pend = parent->position.bgn;
      fbgn = parent->position.end - bhang;
      fend = parent->position.bgn - ahang;
    }

    assert(pbgn < pend);
    assert(fbgn < fend);

    //  Since we don't know the true length of the overlap, if we use just the hangs to place a
    //  fragment, we typically shrink fragments well below their actual length.  In one case, we
    //  shrank a container enough that the containee was placed in the unitig backwards.
    //
    //  We now revert back to placing the end based on the actual length, but will
    //  adjust to maintain a dovetail relationship.
    //
    //  See comments on other instances of this warning.
    //

#warning not knowing the overlap length really hurts.
    fend = fbgn + FI->fragmentLength(frag5.ident);

    //  Make sure that we didn't just make a contained fragment out of a dovetail.  There are two
    //  cases here, either the fragment is before or after the parent.  We'll compare fbgn to the
    //  parent position.  We could use orientations and ends, but this is easier.

    if (fbgn < pbgn) {
      //  Fragment begins before the parent (fragment must therefore be forward since the edge is
      //  off the 5' end) and so fend must be before pend (otherwise fragment would contain parent).
      if (fend >= pend) {
#ifdef DEBUG_PLACE_FRAG
        fprintf(logFile, "RESET5l fend from %d to %d\n", fend, pend - 1);
#endif
        fend = pend - 1;
      }

    } else {
      if (fend <= pend) {
#ifdef DEBUG_PLACE_FRAG
        fprintf(logFile, "RESET5r fend from %d to %d\n", fend, pend + 1);
#endif
        fend = pend + 1;
      }
    }


    //  The new frag is reverse if:
    //    the old frag is forward and we hit its 5' end, or
    //    the old frag is reverse and we hit its 3' end.
    //
    //  The new frag is forward if:
    //    the old frag is forward and we hit its 3' end, or
    //    the old frag is reverse and we hit its 5' end.
    //
    bool flip = (((parent->position.bgn < parent->position.end) && (bestedge5->frag3p() == false)) ||
                 ((parent->position.end < parent->position.bgn) && (bestedge5->frag3p() == true)));

#ifdef DEBUG_PLACE_FRAG
    fprintf(logFile, "bestedge5:  parent iid %d pos %d,%d   b_iid %d ovl %d,%d,%d  pos %d,%d flip %d\n",
            parent->ident, parent->position.bgn, parent->position.end,
            bestedge5->fragId(), bestedge5->frag3p(), bestedge5->ahang(), bestedge5->bhang(), fbgn, fend, flip);
#endif

    frag5.contained    = 0;
    frag5.parent       = bestedge5->fragId();
    frag5.ahang        = ahang;
    frag5.bhang        = bhang;
    frag5.position.bgn = (flip) ? fend : fbgn;
    frag5.position.end = (flip) ? fbgn : fend;
  }


  if ((bestedge3) && (bidx3 != -1)) {
    ufNode *parent = &ufpath[bidx3];

    assert(parent->ident == bestedge3->fragId());

    int ahang = -bestedge3->ahang();
    int bhang = -bestedge3->bhang();

    if (bestedge3->frag3p() == true) {
      ahang = bestedge3->bhang();
      bhang = bestedge3->ahang();
    }

    int  pbgn, pend;
    int  fbgn, fend;

    if (parent->position.bgn < parent->position.end) {
      pbgn = parent->position.bgn;
      pend = parent->position.end;
      fbgn = parent->position.bgn + ahang;
      fend = parent->position.end + bhang;
    } else {
      pbgn = parent->position.end;
      pend = parent->position.bgn;
      fbgn = parent->position.end - bhang;
      fend = parent->position.bgn - ahang;
    }

    assert(pbgn < pend);
    assert(fbgn < fend);

#warning not knowing the overlap length really hurts.
    fend = fbgn + FI->fragmentLength(frag3.ident);

    //  Make sure that we didn't just make a contained fragment out of a dovetail.  There are two
    //  cases here, either the fragment is before or after the parent.  We'll compare fbgn to the
    //  parent position.  We could use orientations and ends, but this is easier.

    if (fbgn < pbgn) {
      //  Fragment begins before the parent (fragment must therefore be reverse since the edge is
      //  off the 3' end) and so fend must be before pend (otherwise fragment would contain parent).
      if (fend >= pend) {
#ifdef DEBUG_PLACE_FRAG
        fprintf(logFile, "RESET3l fend from %d to %d\n", fend, pend-1);
#endif
        fend = pend - 1;
      }

    } else {
      if (fend <= pend) {
#ifdef DEBUG_PLACE_FRAG
        fprintf(logFile, "RESET3r fend from %d to %d\n", fend, pend+1);
#endif
        fend = pend + 1;
      }
    }


    //  The new frag is reverse if:
    //    the old frag is forward and we hit its 3' end, or
    //    the old frag is reverse and we hit its 5' end.
    //
    //  The new frag is forward if:
    //    the old frag is forward and we hit its 5' end, or
    //    the old frag is reverse and we hit its 3' end.
    //
    bool flip = (((parent->position.bgn < parent->position.end) && (bestedge3->frag3p() == true)) ||
                 ((parent->position.end < parent->position.bgn) && (bestedge3->frag3p() == false)));

#ifdef DEBUG_PLACE_FRAG
    fprintf(logFile, "bestedge3:  parent iid %d pos %d,%d   b_iid %d ovl %d,%d,%d  pos %d,%d flip %d\n",
            parent->ident, parent->position.bgn, parent->position.end,
            bestedge3->fragId(), bestedge3->frag3p(), bestedge3->ahang(), bestedge3->bhang(), fbgn, fend, flip);
#endif

    frag3.contained    = 0;
    frag3.parent       = bestedge3->fragId();
    frag3.ahang        = ahang;
    frag3.bhang        = bhang;
    frag3.position.bgn = (flip) ? fend : fbgn;
    frag3.position.end = (flip) ? fbgn : fend;
  }

  return((bidx5 != -1) || (bidx3 != -1));
}



bool
Unitig::placeFrag(ufNode &frag, BestContainment *bestcont) {

  frag.contained         = 0;
  frag.parent            = 0;
  frag.ahang             = 0;
  frag.bhang             = 0;
  frag.position.bgn      = 0;
  frag.position.end      = 0;
  frag.containment_depth = 0;

  ufNode *parent = &ufpath[pathPosition(bestcont->container)];

#if 0
  //  This block is useful for debugging (maybe).  It is usually triggered only during popBubbles(),
  //  when we try to place a contained fragment into a fragment that has not been moved into the new
  //  unitig yet.  It might be useful if pathPosition ever gets messed up.
  //
  if ((parent == NULL) || (parent->ident != bestcont->container)) {
    ufNode *found = parent;

    for (int fi=0; fi<ufpath.size(); fi++)
      if (ufpath[fi].ident == bestcont->container)
        parent = &ufpath[fi];

    if (parent) {
      fprintf(logFile, "WARNING:  Didn't find the correct parent frag (%d) for contained frag %d -- pathPosition screwed up.\n",
              bestcont->container, frag.ident);
      fprintf(logFile, "          Found frag %d instead.\n", (parent == NULL) ? -1 : parent->ident);

      for (int fi=0; fi<ufpath.size(); fi++) {
        ufNode *ix = &ufpath[fi];

        fprintf(logFile, "          path[%4d,%4d] is frag %d %s\n",
                fi, pathPosition(ix->ident),
                ix->ident,
                (ix->ident == bestcont->container) ? " CORRECT PARENT!" : "");
      }
    }
  }
#endif

  if ((parent == NULL) || (parent->ident != bestcont->container)) {
    fprintf(logFile, "Unitig::placeFrag()-- WARNING:  Failed to place frag %d into unitig %d; parent not here.\n",
            frag.ident, id());
    return(false);
  }
  
  //  Adjust orientation.
  //
  //  NOTE!  The hangs are from the (parent) container to the (child)
  //  containee.  This is opposite as to how dovetail edges are stored.

  //  isContained
  //    If true, this is a true BestContainment; this fragment is contained in bestcont->container.
  //
  //    If false, the containment relationship is inverted; this fragment contains
  //    bestcont->contaier.  This is used when placing fragments by overlaps,
  //    and cannot be added to a unitig directly.

  if (bestcont->isContained) {
    assert(bestcont->a_hang >= 0);
    assert(bestcont->b_hang <= 0);
  } else {
    assert(bestcont->a_hang <= 0);
    assert(bestcont->b_hang >= 0);
  }

  if (parent->position.bgn < parent->position.end) {
    //  Container is forward.
    frag.contained = bestcont->container;
    frag.parent    = bestcont->container;
    frag.ahang     = bestcont->a_hang;
    frag.bhang     = bestcont->b_hang;

    if (bestcont->sameOrientation) {
      //  ...and so is containee.
      frag.position.bgn = parent->position.bgn + frag.ahang;
      frag.position.end = parent->position.end + frag.bhang;
    } else {
      //  ...but containee is reverse.
      frag.position.bgn = parent->position.end + frag.bhang;
      frag.position.end = parent->position.bgn + frag.ahang;
    }

  } else {
    //  Container is reverse.
    frag.contained = bestcont->container;
    frag.parent    = bestcont->container;
    frag.ahang     = -bestcont->b_hang;
    frag.bhang     = -bestcont->a_hang;

    if (bestcont->sameOrientation) {
      //  ...and so is containee.
      frag.position.bgn = parent->position.bgn + frag.bhang;
      frag.position.end = parent->position.end + frag.ahang;
    } else {
      //  ...but containee is forward.
      frag.position.bgn = parent->position.end + frag.ahang;
      frag.position.end = parent->position.bgn + frag.bhang;
    }
  }


  if (bestcont->isContained == false)
    //  If we're a false containment overlap, skip the adjustment below.
    return(true);

  //  Containments are particularily painful.  A beautiful example: a fragment of length 253bp is
  //  contained in a fragment of length 251bp (both hangs are zero).  In this case, the
  //  "ahang+length" method fails, placing the contained fragment outside the container (and if
  //  backwards oriented, _BEFORE_ the contained fragment).  The "ahang,bhang" method works here,
  //  but fails on other instances, shrinking deep containments to nothing.
  //
  //  We can use either method first, then adjust using the other method.
  //
  //  We'll use 'ahang,bhang' first (mostly because it was already done, and we need to compute
  //  those values anyway) then reset the end based on the length, limited to maintain a containment
  //  relationship.
  //
#warning not knowing the overlap length really hurts.
  if (frag.position.bgn < frag.position.end) {
    frag.position.end = frag.position.bgn + FI->fragmentLength(frag.ident);
    if (frag.position.end > MAX(parent->position.bgn, parent->position.end))
      frag.position.end = MAX(parent->position.bgn, parent->position.end);
  } else {
    frag.position.bgn = frag.position.end + FI->fragmentLength(frag.ident);
    if (frag.position.bgn > MAX(parent->position.bgn, parent->position.end))
      frag.position.bgn = MAX(parent->position.bgn, parent->position.end);
  }

  //  So we can sort properly, set the depth of this contained fragment.
  frag.containment_depth = parent->containment_depth + 1;

  return(true);
}

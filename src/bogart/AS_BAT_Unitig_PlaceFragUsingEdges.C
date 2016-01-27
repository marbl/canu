
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
 *    src/AS_BAT/AS_BAT_Unitig_PlaceFragUsingEdges.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2015-MAR-06
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

//  This provides low level (and usually too much) detail on placing a read using an edge.
#undef DEBUG_PLACE_FRAG



void
Unitig::placeFrag_computePlacement(ufNode          &frag,
                                   int32           &bidx,
                                   BestEdgeOverlap *bestedge,
                                   bool             bestIs3) {
  ufNode *parent = &ufpath[bidx];

  assert(parent->ident == bestedge->fragId());

  //  Scale the hangs based on the placement of the parent read.  This isn't perfect; we should really only
  //  scale the hang that is into the parent read (either positive A or negative B) and let the other
  //  hang be based on the scaling for this read -- but we don't know the scaling for this read.

  uint32  parentPlacedLen = (parent->position.bgn < parent->position.end) ? (parent->position.end - parent->position.bgn) : (parent->position.bgn - parent->position.end);
  uint32  parentRealLen   = FI->fragmentLength(parent->ident);

  double  intraScale = (double)parentPlacedLen / parentRealLen;  //  Within the parent read overlap
  double  interScale = 1.0;                                      //  Outside the parent read overlap

  //  Overlap is stored using 'node' as the A frag, and we negate the hangs to make them relative
  //  to the 'parent'.  (This is opposite from how containment edges are saved.)  A special case
  //  exists when we overlap to the 5' end of the other fragment; we need to flip the overlap to
  //  ensure the (new) A frag is forward.

  int32 ahang = -bestedge->ahang();
  int32 bhang = -bestedge->bhang();

  if (bestedge->frag3p() == bestIs3) {
    ahang = bestedge->bhang();
    bhang = bestedge->ahang();
  }

  int32  bgnhang = 0;
  int32  endhang = 0;

  int32  pbgn, pend;
  int32  fbgn, fend;

  bool   adjustBgn = false;

  //  Place the new fragment using the overlap.  We don't worry about the orientation of the new
  //  fragment, only the location.  Orientation of the parent fragment matters (1) to know which
  //  coordinate is the lower, and (2) to decide if the overlap needs to be flipped (again).

  if (parent->position.bgn < parent->position.end) {
    pbgn = parent->position.bgn;
    pend = parent->position.end;

    bgnhang = ahang;
    endhang = bhang;

  } else {
    pbgn = parent->position.end;
    pend = parent->position.bgn;

    bgnhang = -bhang;
    endhang = -ahang;
  }

  assert((bgnhang >=0) == (endhang >= 0));

  if (bgnhang > 0) {
    fbgn = pbgn + bgnhang * intraScale;  //  hang is moving low  to the right, inside the parent
    fend = pend + endhang * interScale;  //  hang is moving high to the right, outside the parent
  } else {
    fbgn = pbgn + bgnhang * interScale;  //  hang is moving low  to the left,  outside the parent
    fend = pend + endhang * intraScale;  //  hang is moving high to the left,  inside the parent
  }


  //  Since we don't know the true length of the overlap, if we use just the hangs to place a
  //  fragment, we typically shrink fragments well below their actual length.  In one case, we
  //  shrank a container enough that the containee was placed in the unitig backwards.
  //
  //  We now revert back to placing the end based on the actual length, but will
  //  adjust to maintain a dovetail relationship.
  //
  //  See comments on other instances of this warning.

#warning not knowing the overlap length really hurts.

#if 1

  //  If true, we've moved fend outside the parent range, so that can be adjusted.
  //  If false, the begin point can be adjusted.

  if (bgnhang > 0)
    fend = fbgn + FI->fragmentLength(frag.ident);
  else
    fbgn = fend - FI->fragmentLength(frag.ident);

#else
  //  This was an attempt to adjust position to better capture the length of the read.  It bombed
  //  because it violates hang restrictions - for example, when building the initial unitig, the
  //  path is from positive ahang overlaps, but this change can result in negative hang
  //  postioning.
  //
  //  This typically fails with
  //    AS_BAT_Unitig_AddFrag.C:68:
  //      void Unitig::addFrag(ufNode, int, bool):
  //        Assertion `node.position.end >= 0' failed
  //
  int  fpos = (fbgn + fend) / 2;

  fbgn = fpos - FI->fragmentLength(frag.ident) / 2;
  fend = fpos + FI->fragmentLength(frag.ident) / 2;
#endif

  //  Make sure that we didn't just make a contained fragment out of a dovetail.  There are two
  //  cases here, either the fragment is before or after the parent.  We'll compare fbgn to the
  //  parent position.  We could use orientations and ends, but this is easier.

  if (fbgn < pbgn) {
    if (fend >= pend) {
      fend = pend - 1;
    }

  } else {
    if (fend <= pend) {
      fend = pend + 1;
    }
  }

  if (pbgn >= pend)
    writeLog("placeFrag()-- ERROR: %c' parent placement inconsistent iid=%d %d,%d\n",
             (bestIs3) ? '3' : '5', parent->ident, parent->position.bgn, parent->position.end);
  if (fbgn >= fend)
    writeLog("placeFrag()-- ERROR: %c' placement inconsistent parent=%d %d,%d hang %d,%d this %d %d,%d\n",
             (bestIs3) ? '3' : '5', parent->ident, parent->position.bgn, parent->position.end,
             ahang, bhang,
             frag.ident, fbgn, fend);

  //if ((pbgn >= pend) || (fbgn >= fend))
  //  return(false);

  assert(pbgn < pend);
  assert(fbgn < fend);


  //  The new frag is reverse if:
  //    the old frag is forward and we hit its 5' end, or
  //    the old frag is reverse and we hit its 3' end.
  //
  //  The new frag is forward if:
  //    the old frag is forward and we hit its 3' end, or
  //    the old frag is reverse and we hit its 5' end.
  //
  bool flip = (((parent->position.bgn < parent->position.end) && (bestedge->frag3p() == bestIs3)) ||
               ((parent->position.end < parent->position.bgn) && (bestedge->frag3p() != bestIs3)));

#ifdef DEBUG_PLACE_FRAG
  writeLog("placeFrag()-- bestedge:  parent iid %d pos %d,%d   b_iid %d ovl %d,%d,%d  pos %d,%d flip %d\n",
           parent->ident, parent->position.bgn, parent->position.end,
           bestedge->fragId(), bestedge->frag3p(), bestedge->ahang(), bestedge->bhang(), fbgn, fend, flip);
#endif

  frag.contained    = 0;
  frag.parent       = bestedge->fragId();
  frag.ahang        = ahang;
  frag.bhang        = bhang;
  frag.position.bgn = (flip) ? fend : fbgn;
  frag.position.end = (flip) ? fbgn : fend;
}





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

  if ((bestedge5) && (bidx5 != -1))
    placeFrag_computePlacement(frag5, bidx5, bestedge5, false);

  if ((bestedge3) && (bidx3 != -1))
    placeFrag_computePlacement(frag3, bidx3, bestedge3, true);

  //  Return success if we computed.

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
      writeLog("placeFrag()-- WARNING:  Didn't find the correct parent frag (%d) for contained frag %d -- pathPosition screwed up.\n",
               bestcont->container, frag.ident);
      writeLog("placeFrag()--           Found frag %d instead.\n", (parent == NULL) ? -1 : parent->ident);

      for (int fi=0; fi<ufpath.size(); fi++) {
        ufNode *ix = &ufpath[fi];

        writeLog("placeFrag()--           path[%4d,%4d] is frag %d %s\n",
                 fi, pathPosition(ix->ident),
                 ix->ident,
                 (ix->ident == bestcont->container) ? " CORRECT PARENT!" : "");
      }
    }
  }
#endif

  if ((parent == NULL) || (parent->ident != bestcont->container)) {
#ifdef DEBUG_PLACE_FRAG
    writeLog("placeFrag()-- WARNING:  Failed to place frag %d into unitig %d; parent not here.\n",
             frag.ident, id());
#endif
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

  //  Scale the hangs to the current placed read length.  int32 will overflow for reads > 15 bits, so we use a double.

  if (parent->position.bgn < parent->position.end) {
    //  Container is forward.
    frag.contained = bestcont->container;
    frag.parent    = bestcont->container;
    frag.ahang     = bestcont->a_hang;
    frag.bhang     = bestcont->b_hang;

    double  scale = (double)(parent->position.end - parent->position.bgn) / FI->fragmentLength(parent->ident);

    if ((scale < 0.75) ||
        (1.25  < scale))
      writeLog("placeFrag()-- extreme scaling FWD %d at %d,%d (%d) len %d by %f\n",
               parent->ident, parent->position.bgn, parent->position.end, parent->position.end - parent->position.bgn, FI->fragmentLength(parent->ident), scale);

    if (bestcont->sameOrientation) {
      //  ...and so is containee.
      frag.position.bgn = parent->position.bgn + frag.ahang * scale;
      frag.position.end = parent->position.end + frag.bhang * scale;
    } else {
      //  ...but containee is reverse.
      frag.position.bgn = parent->position.end + frag.bhang * scale;
      frag.position.end = parent->position.bgn + frag.ahang * scale;
    }

  } else {
    //  Container is reverse.
    frag.contained = bestcont->container;
    frag.parent    = bestcont->container;
    frag.ahang     = -bestcont->b_hang;
    frag.bhang     = -bestcont->a_hang;

    double  scale = (double)(parent->position.bgn - parent->position.end) / FI->fragmentLength(parent->ident);

    if ((scale < 0.75) ||
        (1.25  < scale))
      writeLog("placeFrag()-- extreme scaling REV %d at %d,%d (%d) len %d by %f\n",
               parent->ident, parent->position.end, parent->position.bgn, parent->position.bgn - parent->position.end, FI->fragmentLength(parent->ident), scale);

    if (bestcont->sameOrientation) {
      //  ...and so is containee.
      frag.position.bgn = parent->position.bgn + frag.bhang * scale;
      frag.position.end = parent->position.end + frag.ahang * scale;
    } else {
      //  ...but containee is forward.
      frag.position.bgn = parent->position.end + frag.ahang * scale;
      frag.position.end = parent->position.bgn + frag.bhang * scale;
    }
  }


  if (bestcont->isContained == false)
    //  If we're a false containment overlap, skip the adjustment below.
    return(true);


#ifdef DEBUG_PLACE_FRAG
  writeLog("placeFrag()-- contained frag %u at %d,%d -- hangs %d,%d\n",
           frag.ident, frag.position.bgn, frag.position.end, frag.ahang, frag.bhang);
#endif

  //  Reset the position.  Try to accomodate the hangs and full read length.
  //  Note that this CAN break containment relationships.

  int32  fragPos   = (frag.position.bgn + frag.position.end) / 2;
  int32  placedLen = 0;

  if (frag.position.bgn < frag.position.end) {
    int32  placedLen = frag.position.end - frag.position.bgn;
    int32  aveLen    = (placedLen + FI->fragmentLength(frag.ident)) / 2;

    frag.position.bgn = fragPos - aveLen / 2;
    frag.position.end = fragPos + aveLen / 2;

#ifdef DEBUG_PLACE_FRAG
    writeLog("placeFrag()-- contained frag %u at %d,%d fwd from parent frag %u at %d,%d placedLen %d readLen %d aveLen %d\n",
             frag.ident,    frag.position.bgn,    frag.position.end,
             parent->ident, parent->position.bgn, parent->position.end,
             placedLen, FI->fragmentLength(frag.ident), aveLen);
#endif

  } else {
    int32  placedLen = frag.position.bgn - frag.position.end;
    int32  aveLen    = (placedLen + FI->fragmentLength(frag.ident)) / 2;

    frag.position.bgn = fragPos + aveLen / 2;
    frag.position.end = fragPos - aveLen / 2;

#ifdef DEBUG_PLACE_FRAG
    writeLog("placeFrag()-- contained frag %u at %d,%d rev from parent frag %u at %d,%d placedLen %d readLen %d aveLen %d\n",
             frag.ident,    frag.position.bgn,    frag.position.end,
             parent->ident, parent->position.bgn, parent->position.end,
             placedLen, FI->fragmentLength(frag.ident), aveLen);
#endif
  }

  //  If we're pushed outside the container, adjust.

  int32   minParent = MIN(parent->position.bgn, parent->position.end);
  int32   maxParent = MAX(parent->position.bgn, parent->position.end);

  //writeLog("min/max %d %d frag %d %d\n", minParent, maxParent, frag.position.bgn, frag.position.end);

  if (frag.position.bgn < minParent)   frag.position.bgn = minParent;
  if (frag.position.end < minParent)   frag.position.end = minParent;

  if (frag.position.bgn > maxParent)   frag.position.bgn = maxParent;
  if (frag.position.end > maxParent)   frag.position.end = maxParent;

  //writeLog("min/max %d %d frag %d %d\n", minParent, maxParent, frag.position.bgn, frag.position.end);

  assert(frag.position.bgn >= 0);
  assert(frag.position.end >= 0);
  assert(frag.position.bgn <= getLength());
  assert(frag.position.end <= getLength());

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
#if 0
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
#endif

  //  So we can sort properly, set the depth of this contained fragment.
  frag.containment_depth = parent->containment_depth + 1;

  return(true);
}

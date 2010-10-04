
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

static const char *rcsid = "$Id: AS_BOG_Unitig_PlaceFragUsingEdges.cc,v 1.2 2010-10-04 16:37:47 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_Unitig.hh"
#include "AS_BOG_BestOverlapGraph.hh"



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

  //  If we have an incoming edge, AND the fragment for that edge is in this unitig, look up its
  //  index.  Otherwise, discard the edge to prevent placement.
  //
  if ((bestedge5) && (fragIn(bestedge5->frag_b_id) == id())) {
    bidx5 = pathPosition(bestedge5->frag_b_id);
    assert(bestedge5->frag_b_id == ufpath[bidx5].ident);
  } else {
    bestedge5 = NULL;
    bidx5     = -1;
  }

  if ((bestedge3) && (fragIn(bestedge3->frag_b_id) == id())) {
    bidx3 = pathPosition(bestedge3->frag_b_id);;
    assert(bestedge3->frag_b_id == ufpath[bidx3].ident);
  } else {
    bestedge3 = NULL;
    bidx3     = -1;

  }
  //  Now, just compute the placement based on edges that exist.

  if ((bestedge5) && (bidx5 != -1)) {
    ufNode *parent = &ufpath[bidx5];

    assert(parent->ident == bestedge5->frag_b_id);

    //  Overlap is stored using 'node' as the A frag, and we negate the hangs to make them relative
    //  to the 'parent'.  (This is opposite from how containment edges are saved.)  A special case
    //  exists when we overlap to the 5' end of the other fragment; we need to flip the overlap to
    //  ensure the (new) A frag is forward.

    int ahang = -bestedge5->ahang;
    int bhang = -bestedge5->bhang;

    if (bestedge5->bend == FIVE_PRIME) {
      ahang = bestedge5->bhang;
      bhang = bestedge5->ahang;
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
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "RESET5l fend from %d to %d\n", fend, pend - 1);
        fend = pend - 1;
      }

    } else {
      if (fend <= pend) {
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "RESET5r fend from %d to %d\n", fend, pend + 1);
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
    bool flip = (((parent->position.bgn < parent->position.end) && (bestedge5->bend == FIVE_PRIME)) ||
                 ((parent->position.end < parent->position.bgn) && (bestedge5->bend == THREE_PRIME)));

    if (logFileFlagSet(LOG_PLACE_FRAG))
      fprintf(logFile, "bestedge5:  parent iid %d pos %d,%d   b_iid %d ovl %d,%d,%d  pos %d,%d flip %d\n",
              parent->ident, parent->position.bgn, parent->position.end,
              bestedge5->frag_b_id, bestedge5->bend, bestedge5->ahang, bestedge5->bhang, fbgn, fend, flip);

    frag5.contained    = 0;
    frag5.parent       = bestedge5->frag_b_id;
    frag5.ahang        = ahang;
    frag5.bhang        = bhang;
    frag5.position.bgn = (flip) ? fend : fbgn;
    frag5.position.end = (flip) ? fbgn : fend;
  }


  if ((bestedge3) && (bidx3 != -1)) {
    ufNode *parent = &ufpath[bidx3];

    assert(parent->ident == bestedge3->frag_b_id);

    int ahang = -bestedge3->ahang;
    int bhang = -bestedge3->bhang;

    if (bestedge3->bend == THREE_PRIME) {
      ahang = bestedge3->bhang;
      bhang = bestedge3->ahang;
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
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "RESET3l fend from %d to %d\n", fend, pend-1);
        fend = pend - 1;
      }

    } else {
      if (fend <= pend) {
        if (logFileFlagSet(LOG_PLACE_FRAG))
          fprintf(logFile, "RESET3r fend from %d to %d\n", fend, pend+1);
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
    bool flip = (((parent->position.bgn < parent->position.end) && (bestedge3->bend == THREE_PRIME)) ||
                 ((parent->position.end < parent->position.bgn) && (bestedge3->bend == FIVE_PRIME)));

    if (logFileFlagSet(LOG_PLACE_FRAG))
      fprintf(logFile, "bestedge3:  parent iid %d pos %d,%d   b_iid %d ovl %d,%d,%d  pos %d,%d flip %d\n",
              parent->ident, parent->position.bgn, parent->position.end,
              bestedge3->frag_b_id, bestedge3->bend, bestedge3->ahang, bestedge3->bhang, fbgn, fend, flip);

    frag3.contained    = 0;
    frag3.parent       = bestedge3->frag_b_id;
    frag3.ahang        = ahang;
    frag3.bhang        = bhang;
    frag3.position.bgn = (flip) ? fend : fbgn;
    frag3.position.end = (flip) ? fbgn : fend;
  }

  return((bidx5 != -1) || (bidx3 != -1));
}

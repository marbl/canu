
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

static const char *rcsid = "$Id: AS_BOG_Unitig_PlaceFrag.cc,v 1.1 2010-10-04 04:41:33 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_Unitig.hh"
#include "AS_BOG_BestOverlapGraph.hh"



//  Given an implicit fragment, and at least one best edge to a fragment in this unitig, compute the
//  position of the fragment in this unitig.  If both edges are given, both will independently
//  compute a placement, which might disagree.
//
//  If a placement is not found for an edge, the corresponding bidx value is set to -1.  Otherwise,
//  it is set to the position in the fragment list of the fragment in this unitig (from above).
//
//  Returns true if any placement is found, false otherwise.
//  The bidx value is set to -1 if no placement is found for that end.
//
bool
Unitig::placeFrag(ufNode &frag5, int32 &bidx5, BestEdgeOverlap *bestedge5,
                  ufNode &frag3, int32 &bidx3, BestEdgeOverlap *bestedge3) {
  bool  verbose = false;

  //frag5.ident
  frag5.contained    = 0;
  frag5.parent       = 0;
  frag5.ahang        = 0;
  frag5.bhang        = 0;
  frag5.position.bgn = 0;
  frag5.position.end = 0;

  //frag3.ident
  frag3.contained    = 0;
  frag3.parent       = 0;
  frag3.ahang        = 0;
  frag3.bhang        = 0;
  frag3.position.bgn = 0;
  frag3.position.end = 0;

  assert(frag3.ident > 0);
  assert(frag5.ident > 0);

  assert(frag3.ident <= FI->numFragments());
  assert(frag5.ident <= FI->numFragments());

  bidx5              = -1;
  bidx3              = -1;

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

    int  bgn, end;

    //  Place the new fragment using the overlap.  We don't worry about the orientation of the new
    //  fragment, only the location.  Orientation of the parent fragment matters (1) to know which
    //  coordinate is the lower, and (2) to decide if the overlap needs to be flipped (again).

    if (parent->position.bgn < parent->position.end) {
      bgn = parent->position.bgn + ahang;
      end = parent->position.end + bhang;
    } else {
      bgn = parent->position.end - bhang;
      end = parent->position.bgn - ahang;
    }

    assert(bgn < end);

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
    end = bgn + FI->fragmentLength(frag5.ident);
    if (end <= MAX(parent->position.bgn, parent->position.end))
      end = MAX(parent->position.bgn, parent->position.end) + 1;

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

    if (verbose)
      fprintf(logFile, "bestedge5:  parent iid %d pos %d,%d   b_iid %d ovl %d,%d,%d  pos %d,%d  flip %d\n",
              parent->ident, parent->position.bgn, parent->position.end,
              bestedge5->frag_b_id, bestedge5->bend, bestedge5->ahang, bestedge5->bhang, bgn, end, flip);

    frag5.contained    = 0;
    frag5.parent       = bestedge5->frag_b_id;
    frag5.ahang        = ahang;
    frag5.bhang        = bhang;
    frag5.position.bgn = (flip) ? end : bgn;
    frag5.position.end = (flip) ? bgn : end;
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

    int  bgn, end;

    if (parent->position.bgn < parent->position.end) {
      bgn = parent->position.bgn + ahang;
      end = parent->position.end + bhang;
    } else {
      bgn = parent->position.end - bhang;
      end = parent->position.bgn - ahang;
    }

    assert(bgn < end);

#warning not knowing the overlap length really hurts.
    end = bgn + FI->fragmentLength(frag3.ident);
    if (end <= MAX(parent->position.bgn, parent->position.end))
      end = MAX(parent->position.bgn, parent->position.end) + 1;

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

    if (verbose)
      fprintf(logFile, "bestedge3:  parent iid %d %d,%d   b_iid %d ovl %d,%d,%d  pos %d,%d flip %d\n",
              parent->ident, parent->position.bgn, parent->position.end,
              bestedge3->frag_b_id, bestedge3->bend, bestedge3->ahang, bestedge3->bhang, bgn, end, flip);

    frag3.contained    = 0;
    frag3.parent       = bestedge3->frag_b_id;
    frag3.ahang        = ahang;
    frag3.bhang        = bhang;
    frag3.position.bgn = (flip) ? end : bgn;
    frag3.position.end = (flip) ? bgn : end;
  }

  return((bidx5 != -1) || (bidx3 != -1));
}

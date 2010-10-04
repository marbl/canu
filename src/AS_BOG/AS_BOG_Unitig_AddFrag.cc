
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

static const char *rcsid = "$Id: AS_BOG_Unitig_AddFrag.cc,v 1.1 2010-10-04 04:41:33 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_Unitig.hh"
#include "AS_BOG_BestOverlapGraph.hh"




void
Unitig::addFrag(ufNode node, int offset, bool report) {

  node.position.bgn += offset;
  node.position.end += offset;

  assert(node.ident > 0);

  // keep track of the unitig a frag is in
  _inUnitig[node.ident]     = _id;
  _pathPosition[node.ident] = ufpath.size();

  // keep track of max position in unitig
  int32 frgEnd = MAX(node.position.bgn, node.position.end);
  if (frgEnd > _length)
    _length = frgEnd;

  ufpath.push_back(node);

  if ((report) || (node.position.bgn < 0) || (node.position.end < 0)) {
    int32 len = FI->fragmentLength(node.ident);
    int32 pos = (node.position.end > node.position.bgn) ? (node.position.end - node.position.bgn) : (node.position.bgn - node.position.end);

    if (node.contained)
      fprintf(logFile, "Added frag %d (len %d) to unitig %d at %d,%d (idx %d) (lendiff %d) (contained in %d)\n",
              node.ident, len, _id, node.position.bgn, node.position.end,
              ufpath.size() - 1,
              pos - len,
              node.contained);
    else
      fprintf(logFile, "Added frag %d (len %d) to unitig %d at %d,%d (idx %d) (lendiff %d)\n",
              node.ident, len, _id, node.position.bgn, node.position.end,
              ufpath.size() - 1,
              pos - len);
  }

  assert(node.position.bgn >= 0);
  assert(node.position.end >= 0);
}


//  This will add a contained fragment to a unitig, adjusting the position as needed.  It is only
//  needed when moving a contained read from unitig A to unitig B.  It is NOT needed when rebuilding
//  a unitig.
//
bool
Unitig::addContainedFrag(int32 fid, BestContainment *bestcont, bool report) {
  ufNode  frag;
  ufNode *parent = NULL;

  frag.ident        = fid;
  frag.contained    = bestcont->container;
  frag.parent       = bestcont->container;
  frag.ahang        = 0;
  frag.bhang        = 0;
  frag.position.bgn = 0;
  frag.position.end = 0;

  parent = &ufpath[pathPosition(frag.contained)];

#if 0
  //  This block is useful for debugging (maybe).  It is usually triggered only during popBubbles(),
  //  when we try to place a contained fragment into a fragment that has not been moved into the new
  //  unitig yet.  It might be useful if pathPosition ever gets messed up.
  //
  if ((parent == NULL) || (parent->ident != frag.contained)) {
    fprintf(logFile, "WARNING:  Didn't find the correct parent frag (%d) for contained frag %d.\n",
            frag.contained, frag.ident);
    fprintf(logFile, "          Found frag %d instead.\n", (parent == NULL) ? -1 : parent->ident);

    parent = NULL;

    for (int fi=0; fi<ufpath.size(); fi++) {
      ufNode *ix = &ufpath[fi];

      fprintf(logFile, "          path[%4d,%4d] is frag %d %s\n",
              fi, pathPosition(ix->ident),
              ix->ident,
              (ix->ident == frag.contained) ? " CORRECT PARENT!" : "");

      if (ix->ident == frag.contained)
        parent = ix;
    }
  }
#endif

  if ((parent == NULL) || (parent->ident != frag.contained)) {
    fprintf(logFile, "Unitig::addContainedFrag()-- WARNING:  Failed to place frag %d into unitig %d; parent not here.\n",
            fid, id());
    return(false);
  }
  
  //  Adjust orientation.  See comments in AS_BOG_Unitig.cc::placeContains().
  //
  //  NOTE!  Code is duplicated there.
  //
  //  NOTE!  The hangs are from the (parent) container to the (child) containee.  This is opposite
  //  as to how dovetail edges are stored.

  assert(bestcont->a_hang >= 0);
  assert(bestcont->b_hang <= 0);

  if (parent->position.bgn < parent->position.end) {
    //  Container is forward.
    frag.ahang = bestcont->a_hang;
    frag.bhang = bestcont->b_hang;

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
    frag.ahang = -bestcont->b_hang;
    frag.bhang = -bestcont->a_hang;

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


  //  Containments are particularily painful.  A beautiful example: a fragment of length 253bp is
  //  contained in a fragment of length 251bp (both hangs are zero).  In this case, the
  //  "ahang+length" method fails, placing the contained fragment outside the container (and if
  //  backwards oriented, _BEFORE_ the contained fragment).  The "ahang,bhang" method works here,
  //  but fails on other instances, shrinking deep containments to nothing.
  //
  //  We can use either method first, then adjust using the other method.
  //
  //  We'll use 'ahang,bhang' first (mostly because it was already done, and we need to compute
  //  those values anyway) then reset the end based on the length, limited to maintain a
  //  containment relationship.
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

  addFrag(frag, 0, report);

#if 0
  //  Bump that new fragment up to be in the correct spot -- we can't
  //  use the sort() method on Unitig, since we lost the
  //  containPartialOrder.
  //
  int             i = ufpath.size() - 1;
  ufNode   *f = &ufpath.front();

  //  Only needed if the frag we just added (i) begins before the second to last frag.

  if (MIN(f[i].position.bgn, f[i].position.end) < MIN(f[i-1].position.bgn, f[i-1].position.end)) {
    ufNode          containee    = f[i];
    int             containeeMin = MIN(containee.position.bgn, containee.position.end);

    while ((i > 0) &&
           (containee.contained != f[i-1].ident) &&
           (containeeMin < MIN(f[i-1].position.bgn, f[i-1].position.end))) {
      f[i] = f[i-1];
      i--;
    }

    f[i] = containee;
  }
#endif

  return(true);
}


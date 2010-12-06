
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

static const char *rcsid = "$Id: AS_BAT_Unitig_AddFrag.C,v 1.2 2010-12-06 08:03:48 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"




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
      fprintf(logFile, "Added frag %d (len %d) to unitig %d at %d,%d (idx %lu) (lendiff %d) (contained in %d)\n",
              node.ident, len, _id, node.position.bgn, node.position.end,
              ufpath.size() - 1,
              pos - len,
              node.contained);
    else
      fprintf(logFile, "Added frag %d (len %d) to unitig %d at %d,%d (idx %lu) (lendiff %d)\n",
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

  assert(bestcont->isContained);

  frag.ident        = fid;

  if (placeFrag(frag, bestcont) == false) {
    fprintf(logFile, "addContainedFrag()-- Failed to place contained frag %d using bestcont %d (hang %d,%d same orient %d).\n",
            fid, bestcont->container, bestcont->a_hang, bestcont->b_hang, bestcont->sameOrientation);
    return(false);
  }

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


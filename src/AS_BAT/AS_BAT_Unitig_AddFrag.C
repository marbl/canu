
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

static const char *rcsid = "$Id: AS_BAT_Unitig_AddFrag.C,v 1.4 2012-07-30 01:21:01 brianwalenz Exp $";

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
    int32 trulen = FI->fragmentLength(node.ident);
    int32 poslen = (node.position.end > node.position.bgn) ? (node.position.end - node.position.bgn) : (node.position.bgn - node.position.end);

    if (node.contained)
      writeLog("Added frag %d (len %d) to unitig %d at %d,%d (idx %lu) (lendiff %d) (contained in %d)\n",
              node.ident, trulen, _id, node.position.bgn, node.position.end,
              ufpath.size() - 1,
              poslen - trulen,
              node.contained);
    else
      writeLog("Added frag %d (len %d) to unitig %d at %d,%d (idx %lu) (lendiff %d)\n",
              node.ident, trulen, _id, node.position.bgn, node.position.end,
              ufpath.size() - 1,
              poslen - trulen);

    assert(poslen / trulen < 10);
    assert(trulen / poslen < 10);
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
    writeLog("addContainedFrag()-- Failed to place contained frag %d using bestcont %d (hang %d,%d same orient %d).\n",
            fid, bestcont->container, bestcont->a_hang, bestcont->b_hang, bestcont->sameOrientation);
    return(false);
  }

  addFrag(frag, 0, report);

  return(true);
}



//  Percolate the last fragment to the correct spot in the list.
void
Unitig::bubbleSortLastFrag(void) {
  uint32   previd  = ufpath.size() - 2;
  uint32   lastid  = ufpath.size() - 1;

  ufNode   last    = ufpath[lastid];
  uint32   lastbgn = MIN(last.position.bgn, last.position.end);

  while ((lastid > 0) &&
         (lastbgn < MIN(ufpath[previd].position.bgn, ufpath[previd].position.end))) {
    ufpath[lastid] = ufpath[previd];

    _pathPosition[ufpath[lastid].ident] = lastid;

    lastid--;
    previd--;
  }

  _pathPosition[last.ident] = lastid;

  if (lastid < ufpath.size() - 1)
    ufpath[lastid] = last;
}

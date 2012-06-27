
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id: classifyMates-DFS.C,v 1.10 2012-06-27 20:11:51 brianwalenz Exp $";

#include "AS_global.h"

#include "classifyMates.H"
#include "classifyMates-globalData.H"

#include <set>
using namespace std;


//  Attempt to find a path from the 5' end of this fragment to the 5' end of the next iid
//         <---- -path- ---->
//  If we find a path, declare this a PE and not a MP.
//
//  The search is depth first, stopping when we find a path, or when the path gets implausibly long.
//
void
cmGlobalData::doSearchDFS(cmComputation *c,
                          cmThreadData  *t) {
  set<uint32>  visited5p3;
  set<uint32>  visited3p5;

  //  +3 because:
  //  [0] is unused -- it's our termination condition
  //  [1] is the start node
  //  [2] is the first search node (depth 1)
  //  [3] is the next search node (the ones we test)
  //
  if (t->path == NULL) {
    t->pathMax = depthMax + 3;
    t->path    = new searchNode [t->pathMax];
  }

  t->clear();

  t->pathPos = 1;  //  CRITICAL; end of loop we increment t->pathOPos[t->pathPos-1]

  t->path[t->pathPos].pIID = c->fragIID;
  t->path[t->pathPos].p5p3 = c->frag5p3;
  t->path[t->pathPos].pLen = fi[c->fragIID].clearLength;
  t->path[t->pathPos].pDpt = 1;
  t->path[t->pathPos].oMax = bbLen[c->fragIID];
  t->path[t->pathPos].oPos = 0;
  t->path[t->pathPos].oLst = bbPos[c->fragIID];

  //  While we still have paths to follow
  //
  while (t->pathPos > 0) {

    //  Over all overlaps at this depth
    //
    for (;
         t->path[t->pathPos].oPos < t->path[t->pathPos].oMax;
         t->path[t->pathPos].oPos++) {
      overlapInfo  *novl = t->path[t->pathPos].oLst + t->path[t->pathPos].oPos;
      uint32        niid = novl->iid;
      bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
      int32         nlen = 0;

      //  Not true.  novl is a BB overlap, so the iid is the next fragment we'd move to.
      //  pIID is the BB fragment we are at currently.
      //assert(niid == t->path[t->pathPos].pIID);

      set<uint32> &visited = (t->path[t->pathPos].p5p3 == true) ? visited5p3 : visited3p5;
      if (visited.find(niid) != visited.end())
        //  Been here already.
        continue;

      if (fi[niid].isBackbone == false)
        //  Not a backbone read.
        continue;

      computeNextPlacement(t, novl, niid, n5p3, nlen);

      if (nlen < t->path[t->pathPos].pLen)
        //  Went backwards!
        continue;

      //  Went forwards.  Save the extension, unless that would put us over the depth limit.

      t->pathPos++;
      t->searchIter++;

      t->path[t->pathPos].pIID = niid;
      t->path[t->pathPos].p5p3 = n5p3;
      t->path[t->pathPos].pLen = nlen;
      t->path[t->pathPos].pDpt = t->path[t->pathPos-1].pDpt + 1;
      t->path[t->pathPos].oMax = bbLen[niid];
      t->path[t->pathPos].oPos = 0;
      t->path[t->pathPos].oLst = bbPos[niid];

      if (testSearch(c, t, tgPos, tgLen))
        //  If any of the target overlaps are the answer
        return;

      if (t->pathPos >= depthMax)
        //  End of the line.  Do not pass go.  Proceed directly to, ummm, the next overlap.
        t->pathPos--;
    }

    //  We've exhausted the paths for this fragment.  Mark it finished and back up one.

    set<uint32> &visited = (t->path[t->pathPos].p5p3 == true) ? visited5p3 : visited3p5;
    visited.insert(t->path[t->pathPos].pIID);

    t->pathPos--;
  }

  //  Not found.
  assert(c->result.classified == false);
}


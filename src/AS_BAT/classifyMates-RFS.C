
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

static const char *rcsid = "$Id: classifyMates-RFS.C,v 1.8 2011-08-22 16:44:19 brianwalenz Exp $";

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
cmGlobalData::doSearchRFS(cmComputation *c,
                          cmThreadData  *t) {

  t->pathPos    = 0;
  t->pathAdd    = 0;
  t->searchIter = 0;

  if (t->path == NULL) {
    t->pathMax = distMax;  //  Can never have more than one overlap per base if distance
    t->path    = new searchNode [t->pathMax];
  }

  for (uint32 iter=0; iter<pathsMax; iter++) {
    t->pathPos= 0;
    t->searchIter++;

    t->path[t->pathPos].pIID = c->fragIID;
    t->path[t->pathPos].p5p3 = c->frag5p3;
    t->path[t->pathPos].pLen = fi[c->fragIID].clearLength;
    t->path[t->pathPos].oMax = bbLen[c->fragIID];
    t->path[t->pathPos].oPos = 0;
    t->path[t->pathPos].oLst = bbPos[c->fragIID];

    //  Follow random paths until we get too long or too deep.  If we find the answer we immediately
    //  return.  If we don't find the answer we exit the while and do another iteration.

    while (t->path[t->pathPos].pLen < distMax) {
      if (testSearch(c, t, tgPos, tgLen))
        //  If any of the target overlaps are the answer
        return;

      //  Compute which edges extend in the correct direction.
      //
      t->extLen = 0;
      for (uint32 test=0; test<bbLen[t->path[t->pathPos].pIID]; test++) {
        overlapInfo  *novl = t->path[t->pathPos].oLst + test;
        uint32        niid = novl->iid;
        bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
        uint32        nlen = 0;

        if (fi[niid].isBackbone == false)
          //  Not a backbone read
          continue;

        computeNextPlacement(c, t, novl, niid, n5p3, nlen);

        if (nlen > t->path[t->pathPos].pLen)
          t->ext[t->extLen++] = test;
      }

      if (t->extLen == 0)
        //  If there are no backbone overlaps out of here
        break;

      t->path[t->pathPos].oPos = t->ext[lrand48() % t->extLen];

      overlapInfo  *novl = t->path[t->pathPos].oLst + t->path[t->pathPos].oPos;
      uint32        niid = novl->iid;
      bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
      uint32        nlen = 0;

      if (fi[niid].isBackbone == false)
        //  Not a backbone read
        continue;

      computeNextPlacement(c, t, novl, niid, n5p3, nlen);

      //  We're guaranteed to always advance the path (unlike DFS) so do it.

      if (t->pathPos < t->pathMax) {
        t->pathPos++;

        t->path[t->pathPos].pIID = niid;
        t->path[t->pathPos].p5p3 = n5p3;
        t->path[t->pathPos].pLen = nlen;
        t->path[t->pathPos].oMax = bbLen[niid];
        t->path[t->pathPos].oPos = 0;
        t->path[t->pathPos].oLst = bbPos[niid];
      }

      assert(t->path[t->pathPos].pLen > t->path[t->pathPos-1].pLen);
    }

  }  //  Try a bunch of random stabs to find the path

  //  Not found.
  assert(c->result.classified == false);
}

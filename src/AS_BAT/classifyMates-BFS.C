
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

static const char *rcsid = "$Id: classifyMates-BFS.C,v 1.13 2012-07-05 04:07:51 brianwalenz Exp $";

#include "AS_global.h"

#include "classifyMates.H"
#include "classifyMates-globalData.H"

#include <set>
using namespace std;


void
cmGlobalData::doSearchBFS(cmComputation *c,
                          cmThreadData  *t) {

  if (t->path == NULL) {
    t->pathMax        = nodesMax;
    t->path           = new searchNode [t->pathMax];

    t->solutionSet    = new bits(numFrags + 1, "solutionSet");

    t->visited5p3bits = new bits(numFrags + 1, "visited5p3");
    t->visited3p5bits = new bits(numFrags + 1, "visited3p5");

    t->visitedListMax = 16 * 1024 * 1024;
    t->visitedList    = new uint32 [t->visitedListMax];
  }

  t->clear();

  //  Build the map from backbone fragment to solution overlap.  The map goes from a-fragID to
  //  overlap with the mateIID b-frag.  In other words, if this is set, the backbone fragment has an
  //  overlap to the mate fragment we are searching for.

  t->solutionSet->testClear();

  for (uint32 ii=0; ii<gtLen[c->mateIID]; ii++)
    t->solutionSet->set(gtPos[c->mateIID][ii].iid);

  //  Seed the search with the first fragment

  t->path[t->pathAdd].pIID = c->fragIID;
  t->path[t->pathAdd].p5p3 = c->frag5p3;
  t->path[t->pathAdd].pLen = fi[c->fragIID].clearLength;
  t->path[t->pathAdd].pDpt = 1;
  t->path[t->pathAdd].oMax = 0;
  t->path[t->pathAdd].oPos = 0;
  t->path[t->pathAdd].oLst = bbPos[c->fragIID];

  t->pathAdd++;

  uint32  nBB=0, nBACK=0, nFAR=0, nHERE=0;
  uint32  maxLen=0;

  for (;
       ((t->pathPos < t->pathMax) &&
        (t->pathPos < t->pathAdd));
       t->pathPos++) {

    if (maxLen < t->path[t->pathPos].pLen)
      maxLen = t->path[t->pathPos].pLen;

    if ((distMin                  <= t->path[t->pathPos].pLen) &&
        (t->path[t->pathPos].pLen <= distMax) &&
        (testSearch(c, t)))
      //  If any of the target overlaps are the answer
      goto returnBFS;

    if (t->pathAdd >= t->pathMax)
      //  No space for more overlaps, abort adding any.
      continue;

    //  Add more fragments to the search.

    for (uint32 o=0;
         ((t->pathAdd < t->pathMax) &&
          (o < bbLen[t->path[t->pathPos].pIID])); o++) {
      overlapInfo  *novl = t->path[t->pathPos].oLst + o;
      uint32        niid = novl->iid;
      bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
      int32         nlen = 0;

      if (fi[niid].isBackbone == false) {
        //  Not a backbone read
        nBB++;
        continue;
      }

      computeNextPlacement(t, novl, niid, n5p3, nlen);

      if (nlen <= t->path[t->pathPos].pLen) {
        //  Path went backwards.
        nBACK++;
        continue;
      }

      if (nlen > distMax) {
        //  Path too far, don't add
        nFAR++;
        continue;
      }

      bits  *visited = (t->path[t->pathPos].p5p3 == true) ? t->visited5p3bits : t->visited3p5bits;

      if (visited->isSet(niid) == true) {
        //  Been here already.
        nHERE++;
        continue;
      }

      visited->set(niid);

      t->visitedList[t->visitedListLen++] = niid;
      assert(t->visitedListLen < t->visitedListMax);

      t->path[t->pathAdd].pIID = niid;
      t->path[t->pathAdd].p5p3 = n5p3;
      t->path[t->pathAdd].pLen = nlen;
      t->path[t->pathAdd].pDpt = t->path[t->pathPos].pDpt + 1;
      t->path[t->pathAdd].oMax = 0;
      t->path[t->pathAdd].oPos = 0;
      t->path[t->pathAdd].oLst = bbPos[niid];

      t->pathAdd++;
      t->searchIter++;
    }
  }

  if (t->pathAdd >= t->pathMax)
    c->result.limited = true;
  else
    c->result.exhausted = true;

  //  Not found.
  assert(c->result.classified == false);

 returnBFS:

  //fprintf(stderr, "SEARCH:  maxLen %u visited %u BB %u BACK %u FAR %u HERE %u\n",
  //        maxLen, t->pathAdd, nBB, nBACK, nFAR, nHERE);

  for (uint32 ii=0; ii<t->visitedListLen; ii++) {
    t->visited5p3bits->clear(t->visitedList[ii]);
    t->visited3p5bits->clear(t->visitedList[ii]);
  }

  for (uint32 ii=0; ii<gtLen[c->mateIID]; ii++)
    t->solutionSet->clear( gtPos[c->mateIID][ii].iid );

  t->solutionSet->testClear();
}


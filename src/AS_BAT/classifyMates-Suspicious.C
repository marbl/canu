
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

static const char *rcsid = "$Id: classifyMates-Suspicious.C,v 1.3 2012-07-09 06:04:33 brianwalenz Exp $";

#include "AS_global.h"

#include "classifyMates.H"
#include "classifyMates-globalData.H"

#include <set>
using namespace std;

//  Search for suspicious mate overlaps.  This modified from BFS and RFS.
//
//  Suspicious is defined as ANY short path in ANY orientation from the frag to the mate.


//  For BFS, the max distance we we'll search for, and the path depth.  EVERYTHING found below this
//  size is bad.  Including any INNIE-MP found.  DISTMAX must be smaller than the any valid MP.
//
#define DISTMAX  500
#define DEPTHMAX 10

//  The number of iterations of random-path search.  INNIE-MP found here are not bad.
//
#define RFSITERS    100
#define RFSDISTMAX  5000


void
cmGlobalData::testSuspicious(uint32           fragIID,
                             bool             frag5p3,
                             uint32           mateIID,
                             cmThreadData    *t,
                             vector<int32>   &dist3p5,
                             vector<int32>   &dist5p3,
                             bool             onlyNormal) {

  uint32  thisIID = t->path[t->pathPos].pIID;

  //  If not set, then this fragment does not have an overlap to our mate fragment.

  if (t->solutionSet->isSet(thisIID) == false)
    return;

  //  Now search for the overlap to the mate fragment.  When we built the gt overlaps,
  //  we used the source iid, not the mate iid -- we just search for the overlap that
  //  has the same iid as our source fragment.

  overlapInfo  *pos = gtPos[mateIID];
  uint32        len = gtLen[mateIID];
  uint32        iii = 0;

  while ((iii<len) && (thisIID != pos[iii].iid))
    iii++;

  assert(iii < len);

  overlapInfo  *novl = pos + iii;
  uint32        niid = mateIID;
  bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
  int32         nlen = 0;

  computeNextPlacement(t, novl, niid, n5p3, nlen);

  //  Save the orientation and distance to the mate.

  //  We get negative distances here in one particularily nasty case.
  //
  //  bb read ------------------------------------>
  //                                   --A-->
  //                <--B--
  //
  //  The search starts at A, moving onto the BB read.  The BB read has a terminating overlap to the
  //  B read, but it is at a negative distance.  This is nasty because even if we start from the B
  //  read, the same thing happens.  The only fix is to search from both ends of the A read (we do
  //  that; see cmWorker()).  When searching off the beginning of A, we add in the large 5' end of
  //  the BB read, then subtract out the small hang to get back to the B read.
  //
  //  So, long story short, just ignore negative stuff here.

  if (nlen < 0)
    return;

  //fprintf(stderr, "FOUND at len %d depth %u\n", nlen, t->path[t->pathPos].pDpt);

  if ((onlyNormal) && (frag5p3 != n5p3))
    return;

  if (n5p3)
    dist5p3.push_back(nlen);
  else
    dist3p5.push_back(nlen);
}



void
cmGlobalData::doSearchSuspicious(uint32          fragIID,
                                 bool            frag5p3,
                                 uint32          mateIID,
                                 cmThreadData   *t,
                                 vector<int32>  &dist3p5close,
                                 vector<int32>  &dist3p5far,
                                 vector<int32>  &dist5p3close,
                                 vector<int32>  &dist5p3far) {

  if (t->path == NULL) {
    t->pathMax        = nodesMax;
    t->path           = new searchNode [t->pathMax];

    t->solutionSet    = new bits(numFrags + 1, "solutionSet");

    t->visited5p3bits = new bits(numFrags + 1, "visited5p3bits");
    t->visited3p5bits = new bits(numFrags + 1, "visited3p5bits");

    t->visitedListMax = 16 * 1024 * 1024;
    t->visitedList    = new uint32 [t->visitedListMax];
  }

  dist5p3close.clear();
  dist5p3far.clear();
  dist3p5close.clear();
  dist3p5far.clear();

  t->clear();

  //  Build the map from backbone fragment to solution overlap.  The map goes from a-fragID to
  //  overlap with the mateIID b-frag.  In other words, if this is set, the backbone fragment has an
  //  overlap to the mate fragment we are searching for.

  t->solutionSet->testClear();

  for (uint32 ii=0; ii<gtLen[mateIID]; ii++)
    t->solutionSet->set(gtPos[mateIID][ii].iid);

  ////////////////////////////////////////
  //
  //  Pass 1: Breadth first search for short insert crud.
  //
  t->path[t->pathAdd].pIID = fragIID;
  t->path[t->pathAdd].p5p3 = frag5p3;
  t->path[t->pathAdd].pLen = fi[fragIID].clearLength;
  t->path[t->pathAdd].pDpt = 1;
  t->path[t->pathAdd].oMax = 0;
  t->path[t->pathAdd].oPos = 0;
  t->path[t->pathAdd].oLst = bbPos[fragIID];

  t->pathAdd++;

  for (;
       ((t->pathPos < t->pathMax) &&
        (t->pathPos < t->pathAdd));
       t->pathPos++) {

    testSuspicious(fragIID, frag5p3, mateIID, t, dist3p5close, dist5p3close, false);

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

      if (fi[niid].isBackbone == false)
        //  Not a backbone read
        continue;

      computeNextPlacement(t, novl, niid, n5p3, nlen);

      if (nlen <= t->path[t->pathPos].pLen)
        //  Path went backwards.
        continue;

      if (nlen > DISTMAX)
        //  Path too far, don't add
        continue;

      if (t->path[t->pathPos].pDpt + 1 > DEPTHMAX)
        //  Path too deep, don't add
        continue;

      bits  *visited = (t->path[t->pathPos].p5p3 == true) ? t->visited5p3bits : t->visited3p5bits;

      if (visited->isSet(niid) == true)
        //  Been here already.
        continue;

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
    }
  }

  //  Clean up the markings.

  for (uint32 ii=0; ii<t->visitedListLen; ii++) {
    t->visited5p3bits->clear(t->visitedList[ii]);
    t->visited3p5bits->clear(t->visitedList[ii]);
  }

  ////////////////////////////////////////
  //
  //  Pass 2: Random-path search for larger insert normal-oriented crud.
  //
  //  Benchmark - with this enabled at 100 iters, 492/sec.  With it disabled, 673/sec.

  for (uint32 iter=0; iter<RFSITERS; iter++) {
    t->clear();

    t->path[t->pathPos].pIID = fragIID;
    t->path[t->pathPos].p5p3 = frag5p3;
    t->path[t->pathPos].pLen = fi[fragIID].clearLength;
    t->path[t->pathPos].pDpt = 1;
    t->path[t->pathPos].oMax = bbLen[fragIID];
    t->path[t->pathPos].oPos = 0;
    t->path[t->pathPos].oLst = bbPos[fragIID];

    //  Don't know how far to search.  The one sample I have seen, supposedly a 3000 bp library,
    //  tailed off at 3500 bp.

    while (t->path[t->pathPos].pLen < RFSDISTMAX) {

      testSuspicious(fragIID, frag5p3, mateIID, t, dist3p5far, dist5p3far, true);

      //  Compute which edges extend in the correct direction.
      //
      t->extLen = 0;

      for (uint32 test=0; test<bbLen[t->path[t->pathPos].pIID]; test++) {
        overlapInfo  *novl = t->path[t->pathPos].oLst + test;
        uint32        niid = novl->iid;
        bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
        int32         nlen = 0;

        if (fi[niid].isBackbone == false)
          //  Not a backbone read
          continue;

        computeNextPlacement(t, novl, niid, n5p3, nlen);

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
      int32         nlen = 0;

      if (fi[niid].isBackbone == false)
        //  Not a backbone read
        continue;

      computeNextPlacement(t, novl, niid, n5p3, nlen);

      assert(t->path[t->pathPos].pLen < nlen);

      //  We're guaranteed to always advance the path (unlike DFS) so do it.

      if (t->pathPos < t->pathMax) {
        t->pathPos++;

        t->path[t->pathPos].pIID = niid;
        t->path[t->pathPos].p5p3 = n5p3;
        t->path[t->pathPos].pLen = nlen;
        t->path[t->pathPos].pDpt = t->path[t->pathPos-1].pDpt + 1;
        t->path[t->pathPos].oMax = bbLen[niid];
        t->path[t->pathPos].oPos = 0;
        t->path[t->pathPos].oLst = bbPos[niid];
      }

      assert(t->path[t->pathPos].pLen > t->path[t->pathPos-1].pLen);
    }  //  End of extending the path
  }  //  End of random path loop.

  //  Clean up our markings.

  for (uint32 ii=0; ii<gtLen[mateIID]; ii++)
    t->solutionSet->clear(gtPos[mateIID][ii].iid);

  t->solutionSet->testClear();
}

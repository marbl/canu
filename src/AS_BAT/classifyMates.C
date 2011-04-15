
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

const char *mainid = "$Id: classifyMates.C,v 1.4 2011-04-15 05:02:29 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"

#include <set>

#include "sweatShop.H"

using namespace std;

bool   VERBOSE1 = false;
bool   VERBOSE2 = false;
bool   VERBOSE3 = false;
bool   VERBOSE4 = false;
bool   VERBOSE5 = false;

#define NDISTD  6  //  Number of dovetail edges off each end
#define NDISTC  5  //  Number of contains for each end


class fragmentInfo {
public:
  fragmentInfo() {
  };
  ~fragmentInfo() {
  };

  uint64  unused      : 15;
  uint64  isBackbone  : 1;
  uint64  doSearch    : 1;
  uint64  clearLength : 15;
  uint64  mateIID     : 32;
};


class overlapInfo {
public:
  overlapInfo() {
  };
  overlapInfo(OVSoverlap &ovl) {
    unused  = 0;
    flipped = ovl.dat.ovl.flipped;
    ahang   = ovl.dat.ovl.a_hang;
    bhang   = ovl.dat.ovl.b_hang;
    iid     = ovl.b_iid;
  };
  ~overlapInfo() {
  };

  uint32  unused  : 1;
  uint32  flipped : 1;
  int32   ahang   : 15;
  int32   bhang   : 15;

  uint32 iid;
};



class
saveDistance {
public:
  saveDistance() {
  };

  ~saveDistance() {
  };

  //  We save overlaps based on the size of the overlap.
  //
  //  Keep the best and near best overlaps for dovetail overlaps.
  //
  //    ------------------------
  //       ----------------------------- (best)
  //          ---------------------- (near best)
  //                 --------------
  //                     ----------------------------------
  //                          -----------------------------------
  //
  //  Keep containment overlaps that are near the end.
  //
  //                     ---------
  //    ----------------------------- (near the end)
  //        -------------------------------
  //            ----------------------------------
  //                 ------------------------------------ (near the end)
  //                     --------------------------------------- (near the end)
  //
  //  doveDist - Dovetail overlaps with overlap length at least this big are saved.
  //             The 'near best' above might not be informative, but it is still
  //             the second best overlap and is kept.
  //
  //  coneDist - Containee overlaps (A contained in B) are saved if they are
  //             at least this close to the end of the other fragment.
  //
  //  conrDist - Container overlaps (A contains B), likewise.

  int32   doveDist5arr[NDISTD];  //  scratch array for finding the nth largest distance
  int32   coneDist5arr[NDISTC];
  int32   conrDist5arr[NDISTC];

  int32   doveDist3arr[NDISTD];
  int32   coneDist3arr[NDISTC];
  int32   conrDist3arr[NDISTC];

  int32   doveDist5;     //  minimum distance we should be saving
  int32   coneDist5;
  int32   conrDist5;

  int32   doveDist3;
  int32   coneDist3;
  int32   conrDist3;



  //  Save the N largest values - sorted largest to smallest
  void    saveDistMax(int32 *darr, int32  dist) {

    assert(dist >= 0);

    if (dist < darr[NDISTD-1])
      //  Distance less than the smallest distance we want to keep, don't save
      return;

    //  We either want to save a new distance, pushing the last one off of the array,
    //  or notice that we've already saved this distance and leave the array alone.

    for (int32 i=0; i<NDISTD; i++) {
      if (darr[i] == dist)
        //  Saved it already.
        return;

      if (darr[i] < dist) {
        //  Save at i, push i and following down one slot.
        for (int32 j=NDISTD-1; j>i; j--)
          darr[j] = darr[j-1];
        darr[i] = dist;
        return;
      }
    }

    //  Fail, we should never get here.
    assert(0);
  };



  //  Save the N smallest values - sorted smallest to largest
  void    saveDistMin(int32 *darr, int32  dist) {

    assert(dist >= 0);

    if (dist > darr[NDISTC-1])
      //  Distance more than the biggest distance we want to keep, don't save
      return;

    //  We either want to save a new distance, pushing the last one off of the array,
    //  or notice that we've already saved this distance and leave the array alone.

    for (int32 i=0; i<NDISTC; i++) {
      if (darr[i] == dist)
        //  Saved it already.
        return;

      if (darr[i] > dist) {
        //  Save at i, push i and following down one slot.
        for (int32 j=NDISTC-1; j>i; j--)
          darr[j] = darr[j-1];
        darr[i] = dist;
        return;
      }
    }

    //  Fail, we should never get here.
    assert(0);
  };



  void   compute(fragmentInfo     *fi,
                 OVSoverlap       *ovl,
                 uint32            ovlLen) {

    if (ovlLen == 0)
      return;

    for (int32 i=0; i<NDISTD; i++)
      doveDist5arr[i] = doveDist3arr[i] = INT32_MIN;
    for (int32 i=0; i<NDISTC; i++) {
      coneDist5arr[i] = coneDist3arr[i] = INT32_MAX;
      conrDist5arr[i] = conrDist3arr[i] = INT32_MAX;
    }

    for (uint32 i=0; i<ovlLen; i++) {
      int32  ah = ovl[i].dat.ovl.a_hang;
      int32  bh = ovl[i].dat.ovl.b_hang;
      int32  fa = fi[ovl[i].a_iid].clearLength;
      int32  fb = fi[ovl[i].b_iid].clearLength;

      if        (AS_OVS_overlapAEndIs5prime(ovl[i])) {
        //  ah < 0 && bh < 0
        saveDistMax(doveDist5arr, fb - -ah);

      } else if (AS_OVS_overlapAEndIs3prime(ovl[i])) {
        //  ah > 0 && bh > 0
        saveDistMax(doveDist3arr, fb -  bh);

      } else if (AS_OVS_overlapAIsContained(ovl[i])) {
        //  ah <= 0 && bh >= 0
        saveDistMin(coneDist5arr, -ah);
        saveDistMin(coneDist3arr,  bh);
        
      } else if (AS_OVS_overlapAIsContainer(ovl[i])) {
        //  ah >= 0 && bh <= 0
        saveDistMin(conrDist5arr,  ah);
        saveDistMin(conrDist3arr, -bh);

      } else {
        assert(0);
      }
    }

    doveDist5 = doveDist5arr[NDISTD-1];
    coneDist5 = coneDist5arr[NDISTC-1];
    conrDist5 = conrDist5arr[NDISTC-1];

    doveDist3 = doveDist3arr[NDISTD-1];
    coneDist3 = coneDist3arr[NDISTC-1];
    conrDist3 = conrDist3arr[NDISTC-1];

#if 0
    fprintf(stderr, "IID %8d  dove %3d,%3d  cone %3d,%3d  conr %3d,%3d\n",
            ovl[0].a_iid,
            doveDist5, doveDist3,
            coneDist5, coneDist3,
            conrDist5, conrDist3);
#endif
  };
};









class cmThreadData {
public:
  cmThreadData() {
  };
  ~cmThreadData() {
  };

  //  the arrays in cmComputation can be moved here
};


class cmComputation {
public:
  cmComputation(uint32 iid_) {
    iid       = iid_;

    pathDepth = 1;
    pathFound = false;

    memset(pathIID,  1024 * sizeof(uint32), 0);
    memset(path5p3,  1024 * sizeof(uint32), 0);

    memset(pathLen,  1024 * sizeof(uint32), 0);

    memset(pathRoot, 1024 * sizeof(overlapInfo *), 0);
    memset(pathPosn, 1024 * sizeof(uint32), 0);
    memset(pathMaxp, 1024 * sizeof(uint32), 0);
  };
  ~cmComputation() {
  };

public:
  uint32         iid;

  uint32         pathDepth;
  bool           pathFound;

  uint32         pathIID[1024];   //  Fragment here
  uint32         path5p3[1024];   //  Fragment here is 5' to 3'

  uint32         pathLen[1024];   //  The length, in bp, of the path up till now

  overlapInfo   *pathRoot[1024];  //  Root of the list of overlaps for this iid
  uint32         pathPosn[1024];  //  Position we are at in the list of overlaps
  uint32         pathMaxp[1024];  //  Number of ovelerlaps for this iid
};


class cmGlobalData {
public:
  cmGlobalData(FILE    *resultsFile_,
               uint32   pathMin_,
               uint32   pathMax_,
               bool     pathInnie_,
               uint32   depthMax_) {

    resultsFile  = resultsFile_;

    pathMin      = pathMin_;
    pathMax      = pathMax_;
    pathInnie    = pathInnie_;
    depthMax     = depthMax_;

    curFrag      = 0;
    numFrags     = 0;
    fi           = 0L;

    oiPos        = 0L;
    oiLen        = 0L;

    oiStorageMax = 0;
    oiStorageLen = 0;
    oiStorageArr = 0L;
    oiStorage    = 0L;
  };

  ~cmGlobalData() {
    for (uint32 i=0; i<oiStorageMax; i++)
      delete [] oiStorageArr[i];

    delete [] oiPos;
    delete [] oiLen;

    delete [] oiStorageArr;

    delete [] fi;
  };

  void   loadFragments(char    *gkpStorePath,
                       bool     searchLibs,
                       uint32  *searchLib,
                       bool     backboneLibs,
                       uint32  *backboneLib);

  void   loadOverlaps(char   *ovlStorePath);

  void   doSearch(cmComputation  *c,
                  cmThreadData   *t);

  //  Overlap storage.  Each fragment has an entry in oiPos that points to actual data in an array.
  //  There are oiLen overlaps at this location.  The array is NOT contiguous.
  //
  //  oiStorage is a list of allocated memory.  We don't really need to keep the list around, but
  //  we'll play nice and clean up at the end.  Each chunk of allocated memory holds a fixed number
  //  of overlaps -- here, 128 million.  There are 1024 chunks, so we can store 128 billion
  //  overlaps, far far more than any overlap store we've created (about 3.3TB of data).

public:
  FILE             *resultsFile;

  uint32            pathMin;
  uint32            pathMax;
  bool              pathInnie;
  uint32            depthMax;

  uint32            curFrag;
  uint32            numFrags;
  fragmentInfo     *fi;

  overlapInfo     **oiPos;  //  Pointer to start of overlaps for this frag
  uint32           *oiLen;  //  Number of overlaps for this frag

  uint32            oiStorageMax;
  uint32            oiStorageLen;
  overlapInfo     **oiStorageArr;
  overlapInfo      *oiStorage;
};






void
cmGlobalData::loadFragments(char    *gkpStorePath,
                            bool     searchLibs,
                            uint32  *searchLib,
                            bool     backboneLibs,
                            uint32  *backboneLib) {

  fprintf(stderr, "LOADING FRAGMENTS...\n");

  gkStore          *gkpStore  = new gkStore(gkpStorePath, FALSE, FALSE);
  gkStream         *gkpStream = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);
  gkFragment        frg;

  numFrags = gkpStore->gkStore_getNumFragments();
  fi       = new fragmentInfo [numFrags + 1];

  memset(fi, 0, sizeof(fragmentInfo) * (numFrags + 1));

  //  CHANGE CLEAR RANGE FOR PRODUCTION!  We need some way of asking if
  //  a clear range exists.  For now, we'll use a not-so-good method.

  while (gkpStream->next(&frg)) {
    uint32 fid = frg.gkFragment_getReadIID();
    uint32 lib = frg.gkFragment_getLibraryIID();

    uint32  clobt = frg.gkFragment_getClearRegionEnd(AS_READ_CLEAR_OBTCHIMERA);
    uint32  cllat = frg.gkFragment_getClearRegionEnd(AS_READ_CLEAR_LATEST);
    uint32  cllen = 0;

    if ((clobt == 0) && (cllat != 0))
      //  OBT doesn't exist, but latest does.
      cllen = frg.gkFragment_getClearRegionLength(AS_READ_CLEAR_LATEST);
    else
      //  OBT exists.
      cllen = frg.gkFragment_getClearRegionLength(AS_READ_CLEAR_OBTCHIMERA);


    fi[fid].unused      = 0;
    fi[fid].doSearch    = 0;
    fi[fid].clearLength = cllen;
    fi[fid].mateIID     = frg.gkFragment_getMateIID();

    if (searchLibs == false)
      fi[fid].doSearch = gkpStore->gkStore_getLibrary(lib)->doMerBasedTrimming;
    else
      fi[fid].doSearch = searchLib[lib];

    if (backboneLibs == false)
      fi[fid].isBackbone = true;
    else
      fi[fid].isBackbone = backboneLib[lib];
  }

  delete gkpStream;
  delete gkpStore;

  fprintf(stderr, "LOADING FRAGMENTS...%u fragments loaded.\n", numFrags);
}


void
cmGlobalData::loadOverlaps(char  *ovlStorePath) {

  fprintf(stderr, "LOADING OVERLAPS...\n");

  OverlapStore     *ovlStore  = AS_OVS_openOverlapStore(ovlStorePath);

  uint32            ovlMax = 1048576;
  uint32            ovlLen = 1;
  OVSoverlap       *ovl    = new OVSoverlap [ovlMax];


  oiPos = new overlapInfo * [numFrags + 1];  //  Pointer to start of overlaps for this frag
  oiLen = new uint32        [numFrags + 1];  //  Number of overlaps for this frag

  memset(oiPos, 0, sizeof(overlapInfo) * (numFrags + 1));
  memset(oiLen, 0, sizeof(uint32)      * (numFrags + 1));

  oiStorageMax = 1024;
  oiStorageLen = 0;
  oiStorageArr = new overlapInfo * [oiStorageMax];
  oiStorage    = 0L;

  memset(oiStorageArr, 0, sizeof(overlapInfo *) * oiStorageMax);

  uint32            oiStorageBS  = 128 * 1024 * 1024;  //  Fixed block size
  uint32            oiStoragePos = 0;                  //  Position in the current block

  oiStorage = oiStorageArr[0] = new overlapInfo [oiStorageBS];

  uint64            numOvl  = 0;
  uint64            numFilt = 0;

  ovlLen = AS_OVS_readOverlapsFromStore(ovlStore, ovl, ovlMax, AS_OVS_TYPE_OVL);

  saveDistance     *dist = new saveDistance;

  while (ovlLen > 0) {
    //  Filter out overlaps we don't care about.

    uint32  iid = ovl[0].a_iid;

    oiPos[iid] = NULL;
    oiLen[iid] = 0;

    int32  doveDist5 = INT32_MAX, coneDist5 = INT32_MAX, conrDist5 = INT32_MAX;
    int32  doveDist3 = INT32_MAX, coneDist3 = INT32_MAX, conrDist3 = INT32_MAX;

    dist->compute(fi, ovl, ovlLen);

    for (uint32 i=0; i<ovlLen; i++) {
      bool  saveOverlap = false;

      int32  ah = ovl[i].dat.ovl.a_hang;
      int32  bh = ovl[i].dat.ovl.b_hang;
      int32  fa = fi[ovl[i].a_iid].clearLength;
      int32  fb = fi[ovl[i].b_iid].clearLength;

      if        (AS_OVS_overlapAEndIs5prime(ovl[i])) {
        //  ah < 0 && bh < 0
        if (dist->doveDist5 <= fb - -ah) {
          saveOverlap = true;
        }

      } else if (AS_OVS_overlapAEndIs3prime(ovl[i])) {
        //  ah > 0 && bh > 0
        if (dist->doveDist3 <= fb -  bh) {
          saveOverlap = true;
        }

      } else if (AS_OVS_overlapAIsContained(ovl[i])) {
        //  ah <= 0 && bh >= 0
        if (dist->coneDist5 >= -ah) {
          saveOverlap = true;
        }
        if (dist->coneDist3 >=  bh) {
          saveOverlap = true;
        }

      } else if (AS_OVS_overlapAIsContainer(ovl[i])) {
        //  ah >= 0 && bh <= 0
        if (dist->conrDist5 >=  ah) {
          saveOverlap = true;
        }
        if (dist->conrDist3 >= -bh) {
          saveOverlap = true;
        }

      } else {
        assert(0);
      }

      //  If the user specified libraries for the backbone, ignore this overlap if either fragment
      //  isn't in the backbone, and neither fragment is one we are searching.
      //  
      if (((fi[ovl[i].a_iid].isBackbone == false) || (fi[ovl[i].b_iid].isBackbone == false)) &&
          (fi[ovl[i].a_iid].doSearch == false) &&
          (fi[ovl[i].b_iid].doSearch == false))
        saveOverlap = false;

      if (saveOverlap)
        //  Save the overlap, count how many to save so we can get space.
        oiLen[iid]++;
      else
        //  Don't save overlap, remove the a_iid and we'll filter it out.
        ovl[i].a_iid = 0;
    }

    //fprintf(stderr, "Filtered IID %d from %d overlaps to %d overlaps.\n",
    //        iid, ovlLen, oiLen[iid]);

    //  If not enough space in oiStorage, get more space.

    if (oiStoragePos + oiLen[iid] > oiStorageBS) {
      oiStorage = oiStorageArr[oiStorageLen] = new overlapInfo [oiStorageBS];
      oiStorageLen++;
      oiStoragePos = 0;
    }

    //  Add the overlaps.

    oiPos[iid] = oiStorage + oiStoragePos;

    for (uint32 i=0; i<ovlLen; i++) {
      if (ovl[i].a_iid == iid) {
        oiStorage[oiStoragePos++] = overlapInfo(ovl[i]);
        numOvl++;
      } else {
        numFilt++;
      }
    }

    assert(oiStoragePos <= oiStorageBS);

    ovlLen = AS_OVS_readOverlapsFromStore(ovlStore, ovl, ovlMax, AS_OVS_TYPE_OVL);

    if ((iid % 100000) == 0)
      fprintf(stderr, "LOADING OVERLAPS...at IID %u (%06.2f%%): %lu (%06.2f%%) filtered, %lu (%06.2f%%) loaded.\r",
              iid,     100.0 * iid     / numFrags,
              numFilt, 100.0 * numFilt / (numFilt + numOvl),
              numOvl,  100.0 * numOvl  / (numFilt + numOvl));
  }

  AS_OVS_closeOverlapStore(ovlStore);

  delete [] ovl;

  fprintf(stderr, "\nLOADING OVERLAPS...%lu overlaps loaded.\n", numOvl);
}



//  Attempt to find a path from the 5' end of this fragment to the 5' end of the next iid
//         <---- -path- ---->
//  If we find a path, declare this a PE and not a MP.
//
//  The search is depth first, stopping when we find a path, or when the path gets implausibly long.
//
void
cmGlobalData::doSearch(cmComputation *c,
                       cmThreadData  *t) {
  bool         pathFound = false;
  set<uint32>  visited5p3;
  set<uint32>  visited3p5;

  //  If pathInnie == false, then we're attempting to find a path for outtie oriented fragments.
  //  In this case, we start with the first fragment 5p3=false, and need to end with 5p3=true.

  bool    bgn5p3 = (pathInnie == false) ? false : true;
  bool    end5p3 = (pathInnie == false) ? true : false;

  c->pathDepth                  = 1;

  c->pathIID[c->pathDepth]  = c->iid;
  c->path5p3[c->pathDepth]  = bgn5p3;
  c->pathLen[c->pathDepth]  = fi[c->iid].clearLength;
  c->pathRoot[c->pathDepth] = oiPos[c->iid];
  c->pathPosn[c->pathDepth] = 0;
  c->pathMaxp[c->pathDepth] = oiLen[c->iid];

  //if (VERBOSE1)
  //  fprintf(stderr, "FRAG %5u/%s with %5u overlaps (len %u).\n",
  //          iid, (c->path5p3[c->pathDepth] == true) ? "5'3'" : "3'5'", oiLen[iid], c->pathLen[c->pathDepth]);

  while (c->pathDepth > 0) {
    while (c->pathPosn[c->pathDepth] < c->pathMaxp[c->pathDepth]) {

      //
      //  If we've already seen this fragment in this orientation, get out of here.
      //

      set<uint32> &visited = (c->path5p3[c->pathDepth] == true) ? visited5p3 : visited3p5;

      if (visited.find(c->pathIID[c->pathDepth]) != visited.end()) {
        //if (VERBOSE2)
        //  fprintf(stderr, "DONE %5u/%s\n", c->pathIID[c->pathDepth], (c->path5p3[c->pathDepth] == true) ? "5'3'" : "3'5'");

        c->pathPosn[c->pathDepth]++;
        continue;
      }

      //
      //  Are we finished??
      //

      if ((c->pathIID[c->pathDepth] == fi[c->iid].mateIID) &&
          (c->path5p3[c->pathDepth] == end5p3) &&
          (c->pathLen[c->pathDepth] >= pathMin) &&
          (c->pathLen[c->pathDepth] <= pathMax)) {
        c->pathFound = true;
        return;
      }

      //
      //  Try extending into the fragment at this overlap
      //

      overlapInfo  *novl = c->pathRoot[c->pathDepth] + c->pathPosn[c->pathDepth];

      uint32        niid = novl->iid;                                                         //  Next fragment
      bool          n5p3 = (novl->flipped) ? (!c->path5p3[c->pathDepth]) : (c->path5p3[c->pathDepth]);    //  Next fragment orientation;
      uint32        nlen = 0;                                                                 //  New length of the path, if zero, can't extend

      //if (VERBOSE3)
      //  fprintf(stderr, "TRY  %5u/%s -> %5u/%s.\n",
      //          c->pathIID[c->pathDepth], (c->path5p3[c->pathDepth] == true) ? "5'3'" : "3'5'",
      //          niid,               (n5p3               == true) ? "5'3'" : "3'5'");

      //  Frag is forward, overlap is same, hang is positive
      if ((c->path5p3[c->pathDepth] == true) && (novl->flipped == false)) {
        nlen = c->pathLen[c->pathDepth] + novl->bhang;
        assert(n5p3 == true);
      }

      //  Frag is forward, overlap is flipped, hang is positive
      if ((c->path5p3[c->pathDepth] == true) && (novl->flipped == true)) {
        nlen = c->pathLen[c->pathDepth] + novl->bhang;
        assert(n5p3 == false);
      }

      //  Frag is reverse, overlap is same, hang is positive
      if ((c->path5p3[c->pathDepth] == false) && (novl->flipped == false)) {
        nlen = c->pathLen[c->pathDepth] + -novl->ahang;
        assert(n5p3 == false);
      }

      //  Frag is reverse, overlap is flipped, hang is positive
      if ((c->path5p3[c->pathDepth] == false) && (novl->flipped == true)) {
        nlen = c->pathLen[c->pathDepth] + -novl->ahang;
        assert(n5p3 == true);
      }

      //
      //  Are we _now_ finished?  This duplicates the test above, and ALSO tests contained
      //  overlaps that we do not extend into (when nlen <= pathLen).
      //

      if ((niid == fi[c->iid].mateIID) &&
          (n5p3 == end5p3) &&
          (nlen >= pathMin) &&
          (nlen <= pathMax)) {
        c->pathDepth++;

        c->pathIID[c->pathDepth]  = niid;
        c->path5p3[c->pathDepth]  = n5p3;
        c->pathLen[c->pathDepth]  = nlen;
        c->pathRoot[c->pathDepth] = oiPos[niid];
        c->pathPosn[c->pathDepth] = 0;
        c->pathMaxp[c->pathDepth] = oiLen[niid];

        c->pathFound = true;
        return;
      }

      //
      //  If the length is not too large, extend into it, otherwise, try the next fragment.
      //

      if ((nlen      >  c->pathLen[c->pathDepth]) &&
          (nlen      <= pathMax) &&
          (c->pathDepth <= depthMax)) {
        c->pathDepth++;

        c->pathIID[c->pathDepth]  = niid;
        c->path5p3[c->pathDepth]  = n5p3;
        c->pathLen[c->pathDepth]  = nlen;
        c->pathRoot[c->pathDepth] = oiPos[niid];
        c->pathPosn[c->pathDepth] = 0;
        c->pathMaxp[c->pathDepth] = oiLen[niid];

        //if (VERBOSE4)
        //  fprintf(stderr, "WALK %5u/%s with %5u overlaps at depth %2u of length %5d.\n",
        //          c->pathIID[c->pathDepth], (c->path5p3[c->pathDepth] == true) ? "5'3'" : "3'5'",
        //          c->pathMaxp[c->pathDepth], c->pathDepth, c->pathLen[c->pathDepth]);

        continue;
      }

      //if (VERBOSE5)
      //  fprintf(stderr, "ABRT %5u/%s\n", c->pathIID[c->pathDepth], (c->path5p3[c->pathDepth] == true) ? "5'3'" : "3'5'");

      //  Move to the next overlap for this fragment.
      c->pathPosn[c->pathDepth]++;
    }

    //  We've exhausted the paths for this fragment.  Mark it finished and back up one.

    set<uint32> &visited = (c->path5p3[c->pathDepth] == true) ? visited5p3 : visited3p5;
    visited.insert(c->pathIID[c->pathDepth]);

    c->pathDepth--;
    c->pathPosn[c->pathDepth]++;
  }

  c->pathFound = false;
  return;
}






void
cmWorker(void *G, void *T, void *S) {
  cmGlobalData    *g = (cmGlobalData  *)G;
  cmThreadData    *t = (cmThreadData  *)T;
  cmComputation   *s = (cmComputation *)S;

  g->doSearch(s, t);
}


void *
cmReader(void *G) {
  cmGlobalData    *g = (cmGlobalData  *)G;
  cmComputation   *s = NULL;

  for (; g->curFrag < g->numFrags; g->curFrag++) {
    if (g->fi[g->curFrag].doSearch == false)
      continue;

    if (g->fi[g->curFrag].mateIID == 0)
      continue;

    if (g->fi[g->curFrag].mateIID < g->curFrag)
      continue;

    s = new cmComputation(g->curFrag++);

    break;
  }

  return(s);
}



void
cmWriter(void *G, void *S) {
  cmGlobalData    *g = (cmGlobalData  *)G;
  cmComputation   *s = (cmComputation *)S;

  if (s->pathFound == true)
    fprintf(g->resultsFile, "Path from %d/%s to %d/%s found at depth %d of length %d.\n",
            s->iid, (s->path5p3[1] == true) ? "5'3'" : "3'5'", 
            s->pathIID[s->pathDepth], (s->path5p3[s->pathDepth] == true) ? "5'3'" : "3'5'", s->pathDepth, s->pathLen[s->pathDepth]);
  else
    fprintf(g->resultsFile, "Path from %d/%s NOT FOUND.\n",
            s->iid, (s->path5p3[1] == true) ? "5'3'" : "3'5'");

  delete s;
}




int
main(int argc, char **argv) {
  char      *gkpStorePath      = NULL;
  char      *ovlStorePath      = NULL;
  char      *resultsPath       = NULL;

  uint32     pathMin           = 500;    //  100
  uint32     pathMax           = 6000;   //  800
  bool       pathInnie         = false;  //  require mates to be innie

  uint32     depthMax          = 20;     //  100

  bool       searchLibs        = false;
  uint32     searchLib[1024]   = {0};

  bool       backboneLibs      = false;
  uint32     backboneLib[1024] = {0};

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      resultsPath = argv[++arg];

    } else if (strcmp(argv[arg], "-sl") == 0) {
      searchLibs = true;
      searchLib[atoi(argv[++arg])] = 1;

      char *a = argv[arg];
      char *b = strchr(a, '-');

      if (b) {
        uint32 bgn = atoi(a);
        uint32 end = atoi(b+1);

        for (uint32 i=bgn; i<=end; i++)
          searchLib[i] = 1;
      }
      
    } else if (strcmp(argv[arg], "-bl") == 0) {
      backboneLibs = true;
      backboneLib[atoi(argv[++arg])] = 1;

      char *a = argv[arg];
      char *b = strchr(a, '-');

      if (b) {
        uint32 bgn = atoi(a);
        uint32 end = atoi(b+1);

        for (uint32 i=bgn; i<=end; i++)
          backboneLib[i] = 1;
      }

    } else if (strcmp(argv[arg], "-min") == 0) {
      pathMin = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      pathMax = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-innie") == 0) {
      pathInnie = true;

    } else if (strcmp(argv[arg], "-outtie") == 0) {
      pathInnie = false;

    } else if (strcmp(argv[arg], "-depth") == 0) {
      depthMax = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (gkpStorePath == 0L)
    fprintf(stderr, "No gatekeeper store (-G) supplied.\n"), err++;
  if (ovlStorePath == 0L)
    fprintf(stderr, "No overlap store (-O) supplied.\n"), err++;
  if (resultsPath == 0L)
    fprintf(stderr, "No results output (-o) supplied.\n"), err++;
  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore -o resultsFile ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -sl l[-m]    Search for mates in libraries l-m\n");
    fprintf(stderr, "  -bl l[-m]    Use libraries l-m for searching\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -min m       Mates must be at least m bases apart\n");
    fprintf(stderr, "  -max m       Mates must be at most m bases apart\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -innie       Mates must be innie\n");
    fprintf(stderr, "  -outtie      Mates must be outtie\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -depth d     Search to at most d overlaps\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  //  Open this EARLY so we can fail without loading fragments or overlaps.
  errno = 0;
  FILE          *resultsFile = fopen(resultsPath, "w");
  if (errno)
    fprintf(stderr, "Failed to open results file '%s': %s\n", resultsPath, strerror(errno)), exit(1);


  cmGlobalData  *g = new cmGlobalData(resultsFile, pathMin, pathMax, pathInnie, depthMax);

  g->loadFragments(gkpStorePath, searchLibs, searchLib, backboneLibs, backboneLib);
  g->loadOverlaps(ovlStorePath);

  sweatShop *ss = new sweatShop(cmReader, cmWorker, cmWriter);

  ss->setLoaderQueueSize(16384);
  ss->setWriterQueueSize(1024);

  ss->setNumberOfWorkers(4);

  for (u32bit w=0; w<4; w++)
    ss->setThreadData(w, new cmThreadData());  //  these leak

  ss->run(g, true);

  fclose(resultsFile);

  delete g;
}

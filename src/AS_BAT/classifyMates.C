
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

const char *mainid = "$Id: classifyMates.C,v 1.2 2011-02-28 17:18:06 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"

#include <set>

using namespace std;


bool   VERBOSE1 = false;
bool   VERBOSE2 = false;
bool   VERBOSE3 = false;
bool   VERBOSE4 = false;
bool   VERBOSE5 = false;

uint32 MINPATHLENGTH  = 100;
uint32 MAXPATHLENGTH  = 500;
uint32 MAXPATHDEPTH   =   4;


class fragmentInfo {
public:
  fragmentInfo() {
  };
  ~fragmentInfo() {
  };

  uint64  unused      : 16;
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

#define NDISTD  2  //  Number of dovetail edges off each end
#define NDISTC  5  //  Number of contains for each end

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


int
main(int argc, char **argv) {
  char      *gkpStorePath = NULL;
  char      *ovlStorePath = NULL;
  char      *resultsPath  = NULL;

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

    } else {
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
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore -o resultsFile\n", argv[0]);
    exit(1);
  }

  //  Open this EARLY so we can fail without loading fragments or overlaps.
  errno = 0;
  FILE          *resultsFile = fopen(resultsPath, "w");
  if (errno)
    fprintf(stderr, "Failed to open results file '%s': %s\n", resultsPath, strerror(errno)), exit(1);



  fprintf(stderr, "LOADING FRAGMENTS...\n");

  gkStore          *gkpStore  = new gkStore(gkpStorePath, FALSE, FALSE);
  gkStream         *gkpStream = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);
  gkFragment        frg;

  uint32            numFrags = gkpStore->gkStore_getNumFragments();

  fragmentInfo     *fi = new fragmentInfo [numFrags + 1];

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
    fi[fid].doSearch    = gkpStore->gkStore_getLibrary(lib)->doMerBasedTrimming;
    fi[fid].clearLength = cllen;
    fi[fid].mateIID     = frg.gkFragment_getMateIID();
  }

  delete gkpStream;
  delete gkpStore;

  fprintf(stderr, "LOADING FRAGMENTS...%u fragments loaded.\n", numFrags);


  //  Overlap storage.  Each fragment has an entry in oiPos that points to actual data in an array.
  //  There are oiLen overlaps at this location.  The array is NOT contiguous.
  //
  //  oiStorage is a list of allocated memory.  We don't really need to keep the list around, but
  //  we'll play nice and clean up at the end.  Each chunk of allocated memory holds a fixed number
  //  of overlaps -- here, 128 million.  There are 1024 chunks, so we can store 128 billion
  //  overlaps, far far more than any overlap store we've created (about 3.3TB of data).

  fprintf(stderr, "LOADING OVERLAPS...\n");

  OverlapStore     *ovlStore  = AS_OVS_openOverlapStore(ovlStorePath);

  uint32            ovlMax = 1048576;
  uint32            ovlLen = 1;
  OVSoverlap       *ovl    = new OVSoverlap [ovlMax];


  overlapInfo     **oiPos = new overlapInfo * [numFrags + 1];  //  Pointer to start of overlaps for this frag
  uint32           *oiLen = new uint32        [numFrags + 1];  //  Number of overlaps for this frag

  memset(oiPos, 0, sizeof(overlapInfo) * (numFrags + 1));
  memset(oiLen, 0, sizeof(uint32)      * (numFrags + 1));

  uint32            oiStorageMax = 1024;
  uint32            oiStorageLen = 0;
  overlapInfo     **oiStorageArr = new overlapInfo * [oiStorageMax];
  overlapInfo      *oiStorage    = 0L;

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

  fprintf(stderr, "\nLOADING OVERLAPS...%lu overlaps loaded.\n", numOvl);


  uint32         pathIID[1024] = {0};   //  Fragment here
  uint32         path5p3[1024] = {0};   //  Fragment here is 5' to 3'

  uint32         pathLen[1024] = {0};   //  The length, in bp, of the path up till now

  overlapInfo   *pathRoot[1024] = {0};  //  Root of the list of overlaps for this iid
  uint32         pathPosn[1024] = {0};  //  Position we are at in the list of overlaps
  uint32         pathMaxp[1024] = {0};  //  Number of ovelerlaps for this iid

  for (uint32 iid=1; iid<=numFrags; iid++) {
    if (fi[iid].doSearch == false)
      continue;

    if (fi[iid].mateIID < iid)
      continue;

    //fprintf(stderr, "SEARCHING fragment %u\n", iid);

    //  Attempt to find a path from the 5' end of this fragment to the 5' end of the next iid
    //         <---- -path- ---->
    //  If we find a path, declare this a PE and not a MP.
    //
    //  The search is depth first, stopping when we find a path, or when the path gets implausibly long.

    bool         pathFound = false;
    set<uint32>  visited5p3;
    set<uint32>  visited3p5;

    uint32       pathDepth = 1;

    pathIID[pathDepth]  = iid;
    path5p3[pathDepth]  = false;
    pathLen[pathDepth]  = fi[iid].clearLength;
    pathRoot[pathDepth] = oiPos[iid];
    pathPosn[pathDepth] = 0;
    pathMaxp[pathDepth] = oiLen[iid];

    if (VERBOSE1)
      fprintf(stderr, "FRAG %5u/%s with %5u overlaps (len %u).\n",
              iid, (path5p3[pathDepth] == true) ? "5'3'" : "3'5'", oiLen[iid], pathLen[pathDepth]);

    while (pathDepth > 0) {
      while (pathPosn[pathDepth] < pathMaxp[pathDepth]) {

        //
        //  If we've already seen this fragment in this orientation, get out of here.
        //

        set<uint32> &visited = (path5p3[pathDepth] == true) ? visited5p3 : visited3p5;

        if (visited.find(pathIID[pathDepth]) != visited.end()) {
          if (VERBOSE2)
            fprintf(stderr, "DONE %5u/%s\n", pathIID[pathDepth], (path5p3[pathDepth] == true) ? "5'3'" : "3'5'");

          pathPosn[pathDepth]++;
          continue;
        }

        //
        //  Are we finished??
        //

        if ((pathIID[pathDepth] == fi[iid].mateIID) &&
            (path5p3[pathDepth] == true) &&
            (pathLen[pathDepth] >= MINPATHLENGTH) &&
            (pathLen[pathDepth] <= MAXPATHLENGTH)) {
          fprintf(resultsFile, "Path from %d/%s to %d/%s found at depth %d of length %d.\n",
                  iid, (path5p3[1] == true) ? "5'3'" : "3'5'", 
                  pathIID[pathDepth], (path5p3[pathDepth] == true) ? "5'3'" : "3'5'", pathDepth, pathLen[pathDepth]);
          pathFound = true;
          pathDepth = 1;
          break;
        }

        //
        //  Try extending into the fragment at this overlap
        //

        overlapInfo  *novl = pathRoot[pathDepth] + pathPosn[pathDepth];

        uint32        niid = novl->iid;                                                         //  Next fragment
        bool          n5p3 = (novl->flipped) ? (!path5p3[pathDepth]) : (path5p3[pathDepth]);    //  Next fragment orientation;
        uint32        nlen = 0;                                                                 //  New length of the path, if zero, can't extend

        if (VERBOSE3)
          fprintf(stderr, "TRY  %5u/%s -> %5u/%s.\n",
                  pathIID[pathDepth], (path5p3[pathDepth] == true) ? "5'3'" : "3'5'",
                  niid,               (n5p3               == true) ? "5'3'" : "3'5'");

        //  Frag is forward, overlap is same, hang is positive
        if ((path5p3[pathDepth] == true) && (novl->flipped == false)) {
          nlen = pathLen[pathDepth] + novl->bhang;
          assert(n5p3 == true);
        }

        //  Frag is forward, overlap is flipped, hang is positive
        if ((path5p3[pathDepth] == true) && (novl->flipped == true)) {
          nlen = pathLen[pathDepth] + novl->bhang;
          assert(n5p3 == false);
        }

        //  Frag is reverse, overlap is same, hang is positive
        if ((path5p3[pathDepth] == false) && (novl->flipped == false)) {
          nlen = pathLen[pathDepth] + -novl->ahang;
          assert(n5p3 == false);
        }

        //  Frag is reverse, overlap is flipped, hang is positive
        if ((path5p3[pathDepth] == false) && (novl->flipped == true)) {
          nlen = pathLen[pathDepth] + -novl->ahang;
          assert(n5p3 == true);
        }

        //
        //  Are we _now_ finished?  This duplicates the test above, and ALSO tests contained
        //  overlaps that we do not extend into (when nlen <= pathLen).
        //

        if ((niid == fi[iid].mateIID) &&
            (n5p3 == true) &&
            (nlen >= MINPATHLENGTH) &&
            (nlen <= MAXPATHLENGTH)) {
          fprintf(resultsFile, "Path from %d/%s to %d/%s found at depth %d of length %d.\n",
                  iid, (path5p3[1] == true) ? "5'3'" : "3'5'", 
                  niid, (n5p3 == true) ? "5'3'" : "3'5'", pathDepth+1, nlen);
          pathFound = true;
          pathDepth = 1;
          break;
        }

        //
        //  If the length is not too large, extend into it, otherwise, try the next fragment.
        //

        if ((nlen      >  pathLen[pathDepth]) &&
            (nlen      <= MAXPATHLENGTH) &&
            (pathDepth <= MAXPATHDEPTH)) {
          pathDepth++;

          pathIID[pathDepth]  = niid;
          path5p3[pathDepth]  = n5p3;
          pathLen[pathDepth]  = nlen;
          pathRoot[pathDepth] = oiPos[niid];
          pathPosn[pathDepth] = 0;
          pathMaxp[pathDepth] = oiLen[niid];

          if (VERBOSE4)
            fprintf(stderr, "WALK %5u/%s with %5u overlaps at depth %2u of length %5d.\n",
                    pathIID[pathDepth], (path5p3[pathDepth] == true) ? "5'3'" : "3'5'",
                    pathMaxp[pathDepth], pathDepth, pathLen[pathDepth]);

          continue;
        }

        if (VERBOSE5)
          fprintf(stderr, "ABRT %5u/%s\n", pathIID[pathDepth], (path5p3[pathDepth] == true) ? "5'3'" : "3'5'");

        //  Move to the next overlap for this fragment.
        pathPosn[pathDepth]++;
      }

      //  We've exhausted the paths for this fragment.  Mark it finished and back up one.

      set<uint32> &visited = (path5p3[pathDepth] == true) ? visited5p3 : visited3p5;
      visited.insert(pathIID[pathDepth]);

      pathDepth--;
      pathPosn[pathDepth]++;
    }

#if 1
    if (pathFound == false)
      fprintf(resultsFile, "Path from %d/%s NOT FOUND.\n",
              iid, (path5p3[1] == true) ? "5'3'" : "3'5'");
#endif
  }

  fclose(resultsFile);

  for (uint32 i=0; i<oiStorageMax; i++)
    delete [] oiStorageArr[i];

  delete [] oiStorageArr;

  delete [] oiPos;
  delete [] oiLen;

  delete [] ovl;

  delete [] fi;
}

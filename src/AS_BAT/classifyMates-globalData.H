
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

#ifndef CLASSIFYMATES_GLOBALDATA_H
#define CLASSIFYMATES_GLOBALDATA_H

#include "AS_global.H"
#include "AS_OVS_overlapStore.H"
#include "AS_PER_gkpStore.H"

#include "memoryMappedFile.H"

#include "classifyMates-runningTime.H"

#include <vector>
#include <map>
#include <set>
using namespace std;

class fragmentInfo {
public:
  uint64  libIID      : 12;
  uint64  contained   : 1;
  uint64  end5covered : 1;
  uint64  end3covered : 1;
  uint64  isBackbone  : 1;
  uint64  doSearch    : 1;
  uint64  clearLength : 15;
  uint64  mateIID     : 32;
};


//  Save distance needs to know fragmentInfo.  Ugly.
#include "classifyMates-saveDistance.H"

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

#if AS_OVS_HNGBITS < 16
  uint32  unused  : 32 - 1 - AS_OVS_HNGBITS - AS_OVS_HNGBITS;
  uint32  flipped : 1;                //  B-fragment is flipped
  int32   ahang   : AS_OVS_HNGBITS;   //
  int32   bhang   : AS_OVS_HNGBITS;   //
#else
  uint64  unused  : 64 - 1 - AS_OVS_HNGBITS - AS_OVS_HNGBITS;
  uint64  flipped : 1;                //  B-fragment is flipped
  int64   ahang   : AS_OVS_HNGBITS;   //
  int64   bhang   : AS_OVS_HNGBITS;   //
#endif

  uint32 iid;             //  B-fragment ID
};




class searchNode {
public:
  uint32         pIID;   //  Fragment here
  uint32         p5p3;   //  Fragment here is 5' to 3'

  int32          pLen;   //  The length, in bp, of the path up till now
  uint32         pDpt;   //  The depth, in fragments, of the path up till now

  uint32         oMax;  //  Number of ovelerlaps for this iid
  uint32         oPos;  //  Position we are at in the list of overlaps
  overlapInfo   *oLst;  //  Root of the list of overlaps for this iid
};



class bits {
public:
  bits(uint64 num, char *label) {
    nm = num / 32 + 1;
    ba = new uint32 [nm];
    memset(ba, 0, sizeof(uint32) * nm);
    strcpy(la, label);
  };
  ~bits() {
    delete [] ba;
  };

  bool      isSet(uint64 x) {
    uint64 wrd = (x >> 5);
    uint64 bit = (x     ) & 0x000000000000001fllu;

    return((ba[wrd] >> bit) & 0x00000001llu);
  };

  void      set(uint64 x) {
    uint64 wrd = (x >> 5);
    uint64 bit = (x     ) & 0x000000000000001fllu;

    ba[wrd] = ba[wrd] | (0x00000001llu << bit);
  };

  void      clear(uint64 x) {
    uint64 wrd = (x >> 5);
    uint64 bit = (x     ) & 0x000000000000001fllu;

    ba[wrd] = ba[wrd] & ~(0x00000001llu << bit);
  };

#if 1
  void      testClear(void) {};
#else
  void      testClear(void) {
    for (uint32 ii=0; ii<nm; ii++) {
      if (ba[ii] != 0)
        fprintf(stderr, "ERROR: bits::%s %u not zero.\n", la, ii);
      assert(ba[ii] == 0);
    }
  };
#endif

private:
  uint32   *ba;
  uint32    nm;
  char      la[64];
};




class cmThreadData {
public:
  cmThreadData() {
    clear();

    pathMax      = 0;
    path         = NULL;

    solutionSet  = NULL;

    extMax       = 131072;  //  Overlaps out of a node
    ext          = new uint32 [extMax];

    edge3Max     = 131072;  //  Overlaps off end of a read
    edge3        = new uint32 [edge3Max];

    edge5Max     = 131072;
    edge5        = new uint32 [edge5Max];
  };
  ~cmThreadData() {
    delete [] path;
    delete [] solutionSet;
    delete [] ext;
    delete [] edge3;
    delete [] edge5;
  };

  void clear(void) {
    pathPos        = 0;
    pathAdd        = 0;
    //pathMax      = 0;      //  NO!  We reuse this!
    //path         = NULL;

    searchIter     = 0;

    extLen         = 0;

    edge3Len       = 0;
    edge5Len       = 0;

    visitedListLen = 0;

    //if (solutionSet)    solutionSet->testClear();  //  Not valid; we test explicitly.

    if (visited5p3bits) visited5p3bits->testClear();
    if (visited3p5bits) visited3p5bits->testClear();
  };

  uint32         pathPos;  //  Position of the search in DFS and BFS
  uint32         pathAdd;  //  Next free spot to add a fragment in BFS
  uint32         pathMax;  //  Number of nodes we have allocated.
  searchNode    *path;

  bits          *solutionSet;

  bits          *visited5p3bits;
  bits          *visited3p5bits;
  uint32         visitedListMax;
  uint32         visitedListLen;
  uint32        *visitedList;

  uint32         searchIter;
  
  uint32         extLen;   //  Use in RFS, list of edges out of a node
  uint32         extMax;
  uint32        *ext;

  uint32         edge3Len;  //  Use in testChimer
  uint32         edge3Max;
  uint32        *edge3;

  uint32         edge5Len;  //  Use in testChimer
  uint32         edge5Max;
  uint32        *edge5;

  vector<int32>  distInnieClose;  //  Use in testSuspicious
  vector<int32>  distOuttieClose;
  vector<int32>  distNormalClose;
  vector<int32>  distAntiClose;

  vector<int32>  distInnieFar;
  vector<int32>  distOuttieFar;
  vector<int32>  distNormalFar;
  vector<int32>  distAntiFar;
};


class cmComputation {
public:

  //  If pathInnie == false, then we're attempting to find a path for outtie oriented fragments.
  //  In this case, we start with the first fragment 5p3=false, and need to end with 5p3=true.

  cmComputation(uint32 iid, uint32 mid, bool innie) {
    fragIID         = iid;
    frag5p3         = (innie == false) ? false : true;

    mateIID         = mid;
    mate5p3         = (innie == false) ? true : false;

    nInnie          = 0;
    nNormal         = 0;
    nAnti           = 0;
    nOuttie         = 0;

    result.innie    = innie;

    result.readIID  = fragIID;
    result.mateIID  = mateIID;
  };
  ~cmComputation() {
  };

public:
  uint32                fragIID;
  bool                  frag5p3;

  uint32                mateIID;
  bool                  mate5p3;

  uint32                nInnie;
  uint32                nNormal;
  uint32                nAnti;
  uint32                nOuttie;

  classifyMatesResult   result;
};



class cmGlobalData {
public:
  cmGlobalData(char    *resultName_,
               uint32   distMin_,
               uint32   distMax_,
               bool     innie_,
               uint32   nodesMax_,
               uint32   depthMax_,
               uint32   pathsMax_,
               uint64   memoryLimit_,
               bool     doSearchSuspicious_);
  ~cmGlobalData();

  bool   load(set<AS_IID>  &searchLibs,
              set<AS_IID>  &backboneLibs,
              double        maxErrorFraction,
              bool          cacheOnly);
  void   save(void);

  bool   checkBB(void);
  bool   checkTG(void);
  bool   checkGT(void);

  void   buildBBSSmap(gkStore     *gkpStore,
                      set<AS_IID> &searchLibs,
                      set<AS_IID> &backboneLibs);

  void   loadFragments(char         *gkpStoreName,
                       set<AS_IID>  &searchLibs,
                       set<AS_IID>  &backboneLibs);


  void   allocateOverlapPointers(void);

  void   resetOverlapStoreRange(OverlapStore *ovlStore,
                                AS_IID curFragIID, AS_IID &minFragIID, AS_IID &maxFragIID);

  void   loadOverlaps(char   *ovlStoreName,
                      double  maxErrorFraction);
  void   loadOverlaps_invert(void);

  void   computeNextPlacement(cmThreadData  *t,
                              overlapInfo  *&novl,
                              uint32        &niid,
                              bool          &n5p3,
                              int32         &nlen);

  void   testSuspicious(uint32           fragIID,
                        bool             frag5p3,
                        uint32           mateIID,
                        cmThreadData    *t,
                        vector<int32>   &dist3p5,
                        vector<int32>   &dist5p3,
                        bool             onlyNormal);

  bool   testSearch(cmComputation  *c,
                    cmThreadData   *t);

  bool   testSearch(cmComputation *c,
                    cmThreadData  *t,
                    overlapInfo  **pos,
                    uint32        *len);

  bool   testSpur(cmComputation *c, cmThreadData *t);

  bool   testChimer(uint32 iid, cmThreadData *t);
  bool   testChimers(cmComputation *c, cmThreadData *t);

  void   doSearchDFS(cmComputation *c, cmThreadData *t);
  void   doSearchBFS(cmComputation *c, cmThreadData *t);
  void   doSearchRFS(cmComputation *c, cmThreadData *t);

  void   doSearchSuspicious(uint32 fragIID, bool frag5p3,
                            uint32 mateIID,
                            cmThreadData   *t,
                            vector<int32>  &dist5p3close,
                            vector<int32>  &dist5p3far,
                            vector<int32>  &dist3p5close,
                            vector<int32>  &dist3p5far);

public:
  char                      resultPrefix[FILENAME_MAX];
  classifyMatesResultFile  *resultOutput;
  FILE                     *resultSuspicious;

  uint32                    distMin;
  uint32                    distMax;
  bool                      innie;

  uint32                    nodesMax;
  uint32                    depthMax;
  uint32                    pathsMax;

  bool                      suspiciousSearch;

  memoryMappedFile         *cmDat;
  memoryMappedFile         *cmOvl;
  memoryMappedFile         *cmInv;

  bool                      saveCache;

  uint32                    numFrags;
  uint32                    numLibs;
  fragmentInfo             *fi;

  uint32                    curFragIID;

  uint32                   *isBB;
  uint32                   *isSS;

  overlapInfo             **bbPos;  //  Pointer to start of overlaps for this BACKBONE frag
  uint32                   *bbLen;  //  Number of overlaps for this BACKBONE frag

  overlapInfo             **tgPos;  //  Pointer to start of overlaps for this TARGET frag
  uint32                   *tgLen;  //  Number of overlaps for this TARGET frag

  overlapInfo             **gtPos;  //  Same as tgPos, but indexed on the b-frag IID
  uint32                   *gtLen;  //

  uint64                    memoryLimit;
  uint64                    memoryUsed;

  vector<overlapInfo *>     oiStorageArr;  //  List of allocated blocks.

  onlineMeanStdDev          runTime;
};


#endif  //  CLASSIFYMATES_GLOBALDATA_H

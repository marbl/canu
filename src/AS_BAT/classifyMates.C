
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

const char *mainid = "$Id: classifyMates.C,v 1.16 2011-06-03 17:34:19 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"

#include <set>
#include <map>

#include "sweatShop.H"

using namespace std;


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


#include "classifyMates-saveDistance.C"


class searchNode {
public:
  uint32         pIID;   //  Fragment here
  uint32         p5p3;   //  Fragment here is 5' to 3'

  uint32         pLen;   //  The length, in bp, of the path up till now

  uint32         oMax;  //  Number of ovelerlaps for this iid
  uint32         oPos;  //  Position we are at in the list of overlaps
  overlapInfo   *oLst;  //  Root of the list of overlaps for this iid
};


class cmThreadData {
public:
  cmThreadData() {
    pathPos = 0;
    pathAdd = 0;
    pathMax = 0;
    path    = NULL;

    extMax = 1048576;
    extLen = 0;
    ext    = new uint32 [extMax];
  };
  ~cmThreadData() {
    delete [] path;
    delete [] ext;
  };

  //  Use in DFS and BFS.
  //    In DFS, the path to the node we are searching.
  //    In BFS, the list of nodes we have/will search.
  uint32         pathPos;  //  Position of the search in DFS and BFS
  uint32         pathAdd;  //  Next free spot to add a fragment in BFS
  uint32         pathMax;  //  Number of nodes we have allocated.
  searchNode    *path;
  
  //  Use in RFS, list of edges out of a node
  uint32         extLen;
  uint32         extMax;
  uint32        *ext;
};


class cmComputation {
public:

  //  If pathInnie == false, then we're attempting to find a path for outtie oriented fragments.
  //  In this case, we start with the first fragment 5p3=false, and need to end with 5p3=true.

  cmComputation(uint32 iid, uint32 mid, bool innie) {
    fragIID    = iid;
    frag5p3    = (innie == false) ? false : true;

    mateIID    = mid;
    mate5p3    = (innie == false) ? true : false;

    sFound     = false;
    sLimited   = false;
    sExhausted = false;
    sPos       = 0;
    sLen       = 0;
  };
  ~cmComputation() {
  };

public:
  uint32         fragIID;
  bool           frag5p3;

  uint32         mateIID;     //  Fragment here
  uint32         mate5p3;     //  Fragment here is 5' to 3'

  bool           sFound;      //  Did we find an answer?
  bool           sLimited;    //  Search stopped due to CPU limits.
  bool           sExhausted;  //  Search stopped due to no more overlaps.
  uint32         sPos;        //  If answer, the position we found it at.
  uint32         sLen;        //  If answer, the length of the path in bp.
};


class cmGlobalData {
public:
  cmGlobalData(char    *resultsPath_,
               uint32   distMin_,
               uint32   distMax_,
               bool     innie_,
               uint32   nodesMax_,
               uint32   depthMax_,
               uint32   pathsMax_) {

    strcpy(resultsPrefix, resultsPath_);

    //  Open this EARLY so we can fail without loading fragments or overlaps.
    errno = 0;
    resultsFile = fopen(resultsPath_, "w");
    if (errno)
      fprintf(stderr, "Failed to open results file '%s': %s\n", resultsPath_, strerror(errno)), exit(1);

    distMin      = distMin_;
    distMax      = distMax_;
    innie        = innie_;

    nodesMax     = nodesMax_;
    depthMax     = depthMax_;
    pathsMax     = pathsMax_;

    numFrags     = 0;
    fi           = 0L;

    curFragIID   = 0;
    minFragIID   = UINT32_MAX;
    maxFragIID   = 0;

    bbPos        = 0L;
    bbLen        = 0L;

    tgPos        = 0L;
    tgLen        = 0L;

    gtPos        = 0L;
    gtLen        = 0L;

    oiStorageMax = 0;
    oiStorageLen = 0;
    oiStorageArr = 0L;
    oiStorage    = 0L;
  };

  ~cmGlobalData() {

    fclose(resultsFile);

    for (uint32 i=0; i<oiStorageMax; i++)
      delete [] oiStorageArr[i];

    delete [] bbPos;
    delete [] bbLen;

    delete [] tgPos;
    delete [] tgLen;

    delete [] gtPos;
    delete [] gtLen;

    delete [] oiStorageArr;

    delete [] fi;
  };

  void   loadFragments(char    *gkpStorePath,
                       bool     searchLibs,
                       uint32  *searchLib,
                       bool     backboneLibs,
                       uint32  *backboneLib);

  void   loadOverlaps(char   *ovlStorePath);

  void   computeNextPlacement(cmComputation *c,
                              cmThreadData  *t,
                              overlapInfo  *&novl,
                              uint32        &niid,
                              bool          &n5p3,
                              uint32        &nlen);

  bool   testSearch(cmComputation              *c,
                    cmThreadData               *t,
                    map<uint32,overlapInfo*>   &pos);
  bool   testSearch(cmComputation *c,
                    cmThreadData  *t,
                    overlapInfo  **pos,
                    uint32        *len);

  void   doSearchDFS(cmComputation  *c, cmThreadData   *t);
  void   doSearchBFS(cmComputation  *c, cmThreadData   *t);
  void   doSearchRFS(cmComputation  *c, cmThreadData   *t);

public:
  char              resultsPrefix[FILENAME_MAX];
  FILE             *resultsFile;

  uint32            distMin;
  uint32            distMax;
  bool              innie;

  uint32            nodesMax;
  uint32            depthMax;
  uint32            pathsMax;

  uint32            numFrags;
  fragmentInfo     *fi;

  uint32            curFragIID;
  uint32            minFragIID;
  uint32            maxFragIID;

  overlapInfo     **bbPos;  //  Pointer to start of overlaps for this BACKBONE frag
  uint32           *bbLen;  //  Number of overlaps for this BACKBONE frag

  overlapInfo     **tgPos;  //  Pointer to start of overlaps for this TARGET frag
  uint32           *tgLen;  //  Number of overlaps for this TARGET frag

  overlapInfo     **gtPos;  //  Same as tgPos, but indexed on the b-frag IID
  uint32           *gtLen;  //

  uint32            oiStorageMax;  //  Maximum number of blocks we can allocate.
  uint32            oiStorageLen;  //  Actual number of blocks allocated.
  overlapInfo     **oiStorageArr;  //  List of allocated blocks.
  overlapInfo      *oiStorage;     //  The current block -- doesn't need to be in the class.
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

  char  cacheName[FILENAME_MAX];
  sprintf(cacheName, "%s.fi", resultsPrefix);

  if (AS_UTL_fileExists(cacheName, FALSE, TRUE)) {
    fprintf(stderr, "LOADING FRAGMENTS...(from cache)\n");

    FILE *F = fopen(cacheName, "r");
    AS_UTL_safeRead(F, &minFragIID, "minFragIID", sizeof(uint32),       1);
    AS_UTL_safeRead(F, &maxFragIID, "maxFragIID", sizeof(uint32),       1);
    AS_UTL_safeRead(F,  fi,         "fi",         sizeof(fragmentInfo), numFrags+1);
    fclose(F);
  } else {
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
        fi[fid].doSearch = gkpStore->gkStore_getLibrary(lib)->doTrim_initialMerBased;
      else
        fi[fid].doSearch = searchLib[lib];

      if (backboneLibs == false)
        fi[fid].isBackbone = true;
      else
        fi[fid].isBackbone = backboneLib[lib];

      if ((fi[fid].doSearch || fi[fid].isBackbone) && (fid < minFragIID))
        minFragIID = fid;
      if ((fi[fid].doSearch || fi[fid].isBackbone) && (maxFragIID < fid))
        maxFragIID = fid;
    }

    FILE *F = fopen(cacheName, "w");
    AS_UTL_safeWrite(F, &minFragIID, "minFragIID", sizeof(uint32),       1);
    AS_UTL_safeWrite(F, &maxFragIID, "maxFragIID", sizeof(uint32),       1);
    AS_UTL_safeWrite(F,  fi,         "fi",         sizeof(fragmentInfo), numFrags+1);
    fclose(F);
  }

  delete gkpStream;
  delete gkpStore;

  fprintf(stderr, "LOADING FRAGMENTS...%u fragments loaded.\n", numFrags);
}


void
cmGlobalData::loadOverlaps(char  *ovlStorePath) {

  fprintf(stderr, "LOADING OVERLAPS...for fragments %u to %u\n", minFragIID, maxFragIID);

  uint32            ovlMax = 1048576;
  uint32            ovlLen = 1;
  OVSoverlap       *ovl    = new OVSoverlap [ovlMax];
  bool             *ovlBB  = new bool       [ovlMax];
  bool             *ovlTG  = new bool       [ovlMax];

  bbPos = new overlapInfo * [numFrags + 1];  //  Pointer to start of overlaps for this frag
  tgPos = new overlapInfo * [numFrags + 1];
  gtPos = new overlapInfo * [numFrags + 1];

  bbLen = new uint32        [numFrags + 1];  //  Number of overlaps for this frag
  tgLen = new uint32        [numFrags + 1];
  gtLen = new uint32        [numFrags + 1];

  memset(bbPos, 0, sizeof(overlapInfo) * (numFrags + 1));
  memset(tgPos, 0, sizeof(overlapInfo) * (numFrags + 1));
  memset(gtPos, 0, sizeof(overlapInfo) * (numFrags + 1));

  memset(bbLen, 0, sizeof(uint32)      * (numFrags + 1));
  memset(tgLen, 0, sizeof(uint32)      * (numFrags + 1));
  memset(gtLen, 0, sizeof(uint32)      * (numFrags + 1));

  //  Overlap storage.  Each fragment has an entry in bbPos & tgPos that points to actual data in an
  //  array.  There are bbLen & tgPos overlaps at this location.  The array is NOT contiguous.

  oiStorageMax = 1024;                              //  Number blocks of memory
  oiStorageLen = 1;
  oiStorageArr = new overlapInfo * [oiStorageMax];  //  List of allocated memory
  oiStorage    = 0L;

  memset(oiStorageArr, 0, sizeof(overlapInfo *) * oiStorageMax);

  uint32            oiStorageBS  = 128 * 1024 * 1024;  //  Fixed block size
  uint32            oiStoragePos = 0;                  //  Position in the current block

  oiStorage = oiStorageArr[0] = new overlapInfo [oiStorageBS];

  uint64            numTT = 0;  //  Number of overlaps read from disk
  uint64            numUU = 0;  //  Number of overlaps not backbone and not target
  uint64            numBB = 0;  //  Number of BB overlaps loaded
  uint64            numTG = 0;  //  Number of TG overlaps loaded
  uint64            numDD = 0;  //  Number of overlaps discarded

#if NDISTD > 0
  saveDistance     *dist = new saveDistance;
#else
  saveDistance     *dist = NULL;
#endif

  OverlapStore     *ovlStore  = AS_OVS_openOverlapStore(ovlStorePath);

  AS_OVS_setRangeOverlapStore(ovlStore, minFragIID, maxFragIID);

  ovlLen = AS_OVS_readOverlapsFromStore(ovlStore, ovl, ovlMax, AS_OVS_TYPE_OVL);

  uint64 maxError = AS_OVS_encodeQuality(3.0);

  while (ovlLen > 0) {
    numTT += ovlLen;

    //  Filter out overlaps we don't care about.

    uint32  iid = ovl[0].a_iid;

    //  Save if both fragments are in the backbone library.
    //  Save if one fragment is in the backbone library and the other is in the test library.
    //  Otherwise, we don't care about this overlap.
    //
    uint32  c=0;

    for (uint32 i=0; i<ovlLen; i++) {
      if (ovl[i].dat.ovl.orig_erate > maxError) {
        numUU++;
        continue;
      }

      if (((fi[ovl[i].a_iid].isBackbone == true) && (fi[ovl[i].b_iid].isBackbone == true)) ||
          ((fi[ovl[i].a_iid].isBackbone == true) && (fi[ovl[i].b_iid].doSearch   == true)) ||
          ((fi[ovl[i].a_iid].doSearch   == true) && (fi[ovl[i].b_iid].isBackbone == true))) {
        if (c != i)
          ovl[c] = ovl[i];
        c++;
      } else {
        numUU++;
      }
    }

    ovlLen = c;


    if (dist)
      dist->compute(fi, ovl, ovlLen);


    for (uint32 i=0; i<ovlLen; i++) {

      //  By default, keep no BB and no TG overlaps, and delete everything.
      ovlBB[i] = false;
      ovlTG[i] = false;


      //  Overlaps from a search A-fragment to a backbone B-fragment must always be kept.  These are
      //  used to initiating the search.
      //
      if ((fi[ovl[i].a_iid].doSearch   == true) && (fi[ovl[i].b_iid].isBackbone == true))
        ovlBB[i] = true;


      //  Overlaps from a backbone A-fragment to a search B-fragment must always be kept.  These are
      //  used to terminate the search.
      //
      if ((fi[ovl[i].a_iid].isBackbone == true) && (fi[ovl[i].b_iid].doSearch   == true))
        ovlTG[i] = true;


      //  Overlaps from backbone to backbone are kept if they are dovetail.
      //
      if ((fi[ovl[i].a_iid].isBackbone == true) && (fi[ovl[i].b_iid].isBackbone == true)) {
        int32  ah = ovl[i].dat.ovl.a_hang;
        int32  bh = ovl[i].dat.ovl.b_hang;
        int32  fa = fi[ovl[i].a_iid].clearLength;
        int32  fb = fi[ovl[i].b_iid].clearLength;

        if ((AS_OVS_overlapAEndIs5prime(ovl[i])) &&
            ((dist == NULL) ||
             (dist->doveDist5 <= fb - -ah)))
          //  ah < 0 && bh < 0
          ovlBB[i] = true;

        if ((AS_OVS_overlapAEndIs3prime(ovl[i])) &&
            ((dist == NULL) ||
             (dist->doveDist3 <= fb -  bh)))
          //  ah > 0 && bh > 0
          ovlBB[i] = true;
      }

      //  Save or discard the overlap.  We are allowed to save the overlap twice, once as a backbone
      //  overlap and once as a termination overlap.

      if ((ovlBB[i] == false) && (ovlTG[i] == false))
        ovl[i].a_iid = 0;

      if (ovlBB[i] == true)
        bbLen[iid]++;

      if (ovlTG[i] == true)
        tgLen[iid]++;
    }  //  Over all overlaps

    //  If not enough space in oiStorage, get more space.

    if (oiStoragePos + bbLen[iid] + tgLen[iid] > oiStorageBS) {
      oiStorage = oiStorageArr[oiStorageLen] = new overlapInfo [oiStorageBS];
      oiStorageLen++;
      oiStoragePos = 0;
    }

    //  Add the overlaps.

    bbPos[iid] = oiStorage + oiStoragePos;
    tgPos[iid] = oiStorage + oiStoragePos + bbLen[iid];

    bbLen[iid] = 0;
    tgLen[iid] = 0;

    for (uint32 i=0; i<ovlLen; i++) {
      if ((ovlBB[i] == false) && (ovlTG[i] == false)) {
        numDD++;
      }

      if (ovlBB[i] == true) {
        bbPos[iid][bbLen[iid]++] = overlapInfo(ovl[i]);
        numBB++;
      }

      if (ovlTG[i] == true) {
        tgPos[iid][tgLen[iid]++] = overlapInfo(ovl[i]);
        numTG++;
      }
    }

    oiStoragePos += bbLen[iid] + tgLen[iid];

    assert(oiStoragePos <= oiStorageBS);

    ovlLen = AS_OVS_readOverlapsFromStore(ovlStore, ovl, ovlMax, AS_OVS_TYPE_OVL);

    if ((iid % 100000) == 0)
      fprintf(stderr, "LOADING OVERLAPS...at IID %u (%06.2f%%): BB %lu (%06.2f%%) TG %lu (%06.2f%%) DD %lu (%06.2f%%).\r",
              iid,   100.0 * iid   / numFrags,
              numBB, 100.0 * numBB / numTT,
              numTG, 100.0 * numTG / numTT,
              numDD, 100.0 * numDD / numTT);
  }

  AS_OVS_closeOverlapStore(ovlStore);

  delete [] ovl;
  delete [] ovlBB;
  delete [] ovlTG;

  fprintf(stderr, "LOADING OVERLAPS...at IID %u (%06.2f%%): BB %lu (%06.2f%%) TG %lu (%06.2f%%) DD %lu (%06.2f%%).\n",
          numFrags, 100.0,
          numBB, 100.0 * numBB / numTT,
          numTG, 100.0 * numTG / numTT,
          numDD, 100.0 * numDD / numTT);
  fprintf(stderr, "LOADING OVERLAPS...%lu overlaps loaded.\n", numBB + numTG);


  fprintf(stderr, "INVERTING OVERLAPS.\n");

  //  Invert the TG overlaps to make the GT overlaps.  The current list is indexed on a-fragID, we
  //  need to index on b-fragID.

  uint64  numGTfrg = 0;
  uint64  numGTovl = 0;

  //  Count the number of overlaps for each b_iid
  for (uint32 ii=0; ii<numFrags+1; ii++) {
    if (tgPos[ii] == NULL)
      continue;

    for (uint32 jj=0; jj<tgLen[ii]; jj++) {
      gtLen[ tgPos[ii][jj].iid ]++;
      numGTovl++;
    }
  }

  //  Allocate space for them
  oiStorage = oiStorageArr[oiStorageLen] = new overlapInfo [numGTovl];
  oiStorageLen++;
  oiStoragePos = 0;

  //  Set pointers to the space.  We'll recount as the overlaps are added.
  for (uint32 ii=0; ii<numFrags+1; ii++) {
    if (gtLen[ii] == 0)
      continue;

    gtPos[ii]      = oiStorage + oiStoragePos;  //  Point into the space
    oiStoragePos  += gtLen[ii];                 //  Reserve the space
    gtLen[ii]      = 0;                         //  Show we have no overlaps loaded
  }

  //  Finally, copy the overlaps
  for (uint32 ii=0; ii<numFrags+1; ii++) {
    if (tgPos[ii] == NULL)
      continue;

    numGTfrg++;

    for (uint32 jj=0; jj<tgLen[ii]; jj++) {
      uint32  iid = tgPos[ii][jj].iid;

      assert(gtPos[iid] != NULL);

      gtPos[iid][gtLen[iid]]     = tgPos[ii][jj];
      gtPos[iid][gtLen[iid]].iid = ii;

      gtLen[iid]++;
    }
  }

  fprintf(stderr, "INVERTING OVERLAPS....%lu frags found.\n", numGTfrg);
}




void
cmGlobalData::computeNextPlacement(cmComputation *c,
                                   cmThreadData  *t,
                                   overlapInfo  *&novl,
                                   uint32        &niid,
                                   bool          &n5p3,
                                   uint32        &nlen) {

  //  Frag is forward, overlap is same, hang is positive
  if ((t->path[t->pathPos].p5p3 == true) && (novl->flipped == false)) {
    nlen = t->path[t->pathPos].pLen + novl->bhang;
    assert(n5p3 == true);
  }

  //  Frag is forward, overlap is flipped, hang is positive
  if ((t->path[t->pathPos].p5p3 == true) && (novl->flipped == true)) {
    nlen = t->path[t->pathPos].pLen + novl->bhang;
    assert(n5p3 == false);
  }

  //  Frag is reverse, overlap is same, hang is positive
  if ((t->path[t->pathPos].p5p3 == false) && (novl->flipped == false)) {
    nlen = t->path[t->pathPos].pLen + -novl->ahang;
    assert(n5p3 == false);
  }

  //  Frag is reverse, overlap is flipped, hang is positive
  if ((t->path[t->pathPos].p5p3 == false) && (novl->flipped == true)) {
    nlen = t->path[t->pathPos].pLen + -novl->ahang;
    assert(n5p3 == true);
  }
}




bool
cmGlobalData::testSearch(cmComputation              *c,
                         cmThreadData               *t,
                         map<uint32,overlapInfo*>   &pos) {

  uint32  thisIID = t->path[t->pathPos].pIID;

  map<uint32,overlapInfo*>::iterator  it = pos.find(thisIID);

  if (it == pos.end())
    return(false);

  //  We are guaranteed by construction that this overlap is to the mate fragment.
  //  See doSearchBFS.

  overlapInfo  *novl = it->second;
  uint32        niid = c->mateIID;
  bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
  uint32        nlen = 0;

  if (n5p3 != c->mate5p3)
    return(false);

  computeNextPlacement(c, t, novl, niid, n5p3, nlen);

  if ((nlen < distMin) ||
      (nlen > distMax))
    return(false);

  c->sFound = true;
  c->sPos = t->pathPos + 1;
  c->sLen   = nlen;

  return(true);
}

bool
cmGlobalData::testSearch(cmComputation  *c,
                         cmThreadData   *t,
                         overlapInfo   **pos,
                         uint32         *len) {

  uint32  thisIID = t->path[t->pathPos].pIID;

  for (uint32 test=0; test<len[thisIID]; test++) {
    overlapInfo  *novl = pos[thisIID] + test;
    uint32        niid = novl->iid;
    bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
    uint32        nlen = 0;

    if ((niid != c->mateIID) ||
        (n5p3 != c->mate5p3))
      continue;

    computeNextPlacement(c, t, novl, niid, n5p3, nlen);

    if ((nlen >= distMin) &&
        (nlen <= distMax)) {
      c->sFound = true;
      c->sPos = t->pathPos + 1;
      c->sLen   = nlen;

      return(true);
    }
  }
  return(false);
}



#include "classifyMates-DFS.C"
#include "classifyMates-BFS.C"
#include "classifyMates-RFS.C"



void
cmWorker(void *G, void *T, void *S) {
  cmGlobalData    *g = (cmGlobalData  *)G;
  cmThreadData    *t = (cmThreadData  *)T;
  cmComputation   *s = (cmComputation *)S;

  if ((g->nodesMax > 0) && (s->sFound == false))
    g->doSearchBFS(s, t);

  if ((g->depthMax > 0) && (s->sFound == false))
    g->doSearchDFS(s, t);

  if ((g->pathsMax > 0) && (s->sFound == false))
    g->doSearchRFS(s, t);
}



void *
cmReader(void *G) {
  cmGlobalData    *g = (cmGlobalData  *)G;
  cmComputation   *s = NULL;

  for (; g->curFragIID < g->numFrags; g->curFragIID++) {
    if (g->fi[g->curFragIID].doSearch == false)
      continue;

    if (g->fi[g->curFragIID].mateIID == 0)
      continue;

    if (g->fi[g->curFragIID].mateIID < g->curFragIID)
      continue;

    s = new cmComputation(g->curFragIID, g->fi[g->curFragIID].mateIID, g->innie);

    g->curFragIID++;

    break;
  }

  return(s);
}



void
cmWriter(void *G, void *S) {
  cmGlobalData    *g = (cmGlobalData  *)G;
  cmComputation   *c = (cmComputation *)S;

  if (c->sFound == true)
    fprintf(g->resultsFile, "Path from %d/%s to %d/%s found at depth %d of length %d.\n",
            c->fragIID, (c->frag5p3 == true) ? "5'3'" : "3'5'", 
            c->mateIID, (c->mate5p3 == true) ? "5'3'" : "3'5'",
            c->sPos,
            c->sLen);
  else
    fprintf(g->resultsFile, "Path from %d/%s NOT FOUND (%s).\n",
            c->fragIID,
            (c->frag5p3 == true) ? "5'3'" : "3'5'",
            (c->sLimited) ? "limited" : "exhausted");

  delete c;
}




int
main(int argc, char **argv) {
  char      *gkpStorePath      = NULL;
  char      *ovlStorePath      = NULL;
  char      *resultsPath       = NULL;

  uint32     distMin           = 0;
  uint32     distMax           = 0;
  bool       innie             = false;  //  require mates to be innie

  uint32     nodesMax          = 0;
  uint32     depthMax          = 0;
  uint32     pathsMax          = 0;

  bool       searchLibs        = false;
  uint32     searchLib[1024]   = {0};

  bool       backboneLibs      = false;
  uint32     backboneLib[1024] = {0};

  uint32     numThreads        = 4;

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

      char *a = argv[++arg];
      char *b = strchr(a, '-');

      b = (b) ? b+1 : a;

      uint32 bgn = atoi(a);
      uint32 end = atoi(b);

      for (uint32 i=bgn; i<=end; i++) {
        fprintf(stderr, "Search only in library %d\n", i);
        searchLib[i] = 1;
      }
      
    } else if (strcmp(argv[arg], "-bl") == 0) {
      backboneLibs = true;

      char *a = argv[++arg];
      char *b = strchr(a, '-');

      b = (b) ? b+1 : a;

      uint32 bgn = atoi(a);
      uint32 end = atoi(b);

      for (uint32 i=bgn; i<=end; i++) {
        fprintf(stderr, "Search using only library %d\n", i);
        backboneLib[i] = 1;
      }

    } else if (strcmp(argv[arg], "-min") == 0) {
      distMin = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      distMax = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-innie") == 0) {
      innie = true;

    } else if (strcmp(argv[arg], "-outtie") == 0) {
      innie = false;

    } else if (strcmp(argv[arg], "-bfs") == 0) {
      nodesMax = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-dfs") == 0) {
      depthMax = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-rfs") == 0) {
      pathsMax = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {
      numThreads = atoi(argv[++arg]);

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
    fprintf(stderr, "  -bfs N       Use 'breadth-first search'; search at most N fragments\n");
    fprintf(stderr, "  -dfs N       Use 'depth-first search'; search to depth N overlaps\n");
    fprintf(stderr, "  -rfs N       Use 'random-path search'; search at most N paths\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  cmGlobalData  *g = new cmGlobalData(resultsPath, distMin, distMax, innie, nodesMax, depthMax, pathsMax);

  g->loadFragments(gkpStorePath, searchLibs, searchLib, backboneLibs, backboneLib);
  g->loadOverlaps(ovlStorePath);

  sweatShop *ss = new sweatShop(cmReader, cmWorker, cmWriter);

  ss->setLoaderQueueSize(1048576);
  ss->setWriterQueueSize(65536);

  ss->setNumberOfWorkers(numThreads);

  for (u32bit w=0; w<numThreads; w++)
    ss->setThreadData(w, new cmThreadData());  //  these leak

  ss->run(g, true);

  delete g;
}

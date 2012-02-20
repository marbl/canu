
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

static const char *rcsid = "$Id: classifyMates-globalData.C,v 1.14 2012-02-20 01:54:39 brianwalenz Exp $";

#include "AS_global.h"

#include "classifyMates.H"
#include "classifyMates-globalData.H"

#include "memoryMappedFile.H"

#include <set>
using namespace std;



//  Overlap storage.  Each fragment has an entry in bbPos & tgPos that points to actual data in an
//  array.  There are bbLen & tgPos overlaps at this location.  The array is NOT contiguous.


cmGlobalData::cmGlobalData(char    *resultsName_,
                           uint32   distMin_,
                           uint32   distMax_,
                           bool     innie_,
                           uint32   nodesMax_,
                           uint32   depthMax_,
                           uint32   pathsMax_,
                           uint64   memoryLimit_) {

  strcpy(resultsPrefix, resultsName_);

  distMin      = distMin_;
  distMax      = distMax_;
  innie        = innie_;

  nodesMax     = nodesMax_;
  depthMax     = depthMax_;
  pathsMax     = pathsMax_;

  if      (nodesMax > 0)
    resultOutput = new classifyMatesResultFile(resultsName_, 'b', nodesMax, innie, distMin, distMax);
  else if (depthMax > 0)
    resultOutput = new classifyMatesResultFile(resultsName_, 'd', depthMax, innie, distMin, distMax);
  else if (pathsMax > 0)
    resultOutput = new classifyMatesResultFile(resultsName_, 'r', pathsMax, innie, distMin, distMax);
  else
    fprintf(stderr, "ERROR: no algorithm specified.\n"), exit(1);

  cmDat        = NULL;
  cmOvl        = NULL;
  cmInv        = NULL;

  numFrags     = 0;
  numLibs      = 0;
  fi           = 0L;

  curFragIID   = 0;

  isBB         = NULL;
  isSS         = NULL;

  bbPos        = 0L;
  bbLen        = 0L;

  tgPos        = 0L;
  tgLen        = 0L;

  gtPos        = 0L;
  gtLen        = 0L;

  memoryLimit  = memoryLimit_;
  memoryUsed   = 0;

  oiStorageMax = 0;
  oiStorageLen = 0;
  oiStorageArr = 0L;
  oiStorage    = 0L;
};




cmGlobalData::~cmGlobalData() {

  delete resultOutput;

  for (uint32 i=0; i<oiStorageMax; i++)
    delete [] oiStorageArr[i];

  if (cmDat) {
    fi    = NULL;
    bbLen = NULL;
    tgLen = NULL;
    gtLen = NULL;
  }

  delete cmDat;
  delete cmOvl;
  delete cmInv;

  delete [] isBB;
  delete [] isSS;

  delete [] bbPos;
  delete [] bbLen;

  delete [] tgPos;
  delete [] tgLen;

  delete [] gtPos;
  delete [] gtLen;

  delete [] oiStorageArr;

  delete [] fi;
};



bool
cmGlobalData::load(set<AS_IID>  &searchLibs,
                   set<AS_IID>  &backboneLibs,
                   double        maxErrorFraction) {

  char  cacheName[FILENAME_MAX];
  sprintf(cacheName, "%s.cmdat", resultsPrefix);

  if (AS_UTL_fileExists(cacheName, FALSE, FALSE) == false)
    return(false);
    
  fprintf(stderr, "LOADING FRAGMENTS...(from cache)\n");

  errno = 0;
  FILE *F = fopen(cacheName, "r");
  if (errno)
    fprintf(stderr, "Failed to open cache file '%s' for reading: %s\n", cacheName, strerror(errno)), exit(1);

  AS_UTL_safeRead(F, &numFrags,   "numFrags",   sizeof(uint32), 1);
  AS_UTL_safeRead(F, &numLibs,    "numLibs",    sizeof(uint32), 1);

  isBB = new uint32 [numLibs + 1];
  isSS = new uint32 [numLibs + 1];

  AS_UTL_safeRead(F,  isSS,       "isSS",       sizeof(uint32), numLibs + 1);
  AS_UTL_safeRead(F,  isBB,       "isBB",       sizeof(uint32), numLibs + 1);

  fclose(F);

  fprintf(stderr, "LOADING FRAGMENTS...%u fragments loaded.\n", numFrags);

  //  Check for a change in libraries.

  uint32  errs = 0;
  uint32  fixs = 0;

  fprintf(stderr, "CHECKING FRAGMENTS...\n");

  for (uint32 i=0; i<=numLibs; i++) {
    if (isBB[i] != ((backboneLibs.size() == 0) || (backboneLibs.count(i) > 0))) {
      if (isBB[i] == false) {
        //  Not in data file, but requested on command line, uh oh!
        fprintf(stderr, "ERROR: Library "F_U32" is not in the backbone data.\n", i);
        errs++;
      } else {
        //  In the data file, but not requested in the compute.  Fix it!
        fprintf(stderr, "Library "F_U32" will be excluded from the backbone.\n", i);
        fixs++;
        isBB[i] = false;
      }
    }

    if (isSS[i] != ((searchLibs.size() == 0) || (searchLibs.count(i)   > 0))) {
      if (isSS[i] == false) {
        fprintf(stderr, "ERROR: Library "F_U32" is not in the search data.\n", i);
        errs++;
      } else {
        fprintf(stderr, "Library "F_U32" will be excluded from the search.\n", i);
        fixs++;
        isSS[i] = false;
      }
    }
  }

  if (errs) {
    fprintf(stderr, "ERROR: libraries changed; cannot use these overlaps.\n");
    exit(1);
  }

  //  No easy way around this.  Allocate space for the pointers, but not the lengths.

  fprintf(stderr, "ALLOCATING...\n");

  allocateOverlapPointers();

  delete [] bbLen;
  delete [] tgLen;
  delete [] gtLen;

  //  Map in the fragment data...

  fprintf(stderr, "MAPPING FILES...\n");

  sprintf(cacheName, "%s.cmdat", resultsPrefix);
  cmDat = new memoryMappedFile(cacheName, 2 * sizeof(uint32) + 2 * (numLibs + 1) * sizeof(uint32));

  fi    = (fragmentInfo *)cmDat->get((numFrags + 1) * sizeof(fragmentInfo));
  bbLen = (uint32       *)cmDat->get((numFrags + 1) * sizeof(uint32));
  tgLen = (uint32       *)cmDat->get((numFrags + 1) * sizeof(uint32));
  gtLen = (uint32       *)cmDat->get((numFrags + 1) * sizeof(uint32));

  //  Update the fragment flags

  if (fixs) {
    fprintf(stderr, "UPDATING FRAGMENT FLAGS...\n");

    for (uint32 fid=0; fid<=numFrags; fid++) {
      uint32   lib = fi[fid].libIID;

      fi[fid].isBackbone = isBB[lib];
      fi[fid].doSearch   = isSS[lib];
    }
  }

  //  ...and the overlap data

  sprintf(cacheName, "%s.cmovl", resultsPrefix);
  cmOvl = new memoryMappedFile(cacheName);

  overlapInfo *oiRoot = (overlapInfo *)cmOvl->get(0, 0);
  overlapInfo *oiNext = oiRoot;

  sprintf(cacheName, "%s.cminv", resultsPrefix);
  cmInv = new memoryMappedFile(cacheName);

  overlapInfo *gtRoot = (overlapInfo *)cmInv->get(0);
  overlapInfo *gtNext = gtRoot;



  //  Then rebuild the pointers to overlaps

  fprintf(stderr, "LOADING OVERLAPS...(from cache)\n");

  uint64  numBB = 0;
  uint64  numTG = 0;
  uint64  numGT = 0;

  for (uint32 ii=0; ii<numFrags+1; ii++) {
    bbPos[ii] = NULL;
    tgPos[ii] = NULL;
    gtPos[ii] = NULL;

    if (bbLen[ii] > 0) {
      bbPos[ii] = oiNext;
      oiNext += bbLen[ii];
      numBB  += bbLen[ii];
    }

    if (tgLen[ii] > 0) {
      tgPos[ii] = oiNext;
      oiNext += tgLen[ii];
      numTG  += tgLen[ii];
    }

    if (gtLen[ii] > 0) {
      gtPos[ii] = gtNext;
      gtNext += gtLen[ii];
      numGT  += gtLen[ii];
    }

    if ((ii % 100000000) == 0)
      fprintf(stderr, ".."F_U32"\n", ii);
  }

  fprintf(stderr, "LOADING OVERLAPS...found "F_U64" BB overlaps and "F_U64" TG overlaps and "F_U64" GT overlaps.\n",
          numBB, numTG, numGT);



#if 0
  if (errs) {
    fprintf(stderr, "LOADING FRAGMENTS...backbone or search libraries changed; rebuilding map.\n");

    for (uint32 i=0; i<=numLibs; i++) {
      isSS[i] = ((searchLibs.size() == 0)   || (searchLibs.count(i)   > 0));
      isBB[i] = ((backboneLibs.size() == 0) || (backboneLibs.count(i) > 0));
    }
      
    for (uint32 fid=0; fid<=numFrags; fid++) {
      uint32   lib = fi[fid].libIID;

      fi[fid].isBackbone  = isBB[lib];
      fi[fid].doSearch    = isSS[lib];
    }
  }

  //  Useful diagnostic output (report library names used) but we don't have gkpStore open
  for (uint32 fid=0; fid<=numFrags; fid++) {
    uint32   lib = fi[fid].libIID;

    if ((fid > 0) && (fi[fid-1].isBackbone == 0) && (fi[fid].isBackbone == 1))
      fprintf(stderr, "  frag "F_U32" start of lib %s is now backbone.\n", fid, gkpStore->gkStore_getLibrary(lib)->libraryName);
    if ((fid > 0) && (fi[fid-1].isBackbone == 1) && (fi[fid].isBackbone == 0))
      fprintf(stderr, "  frag "F_U32" end of lib %s is last backbone.\n", fid-1, gkpStore->gkStore_getLibrary(lib)->libraryName);

    if ((fid > 0) && (fi[fid-1].doSearch == 0) && (fi[fid].doSearch == 1))
      fprintf(stderr, "  frag "F_U32" start of lib %s is searchable.\n", fid, gkpStore->gkStore_getLibrary(lib)->libraryName);
    if ((fid > 0) && (fi[fid-1].doSearch == 1) && (fi[fid].doSearch == 0))
      fprintf(stderr, "  frag "F_U32" end of lib %s is last searchable.\n", fid-1, gkpStore->gkStore_getLibrary(lib)->libraryName);
  }
#endif

  return(true);
}


void
cmGlobalData::save(void) {
  char  cacheName[FILENAME_MAX];
  sprintf(cacheName, "%s.cmdat", resultsPrefix);

  //  Overlap data is already saved (done inline when loading)

  errno = 0;
  FILE *F = fopen(cacheName, "w");
  if (errno)
    fprintf(stderr, "Failed to open cache file '%s' for writing: %s\n", cacheName, strerror(errno)), exit(1);

  AS_UTL_safeWrite(F, &numFrags,   "numFrags",   sizeof(uint32),       1);
  AS_UTL_safeWrite(F, &numLibs,    "numLibs",    sizeof(uint32),       1);
  AS_UTL_safeWrite(F,  isSS,       "isSS",       sizeof(uint32),       numLibs + 1);
  AS_UTL_safeWrite(F,  isBB,       "isBB",       sizeof(uint32),       numLibs + 1);
  AS_UTL_safeWrite(F,  fi,         "fi",         sizeof(fragmentInfo), numFrags + 1);
  AS_UTL_safeWrite(F,  bbLen,      "bbLen",      sizeof(uint32),       numFrags + 1);
  AS_UTL_safeWrite(F,  tgLen,      "tgLen",      sizeof(uint32),       numFrags + 1);
  AS_UTL_safeWrite(F,  gtLen,      "gtLen",      sizeof(uint32),       numFrags + 1);

  fclose(F);
}




void
cmGlobalData::buildBBSSmap(gkStore     *gkpStore,
                           set<AS_IID> &searchLibs,
                           set<AS_IID> &backboneLibs) {

  if (isBB == NULL) {
    isBB = new uint32 [numLibs + 1];
    isSS = new uint32 [numLibs + 1];
  }

  for (uint32 i=0; i<=numLibs; i++) {
    isSS[i] = ((searchLibs.size() == 0)   || (searchLibs.count(i)   > 0));
    isBB[i] = ((backboneLibs.size() == 0) || (backboneLibs.count(i) > 0));

    if (isSS[i] && isBB[i])
      fprintf(stderr, "library '%s' will be used as evidence AND classified.\n", gkpStore->gkStore_getLibrary(i)->libraryName);
    else if (isSS[i])
      fprintf(stderr, "library '%s' will be classified.\n", gkpStore->gkStore_getLibrary(i)->libraryName);
    else if (isBB[i])
      fprintf(stderr, "library '%s' will be used as evidence.\n", gkpStore->gkStore_getLibrary(i)->libraryName);
  }
}




void
cmGlobalData::loadFragments(char        *gkpStoreName,
                            set<AS_IID> &searchLibs,
                            set<AS_IID> &backboneLibs) {

  fprintf(stderr, "LOADING FRAGMENTS...\n");

  gkStore          *gkpStore  = new gkStore(gkpStoreName, FALSE, FALSE);
  gkStream         *gkpStream = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);
  gkFragment        frg;

  numFrags = gkpStore->gkStore_getNumFragments();
  numLibs  = gkpStore->gkStore_getNumLibraries();

  buildBBSSmap(gkpStore, searchLibs, backboneLibs);

  fi = new fragmentInfo [numFrags + 1];
  memset(fi, 0, sizeof(fragmentInfo) * (numFrags + 1));

  //  CHANGE CLEAR RANGE FOR PRODUCTION!  We need some way of asking if
  //  a clear range exists.  For now, we'll use a not-so-good method.

  while (gkpStream->next(&frg)) {
    uint32  fid = frg.gkFragment_getReadIID();
    uint32  lib = frg.gkFragment_getLibraryIID();

    uint32  clobt = frg.gkFragment_getClearRegionEnd(AS_READ_CLEAR_OBTCHIMERA);
    uint32  cllat = frg.gkFragment_getClearRegionEnd(AS_READ_CLEAR_LATEST);
    uint32  cllen = 0;

    if ((clobt == 0) && (cllat != 0))
      //  OBT doesn't exist, but latest does.
      cllen = frg.gkFragment_getClearRegionLength(AS_READ_CLEAR_LATEST);
    else
      //  OBT exists.
      cllen = frg.gkFragment_getClearRegionLength(AS_READ_CLEAR_OBTCHIMERA);

    fi[fid].libIID      = lib;
    fi[fid].contained   = 0;  //
    fi[fid].end5covered = 0;  //  Set when loading overlaps
    fi[fid].end3covered = 0;  //
    fi[fid].isBackbone  = isBB[lib];
    fi[fid].doSearch    = isSS[lib];
    fi[fid].clearLength = cllen;
    fi[fid].mateIID     = frg.gkFragment_getMateIID();

    assert(fi[fid].libIID      == lib);
    assert(fi[fid].clearLength == cllen);
    assert(fi[fid].mateIID     == frg.gkFragment_getMateIID());

    if ((fid % 10000000) == 0)
      fprintf(stderr, "LOADING FRAGMENTS...at IID "F_U32".\r",
              fid);
  }

  delete gkpStream;
  delete gkpStore;

  fprintf(stderr, "LOADING FRAGMENTS...%u fragments loaded.\n", numFrags);
}



void
cmGlobalData::allocateOverlapPointers(void) {

  memoryUsed += 3 * (numFrags + 1) * sizeof(overlapInfo *);
  memoryUsed += 3 * (numFrags + 1) * sizeof(uint32);

  if (memoryUsed > memoryLimit)
    fprintf(stderr, "\nMemory limit reached, aborting.\n"), exit(1);

  if (bbPos == NULL) {
    bbPos = new overlapInfo * [numFrags + 1];  //  Pointer to start of overlaps for this frag
    bbLen = new uint32        [numFrags + 1];  //  Number of overlaps for this frag
  }

  if (tgPos == NULL) {
    tgPos = new overlapInfo * [numFrags + 1];
    tgLen = new uint32        [numFrags + 1];
  }

  if (gtPos == NULL) {
    gtPos = new overlapInfo * [numFrags + 1];
    gtLen = new uint32        [numFrags + 1];
  }
}



void
cmGlobalData::allocateOverlapStorage(void) {

  if (oiStorageMax > 0)
    return;

  memoryUsed += (oiStorageMax) * sizeof(overlapInfo *);
  if (memoryUsed > memoryLimit)
    fprintf(stderr, "\nMemory limit reached, aborting.\n"), exit(1);

  oiStorageMax = 1024;                              //  Number blocks of memory
  oiStorageLen = 1;
  oiStorageArr = new overlapInfo * [oiStorageMax];  //  List of allocated memory
  oiStorage    = 0L;

  memset(oiStorageArr, 0, sizeof(overlapInfo *) * oiStorageMax);
}


void
cmGlobalData::resetOverlapStoreRange(OverlapStore *ovlStore,
                                     AS_IID curFragIID, AS_IID &minFragIID, AS_IID &maxFragIID) {
  bool   changed = false;

  minFragIID = curFragIID + 1;

  while ((minFragIID < maxFragIID) &&
         (fi[minFragIID].isBackbone == false) &&
         (fi[minFragIID].doSearch   == false)) {
    minFragIID++;
    changed = true;
  }

  while ((minFragIID < maxFragIID) &&
         (fi[maxFragIID-1].isBackbone == false) &&
         (fi[maxFragIID-1].doSearch   == false)) {
    maxFragIID--;
    changed = true;
  }

  if ((minFragIID < maxFragIID) &&
      (changed)) {
    fprintf(stderr, "LOADING OVERLAPS...for fragments %u to %u\n", minFragIID, maxFragIID);
    AS_OVS_setRangeOverlapStore(ovlStore, minFragIID, maxFragIID);
  }
}



void
cmGlobalData::loadOverlaps(char  *ovlStoreName,
                           double maxErrorFraction) {

  uint64            maxError    = AS_OVS_encodeQuality(maxErrorFraction);

  uint32            ovlMax      = 1048576;
  uint32            ovlLen      = 1;
  OVSoverlap       *ovl         = new OVSoverlap [ovlMax];
  bool             *ovlBB       = new bool       [ovlMax];
  bool             *ovlTG       = new bool       [ovlMax];

  allocateOverlapPointers();
  allocateOverlapStorage();

  memset(bbPos, 0, sizeof(overlapInfo *) * (numFrags + 1));
  memset(bbLen, 0, sizeof(uint32)        * (numFrags + 1));

  memset(tgPos, 0, sizeof(overlapInfo *) * (numFrags + 1));
  memset(tgLen, 0, sizeof(uint32)        * (numFrags + 1));

  memset(gtPos, 0, sizeof(overlapInfo *) * (numFrags + 1));
  memset(gtLen, 0, sizeof(uint32)        * (numFrags + 1));

  uint32            oiStorageBS  = 32 * 1024 * 1024;  //  Fixed block size (32 == 256mb, 128 = 1gb)
  uint32            oiStoragePos = 0;                  //  Position in the current block

  char  cacheName[FILENAME_MAX];
  sprintf(cacheName, "%s.cmovl", resultsPrefix);

  errno = 0;
  FILE             *oiFile       = fopen(cacheName, "w");
  if (errno)
    fprintf(stderr, "Failed to open overlap storage file '%s' for writing: %s\n", cacheName, strerror(errno)), exit(1);

  memoryUsed += (oiStorageBS) * sizeof(overlapInfo);
  if (memoryUsed > memoryLimit)
    fprintf(stderr, "\nMemory limit reached, aborting.\n"), exit(1);

  oiStorage = oiStorageArr[0] = new overlapInfo [oiStorageBS];

  uint64            numFG = 0;  //  Number of fragments processed

  uint64            numTT = 0;  //  Number of overlaps read from disk
  uint64            numUU = 0;  //  Number of overlaps not backbone and not target
  uint64            numBB = 0;  //  Number of BB overlaps loaded
  uint64            numTG = 0;  //  Number of TG overlaps loaded
  uint64            numDD = 0;  //  Number of overlaps discarded

  saveDistance     *bbDist  = new saveDistance(0);
  saveDistance     *tgDist  = new saveDistance(2);

  OverlapStore     *ovlStore  = AS_OVS_openOverlapStore(ovlStoreName);

  AS_IID            minFragIID = 1;
  AS_IID            maxFragIID = numFrags + 1;

  resetOverlapStoreRange(ovlStore, 1, minFragIID, maxFragIID);

  ovlLen = AS_OVS_readOverlapsFromStore(ovlStore, ovl, ovlMax, AS_OVS_TYPE_OVL);

  while (ovlLen > 0) {
    AS_IID  curFragIID = ovl[0].a_iid;

    //  Save the overlap if both fragments are in the backbone library, or if one fragment is in the
    //  backbone library and the other is in the test library.

    uint32  c=0;

    for (uint32 i=0; i<ovlLen; i++) {
      assert(curFragIID == ovl[i].a_iid);

      if ((ovl[i].dat.ovl.orig_erate    > maxError) ||
          ((fi[ovl[i].a_iid].isBackbone == false) && (fi[ovl[i].a_iid].doSearch   == false)) ||   //  Don't care about a.
          ((fi[ovl[i].b_iid].isBackbone == false) && (fi[ovl[i].b_iid].doSearch   == false)) ||   //  Don't care about b.
          ((fi[ovl[i].a_iid].isBackbone == false) &&
           (fi[ovl[i].a_iid].doSearch   == true)  &&
           (fi[ovl[i].b_iid].isBackbone == false) &&
           (fi[ovl[i].b_iid].doSearch   == true))) {    //  Don't care about target-target overlaps.
        numUU++;
        continue;
      }

      assert(((fi[ovl[i].a_iid].isBackbone == true) && (fi[ovl[i].b_iid].isBackbone == true)) ||
             ((fi[ovl[i].a_iid].isBackbone == true) && (fi[ovl[i].b_iid].doSearch   == true)) ||
             ((fi[ovl[i].a_iid].doSearch   == true) && (fi[ovl[i].b_iid].isBackbone == true)));

      if (c != i)
        ovl[c] = ovl[i];
      c++;
    }

    numTT  += ovlLen;
    ovlLen  = c;

    //  Figure out where to store this overlap (backbone or target) and if it is long enough to save.

    bbDist->compute(fi, ovl, ovlLen);
    tgDist->compute(fi, ovl, ovlLen);


    for (uint32 i=0; i<ovlLen; i++) {
      int32  ah = ovl[i].dat.ovl.a_hang;
      int32  bh = ovl[i].dat.ovl.b_hang;
      int32  fa = fi[ovl[i].a_iid].clearLength;
      int32  fb = fi[ovl[i].b_iid].clearLength;

      bool   cn = AS_OVS_overlapAIsContained(ovl[i]);
      bool   d5 = AS_OVS_overlapAEndIs5prime(ovl[i]);
      bool   d3 = AS_OVS_overlapAEndIs3prime(ovl[i]);
      bool   cr = AS_OVS_overlapAIsContainer(ovl[i]);

      ovlBB[i] = false;
      ovlTG[i] = false;

      //  Overlaps from a backbone A-fragment to a search B-fragment are kept if a is the container
      //  of b, or the dovetail overlap is long.  These are used to terminate the search.
      //  Documentation calls these "TG (target)"
      //
      if ((fi[ovl[i].a_iid].isBackbone == true) && (fi[ovl[i].b_iid].doSearch   == true))
        if ((cr) ||
            ((d5) && (tgDist->doveDist5 <= fb - -ah)) ||
            ((d3) && (tgDist->doveDist3 <= fb -  bh)))
          ovlTG[i] = true;

      //  Overlaps from backbone to backbone are kept if they are dovetail and long.
      //  Documentation calls these "BB (backbone)"
      //
      if ((fi[ovl[i].a_iid].isBackbone == true) && (fi[ovl[i].b_iid].isBackbone == true))
        if (((d5) && (bbDist->doveDist5 <= fb - -ah)) ||
            ((d3) && (bbDist->doveDist3 <= fb -  bh)))
          ovlBB[i] = true;

      //  Overlaps from a search A-fragment to a backbone B-fragment are kept if a is contained, or
      //  the dovetail overlap is long.  These are used to initiate the search.
      //  Documentation calls these "TB (backbone)"
      //
      if ((fi[ovl[i].a_iid].doSearch   == true) && (fi[ovl[i].b_iid].isBackbone == true))
        if ((cn) ||
            ((d5) && (bbDist->doveDist5 <= fb - -ah)) ||
            ((d3) && (bbDist->doveDist3 <= fb -  bh)))
          ovlBB[i] = true;


      //  Save or discard the overlap.  A fragment can be both backbone and search, and thus we are
      //  allowed to save the overlap twice.

      if ((ovlBB[i] == false) && (ovlTG[i] == false))
        ovl[i].a_iid = 0;

      if (ovlBB[i] == true)
        bbLen[ovl[i].a_iid]++;

      if (ovlTG[i] == true)
        tgLen[ovl[i].a_iid]++;


      //  Remember which end this overlap covered.  If contained, both ends are covered.  Note that if
      //  we don't care about this overlap, the a_iid of the overlap is reset to zero, which isn't a
      //  fragment (so we don't actually remember anything).
      //
      fi[ovl[i].a_iid].contained   |= cn;
      fi[ovl[i].a_iid].end5covered |= d5;
      fi[ovl[i].a_iid].end3covered |= d3;

    }  //  Over all overlaps

    //  If not enough space in oiStorage, get more space.

    if (oiStoragePos + bbLen[curFragIID] + tgLen[curFragIID] > oiStorageBS) {
      if ((oiFile) && (oiStoragePos > 0))
        AS_UTL_safeWrite(oiFile, oiStorage, "oiStorage", sizeof(overlapInfo), oiStoragePos);

      oiStorage = oiStorageArr[oiStorageLen] = new overlapInfo [oiStorageBS];
      oiStorageLen++;
      oiStoragePos = 0;

      memoryUsed += (oiStorageBS) * sizeof(overlapInfo);
      if (memoryUsed > memoryLimit)
        fprintf(stderr, "\nMemory limit reached, aborting.\n"), exit(1);
    }

    //  Add the overlaps to storage.

    bbPos[curFragIID] = oiStorage + oiStoragePos;
    tgPos[curFragIID] = oiStorage + oiStoragePos + bbLen[curFragIID];

    bbLen[curFragIID] = 0;
    tgLen[curFragIID] = 0;

    for (uint32 i=0; i<ovlLen; i++) {
      if ((ovlBB[i] == false) && (ovlTG[i] == false)) {
        numDD++;
      }

      if (ovlBB[i] == true) {
        bbPos[curFragIID][bbLen[curFragIID]++] = overlapInfo(ovl[i]);
        numBB++;
      }

      if (ovlTG[i] == true) {
        tgPos[curFragIID][tgLen[curFragIID]++] = overlapInfo(ovl[i]);
        numTG++;
      }
    }

    oiStoragePos += bbLen[curFragIID] + tgLen[curFragIID];

    assert(oiStoragePos <= oiStorageBS);

    if ((++numFG % 3000000) == 0)
      fprintf(stderr, "LOADING OVERLAPS...at IID "F_U32" (%06.2f%%): BB "F_U64" (%06.2f%%) TG "F_U64" (%06.2f%%) DD "F_U64" (%06.2f%%).\n",
              curFragIID, 100.0 * curFragIID / numFrags,
              numBB,      100.0 * numBB      / numTT,
              numTG,      100.0 * numTG      / numTT,
              numDD,      100.0 * numDD      / numTT);


    resetOverlapStoreRange(ovlStore, curFragIID, minFragIID, maxFragIID);

    ovlLen = AS_OVS_readOverlapsFromStore(ovlStore, ovl, ovlMax, AS_OVS_TYPE_OVL);
  }

  AS_OVS_closeOverlapStore(ovlStore);

  delete [] ovl;
  delete [] ovlBB;
  delete [] ovlTG;

  delete bbDist;
  delete tgDist;

  if ((oiFile) && (oiStoragePos > 0))
    AS_UTL_safeWrite(oiFile, oiStorage, "oiStorage", sizeof(overlapInfo), oiStoragePos);

  if (oiFile)
    fclose(oiFile);

  fprintf(stderr, "LOADING OVERLAPS...at IID "F_U32" (%06.2f%%): BB "F_U64" (%06.2f%%) TG "F_U64" (%06.2f%%) DD "F_U64" (%06.2f%%).\n",
          numFrags, 100.0,
          numBB, 100.0 * numBB / numTT,
          numTG, 100.0 * numTG / numTT,
          numDD, 100.0 * numDD / numTT);
  fprintf(stderr, "LOADING OVERLAPS..."F_U64" overlaps loaded.\n", numBB + numTG);

  loadOverlaps_invert();
}




void
cmGlobalData::loadOverlaps_invert(void) {

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

  overlapInfo *oiStorage = oiStorageArr[oiStorageLen++] = new overlapInfo [numGTovl];
  uint32       oiStoragePos = 0;

  memoryUsed += (numGTovl) * sizeof(overlapInfo);
  if (memoryUsed > memoryLimit)
    fprintf(stderr, "\nMemory limit reached, aborting.\n"), exit(1);

  //  Set pointers to the space.  We'll recount as the overlaps are added.
  for (uint32 ii=0; ii<numFrags+1; ii++) {
    gtPos[ii]      = (gtLen[ii] == 0) ? NULL : oiStorage + oiStoragePos;  //  Point into the space
    oiStoragePos  += gtLen[ii];                                           //  Reserve the space
    gtLen[ii]      = 0;                                                   //  Show we have no overlaps loaded
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

  char  cacheName[FILENAME_MAX];
  sprintf(cacheName, "%s.cminv", resultsPrefix);

  errno = 0;
  FILE *F = fopen(cacheName, "w");
  if (errno)
    fprintf(stderr, "Failed to open overlap storage file '%s' for writing: %s\n", cacheName, strerror(errno)), exit(1);

  AS_UTL_safeWrite(F, oiStorage, "oiStorage", sizeof(overlapInfo), oiStoragePos);

  fclose(F);

  fprintf(stderr, "INVERTING OVERLAPS...."F_U64" frags found.\n", numGTfrg);

  fprintf(stderr, "LOADING OVERLAPS..."F_U64" GB used.\n", memoryUsed >> 30);
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
  } else

  //  Frag is forward, overlap is flipped, hang is positive
  if ((t->path[t->pathPos].p5p3 == true) && (novl->flipped == true)) {
    nlen = t->path[t->pathPos].pLen + novl->bhang;
    assert(n5p3 == false);
  } else

  //  Frag is reverse, overlap is same, hang is positive
  if ((t->path[t->pathPos].p5p3 == false) && (novl->flipped == false)) {
    nlen = t->path[t->pathPos].pLen + -novl->ahang;
    assert(n5p3 == false);
  } else

  //  Frag is reverse, overlap is flipped, hang is positive
  if ((t->path[t->pathPos].p5p3 == false) && (novl->flipped == true)) {
    nlen = t->path[t->pathPos].pLen + -novl->ahang;
    assert(n5p3 == true);
  }
}




bool
cmGlobalData::testSearch(cmComputation              *c,
                         cmThreadData               *t) {

  uint32  thisIID = t->path[t->pathPos].pIID;

  if (t->solutionSet->isSet(thisIID) == false)
    return(false);

  //  We are guaranteed by construction that this overlap is to the mate fragment.
  //  See doSearchBFS.

  overlapInfo  *pos = gtPos[c->mateIID];
  uint32        len = gtLen[c->mateIID];
  uint32        iii = 0;

  while ((iii<len) && (thisIID != pos[iii].iid))
    iii++;

  assert(iii < len);

  overlapInfo  *novl = pos + iii;
  uint32        niid = c->mateIID;
  bool          n5p3 = (novl->flipped) ? (!t->path[t->pathPos].p5p3) : (t->path[t->pathPos].p5p3);
  uint32        nlen = 0;

  if (n5p3 != c->mate5p3)
    return(false);

  computeNextPlacement(c, t, novl, niid, n5p3, nlen);

  if ((nlen < distMin) ||
      (nlen > distMax))
    return(false);

  c->result.classified = true;
  c->result.iteration  = t->searchIter;
  c->result.distance   = nlen;

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
      c->result.classified = true;
      c->result.iteration  = t->searchIter;
      c->result.distance   = nlen;

      return(true);
    }
  }
  return(false);
}



//  Returns true if the frag or the mate look like they are spurs.
bool
cmGlobalData::testSpur(cmComputation  *c, cmThreadData   *t) {

  if ((fi[c->fragIID].contained == false) &&
      ((fi[c->fragIID].end5covered == false) || (fi[c->fragIID].end3covered == false)))
    c->result.fragSpur = true;

  if ((fi[c->mateIID].contained == false) &&
      ((fi[c->mateIID].end5covered == false) || (fi[c->mateIID].end3covered == false)))
    c->result.mateSpur = true;

  return(c->result.fragSpur || c->result.mateSpur);
}



//  Returns true if the frag or mate look like they are chimeric.  This tests
//  that there is a pair of overlapping reads overlapping the read:
//
//      read         -----------------
//      backbone ---------------
//      backbone          --------------------
//
//  Or, as the expected case with long backbone reads, that the reads are contained.
//

uint32 *
reallocEdge(uint32 *edge, uint32 &max, uint32 len) {

  if (max > len)
    return(edge);

  delete [] edge;

  while (max < len)
    max *= 2;

  return(new uint32 [max]);
}


bool
cmGlobalData::testChimer(uint32 iid, cmThreadData *t) {
  OVSoverlap  ovl;
  uint32      containCount = 0;
  uint32      overlapCount = 0;

  t->edge5 = reallocEdge(t->edge5, t->edge5Max, bbLen[iid]);
  t->edge3 = reallocEdge(t->edge3, t->edge3Max, bbLen[iid]);

  t->edge5Len = 0;
  t->edge3Len = 0;

  for (uint32 test=0; test<bbLen[iid]; test++) {
    overlapInfo  *novl = bbPos[iid] + test;
    uint32        niid = novl->iid;

    //  Overlaps here are from target (A) to backbone (B).  The backbone ID will be stored in 'iid'
    //  of the overlap.

    assert(fi[iid].doSearch    == true);
    assert(fi[niid].isBackbone == true);

    assert(niid != iid);

    //  Figure out which end this overlap is on.  We have nice functions to do this for us...if the overlap is
    //  an OVSoverlap.

    ovl.dat.ovl.flipped = novl->flipped;
    ovl.dat.ovl.a_hang  = novl->ahang;
    ovl.dat.ovl.b_hang  = novl->bhang;

    if (AS_OVS_overlapAIsContained(ovl))
      //  Backbone contains Target, not chimer.
      containCount++;

    if (AS_OVS_overlapAEndIs5prime(ovl))
      t->edge5[t->edge5Len++] = niid;

    if (AS_OVS_overlapAEndIs3prime(ovl))
      t->edge3[t->edge3Len++] = niid;
  }

#if 0
  for (uint32 i=0; i<t->edge5Len || i<t->edge3Len; i++)
    fprintf(stderr, F_U32"\t"F_U32"\n",
            (i < t->edge5Len) ? t->edge5[i] : 0,
            (i < t->edge3Len) ? t->edge3[i] : 0);
#endif

  //  The hard part.  Find if any pair in (edge5, edge3) have an overlap in the backbone.  We'll
  //  iterate over all the 5' fragment overlaps, and see if one of those is in the 3' set.
  //
  //  We do NOT check if the overlap makes sense, just that it exists.

  set<uint32>  iid3;

  for (uint32 i=0; i<t->edge3Len; i++)
    iid3.insert(t->edge3[i]);

  for (uint32 i=0; i<t->edge5Len; i++) {
    uint32 tiid = t->edge5[i];

    for (uint32 test=0; test<bbLen[tiid]; test++) {
      overlapInfo  *novl = bbPos[tiid] + test;

      if (iid3.count(novl->iid) > 0)
        overlapCount++;
    }
  }

#if 0
  fprintf(stdout, "IID "F_U32" has "F_U32" overlaps, containCount "F_U32" and overlapCount "F_U32"\n",
          iid, bbLen[iid], containCount, overlapCount);
#endif

  return(containCount + overlapCount < 1);
}


bool
cmGlobalData::testChimers(cmComputation *c, cmThreadData *t) {

  if ((fi[c->fragIID].contained == true) &&
      (fi[c->fragIID].contained == true))
    return(false);

  if (testChimer(c->fragIID, t))
    c->result.fragChimer = true;

  if (testChimer(c->mateIID, t))
    c->result.mateChimer = true;

  return(c->result.fragChimer || c->result.mateChimer);
}

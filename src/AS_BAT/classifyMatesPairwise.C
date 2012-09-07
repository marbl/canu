
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

const char *mainid = "$Id: classifyMatesPairwise.C,v 1.2 2012-09-07 01:52:05 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_decodeRange.H"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"

#include "classifyMates.H"
#include "classifyMates-globalData.H"

#include <set>
using namespace std;


class frgInfo {
public:
  bool    isTrusted     : 1;   //  If set, this can be evidence
  bool    isOverlapping : 1;   //  If set, this can be evidence that is conflicting
  bool    doTesting     : 1;   //  If set, this pair needs to be tested
  bool    isGarbage     : 1;   //  If set, this pair was tested, and looks like PE

  uint64  libIID        : 12;
  uint64  clearLength   : 16;
  uint64  mateIID       : 32;
};

class libInfo {
public:
  bool    isTrusted     : 1;
  bool    doTesting     : 1;

  uint32  numTrusted;          //  Num pairs that are trusted in the input
  uint32  numOverlapping;      //  Num pairs found to be overlapping

  uint32  numTested;           //  Num pairs that we need to test in the input
  uint32  numGarbage;          //  Num pairs found to be PE

  int32   mean;
  int32   stddev;
};


vector<frgInfo>   finf;
vector<libInfo>   linf;


void
loadFragmentData(char         *gkpStoreName,
                 set<AS_IID>   trstLibs,
                 set<AS_IID>   testLibs) {
  gkStore          *gkpStore  = new gkStore(gkpStoreName, FALSE, FALSE);
  gkStream         *gkpStream = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);
  gkFragment        frg;

  uint32 numFrags   = gkpStore->gkStore_getNumFragments();
  uint32 numLibs    = gkpStore->gkStore_getNumLibraries();

  linf.resize(numLibs + 1);
  finf.resize(numFrags + 1);

  for (uint32 ii=1; ii<=numLibs; ii++) {
    linf[ii].isTrusted       = (trstLibs.count(ii) > 0);
    linf[ii].doTesting       = (testLibs.count(ii) > 0);

    linf[ii].numTrusted      = 0;
    linf[ii].numOverlapping  = 0;

    linf[ii].numTested       = 0;
    linf[ii].numGarbage      = 0;

    linf[ii].mean            = gkpStore->gkStore_getLibrary(ii)->mean;
    linf[ii].stddev          = gkpStore->gkStore_getLibrary(ii)->stddev;

    fprintf(stderr, "LIB %d TRUST %d TEST %d %d +- %d\n",
            ii,
            linf[ii].isTrusted,
            linf[ii].doTesting,
            linf[ii].mean, linf[ii].stddev);
  }

  uint32 numTrusted = 0;
  uint32 numTested  = 0;

  while (gkpStream->next(&frg)) {
    uint32  fid = frg.gkFragment_getReadIID();
    uint32  lid = frg.gkFragment_getLibraryIID();

    uint32  clobt = frg.gkFragment_getClearRegionEnd(AS_READ_CLEAR_OBTCHIMERA);
    uint32  cllat = frg.gkFragment_getClearRegionEnd(AS_READ_CLEAR_LATEST);
    uint32  cllen = 0;

    if ((clobt == 0) && (cllat != 0))
      //  OBT doesn't exist, but latest does.
      cllen = frg.gkFragment_getClearRegionLength(AS_READ_CLEAR_LATEST);
    else
      //  OBT exists.
      cllen = frg.gkFragment_getClearRegionLength(AS_READ_CLEAR_OBTCHIMERA);

    bool  isMated = (frg.gkFragment_getMateIID() > 0);

    finf[fid].isTrusted       = (linf[lid].isTrusted && isMated);
    finf[fid].isOverlapping   = false;
    finf[fid].doTesting       = (linf[lid].doTesting && isMated);
    finf[fid].isGarbage       = false;

    finf[fid].libIID          = lid;
    finf[fid].clearLength     = cllen;
    finf[fid].mateIID         = frg.gkFragment_getMateIID();

    if (isMated) {
      linf[lid].numTrusted     += linf[lid].isTrusted;
      linf[lid].numTested      += linf[lid].doTesting;

      numTrusted     += linf[lid].isTrusted;
      numTested      += linf[lid].doTesting;
    }

    assert(finf[fid].libIID      == lid);
    assert(finf[fid].clearLength == cllen);
    assert(finf[fid].mateIID     == frg.gkFragment_getMateIID());

    if ((fid % 10000000) == 0)
      fprintf(stderr, "LOADING FRAGMENTS...at IID "F_U32".\r",
              fid);
  }

  delete gkpStream;
  delete gkpStore;

  fprintf(stderr, "LOADING FRAGMENTS...%u fragments loaded (%u trusted; %u to test).\n",
          numFrags, numTrusted, numTested);
}


bool
testIfTruePEAreOverlapping(OVSoverlap  *avl,
                           uint32       avlLen,
                           OVSoverlap  *bvl,
                           uint32       bvlLen) {
  int32   fsa = 0;
  int32   fsb = 0;

  AS_IID  aiid = avl[0].a_iid;
  AS_IID  biid = bvl[0].a_iid;

  assert(finf[aiid].mateIID == biid);
  assert(finf[biid].mateIID == aiid);

  uint32  ai = 0;
  uint32  bi = 0;

  for (ai=0; ai<avlLen; ai++)
    if (avl[ai].b_iid == biid)
      break;

  for (bi=0; bi<bvlLen; bi++)
    if (bvl[bi].b_iid == aiid)
      break;

  if (ai < avlLen) {
    finf[aiid].isOverlapping = true;
    finf[biid].isOverlapping = true;

    fsa = finf[aiid].clearLength + avl[ai].dat.ovl.b_hang;
  }

  if (bi < bvlLen) {
    finf[aiid].isOverlapping = true;
    finf[biid].isOverlapping = true;

    fsb = finf[biid].clearLength + bvl[bi].dat.ovl.a_hang;
  }

  return(MAX(fsa, fsb));
}



int
main(int argc, char **argv) {
  char       *gkpStoreName      = NULL;
  char       *ovlStoreName      = NULL;
  char       *resultsName       = NULL;

  bool        beVerbose         = true;

  double      maxErrorFraction  = 0.045;

  uint32      distMin           = 0;
  uint32      distMax           = 0;
  bool        innie             = false;  //  require mates to be innie

  set<AS_IID> trstLibs;
  set<AS_IID> testLibs;


  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      maxErrorFraction = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      resultsName = argv[++arg];

    } else if (strcmp(argv[arg], "-trusted") == 0) {
      AS_UTL_decodeRange(argv[++arg], trstLibs);
      
    } else if (strcmp(argv[arg], "-test") == 0) {
      AS_UTL_decodeRange(argv[++arg], testLibs);


    }

    arg++;
  }
  if (gkpStoreName == NULL)
    err++;
  if (ovlStoreName == NULL)
    err++;
  if (resultsName == NULL)
    err++;
  if (trstLibs.size() == 0)
    err++;
  if (testLibs.size() == 0)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore -o outName -trusted X-Y -test A-B\n", argv[0]);

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-G) supplied.\n");
    if (ovlStoreName == NULL)
      fprintf(stderr, "ERROR: no ovlStore (-O) supplied.\n");
    if (resultsName == NULL)
      fprintf(stderr, "ERROR: no output name (-o) supplied.\n");
    if (trstLibs.size() == 0)
      fprintf(stderr, "ERROR: no trusted libraries (-tusted) supplied.\n");
    if (testLibs.size() == 0)
      fprintf(stderr, "ERROR: no test libraries (-test) supplied.\n");
    
    exit(1);
  }

  //  Load the mate pair map, and read lenghts.

  loadFragmentData(gkpStoreName, trstLibs, testLibs);


  char              fileName[FILENAME_MAX];

  sprintf(fileName, "%s.classifyMatesPairwise.gkpEdit", resultsName);
  FILE             *edtFile  = fopen(fileName, "w");

  sprintf(fileName, "%s.classifyMatesPairwise.log", resultsName);
  FILE             *logFile  = fopen(fileName, "w");

  OverlapStore     *ovlStore  = AS_OVS_openOverlapStore(ovlStoreName);

  uint32            avlMax      = 1048576;
  uint32            avlLen      = 0;
  OVSoverlap       *avl         = new OVSoverlap [avlMax];

  uint32            bvlMax      = 1048576;
  uint32            bvlLen      = 0;
  OVSoverlap       *bvl         = new OVSoverlap [bvlMax];

  uint32            tvlMax;
  uint32            tvlLen;
  OVSoverlap       *tvl;

  

  //  The store doesn't (efficiently) support the kind of iteration we want to do.  We want to get
  //  all overlaps for a mated pair.  If one of the reads has no overlaps, the store will return
  //  overlaps for the next id.


  AS_IID            minFragIID = finf.size();
  AS_IID            maxFragIID = 0;

  for (uint32 fi=1; fi<finf.size(); fi++) {
    if (finf[fi].doTesting == false)
      continue;

    if (fi < minFragIID)
      minFragIID = fi;
    if (maxFragIID < fi)
      maxFragIID = fi;
  }

  if (minFragIID > maxFragIID)
    fprintf(stderr, "ERROR: nothing to test.\n"), exit(1);

  fprintf(stderr, "Testing reads from "F_U32" to "F_U32".\n",
          minFragIID, maxFragIID);

  //AS_OVS_setRangeOverlapStore(ovlStore, minFragIID, maxFragIID);






  bool  stillMore = true;

  while (stillMore) {
    if (avlLen == 0)
      avlLen = AS_OVS_readOverlapsFromStore(ovlStore, avl, avlMax, AS_OVS_TYPE_OVL);

    if (bvlLen == 0)
      bvlLen = AS_OVS_readOverlapsFromStore(ovlStore, bvl, bvlMax, AS_OVS_TYPE_OVL);

    if ((avlLen == 0) && (bvlLen == 0))
      break;

    bool  AnoGood = false;

    if (finf[avl[0].a_iid].mateIID == 0)
      //  avl overlaps are not for a mated read.
      AnoGood = true;

    if (finf[avl[0].a_iid].isTrusted == false)
      //  avl overlaps are not for a trusted read.
      AnoGood = true;

    if (finf[avl[0].a_iid].mateIID != bvl[0].a_iid)
      //  avl overlaps are not mated to the bovl overlaps.
      AnoGood = true;

    //  If AnoGood, discard the avl overlaps, move the bvl overlaps to avl, and continue
    //  to read a new batch of bvl overlaps.
    if (AnoGood) {
      tvlMax = avlMax;   tvlLen = avlLen;   tvl = avl;
      avlMax = bvlMax;   avlLen = bvlLen;   avl = bvl;
      bvlMax = tvlMax;   bvlLen = 0;        bvl = tvl;
      continue;
    }

    //  We have overlaps.  They should be mate-consistent, and trusted.

    assert(finf[avl[0].a_iid].mateIID   == bvl[0].a_iid);
    assert(finf[bvl[0].a_iid].mateIID   == avl[0].a_iid);

    assert(finf[avl[0].a_iid].isTrusted == true);
    assert(finf[bvl[0].a_iid].isTrusted == true);

    //  Check if we're overlapping - if the other read in an overlap is the mate

    bool   overlappingPE = false;
    int32  fragSize      = testIfTruePEAreOverlapping(avl, avlLen, bvl, bvlLen);

    if (fragSize == 0) {
      //  Nope.
      fragSize = linf[finf[avl[0].a_iid].libIID].mean;
    } else {
      //  Yup.
      overlappingPE = true;
    }


    //  Examine avl and bvl for any 'testable' pairs.
    //  Both overlaps have Aiid from the PE pair.
    //  We need to find two overlaps such that the Biid is a testable pair.

    for (uint32 ai=0, si=0; ai<avlLen; ai++) {
      AS_IID  taiid = avl[ai].b_iid;

      if ((finf[taiid].mateIID == 0) ||       //  Overlapping read isn't mated, don't care about it.
          (finf[taiid].doTesting == false))   //  Overlapping read isn't marked for testing, don't care about it.
        continue;

      //  Starting position for the next loop.  Advance past the overlaps that are below the read
      //  we're looking for.

      while ((si < bvlLen) &&
             (bvl[si].b_iid < finf[taiid].mateIID))
        si++;

      if (si > 0)  //  Because the loop stops on the one after we want.
        si--;

      //  The ai overlap is to a testable read.  Search for the mate in the
      //  other overlap set.

      for (uint32 bi=si; bi<bvlLen; bi++) {
        AS_IID  tbiid = bvl[bi].b_iid;

        if (finf[taiid].mateIID < tbiid)        //  Overlapping read is after the one we're looking for.  It isn't here.
          break;

        if ((tbiid < finf[taiid].mateIID) ||    //  Overlapping read is before the one we're looking for.  Keep looking.
            (finf[tbiid].doTesting == false))   //  Overlapping read isn't marked for testing, don't care about it.
          continue;

        //  Overlapping read is the one we're looking for.  Test it!

        assert(finf[avl[ai].a_iid].mateIID   == bvl[bi].a_iid);  //  The trusted PE pair
        assert(finf[bvl[bi].a_iid].mateIID   == avl[ai].a_iid);

        assert(finf[avl[ai].b_iid].mateIID   == bvl[bi].b_iid);  //  The tested PE pair
        assert(finf[bvl[bi].b_iid].mateIID   == avl[ai].b_iid);

        //  Compute the expected size of the MP pair (do NOT add to fragSize; it is reused for the next pair in the same overlap set)

        int32  adiff = -avl[ai].dat.ovl.a_hang;
        int32  bdiff = -bvl[bi].dat.ovl.a_hang;

        //  MP pair is plausibly PE If
        //    both are inverted
        //    both are normal and the PE overlaps
        //

        bool plausiblyPE = false;

        if ((avl[ai].dat.ovl.flipped == true) && (avl[ai].dat.ovl.flipped == true))
          plausiblyPE = true;

        if ((avl[ai].dat.ovl.flipped == false) && (avl[ai].dat.ovl.flipped == false) && (overlappingPE == true))
          plausiblyPE = true;

        if (plausiblyPE == false) {
          fprintf(logFile, "A %8u-%8u %c %4ld-%4ld  B %8u-%8u %c %4ld-%4ld -- size %4d %c -- SHORT_MP\n",
                  avl[ai].a_iid, avl[ai].b_iid, avl[ai].dat.ovl.flipped ? 'I' : 'N', avl[ai].dat.ovl.a_hang, avl[ai].dat.ovl.b_hang,
                  bvl[bi].a_iid, bvl[bi].b_iid, bvl[bi].dat.ovl.flipped ? 'I' : 'N', bvl[bi].dat.ovl.a_hang, bvl[bi].dat.ovl.b_hang,
                  fragSize + adiff + bdiff, overlappingPE ? 'O' : '-');

        } else {
          fprintf(logFile, "A %8u-%8u %c %4ld-%4ld  B %8u-%8u %c %4ld-%4ld -- size %4d %c -- PLAUSIBLY_PE\n",
                  avl[ai].a_iid, avl[ai].b_iid, avl[ai].dat.ovl.flipped ? 'I' : 'N', avl[ai].dat.ovl.a_hang, avl[ai].dat.ovl.b_hang,
                  bvl[bi].a_iid, bvl[bi].b_iid, bvl[bi].dat.ovl.flipped ? 'I' : 'N', bvl[bi].dat.ovl.a_hang, bvl[bi].dat.ovl.b_hang,
                  fragSize + adiff + bdiff, overlappingPE ? 'O' : '-');

          fprintf(edtFile, "frg iid %u mateiid 0\n", avl[ai].b_iid);
          fprintf(edtFile, "frg iid %u mateiid 0\n", bvl[bi].b_iid);
        }
      }
    }

    //  Done with this pair of overlaps (OK, 'sets of overlaps', but 'pair of sets of overlaps' is cumbersome).

    avlLen = 0;
    bvlLen = 0;
  }

  fclose(edtFile);
  fclose(logFile);

  exit(0);
}

/**************************************************************************
 * Copyright (C) 2014
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

const char *mainid = "$Id:  $";

#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"

#include "memoryMappedFile.H"

#include "intervalList.H"

#include <map>
#include <set>
#include <list>
#include <vector>
#include <algorithm>

using namespace std;

#ifndef BROKEN_CLANG_OpenMP
#include <omp.h>
#endif

#define ERATE_TOLERANCE 0.03

uint32 blockSize = 1000;

class readErrorEstimate {
public:
  readErrorEstimate() {
    clrBgn = 0;
    clrEnd = 0;
    clrLen = 0;

    errorMeanS      = NULL;  //  Sum of error up to this point
    //errorStdDev     = NULL;
    //errorConfidence = NULL;
    errorMeanU      = NULL;  //  Updated error estimate for this point
  };

  ~readErrorEstimate() {
    delete [] errorMeanS;
    //delete [] errorStdDev;
    //delete [] errorConfidence;
    delete [] errorMeanU;
  };

  uint64     initialize(gkFragment &fr) {
    //uint32      lid = fr.gkFragment_getLibraryIID();
    //gkLibrary  *lb  = gkpStore->gkStore_getLibrary(lid);

    if (fr.gkFragment_getIsDeleted() == 1)
      return(sizeof(readErrorEstimate));

    clrBgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_LATEST);
    clrEnd = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_LATEST);
    clrLen = clrEnd - clrBgn;

    errorMeanS      = new uint32 [clrLen + 1];
    //errorStdDev     = new double [clrLen + 1];
    //errorConfidence = new double [clrLen + 1];
    errorMeanU      = new uint16 [clrLen + 1];

#if 0
    uint32  err40 = AS_OVS_encodeQuality(0.40);

    errorMeanS[0] = err40;
    errorMeanU[0] = 0;

    //  instead of this loop, we can build one full (0-max_read_len) array,
    //  then copy the start N elements to each errorMeanS

    for (uint32 ii=1; ii<clrLen; ii++) {
      errorMeanS[ii] = errorMeanS[ii-1] + err40;
      errorMeanU[ii] = 0;
    }
#endif

#if 0
    //memset(errorMeanS, 0xff, sizeof(uint32) * (readProfile[iid].clrLen + 1));
    //memset(errorMeanU, 0x00, sizeof(uint16) * (readProfile[iid].clrLen + 1));
#endif

    return(sizeof(uint16) * (clrLen + 1 + clrLen + 1) + sizeof(readErrorEstimate));
  };

  uint32     clrBgn;
  uint32     clrEnd;
  uint32     clrLen;

  uint32    *errorMeanS;
  //double    *errorStdDev;
  //double    *errorConfidence;
  uint16    *errorMeanU;
};


class ESToverlap {
public:
  uint32   a_iid      : 23;
  uint32   b_iid_hi   :  9;  //  32 bits

  uint32   b_iid_lo   : 14;
  int32    a_hang     : 17;  //  31 bits

  int32    b_hang     : 17;
  uint32   erate      : 12;
  uint32   flipped    : 1;
  uint32   discarded  : 1;   //  31 bits

  void     populate(OVSoverlap &ovl) {
    a_iid     =  ovl.a_iid;
    b_iid_hi  = (ovl.b_iid >> 14) & 0x000001ff;
    b_iid_lo  = (ovl.b_iid >>  0) & 0x00003fff;
    a_hang    =  ovl.dat.ovl.a_hang;
    b_hang    =  ovl.dat.ovl.b_hang;
    erate     =  ovl.dat.ovl.orig_erate;
    flipped   =  ovl.dat.ovl.flipped;
    discarded =  false;

    assert(ovl.a_iid              ==   a_iid);
    assert(ovl.b_iid              == ((b_iid_hi << 14) | (b_iid_lo)));
    assert(ovl.dat.ovl.a_hang     ==   a_hang);
    assert(ovl.dat.ovl.b_hang     ==   b_hang);
    assert(ovl.dat.ovl.orig_erate ==   erate);
    assert(ovl.dat.ovl.flipped    ==   flipped);
  };

  void     populateOVL(OVSoverlap &ovl) {
    ovl.a_iid              =  a_iid;
    ovl.b_iid              = (b_iid_hi << 14) | (b_iid_lo);

    ovl.dat.ovl.a_hang     = a_hang;
    ovl.dat.ovl.b_hang     = b_hang;

    ovl.dat.ovl.flipped    = flipped;
    ovl.dat.ovl.orig_erate = erate;
    ovl.dat.ovl.corr_erate = erate;
    ovl.dat.ovl.type       = AS_OVS_TYPE_OVL;
  };

  void     populateOBT(OVSoverlap &obt, readErrorEstimate *readProfile, uint32 iidMin) {
    obt.a_iid              =  a_iid;
    obt.b_iid              = (b_iid_hi << 14) | (b_iid_lo);

    int32 clrLenA = readProfile[obt.a_iid - iidMin].clrLen;
    int32 clrLenB = readProfile[obt.b_iid - iidMin].clrLen;

    //  Swiped from AS_OVS_overlap.C

    uint32  abgn    = (a_hang < 0) ? (0)                : (a_hang);
    uint32  aend    = (b_hang < 0) ? (clrLenA + b_hang) : (clrLenA);
    uint32  bbgn    = (a_hang < 0) ? (-a_hang)          : (0);
    uint32  bend    = (b_hang < 0) ? (clrLenB)          : (clrLenB - b_hang);

    obt.dat.obt.type     =  AS_OVS_TYPE_OBT;

    obt.dat.obt.a_beg    =  abgn;
    obt.dat.obt.a_end    =  aend;
    obt.dat.obt.b_beg    =  bbgn;
    obt.dat.obt.b_end_hi = (bend >> 9) & 0x00003fff;
    obt.dat.obt.b_end_lo = (bend)      & 0x000001ff;

    obt.dat.obt.fwd      = (flipped == false);

    obt.dat.obt.erate    =  erate;
  };
};








void
saveProfile(uint32             iid,
            uint32             iter,
            readErrorEstimate *readProfile) {
  char    N[FILENAME_MAX];
  sprintf(N, "erate-%08u-%02u.dat", iid, iter);

  FILE   *F = fopen(N, "w");

  for (uint32 pp=0; pp<readProfile[iid].clrLen; pp++)
    fprintf(F, "%u %7.4f\n", pp, AS_OVS_decodeQuality(readProfile[iid].errorMeanS[pp]));

  fclose(F);

  FILE *P = popen("gnuplot ", "w"); 
  fprintf(P, "set terminal png\n");
  fprintf(P, "set output   'erate-%08d-%02u.png'\n", iid, iter);
  fprintf(P, "plot [] [0.00:0.25] 'erate-%08d-%02u.dat' using 1:2 with lines\n", iid, iter);
  pclose(P);
}



#if 0
void
computeInitialErrorProfile(gkStore           *gkpStore,
                           uint32             iidMin,
                           uint32             numIIDs,
                           uint64            *overlapIndex,
                           ESToverlap        *overlaps,
                           readErrorEstimate *readProfile) {

  fprintf(stderr, "Processing from IID "F_U32" to "F_U32" out of "F_U32" reads.\n",
          iidMin,
          iidMin + numIIDs,
          gkpStore->gkStore_getNumFragments());

  //  PASS 1; compute an initial error rate estimate for each read
  //
  //  This doesn't ned any temporary storage, the result can be
  //  saved directly in errorMeanS[].

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 iid=0; iid<numIIDs; iid++) {
    if (readProfile[iid].clrLen == 0)
      //  Deleted read.
      continue;

    //  Build a list of the overlap intervals with their error rate

    intervalList<uint32,double>   eRateList;

    for (uint64 oo=overlapIndex[iid]; oo<overlapIndex[iid+1]; oo++) {
      OVSoverlap  obt;

      overlaps[oo].populateOBT(obt, readProfile, iidMin);

      assert(obt.a_iid == iid + iidMin);
      assert(obt.dat.obt.a_beg <= obt.dat.obt.a_end);

      eRateList.add(obt.dat.obt.a_beg,
                    obt.dat.obt.a_end - obt.dat.obt.a_beg,
                    AS_OVS_decodeQuality(obt.dat.obt.erate) / 2);
    }

    //  Convert the list to a sum of error rate per base

    intervalList<uint32,double>   eRateMap(eRateList);

    //  Unpack the list into an array of mean error rate per base

    for (uint32 ii=0; ii<eRateMap.numberOfIntervals(); ii++) {
      double eVal = (eRateMap.depth(ii) > 0) ? (eRateMap.value(ii) / eRateMap.depth(ii)) : 0;

      assert(0.0 <= eVal);
      assert(eVal <= 1.0);

      assert(eRateMap.hi(ii) <= readProfile[iid].clrLen);

      uint16  eEnc = AS_OVS_encodeQuality(eVal);

      for (uint32 pp=eRateMap.lo(ii); pp < eRateMap.hi(ii); pp++)
        readProfile[iid].errorMeanU[pp] = eEnc;
    }

    //  Convert the array of mean error per base into an array of summed error per base

    readProfile[iid].errorMeanS[0] = readProfile[iid].errorMeanU[0];

    for (uint32 ii=1; ii<readProfile[iid].clrLen; ii++)
      readProfile[iid].errorMeanS[ii] = readProfile[iid].errorMeanS[ii-1] + readProfile[iid].errorMeanU[ii];

    //  Keep users entertained.

    if ((iid % 1000) == 0)
      fprintf(stderr, "IID %u\r", iid);
  }

  fprintf(stderr, "IID %u\n", numIIDs);
}
#endif


double
computeEstimatedErate(uint32 iidMin, OVSoverlap &obt, readErrorEstimate *readProfile) {
  uint64  estErrorA = 0;
  uint64  estErrorB = 0;
  int32   olapLen   = 0;

  uint32  ab = obt.dat.obt.a_beg;
  uint32  ae = obt.dat.obt.a_end;

  assert(ab < ae);

#if 0
  for (uint32 xx=ab; xx <= ae; xx++)
    estErrorA += readProfile[obt.a_iid  - iidMin].errorMean[xx];
  estErrorA /= (ae - ab);
#else
  estErrorA = ((readProfile[obt.a_iid  - iidMin].errorMeanS[ae]) -
               (readProfile[obt.a_iid  - iidMin].errorMeanS[ab]));
#endif

  uint32  bb = obt.dat.obt.b_beg;
  uint32  be = obt.dat.obt.b_end_hi * 512 + obt.dat.obt.b_end_lo;

  assert(bb < be);

#if 0
  for (uint32 xx=bb; xx <= be; xx++)
    estErrorB += readProfile[obt.b_iid - iidMin].errorMean[xx];
  estErrorB /= (be - bb);
#else
  estErrorB = ((readProfile[obt.b_iid  - iidMin].errorMeanS[be]) -
               (readProfile[obt.b_iid  - iidMin].errorMeanS[bb]));
#endif

  return(AS_OVS_decodeQuality((estErrorA / 2) + (estErrorB / 2)));
}



void
recomputeErrorProfile(gkStore           *gkpStore,
                      uint32             iidMin,
                      uint32             numIIDs,
                      uint64            *overlapIndex,
                      ESToverlap        *overlaps,
                      readErrorEstimate *readProfile,
                      uint32             iter) {
  uint32      nDiscarded   = 0;
  uint32      nDiscard     = 0;
  uint32      nRemain      = 0;

  fprintf(stderr, "Processing from IID "F_U32" to "F_U32" out of "F_U32" reads, iteration %u.\n",
          iidMin,
          iidMin + numIIDs,
          gkpStore->gkStore_getNumFragments(),
          iter);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 iid=0; iid<numIIDs; iid++) {
    if (readProfile[iid].clrLen == 0)
      //  Deleted read.
      continue;

    //  Build a list of the overlap intervals with their error rate.  Unlike the initial estimates,
    //  we are allowed to skip previously discarded overlaps, and we need to compute estimates and
    //  discard high error overlaps.

    intervalList<uint32,double>   eRateList;

    for (uint64 oo=overlapIndex[iid]; oo<overlapIndex[iid+1]; oo++) {
      OVSoverlap  obt;

      if (overlaps[oo].discarded == true) {
        nDiscarded++;
        continue;
      }

      overlaps[oo].populateOBT(obt, readProfile, iidMin);

      assert(obt.a_iid == iid + iidMin);
      assert(obt.dat.obt.a_beg <= obt.dat.obt.a_end);

      //  Compute the expected erate for this overlap based on our estimated error in both reads,
      //  and filter the overlap if it is higher than this.

      double erate    = AS_OVS_decodeQuality(obt.dat.obt.erate);

      if (iter > 0) {
        double estError = computeEstimatedErate(iidMin, obt, readProfile);

        if (estError + ERATE_TOLERANCE < erate) {
          overlaps[oo].discarded = true;

          nDiscard++;
          continue;
        }
      }

      //  Otherwise, add it to the list of intervals.

      nRemain++;

      eRateList.add(obt.dat.obt.a_beg, obt.dat.obt.a_end - obt.dat.obt.a_beg, erate / 2);
    }

    //  Convert the list to a sum of error rate per base

    intervalList<uint32,double>   eRateMap(eRateList);

    //  Unpack the list into an array of mean error rate per base

    memset(readProfile[iid].errorMeanU, 0, sizeof(uint16) * (readProfile[iid].clrLen + 1));

    for (uint32 ii=0; ii<eRateMap.numberOfIntervals(); ii++) {
      double eVal = (eRateMap.depth(ii) > 0) ? (eRateMap.value(ii) / eRateMap.depth(ii)) : 0;

      assert(0.0 <= eVal);
      assert(eVal <= 1.0);

      assert(eRateMap.hi(ii) <= readProfile[iid].clrLen);

      uint16  eEnc = AS_OVS_encodeQuality(eVal);

      for (uint32 pp=eRateMap.lo(ii); pp < eRateMap.hi(ii); pp++)
        readProfile[iid].errorMeanU[pp] = eEnc;
    }

    //  Keep users entertained.

    if ((iid % 1000) == 0)
      fprintf(stderr, "IID %u\r", iid);
  }

  //  All new estimates are computed.  Convert the array of mean error per base into an array of
  //  summed error per base

  for (uint32 iid=0; iid<numIIDs; iid++) {
    if (readProfile[iid].clrLen == 0)
      continue;

    readProfile[iid].errorMeanS[0] = readProfile[iid].errorMeanU[0];

    for (uint32 ii=1; ii<readProfile[iid].clrLen; ii++)
      readProfile[iid].errorMeanS[ii] = readProfile[iid].errorMeanS[ii-1] + readProfile[iid].errorMeanU[ii];
  }

  //  Report stats.

  fprintf(stderr, "\n");
  fprintf(stderr, "nDiscarded %u (in previous iterations)\n", nDiscarded);
  fprintf(stderr, "nDiscard   %u (in this iteration)\n", nDiscard);
  fprintf(stderr, "nRemain    %u\n", nRemain);
}







void
outputOverlaps(gkStore           *gkpStore,
               uint32             iidMin,
               uint32             numIIDs,
               uint64            *overlapIndex,
               ESToverlap        *overlaps,
               readErrorEstimate *readProfile,
               char              *outputName) {
  uint32      nDiscarded   = 0;
  uint32      nRemain      = 0;

  //  Open an output store

  OverlapStore *storeFile = AS_OVS_createOverlapStore(outputName, TRUE);

  fprintf(stderr, "Processing from IID "F_U32" to "F_U32" out of "F_U32" reads.\n",
          iidMin,
          iidMin + numIIDs,
          gkpStore->gkStore_getNumFragments());

  //  Can't thread.  This does sequential output.  Plus, it doesn't compute anything.

  for (uint32 iid=0; iid<numIIDs; iid++) {
    if (readProfile[iid].clrLen == 0)
      //  Deleted read.
      continue;

    for (uint64 oo=overlapIndex[iid]; oo<overlapIndex[iid+1]; oo++) {
      OVSoverlap  obt;

      if (overlaps[oo].discarded == true) {
        nDiscarded++;
        continue;
      }

      nRemain++;

      //  Emit the overlap

      overlaps[oo].populateOVL(obt);

      AS_OVS_writeOverlapToStore(storeFile, &obt);
    }

    if ((iid % 1000) == 0)
      fprintf(stderr, "IID %u\r", iid);
  }

  AS_OVS_closeOverlapStore(storeFile);

  fprintf(stderr, "\n");
  fprintf(stderr, "nDiscarded %u (in previous iterations)\n", nDiscarded);
  fprintf(stderr, "nRemain    %u\n", nRemain);
}








int
main(int argc, char **argv) {
  char             *gkpName        = 0L;
  gkStore          *gkpStore       = 0L;
  char             *ovlStoreName   = 0L;
  OverlapStore     *ovlStore       = 0L;

  char             *ovlCacheName   = 0L;

  uint32            errorRate      = AS_OVS_encodeQuality(0.015);
  double            errorLimit     = 2.5;

  char             *outputPrefix  = NULL;
  char              logName[FILENAME_MAX] = {0};
  char              sumName[FILENAME_MAX] = {0};
  FILE             *logFile = 0L;
  FILE             *sumFile = 0L;

  bool              doModify = true;  //  Make this false for testing

  argc = AS_configure(argc, argv);

  uint32            minEvidenceOverlap  = MIN(500, AS_OVERLAP_MIN_LEN);
  uint32            minEvidenceCoverage = 1;

  uint32            iidMin   = UINT32_MAX;
  uint32            iidMax   = UINT32_MAX;
  uint32            partNum  = 0;
  uint32            partMax  = 1;

  uint32            minLen   = 0;
  double            minErate = 0;

  int arg=1;
  int err=0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      ovlCacheName = argv[++arg];

    } else if (strcmp(argv[arg], "-b") == 0) {
      iidMin  = atoi(argv[++arg]);
      partNum = 0;
      partMax = 0;

    } else if (strcmp(argv[arg], "-e") == 0) {
      iidMax = atoi(argv[++arg]);
      partNum = 0;
      partMax = 0;

    } else if (strcmp(argv[arg], "-p") == 0) {
      partNum = atoi(argv[++arg]) - 1;
      partMax = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-L") == 0) {
      //minLen = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      //minErate = atof(argv[++arg]);

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (err) {
    exit(1);
  }

  //  Can't use dynamic threads, screws up the thread-local memory allocations.

  omp_set_dynamic(false);

  //  Open gatekeeper store

  fprintf(stderr, "Opening '%s'\n", gkpName);
  gkpStore = new gkStore(gkpName, FALSE, doModify);

  gkpStore->gkStore_enableClearRange(AS_READ_CLEAR_OBTMERGE);
  gkpStore->gkStore_metadataCaching(true);

  //  Compute what to compute.

  if (partNum < partMax) {
    uint32  nf = gkpStore->gkStore_getNumFragments();

    iidMin = (partNum + 0) * nf / partMax + 1;
    iidMax = (partNum + 1) * nf / partMax;

    if (partNum + 1 == partMax)
      iidMax = nf;
  }

  if (iidMin == UINT32_MAX)
    iidMin = 1;

  if (iidMax == UINT32_MAX)
    iidMax = gkpStore->gkStore_getNumFragments();

  uint32    numIIDs = (iidMax - iidMin + 1);

  fprintf(stderr, "  iidMin  = %9u\n", iidMin);
  fprintf(stderr, "  iidMax  = %9u numFrags = %9u\n", iidMax, gkpStore->gkStore_getNumFragments());
  fprintf(stderr, "  partNum = %9u\n", partNum);
  fprintf(stderr, "  partMax = %9u\n", partMax);

  //fprintf(stderr, "OVSoverlap %lu\n", sizeof(OVSoverlap));
  //fprintf(stderr, "ESToverlap %lu\n", sizeof(ESToverlap));

  //  Load read metadata, clear ranges, read lengths, and deleted status.

  fprintf(stderr, "Initializing profiles\n");

  uint64              readProfileSize = 0;
  readErrorEstimate  *readProfile     = new readErrorEstimate [numIIDs];

  for (uint32 iid=0; iid<numIIDs; iid++) {
    gkFragment   fr;

    gkpStore->gkStore_getFragment(iid + iidMin, &fr, GKFRAGMENT_INF);

    readProfileSize += readProfile[iid].initialize(fr);

    if ((iid % 10000) == 0)
      fprintf(stderr, "  %u reads\r", iid);
  }

  fprintf(stderr, "  %u reads\n", numIIDs);
  fprintf(stderr, "  %lu GB\n", readProfileSize >> 30);

  //  Open overlap stores

  fprintf(stderr, "Opening '%s'\n", ovlStoreName);
  ovlStore = AS_OVS_openOverlapStore(ovlStoreName);

  fprintf(stderr, "Finding number of overlaps\n");

  AS_OVS_setRangeOverlapStore(ovlStore,   iidMin, iidMax);

  uint64    numOvls = AS_OVS_numOverlapsInRange(ovlStore);
  
  uint64      *overlapIndex = new uint64     [numIIDs + 1];
  uint32      *overlapLen   = AS_OVS_numOverlapsPerFrag(ovlStore);

  overlapIndex[0] = 0;

  for (uint32 iid=0; iid<numIIDs; iid++)
    overlapIndex[iid+1] = overlapIndex[iid] + overlapLen[iid];

  assert(overlapIndex[numIIDs] == numOvls);

  safe_free(overlapLen);

  //  Load overlaps.

  fprintf(stderr, "Loading overlaps\n");
  fprintf(stderr, "  number   %lu overlaps\n",           numOvls);
  fprintf(stderr, "  index    %lu GB\n",                 (sizeof(uint64)     * numIIDs) >> 30);
  fprintf(stderr, "  overlaps %lu GB (previous size)\n", (sizeof(OVSoverlap) * numOvls) >> 30);
  fprintf(stderr, "  overlaps %lu GB\n",                 (sizeof(ESToverlap) * numOvls) >> 30);

  ESToverlap       *overlaps     = NULL;
  memoryMappedFile *overlapsMMF  = NULL;

  if (AS_UTL_fileExists(ovlCacheName, FALSE, FALSE)) {
    fprintf(stderr, "  cache '%s' detected, load averted\n", ovlCacheName);

    overlapsMMF = new memoryMappedFile(ovlCacheName);
    overlaps    = (ESToverlap *)overlapsMMF->get(0);

  } else {
    FILE             *ESTcache     = NULL;
    uint32            overlapblock = 100000000;
    OVSoverlap       *overlapsload = new OVSoverlap [overlapblock];

    overlaps       = new ESToverlap [numOvls];

    if (ovlCacheName) {
      errno = 0;
      ESTcache = fopen(ovlCacheName, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for writing: %s\n", ovlCacheName, strerror(errno)), exit(1);
    }

    for (uint64 no=0; no<numOvls; ) {
      uint64 nLoad  = AS_OVS_readOverlapsFromStore(ovlStore,
                                                   overlapsload,
                                                   overlapblock,
                                                   AS_OVS_TYPE_OVL,
                                                   false);

      for (uint32 xx=0; xx<nLoad; xx++)
        overlaps[no++].populate(overlapsload[xx]);

      if (ESTcache)
        fwrite(overlaps + no - nLoad, sizeof(ESToverlap), nLoad, ESTcache);

      fprintf(stderr, "  loading overlaps: %lu out of %lu (%.4f%%)\r",
              no, numOvls, 100.0 * no / numOvls);
    }

    delete [] overlapsload;

    if (ESTcache)
      fclose(ESTcache);

    fprintf(stderr, "\n");
    fprintf(stderr, "  loaded and cached %lu overlaps.\n", numOvls);
  }

  AS_OVS_closeOverlapStore(ovlStore);

  ovlStore   = NULL;

  //  Allocate space for the result.

  double     *erate5 = new double [numIIDs];
  double     *erate3 = new double [numIIDs];

  //  Compute the initial read profile.

#if 0
  computeInitialErrorProfile(gkpStore, iidMin, numIIDs,
                             overlapIndex,
                             overlaps,
                             readProfile);
#endif

  //  Recompute, using the existing profile to weed out probably false overlaps.

  for (uint32 ii=0; ii<4; ii++)
    recomputeErrorProfile(gkpStore, iidMin, numIIDs,
                          overlapIndex,
                          overlaps,
                          readProfile,
                          ii);


  outputOverlaps(gkpStore, iidMin, numIIDs,
                 overlapIndex,
                 overlaps,
                 readProfile,
                 "TEST.ovlStore");

  if (overlapsMMF) {
    delete    overlapsMMF;
  } else {
    delete [] overlaps;
  }

  exit(0);
}

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
#include "AS_MSG_pmesg.H"
#include "AS_PER_gkpStore.H"

#include "AS_OVS_overlap.H"
#include "AS_OVS_overlapFile.H"
#include "AS_OVS_overlapStore.H"


#include "intervalList.H"

#include <map>
#include <set>
#include <list>
#include <vector>
#include <algorithm>

#ifndef BROKEN_CLANG_OpenMP
#include <omp.h>
#endif

#include "../AS_OBT/AS_OBT_overlaps.C"

using namespace std;


#define ERATE_TOLERANCE 0.03


class readErrorEstimate {
public:
  readErrorEstimate() {
    clrBgn = 0;
    clrEnd = 0;
    clrLen = 0;

    errorMean       = NULL;
    //errorStdDev     = NULL;
    //errorConfidence = NULL;
    errorMeanU      = NULL;
  };

  ~readErrorEstimate() {
    delete [] errorMean;
    //delete [] errorStdDev;
    //delete [] errorConfidence;
    delete [] errorMeanU;
  };

  void       initialize(gkFragment &fr) {
    //uint32      lid = fr.gkFragment_getLibraryIID();
    //gkLibrary  *lb  = gkpStore->gkStore_getLibrary(lid);

    clrBgn = 0;
    clrEnd = 0;
    clrLen = 0;

    if (fr.gkFragment_getIsDeleted() == 1)
      return;

    clrBgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_LATEST);
    clrEnd = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_LATEST);
    clrLen = clrEnd - clrBgn;

    errorMean       = new double [clrLen + 1];
    //errorStdDev     = new double [clrLen + 1];
    //errorConfidence = new double [clrLen + 1];
    errorMeanU      = new double [clrLen + 1];
  };

  uint32     clrBgn;
  uint32     clrEnd;
  uint32     clrLen;

  double    *errorMean;
  //double    *errorStdDev;
  //double    *errorConfidence;
  double    *errorMeanU;
};




uint32
loadOverlaps(uint32         iid,
             OVSoverlap   *&ovl,
             uint32        &ovlLen,
             uint32        &ovlMax,
             OverlapStore  *ovlPrimary,
             OverlapStore  *ovlSecondary);




void
saveProfile(uint32             iid,
            uint32             iter,
            readErrorEstimate *readProfile,
            uint32            *erateCnt) {
  char    N[FILENAME_MAX];
  sprintf(N, "erate-%08u-%02u.dat", iid, iter);

  FILE   *F = fopen(N, "w");

  if (erateCnt)
    for (uint32 pp=0; pp<readProfile[iid].clrLen; pp++)
      fprintf(F, "%u %.4f %u\n", pp, readProfile[iid].errorMean[pp], erateCnt[pp]);
  else
    for (uint32 pp=0; pp<readProfile[iid].clrLen; pp++)
      fprintf(F, "%u %.4f\n", pp, readProfile[iid].errorMean[pp]);

  fclose(F);

  FILE *P = popen("gnuplot ", "w"); 
  fprintf(P, "set terminal png\n");
  fprintf(P, "set output   'erate-%08d-%02u.png'\n", iid, iter);
  fprintf(P, "plot [] [0.00:0.25] 'erate-%08d-%02u.dat' using 1:2 with lines\n", iid, iter);
  pclose(P);
}




void
computeInitialErrorProfile(gkStore           *gkpStore,
                           uint32             iidMin,
                           uint32             iidMax,
                           OverlapStore      *ovlPrimary,
                           OverlapStore      *ovlSecondary,
                           readErrorEstimate *readProfile) {

  //  Allocate space for the compute.

  uint32      ovlLen       = 0;
  uint32      ovlMax       = 64 * 1024;
  OVSoverlap *ovl          = new OVSoverlap [ovlMax];

  memset(ovl, 0, sizeof(OVSoverlap) * ovlMax);

#ifdef SLOW
  double     *erateSum     = new double [AS_READ_MAX_NORMAL_LEN];
  uint32     *erateCnt     = new uint32 [AS_READ_MAX_NORMAL_LEN];
#else
  uint32     *erateCnt     = NULL;
#endif

  AS_OVS_resetRangeOverlapStore(ovlPrimary);
  AS_OVS_resetRangeOverlapStore(ovlSecondary);

  fprintf(stderr, "Processing from IID "F_U32" to "F_U32" out of "F_U32" reads.\n",
          iidMin,
          iidMax,
          gkpStore->gkStore_getNumFragments());

  //  PASS 1; compute an initial error rate estimate for each read

  for (uint32 iid=iidMin; iid<=iidMax; iid++) {
    if (readProfile[iid].clrLen == 0)
      //  Deleted read.
      continue;

    if (loadOverlaps(iid, ovl, ovlLen, ovlMax, ovlPrimary, ovlSecondary) == 0)
      continue;
    
#ifdef SLOW
    for (uint32 pp=0; pp<=readProfile[iid].clrLen; pp++) {
      erateSum[pp] = 0;
      erateCnt[pp] = 0;
    }
#endif

    intervalList<uint32,double>   eRateList;

    for (uint32 oo=0; oo<ovlLen; oo++) {
      assert(ovl[oo].a_iid == iid);
      AS_OVS_convertOVLoverlapToOBToverlap(ovl[oo],
                                           readProfile[ovl[oo].a_iid].clrLen,
                                           readProfile[ovl[oo].b_iid].clrLen);
      assert(ovl[oo].a_iid == iid);

      double               erate = AS_OVS_decodeQuality(ovl[oo].dat.obt.erate) / 2;

#ifdef SLOW
      //  Simple slow version.  Can we use intervalList for this?
      for (uint32 xx=ovl[oo].dat.obt.a_beg; xx <= ovl[oo].dat.obt.a_end; xx++) {
        erateSum[xx] += erate;
        erateCnt[xx] += 1;
      }
#else
      eRateList.add(ovl[oo].dat.obt.a_beg, ovl[oo].dat.obt.a_end - ovl[oo].dat.obt.a_beg, erate);
#endif
    }


#ifdef SLOW
    for (uint32 pp=0; pp <= readProfile[iid].clrLen; pp++) {
      readProfile[iid].errorMean[pp]        = 0;
      //readProfile[iid].errorStdDev[pp]      = 0;
      //readProfile[iid].errorConfidence[pp]  = 0;

      if (erateCnt[pp] > 0)
        readProfile[iid].errorMean[pp] = erateSum[pp] / erateCnt[pp];
    }
#else
    intervalList<uint32,double>   eRateMap(eRateList);


    for (uint32 ii=0; ii<eRateMap.numberOfIntervals(); ii++) {
      double eVal = (eRateMap.depth(ii) > 0) ? (eRateMap.value(ii) / eRateMap.depth(ii)) : 0;

      for (uint32 pp=eRateMap.lo(ii); pp < eRateMap.hi(ii); pp++)
        readProfile[iid].errorMean[pp] = eVal;
    }
#endif

#if 0
    if ((iid == 9) || (iid == 295))
      saveProfile(iid, 0, readProfile, erateCnt);
#endif

    //if ((iid % 1000) == 0)
    //  fprintf(stderr, "IID %u OVLs %u\r", iid, ovlLen);

    //if (ovlLen > 100)
    //  fprintf(stderr, "IID %u OVLs %u\n", iid, ovlLen);
  }

  fprintf(stderr, "\n");

  delete [] ovl;
#ifdef SLOW
  delete [] erateSum;
  delete [] erateCnt;
#endif
}










double
computeEstimatedErate(OVSoverlap *ovl, readErrorEstimate *readProfile) {
  double  estError = 0;
  int32   olapLen  = 0;

  uint32  ab = ovl->dat.obt.a_beg;
  uint32  ae = ovl->dat.obt.a_end;

  assert(ab < ae);

  for (uint32 xx=ab; xx <= ae; xx++)
    estError += readProfile[ovl->a_iid].errorMean[xx];

  uint32  bb = ovl->dat.obt.b_beg;
  uint32  be = ovl->dat.obt.b_end_hi * 512 + ovl->dat.obt.b_end_lo;

  if (be < bb) {
    uint32 tt = bb;
    bb = be;
    be = tt;
  }

  assert(bb < be);

  for (uint32 xx=bb; xx <= be; xx++)
    estError += readProfile[ovl->b_iid].errorMean[xx];

  estError = estError / (((ae - ab) + (be - bb)) / 2);

  return(estError);
}




void
recomputeErrorProfile(gkStore           *gkpStore,
                      uint32             iidMin,
                      uint32             iidMax,
                      OverlapStore      *ovlPrimary,
                      OverlapStore      *ovlSecondary,
                      readErrorEstimate *readProfile,
                      uint32             iter) {

  //  Allocate space for the compute.

  uint32      ovlLen       = 0;
  uint32      ovlMax       = 64 * 1024;
  OVSoverlap *ovl          = new OVSoverlap [ovlMax];

  memset(ovl, 0, sizeof(OVSoverlap) * ovlMax);

#ifdef SLOW
  double     *erateSum     = new double [AS_READ_MAX_NORMAL_LEN];
  uint32     *erateCnt     = new uint32 [AS_READ_MAX_NORMAL_LEN];
#else
  uint32     *erateCnt     = NULL;
#endif

  uint32      nDiscards    = 0;
  uint32      nRemain      = 0;

  AS_OVS_resetRangeOverlapStore(ovlPrimary);
  AS_OVS_resetRangeOverlapStore(ovlSecondary);

  fprintf(stderr, "Processing from IID "F_U32" to "F_U32" out of "F_U32" reads, iteration %u.\n",
          iidMin,
          iidMax,
          gkpStore->gkStore_getNumFragments(),
          iter);

  for (uint32 iid=iidMin; iid<=iidMax; iid++) {
    if (readProfile[iid].clrLen == 0)
      //  Deleted read.
      continue;

    if (loadOverlaps(iid, ovl, ovlLen, ovlMax, ovlPrimary, ovlSecondary) == 0)
      continue;

#ifdef SLOW
    for (uint32 pp=0; pp<=readProfile[iid].clrLen; pp++) {
      erateSum[pp] = 0;
      erateCnt[pp] = 0;
    }
#endif

    intervalList<uint32,double>   eRateList;

    for (uint32 oo=0; oo<ovlLen; oo++) {
      assert(ovl[oo].a_iid == iid);
      AS_OVS_convertOVLoverlapToOBToverlap(ovl[oo],
                                           readProfile[ovl[oo].a_iid].clrLen,
                                           readProfile[ovl[oo].b_iid].clrLen);
      assert(ovl[oo].a_iid == iid);

      //  Compute the expected erate for this overlap based on our estimated error in both reads

      double estError = computeEstimatedErate(ovl + oo, readProfile);
      double erate    = AS_OVS_decodeQuality(ovl[oo].dat.obt.erate);

#ifndef NODEBUG
      if (((ovl[oo].a_iid == 1) && (ovl[oo].b_iid == 6836)) ||
          ((ovl[oo].b_iid == 1) && (ovl[oo].a_iid == 6836))) {
        char ovlstr[1024];

        fprintf(stderr, "%s\n", AS_OVS_toString(ovlstr, ovl[oo]));
        //fprintf(stderr, "%s\n", AS_OVS_toString(ovlstr, obt[oo]));

        fprintf(stderr, "  estError = %f  erate = %f\n", estError, erate);
      }
#endif

#if 0
      if ((iid == 9) || (iid == 295))
        fprintf(stdout, "%8u %8u  %5u %5u  %5u %5u  len %5u %5u est %7.5f tru %7.5f %s\n",
                ovl[oo].a_iid, ovl[oo].b_iid,
                ovl[oo].dat.obt.a_beg,
                ovl[oo].dat.obt.a_end,
                ovl[oo].dat.obt.b_beg,
                ovl[oo].dat.obt.b_end_hi * 512 + ovl[oo].dat.obt.b_end_lo,
                readProfile[ovl[oo].a_iid].clrLen,
                readProfile[ovl[oo].b_iid].clrLen,
                estError,
                erate,
                (estError + ERATE_TOLERANCE < erate) ? "DISCARD" : "SAVE");
#endif

      if (estError + ERATE_TOLERANCE < erate) {
        nDiscards++;
        continue;
      }

      nRemain++;

      erate /= 2;

#ifdef SLOW
      //  Simple slow version.  Can we use intervalList for this?
      for (uint32 xx=ovl[oo].dat.obt.a_beg; xx <= ovl[oo].dat.obt.a_end; xx++) {
        erateSum[xx] += erate;
        erateCnt[xx] += 1;
      }
#else
      eRateList.add(ovl[oo].dat.obt.a_beg, ovl[oo].dat.obt.a_end - ovl[oo].dat.obt.a_beg, erate);
#endif

    }

#ifdef SLOW
    for (uint32 pp=0; pp <= readProfile[iid].clrLen; pp++) {
      readProfile[iid].errorMeanU[pp]       = 0;
      //readProfile[iid].errorStdDev[pp]      = 0;
      //readProfile[iid].errorConfidence[pp]  = 0;

      if (erateCnt[pp] > 0)
        readProfile[iid].errorMeanU[pp] = erateSum[pp] / erateCnt[pp];
    }
#else
    intervalList<uint32,double>   eRateMap(eRateList);

    //  Can have zero if there isn't overlap coverage!
    for (uint32 ii=0; ii<eRateMap.numberOfIntervals(); ii++) {
      double eVal = (eRateMap.depth(ii) > 0) ? (eRateMap.value(ii) / eRateMap.depth(ii)) : 0;

      for (uint32 pp=eRateMap.lo(ii); pp < eRateMap.hi(ii); pp++)
        readProfile[iid].errorMeanU[pp] = eVal;
    }
#endif

#if 0
    if ((iid == 9) || (iid == 295))
      saveProfile(iid, iter, readProfile, erateCnt);
#endif

    //if ((iid % 1000) == 0)
    //  fprintf(stderr, "IID %u OVLs %u\r", iid, ovlLen);
  }

  delete [] ovl;
#ifdef SLOW
  delete [] erateSum;
  delete [] erateCnt;
#endif

  //  Update all the means

  for (uint32 iid=iidMin; iid<=iidMax; iid++) {
    double *t = readProfile[iid].errorMean;
  
    readProfile[iid].errorMean  = readProfile[iid].errorMeanU;
    readProfile[iid].errorMeanU = t;
  }

  //  Report stats.

  fprintf(stderr, "\n");
  fprintf(stderr, "nDiscards %u\n", nDiscards);
  fprintf(stderr, "nRemain   %u\n", nRemain);
}







void
outputOverlaps(gkStore           *gkpStore,
               uint32             iidMin,
               uint32             iidMax,
               OverlapStore      *ovlPrimary,
               OverlapStore      *ovlSecondary,
               readErrorEstimate *readProfile,
               char              *outputName) {

  //  Allocate space for the compute.

  uint32      ovlLen       = 0;
  uint32      ovlMax       = 64 * 1024;
  OVSoverlap *ovl          = new OVSoverlap [ovlMax];
  OVSoverlap *obt          = new OVSoverlap [ovlMax];

  memset(ovl, 0, sizeof(OVSoverlap) * ovlMax);

  uint32      nDiscards    = 0;
  uint32      nRemain      = 0;

  AS_OVS_resetRangeOverlapStore(ovlPrimary);
  AS_OVS_resetRangeOverlapStore(ovlSecondary);

  //  Open an output store

  OverlapStore *storeFile = AS_OVS_createOverlapStore(outputName, TRUE);

  fprintf(stderr, "Processing from IID "F_U32" to "F_U32" out of "F_U32" reads.\n",
          iidMin,
          iidMax,
          gkpStore->gkStore_getNumFragments());

  for (uint32 iid=iidMin; iid<=iidMax; iid++) {
    if (readProfile[iid].clrLen == 0)
      //  Deleted read.
      continue;

    if (loadOverlaps(iid, ovl, ovlLen, ovlMax, ovlPrimary, ovlSecondary) == 0)
      continue;

    memcpy(obt, ovl, sizeof(OVSoverlap) * ovlLen);

    for (uint32 oo=0; oo<ovlLen; oo++) {
      assert(ovl[oo].a_iid == iid);
      AS_OVS_convertOVLoverlapToOBToverlap(obt[oo],
                                           readProfile[obt[oo].a_iid].clrLen,
                                           readProfile[obt[oo].b_iid].clrLen);
      assert(ovl[oo].a_iid == iid);

      double estError = computeEstimatedErate(obt + oo, readProfile);
      double erate    = AS_OVS_decodeQuality(obt[oo].dat.obt.erate);

#ifndef NODEBUG
      if (((ovl[oo].a_iid == 1) && (ovl[oo].b_iid == 6836)) ||
          ((ovl[oo].b_iid == 1) && (ovl[oo].a_iid == 6836))) {
        char ovlstr[1024];

        fprintf(stderr, "%s\n", AS_OVS_toString(ovlstr, ovl[oo]));
        fprintf(stderr, "%s\n", AS_OVS_toString(ovlstr, obt[oo]));

        fprintf(stderr, "  estError = %f  erate = %f\n", estError, erate);
      }
#endif

      if (estError + ERATE_TOLERANCE < erate) {
        nDiscards++;
        continue;
      }

      nRemain++;

      //  Emit the overlap

      AS_OVS_writeOverlapToStore(storeFile, ovl + oo);
    }

    //if ((iid % 1000) == 0)
    //  fprintf(stderr, "IID %u OVLs %u\r", iid, ovlLen);
  }

  AS_OVS_closeOverlapStore(storeFile);

  fprintf(stderr, "\n");
  fprintf(stderr, "nDiscards %u\n", nDiscards);
  fprintf(stderr, "nRemain   %u\n", nRemain);

  delete [] ovl;
  delete [] obt;
}








int
main(int argc, char **argv) {
  char             *gkpName          = 0L;
  gkStore          *gkpStore         = 0L;
  char             *ovlPrimaryName   = 0L;
  OverlapStore     *ovlPrimary       = 0L;
  char             *ovlSecondaryName = 0L;
  OverlapStore     *ovlSecondary     = 0L;

  uint32            errorRate    = AS_OVS_encodeQuality(0.015);
  double            errorLimit   = 2.5;

  char             *outputPrefix  = NULL;
  char              logName[FILENAME_MAX] = {0};
  char              sumName[FILENAME_MAX] = {0};
  FILE             *logFile = 0L;
  FILE             *sumFile = 0L;

  bool              doModify = true;  //  Make this false for testing

  uint32            iidMin = 1;
  uint32            iidMax = UINT32_MAX;

  argc = AS_configure(argc, argv);

  uint32            minEvidenceOverlap  = MIN(500, AS_OVERLAP_MIN_LEN);
  uint32            minEvidenceCoverage = 1;

  int arg=1;
  int err=0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      if (ovlPrimaryName == NULL)
        ovlPrimaryName = argv[++arg];
      else if (ovlSecondary == NULL)
        ovlSecondaryName = argv[++arg];
      else {
        fprintf(stderr, "ERROR:  Only two obtStores (-O) allowed.\n");
        err++;
      }

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (err) {
    exit(1);
  }

  //  Open gatekeeper store

  gkpStore = new gkStore(gkpName, FALSE, doModify);

  gkpStore->gkStore_enableClearRange(AS_READ_CLEAR_OBTMERGE);
  gkpStore->gkStore_metadataCaching(true);

  if (iidMin < 1)
    iidMin = 1;
  if (iidMax > gkpStore->gkStore_getNumFragments())
    iidMax = gkpStore->gkStore_getNumFragments();

  //  Open overlap stores

  if (ovlPrimaryName != NULL)
    ovlPrimary = AS_OVS_openOverlapStore(ovlPrimaryName);

  if (ovlSecondary != NULL)
    ovlSecondary = AS_OVS_openOverlapStore(ovlSecondaryName);
 
  //  Load read metadata, clear ranges, read lengths, and deleted status.

  readErrorEstimate  *readProfile = new readErrorEstimate [iidMax + 1];

  for (uint32 iid=iidMin; iid<=iidMax; iid++) {
    gkFragment   fr;

    gkpStore->gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);

    readProfile[iid].initialize(fr);
  }

  //  Allocate space for the result.

  double     *erate5 = new double [iidMax];
  double     *erate3 = new double [iidMax];


  //  Compute the initial read profile.

  computeInitialErrorProfile(gkpStore, iidMin, iidMax,
                             ovlPrimary,
                             ovlSecondary,
                             readProfile);

  //  Recompute, using the existing profile to weed out probably false overlaps.

  for (uint32 ii=1; ii<4; ii++)
    recomputeErrorProfile(gkpStore, iidMin, iidMax,
                          ovlPrimary,
                          ovlSecondary,
                          readProfile,
                          ii);


  outputOverlaps(gkpStore, iidMin, iidMax,
                 ovlPrimary,
                 ovlSecondary,
                 readProfile,
                 "TEST.ovlStore");

  exit(0);
}

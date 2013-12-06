
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2011, J. Craig Venter Institute.
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

const char *mainid = "$Id$";

#include "AS_global.H"

#include "finalTrim.H"

#include "AS_UTL_decodeRange.H"
#include "AS_OBT_overlaps.H"


bool
trimWithoutOverlaps(OVSoverlap  *ovl,
                    uint32       ovlLen,
                    gkFragment  &fr,
                    uint32       ibgn,
                    uint32       iend,
                    uint32      &fbgn,
                    uint32      &fend,
                    char        *logMsg,
                    uint32       errorRate,
                    double       errorLimit,
                    bool         qvTrimAllowed) {

  fbgn = ibgn;
  fend = iend;

  if (qvTrimAllowed == false) {
    strcpy(logMsg, "\tno overlaps, qv trim not allowed");
    return(true);
  }

  double  minQuality = qual.lookupNumber(20);

  doTrim(&fr, minQuality, fbgn, fend);

  //  Intersect the two trims, or return failure.

  if ((iend < fbgn) ||
      (fend < ibgn)) {
    strcpy(logMsg, "\tno overlaps, qv trim allowed, inconsistent trimming");
    return(false);
  } else {
    fbgn = MAX(fbgn, ibgn);
    fend = MIN(fend, iend);

    strcpy(logMsg, "\tno overlaps, qv trim allowed");
    return(true);
  }
}


//  Enforce any maximum clear range, if it exists (mbgn < mend)
//
//  There are six cases:
//
//       ---MAX-RANGE---
//   ---
//     -------------------
//     -----
//            -----
//                   -----
//                       ---
//
//  If the begin is below the max-bgn, we reset it to max-bgn.
//  If the end   is after the max-end, we reset it to max-end.
//
//  If after the resets we have an invalid clear (bgn > end)
//  the original clear range was completely outside the max range.
//
bool
enforceMaximumClearRange(OVSoverlap  *ovl,
                         uint32       ovlLen,
                         gkFragment  &fr,
                         uint32       ibgn,
                         uint32       iend,
                         uint32      &fbgn,
                         uint32      &fend,
                         char        *logMsg,
                         uint32       errorRate,
                         double       errorLimit,
                         bool         qvTrimAllowed) {

  uint32 mbgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_MAX);
  uint32 mend = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_MAX);

  if (mbgn >= mend)
    //  No maximum clear range set.
    return(true);

  if (fbgn == fend)
    return(true);

  assert(mbgn <  mend);
  assert(fbgn <= fend);

  if ((fend < mbgn) ||
      (mend < fbgn)) {
    //  Final clear not intersecting maximum clear.
    strcat(logMsg, (logMsg[0]) ? " - " : "\t");
    strcat(logMsg, "outside maximum allowed clear range");
    return(false);

  } else if ((fbgn < mbgn) ||
             (mend < fend)) {
    //  Final clear extends outside the maximum clear.
    fbgn = MAX(fbgn, mbgn);
    fend = MIN(fend, mend);

    strcat(logMsg, (logMsg[0]) ? " - " : "\t");
    strcat(logMsg, "adjusted to obey maximum allowed clear range");
    return(true);
  } else {
    //  Final clear already within the maximum clear.
    return(true);
  }
}



int
main(int argc, char **argv) {
  gkStore          *gkpStore     = 0L;
  OverlapStore     *ovlPrimary   = 0L;
  OverlapStore     *ovlSecondary = 0L;

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
      gkpStore = new gkStore(argv[++arg], FALSE, doModify);
      gkpStore->gkStore_enableClearRange(AS_READ_CLEAR_OBTMERGE);
      gkpStore->gkStore_metadataCaching(true);

    } else if (strcmp(argv[arg], "-O") == 0) {
      if (ovlPrimary == NULL)
        ovlPrimary = AS_OVS_openOverlapStore(argv[++arg]);
      else if (ovlSecondary == NULL)
        ovlSecondary = AS_OVS_openOverlapStore(argv[++arg]);
      else {
        fprintf(stderr, "ERROR:  Only two obtStores (-O) allowed.\n");
        err++;
      }

    } else if (strcmp(argv[arg], "-e") == 0) {
      double erate = atof(argv[++arg]);
      if (erate > AS_MAX_ERROR_RATE)
        fprintf(stderr, "ERROR:  Error rate (-e) %s too large; must be 'fraction error' and below %f\n", argv[arg], AS_MAX_ERROR_RATE), exit(1);
      errorRate = AS_OVS_encodeQuality(erate);

    } else if (strcmp(argv[arg], "-E") == 0) {
      errorLimit = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-ol") == 0) {
      minEvidenceOverlap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-oc") == 0) {
      minEvidenceCoverage = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-n") == 0) {
      doModify = false;

    } else if (strcmp(argv[arg], "-t") == 0) {
      AS_UTL_decodeRange(argv[++arg], iidMin, iidMax);

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if ((gkpStore     == NULL) ||
      (ovlPrimary   == NULL) ||
      (outputPrefix == NULL) ||
      (err)) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore [-O ovlStore] -o outputPrefix\n", argv[0]);
    fprintf(stderr, "   -e erate       allow 'erate' percent error\n");
    fprintf(stderr, "   -E elimit      allow 'elimit' errors (only used in 'largestCovered')\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   -ol            for 'largest covered', the minimum evidence overlap length\n");
    fprintf(stderr, "   -oc            for 'largest covered', the minimum evidence overlap coverage\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   -n             do not modify reads in gkpStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   -t bgn-end     limit processing to only reads from bgn to end (inclusive)\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  sprintf(logName, "%s.log",     outputPrefix);
  sprintf(sumName, "%s.summary", outputPrefix);

  logFile = fopen(logName, "w");
  sumFile = fopen(sumName, "w");

  fprintf(logFile, "iid\tinitL\tinitR\tfinalL\tfinalR\tmessage (DEL=deleted NOC=no change MOD=modified)\n");

  fprintf(sumFile, "Overlap error rate     <= %.2f%% error\n", 100.0 * AS_OVS_decodeQuality(errorRate));
  fprintf(sumFile, "Overlap error limit    <= %.2f errors\n", errorLimit);
  fprintf(sumFile, "Overlap min overlap    >= %u base%s (for 'largest covered')\n", minEvidenceOverlap,  (minEvidenceOverlap  == 1) ? "" : "s");
  fprintf(sumFile, "Overlap min coverage   >= %u read%s (for 'largest covered')\n", minEvidenceCoverage, (minEvidenceCoverage == 1) ? "" : "s");


  uint32      ovlLen       = 0;
  uint32      ovlMax       = 64 * 1024;
  OVSoverlap *ovl          = new OVSoverlap [ovlMax];

  memset(ovl, 0, sizeof(OVSoverlap) * ovlMax);

  gkFragment  fr;
  char        logMsg[1024] = {0};

  if (iidMin < 1)
    iidMin = 1;
  if (iidMax > gkpStore->gkStore_getNumFragments())
    iidMax = gkpStore->gkStore_getNumFragments();

  fprintf(stderr, "Processing from IID "F_U32" to "F_U32" out of "F_U32" reads.\n",
          iidMin,
          iidMax,
          gkpStore->gkStore_getNumFragments());

  for (uint32 iid=iidMin; iid<=iidMax; iid++) {

    //  The QLT is needed ONLY for Sanger reads, and is a total waste of bandwidth for all other
    //  read types.  Perhaps this can be moved into initialTrim -- compute the strict QV trim there
    //  too, save it in the OBTMERGE clear range
    //
    gkpStore->gkStore_getFragment(iid, &fr, GKFRAGMENT_QLT);

    logMsg[0] = 0;

    uint32      lid = fr.gkFragment_getLibraryIID();
    gkLibrary  *lb  = gkpStore->gkStore_getLibrary(lid);

    //  If the fragment is deleted, do nothing.  If the fragment was deleted AFTER overlaps were
    //  generated, then the overlaps will be out of sync -- we'll get overlaps for these fragments
    //  we skip.
    //
    if (fr.gkFragment_getIsDeleted() == 1)
      continue;

    //  If it did not request trimming, do nothing.  Similar to the above, we'll get overlaps to
    //  fragments we skip.
    //
    if ((lb->doTrim_finalLargestCovered == false) &&
        (lb->doTrim_finalEvidenceBased  == false) &&
        (lb->doTrim_finalBestEdge       == false))
      continue;

    uint32      ibgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_OBTINITIAL);
    uint32      iend = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_OBTINITIAL);

    //  Is the clear range valid?  If we skip initial trim for this read, there is no clear range defined.
    if ((ibgn == 1) && (iend == 0)) {
      ibgn = 0;
      iend = fr.gkFragment_getSequenceLength();
    }

    uint32      fbgn = ibgn;
    uint32      fend = iend;

    loadOverlaps(iid, ovl, ovlLen, ovlMax, ovlPrimary, ovlSecondary);

    bool isGood = false;

    //  If there are no overlaps for this read, do nothing.
    if ((iid < ovl[0].a_iid) ||
        (ovlLen == 0)) {
      isGood = trimWithoutOverlaps(ovl, ovlLen,
                                   fr,
                                   ibgn, iend, fbgn, fend,
                                   logMsg,
                                   errorRate,
                                   errorLimit,
                                   lb->doTrim_initialQualityBased);
      assert(fbgn <= fend);

    }

    //  Use the largest region covered by overlaps as the trim
    else if        (lb->doTrim_finalLargestCovered == true) {
      isGood = largestCovered(ovl, ovlLen,
                              fr,
                              ibgn, iend, fbgn, fend,
                              logMsg,
                              errorRate,
                              errorLimit,
                              lb->doTrim_initialQualityBased,
                              minEvidenceOverlap,
                              minEvidenceCoverage);
      assert(fbgn <= fend);

    }

    //  Use the largest region covered by overlaps as the trim
    else if        (lb->doTrim_finalBestEdge == true) {
      isGood = bestEdge(ovl, ovlLen,
                        fr,
                        ibgn, iend, fbgn, fend,
                        logMsg,
                        errorRate,
                        errorLimit,
                        lb->doTrim_initialQualityBased,
                        minEvidenceOverlap,
                        minEvidenceCoverage);
      assert(fbgn <= fend);

    }

    //  Do Sanger-style heuristics
    else if (lb->doTrim_finalEvidenceBased == true) {
      isGood = evidenceBased(ovl, ovlLen,
                             fr,
                             ibgn, iend, fbgn, fend,
                             logMsg,
                             errorRate,
                             errorLimit,
                             lb->doTrim_initialQualityBased);
      assert(fbgn <= fend);

    }

    //  Do nothing.  Really shouldn't get here.
    else {
      assert(0);
      continue;
    }

    //  Enforce the maximum clear range

    if (isGood)
      isGood = enforceMaximumClearRange(ovl, ovlLen,
                                        fr,
                                        ibgn, iend, fbgn, fend,
                                        logMsg,
                                        errorRate,
                                        errorLimit,
                                        lb->doTrim_initialQualityBased);

    assert(fbgn <= fend);

    //  Bad trimming?  Invalid clear range?  Too small?

    if ((isGood == false) ||
        (fbgn > fend) ||
        (fend - fbgn < AS_READ_MIN_LEN)) {

      assert(fbgn <= fend);

      if (doModify) {
        fr.gkFragment_setClearRegion(fbgn, fend, AS_READ_CLEAR_OBTMERGE);
        gkpStore->gkStore_setFragment(&fr);
        gkpStore->gkStore_delFragment(iid);
      }

      fprintf(logFile, F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\tDEL%s\n",
              iid,
              ibgn, iend,
              fbgn, fend,
              (logMsg[0] == 0) ? "" : logMsg);
      continue;
    }

    //  Good trimming, just didn't do anything.

    assert(fbgn <= fend);

    if ((ibgn == fbgn) &&
        (iend == fend)) {
      //  Clear range did not change.
      fprintf(logFile, F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\tNOC%s\n",
              iid,
              ibgn, iend,
              fbgn, fend,
              (logMsg[0] == 0) ? "" : logMsg);
      continue;
    }

    //  Clear range changed, and we like it!

    if (doModify) {
      fr.gkFragment_setClearRegion(fbgn, fend, AS_READ_CLEAR_OBTMERGE);
      gkpStore->gkStore_setFragment(&fr);
    }

    fprintf(logFile, F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\tMOD%s\n",
            iid,
            ibgn, iend,
            fbgn, fend,
            (logMsg[0] == 0) ? "" : logMsg);
  }

  delete gkpStore;
  delete ovlPrimary;
  delete ovlSecondary;

  fclose(logFile);
  fclose(sumFile);

  exit(0);
}


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

const char *mainid = "$Id: finalTrim.C,v 1.5 2012/03/07 05:51:54 brianwalenz Exp $";

#include "finalTrim.H"


uint32
loadOverlaps(uint32         iid,
             OVSoverlap   *&ovl,
             uint32        &ovlLen,
             uint32        &ovlMax,
             OverlapStore  *ovlPrimary,
             OverlapStore  *ovlSecondary) {

  //  Allocate initial space

  if (ovl == NULL) {
    ovlLen = 0;
    ovlMax = 65 * 1024;
    ovl    = new OVSoverlap [ovlMax];
  }

  //  Return if the overlaps are for the current read, or some future read.

  if (iid <= ovl[0].a_iid)
    return(ovlLen);

  //  Until we load the correct overlap, repeat.

  do {
    //  Count the number of overlaps to load
    ovlLen  = 0;
    ovlLen += AS_OVS_readOverlapsFromStore(ovlPrimary,   NULL, 0, AS_OVS_TYPE_ANY);
    ovlLen += AS_OVS_readOverlapsFromStore(ovlSecondary, NULL, 0, AS_OVS_TYPE_ANY);

    //  Quit now if there are no overlaps.  This simplifies the rest of the loop.
    if (ovlLen == 0)
      return(0);

    //  Allocate space for these overlaps.
    while (ovlMax < ovlLen) {
      ovlMax *= 2;
      delete [] ovl;
      ovl = new OVSoverlap [ovlMax];
    }

    //  Load the overlaps
    ovlLen  = 0;
    ovlLen += AS_OVS_readOverlapsFromStore(ovlPrimary,   ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
    ovlLen += AS_OVS_readOverlapsFromStore(ovlSecondary, ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);

    //fprintf(stderr, "LOADED %d overlaps for a_iid %d\n", ovlLen, ovl[0].a_iid);

    //  If we read overlaps for a fragment after 'iid', we're done.  The client will properly save
    //  these overlaps until the iid becomes active.  And just in case it doesn't, we return above
    //  if the iid passed in is less than the current overlap.
    //
    if (iid <= ovl[0].a_iid)
      return(ovlLen);

    //  On the otherhand, if we read overlaps for a fragment before 'iid', we can either keep reading
    //  until we find the overlaps for this fragment, or jump to the correct spot to read overlaps.
    //
    //  The rule is simple.  If we're within 50 of the correct IID, keep streaming.  Otherwise, make
    //  a jump.  AS_OVS_setRangeOverlapStore() seems to ALWAYS close and open a file, which is somewhat
    //  expensive, especially if the file doesn't actually change.
    //
    if (50 < iid - ovl[0].a_iid) {
      AS_OVS_setRangeOverlapStore(ovlPrimary,   iid, UINT32_MAX);
      AS_OVS_setRangeOverlapStore(ovlSecondary, iid, UINT32_MAX);
    }
  } while (ovl[0].a_iid < iid);

  return(ovlLen);
}


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

  assert(mbgn < mend);
  assert(fbgn < fend);

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

  argc = AS_configure(argc, argv);

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

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-n") == 0) {
      doModify = false;

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
    fprintf(stderr, "   -n             do not modify\n");
    exit(1);
  }


  sprintf(logName, "%s.log",     outputPrefix);
  sprintf(sumName, "%s.summary", outputPrefix);

  logFile = fopen(logName, "w");
  sumFile = fopen(sumName, "w");

  fprintf(logFile, "iid\tinitL\tinitR\tfinalL\tfinalR\tmessage (DEL=deleted NOC=no change MOD=modified)\n");

  uint32      ovlLen       = 0;
  uint32      ovlMax       = 64 * 1024;
  OVSoverlap *ovl          = new OVSoverlap [ovlMax];
  bool        ovlVal       = false;

  memset(ovl, 0, sizeof(OVSoverlap) * ovlMax);

  gkFragment  fr;
  char        logMsg[1024] = {0};

  for (uint32 iid=1; iid<=gkpStore->gkStore_getNumFragments(); iid++) {

    //  The QLT is needed ONLY for Sanger reads, and is a total waste of bandwidth for all other
    //  read types.  Perhaps this can be moved into initialTrim -- compute the strict QV trim there
    //  too, save it in the OBTMERGE clear range
    //
    gkpStore->gkStore_getFragment(iid, &fr, GKFRAGMENT_QLT);

    logMsg[0] = 0;

    uint32      lid = fr.gkFragment_getLibraryIID();
    gkLibrary  *lb  = gkpStore->gkStore_getLibrary(lid);

    //fprintf(stderr, "FRAG %d starts.\n", iid);

    //  If the fragment is deleted, do nothing.  If the fragment was deleted AFTER overlaps were
    //  generated, then the overlaps will be out of sync -- we'll get overlaps for these fragments
    //  we skip.
    if (fr.gkFragment_getIsDeleted() == 1) {
      //fprintf(stderr, "FRAG %d is deleted.\n", iid);
      continue;
    }

    //  If it did not request trimming, do nothing.  Similar to the above, we'll get overlaps to
    //  fragments we skip.
    if ((lb->doTrim_finalLargestCovered == false) &&
        (lb->doTrim_finalEvidenceBased  == false)) {
      //fprintf(stderr, "FRAG %d didn't request trimming.\n", iid);
      continue;
    }

    uint32      ibgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_OBTINITIAL);
    uint32      iend = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_OBTINITIAL);

    //  Is the clear range valid?  If we skip initial trim for this read, there is no clear range defined.
    if ((ibgn == 1) && (iend == 0)) {
      ibgn = 0;
      iend = fr.gkFragment_getSequenceLength();
    }

    uint32      fbgn = ibgn;
    uint32      fend = iend;


    //  Load overlaps, unless we have already loaded them.
    if (ovl[0].a_iid < iid)
      loadOverlaps(iid, ovl, ovlLen, ovlMax, ovlPrimary, ovlSecondary);

    bool isGood = false;

    //  If there are no overlaps for this fragment, do nothing.
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

    } else if        (lb->doTrim_finalLargestCovered == true) {
      //  Use the largest region covered by overlaps as the trim
      isGood = largestCovered(ovl, ovlLen,
                              fr,
                              ibgn, iend, fbgn, fend,
                              logMsg,
                              errorRate,
                              errorLimit,
                              lb->doTrim_initialQualityBased,
                              AS_OVERLAP_MIN_LEN,
                              2);
      assert(fbgn <= fend);

    } else if (lb->doTrim_finalEvidenceBased == true) {
      //  Do Sanger-style heuristics
      isGood = evidenceBased(ovl, ovlLen,
                             fr,
                             ibgn, iend, fbgn, fend,
                             logMsg,
                             errorRate,
                             errorLimit,
                             lb->doTrim_initialQualityBased);
      assert(fbgn <= fend);

    } else {
      //  Do nothing.  Really shouldn't get here.
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

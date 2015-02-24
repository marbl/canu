
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

#include "finalTrim.H"
#include "clearRangeFile.H"

#include "AS_OBT_overlaps.H"
#include "AS_UTL_decodeRange.H"

bool
trimWithoutOverlaps(ovsOverlap  *ovl,
                    uint32       ovlLen,
                    gkRead      *read,
                    gkReadData  *readData,
                    uint32       ibgn,
                    uint32       iend,
                    uint32      &fbgn,
                    uint32      &fend,
                    char        *logMsg,
                    bool         qvTrimAllowed) {

  fbgn = ibgn;
  fend = iend;

  if (qvTrimAllowed == false) {
    strcpy(logMsg, "\tno overlaps, qv trim not allowed");
    return(true);
  }

  //gkp->gkStore_getReadData(read, readData);

  doTrim(read,
         readData,
         qual.lookupNumber(20),
         fbgn, fend);

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
enforceMaximumClearRange(ovsOverlap      *ovl,
                         uint32           ovlLen,
                         gkRead          *read,
                         uint32           ibgn,
                         uint32           iend,
                         uint32          &fbgn,
                         uint32          &fend,
                         char            *logMsg,
                         clearRangeFile  *maxClr) {

  if (maxClr == NULL)
    return(true);

  if (maxClr->isUndefined(read->gkRead_readID()) == true)
    return(true);

  if (fbgn == fend)
    return(true);

  uint32 mbgn = maxClr->bgn(read->gkRead_readID());
  uint32 mend = maxClr->end(read->gkRead_readID());

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
  char             *gkpName = 0L;
  char             *ovsName = 0L;

  char             *iniClrName = NULL;
  char             *maxClrName = NULL;
  char             *outClrName = NULL;

  uint32            errorRate    = AS_OVS_encodeQuality(0.015);

  char             *outputPrefix  = NULL;
  char              logName[FILENAME_MAX] = {0};
  char              sumName[FILENAME_MAX] = {0};
  FILE             *logFile = 0L;
  FILE             *sumFile = 0L;

  bool              doModify = true;  //  Make this false for testing

  uint32            idMin = 1;
  uint32            idMax = UINT32_MAX;

  argc = AS_configure(argc, argv);

  uint32            minEvidenceOverlap  = MIN(500, AS_OVERLAP_MIN_LEN);
  uint32            minEvidenceCoverage = 1;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovsName = argv[++arg];

    } else if (strcmp(argv[arg], "-Ci") == 0) {
      iniClrName = argv[++arg];
    } else if (strcmp(argv[arg], "-Cm") == 0) {
      maxClrName = argv[++arg];
    } else if (strcmp(argv[arg], "-Co") == 0) {
      outClrName = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      double erate = atof(argv[++arg]);
      if (erate > AS_MAX_ERROR_RATE)
        fprintf(stderr, "ERROR:  Error rate (-e) %s too large; must be 'fraction error' and below %f\n", argv[arg], AS_MAX_ERROR_RATE), exit(1);
      errorRate = AS_OVS_encodeQuality(erate);

    } else if (strcmp(argv[arg], "-ol") == 0) {
      minEvidenceOverlap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-oc") == 0) {
      minEvidenceCoverage = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-n") == 0) {
      doModify = false;

    } else if (strcmp(argv[arg], "-t") == 0) {
      AS_UTL_decodeRange(argv[++arg], idMin, idMax);

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if ((gkpName       == NULL) ||
      (ovsName       == NULL) ||
      (outputPrefix  == NULL) ||
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

  gkStore *gkp = new gkStore(gkpName);
  ovStore *ovs = new ovStore(ovsName);

  clearRangeFile   *iniClr = (iniClrName == NULL) ? NULL : new clearRangeFile(iniClrName, gkp->gkStore_getNumReads());
  clearRangeFile   *maxClr = (maxClrName == NULL) ? NULL : new clearRangeFile(maxClrName, gkp->gkStore_getNumReads());
  clearRangeFile   *outClr = (outClrName == NULL) ? NULL : new clearRangeFile(outClrName, gkp->gkStore_getNumReads());


  sprintf(logName, "%s.log",     outputPrefix);
  sprintf(sumName, "%s.summary", outputPrefix);

  logFile = fopen(logName, "w");
  sumFile = fopen(sumName, "w");

  fprintf(logFile, "id\tinitL\tinitR\tfinalL\tfinalR\tmessage (DEL=deleted NOC=no change MOD=modified)\n");

  fprintf(sumFile, "Overlap error rate     <= %.2f%% error\n", 100.0 * AS_OVS_decodeQuality(errorRate));
  fprintf(sumFile, "Overlap min overlap    >= %u base%s (for 'largest covered')\n", minEvidenceOverlap,  (minEvidenceOverlap  == 1) ? "" : "s");
  fprintf(sumFile, "Overlap min coverage   >= %u read%s (for 'largest covered')\n", minEvidenceCoverage, (minEvidenceCoverage == 1) ? "" : "s");


  uint32      ovlLen       = 0;
  uint32      ovlMax       = 64 * 1024;
  ovsOverlap *ovl          = new ovsOverlap [ovlMax];

  memset(ovl, 0, sizeof(ovsOverlap) * ovlMax);

  char        logMsg[1024] = {0};

  if (idMin < 1)
    idMin = 1;
  if (idMax > gkp->gkStore_getNumReads())
    idMax = gkp->gkStore_getNumReads();

  fprintf(stderr, "Processing from ID "F_U32" to "F_U32" out of "F_U32" reads.\n",
          idMin,
          idMax,
          gkp->gkStore_getNumReads());

  gkReadData  readData;

  for (uint32 id=idMin; id<=idMax; id++) {

    //  The QLT is needed ONLY for Sanger reads, and is a total waste of bandwidth for all other
    //  read types.  Perhaps this can be moved into initialTrim -- compute the strict QV trim there
    //  too, save it in the OBTMERGE clear range
    //
    gkRead *read = gkp->gkStore_getRead(id);

    logMsg[0] = 0;

    gkLibrary  *lb  = gkp->gkStore_getLibrary(read->gkRead_libraryID());

    //  If the fragment is deleted, do nothing.  If the fragment was deleted AFTER overlaps were
    //  generated, then the overlaps will be out of sync -- we'll get overlaps for these fragments
    //  we skip.
    //
    if ((iniClr) && (iniClr->isDeleted(id) == true))
      continue;

    //  If it did not request trimming, do nothing.  Similar to the above, we'll get overlaps to
    //  fragments we skip.
    //
    if ((lb->gkLibrary_finalTrim() == FINALTRIM_LARGEST_COVERED) &&
        //(lb->gkLibrary_finalTrim() == FINALTRIM_EVIDENCE_BASED) &&
        (lb->gkLibrary_finalTrim() == FINALTRIM_BEST_EDGE))
      continue;

    uint32      ibgn = iniClr->bgn(id);
    uint32      iend = iniClr->end(id);

    //  Is the clear range valid?  If we skip initial trim for this read, there is no clear range defined.
    if ((ibgn == 1) && (iend == 0)) {
      ibgn = 0;
      iend = read->gkRead_sequenceLength();
    }

    uint32      fbgn = ibgn;
    uint32      fend = iend;

    loadOverlaps(id, ovl, ovlLen, ovlMax, ovs, NULL);

    bool isGood = false;

    //  If there are no overlaps for this read, do nothing.
    if ((id < ovl[0].a_iid) ||
        (ovlLen == 0)) {
      isGood = trimWithoutOverlaps(ovl, ovlLen,
                                   read,
                                   &readData,
                                   ibgn, iend, fbgn, fend,
                                   logMsg,
                                   lb->gkLibrary_initialTrim() == INITIALTRIM_QUALITY_BASED);
      assert(fbgn <= fend);

    }

    //  Use the largest region covered by overlaps as the trim
    else if        (lb->gkLibrary_finalTrim() == FINALTRIM_LARGEST_COVERED) {
      isGood = largestCovered(gkp,
                              ovl, ovlLen,
                              read,
                              ibgn, iend, fbgn, fend,
                              logMsg,
                              errorRate,
                              lb->gkLibrary_initialTrim() == INITIALTRIM_QUALITY_BASED,
                              minEvidenceOverlap,
                              minEvidenceCoverage);
      assert(fbgn <= fend);

    }

    //  Use the largest region covered by overlaps as the trim
    else if        (lb->gkLibrary_finalTrim() == FINALTRIM_BEST_EDGE) {
      isGood = bestEdge(gkp,
                        ovl, ovlLen,
                        read,
                        ibgn, iend, fbgn, fend,
                        logMsg,
                        errorRate,
                        lb->gkLibrary_initialTrim() == INITIALTRIM_QUALITY_BASED,
                        minEvidenceOverlap,
                        minEvidenceCoverage);
      assert(fbgn <= fend);

    }

#if 0
    //  Do Sanger-style heuristics
    else if (lb->gkLibrary_finalTrim() == FINALTRIM_EVIDENCE_BASED) {
      isGood = evidenceBased(gkp,
                             ovl, ovlLen,
                             read,
                             ibgn, iend, fbgn, fend,
                             logMsg,
                             errorRate,
                             lb->gkLibrary_initialTrim() == INITIALTRIM_QUALITY_BASED);
      assert(fbgn <= fend);
    }
#endif

    //  Do nothing.  Really shouldn't get here.
    else {
      assert(0);
      continue;
    }

    //  Enforce the maximum clear range

    if (isGood)
      isGood = enforceMaximumClearRange(ovl, ovlLen,
                                        read,
                                        ibgn, iend, fbgn, fend,
                                        logMsg,
                                        maxClr);

    assert(fbgn <= fend);

    //  Bad trimming?  Invalid clear range?  Too small?

    if ((isGood == false) ||
        (fbgn > fend) ||
        (fend - fbgn < AS_READ_MIN_LEN)) {

      assert(fbgn <= fend);

      if (doModify) {
        outClr->setbgn(id) = fbgn;
        outClr->setend(id) = fend;
        outClr->setDeleted(id);
      }

      fprintf(logFile, F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\tDEL%s\n",
              id,
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
              id,
              ibgn, iend,
              fbgn, fend,
              (logMsg[0] == 0) ? "" : logMsg);
      continue;
    }

    //  Clear range changed, and we like it!

    if (doModify) {
        outClr->setbgn(id) = fbgn;
        outClr->setend(id) = fend;
    }

    fprintf(logFile, F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\tMOD%s\n",
            id,
            ibgn, iend,
            fbgn, fend,
            (logMsg[0] == 0) ? "" : logMsg);
  }

  delete gkp;
  delete ovs;

  fclose(logFile);
  fclose(sumFile);

  exit(0);
}

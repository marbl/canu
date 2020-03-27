
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "trimReads.H"
#include "trimStat.H"
#include "clearRangeFile.H"

#include "strings.H"




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
enforceMaximumClearRange(uint32           readID,
                         uint32    UNUSED(ibgn),
                         uint32    UNUSED(iend),
                         uint32          &fbgn,
                         uint32          &fend,
                         char            *logMsg,
                         clearRangeFile  *maxClr) {

  if (maxClr == NULL)
    return(true);

  if (fbgn == fend)
    return(true);

  uint32 mbgn = maxClr->bgn(readID);
  uint32 mend = maxClr->end(readID);

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
    fbgn = max(fbgn, mbgn);
    fend = min(fend, mend);

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
  char       *seqName = 0L;
  char       *ovsName = 0L;

  char       *iniClrName = NULL;
  char       *maxClrName = NULL;
  char       *outClrName = NULL;

  uint32      errorValue     = AS_OVS_encodeEvalue(0.015);
  uint32      minAlignLength = 40;
  uint32      minReadLength  = 64;

  char       *outputPrefix  = NULL;
  char        logName[FILENAME_MAX] = {0};
  char        sumName[FILENAME_MAX] = {0};
  FILE       *logFile = 0L;
  FILE       *staFile = 0L;

  uint32      idMin = 1;
  uint32      idMax = UINT32_MAX;

  uint32      minEvidenceOverlap  = 40;
  uint32      minEvidenceCoverage = 1;

  //  Statistics on the trimming

  trimStat    readsIn;      //  Read is eligible for trimming
  trimStat    deletedIn;    //  Read was deleted already
  trimStat    noTrimIn;     //  Read not requesting trimming

  trimStat    readsOut;     //  Read was trimmed to a valid read
  trimStat    noOvlOut;     //  Read was deleted; no ovelaps
  trimStat    deletedOut;   //  Read was deleted; too small after trimming
  trimStat    noChangeOut;  //  Read was untrimmed

  trimStat    trim5;        //  Bases trimmed from the 5' end
  trimStat    trim3;


  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

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
      errorValue = AS_OVS_encodeEvalue(erate);

    } else if (strcmp(argv[arg], "-l") == 0) {
      minAlignLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-minlength") == 0) {
      minReadLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-ol") == 0) {
      minEvidenceOverlap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-oc") == 0) {
      minEvidenceCoverage = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      decodeRange(argv[++arg], idMin, idMax);

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if ((seqName       == NULL) ||
      (ovsName       == NULL) ||
      (outClrName    == NULL) ||
      (outputPrefix  == NULL) ||
      (err)) {
    fprintf(stderr, "usage: %s -S seqStore -O ovlStore -Co output.clearFile -o outputPrefix\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S seqStore    path to read store\n");
    fprintf(stderr, "  -O ovlStore    path to overlap store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o name        output prefix, for logging\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t bgn-end     limit processing to only reads from bgn to end (inclusive)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -Ci clearFile  path to input clear ranges (NOT SUPPORTED)\n");
    //fprintf(stderr, "  -Cm clearFile  path to maximal clear ranges\n");
    fprintf(stderr, "  -Co clearFile  path to ouput clear ranges\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e erate       ignore overlaps with more than 'erate' percent error\n");
    //fprintf(stderr, "  -l length      ignore overlaps shorter than 'l' aligned bases (NOT SUPPORTED)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -ol l          the minimum evidence overlap length\n");
    fprintf(stderr, "  -oc c          the minimum evidence overlap coverage\n");
    fprintf(stderr, "                   evidence overlaps must overlap by 'l' bases to be joined, and\n");
    fprintf(stderr, "                   must be at least 'c' deep to be retained\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -minlength l   reads trimmed below this many bases are deleted\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  sqStore          *seq = new sqStore(seqName);
  ovStore          *ovs = new ovStore(ovsName, seq);

  clearRangeFile   *iniClr = (iniClrName == NULL) ? NULL : new clearRangeFile(iniClrName, seq);
  clearRangeFile   *maxClr = (maxClrName == NULL) ? NULL : new clearRangeFile(maxClrName, seq);
  clearRangeFile   *outClr =                               new clearRangeFile(outClrName, seq);

  if (outClr)
    //  If the outClr file exists, those clear ranges are loaded.  We need to reset them
    //  back to 'untrimmed' for now.
    outClr->reset(seq);

  if (iniClr && outClr)
    //  An iniClr file was supplied, so use those as the initial clear ranges.
    outClr->copy(iniClr);


  if (outputPrefix) {
    snprintf(logName, FILENAME_MAX, "%s.log",   outputPrefix);

    logFile = AS_UTL_openOutputFile(logName);

    fprintf(logFile, "id\tinitL\tinitR\tfinalL\tfinalR\tmessage (DEL=deleted NOC=no change MOD=modified)\n");
  }


  uint32      ovlLen       = 0;
  uint32      ovlMax       = 0;
  ovOverlap  *ovl          = NULL;

  char        logMsg[1024] = {0};

  if (idMin < 1)
    idMin = 1;
  if (idMax > seq->sqStore_lastReadID())
    idMax = seq->sqStore_lastReadID();

  fprintf(stderr, "Processing from ID " F_U32 " to " F_U32 " out of " F_U32 " reads.\n",
          idMin,
          idMax,
          seq->sqStore_lastReadID());


  for (uint32 id=idMin; id<=idMax; id++) {
    sqLibrary  *libr = seq->sqStore_getLibraryForRead(id);

    logMsg[0] = 0;

    //  If the fragment is deleted, do nothing.  If the fragment was deleted AFTER overlaps were
    //  generated, then the overlaps will be out of sync -- we'll get overlaps for these fragments
    //  we skip.
    //
    if ((iniClr) && (iniClr->isDeleted(id) == true)) {
      deletedIn += seq->sqStore_getReadLength(id);
      continue;
    }

    //  If it did not request trimming, do nothing.  Similar to the above, we'll get overlaps to
    //  fragments we skip.
    //
#if 0
    //  (yes, this is nonsense)
    if ((libr->sqLibrary_finalTrim() == SQ_FINALTRIM_LARGEST_COVERED) &&
        (libr->sqLibrary_finalTrim() == SQ_FINALTRIM_BEST_EDGE)) {
      noTrimIn += seq->sqStore_getReadLength(id);
      continue;
    }
#endif

    readsIn += seq->sqStore_getReadLength(id);


    //  Decide on the initial trimming.  We copied any iniClr into outClr above, and if there wasn't
    //  an iniClr, then outClr is the full read.

    uint32      ibgn   = outClr->bgn(id);
    uint32      iend   = outClr->end(id);

    //  Set the, ahem, initial final trimming.

    bool        isGood = false;
    uint32      fbgn   = ibgn;
    uint32      fend   = iend;

    //  Load overlaps.

    ovlLen = ovs->loadOverlapsForRead(id, ovl, ovlMax);

    //  Trim!

    //  No overlaps, so mark it as junk.
    if (ovlLen == 0) {
      isGood = false;
    }

    //  Use the largest region covered by overlaps as the trim
    else {

      assert(ovlLen > 0);
      assert(id == ovl[0].a_iid);

      isGood = largestCovered(ovl, ovlLen,
                              id, seq->sqStore_getReadLength(id),
                              ibgn, iend, fbgn, fend,
                              logMsg,
                              errorValue,
                              minEvidenceOverlap,
                              minEvidenceCoverage,
                              minReadLength);
      assert(fbgn <= fend);
    }

#if 0
    //  Use the largest region covered by overlaps as the trim
    else if (libr->sqLibrary_finalTrim() == SQ_FINALTRIM_BEST_EDGE) {

      assert(ovlLen > 0);
      assert(id == ovl[0].a_iid);

      isGood = bestEdge(ovl, ovlLen,
                        id, seq->sqStore_getReadLength(id),
                        ibgn, iend, fbgn, fend,
                        logMsg,
                        errorValue,
                        minEvidenceOverlap,
                        minEvidenceCoverage,
                        minReadLength);
      assert(fbgn <= fend);
    }

    //  Do nothing.  Really shouldn't get here.
    else {
      assert(0);
      continue;
    }
#endif

    //  Enforce the maximum clear range

    if ((isGood) && (maxClr)) {
      isGood = enforceMaximumClearRange(id,
                                        ibgn, iend, fbgn, fend,
                                        logMsg,
                                        maxClr);
      assert(fbgn <= fend);
    }

    //
    //  Trimmed.  Make sense of the result, write some logs, and update the output.
    //


    //  If bad trimming or too small, write the log and keep going.
    //
    if (ovlLen == 0) {
      noOvlOut += seq->sqStore_getReadLength(id);

      outClr->setbgn(id) = fbgn;
      outClr->setend(id) = fend;
      outClr->setDeleted(id);  //  Gah, just obliterates the clear range.

      fprintf(logFile, F_U32"\t" F_U32 "\t" F_U32 "\t" F_U32 "\t" F_U32 "\tNOV%s\n",
              id,
              ibgn, iend,
              fbgn, fend,
              (logMsg[0] == 0) ? "" : logMsg);
    }

    else if ((isGood == false) || (fend - fbgn < minReadLength)) {
      deletedOut += seq->sqStore_getReadLength(id);

      outClr->setbgn(id) = fbgn;
      outClr->setend(id) = fend;
      outClr->setDeleted(id);  //  Gah, just obliterates the clear range.

      fprintf(logFile, F_U32"\t" F_U32 "\t" F_U32 "\t" F_U32 "\t" F_U32 "\tDEL%s\n",
              id,
              ibgn, iend,
              fbgn, fend,
              (logMsg[0] == 0) ? "" : logMsg);
    }

    //  If we didn't change anything, also write a log.
    //
    else if ((ibgn == fbgn) &&
             (iend == fend)) {
      noChangeOut += seq->sqStore_getReadLength(id);

      fprintf(logFile, F_U32"\t" F_U32 "\t" F_U32 "\t" F_U32 "\t" F_U32 "\tNOC%s\n",
              id,
              ibgn, iend,
              fbgn, fend,
              (logMsg[0] == 0) ? "" : logMsg);
      continue;
    }

    //  Otherwise, we actually did something.

    else {
      readsOut += fend - fbgn;

      outClr->setbgn(id) = fbgn;
      outClr->setend(id) = fend;

      assert(ibgn <= fbgn);
      assert(fend <= iend);

      if (fbgn - ibgn > 0)   trim5 += fbgn - ibgn;
      if (iend - fend > 0)   trim3 += iend - fend;

      fprintf(logFile, F_U32"\t" F_U32 "\t" F_U32 "\t" F_U32 "\t" F_U32 "\tMOD%s\n",
              id,
              ibgn, iend,
              fbgn, fend,
              (logMsg[0] == 0) ? "" : logMsg);
    }
  }

  //  Clean up.

  delete seq;

  delete [] ovl;
  delete    ovs;

  delete    iniClr;
  delete    maxClr;
  delete    outClr;

  AS_UTL_closeFile(logFile, logName);

  //  should fprintf() the numbers directly here so an explanation of each category can be supplied;
  //  simpler for now to have report() do it.

  //  Dump the statistics and plots

  if (outputPrefix) {
    snprintf(sumName, FILENAME_MAX, "%s.stats", outputPrefix);

    staFile = AS_UTL_openOutputFile(sumName);
  }

  if (staFile == NULL)
    staFile = stderr;

  fprintf(staFile, "PARAMETERS:\n");
  fprintf(staFile, "----------\n");
  fprintf(staFile, "%7u    (reads trimmed below this many bases are deleted)\n", minReadLength);
  fprintf(staFile, "%7.4f    (use overlaps at or below this fraction error)\n", AS_OVS_decodeEvalue(errorValue));
  fprintf(staFile, "%7u    (break region if overlap is less than this long, for 'largest covered' algorithm)\n", minEvidenceOverlap);
  fprintf(staFile, "%7u    (break region if overlap coverage is less than this many read%s, for 'largest covered' algorithm)\n", minEvidenceCoverage, (minEvidenceCoverage == 1) ? "" : "s");
  fprintf(staFile, "\n");

  fprintf(staFile, "INPUT READS:\n");
  fprintf(staFile, "-----------\n");
  fprintf(staFile, "%6" F_U32P " reads %12" F_U64P " bases (reads processed)\n", readsIn.nReads,  readsIn.nBases);
  fprintf(staFile, "%6" F_U32P " reads %12" F_U64P " bases (reads not processed, previously deleted)\n", deletedIn.nReads, deletedIn.nBases);
  fprintf(staFile, "%6" F_U32P " reads %12" F_U64P " bases (reads not processed, in a library where trimming isn't allowed)\n", noTrimIn.nReads, noTrimIn.nBases);

  readsIn  .generatePlots(outputPrefix, "inputReads",        250);
  deletedIn.generatePlots(outputPrefix, "inputDeletedReads", 250);
  noTrimIn .generatePlots(outputPrefix, "inputNoTrimReads",  250);

  fprintf(staFile, "\n");
  fprintf(staFile, "OUTPUT READS:\n");
  fprintf(staFile, "------------\n");
  fprintf(staFile, "%6" F_U32P " reads %12" F_U64P " bases (trimmed reads output)\n", readsOut.nReads,    readsOut.nBases);
  fprintf(staFile, "%6" F_U32P " reads %12" F_U64P " bases (reads with no change, kept as is)\n", noChangeOut.nReads, noChangeOut.nBases);
  fprintf(staFile, "%6" F_U32P " reads %12" F_U64P " bases (reads with no overlaps, deleted)\n", noOvlOut.nReads,    noOvlOut.nBases);
  fprintf(staFile, "%6" F_U32P " reads %12" F_U64P " bases (reads with short trimmed length, deleted)\n", deletedOut.nReads,  deletedOut.nBases);

  readsOut   .generatePlots(outputPrefix, "outputTrimmedReads",   250);
  noOvlOut   .generatePlots(outputPrefix, "outputNoOvlReads",     250);
  deletedOut .generatePlots(outputPrefix, "outputDeletedReads",   250);
  noChangeOut.generatePlots(outputPrefix, "outputUnchangedReads", 250);

  fprintf(staFile, "\n");
  fprintf(staFile, "TRIMMING DETAILS:\n");
  fprintf(staFile, "----------------\n");
  fprintf(staFile, "%6" F_U32P " reads %12" F_U64P " bases (bases trimmed from the 5' end of a read)\n", trim5.nReads, trim5.nBases);
  fprintf(staFile, "%6" F_U32P " reads %12" F_U64P " bases (bases trimmed from the 3' end of a read)\n", trim3.nReads, trim3.nBases);

  trim5.generatePlots(outputPrefix, "trim5", 25);
  trim3.generatePlots(outputPrefix, "trim3", 25);

  AS_UTL_closeFile(staFile, sumName);

  //  Buh-bye.

  exit(0);
}

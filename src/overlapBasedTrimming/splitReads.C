
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-APR-15 to 2015-JUN-25
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "splitReads.H"
#include "trimStat.H"
#include "clearRangeFile.H"

#include "AS_UTL_decodeRange.H"


int
main(int argc, char **argv) {
  char     *gkpName = NULL;
  char     *ovsName = NULL;

  char     *finClrName = NULL;
  char     *outClrName = NULL;

  double    errorRate       = 0.06;
  //uint32    minAlignLength  = 40;
  uint32    minReadLength   = 64;

  uint32    idMin = 1;
  uint32    idMax = UINT32_MAX;

  char     *outputPrefix = NULL;
  char      outputName[FILENAME_MAX];

  FILE     *staFile      = NULL;
  FILE     *reportFile   = NULL;
  FILE     *subreadFile  = NULL;

  bool      doSubreadLogging        = true;
  bool      doSubreadLoggingVerbose = false;

  //  Statistics on the trimming - the second set are from the old logging, and don't really apply anymore.

  trimStat  readsIn;                  //  Read is eligible for trimming
  trimStat  deletedIn;                //  Read was deleted already
  trimStat  noTrimIn;                 //  Read not requesting trimming

  trimStat  noOverlaps;               //  no overlaps in store
  trimStat  noCoverage;               //  no coverage after adjusting for trimming done

  trimStat  readsProcChimera;         //  Read was processed for chimera signal
  trimStat  readsProcSpur;            //  Read was processed for spur signal
  trimStat  readsProcSubRead;         //  Read was processed for subread signal

#if 0
  trimStat  badSpur5;
  trimStat  badSpur3;
  trimStat  badChimera;
  trimStat  badSubread;
#endif

  trimStat  readsNoChange;

  trimStat  readsBadSpur5,   basesBadSpur5;
  trimStat  readsBadSpur3,   basesBadSpur3;
  trimStat  readsBadChimera, basesBadChimera;
  trimStat  readsBadSubread, basesBadSubread;

  trimStat  readsTrimmed5;
  trimStat  readsTrimmed3;

#if 0
  trimStat  fullCoverage;             //  fully covered by overlaps
  trimStat  noSignalNoGap;            //  no signal, no gaps
  trimStat  noSignalButGap;           //  no signal, with gaps

  trimStat  bothFixed;                //  both chimera and spur signal trimmed
  trimStat  chimeraFixed;             //  only chimera signal trimmed
  trimStat  spurFixed;                //  only spur signal trimmed

  trimStat  bothDeletedSmall;         //  deleted because of both cimera and spur signals
  trimStat  chimeraDeletedSmall;      //  deleted because of chimera signal
  trimStat  spurDeletedSmall;         //  deleted because of spur signal

  trimStat  spurDetectedNormal;       //  normal spur detected
  trimStat  spurDetectedLinker;       //  linker spur detected

  trimStat  chimeraDetectedInnie;     //  innpue-pair chimera detected
  trimStat  chimeraDetectedOverhang;  //  overhanging chimera detected
  trimStat  chimeraDetectedGap;       //  gap chimera detected
  trimStat  chimeraDetectedLinker;    //  linker chimera detected
#endif

  trimStat  deletedOut;               //  Read was deleted by trimming

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovsName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      AS_UTL_decodeRange(argv[++arg], idMin, idMax);

    } else if (strcmp(argv[arg], "-Ci") == 0) {
      finClrName = argv[++arg];
    } else if (strcmp(argv[arg], "-Co") == 0) {
      outClrName = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      errorRate = atof(argv[++arg]);

    //} else if (strcmp(argv[arg], "-l") == 0) {
    //  minAlignLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-minlength") == 0) {
      minReadLength = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }

  if (errorRate < 0.0)
    err++;

  if ((gkpName == 0L) || (ovsName == 0L) || (outputPrefix == NULL) || (err)) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore -Ci input.clearFile -Co output.clearFile -o outputPrefix]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G gkpStore    path to read store\n");
    fprintf(stderr, "  -O ovlStore    path to overlap store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o name        output prefix, for logging\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t bgn-end     limit processing to only reads from bgn to end (inclusive)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -Ci clearFile  path to input clear ranges (NOT SUPPORTED)\n");
    fprintf(stderr, "  -Co clearFile  path to ouput clear ranges\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e erate       ignore overlaps with more than 'erate' percent error\n");
    //fprintf(stderr, "  -l length      ignore overlaps shorter than 'l' aligned bases (NOT SUPPORTED)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -minlength l   reads trimmed below this many bases are deleted\n");
    fprintf(stderr, "\n");

    if (errorRate < 0.0)
      fprintf(stderr, "ERROR: Error rate (-e) value %f too small; must be 'fraction error' and above 0.0\n", errorRate);

    exit(1);
  }

  gkStore         *gkp = gkStore::gkStore_open(gkpName);
  ovStore         *ovs = new ovStore(ovsName, gkp);

  clearRangeFile  *finClr = new clearRangeFile(finClrName, gkp);
  clearRangeFile  *outClr = new clearRangeFile(outClrName, gkp);

  if (outClr)
    //  If the outClr file exists, those clear ranges are loaded.  We need to reset them
    //  back to 'untrimmed' for now.
    outClr->reset(gkp);

  if (finClr && outClr)
    //  A finClr file was supplied, so use those as the clear ranges.
    outClr->copy(finClr);


  sprintf(outputName, "%s.log",         outputPrefix);
  errno = 0;
  reportFile  = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", outputName, strerror(errno)), exit(1);

  sprintf(outputName, "%s.subread.log", outputPrefix);
  errno = 0;
  subreadFile = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", outputName, strerror(errno)), exit(1);


  uint32      ovlLen = 0;
  uint32      ovlMax = 64 * 1024;
  ovOverlap  *ovl    = ovOverlap::allocateOverlaps(gkp, ovlMax);

  memset(ovl, 0, sizeof(ovOverlap) * ovlMax);

  workUnit *w = new workUnit;


  if (idMin < 1)
    idMin = 1;
  if (idMax > gkp->gkStore_getNumReads())
    idMax = gkp->gkStore_getNumReads();

  fprintf(stderr, "Processing from ID "F_U32" to "F_U32" out of "F_U32" reads, using errorRate = %.2f\n",
          idMin,
          idMax,
          gkp->gkStore_getNumReads(),
          errorRate);

  for (uint32 id=idMin; id<=idMax; id++) {
    gkRead     *read = gkp->gkStore_getRead(id);
    gkLibrary  *libr = gkp->gkStore_getLibrary(read->gkRead_libraryID());

    if (finClr->isDeleted(id)) {
      //  Read already trashed.
      deletedIn += read->gkRead_sequenceLength();
      continue;
    }

    if ((libr->gkLibrary_removeSpurReads()     == false) &&
        (libr->gkLibrary_removeChimericReads() == false) &&
        (libr->gkLibrary_checkForSubReads()    == false)) {
      //  Nothing to do.
      noTrimIn += read->gkRead_sequenceLength();
      continue;
    }

    readsIn += read->gkRead_sequenceLength();


    uint32   nLoaded = ovs->readOverlaps(id, ovl, ovlLen, ovlMax);

    //fprintf(stderr, "read %7u with %7u overlaps\r", id, nLoaded);

    if (nLoaded == 0) {
      //  No overlaps, nothing to check!
      noOverlaps += read->gkRead_sequenceLength();
      continue;
    }

    w->clear(id, finClr->bgn(id), finClr->end(id));
    w->addAndFilterOverlaps(gkp, finClr, errorRate, ovl, ovlLen);

    if (w->adjLen == 0) {
      //  All overlaps trimmed out!
      noCoverage += read->gkRead_sequenceLength();
      continue;
    }

    //  Find bad regions.

    //if (libr->gkLibrary_markBad() == true)
    //  //  From an external file, a list of known bad regions.  If no overlaps span
    //  //  the region with sufficient coverage, mark the region as bad.  This was
    //  //  motivated by the old 454 linker detection.
    //  markBad(gkp, w, subreadFile, doSubreadLoggingVerbose);

    //if (libr->gkLibrary_removeSpurReads() == true) {
    //  readsProcSpur += read->gkRead_sequenceLength();
    //  detectSpur(gkp, w, subreadFile, doSubreadLoggingVerbose);
    //  Get stats on spur region detected - save the length of each region to the trimStats object.
    //}

    //if (libr->gkLibrary_removeChimericReads() == true) {
    //  readsProcChimera += read->gkRead_sequenceLength();
    //  detectChimer(gkp, w, subreadFile, doSubreadLoggingVerbose);
    //  Get stats on chimera region detected - save the length of each region to the trimStats object.
    //}

    if (libr->gkLibrary_checkForSubReads() == true) {
      readsProcSubRead += read->gkRead_sequenceLength();
      detectSubReads(gkp, w, subreadFile, doSubreadLoggingVerbose);
    }

    //  Get stats on the bad regions found.  This kind of duplicates code in trimBadInterval(), but
    //  I don't want to pass all the stats objects into there.

    if (w->blist.size() == 0) {
      readsNoChange += read->gkRead_sequenceLength();
    }

    else {
      uint32  nSpur5   = 0, bSpur5   = 0;
      uint32  nSpur3   = 0, bSpur3   = 0;
      uint32  nChimera = 0, bChimera = 0;
      uint32  nSubread = 0, bSubread = 0;

      for (uint32 bb=0; bb<w->blist.size(); bb++) {
        switch (w->blist[bb].type) {
          case badType_5spur:
            nSpur5        += 1;
            basesBadSpur5 += w->blist[bb].end - w->blist[bb].bgn;
            break;
          case badType_3spur:
            nSpur3        += 1;
            basesBadSpur3 += w->blist[bb].end - w->blist[bb].bgn;
            break;
          case badType_chimera:
            nChimera        += 1;
            basesBadChimera += w->blist[bb].end - w->blist[bb].bgn;
            break;
          case badType_subread:
            nSubread        += 1;
            basesBadSubread += w->blist[bb].end - w->blist[bb].bgn;
            break;
          default:
            break;
        }
      }

      if (nSpur5   > 0)   readsBadSpur5   += nSpur5;
      if (nSpur3   > 0)   readsBadSpur3   += nSpur3;
      if (nChimera > 0)   readsBadChimera += nChimera;
      if (nSubread > 0)   readsBadSubread += nSubread;
    }

    //  Find solution.  This coalesces the list (in 'w') of all the bad regions found, picks out the
    //  largest good region, generates a log of the bad regions that support this decision, and sets
    //  the trim points.

    trimBadInterval(gkp, w, minReadLength, subreadFile, doSubreadLoggingVerbose);

    //  Log the solution.

    AS_UTL_safeWrite(reportFile, w->logMsg, "logMsg", sizeof(char), strlen(w->logMsg));

    //  Save the solution....

    outClr->setbgn(w->id) = w->clrBgn;
    outClr->setend(w->id) = w->clrEnd;

    //  And maybe delete the read.

    if (w->isOK == false) {
      deletedOut += read->gkRead_sequenceLength();

      outClr->setDeleted(w->id);
    }

    //  Update stats on what was trimmed.  The asserts say the clear range didn't expand, and the if
    //  tests if the clear range changed.

    assert(w->clrBgn >= w->iniBgn);
    assert(w->iniEnd >= w->clrEnd);

    if (w->clrBgn > w->iniBgn)
      readsTrimmed5 += w->clrBgn - w->iniBgn;

    if (w->iniEnd > w->clrEnd)
      readsTrimmed3 += w->iniEnd - w->clrEnd;
  }


  delete [] ovl;

  delete    w;

  gkp->gkStore_close();

  delete    finClr;
  delete    outClr;

  //  Close log files

  if (reportFile)
    fclose(reportFile);

  if (subreadFile)
    fclose(subreadFile);

  //  Write the summary

  if (outputPrefix) {
    sprintf(outputName, "%s.stats", outputPrefix);

    errno = 0;
    staFile = fopen(outputName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", outputName, strerror(errno));
  }

  if (staFile == NULL)
    staFile = stdout;

  //  Would like to know number of subreads per read

  fprintf(staFile, "PARAMETERS:\n");
  fprintf(staFile, "----------\n");
  fprintf(staFile, "%7u    (reads trimmed below this many bases are deleted)\n", minReadLength);
  fprintf(staFile, "%7.4f    (use overlaps at or below this fraction error)\n", errorRate);
  //fprintf(staFile, "%7u    (use only overlaps longer than this)\n", minAlignLength);  //  NOT SUPPORTED!
  fprintf(staFile, "INPUT READS:\n");
  fprintf(staFile, "-----------\n");
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (reads processed)\n", readsIn.nReads, readsIn.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (reads not processed, previously deleted)\n", deletedIn.nReads, deletedIn.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (reads not processed, in a library where trimming isn't allowed)\n", noTrimIn.nReads, noTrimIn.nBases);
  fprintf(staFile, "\n");
  fprintf(staFile, "PROCESSED:\n");
  fprintf(staFile, "--------\n");
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (no overlaps)\n", noOverlaps.nReads, noOverlaps.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (no coverage after adjusting for trimming done already)\n", noCoverage.nReads, noCoverage.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (processed for chimera)\n",  readsProcChimera.nReads, readsProcChimera.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (processed for spur)\n",     readsProcSpur.nReads,    readsProcSpur.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (processed for subreads)\n", readsProcSubRead.nReads, readsProcSubRead.nBases);
  fprintf(staFile, "\n");
  fprintf(staFile, "READS WITH SIGNALS:\n");
  fprintf(staFile, "------------------\n");
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" signals (number of 5' spur signal)\n", readsBadSpur5.nReads,   readsBadSpur5.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" signals (number of 3' spur signal)\n", readsBadSpur3.nReads,   readsBadSpur3.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" signals (number of chimera signal)\n", readsBadChimera.nReads, readsBadChimera.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" signals (number of subread signal)\n", readsBadSubread.nReads, readsBadSubread.nBases);
  fprintf(staFile, "\n");
  fprintf(staFile, "SIGNALS:\n");
  fprintf(staFile, "-------\n");
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (size of 5' spur signal)\n", basesBadSpur5.nReads,   basesBadSpur5.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (size of 3' spur signal)\n", basesBadSpur3.nReads,   basesBadSpur3.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (size of chimera signal)\n", basesBadChimera.nReads, basesBadChimera.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (size of subread signal)\n", basesBadSubread.nReads, basesBadSubread.nBases);
  fprintf(staFile, "\n");
  fprintf(staFile, "TRIMMING:\n");
  fprintf(staFile, "--------\n");
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (trimmed from the 5' end of the read)\n", readsTrimmed5.nReads, readsTrimmed5.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (trimmed from the 3' end of the read)\n", readsTrimmed3.nReads, readsTrimmed3.nBases);

#if 0
  fprintf(staFile, "DELETED:\n");
  fprintf(staFile, "-------\n");
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (deleted because of both cimera and spur signals)\n", bothDeletedSmall.nReads, bothDeletedSmall.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (deleted because of chimera signal)\n", chimeraDeletedSmall.nReads, chimeraDeletedSmall.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (deleted because of spur signal)\n", spurDeletedSmall.nReads, spurDeletedSmall.nBases);
  fprintf(staFile, "\n");
  fprintf(staFile, "SPUR TYPES:\n");
  fprintf(staFile, "----------\n");
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (normal spur detected)\n", spurDetectedNormal.nReads, spurDetectedNormal.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (linker spur detected)\n", spurDetectedLinker.nReads, spurDetectedLinker.nBases);
  fprintf(staFile, "\n");
  fprintf(staFile, "CHIMERA TYPES:\n");
  fprintf(staFile, "-------------\n");
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (innie-pair chimera detected)\n", chimeraDetectedInnie.nReads, chimeraDetectedInnie.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (overhanging chimera detected)\n", chimeraDetectedOverhang.nReads, chimeraDetectedOverhang.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (gap chimera detected)\n", chimeraDetectedGap.nReads, chimeraDetectedGap.nBases);
  fprintf(staFile, "%6"F_U32P" reads %12"F_U64P" bases (linker chimera detected)\n", chimeraDetectedLinker.nReads, chimeraDetectedLinker.nBases);
#endif

  //  INPUT READS  = ACCEPTED + TRIMMED + DELETED
  //  SPUR TYPE    = TRIMMED and DELETED spur and both categories
  //  CHIMERA TYPE = TRIMMED and DELETED chimera and both categories

  if (staFile != stdout)
    fclose(staFile);

  exit(0);
}

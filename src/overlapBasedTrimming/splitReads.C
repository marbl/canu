
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "splitReads.H"



int
main(int argc, char **argv) {

  argc = AS_configure(argc, argv);

  char              *gkpName = NULL;
  char              *ovsName = NULL;

  char              *finClrName = NULL;
  char              *outClrName = NULL;

  double             errorRate       = 0.06;
  uint32             minAlignLength  = 40;
  uint32             minReadLength   = 64;

  uint32             idMin = 1;
  uint32             idMax = UINT32_MAX;

  char              *outputPrefix = NULL;
  char               outputName[FILENAME_MAX];

  FILE              *summaryFile  = NULL;
  FILE              *reportFile   = NULL;
  FILE              *subreadFile  = NULL;

  bool               doSubreadLogging        = true;
  bool               doSubreadLoggingVerbose = false;

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

    } else if (strcmp(argv[arg], "-l") == 0) {
      minAlignLength = atoi(argv[++arg]);

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

    if (finClr->isDeleted(id))
      //  Read already trashed.
      continue;

    if ((libr->gkLibrary_removeSpurReads()     == false) &&
        (libr->gkLibrary_removeChimericReads() == false) &&
        (libr->gkLibrary_checkForSubReads()    == false))
      //  Nothing to do.
      continue;

    uint32   nLoaded = ovs->readOverlaps(id, ovl, ovlLen, ovlMax);

    //fprintf(stderr, "read %7u with %7u overlaps\r", id, nLoaded);

    if (nLoaded == 0)
      //  No overlaps, nothing to check!
      continue;

    w->clear(id, finClr->bgn(id), finClr->end(id));
    w->addAndFilterOverlaps(gkp, finClr, errorRate, ovl, ovlLen);

    if (w->adjLen == 0)
      //  All overlaps trimmed out!
      continue;

    //  Find bad regions.

    //if (libr->gkLibrary_markBad() == true)
    //  //  From an external file, a list of known bad regions.  If no overlaps span
    //  //  the region with sufficient coverage, mark the region as bad.  This was
    //  //  motivated by the old 454 linker detection.
    //  markBad(gkp, w, subreadFile, doSubreadLoggingVerbose);

    //if (libr->gkLibrary_removeSpurReads() == true)
    //  detectSpur(gkp, w, subreadFile, doSubreadLoggingVerbose);

    //if (libr->gkLibrary_removeChimericReads() == true)
    //  detectChimer(gkp, w, subreadFile, doSubreadLoggingVerbose);

    if (libr->gkLibrary_checkForSubReads() == true)
      detectSubReads(gkp, w, subreadFile, doSubreadLoggingVerbose);

    //  Find solution.

    trimBadInterval(gkp, w, minReadLength, subreadFile, doSubreadLoggingVerbose);

    //  Report solution.

    AS_UTL_safeWrite(reportFile, w->logMsg, "logMsg", sizeof(char), strlen(w->logMsg));

    outClr->setbgn(w->id) = w->clrBgn;
    outClr->setend(w->id) = w->clrEnd;

    if (w->isOK == false) {
      outClr->setDeleted(w->id);
      fprintf(stderr, "\n");
    }
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

  sprintf(outputName, "%s.summary",     outputPrefix);
  errno = 0;
  summaryFile = fopen(outputName, "w");
  if (errno) {
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", outputName, strerror(errno));
    summaryFile = stdout;
  }

#if 0
  fprintf(summaryFile, "READS (= ACEEPTED + TRIMMED + DELETED)\n");
  fprintf(summaryFile, "  total processed       "F_U32"\n", readsProcessed);
  fprintf(summaryFile, "\n");
  fprintf(summaryFile, "ACCEPTED\n");
  fprintf(summaryFile, "  no coverage           "F_U32"\n", noCoverage);
  fprintf(summaryFile, "  full coverage         "F_U32"\n", fullCoverage);
  fprintf(summaryFile, "  no signal, no gaps    "F_U32"\n", noSignalNoGap);
  fprintf(summaryFile, "  no signal, gaps       "F_U32"\n", noSignalButGap);
  fprintf(summaryFile, "\n");
  fprintf(summaryFile, "TRIMMED\n");
  fprintf(summaryFile, "  both                  "F_U32"\n", bothFixed);
  fprintf(summaryFile, "  chimera               "F_U32"\n", chimeraFixed);
  fprintf(summaryFile, "  spur                  "F_U32"\n", spurFixed);
  fprintf(summaryFile, "\n");
  fprintf(summaryFile, "DELETED\n");
  fprintf(summaryFile, "  both                  "F_U32"\n", bothDeletedSmall);
  fprintf(summaryFile, "  chimera               "F_U32"\n", chimeraDeletedSmall);
  fprintf(summaryFile, "  spur                  "F_U32"\n", spurDeletedSmall);
  fprintf(summaryFile, "\n");
  fprintf(summaryFile, "SPUR TYPES (= TRIMMED/DELETED spur + both)\n");
  fprintf(summaryFile, "  normal                "F_U32"\n", spurDetectedNormal);
  fprintf(summaryFile, "  linker                "F_U32"\n", spurDetectedLinker);
  fprintf(summaryFile, "\n");
  fprintf(summaryFile, "CHIMERA TYPES (= TRIMMED/DELETED chimera + both)\n");
  fprintf(summaryFile, "  innie pair            "F_U32"\n", chimeraDetectedInnie);
  fprintf(summaryFile, "  overhang              "F_U32"\n", chimeraDetectedOverhang);
  fprintf(summaryFile, "  gap                   "F_U32"\n", chimeraDetectedGap);
  fprintf(summaryFile, "  linker                "F_U32"\n", chimeraDetectedLinker);
#endif

  if (summaryFile != stdout)
    fclose(summaryFile);

  exit(0);
}

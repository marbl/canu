
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
 *    Brian P. Walenz from 2015-APR-09 to 2015-SEP-21
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-27
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"
#include "tgStore.H"

#include "outputFalcon.H"

#include "stashContains.H"

#include "splitToWords.H"

#include <set>

using namespace std;

//  Debugging on which reads are filtered, which are used, and which are removed
//  to meet coverage thresholds.  Very big.
#undef DEBUG_LAYOUT


//  Generate a layout for the read in ovl[0].a_iid, using most or all of the overlaps
//  in ovl.

tgTig *
generateLayout(gkStore    *gkpStore,
               uint64     *readScores,
	       bool	   legacyScore,
               uint32      minEvidenceLength,
               double      maxEvidenceErate,
               double      maxEvidenceCoverage,
               ovOverlap *ovl,
               uint32      ovlLen,
               FILE       *flgFile) {

  set<uint32_t> children;
  tgTig  *layout = new tgTig;

  layout->_tigID           = ovl[0].a_iid;
  layout->_coverageStat    = 1.0;  //  Default to just barely unique
  layout->_microhetProb    = 1.0;  //  Default to 100% probability of unique

  layout->_class           = tgTig_noclass;
  layout->_suggestRepeat   = false;
  layout->_suggestCircular = false;

  gkRead  *read = gkpStore->gkStore_getRead(ovl[0].a_iid);

  layout->_layoutLen       = read->gkRead_sequenceLength();

  resizeArray(layout->_children, layout->_childrenLen, layout->_childrenMax, ovlLen, resizeArray_doNothing);

  if (flgFile)
    fprintf(flgFile, "Generate layout for read "F_U32" length "F_U32" using up to "F_U32" overlaps.\n",
            layout->_tigID, layout->_layoutLen, ovlLen);

  for (uint32 oo=0; oo<ovlLen; oo++) {

    //  ovlLength, in filterCorrectionOverlaps, is computed on the a read.  That is now the b read here.
    uint64   ovlLength = ((ovl[oo].b_bgn() < ovl[oo].b_end()) ?
                          ovl[oo].b_end() - ovl[oo].b_bgn() :
                          ovl[oo].b_bgn() - ovl[oo].b_end());
    uint64   ovlScore  = 100 * ovlLength * (1 - ovl[oo].erate());
    if (legacyScore) {
       ovlScore  = ovlLength << AS_MAX_EVALUE_BITS;
       ovlScore |= (AS_MAX_EVALUE - ovl[oo].evalue());
    }

    if (ovlLength > AS_MAX_READLEN) {
      char ovlString[1024];
      fprintf(stderr, "ERROR: bogus overlap '%s'\n", ovl[oo].toString(ovlString, ovOverlapAsCoords, false));
    }
    assert(ovlLength < AS_MAX_READLEN);

    if (ovl[oo].erate() > maxEvidenceErate) {
      if (flgFile)
        fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5u erate %.3f - low quality (threshold %.2f)\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), maxEvidenceErate);
      continue;
    }

    if (ovl[oo].a_end() - ovl[oo].a_bgn() < minEvidenceLength) {
      if (flgFile)
        fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5u erate %.3f - too short (threshold %u)\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), minEvidenceLength);
      continue;
    }

    if ((readScores != NULL) &&
        (ovlScore < readScores[ovl[oo].b_iid])) {
      if (flgFile)
        fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5u erate %.3f - filtered by global filter (threshold "F_U64")\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), readScores[ovl[oo].b_iid]);
      continue;
    }

    if (children.find(ovl[oo].b_iid) != children.end()) {
      if (flgFile)
        fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5u erate %.3f - duplicate\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate());
      continue;
    }

    if (flgFile)
      fprintf(flgFile, "  allow  read %9u at position %6u,%6u length %5u erate %.3f\n",
              ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate());

    tgPosition   *pos = layout->addChild();

    //  Set the read.  Parent is always the read we're building for, hangs and position come from
    //  the overlap.  Easy as pie!

    if (ovl[oo].flipped() == false) {
      pos->set(ovl[oo].b_iid,
               ovl[oo].a_iid,
               ovl[oo].a_hang(),
               ovl[oo].b_hang(),
               ovl[oo].a_bgn(), ovl[oo].a_end());

    } else {
      pos->set(ovl[oo].b_iid,
               ovl[oo].a_iid,
               ovl[oo].a_hang(),
               ovl[oo].b_hang(),
               ovl[oo].a_end(), ovl[oo].a_bgn());
    }

    //  Remember the unaligned bit!

    pos->_askip = ovl[oo].dat.ovl.bhg5;
    pos->_bskip = ovl[oo].dat.ovl.bhg3;

    // record the ID
    children.insert(ovl[oo].b_iid);
  }

  //  Use utgcns's stashContains to get rid of extra coverage; we don't care about it, and
  //  just delete it immediately.

  savedChildren *sc = stashContains(layout, maxEvidenceCoverage);

  if ((flgFile) && (sc))
    sc->reportRemoved(flgFile, layout->tigID());

  if (sc) {
    delete sc->children;
    delete sc;
  }

  //  stashContains also sorts by position, so we're done.

#if 0
  if (flgFile)
    for (uint32 ii=0; ii<layout->numberOfChildren(); ii++)
      fprintf(flgFile, "  read %9u at position %6u,%6u hangs %6d %6d %c unAl %5d %5d\n",
              layout->getChild(ii)->_objID,
              layout->getChild(ii)->_min,
              layout->getChild(ii)->_max,
              layout->getChild(ii)->_ahang,
              layout->getChild(ii)->_bhang,
              layout->getChild(ii)->isForward() ? 'F' : 'R',
              layout->getChild(ii)->_askip,
              layout->getChild(ii)->_bskip);
#endif

  return(layout);
}




int
main(int argc, char **argv) {
  char             *gkpName   = 0L;
  char             *ovlName   = 0L;
  char             *scoreName = 0L;
  char             *tigName   = 0L;

  bool              falconOutput = false;  //  To stdout
  bool              trimToAlign  = false;

  uint32            errorRate = AS_OVS_encodeEvalue(0.015);

  char             *outputPrefix = NULL;
  char              logName[FILENAME_MAX] = {0};
  char              sumName[FILENAME_MAX] = {0};
  char              flgName[FILENAME_MAX] = {0};
  FILE             *logFile = 0L;
  FILE             *sumFile = 0L;
  FILE             *flgFile = 0L;

  uint32            minEvidenceOverlap  = 40;
  uint32            minEvidenceCoverage = 1;

  uint32            iidMin       = 0;
  uint32            iidMax       = UINT32_MAX;
  char             *readListName = NULL;
  set<uint32>       readList;

  uint32            minEvidenceLength   = 0;
  double            maxEvidenceErate    = 1.0;
  double            maxEvidenceCoverage = DBL_MAX;

  uint32            minCorLength        = 0;

  bool              filterCorLength     = false;
  bool		    legacyScore	        = false;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {  //  Input gkpStore
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {  //  Input ovlStore
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-S") == 0) {  //  Input scores
      scoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {  //  Output tigStore
      tigName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {  //  Output directly to falcon, not tigStore
      falconOutput = true;
      trimToAlign  = true;

    } else if (strcmp(argv[arg], "-p") == 0) {  //  Output prefix, just logging and summary
      outputPrefix = argv[++arg];


    } else if (strcmp(argv[arg], "-b") == 0) {  //  Begin read range
      iidMin  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {  //  End read range
      iidMax  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-rl") == 0) {  //  List of reads to correct, will also apply -b/-e range
      readListName = argv[++arg];


    } else if (strcmp(argv[arg], "-L") == 0) {  //  Minimum length of evidence overlap
      minEvidenceLength  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {  //  Max error rate of evidence overlap
      maxEvidenceErate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-C") == 0) {  //  Max coverage of evidence reads to emit.
      maxEvidenceCoverage = atof(argv[++arg]);


    } else if (strcmp(argv[arg], "-M") == 0) {  //  Minimum length of a corrected read
      minCorLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-legacy") == 0) {
      legacyScore = true;

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (gkpName == NULL)
    err++;
  if (ovlName == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore [ -T tigStore | -F ] ...\n", argv[0]);
    fprintf(stderr, "  -G gkpStore   mandatory path to gkpStore\n");
    fprintf(stderr, "  -O ovlStore   mandatory path to ovlStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S file       global score (binary) input file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -T corStore   output layouts to tigStore corStore\n");
    fprintf(stderr, "  -F            output falconsense-style input directly to stdout\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -p  name      output prefix name, for logging and summary\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -b  bgnID     \n");
    fprintf(stderr, "  -e  endID     \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -rl file      \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L  length    minimum length of evidence overlaps\n");
    fprintf(stderr, "  -E  erate     maxerror rate of evidence overlaps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -C  coverage  maximum coverage of evidence reads to emit\n");
    fprintf(stderr, "  -M  length    minimum length of a corrected read\n");
    fprintf(stderr, "\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR: no gkpStore input (-G) supplied.\n");
    if (ovlName == NULL)
      fprintf(stderr, "ERROR: no ovlStore input (-O) supplied.\n");

    exit(1);
  }

  //  Open inputs and output tigStore.

  gkStore  *gkpStore = gkStore::gkStore_open(gkpName);
  ovStore  *ovlStore = new ovStore(ovlName, gkpStore);
  tgStore  *tigStore = (tigName != NULL) ? new tgStore(tigName) : NULL;

  //  Load read scores, if supplied.

  uint64   *readScores = NULL;

  if (scoreName) {
    readScores = new uint64 [gkpStore->gkStore_getNumReads() + 1];

    errno = 0;
    FILE *scoreFile = fopen(scoreName, "r");
    if (errno)
      fprintf(stderr, "failed to open '%s' for reading: %s\n", scoreName, strerror(errno)), exit(1);

    AS_UTL_safeRead(scoreFile, readScores, "scores", sizeof(uint64), gkpStore->gkStore_getNumReads() + 1);

    fclose(scoreFile);
  }

  //  Threshold the range of reads to operate on.

  if (gkpStore->gkStore_getNumReads() < iidMin) {
    fprintf(stderr, "ERROR: only "F_U32" reads in the store (IDs 0-"F_U32" inclusive); can't process requested range -b "F_U32" -e "F_U32"\n",
            gkpStore->gkStore_getNumReads(),
            gkpStore->gkStore_getNumReads()-1,
            iidMin, iidMax);
    exit(1);
  }

  if (gkpStore->gkStore_getNumReads() < iidMax)
      iidMax = gkpStore->gkStore_getNumReads();

  ovlStore->setRange(iidMin, iidMax);

  //  If a readList is supplied, load it, respecting the iidMin/iidMax (only to cut down on the
  //  size).

  if (readListName != NULL) {
    errno = 0;

    char  L[1024];
    FILE *R = fopen(readListName, "r");

    if (errno)
      fprintf(stderr, "Failed to open '%s' for reading: %s\n", readListName, strerror(errno)), exit(1);

    fgets(L, 1024, R);
    while (!feof(R)) {
      splitToWords W(L);
      uint32       id = W(0);

      if ((iidMin <= id) &&
          (id <= iidMax))
        readList.insert(W(0));

      fgets(L, 1024, R);
    }

    fclose(R);
  }

  //  Open logging and summary files

  if (outputPrefix) {
    sprintf(logName, "%s.log",        outputPrefix);
    sprintf(sumName, "%s.summary",    outputPrefix);
    sprintf(flgName, "%s.filter.log", outputPrefix);

    errno = 0;

    logFile = fopen(logName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", logName, strerror(errno)), exit(1);

    sumFile = fopen(sumName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", sumName, strerror(errno)), exit(1);

#ifdef DEBUG_LAYOUT
    flgFile = fopen(flgName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", flgName, strerror(errno)), exit(1);
#endif
  }

  if (logFile)
    fprintf(logFile, "read\torigLen\tnumOlaps\tcorLen\n");

  //  Initialize processing.

  uint32       ovlMax = 1024 * 1024;
  uint32       ovlLen = 0;
  ovOverlap   *ovl    = ovOverlap::allocateOverlaps(gkpStore, ovlMax);

  ovlLen = ovlStore->readOverlaps(ovl, ovlMax, true);

  gkReadData   *readData = new gkReadData;

  //  And process.

  while (ovlLen > 0) {
    bool   skipIt        = false;
    char   skipMsg[1024] = {0};

    tgTig *layout = generateLayout(gkpStore,
                                   readScores,
                                   legacyScore,
                                   minEvidenceLength, maxEvidenceErate, maxEvidenceCoverage,
                                   ovl, ovlLen,
                                   flgFile);

    //  If there was a readList, skip anything not in it.

    if ((readListName != NULL) &&
        (readList.count(layout->tigID()) == 0)) {
      strcat(skipMsg, "\tnot_in_readList");
      skipIt = true;
    }

    //  Possibly filter by the length of the uncorrected read.

    gkRead *read = gkpStore->gkStore_getRead(layout->tigID());

    if (read->gkRead_sequenceLength() < minCorLength) {
      strcat(skipMsg, "\tread_too_short");
      skipIt = true;
    }

    //  Possibly filter by the length of the corrected read.

    uint32  minPos = UINT32_MAX;
    uint32  maxPos = 0;
    uint32  corLen = 0;

    for (uint32 ii=0; ii<layout->numberOfChildren(); ii++) {
      tgPosition *pos = layout->getChild(ii);

      if (pos->_min < minPos)
        minPos = pos->_min;

      if (maxPos < pos->_max)
        maxPos = pos->_max;
    }

    if (minPos != UINT32_MAX)
      corLen = maxPos - minPos;

    if (corLen < minCorLength) {
      strcat(skipMsg, "\tcorrection_too_short");
      skipIt = true;
    }

    //  Filter out empty tigs - these either have no overlaps, or failed the
    //  length check in generateLayout.

    if (layout->numberOfChildren() <= 1) {
      strcat(skipMsg, "\tno_children");
      skipIt = true;
    }

    //  Output, if not skipped.

    if (logFile)
      fprintf(logFile, "%u\t%u\t%u\t%u%s\n",
              layout->tigID(), read->gkRead_sequenceLength(), layout->numberOfChildren(), corLen, skipMsg);

    if ((skipIt == false) && (tigStore != NULL))
      tigStore->insertTig(layout, false);

    if ((skipIt == false) && (falconOutput == true))
      outputFalcon(gkpStore, layout, trimToAlign, stdout, readData);

    delete layout;

    //  Load next batch of overlaps.

    ovlLen = ovlStore->readOverlaps(ovl, ovlMax, true);
  }

  if (falconOutput)
    fprintf(stdout, "- -\n");

  delete readData;

  if (logFile != NULL)
    fclose(logFile);

  if (sumFile != NULL)
    fclose(sumFile);

  if (flgFile != NULL)
    fclose(flgFile);

  delete tigStore;
  delete ovlStore;

  gkpStore->gkStore_close();

  return(0);
}


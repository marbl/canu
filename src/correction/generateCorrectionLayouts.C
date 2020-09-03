
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

#include "runtime.H"
#include "sqStore.H"
#include "ovStore.H"
#include "tgStore.H"

#include "strings.H"
#include "files.H"
#include "intervalList.H"
#include "sequence.H"

#include <set>

using namespace std;

//  Debugging on which reads are filtered, which are used, and which are removed
//  to meet coverage thresholds.  Very big.
#undef DEBUG_LAYOUT

//FILE *flgFile = stderr;
FILE *flgFile = NULL;



uint16 *
loadThresholds(sqStore *seqStore,
               ovStore *ovlStore,
               char    *scoreName,
               uint32   expectedCoverage,
               FILE    *scoFile) {
  uint32   numReads   = seqStore->sqStore_lastReadID();
  uint16  *olapThresh = new uint16 [numReads + 1];

  if (scoreName != NULL)
    AS_UTL_loadFile(scoreName, olapThresh, numReads + 1);

  else {
    ovStoreHistogram  *ovlHisto = ovlStore->getHistogram();

    for (uint32 ii=0; ii<numReads+1; ii++)
      olapThresh[ii] = ovlHisto->overlapScoreEstimate(ii, expectedCoverage, scoFile);

    delete ovlHisto;
  }

  return(olapThresh);
}



static
void
generateLayout(tgTig      *layout,
               uint16     *olapThresh,
               uint32      minEvidenceLength,
               double      maxEvidenceErate,
               double      maxEvidenceCoverage,
               ovOverlap *ovl,
               uint32      ovlLen,
               FILE       *logFile) {

  //  Generate a layout for the read in ovl[0].a_iid, using most or all of the overlaps in ovl.

  layout->allocateChildren(ovlLen);

  if (logFile)
    fprintf(logFile, "Generate layout for read " F_U32 " length " F_U32 " using up to " F_U32 " overlaps.\n",
            layout->_tigID, layout->_layoutLen, ovlLen);

  set<uint32_t>  children;

  for (uint32 oo=0; oo<ovlLen; oo++) {
    uint64   ovlLength = ovl[oo].b_len();
    uint16   ovlScore  = ovl[oo].overlapScore(true);

    if (ovlLength > AS_MAX_READLEN) {
      char ovlString[1024];
      fprintf(stderr, "ERROR: bogus overlap '%s'\n", ovl[oo].toString(ovlString, ovOverlapAsCoords, false));
    }
    assert(ovlLength < AS_MAX_READLEN);

    if (ovl[oo].erate() > maxEvidenceErate) {
      if (logFile)
        fprintf(logFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - low quality (threshold %.2f)\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), maxEvidenceErate);
      continue;
    }

    if (ovl[oo].a_end() - ovl[oo].a_bgn() < minEvidenceLength) {
      if (logFile)
        fprintf(logFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - too short (threshold %u)\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), minEvidenceLength);
      continue;
    }

    if ((olapThresh != NULL) &&
        (ovlScore < olapThresh[ovl[oo].b_iid])) {
      if (logFile)
        fprintf(logFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - filtered by global filter (threshold " F_U16 ")\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), olapThresh[ovl[oo].b_iid]);
      continue;
    }

    if (children.find(ovl[oo].b_iid) != children.end()) {
      if (logFile)
        fprintf(logFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - duplicate\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate());
      continue;
    }

    if (logFile)
      fprintf(logFile, "  allow  read %9u at position %6u,%6u length %5lu erate %.3f\n",
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

    //  Remember we added this read - to filter read with both fwd/rev overlaps.

    children.insert(ovl[oo].b_iid);
  }

  //  Dro excess coverage in the evidence by dropping short reads.  This also
  //  sorts by position.

  layout->dropExcessCoverage(maxEvidenceCoverage);
}





int
main(int argc, char **argv) {
  char             *seqName    = 0L;
  char             *ovlName    = 0L;
  char             *corName    = 0L;
  char             *flgName    = 0L;

  char             *scoreName  = 0L;

  uint32            errorRate  = AS_OVS_encodeEvalue(0.015);

  bool              dumpScores = false;
  bool              doLogging  = false;

  uint32            expectedCoverage    = 40;    //  How many overlaps per read to save, global filter
  uint32            minEvidenceOverlap  = 40;

  uint32            iidMin = 1;
  uint32            iidMax = UINT32_MAX;

  uint32            minEvidenceLength   = 0;
  double            maxEvidenceErate    = 1.0;
  double            maxEvidenceCoverage = DBL_MAX;


  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {   //  INPUTS
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-scores") == 0) {
      scoreName = argv[++arg];


    } else if (strcmp(argv[arg], "-C") == 0) {   //  OUTPUT FORMAT
      corName = argv[++arg];


    } else if (strcmp(argv[arg], "-b") == 0) {   //  READ SELECTION
      iidMin  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      iidMax  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-eL") == 0) {   //  EVIDENCE SELECTION
      minEvidenceLength  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-eE") == 0) {
      maxEvidenceErate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-eC") == 0) {
      maxEvidenceCoverage = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-V") == 0) {
      doLogging = true;

    } else if (strcmp(argv[arg], "-D") == 0) {
      dumpScores = true;


    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (seqName == NULL)
    err++;
  if (corName == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -S seqStore -O ovlStore ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS\n");
    fprintf(stderr, "  -S seqStore      mandatory path to seqStore\n");
    fprintf(stderr, "  -O ovlStore      mandatory path to ovlStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -scores sf       overlap score thresholds (from filterCorrectionOverlaps)\n");
    fprintf(stderr, "                   if not supplied, will be estimated from ovlStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUTS\n");
    fprintf(stderr, "  -C corStore      output layouts to store 'corStore'\n");
    fprintf(stderr, "  -V               write extremely verbose logging to 'corStore.log'\n");
    fprintf(stderr, "  -D               dump the data used to estimate overlap scores to 'corStore.scores'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "READ SELECTION\n");
    fprintf(stderr, "  -b bgnID         process reads starting at bgnID\n");
    fprintf(stderr, "  -e endID         process reads up to but not including endID\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "EVIDENCE SELECTION\n");
    fprintf(stderr, "  -eL length       minimum length of evidence overlaps\n");
    fprintf(stderr, "  -eE erate        maximum error rate of evidence overlaps\n");
    fprintf(stderr, "  -eC coverage     maximum coverage of evidence reads to emit\n");
    fprintf(stderr, "\n");

    if (seqName == NULL)
      fprintf(stderr, "ERROR: no input seqStore (-S) supplied.\n");
    if (corName == NULL)
      fprintf(stderr, "ERROR: no output corStore (-C) supplied.\n");
    exit(1);
  }

  //  Open inputs and output tigStore.

  sqRead_setDefaultVersion(sqRead_raw);

  sqStore  *seqStore = new sqStore(seqName);
  ovStore  *ovlStore = new ovStore(ovlName, seqStore);
  tgStore  *corStore = new tgStore(corName);

  uint32    numReads = seqStore->sqStore_lastReadID();

  //  Threshold the range of reads to operate on.

  if (numReads < iidMin) {
    fprintf(stderr, "ERROR: only " F_U32 " reads in the store (IDs 0-" F_U32 " inclusive); can't process requested range -b " F_U32 " -e " F_U32 "\n",
            numReads,
            numReads-1,
            iidMin, iidMax);
    exit(1);
  }
  if (numReads < iidMax)
    iidMax = numReads;

  ovlStore->setRange(iidMin, iidMax);

  //  Open logging and summary files

  FILE *logFile = AS_UTL_openOutputFile(corName, '.', "log",    doLogging);
  FILE *scoFile = AS_UTL_openOutputFile(corName, '.', "scores", dumpScores);

  //  Load read scores, if supplied.

  uint16   *olapThresh = loadThresholds(seqStore, ovlStore, scoreName, expectedCoverage, scoFile);

  //  Initialize processing.

  uint32             ovlMax    = 0;
  ovOverlap         *ovl       = NULL;

  //  And process.

  for (uint32 rr=1; rr<numReads+1; rr++) {
    uint32 ovlLen = ovlStore->loadOverlapsForRead(rr, ovl, ovlMax);

    if (ovlLen > 0) {
      tgTig   *layout = new tgTig;

      layout->_tigID     = rr;
      layout->_layoutLen = seqStore->sqStore_getReadLength(rr, sqRead_raw);

      generateLayout(layout,
                     olapThresh,
                     minEvidenceLength, maxEvidenceErate, maxEvidenceCoverage,
                     ovl, ovlLen,
                     logFile);

      corStore->insertTig(layout, false);

      delete layout;
    }
  }

  //  Close files and clean up.

  AS_UTL_closeFile(logFile);

  delete [] olapThresh;
  delete [] ovl;
  delete    corStore;
  delete    ovlStore;

  delete seqStore;

  fprintf(stderr, "Bye.\n");

  return(0);
}

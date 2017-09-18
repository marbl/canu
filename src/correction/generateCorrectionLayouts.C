
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

#include "stashContains.H"

#include "splitToWords.H"
#include "intervalList.H"
#include "AS_UTL_reverseComplement.H"
#include "AS_UTL_fasta.H"

#include <set>

using namespace std;

//  Debugging on which reads are filtered, which are used, and which are removed
//  to meet coverage thresholds.  Very big.
#undef DEBUG_LAYOUT



uint16 *
loadThresholds(gkStore *gkpStore,
               ovStore *ovlStore,
               char    *scoreName,
               uint32   expectedCoverage) {
  uint32   numReads   = gkpStore->gkStore_getNumReads();
  uint16  *olapThresh = new uint16 [numReads + 1];

  if (scoreName != NULL) {
    errno = 0;
    FILE *S = fopen(scoreName, "r");
    if (errno)
      fprintf(stderr, "failed to open '%s' for reading: %s\n", scoreName, strerror(errno)), exit(1);

    AS_UTL_safeRead(S, olapThresh, "scores", sizeof(uint16), numReads + 1);

    fclose(S);
  }

  else {
    ovStoreHistogram  *ovlHisto = ovlStore->getHistogram();

    for (uint32 ii=0; ii<numReads+1; ii++)
      olapThresh[ii] = ovlHisto->overlapScoreEstimate(ii, expectedCoverage);

    delete ovlHisto;
  }

  return(olapThresh);
}



tgTig *
generateLayout(tgTig      *layout,
               uint16     *olapThresh,
               uint32      minEvidenceLength,
               double      maxEvidenceErate,
               double      maxEvidenceCoverage,
               ovOverlap *ovl,
               uint32      ovlLen) {

  //  Generate a layout for the read in ovl[0].a_iid, using most or all of the overlaps in ovl.

  resizeArray(layout->_children, layout->_childrenLen, layout->_childrenMax, ovlLen, resizeArray_doNothing);

  //if (flgFile)
  //  fprintf(flgFile, "Generate layout for read " F_U32 " length " F_U32 " using up to " F_U32 " overlaps.\n",
  //          layout->_tigID, layout->_layoutLen, ovlLen);

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
      //if (flgFile)
      //  fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - low quality (threshold %.2f)\n",
      //          ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), maxEvidenceErate);
      continue;
    }

    if (ovl[oo].a_end() - ovl[oo].a_bgn() < minEvidenceLength) {
      //if (flgFile)
      //  fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - too short (threshold %u)\n",
      //          ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), minEvidenceLength);
      continue;
    }

    if ((olapThresh != NULL) &&
        (ovlScore < olapThresh[ovl[oo].b_iid])) {
      //if (flgFile)
      //  fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - filtered by global filter (threshold " F_U16 ")\n",
      //          ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), olapThresh[ovl[oo].b_iid]);
      continue;
    }

    if (children.find(ovl[oo].b_iid) != children.end()) {
      //if (flgFile)
      //  fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - duplicate\n",
      //          ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate());
      continue;
    }

    //if (flgFile)
    //  fprintf(flgFile, "  allow  read %9u at position %6u,%6u length %5lu erate %.3f\n",
    //          ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate());

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

  //  Use utgcns's stashContains to get rid of extra coverage; we don't care about it, and
  //  just delete it immediately.

  savedChildren *sc = stashContains(layout, maxEvidenceCoverage);

  //if ((flgFile) && (sc))
  //  sc->reportRemoved(flgFile, layout->tigID());

  if (sc) {
    delete sc->children;
    delete sc;
  }

  //  stashContains also sorts by position, so we're done.

  return(layout);
}





int
main(int argc, char **argv) {
  char             *gkpName   = 0L;
  char             *ovlName   = 0L;
  char             *corName   = 0L;

  char             *scoreName = 0L;

  uint32            errorRate = AS_OVS_encodeEvalue(0.015);

  char             *outputPrefix = NULL;
  char              logName[FILENAME_MAX] = {0};
  char              sumName[FILENAME_MAX] = {0};
  char              flgName[FILENAME_MAX] = {0};
  FILE             *logFile = 0L;
  FILE             *sumFile = 0L;

  uint32            expectedCoverage    = 40;    //  How many overlaps per read to save, global filter
  uint32            minEvidenceOverlap  = 40;
  uint32            minEvidenceCoverage = 4;

  uint32            iidMin       = 0;
  uint32            iidMax       = UINT32_MAX;

  uint32            minEvidenceLength   = 0;
  double            maxEvidenceErate    = 1.0;
  double            maxEvidenceCoverage = DBL_MAX;

  uint32            minCorLength        = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {   //  INPUTS
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-S") == 0) {
      scoreName = argv[++arg];


    } else if (strcmp(argv[arg], "-C") == 0) {   //  OUTPUT FORMAT
      corName = argv[++arg];

    } else if (strcmp(argv[arg], "-p") == 0) {
      outputPrefix = argv[++arg];


    } else if (strcmp(argv[arg], "-b") == 0) {   //  READ SELECTION
      iidMin  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      iidMax  = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-eL") == 0) {   //  EVIDENCE SELECTION
      minEvidenceLength  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-eE") == 0) {
      maxEvidenceErate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-ec") == 0) {
      minEvidenceCoverage = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-eC") == 0) {
      maxEvidenceCoverage = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-eM") == 0) {
      minCorLength = atoi(argv[++arg]);


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
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS\n");
    fprintf(stderr, "  -G gkpStore      mandatory path to gkpStore\n");
    fprintf(stderr, "  -O ovlStore      mandatory path to ovlStore\n");
    fprintf(stderr, "  -S file          overlap score thresholds (from filterCorrectionOverlaps)\n");
    fprintf(stderr, "                     if not supplied, will be estimated from ovlStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUTS\n");
    fprintf(stderr, "  -C corStore      output layouts to store 'corStore'\n");
    fprintf(stderr, "  -p prefix        output prefix name, for logging and summary report\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "READ SELECTION\n");
    fprintf(stderr, "  -b bgnID         process reads starting at bgnID\n");
    fprintf(stderr, "  -e endID         process reads up to but not including endID\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "EVIDENCE SELECTION\n");
    fprintf(stderr, "  -eL length       minimum length of evidence overlaps\n");
    fprintf(stderr, "  -eE erate        maximum error rate of evidence overlaps\n");
    fprintf(stderr, "  -ec coverage     minimum coverage needed in evidence reads\n");       //  not used in canu
    fprintf(stderr, "  -eC coverage     maximum coverage of evidence reads to emit\n");
    fprintf(stderr, "  -eM length       minimum length of a corrected read\n");              //  not used in canu
    fprintf(stderr, "\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR: no input gkpStore (-G) supplied.\n");
    if (ovlName == NULL)
      fprintf(stderr, "ERROR: no input ovlStore (-O) supplied.\n");
    if (corName == NULL)
      fprintf(stderr, "ERROR: no output corStore (-C) supplied.\n");
    exit(1);
  }

  //  Open inputs and output tigStore.

  gkStore  *gkpStore = gkStore::gkStore_open(gkpName);
  ovStore  *ovlStore = new ovStore(ovlName, gkpStore);

  tgStore  *corStore = new tgStore(corName);

  uint32    numReads = gkpStore->gkStore_getNumReads();

  //  Load read scores, if supplied.

  uint16   *olapThresh = loadThresholds(gkpStore, ovlStore, scoreName, expectedCoverage);

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

  logFile = AS_UTL_openOutputFile(outputPrefix, "log");
  sumFile = AS_UTL_openOutputFile(outputPrefix, "summary",    false);    //  Never used!

  //  Initialize processing.

  uint32             ovlMax    = 1024 * 1024;
  ovOverlap         *ovl       = ovOverlap::allocateOverlaps(gkpStore, ovlMax);
  uint32             ovlLen    = ovlStore->readOverlaps(ovl, ovlMax, true);

  gkReadData        *readData  = new gkReadData;

  //  And process.

  for (uint32 ii=0; ii<numReads+1; ii++) {
    uint32   readID = (ovlLen > 0) ? ovl[0].a_iid : UINT32_MAX;   //  Read ID of overlaps, or maximum ID if no overlaps.
    tgTig   *layout = new tgTig;

    layout->_tigID           = ii;
    layout->_layoutLen       = gkpStore->gkStore_getRead(readID)->gkRead_sequenceLength();

    assert(ii <= readID);

    //  If ii is below readID, there are no overlaps for this read.  Make an empty placeholder tig for it.
    //  But if ii is readID, we have overlaps, so process them, then load more.

    if (ii == readID) {
      layout = generateLayout(layout,
                              olapThresh,
                              minEvidenceLength, maxEvidenceErate, maxEvidenceCoverage,
                              ovl, ovlLen);
      ovlLen = ovlStore->readOverlaps(ovl, ovlMax, true);
    }

    //  And save the layout into the corStore.

    corStore->insertTig(layout, false);

    delete layout;
  }

  //  Close files and clean up.

  if (logFile != NULL)   fclose(logFile);
  if (sumFile != NULL)   fclose(sumFile);

  delete [] olapThresh;
  delete    readData;
  delete [] ovl;
  delete    corStore;
  delete    ovlStore;

  gkpStore->gkStore_close();

  fprintf(stderr, "Bye.\n");

  return(0);
}


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
#include "intervalList.H"
#include "AS_UTL_reverseComplement.H"
#include "AS_UTL_fasta.H"

#include "falconConsensus.H"

#include <set>

using namespace std;

//  Debugging on which reads are filtered, which are used, and which are removed
//  to meet coverage thresholds.  Very big.
#undef DEBUG_LAYOUT



class readStatus {
public:
  readStatus() {
    numOlaps   = 0;

    origLength = 0;
    corrLength = 0;

    message = NULL;

    usedAsEvidence  = false;
  };
  ~readStatus() {
    delete [] message;
  };

  uint32   numOlaps;

  uint32   origLength;
  uint32   corrLength;

  char    *message;

  bool     usedAsEvidence;
};





uint16 *
loadThresholds(char *scoreName, uint32 numReads) {

  if (scoreName == NULL)
    return(NULL);

  uint16  *olapThresh = new uint16 [numReads + 1];

  errno = 0;
  FILE *S = fopen(scoreName, "r");
  if (errno)
    fprintf(stderr, "failed to open '%s' for reading: %s\n", scoreName, strerror(errno)), exit(1);

  AS_UTL_safeRead(S, olapThresh, "scores", sizeof(uint16), numReads + 1);

  fclose(S);

  return(olapThresh);
}



void
loadReadList(char *readListName, uint32 iidMin, uint32 iidMax, set<uint32> &readList) {
  char  L[1024];

  if (readListName == NULL)
    return;

  errno = 0;
  FILE *R = fopen(readListName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", readListName, strerror(errno)), exit(1);

  fgets(L, 1024, R);
  while (!feof(R)) {
    splitToWords W(L);
    uint32       id = W(0);

    if ((iidMin <= id) &&
        (id     <= iidMax))
      readList.insert(W(0));

    fgets(L, 1024, R);
  }

  fclose(R);
}



FILE *
openOutputFile(char *outputPrefix,
               char *outputSuffix,
               bool  doOpen = true) {
  char   outputName[FILENAME_MAX];

  if (outputPrefix == NULL)
    return(NULL);

  if (doOpen == false)
    return(NULL);

  snprintf(outputName, FILENAME_MAX, "%s.%s", outputPrefix, outputSuffix);

  errno = 0;

  FILE *F = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", outputName, strerror(errno)), exit(1);

  return(F);
}



//  Generate a layout for the read in ovl[0].a_iid, using most or all of the overlaps
//  in ovl.

tgTig *
generateLayout(readStatus *status,
               uint16     *olapThresh,
               uint32      minEvidenceLength,
               double      maxEvidenceErate,
               double      maxEvidenceCoverage,
               ovOverlap *ovl,
               uint32      ovlLen,
               FILE       *flgFile) {
  tgTig         *layout = new tgTig;

  set<uint32_t>  children;

  layout->_tigID           = ovl[0].a_iid;
  layout->_coverageStat    = 1.0;  //  Default to just barely unique
  layout->_microhetProb    = 1.0;  //  Default to 100% probability of unique

  layout->_class           = tgTig_noclass;
  layout->_suggestRepeat   = false;
  layout->_suggestCircular = false;

  layout->_layoutLen       = status[layout->_tigID].origLength;

  resizeArray(layout->_children, layout->_childrenLen, layout->_childrenMax, ovlLen, resizeArray_doNothing);

  if (flgFile)
    fprintf(flgFile, "Generate layout for read " F_U32 " length " F_U32 " using up to " F_U32 " overlaps.\n",
            layout->_tigID, layout->_layoutLen, ovlLen);

  for (uint32 oo=0; oo<ovlLen; oo++) {
    uint64   ovlLength = ovl[oo].b_len();
    uint16   ovlScore  = ovl[oo].overlapScore(true);

    if (ovlLength > AS_MAX_READLEN) {
      char ovlString[1024];
      fprintf(stderr, "ERROR: bogus overlap '%s'\n", ovl[oo].toString(ovlString, ovOverlapAsCoords, false));
    }
    assert(ovlLength < AS_MAX_READLEN);

    if (ovl[oo].erate() > maxEvidenceErate) {
      if (flgFile)
        fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - low quality (threshold %.2f)\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), maxEvidenceErate);
      continue;
    }

    if (ovl[oo].a_end() - ovl[oo].a_bgn() < minEvidenceLength) {
      if (flgFile)
        fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - too short (threshold %u)\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), minEvidenceLength);
      continue;
    }

    if ((olapThresh != NULL) &&
        (ovlScore < olapThresh[ovl[oo].b_iid])) {
      if (flgFile)
        fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - filtered by global filter (threshold " F_U16 ")\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate(), olapThresh[ovl[oo].b_iid]);
      continue;
    }

    if (children.find(ovl[oo].b_iid) != children.end()) {
      if (flgFile)
        fprintf(flgFile, "  filter read %9u at position %6u,%6u length %5lu erate %.3f - duplicate\n",
                ovl[oo].b_iid, ovl[oo].a_bgn(), ovl[oo].a_end(), ovlLength, ovl[oo].erate());
      continue;
    }

    if (flgFile)
      fprintf(flgFile, "  allow  read %9u at position %6u,%6u length %5lu erate %.3f\n",
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



bool
analyzeExpectedConsensus(tgTig      *layout,
                         uint32      minCorLength,
                         uint32      minEvidenceCoverage,
                         readStatus *status) {
  bool    skipIt        = false;
  char    skipMsg[1024] = {0};

  //  Possibly filter by the length of the uncorrected read.

  uint32  readID        = layout->tigID();

  if (status[readID].origLength < minCorLength) {
    strcat(skipMsg, "\tread_too_short");
    skipIt = true;
  }

  //  Possibly filter by the length of the corrected read, taking into account depth of coverage.

  intervalList<int32>   coverage;

  for (uint32 ii=0; ii<layout->numberOfChildren(); ii++) {
    tgPosition *pos = layout->getChild(ii);

    coverage.add(pos->_min, pos->_max - pos->_min);
  }

  intervalList<int32>   depth(coverage);

  int32    bgn       = INT32_MAX;
  int32    corLen    = 0;

  for (uint32 dd=0; dd<depth.numberOfIntervals(); dd++) {
    if (depth.depth(dd) < minEvidenceCoverage) {
      bgn = INT32_MAX;
      continue;
    }

    if (bgn == INT32_MAX)
      bgn = depth.lo(dd);

    if (corLen < depth.hi(dd) - bgn)
      corLen = depth.hi(dd) - bgn;
  }

  if (corLen < minCorLength) {
    strcat(skipMsg, "\tcorrection_too_short");
    skipIt = true;
  }

  //  Filter out empty tigs - these either have no overlaps, or failed the
  //  length check in generateLayout.

  if (layout->numberOfChildren() == 0) {
    strcat(skipMsg, "\tno_children");
    skipIt = true;
  }

  status[readID].numOlaps     = layout->numberOfChildren();
  status[readID].corrLength   = corLen;

  if (skipMsg[0] != 0)
    status[readID].message = duplicateString(skipMsg);

  for (uint32 ii=0; ii<layout->numberOfChildren(); ii++)
    status[layout->getChild(ii)->ident()].usedAsEvidence = true;

  return(skipIt == false);
}



//  A mash up of falcon_sense.C and outputFalcon.C

void
generateFalconConsensus(falconConsensus   *fc,
                        gkStore           *gkpStore,
                        tgTig             *tig,
                        bool               trimToAlign,
                        FILE              *F,
                        gkReadData        *readData,
                        uint32             minOutputLength) {

  //  Grab and save the raw read for the template.

  fprintf(stderr, "Processing read %u of length %u with %u evidence reads.\n",
          tig->tigID(), tig->length(), tig->numberOfChildren());

  gkpStore->gkStore_loadReadData(tig->tigID(), readData);

  //  Now parse the layout and push all the sequences onto our seqs vector.

  falconInput   *evidence = new falconInput [tig->numberOfChildren() + 1];

  evidence[0].addInput(tig->tigID(),
                       readData->gkReadData_getSequence(),
                       readData->gkReadData_getRead()->gkRead_sequenceLength(),
                       0,
                       readData->gkReadData_getRead()->gkRead_sequenceLength());


  for (uint32 cc=0; cc<tig->numberOfChildren(); cc++) {
    tgPosition  *child = tig->getChild(cc);

    gkpStore->gkStore_loadReadData(child->ident(), readData);

    if (child->isReverse())
      reverseComplementSequence(readData->gkReadData_getSequence(),
                                readData->gkReadData_getRead()->gkRead_sequenceLength());

    //  For debugging/testing, skip one orientation of overlap.
    //
    //if (child->isReverse() == false)
    //  continue;
    //if (child->isReverse() == true)
    //  continue;

    //  Trim the read to the aligned bit
    char   *seq    = readData->gkReadData_getSequence();
    uint32  seqLen = readData->gkReadData_getRead()->gkRead_sequenceLength();

    if (trimToAlign) {
      seq    += child->askip();
      seqLen -= child->askip() + child->bskip();

      seq[seqLen] = 0;
    }

    //  Used to skip if read length was less or equal to min_ovl_len

    evidence[cc+1].addInput(child->ident(), seq, seqLen, child->min(), child->max());
  }

  //  Loaded all reads, build consensus.

  uint32 splitSeqID = 0;

  //FConsensus::consensus_data *consensus_data_ptr = FConsensus::generate_consensus( seqs, min_cov, min_idt, min_ovl_len, max_read_len );

  falconData  *fd = fc->generateConsensus(evidence,
                                          tig->numberOfChildren() + 1);

#ifdef TRACK_POSITIONS
  //const std::string& sequenceToCorrect = seqs.at(0);
  char * originalStringPointer = consensus_data_ptr->sequence;
#endif

  char * split = strtok(fd->seq, "acgt");

  while (split != NULL) {
    if (strlen(split) > minOutputLength) {
      fprintf(stderr, "Generated read %u_%u of length %lu.\n",
              tig->tigID(), splitSeqID, strlen(split));

      AS_UTL_writeFastA(F, split, strlen(split), 60, ">read%u_%d\n", tig->tigID(), splitSeqID);

      splitSeqID++;

#ifdef TRACK_POSITIONS
      int distance_from_beginning = split - originalStringPointer;
      std::vector<int> relevantOriginalPositions(consensus_data_ptr->originalPos.begin() + distance_from_beginning, consensus_data_ptr->originalPos.begin() + distance_from_beginning + strlen(split));
      int firstRelevantPosition = relevantOriginalPositions.front();
      int lastRelevantPosition = relevantOriginalPositions.back();

      std::string relevantOriginalTemplate = seqs.at(0).substr(firstRelevantPosition, lastRelevantPosition - firstRelevantPosition + 1);

      // store relevantOriginalTemplate along with corrected read - not implemented
#endif
    }

    split = strtok(NULL, "acgt");
  }

  delete fd;  //FConsensus::free_consensus_data( consensus_data_ptr );
  delete [] evidence;
}



void
estimateMemoryUsage(gkStore       *gkpStore,
                    uint32         iidMin,
                    uint32         iidMax,
                    set<uint32>   &readList,
                    ovStore       *ovlStore) {

  uint32   maxReadID  = 0;    //  Find the longest read and use that to estimate
  uint32   maxReadLen = 0;    //  the maximum memory usage.

  for (uint32 rr=iidMin; rr<iidMax; rr++) {
    if ((readList.size() > 0) &&
        (readList.count(rr) == 0))
      continue;

    gkRead *read = gkpStore->gkStore_getRead(rr);
    uint32  rLen = read->gkRead_sequenceLength();

    if (maxReadLen < rLen) {
      maxReadID  = rr;
      maxReadLen = rLen;
    }
  }

  //  Estimate the number of evidence reads it'll have.  Just assume all overlaps are used and
  //  that every error is an indel.

  ovlStore->setRange(maxReadID, maxReadID);

  ovOverlap   *overlaps    = NULL;
  uint32       overlapsLen = 0;
  uint32       overlapsMax = 0;

  uint64       nOvl       = ovlStore->numOverlapsInRange();
  uint32       nOvlLoaded = ovlStore->readOverlaps(maxReadID, overlaps, overlapsLen, overlapsMax);

  assert(nOvl == nOvlLoaded);

  uint64       nBasesInOlaps = 0;

  for (uint32 oo=0; oo<nOvl; oo++) {
    assert(overlaps[oo].a_bgn() < overlaps[oo].a_end());

    nBasesInOlaps += (overlaps[oo].a_end() - overlaps[oo].a_bgn()) * (1 + overlaps[oo].erate());
  }

  //  Throw that at falconConsensus and let it decide memory usage.

  falconConsensus *fc = new falconConsensus(0, 0.0, 0);

  uint64 mem = fc->estimateMemoryUsage(nOvl, nBasesInOlaps, maxReadLen);

  fprintf(stdout, "Based on read %u of length %u with %lu overlaps covering %lu bases, expecting to use %lu bytes, %.3f GB memory usage\n",
          maxReadID, maxReadLen, nOvl, nBasesInOlaps,
          mem, mem / 1024.0 / 1024.0 / 1024.0);
}



int
main(int argc, char **argv) {
  char             *gkpName   = 0L;
  char             *ovlName   = 0L;
  char             *scoreName = 0L;
  char             *tigName   = 0L;

  bool              falconOutput    = false;  //  To stdout
  bool              consensusOutput = false;
  bool              trimToAlign     = true;

  bool              estimateMemory  = false;

  uint32            errorRate = AS_OVS_encodeEvalue(0.015);

  char             *outputPrefix = NULL;
  char              logName[FILENAME_MAX] = {0};
  char              sumName[FILENAME_MAX] = {0};
  char              flgName[FILENAME_MAX] = {0};
  FILE             *logFile = 0L;
  FILE             *sumFile = 0L;
  FILE             *flgFile = 0L;

  uint32            minEvidenceOverlap  = 40;
  uint32            minEvidenceCoverage = 4;

  uint32            iidMin       = 0;
  uint32            iidMax       = UINT32_MAX;
  char             *readListName = NULL;
  set<uint32>       readList;

  uint32            minEvidenceLength   = 0;
  double            maxEvidenceErate    = 1.0;
  double            maxEvidenceCoverage = DBL_MAX;

  uint32            minCorLength        = 0;

  //  Consensus parameters

  uint32            numThreads         = 1;
  uint32            minAllowedCoverage = 4;
  double            minIdentity        = 0.5;
  uint32            minOutputLength    = 500;

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

    } else if (strcmp(argv[arg], "-p") == 0) {
      outputPrefix = argv[++arg];


    } else if (strcmp(argv[arg], "-t") == 0) {   //  COMPUTE RESOURCES
      numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-M") == 0) {
      estimateMemory = true;


    } else if (strcmp(argv[arg], "-T") == 0) {   //  OUTPUT FORMAT
      tigName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      falconOutput = true;

    } else if (strcmp(argv[arg], "-C") == 0) {
      consensusOutput = true;


    } else if (strcmp(argv[arg], "-b") == 0) {   //  READ SELECTION
      iidMin  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      iidMax  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-r") == 0) {
      readListName = argv[++arg];


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


    } else if (strcmp(argv[arg], "-cc") == 0) {   //  CONSENSUS
      minAllowedCoverage = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-cl") == 0) {
      minOutputLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-ci") == 0) {
      minIdentity = atoi(argv[++arg]);


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
    fprintf(stderr, "INPUTS (all mandatory)\n");
    fprintf(stderr, "  -G gkpStore      mandatory path to gkpStore\n");
    fprintf(stderr, "  -O ovlStore      mandatory path to ovlStore\n");
    fprintf(stderr, "  -S file          overlap score thresholds (from filterCorrectionOverlaps)\n");
    fprintf(stderr, "  -p prefix        output prefix name, for logging and summary report\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "RESOURCE PARAMETERS\n");
    fprintf(stderr, "  -t numThreads    number of compute threads to use\n");
    fprintf(stderr, "  -M               estimate memory requirements for generating corrected reads (-C)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUT FORMAT\n");
    fprintf(stderr, "  -T store         output layouts to tigStore 'store'\n");
    fprintf(stderr, "  -F               output falconsense-style input directly to stdout (OBSOLETE)\n");
    fprintf(stderr, "  -C               output corrected read consensus sequences to stdout\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "READ SELECTION\n");
    fprintf(stderr, "  -b bgnID         process reads starting at bgnID\n");
    fprintf(stderr, "  -e endID         process reads up to but not including endID\n");
    fprintf(stderr, "  -r file          only process reads listed (must also be bgnID <= ID < endOD)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "EVIDENCE SELECTION\n");
    fprintf(stderr, "  -eL length       minimum length of evidence overlaps\n");
    fprintf(stderr, "  -eE erate        maximum error rate of evidence overlaps\n");
    fprintf(stderr, "  -ec coverage     minimum coverage needed in evidence reads\n");       //  not used in canu
    fprintf(stderr, "  -eC coverage     maximum coverage of evidence reads to emit\n");
    fprintf(stderr, "  -eM length       minimum length of a corrected read\n");              //  not used in canu
    fprintf(stderr, "\n");
    fprintf(stderr, "CONSENSUS PARAMETERS\n");
    fprintf(stderr, "  -cc coverage     minimum consensus coverage to output corrected base\n");
    fprintf(stderr, "  -cl length       minimum length of corrected region\n");
    fprintf(stderr, "  -ci identity     minimum identity of an aligned evidence read\n");
    fprintf(stderr, "\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR: no gkpStore input (-G) supplied.\n");
    if (ovlName == NULL)
      fprintf(stderr, "ERROR: no ovlStore input (-O) supplied.\n");

    exit(1);
  }


  omp_set_num_threads(numThreads);



  //  Open inputs and output tigStore.

  gkStore  *gkpStore = gkStore::gkStore_open(gkpName);
  ovStore  *ovlStore = new ovStore(ovlName, gkpStore);
  tgStore  *tigStore = (tigName != NULL) ? new tgStore(tigName) : NULL;

  uint32    numReads = gkpStore->gkStore_getNumReads();

  //  Load read scores, if supplied.

  uint16   *olapThresh = loadThresholds(scoreName, numReads);

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

  //  If a readList is supplied, load it, respecting the iidMin/iidMax (only to cut down on the
  //  size).

  loadReadList(readListName, iidMin, iidMax, readList);

  //  Open logging and summary files

  logFile = openOutputFile(outputPrefix, "log");
  sumFile = openOutputFile(outputPrefix, "summary",    false);    //  Never used!
  flgFile = openOutputFile(outputPrefix, "filter.log", false);    //  Too much crud!

  //  Estimate memory?

  if (estimateMemory == true) {
    estimateMemoryUsage(gkpStore, iidMin, iidMax, readList, ovlStore);
    exit(0);  //  And exit rather ungracefully; just too many parameters to functionify the rest of main().
  }

  //  Initialize processing.

  uint32             ovlMax = 1024 * 1024;
  uint32             ovlLen = 0;
  ovOverlap         *ovl    = ovOverlap::allocateOverlaps(gkpStore, ovlMax);

  gkReadData        *readData = new gkReadData;

  falconConsensus   *fc = (consensusOutput == false) ? NULL : new falconConsensus(minAllowedCoverage, minIdentity, minOutputLength);

  readStatus        *status = new readStatus [numReads + 1];

  for (uint32 ii=1; ii<numReads; ii++)
    status[ii].origLength = gkpStore->gkStore_getRead(ii)->gkRead_sequenceLength();

  //  And process.

  for (ovlLen = ovlStore->readOverlaps(ovl, ovlMax, true);
       ovlLen > 0;
       ovlLen = ovlStore->readOverlaps(ovl, ovlMax, true)) {

    if ((readList.size() > 0) &&                     //  Skip reads not on the read list.  We need
        (readList.count(ovl[0].a_iid) == 0))         //  to read the olaps, but don't need to make
      continue;                                      //  or process the layout.

    tgTig *layout = generateLayout(status,
                                   olapThresh,
                                   minEvidenceLength, maxEvidenceErate, maxEvidenceCoverage,
                                   ovl, ovlLen,
                                   flgFile);

    if (analyzeExpectedConsensus(layout,             //  If the analysis is good, write outputs.
                                 minCorLength,       //  Otherwise, don't.
                                 minEvidenceCoverage,
                                 status) == true) {
      if (tigStore != NULL)
        tigStore->insertTig(layout, false);

      if (falconOutput == true)
        outputFalcon(gkpStore, layout, trimToAlign, stdout, readData);

      if (consensusOutput == true)
        generateFalconConsensus(fc, gkpStore, layout, trimToAlign, stdout, readData, minOutputLength);
    }

    delete layout;
  }

  //  Terminate output files and write the logging.

  if (falconOutput) {
    fprintf(stdout, "- -\n");
  }

  if (logFile) {
    fprintf(logFile, "read\torigLen\tevidence\tcorrLen\tused\n");

    for (uint32 ii=1; ii<numReads; ii++)
      fprintf(logFile, "%u\t%u\t%u\t%u\t%u%s\n",
              ii,
              status[ii].origLength,
              status[ii].numOlaps,
              status[ii].corrLength,
              status[ii].usedAsEvidence,
              (status[ii].message) ? status[ii].message : "");
  }

  //  Close files and clean up.

  if (logFile != NULL)   fclose(logFile);
  if (sumFile != NULL)   fclose(sumFile);
  if (flgFile != NULL)   fclose(flgFile);

  delete [] status;
  delete    fc;
  delete [] olapThresh;
  delete    readData;
  delete [] ovl;
  delete    tigStore;
  delete    ovlStore;

  gkpStore->gkStore_close();

  return(0);
}

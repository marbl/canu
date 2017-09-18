
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
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"

#include "falconConsensus.H"
//#include "computeGlobalScore.H"

#include "intervalList.H"

#include <vector>
#include <algorithm>

using namespace std;




class readStatus {
public:
  readStatus() {
    readID             = 0;

    numOlaps           = 0;

    origLength         = UINT32_MAX;
    corrLength         = UINT32_MAX;   //  Makes zeroth read first when sorted by length.

    memoryRequired     = 0;

    usedForEvidence    = false;
    usedForCorrection  = false;
    rescued            = false;
  };
  ~readStatus() {
  };

  uint32   readID;

  uint32   numOlaps;

  uint32   origLength;
  uint32   corrLength;

  uint32   memoryRequired;

  bool     usedForEvidence;
  bool     usedForCorrection;
  bool     rescued;
};

bool
sortByCorLength(readStatus const &a, readStatus const &b) {
  return(a.corrLength > b.corrLength);
}

bool
sortByReadID(readStatus const &a, readStatus const &b) {
  return(a.readID < b.readID);
}





class lengthStats {
public:
  lengthStats(uint32 maxN) {
    L = new uint32 [maxN];
    N = 0;

    nReads  = 0;
    nBases  = 0;
    median  = 0;
    mean    = 0;
    n50     = 0;
    minimum = 0;
    maximum = 0;
  };
  ~lengthStats(void) {
    delete [] L;
  };

  void     add(uint32 length) {
    L[N++] = length;
  };

  void     compute(void) {

    if (N == 0)
      return;

    sort(L, L+N);

    nReads  = N;
    nBases  = 0;
    median  = L[N/2];
    mean    = 0;
    n50     = 0;
    minimum = L[0];
    maximum = L[N-1];

    for (uint32 ii=0; ii<N; ii++)
      nBases += L[ii];

    mean = nBases / N;

    uint64  ss = 0;
    uint32  ii = 0;

    while (ss < nBases/2)
      ss += L[ii++];

    n50 = L[ii];
  };


  uint32  *L;
  uint32   N;

  uint32   nReads;
  uint64   nBases;
  uint32   median;
  uint32   mean;
  uint32   n50;
  uint32   minimum;
  uint32   maximum;
};




void
analyzeLength(tgTig            *layout,
              uint32     UNUSED(minCorLength),
              uint32            minEvidenceCoverage,
              readStatus       *status,
              falconConsensus  *fc) {
  uint32  readID        = layout->tigID();

  //  Estimate the length of the corrected read, using overlap position and depth.

  intervalList<int32>   coverage;
  uint64                basesInOlaps = 0;

  for (uint32 ii=0; ii<layout->numberOfChildren(); ii++) {
    tgPosition *pos = layout->getChild(ii);

    coverage.add(pos->_min, pos->_max - pos->_min);

    basesInOlaps += pos->_max - pos->_min;
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

  //  Save our results.

  status[readID].readID         = readID;

  status[readID].numOlaps       = layout->numberOfChildren();

  status[readID].origLength     = layout->length();
  status[readID].corrLength     = corLen;

  status[readID].memoryRequired = fc->estimateMemoryUsage(layout->numberOfChildren(),
                                                          basesInOlaps,
                                                          layout->length());
}



void
markEvidence(tgTig            *layout,
             readStatus       *status) {
  uint32  readID        = layout->tigID();

  //  If not used for correction, don't flag the evidence!

  if (status[readID].usedForCorrection == false)
    return;

  //  Otherwise, used for correction, so flag the evidence.

  for (uint32 ii=0; ii<layout->numberOfChildren(); ii++)
    status[layout->getChild(ii)->ident()].usedForEvidence = true;
}



void
dumpStatistics(FILE *F, readStatus *status, uint32 numReads) {

  lengthStats  lenRaw(numReads+1);        //  All raw reads

  lengthStats  lenNoOlaps(numReads+1);    //  No overlaps at all.

  lengthStats  lenEvidence(numReads+1);   //  Reads used as evidence

  lengthStats  lenCorrOrig(numReads+1);   //  Reads corrected, both original and ocrrected length
  lengthStats  lenCorr(numReads+1);

  lengthStats  lenRescOrig(numReads+1);   //  Reads rescued, both original and ocrrected length
  lengthStats  lenResc(numReads+1);

  lengthStats  lenNoCoOrig(numReads+1);   //  Reads not corrected, both original and ocrrected length
  lengthStats  lenNoCo(numReads+1);

  uint64       maxMemory = 0;

  for (uint32 ti=1; ti<numReads+1; ti++) {
    uint32  O = status[ti].origLength;
    uint32  C = status[ti].corrLength;

    if (status[ti].numOlaps > 0)
      lenRaw.add(O);
    else
      lenNoOlaps.add(O);

    if (status[ti].usedForEvidence == true) {
      lenEvidence.add(O);
    }

    if (status[ti].usedForCorrection == true) {
      lenCorrOrig.add(O);
      lenCorr.add(C);
    }

    if (status[ti].rescued == true) {
      lenRescOrig.add(O);
      lenResc.add(C);
    }

    if ((status[ti].usedForCorrection == false) &&
        (status[ti].rescued == false)) {
      lenNoCoOrig.add(O);
      lenNoCo.add(C);
    }

    if (maxMemory < status[ti].memoryRequired)
      maxMemory = status[ti].memoryRequired;
  }

  lenRaw.compute();

  lenNoOlaps.compute();

  lenEvidence.compute();

  lenCorrOrig.compute();
  lenCorr.compute();

  lenRescOrig.compute();
  lenResc.compute();

  lenNoCoOrig.compute();
  lenNoCo.compute();


#define U32FORMAT  "%13" F_U32P
#define U64FORMAT  "%13" F_U64P

#if 0
  fprintf(F, "                          original      original                 --------corrected---------  ----------rescued----------  --------uncorrected--------\n");
  fprintf(F, "                               raw       with no      evidence                     expected                     expected                     expected\n");
  fprintf(F, "category                     reads      overlaps         reads            raw     corrected            raw     corrected            raw     corrected\n");
  fprintf(F, "-------------------- ------------- ------------- -------------  ------------- -------------  ------------- -------------  ------------- -------------\n");
  fprintf(F, "Number of Reads      " U32FORMAT " " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenRaw.nReads,  lenNoOlaps.nReads,  lenEvidence.nReads,  lenCorrOrig.nReads,  lenCorr.nReads,  lenRescOrig.nReads,  lenResc.nReads,  lenNoCoOrig.nReads,  lenNoCo.nReads);
  fprintf(F, "Number of Bases      " U64FORMAT " " U64FORMAT " " U64FORMAT "  " U64FORMAT " " U64FORMAT "  " U64FORMAT " " U64FORMAT "  " U64FORMAT " " U64FORMAT "\n", lenRaw.nBases,  lenNoOlaps.nBases,  lenEvidence.nBases,  lenCorrOrig.nBases,  lenCorr.nBases,  lenRescOrig.nBases,  lenResc.nBases,  lenNoCoOrig.nBases,  lenNoCo.nBases);
  fprintf(F, "Median               " U32FORMAT " " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenRaw.median,  lenNoOlaps.median,  lenEvidence.median,  lenCorrOrig.median,  lenCorr.median,  lenRescOrig.median,  lenResc.median,  lenNoCoOrig.median,  lenNoCo.median);
  fprintf(F, "Mean                 " U32FORMAT " " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenRaw.mean,    lenNoOlaps.mean,    lenEvidence.mean,    lenCorrOrig.mean,    lenCorr.mean,    lenRescOrig.mean,    lenResc.mean,    lenNoCoOrig.mean,    lenNoCo.mean);
  fprintf(F, "N50                  " U32FORMAT " " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenRaw.n50,     lenNoOlaps.n50,     lenEvidence.n50,     lenCorrOrig.n50,     lenCorr.n50,     lenRescOrig.n50,     lenResc.n50,     lenNoCoOrig.n50,     lenNoCo.n50);
  fprintf(F, "Minimum              " U32FORMAT " " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenRaw.minimum, lenNoOlaps.minimum, lenEvidence.minimum, lenCorrOrig.minimum, lenCorr.minimum, lenRescOrig.minimum, lenResc.minimum, lenNoCoOrig.minimum, lenNoCo.minimum);
  fprintf(F, "Maximum              " U32FORMAT " " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenRaw.maximum, lenNoOlaps.maximum, lenEvidence.maximum, lenCorrOrig.maximum, lenCorr.maximum, lenRescOrig.maximum, lenResc.maximum, lenNoCoOrig.maximum, lenNoCo.maximum);
  fprintf(F, "\n");
  fprintf(F, "Maximum Memory       " U64FORMAT "\n", maxMemory);
#endif

  fprintf(F, "                          original      original\n");
  fprintf(F, "                         raw reads     raw reads\n");
  fprintf(F, "category                w/overlaps  w/o/overlaps\n");
  fprintf(F, "-------------------- ------------- -------------\n");
  fprintf(F, "Number of Reads      " U32FORMAT " " U32FORMAT "\n", lenRaw.nReads,  lenNoOlaps.nReads);
  fprintf(F, "Number of Bases      " U64FORMAT " " U64FORMAT "\n", lenRaw.nBases,  lenNoOlaps.nBases);
  fprintf(F, "Median               " U32FORMAT " " U32FORMAT "\n", lenRaw.median,  lenNoOlaps.median);
  fprintf(F, "Mean                 " U32FORMAT " " U32FORMAT "\n", lenRaw.mean,    lenNoOlaps.mean);
  fprintf(F, "N50                  " U32FORMAT " " U32FORMAT "\n", lenRaw.n50,     lenNoOlaps.n50);
  fprintf(F, "Minimum              " U32FORMAT " " U32FORMAT "\n", lenRaw.minimum, lenNoOlaps.minimum);
  fprintf(F, "Maximum              " U32FORMAT " " U32FORMAT "\n", lenRaw.maximum, lenNoOlaps.maximum);
  fprintf(F, "\n");
  fprintf(F, "                                     --------corrected---------  ----------rescued----------\n");
  fprintf(F, "                          evidence                     expected                     expected\n");
  fprintf(F, "category                     reads            raw     corrected            raw     corrected\n");
  fprintf(F, "-------------------- -------------  ------------- -------------  ------------- -------------\n");
  fprintf(F, "Number of Reads      " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenEvidence.nReads,  lenCorrOrig.nReads,  lenCorr.nReads,  lenRescOrig.nReads,  lenResc.nReads);
  fprintf(F, "Number of Bases      " U64FORMAT "  " U64FORMAT " " U64FORMAT "  " U64FORMAT " " U64FORMAT "\n", lenEvidence.nBases,  lenCorrOrig.nBases,  lenCorr.nBases,  lenRescOrig.nBases,  lenResc.nBases);
  fprintf(F, "Median               " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenEvidence.median,  lenCorrOrig.median,  lenCorr.median,  lenRescOrig.median,  lenResc.median);
  fprintf(F, "Mean                 " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenEvidence.mean,    lenCorrOrig.mean,    lenCorr.mean,    lenRescOrig.mean,    lenResc.mean);
  fprintf(F, "N50                  " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenEvidence.n50,     lenCorrOrig.n50,     lenCorr.n50,     lenRescOrig.n50,     lenResc.n50);
  fprintf(F, "Minimum              " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenEvidence.minimum, lenCorrOrig.minimum, lenCorr.minimum, lenRescOrig.minimum, lenResc.minimum);
  fprintf(F, "Maximum              " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenEvidence.maximum, lenCorrOrig.maximum, lenCorr.maximum, lenRescOrig.maximum, lenResc.maximum);
  fprintf(F, "\n");
  fprintf(F, "                     --------uncorrected--------\n");
  fprintf(F, "                                        expected\n");
  fprintf(F, "category                       raw     corrected\n");
  fprintf(F, "-------------------- ------------- -------------\n");
  fprintf(F, "Number of Reads      " U32FORMAT " " U32FORMAT "\n", lenNoCoOrig.nReads,  lenNoCo.nReads);
  fprintf(F, "Number of Bases      " U64FORMAT " " U64FORMAT "\n", lenNoCoOrig.nBases,  lenNoCo.nBases);
  fprintf(F, "Median               " U32FORMAT " " U32FORMAT "\n", lenNoCoOrig.median,  lenNoCo.median);
  fprintf(F, "Mean                 " U32FORMAT " " U32FORMAT "\n", lenNoCoOrig.mean,    lenNoCo.mean);
  fprintf(F, "N50                  " U32FORMAT " " U32FORMAT "\n", lenNoCoOrig.n50,     lenNoCo.n50);
  fprintf(F, "Minimum              " U32FORMAT " " U32FORMAT "\n", lenNoCoOrig.minimum, lenNoCo.minimum);
  fprintf(F, "Maximum              " U32FORMAT " " U32FORMAT "\n", lenNoCoOrig.maximum, lenNoCo.maximum);
  fprintf(F, "\n");
  fprintf(F, "Maximum Memory       " U64FORMAT "\n", maxMemory);
}



void
dumpLog(FILE *F, readStatus *status, uint32 numReads) {

  fprintf(F, "readID          numOlaps    origLength    corrLength        memory  used\n");
  fprintf(F, "---------- ------------- ------------- ------------- ------------- -----\n");
  for (uint32 ti=1; ti<numReads+1; ti++)
    fprintf(F, "%-10" F_U32P " " U32FORMAT " " U32FORMAT " " U32FORMAT " " U32FORMAT " %c %c %c\n",
            status[ti].readID,
            status[ti].numOlaps,
            status[ti].origLength,
            status[ti].corrLength,
            status[ti].memoryRequired,
            status[ti].usedForEvidence   ? 'e' : '-',
            status[ti].usedForCorrection ? 'c' : '-',
            status[ti].rescued           ? 'r' : '-');
}




int
main(int argc, char **argv) {
  char           *gkpStoreName     = NULL;
  char           *corStoreName     = NULL;
  char           *outName          = NULL;

  bool            filterNone       = false;
  bool            filterStandard   = true;

  uint32          evidenceCoverage = 4;
  uint32          minLength        = 500;

  uint64          genomeSize       = 0;
  uint32          outCoverage      = 40;

  argc = AS_configure(argc, argv);

  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      corStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-R") == 0) {
      outName = argv[++arg];


    } else if (strcmp(argv[arg], "-all") == 0) {
      filterNone     = true;
      filterStandard = false;

    } else if (strcmp(argv[arg], "-normal") == 0) {
      filterNone     = false;
      filterStandard = true;


    } else if (strcmp(argv[arg], "-cc") == 0) {
      evidenceCoverage = strtoul(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-cl") == 0) {
      minLength = strtoul(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-g") == 0) {
      genomeSize = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-c") == 0) {
      outCoverage = strtoul(argv[++arg], NULL, 10);


    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if (gkpStoreName == NULL)
    err++;
  if (corStoreName == NULL)
    err++;
  if (outName == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [options]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Examines correction layouts in corStore, decides which ones to correct.\n");
    fprintf(stderr, "Writes output (-r) *.readsToCorrect.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS and OUTPUTS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G gkpStore              input reads\n");
    fprintf(stderr, "  -C corStore              input correction layouts\n");
    fprintf(stderr, "  -R asm.readsToCorrect    output ascii list of read IDs to correct\n");
    fprintf(stderr, "                           also creates\n");
    fprintf(stderr, "                             asm.readsToCorrect.stats and\n");
    fprintf(stderr, "                             asm.readsToCorrect.log\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "FILTERING STRATEGY and PARAMETERS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -all                     no filtering, correct all reads (NOT IMPLEMENTED)\n");
    fprintf(stderr, "  -normal                  correct longest expected corrected reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -cc                      minimum coverage of evidence reads\n");
    fprintf(stderr, "  -cl                      minimum length of a corrected read\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -g                       estimated genome size\n");
    fprintf(stderr, "  -c                       desired coverage in corrected reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "RESCUE\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -rescue                  enable rescue - if read not used as evidence\n");
    fprintf(stderr, "                           force it to be corrected\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gatekeeper store (-G) supplied.\n");
    if (corStoreName == NULL)
      fprintf(stderr, "ERROR: no correction store (-C) supplied.\n");
    if (outName == NULL)
      fprintf(stderr, "ERROR: no output (-R) supplied.\n");

    exit(1);
  }

  gkStore          *gkpStore = gkStore::gkStore_open(gkpStoreName);
  tgStore          *corStore = new tgStore(corStoreName, 1);

  falconConsensus  *fc       = new falconConsensus(0, 0, 0);  //  For memory estimtes

  uint32            numReads = gkpStore->gkStore_getNumReads();

  readStatus       *status   = new readStatus [numReads + 1];

  FILE             *roc      = AS_UTL_openOutputFile(outName);
  FILE             *stats    = AS_UTL_openOutputFile(outName, "stats");
  FILE             *log      = AS_UTL_openOutputFile(outName, "log");

  //  Scan the tigs, computing expected corrected length.

  for (uint32 ti=1; ti<corStore->numTigs(); ti++) {
    tgTig  *layout = corStore->loadTig(ti);

    analyzeLength(layout, minLength, evidenceCoverage, status, fc);

    corStore->unloadTig(layout->tigID());
  }

  //  Sort by expected corrected length, then mark reads for correction until we get the desired
  //  outCoverage.  Zeroth read has max corrected length, so remains first in sorted list.

  sort(status, status + numReads+1, sortByCorLength);

  uint64   desiredLength = genomeSize * outCoverage;
  uint64   corrLength    = 0;

  fprintf(stderr, "for genomeSize %lu and coverage %u -> %lu bases\n", genomeSize, outCoverage, desiredLength);

  for (uint32 ii=1; ii<numReads+1; ii++) {
    if (status[ii].corrLength == 0)
      continue;

    status[ii].usedForCorrection = (corrLength < desiredLength);

    corrLength += status[ii].corrLength;
  }

  //  Sort by ID.

  sort(status, status + numReads+1, sortByReadID);

  //  Scan the tigs again, this time marking reads used as evidence in the corrected reads.

  for (uint32 ti=1; ti<corStore->numTigs(); ti++) {
    tgTig  *layout = corStore->loadTig(ti);

    markEvidence(layout, status);

    corStore->unloadTig(layout->tigID());
  }

  //  And finally, flag any read for correction if it isn't already used as evidence or being corrected.

  for (uint32 ti=1; ti<numReads+1; ti++) {
    if ((status[ti].usedForCorrection == false) &&
        (status[ti].usedForEvidence   == false) &&
        (status[ti].corrLength         > minLength))
      status[ti].rescued = true;
  }

  //  Output the list of reads to correct.

  for (uint32 ti=1; ti<numReads+1; ti++)
    if ((status[ti].usedForCorrection == true) ||
        (status[ti].rescued           == true))
      fprintf(roc, "%u\n", status[ti].readID);

  fclose(roc);

  //  Write some statistics and logs.

  dumpStatistics(stats, status, numReads);
  dumpLog(log, status, numReads);

  fclose(stats);
  fclose(log);

  //  And say goodbye.

  delete [] status;

  fprintf(stderr, "Bye.\n");

  exit(0);
}

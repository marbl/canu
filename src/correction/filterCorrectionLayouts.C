
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


#include "sqStore.H"
#include "tgStore.H"

#include "falconConsensus.H"
//#include "computeGlobalScore.H"

#include "intervals.H"

#include <vector>
#include <algorithm>


class readStatus {
public:
  uint32   readID            = 0;

  uint32   numOlaps          = 0;

  uint32   origLength        = 0;
  uint32   corrLength        = 0;

  uint64   memoryRequired    = 0;

  bool     usedForEvidence   = false;
  bool     usedForCorrection = false;
  bool     rescued           = false;
};


class lengthStats {
public:
  lengthStats(uint32 maxN)    {  L = new uint32 [maxN];  };
  ~lengthStats(void)          {  delete [] L;            };

  void     add(uint32 length) {  L[N++] = length;        };

  void     compute(uint64 genomeSize) {

    if (N == 0)
      return;

    std::sort(L, L+N);

    nReads  = N;
    nBases  = 0;

    for (uint32 ii=0; ii<N; ii++)
      nBases += L[ii];

    cov    = nBases / (double)genomeSize;
    median = L[N/2];
    mean   = nBases / N;

    uint64  ss = 0;
    uint32  ii = 0;

    while (ss < nBases/2)
      ss += L[ii++];

    n50     = L[ii];
    minimum = L[0];
    maximum = L[N-1];
  };


  uint32  *L        = nullptr;
  uint32   N        = 0;

  uint32   nReads  = 0;
  uint64   nBases  = 0;
  double   cov     = 0;
  uint32   median  = 0;
  uint32   mean    = 0;
  uint32   n50     = 0;
  uint32   minimum = 0;
  uint32   maximum = 0;
};



void
dumpStatistics(FILE *F, readStatus *status, uint32 numReads, uint64 genomeSize) {

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

  lenRaw.compute(genomeSize);

  lenNoOlaps.compute(genomeSize);

  lenEvidence.compute(genomeSize);

  lenCorrOrig.compute(genomeSize);
  lenCorr.compute(genomeSize);

  lenRescOrig.compute(genomeSize);
  lenResc.compute(genomeSize);

  lenNoCoOrig.compute(genomeSize);
  lenNoCo.compute(genomeSize);


#define U32FORMAT  "%13" F_U32P
#define U64FORMAT  "%13" F_U64P
#define FLTFORMAT  "%13.3f"
#if 0
  fprintf(F, "                          original      original                 --------corrected---------  ----------rescued----------  --------uncorrected--------\n");
  fprintf(F, "                               raw       with no      evidence                     expected                     expected                     expected\n");
  fprintf(F, "category                     reads      overlaps         reads            raw     corrected            raw     corrected            raw     corrected\n");
  fprintf(F, "-------------------- ------------- ------------- -------------  ------------- -------------  ------------- -------------  ------------- -------------\n");
  fprintf(F, "Number of Reads      " U32FORMAT " " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "  " U32FORMAT " " U32FORMAT "\n", lenRaw.nReads,  lenNoOlaps.nReads,  lenEvidence.nReads,  lenCorrOrig.nReads,  lenCorr.nReads,  lenRescOrig.nReads,  lenResc.nReads,  lenNoCoOrig.nReads,  lenNoCo.nReads);
  fprintf(F, "Number of Bases      " U64FORMAT " " U64FORMAT " " U64FORMAT "  " U64FORMAT " " U64FORMAT "  " U64FORMAT " " U64FORMAT "  " U64FORMAT " " U64FORMAT "\n", lenRaw.nBases,  lenNoOlaps.nBases,  lenEvidence.nBases,  lenCorrOrig.nBases,  lenCorr.nBases,  lenRescOrig.nBases,  lenResc.nBases,  lenNoCoOrig.nBases,  lenNoCo.nBases);
  fprintf(F, "Coverage             " FLTFORMAT " " FLTFORMAT " " FLTFORMAT "  " FLTFORMAT " " FLTFORMAT "  " FLTFORMAT " " FLTFORMAT "  " FLTFORMAT " " FLTFORMAT "\n", lenRaw.cov,     lenNoOlaps.cov,     lenEvidence.cov,     lenCorrOrig.cov,     lenCorr.cov,     lenRescOrig.cov,     lenResc.cov,     lenNoCoOrig.cov,     lenNoCo.coverage);
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
  fprintf(F, "Coverage             " FLTFORMAT " " FLTFORMAT "\n", lenRaw.cov,     lenNoOlaps.cov);
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
  fprintf(F, "Coverage             " FLTFORMAT "  " FLTFORMAT " " FLTFORMAT "  " FLTFORMAT " " FLTFORMAT "\n", lenEvidence.cov,     lenCorrOrig.cov,     lenCorr.cov,     lenRescOrig.cov,     lenResc.cov);
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
  fprintf(F, "Coverage             " FLTFORMAT " " FLTFORMAT "\n", lenNoCoOrig.cov,     lenNoCo.cov);
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
    fprintf(F, "%-10" F_U32P " " U32FORMAT " " U32FORMAT " " U32FORMAT " " U64FORMAT " %c %c %c\n",
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
  char           *seqStoreName      = NULL;
  char           *corStoreName      = NULL;
  char           *outName           = NULL;

  uint32          minOutputCoverage = 4;
  uint32          minOutputLength   = 500;

  uint64          genomeSize        = 0;
  uint32          outCoverage       = 40;

  argc = AS_configure(argc, argv, 1);

  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      corStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-R") == 0) {
      outName = argv[++arg];

    } else if (strcmp(argv[arg], "-cc") == 0) {
      minOutputCoverage = strtoul(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-cl") == 0) {
      minOutputLength = strtoul(argv[++arg], NULL, 10);

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

  if (seqStoreName == NULL)
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
    fprintf(stderr, "  -S seqStore              input reads\n");
    fprintf(stderr, "  -C corStore              input correction layouts\n");
    fprintf(stderr, "  -R asm.readsToCorrect    output ascii list of read IDs to correct\n");
    fprintf(stderr, "                           also creates\n");
    fprintf(stderr, "                             asm.readsToCorrect.stats and\n");
    fprintf(stderr, "                             asm.readsToCorrect.log\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "FILTERING STRATEGY and PARAMETERS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -cc                      minimum coverage of evidence reads\n");
    fprintf(stderr, "  -cl                      minimum length of a corrected read\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -g                       estimated genome size\n");
    fprintf(stderr, "  -c                       desired coverage in corrected reads\n");
    fprintf(stderr, "\n");

    if (seqStoreName == NULL)
      fprintf(stderr, "ERROR: no sequence store (-S) supplied.\n");
    if (corStoreName == NULL)
      fprintf(stderr, "ERROR: no corStore store (-C) supplied.\n");
    if (outName == NULL)
      fprintf(stderr, "ERROR: no output (-R) supplied.\n");

    exit(1);
  }

  sqRead_setDefaultVersion(sqRead_raw);

  sqStore          *seqStore = new sqStore(seqStoreName);
  tgStore          *corStore = new tgStore(corStoreName, 1);

  uint32            numReads = seqStore->sqStore_lastReadID();

  readStatus       *status   = new readStatus [numReads + 1];

  for (uint32 rr=0; rr<numReads+1; rr++)    //  Initialize read status with the readID
    status[rr].readID = rr;                 //  so sort byReadID works correctly.

  //  Scan the tigs, computing expected corrected length.

  {
    falconConsensus  fc(minOutputCoverage, minOutputLength, 0, 0);  //  For memory estimtes

    for (uint32 ti=1; ti<corStore->numTigs(); ti++) {
      tgTig  *layout = corStore->loadTig(ti);

      if (layout) {
        assert(status[ti].readID == layout->tigID());
        status[ti].readID         = layout->tigID();

        status[ti].numOlaps       = layout->numberOfChildren();

        status[ti].origLength     = layout->length();
        status[ti].corrLength     = 0;

        status[ti].memoryRequired = 0;

        fc.analyzeLength(layout,
                         status[ti].corrLength,         //  output
                         status[ti].memoryRequired);    //  output
      }

      corStore->unloadTig(ti);
    }
  }

  //  Sort by expected corrected length, then mark reads for correction until
  //  we get the desired outCoverage.  The unused zeroth read has max
  //  corrected length, so remains first in sorted list and we start at
  //  readID 1.  When done, sort back to readID order.

  auto  byCorLength = [](readStatus const &a, readStatus const &b)   { return(a.corrLength > b.corrLength); };
  auto  byReadID    = [](readStatus const &a, readStatus const &b)   { return(a.readID     < b.readID);     };

  uint64   desiredLength = genomeSize * outCoverage;
  uint64   corrLength    = 0;

  std::sort(status, status + numReads+1, byCorLength);

  for (uint32 rr=1; rr<numReads+1; rr++) {
    if (status[rr].corrLength == 0)
      continue;

    status[rr].usedForCorrection = (corrLength < desiredLength);

    corrLength += status[rr].corrLength;
  }

  std::sort(status, status + numReads+1, byReadID);

  //  Scan the tigs again, this time marking reads used as evidence in the
  //  corrected reads.

  for (uint32 ti=1; ti<corStore->numTigs(); ti++) {
    tgTig  *layout = corStore->loadTig(ti);

    if ((layout != nullptr) &&
        (status[layout->tigID()].usedForCorrection == true))
      for (uint32 ii=0; ii<layout->numberOfChildren(); ii++)
        status[layout->getChild(ii)->ident()].usedForEvidence = true;

    corStore->unloadTig(ti);
  }

  //  And finally, flag any read for correction if it isn't already used as
  //  evidence or being corrected.

  for (uint32 rr=1; rr<numReads+1; rr++)
    if ((status[rr].usedForCorrection == false) &&
        (status[rr].usedForEvidence   == false) &&
        (status[rr].corrLength         > minOutputLength))
      status[rr].rescued = true;

  //  Output the list of reads to correct.

  FILE  *roc      = merylutil::openOutputFile(outName);

  for (uint32 rr=1; rr<numReads+1; rr++)
    if ((status[rr].usedForCorrection == true) ||
        (status[rr].rescued           == true))
      fprintf(roc, "%u\n", status[rr].readID);

  merylutil::closeFile(roc);

  //  Write some statistics and logs.  These are opened just before they're
  //  written so we can use them for completion tests.

  FILE  *stats    = merylutil::openOutputFile(outName, '.', "stats");
  FILE  *log      = merylutil::openOutputFile(outName, '.', "log");

  dumpStatistics(stats, status, numReads, genomeSize);
  dumpLog(log, status, numReads);

  merylutil::closeFile(stats);
  merylutil::closeFile(log);

  //  And say goodbye.

  delete [] status;
  delete    corStore;
  delete    seqStore;

  fprintf(stderr, "Bye.\n");

  exit(0);
}

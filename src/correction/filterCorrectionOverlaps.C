
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
#include "ovStore.H"
#include "strings.H"

#include "computeGlobalScore.H"

#include <vector>
#include <algorithm>

using namespace std;




FILE *
openOutput(char *fileName, bool doOpen) {
  errno = 0;

  if ((doOpen == false) || (fileName == NULL))
    return(NULL);

  FILE *F = fopen(fileName, "w");

  if (errno)
    fprintf(stderr, "ERROR: failed to open '%s' for writing: %s\n", fileName, strerror(errno)), exit(1);

  return(F);
}



int
main(int argc, char **argv) {
  char           *seqStoreName     = NULL;
  char           *ovlStoreName     = NULL;
  char           *scoreFileName    = NULL;
  char            logFileName[FILENAME_MAX];
  char            statsFileName[FILENAME_MAX];

  bool            doEstimate       = true;
  bool            doExact          = false;
  bool            doCompare        = false;

  bool            noLog            = false;
  bool            noStats          = false;

  uint32          expectedCoverage = 25;

  uint32          minOvlLength     = 500;
  uint32          maxOvlLength     = AS_MAX_READLEN;

  double          maxErate         = 1.0;
  double          minErate         = 1.0;

  argc = AS_configure(argc, argv);

  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-scores") == 0) {
      scoreFileName = argv[++arg];


    } else if (strcmp(argv[arg], "-estimate") == 0) {
      doEstimate = true;
      doExact    = false;

    } else if (strcmp(argv[arg], "-exact") == 0) {
      doEstimate = false;
      doExact    = true;

    } else if (strcmp(argv[arg], "-compare") == 0) {
      doEstimate = true;
      doExact    = true;
      doCompare  = true;


    } else if (strcmp(argv[arg], "-c") == 0) {
      expectedCoverage = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-l") == 0) {
      minOvlLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      decodeRange(argv[++arg], minErate, maxErate);


    } else if (strcmp(argv[arg], "-nolog") == 0) {
      noLog = true;

    } else if (strcmp(argv[arg], "-nostats") == 0) {
      noStats = true;

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if (seqStoreName == NULL)
    err++;
  if (ovlStoreName == NULL)
    err++;
  if (scoreFileName == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [options]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Rewrites an ovlStore, filtering overlaps that shouldn't be used for correcting reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S seqStore     input reads\n");
    fprintf(stderr, "  -O ovlStore     input overlaps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -scores sf      output scores for each read, binary file, to file 'sf'\n");
    fprintf(stderr, "                  per-read logging to 'sf.log' (see -nolog)\n");
    fprintf(stderr, "                  summary statistics to 'sf.stats' (see -nostats)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -estimate       estimate the cutoff from precomputed scores\n");
    fprintf(stderr, "  -exact          compute an exact cutoff by reading all overlaps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -compare        output a comparison of estimated vs exact scores\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c coverage     retain at most this many overlaps per read\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -l length       filter overlaps shorter than this length\n");
    fprintf(stderr, "  -e (min-)max    filter overlaps outside this range of fraction error\n");
    fprintf(stderr, "                    example:  -e 0.20          filter overlaps above 20%% error\n");
    fprintf(stderr, "                    example:  -e 0.05-0.20     filter overlaps below 5%% error\n");
    fprintf(stderr, "                                                            or above 20%% error\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Length and Fraction Error filtering NOT SUPPORTED with -estimate.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nolog          don't create 'scoreFile.log'\n");
    fprintf(stderr, "  -nostats        don't create 'scoreFile.stats'\n");

    if (seqStoreName == NULL)
      fprintf(stderr, "ERROR: no sequence store (-S) supplied.\n");
    if (ovlStoreName == NULL)
      fprintf(stderr, "ERROR: no overlap store (-O) supplied.\n");
    if (scoreFileName == NULL)
      fprintf(stderr, "ERROR: no output scoreFile (-scores) supplied.\n");

    exit(1);
  }

  //  If only one value supplied to -e, both minErate and maxErate are set to it.
  //  Change the max to the value supplied, and set min to zero.  Otherwise, two values
  //  were supplied, and we're all set.  No checking is made that the user isn't an idiot.

  if (minErate == maxErate) {
    maxErate = minErate;
    minErate = 0.0;
  }

  sqRead_setDefaultVersion(sqRead_raw);

  sqStore           *seqStore    = new sqStore(seqStoreName);

  ovStore           *ovlStore    = new ovStore(ovlStoreName, seqStore);
  ovStoreHistogram  *ovlHisto    = ovlStore->getHistogram();

  uint32             *numOlaps   = ovlStore->numOverlapsPerRead();

  uint32              ovlLen     = 0;
  uint32              ovlMax     = 0;
  ovOverlap          *ovl        = NULL;

  uint16             *scores     = new uint16 [seqStore->sqStore_lastReadID() + 1];
  uint16              scoreExact = 0;
  uint16              scoreEstim = 0;

  snprintf(logFileName,   FILENAME_MAX, "%s.log",   scoreFileName);
  snprintf(statsFileName, FILENAME_MAX, "%s.stats", scoreFileName);

  FILE               *scoreFile = openOutput(scoreFileName, true);
  FILE               *logFile   = openOutput(logFileName,   (noLog == false));

  globalScore         *gs       = new globalScore(minOvlLength, maxOvlLength, minErate, maxErate, logFile, (noStats == false));

  uint64              readsNoOlaps = 0;

  if (doCompare) {
    fprintf(stdout, "  readID  exact  estim\n");
    //fprintf(stdout, "-------- ------ ------\n");
  }

  for (uint32 id=0; id <= seqStore->sqStore_lastReadID(); id++) {
    scores[id] = UINT16_MAX;

    if (numOlaps[id] == 0) {
      readsNoOlaps++;
      continue;
    }

    if (doEstimate == true) {
      scores[id] = scoreEstim = ovlHisto->overlapScoreEstimate(id, expectedCoverage);

      gs->estimate(numOlaps[id], expectedCoverage);     //  Just for stats collection
    }

    if (doExact == true) {
      ovlLen = ovlStore->loadOverlapsForRead(id, ovl, ovlMax);

      if (ovlLen > 0) {
        assert(ovlLen == numOlaps[id]);
        assert(ovl[0].a_iid == id);

        scores[id] = scoreExact = gs->compute(ovlLen, ovl, expectedCoverage, 0, NULL);
      }
    }

    if (doCompare) {
      fprintf(stdout, "%8u %6u %6u\n", id, scoreExact, scoreEstim);
    }
  }

  if (scoreFile)
    writeToFile(scores, "scores", seqStore->sqStore_lastReadID() + 1, scoreFile);

  AS_UTL_closeFile(scoreFile, scoreFileName);
  AS_UTL_closeFile(logFile,   logFileName);

  delete [] scores;

  delete [] ovl;
  delete [] numOlaps;
  delete    ovlHisto;
  delete    ovlStore;

  delete seqStore;

  if (noStats == true)
    exit(0);

  FILE  *statsFile = AS_UTL_openOutputFile(statsFileName);

  fprintf(statsFile, "PARAMETERS:\n");
  fprintf(statsFile, "----------\n");
  fprintf(statsFile, "\n");
  fprintf(statsFile, "%7" F_U32P " (expected coverage)\n", expectedCoverage);
  fprintf(statsFile, "%7" F_U32P " (don't use overlaps shorter than this)\n", minOvlLength);
  fprintf(statsFile, "%7.3f (don't use overlaps with erate less than this)\n", minErate);
  fprintf(statsFile, "%7.3f (don't use overlaps with erate more than this)\n", maxErate);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "OVERLAPS:\n");
  fprintf(statsFile, "--------\n");
  fprintf(statsFile, "\n");
  fprintf(statsFile, "IGNORED:\n");
  fprintf(statsFile, "\n");
  fprintf(statsFile, "%12" F_U64P " (< %6.4f fraction error)\n", gs->lowErate(),  minErate);
  fprintf(statsFile, "%12" F_U64P " (> %6.4f fraction error)\n", gs->highErate(), maxErate);
  fprintf(statsFile, "%12" F_U64P " (< %u bases long)\n", gs->tooShort(),  minOvlLength);
  fprintf(statsFile, "%12" F_U64P " (> %u bases long)\n", gs->tooLong(),   maxOvlLength);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "FILTERED:%s\n", (doEstimate == true) ? " (estimated)" : "");
  fprintf(statsFile, "\n");
  fprintf(statsFile, "%12" F_U64P " (too many overlaps, discard these shortest ones)\n", gs->belowCutoff());
  fprintf(statsFile, "\n");
  fprintf(statsFile, "EVIDENCE:%s\n", (doEstimate == true) ? " (estimated)" : "");
  fprintf(statsFile, "\n");
  fprintf(statsFile, "%12" F_U64P " (longest overlaps)\n", gs->retained());
  fprintf(statsFile, "\n");
  fprintf(statsFile, "TOTAL:\n");
  fprintf(statsFile, "\n");
  fprintf(statsFile, "%12" F_U64P " (all overlaps)\n", gs->totalOverlaps());
  fprintf(statsFile, "\n");
  fprintf(statsFile, "READS:%s\n", (doEstimate == true) ? " (estimated)" : "");
  fprintf(statsFile, "-----\n");
  fprintf(statsFile, "\n");
  fprintf(statsFile, "%12" F_U64P " (no overlaps)\n", readsNoOlaps);
  fprintf(statsFile, "%12" F_U64P " (no overlaps filtered)\n",      gs->reads00OlapsFiltered());
  fprintf(statsFile, "%12" F_U64P " (<=  50%% overlaps filtered)\n", gs->reads50OlapsFiltered());
  fprintf(statsFile, "%12" F_U64P " (<=  80%% overlaps filtered)\n", gs->reads80OlapsFiltered());
  fprintf(statsFile, "%12" F_U64P " (<=  95%% overlaps filtered)\n", gs->reads95OlapsFiltered());
  fprintf(statsFile, "%12" F_U64P " (<= 100%% overlaps filtered)\n", gs->reads99OlapsFiltered());
  fprintf(statsFile, "\n");

  AS_UTL_closeFile(statsFile, statsFileName);

  delete gs;

  //  Histogram of overlaps per read
  //  Histogram of overlaps filtered per read
  //  Histogram of overlaps used for correction (number per read, lengths?)

  exit(0);
}

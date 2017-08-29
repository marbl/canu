
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
 *    Brian P. Walenz from 2015-MAY-28 to 2015-JUN-25
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-24
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "ovStore.H"
#include "splitToWords.H"

#include "AS_UTL_decodeRange.H"

#include "computeGlobalScore.H"

#include <vector>
#include <algorithm>

using namespace std;

int
main(int argc, char **argv) {
  char           *gkpStoreName     = NULL;
  char           *ovlStoreName     = NULL;
  char           *scoreFileName    = NULL;
  char            logFileName[FILENAME_MAX];
  char            statsFileName[FILENAME_MAX];

  bool            noLog            = false;
  bool            noStats          = false;

  uint32          expectedCoverage = 25;

  uint32          minOvlLength     = 500;
  uint32          maxOvlLength     = AS_MAX_READLEN;

  double          maxErate         = 1.0;
  double          minErate         = 1.0;

  bool		  legacyScore	   = false;

  argc = AS_configure(argc, argv);

  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-S") == 0) {
      scoreFileName = argv[++arg];


    } else if (strcmp(argv[arg], "-c") == 0) {
      expectedCoverage = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-l") == 0) {
      minOvlLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      AS_UTL_decodeRange(argv[++arg], minErate, maxErate);


    } else if (strcmp(argv[arg], "-nolog") == 0) {
      noLog = true;

    } else if (strcmp(argv[arg], "-nostats") == 0) {
      noStats = true;

    } else if (strcmp(argv[arg], "-legacy") == 0) {
      legacyScore = true;

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if (gkpStoreName == NULL)
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
    fprintf(stderr, "  -G gkpStore     input reads\n");
    fprintf(stderr, "  -O ovlStore     input overlaps\n");
    fprintf(stderr, "  -S scoreFile    output scores for each read, binary file, to 'scoreFile'\n");
    fprintf(stderr, "                  per-read logging to 'scoreFile.log' (see -nolog)\n");
    fprintf(stderr, "                  summary statistics to 'scoreFile.stats' (see -nostats)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c coverage     retain at most this many overlaps per read\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -l length       filter overlaps shorter than this length\n");
    fprintf(stderr, "  -e (min-)max    filter overlaps outside this range of fraction error\n");
    fprintf(stderr, "                    example:  -e 0.20          filter overlaps above 20%% error\n");
    fprintf(stderr, "                    example:  -e 0.05-0.20     filter overlaps below 5%% error\n");
    fprintf(stderr, "                                                            or above 20%% error\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nolog          don't create 'scoreFile.log'\n");
    fprintf(stderr, "  -nostats        don't create 'scoreFile.stats'\n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gatekeeper store (-G) supplied.\n");
    if (ovlStoreName == NULL)
      fprintf(stderr, "ERROR: no overlap store (-O) supplied.\n");
    if (scoreFileName == NULL)
      fprintf(stderr, "ERROR: no output scoreFile (-S) supplied.\n");

    exit(1);
  }

  //  If only one value supplied to -e, both minErate and maxErate are set to it.
  //  Change the max to the value supplied, and set min to zero.  Otherwise, two values
  //  were supplied, and we're all set.  No checking is made that the user isn't an idiot.

  if (minErate == maxErate) {
    maxErate = minErate;
    minErate = 0.0;
  }

  uint32    maxEvalue = AS_OVS_encodeEvalue(maxErate);
  uint32    minEvalue = AS_OVS_encodeEvalue(minErate);;

  gkStore  *gkpStore  = gkStore::gkStore_open(gkpStoreName);

  ovStore  *inpStore  = new ovStore(ovlStoreName, gkpStore);

  uint64   *scores    = new uint64 [gkpStore->gkStore_getNumReads() + 1];


  snprintf(logFileName, FILENAME_MAX, "%s.log", scoreFileName);
  snprintf(statsFileName, FILENAME_MAX, "%s.stats", scoreFileName);

  errno = 0;
  FILE     *scoreFile   = (scoreFileName == NULL) ? NULL : fopen(scoreFileName, "w");
  if (errno)
    fprintf(stderr, "ERROR: failed to open '%s' for writing: %s\n", scoreFileName, strerror(errno)), exit(1);


  errno = 0;
  FILE     *logFile     = (noLog == true) ? NULL : fopen(logFileName, "w");
  if (errno)
    fprintf(stderr, "ERROR: failed to open '%s' for writing: %s\n", logFileName, strerror(errno)), exit(1);


  uint32       ovlLen = 0;
  uint32       ovlMax = 131072;
  ovOverlap   *ovl    = ovOverlap::allocateOverlaps(gkpStore, ovlMax);

  uint64      readsNoOlaps          = 0;

  globalScore *gs     = new globalScore(minOvlLength, maxOvlLength, minErate, maxErate, logFile, (noStats == false));

  for (uint32 id=1; id <= gkpStore->gkStore_getNumReads(); id++) {
    scores[id] = UINT64_MAX;

    inpStore->readOverlaps(id, ovl, ovlLen, ovlMax);

    if ((ovlLen       == 0) ||
        (ovl[0].a_iid != id)) {
      readsNoOlaps++;
      continue;
    }

    scores[id] = gs->compute(ovlLen, ovl, expectedCoverage, 0, NULL);
  }

  if (scoreFile) {
    AS_UTL_safeWrite(scoreFile, scores, "scores", sizeof(uint64), gkpStore->gkStore_getNumReads() + 1);
    fclose(scoreFile);
  }

  if (logFile)
    fclose(logFile);

  delete [] scores;

  delete [] ovl;
  delete    inpStore;

  gkpStore->gkStore_close();

  if (noStats == true)
    exit(0);

  errno = 0;
  FILE     *statsFile = fopen(statsFileName, "w");
  if (errno)
    fprintf(stderr, "ERROR: failed to open '%s' for writing: %s\n", statsFileName, strerror(errno)), exit(1);

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
  fprintf(statsFile, "%12" F_U64P " (< %6.4f fraction error)\n", gs->lowErate(),  AS_OVS_decodeEvalue(minEvalue));
  fprintf(statsFile, "%12" F_U64P " (> %6.4f fraction error)\n", gs->highErate(), AS_OVS_decodeEvalue(maxEvalue));
  fprintf(statsFile, "%12" F_U64P " (< %u bases long)\n", gs->tooShort(),  minOvlLength);
  fprintf(statsFile, "%12" F_U64P " (> %u bases long)\n", gs->tooLong(),   maxOvlLength);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "FILTERED:\n");
  fprintf(statsFile, "\n");
  fprintf(statsFile, "%12" F_U64P " (too many overlaps, discard these shortest ones)\n", gs->belowCutoff());
  fprintf(statsFile, "\n");
  fprintf(statsFile, "EVIDENCE:\n");
  fprintf(statsFile, "\n");
  fprintf(statsFile, "%12" F_U64P " (longest overlaps)\n", gs->retained());
  fprintf(statsFile, "\n");
  fprintf(statsFile, "TOTAL:\n");
  fprintf(statsFile, "\n");
  fprintf(statsFile, "%12" F_U64P " (all overlaps)\n", gs->totalOverlaps());
  fprintf(statsFile, "\n");
  fprintf(statsFile, "READS:\n");
  fprintf(statsFile, "-----\n");
  fprintf(statsFile, "\n");
  fprintf(statsFile, "%12" F_U64P " (no overlaps)\n", readsNoOlaps);
  fprintf(statsFile, "%12" F_U64P " (no overlaps filtered)\n",      gs->reads00OlapsFiltered());
  fprintf(statsFile, "%12" F_U64P " (<  50%% overlaps filtered)\n", gs->reads50OlapsFiltered());
  fprintf(statsFile, "%12" F_U64P " (<  80%% overlaps filtered)\n", gs->reads80OlapsFiltered());
  fprintf(statsFile, "%12" F_U64P " (<  95%% overlaps filtered)\n", gs->reads95OlapsFiltered());
  fprintf(statsFile, "%12" F_U64P " (< 100%% overlaps filtered)\n", gs->reads99OlapsFiltered());
  fprintf(statsFile, "\n");
  fclose(statsFile);

  delete gs;

  //  Histogram of overlaps per read
  //  Histogram of overlaps filtered per read
  //  Histogram of overlaps used for correction (number per read, lengths?)

  exit(0);
}

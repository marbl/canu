


#include "AS_global.H"
#include "ovStore.H"
#include "splitToWords.H"

#include "AS_UTL_decodeRange.H"

#include <vector>
#include <algorithm>

using namespace std;


int
main(int argc, char **argv) {
  char           *gkpStoreName     = NULL;
  char           *ovlStoreName     = NULL;
  char           *scoreFileName    = NULL;
  char           *logFileName      = NULL;

  uint32          expectedCoverage = 25;

  uint32          minOvlLength     = 500;
  uint32          maxOvlLength     = UINT32_MAX;

  double          maxErate         = DBL_MAX;
  double          minErate         = DBL_MAX;

  bool            noContain        = false;
  bool            noDovetail       = false;

  uint32          maxHang          = UINT32_MAX;

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


    } else if (strcmp(argv[arg], "-logfile") == 0) {
      logFileName = argv[++arg];


    } else if (strcmp(argv[arg], "-nocontain") == 0) {
      noContain = true;

    } else if (strcmp(argv[arg], "-nodovetail") == 0) {
      noDovetail = true;


    } else if (strcmp(argv[arg], "-maxHang") == 0) {
      maxHang = atoi(argv[++arg]);


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
    fprintf(stderr, "  -S scoreFile    output scores for each read, binary file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c coverage     retain at most this many overlaps per read\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -l length       filter overlaps shorter than this length\n");
    fprintf(stderr, "  -e (min-)max    filter overlaps outside this range of fraction error\n");
    fprintf(stderr, "                    example:  -e 0.20          filter overlaps above 20%% error\n");
    fprintf(stderr, "                    example:  -e 0.05-0.20     filter overlaps below 5%% error\n");
    fprintf(stderr, "                                                            or above 20%% error\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -logfile L      write detailed per-read logging to file L\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The following are not implemented:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nocontain      filter overlaps that are contained in the target read\n");
    fprintf(stderr, "  -nodovetail     filter dovetail overlaps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -maxhang h      filter overlaps with more than 'h' bases unaligned on both the\n");
    fprintf(stderr, "                  target and evidence read, on either end\n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gatekeeper store (-G) supplied.\n");
    if (ovlStoreName == NULL)
      fprintf(stderr, "ERROR: no overlap store (-O) supplied.\n");
    if (scoreFileName == NULL)
      fprintf(stderr, "ERROR: no log file name (-f) supplied.\n");

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

  gkStore  *gkpStore  = new gkStore(gkpStoreName);

  ovStore  *inpStore  = new ovStore(ovlStoreName);

  uint32   *scores    = new uint32 [gkpStore->gkStore_getNumReads() + 1];


  errno = 0;
  FILE     *scoreFile   = (scoreFileName == NULL) ? NULL : fopen(scoreFileName, "w");
  if (errno)
    fprintf(stderr, "ERROR: failed to open '%s' for writing: %s\n", scoreFileName, strerror(errno)), exit(1);


  errno = 0;
  FILE     *logFile     = (logFileName == NULL) ? NULL : fopen(logFileName, "w");
  if (errno)
    fprintf(stderr, "ERROR: failed to open '%s' for writing: %s\n", logFileName, strerror(errno)), exit(1);


  uint32      ovlLen = 0;
  uint32      ovlMax = 131072;
  ovsOverlap *ovl    = new ovsOverlap [ovlMax];
  ovsOverlap  swapped;

  uint32      histLen = 0;
  uint32      histMax = ovlMax;
  uint32     *hist    = new uint32 [histMax];

  uint64      totalOverlaps = 0;
  uint64      lowErate     = 0;
  uint64      highErate    = 0;
  uint64      tooShort     = 0;
  uint64      tooLong      = 0;
  uint64      isContain    = 0;
  uint64      isDovetail   = 0;
  uint64      belowCutoff  = 0;
  uint64      retained     = 0;

  for (uint32 id=1; id <= gkpStore->gkStore_getNumReads(); id++) {
    scores[id] = UINT32_MAX;

    inpStore->readOverlaps(id, ovl, ovlLen, ovlMax);

    if (ovlLen == 0)
      continue;

    if (ovl[0].a_iid != id)
      continue;

    histLen = 0;

    if (histMax < ovlMax) {
      delete [] hist;

      histMax = ovlMax;
      hist    = new uint32 [ovlMax];
    }

    //  Figure out which overlaps are good enough to consider and save their length.

    for (uint32 oo=0; oo<ovlLen; oo++) {
      uint32  ovlLength  = ovl[oo].a_end(gkpStore) - ovl[oo].a_bgn(gkpStore);
      uint32  ovlScore   = 100 * ovlLength * (1 - ovl[oo].erate());
      bool    isC        = false;
      bool    isD        = false;

      if ((ovl[oo].evalue() < minEvalue)        ||
          (maxEvalue        < ovl[oo].evalue()) ||
          (ovlLength        < minOvlLength)     ||
          (maxOvlLength     < ovlLength)        ||
          ((isC == true) && (noContain  == true)) ||
          ((isD == true) && (noDovetail == true)))
        continue;

      hist[histLen++] = ovlScore;
    }

    //  Sort the lengths of overlaps we would save.

    sort(hist, hist + histLen);

    //  Figure out our threshold score.  Any overlap with score below this should be filtered.

    if (expectedCoverage <= histLen)
      scores[id] = hist[histLen - expectedCoverage];
    else
      scores[id] = 0;

    //  One more pass, just to gather statistics

    uint32 belowCutoffLocal = 0;

    for (uint32 oo=0; oo<ovlLen; oo++) {
      uint32  ovlLength  = ovl[oo].a_end(gkpStore) - ovl[oo].a_bgn(gkpStore);
      uint32  ovlScore   = 100 * ovlLength * (1 - ovl[oo].erate());
      bool    isC        = false;
      bool    isD        = false;
      bool    skipIt     = false;

      totalOverlaps++;

      //  First, count the filtering done above.

      if (ovl[oo].evalue() < minEvalue) {
        lowErate++;
        skipIt = true;
      }

      if (maxEvalue < ovl[oo].evalue()) {
        highErate++;
        skipIt = true;
      }

      if (ovlLength < minOvlLength) {
        tooShort++;
        skipIt = true;
      }

      if (maxOvlLength < ovlLength) {
        tooLong++;
        skipIt = true;
      }

      if ((isC == true) && (noContain  == true)) {
        isContain++;
        skipIt = true;
      }

      if ((isD == true) && (noDovetail == true)) {
        isDovetail++;
        skipIt = true;
      }

      //  Now, apply the global filter cutoff, only if the overlap wasn't already tossed out.

      if ((skipIt == false) &&
          (ovlScore < scores[id])) {
        belowCutoff++;
        belowCutoffLocal++;
        skipIt = true;
      }

      if (skipIt)
        continue;

      retained++;
    }

    if (logFile) {
      if (histLen <= expectedCoverage)
        fprintf(logFile, "%9u - %6u overlaps - %6u scored - %6u filtered - %4u saved (no filtering)\n",
                id, ovlLen, histLen, 0, histLen);
      else
        fprintf(logFile, "%9u - %6u overlaps - %6u scored - %6u filtered - %4u saved (length * erate cutoff %.2f)\n",
                id, ovlLen, histLen, belowCutoffLocal, histLen - belowCutoffLocal, scores[id] / 100.0);
    }
  }

  if (scoreFile)
    AS_UTL_safeWrite(scoreFile, scores, "scores", sizeof(uint32), gkpStore->gkStore_getNumReads() + 1);

  if (scoreFile)
    fclose(scoreFile);

  if (logFile)
    fclose(logFile);

  delete [] scores;

  delete inpStore;
  delete gkpStore;



  fprintf(stderr, "\n");
  fprintf(stderr, "Processed  %lu overlaps.\n", totalOverlaps);
  fprintf(stderr, "\n");
  fprintf(stderr, "lowErate   %lu (< %6.4f fraction error)\n", lowErate,  AS_OVS_decodeEvalue(minEvalue));
  fprintf(stderr, "highErate  %lu (> %6.4f fraction error)\n", highErate, AS_OVS_decodeEvalue(maxEvalue));
  fprintf(stderr, "tooShort   %lu (< %u bases)\n", tooShort,  minOvlLength);
  fprintf(stderr, "tooLong    %lu (> %u bases)\n", tooLong,   maxOvlLength);
  fprintf(stderr, "\n");
  fprintf(stderr, "isContain  %lu (evidence contained in target)\n", isContain);
  fprintf(stderr, "isDovetail %lu (evidence dovetail to target)\n", isDovetail);
  fprintf(stderr, "\n");
  fprintf(stderr, "filtered   %lu (shortest overlaps, not filtered above)\n", belowCutoff);
  fprintf(stderr, "retained   %lu (longest overlaps)\n",  retained);
  fprintf(stderr, "\n");

  exit(0);
}


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
 *  This file is derived from:
 *
 *    src/correction/filterCorrectionOverlaps.C
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

#include "computeGlobalScore.H"

#include <set>
#include <algorithm>

using namespace std;


uint16
globalScore::compute(uint32             ovlLen,
                     ovOverlap         *ovl,
                     uint32             expectedCoverage,
                     uint32             thresholdsLen,
                     uint16            *thresholds) {

  //  Build a list of all the overlap scores.  Ignore
  //  overlaps that are too bad/good or too short/long.

  histLen = 0;

  resizeArray(hist, histLen, histMax, ovlLen, resizeArray_doNothing);

  for (uint32 oo=0; oo<ovlLen; oo++) {
    uint32  ovlLength  = ovl[oo].a_end() - ovl[oo].a_bgn();

    if ((ovl[oo].evalue() < minEvalue)    || (maxEvalue    < ovl[oo].evalue()) ||
        (ovlLength        < minOvlLength) || (maxOvlLength < ovlLength))
      continue;

    hist[histLen++] = ovl[oo].overlapScore();
  }

  //  Sort, reversely.

#ifdef _GLIBCXX_PARALLEL
  __gnu_sequential::
#endif
  sort(hist, hist + histLen, std::greater<uint64>());

  //  Figure out our threshold score.  Any overlap with score below this should be filtered.

  uint16 threshold = (expectedCoverage < histLen) ? hist[expectedCoverage] : 0;

  for (uint32 ii=0; ii<thresholdsLen; ii++) {
    if (ii * 10 < histLen)
      thresholds[ii] = hist[ii * 10];
    else
      thresholds[ii] = 0;
  }

  if (stats == NULL)
    return(threshold);

  //  One more pass, just to gather statistics now that we know the threshold score.

  uint32 belowCutoffLocal = 0;

  for (uint32 oo=0; oo<ovlLen; oo++) {
    uint64  ovlLength  = ovl[oo].a_len();
    uint16  ovlScore   = ovl[oo].overlapScore();

    bool    isC        = false;
    bool    isD        = false;
    bool    skipIt     = false;

    stats->totalOverlaps++;

    //  First, count the filtering done above.

    if (ovl[oo].evalue() < minEvalue) {
      stats->lowErate++;
      skipIt = true;
    }

    if (maxEvalue < ovl[oo].evalue()) {
      stats->highErate++;
      skipIt = true;
    }

    if (ovlLength < minOvlLength) {
      stats->tooShort++;
      skipIt = true;
    }

    if (maxOvlLength < ovlLength) {
      stats->tooLong++;
      skipIt = true;
    }

    //  Now, apply the global filter cutoff, only if the overlap wasn't already tossed out.

    if ((skipIt == false) &&
        (ovlScore < threshold)) {
      stats->belowCutoff++;
      belowCutoffLocal++;
      skipIt = true;
    }

    if (skipIt)
      continue;

    stats->retained++;
  }  //  Over all overlaps

  double  fractionFiltered = (double)belowCutoffLocal / histLen;

  if (fractionFiltered <= 0.00)   stats->reads00OlapsFiltered++;
  if (fractionFiltered <= 0.50)   stats->reads50OlapsFiltered++;
  if (fractionFiltered <= 0.80)   stats->reads80OlapsFiltered++;
  if (fractionFiltered <= 0.95)   stats->reads95OlapsFiltered++;
  if (fractionFiltered <= 1.00)   stats->reads99OlapsFiltered++;

  if (logFile)
    if (histLen <= expectedCoverage)
      fprintf(logFile, "%9u - %6u overlaps - %6u scored - %6u filtered - %4u saved (no filtering)\n",
              ovl[0].a_iid, ovlLen, histLen, 0, histLen);
    else
      fprintf(logFile, "%9u - %6u overlaps - %6u scored - %6u filtered - %4u saved (threshold %u)\n",
              ovl[0].a_iid, ovlLen, histLen, belowCutoffLocal, histLen - belowCutoffLocal, threshold);

  return(threshold);
}



void
globalScore::estimate(uint32            ovlLen,
                      uint32            expectedCoverage) {
  double    fractionFiltered = 0;

  stats->totalOverlaps += ovlLen;

  if (ovlLen < expectedCoverage) {
    stats->retained    += ovlLen;

  } else {
    fractionFiltered = (double)(ovlLen - expectedCoverage) / ovlLen;

    stats->retained    +=          expectedCoverage;
    stats->belowCutoff += ovlLen - expectedCoverage;
  }

  if (fractionFiltered <= 0.00)   stats->reads00OlapsFiltered++;
  if (fractionFiltered <= 0.50)   stats->reads50OlapsFiltered++;
  if (fractionFiltered <= 0.80)   stats->reads80OlapsFiltered++;
  if (fractionFiltered <= 0.95)   stats->reads95OlapsFiltered++;
  if (fractionFiltered <= 1.00)   stats->reads99OlapsFiltered++;
}



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
 *    Brian P. Walenz beginning on 2015-JUL-23
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef STDDEV_H
#define STDDEV_H

#include "AS_global.H"

#include <vector>

using namespace std;



template<typename TT>
void
computeStdDev(vector<TT> dist, TT &mean, TT &median, TT &mode, double &stddev) {
  mean   = 0;
  median = 0;
  mode   = 0;
  stddev = 0;

  if (dist.size() == 0)
    return;

  //  Sort the values.  Lets us approximate the stddev for filtering out outliers,
  //  and makes finding the mode simple.

  sort(dist.begin(), dist.end());

  //  Approximate the stddev to filter out outliers.  This is done by assuming we're normally
  //  distributed, finding the values that would represent 1 standard deviation (about 68.27% of the
  //  data), and using that to find the 5 std.dev. limits.

  median        = dist[dist.size() * 1 / 2];

  TT oneThird   = dist[dist.size() * 1 / 3];
  TT twoThird   = dist[dist.size() * 2 / 3];

  TT approxStd  = max(median - oneThird, twoThird - median);

  TT biggest    = median + approxStd * 5;
  TT smallest   = median - approxStd * 5;

  //  Now, compute the number of samples within our bounts.  And find the mean, too.

  size_t  numSamples = 0;

  for (size_t x=0; x<dist.size(); x++)
    if ((smallest  <= dist[x]) &&
        (dist[x]   <= biggest)) {
      numSamples += 1;
      mean       += dist[x];
    }

  if (numSamples == 0)
    return;

  mean   = mean / numSamples;

  //  Use the standard std.dev. algorithm, tossing out the outliers.

  for (uint64 x=0; x<dist.size(); x++)
    if ((smallest  <= dist[x]) &&
        (dist[x]   <= biggest))
      stddev += (dist[x] - mean) * (dist[x] - mean);

  if (numSamples > 1)
    stddev = sqrt(stddev / (numSamples - 1));

  //  And the mode.  The values are sorted, so we just need to scan the list and remember the
  //  most common value.

  uint32  modeCnt = 0;
  TT      modeVal = 0;

  uint32  modeTmpCnt = 0;
  TT      modeTmpVal = 0;

  for (uint64 x=0; x<dist.size(); x++) {
    if (dist[x] != modeTmpVal) {
      if (modeCnt < modeTmpCnt) {
        modeCnt = modeTmpCnt;
        modeVal = modeTmpVal;
      }

      modeTmpCnt = 1;
      modeTmpVal = dist[x];
    }

    modeTmpCnt++;
  }

  if (modeCnt < modeTmpCnt) {
    modeCnt = modeTmpCnt;
    modeVal = modeTmpVal;
  }

  mode = modeVal;
};



template<typename TT>
TT
computeExponentialMovingAverage(TT alpha, TT ema, TT value) {
  assert(0.0   <= alpha);
  assert(alpha <= 1.0);

  return(alpha * value + (1 - alpha) * ema);
};



#endif  //  STDDEV_H

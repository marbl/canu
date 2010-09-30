
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static const char *rcsid = "$Id: AS_BOG_InsertSizes.cc,v 1.2 2010-09-30 05:50:17 brianwalenz Exp $";

#include "AS_BOG_InsertSizes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_Unitig.hh"


void
InsertSizes::accumulateLibraryStats(Unitig *utg) {

  for (uint32 fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
    DoveTailNode  *frag = &(*utg->dovetail_path_ptr)[fi];

    if (FI->mateIID(frag->ident) == 0)
      //  Unmated fragment.
      continue;

    if (utg->id() != utg->fragIn(FI->mateIID(frag->ident)))
      //  Mate in a different unitig.
      continue;

    uint32        mi   = utg->pathPosition(FI->mateIID(frag->ident));
    DoveTailNode *mate = &(*utg->dovetail_path_ptr)[mi];

    if (frag->ident < mate->ident)
      //  Only do this once for each mate pair.
      continue;

#warning assumes innie mate pairs
    if (isReverse(frag->position) == isReverse(mate->position))
      //  Same orient mates, not a good mate pair.
      continue;

    //  Compute the insert size.

    int32  insertSize = 0;

    if (isReverse(frag->position)) {
      if (frag->position.end <= mate->position.bgn)
        //  Misordered, not a good mate pair.  NOTE: this is specifically allowing mates that
        //  overlap, i.e., insert size is less than the sum of fragment lengths.
        continue;

      //  We must have, at least, the following picture, where the relationship on the left is
      //  strict.  'frag' is allowed to move to the right, and 'mate' is allowed to move to the
      //  left, but moving either in the other direction turns this into an outtie relationship.
      //
      //  (end)  <----- (bgn)  frag
      //  (bgn) ----->  (end) mate
      //
      assert(mate->position.bgn < frag->position.end);

      //  However, we aren't guaranteed that the right side is ordered properly (i.e., misordered by
      //  one fragment contained in the other fragment).  We hope this doesn't happen, but if it
      //  does, we'll catch it.  The insertSize defaults to zero, which is then ignored below.
      //
      if (mate->position.bgn < frag->position.bgn)
        insertSize = frag->position.bgn - mate->position.bgn;

    } else {
      //  (See comments above)
      if (mate->position.end <= frag->position.bgn)
        continue;

      assert(frag->position.bgn < mate->position.end);

      if (frag->position.bgn < mate->position.bgn)
        insertSize = mate->position.bgn - frag->position.bgn;
    }

    assert(insertSize > 0);

    uint32 di = FI->libraryIID(frag->ident);

    if (_distMax[di] <= _distLen[di]) {
      _distMax[di] *= 2;
      int32 *d = new int32 [_distMax[di]];
      memcpy(d, _dist[di], sizeof(int32) * _distLen[di]);
      delete [] _dist[di];
      _dist[di] = d;
    }

    _dist[di][_distLen[di]++] = insertSize;
  }
}



InsertSizes::InsertSizes() {

  _numLibs        = FI->numLibraries() + 1;

  _dist    = new int32 * [_numLibs + 1];
  _distLen = new int32   [_numLibs + 1];
  _distMax = new int32   [_numLibs + 1];

  _mean           = new int32   [_numLibs + 1];
  _stddev         = new int32   [_numLibs + 1];
  _samples        = new int32   [_numLibs + 1];

  _distLen[0] = 0;
  _distMax[0] = 0;
  _dist[0]    = NULL;

  for (uint32 i=1; i<_numLibs + 1; i++) {
    _distLen[i] = 0;
    _distMax[i] = 1048576;
    _dist[i]    = new int32 [_distMax[i]];

    _mean[i]     = FI->mean(i);
    _stddev[i]   = FI->stddev(i);
    _samples[i]  = FI->numMatesInLib(i);
  }

  for (uint32 ti=0; ti<UG->unitigs.size(); ti++) {
    Unitig        *utg = UG->unitigs[ti];

    if ((utg == NULL) ||
        (utg->dovetail_path_ptr->size() < 2))
      continue;

    accumulateLibraryStats(utg);
  }

  for (uint32 i=1; i<_numLibs + 1; i++)
    sort(_dist[i], _dist[i] + _distLen[i]);

  //  Disregard outliers (those outside 5 (estimated) stddevs) and recompute global stddev

  for (uint32 i=1; i<_numLibs + 1; i++) {
    int32     median     = _dist[i][_distLen[i] * 1 / 2];
    int32     oneThird   = _dist[i][_distLen[i] * 1 / 3];
    int32     twoThird   = _dist[i][_distLen[i] * 2 / 3];

    int32     aproxStd   = MAX(median - oneThird, twoThird - median);

    int32     biggest    = median + aproxStd * 5;
    int32     smallest   = median - aproxStd * 5;

    uint32    numPairs   = 0;
    double    sum_Dists  = 0.0;
    double    sumSquares = 0.0;

    for (uint32 d=0; d<_distLen[i]; d++)
      if ((smallest    <= _dist[i][d]) &&
          (_dist[i][d] <= biggest)) {
        numPairs++;
        sum_Dists += _dist[i][d];
      }

    _samples[i] = numPairs;
    _mean[i]    = (numPairs > 0) ? sum_Dists / numPairs : 0;

    for (uint32 d=0; d<_distLen[i]; d++)
      if ((smallest    <= _dist[i][d]) &&
          (_dist[i][d] <= biggest))
        sumSquares += ((_dist[i][d] - _mean[i]) *
                       (_dist[i][d] - _mean[i]));

    _stddev[i]  = (numPairs > 1) ? sqrt(sumSquares / (numPairs - 1)) : 0.0;
  }

#if 0
  fprintf(logFile, "LIB\tn_Dist\tmedian\t1/3rd\t2/3rd\tmaxDiff\tmin\tmax\tnumGood\tmean\tstddev\n");
  for (uint32 i=1; i<_numLibs + 1; i++)
    fprintf(logFile, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.1f\n",
            i, _distLen[i], median, oneThird, twoThird, aproxStd, smallest, biggest,
            numPairs, _mean[i], _stddev[i]);
#endif

  for (uint32 i=0; i<_numLibs + 1; i++)
    delete [] _dist[i];

  delete [] _dist;
  delete [] _distLen;
  delete [] _distMax;

  delete [] _samples;
};



InsertSizes::~InsertSizes() {
  delete [] _mean;
  delete [] _stddev;
};











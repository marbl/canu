
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
 *    src/AS_BAT/AS_BAT_InsertSizes.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static const char *rcsid = "$Id$";

#include "AS_BAT_InsertSizes.H"
#include "AS_BAT_Unitig.H"

void
InsertSizes::accumulateLibraryStats(Unitig *utg) {

  for (uint32 fi=0; fi<utg->ufpath.size(); fi++) {
    ufNode  *frag = &utg->ufpath[fi];

    if (FI->mateIID(frag->ident) == 0)
      //  Unmated fragment.
      continue;

    if (utg->id() != utg->fragIn(FI->mateIID(frag->ident)))
      //  Mate in a different unitig.
      continue;

    uint32        mi   = utg->pathPosition(FI->mateIID(frag->ident));
    ufNode *mate = &utg->ufpath[mi];

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



InsertSizes::InsertSizes(UnitigVector &unitigs) {

  _numLibs        = FI->numLibraries();

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

    _mean[i]     = (int32)FI->mean(i);
    _stddev[i]   = (int32)FI->stddev(i);
    _samples[i]  = FI->numMatesInLib(i);
  }

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig        *utg = unitigs[ti];

    if ((utg == NULL) ||
        (utg->ufpath.size() < 2))
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

    for (int32 d=0; d<_distLen[i]; d++)
      if ((smallest    <= _dist[i][d]) &&
          (_dist[i][d] <= biggest)) {
        numPairs++;
        sum_Dists += _dist[i][d];
      }

    _samples[i] = numPairs;
    _mean[i]    = (numPairs > 0) ? sum_Dists / numPairs : 0;

    for (int32 d=0; d<_distLen[i]; d++)
      if ((smallest    <= _dist[i][d]) &&
          (_dist[i][d] <= biggest))
        sumSquares += ((double)(_dist[i][d] - _mean[i]) *
                       (double)(_dist[i][d] - _mean[i]));

    _stddev[i]  = (numPairs > 1) ? sqrt(sumSquares / (numPairs - 1)) : 0.0;

    writeLog("InsertSizes()-- lib %d mean %d stddev %d samples %d\n", i, _mean[i], _stddev[i], _samples[i]);
  }

  for (uint32 i=0; i<_numLibs + 1; i++)
    delete [] _dist[i];

  delete [] _dist;     _dist    = NULL;
  delete [] _distLen;  _distLen = NULL;
  delete [] _distMax;  _distMax = NULL;
}



InsertSizes::~InsertSizes() {
  delete [] _mean;
  delete [] _stddev;
  delete [] _samples;
}


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
 *    Brian P. Walenz from 2003-JAN-02 to 2003-MAY-06
 *      are Copyright 2003 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-APR-21 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2006-JUL-07 to 2014-APR-11
 *      are Copyright 2006-2008,2011,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "positionDB.H"
#include "bio++.H"


void
adjustHeap(uint64 *C,
           uint64 *P, int64 i, int64 n) {
  uint64  c = C[i];
  uint64  p = P[i];
  int64  j = (i << 1) + 1;  //  let j be the left child

  while (j < n) {
    if (j<n-1 && C[j] < C[j+1])
      j++;                   //  j is the larger child

    if (c >= C[j])           //  a position for M[i] has been found
      break;

    C[(j-1)/2] = C[j];       //  Move larger child up a level
    P[(j-1)/2] = P[j];

    j = (j << 1) + 1;
  }

  C[(j-1)/2] = c;
  P[(j-1)/2] = p;
}


void
positionDB::sortAndRepackBucket(uint64 b) {
  uint64 st = _bucketSizes[b];
  uint64 ed = _bucketSizes[b+1];
  uint32 le = (uint32)(ed - st);

  if (ed < st)
    fprintf(stdout, "ERROR: Bucket "uint64FMT" starts at "uint64FMT" ends at "uint64FMT"?\n", b, st, ed);

  if (le == 0)
    return;

  //  One mer in the list?  It's distinct and unique!  (and doesn't
  //  contribute to the position list space count)
  //
  if (le == 1) {
    _numberOfDistinct++;
    _numberOfUnique++;
    return;
  }

  //  Allocate more space, if we need to.
  //
  if (_sortedMax <= le) {
    _sortedMax = le + 1024;
    delete [] _sortedChck;
    delete [] _sortedPosn;
    _sortedChck = new uint64 [_sortedMax];
    _sortedPosn = new uint64 [_sortedMax];
  }

  //  Unpack the bucket
  //
  uint64   lens[3] = {_chckWidth, _posnWidth, 1 + _sizeWidth};
  uint64   vals[3] = {0};
  for (uint64 i=st, J=st * _wCnt; i<ed; i++, J += _wCnt) {
    getDecodedValues(_countingBuckets, J, 2, lens, vals);
    _sortedChck[i-st] = vals[0];
    _sortedPosn[i-st] = vals[1];
  }

  //  Create the heap of lines.
  //
  int unsetBucket = 0;

  for (int64 t=(le-2)/2; t>=0; t--) {
    if (_sortedPosn[t] == uint64MASK(_posnWidth)) {
      unsetBucket = 1;
      fprintf(stdout, "ERROR: unset posn bucket="uint64FMT" t="int64FMT" le="uint32FMT"\n", b, t, le);
    }

    adjustHeap(_sortedChck, _sortedPosn, t, le);
  }

  if (unsetBucket)
    for (uint32 t=0; t<le; t++)
      fprintf(stdout, uint32FMTW(4)"] chck="uint64HEX" posn="uint64FMT"\n", t, _sortedChck[t], _sortedPosn[t]);

  //  Interchange the new maximum with the element at the end of the tree
  //
  for (int64 t=le-1; t>0; t--) {
    uint64           tc = _sortedChck[t];
    uint64           tp = _sortedPosn[t];

    _sortedChck[t]      = _sortedChck[0];
    _sortedPosn[t]      = _sortedPosn[0];

    _sortedChck[0]      = tc;
    _sortedPosn[0]      = tp;

    adjustHeap(_sortedChck, _sortedPosn, 0, t);
  }

  //  Scan the list of sorted mers, counting the number of distinct and unique,
  //  and the space needed in the position list.

  uint64   entries = 1;  //  For t=0

  for (uint32 t=1; t<le; t++) {
    if (_sortedChck[t-1] > _sortedChck[t])
      fprintf(stdout, "ERROR: bucket="uint64FMT" t="uint32FMT" le="uint32FMT": "uint64HEX" > "uint64HEX"\n",
              b, t, le, _sortedChck[t-1], _sortedChck[t]);

    if (_sortedChck[t-1] != _sortedChck[t]) {
      _numberOfDistinct++;

      if (_maximumEntries < entries)
        _maximumEntries = entries;

      if (entries == 1)
        _numberOfUnique++;
      else
        _numberOfEntries += entries + 1;  //  +1 for the length

      entries = 0;
    }

    entries++;
  }

  //  Don't forget the last mer!
  //
  _numberOfDistinct++;
  if (_maximumEntries < entries)
    _maximumEntries = entries;
  if (entries == 1)
    _numberOfUnique++;
  else
    _numberOfEntries += entries + 1;


  //  Repack the sorted entries
  //
  for (uint64 i=st, J=st * _wCnt; i<ed; i++, J += _wCnt) {
    vals[0] = _sortedChck[i-st];
    vals[1] = _sortedPosn[i-st];
    vals[2] = 0;
    setDecodedValues(_countingBuckets, J, 3, lens, vals);
  }
}


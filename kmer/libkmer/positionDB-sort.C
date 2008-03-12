#include "positionDB.H"
#include "bio++.H"


void
adjustHeap(u64bit *C,
           u64bit *P, s64bit i, s64bit n) {
  u64bit  c = C[i];
  u64bit  p = P[i];
  s64bit  j = (i << 1) + 1;  //  let j be the left child

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
positionDB::sortAndRepackBucket(u64bit b) {
  u64bit st = _bucketSizes[b];
  u64bit ed = _bucketSizes[b+1];
  u32bit le = (u32bit)(ed - st);

  if (ed < st)
    fprintf(stdout, "ERROR: Bucket %10lu starts at %10lu ends at %10lu?\n", b, st, ed);

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
    _sortedChck = new u64bit [_sortedMax];
    _sortedPosn = new u64bit [_sortedMax];
  }

  //  Unpack the bucket
  //
  u64bit   lens[3] = {_chckWidth, _posnWidth, 1 + _sizeWidth};
  u64bit   vals[3] = {0};
  for (u64bit i=st, J=st * _wCnt; i<ed; i++, J += _wCnt) {
    getDecodedValues(_countingBuckets, J, 2, lens, vals);
    _sortedChck[i-st] = vals[0];
    _sortedPosn[i-st] = vals[1];
  }

  //  Create the heap of lines.
  //
  int unsetBucket = 0;

  for (s64bit t=(le-2)/2; t>=0; t--) {
    if (_sortedPosn[t] == u64bitMASK(_posnWidth)) {
      unsetBucket = 1;
      fprintf(stdout, "ERROR: unset posn bucket="u64bitFMT" t="s64bitFMT" le="u32bitFMT"\n", b, t, le);
    }

    adjustHeap(_sortedChck, _sortedPosn, t, le);
  }

  if (unsetBucket)
    for (u32bit t=0; t<le; t++)
      fprintf(stdout, u32bitFMTW(4)"] chck="u64bitHEX" posn="u64bitFMT"\n", t, _sortedChck[t], _sortedPosn[t]);

  //  Interchange the new maximum with the element at the end of the tree
  //
  for (s64bit t=le-1; t>0; t--) {
    u64bit           tc = _sortedChck[t];
    u64bit           tp = _sortedPosn[t];

    _sortedChck[t]      = _sortedChck[0];
    _sortedPosn[t]      = _sortedPosn[0];

    _sortedChck[0]      = tc;
    _sortedPosn[0]      = tp;

    adjustHeap(_sortedChck, _sortedPosn, 0, t);
  }

  //  Scan the list of sorted mers, counting the number of distinct and unique,
  //  and the space needed in the position list.

  u64bit   entries = 1;  //  For t=0

  for (u32bit t=1; t<le; t++) {
    if (_sortedChck[t-1] > _sortedChck[t])
      fprintf(stdout, "ERROR: bucket="u64bitFMT" t="u32bitFMT" le="u32bitFMT": "u64bitHEX" > "u64bitHEX"\n",
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
  for (u64bit i=st, J=st * _wCnt; i<ed; i++, J += _wCnt) {
    vals[0] = _sortedChck[i-st];
    vals[1] = _sortedPosn[i-st];
    vals[2] = 0;
    setDecodedValues(_countingBuckets, J, 3, lens, vals);
  }
}


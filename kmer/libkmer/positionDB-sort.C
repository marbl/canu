#include "positionDB.H"
#include "bio++.H"

//#define ERROR_CHECK
//#define SORT_CHECK
//#define CHECK_FOR_EQUALITY

void
adjustHeap(heapbit *M, s64bit i, s64bit n) {
  heapbit   m = M[i];
  s64bit    j = (i << 1) + 1;  //  let j be the left child

  while (j < n) {
    if (j<n-1 && M[j] < M[j+1])
      j++;                   //  j is the larger child

    if (m >= M[j])           //  a position for M[i] has been found
      break;

    M[(j-1)/2] = M[j];       //  Move larger child up a level

    j = (j << 1) + 1;
  }

  M[(j-1)/2] = m;
}

#ifdef SORT_CHECK
int
cmpheap(void const *a, void const *b) {
  heapbit A = *((heapbit const *)a);
  heapbit B = *((heapbit const *)b);

  if (A < B) return(-1);
  if (A > B) return(1);
  return(0);
}
#endif

void
positionDB::sortAndRepackBucket(u64bit b) {
  u64bit st = _bucketSizes[b];
  u64bit ed = _bucketSizes[b+1];

#ifdef ERROR_CHECK
  //  This passes!
  //
  if (ed < st)
    fprintf(stdout, "ERROR: Bucket %10lu starts at %10lu ends at %10lu?\n", b, st, ed);
#endif

  _sortedListLen = (u32bit)(ed - st);

  //  No mers in the list?  We're done.
  //
  if (_sortedListLen == 0)
    return;

  //  One mer in the list?  It's distinct and unique!  (and doesn't contribute
  //  to the position list space count)
  //
  if (_sortedListLen == 1) {
    _numberOfDistinct++;
    _numberOfUnique++;
    return;
  }

  //  Allocate more space, if we need to.
  //
  if (_sortedListLen > _sortedListMax) {
    delete [] _sortedList;
    _sortedList    = new heapbit [_sortedListLen + 1];
    _sortedListMax = _sortedListLen;
  }

  //  Unpack the check values
  // 
  for (u64bit i=st, J=st * _wCnt; i<ed; i++, J += _wCnt)
   _sortedList[i-st] = getDecodedValue(_countingBuckets, J, _wCnt);


#ifdef CHECK_FOR_EQUALITY
  //  Check for duplicate entries in the _sortedList.
  //
  u32bit duplicates = 0;

  for (u32bit i=0; i<_sortedListLen; i++)
    for (u32bit j=i+1; j<_sortedListLen; j++)
      if (_sortedList[i] == _sortedList[j])
        duplicates++;

  if (duplicates > 0)
    fprintf(stdout, "Found %u (pairs of) duplicates!\n", duplicates);
#endif



#ifdef ERROR_CHECK
  int unsetBucket = 0;
  for (u32bit t=0; t<_sortedListLen; t++)
    if ((_sortedList[t] & _posnMask) == _posnMask) {
      unsetBucket = 1;
      fprintf(stdout, "WARNING!  Unset countingBucket --i=%lu len=%lu 0x%016lx\n", t, _sortedListLen, _sortedList[t]);
      fprintf(stdout, "          Bucket %10lu starts at %10lu ends at %10lu\n\n", b, st, ed);
    }
#endif


#ifdef ERROR_CHECK
  if (unsetBucket) {
    fprintf(stdout, "entries are BEFORE:\n");
    for (u32bit t=0; t<_sortedListLen; t++)
      fprintf(stdout, "%4u] 0x%016lx\n", t, _sortedList[t]);
    fprintf(stdout, "\n");
  }
#endif



#ifdef SORT_CHECK
  //  Error check the sorting.  Sort once using qsort (which we
  //  assume is correct), then sort it again with heapsort.
  //  compare the two outputs.
  //

  heapbit *sortTest = new heapbit [_sortedListLen];
  for (u32bit t=0; t<_sortedListLen; t++)
    sortTest[t] = _sortedList[t];

  qsort(sortTest, _sortedListLen, sizeof(heapbit), cmpheap);
#endif



  //  Create the heap of lines.
  //
  for (s64bit t=(_sortedListLen-2)/2; t>=0; t--)
    adjustHeap(_sortedList, t, _sortedListLen);

  //  Interchange the new maximum with the element at the end of the tree
  //
  for (s64bit t=_sortedListLen-1; t>0; t--) {
    heapbit          tv = _sortedList[t];
    _sortedList[t]      = _sortedList[0];
    _sortedList[0]      = tv;

    adjustHeap(_sortedList, 0, t);
  }


#if 0
  //  Define this (and SORT_CHECK) to use qsort() exclusively
  for (u32bit t=0; t<_sortedListLen; t++)
    _sortedList[t] = sortTest[t];

#endif


#ifdef SORT_CHECK
  for (u32bit t=0; t<_sortedListLen; t++)
    if (sortTest[t] != _sortedList[t])
      fprintf(stderr, "ERROR with sort: %5u/%5u\n", t, _sortedListLen);
  delete [] sortTest;
#endif

#ifdef ERROR_CHECK
  if (unsetBucket) {
    fprintf(stdout, "entries are AFTER:\n");
    for (u32bit t=0; t<_sortedListLen; t++)
      fprintf(stdout, "%4u] 0x%016lx\n", t, _sortedList[t]);
    fprintf(stdout, "\n");
  }
#endif

#ifdef ERROR_CHECK
  //  Are we actually sorted?
  //
  //  We get EQUAL (on real data), and lots of:
  //    ERROR:    38/ 174: 0x0001ffffffffffff == 0x0001ffffffffffff
  //
  //  ON BOTH versions.
  //
  for (u32bit t=1; t<_sortedListLen; t++) {
    if (_sortedList[t-1] == _sortedList[t])
      fprintf(stdout, "ERROR: %6u %4u/%4u: 0x%016llx == 0x%016llx\n", b, t, _sortedListLen, _sortedList[t-1], _sortedList[t]);
    if (_sortedList[t-1] > _sortedList[t])
      fprintf(stdout, "ERROR: %6u %4u/%4u: 0x%016llx  > 0x%016llx\n", b, t, _sortedListLen, _sortedList[t-1], _sortedList[t]);
  }
#endif



  //  First was original, second has the same effect
  //
  //u64bit  checkMask = ~_posnMask;
  u64bit  checkMask = _mask2 << _posnWidth;

  //
  //  Scan the list of sorted mers, counting the number of distinct and unique,
  //  and the space needed in the position list.
  //

  //  Count the first mer.
  //
  u64bit   entries = 1;

  for (u32bit t=1; t<_sortedListLen; t++) {

    //  If the current check is not the last check, then we have a new
    //  mer.  Update the counts.
    //
    if ((_sortedList[t-1] & checkMask) != (_sortedList[t] & checkMask)) {
      _numberOfDistinct++;

      if (_maximumEntries < entries)
        _maximumEntries = entries;

      if (entries == 1)
        _numberOfUnique++;
      else
        _numberOfEntries += entries + 1;

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
  for (u64bit i=st, J=st * _wCnt; i<ed; i++, J += _wCnt)
    setDecodedValue(_countingBuckets, J, _wCnt, _sortedList[i-st]);
}


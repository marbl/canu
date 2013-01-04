#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "existDB.H"
#include "bio++.H"
#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

bool
existDB::createFromSequence(char const  *sequence,
                            u32bit       merSize,
                            u32bit       flags) {

  bool               beVerbose  = false;
  bool               rebuilding = false;

  _hashTable = 0L;
  _buckets   = 0L;
  _counts    = 0L;

  _merSizeInBases        = merSize;

  _searchForDupe = true;

  if ((flags & existDBcompressHash) ||
      (flags & existDBcompressBuckets) ||
      (flags & existDBcompressCounts))
    fprintf(stderr, "existDB::createFromSequence: compression not supported.\n"), exit(1);

  //  This (at =22) eats up 16MB, and should allow a lot of mers at big sizes.  Unfortunately, we
  //  know nothing about how man mers are going to be in the input.
  //
  //  Setting this too high drastically reduces performance, suspected because of cache misses.
  //  Setting this too low will also reduce performance, by increasing the search time in a bucket.
  //
  u32bit tblBits = logBaseTwo64(strlen(sequence));

 rebuild:
  fprintf(stderr, "tblBits %d seqlen %d\n", tblBits, strlen(sequence));

  _shift1                = 2 * _merSizeInBases - tblBits;
  _shift2                = _shift1 / 2;
  _mask1                 = u64bitMASK(tblBits);
  _mask2                 = u64bitMASK(_shift1);

  _hshWidth              = u32bitZERO;
  _chkWidth              = 2 * merSize - tblBits;
  _cntWidth              = 16;

  u64bit  tableSizeInEntries = u64bitONE << tblBits;
  u64bit  numberOfMers       = u64bitZERO;
  u64bit *countingTable      = new u64bit [tableSizeInEntries + 1];

  for (u64bit i=tableSizeInEntries+1; i--; )
    countingTable[i] = 0;

  _isCanonical = flags & existDBcanonical;
  _isForward   = flags & existDBforward;

  assert(_isCanonical + _isForward == 1);

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  1)  Count bucket sizes
  //
  merStream    *M = new merStream(new kMerBuilder(_merSizeInBases),
                                  new seqStream(sequence, strlen(sequence)),
                                  true, true);

  while (M->nextMer()) {
    if (_isForward) {
      countingTable[ HASH(M->theFMer()) ]++;
      numberOfMers++;
    }

    if (_isCanonical) {
      countingTable[ HASH(M->theCMer()) ]++;
      numberOfMers++;
    }
  }

  delete M;

#ifdef STATS
  u64bit  dist[32] = {0};
  u64bit  maxcnt = 0;
  for (u64bit i=tableSizeInEntries+1; i--; ) {
    if (countingTable[i] > maxcnt)
      maxcnt = countingTable[i];

    if (countingTable[i] < 32)
      dist[countingTable[i]]++;
  }

  for(u64bit i=0; i<32; i++)
    fprintf(stderr, "existDB::usage[%2d] = %d\n", i, dist[i]);
  fprintf(stderr, "existDB::maxcnt    = %d\n", maxcnt);
#endif


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Determine how many bits we need to hold the value
  //  numberOfMers.....then....
  //
  //  This is numberOfMers+1 because we need to store the
  //  first position after the last mer.  That is, if there are two
  //  mers, we will store that the first mer is at position 0, the
  //  second mer is at position 1, and the end of the second mer is at
  //  position 2.
  //
  if (_compressedHash) {
    _hshWidth = 1;
    while ((numberOfMers+1) > (u64bitONE << _hshWidth))
      _hshWidth++;
  }


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  2)  Allocate a hash table and some mer storage buckets.
  //
  _hashTableWords = tableSizeInEntries + 2;
  if (_compressedHash)
    _hashTableWords = _hashTableWords * _hshWidth / 64 + 1;

  _bucketsWords = numberOfMers + 2;
  if (_compressedBucket)
    _bucketsWords = _bucketsWords * _chkWidth / 64 + 1;

  _countsWords = numberOfMers + 2;
  if (_compressedCounts)
    _countsWords = _countsWords * _cntWidth / 64 + 1;

  if (beVerbose) {
    fprintf(stderr, "existDB::createFromSequence()-- hashTable is "u64bitFMT"MB\n", _hashTableWords >> 17);
    fprintf(stderr, "existDB::createFromSequence()-- buckets is "u64bitFMT"MB\n", _bucketsWords >> 17);
    if (flags & existDBcounts)
      fprintf(stderr, "existDB::createFromSequence()-- counts is "u64bitFMT"MB\n", _countsWords >> 17);
  }

  _hashTable   = new u64bit [_hashTableWords];
  _buckets     = new u64bit [_bucketsWords];
  _countsWords = (flags & existDBcounts) ?             _countsWords  : 0;
  _counts      = (flags & existDBcounts) ? new u64bit [_countsWords] : 0L;

  //  These aren't strictly needed.  _buckets is cleared as it is initialied.  _hashTable
  //  is also cleared as it is initialized, but in the _compressedHash case, the last
  //  few words might be uninitialized.  They're unused.

  //memset(_hashTable, 0, sizeof(u64bit) * _hashTableWords);
  //memset(_buckets,   0, sizeof(u64bit) * _bucketsWords);  //  buckets is cleared as it is built
  //memset(_counts,    0, sizeof(u64bit) * _countsWords);

  _hashTable[_hashTableWords-1] = 0;
  _hashTable[_hashTableWords-2] = 0;
  _hashTable[_hashTableWords-3] = 0;
  _hashTable[_hashTableWords-4] = 0;

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Make the hash table point to the start of the bucket, and reset
  //  the counting table -- we're going to use it to fill the buckets.
  //
  u64bit  tmpPosition = 0;
  u64bit  begPosition = 0;
  u64bit  ptr         = 0;

  if (_compressedHash) {
    for (u64bit i=0; i<tableSizeInEntries; i++) {
      tmpPosition    = countingTable[i];
      countingTable[i] = begPosition;

      setDecodedValue(_hashTable, ptr, _hshWidth, begPosition);
      ptr         += _hshWidth;

      begPosition += tmpPosition;
    }

    setDecodedValue(_hashTable, ptr, _hshWidth, begPosition);
  } else {
    for (u64bit i=0; i<tableSizeInEntries; i++) {
      tmpPosition    = countingTable[i];
      countingTable[i] = begPosition;

      _hashTable[i] = begPosition;

      begPosition += tmpPosition;
    }

    //  Set the last position in the hash, but we don't care about
    //  the temporary counting table.
    //
    _hashTable[tableSizeInEntries] = begPosition;
  }




  ////////////////////////////////////////////////////////////////////////////////
  //
  //  3)  Build list of mers, placed into buckets
  //
  M  = new merStream(new kMerBuilder(_merSizeInBases),
                     new seqStream(sequence, strlen(sequence)),
                     true, true);

  while (M->nextMer()) {
    if (_isForward)
      insertMer(HASH(M->theFMer()), CHECK(M->theFMer()), 1, countingTable);

    if (_isCanonical)
      insertMer(HASH(M->theCMer()), CHECK(M->theCMer()), 1, countingTable);
  }

  delete M;

  //  Compress out the gaps we have from redundant kmers.

  u64bit  pos = 0;
  u64bit  frm = 0;
  u64bit  len = 0;

  for (u64bit i=0; i<tableSizeInEntries; i++) {
    frm = _hashTable[i];
    len = countingTable[i] - _hashTable[i];

    _hashTable[i] = pos;

    for (u64bit j=0; j<len; j++) {
      if (_counts)
        _counts[pos]  = _counts[frm];

      _buckets[pos++] = _buckets[frm++];
    }
  }

  fprintf(stderr, "Compressed from "u64bitFMT" to "u64bitFMT" ("u64bitFMT" bits)\n",
          _hashTable[tableSizeInEntries], pos, logBaseTwo64(pos));

  while (pos < _bucketsWords)
    _buckets[pos++] = 0;

  _hashTable[tableSizeInEntries] = pos;

  //  All done.  Delete temporary stuff
  //
  delete [] countingTable;

  //  But if we horribly screwed up the estimate of tblBits, reset and recompute

  if ((logBaseTwo64(pos) < tblBits) &&
      (rebuilding == false)) {
    rebuilding = true;

    delete [] _hashTable;
    delete [] _buckets;
    delete [] _counts;

    _hashTable = 0L;
    _buckets   = 0L;
    _counts    = 0L;

    tblBits = logBaseTwo64(pos);

    goto rebuild;
  }

  return(true);
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "existDB.H"
#include "libmeryl.H"


bool
existDB::createFromMeryl(char const  *prefix,
                         u32bit       lo,
                         u32bit       hi,
                         u32bit       flags) {

  merylStreamReader *M = new merylStreamReader(prefix);

  bool               beVerbose = true;

  _hashTable = 0L;
  _buckets   = 0L;
  _counts    = 0L;

  _merSizeInBases        = M->merSize();

  //  We can set this exactly, but not memory optimal (see
  //  meryl/estimate.C:optimalNumberOfBuckets()).  Instead,
  //  we just blindly use whatever meryl used.
  //
  //u32bit tblBits = logBaseTwo64(M->numberOfDistinctMers() + 1);
  //
  u32bit tblBits = M->prefixSize();

  _shift1                = 2 * _merSizeInBases - tblBits;
  _shift2                = _shift1 / 2;
  _mask1                 = u64bitMASK(tblBits);
  _mask2                 = u64bitMASK(_shift1);

  _hshWidth              = u32bitZERO;
  _chkWidth              = 2 * _merSizeInBases - tblBits;
  _cntWidth              = 16;

  u64bit  tableSizeInEntries = u64bitONE << tblBits;
  u64bit  numberOfMers       = u64bitZERO;
  u64bit  maxCount           = u64bitZERO;
  u64bit *countingTable      = new u64bit [tableSizeInEntries + 1];

  for (u64bit i=tableSizeInEntries+1; i--; )
    countingTable[i] = 0;

  _isCanonical = flags & existDBcanonical;
  _isForward   = flags & existDBforward;

  //  1) Count bucket sizes
  //     While we don't know the bucket sizes right now, but we do know
  //     how many buckets and how many mers.
  //
  //  Because we could be inserting both forward and reverse, we can't
  //  really move the direction testing outside the loop, unless we
  //  want to do two iterations over M.
  //
  speedCounter  *C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);

  while (M->nextMer()) {
    if ((lo <= M->theCount()) && (M->theCount() <= hi)) {
      if (_isForward) {
        countingTable[ HASH(M->theFMer()) ]++;
        numberOfMers++;
      }

      if (_isCanonical) {
        kMer  r = M->theFMer();
        r.reverseComplement();

        if (M->theFMer() < r)
          countingTable[ HASH(M->theFMer()) ]++;
        else
          countingTable[ HASH(r) ]++;
        numberOfMers++;
      }

      if (maxCount < M->theCount())
        maxCount = M->theCount();

      C->tick();
    }
  }

  delete C;
  delete M;

  if (_compressedHash) {
    _hshWidth = 1;
    while ((numberOfMers+1) > (u64bitONE << _hshWidth))
      _hshWidth++;
  }

  if (beVerbose) {
    fprintf(stderr, "existDB::createFromMeryl()-- Found "u64bitFMT" mers between count of "u32bitFMT" and "u32bitFMT"\n",
            numberOfMers, lo, hi);
  }

  //  2) Allocate hash table, mer storage buckets
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
    fprintf(stderr, "existDB::createFromMeryl()-- hashTable is "u64bitFMT"MB\n", _hashTableWords >> 17);
    fprintf(stderr, "existDB::createFromMeryl()-- buckets is "u64bitFMT"MB\n", _bucketsWords >> 17);
    if (flags & existDBcounts)
      fprintf(stderr, "existDB::createFromMeryl()-- counts is "u64bitFMT"MB\n", _countsWords >> 17);
  }

  _hashTable = new u64bit [_hashTableWords];
  _buckets   = new u64bit [_bucketsWords];
  _counts    = (flags & existDBcounts) ? new u64bit [_countsWords] : 0L;

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


  ///////////////////////////////////////////////////////////////////////////////
  //
  //  3)  Build list of mers, placed into buckets
  //
  M = new merylStreamReader(prefix);
  C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);

  while (M->nextMer()) {
    if ((lo <= M->theCount()) && (M->theCount() <= hi)) {
      if (_isForward)
        insertMer(HASH(M->theFMer()), CHECK(M->theFMer()), M->theCount(), countingTable);

      if (_isCanonical) {
        kMer  r = M->theFMer();
        r.reverseComplement();

        if (M->theFMer() < r)
          insertMer(HASH(M->theFMer()), CHECK(M->theFMer()), M->theCount(), countingTable);
        else
          insertMer(HASH(r), CHECK(r), M->theCount(), countingTable);
        numberOfMers++;
      }


      C->tick();
    }
  }

  delete C;
  delete M;
  delete [] countingTable;

  return(true);
}

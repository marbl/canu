#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "existDB.H"
#include "positionDB.H"
#include "libmeryl.H"


bool
existDB::createFromMeryl(char const  *prefix,
                         u32bit       lo,
                         u32bit       hi,
                         u32bit       tblBits) {

  merylStreamReader *M = new merylStreamReader(prefix);

  _hashTable = 0L;
  _buckets   = 0L;

  _merSizeInBases        = M->merSize();

  _shift1                = 2 * _merSizeInBases - tblBits;
  _shift2                = _shift1 / 2;
  _mask1                 = u64bitMASK(tblBits);
  _mask2                 = u64bitMASK(_shift1);

  _hashWidth             = u32bitZERO;
  _chckWidth             = 2 * _merSizeInBases - tblBits;

  _hashMask              = u64bitMASK(tblBits);
  _chckMask              = u64bitMASK(2 * _merSizeInBases - tblBits);

  u64bit  tableSizeInEntries = u64bitONE << tblBits;
  u64bit  numberOfMers       = u64bitZERO;
  u64bit *countingTable      = new u64bit [tableSizeInEntries + 1];

#if 1
  fprintf(stderr, "existDB::createFromMeryl()-- countingTable is "u64bitFMT"MB\n",
          tableSizeInEntries >> 17);
#endif

  for (u64bit i=tableSizeInEntries+1; i--; )
    countingTable[i] = 0;


  //  1) Count bucket sizes
  //     While we don't know the bucket sizes right now, but we do know
  //     how many buckets and how many mers.  Unfortunately, we still
  //     need to HASH() each of the mers, so that we can use the existing
  //     exists() method.
  //

  while (M->nextMer()) {
    if ((lo <= M->theCount()) && (M->theCount() <= hi)) {
      countingTable[ HASH(M->theFMer()) ]++;
      numberOfMers++;
    }
  }

  delete M;
  
  if (_compressedHash) {
    _hashWidth = 1;
    while ((numberOfMers+1) > (u64bitONE << _hashWidth))
      _hashWidth++;
  }

  fprintf(stderr, "existDB: Found "u64bitFMT" mers between count of "u32bitFMT" and "u32bitFMT"\n",
          numberOfMers, lo, hi);


  //  2) Allocate hash table, mer storage buckets
  //
  _hashTableWords = tableSizeInEntries + 2;
  if (_compressedHash)
    _hashTableWords = _hashTableWords * _hashWidth / 64 + 1;

  _bucketsWords = numberOfMers + 2;
  if (_compressedBucket)
    _bucketsWords = _bucketsWords * _chckWidth / 64 + 1;

  _hashTable = new u64bit [_hashTableWords];
  _buckets   = new u64bit [_bucketsWords];

#if 1
  fprintf(stderr, "existDB::createFromMeryl()-- hashTable is "u64bitFMT"MB\n",
          _hashTableWords >> 17);
  fprintf(stderr, "existDB::createFromMeryl()-- buckets is "u64bitFMT"MB\n",
          _bucketsWords >> 17);
#endif

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

      setDecodedValue(_hashTable, ptr, _hashWidth, begPosition);
      ptr         += _hashWidth;

      begPosition += tmpPosition;
    }

    setDecodedValue(_hashTable, ptr, _hashWidth, begPosition);
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

  u64bit  h;

  while (M->nextMer()) {
    if ((lo <= M->theCount()) && (M->theCount() <= hi)) {
      h = HASH(M->theFMer());

      if (_compressedBucket)
        setDecodedValue(_buckets,
                        countingTable[h] * _chckWidth,
                        _chckWidth,
                        CHECK(M->theFMer()));
      else
        _buckets[countingTable[h]] = CHECK(M->theFMer());

      countingTable[h]++;
    }
  }

  delete M;
  delete [] countingTable;

  return(true);
}

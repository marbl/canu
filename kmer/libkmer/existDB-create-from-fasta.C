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
existDB::createFromFastA(char const  *filename,
                         u32bit       merSize,
                         u32bit       flags) {

  _hashTable = 0L;
  _buckets   = 0L;
  _counts    = 0L;

  _merSizeInBases        = merSize;

  //  This eats up 16MB, and should allow a lot of mers at big sizes.
  //  Unfortunately, we know nothing about how man mers are going to
  //  be in the input.
  //
  u32bit tblBits = 22;

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

  if (flags & existDBcounts) {
    fprintf(stderr, "existDB:createFromFastA()--  ERROR!  I don't support existDBcounts.\n");
    exit(1);
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  1)  Count bucket sizes
  //
  merStream    *M = new merStream(new kMerBuilder(_merSizeInBases),
                                  new seqStream(filename),
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

  _hashTable = new u64bit [_hashTableWords];
  _buckets   = new u64bit [_bucketsWords];

#ifdef STATS
  fprintf(stderr, "existDB::allocated %lu bytes of storage.\n", 8 * (_hashTableWords + _bucketsWords));
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
                     new seqStream(filename),
                     true, true);

  while (M->nextMer()) {
    if (_isForward)
      insertMer(HASH(M->theFMer()), CHECK(M->theFMer()), 0, countingTable);

    if (_isCanonical)
      insertMer(HASH(M->theCMer()), CHECK(M->theCMer()), 0, countingTable);
  }

  delete M;

  //  All done.  Delete temporary stuff
  //
  delete [] countingTable;

  return(true);
}

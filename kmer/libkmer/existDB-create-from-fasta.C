#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "existDB.H"
#include "positionDB.H"
#include "bio++.H"


bool
existDB::createFromFastA(char const  *filename,
                         u32bit       merSize,
                         u32bit       tblBits,
                         u32bit       flags) {

  _hashTable = 0L;
  _buckets   = 0L;

  _merSizeInBases        = merSize;

  _shift1                = 2 * _merSizeInBases - tblBits;
  _shift2                = _shift1 / 2;
  _mask1                 = u64bitMASK(tblBits);
  _mask2                 = u64bitMASK(_shift1);

  _hashWidth             = u32bitZERO;
  _chckWidth             = 2 * merSize - tblBits;

  _hashMask              = u64bitMASK(tblBits);
  _chckMask              = u64bitMASK(2 * merSize - tblBits);

  u64bit  tableSizeInEntries = u64bitONE << tblBits;
  u64bit  numberOfMers       = u64bitZERO;
  u64bit *countingTable      = new u64bit [tableSizeInEntries + 1];

  if (_beVerbose)
    fprintf(stderr, "existDB::createFromFastA()-- countingTable is "u64bitFMT"MB\n",
            tableSizeInEntries >> 17);

  for (u64bit i=tableSizeInEntries+1; i--; )
    countingTable[i] = 0;

  bool  doCanonical = flags & existDBcanonical;
  bool  doForward   = flags & existDBforward;
  bool  doReverse   = flags & existDBreverse;

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  1)  Count bucket sizes
  //
  FastAstream   *F = new FastAstream(filename);
  merStream     *M = new merStream(_merSizeInBases, F);

  while (M->nextMer()) {
    if (doForward) {
      countingTable[ HASH(M->theFMer()) ]++;
      numberOfMers++;
    }

    if (doReverse) {
      countingTable[ HASH(M->theRMer()) ]++;
      numberOfMers++;
    }

    if (doCanonical) {
      fprintf(stderr, "ERROR:  canonical mers in existDB not implemented.\n");
      exit(1);
    }
  }

  delete M;
  delete F;

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
    _hashWidth = 1;
    while ((numberOfMers+1) > (u64bitONE << _hashWidth))
      _hashWidth++;
  }


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  2)  Allocate a hash table and some mer storage buckets.
  //
  _hashTableWords = tableSizeInEntries + 2;
  if (_compressedHash)
    _hashTableWords = _hashTableWords * _hashWidth / 64 + 1;

  _bucketsWords = numberOfMers + 2;
  if (_compressedBucket)
    _bucketsWords = _bucketsWords * _chckWidth / 64 + 1;

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




  ////////////////////////////////////////////////////////////////////////////////
  //
  //  3)  Build list of mers, placed into buckets
  //
  F = new FastAstream(filename);
  M = new merStream(_merSizeInBases, F);

  while (M->nextMer()) {
    if (doForward)
      INSERT(HASH(M->theFMer()), CHECK(M->theFMer()), countingTable);

    if (doReverse)
      INSERT(HASH(M->theRMer()), CHECK(M->theRMer()), countingTable);

    if (doCanonical)
      ;  //  Not implemented, caught above
  }

  delete M;
  delete F;

  //  All done.  Delete temporary stuff
  //
  delete [] countingTable;

  return(true);
}

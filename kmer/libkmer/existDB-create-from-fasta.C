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
                         positionDB  *posDB) {

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

  for (u64bit i=tableSizeInEntries+1; i--; )
    countingTable[i] = 0;



  ////////////////////////////////////////////////////////////////////////////////
  //
  //  1)  Count bucket sizes
  //
  FastAstream   *F = new FastAstream(filename);
  merStream     *M = new merStream(_merSizeInBases, F);

  if (posDB) {
    while (M->nextMer()) {
      if (posDB->exists(M->theFMer())) {
        countingTable[ HASH(M->theFMer()) ]++;
        numberOfMers++;
      }
      if (posDB->exists(M->theRMer())) {
        countingTable[ HASH(M->theRMer()) ]++;
        numberOfMers++;
      }
    }
  } else {
    while (M->nextMer()) {
      countingTable[ HASH(M->theFMer()) ]++;
      numberOfMers++;
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

  //  XXX:  Pretty big code bloat here
  //
  //  I don't want to bury a test in the posDB on every mer,
  //  especially if we are never using the capability, as ESTmapper
  //  searchGENOME does.
  //
  //  So, I bloated the code (a little bit).
  //
  if (posDB) {
    u64bit  fmer;
    u64bit  rmer;
    u64bit  h;

    while (M->nextMer()) {
      fmer = M->theFMer();
      rmer = M->theRMer();

      if (posDB->exists(fmer)) {
        h = HASH(fmer);

        if (_compressedBucket)
          setDecodedValue(_buckets,
                          countingTable[h] * _chckWidth,
                          _chckWidth,
                          CHECK(fmer));
        else
          _buckets[countingTable[h]] = CHECK(fmer);

        countingTable[h]++;
      }
      if (posDB->exists(rmer)) {
        h = HASH(rmer);

        if (_compressedBucket)
          setDecodedValue(_buckets,
                          countingTable[h] * _chckWidth,
                          _chckWidth,
                          CHECK(rmer));
        else
          _buckets[countingTable[h]] = CHECK(rmer);

        countingTable[h]++;
      }
    }
  } else {
    u64bit  h;

    while (M->nextMer()) {
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
  delete F;

  //  All done.  Delete temporary stuff
  //
  delete [] countingTable;

  return(true);
}

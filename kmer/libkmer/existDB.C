#include "existDB.H"
#include "positionDB.H"
#include "merstream.H"

//  Print some statistics
//
//#define STATS

//
//  2x faster and 2x bigger if not compressed (for testexist.C)
//  8MB for 782439 mers (14-mers and 12-bit table)
//

#include "bit-packing.H"

existDB::existDB(char        *filename,
                 u32bit       merSize,
                 u32bit       tblBits,
                 positionDB  *posDB) {

  _hashTable = 0L;
  _buckets   = 0L;

  _merSizeInBases        = merSize;
  //_merSizeInBases_div2   = merSize >> 1;

  _shift1                = 2 * _merSizeInBases - tblBits;
  _shift2                = _shift1 / 2;
  _mask1                 = u64bitMASK(tblBits);
  _mask2                 = u64bitMASK(_shift1);

#ifdef COMPRESSED_HASH
  _hashWidth             = u32bitZERO;
#endif
#ifdef COMPRESSED_BUCKET
  _chckWidth             = 2 * merSize - tblBits;
#endif
  _hashMask              = u64bitMASK(tblBits);
  _chckMask              = u64bitMASK(2 * merSize - tblBits);

  u64bit  tableSizeInEntries = u64bitONE << tblBits;
  u64bit  numberOfMers       = u64bitZERO;
  u64bit *countingTable      = 0L;

#if 0
  fprintf(stderr, "merSize:    %u\n", _merSizeInBases);
  fprintf(stderr, "chckWidth:  %u\n", 2 * merSize - tblBits);
#endif

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  1)  Count bucket sizes
  //
  countingTable = new u64bit [tableSizeInEntries + 1];
  for (u64bit i=tableSizeInEntries+1; i--; )
    countingTable[i] = 0;

  merStream     *M = new merStream(_merSizeInBases, filename);

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
#ifdef COMPRESSED_HASH
  _hashWidth = 1;
  while ((numberOfMers+1) > (u64bitONE << _hashWidth))
    _hashWidth++;
#endif


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  2)  Allocate a hash table and some mer storage buckets.
  //
#ifdef COMPRESSED_HASH
  u64bit   hs = (tableSizeInEntries + 1) * _hashWidth / 64 + 1;
#else
  u64bit   hs = tableSizeInEntries + 2;
#endif

#ifdef COMPRESSED_BUCKET
  u64bit   bs = (numberOfMers + 1) * _chckWidth / 64 + 1;
#else
  u64bit   bs = numberOfMers + 2;
#endif

  _hashTable = new u64bit [hs];
  _buckets   = new u64bit [bs];


#ifdef STATS
  fprintf(stderr, "existDB::allocated %lu bytes of storage.\n", 8 * (hs + bs));
#endif

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Make the hash table point to the start of the bucket, and reset
  //  the counting table -- we're going to use it to fill the buckets.
  //
  u64bit  tmpPosition = 0;
  u64bit  begPosition = 0;
#ifdef COMPRESSED_HASH
  u64bit  ptr         = 0;
#endif

  for (u64bit i=0; i<tableSizeInEntries; i++) {
    tmpPosition    = countingTable[i];
    countingTable[i] = begPosition;

#ifdef COMPRESSED_HASH
    setDecodedValue(_hashTable, ptr, _hashWidth, begPosition);
    ptr         += _hashWidth;
#else
    _hashTable[i] = begPosition;
#endif

    begPosition += tmpPosition;
  }

  //  Set the last position in the hash, but we don't care about
  //  the temporary counting table.
  //
#ifdef COMPRESSED_HASH
  setDecodedValue(_hashTable, ptr, _hashWidth, begPosition);
#else
  _hashTable[tableSizeInEntries] = begPosition;
#endif



  ////////////////////////////////////////////////////////////////////////////////
  //
  //  3)  Build list of mers, placed into buckets
  //
  M = new merStream(_merSizeInBases, filename);

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

#ifdef COMPRESSED_BUCKET
        setDecodedValue(_buckets,
                        countingTable[h] * _chckWidth,
                        _chckWidth,
                        CHECK(fmer));
#else
        _buckets[countingTable[h]] = CHECK(fmer);
#endif

        countingTable[h]++;
      }
      if (posDB->exists(rmer)) {
        h = HASH(rmer);

#ifdef COMPRESSED_BUCKET
        setDecodedValue(_buckets,
                        countingTable[h] * _chckWidth,
                        _chckWidth,
                        CHECK(rmer));
#else
        _buckets[countingTable[h]] = CHECK(rmer);
#endif

        countingTable[h]++;
      }
    }
  } else {
    u64bit  h;

    while (M->nextMer()) {
      h = HASH(M->theFMer());

#ifdef COMPRESSED_BUCKET
      setDecodedValue(_buckets,
                      countingTable[h] * _chckWidth,
                      _chckWidth,
                      CHECK(M->theFMer()));
#else
      _buckets[countingTable[h]] = CHECK(M->theFMer());
#endif

      countingTable[h]++;
    }
  }

  delete M;

  //  All done.  Delete temporary stuff
  //
  delete [] countingTable;
}


existDB::~existDB() {
  delete [] _hashTable;
  delete [] _buckets;
}



bool
existDB::exists(u64bit mer) {
#ifdef COMPRESSED_HASH
  u64bit  h = HASH(mer) * _hashWidth;
  u64bit st = getDecodedValue(_hashTable, h,              _hashWidth);
  u64bit ed = getDecodedValue(_hashTable, h + _hashWidth, _hashWidth);
#else
  u64bit  h = HASH(mer);
  u64bit st = _hashTable[h];
  u64bit ed = _hashTable[h+1];
#endif

  if (st == ed)
    return(false);

#ifdef COMPRESSED_BUCKET
  st *= _chckWidth;
  ed *= _chckWidth;
#endif

  u64bit  c = CHECK(mer);


#ifdef COMPRESSED_BUCKET
  for (; st<ed; st += _chckWidth) {
    if (getDecodedValue(_buckets, st, _chckWidth) == c)
      return(true);
  }
#else
  for (; st<ed; st++) {
    if (_buckets[st] == c)
      return(true);
  }
#endif

  return(false);
}


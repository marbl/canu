#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
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

char     magic[16] = { 'e', 'x', 'i', 's', 't', 'D', 'B', '1', 
                       ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '  };


existDB::existDB(char const  *filename,
                 bool         loadData) {
  if (loadState(filename, true, loadData) == false) {
    fprintf(stderr, "existDB::existDB()-- Tried to read state from '%s', but failed.\n", filename);
    exit(1);
  }
}

existDB::existDB(char const  *filename,
                 u32bit       merSize,
                 u32bit       tblBits,
                 positionDB  *posDB) {

  //  Try to read state from the filename.  If successful, make sure that
  //  the merSize and tblBits are correct.
  //
  if (loadState(filename)) {
    bool fail = false;

    if (_merSizeInBases != merSize) {
      fprintf(stderr, "existDB::existDB()-- Read state from '%s', but got different mer sizes\n", filename);
      fprintf(stderr, "existDB::existDB()-- Got %d, expected %d\n", _merSizeInBases, merSize);
      fail = true;
    }
    if (_mask1 != u64bitMASK(tblBits)) {
      fprintf(stderr, "existDB::existDB()-- Read state from '%s', but got different table sizes\n", filename);
      fprintf(stderr, "existDB::existDB()-- Got %d, expected %d\n", 2 * _merSizeInBases - _shift1, tblBits);
      fail = true;
    }

    if (fail)
      exit(1);

    return;
  }

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
  _hashTableWords = (tableSizeInEntries + 1) * _hashWidth / 64 + 1;
#else
  _hashTableWords = tableSizeInEntries + 2;
#endif

#ifdef COMPRESSED_BUCKET
  _bucketsWords = (numberOfMers + 1) * _chckWidth / 64 + 1;
#else
  _bucketsWords = numberOfMers + 2;
#endif

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



void
existDB::saveState(char const *filename) {

  errno = 0;
  FILE *F = fopen(filename, "wb");
  if (errno) {
    fprintf(stderr, "Can't open '%s' for writing\n%s\n", strerror(errno));
    exit(1);
  }

#ifdef COMPRESSED_HASH
  magic[8] = 'h';
#endif

#ifdef COMPRESSED_BUCKET
  magic[9] = 'b';
#endif

  fwrite(magic, sizeof(char), 16, F);

  fwrite(&_merSizeInBases, sizeof(u32bit), 1, F);
  fwrite(&_shift1, sizeof(u32bit), 1, F);
  fwrite(&_shift2, sizeof(u32bit), 1, F);
  fwrite(&_mask1, sizeof(u64bit), 1, F);
  fwrite(&_mask2, sizeof(u64bit), 1, F);

#ifdef COMPRESSED_HASH
  fwrite(&_hashWidth, sizeof(u32bit), 1, F);
#endif

#ifdef COMPRESSED_BUCKET
  fwrite(&_chckWidth, sizeof(u32bit), 1, F);
#endif

  fwrite(&_hashMask, sizeof(u64bit), 1, F);
  fwrite(&_chckMask, sizeof(u64bit), 1, F);

  fwrite(&_hashTableWords, sizeof(u64bit), 1, F);
  fwrite(&_bucketsWords,   sizeof(u64bit), 1, F);

  fwrite(_hashTable, sizeof(u64bit), _hashTableWords, F);
  fwrite(_buckets,   sizeof(u64bit), _bucketsWords,   F);

  fclose(F);

  if (errno) {
    fprintf(stderr, "existDB::saveState()-- Write failure.\n%s\n", strerror(errno));
    exit(1);
  }
}



bool
existDB::loadState(char const *filename,
                   bool        beNoisy,
                   bool        loadData) {
  char     cigam[16];

  errno = 0;
  FILE *F = fopen(filename, "rb");
  if (errno) {
    fprintf(stderr, "Can't open '%s' for reading pre-built existDB\n%s\n", strerror(errno));
    return(false);
  }

#ifdef COMPRESSED_HASH
  magic[8] = 'h';
#endif

#ifdef COMPRESSED_BUCKET
  magic[9] = 'b';
#endif

  fread(cigam, sizeof(char), 16, F);

  if (strncmp(magic, cigam, 16) != 0) {
    if (beNoisy) {
      fprintf(stderr, "existDB::loadState()-- Not an existDB binary file, maybe a sequence file?\n");
      fprintf(stderr, "existDB::loadState()-- Read     '%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c'\n",
              cigam[0],  cigam[1],  cigam[2],  cigam[3],
              cigam[4],  cigam[5],  cigam[6],  cigam[7],
              cigam[8],  cigam[9],  cigam[10], cigam[11],
              cigam[12], cigam[13], cigam[14], cigam[15]);
      fprintf(stderr, "existDB::loadState()-- Expected '%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c'\n",
              magic[0],  magic[1],  magic[2],  magic[3],
              magic[4],  magic[5],  magic[6],  magic[7],
              magic[8],  magic[9],  magic[10], magic[11],
              magic[12], magic[13], magic[14], magic[15]);
    }

    fclose(F);
    return(false);
  }

  fread(&_merSizeInBases, sizeof(u32bit), 1, F);
  fread(&_shift1, sizeof(u32bit), 1, F);
  fread(&_shift2, sizeof(u32bit), 1, F);
  fread(&_mask1, sizeof(u64bit), 1, F);
  fread(&_mask2, sizeof(u64bit), 1, F);

#ifdef COMPRESSED_HASH
  fread(&_hashWidth, sizeof(u32bit), 1, F);
#endif

#ifdef COMPRESSED_BUCKET
  fread(&_chckWidth, sizeof(u32bit), 1, F);
#endif

  fread(&_hashMask, sizeof(u64bit), 1, F);
  fread(&_chckMask, sizeof(u64bit), 1, F);

  fread(&_hashTableWords, sizeof(u64bit), 1, F);
  fread(&_bucketsWords,   sizeof(u64bit), 1, F);

  _hashTable = 0L;
  _buckets   = 0L;

  if (loadData) {
    _hashTable = new u64bit [_hashTableWords];
    _buckets   = new u64bit [_bucketsWords];

    fread(_hashTable, sizeof(u64bit), _hashTableWords, F);
    fread(_buckets,   sizeof(u64bit), _bucketsWords,   F);
  }

  fclose(F);

  if (errno) {
    fprintf(stderr, "existDB::loadState()-- Read failure.\n%s\n", strerror(errno));
    exit(1);
  }

  return(true);
}


void
existDB::printState(FILE *stream) {

  fprintf(stream, "merSizeInBases:   %u\n", _merSizeInBases);
  fprintf(stream, "tableBits         %u\n", 2 * _merSizeInBases - _shift1);
  fprintf(stream, "-----------------\n");
  fprintf(stream, "_hashTableWords   %lu (%lu KB)\n", _hashTableWords, _hashTableWords >> 7);
  fprintf(stream, "_bucketsWords     %lu (%lu KB)\n", _bucketsWords, _bucketsWords >> 7);
  fprintf(stream, "-----------------\n");
  fprintf(stream, "_shift1:          %u\n", _shift1);
  fprintf(stream, "_shift2           %u\n", _shift2);
  fprintf(stream, "_mask1            0x%016lx\n", _mask1);
  fprintf(stream, "_mask2            0x%016lx\n", _mask2);

#ifdef COMPRESSED_HASH
  fprintf(stream, "COMPRESSED_HASH   enabled\n");
  fprintf(stream, "_hashWidth        %u\n", _hashWidth);
#endif

#ifdef COMPRESSED_BUCKET
  fprintf(stream, "COMPRESSED_BUCKET enabled\n");
  fprintf(stream, "_chckWidth        %u\n", _chckWidth);
#endif

  fprintf(stream, "_hashMask         0x%016lx\n", _hashMask);
  fprintf(stream, "_chckMask         0x%016lx\n", _chckMask);
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


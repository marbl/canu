#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "existDB.H"
#include "libmeryl.H"


bool
existDB::createFromMeryl(char const  *prefix,
                         u32bit       lo,
                         u32bit       hi) {

  fprintf(stderr, "Reading mers from meryl stream %s\n", prefix);

  merStreamFromMeryl *M = new merStreamFromMeryl(prefix);

  _hashTable = 0L;
  _buckets   = 0L;

  _merSizeInBases        = M->mcd()._merSizeInBases;

  //  XXX:  Should probably be a parameter.

  u32bit tblBits = 19;

  _shift1                = 2 * _merSizeInBases - tblBits;
  _shift2                = _shift1 / 2;
  _mask1                 = u64bitMASK(tblBits);
  _mask2                 = u64bitMASK(_shift1);

#ifdef COMPRESSED_HASH
  _hashWidth             = u32bitZERO;
#endif
#ifdef COMPRESSED_BUCKET
  _chckWidth             = 2 * _merSizeInBases - tblBits;
#endif
  _hashMask              = u64bitMASK(tblBits);
  _chckMask              = u64bitMASK(2 * _merSizeInBases - tblBits);

  u64bit  tableSizeInEntries = u64bitONE << tblBits;
  u64bit  numberOfMers       = u64bitZERO;
  u64bit *countingTable      = new u64bit [tableSizeInEntries + 1];

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
      countingTable[ HASH(M->theMer()) ]++;
      numberOfMers++;
    }
  }

  delete M;
  
#ifdef COMPRESSED_HASH
  _hashWidth = 1;
  while ((numberOfMers+1) > (u64bitONE << _hashWidth))
    _hashWidth++;
#endif

  fprintf(stderr, "Found %u mers between count of %u and %u\n", numberOfMers, lo, hi);


  //  2) Allocate hash table, mer storage buckets
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



  ///////////////////////////////////////////////////////////////////////////////
  //
  //  3)  Build list of mers, placed into buckets
  //
  M = new merStreamFromMeryl(prefix);

  u64bit  h;
    
  while (M->nextMer()) {
    if ((lo <= M->theCount()) && (M->theCount() <= hi)) {
      h = HASH(M->theMer());

#ifdef COMPRESSED_BUCKET
      setDecodedValue(_buckets,
                      countingTable[h] * _chckWidth,
                      _chckWidth,
                      CHECK(M->theMer()));
#else
      _buckets[countingTable[h]] = CHECK(M->theMer());
#endif

      countingTable[h]++;
    }
  }

  delete M;
  delete [] countingTable;

  fprintf(stderr, "All done.\n");

  return(true);
}

#ifndef EXISTDB_H
#define EXISTDB_H

#include "libbri.H"

//
//  Takes as input a list of mers (in a file) and builds
//  a searchable structure listing those mers.  Doesn't
//  remove duplicate mers.
//
//  Modified from existDB (04nov2002) to insert either
//  the forward of the reverse mer only if it exists
//  in the positionDB
//

//#define COMPRESSED_HASH
//#define COMPRESSED_BUCKET

class positionDB;

class existDB {
public:
  existDB(char        *filename,
          u32bit       merSize,
          u32bit       tblBits,
          positionDB  *posDB = 0L);
  ~existDB();

  //  Inlining this is a loss on x86.
  //
  bool        exists(u64bit mer);
private:

  u64bit       HASH(u64bit k) {
    return(((k >> _shift1) ^ (k >> _shift2) ^ k) & _mask1);
  };

  u64bit       CHECK(u64bit k) {
    return(k & _mask2);
  };

  u64bit     *_hashTable;
  u64bit     *_buckets;

  u32bit      _merSizeInBases;
  //u32bit      _merSizeInBases_div2;

  u32bit      _shift1;
  u32bit      _shift2;
  u64bit      _mask1;
  u64bit      _mask2;

#ifdef COMPRESSED_HASH
  u32bit      _hashWidth;
#endif

#ifdef COMPRESSED_BUCKET
  u32bit      _chckWidth;
#endif

  u64bit      _hashMask;
  u64bit      _chckMask;
};

#endif  //  EXISTDB_H


#ifndef EXISTDB_H
#define EXISTDB_H

#include "libbri.H"

//
//  Takes as input a list of mers (in a file) and builds a searchable
//  structure listing those mers.  Doesn't remove duplicate mers.
//
//  Insert either the forward of the reverse mer only if it exists in
//  the positionDB.
//

//  If both of these are defined, the tables are half as big and half as fast.
//
//#define COMPRESSED_HASH
//#define COMPRESSED_BUCKET

#if defined COMPRESSED_HASH || defined COMPRESSED_BUCKET
#include "bit-packing.H"
#endif

class positionDB;

class existDB {
public:

  //  Read state from an existDB file
  existDB(char const  *filename,
          bool         loadData=true);

  //  Load mers from an existing existDB file, a fastafile, or a meryl database
  existDB(char const  *filename,
          u32bit       merSize,
          u32bit       tblBits,
          u32bit       lo = 0,
          u32bit       hi = ~u32bitZERO,
          positionDB  *posDB = 0L);
  ~existDB();

  void        saveState(char const *filename);

  void        printState(FILE *stream);

  //  Inlining this is a loss on x86.
  //
  bool        exists(u64bit mer);
private:
  bool        loadState(char const *filename, bool beNoisy=false, bool loadData=true);
  bool        createFromFastA(char const  *filename,
                              u32bit       merSize,
                              u32bit       tblBits,
                              positionDB  *posDB);
  bool        createFromMeryl(char const  *filename,
                              u32bit       merSize,
                              u32bit       tblBits,
                              positionDB  *posDB);

  u64bit       HASH(u64bit k) {
    return(((k >> _shift1) ^ (k >> _shift2) ^ k) & _mask1);
  };

  u64bit       CHECK(u64bit k) {
    return(k & _mask2);
  };

  u32bit      _merSizeInBases;

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

  u64bit      _hashTableWords;
  u64bit     *_hashTable;
  u64bit      _bucketsWords;
  u64bit     *_buckets;
};

#endif  //  EXISTDB_H

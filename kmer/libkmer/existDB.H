#ifndef EXISTDB_H
#define EXISTDB_H

#include "bio++.H"

//
//  Takes as input a list of mers (in a file) and builds a searchable
//  structure listing those mers.  Doesn't remove duplicate mers.
//
//  Insert either the forward of the reverse mer only if it exists in
//  the positionDB.
//

//  Extend to allow runtime selection of COMPRESSED_HASH and COMPRESSED_BUCKET.
//  Extend to allow runtime enable of VERBOSE (currently always not verbose)
//  Extend to allow optimal memory sizing ala meryl (only if we are given meryl as input)

//  Using a compressed hash is extreme -- really slow and doesn't save that much.

//  Print some statistics
//
//#define STATS


typedef u32bit existDBflags;
const existDBflags  existDBcompressHash    = 0x0001;
const existDBflags  existDBcompressBuckets = 0x0002;
const existDBflags  existDBverbose         = 0x0004;
const existDBflags  existDBcanonical       = 0x0008;
const existDBflags  existDBforward         = 0x0010;
const existDBflags  existDBreverse         = 0x0020;

class existDB {
public:

  //  Read state from an existDB file
  existDB(char const  *filename,
          bool         loadData=true);

  //  Load mers from an existing existDB file, a fastafile, or a meryl database
  existDB(char const    *filename,
          u32bit         merSize,
          u32bit         tblBits,
          u32bit         lo = 0,
          u32bit         hi = ~u32bitZERO,
          existDBflags   flags = 0x0000);
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
                              u32bit       tblBits);
  bool        createFromMeryl(char const  *filename,
                              u32bit       lo,
                              u32bit       hi,
                              u32bit       tblBits);

  u64bit       HASH(u64bit k) {
    return(((k >> _shift1) ^ (k >> _shift2) ^ k) & _mask1);
  };

  u64bit       CHECK(u64bit k) {
    return(k & _mask2);
  };

  bool        _compressedHash;
  bool        _compressedBucket;
  bool        _beVerbose;

  u32bit      _merSizeInBases;

  u32bit      _shift1;
  u32bit      _shift2;
  u64bit      _mask1;
  u64bit      _mask2;

  u32bit      _hashWidth;  //  Only for the compressed hash
  u32bit      _chckWidth;  //  Only for the compressed bucket

  u64bit      _hashMask;
  u64bit      _chckMask;

  u64bit      _hashTableWords;
  u64bit      _bucketsWords;
  u64bit     *_hashTable;
  u64bit     *_buckets;
};

#endif  //  EXISTDB_H

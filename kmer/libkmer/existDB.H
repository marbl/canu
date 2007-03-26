#ifndef EXISTDB_H
#define EXISTDB_H

#include "bio++.H"

//  Takes as input a list of mers (in a file) and builds a searchable
//  structure listing those mers.  Duplicate mers are not removed and
//  will be stored multiple times.
//
//  Using a compressed hash is allowed, but somewhat useless -- it is
//  really slow and doesn't save that much.
//
//  If existDBcanonical is requested, this will store only the
//  canonical mer.  It is up to the client to be sure that is
//  appropriate!  See positionDB.H for more.

//#define STATS

typedef u32bit existDBflags;
const existDBflags  existDBcompressHash    = 0x0001;
const existDBflags  existDBcompressBuckets = 0x0002;
const existDBflags  existDBcanonical       = 0x0008;
const existDBflags  existDBforward         = 0x0010;

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

  bool        isForward(void)    { return(_isForward);   };
  bool        isCanonical(void)  { return(_isCanonical); };

  bool        exists(u64bit mer);
  bool        exists(u64bit b, u64bit c);

private:
  bool        loadState(char const *filename, bool beNoisy=false, bool loadData=true);
  bool        createFromFastA(char const  *filename,
                              u32bit       merSize,
                              u32bit       tblBits,
                              u32bit       flags);
  bool        createFromMeryl(char const  *filename,
                              u32bit       lo,
                              u32bit       hi,
                              u32bit       tblBits,
                              u32bit       flags);

  u64bit       HASH(u64bit k) {
    return(((k >> _shift1) ^ (k >> _shift2) ^ k) & _mask1);
  };

  u64bit       CHECK(u64bit k) {
    return(k & _mask2);
  };

  void         insertMer(u64bit h, u64bit c, u64bit *countingTable) {
    if (_compressedBucket)
      setDecodedValue(_buckets, countingTable[h] * _chckWidth, _chckWidth, c);
    else
      _buckets[countingTable[h]] = c;
    countingTable[h]++;
  };

  bool        _compressedHash;
  bool        _compressedBucket;
  bool        _isForward;
  bool        _isCanonical;

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

  //  Used by the streaming exists interface
#if 0
  u64bit      _es_bucket;
  u64bit      _es_st;
  u64bit      _es_ed;
#endif

  void clear(void) {
#if 0
    _es_bucket = ~u64bitZERO;
    _es_st     = 0;
    _es_ed     = 0;
#endif
  };
};

#endif  //  EXISTDB_H

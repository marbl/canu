#ifndef POSITIONDB_H
#define POSITIONDB_H

#include "libbri.H"

typedef u64bit heapbit;

class positionDB {
public:
  positionDB(unsigned char const  *seq,
             char const           *filename,
             u32bit                merSize,
             u32bit                merSkip,    //  mers to skip between mers
             u32bit                tblBits,
             bool                  beVerbose);
  ~positionDB();

  void        dump(u64bit mer);
  void        dumpTable(void);

  //  Returns position in posn, resizing it if needed.  Space is
  //  allocated if none supplied.  The following is valid:
  //
  //    u64bit  *posn    = 0L;
  //    u64bit   posnMax = 0;
  //    u64bit   posnLen = 0;
  //    if (get(somemer, posn, posnMax, posnLen)) {
  //      do something with the positions
  //    }
  //
  bool        get(u64bit   mer,
                  u64bit*& posn,
                  u64bit&  posnMax,
                  u64bit&  posnLen);

  //  True/False on existence of mer.
  //
  bool        exists(u64bit   mer);
private:

  u64bit       HASH(u64bit k) {
    return(((k >> _shift1) ^ (k >> _shift2) ^ k) & _mask1);
  };

  u64bit       CHECK(u64bit k) {
    return(k & _mask2);
  };

  void         sortAndRepackBucket(u64bit b);

  u32bit     *_bucketSizes;
  u64bit     *_countingBuckets;
  u64bit     *_hashTable;
  u64bit     *_buckets;
  u64bit     *_positions;

  u32bit      _merSizeInBases;
  u32bit      _merSizeInBits;

  u32bit      _merSkipInBases;

  u32bit      _tableSizeInBits;
  u64bit      _tableSizeInEntries;

  u32bit      _hashWidth;
  u64bit      _hashMask;

  u32bit      _chckWidth;
  u64bit      _chckMask;

  u32bit      _posnWidth;
  u64bit      _posnMask;

  u32bit      _wCnt;
  u32bit      _wFin;

  u32bit      _shift1;
  u32bit      _shift2;
  u64bit      _mask1;
  u64bit      _mask2;

  u64bit      _numberOfMers;
  u64bit      _numberOfPositions;
  u64bit      _numberOfDistinct;
  u64bit      _numberOfUnique;
  u64bit      _numberOfEntries;
  u64bit      _maximumEntries;

  //  For sorting the mers
  //
  u32bit   _sortedListMax;
  u32bit   _sortedListLen;
  heapbit *_sortedList;
};

#endif  //  POSITIONDB_H

#ifndef MCBUCKET_H
#define MCBUCKET_H

#include <stdio.h>
#include "mcDescription.H"
#include "libbri.H"



//  A mer that we wish to get the count for.
//
//  mer:       the bit-packed representation of the mer.
//  sequence:  the sequence ID the mer is from.
//  position:  the location in the sequence of the mer.
//  count:     the count that we are trying to find.
//
class mcMer {
public:
  u64bit    mer;
  u32bit    sequence;
  u32bit    position;
  u32bit    count;
};


//  A single bucket read from the mc files.
//
class mcBucket {
private:
public:
  bitPackedFileReader   *_IDX;
  bitPackedFileReader   *_DAT;

  u32bit                 _chckBits;
  u64bit                 _chckMask;
  u64bit                 _firstBit;

  u64bit                 _bucketID;
  u64bit                 _lastBucketID;

  u64bit                 _items;
  u64bit                 _itemsMax;
  u64bit                *_checks;
  u64bit                *_counts;

  u64bit                 _bitsRead;

public:
  mcBucket(bitPackedFileReader *IDX, bitPackedFileReader *DAT, mcDescription *mcd) {
    _IDX = IDX;
    _DAT = DAT;

    _chckBits     = mcd->_chckBits;
    _chckMask     = u64bitMASK(_chckBits);
    _firstBit     = u64bitONE << _chckBits;

    _bucketID     = 0;
    _lastBucketID = mcd->_tableSizeInEntries;

    //  Allocate some space
    //
    _items    = 0;
    _itemsMax = 131072;
    _checks   = new u64bit [_itemsMax];
    _counts   = new u64bit [_itemsMax];

    _bitsRead = 0;

    //  Read the first bucket
    //
    readBucket();
    _bucketID = 0;
  };
  ~mcBucket() {
    delete [] _checks;
    delete [] _counts;
  };

  bool     readBucket(void);

  void     read(mcMer *mer);
  void     scan(mcMer *mer);
  void     dump(FILE *F);

  u64bit   bitsRead(void) {
    return(_bitsRead);
  }

  void     reallocate(void) {
    if (_items >= _itemsMax) {
      _itemsMax <<= 1;

      u64bit *ch = new u64bit [_itemsMax];
      u64bit *co = new u64bit [_itemsMax];

      for (u64bit i=_items; i--; ) {
        ch[i] = _checks[i];
        co[i] = _counts[i];
      }

      delete [] _checks;
      delete [] _counts;

      _checks = ch;
      _counts = co;
    }
  };
};


#endif  //  MCBUCKET_H

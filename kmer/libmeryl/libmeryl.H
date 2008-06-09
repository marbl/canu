#ifndef LIBMERYL_H
#define LIBMERYL_H

#include "bio++.H"

//  A merStream reader/writer for meryl mercount data.
//
//  merSize is used to check that the meryl file is the correct size.
//  If it isn't the code fails.
//
//  The reader returns mers in lexicographic order.  No random access.
//  The writer assumes that mers come in sorted increasingly.
//
//  numUnique    the total number of mers with count of one
//  numDistinct  the total number of distinct mers in this file
//  numTotal     the total number of mers in this file


class merylStreamReader {
public:
  merylStreamReader(const char *fn, u32bit ms=0);
  ~merylStreamReader();

  kMer           &theFMer(void)      { return(_thisMer);          };
  u64bit          theCount(void)     { return(_thisMerCount);     };

  bool            hasPositions(void)    { return(_POS != 0L); };
  u32bit         *thePositions(void)    { return(_thisMerPositions); };
  u32bit          getPosition(u32bit i) { return(((_POS) && (i < _thisMerCount)) ? _thisMerPositions[i] : ~u32bitZERO); };

  u32bit          merSize(void)         { return(_merSizeInBits >> 1); };
  u32bit          merCompression(void)  { return(_merCompression); };

  u32bit          prefixSize(void) { return(_prefixSize); };

  u64bit          numberOfUniqueMers(void)   { return(_numUnique); };
  u64bit          numberOfDistinctMers(void) { return(_numDistinct); };
  u64bit          numberOfTotalMers(void)    { return(_numTotal); };

  u64bit          histogram(u32bit i)         { return((i < _histogramLen) ? _histogram[i] : ~u64bitZERO); };
  u64bit          histogramLength(void)       { return(_histogramLen); };
  u64bit          histogramHuge(void)         { return(_histogramHuge); };
  u64bit          histogramMaximumCount(void) { return(_histogramMaxValue); };

  bool            nextMer(void);
  bool            validMer(void) { return(_validMer); };
private:
  bitPackedFile         *_IDX;
  bitPackedFile         *_DAT;
  bitPackedFile         *_POS;

  u64bit                  getIDXnumber(void) {
    u64bit n = 1;

    if (_idxIsPacked)
      n = _IDX->getNumber();
    else
      n = _IDX->getBits(32);

    return(n);
  };
  u64bit                  getDATnumber(void) {
    u64bit n = 1;

    if (_datIsPacked) {
      if (_DAT->getBits(1))
        n = _DAT->getNumber() + 2;
    } else {
      n = _DAT->getBits(32);
    }

    return(n);
  };

  //  Why not bool?  Seems like the bitPackedFile is incompatible
  //  with bools.
  u32bit                 _idxIsPacked;
  u32bit                 _datIsPacked;
  u32bit                 _posIsPacked;

  u32bit                 _merSizeInBits;
  u32bit                 _merCompression;
  u32bit                 _prefixSize;
  u32bit                 _merDataSize;
  u64bit                 _thisBucket;
  u64bit                 _thisBucketSize;
  u64bit                 _numBuckets;

  kMer                   _thisMer;
  u64bit                 _thisMerCount;

  u32bit                 _thisMerPositionsMax;
  u32bit                *_thisMerPositions;

  u64bit                 _numUnique;
  u64bit                 _numDistinct;
  u64bit                 _numTotal;

  u64bit                 _histogramHuge;       // number that are bigger than Len
  u64bit                 _histogramLen;        // number of entries in the histo
  u64bit                 _histogramMaxValue;   // highest count ever seen
  u64bit                *_histogram;

  bool                   _validMer;
};


class merylStreamWriter {
public:
  merylStreamWriter(const char *filePrefix,
                    u32bit merSize,          //  In bases
                    u32bit merComp,          //  A length, bases
                    u32bit prefixSize,       //  In bits
                    bool   positionsEnabled);
  ~merylStreamWriter();

  void                    addMer(kMer &mer, u32bit count=1, u32bit *positions=0L);
private:
  void                    writeMer(void);

  void                    setIDXnumber(u64bit n) {
    if (_idxIsPacked)
      _IDX->putNumber(n);
    else
      _IDX->putBits(n, 32);
  };
  void                    setDATnumber(u64bit n) {
    if (_datIsPacked) {
      if (n == 1) {
        _DAT->putBits(u64bitZERO, 1);
      } else {
        _DAT->putBits(u64bitONE, 1);
        _DAT->putNumber(n-2);
      }
    } else {
      _DAT->putBits(n, 32);
    }
  };


  bitPackedFile         *_IDX;
  bitPackedFile         *_DAT;
  bitPackedFile         *_POS;

  u32bit                 _idxIsPacked;
  u32bit                 _datIsPacked;
  u32bit                 _posIsPacked;

  u32bit                 _merSizeInBits;
  u32bit                 _merCompression;
  u32bit                 _prefixSize;
  u32bit                 _merDataSize;
  u64bit                 _thisBucket;
  u64bit                 _thisBucketSize;
  u64bit                 _numBuckets;

  u64bit                 _numUnique;
  u64bit                 _numDistinct;
  u64bit                 _numTotal;

  u64bit                 _histogramHuge;       // number that are bigger than Len
  u64bit                 _histogramLen;        // number of entries in the histo
  u64bit                 _histogramMaxValue;   // highest count ever seen
  u64bit                *_histogram;

  kMer                   _thisMer;
  u64bit                 _thisMerCount;
};

#endif  //  LIBMERYL_H

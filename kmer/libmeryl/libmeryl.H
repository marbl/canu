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

#define LIBMERYL_HISTOGRAM_MAX  1048576

class merylStreamReader {
public:
  merylStreamReader(const char *fn, u32bit ms=0);
  ~merylStreamReader();

  kMer           &theFMer(void)  { return(_thisMer); };

  u64bit          theCount(void) { return(_thisMerCount); };

  u32bit          getPosition(u32bit i);

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


inline
u32bit
merylStreamReader::getPosition(u32bit i) {
  if ((i < _thisMerCount) && (_POS))
    return(_thisMerPositions[i]);
  return(~u32bitZERO);
}


inline
bool
merylStreamReader::nextMer(void) {

  //  Use a while here, so that we skip buckets that are empty
  //
  while ((_thisBucketSize == 0) && (_thisBucket < _numBuckets)) {
    _thisBucketSize = _IDX->getNumber();
    _thisBucket++;
  }

  if (_thisBucket >= _numBuckets)
    return(_validMer = false);

  //  Before you get rid of the clear() -- if, say, the list of mers
  //  is sorted and we can shift the mer to make space for the new
  //  stuff -- make sure that nobody is calling reverseComplement()!
  //
  _thisMer.clear();
  _thisMer.readFromBitPackedFile(_DAT, _merDataSize);
  _thisMer.setBits(_merDataSize, _prefixSize, _thisBucket);

  _thisMerCount = 1;

  if (_DAT->getBits(1))
    _thisMerCount = _DAT->getNumber() + 2;

  _thisBucketSize--;

  if (_POS) {
    if (_thisMerPositionsMax < _thisMerCount) {
      delete [] _thisMerPositions;
      _thisMerPositionsMax = _thisMerCount + 1024;
      _thisMerPositions    = new u32bit [_thisMerPositionsMax];
    }
    for (u32bit i=0; i<_thisMerCount; i++) {
      _thisMerPositions[i] = _POS->getBits(32);
    }
  }

  return(true);
}


class merylStreamWriter {
public:
  merylStreamWriter(const char *filePrefix,
                    u32bit merSize,          //  In bases
                    u32bit merComp,          //  A length, bases
                    u32bit prefixSize);      //  In bits

  ~merylStreamWriter();

  void                    addMer(kMer &mer, u32bit count=1, u32bit position=~u32bitZERO);
private:
  void                    writeMer(void);

  bitPackedFile         *_IDX;
  bitPackedFile         *_DAT;
  bitPackedFile         *_POS;

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

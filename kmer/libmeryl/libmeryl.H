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
//

class merylStreamReader {
public:
  merylStreamReader(const char *fn, u32bit ms=0);
  ~merylStreamReader();

  //  theRMer() is used in existDB-create-from-meryl.C:61, 62, 166
  kMer           &theFMer(void)  { return(_thisMer); };

  u64bit          theCount(void) { return(_thisMerCount); };

  u32bit          merSize(void)    { return(_merSizeInBits >> 1); };
  u32bit          prefixSize(void) { return(_prefixSize); };

  u64bit          numberOfUniqueMers(void)   { return(_numUnique); };
  u64bit          numberOfDistinctMers(void) { return(_numDistinct); };
  u64bit          numberOfTotalMers(void)    { return(_numTotal); };

  bool            nextMer(void);
  bool            validMer(void) { return(_validMer); };
private:
  bitPackedFile         *_IDX;
  bitPackedFile         *_DAT;

  u32bit                 _merSizeInBits;
  u32bit                 _prefixSize;
  u32bit                 _merDataSize;
  u64bit                 _thisBucket;
  u64bit                 _thisBucketSize;
  u64bit                 _numBuckets;

  kMer                   _thisMer;
  u64bit                 _thisMerCount;

  u64bit                 _numUnique;
  u64bit                 _numDistinct;
  u64bit                 _numTotal;

  bool                   _validMer;
};


inline
bool
merylStreamReader::nextMer(void) {

  //  Use a while here, so that we skip buckets that are empty
  //
  while ((_thisBucketSize == 0) && (_thisBucket < _numBuckets)) {
    _thisBucketSize = _IDX->getNumber();
    _thisBucket++;
  }

#if 0
  fprintf(stderr, "thisBucket: "u64bitFMT" (max "u64bitFMT") size "u64bitFMT"\n",
          _thisBucket, _numBuckets, _thisBucketSize);
#endif

  if (_thisBucket >= _numBuckets)
    return(_validMer = false);

  _thisMer.clear();
  _thisMer.readFromBitPackedFile(_DAT, _merDataSize);
  _thisMer.setBits(_merDataSize, _prefixSize, _thisBucket);

  _thisMerCount = 1;

  if (_DAT->getBits(1))
    _thisMerCount = _DAT->getNumber() + 2;

  _thisBucketSize--;

  return(true);
}







class merylStreamWriter {
public:
  merylStreamWriter(const char *filePrefix,
                    u32bit merSize,          //  In bases
                    u32bit prefixSize);      //  In bits

  ~merylStreamWriter();

  void                    addMer(kMer &mer, u32bit count=1);
private:
  void                    writeMer(void);

  bitPackedFile         *_IDX;
  bitPackedFile         *_DAT;

  u32bit                 _merSizeInBits;
  u32bit                 _prefixSize;
  u32bit                 _merDataSize;
  u64bit                 _thisBucket;
  u64bit                 _thisBucketSize;
  u64bit                 _numBuckets;

  u64bit                 _numUnique;
  u64bit                 _numDistinct;
  u64bit                 _numTotal;

  kMer                   _thisMer;
  u64bit                 _thisMerCount;
};

#endif  //  LIBMERYL_H

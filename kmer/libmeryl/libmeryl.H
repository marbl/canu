#ifndef LIBMERYL_H
#define LIBMERYL_H

#include "libbri.H"

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


  u64bit          theFMer(void)  { return(_thisMer); };
  u64bit          theRMer(void)  { return(reverseComplementMer(_merSizeInBits >> 1, _thisMer)); };

  u64bit          theCount(void) { return(_thisMerCount); };

  u32bit          merSize(void)    { return(_merSizeInBits >> 1); };
  u32bit          prefixSize(void) { return(_prefixSize); };

  u64bit          numberOfUniqueMers(void)   { return(_numUnique); };
  u64bit          numberOfDistinctMers(void) { return(_numDistinct); };
  u64bit          numberOfTotalMers(void)    { return(_numTotal); };

  bool            nextMer(void);
  bool            validMer(void) { return(_validMer); };
private:
  bitPackedFileReader   *_IDX;
  bitPackedFileReader   *_DAT;

  u32bit                 _merSizeInBits;
  u32bit                 _prefixSize;
  u64bit                 _prefixMask;
  u32bit                 _merDataSize;
  u64bit                 _merDataMask;

  u64bit                 _thisBucket;
  u64bit                 _thisBucketSize;

  u64bit                 _thisMer;
  u64bit                 _thisMerCount;

  u64bit                 _numUnique;
  u64bit                 _numDistinct;
  u64bit                 _numTotal;

  bool                   _validMer;
};



class merylStreamWriter {
public:
  merylStreamWriter(const char *filePrefix,
                    u32bit merSize,          //  In bases
                    u32bit prefixSize,       //  In bits
                    u64bit numUnique=0,
                    u64bit numDistinct=0,
                    u64bit numTotal=0);
  ~merylStreamWriter();

  void                    addMer(u64bit mer, u32bit count=1);
private:
  bitPackedFileWriter   *_IDX;
  bitPackedFileWriter   *_DAT;

  u32bit                 _merSizeInBits;
  u32bit                 _prefixSize;
  u64bit                 _prefixMask;
  u32bit                 _merDataSize;
  u64bit                 _merDataMask;

  u64bit                 _thisBucket;
  u64bit                 _thisBucketSize;

  u64bit                 _thisMer;
  u64bit                 _thisMerCount;
};

#endif  //  LIBMERYL_H

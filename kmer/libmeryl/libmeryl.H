#ifndef LIBMERYL_H
#define LIBMERYL_H

#include "libbri.H"
#include "mcBucket.H"
#include "mcDescription.H"

//
//  Implements a merStream that gets its mers from a meryl database
//
//  merSize is used to check that the meryl file is the correct size.
//  If it isn't the code fails.
//

class merStreamFromMeryl {
public:
  merStreamFromMeryl(const char *prefix, u32bit merSize=0);
  ~merStreamFromMeryl();

  u64bit          theMer(void)         { return(_mer); };
  u64bit          theCount(void)       { return(_cnt); };

  u32bit          merSize(void)        { return(_mcd._merSizeInBases); };

  mcDescription  &mcd(void)            { return(_mcd); };

  bool            nextMer(void);
private:
  mcDescription         _mcd;
  bitPackedFileReader  *_IDX;
  bitPackedFileReader  *_DAT;

  mcBucket             *_B;
  u32bit                _tPos;  //  Position in the table of buckets
  u32bit                _bPos;  //  Position in the bucket

  u64bit                _mer;
  u64bit                _cnt;
};

inline
bool
merStreamFromMeryl::nextMer(void) {

  //  Use a while here, so that we skip buckets that are empty
  //
  while ((_bPos >= _B->_items) &&
         (_tPos <  _mcd._tableSizeInEntries)) {
    _tPos++;
    _bPos = 0;
    _B->readBucket();
  }

  if (_tPos >= _mcd._tableSizeInEntries)
    return(false);

  _mer = _B->_bucketID << _B->_chckBits | _B->_checks[_bPos];
  _cnt = _B->_counts[_bPos];

  _bPos++;

  return(true);
}





//  this class is tricky, because it needs to have the mers fed into it in
//  a sorted order.

#if 0

class merylFileWriter {
public:
  merylFileWriter(const char *prefix, mcDescription *mcd);
  ~merylFileWriter();

  put(u64bit mer, u64bit cnt);
private:
};

#endif



#endif  //  LIBMERYL_H

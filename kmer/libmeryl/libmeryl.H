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

//  Define this to test the Fmer -> Rmer computation.
//#define TEST_RMER

class merStreamFromMeryl {
public:
  merStreamFromMeryl(const char *prefix, u32bit merSize=0);
  ~merStreamFromMeryl();

  //  Dang.  Meryl databases only have the forward mer, but we need to
  //  examine both here.  So, we painfully construct the reverse mer.
  //
  //  XXX:  This abuses the fact that a=00, c=01, g=10, t=11, and
  //  so we can reverse complement with xor.

  u64bit          theFMer(void) {
    return(_mer);
  };
  u64bit          theRMer(void) {
    u64bit  fmer = _mer;
    u64bit  rmer = u64bitZERO;

    //  Complement everything
    fmer ^= ~u64bitZERO;

    //  Reverse everything -- we let the compiler unroll this, at expense
    //  of maybe 50% more iterations
    //
    for (u32bit i=0; i<32; i++) {
      rmer <<= 2;
      rmer  |= fmer & 0x3LLU;
      fmer >>= 2;
    }

#ifdef TEST_RMER
    //  We can easily test that this works by just doing it again!
    //
    u64bit tr = rmer;
    u64bit tf = u64bitZERO;

    tr ^= ~u64bitZERO;

    //  Reverse everything -- we let the compiler unroll this, at expense
    //  of maybe 50% more iterations
    //
    for (u32bit i=0; i<32; i++) {
      tf <<= 2;
      tf  |= tr & 0x3LLU;
      tr >>= 2;
    }

    if (tf != _mer)
      fprintf(stderr, "ERROR -- reverse complement didn't work for fmer=0x%016lx (got 0x%016lx)!\n", _mer, tf);
#endif

    return(rmer);
  };

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

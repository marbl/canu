#ifndef MERSTREAM_H
#define MERSTREAM_H

#include "libbritypes.h"
#include "fastastream.H"

//
//  This is the string-based merStreamer.  It gets sequence from
//  a character string.
//

class merStream {
public:
  merStream(u32bit merSize, const char *filename);
  merStream(u32bit merSize, char const *seq, u32bit len);
  ~merStream();

  u64bit       theFMer(void)        { return(_theFMer); };
  u64bit       theRMer(void)        { return(_theRMer); };

  char const  *theFMerString(void);
  char const  *theRMerString(void);

  u64bit       thePosition(void) {
    if (_theFile)
      return(_theFile->thePosition() - _merSize);
    else
      return(_thePos - _merSize);
  };
  u64bit       theSequenceNumber(void) {
    if (_theFile)
      return(_theFile->theSequenceNumber());
    else
      return(u64bitZERO);
  };
  char const  *theDefLine(void) {
    if (_theFile)
      return(_theFile->theDefLine());
    else
      return(">");
  };

  bool               nextMer(u32bit skip=0);
private:
  void               loadMer(u32bit s);
  char               nextSymbol(void) {
    if (_theFile)
      return(_theFile->nextSymbol());
    else
      return(_theSeq[_thePos++]);
  };

  FastAstream          *_theFile;
  char const           *_theSeq;
  u32bit                _theLen;
  u32bit                _thePos;

  u32bit                _merSize;
  s32bit                _timeUntilValid;
  u64bit                _theMerMask;

  u64bit                _theFMer;
  u64bit                _theRMer;

  u32bit                _theRMerShift;

  char                  _theMerString[33];
};




//  We have a problem; at the start of a stream, we want to
//  initialize the mer with merSize-1 bases.  In the middle of the
//  stream, we need to load the mer with merSize bases.
//
//  loadMer will push s bases onto the mer, restarting if it hits a
//  mer break.
//
//  No masking of the mer is performed.
//
inline
void
merStream::loadMer(u32bit s) {
  u64bit   ch = 255;
  u64bit   cf = 0;
  u64bit   cr = 0;

  _timeUntilValid = s;

  //  While we are invalid, and still in the sequence
  //  push characters onto the mer.  The valid
  //  count is updated if we hit an invalid base.
  //
  while ((_timeUntilValid > 0) && (ch != 0)) {
    ch = nextSymbol();
    cf = validCompressedSymbol[ch];

    //  Rather than take the chance of generating a cache miss accessing
    //  another array, we just mask out the upper bits of validCompressedSymbol[];
    //  this is exactly the same as using compressSymbol[].
    //
    //  We need to mask the upper bits for reverse (but not for forward)
    //  because, in reverse, we shift to the right.  If we don't mask these
    //  out, we will have extra bits set in the mer.
    //
    //  Example: Consider placing 255 (== 11111111, the invalid symbol returned
    //  from validCompressedSymbol) at the fourth base in a four mer:
    //    000000**xxxxxx -- the ** are is the fourth base.
    //
    //  Without masking, we'd set all the 0's and all the *'s to one, in effect,
    //  preloading the next three bases of the mer.
    //
    cr = validCompressedSymbol[complementSymbol[ch]] & 0x03;

    _timeUntilValid--;

    if (cf & 0xfc)
      _timeUntilValid = s;

    //  If the ch is valid, we can obviously add it to the mer.
    //  If it is invalid, we don't care; by the time the mer is valid
    //  all bits of any invalid mers are removed.  Sure, we could
    //  put this into the else, but I suspect that it will be faster outside.
    //
    _theFMer <<= 2;
    _theRMer >>= 2;

    _theFMer |= cf;
    _theRMer |= cr << _theRMerShift;
  }
}

inline
bool
merStream::nextMer(u32bit skip) {

 again:
  u64bit  ch = nextSymbol();
  u64bit  cf = validCompressedSymbol[ch];
  u64bit  cr = validCompressedSymbol[complementSymbol[ch]] & 0x03;

  //  EOF?
  if (ch == 0)
    return(false);

  //  Push the ch onto the mer.
  //
  _theFMer <<= 2;
  _theRMer >>= 2;

  _theFMer |= cf;
  _theRMer |= cr << _theRMerShift;

  //  If the ch is invalid, we need to make a whole new mer.
  //  This also implies that we're done skipping.
  //
  if (cf & (unsigned char)0xfc) {
    loadMer(_merSize);
    skip = 0;
  }

  if (skip--)
    goto again;

  _theFMer &= _theMerMask;
  _theRMer &= _theMerMask;

  //  Still need to check if we are valid -- we could
  //  have run off the end of the sequence before a valid
  //  mer was created.
  //
  return(_timeUntilValid <= 0);
}

#endif  //  MERSTREAM_H

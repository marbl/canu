#ifndef FASTA_ACCESSOR_H
#define FASTA_ACCESSOR_H

#include "bio++.H"


//  Define this to do bounds checking
//
#if 1
#define SANITY(NAM, POS) \
  if ((POS) > _len) { \
    fprintf(stderr, "%s-- position "uint32FMT" larger than length "uint32FMT"\n", \
            NAM, (POS), _len); \
    assert((POS) <= _len); \
  }
#else
#define SANITY(NAM, POS)
#endif


//  Externally, we show the coordinate in the forward strand.
//  Internally, we represent the coordinate on the physical string.
//  The user is required to maintain the range that the
//  reverse-complement string is valid over.  It is not possible to
//  randomly access sequence outside the range, but it is possible to
//  iterate over it (but then you cannot get the coordinate of where
//  you are at!)


class FastAAccessor {
private:
public:
  char    *_seq;
  uint32   _pos;
  uint32   _len;

  uint32   _rcBase;
  uint32   _rcLen;

  bool     _doRevComp;

private:
  void FastAAccessorInit(char *S, uint32 length, bool revcomp) {
    _seq = S;
    _pos = 0;
    _len = length;
    if (length == 0)
      _len = (uint32)strlen(S);

    _rcBase = 0;
    _rcLen  = _len;

    _doRevComp = revcomp;

    if (_doRevComp)
      _pos = _len-1;
  };

public:
  FastAAccessor(seqInCore *S, bool revcomp=false) {
    FastAAccessorInit(S->sequence(), S->sequenceLength(), revcomp);
  };

  FastAAccessor(char *S, uint32 length=0, bool revcomp=false) {
    FastAAccessorInit(S, length, revcomp);
  };

private:

  //  Given a range in the forward string, we can reverse-complement
  //  just that range.  This amounts to translating the forward string
  //  to make the beginning of the range be the origin, then
  //  reverseing the range, then translating the sequence back to the
  //  original origin.
  //
  uint32  rc(uint32 p) const {
    return(_rcBase + _rcLen - (p - _rcBase) - 1);
  };

public:

  //  For iterating over reverse complement regions of a forward
  //  sequence.
  //
  //  e.g.,  (500, 250) would be:
  //
  //  |-----|-----|XXXXX|-----|
  //  0    250   500   750 1000
  //
  //  Set both to zero (also the default) to unset the range
  //
  //  The physical location (_pos) doesn't change, but this will
  //  change the value of the corresponding forward coordinate, but
  //  not the meaning.
  //
  bool  setRange(uint32 base=0, uint32 length=0) {
    if ((base < _len) && (base + length <= _len)) {
      if ((base == 0) && (length == 0)) {
        _rcBase = 0;
        _rcLen  = _len;
      } else {
        _rcBase = base;
        _rcLen  = length;
      }
      return(true);
    } else {
      fprintf(stderr, "FastAAccessor::setRange()-- base="uint32FMT" and length="uint32FMT" exceed sequence length of "uint32FMT"\n",
              base, length, _len);
      assert(0);
      return(false);
    }
  };


  //  True if this physical location is valid.
  //
  bool isValid(void) {
    return(_pos < _len);
  };


  char operator[](uint32 p) const {
    if ((_rcBase <= p) && (p < _rcBase + _rcLen)) {
      if (_doRevComp) return(complementSymbol[_seq[rc(p)]]);
      else            return(_seq[p]);
    } else {
      fprintf(stderr, "operator[]-- Tried to access to "uint32FMT", but range is "uint32FMT"-"uint32FMT"\n",
              p, _rcBase, _rcBase + _rcLen);
      assert(0);
      return(0);
    }
  };



  //  Set the accessor to some position.
  //
  bool setPosition(uint32 p) {
    if ((_rcBase <= p) && (p < _rcBase + _rcLen)) {
      if (_doRevComp) _pos = rc(p);
      else            _pos = p;
      return(true);
    } else {
      fprintf(stderr, "setPosition()-- Tried to set to "uint32FMT", but range is "uint32FMT"-"uint32FMT".\n",
              p, _rcBase, _rcBase + _rcLen);
      assert(0);
      return(false);
    }
  };

  uint32 getPosition(void) {
    if (_doRevComp) return(rc(_pos));
    else            return(_pos);
  };

  uint32 getRangeBegin(void)  { return(_rcBase); };
  uint32 getRangeEnd(void)    { return(_rcBase + _rcLen); };
  uint32 getRangeLength(void) { return(_rcLen); };

  bool   extendLeft(int32 x) {
    _rcLen  += x;
    if (_doRevComp == false)
      _rcBase -= x;

    if ((_rcBase > _len) || (_rcBase + _rcLen > _len)) {
      fprintf(stderr, "FastAAccessor::extendLeft()-- extend by "int32FMT" makes invalid: length is "uint32FMT", new range is "uint32FMT"-"uint32FMT"\n",
              x, _len, _rcBase, _rcBase + _rcLen);
      assert(0);
      return(false);
    }
    return(true);
  };

  bool   extendRight(int32 x)   {
    _rcLen  += x;
    if (_doRevComp == true)
      _rcBase -= x;

    if ((_rcBase > _len) || (_rcBase + _rcLen > _len)) {
      fprintf(stderr, "FastAAccessor::extendRight()-- extend by "int32FMT" makes invalid: length is "uint32FMT", new range is "uint32FMT"-"uint32FMT"\n",
              x, _len, _rcBase, _rcBase + _rcLen);
      assert(0);
      return(false);
    }
    return(true);
  };

  char operator*(void) const {
    SANITY("FastAAccessor::operator*()", _pos);
    if (_doRevComp) return(complementSymbol[_seq[_pos]]);
    else            return(_seq[_pos]);
  };

  char get(void) const {
    SANITY("FastAAccessor::get()", _pos);
    if (_doRevComp) return(complementSymbol[_seq[_pos]]);
    else            return(_seq[_pos]);
  };

  FastAAccessor &operator--(void) {
    if (_doRevComp) _pos++;
    else            _pos--;
    return(*this);
  };

  FastAAccessor &operator++(void) {
    if (_doRevComp) _pos--;
    else            _pos++;
    return(*this);
  };
};

#endif  //  FASTA_ACCESSOR_H

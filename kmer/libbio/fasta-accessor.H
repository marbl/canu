#ifndef FASTA_ACCESSOR_H
#define FASTA_ACCESSOR_H

#include "bio++.H"


//  Define this to do bounds checking
//
#if 1
#define SANITY(NAM, POS) \
  if ((POS) > _len) { \
    fprintf(stderr, "%s-- position "u32bitFMT" larger than length "u32bitFMT"\n", \
            NAM, (POS), _len); \
    abort(); \
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
  u32bit   _pos;
  u32bit   _len;

  u32bit   _rcBase;
  u32bit   _rcLen;

  bool     _doRevComp;

private:
  void FastAAccessorInit(char *S, u32bit length, bool revcomp) {
    _seq = S;
    _pos = 0;
    _len = length;
    if (length == 0)
      _len = (u32bit)strlen(S);

    _rcBase = 0;
    _rcLen  = _len;

    _doRevComp = revcomp;

    if (_doRevComp)
      _pos = _len-1;
  };

public:
  FastAAccessor(FastASequenceInCore *S, bool revcomp=false) {
    FastAAccessorInit(S->sequence(), S->sequenceLength(), revcomp);
  };

  FastAAccessor(char *S, u32bit length=0, bool revcomp=false) {
    FastAAccessorInit(S, length, revcomp);
  };

private:

  //  Given a range in the forward string, we can reverse-complement
  //  just that range.  This amounts to translating the forward string
  //  to make the beginning of the range be the origin, then
  //  reverseing the range, then translating the sequence back to the
  //  original origin.
  //
  u32bit  rc(u32bit p) const {
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
  bool  setRange(u32bit base=0, u32bit length=0) {
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
      fprintf(stderr, "FastAAccessor::setReverseComplementRange()-- base="u32bitFMT" and length="u32bitFMT" exceed sequence length of "u32bitFMT"\n",
              base, length, _len);
      abort();
      return(false);
    }
  };


  //  True if this physical location is valid.
  //
  bool isValid(void) {
    return(_pos < _len);
  };


  char operator[](u32bit p) const {
    if ((_rcBase <= p) && (p < _rcBase + _rcLen)) {
      if (_doRevComp) return(complementSymbol[(int)_seq[rc(p)]]);
      else            return(_seq[p]);
    } else {
      fprintf(stderr, "operator[]-- Tried to access to "u32bitFMT", but range is "u32bitFMT"-"u32bitFMT"\n",
              p, _rcBase, _rcBase + _rcLen);
      abort();
      return(0);
    }
  };



  //  Set the accessor to some position.
  //
  bool setPosition(u32bit p) {
    if ((_rcBase <= p) && (p < _rcBase + _rcLen)) {
      if (_doRevComp) _pos = rc(p);
      else            _pos = p;
      return(true);
    } else {
      fprintf(stderr, "setPosition()-- Tried to set to "u32bitFMT", but range is "u32bitFMT"-"u32bitFMT".\n",
              p, _rcBase, _rcBase + _rcLen);
      abort();
      return(false);
    }
  };

  u32bit getPosition(void) {
    if (_doRevComp) return(rc(_pos));
    else            return(_pos);
  };

  u32bit getRangeBegin(void)  { return(_rcBase); };
  u32bit getRangeEnd(void)    { return(_rcBase + _rcLen); };
  u32bit getRangeLength(void) { return(_rcLen); };

  bool   extendLeft(s32bit x) {
    _rcLen  += x;
    if (_doRevComp == false)
      _rcBase -= x;

    if ((_rcBase < _len) && (_rcBase + _rcLen <= _len)) {
      fprintf(stderr, "FastAAccessor::extendLeft()-- extend by "s32bitFMT" makes invalid: length is "u32bitFMT", new range is "u32bitFMT"-"u32bitFMT"\n",
              x, _len, _rcBase, _rcBase + _rcLen);
      abort();
      return(false);
    }
    return(true);
  };

  bool   extendRight(s32bit x)   {
    _rcLen  += x;
    if (_doRevComp == true)
      _rcBase -= x;

    if ((_rcBase < _len) && (_rcBase + _rcLen <= _len)) {
      fprintf(stderr, "FastAAccessor::extendRight()-- extend by "s32bitFMT" makes invalid: length is "u32bitFMT", new range is "u32bitFMT"-"u32bitFMT"\n",
              x, _len, _rcBase, _rcBase + _rcLen);
      abort();
      return(false);
    }
    return(true);
  };

  char operator*(void) const {
    SANITY("FastAAccessor::operator*()", _pos);
    if (_doRevComp) return(complementSymbol[(int)_seq[_pos]]);
    else            return(_seq[_pos]);
  };

  char get(void) const {
    SANITY("FastAAccessor::get()", _pos);
    if (_doRevComp) return(complementSymbol[(int)_seq[_pos]]);
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

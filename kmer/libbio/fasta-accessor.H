#ifndef FASTA_ACCESSOR_H
#define FASTA_ACCESSOR_H

#include "bio++.H"


//  Define this to do bounds checking
//
#if 1
#define SANITY(NAM, POS) \
  if (((_doRevComp) ? rc((POS)) : (POS)) > _len) \
    fprintf(stderr, "%s-- position "u32bitFMT" larger than length "u32bitFMT"\n", NAM, ((_doRevComp) ? rc((POS)) : (POS)), _len)
#else
#define SANITY(NAM, POS)
#endif


class FastAAccessor {
public:
  FastAAccessor(FastASequenceInCore *S, bool revcomp=false) {
    _seq = S->sequence();
    _len = S->sequenceLength();
    _pos = 0;

    _doRevComp = revcomp;

    _rcBase = 0;
    _rcLen  = _len;

    if (_doRevComp)
      _pos = _len-1;
  };

  FastAAccessor(char *S, u32bit length=0, bool revcomp=false) {
    _seq = S;
    _len = length;
    if (length == 0)
      _len = (u32bit)strlen(S);
    _pos = 0;

    _doRevComp = revcomp;

    _rcBase = 0;
    _rcLen  = _len;

    if (_doRevComp)
      _pos = _len-1;
  };

  ~FastAAccessor() {
  };

private:
  u32bit  rc(u32bit p) const {
    return(_rcBase + _rcLen - (p - _rcBase) - 1);
  };

public:

  //  For iterating over reverse complement regions of a sequence.
  //  The base is the start of the range, and the length is, the
  //  length.  In the forward coordinate system.
  //
  //  e.g.,  (500, 250) would be:
  //
  //  |-----|-----|XXXXX|-----|
  //  0    250   500   750 1000
  //
  //  Set both to zero (also the default) to unset the range
  //
  void  setRange(u32bit base=0, u32bit length=0) {
    SANITY("FastAAccessor::setReverseComplementRange()", base);
    SANITY("FastAAccessor::setReverseComplementRange()", base+length);
    u32bit p = rc(_pos);
    if ((base == 0) && (length == 0)) {
      _rcBase = 0;
      _rcLen  = _len;
    } else {
      _rcBase = base;
      _rcLen  = length;
    }
    if (_doRevComp)
      _pos = rc(p);
  };

  u32bit getRangeBegin(void)  { return(_rcBase); };
  u32bit getRangeEnd(void)    { return(_rcBase + _rcLen); };

  void   extendLeft(s32bit x) {
    fprintf(stderr, "extendLeft()-- by %d\n", (int)x);
    if (_doRevComp) {
      u32bit p = rc(_pos);
      _rcLen  += x;
      _pos = rc(p);
    } else {
      _rcBase -= x;
      _rcLen  += x;
    }
  };
  void   extendRight(s32bit x)   {
    fprintf(stderr, "extendRight()-- by %d\n", (int)x);
    if (_doRevComp) {
      u32bit p = rc(_pos);
      _rcBase -= x;
      _rcLen  += x;
      _pos = rc(p);
    } else {
      _rcLen  += x;
    }
  };

  char operator[](u32bit p) const {
    SANITY("FastAAccessor::operator[]", p);
  
    if (_doRevComp) return(complementSymbol[(int)_seq[rc(p)]]);
    else            return(_seq[p]);
  };


  char operator*(void) const {
    SANITY("FastAAccessor::operator()", _pos);
    if (_doRevComp) return(complementSymbol[(int)_seq[_pos]]);
    else            return(_seq[_pos]);
  };

  char get(void) const {
    SANITY("FastAAccessor::get()", _pos);
    if (_doRevComp) return(complementSymbol[(int)_seq[_pos]]);
    else            return(_seq[_pos]);
  };

  //  Prefix operators -- similar to setPosition, we let these be
  //  invalid without sanity checking, because often we want to do
  //  "for (i=0; i<len; i++) { .... ++ACCESSOR; } which would
  //  increment one past the actual end, but we never access the base
  //  there.
  //
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

  //bool isValid(s32bit offset=0, bool report=false) {
  bool isValid(void) {
#if 0
    if (_doRevComp) return( (-offset <= rc(_pos)) && ((rc(_pos) + offset) < _len) );
    else            return( (-offset <=    _pos)  && ((   _pos  + offset) < _len) );
#else
    //  True if this physical location is valid.
    return(_pos < _len);
#endif
  };

  bool setPosition(u32bit p) {
    //  We don't check sanity here, because this is also useful as a
    //  test, e.g., set the position somewhere, if valid it returns
    //  true.
    if (_doRevComp) _pos = rc(p);
    else            _pos = p;
    return(_pos < _len);
  };

  u32bit getPosition(void) {
    if (_doRevComp) return(rc(_pos));  //  oddly symmetric
    else            return(_pos);
  };

  u32bit getLength(void) {
    return(_len);
  };

private:
public:
  char    *_seq;
  u32bit   _pos;
  u32bit   _len;

  u32bit   _rcBase;
  u32bit   _rcLen;

  bool     _doRevComp;
};

#endif  //  FASTA_ACCESSOR_H

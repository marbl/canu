#ifndef FASTA_ACCESSOR_H
#define FASTA_ACCESSOR_H

#include "bio++.H"


//  Define this to do bounds checking
//

#if 0
//  The old sanity check
  if (((_doRevComp) ? rc((POS)) : (POS)) > _len) \
    fprintf(stderr, "%s-- position "u32bitFMT" larger than length "u32bitFMT"\n", NAM, ((_doRevComp) ? rc((POS)) : (POS)), _len)
#endif


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


class FastAAccessor {
private:
public:
  char    *_seq;
  u32bit   _pos;
  u32bit   _len;

  u32bit   _rcBase;
  u32bit   _rcLen;

  bool     _doRevComp;

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
  //  The physical location (_pos) doesn't change, but this
  //  will change the mapping back to a forward coordinate.
  //
  void  setRange(u32bit base=0, u32bit length=0) {
    SANITY("FastAAccessor::setReverseComplementRange()", base);
    SANITY("FastAAccessor::setReverseComplementRange()", base+length);

    if ((base == 0) && (length == 0)) {
      _rcBase = 0;
      _rcLen  = _len;
    } else {
      _rcBase = base;
      _rcLen  = length;
    }
  };


#if 0
  char operator[](u32bit p) const {
    SANITY("FastAAccessor::operator[]", p);
  
    if (_doRevComp) return(complementSymbol[(int)_seq[rc(p)]]);
    else            return(_seq[p]);
  };
#endif


  //  True if this physical location is valid.
  bool isValid(void) {
    return(_pos < _len);
  };

  bool setPosition(u32bit p) {
    if ((_rcBase <= p) && (p < _rcBase + _rcLen)) {
      if (_doRevComp) _pos = rc(p);
      else            _pos = p;
      return(true);
    } else {
      fprintf(stderr, "setPosition()-- Tried to set to %d, but range is %d-%d.\n",
              p, _rcBase, _rcBase + _rcLen);
      abort();
      return(false);
    }
  };

  u32bit getPosition(void) {
    if (_doRevComp) return(rc(_pos));  //  oddly symmetric
    else            return(_pos);
  };




  ////////////////////////////////////////
  //
  //  Stuff below here is probably correct.
  //

  u32bit getRangeBegin(void)  { return(_rcBase); };
  u32bit getRangeEnd(void)    { return(_rcBase + _rcLen); };
  u32bit getRangeLength(void) { return(_rcLen); };

  void   extendLeft(s32bit x) {
    //fprintf(stderr, "extendLeft()-- by %d\n", (int)x);
    if (_doRevComp) {
      _rcLen  += x;
    } else {
      _rcBase -= x;
      _rcLen  += x;
    }
  };

  void   extendRight(s32bit x)   {
    //fprintf(stderr, "extendRight()-- by %d\n", (int)x);
    if (_doRevComp) {
      _rcBase -= x;
      _rcLen  += x;
    } else {
      _rcLen  += x;
    }
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

  u32bit getLength(void) {
    return(_len);
  };
};

#endif  //  FASTA_ACCESSOR_H

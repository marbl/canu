#ifndef FASTA_ACCESSOR_H
#define FASTA_ACCESSOR_H

#include "bio++.H"

//  Define this to do bounds checking
#define SANITY

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

  FastAAccessor(FastASequenceOnDisk *S, bool revcomp=false) {
    fprintf(stderr, "SORRY!  FastAAccessor() isn't implemented for FastASequenceOnDisk.\n");
    exit(1);
  };

  FastAAccessor(char *S, u32bit length=0, bool revcomp=false) {
    _seq = S;
    _len = strlen(S);
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
    return(_rcBase + _rcBase + _rcLen - p - 1);
  };

public:

  //  For interating over reverse complement regions of a sequence.
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
  void  setReverseComplementRange(u32bit base=0, u32bit length=0) {
#ifdef SANITY
    if (base > _len)
      fprintf(stderr, "FastAAccessor::setReverseComplementRange()-- position "u32bitFMT" larger than length "u32bitFMT"\n", base, _len);
#endif
    if ((base == 0) && (length == 0)) {
      _rcBase = 0;
      _rcLen  = _len;
    } else {
      _rcBase = base;
      _rcLen  = length;
    }
  };

  char operator[](u32bit p) const {
#ifdef SANITY
    if (p > _len)
      fprintf(stderr, "FastAAccessor::operator[]-- position "u32bitFMT" larger than length "u32bitFMT"\n", p, _len);
#endif
    if (_doRevComp) return(complementSymbol[(int)_seq[rc(p)]]);
    else            return(_seq[p]);
  };


  char operator*(void) const {
#ifdef SANITY
    if (_pos > _len)
      fprintf(stderr, "FastAAccessor::operator*()-- position "u32bitFMT" larger than length "u32bitFMT"\n", _pos, _len);
#endif
    if (_doRevComp) return(complementSymbol[(int)_seq[_pos]]);
    else            return(_seq[_pos]);
  };

  char get(void) const {
#ifdef SANITY
    if (_pos > _len)
      fprintf(stderr, "FastAAccessor::get()-- position "u32bitFMT" larger than length "u32bitFMT" rc=%d\n", _pos, _len, (int)_doRevComp);
#endif
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

  bool isValid(void) {
    return(_pos < _len);
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


#if 0
  //  Postfix operators -- a little expensive, we would need to
  //  copy construct the return value.
  FastAAccessor const operator--(int x) {
  };
  FastAAccessor const operator++(int x) {
  };
#endif

private:
  char    *_seq;
  u32bit   _pos;
  u32bit   _len;

  u32bit   _rcBase;
  u32bit   _rcLen;

  bool     _doRevComp;
};

#endif  //  FASTA_ACCESSOR_H

#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include "seatac.H"
#include "bio++.H"

encodedQuery::encodedQuery(char const    *seq,
                           u32bit         seqLen,
                           u32bit         k,
                           bool           rc) {
  _seq            = seq;
  _seqLen         = seqLen;
  _merSize        = k;
  _rc             = rc;

  _seqPos         = 0;

  _substring      = u64bitZERO;
  _mermask        = u64bitMASK(2 * _merSize);
  _timeUntilValid = _merSize;
}


bool
encodedQuery::getMer(u64bit &mer, u32bit &pos) {
  bool  found = false;

  mer = u64bitZERO;
  pos = u32bitZERO;

  if (_rc) {

    while (!found && (_seqPos < _seqLen)) {
      _substring <<= 2;
      _substring  &= _mermask;

      if (letterToBits[_seq[_seqLen - 1 - _seqPos]] != 0xff) {
        _substring |= letterToBits[ complementSymbol[ _seq[_seqLen - 1 - _seqPos] ]];
        _timeUntilValid--;
      } else {
        _timeUntilValid = _merSize;
      }

      _seqPos++;

      if (_seqPos >= _merSize) {
        mer   = _substring;
        pos   = _seqPos - _merSize;
        found = _timeUntilValid <= 0;
      }
    }

  } else {

    while (!found && (_seqPos < _seqLen)) {
      _substring <<= 2;
      _substring  &= _mermask;

      if (letterToBits[_seq[_seqPos]] != 0xff) {
        _substring |= letterToBits[_seq[_seqPos]];
        _timeUntilValid--;
      } else {
        _timeUntilValid = _merSize;
      }

      _seqPos++;

      if (_seqPos >= _merSize) {
        mer   = _substring;
        pos   = _seqPos - _merSize;
        found = _timeUntilValid <= 0;
      }
    }

  }

  return(found);
}




encodedQuery::~encodedQuery() {
}

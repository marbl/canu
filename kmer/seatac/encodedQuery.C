#include <stdio.h>
#include <stdlib.h>
#include "seatac.H"
#include "bio++.H"

encodedQuery::encodedQuery(char const    *seq,
                           uint32         seqLen,
                           uint32         k,
                           bool           rc) {
  _seq            = seq;
  _seqLen         = seqLen;
  _merSize        = k;
  _rc             = rc;

  _seqPos         = 0;

  _substring      = uint64ZERO;
  _mermask        = uint64MASK(2 * _merSize);
  _timeUntilValid = _merSize;
}


bool
encodedQuery::getMer(uint64 &mer, uint32 &pos) {
  bool  found = false;

  mer = uint64ZERO;
  pos = uint32ZERO;

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

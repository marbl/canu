#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include "searchGENOME.H"
#include "libbri.H"

encodedQuery::encodedQuery(char const           *seq,
                           u32bit                seqLen,
                           u32bit                k,
                           bool                  rc) {

  _mersTotal = 0;
  _mersAvail = 0;
  _mers      = 0L;
  _skip      = 0L;

  if (k > seqLen)
    return;

  _mersTotal = seqLen - k + 1;
  _mersAvail = 0;
  _mers      = new u64bit [_mersTotal];
  _skip      = new u8bit  [_mersTotal];

  u64bit   substring      = u64bitZERO;
  u64bit   mermask        = u64bitMASK(2 * k);
  s32bit   timeUntilValid = k;

  if (rc) {
    for (u32bit i=0; i<seqLen; ) {
      substring <<= 2;
      substring  &= mermask;

      if (validSymbol[seq[seqLen - 1 - i]]) {
        substring |= compressSymbol[ complementSymbol[ seq[seqLen - 1 - i] ]];
        timeUntilValid--;
      } else {
        timeUntilValid = k;
      }

      i++;

      if (i >= k) {
        _mers[i-k]  = substring;
        _skip[i-k]  = timeUntilValid > 0;
        _mersAvail += 1 - _skip[i-k];
      }
    }
  } else {
    for (u32bit i=0; i<seqLen; ) {
      substring <<= 2;
      substring  &= mermask;

      if (validSymbol[seq[i]]) {
        substring |= compressSymbol[ seq[i] ];
        timeUntilValid--;
      } else {
        timeUntilValid = k;
      }

      i++;

      if (i >= k) {
        _mers[i-k]  = substring;
        _skip[i-k]  = timeUntilValid > 0;
        _mersAvail += 1 - _skip[i-k];
      }
    }
  }
}

encodedQuery::~encodedQuery() {
  delete [] _mers;
  delete [] _skip;
}

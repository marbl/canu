#include "posix.H"
#include <stdio.h>
#include <stdlib.h>

#include "searchGENOME.H"
#include "libbri.H"

encodedQuery::encodedQuery(unsigned char const  *seq,
                           u32bit                seqLen,
                           u32bit                k,
                           bool                  rc) {
  _mers = 0L;
  _skip = 0L;

  if (k > seqLen) {
    _mersLen = 0;
    return;
  }

  _mersLen   = seqLen - k + 1;
  _mers      = new u64bit [_mersLen];
  _skip      = new u32bit [_mersLen];

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
        _mers[i-k] = substring;
        _skip[i-k] = timeUntilValid > 0;
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
        _mers[i-k] = substring;
        _skip[i-k] = timeUntilValid > 0;
      }
    }
  }
}

encodedQuery::~encodedQuery() {
  delete [] _mers;
  delete [] _skip;
}

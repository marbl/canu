#ifndef UNARY_ENCODING_H
#define UNARY_ENCODING_H

#include "bitPacking.h"


//  Routines to store and retrieve a unary encoded number to/from a
//  bit packed word array based at 'ptr' and currently at location
//  'pos'.  Both routines return the size of the encoded number in
//  'siz'.



//  The usual unary encoding.  Store the number n as n 0 bits followed
//  by a single 1 bit.
//
//  0 -> 1
//  1 -> 01
//  2 -> 001
//  3 -> 0001
//  4 -> 00001
//
//  See the decoder as to why we use 0 instead of 1 for the count.


inline
void
setUnaryEncodedNumber(uint64 *ptr,
                      uint64  pos,
                      uint64 *siz,
                      uint64  val) {

  *siz = val + 1;

  while (val >= 64) {
    setDecodedValue(ptr, pos, 64, uint64ZERO);
    pos += 64;
    val -= 64;
    siz += 64;
  }

  setDecodedValue(ptr, pos, val + 1, uint64ONE);
  pos += val + 1;
}



inline
uint64
getUnaryEncodedNumber(uint64 *ptr,
                      uint64  pos,
                      uint64 *siz) {
  uint64 val = uint64ZERO;
  uint64 enc = uint64ZERO;

  //  How many whole words are zero?
  //
  enc = getDecodedValue(ptr, pos, 64);
  while (enc == uint64ZERO) {
    val += 64;
    pos += 64;
    enc  = getDecodedValue(ptr, pos, 64);
  }

  //  This word isn't zero.  Count how many bits are zero (see, the
  //  choice of 0 or 1 for the encoding wasn't arbitrary!)
  //
  val += 64 - logBaseTwo64(enc);

  *siz = val + 1;

  return(val);
}


#endif  //  UNARY_ENCODING_H

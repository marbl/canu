#ifndef ELIAS_DELTA_ENCODING_H
#define ELIAS_DELTA_ENCODING_H

#include "bitPacking.h"

inline
void
setEliasDeltaEncodedNumber(uint64 *ptr,
                           uint64  pos,
                           uint64 *siz,
                           uint64  val) {
  uint64 b = logBaseTwo64(val);
  setEliasGammaEncodedNumber(ptr, pos, siz, b);
  pos += *siz;
  setDecodedValue(ptr, pos, b-1, val);
  *siz += b-1;
}


inline
uint64
getEliasDeltaEncodedNumber(uint64 *ptr,
                           uint64  pos,
                           uint64 *siz) {
  uint64 b = getEliasGammaEncodedNumber(ptr, pos, siz) - 1;
  pos  += *siz;
  *siz += b;
  return(uint64ONE << b | getDecodedValue(ptr, pos, b));
}



#endif  //  ELIAS_DELTA_ENCODING_H

#ifndef ELIAS_DELTA_ENCODING_H
#define ELIAS_DELTA_ENCODING_H

#include "bitPacking.h"

inline
void
setEliasDeltaEncodedNumber(u64bit *ptr,
                           u64bit  pos,
                           u64bit *siz,
                           u64bit  val) {
  u64bit b = logBaseTwo64(val);
  setEliasGammaEncodedNumber(ptr, pos, siz, b);
  pos += *siz;
  setDecodedValue(ptr, pos, b-1, val);
  *siz += b-1;
}


inline
u64bit
getEliasDeltaEncodedNumber(u64bit *ptr,
                           u64bit  pos,
                           u64bit *siz) {
  u64bit b = getEliasGammaEncodedNumber(ptr, pos, siz) - 1;
  pos  += *siz;
  *siz += b;
  return(u64bitONE << b | getDecodedValue(ptr, pos, b));
}



#endif  //  ELIAS_DELTA_ENCODING_H

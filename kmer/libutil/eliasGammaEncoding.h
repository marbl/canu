#ifndef ELIAS_GAMMA_ENCODING_H
#define ELIAS_GAMMA_ENCODING_H

#include "bitPacking.h"

inline
void
setEliasGammaEncodedNumber(uint64 *ptr,
                           uint64  pos,
                           uint64 *siz,
                           uint64  val) {
  uint64 b = logBaseTwo64(val);
  setUnaryEncodedNumber(ptr, pos, siz, b);
  pos += *siz;
  setDecodedValue(ptr, pos, b, val);
  *siz += b;
}


inline
uint64
getEliasGammaEncodedNumber(uint64 *ptr,
                           uint64  pos,
                           uint64 *siz) {
  uint64 b = getUnaryEncodedNumber(ptr, pos, siz);
  pos  += *siz;
  *siz += b;
  return(getDecodedValue(ptr, pos, b));
}



#endif  //  ELIAS_GAMMA_ENCODING_H

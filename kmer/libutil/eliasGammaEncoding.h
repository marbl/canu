#ifndef ELIAS_GAMMA_ENCODING_H
#define ELIAS_GAMMA_ENCODING_H

#include "bitPacking.h"

inline
void
setEliasGammaEncodedNumber(u64bit *ptr,
                           u64bit  pos,
                           u64bit *siz,
                           u64bit  val) {
  u64bit b = logBaseTwo64(val);
  setUnaryEncodedNumber(ptr, pos, siz, b);
  pos += *siz;
  setDecodedValue(ptr, pos, b, val);
  *siz += b;
}


inline
u64bit
getEliasGammaEncodedNumber(u64bit *ptr,
                           u64bit  pos,
                           u64bit *siz) {
  u64bit b = getUnaryEncodedNumber(ptr, pos, siz);
  pos  += *siz;
  *siz += b;
  return(getDecodedValue(ptr, pos, b));
}



#endif  //  ELIAS_GAMMA_ENCODING_H

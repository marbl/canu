#ifndef BIT_PACKING_H
#define BIT_PACKING_H

#include "libbritypes.h"

inline
u64bit
getDecodedValue(u64bit *ptr,
                u64bit  pos,
                u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1
  u64bit ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);
  } else {
    ret  = (ptr[wrd] & u64bitMASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & u64bitMASK(b2);
  }

  ret &= u64bitMASK(siz);

  return(ret);
}

inline
void
setDecodedValue(u64bit *ptr,
                u64bit  pos,
                u64bit  siz,
                u64bit  val) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1

  val &= u64bitMASK(siz);

  if (b1 >= siz) {
    ptr[wrd] &= ~( u64bitMASK(siz) << (b1-siz) );
    ptr[wrd] |= val << (b1-siz);
  } else {
    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (val & (u64bitMASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(u64bitMASK(b2) << (64-b2));
    ptr[wrd+1] |= (val & (u64bitMASK(b2))) << (64-b2);
  }
}

//  Special case for storing 32-bit values
//
#if 0
#error Hey!  getDecodedValue32() and setDecodedValue32() are untested!
inline
u64bit
getDecodedValue32(u64bit *ptr,
                  u64bit  pos) {
  u32bit s = (pos % u64bitONE) << 5;

  return((ptr[pos >> 1] >> s) & 0x00000000ffffffffllu);
}

inline
void
setDecodedValue32(u64bit *ptr,
                  u64bit  pos,
                  u64bit  val) {
  u32bit s = (pos % u64bitONE) << 5;

  ptr[pos >> 1] &= 0x00000000ffffffffllu << s;
  ptr[pos >> 1] |= val << s;
}
#endif

inline
u64bit
preIncrementDecodedValue(u64bit *ptr,
                         u64bit  pos,
                         u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1
  u64bit ret  = 0;

  if (b1 >= siz) {
    ret  = ptr[wrd] >> (b1 - siz);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1-siz) );
    ptr[wrd] |= ret << (b1-siz);
  } else {
    ret  = (ptr[wrd] & u64bitMASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & u64bitMASK(b2);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(u64bitMASK(b2) << (64-b2));
    ptr[wrd+1] |= (ret & (u64bitMASK(b2))) << (64-b2);
  }

  return(ret);
}

inline
u64bit
preDecrementDecodedValue(u64bit *ptr,
                         u64bit  pos,
                         u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1
  u64bit ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1-siz) );
    ptr[wrd] |= ret << (b1-siz);
  } else {
    ret  = (ptr[wrd] & u64bitMASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & u64bitMASK(b2);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(u64bitMASK(b2) << (64-b2));
    ptr[wrd+1] |= (ret & (u64bitMASK(b2))) << (64-b2);
  }

  return(ret);
}




inline
u64bit
postIncrementDecodedValue(u64bit *ptr,
                          u64bit  pos,
                          u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1
  u64bit ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1-siz) );
    ptr[wrd] |= ret << (b1-siz);
  } else {
    ret  = (ptr[wrd] & u64bitMASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & u64bitMASK(b2);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(u64bitMASK(b2) << (64-b2));
    ptr[wrd+1] |= (ret & (u64bitMASK(b2))) << (64-b2);
  }

  ret--;
  ret &= u64bitMASK(siz);

  return(ret);
}

inline
u64bit
postDecrementDecodedValue(u64bit *ptr,
                          u64bit  pos,
                          u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1
  u64bit ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1-siz) );
    ptr[wrd] |= ret << (b1-siz);
  } else {
    ret  = (ptr[wrd] & u64bitMASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & u64bitMASK(b2);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(u64bitMASK(b2) << (64-b2));
    ptr[wrd+1] |= (ret & (u64bitMASK(b2))) << (64-b2);
  }

  ret++;
  ret &= u64bitMASK(siz);

  return(ret);
}

inline
u64bit
sumDecodedValue(u64bit *ptr,
                u64bit  pos,
                u64bit  siz,
                u64bit  val) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit b2  = siz - b1;  //  Only used if siz > b1
  u64bit ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret += val;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1-siz) );
    ptr[wrd] |= ret << (b1-siz);
  } else {
    ret  = (ptr[wrd] & u64bitMASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & u64bitMASK(b2);

    ret += val;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(u64bitMASK(b2) << (64-b2));
    ptr[wrd+1] |= (ret & (u64bitMASK(b2))) << (64-b2);
  }

  return(ret);
}

#endif  //  BIT_PACKING_H

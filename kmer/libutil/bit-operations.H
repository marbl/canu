#ifndef BIT_OPERATIONS_H
#define BIT_OPERATIONS_H

#include "libbritypes.h"

//  countNumberOfSetBits:  returns the number of set bits.
//  The super-slow way to do this is:
//    u32bit r = u32bitZERO;
//    while (x) {
//      r++;
//      x >>= 1;
//    }
//    return(r);
//  But we have a super-fast way.
//  Credit:  http://remus.rutgers.edu/~rhoads/Code/bitcount3.c

//  logBaseTwo:  Returns the 'index' of the highest bit set.

//  reverseBits:  Reverse the bits in a word.
//  Credit:  The BSD fortunes file.


inline
u32bit
countNumberOfSetBits(u32bit x) {
  x -= (0xaaaaaaaa & x) >> 1;
  x  = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  x += x >> 4;
  x &= 0x0f0f0f0f;
  x += x >> 8;
  x += x >> 16;
  x &= 0x000000ff;

  return(x);
}

inline
u32bit
logBaseTwo(u32bit x) {

  //  Step one, set all the bits below the highest bit set.
  //  Credit:  Ross Lippert
  //  
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;

  //  Step two, count the number of bits set.
  //
  return(countNumberOfSetBits(x));
}


inline
u32bit
logBaseTwo_64(u64bit x) {
  u32bit r = 64;
  u64bit m = u64bitONE << 63;

  while ((x & m) == 0) {
    r--;
    m = u64bitONE << (r-1);
  }

  return(r);
}


inline
u32bit
reverseBits(u32bit x) {
  x = ((x >>  1) & 0x55555555) | ((x <<  1) & 0xaaaaaaaa);
  x = ((x >>  2) & 0x33333333) | ((x <<  2) & 0xcccccccc);
  x = ((x >>  4) & 0x0f0f0f0f) | ((x <<  4) & 0xf0f0f0f0);
  x = ((x >>  8) & 0x00ff00ff) | ((x <<  8) & 0xff00ff00);
  x = ((x >> 16) & 0x0000ffff) | ((x << 16) & 0xffff0000);

  return(x);
}


inline
u64bit
reverseComplementMer(u32bit ms, u64bit fmer) {

  //  Reverse the mer
  //
  u64bit rmer = fmer;
  rmer = ((rmer >>  2) & 0x3333333333333333llu) | ((rmer <<  2) & 0xccccccccccccccccllu);
  rmer = ((rmer >>  4) & 0x0f0f0f0f0f0f0f0fllu) | ((rmer <<  4) & 0xf0f0f0f0f0f0f0f0llu);
  rmer = ((rmer >>  8) & 0x00ff00ff00ff00ffllu) | ((rmer <<  8) & 0xff00ff00ff00ff00llu);
  rmer = ((rmer >> 16) & 0x0000ffff0000ffffllu) | ((rmer << 16) & 0xffff0000ffff0000llu);
  rmer = ((rmer >> 32) & 0x00000000ffffffffllu) | ((rmer << 32) & 0xffffffff00000000llu);

  //  Complement the bases
  //
  rmer ^= 0xffffffffffffffffllu;

  //  Shift and mask
  //
  rmer >>= 64 - ms * 2;
  rmer  &= u64bitMASK(ms * 2);
  return(rmer);
}



#endif  //  BIT_OPERATIONS_H

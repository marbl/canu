#ifndef BRI_BITS_H
#define BRI_BITS_H

//  For dealing with the bits in bytes.

//  I wish I could claim these.
//
//  Freed, Edwin E. 1983. "Binary Magic Number" Dr. Dobbs Journal
//  Vol. 78 (April) pp. 24-37
//
//  Supposedly tells us how to reverse the bits in a word, count the number
//  of set bits in a words and more.
//
//  A bit of verbage on counting the number of set bits.  The naive way
//  is to loop and shift:
//
//      u32bit r = u32bitZERO;
//      while (x) {
//        r++;
//        x >>= 1;
//      }
//      return(r);
//
//  http://remus.rutgers.edu/~rhoads/Code/bitcount3.c has an optimized
//  method:
//
//      x -= (0xaaaaaaaa & x) >> 1;
//      x  = (x & 0x33333333) + ((x >> 2) & 0x33333333);
//      x += x >> 4;
//      x &= 0x0f0f0f0f;
//      x += x >> 8;
//      x += x >> 16;
//      x &= 0x000000ff;
//      return(x);
//
//  No loops!
//
//  Freed's methods are easier to understand, and just as fast.
//
//  Using our bit counting routines, Ross Lippert suggested a nice
//  way of computing log2 -- use log2 shifts to fill up the lower
//  bits, then count bits.  See logBaseTwo*()
//


inline
u32bit
reverseBits32(u32bit x) {
  x = ((x >>  1) & 0x55555555LU) | ((x <<  1) & 0xaaaaaaaaLU);
  x = ((x >>  2) & 0x33333333LU) | ((x <<  2) & 0xccccccccLU);
  x = ((x >>  4) & 0x0f0f0f0fLU) | ((x <<  4) & 0xf0f0f0f0LU);
  x = ((x >>  8) & 0x00ff00ffLU) | ((x <<  8) & 0xff00ff00LU);
  x = ((x >> 16) & 0x0000ffffLU) | ((x << 16) & 0xffff0000LU);
  return(x);
}

inline
u64bit
reverseBits64(u64bit x) {
  x = ((x >>  1) & 0x5555555555555555LLU) | ((x <<  1) & 0xaaaaaaaaaaaaaaaaLLU);
  x = ((x >>  2) & 0x3333333333333333LLU) | ((x <<  2) & 0xccccccccccccccccLLU);
  x = ((x >>  4) & 0x0f0f0f0f0f0f0f0fLLU) | ((x <<  4) & 0xf0f0f0f0f0f0f0f0LLU);
  x = ((x >>  8) & 0x00ff00ff00ff00ffLLU) | ((x <<  8) & 0xff00ff00ff00ff00LLU);
  x = ((x >> 16) & 0x0000ffff0000ffffLLU) | ((x << 16) & 0xffff0000ffff0000LLU);
  x = ((x >> 32) & 0x00000000ffffffffLLU) | ((x << 32) & 0xffffffff00000000LLU);
  return(x);
}



inline
u32bit
countNumberOfSetBits32(u32bit x) {
  x = ((x >>  1) & 0x55555555LU) + (x & 0x55555555LU);
  x = ((x >>  2) & 0x33333333LU) + (x & 0x33333333LU);          
  x = ((x >>  4) & 0x0f0f0f0fLU) + (x & 0x0f0f0f0fLU);
  x = ((x >>  8) & 0x00ff00ffLU) + (x & 0x00ff00ffLU);              
  x = ((x >> 16) & 0x0000ffffLU) + (x & 0x0000ffffLU);
  return(x);
}

inline
u64bit
countNumberOfSetBits64(u64bit x) {
  x = ((x >>  1) & 0x5555555555555555LLU) + (x & 0x5555555555555555LLU);
  x = ((x >>  2) & 0x3333333333333333LLU) + (x & 0x3333333333333333LLU);          
  x = ((x >>  4) & 0x0f0f0f0f0f0f0f0fLLU) + (x & 0x0f0f0f0f0f0f0f0fLLU);
  x = ((x >>  8) & 0x00ff00ff00ff00ffLLU) + (x & 0x00ff00ff00ff00ffLLU);              
  x = ((x >> 16) & 0x0000ffff0000ffffLLU) + (x & 0x0000ffff0000ffffLLU);
  x = ((x >> 32) & 0x00000000ffffffffLLU) + (x & 0x00000000ffffffffLLU);
  return(x);
}



inline
u32bit
logBaseTwo32(u32bit x) {
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return(countNumberOfSetBits32(x));
}

inline
u64bit
logBaseTwo64(u64bit x) {
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  return(countNumberOfSetBits64(x));
}




#endif  //  BRI_BITS_H

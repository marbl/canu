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
  x = ((x >>  1) & u32bitNUMBER(0x55555555)) | ((x <<  1) & u32bitNUMBER(0xaaaaaaaa));
  x = ((x >>  2) & u32bitNUMBER(0x33333333)) | ((x <<  2) & u32bitNUMBER(0xcccccccc));
  x = ((x >>  4) & u32bitNUMBER(0x0f0f0f0f)) | ((x <<  4) & u32bitNUMBER(0xf0f0f0f0));
  x = ((x >>  8) & u32bitNUMBER(0x00ff00ff)) | ((x <<  8) & u32bitNUMBER(0xff00ff00));
  x = ((x >> 16) & u32bitNUMBER(0x0000ffff)) | ((x << 16) & u32bitNUMBER(0xffff0000));
  return(x);
}

inline
u64bit
reverseBits64(u64bit x) {
  x = ((x >>  1) & u64bitNUMBER(0x5555555555555555)) | ((x <<  1) & u64bitNUMBER(0xaaaaaaaaaaaaaaaa));
  x = ((x >>  2) & u64bitNUMBER(0x3333333333333333)) | ((x <<  2) & u64bitNUMBER(0xcccccccccccccccc));
  x = ((x >>  4) & u64bitNUMBER(0x0f0f0f0f0f0f0f0f)) | ((x <<  4) & u64bitNUMBER(0xf0f0f0f0f0f0f0f0));
  x = ((x >>  8) & u64bitNUMBER(0x00ff00ff00ff00ff)) | ((x <<  8) & u64bitNUMBER(0xff00ff00ff00ff00));
  x = ((x >> 16) & u64bitNUMBER(0x0000ffff0000ffff)) | ((x << 16) & u64bitNUMBER(0xffff0000ffff0000));
  x = ((x >> 32) & u64bitNUMBER(0x00000000ffffffff)) | ((x << 32) & u64bitNUMBER(0xffffffff00000000));
  return(x);
}


#if (__GNUC__ > 3) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#define PREFETCH(x) __builtin_prefetch((x), 0, 0)
#else
#define PREFETCH(x)
#endif




//  Amazingingly, this is slower.  From what I can google, the builtin
//  is using the 2^16 lookup table method - so a 64-bit popcount does
//  4 lookups in the table and sums.  Bad cache performance in codes
//  that already have bad cache performance, I'd guess.
//
//#if (__GNUC__ > 3) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
//#define BUILTIN_POPCOUNT
//#endif

#ifdef BUILTIN_POPCOUNT

inline
u32bit
countNumberOfSetBits32(u32bit x) {
  return(__builtin_popcount(x));
}

inline
u64bit
countNumberOfSetBits64(u64bit x) {
  return(__builtin_popcountll(x));
}

#else

inline
u32bit
countNumberOfSetBits32(u32bit x) {
  x = ((x >>  1) & u32bitNUMBER(0x55555555)) + (x & u32bitNUMBER(0x55555555));
  x = ((x >>  2) & u32bitNUMBER(0x33333333)) + (x & u32bitNUMBER(0x33333333));
  x = ((x >>  4) & u32bitNUMBER(0x0f0f0f0f)) + (x & u32bitNUMBER(0x0f0f0f0f));
  x = ((x >>  8) & u32bitNUMBER(0x00ff00ff)) + (x & u32bitNUMBER(0x00ff00ff));
  x = ((x >> 16) & u32bitNUMBER(0x0000ffff)) + (x & u32bitNUMBER(0x0000ffff));
  return(x);
}

inline
u64bit
countNumberOfSetBits64(u64bit x) {
  x = ((x >>  1) & u64bitNUMBER(0x5555555555555555)) + (x & u64bitNUMBER(0x5555555555555555));
  x = ((x >>  2) & u64bitNUMBER(0x3333333333333333)) + (x & u64bitNUMBER(0x3333333333333333));
  x = ((x >>  4) & u64bitNUMBER(0x0f0f0f0f0f0f0f0f)) + (x & u64bitNUMBER(0x0f0f0f0f0f0f0f0f));
  x = ((x >>  8) & u64bitNUMBER(0x00ff00ff00ff00ff)) + (x & u64bitNUMBER(0x00ff00ff00ff00ff));
  x = ((x >> 16) & u64bitNUMBER(0x0000ffff0000ffff)) + (x & u64bitNUMBER(0x0000ffff0000ffff));
  x = ((x >> 32) & u64bitNUMBER(0x00000000ffffffff)) + (x & u64bitNUMBER(0x00000000ffffffff));
  return(x);
}

#endif



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

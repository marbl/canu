#ifndef FIBONACCI_ENCODING_H
#define FIBONACCI_ENCODING_H

#include "bitPacking.h"


extern u32bit fibonacciValuesLen;
extern u64bit fibonacciValues[92];



//  Routines to store and retrieve a Fibonacci encoded number to/from a
//  bit packed word array based at 'ptr' and currently at location
//  'pos'.  Both routines return the size of the encoded number in
//  'siz'.
//
void
setFibonacciEncodedNumber(u64bit *ptr,
                          u64bit  pos,
                          u64bit *siz,
                          u64bit  val);

u64bit
getFibonacciEncodedNumber(u64bit *ptr,
                          u64bit  pos,
                          u64bit *siz);





//  FibEncoding can store values up to 17,167,680,177,565 (slightly
//  below 2^45, so at most a 44-bit number) in a 64-bit quantity.
//
//  93 bits (92 + 1) are needed to store up to 64-bit values.
//
//  1000 0's is 256 bytes = 2048 bits -> 2 bits each
//  1000 1's is 376 bytes = 3008 bits -> 3 bits each
//  1000 2's is 504 bytes = 4032 bits -> 4 bits each
//  1000 3's is 504 bytes = 4032 bits -> 4 bits each
//  1000 4's is 632 bytes = 5056 bits -> 5 bits each
//


inline
void
setFibonacciEncodedNumber(u64bit *ptr,
                          u64bit  pos,
                          u64bit *siz,
                          u64bit  val) {
  u64bit  out1   = u64bitZERO;
  u64bit  out2   = u64bitZERO;
  u32bit  fib    = fibonacciValuesLen;
  u32bit  fibmax = u64bitZERO;

  //  We cannot store zero as a fibonacci number, so we simply
  //  increase everything by one.
  //
  val++;


  //  Estimate a starting point for our search; we need a function
  //  that is always slightly more than fib()
  //
  //  Find the highest bit set, do a lookup
  //
  //  XXX: Still need this!


  while (fib-- > 0) {
    if (val >= fibonacciValues[fib]) {
      if (fib >= 64)
        out2 |= u64bitONE << (127 - fib);
      else
        out1 |= u64bitONE << (63  - fib);

      val -= fibonacciValues[fib];

      if (fibmax == u64bitZERO) {
        fibmax = fib + 1;
        if (fibmax >= 64)
          out2 |= u64bitONE << (127 - fibmax);
        else
          out1 |= u64bitONE << (63  - fibmax);
      }
    }
  }

  fibmax++;

  //  Write the encoded numbers to the stream
  //
  if (fibmax > 64) {
    setDecodedValue(ptr, pos,          64, out1);
    pos += 64;
    out2 >>= (128 - fibmax);
    setDecodedValue(ptr, pos, fibmax - 64, out2);
  } else {
    out1 >>= (64 - fibmax);
    setDecodedValue(ptr, pos, fibmax,      out1);
  }

  *siz = fibmax;
}





inline
u64bit
getFibonacciEncodedNumber(u64bit *ptr,
                          u64bit  pos,
                          u64bit *siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit sft = 0x8000000000000000llu >> (pos & 0x000000000000003fllu);
  u64bit val = 0;
  u32bit fib = 0;
  u64bit newbit;
  u64bit oldbit;

  oldbit = ptr[wrd] & sft;
  sft >>= 1;
  if (sft == u64bitZERO) {
    wrd++;
    sft = 0x8000000000000000llu;
  }

  newbit = ptr[wrd] & sft;
  sft >>= 1;
  if (sft == u64bitZERO) {
    wrd++;
    sft = 0x8000000000000000llu;
  }

  while (!oldbit || !newbit) {
    if (oldbit)
      val += fibonacciValues[fib];

    fib++;

    oldbit = newbit;
    newbit = ptr[wrd] & sft;
    sft >>= 1;
    if (sft == u64bitZERO) {
      wrd++;
      sft = 0x8000000000000000llu;
    }
  }

  val += fibonacciValues[fib];

  (*siz) = fib + 2;

  //  We stored val+1, remember?  Probably not, because the encoder is
  //  next.
  //
  return(val - 1);
}


#endif  //  FIBONACCI_ENCODING_H

#ifndef BRI_BITPACKING_H
#define BRI_BITPACKING_H


extern u32bit fibonacciValuesLen;
extern u64bit fibonacciValues[92];


//  NOTE:
//
//  I assume the addresses of bits in two consectuve words would be:
//
//    [ 3 2 1 0 ][ 3 2 1 0 ]
//
//  However, setDecodedValue() works as:
//
//    bitstream = [xxxx][xx00][0000]
//    value = 0cba (three bits have information. the fourth is empty)
//
//    setDecodedValue(bitstream, value, 3);
//
//    bitstream = [xxxx][xxcb][a000]
//
//  So, if you're reading the bitstream bit-by-bit, you get the values
//  high-order first.
//



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

  //fprintf(stderr, "wrd=%lu\n", wrd);

  oldbit = ptr[wrd] & sft;
  //fprintf(stderr, "got 0x%016lx at bit=0x%016lx\n", oldbit, sft);
  sft >>= 1;
  if (sft == u64bitZERO) {
    wrd++;
    sft = 0x8000000000000000llu;
  }

  newbit = ptr[wrd] & sft;
  //fprintf(stderr, "got 0x%016lx at bit=0x%016lx\n", newbit, sft);
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
    //fprintf(stderr, "got 0x%016lx at bit=0x%016lx\n", newbit, sft);
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
      if (fib >= 64) {
        out2 |= u64bitONE << (127 - fib);
      } else {
        out1 |= u64bitONE << (63  - fib);
      }

      val -= fibonacciValues[fib];

      if (fibmax == u64bitZERO) {
        fibmax = fib + 1;
        //fprintf(stderr, "set term %lu\n", fibmax);
        if (fibmax >= 64) {
          out2 |= u64bitONE << (127 - fibmax);
        } else {
          out1 |= u64bitONE << (63  - fibmax);
        }
      }

      //fprintf(stderr, "set bit %lu\n", fib);
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
    //fprintf(stderr, "encoding 0x%016lx shifted to ", out1);
    out1 >>= (64 - fibmax);
    //fprintf(stderr, "0x%016lx\n", out1);
    setDecodedValue(ptr, pos, fibmax,      out1);
  }

  *siz = fibmax;
}



#endif  //  BRI_BITPACKING_H

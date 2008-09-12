#ifndef BRI_BITPACKING_H
#define BRI_BITPACKING_H

#include <stdio.h>
#include <assert.h>

//  Routines used for stuffing bits into a word array.

//  Define this to enable testing that the width of the data element
//  is greater than zero.  The u64bitMASK() macro (bri.h) does not
//  generate a mask for 0.  Compiler warnings are issued, because you
//  shouldn't use this in production code.
//
//#define CHECK_WIDTH

//  As CHECK_WIDTH is kind of expensive, we'll warn.
#ifdef CHECK_WIDTH
#warning libutil/bitPacking.h defined CHECK_WIDTH
#endif

//  Returns 'siz' bits from the stream based at 'ptr' and currently at
//  location 'pos'.  The position of the stream is not changed.
//
//  Retrieves a collection of values; the number of bits advanced in
//  the stream is returned.
//
//  Copies the lowest 'siz' bits in 'val' to the stream based at 'ptr'
//  and currently at 'pos'.  The position of the stream is not
//  changed.
//
//  Sets a collection of values; the number of bits advanced in the
//  stream is returned.
//
u64bit getDecodedValue (u64bit *ptr, u64bit  pos, u64bit  siz);
u64bit getDecodedValues(u64bit *ptr, u64bit  pos, u64bit  num, u64bit *sizs, u64bit *vals);
void   setDecodedValue (u64bit *ptr, u64bit  pos, u64bit  siz, u64bit  val);
u64bit setDecodedValues(u64bit *ptr, u64bit  pos, u64bit  num, u64bit *sizs, u64bit *vals);


//  Like getDecodedValue() but will pre/post increment/decrement the
//  value stored in the stream before in addition to returning the
//  value.
//
//  preIncrementDecodedValue(ptr, pos, siz) === x = getDecodedValue(ptr, pos, siz) + 1;
//                                              setDecodedValue(ptr, pos, siz, x);
//
//  preDecrementDecodedValue(ptr, pos, siz) === x = getDecodedValue(ptr, pos, siz) - 1;
//                                              setDecodedValue(ptr, pos, siz, x);
//
//  postIncrementDecodedValue(ptr, pos, siz) === x = getDecodedValue(ptr, pos, siz);
//                                               setDecodedValue(ptr, pos, siz, x + 1);
//
//  postDecrementDecodedValue(ptr, pos, siz) === x = getDecodedValue(ptr, pos, siz);
//                                               setDecodedValue(ptr, pos, siz, x - 1);
//
u64bit preIncrementDecodedValue(u64bit *ptr, u64bit  pos, u64bit  siz);
u64bit preDecrementDecodedValue(u64bit *ptr, u64bit  pos, u64bit  siz);
u64bit postIncrementDecodedValue(u64bit *ptr, u64bit  pos, u64bit  siz); 
u64bit postDecrementDecodedValue(u64bit *ptr, u64bit  pos, u64bit  siz);



//  N.B. - I assume the bits in words are big-endian, which is
//  backwards from the way we shift things around.
//
//  I define the "addresses" of bits in two consectuve words as
//  [0123][0123].  When adding words to the bit array, they're added
//  from left to right:
//
//  setDecodedValue(bitstream, %0abc, 3)
//  setDecodedValue(bitstream, %0def, 3)
//
//  results in [abcd][ef00]
//
//  But when shifting things around, we typically do it from the right
//  side, since that is where the machine places numbers.
//
//  A picture or two might help.
//
//
//         |----b1-----|
//  |-bit-||-sz-|
//         XXXXXX     
//  [0---------------63]
//         ^
//        pos
//
//
//  If the bits span two words, it'll look like this; b1 is smaller
//  than siz, and we update bit to be the "uncovered" piece of XXX
//  (all the stuff in word2).  The first word is masked, then those
//  bits are shifted onto the result in the correct place.  The second
//  word has the correct bits shifted to the right, then those are
//  appended to the result.
//
//                 |b1-|
//  |-----bit-----||---sz---|
//                 XXXXXXXXXX              
//  [0------------word1][0-------------word2]
//                 ^
//                pos
//


inline
u64bit
getDecodedValue(u64bit *ptr,
                u64bit  pos,
                u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  //PREFETCH(ptr + wrd);  makes it worse
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit ret = 0;

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: getDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("u64bitFMT")!\n", siz);
    abort();
  }
#endif

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);
  } else {
    bit  = siz - b1;
    ret  = (ptr[wrd] & u64bitMASK(b1)) << bit;
    wrd++;
    ret |= (ptr[wrd] >> (64 - bit)) & u64bitMASK(bit);
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

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: setDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("u64bitFMT")!\n", siz);
    abort();
  }
#endif

  val &= u64bitMASK(siz);

  if (b1 >= siz) {
    ptr[wrd] &= ~( u64bitMASK(siz) << (b1 - siz) );
    ptr[wrd] |= val << (b1 - siz);
  } else {
    bit = siz - b1;
    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (val & (u64bitMASK(b1) << (bit))) >> (bit);
    wrd++;
    ptr[wrd] &= ~(u64bitMASK(bit) << (64 - bit));
    ptr[wrd] |= (val & (u64bitMASK(bit))) << (64 - bit);
  }
}


inline
u64bit
getDecodedValues(u64bit *ptr,
                 u64bit  pos,
                 u64bit  num,
                 u64bit *sizs,
                 u64bit *vals) {

  //  compute the location of the start of the encoded words, then
  //  just walk through to get the remaining words.

  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  //PREFETCH(ptr + wrd);  makes it worse
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 0;

  for (u64bit i=0; i<num; i++) {
    b1 = 64 - bit;

#ifdef CHECK_WIDTH
    if (siz[i] == 0) {
      fprintf(stderr, "ERROR: postDecrementDecodedValue() called with zero size!\n");
      abort();
    }
    if (siz[i] > 64) {
      fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("u64bitFMT")!\n", siz);
      abort();
    }
#endif

    if (b1 >= sizs[i]) {
      //fprintf(stderr, "get-single pos=%d b1=%d bit=%d wrd=%d\n", pos, b1, bit, wrd);
      vals[i] = ptr[wrd] >> (b1 - sizs[i]);
      bit += sizs[i];
    } else {
      //fprintf(stderr, "get-double pos=%d b1=%d bit=%d wrd=%d bitafter=%d\n", pos, b1, bit, wrd, sizs[i]-b1);
      bit = sizs[i] - b1;
      vals[i]  = (ptr[wrd] & u64bitMASK(b1)) << bit;
      wrd++;
      vals[i] |= (ptr[wrd] >> (64 - bit)) & u64bitMASK(bit);
    }

    if (bit == 64) {
      wrd++;
      bit = 0;
    }

    assert(bit < 64);

    vals[i] &= u64bitMASK(sizs[i]);
    pos     += sizs[i];
  }

  return(pos);
}


inline
u64bit
setDecodedValues(u64bit *ptr,
                 u64bit  pos,
                 u64bit  num,
                 u64bit *sizs,
                 u64bit *vals) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 0;

  for (u64bit i=0; i<num; i++) {
    vals[i] &= u64bitMASK(sizs[i]);

    b1 = 64 - bit;

#ifdef CHECK_WIDTH
    if (siz[i] == 0) {
      fprintf(stderr, "ERROR: postDecrementDecodedValue() called with zero size!\n");
      abort();
    }
    if (siz[i] > 64) {
      fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("u64bitFMT")!\n", siz);
      abort();
    }
#endif

    if (b1 >= sizs[i]) {
      //fprintf(stderr, "set-single pos=%d b1=%d bit=%d wrd=%d\n", pos, b1, bit, wrd);
      ptr[wrd] &= ~( u64bitMASK(sizs[i]) << (b1 - sizs[i]) );
      ptr[wrd] |= vals[i] << (b1 - sizs[i]);
      bit += sizs[i];
    } else {
      //fprintf(stderr, "set-double pos=%d b1=%d bit=%d wrd=%d bitafter=%d\n", pos, b1, bit, wrd, sizs[i]-b1);
      bit = sizs[i] - b1;
      ptr[wrd] &= ~u64bitMASK(b1);
      ptr[wrd] |= (vals[i] & (u64bitMASK(b1) << (bit))) >> (bit);
      wrd++;
      ptr[wrd] &= ~(u64bitMASK(bit) << (64 - bit));
      ptr[wrd] |= (vals[i] & (u64bitMASK(bit))) << (64 - bit);
    }

    if (bit == 64) {
      wrd++;
      bit = 0;
    }

    assert(bit < 64);

    pos += sizs[i];
  }

  return(pos);
}












inline
u64bit
preIncrementDecodedValue(u64bit *ptr,
                         u64bit  pos,
                         u64bit  siz) {
  u64bit wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (pos     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bit;
  u64bit ret = 0;

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: preIncrementDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("u64bitFMT")!\n", siz);
    abort();
  }
#endif

  if (b1 >= siz) {
    ret  = ptr[wrd] >> (b1 - siz);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1 - siz) );
    ptr[wrd] |= ret << (b1 - siz);
  } else {
    bit  = siz - b1;

    ret  = (ptr[wrd] & u64bitMASK(b1)) << bit;
    ret |= (ptr[wrd+1] >> (64 - bit)) & u64bitMASK(bit);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (bit))) >> (bit);
    wrd++;
    ptr[wrd] &= ~(u64bitMASK(bit) << (64 - bit));
    ptr[wrd] |= (ret & (u64bitMASK(bit))) << (64 - bit);
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
  u64bit ret = 0;

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: preDecrementDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("u64bitFMT")!\n", siz);
    abort();
  }
#endif

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1 - siz) );
    ptr[wrd] |= ret << (b1 - siz);
  } else {
    bit  = siz - b1;

    ret  = (ptr[wrd] & u64bitMASK(b1)) << bit;
    ret |= (ptr[wrd+1] >> (64 - bit)) & u64bitMASK(bit);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (bit))) >> (bit);
    wrd++;
    ptr[wrd] &= ~(u64bitMASK(bit) << (64 - bit));
    ptr[wrd] |= (ret & (u64bitMASK(bit))) << (64 - bit);
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
  u64bit ret = 0;

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: postIncrementDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("u64bitFMT")!\n", siz);
    abort();
  }
#endif

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1 - siz) );
    ptr[wrd] |= ret << (b1 - siz);
  } else {
    bit  = siz - b1;

    ret  = (ptr[wrd] & u64bitMASK(b1)) << bit;
    ret |= (ptr[wrd+1] >> (64 - bit)) & u64bitMASK(bit);

    ret++;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (bit))) >> (bit);
    wrd++;
    ptr[wrd] &= ~(u64bitMASK(bit) << (64 - bit));
    ptr[wrd] |= (ret & (u64bitMASK(bit))) << (64 - bit);
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
  u64bit ret = 0;

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: postDecrementDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("u64bitFMT")!\n", siz);
    abort();
  }
#endif

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~( u64bitMASK(siz) << (b1 - siz) );
    ptr[wrd] |= ret << (b1 - siz);
  } else {
    bit  = siz - b1;

    ret  = (ptr[wrd] & u64bitMASK(b1)) << bit;
    ret |= (ptr[wrd+1] >> (64 - bit)) & u64bitMASK(bit);

    ret--;
    ret &= u64bitMASK(siz);

    ptr[wrd] &= ~u64bitMASK(b1);
    ptr[wrd] |= (ret & (u64bitMASK(b1) << (bit))) >> (bit);
    wrd++;
    ptr[wrd] &= ~(u64bitMASK(bit) << (64 - bit));
    ptr[wrd] |= (ret & (u64bitMASK(bit))) << (64 - bit);
  }

  ret++;
  ret &= u64bitMASK(siz);

  return(ret);
}



#endif  //  BRI_BITPACKING_H

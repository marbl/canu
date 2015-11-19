
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    kmer/libutil/bitPacking.h
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2004-APR-27 to 2004-JUN-23
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2004-MAY-06
 *      are Copyright 2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2007-DEC-05 to 2014-APR-11
 *      are Copyright 2007-2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef BRI_BITPACKING_H
#define BRI_BITPACKING_H

#include <stdio.h>
#include <assert.h>

//  Routines used for stuffing bits into a word array.

//  Define this to enable testing that the width of the data element
//  is greater than zero.  The uint64MASK() macro (bri.h) does not
//  generate a mask for 0.  Compiler warnings are issued, because you
//  shouldn't use this in production code.
//
//#define CHECK_WIDTH

//  As CHECK_WIDTH is kind of expensive, we'll warn.
#ifdef CHECK_WIDTH
#warning CHECK_WIDTH is EXPENSIVE
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
uint64 getDecodedValue (uint64 *ptr, uint64  pos, uint64  siz);
uint64 getDecodedValues(uint64 *ptr, uint64  pos, uint64  num, uint64 *sizs, uint64 *vals);
void   setDecodedValue (uint64 *ptr, uint64  pos, uint64  siz, uint64  val);
uint64 setDecodedValues(uint64 *ptr, uint64  pos, uint64  num, uint64 *sizs, uint64 *vals);


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
uint64 preIncrementDecodedValue(uint64 *ptr, uint64  pos, uint64  siz);
uint64 preDecrementDecodedValue(uint64 *ptr, uint64  pos, uint64  siz);
uint64 postIncrementDecodedValue(uint64 *ptr, uint64  pos, uint64  siz);
uint64 postDecrementDecodedValue(uint64 *ptr, uint64  pos, uint64  siz);



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
uint64
getDecodedValue(uint64 *ptr,
                uint64  pos,
                uint64  siz) {
  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  //PREFETCH(ptr + wrd);  makes it worse
  uint64 bit = (pos     ) & 0x000000000000003fllu;
  uint64 b1  = 64 - bit;
  uint64 ret = 0;

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: getDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("uint64FMT")!\n", siz);
    abort();
  }
#endif

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);
  } else {
    bit  = siz - b1;
    ret  = (ptr[wrd] & uint64MASK(b1)) << bit;
    wrd++;
    ret |= (ptr[wrd] >> (64 - bit)) & uint64MASK(bit);
  }

  ret &= uint64MASK(siz);

  return(ret);
}


inline
void
setDecodedValue(uint64 *ptr,
                uint64  pos,
                uint64  siz,
                uint64  val) {
  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (pos     ) & 0x000000000000003fllu;
  uint64 b1  = 64 - bit;

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: setDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("uint64FMT")!\n", siz);
    abort();
  }
#endif

  val &= uint64MASK(siz);

  if (b1 >= siz) {
    ptr[wrd] &= ~( uint64MASK(siz) << (b1 - siz) );
    ptr[wrd] |= val << (b1 - siz);
  } else {
    bit = siz - b1;
    ptr[wrd] &= ~uint64MASK(b1);
    ptr[wrd] |= (val & (uint64MASK(b1) << (bit))) >> (bit);
    wrd++;
    ptr[wrd] &= ~(uint64MASK(bit) << (64 - bit));
    ptr[wrd] |= (val & (uint64MASK(bit))) << (64 - bit);
  }
}


inline
uint64
getDecodedValues(uint64 *ptr,
                 uint64  pos,
                 uint64  num,
                 uint64 *sizs,
                 uint64 *vals) {

  //  compute the location of the start of the encoded words, then
  //  just walk through to get the remaining words.

  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  //PREFETCH(ptr + wrd);  makes it worse
  uint64 bit = (pos     ) & 0x000000000000003fllu;
  uint64 b1  = 0;

  for (uint64 i=0; i<num; i++) {
    b1 = 64 - bit;

#ifdef CHECK_WIDTH
    if (siz[i] == 0) {
      fprintf(stderr, "ERROR: postDecrementDecodedValue() called with zero size!\n");
      abort();
    }
    if (siz[i] > 64) {
      fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("uint64FMT")!\n", siz);
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
      vals[i]  = (ptr[wrd] & uint64MASK(b1)) << bit;
      wrd++;
      vals[i] |= (ptr[wrd] >> (64 - bit)) & uint64MASK(bit);
    }

    if (bit == 64) {
      wrd++;
      bit = 0;
    }

    assert(bit < 64);

    vals[i] &= uint64MASK(sizs[i]);
    pos     += sizs[i];
  }

  return(pos);
}


inline
uint64
setDecodedValues(uint64 *ptr,
                 uint64  pos,
                 uint64  num,
                 uint64 *sizs,
                 uint64 *vals) {
  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (pos     ) & 0x000000000000003fllu;
  uint64 b1  = 0;

  for (uint64 i=0; i<num; i++) {
    vals[i] &= uint64MASK(sizs[i]);

    b1 = 64 - bit;

#ifdef CHECK_WIDTH
    if (siz[i] == 0) {
      fprintf(stderr, "ERROR: postDecrementDecodedValue() called with zero size!\n");
      abort();
    }
    if (siz[i] > 64) {
      fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("uint64FMT")!\n", siz);
      abort();
    }
#endif

    if (b1 >= sizs[i]) {
      //fprintf(stderr, "set-single pos=%d b1=%d bit=%d wrd=%d\n", pos, b1, bit, wrd);
      ptr[wrd] &= ~( uint64MASK(sizs[i]) << (b1 - sizs[i]) );
      ptr[wrd] |= vals[i] << (b1 - sizs[i]);
      bit += sizs[i];
    } else {
      //fprintf(stderr, "set-double pos=%d b1=%d bit=%d wrd=%d bitafter=%d\n", pos, b1, bit, wrd, sizs[i]-b1);
      bit = sizs[i] - b1;
      ptr[wrd] &= ~uint64MASK(b1);
      ptr[wrd] |= (vals[i] & (uint64MASK(b1) << (bit))) >> (bit);
      wrd++;
      ptr[wrd] &= ~(uint64MASK(bit) << (64 - bit));
      ptr[wrd] |= (vals[i] & (uint64MASK(bit))) << (64 - bit);
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
uint64
preIncrementDecodedValue(uint64 *ptr,
                         uint64  pos,
                         uint64  siz) {
  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (pos     ) & 0x000000000000003fllu;
  uint64 b1  = 64 - bit;
  uint64 ret = 0;

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: preIncrementDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("uint64FMT")!\n", siz);
    abort();
  }
#endif

  if (b1 >= siz) {
    ret  = ptr[wrd] >> (b1 - siz);

    ret++;
    ret &= uint64MASK(siz);

    ptr[wrd] &= ~( uint64MASK(siz) << (b1 - siz) );
    ptr[wrd] |= ret << (b1 - siz);
  } else {
    bit  = siz - b1;

    ret  = (ptr[wrd] & uint64MASK(b1)) << bit;
    ret |= (ptr[wrd+1] >> (64 - bit)) & uint64MASK(bit);

    ret++;
    ret &= uint64MASK(siz);

    ptr[wrd] &= ~uint64MASK(b1);
    ptr[wrd] |= (ret & (uint64MASK(b1) << (bit))) >> (bit);
    wrd++;
    ptr[wrd] &= ~(uint64MASK(bit) << (64 - bit));
    ptr[wrd] |= (ret & (uint64MASK(bit))) << (64 - bit);
  }

  return(ret);
}



inline
uint64
preDecrementDecodedValue(uint64 *ptr,
                         uint64  pos,
                         uint64  siz) {
  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (pos     ) & 0x000000000000003fllu;
  uint64 b1  = 64 - bit;
  uint64 ret = 0;

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: preDecrementDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("uint64FMT")!\n", siz);
    abort();
  }
#endif

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret--;
    ret &= uint64MASK(siz);

    ptr[wrd] &= ~( uint64MASK(siz) << (b1 - siz) );
    ptr[wrd] |= ret << (b1 - siz);
  } else {
    bit  = siz - b1;

    ret  = (ptr[wrd] & uint64MASK(b1)) << bit;
    ret |= (ptr[wrd+1] >> (64 - bit)) & uint64MASK(bit);

    ret--;
    ret &= uint64MASK(siz);

    ptr[wrd] &= ~uint64MASK(b1);
    ptr[wrd] |= (ret & (uint64MASK(b1) << (bit))) >> (bit);
    wrd++;
    ptr[wrd] &= ~(uint64MASK(bit) << (64 - bit));
    ptr[wrd] |= (ret & (uint64MASK(bit))) << (64 - bit);
  }

  return(ret);
}



inline
uint64
postIncrementDecodedValue(uint64 *ptr,
                          uint64  pos,
                          uint64  siz) {
  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (pos     ) & 0x000000000000003fllu;
  uint64 b1  = 64 - bit;
  uint64 ret = 0;

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: postIncrementDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("uint64FMT")!\n", siz);
    abort();
  }
#endif

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret++;
    ret &= uint64MASK(siz);

    ptr[wrd] &= ~( uint64MASK(siz) << (b1 - siz) );
    ptr[wrd] |= ret << (b1 - siz);
  } else {
    bit  = siz - b1;

    ret  = (ptr[wrd] & uint64MASK(b1)) << bit;
    ret |= (ptr[wrd+1] >> (64 - bit)) & uint64MASK(bit);

    ret++;
    ret &= uint64MASK(siz);

    ptr[wrd] &= ~uint64MASK(b1);
    ptr[wrd] |= (ret & (uint64MASK(b1) << (bit))) >> (bit);
    wrd++;
    ptr[wrd] &= ~(uint64MASK(bit) << (64 - bit));
    ptr[wrd] |= (ret & (uint64MASK(bit))) << (64 - bit);
  }

  ret--;
  ret &= uint64MASK(siz);

  return(ret);
}





inline
uint64
postDecrementDecodedValue(uint64 *ptr,
                          uint64  pos,
                          uint64  siz) {
  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (pos     ) & 0x000000000000003fllu;
  uint64 b1  = 64 - bit;
  uint64 ret = 0;

#ifdef CHECK_WIDTH
  if (siz == 0) {
    fprintf(stderr, "ERROR: postDecrementDecodedValue() called with zero size!\n");
    abort();
  }
  if (siz > 64) {
    fprintf(stderr, "ERROR: getDecodedValue() called with huge size ("uint64FMT")!\n", siz);
    abort();
  }
#endif

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret--;
    ret &= uint64MASK(siz);

    ptr[wrd] &= ~( uint64MASK(siz) << (b1 - siz) );
    ptr[wrd] |= ret << (b1 - siz);
  } else {
    bit  = siz - b1;

    ret  = (ptr[wrd] & uint64MASK(b1)) << bit;
    ret |= (ptr[wrd+1] >> (64 - bit)) & uint64MASK(bit);

    ret--;
    ret &= uint64MASK(siz);

    ptr[wrd] &= ~uint64MASK(b1);
    ptr[wrd] |= (ret & (uint64MASK(b1) << (bit))) >> (bit);
    wrd++;
    ptr[wrd] &= ~(uint64MASK(bit) << (64 - bit));
    ptr[wrd] |= (ret & (uint64MASK(bit))) << (64 - bit);
  }

  ret++;
  ret &= uint64MASK(siz);

  return(ret);
}



#endif  //  BRI_BITPACKING_H

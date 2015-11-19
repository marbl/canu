
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
 *  Modifications by:
 *
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "bitPacking.H"
#include "bitOperations.H"

//  Methods to encode and decode numbers from an array of bits.
//
//  The bit packed array is rooted at location 'ptr' and currently at location
//  'pos'.  The size of the encoded number (in bits) is returned in 'siz'.
//
//  Came from kmer/libutil/
//    eliasDeltaEncoding.h
//    eliasGammaEncoding.h
//    fibonacciEncoding.h
//    generalizedUnaryEncoding.h
//    unaryEncoding.h
//    unaryEncodingTester.C


//
//  Unary Encoding
//
//  Routines to store and retrieve a unary encoded number to/from a
//  bit packed word array based at 'ptr' and currently at location
//  'pos'.  Both routines return the size of the encoded number in
//  'siz'.
//
//  The usual unary encoding.  Store the number n as n 0 bits followed
//  by a single 1 bit.
//
//  0 -> 1
//  1 -> 01
//  2 -> 001
//  3 -> 0001
//  4 -> 00001
//
//  See the decoder as to why we use 0 instead of 1 for the count.
//

inline
void
setUnaryEncodedNumber(uint64 *ptr,
                      uint64  pos,
                      uint64 *siz,
                      uint64  val) {

  *siz = val + 1;

  while (val >= 64) {
    setDecodedValue(ptr, pos, 64, uint64ZERO);
    pos += 64;
    val -= 64;
    siz += 64;
  }

  setDecodedValue(ptr, pos, val + 1, uint64ONE);
  pos += val + 1;
}



inline
uint64
getUnaryEncodedNumber(uint64 *ptr,
                      uint64  pos,
                      uint64 *siz) {
  uint64 val = uint64ZERO;
  uint64 enc = uint64ZERO;

  //  How many whole words are zero?
  //
  enc = getDecodedValue(ptr, pos, 64);
  while (enc == uint64ZERO) {
    val += 64;
    pos += 64;
    enc  = getDecodedValue(ptr, pos, 64);
  }

  //  This word isn't zero.  Count how many bits are zero (see, the
  //  choice of 0 or 1 for the encoding wasn't arbitrary!)

  val += 64 - logBaseTwo64(enc);

  *siz = val + 1;

  return(val);
}



//
//  Generalized Unary Encoding
//
//  Generalized unary encodings.  Defined by (start, step, stop).
//  This implementation uses stop=infinity to encode all possible
//  numbers.  If you know the highest number possible, you'll get a
//  slight decrease in space used ...
//
//  The method:
//
//  The mth code word consists of 'm' unary encoded, followed by w =
//  start + m * step binary encoded bits.  If a == stop, then the
//  terminator in the unary code is dropped.
//
//  Encoding is tricky.  Take the 3,2,9 example:
//    m  w  template    # vals     #'s
//    0  3  1xxx             8    0-  7
//    1  5  01xxxxx         32    8- 39
//    2  7  001xxxxxxx     128   40-167
//    3  9  000xxxxxxxxx   512  168-679
//
//  I don't see a nice way of mapping our number n to the prefix m,
//  short of some sort of search.  The implementation below is
//  probably very slow.
//
//  On the bright side, decoding is trivial.  Read back the unary
//  encoded number, then read that many bits to get the value.
//


static const uint64 _genunary_start = 3;
static const uint64 _genunary_step  = 2;


inline
void
setGeneralizedUnaryEncodedNumber(uint64 *ptr,
                                 uint64  pos,
                                 uint64 *siz,
                                 uint64  val) {
  uint64 m = uint64ZERO;
  uint64 w = _genunary_start;
  uint64 n = uint64ONE << w;

  //  Search for the prefix m, given our number 'val'.
  //  While doing this, we get rid of all the implicitly stored values from 'val'.
  //
#ifdef DEBUG_GENERALIZEDUNARYENCODING
  fprintf(stderr, "  val="uint64FMT" try n="uint64FMT" for m="uint64FMT"\n", val, n, m);
#endif

  while (n <= val) {
    val -= n;
    w   += _genunary_step;
    n    = uint64ONE << w;
    m++;
#ifdef DEBUG_GENERALIZEDUNARYENCODING
    fprintf(stderr, "  val="uint64FMT" try n="uint64FMT" for m="uint64FMT"\n", val, n, m);
#endif
  }

#ifdef DEBUG_GENERALIZEDUNARYENCODING
  fprintf(stderr, "val="uint64FMT" found m="uint64FMT"\n", val, m);
#endif

  //  Now just encode the number
  //    m    - the unary encoded prefix
  //    w    - the size of the binary encoded number

  setUnaryEncodedNumber(ptr, pos, siz, m);
  setDecodedValue(ptr, pos+*siz, w, val);
  *siz = m + 1 + w;
}



inline
uint64
getGeneralizedUnaryEncodedNumber(uint64 *ptr,
                                 uint64  pos,
                                 uint64 *siz) {
  uint64 val = uint64ZERO;
  uint64 m   = uint64ZERO;
  uint64 w   = uint64ZERO;

  //  Comments in the encoder apply here too.

  m    = getUnaryEncodedNumber(ptr, pos, siz);
  w    = _genunary_start + m * _genunary_step;
  val  = getDecodedValue(ptr, pos + *siz, w);
  *siz = m + 1 + w;

#ifdef DEBUG_GENERALIZEDUNARYENCODING
  fprintf(stderr, "m="uint64FMT" w="uint64FMT" val="uint64FMT"\n", m, w, val);
#endif

  //  Add in the implcitly stored pieces of the number

  while (m--) {
    w -= _genunary_step;
    val += uint64ONE << w;
  }

  return(val);
}


//
//  Elias' Gamma Encoding
//

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

//
//  Elias' Delta Encoding
//

inline
void
setEliasDeltaEncodedNumber(uint64 *ptr,
                           uint64  pos,
                           uint64 *siz,
                           uint64  val) {
  uint64 b = logBaseTwo64(val);
  setEliasGammaEncodedNumber(ptr, pos, siz, b);
  pos += *siz;
  setDecodedValue(ptr, pos, b-1, val);
  *siz += b-1;
}


inline
uint64
getEliasDeltaEncodedNumber(uint64 *ptr,
                           uint64  pos,
                           uint64 *siz) {
  uint64 b = getEliasGammaEncodedNumber(ptr, pos, siz) - 1;
  pos  += *siz;
  *siz += b;
  return(uint64ONE << b | getDecodedValue(ptr, pos, b));
}




//
//  Fibonacci Encoding
//
//  FibEncoding can store values up to 17,167,680,177,565 (slightly
//  below 2^45, so at most a 44-bit number) in a 64-bit quantity.
//  93 bits (92 + 1) are needed to store up to 64-bit values.
//  Since zero cannot be stored, the value is incremented by one before storing.
//  The actual space used is:
//
//    ####  bits
//       0  2
//       1  3
//       2  4
//       3  4
//       4  5
//       5  5
//       6  5
//       7  6
//       8  6
//       9  6
//      10  6
//      11  6
//      12  7
//      20  8
//      33  9
//      54  10
//      88  11
//     143  12
//     232  13
//     376  14
//     609  15
//     986  16
//    1596  17
//    2583  18
//    4180  19
//    6764  20
//   10945  21
//   17710  22
//   28656  23
//   46387  24
//   75024  25
//  121392  26

extern const uint32 fibonacciValuesLen;
extern const uint64 fibonacciValues[92];

inline
void
setFibonacciEncodedNumber(uint64 *ptr,
                          uint64  pos,
                          uint64 *siz,
                          uint64  val) {
  uint64  out1   = uint64ZERO;
  uint64  out2   = uint64ZERO;
  uint32  fib    = fibonacciValuesLen;
  uint32  fibmax = uint64ZERO;

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
        out2 |= uint64ONE << (127 - fib);
      else
        out1 |= uint64ONE << (63  - fib);

      val -= fibonacciValues[fib];

      if (fibmax == uint64ZERO) {
        fibmax = fib + 1;
        if (fibmax >= 64)
          out2 |= uint64ONE << (127 - fibmax);
        else
          out1 |= uint64ONE << (63  - fibmax);
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
uint64
getFibonacciEncodedNumber(uint64 *ptr,
                          uint64  pos,
                          uint64 *siz) {
  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  uint64 sft = 0x8000000000000000llu >> (pos & 0x000000000000003fllu);
  uint64 val = 0;
  uint32 fib = 0;
  uint64 newbit;
  uint64 oldbit;

  oldbit = ptr[wrd] & sft;
  sft >>= 1;
  if (sft == uint64ZERO) {
    wrd++;
    sft = 0x8000000000000000llu;
  }

  newbit = ptr[wrd] & sft;
  sft >>= 1;
  if (sft == uint64ZERO) {
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
    if (sft == uint64ZERO) {
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

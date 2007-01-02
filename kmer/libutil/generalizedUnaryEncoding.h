#ifndef GENERALIZED_UNARY_ENCODING_H
#define GENERALIZED_UNARY_ENCODING_H

#include "bitPacking.h"

//  Lots and lots of semi-useless debugging information
//#define DEBUG_GENERALIZEDUNARYENCODING


//  Generalized unary encodings.  Defined by (start, step, stop).
//  This implementation uses stop=infinity to encode all possible
//  numbers.  If you know the highest number possible, you'll get a
//  slight decrease in space used ...

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

static const u64bit _genunary_start = 3;
static const u64bit _genunary_step  = 2;
//static const u64bit _genunary_stop  = ~u64bitZERO;


inline
void
setGeneralizedUnaryEncodedNumber(u64bit *ptr,
                                 u64bit  pos,
                                 u64bit *siz,
                                 u64bit  val) {
  u64bit m = u64bitZERO;
  u64bit w = _genunary_start;
  u64bit n = u64bitONE << w;

  //  Search for the prefix m, given our number 'val'.
  //  While doing this, we get rid of all the implicitly stored values from 'val'.
  //
#ifdef DEBUG_GENERALIZEDUNARYENCODING
  fprintf(stderr, "  val="u64bitFMT" try n="u64bitFMT" for m="u64bitFMT"\n", val, n, m);
#endif

  while (n <= val) {
    val -= n;
    w   += _genunary_step;
    n    = u64bitONE << w;
    m++;
#ifdef DEBUG_GENERALIZEDUNARYENCODING
    fprintf(stderr, "  val="u64bitFMT" try n="u64bitFMT" for m="u64bitFMT"\n", val, n, m);
#endif
  }

#ifdef DEBUG_GENERALIZEDUNARYENCODING
  fprintf(stderr, "val="u64bitFMT" found m="u64bitFMT"\n", val, m);
#endif

  //  Now just encode the number
  //    m    - the unary encoded prefix
  //    w    - the size of the binary encoded number

  setUnaryEncodedNumber(ptr, pos, siz, m);
  setDecodedValue(ptr, pos+*siz, w, val);
  *siz = m + 1 + w;
}



inline
u64bit
getGeneralizedUnaryEncodedNumber(u64bit *ptr,
                                 u64bit  pos,
                                 u64bit *siz) {
  u64bit val = u64bitZERO;
  u64bit m   = u64bitZERO;
  u64bit w   = u64bitZERO;

  //  Comments in the encoder apply here too.

  m    = getUnaryEncodedNumber(ptr, pos, siz);
  w    = _genunary_start + m * _genunary_step;
  val  = getDecodedValue(ptr, pos + *siz, w);
  *siz = m + 1 + w;

#ifdef DEBUG_GENERALIZEDUNARYENCODING
  fprintf(stderr, "m="u64bitFMT" w="u64bitFMT" val="u64bitFMT"\n", m, w, val);
#endif

  //  Add in the implcitly stored pieces of the number
  //
  while (m--) {
    w -= _genunary_step;
    val += u64bitONE << w;
  }

  return(val);
}




#endif  //  GENERALIZED_UNARY_ENCODING_H

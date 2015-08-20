
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
 *    Brian P. Walenz on 2004-APR-27
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2007-JAN-02 to 2014-APR-11
 *      are Copyright 2007,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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

static const uint64 _genunary_start = 3;
static const uint64 _genunary_step  = 2;
//static const uint64 _genunary_stop  = ~uint64ZERO;


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
  //
  while (m--) {
    w -= _genunary_step;
    val += uint64ONE << w;
  }

  return(val);
}




#endif  //  GENERALIZED_UNARY_ENCODING_H

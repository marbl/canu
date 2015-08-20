
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
 *    Brian P. Walenz from 2008-JUN-09 to 2014-APR-11
 *      are Copyright 2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef FIBONACCI_ENCODING_H
#define FIBONACCI_ENCODING_H

#include "bitPacking.h"

//  Routines to store and retrieve a Fibonacci encoded number to/from a
//  bit packed word array based at 'ptr' and currently at location
//  'pos'.  Both routines return the size of the encoded number in
//  'siz'.
//
//  FibEncoding can store values up to 17,167,680,177,565 (slightly
//  below 2^45, so at most a 44-bit number) in a 64-bit quantity.
//
//  93 bits (92 + 1) are needed to store up to 64-bit values.
//
//  Remember that since we can't store 0, we increment all incoming
//  values, so the actual space used is:
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

extern uint32 fibonacciValuesLen;
extern uint64 fibonacciValues[92];

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


#endif  //  FIBONACCI_ENCODING_H

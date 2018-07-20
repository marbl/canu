
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
 *    kmer/libutil/bitOperations.h
 *
 *  Modifications by:
 *
 *    Brian P. Walenz on 2004-APR-27
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-FEB-23 to 2014-APR-11
 *      are Copyright 2005-2006,2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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
//      uint32 r = uint32ZERO;
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
uint64
reverseBits64(uint64 x) {
  x = ((x >>  1) & uint64NUMBER(0x5555555555555555)) | ((x <<  1) & uint64NUMBER(0xaaaaaaaaaaaaaaaa));
  x = ((x >>  2) & uint64NUMBER(0x3333333333333333)) | ((x <<  2) & uint64NUMBER(0xcccccccccccccccc));
  x = ((x >>  4) & uint64NUMBER(0x0f0f0f0f0f0f0f0f)) | ((x <<  4) & uint64NUMBER(0xf0f0f0f0f0f0f0f0));
  x = ((x >>  8) & uint64NUMBER(0x00ff00ff00ff00ff)) | ((x <<  8) & uint64NUMBER(0xff00ff00ff00ff00));
  x = ((x >> 16) & uint64NUMBER(0x0000ffff0000ffff)) | ((x << 16) & uint64NUMBER(0xffff0000ffff0000));
  x = ((x >> 32) & uint64NUMBER(0x00000000ffffffff)) | ((x << 32) & uint64NUMBER(0xffffffff00000000));
  return(x);
}

inline
uint32
reverseBits32(uint32 x) {
  x = ((x >>  1) & uint32NUMBER(0x55555555)) | ((x <<  1) & uint32NUMBER(0xaaaaaaaa));
  x = ((x >>  2) & uint32NUMBER(0x33333333)) | ((x <<  2) & uint32NUMBER(0xcccccccc));
  x = ((x >>  4) & uint32NUMBER(0x0f0f0f0f)) | ((x <<  4) & uint32NUMBER(0xf0f0f0f0));
  x = ((x >>  8) & uint32NUMBER(0x00ff00ff)) | ((x <<  8) & uint32NUMBER(0xff00ff00));
  x = ((x >> 16) & uint32NUMBER(0x0000ffff)) | ((x << 16) & uint32NUMBER(0xffff0000));
  return(x);
}



inline
uint64
uint64Swap(uint64 x) {
  x = ((x >>  8) & uint64NUMBER(0x00ff00ff00ff00ff)) | ((x <<  8) & uint64NUMBER(0xff00ff00ff00ff00));
  x = ((x >> 16) & uint64NUMBER(0x0000ffff0000ffff)) | ((x << 16) & uint64NUMBER(0xffff0000ffff0000));
  x = ((x >> 32) & uint64NUMBER(0x00000000ffffffff)) | ((x << 32) & uint64NUMBER(0xffffffff00000000));
  return(x);
}

inline
uint32
uint32Swap(uint32 x) {
  x = ((x >>  8) & uint32NUMBER(0x00ff00ff)) | ((x <<  8) & uint32NUMBER(0xff00ff00));
  x = ((x >> 16) & uint32NUMBER(0x0000ffff)) | ((x << 16) & uint32NUMBER(0xffff0000));
  return(x);
}

inline
uint16
uint16Swap(uint16 x) {
  x = ((x >>  8) & 0x00ff) | ((x <<  8) & 0xff00);
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
uint32
countNumberOfSetBits32(uint32 x) {
  return(__builtin_popcount(x));
}

inline
uint64
countNumberOfSetBits64(uint64 x) {
  return(__builtin_popcountll(x));
}

#else

inline
uint32
countNumberOfSetBits32(uint32 x) {
  x = ((x >>  1) & uint32NUMBER(0x55555555)) + (x & uint32NUMBER(0x55555555));
  x = ((x >>  2) & uint32NUMBER(0x33333333)) + (x & uint32NUMBER(0x33333333));
  x = ((x >>  4) & uint32NUMBER(0x0f0f0f0f)) + (x & uint32NUMBER(0x0f0f0f0f));
  x = ((x >>  8) & uint32NUMBER(0x00ff00ff)) + (x & uint32NUMBER(0x00ff00ff));
  x = ((x >> 16) & uint32NUMBER(0x0000ffff)) + (x & uint32NUMBER(0x0000ffff));
  return(x);
}

inline
uint64
countNumberOfSetBits64(uint64 x) {
  x = ((x >>  1) & uint64NUMBER(0x5555555555555555)) + (x & uint64NUMBER(0x5555555555555555));
  x = ((x >>  2) & uint64NUMBER(0x3333333333333333)) + (x & uint64NUMBER(0x3333333333333333));
  x = ((x >>  4) & uint64NUMBER(0x0f0f0f0f0f0f0f0f)) + (x & uint64NUMBER(0x0f0f0f0f0f0f0f0f));
  x = ((x >>  8) & uint64NUMBER(0x00ff00ff00ff00ff)) + (x & uint64NUMBER(0x00ff00ff00ff00ff));
  x = ((x >> 16) & uint64NUMBER(0x0000ffff0000ffff)) + (x & uint64NUMBER(0x0000ffff0000ffff));
  x = ((x >> 32) & uint64NUMBER(0x00000000ffffffff)) + (x & uint64NUMBER(0x00000000ffffffff));
  return(x);
}

#endif



inline
uint32
logBaseTwo32(uint32 x) {
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return(countNumberOfSetBits32(x));
}

inline
uint64
logBaseTwo64(uint64 x) {
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  return(countNumberOfSetBits64(x));
}




#endif  //  BRI_BITS_H

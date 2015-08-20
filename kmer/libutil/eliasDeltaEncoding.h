
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
 *    Brian P. Walenz on 2014-APR-11
 *      are Copyright 2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef ELIAS_DELTA_ENCODING_H
#define ELIAS_DELTA_ENCODING_H

#include "bitPacking.h"

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



#endif  //  ELIAS_DELTA_ENCODING_H

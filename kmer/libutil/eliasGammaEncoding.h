
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

#ifndef ELIAS_GAMMA_ENCODING_H
#define ELIAS_GAMMA_ENCODING_H

#include "bitPacking.h"

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



#endif  //  ELIAS_GAMMA_ENCODING_H

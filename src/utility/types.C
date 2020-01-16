
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "types.H"

char   hex[16]     = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };
char   str[16][33] = { 0 };
uint32 pos         =   0;

char const *
toHex(uint128 v) {
  char *ret = str[pos++];

  if (pos >= 16)
    pos = 0;

  uint32  p = 32;
  uint32  s = 0;

  while (p > 0) {
    p -= 1;
    ret[p] = hex[ (v >> s) & 0xf ];
    s += 4;
  }

  ret[32] = 0;

  return(ret);
}

char const *
toHex(uint64  v) {
  char *ret = str[pos++];

  if (pos >= 16)
    pos = 0;

  uint32  p = 16;
  uint32  s = 0;

  while (p > 0) {
    p -= 1;
    ret[p] = hex[ (v >> s) & 0xf ];
    s += 4;
  }

  ret[16] = 0;

  return(ret);
}

char const *
toHex(uint32  v) {
  char *ret = str[pos++];

  if (pos >= 16)
    pos = 0;

  uint32  p = 8;
  uint32  s = 0;

  while (p > 0) {
    p -= 1;
    ret[p] = hex[ (v >> s) & 0xf ];
    s += 4;
  }

  ret[8] = 0;

  return(ret);
}


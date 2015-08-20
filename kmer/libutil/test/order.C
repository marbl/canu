
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
 *    Brian P. Walenz from 2005-DEC-04 to 2014-APR-11
 *      are Copyright 2005,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>

#include "../util.h"

//#include <sys/param.h>

union u64 {
  uint64          u;
  unsigned char   c[8];
};

union u32 {
  uint32          u;
  unsigned char   c[4];
};

union u16 {
  uint16          u;
  unsigned char   c[2];
};


uint64
uint64Swap(uint64 x) {
  x = ((x >>  8) & uint64NUMBER(0x00ff00ff00ff00ff)) | ((x <<  8) & uint64NUMBER(0xff00ff00ff00ff00));
  x = ((x >> 16) & uint64NUMBER(0x0000ffff0000ffff)) | ((x << 16) & uint64NUMBER(0xffff0000ffff0000));
  x = ((x >> 32) & uint64NUMBER(0x00000000ffffffff)) | ((x << 32) & uint64NUMBER(0xffffffff00000000));
  return(x);
}

uint32
uint32Swap(uint32 x) {
  x = ((x >>  8) & uint32NUMBER(0x00ff00ff)) | ((x <<  8) & uint32NUMBER(0xff00ff00));
  x = ((x >> 16) & uint32NUMBER(0x0000ffff)) | ((x << 16) & uint32NUMBER(0xffff0000));
  return(x);
}

uint16
uint16Swap(uint16 x) {
  x = ((x >>  8) & 0x00ff) | ((x <<  8) & 0xff00);
  return(x);
}



int
main(int argc, char **argv) {
  u64  u64v;
  u32  u32v;
  u16  u16v;

  u64v.u = 0x1234567890abcdefLLU;
  u32v.u = 0x12345678;
  u16v.u = 0x1234;

  for (int i=0; i<8; i++)
    fprintf(stderr, "%02x", u64v.c[i]);
  fprintf(stderr, "\n");

  for (int i=0; i<4; i++)
    fprintf(stderr, "%02x", u32v.c[i]);
  fprintf(stderr, "\n");

  for (int i=0; i<2; i++)
    fprintf(stderr, "%02x", u16v.c[i]);
  fprintf(stderr, "\n");

  u64v.u = uint64Swap(u64v.u);
  u32v.u = uint32Swap(u32v.u);
  u16v.u = uint16Swap(u16v.u);

  for (int i=0; i<8; i++)
    fprintf(stderr, "%02x", u64v.c[i]);
  fprintf(stderr, "\n");

  for (int i=0; i<4; i++)
    fprintf(stderr, "%02x", u32v.c[i]);
  fprintf(stderr, "\n");

  for (int i=0; i<2; i++)
    fprintf(stderr, "%02x", u16v.c[i]);
  fprintf(stderr, "\n");
}


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
 *    Brian P. Walenz beginning on 2016-MAR-10
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

typedef int8_t  int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t  uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

#include "intervalList.H"

int
main(int argc, char **argv) {

  intervalList<int32>  t1;

  t1.add(0, 10);
  t1.add(11,7);
  t1.add(20, 8);

  fprintf(stderr, "BEFORE:\n");
  for (uint32 ii=0; ii<t1.numberOfIntervals(); ii++)
    fprintf(stderr, "%2d %3d-%3d\n", ii, t1.lo(ii), t1.hi(ii));

  t1.merge(-1);

  fprintf(stderr, "AFTER:\n");
  for (uint32 ii=0; ii<t1.numberOfIntervals(); ii++)
    fprintf(stderr, "%2d %3d-%3d\n", ii, t1.lo(ii), t1.hi(ii));

  exit(0);
}

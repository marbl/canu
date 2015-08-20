
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
 *    Brian P. Walenz on 2004-MAY-06
 *      are Copyright 2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2004-OCT-10
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

#include <stdio.h>

#include "util.h"

int
main(void) {
  uint32   errors = 0;
  uint32   u3 = -1;
  int32   s3 = -1;
  uint64   u6 = -1;
  int64   s6 = -1;

  if (sizeof(uint32) != 4)
    fprintf(stderr, "uint32 has %d bytes (should be 4)!\n", (int)sizeof(uint32)), errors++;

  if (sizeof(uint64) != 8)
    fprintf(stderr, "uint64 has %d bytes (should be 8)!\n", (int)sizeof(uint64)), errors++;

  if (u3 < 0)
    fprintf(stderr, "uint32 is signed (should be unsigned)!\n"), errors++;

  if (s3 > 0)
    fprintf(stderr, "int32 is unsigned (should be signed)!\n"), errors++;

  if (u6 < 0)
    fprintf(stderr, "uint64 is signed (should be unsigned)!\n"), errors++;

  if (s6 > 0)
    fprintf(stderr, "int64 is unsigned (should be signed)!\n"), errors++;

  return(errors);
}



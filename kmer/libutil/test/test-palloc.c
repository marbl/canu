
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>

#include "util.h"

int
main(int argc, char **argv) {
  mt_s   *mtctx;
  int     i;

  psetdebug(2);
  psetblocksize(1024);

  palloc(2048);
  palloc(128);
  palloc(999);
  palloc(1);
  palloc(2);
  palloc(3);
  palloc(4);
  palloc(2056);
  palloc(8);
  palloc(2064);
  palloc(8);
  palloc(2072);
  palloc(8);

  pdumppalloc();

  pfree();

  fprintf(stderr, "----------------------------------------\n");

  psetblocksize(10240);

  palloc(2048);
  palloc(128);
  palloc(999);
  palloc(8);
  palloc(8);
  palloc(8);
  palloc(8);
  palloc(2056);
  palloc(8);
  palloc(2064);
  palloc(8);
  palloc(2072);
  palloc(8);

  pdumppalloc();

  pfree();

  psetdebug(0);
  psetblocksize(16 * 1024 * 1024);

  mtctx = mtInit(time(NULL));
  for (i=0; i<512 * 1024; i++)
    palloc(mtRandom32(mtctx) & 0xfff);
  psetdebug(1);
  pfree();

  return(0);
}


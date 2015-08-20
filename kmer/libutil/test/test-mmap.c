
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
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

#include "util.h"

//  Does a quick test of memory mapped files.  First, writes a small
//  file, then it reads it back, checking the data.
//
//  Takes one optional argument, the size in MB of the file to map.

int
main(int argc, char **argv) {
  size_t   lw;
  uint32  *ww = 0L;
  uint32   idx = 0;
  uint32   err = 0;
  FILE    *out;
  uint32   blockSize = 1048576;
  uint32   numBlocks = 32;

  if (argc == 2)
    numBlocks = strtouint32(argv[1], 0L);

  //  The file must exist, and it must be large enough to contain all
  //  that we want to write.  So, we create the file and fill it with
  //  junk.
  //
  ww  = (uint32 *)malloc(sizeof(uint32) * blockSize);
  if (ww == NULL) {
    fprintf(stderr, "can't allocate %d uint32's for clearing the file.\n", blockSize);
    exit(1);
  }
  errno = 0;
  out = fopen("mmap.test.junk", "w");
  if (errno) {
    fprintf(stderr, "can't open 'mmap.test.junk' to fill with junk: %s\n", strerror(errno));
    exit(1);
  }
  for (idx=0; idx<numBlocks; idx++) {
    fprintf(stderr, "Writing initial blocks: "uint32FMT"/"uint32FMT"\r", idx, numBlocks), fflush(stderr);
    fwrite(ww, sizeof(uint32), 1048576, out);
    if (errno) {
      fprintf(stderr, "can't write to 'mmap.test.junk': %s\n", strerror(errno));
      exit(1);
    }
  }
  fclose(out);
  free(ww);
  fprintf(stderr, "\n");

  //  Now, map it, and fill it with real data.
  //
  ww = (uint32 *)mapFile("mmap.test.junk", &lw, 'w');
  for (idx=0; idx<numBlocks * blockSize; idx++) {
    if ((idx & 0xfff) == 0)
      fprintf(stderr, "Writing: "uint32FMT"/"uint32FMT"\r", idx, numBlocks * blockSize), fflush(stderr);
    ww[idx] = idx;
  }
  unmapFile(ww, lw);
  fprintf(stderr, "\n");

  //  Map again, and check the data.
  //
  ww = mapFile("mmap.test.junk", &lw, 'r');
  for (idx=0; idx<numBlocks * blockSize; idx++) {
    if ((idx & 0xfff) == 0)
      fprintf(stderr, "Verifying: "uint32FMT"/"uint32FMT"\r", idx, numBlocks * blockSize), fflush(stderr);
    if (ww[idx] != idx)
      err++;
  }
  unmapFile(ww, lw);
  fprintf(stderr, "\n");

  unlink("mmap.test.junk");

  return (err != 0);
}

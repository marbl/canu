
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
 *    src/utility/splitToWordsTest.C
 *
 *  Modifications by:
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "strings.H"

int
main(int argc, char **argv) {
  splitToWords  W;
  splitType     type = splitWords;

  for (uint32 arg=1; arg<argc; arg++) {
    if (strcmp(argv[arg], "-p") == 0) {
      type = splitPaths;
      continue;
    }

    if (strcmp(argv[arg], "-w") == 0) {
      type = splitWords;
      continue;
    }

    W.split(argv[arg], type);

    fprintf(stderr, "'%s'\n", argv[arg]);

    for (uint32 ii=0; ii<W.numWords(); ii++)
      fprintf(stderr, "%02u - '%s'\n", ii, W[ii]);
  }

  exit(0);
}

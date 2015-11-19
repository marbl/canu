
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
 *    kmer/meryl/meryl.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-JAN-02 to 2004-APR-07
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-MAR-25 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAY-23 to 2009-AUG-07
 *      are Copyright 2005,2007-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "meryl.H"

int
main(int argc, char **argv) {
  merylArgs   *args = new merylArgs(argc, argv);

  switch (args->personality) {
    case 'P':
      estimate(args);
      break;

    case 'B':
      build(args);
      break;

    case 'd':
      dumpDistanceBetweenMers(args);
      break;
    case 't':
      dumpThreshold(args);
      break;
    case 'p':
      dumpPositions(args);
      break;
    case 'c':
      countUnique(args);
      break;
    case 'h':
      plotHistogram(args);
      break;

    case PERSONALITY_MIN:
    case PERSONALITY_MINEXIST:
    case PERSONALITY_MAX:
    case PERSONALITY_MAXEXIST:
    case PERSONALITY_ADD:
    case PERSONALITY_AND:
    case PERSONALITY_NAND:
    case PERSONALITY_OR:
    case PERSONALITY_XOR:
      multipleOperations(args);
      break;

    case PERSONALITY_SUB:
    case PERSONALITY_ABS:
    case PERSONALITY_DIVIDE:
      binaryOperations(args);
      break;

    case PERSONALITY_LEQ:
    case PERSONALITY_GEQ:
    case PERSONALITY_EQ:
      unaryOperations(args);
      break;

    default:
      args->usage();
      fprintf(stderr, "%s: unknown personality.  Specify -P, -B, -S or -M!\n", args->execName);
      exit(1);
      break;
  }

  delete args;

  return(0);
}

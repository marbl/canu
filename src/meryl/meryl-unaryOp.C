
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
 *    kmer/meryl/unaryOp.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-JAN-02 to 2004-APR-07
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2007-OCT-12 to 2008-JUN-09
 *      are Copyright 2007-2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "meryl.H"
#include "libmeryl.H"


void
unaryOperations(merylArgs *args) {

  if (args->mergeFilesLen != 1) {
    fprintf(stderr, "ERROR - must have exactly one file!\n");
    exit(1);
  }
  if (args->outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }
  if ((args->personality != PERSONALITY_LEQ) &&
      (args->personality != PERSONALITY_GEQ) &&
      (args->personality != PERSONALITY_EQ)) {
    fprintf(stderr, "ERROR - only personalities lessthan, lessthanorequal,\n");
    fprintf(stderr, "ERROR - greaterthan, greaterthanorequal, and equal\n");
    fprintf(stderr, "ERROR - are supported in unaryOperations().\n");
    fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
    exit(1);
  }

  //  Open the input and output files -- we don't know the number
  //  unique, distinct, and total until after the operation, so we
  //  leave them zero.
  //
  merylStreamReader   *R = new merylStreamReader(args->mergeFiles[0]);
  merylStreamWriter   *W = new merylStreamWriter(args->outputFile, R->merSize(), R->merCompression(), R->prefixSize(), R->hasPositions());

  switch (args->personality) {
    case PERSONALITY_LEQ:
      while (R->nextMer())
        if (R->theCount() <= args->desiredCount)
          W->addMer(R->theFMer(), R->theCount(), R->thePositions());
      break;

    case PERSONALITY_GEQ:
      while (R->nextMer())
        if (R->theCount() >= args->desiredCount)
          W->addMer(R->theFMer(), R->theCount(), R->thePositions());
      break;

    case PERSONALITY_EQ:
      while (R->nextMer())
        if (R->theCount() == args->desiredCount)
          W->addMer(R->theFMer(), R->theCount(), R->thePositions());
      break;
  }

  delete R;
  delete W;
}

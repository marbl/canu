
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
 *    kmer/meryl/binaryOp.C
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
 *    Brian P. Walenz from 2005-MAY-23 to 2014-APR-11
 *      are Copyright 2005,2007-2009,2014 J. Craig Venter Institute, and
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
#include <math.h>

#include "meryl.H"
#include "libmeryl.H"


void
binaryOperations(merylArgs *args) {

  if (args->mergeFilesLen != 2) {
    fprintf(stderr, "ERROR - must have exactly two files!\n");
    exit(1);
  }
  if (args->outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }
  if ((args->personality != PERSONALITY_SUB) &&
      (args->personality != PERSONALITY_ABS) &&
      (args->personality != PERSONALITY_DIVIDE)) {
    fprintf(stderr, "ERROR - only personalities sub and abs\n");
    fprintf(stderr, "ERROR - are supported in binaryOperations().\n");
    fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
    exit(1);
  }

  //  Open the input files, read in the first mer
  //
  merylStreamReader *A = new merylStreamReader(args->mergeFiles[0]);
  merylStreamReader *B = new merylStreamReader(args->mergeFiles[1]);

  A->nextMer();
  B->nextMer();

  //  Make sure that the mersizes agree, and pick a prefix size for
  //  the output
  //
  if (A->merSize() != B->merSize()) {
    fprintf(stderr, "ERROR - mersizes are different!\n");
    fprintf(stderr, "ERROR - mersize of '%s' is "F_U32"\n", args->mergeFiles[0], A->merSize());
    fprintf(stderr, "ERROR - mersize of '%s' is "F_U32"\n", args->mergeFiles[1], B->merSize());
    exit(1);
  }

  //  Open the output file, using the larger of the two prefix sizes
  //
  merylStreamWriter *W = new merylStreamWriter(args->outputFile,
                                               A->merSize(),
                                               A->merCompression(),
                                               (A->prefixSize() > B->prefixSize()) ? A->prefixSize() : B->prefixSize(),
                                               A->hasPositions());




  //  SUB - report A - B
  //  ABS - report the absolute difference between the two files
  //
  //  These two operations are very similar (SUB was derived from ABS), so
  //  any bug found in one is probably in the other.
  //
  kMer    Amer;
  uint32  Acnt = uint32ZERO;
  kMer    Bmer;
  uint32  Bcnt = uint32ZERO;

  switch (args->personality) {
    case PERSONALITY_SUB:
      while (A->validMer() || B->validMer()) {
        Amer = A->theFMer();
        Acnt = A->theCount();
        Bmer = B->theFMer();
        Bcnt = B->theCount();

        //  If the A stream is all out of mers, set Amer to be the
        //  same as Bmer, and set Acnt to zero.  Similar for B.
        //
        if (!A->validMer()) {
          Amer = Bmer;
          Acnt = uint32ZERO;
        }
        if (!B->validMer()) {
          Bmer = Amer;
          Bcnt = uint32ZERO;
        }

        //fprintf(stderr, "sub A="uint64HEX" B="uint64HEX"\n", Amer, Bmer);

        if (Amer == Bmer) {
          W->addMer(Amer, (Acnt > Bcnt) ? Acnt - Bcnt : 0);
          A->nextMer();
          B->nextMer();
        } else if (Amer < Bmer) {
          W->addMer(Amer, Acnt);
          A->nextMer();
        } else {
          B->nextMer();
        }
      }
      break;
    case PERSONALITY_ABS:
      while (A->validMer() || B->validMer()) {
        Amer = A->theFMer();
        Acnt = A->theCount();
        Bmer = B->theFMer();
        Bcnt = B->theCount();

        //  If the A stream is all out of mers, set Amer to be the
        //  same as Bmer, and set Acnt to zero.  Similar for B.
        //
        if (!A->validMer()) {
          Amer = Bmer;
          Acnt = uint32ZERO;
        }
        if (!B->validMer()) {
          Bmer = Amer;
          Bcnt = uint32ZERO;
        }

        if (Amer == Bmer) {
          W->addMer(Amer, (Acnt > Bcnt) ? Acnt - Bcnt : Bcnt - Acnt);
          A->nextMer();
          B->nextMer();
        } else if (Amer < Bmer) {
          W->addMer(Amer, Acnt);
          A->nextMer();
        } else {
          W->addMer(Bmer, Bcnt);
          B->nextMer();
        }
      }
      break;
    case PERSONALITY_DIVIDE:
      while (A->validMer() || B->validMer()) {
        Amer = A->theFMer();
        Acnt = A->theCount();
        Bmer = B->theFMer();
        Bcnt = B->theCount();

        //  If the A stream is all out of mers, set Amer to be the
        //  same as Bmer, and set Acnt to zero.  Similar for B.
        //
        if (!A->validMer()) {
          Amer = Bmer;
          Acnt = uint32ZERO;
        }
        if (!B->validMer()) {
          Bmer = Amer;
          Bcnt = uint32ZERO;
        }

        if (Amer == Bmer) {
          if ((Acnt > 0) && (Bcnt > 0)) {
            double d = 1000.0 * (double)Acnt / (double)Bcnt;
            if (d > 4096.0 * 1024.0 * 1024.0)
              d = 4096.0 * 1024.0 * 1024.0;
            W->addMer(Amer, (uint32)floor(d));
          }
          A->nextMer();
          B->nextMer();
        } else if (Amer < Bmer) {
          A->nextMer();
        } else {
          B->nextMer();
        }
      }
      break;
  }

  delete A;
  delete B;
  delete W;
}


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
 *    kmer/meryl/merge.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-JAN-02 to 2004-APR-08
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-APR-09 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAY-23 to 2014-APR-11
 *      are Copyright 2005-2009,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "meryl.H"
#include "libmeryl.H"



void
multipleOperations(merylArgs *args) {

  if (args->mergeFilesLen < 2) {
    fprintf(stderr, "ERROR - must have at least two databases (you gave "F_U32")!\n", args->mergeFilesLen);
    exit(1);
  }
  if (args->outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }
  if ((args->personality != PERSONALITY_MERGE) &&
      (args->personality != PERSONALITY_MIN) &&
      (args->personality != PERSONALITY_MINEXIST) &&
      (args->personality != PERSONALITY_MAX) &&
      (args->personality != PERSONALITY_MAXEXIST) &&
      (args->personality != PERSONALITY_ADD) &&
      (args->personality != PERSONALITY_AND) &&
      (args->personality != PERSONALITY_NAND) &&
      (args->personality != PERSONALITY_OR) &&
      (args->personality != PERSONALITY_XOR)) {
    fprintf(stderr, "ERROR - only personalities min, minexist, max, maxexist, add, and, nand, or, xor\n");
    fprintf(stderr, "ERROR - are supported in multipleOperations().  (%d)\n", args->personality);
    fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
    exit(1);
  }

  merylStreamReader  **R = new merylStreamReader* [args->mergeFilesLen];
  merylStreamWriter   *W = 0L;

  //  Open the input files, read in the first mer
  //
  for (uint32 i=0; i<args->mergeFilesLen; i++) {
    R[i] = new merylStreamReader(args->mergeFiles[i]);
    R[i]->nextMer();
  }

  //  Verify that the mersizes are all the same
  //
  bool    fail       = false;
  uint32  merSize    = R[0]->merSize();
  uint32  merComp    = R[0]->merCompression();

  for (uint32 i=0; i<args->mergeFilesLen; i++) {
    fail |= (merSize != R[i]->merSize());
    fail |= (merComp != R[i]->merCompression());
  }

  if (fail)
    fprintf(stderr, "ERROR:  mer sizes (or compression level) differ.\n"), exit(1);

  //  Open the output file, using the largest prefix size found in the
  //  input/mask files.
  //
  uint32  prefixSize = 0;
  for (uint32 i=0; i<args->mergeFilesLen; i++)
    if (prefixSize < R[i]->prefixSize())
      prefixSize = R[i]->prefixSize();

  W = new merylStreamWriter(args->outputFile, merSize, merComp, prefixSize, args->positionsEnabled);

  //  We will find the smallest mer in any file, and count the number of times
  //  it is present in the input files.

  bool     moreInput        = true;

  kMer     currentMer;                      //  The current mer we're operating on
  uint32   currentCount     =  uint32ZERO;  //  The count (operation dependent) of this mer
  uint32   currentTimes     =  uint32ZERO;  //  Number of files it's in

  uint32   currentPositionsMax =  0;
  uint32  *currentPositions    =  0L;

  kMer     thisMer;                         //  The mer we just read
  uint32   thisFile         = ~uint32ZERO;  //  The file we read it from
  uint32   thisCount        =  uint32ZERO;  //  The count of the mer we just read

  speedCounter *C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);

  currentMer.setMerSize(merSize);
  thisMer.setMerSize(merSize);

  while (moreInput) {

    //  Find the smallest mer present in any input file.
    //
    moreInput     = false;
    thisMer.clear();
    thisFile      = ~uint32ZERO;
    thisCount     =  uint32ZERO;

    //  Load thisMer with the first valid mer
    for (uint32 i=0; i<args->mergeFilesLen && !moreInput; i++)
      if (R[i]->validMer()) {
        moreInput = true;
        thisCount = R[i]->theCount();
        thisFile  = i;
        thisMer   = R[i]->theFMer();
      }

    //  Now find the smallest one
    if (moreInput) {
      for (uint32 i=thisFile+1; i<args->mergeFilesLen; i++)
        if ((R[i]->validMer()) && (R[i]->theFMer()) < thisMer) {
          moreInput = true;
          thisCount = R[i]->theCount();
          thisFile  = i;
          thisMer   = R[i]->theFMer();
        }
    }

    //  If we've hit a different mer, write out the last one
    if ((moreInput == false) || (thisMer != currentMer)) {
      switch (args->personality) {
        case PERSONALITY_MIN:
        case PERSONALITY_MAX:
          if (currentTimes == args->mergeFilesLen)
            W->addMer(currentMer, currentCount);
          break;
        case PERSONALITY_MERGE:
        case PERSONALITY_MINEXIST:
        case PERSONALITY_MAXEXIST:
        case PERSONALITY_ADD:
          W->addMer(currentMer, currentCount, currentPositions);
          break;
        case PERSONALITY_AND:
          if (currentTimes == args->mergeFilesLen)
            W->addMer(currentMer, currentCount);
          break;
        case PERSONALITY_NAND:
          if (currentTimes != args->mergeFilesLen)
            W->addMer(currentMer, currentCount);
          break;
        case PERSONALITY_OR:
          W->addMer(currentMer, currentCount);
          break;
        case PERSONALITY_XOR:
          if ((currentTimes % 2) == 1)
            W->addMer(currentMer, currentCount);
          break;
        default:
          fprintf(stderr, "ERROR - invalid personality in multipleOperations::write\n");
          fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
          exit(1);
          break;
      }

      currentMer = thisMer;

      currentCount = uint32ZERO;
      currentTimes = uint32ZERO;

      C->tick();
    }

    //  All done?  Exit.
    if (moreInput == false)
      continue;

    //  Perform the operation
    switch (args->personality) {
      case PERSONALITY_MERGE:
        if (R[thisFile]->thePositions()) {
          if (currentPositionsMax == 0) {
            currentPositionsMax = 1048576;
            currentPositions    = new uint32 [currentPositionsMax];
          }

          if (currentPositionsMax < currentCount + thisCount) {
            while (currentPositionsMax < currentCount + thisCount)
              currentPositionsMax *= 2;

            uint32 *t = new uint32 [currentPositionsMax];
            memcpy(t, currentPositions, sizeof(uint32) * currentCount);
            delete [] currentPositions;
            currentPositions = t;
          }

          if (thisCount < 16) {
            uint32 *p = R[thisFile]->thePositions();
            for (uint32 i=0; i<thisCount; i++)
              currentPositions[currentCount + i] = p[i];
          } else {
            memcpy(currentPositions + currentCount, R[thisFile]->thePositions(), sizeof(uint32) * thisCount);
          }
        }
        //  Otherwise, we're the same as ADD.
        currentCount += thisCount;
        break;
      case PERSONALITY_MIN:
      case PERSONALITY_MINEXIST:
        if (currentTimes == 0) {
          currentCount = thisCount;
        } else {
          if (currentCount > thisCount)
            currentCount = thisCount;
        }
        break;
      case PERSONALITY_MAX:
      case PERSONALITY_MAXEXIST:
        if (currentCount < thisCount)
          currentCount = thisCount;
        break;
      case PERSONALITY_ADD:
        currentCount += thisCount;
        break;
      case PERSONALITY_AND:
      case PERSONALITY_NAND:
      case PERSONALITY_OR:
      case PERSONALITY_XOR:
        currentCount = 1;
        break;
      default:
        fprintf(stderr, "ERROR - invalid personality in multipleOperations::operate\n");
        fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
        exit(1);
        break;
    }

    currentTimes++;

    //  Move the file we just read from to the next mer
    R[thisFile]->nextMer();
  }

  for (uint32 i=0; i<args->mergeFilesLen; i++)
    delete R[i];
  delete R;
  delete W;
  delete C;
}

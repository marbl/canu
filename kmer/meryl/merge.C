#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "meryl.H"
#include "libmeryl.H"



void
multipleOperations(merylArgs *args) {

  if (args->mergeFilesLen < 2) {
    fprintf(stderr, "ERROR - must have at least two databases!\n");
    exit(1);
  }
  if (args->outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }
  if ((args->personality != PERSONALITY_MIN) &&
      (args->personality != PERSONALITY_MINEXIST) &&
      (args->personality != PERSONALITY_MAX) &&
      (args->personality != PERSONALITY_ADD) &&
      (args->personality != PERSONALITY_AND) &&
      (args->personality != PERSONALITY_NAND) &&
      (args->personality != PERSONALITY_OR) &&
      (args->personality != PERSONALITY_XOR)) {
    fprintf(stderr, "ERROR - only personalities min, minexist, max, add, and, nand, or, xor\n");
    fprintf(stderr, "ERROR - are supported in multipleOperations().\n");
    fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
    exit(1);
  }



  merylStreamReader  **R = new merylStreamReader* [args->mergeFilesLen];
  merylStreamWriter   *W = 0L;

  //  Open the input files, read in the first mer
  //
  for (u32bit i=0; i<args->mergeFilesLen; i++) {
    R[i] = new merylStreamReader(args->mergeFiles[i]);
    R[i]->nextMer();
  }

  //  Verify that the mersizes are all the same
  //
  bool    fail       = false;
  u32bit  merSize    = R[0]->merSize();

  for (u32bit i=0; i<args->mergeFilesLen; i++)
    fail |= (merSize != R[i]->merSize());

  if (fail) {
    fprintf(stderr, "ERROR:  mer sizes are different.\n");
    exit(1);
  }


  //  Open the output file, using the largest prefix size found in the
  //  input/mask files.
  //
  u32bit  prefixSize = 0;
  for (u32bit i=0; i<args->mergeFilesLen; i++)
    if (prefixSize < R[i]->prefixSize())
      prefixSize = R[i]->prefixSize();

  W = new merylStreamWriter(args->outputFile,
                            merSize,
                            prefixSize);



  //  We will find the smallest mer in any file, and count the number of times
  //  it is present in the input files.

  bool     moreInput        = true;

  u64bit   currentMer       =  u64bitZERO;  //  The current mer we're operating on
  u32bit   currentCount     =  u32bitZERO;  //  The count (operation dependent) of this mer
  u32bit   currentTimes     =  u32bitZERO;  //  Number of files it's in

  u64bit   thisMer          = ~u64bitZERO;  //  The mer we just read
  u32bit   thisFile         = ~u32bitZERO;  //  The file we read it from
  u32bit   thisCount        =  u32bitZERO;  //  The count of the mer we just read

  speedCounter *C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);

  while (moreInput) {

    //  Find the smallest mer present in any input file.
    //
    moreInput     = false;
    thisMer       = ~u64bitZERO;
    thisFile      = ~u32bitZERO;
    thisCount     =  u32bitZERO;

    for (u32bit i=0; i<args->mergeFilesLen; i++)
      if ((R[i]->validMer()) && (R[i]->theFMer()) < thisMer) {
        moreInput     = true;
        thisMer       = R[i]->theFMer();
        thisCount     = R[i]->theCount();
        thisFile      = i;
      }

    //  If we've hit a different mer, write out the last one
    //
    if ((moreInput == false) || (thisMer != currentMer)) {
      switch (args->personality) {
        case PERSONALITY_MIN:
          if (currentTimes == args->mergeFilesLen)
            W->addMer(currentMer, currentCount);
          break;
        case PERSONALITY_MINEXIST:
        case PERSONALITY_MAX:
        case PERSONALITY_ADD:
          W->addMer(currentMer, currentCount);
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

      currentMer   = thisMer;
      currentCount = u32bitZERO;
      currentTimes = u32bitZERO;

      C->tick();
    }


    //  All done?  Exit.
    if (moreInput == false)
      continue;


    //  Move the file we just read from to the next mer
    //
    R[thisFile]->nextMer();


    //  Perform the operation
    //
    switch (args->personality) {
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
  }

  if (args->beVerbose)
    fprintf(stderr, "\n");

  for (u32bit i=0; i<args->mergeFilesLen; i++)
    delete R[i];
  delete R;
  delete W;
  delete C;
}

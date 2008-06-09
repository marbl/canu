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

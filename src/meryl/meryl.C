#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "bio++.H"
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

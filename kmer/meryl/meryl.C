#include <stdio.h>
#include <stdlib.h>
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
    case 'S':
      scan(args);
      break;

    case 't':
      dumpThreshold(args);
      break;
    case 'c':
      countUnique(args);
      break;
    case 'p':
      plotDistanceBetweenMers(args);
      break;
    case 'h':
      plotHistogram(args);
      break;

    case PERSONALITY_MIN:
    case PERSONALITY_MINEXIST:
    case PERSONALITY_MAX:
    case PERSONALITY_ADD:
    case PERSONALITY_AND:
    case PERSONALITY_NAND:
    case PERSONALITY_OR:
    case PERSONALITY_XOR:
      multipleOperations(args);
      break;

    case PERSONALITY_SUB:
    case PERSONALITY_ABS:
      binaryOperations(args);
      break;

    case PERSONALITY_LEQ:
    case PERSONALITY_GEQ:
    case PERSONALITY_EQ:
      unaryOperations(args);
      break;

    case PERSONALITY_NOR:
    case PERSONALITY_NOT:
      //negate(args);
      break;

    default:
      args->usage();
      fprintf(stderr, "\nERROR:  Unknown personality.  Specify -P, -B, -S or -M!\n");
      exit(1);
      break;
  }

  if (args->statsFile) {
    errno = 0;
    FILE *F = fopen(args->statsFile, "w");
    if (errno) {
      fprintf(stderr, "WARNING: Failed to open stats file '%s'\n%s\n", args->statsFile, strerror(errno));
    } else {
      write_rusage(F);
      fclose(F);
    }
  }
}

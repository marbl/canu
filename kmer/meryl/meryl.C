#include <stdio.h>
#include <stdlib.h>
#include "meryl.H"
#include "libbri.H"
#include "buildinfo-meryl.h"
#include "buildinfo-libbri.h"

extern const char *usage;

int
main(int argc, char **argv) {
  bool              doForward        = true;
  bool              doReverse        = false;
  bool              doCanonical      = false;
  u32bit            merSize          = 20;
  u32bit            tblSize          = 24;
  u32bit            hashSize         = 0;
  char             *inputFile        = 0L;
  char             *outputFile       = 0L;
  char             *queryFile        = 0L;
  char             *maskFile         = 0L;
  bool              beVerbose        = false;

  u32bit            lowCount         = 0;
  u32bit            highCount        = ~lowCount;
  u32bit            desiredCount     = 0;

  bool              outputCount      = false;
  bool              outputAll        = false;
  bool              outputPosition   = false;

  bool              includeDefLine   = false;
  bool              includeMer       = false;

  u32bit            mergeFilesMax    = 0;
  u32bit            mergeFilesLen    = 0;
  char            **mergeFiles       = 0L;

  u64bit            estimatedNumMers = 0;

  char              personality      = 0;

  char             *statsFile        = 0L;

  if (argc == 1) {
    fprintf(stderr, usage, argv[0], argv[0], argv[0]);
    exit(1);
  }

  //  Count how many '-s' switches there are, then allocate space
  //  for them in mergeFiles.
  //
  for (int arg=1; arg < argc; arg++) {
    if (strcmp(argv[arg], "-s") == 0)
      mergeFilesMax++;
  }
  mergeFiles = new char* [mergeFilesMax];


  //  Parse the options
  //
  for (int arg=1; arg < argc; arg++) {
    if        (strncmp(argv[arg], "--buildinfo", 3) == 0) {
      buildinfo_meryl(stderr);
      buildinfo_libbri(stderr);
      exit(1);
    } else if (strcmp(argv[arg], "-m") == 0) {
      arg++;
      merSize = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-mask") == 0) {
      arg++;
      maskFile = argv[arg];
    } else if (strcmp(argv[arg], "-s") == 0) {
      arg++;
      mergeFiles[mergeFilesLen++] = inputFile = argv[arg];
    } else if (strcmp(argv[arg], "-stats") == 0) {
      arg++;
      statsFile = argv[arg];
    } else if (strcmp(argv[arg], "-n") == 0) {
      arg++;
      estimatedNumMers = atol(argv[arg]);
    } else if (strcmp(argv[arg], "-f") == 0) {
      doForward   = true;
      doReverse   = false;
      doCanonical = false;
    } else if (strcmp(argv[arg], "-4") == 0) {
      doForward   = false;
      doReverse   = true;
      doCanonical = false;
    } else if (strcmp(argv[arg], "-C") == 0) {
      doForward   = false;
      doReverse   = false;
      doCanonical = true;
    } else if (strcmp(argv[arg], "-L") == 0) {
      arg++;
      lowCount = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-U") == 0) {
      arg++;
      highCount = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-t") == 0) {
      arg++;
      tblSize = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-H") == 0) {
      arg++;
      hashSize = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-o") == 0) {
      arg++;
      outputFile = argv[arg];
    } else if (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;
    } else if (strcmp(argv[arg], "-q") == 0) {
      arg++;
      queryFile = argv[arg];
    } else if (strcmp(argv[arg], "-d") == 0) {
      includeDefLine = true;
    } else if (strcmp(argv[arg], "-e") == 0) {
      includeMer = true;
    } else if (strcmp(argv[arg], "-c") == 0) {
      outputCount    = true;
      outputAll      = false;
      outputPosition = false;
    } else if (strcmp(argv[arg], "-a") == 0) {
      outputCount    = false;
      outputAll      = true;
      outputPosition = false;
    } else if (strcmp(argv[arg], "-p") == 0) {
      outputCount    = false;
      outputAll      = false;
      outputPosition = true;
    } else if (strcmp(argv[arg], "-P") == 0) {
      personality = 'P';
    } else if (strcmp(argv[arg], "-B") == 0) {
      personality = 'B';
    } else if (strcmp(argv[arg], "-S") == 0) {
      personality = 'S';
    } else if (strcmp(argv[arg], "-M") == 0) {
      arg++;
      if        (strcmp(argv[arg], "min") == 0) {
        personality = PERSONALITY_MIN;
      } else if (strcmp(argv[arg], "minexist") == 0) {
        personality = PERSONALITY_MINEXIST;
      } else if (strcmp(argv[arg], "max") == 0) {
        personality = PERSONALITY_MAX;
      } else if (strcmp(argv[arg], "add") == 0) {
        personality = PERSONALITY_ADD;
      } else if (strcmp(argv[arg], "sub") == 0) {
        personality = PERSONALITY_SUB;
      } else if (strcmp(argv[arg], "abs") == 0) {
        personality = PERSONALITY_ABS;
      } else if (strcmp(argv[arg], "and") == 0) {
        personality = PERSONALITY_AND;
      } else if (strcmp(argv[arg], "nand") == 0) {
        personality = PERSONALITY_NAND;
      } else if (strcmp(argv[arg], "or") == 0) {
        personality = PERSONALITY_OR;
      } else if (strcmp(argv[arg], "nor") == 0) {
        personality = PERSONALITY_NOR;
      } else if (strcmp(argv[arg], "not") == 0) {
        personality = PERSONALITY_NOT;
      } else if (strcmp(argv[arg], "xor") == 0) {
        personality = PERSONALITY_XOR;
      } else if (strcmp(argv[arg], "lessthan") == 0) {
        personality = PERSONALITY_LEQ;
        arg++;
        desiredCount = atoi(argv[arg]) - 1;
      } else if (strcmp(argv[arg], "lessthanorequal") == 0) {
        personality = PERSONALITY_LEQ;
        arg++;
        desiredCount = atoi(argv[arg]);
      } else if (strcmp(argv[arg], "greaterthan") == 0) {
        personality = PERSONALITY_GEQ;
        arg++;
        desiredCount = atoi(argv[arg]) + 1;
      } else if (strcmp(argv[arg], "greaterthanorequal") == 0) {
        personality = PERSONALITY_GEQ;
        arg++;
        desiredCount = atoi(argv[arg]);
      } else if (strcmp(argv[arg], "equal") == 0) {
        personality = PERSONALITY_EQ;
        arg++;
        desiredCount = atoi(argv[arg]);
      } else {
        fprintf(stderr, "ERROR: unknown math personality %s\n", argv[arg]);
        exit(1);
      }
    } else if (strcmp(argv[arg], "-Dc") == 0) {
      personality = 'c';
    } else if (strcmp(argv[arg], "-Dp") == 0) {
      personality = 'p';
    } else if (strcmp(argv[arg], "-Dt") == 0) {
      personality = 't';
    } else if (strcmp(argv[arg], "-Dh") == 0) {
      personality = 'h';
    } else if (strcmp(argv[arg], "-D2") == 0) {
      personality = '2';
    } else {
      fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
    }
  }

  switch (personality) {
    case 'P':
      estimate(inputFile, merSize, estimatedNumMers, beVerbose);
      break;
    case 'B':
      build(inputFile, outputFile, merSize, tblSize, hashSize, lowCount, highCount, doForward, doReverse, doCanonical, beVerbose);
      break;
    case 'S':
      scan(queryFile, inputFile, outputFile,
           includeDefLine, includeMer, doForward, doReverse, doCanonical,
           outputCount, outputAll, outputPosition, beVerbose);
      break;
    case 't':
      dumpThreshold(inputFile, (u32bit)estimatedNumMers);
      break;
    case 'c':
      countUnique(inputFile);
      break;
    case 'p':
      plotDistanceBetweenMers(inputFile);
      break;
    case 'h':
      plotHistogram(inputFile);
      break;

    case PERSONALITY_MIN:
    case PERSONALITY_MINEXIST:
    case PERSONALITY_MAX:
    case PERSONALITY_ADD:
    case PERSONALITY_AND:
    case PERSONALITY_NAND:
    case PERSONALITY_OR:
    case PERSONALITY_XOR:
      multipleOperations(personality, mergeFiles, mergeFilesLen, maskFile, outputFile, beVerbose);
      break;

    case PERSONALITY_SUB:
    case PERSONALITY_ABS:
      binaryOperations(personality, mergeFiles, mergeFilesLen, maskFile, outputFile, beVerbose);
      break;

    case PERSONALITY_LEQ:
    case PERSONALITY_GEQ:
    case PERSONALITY_EQ:
      unaryOperations(personality, mergeFiles, mergeFilesLen, desiredCount, outputFile, beVerbose);
      break;

    case PERSONALITY_NOR:
    case PERSONALITY_NOT:
      //negate(mergeFiles, mergeFilesLen, outputFile, beVerbose);
      break;

    default:
      fprintf(stderr, usage, argv[0], argv[0]);
      fprintf(stderr, "\nERROR:  Unknown personality.  Specify -P, -B, -S or -M!\n");
      exit(1);
      break;
  }

  if (statsFile) {
    errno = 0;
    FILE *F = fopen(statsFile, "w");
    if (errno) {
      fprintf(stderr, "WARNING: Failed to open stats file '%s'\n%s\n", statsFile, strerror(errno));
    } else {
      write_rusage(F);
      fclose(F);
    }
  }
}

#ifndef MERYL_H
#define MERYL_H

#include "libbri.H"
#include "mcBucket.H"
#include "mcDescription.H"


#define PERSONALITY_MIN           0x01
#define PERSONALITY_MINEXIST      0x02
#define PERSONALITY_MAX           0x03
#define PERSONALITY_ADD           0x04
#define PERSONALITY_SUB           0x05
#define PERSONALITY_ABS           0x06
#define PERSONALITY_AND           0x07
#define PERSONALITY_NAND          0x08
#define PERSONALITY_OR            0x09
#define PERSONALITY_NOR           0x0a
#define PERSONALITY_NOT           0x0b
#define PERSONALITY_XOR           0x0c
#define PERSONALITY_LEQ           0x0d
#define PERSONALITY_GEQ           0x0e
#define PERSONALITY_EQ            0x0f


void
estimate(char   *inputFile,
         u32bit  merSize,
         u64bit  numMers,
         bool    beVerbose);
         
void
build(char   *inputFile,
      char   *outputFile,
      u32bit  merSize,
      u32bit  tblSize,
      u32bit  hashSize,
      u32bit  lowCount,
      u32bit  highCount,
      bool    doForward,
      bool    doReverse,
      bool    doCanonical,
      bool    beVerbose);

void
scan(char   *queryFile,
     char   *inputFile,
     char   *outputFile,
     bool    includeDefLine,
     bool    includeMer,
     bool    doForward,
     bool    doReverse,
     bool    doCanonical,
     bool    outputCount,
     bool    outputAll,
     bool    outputPosition,
     bool    beVerbose);

bool
checkSingleDescription(mcDescription *a, char *A,
                       mcDescription *b, char *B);

void
multipleOperations(char      personality,
                   char    **mergeFiles,
                   u32bit    mergeFilesLen,
		   char     *maskFile,
                   char     *outputFile,
                   bool      beVerbose);

void
binaryOperations(char     personality,
                 char   **mergeFiles,
                 u32bit   mergeFilesLen,
		 char    *maskFile,
                 char    *outputFile,
                 bool     beVerbose);

void
unaryOperations(char     personality,
                char   **mergeFiles,
                u32bit   mergeFilesLen,
                u32bit   desiredCount,
                char    *outputFile,
                bool     beVerbose);



//  Debugging stuff
//
void dump(char *inputFile);
void dumpThreshold(char *inputFile, u32bit threshold);
void countUnique(char *inputFile);
void plotDistanceBetweenMers(char *inputFile);
void plotHistogram(char *inputFile);



#endif  //  MERYL_H

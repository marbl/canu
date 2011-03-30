#ifndef MERYL_H
#define MERYL_H

#include "bio++.H"

#define PERSONALITY_MERGE         0xff

#define PERSONALITY_MIN           0x01
#define PERSONALITY_MINEXIST      0x02
#define PERSONALITY_MAX           0x03
#define PERSONALITY_MAXEXIST      0x04
#define PERSONALITY_ADD           0x05
#define PERSONALITY_SUB           0x06
#define PERSONALITY_DIVIDE        0x07
#define PERSONALITY_ABS           0x08

#define PERSONALITY_AND           0x10
#define PERSONALITY_NAND          0x11
#define PERSONALITY_OR            0x12
#define PERSONALITY_XOR           0x13
#define PERSONALITY_LEQ           0x14
#define PERSONALITY_GEQ           0x15
#define PERSONALITY_EQ            0x16

class merylArgs {
public:
  merylArgs(int argc, char **argv);
  merylArgs(const char *prefix);
  ~merylArgs();

  void              usage(void);
  void              clear(void);

  u64bit            hash(kMer const &mer) {
    return(mer.startOfMer(numBuckets_log2));
  };

  bool              writeConfig(void);
  bool              readConfig(const char *prefix);
public:
  char             *execName;
  char             *options;

  bool              beVerbose;
  bool              doForward;
  bool              doReverse;
  bool              doCanonical;

  char             *inputFile;
  char             *outputFile;
  char             *queryFile;

  u32bit            merSize;
  u32bit            merComp;
  bool              positionsEnabled;

  u64bit            numMersEstimated;
  u64bit            numMersActual;

  u64bit            numBasesActual;

  u64bit            mersPerBatch;
  u64bit            basesPerBatch;

  u64bit            numBuckets;
  u32bit            numBuckets_log2;
  u32bit            merDataWidth;
  u64bit            merDataMask;
  u32bit            bucketPointerWidth;

  u32bit            numThreads;
  u64bit            memoryLimit;
  u64bit            segmentLimit;
  bool              configBatch;
  bool              countBatch;
  bool              mergeBatch;
  u32bit            batchNumber;

  char             *sgeJobName;
  char             *sgeBuildOpt;
  char             *sgeMergeOpt;
  bool              isOnGrid;

  u32bit            lowCount;
  u32bit            highCount;
  u32bit            desiredCount;

  bool              outputCount;
  bool              outputAll;
  bool              outputPosition;

  bool              includeDefLine;
  bool              includeMer;

  u32bit            mergeFilesMax;
  u32bit            mergeFilesLen;
  char            **mergeFiles;

  u32bit            personality;
};


u64bit
estimateNumMersInMemorySize(u32bit merSize,
                            u32bit mem,
                            bool   positionsEnabled,
                            bool   beVerbose);

u32bit
optimalNumberOfBuckets(u32bit merSize,
                       u64bit numMers,
                       bool   positionsEnabled);

void estimate(merylArgs *args);
void build(merylArgs *args);

void multipleOperations(merylArgs *args);
void binaryOperations(merylArgs *args);
void unaryOperations(merylArgs *args);

void dump(merylArgs *args);
void dumpThreshold(merylArgs *args);
void dumpPositions(merylArgs *args);
void countUnique(merylArgs *args);
void dumpDistanceBetweenMers(merylArgs *args);
void plotHistogram(merylArgs *args);

#endif  //  MERYL_H

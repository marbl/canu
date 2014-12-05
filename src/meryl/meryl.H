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

  uint64            hash(kMer const &mer) {
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

  uint32            merSize;
  uint32            merComp;
  bool              positionsEnabled;

  uint64            numMersEstimated;
  uint64            numMersActual;

  uint64            numBasesActual;

  uint64            mersPerBatch;
  uint64            basesPerBatch;

  uint64            numBuckets;
  uint32            numBuckets_log2;
  uint32            merDataWidth;
  uint64            merDataMask;
  uint32            bucketPointerWidth;

  uint32            numThreads;
  uint64            memoryLimit;
  uint64            segmentLimit;
  bool              configBatch;
  bool              countBatch;
  bool              mergeBatch;
  uint32            batchNumber;

  char             *sgeJobName;
  char             *sgeBuildOpt;
  char             *sgeMergeOpt;
  bool              isOnGrid;

  uint32            lowCount;
  uint32            highCount;
  uint32            desiredCount;

  bool              outputCount;
  bool              outputAll;
  bool              outputPosition;

  bool              includeDefLine;
  bool              includeMer;

  uint32            mergeFilesMax;
  uint32            mergeFilesLen;
  char            **mergeFiles;

  uint32            personality;
};


uint64
estimateNumMersInMemorySize(uint32 merSize,
                            uint32 mem,
                            bool   positionsEnabled,
                            bool   beVerbose);

uint32
optimalNumberOfBuckets(uint32 merSize,
                       uint64 numMers,
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

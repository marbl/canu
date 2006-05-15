#ifndef MERYL_H
#define MERYL_H

#include "bio++.H"

#define PERSONALITY_MIN           0x01
#define PERSONALITY_MINEXIST      0x02
#define PERSONALITY_MAX           0x03
#define PERSONALITY_ADD           0x04
#define PERSONALITY_SUB           0x05
#define PERSONALITY_ABS           0x06
#define PERSONALITY_AND           0x07
#define PERSONALITY_NAND          0x08
#define PERSONALITY_OR            0x09
#define PERSONALITY_XOR           0x0c
#define PERSONALITY_LEQ           0x0d
#define PERSONALITY_GEQ           0x0e
#define PERSONALITY_EQ            0x0f

class merylArgs {
public:
  merylArgs(int argc, char **argv);
  merylArgs(const char *prefix);
  ~merylArgs();

  void              usage(void);

  u64bit            hash(kMer const &mer) {
    return(mer.startOfMer(numBuckets_log2));
    //return((mer >> merDataWidth) & bucketPointerMask);
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
  char             *statsFile;

  u32bit            merSize;

  u64bit            numMersEstimated;
  u64bit            numMersActual;

  u64bit            mersPerBatch;
  u64bit            numBuckets;
  u32bit            numBuckets_log2;
  u32bit            merDataWidth;
  u64bit            merDataMask;
  //u64bit            bucketPointerMask;
  u32bit            bucketPointerWidth;

  u32bit            numThreads;
  u64bit            memoryLimit;
  u64bit            segmentLimit;
  bool              configBatch;
  bool              countBatch;
  bool              mergeBatch;
  u32bit            batchNumber;

  char             *sgeJobName;
  char             *sgeOptions;
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

  char              personality;
#ifdef __APPLE__
  char              pad1[3];
#endif

};


u64bit
estimateNumMersInMemorySize(u32bit merSize,
                            u32bit mem,
                            bool   beVerbose);

u32bit
optimalNumberOfBuckets(u32bit merSize,
                       u64bit numMers);

void estimate(merylArgs *args);
void build(merylArgs *args);

void multipleOperations(merylArgs *args);
void binaryOperations(merylArgs *args);
void unaryOperations(merylArgs *args);

//  Debugging stuff
//
void dump(merylArgs *args);
void dumpThreshold(merylArgs *args);
void countUnique(merylArgs *args);
void plotDistanceBetweenMers(merylArgs *args);
void plotHistogram(merylArgs *args);

#endif  //  MERYL_H

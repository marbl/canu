#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libbri.H"
#include "mcBucket.H"
#include "mcDescription.H"
#include "libmeryl.H"

#ifdef TRUE64BIT
const char *str_count      = ">%lu\n%s\n";
const char *str_uniq       = "Found %lu mers.\nFound %lu distinct mers.\nFound %lu unique mers.\n";
const char *str_u32bit     = "%lu\n";
const char *str_hist       = "Found %lu mers.\nFound %lu distinct mers.\nFound %lu unique mers.\nFound %lu huge mers (count >= %lu).\nLargest mercount is %lu.\n";
#else
const char *str_count      = ">%llu\n%s\n";
const char *str_uniq       = "Found %llu mers.\nFound %llu distinct mers.\nFound %llu unique mers.\n";
const char *str_u32bit     = "%llu\n";
const char *str_hist       = "Found %llu mers.\nFound %llu distinct mers.\nFound %llu unique mers.\nFound %llu huge mers (count >= %llu).\nLargest mercount is %llu.\n";
#endif


void
dumpThreshold(char *inputFile, u32bit threshold) {
  merStreamFromMeryl  *M = new merStreamFromMeryl(inputFile);
  char                 ms[33];

  while (M->nextMer()) {
    if (M->theCount() >= threshold) {
      for (u32bit z=0; z<M->merSize(); z++)
        ms[M->merSize()-z-1] = decompressSymbol[(M->theMer() >> (2*z)) & 0x03];
      ms[M->merSize()] = 0;

      fprintf(stdout, str_count, M->theCount(), ms);
    }
  }

  delete M;
}


void
countUnique(char *inputFile) {
  merStreamFromMeryl  *M = new merStreamFromMeryl(inputFile);

  u64bit numDistinct = 0;
  u64bit numUnique   = 0;
  u64bit numMers     = 0;
  u64bit c           = 0;

  while (M->nextMer()) {
    c = M->theCount();

    numDistinct++;
    if (c == 1)
      numUnique++;
    numMers += c;
  }

  delete M;

  fprintf(stderr, str_uniq, numMers, numDistinct, numUnique);
}


void
plotHistogram(char *inputFile) {
  merStreamFromMeryl  *M = new merStreamFromMeryl(inputFile);

  u64bit   numMers     = 0;
  u64bit   numUnique   = 0;
  u64bit   numDistinct = 0;
  u64bit   numHuge     = 0;
  u64bit   maxCount    = 0;
  u64bit   hugeCount   = 64 * 1024 * 1024;
  u32bit  *H           = new u32bit [hugeCount];
  u64bit   c           = 0;

  for (u32bit i=0; i<hugeCount; i++)
    H[i] = 0;

  while (M->nextMer()) {
    c = M->theCount();

    numDistinct++;
    if (c == 1)
      numUnique++;
    numMers += c;

    if (c < hugeCount) {
      H[c]++;
      if (maxCount < c)
        maxCount = c;
    } else {
      numHuge++;
    }
  }

  fprintf(stderr, str_hist, numMers, numDistinct, numUnique, numHuge, hugeCount, maxCount);

  for (u32bit i=0; i<=maxCount; i++)
    fprintf(stdout, str_u32bit, H[i]);

  delete    M;
  delete [] H;
}


void
plotDistanceBetweenMers(char *inputFile) {
  merStreamFromMeryl  *M = new merStreamFromMeryl(inputFile);
  u32bit               hugeD = 16 * 1024 * 1024;
  u32bit              *D = new u32bit [hugeD];

  for (u32bit i=0; i<hugeD; i++)
    D[i] = 0;

  u64bit  lastMer = 0;
  u64bit  thisMer = 0;

  while (M->nextMer()) {
    thisMer = M->theMer();

    if ((thisMer - lastMer) < hugeD)
      D[thisMer - lastMer]++;
    else
      fprintf(stderr, "Too large!  %lu\n", thisMer - lastMer);

    lastMer = thisMer;
  }

  for (u32bit i=0; i< hugeD; i++)
    fprintf(stdout, str_u32bit, D[i]);

  delete    M;
  delete [] D;
}

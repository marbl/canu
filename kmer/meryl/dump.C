#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "meryl.H"
#include "libmeryl.H"


void
dumpThreshold(merylArgs *args) {
  merylStreamReader   *M = new merylStreamReader(args->inputFile);
  char                 str[1025];

  while (M->nextMer()) {
    if (M->theCount() >= args->numMersEstimated)
      fprintf(stdout, ">"u64bitFMT"\n%s\n",
              M->theCount(),
              M->theFMer().merToString(str));
  }

  delete M;
}


void
countUnique(merylArgs *args) {
  u64bit numDistinct     = 0;
  u64bit numUnique       = 0;
  u64bit numMers         = 0;
  u64bit c               = 0;

  merylStreamReader   *M = new merylStreamReader(args->inputFile);

  while (M->nextMer()) {
    c = M->theCount();

    numDistinct++;
    if (c == 1)
      numUnique++;
    numMers += c;
  }

  delete M;

  fprintf(stdout, "Found "u64bitFMT" mers.\n", numMers);
  fprintf(stdout, "Found "u64bitFMT" distinct mers.\n", numDistinct);
  fprintf(stdout, "Found "u64bitFMT" unique mers.\n", numUnique);
}


void
plotHistogram(merylArgs *args) {
  u64bit   numMers       = 0;
  u64bit   numUnique     = 0;
  u64bit   numDistinct   = 0;
  u64bit   numHuge       = 0;
  u64bit   maxCount      = 0;
  u64bit   hugeCount     = 64 * 1024 * 1024;
  u64bit  *hist          = new u64bit [hugeCount];
  u64bit   count         = 0;

  for (u32bit i=0; i<hugeCount; i++)
    hist[i] = 0;

  merylStreamReader   *M = new merylStreamReader(args->inputFile);

  while (M->nextMer()) {
    count = M->theCount();

    numDistinct++;
    if (count == 1)
      numUnique++;
    numMers += count;

    if (count < hugeCount)
      hist[count]++;
    else
      numHuge++;

    if (maxCount < count)
      maxCount = count;
  }

  if (hugeCount < maxCount)
    maxCount = hugeCount;

  fprintf(stderr, "Found "u64bitFMT" mers.\n", numMers);
  fprintf(stderr, "Found "u64bitFMT" distinct mers.\n", numDistinct);
  fprintf(stderr, "Found "u64bitFMT" unique mers.\n", numUnique);
  fprintf(stderr, "Found "u64bitFMT" huge mers (count >= "u64bitFMT").\n", numHuge, hugeCount);
  fprintf(stderr, "Largest mercount is "u64bitFMT".\n", maxCount);

  u64bit  distinct = 0;
  u64bit  total    = 0;

  for (u32bit i=1; i<=maxCount; i++) {
    if (hist[i] > 0) {
      distinct += hist[i];
      total    += hist[i] * i;

      fprintf(stdout, u32bitFMT"\t"u64bitFMT"\t%.4f\t%.4f\n",
              i,
              hist[i],
              distinct / (double)numDistinct,
              total / (double)numMers);
    }
  }

  delete    M;
  delete [] hist;
}

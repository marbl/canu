#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "meryl.H"
#include "libmeryl.H"


void
dumpThreshold(merylArgs *args) {
  merylStreamReader   *M = new merylStreamReader(args->inputFile);
  char                 ms[33];

  while (M->nextMer()) {
    if (M->theCount() >= args->numMersEstimated) {
      for (u32bit z=0; z<M->merSize(); z++)
        ms[M->merSize()-z-1] = decompressSymbol[(M->theFMer() >> (2*z)) & 0x03];
      ms[M->merSize()] = 0;

      fprintf(stdout, ">"u64bitFMT"\n%s\n", M->theCount(), ms);
    }
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

  fprintf(stderr, "Found "u64bitFMT" mers.\n", numMers);
  fprintf(stderr, "Found "u64bitFMT" distinct mers.\n", numDistinct);
  fprintf(stderr, "Found "u64bitFMT" unique mers.\n", numUnique);
}


void
plotHistogram(merylArgs *args) {
  u64bit   numMers       = 0;
  u64bit   numUnique     = 0;
  u64bit   numDistinct   = 0;
  u64bit   numHuge       = 0;
  u64bit   maxCount      = 0;
  u64bit   hugeCount     = 64 * 1024 * 1024;
  u32bit  *H             = new u32bit [hugeCount];
  u64bit   c             = 0;

  for (u32bit i=0; i<hugeCount; i++)
    H[i] = 0;

  merylStreamReader   *M = new merylStreamReader(args->inputFile);

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

  fprintf(stderr, "Found "u64bitFMT" mers.\n", numMers);
  fprintf(stderr, "Found "u64bitFMT" distinct mers.\n", numDistinct);
  fprintf(stderr, "Found "u64bitFMT" unique mers.\n", numUnique);
  fprintf(stderr, "Found "u64bitFMT" huge mers (count >= "u64bitFMT").\n", numHuge, hugeCount);
  fprintf(stderr, "Largest mercount is "u64bitFMT".\n", maxCount);

  for (u32bit i=0; i<=maxCount; i++)
    fprintf(stdout, u32bitFMT"\n", H[i]);

  delete    M;
  delete [] H;
}


void
plotDistanceBetweenMers(merylArgs *args) {
  u32bit               hugeD = 16 * 1024 * 1024;
  u32bit              *D = new u32bit [hugeD];
  u32bit               numHuge = u32bitZERO;
  u32bit               maxDist = u32bitZERO;
  u64bit               lastMer = u64bitZERO;
  u64bit               thisMer = u64bitZERO;
  merylStreamReader   *M = new merylStreamReader(args->inputFile);

  for (u32bit i=0; i<hugeD; i++)
    D[i] = 0;

  while (M->nextMer()) {
    thisMer = M->theFMer();

    if ((thisMer - lastMer) < hugeD) {
      D[thisMer - lastMer]++;
    } else {
      numHuge++;
      if (maxDist < thisMer - lastMer)
        maxDist = thisMer - lastMer;
    }

    lastMer = thisMer;
  }

  fprintf(stdout, u32bitFMT" distances at least "u32bitFMT", largest is "u32bitFMT"\n",
          numHuge,
          hugeD, maxDist);
  for (u32bit i=0; i< hugeD; i++)
    fprintf(stdout, u32bitFMT"\n", D[i]);

  delete    M;
  delete [] D;
}

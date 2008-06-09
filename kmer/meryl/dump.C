#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "meryl.H"
#include "libmeryl.H"

#include <algorithm>

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
dumpPositions(merylArgs *args) {
  merylStreamReader   *M = new merylStreamReader(args->inputFile);
  char                 str[1025];

  if (M->hasPositions() == false) {
    fprintf(stderr, "File '%s' contains no position information.\n", args->inputFile);
  } else {
    while (M->nextMer()) {
      fprintf(stdout, ">"u64bitFMT, M->theCount());
      for (u32bit i=0; i<M->theCount(); i++)
        fprintf(stdout, " "u32bitFMT, M->getPosition(i));
      fprintf(stdout, "\n%s\n", M->theFMer().merToString(str));
    }
  }

  delete M;
}


void
countUnique(merylArgs *args) {
  merylStreamReader   *M = new merylStreamReader(args->inputFile);

#warning make this a test
#if 0
  u64bit numDistinct     = 0;
  u64bit numUnique       = 0;
  u64bit numMers         = 0;
  u64bit c               = 0;

  while (M->nextMer()) {
    c = M->theCount();

    numDistinct++;
    if (c == 1)
      numUnique++;
    numMers += c;
  }

  assert(numMers     == M->numberOfTotalMers());
  assert(numDistinct == M->numberOfDistinctMers());
  assert(numUnique   == M->numberOfUniqueMers());
  fprintf(stderr, "OK\n");
#endif

  fprintf(stdout, "Found "u64bitFMT" mers.\n",          M->numberOfTotalMers());
  fprintf(stdout, "Found "u64bitFMT" distinct mers.\n", M->numberOfDistinctMers());
  fprintf(stdout, "Found "u64bitFMT" unique mers.\n",   M->numberOfUniqueMers());

  delete M;
}


void
plotHistogram(merylArgs *args) {
  u64bit  distinct = 0;
  u64bit  total    = 0;

  merylStreamReader   *M = new merylStreamReader(args->inputFile);

  fprintf(stderr, "Found "u64bitFMT" mers.\n",          M->numberOfTotalMers());
  fprintf(stderr, "Found "u64bitFMT" distinct mers.\n", M->numberOfDistinctMers());
  fprintf(stderr, "Found "u64bitFMT" unique mers.\n",   M->numberOfUniqueMers());

  fprintf(stderr, "Largest mercount is "u64bitFMT"; "u64bitFMT" mers are too big for histogram.\n",
          M->histogramMaximumCount(), M->histogramHuge());

  for (u32bit i=1; i<M->histogramLength(); i++) {
    u64bit hist = M->histogram(i);

    if (hist > 0) {
      distinct += hist;
      total    += hist * i;

      fprintf(stdout, u32bitFMT"\t"u64bitFMT"\t%.4f\t%.4f\n",
              i,
              hist,
              distinct / (double)M->numberOfDistinctMers(),
              total    / (double)M->numberOfTotalMers());
    }
  }

  delete    M;
}



void
dumpDistanceBetweenMers(merylArgs *args) {
  merylStreamReader   *M = new merylStreamReader(args->inputFile);

  //  This is now tough because we don't know where the sequences end,
  //  and our positions encode position in the chain.

  u32bit  histMax  = 64 * 1024 * 1024;
  u64bit *hist     = new u64bit [histMax];
  u64bit  histHuge = 0;

  if (M->hasPositions() == false) {
    fprintf(stderr, "File '%s' contains no position information.\n", args->inputFile);
  } else {
    while (M->nextMer()) {
      std::sort(M->thePositions(), M->thePositions() + M->theCount());

      for (u32bit i=1; i<M->theCount(); i++) {
        u32bit d = M->getPosition(i) - M->getPosition(i-1);
        if (d < histMax)
          hist[d]++;
        else
          histHuge++;
      }
    }

    u32bit maxd = 0;

    for (u32bit d=0; d<histMax; d++)
      if (hist[d])
        maxd = d+1;

    for (u32bit d=0; d<maxd; d++)
      if (hist[d])
        fprintf(stderr, u32bitFMT"\t"u64bitFMT"\n", d, hist[d]);

    if (histHuge)
      fprintf(stderr, "huge\t"u64bitFMT"\n", histHuge);
  }

  delete [] hist;
  delete    M;
}

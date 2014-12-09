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
      fprintf(stdout, ">"F_U64"\n%s\n",
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
      fprintf(stdout, ">"F_U64, M->theCount());
      for (uint32 i=0; i<M->theCount(); i++)
        fprintf(stdout, " "F_U32, M->getPosition(i));
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
  uint64 numDistinct     = 0;
  uint64 numUnique       = 0;
  uint64 numMers         = 0;
  uint64 c               = 0;

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

  fprintf(stdout, "Found "F_U64" mers.\n",          M->numberOfTotalMers());
  fprintf(stdout, "Found "F_U64" distinct mers.\n", M->numberOfDistinctMers());
  fprintf(stdout, "Found "F_U64" unique mers.\n",   M->numberOfUniqueMers());

  delete M;
}


void
plotHistogram(merylArgs *args) {
  uint64  distinct = 0;
  uint64  total    = 0;

  merylStreamReader   *M = new merylStreamReader(args->inputFile);

  fprintf(stderr, "Found "F_U64" mers.\n",          M->numberOfTotalMers());
  fprintf(stderr, "Found "F_U64" distinct mers.\n", M->numberOfDistinctMers());
  fprintf(stderr, "Found "F_U64" unique mers.\n",   M->numberOfUniqueMers());

  fprintf(stderr, "Largest mercount is "F_U64"; "F_U64" mers are too big for histogram.\n",
          M->histogramMaximumCount(), M->histogramHuge());

  for (uint32 i=1; i<M->histogramLength(); i++) {
    uint64 hist = M->histogram(i);

    if (hist > 0) {
      distinct += hist;
      total    += hist * i;

      fprintf(stdout, F_U32"\t"F_U64"\t%.4f\t%.4f\n",
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

  uint32  histMax  = 64 * 1024 * 1024;
  uint64 *hist     = new uint64 [histMax];
  uint64  histHuge = 0;

  if (M->hasPositions() == false) {
    fprintf(stderr, "File '%s' contains no position information.\n", args->inputFile);
  } else {
    while (M->nextMer()) {
      std::sort(M->thePositions(), M->thePositions() + M->theCount());

      for (uint32 i=1; i<M->theCount(); i++) {
        uint32 d = M->getPosition(i) - M->getPosition(i-1);
        if (d < histMax)
          hist[d]++;
        else
          histHuge++;
      }
    }

    uint32 maxd = 0;

    for (uint32 d=0; d<histMax; d++)
      if (hist[d])
        maxd = d+1;

    for (uint32 d=0; d<maxd; d++)
      if (hist[d])
        fprintf(stderr, F_U32"\t"F_U64"\n", d, hist[d]);

    if (histHuge)
      fprintf(stderr, "huge\t"F_U64"\n", histHuge);
  }

  delete [] hist;
  delete    M;
}

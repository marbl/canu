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
dumpPositions(merylArgs *args) {
  merylStreamReader   *M = new merylStreamReader(args->inputFile);
  char                 str[1025];

#ifdef WITH_POSITIONS
  while (M->nextMer()) {
    fprintf(stdout, ">"u64bitFMT, M->theCount());
    for (u32bit i=0; i<M->theCount(); i++)
      fprintf(stdout, " "u32bitFMT, M->getPosition(i));
    fprintf(stdout, "\n%s\n", M->theFMer().merToString(str));
  }
#else
  fprintf(stderr, "This meryl not compiled with position support.\n");
#endif

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

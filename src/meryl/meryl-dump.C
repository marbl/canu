
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    kmer/meryl/dump.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-JAN-02 to 2004-APR-07
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-APR-12 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAY-23 to 2014-APR-11
 *      are Copyright 2005,2007-2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-FEB-25
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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

  fprintf(stderr, "Largest mercount is "F_U64".\n",
          M->histogramMaximumCount());

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

  memset(hist, 0, sizeof(uint64) * histMax);

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

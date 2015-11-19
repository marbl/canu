
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
 *    kmer/meryl/simple.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-AUG-31 to 2014-APR-11
 *      are Copyright 2010-2011,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2015-APR-24
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "AS_UTL_fileIO.H"

#include "kMer.H"
#include "bitPackedFile.H"
#include "libmeryl.H"

#include "seqStream.H"
#include "merStream.H"
#include "speedCounter.H"

#include "speedCounter.H"
#include "timeAndSize.H"

#include <algorithm>

using namespace std;


#if 0
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include <algorithm>

#include "bio++.H"
#include "meryl.H"

#include "libmeryl.H"
#include "seqStream.H"
#include "merStream.H"

#endif


//  A very simple mer counter.  Allocates a gigantic 32-bit array,
//  populates the array with mers, sorts, writes output.

int
main(int argc, char **argv) {
  char    *inName         = 0L;
  char    *otName         = 0L;
  uint32   merSize        = 22;
  uint32   merCompression = 1;

  bool     doForward      = false;
  bool     doReverse      = false;
  bool     doCanonical    = false;

  speedCounter        *C = 0L;
  merStream           *M = 0L;
  merylStreamWriter   *W = 0L;

  uint64         numMers = 0;

  bool           algSorting  = false;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-i") == 0) {
      inName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      otName = argv[++arg];

    } else if (strcmp(argv[arg], "-m") == 0) {
      merSize = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-f") == 0) {
      doForward   = true;

    } else if (strcmp(argv[arg], "-r") == 0) {
      doReverse   = true;

    } else if (strcmp(argv[arg], "-C") == 0) {
      doCanonical = true;

    } else if (strcmp(argv[arg], "-c") == 0) {
      merCompression = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-sort") == 0) {
      algSorting = true;

    } else if (strcmp(argv[arg], "-direct") == 0) {
      algSorting = false;

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (inName == 0L) {
    fprintf(stderr, "no input given with '-i'\n");
    err++;
  }
  if (otName == 0L) {
    fprintf(stderr, "no output given with '-o'\n");
    err++;
  }
  if (err)
    exit(1);


  M = new merStream(new kMerBuilder(merSize, merCompression),
                    new seqStream(inName),
                    true, true);
  numMers = M->approximateNumberOfMers();
  delete M;

  fprintf(stderr, "Guessing "F_U64" mers in input '%s'\n", numMers, inName);

  M = new merStream(new kMerBuilder(merSize, merCompression),
                    new seqStream(inName),
                    true, true);





  if (algSorting) {
    uint64  theMersLen = 0;
    uint64  theMersMax = 2 * numMers;  //  for allowing both -f and -r
    uint32 *theMers    = new uint32 [theMersMax];

    fprintf(stderr, "Allocating "F_U64"MB for mer storage.\n", numMers * sizeof(uint64) >> 20);

    C = new speedCounter(" Filling mer list: %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, 1);

    while (M->nextMer()) {
      if (doForward)
        theMers[theMersLen++] = M->theFMer();

      if (doReverse)
        theMers[theMersLen++] = M->theRMer();

      if (doCanonical)
        theMers[theMersLen++] = (M->theFMer() <= M->theRMer()) ? M->theFMer() : M->theRMer();

      C->tick();
    }

    delete C;
    delete M;

    fprintf(stderr, "Found "F_U64" mers in input '%s'\n", theMersLen, inName);

    if (theMersLen > theMersMax)
      fprintf(stderr, "ERROR:  too many mers in input!\n"), exit(1);

    sort(theMers, theMers + theMersLen);

    C = new speedCounter(" Writing output:           %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, 1);
    W = new merylStreamWriter(otName,
                              merSize, merCompression,
                              16,
                              false);

    kMer mer(merSize);

    for (uint64 i=0; i<theMersLen; i++) {
      mer.setWord(0, theMers[i]);
      W->addMer(mer, 1, 0L);
      C->tick();
    }

    delete C;
    delete W;

    delete [] theMers;
  }


  else {
    uint64  numCounts = ((uint64)1) << (2 * merSize);
    uint32 *theCounts = new uint32 [numCounts];

    fprintf(stderr, "Allocating "F_U64"MB for count storage.\n", numCounts * sizeof(uint32) >> 20);

    memset(theCounts, 0, sizeof(uint32) * numCounts);

    C = new speedCounter(" Filling mer counts: %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, 1);

    while (M->nextMer()) {
      if (doForward)
        theCounts[M->theFMer()]++;

      if (doReverse)
        theCounts[M->theRMer()]++;

      if (doCanonical)
        theCounts[(M->theFMer() <= M->theRMer()) ? M->theFMer() : M->theRMer()]++;

      C->tick();
    }

    delete C;
    delete M;

    C = new speedCounter(" Writing output:           %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, 1);
    W = new merylStreamWriter(otName,
                              merSize, merCompression,
                              16,
                              false);

    kMer mer(merSize);

    for (uint64 i=0; i<numCounts; i++) {
      if (theCounts[i] > 0) {
        mer.setWord(0, i);
        W->addMer(mer, theCounts[i], 0L);
        C->tick();
      }
    }

    delete C;
    delete W;

    delete [] theCounts;
  }


  exit(0);
}

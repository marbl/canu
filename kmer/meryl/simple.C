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

using namespace std;

//  A very simple mer counter.  Allocates a gigantic 32-bit array,
//  populates the array with mers, sorts, writes output.

int
main(int argc, char **argv) {
  char    *inName         = 0L;
  char    *otName         = 0L;
  u32bit   merSize        = 22;
  u32bit   merCompression = 1;

  bool     doForward      = true;
  bool     doReverse      = false;
  bool     doCanonical    = false;

  speedCounter        *C = 0L;
  merStream           *M = 0L;
  merylStreamWriter   *W = 0L;

  u64bit         numMers = 0;

  u64bit        *theMers    = 0L;
  u64bit         theMersMax = 0;
  u64bit         theMersLen = 0;

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
      doReverse   = false;
      doCanonical = false;

    } else if (strcmp(argv[arg], "-r") == 0) {
      doForward   = false;
      doReverse   = true;
      doCanonical = false;

    } else if (strcmp(argv[arg], "-C") == 0) {
      doForward   = false;
      doReverse   = false;
      doCanonical = true;

    } else if (strcmp(argv[arg], "-c") == 0) {
      merCompression = atoi(argv[++arg]);

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


  {
    M = new merStream(new kMerBuilder(merSize, merCompression),
                      new seqStream(inName),
                      true, true);
    numMers = M->approximateNumberOfMers();
    delete M;
  }

  fprintf(stderr, "Guessing "u64bitFMT" mers in input '%s'\n", numMers, inName);
  fprintf(stderr, "Allocating "u64bitFMT"MB for mer storage.\n", numMers * 8 >> 20);
  
  theMers    = new u64bit [numMers];
  theMersLen = 0;
  theMersMax = numMers;

  ////////////////////////////////////////

  C = new speedCounter(" Counting mers in buckets: %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, 1);
  M = new merStream(new kMerBuilder(merSize, merCompression),
                    new seqStream(inName),
                    true, true);
  //M->setRange(args->mersPerBatch * segment, args->mersPerBatch * segment + args->mersPerBatch);

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

  fprintf(stderr, "Found "u64bitFMT" mers in input '%s'\n", theMersLen, inName);

  if (theMersLen > theMersMax)
    fprintf(stderr, "ERROR:  too many mers in input!\n"), exit(1);

  ////////////////////////////////////////

  fprintf(stderr, "sorting\n");

  sort(theMers, theMers + theMersLen);

  ////////////////////////////////////////

  C = new speedCounter(" Writing output:           %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, 1);
  W = new merylStreamWriter(otName,
                            merSize, merCompression,
                            16,
                            false);

  kMer mer(merSize);

  for (u64bit i=0; i<theMersLen; i++) {
    mer.setWord(0, theMers[i]);
    W->addMer(mer, 1, 0L);
    C->tick();
  }

  delete C;
  delete W;

  ////////////////////////////////////////

  delete [] theMers;

  exit(0);
}

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "meryl.H"
#include "libbri.H"
#include "britime.H"
#include "buildinfo-merle.h"

//
//  Very simple mer counter, works only for small mers.
//
//  Only outputs the histogram, not the full meryl output format.
//

int
main(int argc, char **argv) {
  u32bit   merSizeInBases = 0;
  u32bit   merSizeInBits  = 0;
  u32bit   numMers        = 0;
  char    *inputName      = 0L;
  u32bit  *counts         = 0L;
  bool     forward        = false;
  bool     reverse        = false;
  bool     canonical      = true;

  int arg = 1;
  while (arg < argc) {
    if (strncmp(argv[arg], "-mersize", 2) == 0) {
      merSizeInBases = atoi(argv[++arg]);
      merSizeInBits  = merSizeInBases << 1;
      numMers        = u32bitONE << merSizeInBits;
    } else if (strncmp(argv[arg], "-sequence", 2) == 0) {
      inputName = argv[++arg];
    } else if (strncmp(argv[arg], "-f", 2) == 0) {
      forward   = true;
      reverse   = false;
      canonical = false;
    } else if (strncmp(argv[arg], "-r", 2) == 0) {
      forward   = false;
      reverse   = true;
      canonical = false;
    } else if (strncmp(argv[arg], "-c", 2) == 0) {
      forward   = false;
      reverse   = false;
      canonical = true;
    }
    arg++;
  }

  if ((inputName == 0L) || (merSizeInBits == 0)) {
    fprintf(stderr, "usage: %s [-f | -r | -c] -mersize n -sequence inputname\n", argv[0]);
    fprintf(stderr, "           n <= 14\n");
    exit(1);
  }

  if (merSizeInBits > 28) {
    fprintf(stderr, "ERROR: mersize too large.  Use meryl for mersize > 14.\n");
    exit(1);
  }

  fprintf(stderr, "Using mersize of %u bases (%u bits).\n", merSizeInBases, merSizeInBits);
  fprintf(stderr, "This mersize has %u mers.\n", numMers);
  fprintf(stderr, "Attempting to allocate %u bytes for counts.\n", numMers * sizeof(u32bit));
  counts = new u32bit [numMers];

  fprintf(stderr, "Clearing.\n");
  bzero(counts, sizeof(u32bit) * numMers);

  merStream    M(merSizeInBases, inputName);
  speedCounter C("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, true);

  if (forward) {
    while (M.nextMer()) {
      counts[M.theFMer()]++;
      C.tick();
    }
  }

  if (reverse) {
    while (M.nextMer()) {
      counts[M.theRMer()]++;
      C.tick();
    }
  }

  if (canonical) {
    while (M.nextMer()) {
      if (M.theFMer() < M.theRMer())
        counts[M.theFMer()]++;
      else
        counts[M.theRMer()]++;
      C.tick();
    }
  }


  for (u32bit mer=0; mer<numMers; mer++) {
    char  ms[33];

    if (counts[mer] > 0) {
      for (u32bit z=0; z<merSizeInBases; z++)
        ms[merSizeInBases-z-1] = decompressSymbol[(mer >> (2*z)) & 0x03];
      ms[merSizeInBases] = 0;

      fprintf(stdout, ">%u\n%s\n", counts[mer], ms);
    }
  }
}

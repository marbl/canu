#include <stdio.h>
#include <stdlib.h>
#include "bio++.H"

//  Prints the number of mers in a merStreamFile.

int
main(int argc, char **argv) {

#if 0
  if (argc < 2) {
    fprintf(stderr, "usage: %s some.merStreamFile\n", argv[0]);
    exit(1);
  }

  merStreamFileReader *MSFR = new merStreamFileReader(argv[1]);
  fprintf(stdout, u64bitFMT"\n", MSFR->numberOfMers());
  delete MSFR;
#endif

#warning I AM BROKEN

  exit(0);
}

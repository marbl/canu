#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "existDB.H"

//  Driver for the existDB creation.  Reads a sequence.fasta, builds
//  an existDB for the mers in the file, and then writes the internal
//  structures to disk.
//
//  The existDB constructor is smart enough to read either a pre-built
//  image or a regular multi-fasta file.
//

int
main(int argc, char **argv) {
  u32bit    mersize = 20;
  u32bit    tblsize = 19;

  if (argc < 3) {
    fprintf(stderr, "usage: %s [stuff]\n", argv[0]);
    fprintf(stderr, "       To create an image:\n");
    fprintf(stderr, "         [-m mersize] [-t tablesize] [-C] <sequence.fasta> <datafile>\n", argv[0]);
    fprintf(stderr, "           -C   Insert both the forward and reverse mer into the table.\n");
    fprintf(stderr, "           defaults: -m 20 -t 19 -C (suitable for searchGENOME*)\n");
    fprintf(stderr, "       To dump information about an image:\n");
    fprintf(stderr, "         -d datafile\n");
    exit(1);
  }

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-m", 2) == 0) {
      arg++;
      mersize = atoi(argv[arg]);
    } else if (strncmp(argv[arg], "-t", 2) == 0) {
      arg++;
      tblsize = atoi(argv[arg]);
    } else if (strncmp(argv[arg], "-d", 2) == 0) {
      existDB *e = new existDB(argv[argc-1], false);
      e->printState(stdout);
      delete e;
      exit(0);
    }
    arg++;
  }

  existDB  *e = new existDB(argv[argc-2], mersize, tblsize);
  e->saveState(argv[argc-1]);
  delete e;
  exit(0);
}

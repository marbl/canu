#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

//#include "bio++.H"
#include "sim4.H"

//  Writes n polishes from stdin to stdout, default 1.

int
main(int argc, char **argv) {
  u32bit     numToPrint   = 1;
  FILE      *F            = stdin;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-h", 2) == 0) {
      fprintf(stderr, "usage: %s [-h] [-numPolishes] [polishes-file]\n", argv[0]);
      exit(1);
    } else if (strncmp(argv[arg], "-", 1) == 0) {
      numToPrint = atoi(argv[arg] + 1);
    } else {
      F = fopen(argv[arg], "r");
      if (F == 0L) {
        fprintf(stderr, "Failed to open '%s'.\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }
    }
    arg++;
  }

  sim4polish *p = new sim4polish(F);
  while ((numToPrint--) && (p->_numExons > 0)) {
    p->s4p_printPolish(stdout, S4P_PRINTPOLISH_FULL);
    delete p;
    p = new sim4polish(F);
  }

  return(0);
}


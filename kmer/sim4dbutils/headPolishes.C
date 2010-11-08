#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "sim4.H"

//  Writes n polishes from stdin to stdout, default 1.

int
main(int argc, char **argv) {
  u32bit              numToPrint   = 1;
  sim4polishReader   *R = 0L;
  sim4polishWriter   *W = 0L;

  sim4polishStyle     style = sim4polishStyleDefault;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-h", 2) == 0) {
      err++;

    } else if (strncmp(argv[arg], "-n", 2) == 0) {
      numToPrint = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-gff3") == 0) {
      style = sim4polishGFF3;

    } else if (strncmp(argv[arg], "-", 1) == 0) {
      numToPrint = atoi(argv[arg] + 1);

    } else {
      R = new sim4polishReader(argv[arg]);
    }

    arg++;
  }
  if ((err) || ((R == 0L) && (isatty(fileno(stdin))))) {
    fprintf(stderr, "usage: %s [-h] [-# | -n #] [-gff3] [polishes-file]\n", argv[0]);
    exit(1);
  }

  if (R == 0L)
    R = new sim4polishReader("-");

  if (W == 0L)
    W = new sim4polishWriter("-", style);

  if (R->getsim4polishStyle() != style)
    fprintf(stderr, "warning: Input format and output format differ.\n");

  sim4polish *p = 0L;

  while ((numToPrint--) && (R->nextAlignment(p)))
    W->writeAlignment(p);

  delete W;
  delete R;

  return(0);
}


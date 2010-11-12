#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "bio.h"
#include "sim4.H"

int
main(int argc, char ** argv) {
  sim4polishWriter   *GOOD       = 0L;
  sim4polishStyle     in_style, out_style;

  //  We limit scaffolds to be below the number of open files per
  //  process.
  //

  if (argc != 1) {
    fprintf(stderr, "S4DB to GFF3 format converter.\nUsage: %s < input_file > output_file\n", argv[0]);
    exit(1);
  }

  sim4polishReader *R = new sim4polishReader("-");
  sim4polish       *p = 0L;
 
  in_style = R->getsim4polishStyle();

  if (in_style == sim4polishS4DB)
    out_style = sim4polishGFF3;
  else if (in_style == sim4polishGFF3) 
    out_style = sim4polishS4DB;
  else {
    fprintf(stderr, "ERROR: Unrecognized or unsupported polishes format. Aborting.\n"); exit(1);
  }

  if (GOOD == 0L)
    GOOD = new sim4polishWriter("-", out_style);

  while (R->nextAlignment(p)) {

#if 0
    if (noDefLines)
      p->s4p_removeDefLines();
    if (noAlignments)
      p->s4p_removeAlignments();
#endif

    GOOD->writeAlignment(p);
  }

  delete R;

  delete GOOD;

  return(0);
}

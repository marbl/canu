#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "sim4polish.h"

//  Reads a file of matches, strips out the defline and alignments,
//  and writes a new file.

const char *usage =
"usage: %s [opts]\n"
"  -nodeflines    Strip out deflines\n"
"  -noalignments  Strip out alignments\n"
"  -normalized    Strip out the genomic region (make the polish relative\n"
"                 to the start of the sequence)\n"
"\n";

int
main(int argc, char **argv) {
  sim4polish  *p = 0L;
  u32bit       f = S4P_PRINTPOLISH_NOTVALUABLE;
  int          i  = 0;

  int  arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-nodeflines", 4) == 0) {
      f |= S4P_PRINTPOLISH_NODEFS;
    } else if (strncmp(argv[arg], "-noalignments", 4) == 0) {
      f |= S4P_PRINTPOLISH_NOALIGNS;
    } else if (strncmp(argv[arg], "-normalized", 4) == 0) {
      f |= S4P_PRINTPOLISH_NORMALIZED;
    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }

  if ((isatty(fileno(stdin)) || isatty(fileno(stdout)))) {
    fprintf(stderr, usage, argv[0]);

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");

    if (isatty(fileno(stdout)))
      fprintf(stderr, "error: Please redirect the output to a file.\n\n");

    exit(1);
  }

  while (!feof(stdin)) {
    p = s4p_readPolish(stdin);
    s4p_printPolish(stdout, p, f);
    s4p_destroyPolish(p);
  }

  return(0);
}

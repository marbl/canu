#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "sim4reader.h"


char const *usage =
"usage: %s [-uniq | -dupl] < file > file\n"
"\n";


void
pickBest(sim4polish **p, int pNum, int uniq) {
  int i;

  if (pNum == 1) {
    if (uniq)
      printPolish(stdout, p[0]);
  } else {
    if (!uniq)
      for (i=0; i<pNum; i++)
        printPolish(stdout, p[i]);
  }

  for (i=0; i<pNum; i++)
    destroyPolish(p[i]);
}




int
main(int argc, char **argv) {
  int          pNum   = 0;
  int          pAlloc = 8388608;
  sim4polish **p      = 0L;
  sim4polish  *q      = 0L;
  int          estID  = ~0;

  int          uniq   = 1;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-uniq", 2) == 0) {
      uniq = 1;
    } else if (strncmp(argv[arg], "-dupl", 2) == 0) {
      uniq = 0;
    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }

  if (isatty(fileno(stdin))) {
    fprintf(stderr, usage, argv[0]);

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");

    exit(1);
  }

  //  Read polishes, picking the best when we see a change in
  //  the estID.

  p = (sim4polish **)malloc(sizeof(sim4polish *) * pAlloc);

  while ((q = readPolish(stdin)) != 0L) {
    if ((q->estID != estID) && (pNum > 0)) {
      pickBest(p, pNum, uniq);
      pNum  = 0;
    }

    //  Reallocate pointers?
    //
    if (pNum >= pAlloc) {
      p = (sim4polish **)realloc(p, sizeof(sim4polish *) * (pAlloc *= 2));
      if (p == 0L) {
        fprintf(stderr, "Out of memory: Couldn't allocate space for polish pointers.\n");
        exit(1);
      }
    }

    p[pNum++] = q;
    estID     = q->estID;
  }

  if (pNum > 0)
    pickBest(p, pNum, uniq);

  return(0);
}


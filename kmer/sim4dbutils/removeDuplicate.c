#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "sim4.H"

char const *usage =
"usage: %s < file > file\n"
"\n";

//  Input matches should be sorted by cDNA, and ran through pickBest.
//  This code will remove all matches that have the same genomic span,
//  and warn when two matches have nearly the same genomic span.

void
pickBest(sim4polish **p, int pNum) {
  int i, j;

  //fprintf(stderr, "Pick for %d with %d things\n", p[0]->estID, pNum);

  for (i=0; i<pNum; i++) {
    for (j=i+1; j<pNum; j++) {
      if ((p[i]) &&
          (p[j]) &&
          (p[i]->numExons == p[j]->numExons) &&
          (p[i]->genID == p[j]->genID)) {
        int a, b;
        int sd = 666;
        int ed = 666;

        a = p[i]->exons[0].genFrom;
        b = p[j]->exons[0].genFrom;
        if (a < b)
          sd = b - a;
        else
          sd = a - b;

        a = p[i]->exons[p[i]->numExons-1].genTo;
        b = p[j]->exons[p[j]->numExons-1].genTo;
        if (a < b)
          ed = b - a;
        else
          ed = a - b;

        if ((sd == 0) && (ed == 0)) {
          //fprintf(stderr, "%d and %d are exact; %d removed.\n", i, j, j);
          s4p_destroyPolish(p[j]);
          p[j] = 0L;
        } else if ((sd < 10) && (ed < 10)) {
          fprintf(stderr, "----------------------------------------\n");
          fprintf(stderr, "Warning: %d and %d are similar.\n", i, j);
          s4p_printPolish(stderr, p[i], S4P_PRINTPOLISH_FULL);
          s4p_printPolish(stderr, p[j], S4P_PRINTPOLISH_FULL);
          fprintf(stderr, "----------------------------------------\n");
        }
      }
    }
  }

  for (i=0; i<pNum; i++) {
    if (p[i]) {
      s4p_printPolish(stdout, p[i], S4P_PRINTPOLISH_FULL);
      s4p_destroyPolish(p[i]);
    }
  }
}

int
main(int argc, char **argv) {
  int          pNum   = 0;
  int          pAlloc = 8388608;
  sim4polish **p      = 0L;
  sim4polish  *q      = 0L;
  int          estID  = ~0;

  if (isatty(fileno(stdin))) {
    fprintf(stderr, usage, argv[0]);

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");

    exit(1);
  }

  //  Read polishes, picking the best when we see a change in
  //  the estID.

  p = (sim4polish **)malloc(sizeof(sim4polish *) * pAlloc);

  while ((q = s4p_readPolish(stdin)) != 0L) {
    if ((q->estID != estID) && (pNum > 0)) {
      pickBest(p, pNum);
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
    pickBest(p, pNum);

  return(0);
}


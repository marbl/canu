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

  //fprintf(stderr, "Pick for %d with %d things\n", p[0]->_estID, pNum);

  for (i=0; i<pNum; i++) {
    for (j=i+1; j<pNum; j++) {
      if ((p[i]) &&
          (p[j]) &&
          (p[i]->_numExons == p[j]->_numExons) &&
          (p[i]->_genID == p[j]->_genID)) {
        int a, b;
        int sd = 666;
        int ed = 666;

        a = p[i]->_exons[0]._genFrom;
        b = p[j]->_exons[0]._genFrom;
        if (a < b)
          sd = b - a;
        else
          sd = a - b;

        a = p[i]->_exons[p[i]->_numExons-1]._genTo;
        b = p[j]->_exons[p[j]->_numExons-1]._genTo;
        if (a < b)
          ed = b - a;
        else
          ed = a - b;

        if ((sd == 0) && (ed == 0)) {
          //fprintf(stderr, "%d and %d are exact; %d removed.\n", i, j, j);
          delete p[j];
          p[j] = 0L;
        } else if ((sd < 10) && (ed < 10)) {
          fprintf(stderr, "----------------------------------------\n");
          fprintf(stderr, "Warning: %d and %d are similar.\n", i, j);
          p[i]->s4p_printPolish(stderr, S4P_PRINTPOLISH_FULL);
          p[j]->s4p_printPolish(stderr, S4P_PRINTPOLISH_FULL);
          fprintf(stderr, "----------------------------------------\n");
        }
      }
    }
  }

  for (i=0; i<pNum; i++) {
    if (p[i]) {
      p[i]->s4p_printPolish(stdout, S4P_PRINTPOLISH_FULL);
      delete p[i];
    }
  }
}

int
main(int argc, char **argv) {
  u32bit   pNum   = 0;
  u32bit   pAlloc = 8388608;
  u32bit   estID  = ~u32bitZERO;

  if (isatty(fileno(stdin))) {
    fprintf(stderr, usage, argv[0]);

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");

    exit(1);
  }

  //  Read polishes, picking the best when we see a change in
  //  the estID.

  sim4polish **p = new sim4polish * [pAlloc];
  sim4polish  *q = new sim4polish (stdin);

  while (q->_numExons > 0) {
    if ((q->_estID != estID) && (pNum > 0)) {
      pickBest(p, pNum);
      pNum  = 0;
    }

    if (pNum >= pAlloc) {
      sim4polish **P = new sim4polish * [pAlloc * 2];
      memcpy(p, P, sizeof(sim4polish *) * pAlloc);
      delete [] p;
      p = P;
      pAlloc *= 2;
    }

    p[pNum++] = q;
    estID     = q->_estID;

    q = new sim4polish (stdin);
  }

  if (pNum > 0)
    pickBest(p, pNum);

  delete [] p;

  return(0);
}


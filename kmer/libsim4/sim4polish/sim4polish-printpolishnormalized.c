#include "sim4polish.h"

void
s4p_printPolishNormalized(FILE *O, sim4polish *p) {
  int  i;
  int  genLo = p->genLo;

  p->genLo = 0;

  for (i=0; i<p->numExons; i++) {
    p->exons[i].genFrom += genLo;
    p->exons[i].genTo   += genLo;
  }

  s4p_printPolish(O, p);

  for (i=0; i<p->numExons; i++) {
    p->exons[i].genFrom -= genLo;
    p->exons[i].genTo   -= genLo;
  }

  p->genLo = genLo;
}


#include "sim4polish.h"


void
s4p_removeAlignments(sim4polish *p) {
  int i;

  if (p) {
    for (i=0; i<p->numExons; i++) {
      free(p->exons[i].estAlignment);
      p->exons[i].estAlignment = 0L;

      free(p->exons[i].genAlignment);
      p->exons[i].genAlignment = 0L;
    }
  }  
}


void
s4p_removeDeflines(sim4polish *p) {
  int i;

  if (p) {
    free(p->estDefLine);
    p->estDefLine = 0L;

    free(p->genDefLine);
    p->genDefLine = 0L;
  }
}


void
s4p_destroyPolish(sim4polish *p) {
  int i;

  if (p) {
    for (i=0; i<p->numExons; i++) {
      free(p->exons[i].estAlignment);
      free(p->exons[i].genAlignment);
    }
    free(p->exons);
    free(p->comment);
    free(p->estDefLine);
    free(p->genDefLine);
    free(p);
  }
}

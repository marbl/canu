#include "sim4polish.h"


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

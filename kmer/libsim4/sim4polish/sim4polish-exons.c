#include "sim4polish.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void
s4p_swapExons(sim4polish *p, int a, int b) {
  sim4polishExon  copyofa;

  //  Save the a exon into copyofa; copy b into a; copy copyofa into b
  //
  memcpy(&copyofa,   p->exons+a, sizeof(sim4polishExon));
  memcpy(p->exons+a, p->exons+b, sizeof(sim4polishExon));
  memcpy(p->exons+b, &copyofa,   sizeof(sim4polishExon));
}

// p[a] --> e
void
s4p_copyExon(sim4polish *p, int a, sim4polishExon *e) {
  memcpy(e, p->exons+a, sizeof(sim4polishExon));
}


// e --> p[a]
void
s4p_overwriteExon(sim4polish *p, sim4polishExon *e, int a) {
  memcpy(p->exons+a, e, sizeof(sim4polishExon));
}


void
s4p_insertExon(sim4polish *p, int a, sim4polishExon *e) {
  fprintf(stderr, "s4p_insertExon() not implemented!\n");
}




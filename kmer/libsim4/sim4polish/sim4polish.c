#include "sim4polish.h"


int
s4p_makeForward(sim4polish *p) {
  if (p->matchOrientation == SIM4_MATCH_COMPLEMENT) {
    int e, t;
    for (e=0; e < p->numExons; e++) {
      t = p->estLen - p->exons[e].estFrom + 1;
      p->exons[e].estFrom = p->estLen - p->exons[e].estTo   + 1;
      p->exons[e].estTo = t;
    }
    p->matchOrientation = SIM4_MATCH_FORWARD;
    return(1);
  } else {
    return(0);
  }
}


int
s4p_makeReverse(sim4polish *p) {
  if (p->matchOrientation == SIM4_MATCH_FORWARD) {
    int e, t;
    for (e=0; e < p->numExons; e++) {
      t = p->estLen - p->exons[e].estFrom + 1;
      p->exons[e].estFrom = p->estLen - p->exons[e].estTo   + 1;
      p->exons[e].estTo = t;
    }
    p->matchOrientation = SIM4_MATCH_COMPLEMENT;
    return(1);
  } else {
    return(0);
  }
}

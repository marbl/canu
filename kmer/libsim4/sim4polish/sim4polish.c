#include "sim4polish.h"


//  Bunch of misc functions that should have a better home


void
s4p_reverseComplement(sim4polish *p) {
  int e, t;

  if (p->matchOrientation == SIM4_MATCH_COMPLEMENT) {
    p->matchOrientation = SIM4_MATCH_FORWARD;

    for (e=0; e < p->numExons; e++) {
      t = p->estLen - p->exons[e].estFrom + 1;
      p->exons[e].estFrom = p->estLen - p->exons[e].estTo   + 1;
      p->exons[e].estTo = t;
    }
  } else {
    p->matchOrientation = SIM4_MATCH_COMPLEMENT;
  }
}


int
s4p_makeForward(sim4polish *p) {
  if (p->matchOrientation == SIM4_MATCH_COMPLEMENT) {
    s4p_reverseComplement(p);
    return(1);
  } else {
    return(0);
  }
}


int
s4p_makeReverse(sim4polish *p) {
  if (p->matchOrientation == SIM4_MATCH_FORWARD) {
    s4p_reverseComplement(p);
    return(1);
  } else {
    return(0);
  }
}


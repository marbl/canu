#include "sim4polish.h"
#include "memory.h"

#include <errno.h>
#include <string.h>

sim4polish *
s4p_copyPolish(sim4polish *orig) {
  int         i;
  sim4polish *copy;

  copy = (sim4polish *)memdup(orig, sizeof(sim4polish));

  copy->estDefLine = (char *)memdup(orig->estDefLine, sizeof(char) * (strlen(orig->estDefLine) + 1));
  copy->genDefLine = (char *)memdup(orig->genDefLine, sizeof(char) * (strlen(orig->genDefLine) + 1));

  copy->exons = (sim4polishExon *)memdup(orig->exons, sizeof(sim4polishExon) * copy->numExons);

  for (i=0; i<copy->numExons; i++) {
    copy->exons[i].estAlignment = (char *)memdup(orig->exons[i].estAlignment, sizeof(char) * (strlen(copy->exons[i].estAlignment) + 1));
    copy->exons[i].genAlignment = (char *)memdup(orig->exons[i].genAlignment, sizeof(char) * (strlen(copy->exons[i].genAlignment) + 1));
  }

  return(copy);
}

#include "sim4polish.h"
#include "memory.h"

#include <errno.h>
#include <string.h>


sim4polish *
s4p_copyPolish(sim4polish *orig) {
  int         i;
  sim4polish *copy = 0L;

  if (orig == 0L)
    return(copy);;

  copy = (sim4polish *)memdup(orig, sizeof(sim4polish));

  if (orig->estDefLine)
    copy->estDefLine = (char *)memdup(orig->estDefLine, sizeof(char) * (strlen(orig->estDefLine) + 1));

  if (orig->genDefLine)
    copy->genDefLine = (char *)memdup(orig->genDefLine, sizeof(char) * (strlen(orig->genDefLine) + 1));

  copy->exons = (sim4polishExon *)memdup(orig->exons, sizeof(sim4polishExon) * copy->numExons);

  for (i=0; i<copy->numExons; i++) {
    if (orig->exons[i].estAlignment)
      copy->exons[i].estAlignment = (char *)memdup(orig->exons[i].estAlignment, sizeof(char) * (strlen(copy->exons[i].estAlignment) + 1));

    if (orig->exons[i].genAlignment)
      copy->exons[i].genAlignment = (char *)memdup(orig->exons[i].genAlignment, sizeof(char) * (strlen(copy->exons[i].genAlignment) + 1));
  }

  return(copy);
}


sim4polish *
s4p_copyPolish_OneExon(sim4polish *orig, int exon) {
  sim4polish *copy = 0L;

  if (orig == 0L)
    return(copy);;

  copy = (sim4polish *)memdup(orig, sizeof(sim4polish));

  if (orig->estDefLine)
    copy->estDefLine = (char *)memdup(orig->estDefLine, sizeof(char) * (strlen(orig->estDefLine) + 1));

  if (orig->genDefLine)
    copy->genDefLine = (char *)memdup(orig->genDefLine, sizeof(char) * (strlen(orig->genDefLine) + 1));

  copy->exons    = (sim4polishExon *)memdup(orig->exons + exon, sizeof(sim4polishExon));
  copy->numExons = 1;

  if (orig->exons[exon].estAlignment)
    copy->exons[0].estAlignment = (char *)memdup(orig->exons[exon].estAlignment, sizeof(char) * (strlen(orig->exons[exon].estAlignment) + 1));

  if (orig->exons[exon].genAlignment)
    copy->exons[0].genAlignment = (char *)memdup(orig->exons[exon].genAlignment, sizeof(char) * (strlen(orig->exons[exon].genAlignment) + 1));

  //  Rebuild stats
  copy->numMatches       = copy->exons[0].numMatches;
  copy->numMatchesN      = copy->exons[0].numMatchesN;
  copy->numCovered       = copy->exons[0].genTo - copy->exons[0].genFrom + 1;
  copy->percentIdentity  = copy->exons[0].percentIdentity;
  copy->querySeqIdentity = s4p_percentCoverageApprox(copy);

  return(copy);
}

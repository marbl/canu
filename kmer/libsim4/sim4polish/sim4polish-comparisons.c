#include "sim4polish.h"



int
s4p_estIDcompare(const void *a, const void *b) {
  sim4polish *A = (*(sim4polish **)a);
  sim4polish *B = (*(sim4polish **)b);

  if (A == 0L)  return(1);
  if (B == 0L)  return(-1);

  if (A->estID < B->estID) return(-1);
  if (A->estID > B->estID) return(1);
  if (A->genID < B->genID) return(-1);
  if (A->genID > B->genID) return(1);
  if (A->genLo < B->genLo) return(-1);
  if (A->genLo > B->genLo) return(1);
  if (A->exons[0].genFrom < B->exons[0].genFrom) return(-1);
  if (A->exons[0].genFrom > B->exons[0].genFrom) return(1);

  return(0);
}



int
s4p_genIDcompare(const void *a, const void *b) {
  sim4polish *A = (*(sim4polish **)a);
  sim4polish *B = (*(sim4polish **)b);

  if (A == 0L)  return(1);
  if (B == 0L)  return(-1);

  if (A->genID < B->genID) return(-1);
  if (A->genID > B->genID) return(1);
  if (A->estID < B->estID) return(-1);
  if (A->estID > B->estID) return(1);
  if (A->genLo < B->genLo) return(-1);
  if (A->genLo > B->genLo) return(1);
  if (A->exons[0].genFrom < B->exons[0].genFrom) return(-1);
  if (A->exons[0].genFrom > B->exons[0].genFrom) return(1);

  return(0);
}


#include "util++.H"
#include "sim4.H"

//  Build an interval list with all exons (from both guys), merge
//  overlapping regions, compute the length, subtract from the total.
//  Result: the number of bp that the two matches overlap in the
//  genomic.
//
u32bit
findOverlap(sim4polish *A, sim4polish *B) {

  if ((A->genID != B->genID) || (A->matchOrientation != B->matchOrientation))
    return(0);

  u32bit        length = 0;
  u32bit        total  = 0;
  intervalList  IL;

  for (u32bit i=0; i<A->numExons; i++) {
    length = A->exons[i].genTo - A->exons[i].genFrom + 1;
    total  += length;
    IL.add(A->genLo + A->exons[i].genFrom, length);
  }

  for (u32bit i=0; i<B->numExons; i++) {
    length = B->exons[i].genTo - B->exons[i].genFrom + 1;
    total  += length;
    IL.add(B->genLo + B->exons[i].genFrom, length);
  }

  IL.merge();

  return(total - IL.sumOfLengths());
}

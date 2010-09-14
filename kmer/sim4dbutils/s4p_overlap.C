#include "util++.H"
#include "sim4.H"

//  Build an interval list with all exons (from both guys), merge
//  overlapping regions, compute the length, subtract from the total.
//  Result: the number of bp that the two matches overlap in the
//  genomic.
//
u32bit
findOverlap(sim4polish *A, sim4polish *B) {

  if ((A->_genID != B->_genID) || (A->_matchOrientation != B->_matchOrientation))
    return(0);

  u32bit        length = 0;
  u32bit        total  = 0;
  intervalList  IL;

  for (u32bit i=0; i<A->_numExons; i++) {
    length = A->_exons[i]._genTo - A->_exons[i]._genFrom + 1;
    total  += length;
    IL.add(A->_exons[i]._genFrom, length);
  }

  for (u32bit i=0; i<B->_numExons; i++) {
    length = B->_exons[i]._genTo - B->_exons[i]._genFrom + 1;
    total  += length;
    IL.add(B->_exons[i]._genFrom, length);
  }

  IL.merge();

#ifdef OLAP_IS_SHORT
  if (total - IL.sumOfLengths() > 65536) {
    fprintf(stderr, "findOverlap()-- ERROR!  The overlap is larger than the return type!\n");
    fprintf(stderr, "findOverlap()--         Switch to 32-bit ints in s4p_overlap.H.\n");
  }
#endif

  return(total - IL.sumOfLengths());
}

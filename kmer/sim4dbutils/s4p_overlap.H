#ifndef S4P_OVERLAP_H
#define S4P_OVERLAP_H

//  Using 16-bit ints for storing the amount overlapped gives a big
//  memory reduction, but will fail for long sequences (mRNA,
//  probably).  findOverlap() checks for overflow.
//
#define OLAP_IS_SHORT

#ifdef OLAP_IS_SHORT
typedef  uint16    olap_t;
#define  OLAPTFMT  uint16FMT
#else
typedef  uint32    olap_t;
#define  OLAPTFMT  uint32FMT
#endif

olap_t   findOverlap(sim4polish *A, sim4polish *B);

#endif  //  S4P_OVERLAP_H

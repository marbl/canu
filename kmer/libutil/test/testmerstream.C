#include <stdio.h>
#include <stdlib.h>

#include "merstream.H"

//  This is the input file
//
#if 0
>(one)
TTTTTTTTTTAAAAAAAAAACGNTTTTTTTTTTNGGGGGGGGGGNAAAAAAAAANTTTTTTTTTT
>(two)
TTTTTTTTAA
>(thr)
TTTTTTTTAC
>(fou)
T       T       T       T       T
T       T       T       A       G
>(fiv)
T T T T T T T T A T

>(six)

T
T
T
T
T
T
T
T
T
T
#endif

//  This is the correct output
//
u64bit correct[21] = {
  0x00000000000fffff,
  0x00000000000ffffc,
  0x00000000000ffff0,
  0x00000000000fffc0,
  0x00000000000fff00,
  0x00000000000ffc00,
  0x00000000000ff000,
  0x00000000000fc000,
  0x00000000000f0000,
  0x00000000000c0000,
  0x0000000000000000,
  0x0000000000000001,
  0x0000000000000006,
  0x00000000000fffff,
  0x00000000000aaaaa,
  0x00000000000fffff,
  0x00000000000ffff0,
  0x00000000000ffff1,
  0x00000000000ffff2,
  0x00000000000ffff3,
  0x00000000000fffff
};


int
main(int argc, char **argv) {
  merStream  *S = new merStream(10, "1.fasta");
  u32bit      c = 0;

  while (S->nextMer()) {
    fprintf(stderr, "0x%016lx '%s' pos=%4lu seq=%4lu\n",
            S->theFMer(), S->theDefLine(), S->thePosition(), S->theSequenceNumber());
    if (S->theFMer() != correct[c])
      fprintf(stderr, "0x%016lx != 0x%016lx at position %2u.\n", S->theFMer(), correct[c], c);
    c++;
  }

  if (c != 21)
    fprintf(stderr, "Wrong number.  Should be 21 was %u\n", c);
}


#include <stdio.h>
#include "intervalList.H"

int
main(int argc, char **argv) {
  intervalList  I;

  I.add( 5,  3);
  I.add( 5,  3);
  I.add( 5,  3);
  I.add(10,  5);
  I.add(15,  5);
  I.add(20, 10);
  I.add(30, 10);
  I.add(32,  5);
  I.add(50, 10);
  I.add(55, 10);
  I.add(70,  3);
  I.add(71,  3);
  I.add(72,  3);
  I.add(73,  3);

  for (u32bit i=0; i<I.numberOfIntervals(); i++)
    fprintf(stderr, "%2d] %2d %2d\n", i, I.lo(i), I.hi(i));

  fprintf(stderr, "sorting\n");
  I.sort();
  for (u32bit i=0; i<I.numberOfIntervals(); i++)
    fprintf(stderr, "%2d] %2d %2d\n", i, I.lo(i), I.hi(i));

  fprintf(stderr, "merging\n");
  I.merge();
  for (u32bit i=0; i<I.numberOfIntervals(); i++)
    fprintf(stderr, "%2d] %2d %2d\n", i, I.lo(i), I.hi(i));

  fprintf(stderr, "Length: %u\n", I.sumOfLengths());
}

#include <stdio.h>

#include "bri++.H"

int
main(int argc, char **argv) {
  int           e = 0;
  intervalList  I;

  I.add(71,  3);
  I.add( 5,  3);
  I.add(32,  5);
  I.add(73,  3);
  I.add(55, 10);
  I.add( 5,  3);
  I.add(10,  5);
  I.add(20, 10);
  I.add(30, 10);
  I.add(50, 10);
  I.add(70,  3);
  I.add(72,  3);
  I.add( 5,  3);
  I.add(15,  5);

#if 0
  for (u32bit i=0; i<I.numberOfIntervals(); i++)
    fprintf(stderr, "%2d] %2d %2d\n", i, I.lo(i), I.hi(i));
#endif

  I.sort();
  I.merge();

  if (I.sumOfLengths() != 54)
    fprintf(stderr, "Lengths don't add up.\n"), e++;

  if (I.numberOfIntervals() != 4)
    fprintf(stderr, "Wrong number of intervals.\n"), e++;

  if ((I.lo(0) != 5) || (I.hi(0) != 8))
    fprintf(stderr, "Interval 0 is wrong.\n"), e++;
  if ((I.lo(1) != 10) || (I.hi(1) != 40))
    fprintf(stderr, "Interval 1 is wrong.\n"), e++;
  if ((I.lo(2) != 50) || (I.hi(2) != 65))
    fprintf(stderr, "Interval 2 is wrong.\n"), e++;
  if ((I.lo(3) != 70) || (I.hi(3) != 76))
    fprintf(stderr, "Interval 3 is wrong.\n"), e++;

  exit(e);
}

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#define SIZE (1073741824)
#define ITER (32 * 1024 * 1024)

int
main(int argc, char **argv) {
  char  *bits = new char [SIZE];
  int    coll = 0;

  bzero(bits, SIZE);
  srand48(time(NULL));

  for (int i=0; i<ITER; i++) {
    int  p = lrand48() % SIZE;
    if (bits[p])
      coll++;
    bits[p] = 1;
  }

  fprintf(stderr, "collisions: %d\n", coll);

  double expected = ITER / 2;
  expected *= (ITER-1);
  expected /= SIZE;

  fprintf(stderr, "expected:   %f\n", expected); 
}

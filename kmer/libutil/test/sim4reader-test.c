#include <stdio.h>
#include <stdlib.h>
#include "sim4reader.h"

int
main(int argc, char **argv) {
  sim4polish  *p;

  while (!feof(stdin)) {
    p = readPolish(stdin);

    if (p)
      if (argc > 1)
        printPolishColumn(stdout, p);
      else
        printPolish(stdout, p);

    destroyPolish(p);
  }

  return(0);
}

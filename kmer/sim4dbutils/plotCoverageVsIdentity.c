#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "sim4reader.h"

int
main(int argc, char ** argv) {
  sim4polish  *p;
  int          c[101] = {0};
  int          i[101] = {0};
  int          x;
  FILE        *C;
  FILE        *I;
  FILE        *S;

  if (isatty(fileno(stdin))) {
    fprintf(stderr, "creates three files:\n");
    fprintf(stderr, "  coverage.histogram\n");
    fprintf(stderr, "  identity.histogram\n");
    fprintf(stderr, "  c-vs-i.scatter\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");
    exit(1);
  }

  C = fopen("coverage.histogram", "w");
  I = fopen("identity.histogram", "w");
  S = fopen("c-vs-i.scatter", "w");

  while ((p = readPolish(stdin)) != 0L) {
    fprintf(S, "%d %d\n", p->percentIdentity, p->querySeqIdentity);

    i[p->percentIdentity]++;
    c[p->querySeqIdentity]++;

    destroyPolish(p);
  }

  for (x=0; x<101; x++) {
    fprintf(C, "%d\n", c[x]);
    fprintf(I, "%d\n", i[x]);
  }

  fclose(C);
  fclose(I);
  fclose(S);

  return(0);
}

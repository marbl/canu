#include <stdio.h>

#include "outputString.H"

int
main(int argc, char **argv) {
  FILE          *F;
  char           c[129];
  int            l;
  outputString   O;

  if (argc != 2) {
    fprintf(stderr, "usage: %s <file>\n", argv[0]);
    exit(1);
  }

  F = fopen(argv[1], "r");
  while (!feof(F)) {
    l = fread(c, sizeof(char), 128, F);

    //fprintf(stderr, "l=%d\n", l);

    c[l] = 0;

    O.appendString(c);
  }
  fclose(F);

  fputs(O.getString(), stdout);
}

#include <stdio.h>
#include <stdlib.h>

#include "merstream.H"
#include "util.H"

int
main(int argc, char **argv) {
  merStreamFromFile  *S = new merStreamFromFile(10, argv[1]);
  double              s = getTime();
  u64bit              c = 0;

  s = getTime();
  c = 0;
  FILE *F = fopen(argv[1], "r");
  while (fgetc(F) != EOF) {
    if ((++c & (u64bit)0xffffff) == u64bitZERO) {
      fprintf(stderr, " %4ld Mbases -- %8.5f Mb per second\r",
              c >> 20,
              c / (getTime() - s) / (1000000.0));
      fflush(stderr);
    }
  }
  fprintf(stderr, "\n");

  s = getTime();
  c = 0;
  FastAstream   *_theFile = new FastAstream(argv[1]);
  while (_theFile->nextSymbol() != 255) {
    if ((++c & (u64bit)0xffffff) == u64bitZERO) {
      fprintf(stderr, " %4ld Mbases -- %8.5f Mb per second\r",
              c >> 20,
              c / (getTime() - s) / (1000000.0));
      fflush(stderr);
    }
  }
  fprintf(stderr, "\n");

  s = getTime();
  c = 0;
  while (S->nextMer()) {
    if ((++c & (u64bit)0xffffff) == u64bitZERO) {
      fprintf(stderr, " %4ld Mbases -- %8.5f Mb per second\r",
              c >> 20,
              c / (getTime() - s) / (1000000.0));
      fflush(stderr);
    }
  }
  fprintf(stderr, "\n");
}


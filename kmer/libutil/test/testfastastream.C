#include <stdio.h>
#include <stdlib.h>

#include "fastastream.H"

int
main(int argc, char **argv) {
  FastAstream  *F = new FastAstream("1.fasta");

  unsigned char   ch;

  for (ch = F->nextSymbol(); ch != 255; ch = F->nextSymbol()) {
    if (ch == 254)
      fprintf(stderr, "\n%s\n", F->theDefLine());

    if (!(ch & 0x80))
      fprintf(stderr, "%c", ch);
  }
  fprintf(stderr, "\n");
}

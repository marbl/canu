#include <stdio.h>
#include <stdlib.h>
#include "bio++.H"

//  Constructs a seqStore, and counts the number of ACGT in it.

int
main(int argc, char **argv) {

  if (argc < 2) {
    fprintf(stderr, "usage: %s some.fasta\n", argv[0]);
    exit(1);
  }

  seqStream *ST = new seqStream(argv[1], true);
  seqStore  *SS = new seqStore(argv[1], ST);

  fprintf(stdout, "Found "u64bitFMT" ACGT in '%s'\n",
          SS->numberOfACGT(), argv[1]);

  delete SS;
  delete ST;

  exit(0);
}

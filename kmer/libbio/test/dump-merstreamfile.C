#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bio++.H"

//  Dumps each mer in a merStreamFile.

#define BUILD_SIZE 10
#define TEST_SIZE  10

int
main(int argc, char **argv) {
  char   str[1024];

  if (argc != 2) {
    fprintf(stderr, "usage: %s some.fasta\n", argv[0]);
    fprintf(stderr, "       Builds a merStreamFile, writes the contents to stdout.\n");
    exit(1);
  }

  merStreamFileBuilder   *B = new merStreamFileBuilder(BUILD_SIZE, argv[1], "junk");
  u64bit numMers = B->build(true);
  delete B;

  fprintf(stderr, "Found "u64bitFMT" mers in %s at mersize %d\n", numMers, argv[1], BUILD_SIZE);

  merStreamFileReader *R = new merStreamFileReader("junk", TEST_SIZE);
  merStream           *M = new merStream(R);

  while (M->nextMer())
    fprintf(stdout, "'%s' seqNum:"u64bitFMT" posInSeq:"u64bitFMT" posInStream:"u64bitFMT" '%s'\n",
            M->theFMer().merToString(str),
            M->theSequenceNumber(),
            M->thePositionInSequence(),
            M->thePositionInStream(),
            M->theDefLine());

  exit(0);
}

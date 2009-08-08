#include <stdio.h>
#include <stdlib.h>

#include "bio++.H"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

int
main(int argc, char **argv) {
  speedCounter          *C    = 0L;
  FILE                  *F    = 0L;
  seqStream             *S    = 0L;
  merStream             *M    = 0L;

  if (argc != 2) {
    fprintf(stderr, "usage: %s some.fasta\n", argv[0]);
    fprintf(stderr, "Reads some.fasta using fgetc(), the seqStream and the merStream,\n");
    fprintf(stderr, "reporting the speed of each method.\n");
    exit(1);
  }

  ////////////////////////////////////////
  F = fopen(argv[1], "r");
  C = new speedCounter("fgetc():        %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (!feof(F))
    fgetc(F), C->tick();
  delete C;
  fclose(F);

  ////////////////////////////////////////
  S = new seqStream(argv[1]);
  C = new speedCounter("seqStream:      %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (S->get())
    C->tick();
  delete C;
  delete S;

  ////////////////////////////////////////
  M = new merStream(new kMerBuilder(20),
                    new seqStream(argv[1]),
                    true, true);
  C = new speedCounter("seqStream -> merStream:   %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (M->nextMer())
    C->tick();
  delete C;
  delete M;

  exit(0);
}


#include <stdio.h>
#include <stdlib.h>

#include "bri++.H"

int
main(int argc, char **argv) {
  speedCounter          *C = 0L;
  FILE                  *F = 0L;
  FastAstream           *S = 0L;
  merStream             *M = 0L;
  merStreamFileBuilder  *msfB;
  merStreamFileReader   *msfR;


  if (argc != 2) {
    fprintf(stderr, "usage: %s some.fasta\n", argv[0]);
    fprintf(stderr, "Reads some.fasta using fgetc(), the FastAstream and the merStream,\n");
    fprintf(stderr, "reporting the speed of each method.\n");
    exit(1);
  }

  F = fopen(argv[1], "r");
  C = new speedCounter("fgetc():       %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (!feof(F))
    fgetc(F), C->tick();
  delete C;
  fclose(F);

  S = new FastAstream(argv[1]);
  C = new speedCounter("FastAstream:   %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (S->nextSymbol())
    C->tick();
  delete C;
  delete S;

  M = new merStream(20, argv[1]);
  C = new speedCounter("merStream:     %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (M->nextMer())
    C->tick();
  delete C;
  delete M;

  msfB = new merStreamFileBuilder(20, argv[1], "merStreamFileTest");
  msfB->build(true);
  delete msfB;

  msfR = new merStreamFileReader("merStreamFileTest");
  C    = new speedCounter("merStreamFile: %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (msfR->nextMer())
    C->tick();
  delete C;
  delete msfR;

  exit(0);
}


#include <stdio.h>
#include <stdlib.h>

#include "bio++.H"

int
main(int argc, char **argv) {
  speedCounter          *C    = 0L;
  FILE                  *F    = 0L;
  seqStream             *S    = 0L;
  merStream             *M    = 0L;
  merStreamFileBuilder  *msfB = 0L;
  merStreamFileReader   *msfR = 0L;
  chainedSequence       *cs   = 0L;

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
  C = new speedCounter("FS:             %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (S->nextSymbol())
    C->tick();
  delete C;
  delete S;

  ////////////////////////////////////////
  msfB = new merStreamFileBuilder(20, argv[1], "junk");
  msfB->build(true);
  delete msfB;

  msfR = new merStreamFileReader("junk");
  C    = new speedCounter("MSF:           %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (msfR->nextMer())
    C->tick();
  delete C;
  delete msfR;

  ////////////////////////////////////////
  msfB = new merStreamFileBuilder(20, argv[1], "junk");
  msfB->build(true);
  delete msfB;

  msfR = new merStreamFileReader("junk");
  M = new merStream(msfR);
  C    = new speedCounter("MSF -> MS:      %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (M->nextMer())
    C->tick();
  delete C;
  delete msfR;
  delete M;

  ////////////////////////////////////////
  cs = new chainedSequence();
  cs->setSource(argv[1]);
  cs->finish();

  C = new speedCounter("CS:             %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (cs->get())
    C->tick();
  delete C;
  delete cs;

  ////////////////////////////////////////
  cs = new chainedSequence();
  cs->setSource(argv[1]);
  cs->finish();

  M = new merStream(20, cs);
  C = new speedCounter("CS -> MS:       %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (M->nextMer())
    C->tick();
  delete C;
  delete M;
  delete cs;


  exit(0);
}


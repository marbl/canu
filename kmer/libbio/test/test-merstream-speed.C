#include <stdio.h>
#include <stdlib.h>

#include "bio++.H"

//  Results on oi not so great.
//
//  fgetc():         55.06 Mthings --  7.07 Mthings/second
//  FastAstream:     55.00 Mthings --  9.80 Mthings/second
//  merStream:       54.80 Mthings --  4.05 Mthings/second
//  merStreamFile:   54.80 Mthings --  3.27 Mthings/second
//  chainedSeq -> FastAstream -> merStream:    54.99 Mthings --  3.14 Mthings/second
//  chainedSeq -> FastAstream -> merStream:    54.99 Mthings --  3.42 Mthings/second
//
//  On highfive better (not the same code version as above)
//
//  fgetc():          55.25 Mthings -- 24.75 Mthings/second
//  FS:               55.19 Mthings -- 26.34 Mthings/second
//  FS -> MS:         54.99 Mthings -- 10.65 Mthings/second
//  CS:               55.29 Mthings -- 25.07 Mthings/second
//  CS -> FS:         55.29 Mthings -- 17.26 Mthings/second
//  CS -> MS:         54.99 Mthings -- 10.62 Mthings/second
//  CS -> FS -> MS:   54.99 Mthings --  8.66 Mthings/second
//

int
main(int argc, char **argv) {
  speedCounter          *C    = 0L;
  FILE                  *F    = 0L;
  FastAstream           *S    = 0L;
  merStream             *M    = 0L;
  merStreamFileBuilder  *msfB = 0L;
  merStreamFileReader   *msfR = 0L;
  chainedSequence       *cs   = 0L;


  if (argc != 2) {
    fprintf(stderr, "usage: %s some.fasta\n", argv[0]);
    fprintf(stderr, "Reads some.fasta using fgetc(), the FastAstream and the merStream,\n");
    fprintf(stderr, "reporting the speed of each method.\n");
    exit(1);
  }


  F = fopen(argv[1], "r");
  C = new speedCounter("fgetc():        %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (!feof(F))
    fgetc(F), C->tick();
  delete C;
  fclose(F);


  S = new FastAstream(argv[1]);
  C = new speedCounter("FS:             %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (S->nextSymbol())
    C->tick();
  delete C;
  delete S;


  S = new FastAstream(argv[1]);
  M = new merStream(20, S);
  C = new speedCounter("FS -> MS:       %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (M->nextMer())
    C->tick();
  delete C;
  delete M;
  delete S;

#if 1
  msfB = new merStreamFileBuilder(20, argv[1], "junk");
  msfB->build(true);
  delete msfB;

  msfR = new merStreamFileReader("junk");
  C    = new speedCounter("MSF:         %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (msfR->nextMer())
    C->tick();
  delete C;
  delete msfR;
#endif


  cs = new chainedSequence();
  cs->setSource(argv[1]);
  cs->finish();

  C = new speedCounter("CS:             %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (cs->get())
    C->tick();
  delete C;
  delete cs;

  cs = new chainedSequence();
  cs->setSource(argv[1]);
  cs->finish();

  S = new FastAstream(cs);
  C = new speedCounter("CS -> FS:       %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (S->nextSymbol())
    C->tick();
  delete C;
  delete S;
  delete cs;


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


  cs = new chainedSequence();
  cs->setSource(argv[1]);
  cs->finish();

  S = new FastAstream(cs);
  M = new merStream(20, S);
  C = new speedCounter("CS -> FS -> MS: %7.2f Mthings -- %5.2f Mthings/second\r", 1000000.0, 0x3fffff, true);
  while (M->nextMer())
    C->tick();
  delete C;
  delete M;
  delete S;
  delete cs;


  exit(0);
}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bio++.H"

#if 0

const char *sequence =
">(one)\n"
"TTTTTTTTTT\n"
"AAAAAAAAAACGNTTTTTTTTTTNGGGGGGGGGGNAAAAAAAAANTTTTTTTTTT\n"
">(two)\n"
"TTTTTTTTAA\n"
">(thr)\n"
"TTTTTTTTAC\n"
">(fou)\n"
"T       T       T       T       T\n"
"T       T       T       A       G\n"
">(fiv)\n"
"T T T T T T T T A T C\n"
"\n"
">(six)\n"
"\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n";

#else

const char *sequence =
">one\n"
"AAAAAAAAAAAAAAAAAAAA\n"
">two\n"
"CCCCCCCCCCCCCCCCCCCC\n"
">three\n"
"GGGGGGGGGGGGGGGGGGGG\n";

#endif



int
main(int argc, char **argv) {
  chainedSequence *CS  = 0L;
  FastAstream     *F1  = 0L;
  FastAstream     *F2  = 0L;
  unsigned char    a   = 1;
  unsigned char    b   = 1;
  int              err = 0;

  FILE *F = fopen("test.junk", "w");
  fwrite(sequence, sizeof(char), strlen(sequence), F);
  fclose(F);

#if 1
  CS = new chainedSequence();
  CS->setSource("test.junk");
  CS->finish();
  
  char x = CS->get();
  while (x) {
    fprintf(stderr, "%c", x);
    x = CS->get();
  }

  fprintf(stderr, "\n");

  delete CS;
#endif

  CS = new chainedSequence();
  CS->setSource("test.junk");
  CS->finish();

  F1 = new FastAstream("test.junk");
  F2 = new FastAstream(CS);

  a = F1->nextSymbol();
  b = F2->nextSymbol();
  while ((a != 0) && (b != 0)) {

    //  Skip the _sequence_ _breaks_ in the fasta sourced stream
    if (a == 254)
      a = F1->nextSymbol();

    //  Skip the _sequence_ _gaps_ in the chainedSequence sourced stream
    while (b == 253)
      b = F2->nextSymbol();

    if (a != b)
      fprintf(stderr, "%c != %c\n", a, b);

    a = F1->nextSymbol();
    b = F2->nextSymbol();
  }

  if (a)
    fprintf(stderr, "FastAstream(filename) has more stuff (%d %c)!\n", a, a), err=1;
  if (b)
    fprintf(stderr, "FastAstream(chainedSequence) has more stuff (%d %c)!\n", b, b), err=1;

  delete F2;
  delete F1;
  delete CS;

  unlink("test.junk");
  unlink("test.junkidx");

  if (err)
    fprintf(stderr, "Test failed.\n");
  else
    fprintf(stderr, "Test passed.\n");

  return(err);
}


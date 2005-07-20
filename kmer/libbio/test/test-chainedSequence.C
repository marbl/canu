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
  chainedSequence *CT  = 0L;
  FastAstream     *F1  = 0L;
  FastAstream     *F2  = 0L;
  FastAstream     *F3  = 0L;
  unsigned char    a   = 1;
  unsigned char    b   = 1;
  unsigned char    c   = 1;
  int              err = 0;
  char const      *filename = 0L;

  if (argc == 1) {
    FILE *F = fopen("junk.fasta", "w");
    fwrite(sequence, sizeof(char), strlen(sequence), F);
    fclose(F);

    filename = "junk.fasta";
  } else if (argc == 2) {
    filename = argv[1];
  } else {
    fprintf(stderr, "usage: %s <some.fasta>\n", argv[0]);
    fprintf(stderr, "  if <some.fasta>, use that for testing, else use internal sequence\n");
    exit(1);
  }

#if 0
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
  CS->setSource(filename);
  CS->finish();

  CS->saveState("junk.chainedSequence");

  CT = new chainedSequence();
  CT->loadState("junk.chainedSequence");
  CT->saveState("junk.copied.chainedSequence");

  F1 = new FastAstream(filename);
  F2 = new FastAstream(CS);
  F3 = new FastAstream(CT);

  a = F1->nextSymbol();
  b = F2->nextSymbol();
  c = F3->nextSymbol();

  while ((a != 0) && (b != 0) && (c != 0)) {

    //  Skip the _sequence_ _breaks_ in the fasta sourced stream
    if (a == 254)
      a = F1->nextSymbol();

    //  Skip the _sequence_ _gaps_ in the chainedSequence sourced stream
    while (b == 253)
      b = F2->nextSymbol();
    while (c == 253)
      c = F3->nextSymbol();

    if ((a != b) || (a != c) || (b != c))
      fprintf(stderr, "FS:%c != CSoriginal:%c != CSrestored:%c\n", a, b, c);

    a = F1->nextSymbol();
    b = F2->nextSymbol();
    c = F3->nextSymbol();
  }

  if (a)
    fprintf(stderr, "FastAstream(filename) has more stuff (%d %c)!\n", a, a), err=1;
  if (b)
    fprintf(stderr, "FastAstream(chainedSequence) original has more stuff (%d %c)!\n", b, b), err=1;
  if (c)
    fprintf(stderr, "FastAstream(chainedSequence) restored has more stuff (%d %c)!\n", c, c), err=1;

  delete F2;
  delete F1;
  delete CS;
  delete CT;

  if (err)
    fprintf(stderr, "Test failed.\n");
  else
    fprintf(stderr, "Test passed.\n");

  return(err);
}


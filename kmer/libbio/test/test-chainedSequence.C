#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
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


#if 0  //  redundant, nothing here



int
main(int argc, char **argv) {
  seqStream       *CS  = 0L;
  seqStream       *CT  = 0L;
  seqStream       *F1  = 0L;
  seqStream       *F2  = 0L;
  seqStream       *F3  = 0L;
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
    CS = new seqStream();
    CS->setFile("test.junk", true);
  
    char x = CS->get();
    while (x) {
      fprintf(stderr, "%c", x);
      x = CS->get();
    }

    fprintf(stderr, "\n");

    delete CS;
#endif

#if 0
  //  We don't do saveState/loadState anymore
  CS = new seqStream(filename, true);

  CS->saveState("junk.seqStream");

  CT = new seqStream();
  CT->loadState("junk.seqStream");
  CT->saveState("junk.copied.seqStream");
#else
  //  Which makes this test kind of pointless
  CS = new seqStream(filename, true);
  CT = new seqStream(filename, true);
#endif

  F1 = new seqStream(filename, true);
  F2 = new seqStream(CS);
  F3 = new seqStream(CT);

  a = F1->nextSymbol();
  b = F2->nextSymbol();
  c = F3->nextSymbol();

  while ((a != 0) && (b != 0) && (c != 0)) {
    while ((a == 254) || (a == 253))
      a = F1->nextSymbol();
    while ((b == 254) || (b == 253))
      b = F2->nextSymbol();
    while ((c == 254) || (c == 253))
      c = F3->nextSymbol();

    if ((a != b) || (a != c) || (b != c))
      fprintf(stderr, "FS:%c(%3d) != CSoriginal:%c(%3d) != CSrestored:%c(%3d)\n", a, a, b, b, c, c);

    a = F1->nextSymbol();
    b = F2->nextSymbol();
    c = F3->nextSymbol();
  }

  if (a)
    fprintf(stderr, "seqStream(filename) has more stuff (%d %c)!\n", a, a), err=1;
  if (b)
    fprintf(stderr, "seqSstream(seqStream) original has more stuff (%d %c)!\n", b, b), err=1;
  if (c)
    fprintf(stderr, "seqStream(seqStream) restored has more stuff (%d %c)!\n", c, c), err=1;

  delete F2;
  delete F1;
  delete CS;
  delete CT;

  if (err)
    fprintf(stderr, "Test failed.\n");
  else
    fprintf(stderr, "Test passed.\n");

  unlink("junk.seqStream");
  unlink("junk.copied.seqStream");

  return(err);
}


#else  //  redundant, nothing here

int main(int argc, char **argv) {
  return(0);
}

#endif  //  redundant, nothing here

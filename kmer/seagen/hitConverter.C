#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include "aHit.H"

//  $Id$

#ifdef TRUE64BIT
char const *hitStr = "%u hits.\r";
#else
char const *hitStr = "%lu hits.\r";
#endif

void
bin2asc(FILE *I, FILE *O) {
  u32bit  i = 0;
  aHit    a;

  fprintf(stderr, "Converting BINARY to ASCII.\n");

  while (!feof(I)) {
    ahit_readBinary(&a, I);

    if (!feof(I)) {
      ahit_printASCII(&a, O);

      if ((++i & 0xffff) == 0) {
        fprintf(stderr, hitStr, i);
        fflush(stderr);
      }
    }
  }

  fprintf(stderr, hitStr, i);
  fprintf(stderr, "\n");
}


void
asc2bin(FILE *I, FILE *O) {
  u32bit  i = 0;
  aHit    a;
  char    b[1025];

  fprintf(stderr, "Converting ASCII to BINARY.\n");

  while (!feof(I)) {
    fgets(b, 1024, I);

    if (!feof(I)) {
      ahit_parseString(&a, b);
      ahit_writeBinary(&a, O);

      if ((++i & 0xffff) == 0) {
        fprintf(stderr, hitStr, i);
        fflush(stderr);
      }
    }
  }

  fprintf(stderr, hitStr, i);
  fprintf(stderr, "\n");
}


int
main(int argc, char **argv) {

  if (argc != 1) {
    fprintf(stderr, "%s: I only read stdin and write stdout.\n", argv[0]);
    exit(1);
  }

  //  If the first character in the stream is ascii, do ASCII -> BINARY.
  //  Else, do BINARY -> ASCII.
  //
  char x = (char)fgetc(stdin);
  ungetc(x, stdin);

  if (x == '-')
    asc2bin(stdin, stdout);
  else
    bin2asc(stdin, stdout);

  return(0);
}

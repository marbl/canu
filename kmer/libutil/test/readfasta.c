#include <stdio.h>
#include <stdlib.h>

#include "fasta-c.h"

//  You need to compile your source files with:
//    cc -c -o somefile.o somefile.c -I/home/walenzbp/projects/libbri
//
//  The executable must be linked by the C++ compiler:
//    cxx -o program $(OBJS) -L/home/walenzbp/projects/libbri -lbri
//
//  cc -c -o readfasta.o readfasta.c
//  cxx -o readfasta readfasta.o -L. -lbri
//
void
main(int argc, char **argv) {
  void  *f;
  int    i;

  if (argc != 3) {
    fprintf(stderr, "usage: %s seq.fasta number-of-seq-to-get\n", argv[0]);
    exit(1);
  }

  i = atoi(argv[2]);

  //  this creates the FastA reader "object".  Do it ONCE.
  //
  f = createFastA(argv[1]);

  //  These will extract a single sequence, length and header from the
  //  fasta
  //
  fprintf(stdout, "%d\n", getFastAsequenceLength(f, i));
  fprintf(stdout, "%s\n", getFastAheader(f, i));
  fprintf(stdout, "%s\n", getFastAsequence(f, i));

  //  When you're all done reading sequences, destroy the FastA
  //  reader.
  //
  destroyFastA(f);
}

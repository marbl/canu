#include <stdio.h>
#include "/home/walenzbp/projects/fasta/fasta-c.h"

//
//  Test program for the C wrapper of FastAWrapper.
//

//  cxx -c -o fasta-c-test.o fasta-c-test.c
//  cxx -o fasta-c-test fasta-c-test.o -L/home/walenzbp/projects/fasta -L/home/walenzbp/projects/libbri -lfasta -lbri

void
main(int argc, char **argv) {
  void *fastafile;

  fastafile = createFastA("test-50bp.fasta");

  fprintf(stderr, "length =  %u\n",  getFastAsequenceLength(fastafile, 10));
  fprintf(stderr, "header = '%s'\n", getFastAheader(fastafile, 10));
  fprintf(stderr, "seq    = '%s'\n", getFastAsequence(fastafile, 10));

  fprintf(stderr, "length =  %u\n",  getFastAsequenceLength(fastafile, 4));
  fprintf(stderr, "header = '%s'\n", getFastAheader(fastafile, 4));
  fprintf(stderr, "seq    = '%s'\n", getFastAsequence(fastafile, 4));

  destroyFastA(fastafile);
}

#include <stdio.h>
#include <stdlib.h>

#include "util++.H"
#include "bio++.H"

//  Sorts 100-mers in an input sequence, using qsort

#define MAXBASES (1024 * 1024 * 1024)

char *SA;

int
compare100(const void *a, const void *b) {
  u32bit A = *(u32bit *)a;
  u32bit B = *(u32bit *)b;
  return(strncmp(SA+A, SA+B, 100));
}


int
main(int argc, char **argv) {
  
  if (argc != 2) {
    fprintf(stderr, "usage: %s seq.fasta\n", argv[0]);
    exit(1);
  }

  //  Read all sequence, put into a single string.  Separate sequences
  //  with a dot.  Hopefully, use the chainedSequence, if that
  //  actually keeps the sequence in memory.  Or not.  Instead
  //  use the fastastream, which does separate sequences, and automagically
  //  filters out N!
  //
  FastAstream  *FS = new FastAstream(argv[1]);

  fprintf(stderr, "Copying string.\n");

  //  Guess the string length, so we can copy the FS into a single character
  //  array.
  //
  SA = new char [MAXBASES + 2 * 1024 * 1024];

  //  Then copy the string into memory.
  u32bit  numBases = 0;
  u32bit  numMers  = 0;

  char  x = FS->nextSymbol();
  while ((x) && (numBases < MAXBASES)) {
    if (validSymbol[x])
      SA[numBases++] = toUpper[x];
    else
      SA[numBases++] = 0;

    x = FS->nextSymbol();
  }
  SA[numBases] = 0;

  delete FS;

  fprintf(stderr, "Found "u32bitFMT" bases.\n", numBases);


#if 0
  //  Dump the sequence for meryl
  for (u32bit i=0; i<numBases; i++)
    if (SA[i] == 0)
      SA[i] = 'N';
  fprintf(stdout, ">testsequence\n");
  fwrite(SA, sizeof(char), numBases, stdout);
  fprintf(stdout, "\n");
  exit(0);
#endif



  fprintf(stderr, "Setup offsets.\n");

  //  Allocate offsets into FS, this is what we sort.
  //
  u32bit  *PS = new u32bit [numBases];

  //  Set the offsets
  //
  for (u32bit i=0; i<numBases; i++)
    if (SA[i])
      PS[numMers++] = i;

  fprintf(stderr, "Found "u32bitFMT" mers.\n", numMers);
  fprintf(stderr, "Sort.\n");

  qsort(PS, numMers-1, sizeof(u32bit), compare100);

  //  Now just dump!
  //
#if 0
  char mer[1024];
  for (u32bit i=0; i<numMers; i++) {
    strncpy(mer, SA + PS[i], 100);
    mer[100] = 0;
    fprintf(stdout, u32bitFMTW(6)" "u32bitFMTW(8)" %s\n", i, PS[i], mer);
  }
#endif

  return(0);
}

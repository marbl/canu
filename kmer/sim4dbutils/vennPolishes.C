#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "bio++.H"
#include "sim4.H"


const char *usage = 
"usage: %s [options] <polishes-file-1> <polishes-file-2> ...\n"
"\n"
"  Given n sets of sim4 polishes (of the same set of cDNA to the same\n"
"  set of genomic, but this isn't enforced) this code will generate a\n"
"  Venn diagram of how the sequences map.\n"
"\n"
"  -n <num-seqs>     there are <num-seqs> in the input set\n"
"  -i <min-ident>    filter matches to be >= <min-ident> identity\n"
"                    default = 95\n"
"  -c <min-cover>    filter matches to be >= <min-cover> coverage\n"
"                    default = 50\n"
"  -d <class-id>     dump the sequence IIDs in <class-id> to stdout\n"
"\n"
"  -plot             write a plot-able datafile of the venn diagram\n"
"                    for percent identity <min-idenit> to 100 (inclusive)\n"
"                    and <min-cover> coverage.\n";


//  Yes, yes.  Tell me all about how bad globals are.
u32bit   minI      = 95;
u32bit   minC      = 50;
u32bit   foundMax  = 100000;
u32bit   dumpIID   = ~u32bitZERO;
int      numArgs   = 0;
bool     plot      = false;
u32bit   numFiles = 0;
u32bit **found    = 0L;
u32bit   indexMax = 0;
u32bit  *counts   = 0L;
u32bit  *sizes    = 0L;


void
doVenn(u32bit minI, u32bit minC) {

  //  Count how many elements are in each set
  for (u32bit i=0; i<numFiles; i++) {
    sizes[i] = 0;
    for (u32bit j=0; j<foundMax; j++)
      if (found[i][j] >= minI)
        sizes[i]++;
  }


  for (u32bit i=0; i<indexMax; i++)
    counts[i] = 0;


  //  For each guy in the datasets
  //
  for (u32bit thisguy=0; thisguy < foundMax; thisguy++) {

    //  Compute which class he is in.  'class' is a reserved word,
    //  so we use 'membership' instead.
    //
    u32bit   membership = 0;

    for (u32bit dataset=0; dataset < numFiles; dataset++) {
      if (found[dataset][thisguy] >= minI)
        membership |= 1 << dataset;
    }

    if (membership == dumpIID)
      fprintf(stdout, u32bitFMT"\n", thisguy);

    counts[membership]++;
  }
}



int
main(int argc, char **argv) {

  if ((argc < 5)) {
    fprintf(stderr, usage, argv[0]);
    exit(1);
  }

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-n") == 0) {
      foundMax = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-i") == 0) {
      minI = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      minC = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-d") == 0) {
      dumpIID = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-plot") == 0) {
      plot = true;
    } else {
      //  Assume we got all the options, and we are at a file.
      //
      numArgs = arg;
      arg = argc;
    }

    arg++;
  }

  numFiles = argc - numArgs;
  found    = new u32bit * [numFiles];

  if (numFiles > 16) {
    fprintf(stderr, "WARNING:  You gave me "u32bitFMT" files!  That's pretty big.  I don't know\n", numFiles);
    fprintf(stderr, "          if I'm up to it.  Fasten seat belts and hang on!\n");
  }

  for (int arg=numArgs; arg<argc; arg++) {
    fprintf(stderr, "Reading '%s'\n", argv[arg]);

    found[arg-numArgs]    = new u32bit [foundMax];

    for (u32bit i=0; i<foundMax; i++)
      found[arg-numArgs][i] = 0;

    errno = 0;
    FILE *F = fopen(argv[arg], "r");
    if (errno) {
      fprintf(stderr, "Failed to open '%s': %s\n", argv[arg], strerror(errno));
      exit(1);
    }

    sim4polish *p = new sim4polish(F);
    while (p->_numExons > 0) {
      if ((p->_percentIdentity  >= minI) &&
          (p->_querySeqIdentity >= minC)) {

        if (p->_estID >= foundMax) {
          fprintf(stderr, "Please increase foundMax, or make me reallocate storage.\n");
          exit(1);
        }

        if (found[arg-numArgs][p->_estID] < p->_percentIdentity)
          found[arg-numArgs][p->_estID] = p->_percentIdentity;
      }

      delete p;
      p = new sim4polish(F);
    }
  }


  //  There are 2^n categories for n files.
  //
  //  If A and B, then there is
  //
  //  A B
  //  0 0 - neither (we can't compute this)
  //  0 1 - only B
  //  1 0 - only A
  //  1 1 - both A and B
  //
  //  So, we make an array of size 2^n that holds the ocunts of each
  //  class.  It's indexed by a bit vector.
  //
  indexMax = 1 << numFiles;
  counts   = new u32bit [indexMax];
  sizes    = new u32bit [numFiles];

  if        (dumpIID != ~u32bitZERO) {
    doVenn(minI, minC);
  } else if (plot) {
    for (u32bit id=minI; id <= 100; id++) {
      doVenn(id, minC);

      fprintf(stdout, u32bitFMTW(3)" ", id);
      for (u32bit i=0; i<numFiles; i++)
        fprintf(stdout, u32bitFMTW(8)" ", sizes[i]);
      for (u32bit index=0; index < indexMax; index++) {
        for (u32bit dataset=0; dataset < numFiles; dataset++)
          fprintf(stdout, "%c", (index & (1 << dataset)) ? 'A' + (char)dataset : '-');
        fprintf(stdout, " "u32bitFMTW(8)" ", counts[index]);
      }
      fprintf(stdout, "\n");
    }
  } else {
    doVenn(minI, minC);

    for (u32bit i=0; i<numFiles; i++)
      fprintf(stdout, "%c = ("u32bitFMTW(8)" total) %s\n", 'A' + (char)i, sizes[i], argv[i+numArgs]);

    for (u32bit index=0; index < indexMax; index++) {
      fprintf(stdout, u32bitFMTW(4)" [", index);
      for (u32bit dataset=0; dataset < numFiles; dataset++)
        fprintf(stdout, "%c", (index & (1 << dataset)) ? 'A' + (char)dataset : '-');
      fprintf(stdout, "]  "u32bitFMT"\n", counts[index]);
    }
  }
}

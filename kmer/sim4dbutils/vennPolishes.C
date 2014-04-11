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
uint32   minI      = 95;
uint32   minC      = 50;
uint32   foundMax  = 100000;
uint32   dumpIID   = ~uint32ZERO;
int      numArgs   = 0;
bool     plot      = false;
uint32   numFiles = 0;
uint32 **found    = 0L;
uint32   indexMax = 0;
uint32  *counts   = 0L;
uint32  *sizes    = 0L;


void
doVenn(uint32 minI, uint32 minC) {

  //  Count how many elements are in each set
  for (uint32 i=0; i<numFiles; i++) {
    sizes[i] = 0;
    for (uint32 j=0; j<foundMax; j++)
      if (found[i][j] >= minI)
        sizes[i]++;
  }


  for (uint32 i=0; i<indexMax; i++)
    counts[i] = 0;


  //  For each guy in the datasets
  //
  for (uint32 thisguy=0; thisguy < foundMax; thisguy++) {

    //  Compute which class he is in.  'class' is a reserved word,
    //  so we use 'membership' instead.
    //
    uint32   membership = 0;

    for (uint32 dataset=0; dataset < numFiles; dataset++) {
      if (found[dataset][thisguy] >= minI)
        membership |= 1 << dataset;
    }

    if (membership == dumpIID)
      fprintf(stdout, uint32FMT"\n", thisguy);

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
  found    = new uint32 * [numFiles];

  if (numFiles > 16) {
    fprintf(stderr, "WARNING:  You gave me "uint32FMT" files!  That's pretty big.  I don't know\n", numFiles);
    fprintf(stderr, "          if I'm up to it.  Fasten seat belts and hang on!\n");
  }

  for (int arg=numArgs; arg<argc; arg++) {
    fprintf(stderr, "Reading '%s'\n", argv[arg]);

    found[arg-numArgs]    = new uint32 [foundMax];

    for (uint32 i=0; i<foundMax; i++)
      found[arg-numArgs][i] = 0;

    sim4polishReader *R = new sim4polishReader("-");
    sim4polish       *p = 0L;

    while (R->nextAlignment(p)) {
      if ((p->_percentIdentity  >= minI) &&
          (p->_querySeqIdentity >= minC)) {

        if (p->_estID >= foundMax) {
          fprintf(stderr, "Please increase foundMax, or make me reallocate storage.\n");
          exit(1);
        }

        if (found[arg-numArgs][p->_estID] < p->_percentIdentity)
          found[arg-numArgs][p->_estID] = p->_percentIdentity;
      }
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
  counts   = new uint32 [indexMax];
  sizes    = new uint32 [numFiles];

  if        (dumpIID != ~uint32ZERO) {
    doVenn(minI, minC);
  } else if (plot) {
    for (uint32 id=minI; id <= 100; id++) {
      doVenn(id, minC);

      fprintf(stdout, uint32FMTW(3)" ", id);
      for (uint32 i=0; i<numFiles; i++)
        fprintf(stdout, uint32FMTW(8)" ", sizes[i]);
      for (uint32 index=0; index < indexMax; index++) {
        for (uint32 dataset=0; dataset < numFiles; dataset++)
          fprintf(stdout, "%c", (index & (1 << dataset)) ? 'A' + (char)dataset : '-');
        fprintf(stdout, " "uint32FMTW(8)" ", counts[index]);
      }
      fprintf(stdout, "\n");
    }
  } else {
    doVenn(minI, minC);

    for (uint32 i=0; i<numFiles; i++)
      fprintf(stdout, "%c = ("uint32FMTW(8)" total) %s\n", 'A' + (char)i, sizes[i], argv[i+numArgs]);

    for (uint32 index=0; index < indexMax; index++) {
      fprintf(stdout, uint32FMTW(4)" [", index);
      for (uint32 dataset=0; dataset < numFiles; dataset++)
        fprintf(stdout, "%c", (index & (1 << dataset)) ? 'A' + (char)dataset : '-');
      fprintf(stdout, "]  "uint32FMT"\n", counts[index]);
    }
  }
}

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "libbri.H"
#include "fasta.H"
#include "sim4polish.h"

//
//  usage: mergeInput -m match1 cdna1 -m match2 cdna2 -m ... -o match cdna
//
//  Merges the results from two ESTmapper runs.  The runs MUST be on
//  the same genomic sequence using DIFFERENT cDNA inputs.
//

int
main(int argc, char **argv) {
  char       **inMatchName = new char * [argc];
  char       **inSeqName   = new char * [argc];
  char        *otMatchName = 0L;
  char        *otSeqName   = 0L;

  FILE       **inMatch = new FILE * [argc];
  FILE        *otMatch = 0L;
  int         *numSeqs = new int [argc];

  sim4polish **polishes = new sim4polish * [argc];

  int          numIn = 0;

  int arg = 1;
  while (arg < argc) {
    if (strcmp(argv[arg], "-m") == 0) {
      arg++;
      inMatchName[numIn] = (char *)argv[arg++];
      inSeqName[numIn]   = (char *)argv[arg++];

      errno = 0;
      inMatch[numIn] = fopen(inMatchName[numIn], "r");
      if (inMatch[numIn] == 0L) {
        fprintf(stderr, "Couldn't open '%s'\n%s\n", inMatchName[numIn], strerror(errno));
        exit(1);
      }
      numIn++;
    } else if (strcmp(argv[arg], "-o") == 0) {
      arg++;
      otMatchName = (char *)argv[arg++];
      otSeqName   = (char *)argv[arg++];

      otMatch = fopen(otMatchName, "w");
    }
  }

  if ((numIn < 1) || (otMatch == 0L)) {
    fprintf(stderr, "usage: %s -o match cdna -m match1 cdna1 -m match2 cdna2 -m ...\n", argv[0]);
    exit(1);
  }


  //
  //  Merge the input sequences into the output sequence.  We also count
  //  the number of sequences here, so we don't need random-access
  //  of the input.
  //
  fprintf(stderr, "Merging sequences.\n");

  FILE         *O = fopen(otSeqName, "w");
  for (int i=0; i<numIn; i++) {
    FastAWrapper  *I = new FastAWrapper(inSeqName[i]);

    numSeqs[i] = 0;

    while (!I->eof()) {
      FastASequenceInCore *B = I->getSequence();

      fprintf(O, "%s\n%s\n", B->header(), B->sequence());
      numSeqs[i]++;

      delete B;
    }

    delete I;
  }
  fclose(O);


  //  Make numSeqs[] be the offset needed to convert a polish
  //  in each inMatch[] file into a polish in the merged file.
  //
  int o = 0;
  int s = 0;
  for (int i=0; i<numIn; i++) {
    o  = numSeqs[i];
    numSeqs[i] = s;
    s += o;
  }



  //
  //  Load the initial polishes
  //
  fprintf(stderr, "Loading initial.\n");

  for (int i=0; i<numIn; i++) {
    polishes[i]           = s4p_readPolish(inMatch[i]);
    if (polishes[i])
      polishes[i]->estID += numSeqs[i];

    if (polishes[i] == 0L)
      fprintf(stderr, "Error: no matches found in %s\n", inMatchName[i]);
  }


  //
  //  Do a merge, until no more input exists
  //
  fprintf(stderr, "Merging polishes.\n");

  bool keepGoing = true;
  while (keepGoing) {

    //  Find the lowest polish -- we want to keep these sorted by scaffold,
    //  then by cDNA.
    //
    int first = 0;
    while ((polishes[first] == 0L) && (first < numIn))
      first++;

    if (polishes[first]) {
      for (int i=first+1; i<numIn; i++)
        if ((polishes[i]) &&
            (s4p_genIDcompare(polishes + first, polishes + i) > 0))
          first = i;


      //  Dump the 'first', read in a new polish
      //
      s4p_printPolish(otMatch, polishes[first], S4P_PRINTPOLISH_FULL);

      s4p_destroyPolish(polishes[first]);

      polishes[first]           = s4p_readPolish(inMatch[first]);
      if (polishes[first])
        polishes[first]->estID += numSeqs[first];
    } else {
      keepGoing = false;
    }
  }

  fclose(otMatch);
}

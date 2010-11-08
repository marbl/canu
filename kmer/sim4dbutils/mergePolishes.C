#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "bio++.H"
//#include "fasta.H"
#include "sim4.H"

//  usage: mergeInput -m match1 cdna1 -m match2 cdna2 -m ... -o match cdna [-gff3]
//
//  Merges the results from two ESTmapper runs.  The runs MUST be on
//  the same genomic sequence using DIFFERENT cDNA inputs.

static
void
loadNext(u32bit idx, sim4polish **polishes, sim4polishReader **inMatch, u32bit *numSeqs) {
  if (inMatch[idx]->nextAlignment(polishes[idx]))
    polishes[idx]->_estID += numSeqs[idx];
}

int
main(int argc, char **argv) {
  char              **inMatchName = new char * [argc];
  char              **inSeqName   = new char * [argc];
  char               *otMatchName = 0L;
  char               *otSeqName   = 0L;

  sim4polishReader  **inMatch  = new sim4polishReader * [argc];
  sim4polish        **polishes = new sim4polish * [argc];
  sim4polishWriter   *otMatch  = 0L;
  u32bit             *numSeqs  = new u32bit [argc];

  u32bit              numIn = 0;

  sim4polishStyle    style = sim4polishStyleDefault;

  int arg = 1;
  while (arg < argc) {
    if (strcmp(argv[arg], "-m") == 0) {
      arg++;
      inMatchName[numIn] = (char *)argv[arg++];
      inSeqName[numIn]   = (char *)argv[arg++];

      inMatch[numIn] = new sim4polishReader(inMatchName[numIn]);
      numIn++;

    } else if (strcmp(argv[arg], "-o") == 0) {
      arg++;
      otMatchName = (char *)argv[arg++];
      otSeqName   = (char *)argv[arg++];

    } else if (strcmp(argv[arg], "-gff3") == 0) {
      style = sim4polishGFF3;

    }
  }

  if ((numIn < 1) || (otMatch == 0L)) {
    fprintf(stderr, "usage: %s -o match cdna -m match1 cdna1 -m match2 cdna2 -m ... [-gff3]\n", argv[0]);
    exit(1);
  }

  otMatch = new sim4polishWriter(otMatchName, style);

  for (u32bit i=0; i<numIn; i++)
    if (inMatch[i]->getsim4polishStyle() != style) {
      fprintf(stderr, "warning: input format and output format may differ.\n");
      break;
    }

  //  Merge the input sequences into the output sequence.  We also count the number of sequences
  //  here, so we don't need random-access of the input.
  //
  fprintf(stderr, "Merging sequences.\n");

  FILE         *O = fopen(otSeqName, "w");
  for (u32bit i=0; i<numIn; i++) {
    seqCache  *I = new seqCache(inSeqName[i]);
    seqInCore *B = I->getSequenceInCore();

    numSeqs[i] = 0;

    while (B) {
      fprintf(O, ">%s\n%s\n", B->header(), B->sequence());
      numSeqs[i]++;

      delete B;
      B = I->getSequenceInCore();
    }

    delete I;
  }
  fclose(O);


  //  Make numSeqs[] be the offset needed to convert a polish in each inMatch[] file into a polish
  //  in the merged file.
  //
  u32bit o = 0;
  u32bit s = 0;
  for (u32bit i=0; i<numIn; i++) {
    o  = numSeqs[i];
    numSeqs[i] = s;
    s += o;
  }

  //  Load the initial polishes
  //
  for (u32bit i=0; i<numIn; i++)
    loadNext(i, polishes, inMatch, numSeqs);

  //  Merge, until no more input is left.  Each round we scan the list of loaded polishes[] and
  //  remember the lowest, which is then output and a new polish is loaded in its place.
  //
  bool keepGoing = true;
  while (keepGoing) {

    u32bit first = 0;
    while ((polishes[first] == 0L) && (first < numIn))
      first++;

    if (polishes[first] == 0L) {
      keepGoing = 0L;
      continue;
    }

    for (u32bit i=first+1; i<numIn; i++)
      if ((polishes[i]) &&
          (s4p_genIDcompare(polishes + first, polishes + i) > 0))
        first = i;

    otMatch->writeAlignment(polishes[first]);

    loadNext(first, polishes, inMatch, numSeqs);
  }

  delete [] polishes;

  delete inMatch;
  delete otMatch;
}




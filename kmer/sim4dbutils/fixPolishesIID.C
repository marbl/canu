#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "bio.h"
#include "sim4.H"

#include <map>
#include <string>

using namespace std;

//  Updates the IID's in a set of polishes.  If a file of deflines (or
//  fasta file) is supplied, the IIDs will match those, otherwise,
//  they remain the same.

void
addToDict(map<string, u64bit> &d, char *n) {

  if (n == 0L)
    return;

  seqCache  *F = new seqCache(n);
  seqInCore *S = F->getSequenceInCore();

  while (S) {
    string   s = S->header();

    d[s] = S->getIID();

    delete S;
    S = F->getSequenceInCore();
  }

  delete F;
}



int
main(int argc, char **argv) {
  char  *cDeflines = 0L;
  char  *gDeflines = 0L;

  sim4polishStyle style = sim4polishStyleDefault;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-c") == 0) {
      cDeflines = argv[++arg];

    } else if (strcmp(argv[arg], "-g") == 0) {
      gDeflines = argv[++arg];

    } else if (strcmp(argv[arg], "-gff3") == 0) {
      style = sim4polishGFF3;

    } else {
      fprintf(stderr, "Unknown arg: %s\n", argv[arg]);
    }
    arg++;
  }
  if (isatty(fileno(stdin))) {
    fprintf(stderr, "usage: %s [-c c.fasta] [-g g.fasta] [-gff3] < polishes > polishes\n", argv[0]);
    fprintf(stderr, "     -c c.fasta   Read cDNA deflines from c.fasta\n");
    fprintf(stderr, "     -g g.fasta   Read genomic deflines from g.fasta\n");
    fprintf(stderr, "     -gff3        Write output as GFF3\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Rewrites the input polishes, updating the sequence index to match\n");
    fprintf(stderr, "  that of the associated fasta file.  One or both of -c and -g may be used.\n");
    fprintf(stderr, "  Polishes that refer to a sequence not present in the input fasta file are\n");
    fprintf(stderr, "  not output.\n");
    exit(1);
  }

  //  We parse args, then build the dictionaries, so we can do
  //  any quick error detection first.

  map<string, u64bit>   g;
  map<string, u64bit>   c;

  if (gDeflines) {
    fprintf(stderr, "Reading genomic deflines from '%s'\n", gDeflines);
    addToDict(g, gDeflines);
  }

  if (cDeflines) {
    fprintf(stderr, "Reading genomic deflines from '%s'\n", cDeflines);
    addToDict(c, cDeflines);
  }

  //  Read all the matches, changing IIDs.  If we find a defline
  //  with no IID, holler and die.

  sim4polishWriter *W = new sim4polishWriter("-", style);
  sim4polishReader *R = new sim4polishReader("-");
  sim4polish       *p = 0L;

  if (R->getsim4polishStyle() != style) 
    fprintf(stderr, "warning: input format and output format differ.\n");

  fprintf(stderr, "Filtering polishes.\n");

  while (R->nextAlignment(p)) {
    string cd = p->_estDefLine;
    string gd = p->_genDefLine;

    if (cDeflines != 0L) {
      if (c.find(cd) == c.end())
        //  EST defline not in the input sequences, don't output.
        continue;
      p->_estID = c[cd];
    }

    if (gDeflines != 0L) {
      if (g.find(gd) == g.end())
        //  Genomic defline not in the input sequences, don't output.
        continue;
      p->_genID = g[gd];
    }

    W->writeAlignment(p);
  }

  delete R;
  delete W;
}

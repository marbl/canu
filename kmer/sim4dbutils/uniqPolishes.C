#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "sim4.H"


void
pickBest(sim4polishWriter *W, sim4polish **p, int pNum, int uniq) {
  int i;

  if (pNum == 1) {
    if (uniq)
      W->writeAlignment(p[0]);
  } else {
    if (!uniq)
      for (i=0; i<pNum; i++)
        W->writeAlignment(p[0]);
  }

  for (i=0; i<pNum; i++)
    delete p[i];
}




int
main(int argc, char **argv) {
  u32bit       pNum   = 0;
  u32bit       pAlloc = 8388608;
  u32bit       estID  = ~u32bitZERO;

  u32bit       uniq   = 1;

  sim4polishStyle style = sim4polishStyleDefault;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-uniq", 2) == 0) {
      uniq = 1;
    } else if (strncmp(argv[arg], "-dupl", 2) == 0) {
      uniq = 0;
    } else if (strcmp(argv[arg], "-gff3") == 0) {
      style = sim4polishGFF3;
    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }

  if (isatty(fileno(stdin))) {
    fprintf(stderr, "usage: %s [-uniq | -dupl] [-gff3] < file > file\n", argv[0]);

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");

    exit(1);
  }

  //  Read polishes, picking the best when we see a change in
  //  the estID.

  sim4polishWriter  *W = new sim4polishWriter("-", style);
  sim4polishReader  *R = new sim4polishReader("-");
  sim4polish       **p = new sim4polish * [pAlloc];
  sim4polish        *q = 0L;

  if (R->getsim4polishStyle() != style) 
    fprintf(stderr, "warning: input format and output format differ.\n");

  while (R->nextAlignment(q)) {
    if ((q->_estID != estID) && (pNum > 0)) {
      pickBest(W, p, pNum, uniq);
      pNum  = 0;
    }

    if (pNum >= pAlloc) {
      sim4polish **P = new sim4polish * [pAlloc * 2];
      memcpy(p, P, sizeof(sim4polish *) * pAlloc);
      delete [] p;
      p = P;
      pAlloc *= 2;
    }

    p[pNum++] = q;
    estID     = q->_estID;

    q = 0L;  //  Else we'll delete the polish we just saved!
  }

  if (pNum > 0)
    pickBest(W, p, pNum, uniq);

  delete [] p;
  delete    R;
  delete    W;

  return(0);
}


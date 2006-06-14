#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "bio.h"
#include "sim4.H"

//  Reports the amount of sequence covered by ALL matches for that
//  sequence.  Example: if sequence iid 4 has two matches, one
//  covering the first 30% and the second covering the last 30%, this
//  will report that sequence iid 4 is covered 60%.
//
//  Takes no options, reads from stdin, writes to stdout.

int
main(int argc, char **argv) {
  sim4polish     *p;
  u32bit          covMax = 8 * 1024 * 1024;
  intervalList  **cov    = new intervalList * [covMax];
  u32bit         *len    = new u32bit [covMax];


  for (u32bit i=0; i<covMax; i++)
    cov[i] = 0L;


  if (isatty(fileno(stdin))) {
    fprintf(stderr, "usage: %s < file > file\n", argv[0]);

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");

    exit(1);
  }



  while ((p = s4p_readPolish(stdin)) != 0L) {
    if (p->estID > covMax)
      fprintf(stderr, "DIE!  I should make my space bigger.\n"), exit(1);

    if (cov[p->estID] == 0L) {
      cov[p->estID] = new intervalList;
      len[p->estID] = p->estLen;
    }

    for (u32bit e=0; e<p->numExons; e++)
      cov[p->estID]->add(p->exons[e].estFrom, p->exons[e].estTo - p->exons[e].estFrom + 1);

    s4p_destroyPolish(p);
  }


  //  Scan the list of intervalLists, compute the amount covered, print.
  //
  for (u32bit i=0; i<covMax; i++) {
    if (cov[i]) {
      u32bit  numRegions = cov[i]->numberOfIntervals();
      u32bit  sumLengths = cov[i]->sumOfLengths();

      cov[i]->merge();

      fprintf(stdout, u32bitFMT"\t"u32bitFMT"\t%5.3f\t"u32bitFMT"\t"u32bitFMT"\n",
              i,
              len[i],
              cov[i]->sumOfLengths() / (double)len[i],
              numRegions,
              sumLengths);
    }
  }
}

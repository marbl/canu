#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "sim4reader.h"

//
//  Outputs some statistics on the matches
//

#define HISTBIN (1000)
#define HISTMAX (300000000 / HISTBIN)

int
main(int argc, char ** argv) {
  u32bit       dumpSize = 0;
  u32bit      *hist;
  FILE        *all;
  FILE        *big;
  int          i, j;

  if (isatty(fileno(stdin))) {
    fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");
    exit(1);
  }

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-dump", 2) == 0) {
      dumpSize = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-all", 2) == 0) {
      all = fopen(argv[++arg], "w");
      if (all == 0L) {
        fprintf(stderr, "Can't open '%s' for writing\n", argv[arg]);
        exit(1);
      }
    } else if (strncmp(argv[arg], "-big", 2) == 0) {
      big = fopen(argv[++arg], "w");
      if (big == 0L) {
        fprintf(stderr, "Can't open '%s' for writing\n", argv[arg]);
        exit(1);
      }
    } else {
      fprintf(stderr, "Unknown option: '%s'\n", argv[arg]);
    }

    arg++;
  }

  if (all || big) {
    hist = new u32bit [HISTMAX];
    memset(hist, 0, sizeof(u32bit) * HISTMAX);
  }

  sim4polish  *p = new sim4polish(stdin);
  while (p->_numExons > 0) {
    if (p->numExons > 1) {
      int   exA;
      int   exB;
      int   biggestIntron = 0;

      for (exA=0, exB=1; exB < p->numExons; exA++, exB++) {
        int dist = p->exons[exB].genFrom - p->exons[exA].genTo + 1;
        if (dist > biggestIntron)
          biggestIntron = dist;
        if (all)
          hist[dist / HISTBIN]++;
      }

      if (big)
        hist[biggestIntron / HISTBIN]++;

      //fprintf(stdout, "%d\n", biggestIntron);

      if ((dumpSize > 0) && (biggestIntron > dumpSize))
        printPolish(stdout, p);
    }

    destroyPolish(p);
  }


  if (all) {
    for (j=HISTMAX-1; hist[j]==0 && j>=0; j--)
      ;
    for (i=0; i<j; i++)
      fprintf(all, "%d\n", hist[i]);
    fclose(all);
  }

  if (big) {
    for (j=HISTMAX-1; hist[j]==0 && j>=0; j--)
      ;
    for (i=0; i<j; i++)
      fprintf(big, "%d\n", hist[i]);
    fclose(big);
  }

  delete [] hist;

  return(0);
}

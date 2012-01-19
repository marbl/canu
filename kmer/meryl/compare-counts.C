#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libmeryl.H"


int
main(int argc, char **argv) {
  merylStreamReader  *A = 0L;
  merylStreamReader  *B = 0L;
  bool                beVerbose = true;

  int arg = 1;
  int err = 0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-a") == 0) {
      A = new merylStreamReader(argv[++arg]);

    } else if (strcmp(argv[arg], "-b") == 0) {
      B = new merylStreamReader(argv[++arg]);

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((A == 0L) || (B == 0L) || (err)) {
    fprintf(stderr, "usage: %s [opts] [-v] [-dump prefix]\n", argv[0]);
    fprintf(stderr, "At least one fragcounts and one contigcounts are needed.\n");
    fprintf(stderr, "          -af | -tf        fragcounts\n");
    fprintf(stderr, "          -ac | -dc | -co  contigcounts \n");
    fprintf(stderr, "Dumping is probably only useful with exactly one frag and\n");
    fprintf(stderr, "one contig, but I'll let you do it with any number.\n");
    exit(1);
  }

  speedCounter *C = new speedCounter(" Examining: %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1ffffff, beVerbose);

#define MAXA 150
#define MAXB 150

  double   heatraw[MAXA][MAXB];
  double   heatsca[MAXA][MAXB];

  for (u32bit i=0; i<MAXA; i++)
    for (u32bit j=0; j<MAXB; j++)
      heatraw[i][j] = heatsca[i][j] = 0;

  A->nextMer();
  B->nextMer();

  while ((A->validMer()) ||
         (B->validMer())) {
    kMer  &a = A->theFMer();
    kMer  &b = B->theFMer();

    u32bit ac = A->theCount();
    u32bit bc = B->theCount();

    if (ac >= MAXA)
      ac = MAXA-1;

    if (bc >= MAXB)
      bc = MAXB-1;

    if (A->validMer() == false) {
      ac = 0;
      heatraw[ac][bc]++;
      B->nextMer();
      continue;
    }

    if (B->validMer() == false) {
      bc = 0;
      heatraw[ac][bc]++;
      A->nextMer();
      continue;
    }

    if (a == b) {
      heatraw[ac][bc]++;
      A->nextMer();
      B->nextMer();

    } else if (a < b) {
      heatraw[ac][0]++;
      A->nextMer();

    } else {
      heatraw[0][bc]++;
      B->nextMer();
    }

    C->tick();
  }

  delete C;
  delete A;
  delete B;

  //  Scale each row to be between 0 and 1

#if 0
  for (u32bit j=0; j<MAXB; j++) {
    double  mina = heatraw[0][j];
    double  maxa = heatraw[0][j];

    for (u32bit ii=0; ii<MAXA; ii++) {
      if (maxa < heatraw[ii][j])
        maxa = heatraw[ii][j];
      if (heatraw[ii][j] < mina)
        mina = heatraw[ii][j];
    }

    for (u32bit i=0; i<MAXA; i++)
      heatsca[i][j] = (heatraw[i][j] - mina) / (maxa - mina);
  }
#endif


  for (u32bit i=0; i<MAXA; i++)
    for (u32bit j=0; j<MAXB; j++)
      fprintf(stdout, u32bitFMT"\t"u32bitFMT"\t%f\n", i, j, log(heatraw[i][j]));
}

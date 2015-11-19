
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    kmer/meryl/compare-counts.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-AUG-28 to 2014-APR-11
 *      are Copyright 2010,2012,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libmeryl.H"




#if 0
void
heatMap() {
  speedCounter *C = new speedCounter(" Examining: %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1ffffff, false);

#define MAXA 150
#define MAXB 150

  double   heatraw[MAXA][MAXB];
  double   heatsca[MAXA][MAXB];

  for (uint32 i=0; i<MAXA; i++)
    for (uint32 j=0; j<MAXB; j++)
      heatraw[i][j] = heatsca[i][j] = 0;

  A->nextMer();
  B->nextMer();

  while ((A->validMer()) ||
         (B->validMer())) {
    kMer  &a = A->theFMer();
    kMer  &b = B->theFMer();

    uint32 ac = A->theCount();
    uint32 bc = B->theCount();

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
  for (uint32 j=0; j<MAXB; j++) {
    double  mina = heatraw[0][j];
    double  maxa = heatraw[0][j];

    for (uint32 ii=0; ii<MAXA; ii++) {
      if (maxa < heatraw[ii][j])
        maxa = heatraw[ii][j];
      if (heatraw[ii][j] < mina)
        mina = heatraw[ii][j];
    }

    for (uint32 i=0; i<MAXA; i++)
      heatsca[i][j] = (heatraw[i][j] - mina) / (maxa - mina);
  }
#endif


  for (uint32 i=0; i<MAXA; i++)
    for (uint32 j=0; j<MAXB; j++)
      fprintf(stdout, uint32FMT"\t"uint32FMT"\t%f\n", i, j, log(heatraw[i][j]));
}
#endif





int
main(int argc, char **argv) {
  merylStreamReader  *T = 0L;
  merylStreamReader  *S = 0L;
  char               *outputPrefix = NULL;
  char               *plotTitle    = NULL;
  int arg = 1;
  int err = 0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-truth") == 0) {
      T = new merylStreamReader(argv[++arg]);

    } else if (strcmp(argv[arg], "-sample") == 0) {
      S = new merylStreamReader(argv[++arg]);

    } else if (strcmp(argv[arg], "-output") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-title") == 0) {
      plotTitle = argv[++arg];

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((T == 0L) || (S == 0L) || (outputPrefix == 0L) || (plotTitle == 0L) || (err)) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "  -truth     k-mers from reference\n");
    fprintf(stderr, "  -sample    k-mers from sample\n");
    fprintf(stderr, "  -output    output prefix\n");
    fprintf(stderr, "  -title     plot label\n");
    exit(1);
  }

  uint32   kmerSize = T->merSize();

#define HMAX  64 * 1024

  uint32  *Htrue  = new uint32 [HMAX];
  uint32  *Hnoise = new uint32 [HMAX];

  for (uint32 i=0; i<HMAX; i++)
    Htrue[i] = Hnoise[i] = 0;

  T->nextMer();
  S->nextMer();

  while ((T->validMer()) ||
         (S->validMer())) {
    kMer  &t = T->theFMer();
    kMer  &s = S->theFMer();

    uint32 tc = T->theCount();
    uint32 sc = S->theCount();

    if (tc >= HMAX)   tc = HMAX-1;
    if (sc >= HMAX)   sc = HMAX-1;

    //  If we're out of truth kmers, the sample is noise.
    if (T->validMer() == false) {
      Hnoise[sc]++;
      S->nextMer();
      continue;
    }

    //  If we're out of sample kmers, do nothing but go to the next truth kmer.
    if (S->validMer() == false) {
      T->nextMer();
      continue;
    }

    //  If the kmers are equal, this is a true kmer
    if (t == s) {
      Htrue[sc]++;
      T->nextMer();
      S->nextMer();
    }

    //  If the truth kmer is the lesser, get the next truth.
    else if (t < s) {
      T->nextMer();
    }

    //  Else the sample kmer is smaller, add it to the noise pile, and get the next.
    else {
      Hnoise[sc]++;
      S->nextMer();
    }
  }

  delete T;
  delete S;

  char  outputName[FILENAME_MAX];

  sprintf(outputName, "%s.gp", outputPrefix);
  FILE *outputGP  = fopen(outputName, "w");

  sprintf(outputName, "%s.dat", outputPrefix);
  FILE *outputDAT = fopen(outputName, "w");

  fprintf(outputGP, "set terminal png\n");
  fprintf(outputGP, "set output \"%s.png\"\n", outputPrefix);
  fprintf(outputGP, "set title \"%s true/false %d-mers\"\n", plotTitle, kmerSize);
  fprintf(outputGP, "set xlabel \"k-mer count\"\n");
  fprintf(outputGP, "set ylabel \"number of kmers\"\n");
  fprintf(outputGP, "plot [0:100] [0:1000000] \"%s.dat\" using 1:2 with lines title \"true\", \"%s.dat\" using 1:3 with lines title \"false\"\n",
          outputPrefix, outputPrefix);

  fclose(outputGP);

  for (uint32 i=0; i<HMAX; i++)
    fprintf(outputDAT, uint32FMT"\t"uint32FMT"\t"uint32FMT"\n", i, Htrue[i], Hnoise[i]);

  fclose(outputDAT);

  sprintf(outputName, "gnuplot < %s.gp", outputPrefix);
  system(outputName);

  exit(0);
}

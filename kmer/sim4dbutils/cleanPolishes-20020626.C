#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "sim4reader.h"

#define SHOWTRIMMING

char const *usage =
"usage: %s [-save splitFile] [-threshold t]\n"
"  -threshold    Introns bigger than this are split into two matches (default = 150000).\n"
"  -savesplits   Saves a before/after of each split match.\n"
"                  All matches are printed to stdout (untrimmed and trimmed).\n"
"\n";



bool
lowComplexityExon(char *s) {
  int    cnt[5][5] = {0};
  int    map[256]  = {0};
  int    i, j, len = 0;
  int    a=0, b=0, c=0;
  double qual = 0.0;

  if (s == 0L)
    return(false);

  map['A'] = map['a'] = 1;
  map['C'] = map['c'] = 2;
  map['G'] = map['g'] = 3;
  map['T'] = map['t'] = 4;

  for (i=0; i<5; i++)
    for (j=0; j<5; j++)
      cnt[i][j] = 0;

  for (i=0, j=1; s[j]; i++, j++) {
    cnt[map[s[i]]][map[s[j]]]++;
    len++;
  }

  for (i=0; i<5; i++) {
    for (j=0; j<5; j++) {
      if (a < cnt[i][j]) {
        c = b;
        b = a;
        a = cnt[i][j];
      } else if (b < cnt[i][j]) {
        c = b;
        b = cnt[i][j];
      } else if (c < cnt[i][j]) {
        c = cnt[i][j];
      }
    }
  }

  qual = (double)(a+b+c) / (double)(len);

  if (len > 50)
    qual = 0.0;

  //if (qual > 0.75)
  //fprintf(stdout, "%8.5f:\t%s\n", qual, s);

  return(qual > 0.75);
}




int
main(int argc, char ** argv) {
  int          arg  = 1;
  FILE        *splitFile = 0L;
  int          intronLimit = 150000;
  sim4polish  *p;

#if 0
  if (isatty(fileno(stdin)) || isatty(fileno(stdout))) {
    fprintf(stderr, usage, argv[0]);

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");

    if (isatty(fileno(stdout)))
      fprintf(stderr, "error: Please redirect the polishes to a file.\n      (They are on stdout)\n\n");

    exit(1);
  }
#endif

  arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-savesplits", 2) == 0) {
      arg++;
      errno=0;
      splitFile = fopen(argv[arg], "w");
      if (errno) {
        fprintf(stderr, "Can't open '%s' for writing\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }
    } else if (strncmp(argv[arg], "-threshold", 2) == 0) {
      intronLimit = atoi(argv[++arg]);
    }

    arg++;
  }


  //  Statistics on the splitting quality / frequency
  int  totMatches = 0;
  int  oneExon    = 0;
  int  smaIntron  = 0;
  int  junkFirst  = 0;
  int  junkLast   = 0;
  int  junkBoth   = 0;
  int  splitOnGap = 0;
  int  goodQual   = 0;
  int  flanking   = 0;


  FILE *junkF   = fopen("spl.junkfirst", "w");
  FILE *junkL   = fopen("spl.junklast", "w");
  FILE *junkB   = fopen("spl.junkboth", "w");
  FILE *splGap  = fopen("spl.splitGap", "w");
  FILE *good    = fopen("spl.good", "w");
  FILE *flank   = fopen("spl.flanking", "w");

  

  while ((p = readPolish(stdin)) != 0L) {
    int exA;
    int exB;


    if (p->numExons == 1) {
      oneExon++;
    } else {

      //  Find the big intron.  We assume there is only one big intron.
      //
      int biggestIntron = 0;
      int intronSplit   = 0;
      int intronOri     = 0;

      for (exA=0, exB=1; exB < p->numExons; exA++, exB++) {
        int dist = p->exons[exB].genFrom - p->exons[exA].genTo + 1;
        if (dist > biggestIntron) {
          biggestIntron = dist;
          intronSplit   = exB;
          intronOri     = p->exons[exA].intronOrientation;
        }
      }

      if (intronOri == 0) {
        fprintf(stderr, "didn't find the largest intron? (got zero)?\n");
        exit(1);
      }

      if (intronOri == INTRON_NONE) {
        fprintf(stderr, "biggest intron isn't an intron? (got none)?\n");
        exit(1);
      }

      if (biggestIntron < 100000) {
        smaIntron++;
      } else {


        //  Declare the split obvious if all exons on either side are
        //  below 30bp, difficult otherwise.
        //
        bool  killFirst = true;
        bool  killLast  = true;

        for (int i=0; i<intronSplit; i++)
          if ((p->exons[i].estTo - p->exons[i].estFrom + 1 >= 50) &&
              (p->exons[i].percentIdentity >= 88) &&
              (lowComplexityExon(p->exons[i].estAlignment) == false))
            killFirst = false;

        for (int i=intronSplit; i<p->numExons; i++)
          if ((p->exons[i].estTo - p->exons[i].estFrom + 1 >= 50) &&
              (p->exons[i].percentIdentity >= 88) &&
              (lowComplexityExon(p->exons[i].estAlignment) == false))
            killLast = false;


        //  We shouldn't ever want to kill both sides.
        //
        if ((killFirst == true) && (killLast == true)) {
          junkBoth++;
          fprintf(junkB, "==============================JUNK FIRST AND LAST?\n");
          printPolish(junkB, p);
        }

        if ((killFirst == true) && (killLast == false)) {
          junkFirst++;
          printPolish(junkF, p);
          fprintf(junkF, "==============================\n");
        }

        if ((killFirst == false) && (killLast == true)) {
          junkLast++;
          printPolish(junkL, p);
          fprintf(junkL, "==============================\n");
        }

        if ((killFirst == false) && (killLast == false)) {
          if (intronOri == INTRON_GAP) {
            splitOnGap++;
            printPolish(splGap, p);
            fprintf(splGap, "==============================\n");
          } else {

            //  If there is a valid strand prediction and
            //      a) all exons >= 90%
            //      b) all exons >= 95%
            //      c) all exons >= 95%, except first and last, which can be >= 90%
            //  save the match as is.
            //
            bool  validStrand = false;
            if ((p->strandOrientation == STRAND_POSITIVE) ||
                (p->strandOrientation == STRAND_NEGATIVE))
              validStrand = true;

#if 0
            bool  qualIsA = true;
            for (exA=0; exA < p->numExons; exA++)
              if (p->exons[exA].percentIdentity < 90)
                qualIsA = false;

            bool  qualIsB = true;
            for (exA=0; exA < p->numExons; exA++)
              if (p->exons[exA].percentIdentity < 95)
                qualIsB = false;
#endif

            bool  qualIsC = true;
            if (p->exons[0].percentIdentity < 90)
              qualIsC = false;
            if (p->exons[p->numExons-1].percentIdentity < 90)
              qualIsC = false;
            for (exA=1; exA < p->numExons-1; exA++)
              if (p->exons[exA].percentIdentity < 95)
                qualIsC = false;

            //  If the match looks good, but just has a large intron, keep it.
            //
            if (validStrand && qualIsC) {
              printPolish(good, p);
              fprintf(good, "==============================\n");
              goodQual++;
            } else {
              flanking++;
              printPolish(flank, p);
              fprintf(flank, "==============================\n");
            }
          }
        }

      }  //  Has a big intron
    }  //  More than one exon

    totMatches++;
    if ((totMatches % 3759) == 0) {
      fprintf(stderr, "tot: %7d ", totMatches);
      fprintf(stderr, "one: %7d ", oneExon);
      fprintf(stderr, "sma: %7d ", smaIntron);
      fprintf(stderr, "jnkF: %7d ", junkFirst);
      fprintf(stderr, "jnkL: %7d ", junkLast);
      fprintf(stderr, "jnkB: %7d ", junkBoth);
      fprintf(stderr, "onGap: %7d ", splitOnGap);
      fprintf(stderr, "good: %7d ", goodQual);
      fprintf(stderr, "flank: %7d\r", flanking);
    }

    destroyPolish(p);
  }

  fclose(junkF);
  fclose(junkL);
  fclose(junkB);
  fclose(splGap);
  fclose(good);
  fclose(flank);

  fprintf(stderr, "tot: %7d ", totMatches);
  fprintf(stderr, "one: %7d ", oneExon);
  fprintf(stderr, "sma: %7d ", smaIntron);
  fprintf(stderr, "jnkF: %7d ", junkFirst);
  fprintf(stderr, "jnkL: %7d ", junkLast);
  fprintf(stderr, "jnkB: %7d ", junkBoth);
  fprintf(stderr, "onGap: %7d ", splitOnGap);
  fprintf(stderr, "good: %7d ", goodQual);
  fprintf(stderr, "flank: %7d\n", flanking);

  return(0);
}

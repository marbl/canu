#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "sim4reader.h"

#define SHOWTRIMMING

char const *usage =
"usage: %s [-save trimmedFile]\n"
"  -savetrimming   Saves a before/after of each trimmed match.\n"
"                  All matches are printed to stdout (untrimmed and trimmed).\n"
"\n";

int
main(int argc, char ** argv) {
  int          arg  = 1;
  FILE        *trimmedFile = 0L;
  int          beVerbose = 0;
  sim4polish  *p;
  int          polishesProcessed = 0;
  int          polishesTrimmed   = 0;

  if (isatty(fileno(stdin)) || isatty(fileno(stdout))) {
    fprintf(stderr, usage, argv[0]);

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");

    if (isatty(fileno(stdout)))
      fprintf(stderr, "error: Please redirect the polishes to a file.\n      (They are on stdout)\n\n");

    exit(1);
  }

  arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-savetrimming", 2) == 0) {
      arg++;
      errno=0;
      trimmedFile = fopen(argv[arg], "w");
      if (errno) {
        fprintf(stderr, "Can't open '%s' for writing\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }
    } else if (strncmp(argv[arg], "-verbose", 2) == 0) {
      beVerbose = 1;
    }

    arg++;
  }

  while ((p = readPolish(stdin)) != 0L) {
    int trimFirst = 0;
    int trimLast  = 0;

    /*  Decide if we need to trim anything
     */
    if (p->numExons > 1) {
      int   exA;
      int   exB;
      int   dist;
      int   qual;
      int   size;

      exA = 0;  //  First exon
      exB = 1;  //  Second exon
      dist = p->exons[exB].genFrom - p->exons[exA].genTo + 1;
      qual = p->exons[exA].percentIdentity;
      size = p->exons[exA].estTo - p->exons[exA].estFrom + 1;

      trimFirst = 1;

      if (dist < 100000)
        trimFirst = 0;

      if (size >= 50)
        trimFirst = 0;

      if (size >= 25 + (int)((dist - 100000) * 25.0 / 900000.0))
        trimFirst = 0;

      if ((qual >= 98) &&
          (size >= 25 + (int)((dist - 100000) * 25.0 / 1400000.0)))
        trimFirst = 0;

      //  Reverse our decision if the first exon is of low quality.
      //
      if ((qual < 85) && (dist >= 10000)) {
        if (trimFirst == 0)
          fprintf(trimmedFile, "Trimming frist exon based only on percent ID\n");
        trimFirst = 1;
      }

      exA = p->numExons - 1;  //  Last exon
      exB = p->numExons - 2;  //  Second to last
      dist = p->exons[exA].genFrom - p->exons[exB].genTo + 1;
      qual = p->exons[exA].percentIdentity;
      size = p->exons[exA].estTo - p->exons[exA].estFrom + 1;

      trimLast  = 1;

      if (dist < 100000)
        trimLast = 0;

      if (size >= 50)
        trimLast = 0;

      if (size >= 25 + (int)((dist - 100000) * 25.0 / 900000.0))
        trimLast = 0;

      if ((qual >= 98) &&
          (size >= 25 + (int)((dist - 100000) * 25.0 / 1400000.0)))
        trimLast = 0;

      //  Reverse our decision if the first exon is of low quality.
      //
      if ((qual < 85) && (dist >= 10000)) {
        if (trimLast == 0)
          fprintf(trimmedFile, "Trimming last exon based only on percent ID\n");
        trimLast = 1;
      }
    }

    if (trimmedFile && (trimFirst || trimLast)) {
      fprintf(trimmedFile, "------------------------------------------------------------BEFORE\n");
      printPolish(trimmedFile, p);
    }


    if (beVerbose) {
      polishesProcessed++;
      if (trimFirst || trimLast)
        polishesTrimmed++;
      if ((polishesProcessed % 10000) == 0) {
        fprintf(stderr, " %d processed, %d trimmed (%8.5f%%)\r",
                polishesProcessed, polishesTrimmed,
                100.0 * (double)polishesTrimmed / (double)polishesProcessed);
        fflush(stderr);
      }
    }


    //  If there is one intron, and we've been asked to remove
    //  either the first or the last (it should say to remove
    //  both), then remove the shorter of the two.
    //
    if ((trimFirst || trimLast) && (p->numExons == 2)) {
      trimFirst = 0;
      trimLast  = 0;

      if ((p->exons[0].estTo - p->exons[0].estFrom) > (p->exons[1].estTo - p->exons[1].estFrom))
        trimLast = 1;
      else
        trimFirst = 1;
    }



    //  Remove the first exon, by circularly shifting the list of
    //  exons.  The exon trimmed from the start is moved to the end of
    //  the exon list.
    //
    if (trimFirst) {
      int              i;
      sim4polishExon   save;

      memcpy(&save, p->exons+0, sizeof(sim4polishExon));

      for (i=1; i<p->numExons; i++)
        memcpy(p->exons+i-1, p->exons+i, sizeof(sim4polishExon));

      memcpy(p->exons+p->numExons-1, &save, sizeof(sim4polishExon));
             
      p->numExons--;
    }


    //  Trimming the last exon is easy; just decrement the size of the
    //  list.
    //
    if (trimLast) {
      p->numExons--;

      //  We also need to clear the intron orientation flag in the new
      //  last exon
      //
      p->exons[p->numExons-1].intronOrientation = INTRON_NONE;
    }

    if (trimmedFile && (trimFirst || trimLast)) {
      fprintf(trimmedFile, "------------------------------------------------------------AFTER\n");
      printPolish(trimmedFile, p);
      fprintf(trimmedFile, "============================================================EOP\n");
    }

    printPolish(stdout, p);

    //  Insert the exons back in, so they will be destroyed properly.
    //
    if (trimFirst)  p->numExons++;
    if (trimLast)   p->numExons++;

    destroyPolish(p);
  }

  return(0);
}

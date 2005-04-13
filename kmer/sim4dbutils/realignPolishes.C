#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "bio++.H"
#include "sim4.H"

//  This code takes basic sim4db format polishes and recomputes the
//  alignments and scores.  Required in the input polishes are the EST
//  id, genomic id, exon coordinates and an orientation.

int
main(int argc, char **argv) {

  //  Load all the sequences.  We really do need them all the ESTs in
  //  core, since they probably aren't in a useful sorted order.  You
  //  can probably figure out a way to get rid of the FastACache for
  //  the GEN.  Doing so will reduce memory usage by about 50%.

  FastACache *EST = 0L;
  FastACache *GEN = 0L;
  int         mergeTolerancePerc = 0;
  int         mergeToleranceBase = 0;
  int         statsOnly          = 0;

  //  Statistics on the exon merge

  int   mergedExons   = 0;
  int   mergedMatches = 0;

  int   numcdnagaps        = 0;
  int   nummatcheswithgaps = 0;

  FILE *mergeLog      = 0L;

  int     arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-merge", 2) == 0) {
      mergeTolerancePerc = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-b", 2) == 0) {
      mergeToleranceBase = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-M", 2) == 0) {
      mergeLog = fopen(argv[++arg], "w");
    } else if (strncmp(argv[arg], "-e", 2) == 0) {
      if (statsOnly)
        EST = new FastACache(argv[++arg], 1000, false, false);  //  debugging only!
      else 
        //EST = new FastACache(argv[++arg],    0, true);
        EST = new FastACache(argv[++arg], 1000, false, false);  //  debugging only!
    } else if (strncmp(argv[arg], "-g", 2) == 0) {
      GEN = new FastACache(argv[++arg],    1, false, true);
    } else if (strncmp(argv[arg], "-q", 2) == 0) {
      statsOnly = 1;
    }
    arg++;
  }

  if ((statsOnly == 0) && (!EST || !GEN)) {
    fprintf(stderr, "usage: %s [-merge percent-tolerance] [-M merge-log] [-q] -e est.fasta -g genome.fasta < polishes > somewhere\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "       Polishes _MUST_ be sorted by genomic index.\n");
    fprintf(stderr, "       If not, performance will be worse than atrocious.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       percent-tolerance -- merge exons separated by gap if\n");
    fprintf(stderr, "       the cDNA and genomic gaps differ by less than p percent.\n");
    fprintf(stderr, "       A value of 5 means 5%%\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -q: Don't actually do the work, just count the statistics\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    exit(1);
  }


  char   *s1 = new char [16 * 1024 * 1024];
  char   *s2 = new char [16 * 1024 * 1024];
  int     l1 = 0;
  int     l2 = 0;

  speedCounter  *C = new speedCounter("%12.0f polishes -- %12.0f polishes/second\r",
                                      1.0, 0xf, true);

  sim4polish *p = 0L;
  while ((p = s4p_readPolish(stdin)) != 0L) {

    //  If we have a mergeTolerance, merge adjacent exons that are
    //  separated my approximately equal sized cDNA and genomic gaps.
    //
    //  Possible a better way to do this is to check if the identity
    //  of the missing region is decent, too.

    //  Remember the id/cv of this guy for the log
    //
    double id = 0.0;
    double cv = 0.0;
    if (mergeLog) {
      id = s4p_percentIdentity(p);
      cv = s4p_percentCoverage(p);
    }

    int   merged = 0;
    int   gapped = 0;

    if ((mergeTolerancePerc > 0) || (mergeToleranceBase > 0)) {
      for (u32bit i=1; i<p->numExons; i++) {
        int cgap = p->exons[i].estFrom - p->exons[i-1].estTo;
        int ggap = p->exons[i].genFrom - p->exons[i-1].genTo;

        bool  mergeGap = false;

        //  New method -- check if the gaps are within 20bp of each other
        //
        int diff = cgap - ggap;
        if (diff < 0)
          diff = -diff;

        if (diff < mergeToleranceBase)
          mergeGap = true;


        //  Original method -- cehck if the gaps are within 10% of each other
        //
        int ctol = cgap * (100 + mergeTolerancePerc);
        int gtol = ggap * (100 + mergeTolerancePerc);

        cgap *= 100;
        ggap *= 100;

        if (((cgap < ggap) && (ctol > ggap)) ||
            ((ggap < cgap) && (gtol > cgap)))
          mergeGap = true;

        if (cgap > 1) {
          numcdnagaps++;
          gapped++;
        }

        if ((cgap > 1) && (mergeGap)) {

          //  Merge i and i-1 if adding in the tolerance makes either
          //  the cgap or the ggap longer than the other gap.  i.e., the
          //  cgap was shorter, but including the tolerance makes it
          //  longer, so they're about the same size.

          if (mergeLog)
            fprintf(mergeLog,
                    "MERGE: "u32bitFMTW(4)"-"u32bitFMTW(4)" (%6.2f,%6.2f) "u32bitFMTW(4)"-"u32bitFMTW(4)
                    " and "u32bitFMTW(8)"-"u32bitFMTW(8)" (%6.2f,%6.2f) "u32bitFMTW(8)"-"u32bitFMTW(8)"\n",
                    p->exons[i-1].estFrom, p->exons[i-1].estTo,
                    cgap / 100.0, ctol / 100.0,
                    p->exons[i].estFrom, p->exons[i].estTo,
                    p->exons[i-1].genFrom, p->exons[i-1].genTo,
                    ggap / 100.0, gtol / 100.0,
                    p->exons[i].genFrom, p->exons[i].genTo);

          //  merge exons
          p->exons[i-1].estTo = p->exons[i].estTo;
          p->exons[i-1].genTo = p->exons[i].genTo;

          //  delete this exon
          s4p_deleteExon(p, i);

          //  Do it again!
          i--;

          merged++;
          mergedExons++;
        }
      }

      if (merged)
        mergedMatches++;
      if (gapped)
        nummatcheswithgaps++;
    }


    //  For each exon, generate an alignment


    if (statsOnly == 0) {
      p->estLen   = EST->getSequence(p->estID)->sequenceLength();
      p->estPolyA = 0;
      p->estPolyT = 0;

      for (u32bit i=0; i<p->numExons; i++) {
        l1 = p->exons[i].estTo - p->exons[i].estFrom;
        l2 = p->exons[i].genTo - p->exons[i].genFrom;

        strncpy(s1, EST->getSequence(p->estID)->sequence() + p->exons[i].estFrom, l1);
        strncpy(s2, GEN->getSequence(p->genID)->sequence() + p->exons[i].genFrom, l2);

        if (p->matchOrientation == SIM4_MATCH_COMPLEMENT) {
          strncpy(s1, EST->getSequence(p->estID)->sequence() + p->estLen - p->exons[i].estTo, l1);
          reverseComplementSequence(s1, l1);
        }

        free(p->exons[i].estAlignment);
        free(p->exons[i].genAlignment);

        p->exons[i].estAlignment = (char *)malloc(sizeof(char) * (l1+l2+1));
        p->exons[i].genAlignment = (char *)malloc(sizeof(char) * (l1+l2+1));

        halign(s1, s2,
               l1, l2,
               p->exons[i].estAlignment,
               p->exons[i].genAlignment);
      }

      //  There isn't an intron after the last exon.  Force it.
      //
      p->exons[p->numExons-1].intronOrientation = SIM4_INTRON_NONE;

      s4p_updateAlignmentScores(p);
      s4p_printPolish(stdout, p, 0);
    }

    if (merged) {
      fprintf(mergeLog, "MERGED\tEST\t"u32bitFMT"\tfrom\t%8.3f\t%8.3f\tto\t%8.3f\t%8.3f\n",
              p->estID, id, cv, s4p_percentIdentity(p), s4p_percentCoverage(p));
    }

    s4p_destroyPolish(p);

    C->tick();
  }

  if ((mergeTolerancePerc > 0) || (mergeToleranceBase > 0)) {
    fprintf(stderr, "FOUND:  %d gaps in %d matches.\n", numcdnagaps, nummatcheswithgaps);
    fprintf(stderr, "MERGED: %d gaps in %d matches.\n", mergedExons, mergedMatches);
  }

  delete GEN;
  delete EST;

  return(0);
}

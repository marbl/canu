#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "bri++.H"
#include "fasta.H"
#include "fasta-cache.H"
#include "sim4polish.h"

#include "align.h"


void
align(const char *string1,
      const char *string2,
      const int len1,
      const int len2,
      char *alnline1,
      char *alnline2);


//
//  Converts blat's sim4-like output into sim4db format, including
//  alignments.  There is a helper perl-script that converts blat
//  sim4-like to a basic sim4db.  This code takes the basic sim4db and
//  recomputes the alignments and scores.
//

int
main(int argc, char **argv) {


  //  Usage:
  //     %s est.fasta genome.fasta

  //  Load all the sequences.  We really do need them all in core,
  //  since they probably aren't in a useful sorted order.  Maybe ESTs
  //  are.
  //
  fprintf(stderr, "Loading ESTs\n");
  FastACache *EST = new FastACache(argv[1],    0, true);

  fprintf(stderr, "Loading GENs\n");
  FastACache *GEN = new FastACache(argv[2], 1000, false);

  fprintf(stderr, "Begin\n");

  char   *s1 = new char [16 * 1024 * 1024];
  char   *s2 = new char [16 * 1024 * 1024];
  int     l1 = 0;
  int     l2 = 0;

  sim4polish *p = 0L;
  while ((p = s4p_readPolish(stdin)) != 0L) {
    s4p_printPolish(stderr, p, 0);

    p->estLen   = EST->getSequence(p->estID)->sequenceLength();
    p->estPolyA = 0;
    p->estPolyT = 0;

    //  For each exon, generate an alignment

    for (int i=0; i<p->numExons; i++) {
#if 0
      fprintf(stderr, "align exon %d %d-%d to %d-%d %c\n",
              i,
              p->exons[i].estFrom, p->exons[i].estTo,
              p->exons[i].genFrom, p->exons[i].genTo,
              (p->matchOrientation == SIM4_MATCH_FORWARD) ? 'f' : 'r');
#endif

      l1 = p->exons[i].estTo - p->exons[i].estFrom;
      l2 = p->exons[i].genTo - p->exons[i].genFrom;

      strncpy(s1, EST->getSequence(p->estID)->sequence() + p->exons[i].estFrom, l1);
      strncpy(s2, GEN->getSequence(p->genID)->sequence() + p->exons[i].genFrom, l2);

      if (p->matchOrientation == SIM4_MATCH_COMPLEMENT) {
        strncpy(s1, EST->getSequence(p->estID)->sequence() + p->estLen - p->exons[i].estTo, l1);
        reverseComplementSequence(s1, l1);
      }

      p->exons[i].estAlignment = (char *)malloc(sizeof(char) * (l1+l2+1));
      p->exons[i].genAlignment = (char *)malloc(sizeof(char) * (l1+l2+1));

      align(s1, s2,
            l1, l2,
            p->exons[i].estAlignment,
            p->exons[i].genAlignment);

#if 0
      fprintf(stderr, "a1 = %s\n", a1);
      fprintf(stderr, "a2 = %s\n", a2);
      fprintf(stderr, "done!\n");
#endif
    }

    s4p_updateAlignmentScores(p);
    s4p_printPolish(stdout, p, 0);
    s4p_destroyPolish(p);
  }


  delete GEN;
  delete EST;
}


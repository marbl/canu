#ifndef EXON_H
#define EXON_H

//  Define this only if Liliana is watching you.
//#define SPLSCORE

typedef struct exon {
  struct exon *next_exon;

  int     frGEN, toGEN;  //  Genomic coords
  int     frEST, toEST;  //  cDNA coords

  char    ori;
  int     length;

  // - 'flag' controls whether the second blast pass should be run 
  // or not on the adjacent (unmatched) fragment.
  // - 'edist' is an approximation for the error rate within the exon
  //
  int     flag;
  int     edist;

  int     percentID;
  int     alignmentLength;
  int     numMatches;
  int     numNs;
  int     numInDel;
  int     numEdits;
#ifdef SPLSCORE
  double  splScore;
#endif
} *Exon_ptr;

#endif  //  EXON_H

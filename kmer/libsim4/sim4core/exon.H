#ifndef EXON_H
#define EXON_H

#include "ckalloc.H"

//  Define this only if Liliana is watching you.
//#define SPLSCORE

typedef struct exon {
  struct exon *next_exon;

  int     frGEN, toGEN;  //  Genomic coords
  int     frEST, toEST;  //  cDNA coords

  int     ori;
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

typedef struct exon   Exon;


inline
Exon_ptr
new_exon(int f1, int f2, int t1, int t2, int len, int edist, int flag, Exon_ptr next) {
  Exon_ptr newthing = (Exon_ptr )ckalloc(sizeof(Exon));

  newthing->next_exon = next;

  newthing->frGEN = f1;
  newthing->frEST = f2;
  newthing->toGEN = t1;
  newthing->toEST = t2;

  newthing->ori = 'U';
  newthing->length = (len < 0) ? (t2-f2+1) : len;

  newthing->edist = edist;
  newthing->flag = flag;

  newthing->percentID = 0;
  newthing->alignmentLength = 0;
  newthing->numMatches = 0;
  newthing->numNs = 0;
  newthing->numInDel = 0;
  newthing->numEdits = 0;

#ifdef SPLSCORE
  newthing->splScore = -999999;
#endif

  return newthing;
}


#endif  //  EXON_H

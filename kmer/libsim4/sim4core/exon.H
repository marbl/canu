#ifndef EXON_H
#define EXON_H

#include <stdio.h>

//  Define this only if Liliana is watching you.
//#define SPLSCORE


class Exon {
public:
  Exon(int f1, int f2, int t1, int t2, int len, int edist, int flag, Exon *next);
  ~Exon();

  void  printList(char *label) {
    Exon *l = this;

    fprintf(stdout, label);
    while (l) {
      fprintf(stdout, "GEN f=%8d t=%8d  EST f=%8d t=%8d   flag=%d\n",
              l->frGEN,
              l->toGEN,
              l->frEST,
              l->toEST,
              l->flag);
    
      l = l->next_exon;
    }
    fprintf(stdout, "----------------------------------------\n");
    fflush(stdout);
  };

  Exon   *next_exon;

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
};

inline
Exon::Exon(int f1, int f2, int t1, int t2, int len, int edost, int flog, Exon *next) {
  next_exon = next;

  frGEN = f1;
  frEST = f2;
  toGEN = t1;
  toEST = t2;

  ori = 'U';
  length = (len < 0) ? (t2-f2+1) : len;

  edist = edost;
  flag = flog;

  percentID       = 0;
  alignmentLength = 0;
  numMatches      = 0;
  numNs           = 0;
  numInDel        = 0;
  numEdits        = 0;

#ifdef SPLSCORE
  splScore = -999999;
#endif
}

inline
Exon::~Exon() {
}


inline
void
freeExon(Exon *&l) {
  delete l;
  l = 0L;
}

inline
void 
freeExonList(Exon *&l) {
  Exon *t;

  while (l) {
    t = l->next_exon;
    freeExon(l);
    l = t;
  }

  l = 0L;
}    



#endif  //  EXON_H

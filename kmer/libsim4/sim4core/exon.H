#ifndef EXON_H
#define EXON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//  Define this only if Liliana is watching you.
//#define SPLSCORE

class Exon {
private:
  Exon() {
    next_exon = 0L;
    frGEN = 0;
    frEST = 0;
    toGEN = 0;
    toEST = 0;
    ori    = 'U';
    length = 0;

    edist = 0;
    flag  = 0;

    percentID       = 0;
    alignmentLength = 0;
    numMatches      = 0;
    numNs           = 0;
    numInDel        = 0;
    numEdits        = 0;

#ifdef SPLSCORE
    splScore = -999999;
#endif
  };
  ~Exon() {};

  void    init(int f1, int f2, int t1, int t2, int len, int edost, int flog, Exon *next) {
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
  };

  friend class exonManager;

public:
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






//  Just for allocating exons.  It takes care of deleting exons.  You
//  are not allowed to delete them.
//
class exonManager {
  
public:
  exonManager() {
    _listLen = 1;
    _listMax = 8;
    _list    = new Exon* [_listMax];

    _curLen  = 0;
    _curMax  = 256;
    _list[0] = new Exon [_curMax];
  };
  ~exonManager() {
    for (u32bit i=0; i<_listLen; i++)
      delete [] _list[i];
    delete [] _list;
  };

  Exon  *newExon(int f1, int f2, int t1, int t2, int len, int edist, int flag, Exon *next) {

    //  If the current list is full, move to the next one, or allocate
    //  more space.
    //
    if (_curLen >= _curMax) {
      if (_listLen >= _listMax) {
        _listMax *= 2;
        Exon **l = new Exon* [_listMax];
        memcpy(l, _list, sizeof(Exon**) * _listLen);
        delete [] _list;
        _list = l;
      }

      _list[_listLen++] = new Exon [_curMax];
      _curLen = 0;
    }

    Exon *e = _list[_listLen - 1] + _curLen;
    _curLen++;
    e->init(f1, f2, t1, t2, len, edist, flag, next);

    //fprintf(stderr, "exonManager::newExon()--  return exon at %p\n", e);

    return(e);;
  };

private:
  //  Exon pointers are valid throughout the lifetime of execution,
  //  so we can't use realloc here.  Instead, we keep a list of arrays
  //  of exons.
  //
  u32bit   _listLen;
  u32bit   _listMax;
  Exon   **_list;

  u32bit   _curLen;
  u32bit   _curMax;
};


#endif  //  EXON_H

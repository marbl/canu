


#ifndef EXONMANAGER_H
#define EXONMANAGER_H

struct exon {
  int          from1, from2, to1, to2;
  int          min_diag, max_diag;
  int          match;
  char         ori;
  int          length;
  int          flag;
  int          ematches;
  int          edist;
  struct exon *next_exon;
};



class exonManager {
private:
  unsigned int   _allocExons;
  unsigned int   _numExons;
  exon          *_allExons;
public:
  exonManager() {
    _allocExons   = 16384;
    _numExons     = 0;
    _allExons     = new msp [_allocExons];
  };

  ~exonManager() {
    delete [] _allExons;
  };


  void   addExon(int l, int p1, int p2, int sc) {

    //  Allocate more Exons, if we need to.
    //
    if (_numExons >= _allocExons) {
      _allocExons *= 2;

      msp *n = new msp [_allocExons];

      for (unsigned int i=0; i<_numExons; i++) {
        n[i].len          = _allExons[i].len;
        n[i].pos1         = _allExons[i].pos1;
        n[i].pos2         = _allExons[i].pos2;
        n[i].score        = _allExons[i].score;
        n[i].linkingScore = _allExons[i].linkingScore;
        n[i].prev         = _allExons[i].prev;
      }

      delete _allExons;
      _allExons = n;
    }

    _allExons[_numExons].len          = l;
    _allExons[_numExons].pos1         = p1;
    _allExons[_numExons].pos2         = p2;
    _allExons[_numExons].score        = sc;
    _allExons[_numExons].linkingScore = 0;
    _allExons[_numExons].prev         = 0;

    _numExons++;
  };


  //  Clear the msp list -- doesn't delete them, just puts them on
  //  a free list.
  //
  void      clear(void) {
    _numExons = 0;
  };

};




#endif  //  EXONMANAGER_H

#ifndef HITMATRIX_H
#define HITMATRIX_H

#include "libbri.H"
#include "positionDB.H"
#include "intervalList.H"

//  $Id$

struct diagonalLine {
  u32bit   _qsPos;
  u32bit   _dsPos;
  u32bit   _diagonalID;
};

class trapMatch {
public:
  u32bit         _qsLo, _qsHi;
  u32bit         _dsLo, _dsHi;
  u32bit         _numMatches;
  char           _direction;
  intervalList  *_IL;
  trapMatch     *_next;

  trapMatch(u32bit         qsLo,
            u32bit         qsHi,
            u32bit         dsLo,
            u32bit         dsHi,
            intervalList  *IL,
            char           direction) {
    _qsLo       = qsLo;
    _qsHi       = qsHi;
    _dsLo       = dsLo;
    _dsHi       = dsHi;
    _IL         = IL;
    _direction  = direction;
    _next       = 0L;
  };
};

class hitMatrix {
private:
  u32bit                        _qsLen;   //  Seq Len of Q
  u32bit                        _qsMers;  //  Valid mers in Q
  u32bit                        _qsIdx;   //  Index of Q in the FastA

  //  Instead of building the lines during add(), we store
  //  the information used to build lines, and then build them
  //  in chain().  This was done to reduce simultaneous memory
  //  usage, as the lineArrayMap and etc take up considerable space.
  //
  u32bit                        _hitsLen;
  u32bit                        _hitsMax;
  diagonalLine                 *_hits;


  //  Making sense of the raw output from the search is not a trivial
  //  task for perl.  SMALL searches (dbEST vs 0.5MB sequence) used more
  //  than 4GB of memory in perl.
  //
  //  So, we bite the bullet and do it here.
  //
  //  _matches is a sorted linked list of the regions we have found.
  //  The list is kept in REVERSE order, as we usually add regions
  //  in the correct order (correct reverse order), occasionally
  //  we need to swap the last two.
  //
  //  The list is deleted in filter()
  //
  trapMatch                    *_matches;

  void addMatch(u32bit         qsLo,
                u32bit         qsHi,
                u32bit         dsLo,
                u32bit         dsHi,
                intervalList  *IL,
                char           direction);

public:
  hitMatrix(u32bit qsLen, u32bit q, u32bit qsIdx);
  ~hitMatrix();

  void    addHit(u32bit qi,
                 u32bit ps);
  void    addHits(u32bit  qi,
                  u64bit *ps,
                  u64bit  cn);

  void    sort_diagonal(void);
  void    sort_dsPos(void);

  void    filter(char direction, char *&theOutput, u32bit &theOutputPos, u32bit &theOutputMax);
};


inline
void
hitMatrix::addHits(u32bit  qi,
                   u64bit *ps,
                   u64bit  cn) {

  if ((_hitsLen + cn) >= _hitsMax) {
    _hitsMax = _hitsMax + _hitsMax + (u32bit)cn;

    diagonalLine  *h = new diagonalLine [_hitsMax];

    for (u32bit z=_hitsLen; z--; ) {
      h[z]._qsPos      = _hits[z]._qsPos;
      h[z]._dsPos      = _hits[z]._dsPos;
      h[z]._diagonalID = _hits[z]._diagonalID;
    }
          
    delete [] _hits;

    _hits = h;
  }

  for (u64bit i=0; i<cn; i++) {
    //fprintf(stderr, "Hit: %u %lu\n", qi, ps[i]);

    _hits[_hitsLen]._qsPos      = (u32bit)(qi);
    _hits[_hitsLen]._dsPos      = (u32bit)(ps[i]);
    _hits[_hitsLen]._diagonalID = (u32bit)(_qsLen - qi - 1 + ps[i]);
    _hitsLen++;
  }
}


inline
void
hitMatrix::addHit(u32bit qi,
                  u32bit ps) {

  if (_hitsLen >= _hitsMax) {
    _hitsMax <<= 1;

    diagonalLine  *h = new diagonalLine [_hitsMax];

    for (u32bit z=_hitsLen; z--; ) {
      h[z]._qsPos      = _hits[z]._qsPos;
      h[z]._dsPos      = _hits[z]._dsPos;
      h[z]._diagonalID = _hits[z]._diagonalID;
    }
          
    delete [] _hits;

    _hits = h;
  }

  _hits[_hitsLen]._qsPos      = qi;
  _hits[_hitsLen]._dsPos      = ps;
  _hits[_hitsLen]._diagonalID = _qsLen - qi - 1 + ps;
        
  _hitsLen++;
}



#endif  //  HITMATRIX_H

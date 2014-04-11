#ifndef HITMATRIX_H
#define HITMATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "bio++.H"
#include "positionDB.H"
#include "encodedQuery.H"

//  Define this to cut the space required for storing hits by 1/3 --
//  from 12 byyes to 8 bytes -- at a slight computational expense --
//  negligible on real hardware, I hope.
//
//  The original definition of diagonalID was
//      qsLen - qsPos - 1 + dsPos
//  but qsLen is fixed for everyone, so we could reduce it to
//      dsPos - qsPos
//  but that's not unsigned.
//
//  Results: on a human mapping, using chromosomes as the stream and
//  the whole human as the table (so we need to actually store a large
//  number of hits), we see a savings of 2GB and a small drop in
//  runtime.  Process size went from 20.7GB to 18.7GB, CPU time from
//  20578 to 20193 seconds (833MHz EV6.8AL (21264B)).
//
#define WITHOUT_DIAGONALID


struct diagonalLine {
  uint32   _qsPos;
  uint32   _dsPos;
#ifndef WITHOUT_DIAGONALID
  uint32   _diagonalID;
#endif
};

class trapMatch {
public:
  uint32         _qsLo, _qsHi;
  uint32         _dsLo, _dsHi;
  merCovering   *_IL;
  trapMatch     *_next;

  trapMatch(uint32         qsLo,
            uint32         qsHi,
            uint32         dsLo,
            uint32         dsHi,
            merCovering   *IL) {
    _qsLo       = qsLo;
    _qsHi       = qsHi;
    _dsLo       = dsLo;
    _dsHi       = dsHi;
    _IL         = IL;
    _next       = 0L;
  };
};

class hitMatrix {
public:
  hitMatrix(uint32 qsLen,
            uint32 qsMers,
            uint32 qsIdx);
  ~hitMatrix();

  void    addHits(uint32  qi,
                  uint64 *ps,
                  uint64  cn);

  void    sort_diagonal(void);
  void    sort_dsPos(void);

  void    filter(encodedQuery *query, bool isReverse);
private:
  uint32                        _qsLen;   //  Seq Len of Q
  uint32                        _qsMers;  //  Valid mers in Q
  uint32                        _qsIdx;   //  Index of Q in the FastA

  //  Instead of building the lines during add(), we store
  //  the information used to build lines, and then build them
  //  in chain().  This was done to reduce simultaneous memory
  //  usage, as the lineArrayMap and etc take up considerable space.
  //
  uint32                        _hitsLen;
  uint32                        _hitsMax;
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

  void addMatch(uint32         qsLo,
                uint32         qsHi,
                uint32         dsLo,
                uint32         dsHi,
                merCovering   *IL);

};


inline
void
hitMatrix::addHits(uint32  qi,
                   uint64 *ps,
                   uint64  cn) {

  if ((_hitsLen + cn) >= _hitsMax) {
    _hitsMax = _hitsMax + _hitsMax + (uint32)cn;

    diagonalLine *h;
    try {
      h = new diagonalLine [_hitsMax];
    } catch (std::bad_alloc) {
      fprintf(stderr, "hitMatrix::addHits()-- caught std::bad_alloc in %s at line %d.\n", __FILE__, __LINE__);
      fprintf(stderr, "hitMatrix::addHits()-- have "uint32FMT" hits, tried to add "uint64FMT" more\n", _hitsLen, cn);
      exit(1);
    }

    for (uint32 z=_hitsLen; z--; ) {
      h[z]._qsPos      = _hits[z]._qsPos;
      h[z]._dsPos      = _hits[z]._dsPos;
#ifndef WITHOUT_DIAGONALID
      h[z]._diagonalID = _hits[z]._diagonalID;
#endif
    }
          
    delete [] _hits;

    _hits = h;
  }

  for (uint64 i=0; i<cn; i++) {
    _hits[_hitsLen]._qsPos      = (uint32)(qi);
    _hits[_hitsLen]._dsPos      = (uint32)(ps[i]);
#ifndef WITHOUT_DIAGONALID
    _hits[_hitsLen]._diagonalID = (uint32)(_qsLen - qi - 1 + ps[i]);
#endif
    _hitsLen++;
  }
}


#endif  //  HITMATRIX_H

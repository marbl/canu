#ifndef HITMATRIX_H
#define HITMATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "bio++.H"
#include "positionDB.H"
#include "filterObj.H"

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
  u32bit   _qsPos;
  u32bit   _dsPos;
#ifndef WITHOUT_DIAGONALID
  u32bit   _diagonalID;
#endif
};


class hitMatrix {
public:
  hitMatrix(u32bit qsLen,
            u32bit qsIdx);
  ~hitMatrix();

  void    addHits(u32bit  qi,
                  u64bit *ps,
                  u64bit  cn);

  void    sort_diagonal(void);
  void    sort_dsPos(void);

  void    processMatrix(char direction, filterObj *FO);

private:
  u32bit                        _qsLen;   //  Seq Len of Q
  u32bit                        _qsIdx;   //  Index of Q in the FastA

  //  Instead of building the lines during add(), we store
  //  the information used to build lines, and then build them
  //  in chain().  This was done to reduce simultaneous memory
  //  usage, as the lineArrayMap and etc take up considerable space.
  //
  u32bit                        _hitsLen;
  u32bit                        _hitsMax;
  diagonalLine                 *_hits;
};


inline
void
hitMatrix::addHits(u32bit  qi,
                   u64bit *ps,
                   u64bit  cn) {

  if ((_hitsLen + cn) >= _hitsMax) {
    _hitsMax = _hitsMax + _hitsMax + (u32bit)cn;

    diagonalLine *h;
    try {
      h = new diagonalLine [_hitsMax];
    } catch (std::bad_alloc) {
      fprintf(stderr, "hitMatrix::addHits()-- caught std::bad_alloc in %s at line %d.\n", __FILE__, __LINE__);
      fprintf(stderr, "hitMatrix::addHits()-- have "u32bitFMT" hits, tried to add "u64bitFMT" more\n", _hitsLen, cn);
      exit(1);
    }

    for (u32bit z=_hitsLen; z--; ) {
      h[z]._qsPos      = _hits[z]._qsPos;
      h[z]._dsPos      = _hits[z]._dsPos;
#ifndef WITHOUT_DIAGONALID
      h[z]._diagonalID = _hits[z]._diagonalID;
#endif
    }
          
    delete [] _hits;

    _hits = h;
  }

  for (u64bit i=0; i<cn; i++) {
    _hits[_hitsLen]._qsPos      = (u32bit)(qi);
    _hits[_hitsLen]._dsPos      = (u32bit)(ps[i]);
#ifndef WITHOUT_DIAGONALID
    _hits[_hitsLen]._diagonalID = (u32bit)(_qsLen - qi - 1 + ps[i]);
#endif
    _hitsLen++;
  }
}


#endif  //  HITMATRIX_H

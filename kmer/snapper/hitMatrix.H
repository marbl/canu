#ifndef HITMATRIX_H
#define HITMATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "bio++.H"
#include "positionDB.H"


//  aHit -- storing the internal hits
//
//  _status
//  & 0x00000001 -- direction, forward if set, otherwise reverse
//    0x00000002 -- not filtered, if set, polish the hit
//    0x00000004 -- not filtered, if set, polish the hit because it has something unique-ish
//    0x00000008 -- match verified via polishng
//    0x00000010 -- match discarded via hit refinement
//    0x0000fff8 -- unused
//    0x00ff0000 -- percent idendity of match
//    0xff000000 -- percent coverage of match
//
#define AHIT_DIRECTION_MASK  0x00000001
#define AHIT_POLISHABLE      0x00000002
#define AHIT_HAS_UNIQUE      0x00000004
#define AHIT_VERIFIED        0x00000008
#define AHIT_DISCARDED       0x00000010


struct diagonalLine {
  union {
    u64bit  all;
    struct {
      u64bit     uniq : 10;  //  uniqueness score for the mer here
      u64bit     qPos : 22;  //  position in the query, 0 to 4M
      u64bit     dPos : 32;  //  position in the genome, 0 to 4G
    } val;
  };
};


struct aHit {
  u32bit     _status;
  u32bit     _qsIdx;
  u32bit     _dsIdx;
  u32bit     _dsLo;
  u32bit     _dsHi;
  u32bit     _covered;
  u32bit     _matched;
  u32bit     _numMers;
  merList   *_ML;
};


class trapMatch {
public:
  u32bit         _unique;
  u32bit         _qsLo;
  u32bit         _qsHi;
  u32bit         _dsLo;
  u32bit         _dsHi;

  merCovering   *_IL;
  merList       *_ML;

  trapMatch     *_next;

  trapMatch(u32bit         isunique,
            u32bit         qsLo,
            u32bit         qsHi,
            u32bit         dsLo,
            u32bit         dsHi,
            merCovering   *IL,
            merList       *ML) {
    _unique     = isunique;
    _qsLo       = qsLo;
    _qsHi       = qsHi;
    _dsLo       = dsLo;
    _dsHi       = dsHi;
    _IL         = IL;
    _ML         = ML;
    _next       = 0L;
  };
};




class hitMatrix {
public:
  hitMatrix(u32bit qsLen, u32bit q, u32bit qsIdx);
  ~hitMatrix();

  void    addHits(u32bit  qi,
                  u64bit *ps,
                  u64bit  cn,
                  u64bit  ad=0);

  void    sort_diagonal(void);
  void    sort_dsPos(void);

  void    filter(char      direction,
                 double    minHitCoverage,
                 u32bit    minHitLength,
                 aHit    *&theOutput,
                 u32bit   &theOutputPos,
                 u32bit   &theOutputMax);
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

  void addMatch(u32bit         isunique,
                u32bit         qsLo,
                u32bit         qsHi,
                u32bit         dsLo,
                u32bit         dsHi,
                merCovering   *IL,
                merList       *ML);

};


inline
void
hitMatrix::addHits(u32bit  qi,
                   u64bit *ps,
                   u64bit  cn,
                   u64bit  ad) {

  if ((_hitsLen + cn) >= _hitsMax) {
    _hitsMax = _hitsMax + _hitsMax + (u32bit)cn;

    diagonalLine *h;
    try {
      h = new diagonalLine [_hitsMax];
    } catch (std::bad_alloc) {
      fprintf(stderr, "hitMatrix::addHits()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
      fprintf(stderr, "hitMatrix::addHits()-- have "u32bitFMT" hits, tried to add "u64bitFMT" more\n", _hitsLen, cn);
      exit(1);
    }

    memcpy(h, _hits, sizeof(u64bit) * _hitsLen);
    delete [] _hits;
    _hits = h;
  }

  u64bit  uniq = cn;
  if (ad > 0)
    uniq = ad;
  if (uniq > 0x000003ff)
    uniq = 0x000003ff;

  for (u64bit i=0; i<cn; i++) {
    _hits[_hitsLen].val.uniq = uniq;
    _hits[_hitsLen].val.qPos = qi;
    _hits[_hitsLen].val.dPos = ps[i];
    _hitsLen++;
  }
}


#endif  //  HITMATRIX_H

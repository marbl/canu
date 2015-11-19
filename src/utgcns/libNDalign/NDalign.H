
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/overlapInCore/overlapAlign.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-JUN-17 to 2015-AUG-25
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef NDALIGN_H
#define NDALIGN_H

#include "AS_global.H"

#include "ovStore.H"

#include "NDalgorithm.H"

#include "AS_UTL_reverseComplement.H"

#include <map>
#include <vector>
#include <algorithm>

using namespace std;

class NDalignResult {
public:
  NDalignResult() {
    clear();

    _deltaMax = 0;
    _delta    = NULL;
  };
  ~NDalignResult() {
    delete [] _delta;
  };

  //  Accessors

  double  score(void)    { return(_olapScore); };

  int32   deltaLen(void) { return(_deltaLen); };
  int32  *delta(void)    { return(_delta);    };

  //  Setters

  void       clear(void) {
    _aLo       = 0;
    _aHi       = 0;

    _bLo       = 0;
    _bHi       = 0;

    _olapLen   = 0;
    _olapScore = DBL_MIN;
    _olapErate = 0;
    _olapType  = pedBothBranch;

    _deltaLen = 0;
  };

  void       save(int32 aLo, int32 aHi, int32 bLo, int32 bHi, int32 score, pedOverlapType type, int32 dl, int32 *d) {
    _aLo       = aLo;
    _aHi       = aHi;

    _bLo       = bLo;
    _bHi       = bHi;

    _olapLen   = ((aHi - aLo) + (bHi - bLo) + dl) / 2;
    _olapScore = score;
    _olapErate = DBL_MAX;
    _olapType  = type;

    _deltaLen = dl;

    if (_deltaMax < _deltaLen) {
      delete [] _delta;
      _deltaMax = _deltaLen + 4096;
      _delta    = new int32 [_deltaMax];
    }

    memcpy(_delta, d, sizeof(int32) * _deltaLen);
  };

  void       setErate(double olapErate) {
    _olapErate = olapErate;
  };

  int32               _aLo;
  int32               _aHi;

  int32               _bLo;
  int32               _bHi;

  int32               _olapLen;
  double              _olapScore;
  double              _olapErate;
  pedOverlapType      _olapType;

  int32               _deltaMax;
  int32               _deltaLen;
  int32              *_delta;
};



class NDalign {
public:
  NDalign(pedAlignType  alignType,
               double        maxErate,
               int32         merSize);
  ~NDalign();

  void             initialize(uint32 aID, char *aStr, int32 aLen, int32 aLo, int32 aHi,
                              uint32 bID, char *bStr, int32 bLen, int32 bLo, int32 bHi, bool bFlipped);

  //  Algorithm

  bool             findMinMaxDiagonal(int32 minLength,
                                      uint32 AbgnAdj=0, uint32 AendAdj=0,
                                      uint32 BbgnAdj=0, uint32 BendAdj=0);
  bool             findSeeds(bool dupIgnore);
  bool             findHits(void);
  bool             chainHits(void);

  bool             makeNullHit(void);

  bool             processHits(void);

  //  Diagnostic

  bool             scanDeltaForBadness(bool verbose=false, bool showAlign=false);
  void             realignForward(bool verbose=false, bool displayAlign=false);
  void             realignBackward(bool verbose=false, bool displayAlign=false);

  void             display(char    *prefix,
                           int32    aLo,   int32   aHi,
                           int32    bLo,   int32   bHi,
                           int32   *delta, int32   deltaLen,
                           bool     displayIt = true,
                           bool     saveStats = false);

  void             display(char    *prefix,
                           bool     displayIt = true);

  //  Result reporting

  char            *astr(void)     { return(_aStr); };
  char            *bstr(void)     { return(_bStr); };

  int32            abgn(void)     { return(_bestResult._aLo); };
  int32            aend(void)     { return(_bestResult._aHi); };
  int32            bbgn(void)     { return((_bFlipped == false) ? (_bestResult._bLo) : (_bestResult._bLo)); };
  int32            bend(void)     { return((_bFlipped == false) ? (_bestResult._bHi) : (_bestResult._bHi)); };

  int32            ahg5(void)     { return(        _bestResult._aLo); };
  int32            ahg3(void)     { return(_aLen - _bestResult._aHi); };
  int32            bhg5(void)     { return((_bFlipped == false) ? (        _bestResult._bLo) : (_bLen - _bestResult._bLo)); };
  int32            bhg3(void)     { return((_bFlipped == false) ? (_bLen - _bestResult._bHi) : (        _bestResult._bHi)); };

  int32            length(void)   { return(_bestResult._olapLen);  };
  double           erate(void)    { return(_bestResult._olapErate); };
  int32            score(void)    { return(_bestResult._olapScore); };

  pedOverlapType   type(void)     { return(_bestResult._olapType); };

  int32            deltaLen(void) { return(_bestResult.deltaLen()); };
  int32           *delta(void)    { return(_bestResult.delta());    };

private:

  class exactMatch {
  public:
    exactMatch(int32 a, int32 b, int32 l) {
      aBgn = a;
      bBgn = b;
      tLen = l;
    };

    int32  aBgn;  //  Signed to allow for easy compute of diagonal
    int32  bBgn;
    int32  tLen;

    bool operator<(exactMatch const that) const {
      if (tLen == that.tLen)
        return(aBgn < that.aBgn);

      return(tLen > that.tLen);
    };
  };

  //  Parameters of the alignment

  pedAlignType        _alignType;
  double              _maxErate;
  int32               _merSizeInitial;

  //  From the parameters.

  uint32              _aID;
  char               *_aStr;
  int32               _aLen;
  int32               _aLoOrig;
  int32               _aHiOrig;

  uint32              _bID;
  char               *_bStr;
  int32               _bLen;
  int32               _bLoOrig;
  int32               _bHiOrig;

  //  If flipped, we need to save the reverse complement somewhere

  bool                _bFlipped;
  uint32              _bRevMax;
  char               *_bRev;

  //  Computed stuff, for each alignment.

  NDalgorithm        *_editDist;

  int32               _minDiag;
  int32               _maxDiag;

  int32               _merSize;

  map<uint64,int32>   _aMap;  //  Signed, to allow for easy compute of diagonal
  map<uint64,int32>   _bMap;

  vector<exactMatch>  _rawhits;
  vector<exactMatch>  _hits;

  uint32              _hitr;

  //  The result.

  NDalignResult       _bestResult;

  //  For displaying and stats.

  char               *_topDisplay;
  char               *_botDisplay;
  char               *_resDisplay;

  uint32              _matches;
  uint32              _errors;
  uint32              _gapmatches;
  uint32              _freegaps;

  //  Constants

  uint64  acgtToBit[256];
  uint64  acgtToVal[256];
  uint64  merMask[33];

  void    fastFindMersA(bool dupIgnore);
  void    fastFindMersB(bool dupIgnore);
};

#endif  //  NDALIGN_H

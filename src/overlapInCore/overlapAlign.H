
#ifndef OVERLAPALIGN_H
#define OVERLAPALIGN_H

#include "AS_global.H"

#include "ovStore.H"

#include "prefixEditDistance.H"

#include "AS_UTL_reverseComplement.H"

#include <map>
#include <vector>
#include <algorithm>

using namespace std;

class overlapAlignResult {
public:
  overlapAlignResult() {
    clear();

    _deltaMax = 0;
    _delta    = NULL;
  };
  ~overlapAlignResult() {
    delete [] _delta;
  };

  //  Accessors

  double  score(void)    { return(_olapLen * (1 - _olapQual)); };

  int32   deltaLen(void) { return(_deltaLen); };
  int32  *delta(void)    { return(_delta);    };

  //  Setters

  void       clear(void) {
    _aLo      = 0;
    _aHi      = 0;

    _bLo      = 0;
    _bHi      = 0;

    _olapLen  = 0;
    _olapQual = 0;

    _deltaLen = 0;
  };

  void       saveCoords(int32 aLo, int32 aHi, int32 bLo, int32 bHi, int32 olapLen, double olapQual) {
    _aLo      = aLo;
    _aHi      = aHi;

    _bLo      = bLo;
    _bHi      = bHi;

    _olapLen  = olapLen;
    _olapQual = olapQual;
  };

  void       saveDelta(int32 dl, int32 *d) {
    _deltaLen = dl;

    if (_deltaMax < _deltaLen) {
      delete [] _delta;
      _deltaMax = _deltaLen + 4096;
      _delta    = new int32 [_deltaMax];
    }

    memcpy(_delta, d, sizeof(int32) * _deltaLen);
  }
  
  int32               _aLo;
  int32               _aHi;

  int32               _bLo;
  int32               _bHi;

  int32               _olapLen;
  double              _olapQual;

  int32               _deltaMax;
  int32               _deltaLen;
  int32              *_delta;
};



class overlapAlign {
public:
  overlapAlign(bool    partialOverlaps,
               double  maxErate,
               int32   merSize);
  ~overlapAlign();

  void  initialize(uint32 aID, char *aStr, int32 aLen, int32 aLo, int32 aHi,
                   uint32 bID, char *bStr, int32 bLen, int32 bLo, int32 bHi, bool bFlipped);

  //  Algorithm

  bool  findMinMaxDiagonal(int32 minLength);
  bool  findSeeds(bool dupIgnore);
  bool  findHits(void);
  bool  chainHits(void);
  bool  processHits(void);

  //  Diagnostic

  void    display(bool withAlign=true);

  //  Result reporting

  char   *astr(void)    { return(_aStr); };
  char   *bstr(void)    { return(_bStr); };

  int32   abgn(void)    { return(_bestResult._aLo); };
  int32   aend(void)    { return(_bestResult._aHi); };
  int32   bbgn(void)    { return((_bFlipped == false) ? (_bestResult._bLo) : (_bestResult._bLo)); };
  int32   bend(void)    { return((_bFlipped == false) ? (_bestResult._bHi) : (_bestResult._bHi)); };

  int32   ahg5(void)    { return(        _bestResult._aLo); };
  int32   ahg3(void)    { return(_aLen - _bestResult._aHi); };
  int32   bhg5(void)    { return((_bFlipped == false) ? (        _bestResult._bLo) : (_bLen - _bestResult._bLo)); };
  int32   bhg3(void)    { return((_bFlipped == false) ? (_bLen - _bestResult._bHi) : (        _bestResult._bHi)); };

  int32   length(void)   { return(_bestResult._olapLen);  };
  double  erate(void)    { return(_bestResult._olapQual); };

  int32   deltaLen(void) { return(_bestResult.deltaLen()); };
  int32  *delta(void)    { return(_bestResult.delta());    };

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

  bool                _partialOverlaps;
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

  prefixEditDistance *_editDist;

  int32               _minDiag;
  int32               _maxDiag;

  int32               _merSize;

  map<uint64,int32>   _aMap;  //  Signed, to allow for easy compute of diagonal
  map<uint64,int32>   _bMap;

  vector<exactMatch>  _rawhits;
  vector<exactMatch>  _hits;

  //  The result.

  overlapAlignResult  _bestResult;

  //  Constants

  uint64  acgtToBit[256];
  uint64  acgtToVal[256];
  uint64  merMask[33];

  void    fastFindMersA(bool dupIgnore);
  void    fastFindMersB(bool dupIgnore);
};

#endif  //  OVERLAPALIGN_H


#include "AS_global.H"

#include "ovStore.H"

#include "prefixEditDistance.H"
//#include "Binomial_Bound.H"  //  liboverlap
#include "AS_UTL_reverseComplement.H"

#include <map>
#include <vector>
#include <algorithm>

using namespace std;


class overlapAlign {
public:
  overlapAlign(bool    partialOverlaps,
              double  maxErate,
              int32   merSize) {
    _partialOverlaps  = partialOverlaps;
    _maxErate         = maxErate;
    _merSizeInitial   = merSize;
    _editDist         = new prefixEditDistance(_partialOverlaps, _maxErate);

    _deltaMax         = 16384;
    _deltaLen         = 0;
    _delta            = new int32 [_deltaMax];
  };

  ~overlapAlign() {
    delete    _editDist;
    delete [] _delta;
  };



  void  initialize(uint32 aID, char *aStr, int32 aLen, int32 aLo, int32 aHi,
                   uint32 bID, char *bStr, int32 bLen, int32 bLo, int32 bHi, bool bFlipped) {

    _aID  = aID;
    _aStr = aStr;
    _aLen = _aLenOrig = aLen;
    _aLo  = _aLoOrig  = aLo;
    _aHi  = _aHiOrig  = aHi;

    _bID  = bID;
    _bStr = bStr;
    _bLen = _bLenOrig = bLen;
    _bLo  = _bLoOrig  = bLo;
    _bHi  = _bHiOrig  = bHi;

    _bFlipped  = bFlipped;

    if (_bFlipped == true) {
      if (_bRevMax < _bLen) {
        delete _bRev;
        _bRevMax = _bLen + 1000;
        _bRev    = new char [_bRevMax];
      }

      memcpy(_bRev, bStr, sizeof(char) * (_bLen + 1));

      reverseComplementSequence(_bRev, _bLen);

      _bStr = _bRev;

      _bLo = _bLen - bLo;  //  Now correct for the reverse complemented sequence
      _bHi = _bLen - bHi;
    }

    assert(_aLo < _aHi);
    assert(_bLo < _bHi);

    //_editDist doesn't need to be cleared.

    _minDiag = 0;
    _maxDiag = 0;

    _merSize = 0;

    _aMap.clear();
    _bMap.clear();

    _rawhits.clear();
    _hits.clear();
  };

  //  Algorithm

  bool  findMinMaxDiagonal(int32 minLength);
  bool  findSeeds(bool dupIgnore);
  void  findHits(void);
  void  chainHits(void);
  bool  processHits(ovOverlap *result);

  //  Result reporting

  int32   abgn(void)    { return(_aLo); };
  int32   aend(void)    { return(_aHi); };
  int32   bbgn(void)    { return((_bFlipped == false) ? (_bLo) : (_bLo)); };
  int32   bend(void)    { return((_bFlipped == false) ? (_bHi) : (_bHi)); };

  int32   ahg5(void)    { return(        _aLo); };
  int32   ahg3(void)    { return(_aLen - _aHi); };
  int32   bhg5(void)    { return((_bFlipped == false) ? (        _bLo) : (_bLen - _bLo)); };
  int32   bhg3(void)    { return((_bFlipped == false) ? (_bLen - _bHi) : (        _bHi)); };

  //result->dat.ovl.ahg5   = (_bFlipped == false) ? (        _aLo    ) : (        _aLo    );
  //result->dat.ovl.ahg3   = (_bFlipped == false) ? (_aLen - _aHi - 1) : (_aLen - _aHi - 1);
  //result->dat.ovl.bhg5   = (_bFlipped == false) ? (        _bLo    ) : (_bLen - _bLo - 1);
  //result->dat.ovl.bhg3   = (_bFlipped == false) ? (_bLen - _bHi - 1) : (        _bHi    );

  int32   length(void)  { return(_olapLen);  };
  double  erate(void)   { return(_olapQual); };

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

  //  From the parameters.  Updated with the result.

  uint32              _aID;
  char               *_aStr;
  int32               _aLen,  _aLenOrig;
  int32               _aLo,   _aLoOrig;
  int32               _aHi,   _aHiOrig;

  uint32              _bID;
  char               *_bStr;
  int32               _bLen,  _bLenOrig;
  int32               _bLo,   _bLoOrig;
  int32               _bHi,   _bHiOrig;

  int32               _olapLen;
  double              _olapQual;

  int32               _deltaMax;
  int32               _deltaLen;
  int32              *_delta;

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

  //  Constants

  uint64  acgtToBit[256];
  uint64  acgtToVal[256];
  uint64  merMask[33];

  void    initializeConstants(void);

  void    fastFindMersA(bool dupIgnore);
  void    fastFindMersB(bool dupIgnore);
};




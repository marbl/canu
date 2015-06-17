
#include "AS_global.H"

#include "ovStore.H"

#include "prefixEditDistance.H"
//#include "Binomial_Bound.H"  //  liboverlap

#include <map>
#include <vector>
#include <algorithm>

using namespace std;


class overlapAlign {
public:
  overlapAlign(bool    partialOverlaps,
              double  maxErate,
              uint32  merSize) {
    _partialOverlaps  = partialOverlaps;
    _maxErate         = maxErate;
    _merSizeInitial   = merSize;
    _editDist         = new prefixEditDistance(_partialOverlaps, _maxErate);
  };

  ~overlapAlign() {
    delete _editDist;
  };



  void  initialize(char *aStr, uint32 aLen, int32 aLo, int32 aHi,
                   char *bStr, uint32 bLen, int32 bLo, int32 bHi) {

    _aStr    = aStr;
    _aLen    = aLen;
    _aLo     = aLo;
    _aHi     = aHi;

    _bStr    = bStr;
    _bLen    = bLen;
    _bLo     = bLo;
    _bHi     = bHi;

    //_editDist doesn't need to be cleared.

    _minDiag = 0;
    _maxDiag = 0;

    _merSize = 0;

    _aMap.clear();
    _bMap.clear();

    _rawhits.clear();
    _hits.clear();
  };


  bool  findMinMaxDiagonal(int32 minLength);
  bool  findSeeds(bool dupIgnore);
  void  findHits(void);
  void  chainHits(void);
  bool  processHits(ovOverlap *result);

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
  uint32              _merSizeInitial;

  //  From the parameters.

  char               *_aStr;
  uint32              _aLen;
  int32               _aLo;
  int32               _aHi;

  char               *_bStr;
  uint32              _bLen;
  int32               _bLo;
  int32               _bHi;

  //  Computed stuff, for each alignment.

  prefixEditDistance *_editDist;

  int32               _minDiag;
  int32               _maxDiag;

  uint32              _merSize;

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




#ifndef MSP_MANAGER_H
#define MSP_MANAGER_H

#include "libbri.H"
#include "exon.H"

struct msp {
  int         len;
  int         pos1;
  int         pos2;
  int         score;
  int         linkingScore;
  int         prev;
};


//
//  How to handle memory allocation?
//
//  Just use an array of msp objects, reallocate
//  when needed.  Allocate a large number of these
//  initially.
//
class mspManager {
public:
  mspManager() {
    _allocMSPs        = 16384;
    _numMSPs          = 0;
    _allMSPs          = new msp [_allocMSPs];

    //  The following four variables are for aborting expensive
    //  polishes -- ones that have proven to be large chunks of
    //  genomic labeled as cDNA, and that have (ESTmapper) signals
    //  across entire scafflds.
    //
    _tooManyMSPs      = false;
    _cDNALength       = 0;
    _mspLimitAbsolute = 0;
    _mspLimitPercent  = 0;

    //  These need to be reset with setParameters.  The code will die
    //  during link() if they are not set.
    //
    _match            = 0;
    _matchdiff        = 0;
    _percentError     = 0.0;
  };

  ~mspManager() {
    delete [] _allMSPs;
  };

  void   setParameters(int match, int mismatch, double percenterror) {
    _match            = match;
    _matchdiff        = match - mismatch;
    _percentError     = percenterror;
  }

  void   setLimits(u32bit a, double p) {
    _mspLimitAbsolute = a;
    _mspLimitPercent  = p;
  };

  void   addMSP(int l, int p1, int p2, int sc) {

    //  Allocate more MSPs, if we need to.
    //
    if (_numMSPs >= _allocMSPs) {
      _allocMSPs *= 2;

      msp *n = new msp [_allocMSPs];

      for (u32bit i=0; i<_numMSPs; i++) {
        n[i].len          = _allMSPs[i].len;
        n[i].pos1         = _allMSPs[i].pos1;
        n[i].pos2         = _allMSPs[i].pos2;
        n[i].score        = _allMSPs[i].score;
        n[i].linkingScore = _allMSPs[i].linkingScore;
        n[i].prev         = _allMSPs[i].prev;
      }

      delete _allMSPs;
      _allMSPs = n;
    }

    _allMSPs[_numMSPs].len          = l;
    _allMSPs[_numMSPs].pos1         = p1;
    _allMSPs[_numMSPs].pos2         = p2;
    _allMSPs[_numMSPs].score        = sc;
    _allMSPs[_numMSPs].linkingScore = 0;
    _allMSPs[_numMSPs].prev         = 0;

    _numMSPs++;
  };

  //  Clear the msp list -- doesn't delete them, just puts them on
  //  a free list.
  //
  void      clear(void) {
    _numMSPs = 0;
  };

  bool          tooManyMSPs(void) {
    return(_tooManyMSPs);
  };
  u32bit        numberOfMSPs(void) {
    return(_numMSPs);
  };
  void          setLength(int l) {
    _cDNALength = l;
  };

  //
  //  What access methods for linking?
  //

  void      sort(void);
  Exon_ptr  link(int weight, int drange,
                 int offset1, int offset2,
                 int flag, int relinkFlag,
                 char *s1, char *s2);

  void      swapMSPpositions(void) {
    register int   p, q, i;

    for (i=_numMSPs-1; i>=0; i--) {
      p = _allMSPs[i].pos1;
      q = _allMSPs[i].pos2;

      _allMSPs[i].pos1 = q;
      _allMSPs[i].pos2 = p;
    }
  };
private:

  //  See the constructor for descriptions

  u32bit         _allocMSPs;
  u32bit         _numMSPs;
  msp           *_allMSPs;

  bool           _tooManyMSPs;
  int            _cDNALength;
  u32bit         _mspLimitAbsolute;
  double         _mspLimitPercent;

  int            _match;
  int            _matchdiff;
  double         _percentError;
};


#endif  //  MSP_MANAGER_H

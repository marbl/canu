#include "exon.H"

#define ABORT_EXPENSIVE 1

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
private:
  unsigned int   _allocMSPs;
  unsigned int   _numMSPs;
  msp           *_allMSPs;
#if ABORT_EXPENSIVE
  bool           _tooManyMSPs;
  int            _cDNALength;
#endif
public:
  mspManager() {
    _allocMSPs   = 16384;
    _numMSPs     = 0;
    _allMSPs     = new msp [_allocMSPs];
#if ABORT_EXPENSIVE
    _tooManyMSPs = false;
    _cDNALength  = 0;
#endif
  };

  ~mspManager() {
    delete [] _allMSPs;
  };

  void   addMSP(int l, int p1, int p2, int sc) {

    //  Allocate more MSPs, if we need to.
    //
    if (_numMSPs >= _allocMSPs) {
      _allocMSPs *= 2;

      msp *n = new msp [_allocMSPs];

      for (unsigned int i=0; i<_numMSPs; i++) {
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

#if ABORT_EXPENSIVE
  bool          tooManyMSPs(void) {
    return(_tooManyMSPs);
  };
  unsigned int  numberOfMSPs(void) {
    return(_numMSPs);
  };
  void          setLength(int l) {
    _cDNALength = l;
  };
#endif

  //
  //  What access methods for linking?
  //

  void      sort(void);
  Exon_ptr  link(int weight, int drange,
                 int offset1, int offset2,
                 int flag, int relinkFlag,
                 unsigned char *s1, unsigned char *s2);

  void      swapMSPpositions(void) {
    register int   p, q, i;

    for (i=_numMSPs-1; i>=0; i--) {
      p = _allMSPs[i].pos1;
      q = _allMSPs[i].pos2;

      _allMSPs[i].pos1 = q;
      _allMSPs[i].pos2 = p;
    }
  };
};

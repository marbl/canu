#ifndef MSP_MANAGER_H
#define MSP_MANAGER_H

#include "bri++.H"
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
  mspManager();
  ~mspManager();

  void   setParameters(int match, int mismatch, double percenterror, int wordextallow) {
    _match            = match;
    _mismatch         = mismatch;
    _matchdiff        = match - mismatch;
    _percentError     = percenterror;
    _wordExtAllow     = wordextallow;
    //int MATCH    = globalParams->_match;
    //int MISMATCH = globalParams->_mismatch;
  };

  void   setLimits(u32bit a, double p) {
    _mspLimitAbsolute = a;
    _mspLimitPercent  = p;
  };

  bool          tooManyMSPs(void)  { return(_tooManyMSPs); };
  u32bit        numberOfMSPs(void) { return(_numMSPs);     };

  void          setLength(int l)   { _cDNALength = l;      };

  void          clear(void);


  //  add an extended MSP to the list
  //
  void   addMSP(int l, int p1, int p2, int sc);

  //  add a single unextended hit to the list.
  //  this will do extensions if we haven't already extended through it.
  //
  void  addHit(char *s1, char *s2,
               int   l1, int   l2,
               int   p1, int   p2,
               int   W,
               int   K);


  Exon_ptr  doLinking(int    weight,
                      int    drange,
                      int    offset1,
                      int    offset2,
                      int    flag,
                      int    relinkFlag,
                      char  *s1,
                      char  *s2);

private:
  bool           _sorted;

  u32bit         _allocMSPs;
  u32bit         _numMSPs;
  msp           *_allMSPs;

  bool           _tooManyMSPs;
  int            _cDNALength;

  u32bit         _mspLimitAbsolute;
  double         _mspLimitPercent;

  int            _match;
  int            _mismatch;
  int            _matchdiff;
  double         _percentError;
  int            _wordExtAllow;

  int            _diagMax;
  int           *_diagExt;
};




inline
void
mspManager::addMSP(int l, int p1, int p2, int sc) {

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

#if 0
  fprintf(stdout, "ADDMSP: p1=%8d p2=%8d l=%8d sc=%8d\n",
          p1, p2, l, sc);
#endif

  _allMSPs[_numMSPs].len          = l;
  _allMSPs[_numMSPs].pos1         = p1;
  _allMSPs[_numMSPs].pos2         = p2;
  _allMSPs[_numMSPs].score        = sc;
  _allMSPs[_numMSPs].linkingScore = 0;
  _allMSPs[_numMSPs].prev         = 0;

  _numMSPs++;

  _sorted = false;
}


#endif  //  MSP_MANAGER_H

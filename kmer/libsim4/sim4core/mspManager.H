#ifndef MSP_MANAGER_H
#define MSP_MANAGER_H

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "util++.H"
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
    _matchdiff        = match - mismatch;
    _percentError     = percenterror;
    _wordExtAllow     = wordextallow;
  };

  void   setLimits(u32bit a, double p) {
    _mspLimitAbsolute = a;
    _mspLimitPercent  = p;
  };

  void   setScoreThreshold(int K, int interspecies);

  bool          tooManyMSPs(void)  { return(_tooManyMSPs); };
  u32bit        numberOfMSPs(void) { return(_numMSPs);     };

  void          setLength(int l)   {
    _cDNALength = l;
  };

  void          clear(void) { _numMSPs = 0; };

  void          clearDiagonal(int genlen, int estlen) {

    //  XXX: These aren't always the real EST and GENOMIC lengths.  If
    //  we are working in a subsequence of the whole sequence they
    //  will be the length of the subsequence.
    //
    _GENlen = genlen;
    _ESTlen = estlen;

    //  Allocate more space, if needed.
    //
    if (_GENlen + _ESTlen + 1 > _diagMax) {
      delete _diagExt;

      //  Allocate space for the list of extension ends.  Each diagonal
      //  remembers the lowest position that it has been extended to.  This
      //  lets us throw out new hits without extending or merging in new
      //  extensions.  Assumes that hits are added in order.
      //
      _diagMax = _GENlen + _ESTlen + 1;
      _diagExt = new int [_diagMax];
    }

    //  Using the obvious for loop for this hurts.  Don't do it.
    //
    bzero(_diagExt, sizeof(int) * (_GENlen + _ESTlen + 1));
  };

  //  add an extended MSP to the list
  //
  void   addMSP(int l, int p1, int p2, int sc);

  //  add a single unextended hit to the list.
  //  this will do extensions if we haven't already extended through it.
  //
  void  addHit(char *genSeq, char *estSeq,
               int   genLen, int   estLen,
               int   genPos, int   estPos,
               int   W) {

#if 0
    //  We'd like to tie this into DEBUG_EXTENSION, but I want
    //  to keep those defines in the source file.  Oh well.
    //
    fprintf(stderr, "mspManager::addHit()-- adding hit from GEN %d to %d and EST %d to %d (length = %d)\n", 
            genPos-W, genPos, estPos-W, estPos, W);
#endif

    if (_diagExt[estLen + genPos - estPos - 1] <= genPos)
      addHit_(genSeq, estSeq, genLen, estLen, genPos, estPos, W);
  };


  Exon     *doLinking(int    weight,
                      int    drange,
                      int    offset1,
                      int    offset2,
                      int    flag,
                      int    relinkFlag,
                      char  *s1,
                      char  *s2);

private:
  void  addHit_(char *genSeq, char *estSeq,
                int   genLen, int   estLen,
                int   genPos, int   estPos,
                int   W);

  bool           _sorted;

  int            _ESTlen;
  int            _GENlen;

  u32bit         _allocMSPs;
  u32bit         _numMSPs;
  msp           *_allMSPs;

  bool           _tooManyMSPs;
  int            _cDNALength;

  double         _mspLimitPercent;
  u32bit         _mspLimitAbsolute;

  int            _match;
  int            _matchdiff;

#ifdef __APPLE__
  int            _pad1;
#endif

  double         _percentError;
  int            _wordExtAllow;

  int            _minMSPScore;

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

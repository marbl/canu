#include "seqCache.H"
#include "seqFactory.H"
#include "alphabet.h"

#undef DEBUG


seqCache::seqCache(const char *filename, u32bit cachesize, bool verbose) {

  _fb                  = openSeqFile(filename);
  _idToGetNext         = 0;
  _allSequencesLoaded  = false;
  _reportLoading       = verbose;

  _cacheMap            = 0L;
  _cacheSize           = 0;
  _cacheNext           = 0;
  _cache               = 0L;

  setCacheSize(cachesize);
}



seqCache::~seqCache() {
  flushCache();
  delete    _fb;
  delete [] _cacheMap;
  delete [] _cache;
}



u32bit
seqCache::getSequenceIID(char *name) {
  u32bit iid = _fb->find(name);

  //  Nothing?  If the name is all integers, return that, otherwise,
  //  fail.
  //
  if (iid == ~u32bitZERO) {
    bool  isInt = true;
    char *x = name;

    while (*x)
      if ((*x < '0') || ('9' < *x))
        isInt = false;

    if (isInt)
      iid = strtou32bit(name, 0L);
  }

#ifdef DEBUG
  fprintf(stderr, "seqCache::getSequenceIID()-- '%s' -> "u32bitFMT"\n", name, iid);
#endif

  return(iid);
}



seqInCore *
seqCache::getSequenceInCore(u32bit iid) {

#ifdef DEBUG
  fprintf(stderr, "seqCache::getSequenceInCore(iid)-- "u32bitFMT"\n", iid);
#endif

  if (iid >= _fb->getNumberOfSequences()) {
#ifdef DEBUG
    fprintf(stderr, "seqCache::getSequenceInCore(iid)-- iid="u32bitFMT" >= numSeq="u32bitFMT"\n", iid, _fb->getNumberOfSequences());
#endif
    return(0L);
  }

  if (_allSequencesLoaded == true)
    return(_cache[iid]);

  if ((_cacheSize > 0) && (_cacheMap[iid] != ~u32bitZERO))
    return(_cache[_cacheMap[iid]]);

  u32bit  hLen=0, hMax=0, sLen=0, sMax=0;
  char   *h=0L, *s=0L;

  if (_fb->getSequence(iid, h, hLen, hMax, s, sLen, sMax) == false) {
#ifdef DEBUG
    fprintf(stderr, "seqCache::getSequenceInCore(iid)-- failed getSequence().");
#endif
    return(0L);
  }

  seqInCore *sc = new seqInCore(iid, h, hLen, s, sLen);

  //  Remove and old cached sequence

  if (_cache) {
    if (_cache[_cacheNext]) {
      _cacheMap[_cache[_cacheNext]->getIID()] = ~u32bitZERO;
      delete _cache[_cacheNext];
      _cache[_cacheNext] = 0L;
    }

    //  Store the new one in the cache

    _cache[_cacheNext] = sc;
    _cacheMap[iid]     = _cacheNext;

    _cacheNext = (_cacheNext + 1) % _cacheSize;
  }

  return(sc);
}



void
seqCache::setCacheSize(u32bit cachesize) {

  flushCache();

  if (cachesize == 0) {
    _cacheMap   = 0L;
    _cacheSize  = 0;
    _cacheNext  = 0;
    _cache      = 0L;
    return;
  }

  _cacheMap   = new u32bit [_fb->getNumberOfSequences()];
  _cacheSize  = cachesize;
  _cacheNext  = 0;
  _cache      = new seqInCore * [_cacheSize];

  for (u32bit i=_fb->getNumberOfSequences()-1; i--; )
    _cacheMap[i] = ~u32bitZERO;

  for (u32bit i=_cacheSize-1; i--; )
    _cache[i] = 0L;
}



void
seqCache::loadAllSequences(void) {

  if (_allSequencesLoaded)
    return;

  flushCache();

  delete [] _cacheMap;
  delete [] _cache;

  _cacheMap   = 0L;
  _cacheSize  = _fb->getNumberOfSequences();
  _cacheNext  = 0;
  _cache      = new seqInCore * [_cacheSize];


  for (u32bit iid=0; iid<_cacheSize; iid++) {
    u32bit  hLen=0, hMax=0, sLen=0, sMax=0;
    char   *h=0L, *s=0L;

    if (_fb->getSequence(iid, h, hLen, hMax, s, sLen, sMax) == false)
      fprintf(stderr, "seqCache::loadAllSequences()-- Failed to load iid "u32bitFMT".\n",
              iid), exit(1);

    _cache[iid] = new seqInCore(iid, h, hLen, s, sLen);
  }

  _allSequencesLoaded = true;
}

void
seqCache::flushCache(void) {
  u32bit ns = _fb->getNumberOfSequences();

  if (_cacheMap)
    for (u32bit i=0; i<ns; i++)
      _cacheMap[i] = ~u32bitZERO;

  if (_cache)
    for (u32bit i=0; i<_cacheSize; i++) {
      delete _cache[i];
      _cache[i] = 0L;
    }
}

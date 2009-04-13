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
  u32bit iid = ~u32bitZERO;

  //  If the name is all integers, AND below the number of sequences
  //  we have, return that, otherwise, look it up.
  //
  bool  isInt = true;
  char *x     = name;

  while (*x) {
    if ((*x < '0') || ('9' < *x))
      isInt = false;
    x++;
  }

  if (isInt)
    iid = strtou32bit(name, 0L);

  if (iid >= _fb->getNumberOfSequences())
    iid = _fb->find(name);

#ifdef DEBUG
  fprintf(stderr, "seqCache::getSequenceIID()-- '%s' -> "u32bitFMT"\n", name, iid);
#endif

  return(iid);
}



seqInCore *
seqCache::getSequenceInCore(u32bit iid) {
  u32bit       cacheID = ~u32bitZERO;
  seqInCore   *retSeq  = 0L;

  if (iid >= _fb->getNumberOfSequences())
    return(0L);


  if (_allSequencesLoaded == true) {
    cacheID = iid;

  } else if ((_cacheSize > 0) && (_cacheMap[iid] != ~u32bitZERO)) {
    cacheID = _cacheMap[iid];

  } else {
    u32bit  hLen=0, hMax=0, sLen=0, sMax=0;
    char   *h=0L, *s=0L;

    if (_fb->getSequence(iid, h, hLen, hMax, s, sLen, sMax) == false)
      return(0L);

    retSeq = new seqInCore(iid, h, hLen, s, sLen, true);

    //  Remove any old cached sequence, then store the one we just made

    if (_cache) {
      if (_cache[_cacheNext]) {
        _cacheMap[_cache[_cacheNext]->getIID()] = ~u32bitZERO;
        delete _cache[_cacheNext];
      }

      _cache[_cacheNext] = retSeq;
      _cacheMap[iid]     = _cacheNext;

      cacheID = _cacheNext;
      retSeq  = 0L;

      _cacheNext = (_cacheNext + 1) % _cacheSize;
    }
  }

  //  If no retSeq set, make a copy of the one we have in the cache.

  if ((retSeq == 0L) && (cacheID != ~u32bitZERO))
    retSeq  = new seqInCore(iid,
                            _cache[cacheID]->header(),   _cache[cacheID]->headerLength(),
                            _cache[cacheID]->sequence(), _cache[cacheID]->sequenceLength(),
                            false);

  return(retSeq);
}



void
seqCache::setCacheSize(u32bit cachesize) {
  u32bit ns = _fb->getNumberOfSequences();

  flushCache();

  if (cachesize == 0) {
    _cacheMap   = 0L;
    _cacheSize  = 0;
    _cacheNext  = 0;
    _cache      = 0L;
    return;
  }

  _cacheMap   = new u32bit [ns];
  _cacheSize  = cachesize;
  _cacheNext  = 0;
  _cache      = new seqInCore * [_cacheSize];

  for (u32bit i=0; i<ns; i++)
    _cacheMap[i] = ~u32bitZERO;

  for (u32bit i=0; i<_cacheSize; i++)
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

    _cache[iid] = new seqInCore(iid, h, hLen, s, sLen, true);
  }

  _allSequencesLoaded = true;
}

void
seqCache::flushCache(void) {

  if (_fb == 0L)
    return;

  if (_cacheMap) {
    u32bit ns = _fb->getNumberOfSequences();
    for (u32bit i=0; i<ns; i++)
      _cacheMap[i] = ~u32bitZERO;
  }

  if (_cache)
    for (u32bit i=0; i<_cacheSize; i++) {
      delete _cache[i];
      _cache[i] = 0L;
    }
}

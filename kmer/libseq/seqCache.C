#include "seqCache.H"
#include "seqFactory.H"
#include "alphabet.h"

#undef DEBUG


seqCache::seqCache(const char *filename, uint32 cachesize, bool verbose) {

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



uint32
seqCache::getSequenceIID(char *name) {
  uint32 iid = ~uint32ZERO;

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
    iid = strtouint32(name, 0L);

  if (iid >= _fb->getNumberOfSequences())
    iid = _fb->find(name);

#ifdef DEBUG
  fprintf(stderr, "seqCache::getSequenceIID()-- '%s' -> "uint32FMT"\n", name, iid);
#endif

  return(iid);
}



seqInCore *
seqCache::getSequenceInCore(uint32 iid) {
  uint32       cacheID = ~uint32ZERO;
  seqInCore   *retSeq  = 0L;

  if ((_fb->randomAccessSupported() == true) &&
      (iid >= _fb->getNumberOfSequences()))
    return(0L);

  if (_allSequencesLoaded == true) {
    cacheID = iid;

  } else if ((_cacheSize > 0) && (_cacheMap[iid] != ~uint32ZERO)) {
    cacheID = _cacheMap[iid];

  } else {
    uint32  hLen=0, hMax=0, sLen=0, sMax=0;
    char   *h=0L, *s=0L;

    if (_fb->getSequence(iid, h, hLen, hMax, s, sLen, sMax) == false)
      return(0L);

    retSeq = new seqInCore(iid, h, hLen, s, sLen, true);

    //  Remove any old cached sequence, then store the one we just made

    if (_cache) {
      if (_cache[_cacheNext]) {
        _cacheMap[_cache[_cacheNext]->getIID()] = ~uint32ZERO;
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

  if ((retSeq == 0L) && (cacheID != ~uint32ZERO))
    retSeq  = new seqInCore(iid,
                            _cache[cacheID]->header(),   _cache[cacheID]->headerLength(),
                            _cache[cacheID]->sequence(), _cache[cacheID]->sequenceLength(),
                            false);

  return(retSeq);
}



void
seqCache::setCacheSize(uint32 cachesize) {
  uint32 ns = _fb->getNumberOfSequences();

  flushCache();

  if (cachesize == 0) {
    _cacheMap   = 0L;
    _cacheSize  = 0;
    _cacheNext  = 0;
    _cache      = 0L;
    return;
  }

  _cacheMap   = new uint32 [ns];
  _cacheSize  = cachesize;
  _cacheNext  = 0;
  _cache      = new seqInCore * [_cacheSize];

  for (uint32 i=0; i<ns; i++)
    _cacheMap[i] = ~uint32ZERO;

  for (uint32 i=0; i<_cacheSize; i++)
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


  for (uint32 iid=0; iid<_cacheSize; iid++) {
    uint32  hLen=0, hMax=0, sLen=0, sMax=0;
    char   *h=0L, *s=0L;

    if (_fb->getSequence(iid, h, hLen, hMax, s, sLen, sMax) == false)
      fprintf(stderr, "seqCache::loadAllSequences()-- Failed to load iid "uint32FMT".\n",
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
    uint32 ns = _fb->getNumberOfSequences();
    for (uint32 i=0; i<ns; i++)
      _cacheMap[i] = ~uint32ZERO;
  }

  if (_cache)
    for (uint32 i=0; i<_cacheSize; i++) {
      delete _cache[i];
      _cache[i] = 0L;
    }
}

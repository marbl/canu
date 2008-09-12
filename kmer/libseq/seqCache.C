#include "seqCache.H"
#include "seqFactory.H"
#include "alphabet.h"

#undef DEBUG

seqOnDisk::seqOnDisk(char const *filename,
                     u32bit iid,
                     u64bit hdrstart, u32bit hdrlen,
                     u64bit seqstart, u32bit seqlen) {

  _idx               = iid;

  _headerLength      = hdrlen;
  _headerStart       = hdrstart;

  _sequenceLength    = seqlen;
  _sequenceStart     = seqstart;

  _rb        = new readBuffer(filename, MIN(1024 * 1024, _headerLength + _sequenceLength + 16));

  _rb->seek(_headerStart);

  _header = new char [_headerLength + 1];
  _rb->read(_header, _headerLength);
  _header[_headerLength] = 0;

  _sequence = 0L;

  _rb->seek(_sequenceStart);

  _sequencePosition  = 0;
}



seqOnDisk::seqOnDisk(u32bit iid,
            char *hdr, u32bit hdrlen,
                     char *seq, u32bit seqlen) {
  _rb       = 0L;;
  _idx              = iid;
  _headerLength     = hdrlen;
  _sequenceLength   = seqlen;
  _headerStart      = ~u64bitZERO;
  _sequenceStart    = ~u64bitZERO;
  _header           = hdr;
  _sequence         = seq;
  _sequencePosition = 0;
}



seqOnDisk::~seqOnDisk() {
  delete    _rb;
  delete [] _header;
  delete [] _sequence;
}



char
seqOnDisk::read(void) {

  if (_sequencePosition >= _sequenceLength)
    return(0);

  if (_sequence)
    return(_sequence[_sequencePosition++]);

  //  Assumptions; there is another non-space letter out there,
  //  otherwise, _sequencePosition would equal _sequenceLength.

  char   x = _rb->read();

  while (whitespaceSymbol[x]) {
    if (_rb->eof()) {
      fprintf(stderr, "seqOnDisk::read()--  WARNING hit unexpected eof.\n");
      return(0);
    }

    x = _rb->read();
  }

  _sequencePosition++;

  return(x);
}




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
  u32bit  iid = ~u32bitZERO;

#ifdef DEBUG
  fprintf(stderr, "seqCache::getSequenceIID()-- '%s'\n", name);
#endif

#warning mostly unimplemented
  iid = strtou32bit(name, 0L);

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



seqOnDisk *
seqCache::getSequenceOnDisk(u32bit iid) {

#ifdef DEBUG
  fprintf(stderr, "seqCache::getSequenceOnDisk()-- "u32bitFMT"\n", iid);
#endif

  if (iid >= _fb->getNumberOfSequences())
    return(0L);

#warning unimplemented
  fprintf(stderr, "seqCache::getSequenceOnDisk()--  Not implemented.\n");
  exit(1);
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

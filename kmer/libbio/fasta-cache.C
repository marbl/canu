#include "fasta-cache.H"


FastACache::FastACache(char *filename, u32bit cachesize, bool loadall) {
  _fastawrapper = new FastAWrapper(filename);
  _fastawrapper->openIndex();

  if (loadall == false) {
    _allSequencesLoaded = false;
  
    _cacheMap  = new u32bit [_fastawrapper->getNumberOfSequences()];

    for (u32bit i=_fastawrapper->getNumberOfSequences(); i-- > 0; )
      _cacheMap[i] = ~u32bitZERO;

    _cacheSize = cachesize;
    _cacheNext = 0;
    _cache     = new FastASequenceInCore* [_cacheSize];

    for (u32bit i=0; i<_cacheSize; i++)
      _cache[i] = 0L;
  } else {
    _allSequencesLoaded = true;

    _cacheMap  = 0L;
    _cacheSize = _fastawrapper->getNumberOfSequences();
    _cacheNext = 0;
    _cache     = new FastASequenceInCore* [_cacheSize];

    fprintf(stderr, "Loading %u sequences from '%s'\n",
            _cacheSize,
            filename);

    for (u32bit i=0; i<_cacheSize; i++)
      _cache[i] = _fastawrapper->getSequence();
  }
}


FastACache::~FastACache() {
  for (u32bit i=0; i<_cacheSize; i++)
    delete _cache[i];
  delete [] _cache;
  delete [] _cacheMap;
  delete    _fastawrapper;
}



FastASequenceInCore*
FastACache::getSequence(u32bit iid)  {

  if (_allSequencesLoaded) {
    if (iid < _cacheSize) {
      return(_cache[iid]);
    } else {
      fprintf(stderr, "ERROR: FastACache of '%s' was asked for iid=%u, but only %u available.\n",
              _fastawrapper->getSourceName(), iid, _cacheSize);
      return(0L);
    }
  }

  //  Not all sequences loaded.

  if (_cacheMap[iid] != ~u32bitZERO)
    return(_cache[_cacheMap[iid]]);

  //  The sequence we are looking for isn't loaded.  If there
  //  isn't space in the cache, make space.

  //fprintf(stderr, "FastACache::getSequence()-- "u32bitFMT" isn't loaded -- loading.\n", iid);

  if (_cache[_cacheNext]) {
    _cacheMap[ _cache[_cacheNext]->getIID() ] = ~u32bitZERO;
    delete _cache[_cacheNext];
  }

  //  Load the sequence into the cache

  _fastawrapper->find(iid);
  _cache[_cacheNext] = _fastawrapper->getSequence();
  _cacheMap[iid] = _cacheNext;
  _cacheNext++;

  if (_cacheNext >= _cacheSize)
    _cacheNext = 0;

  //  Finally, return the sequence.

  return(_cache[_cacheNext-1]);
}


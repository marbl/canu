#ifndef FASTA_CACHE_H
#define FASTA_CACHE_H

#include "bio++.H"

//  Keeps a cache of the last N sequences accessed.  If N==0, then the
//  cache is unlimited in size.

class FastACache {
public:
  FastACache(char *filename, u32bit cachesize=32, bool loadall=false, bool report=false);
  ~FastACache();

  //  Returns sequence iid.  If the sequence is not in the cache,
  //  the sequence is loaded.
  //
  FastASequenceInCore    *getSequence(u32bit iid);

  //  Returns a pointer to the fasta wrapper we are using.
  //
  FastAWrapper           *fasta(void) {
    return(_fastawrapper);
  };

private:
  FastAWrapper          *_fastawrapper;

  bool                   _allSequencesLoaded;
  bool                   _reportLoading;

  //  cacheSize:  The size of the cache.  If non-zero, this is the size
  //  of _cacheMap and _seqs.
  //
  u32bit                *_cacheMap;
  u32bit                 _cacheSize;
  u32bit                 _cacheNext;
  FastASequenceInCore  **_cache;
};

#endif  //  FASTA_CACHE_H

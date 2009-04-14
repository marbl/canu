#ifndef SEQCACHE_H
#define SEQCACHE_H

#include "util++.H"
#include "seqFile.H"


class seqInCore {
private:
  seqInCore(u32bit iid, char *hdr, u32bit hdrlen, char *seq, u32bit seqlen, bool deletable) {
    _idx       = iid;

    _deletable = deletable;

    _headerLen = hdrlen;
    _header    = hdr;

    _seqLen = seqlen;
    _seq    = seq;
  };

  friend class seqCache;

public:
  ~seqInCore() {
    if (_deletable) {
      delete [] _header;  _header = 0L;
      delete [] _seq;     _seq    = 0L;
    }
  };

  char           *header(void) const         { return(_header); };
  u32bit          headerLength(void) const   { return(_headerLen); };

  char           *sequence(void) const       { return(_seq); };
  u32bit          sequenceLength(void) const { return(_seqLen); };

  u32bit          getIID(void) const         { return(_idx); };

  //  Used only by searchGENOME (as far as I know)
  seqInCore      *copy(void) {
    char *h = new char [_headerLen + 1];
    char *s = new char [_seqLen    + 1];

    memcpy(h, _header, _headerLen + 1);
    memcpy(s, _seq,    _seqLen    + 1);

    return(new seqInCore(_idx, h, _headerLen, s, _seqLen, true));
  };

private:
  u32bit         _idx;

  bool           _deletable;

  u32bit         _headerLen;
  char          *_header;

  u32bit         _seqLen;
  char          *_seq;
};



class seqCache {
public:
  seqCache(const char *filename, u32bit cachesize=0, bool verbose=false);
  ~seqCache();

  //  Returns IID for a name, either the first word on the defline, or
  //  the ascii IID.
  u32bit                  getSequenceIID(char *name);

  seqInCore              *getSequenceInCore(u32bit  iid);
  seqInCore              *getSequenceInCore(char   *name) { return(getSequenceInCore(getSequenceIID(name))); };
  seqInCore              *getSequenceInCore(void)         { return(getSequenceInCore(_idToGetNext++)); };

  const char             *getSourceName(void)    { return(_fb->getSourceName()); };
  const char             *getFileTypeName(void)  { return(_fb->getFileTypeName()); };

  u32bit                  getNumberOfSequences(void) { return(_fb->getNumberOfSequences()); };

  u32bit                  getSequenceLength(u32bit iid) { return(_fb->getSequenceLength(iid)); };

  void                    setCacheSize(u32bit cachesize);

  void                    loadAllSequences(void);
  void                    flushCache(void);

private:
  seqFile               *_fb;
  u32bit                 _idToGetNext;

  bool                   _allSequencesLoaded;
  bool                   _reportLoading;

  u32bit                *_cacheMap;   //  Maps ID to cache entry
  u32bit                 _cacheSize;  //  Size of cache
  u32bit                 _cacheNext;  //  Next cache spot to use
  seqInCore            **_cache;      //  Cache of sequences
};


#endif  //  SEQCACHE_H

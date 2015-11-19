
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef SEQCACHE_H
#define SEQCACHE_H

#include "seqFile.H"


class seqInCore {
private:
  seqInCore(uint32 iid, char *hdr, uint32 hdrlen, char *seq, uint32 seqlen, bool deletable) {
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
  uint32          headerLength(void) const   { return(_headerLen); };

  char           *sequence(void) const       { return(_seq); };
  uint32          sequenceLength(void) const { return(_seqLen); };

  uint32          getIID(void) const         { return(_idx); };

  //  Used only by searchGENOME (as far as I know)
  seqInCore      *copy(void) {
    char *h = new char [_headerLen + 1];
    char *s = new char [_seqLen    + 1];

    memcpy(h, _header, _headerLen + 1);
    memcpy(s, _seq,    _seqLen    + 1);

    return(new seqInCore(_idx, h, _headerLen, s, _seqLen, true));
  };

private:
  uint32         _idx;

  bool           _deletable;

  uint32         _headerLen;
  char          *_header;

  uint32         _seqLen;
  char          *_seq;
};



class seqCache {
public:
  seqCache(const char *filename, uint32 cachesize=0, bool verbose=false);
  ~seqCache();

  //  Returns IID for a name, either the first word on the defline, or
  //  the ascii IID.
  uint32                  getSequenceIID(char *name);

  seqInCore              *getSequenceInCore(uint32  iid);
  seqInCore              *getSequenceInCore(char   *name) { return(getSequenceInCore(getSequenceIID(name))); };
  seqInCore              *getSequenceInCore(void)         { return(getSequenceInCore(_idToGetNext++)); };

  const char             *getSourceName(void)    { return(_fb->getSourceName()); };
  const char             *getFileTypeName(void)  { return(_fb->getFileTypeName()); };

  bool                    randomAccessSupported(void) { return(_fb->randomAccessSupported()); };

  uint32                  getNumberOfSequences(void) { return(_fb->getNumberOfSequences()); };

  uint32                  getSequenceLength(uint32 iid) { return(_fb->getSequenceLength(iid)); };

  void                    setCacheSize(uint32 cachesize);

  void                    loadAllSequences(void);
  void                    flushCache(void);

private:
  seqFile               *_fb;
  uint32                 _idToGetNext;

  bool                   _allSequencesLoaded;
  bool                   _reportLoading;

  uint32                *_cacheMap;   //  Maps ID to cache entry
  uint32                 _cacheSize;  //  Size of cache
  uint32                 _cacheNext;  //  Next cache spot to use
  seqInCore            **_cache;      //  Cache of sequences
};


#endif  //  SEQCACHE_H

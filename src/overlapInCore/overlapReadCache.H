
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"
#include "sqStore.H"
#include "ovStore.H"
#include "tgStore.H"

#include <set>
using namespace std;

class overlapReadCache {
public:
  overlapReadCache(sqStore *seqStore_, uint64 memLimit);
  ~overlapReadCache();

private:
  void         loadRead(uint32 id);
  void         loadReads(set<uint32> reads);
  void         markForLoading(set<uint32> &reads, uint32 id);

public:
  void         loadReads(ovOverlap *ovl, uint32 nOvl);
  void         loadReads(tgTig *tig);

  void         purgeReads(void);

  char        *getRead(uint32 id) {
    assert(readLen[id] > 0);
    return(readSeqFwd[id]);
  };

  uint32       getLength(uint32 id) {
    assert(readLen[id] > 0);
    return(readLen[id]);
  };

private:
  sqStore     *seqStore;
  uint32       nReads;

  uint32      *readAge;
  uint32      *readLen;
  char       **readSeqFwd;

  sqRead       read;

  uint64       memoryLimit;
};




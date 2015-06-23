
#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"
#include "tgStore.H"

class overlapReadCache {
public:
  overlapReadCache(gkStore *gkpStore_, uint64 memLimit);
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
  gkStore     *gkpStore;
  uint32       nReads;

  uint32      *readAge;
  uint32      *readLen;
  char       **readSeqFwd;
  //char       **readSeqRev;  //  Save it, or recompute?

  gkReadData   readdata;

  uint64       memoryLimit;
};




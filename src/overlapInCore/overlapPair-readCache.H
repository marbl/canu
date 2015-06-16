

#include "AS_global.H"
#include "gkStore.H"


class readCache {
public:
  readCache(gkStore *gkpStore_, uint64 memLimit) {
    gkpStore    = gkpStore_;
    nReads      = gkpStore->gkStore_getNumReads();

    readAge     = new uint32 [nReads + 1];
    readLen     = new uint32 [nReads + 1];

    memset(readAge, 0, sizeof(uint32) * (nReads + 1));
    memset(readLen, 0, sizeof(uint32) * (nReads + 1));

    readSeqFwd  = new char * [nReads + 1];
    //readSeqRev  = new char * [nReads + 1];

    memset(readSeqFwd, 0, sizeof(char *) * (nReads + 1));
    //memset(readSeqRev, 0, sizeof(char *) * (nReads + 1));

    memoryLimit = memLimit * 1024 * 1024 * 1024;
  };

  ~readCache() {
    delete [] readAge;
    delete [] readLen;

    for (uint32 rr=0; rr<=nReads; rr++) {
      delete [] readSeqFwd[rr];
      //delete [] readSeqRev[rr];
    }

    delete [] readSeqFwd;
    //delete [] readSeqRev;
  };


  void         loadRead(uint32 id) {
    gkRead *read = gkpStore->gkStore_getRead(id);
    gkpStore->gkStore_loadReadData(read, &readdata);

    readLen[id] = read->gkRead_sequenceLength();

    readSeqFwd[id] = new char [readLen[id] + 1];
    //readSeqRev[id] = new char [readLen[id] + 1];

    memcpy(readSeqFwd[id], readdata.gkReadData_getSequence(), sizeof(char) * readLen[id]);

    readSeqFwd[id][readLen[id]] = 0;
  };


  void         loadReads(ovsOverlap *ovl, uint32 nOvl) {
    uint32  nLoaded  = 0;
    uint32  nUpdated = 0;
    uint64  memUsed  = 0;

    //fprintf(stderr, "loadReads()--  begins\n");

    for (uint32 oo=0; oo<nOvl; oo++) {
      uint32  aid = ovl[oo].a_iid;
      uint32  bid = ovl[oo].b_iid;

      if (readLen[aid] == 0) {
        nLoaded++;
        loadRead(aid);
      } else {
        nUpdated++;
      }

      if (readLen[bid] == 0) {
        nLoaded++;
        loadRead(bid);
      } else {
        nUpdated++;
      }

      readAge[aid] = 0;
      readAge[bid] = 0;
    }

    //  Age all the reads (also count the space used)

    for (uint32 id=0; id<nReads; id++) {
      readAge[id]++;
      memUsed += readLen[id];
    }

    if (nLoaded + nUpdated > 0)
      fprintf(stderr, "loadReads()--  loaded %u updated %u -- %f.3 GB used\n",
              nLoaded, nUpdated, memUsed / 1024.0 / 1024.0 / 1024.0);
  };


  void         purgeReads(void) {
    uint32  maxAge     = 0;
    uint64  memoryUsed = 0;

    //  Find maxAge, and sum memory used

    for (uint32 rr=0; rr<=nReads; rr++) {
      if (maxAge < readAge[rr])
        maxAge = readAge[rr];

      memoryUsed += readLen[rr];
    }

    //  Purge oldest until memory is below watermark

    while (memoryLimit < memoryUsed) {
      for (uint32 rr=0; rr<=nReads; rr++) {
        if (maxAge == readAge[rr]) {
          memoryUsed -= readLen[rr];

          delete [] readSeqFwd[rr];  readSeqFwd[rr] = NULL;
          //delete [] readSeqRev[rr];  readSeqRev[rr] = NULL;

          readLen[rr] = 0;
          readAge[rr] = 0;
        }
      }

      maxAge--;
    }
  };


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




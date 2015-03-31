
const char *mainid = "$Id:  $";

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include "overlapInCore.H"

#include "AS_UTL_reverseComplement.H"

#include "Binomial_Bound.H"  //  liboverlap

#include "kMer.H"
#include "merStream.H"

#include <map>
#include <vector>
#include <algorithm>

using namespace std;



class exactMatch {
public:
  exactMatch(int32 a, int32 b, int32 l) {
    aBgn = a;
    bBgn = b;
    tLen = l;
  };

  int32  aBgn;  //  Signed to allow for easy compute of diagonal
  int32  bBgn;
  int32  tLen;

  bool operator<(exactMatch const that) const {
    if (tLen == that.tLen)
      return(aBgn < that.aBgn);

    return(tLen > that.tLen);
  };
};



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

    fprintf(stderr, "loadReads()--\n");

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








Overlap_t
Extend_Alignment(Match_Node_t  *Match,
                 char          *S,     int32   S_Len,
                 char          *T,     int32   T_Len,
                 int32         &S_Lo,  int32   &S_Hi,
                 int32         &T_Lo,  int32   &T_Hi,
                 int32         &Errors,
                 Work_Area_t   *WA);




//uint32        *overlapsLen  = NULL;
//ovsOverlap    *overlaps     = NULL;

readCache        *cache        = NULL;

oicParameters     G;

pthread_mutex_t   balanceMutex;
uint32            batchBgnID = 0;
uint32            batchEndID = 0;

#define THREAD_SIZE  512
#define BATCH_SIZE   512 * 256


bool
getRange(uint32 &bgnID, uint32 &endID) {

  pthread_mutex_lock(&balanceMutex);

  bgnID       = batchBgnID;
  batchBgnID += THREAD_SIZE;         //  Supposed to overflow.
  endID       = batchBgnID;

  if (endID > batchEndID)
    endID = batchEndID;

  pthread_mutex_unlock(&balanceMutex);

  //  If we're out of overlaps, batchBgnID is more than batchEndID (from the last call to this
  //  function), which makes bgnID > endID (in this call).

  return(bgnID < endID);
}



void *
recomputeOverlaps(void *ptr) {
  Work_Area_t  *WA = (Work_Area_t *)ptr;

  char         *bRev = new char [AS_MAX_READLEN];

  uint32        bgnID = 0;
  uint32        endID = 0;

  while (getRange(bgnID, endID)) {
    fprintf(stderr, "Thread %u computes range %u - %u\n", WA->thread_id, bgnID, endID);

    for (uint32 oo=bgnID; oo<endID; oo++) {
      ovsOverlap  *ovl = WA->overlaps + oo;

      //  Load A.

      char   *aStr = cache->getRead  (ovl->a_iid);
      uint32  aLen = cache->getLength(ovl->a_iid);

      int32   aLo = ovl->a_bgn(WA->gkpStore);
      int32   aHi = ovl->a_end(WA->gkpStore);

      assert(aLo < aHi);

      //  Load B.

      char   *bStr = cache->getRead  (ovl->b_iid);
      uint32  bLen = cache->getLength(ovl->b_iid);

      int32   bLo = ovl->b_bgn(WA->gkpStore);
      int32   bHi = ovl->b_end(WA->gkpStore);

      if (ovl->flipped() == true) {
        memcpy(bRev, bStr, sizeof(char) * (bLen + 1));

        reverseComplementSequence(bRev, bLen);

        bStr = bRev;

        bLo = bLen - ovl->b_bgn(WA->gkpStore);
        bHi = bLen - ovl->b_end(WA->gkpStore);
      }

      assert(bLo < bHi);

      //  Set the min/max diagonal we will accept seeds for.  It's just the min/max diagonal for the
      //  two endpoints extended by half the erate times the align length.

      int32  minDiag = 0;
      int32  maxDiag = 0;

      {
        int32  aALen = aHi - aLo;
        int32  bALen = bHi - bLo;

        int32  alignLen = (aALen < bALen) ? bALen : aALen;

        int32  bgnDiag = aLo - bLo;
        int32  endDiag = aHi - bHi;

        if (bgnDiag < endDiag) {
          minDiag = bgnDiag - G.maxErate * alignLen / 2;
          maxDiag = endDiag + G.maxErate * alignLen / 2;

        } else {
          minDiag = endDiag - G.maxErate * alignLen / 2;
          maxDiag = bgnDiag + G.maxErate * alignLen / 2;
        }

        assert(minDiag < maxDiag);
      }


      //  Find seeds - hash the kmer and position from the first read, then lookup
      //  each kmer in the second read.  For unique hits, save the diagonal.  Then what?
      //  If the diagonal is too far from the expected diagonal (based on the overlap),
      //  ignore the seed.

      map<uint64,int32>   aMap;  //  Signed, to allow for easy compute of diagonal
      map<uint64,int32>   bMap;

      uint32 merSize   = 18;
      bool   dupIgnore = false;  //  Ignore dups?  Otherwise, use the first occurrence

      //  Find mers in A
    findMersAgain:

      //fprintf(stderr, "Finding mers of size %u; dupIgnore=%s\n", merSize, dupIgnore ? "true" : "false");

      {
        kMerBuilder *kb = new kMerBuilder(merSize);
        seqStream   *ss = new seqStream(aStr, aLen);
        merStream   *ms = new merStream(kb, ss, true, true);

        while (ms->nextMer() == true) {
          uint64  kmer = ms->theFMer();

          if (aMap.find(kmer) != aMap.end()) {
            if (dupIgnore == true)
              aMap[kmer] = INT32_MAX;  //  Duplicate mer, now ignored!
          } else {
            aMap[kmer] = (int32)ms->thePositionInSequence();
          }
        }

        delete ms;
      }

      if (aMap.size() == 0) {
        aMap.clear();
        bMap.clear();

        merSize--;

        if ((merSize < 8) && (dupIgnore == true)) {
          merSize   = 20;
          dupIgnore =  false;
        }

        if (merSize >= 8)
          goto findMersAgain;
      }

      //  Find mers in B

      {
        kMerBuilder *kb = new kMerBuilder(merSize);
        seqStream   *ss = new seqStream(bStr, bLen);
        merStream   *ms = new merStream(kb, ss, true, true);

        while (ms->nextMer() == true) {
          uint64  kmer = ms->theFMer();

          if (aMap.find(kmer) == aMap.end())
            //  Doesn't exist in aSeq, don't care about it.
            continue;

          int32  apos = aMap[kmer];
          int32  bpos = (int32)ms->thePositionInSequence();

          if (apos == INT32_MAX)
            //  Exists too many times in aSeq, don't care about it.
            continue;

          if ((apos - bpos < minDiag) ||
              (apos - bpos > maxDiag))
            //  Too different.
            continue;

          if (bMap.find(kmer) != bMap.end()) {
            if (dupIgnore == true)
              bMap[kmer] = INT32_MAX;  //  Duplicate mer, now ignored!
          } else {
            bMap[kmer] = bpos;
          }
        }

        delete ms;
      }

      if (bMap.size() == 0) {
        aMap.clear();
        bMap.clear();

        merSize--;

        if ((merSize < 8) && (dupIgnore == true)) {
          merSize   = 20;
          dupIgnore =  false;
        }

        if (merSize >= 8)
          goto findMersAgain;
      }

      //  Still zero?  Didn't find any unique seeds anywhere.

      if (bMap.size() == 0) {
        fprintf(stderr, "A %6u %5d-%5d ->   B %6u %5d-%5d %s NO SEEDS.\n",
                ovl->a_iid, ovl->a_bgn(WA->gkpStore), ovl->a_end(WA->gkpStore),
                ovl->b_iid, ovl->b_bgn(WA->gkpStore), ovl->b_end(WA->gkpStore),
                ovl->flipped() ? "<-" : "->");
        continue;
      }

      //  For unique mers in B, if the mer is also unique in A, add a hit.

      vector<exactMatch>   rawhits;
      vector<exactMatch>   hits;

      {
        for (map<uint64,int32>::iterator bit=bMap.begin(); bit != bMap.end(); bit++) {
          uint64  kmer = bit->first;
          int32   bpos = bit->second;

          if (bpos == INT32_MAX)
            //  Exists too many times in bSeq, don't care about it.
            continue;

          int32   apos = aMap[kmer];

          assert(apos != INT32_MAX);       //  Should never get a bMap for these

          assert(apos - bpos >= minDiag);  //  ...these too.
          assert(apos - bpos <= maxDiag);

          rawhits.push_back(exactMatch(apos, bpos, merSize));
        }
      }

      //  Sort by aPos (actually by length, then by aPos, but length is constant here).

      sort(rawhits.begin(), rawhits.end());

#if 0
      for (uint32 rr=0; rr<rawhits.size(); rr++)
        fprintf(stderr, "HIT: %d - %d diag %d\n", rawhits[rr].aBgn, rawhits[rr].bBgn, rawhits[rr].aBgn - rawhits[rr].bBgn);
#endif

      //  Chain the hits.

      hits.push_back(rawhits[0]);

      for (uint32 rr=1; rr<rawhits.size(); rr++) {
        uint32  hh = hits.size() - 1;

        int32   da = rawhits[rr].aBgn - hits[hh].aBgn;
        int32   db = rawhits[rr].bBgn - hits[hh].bBgn;

        if ((da > 0) && (da < 2 * merSize) && (da == db))
          hits[hh].tLen += da;
        else
          hits.push_back(rawhits[rr]);
      }

      //  Sort by longest

      sort(hits.begin(), hits.end());

#if 0
      for (uint32 hh=0; hh<hits.size(); hh++) {
        fprintf(stderr, "hit %02u %5d-%5d diag %d len %3u\n",
                hh,
                hits[hh].aBgn, hits[hh].bBgn,
                hits[hh].aBgn - hits[hh].bBgn,
                hits[hh].tLen);
      }
#endif

      //  Recompute.  First, mark the overlap as junk.  If we find a new overlap, it is marked as good.

      ovl->dat.ovl.forOBT = false;
      ovl->dat.ovl.forDUP = false;
      ovl->dat.ovl.forUTG = false;

      for (uint32 hh=0; hh<hits.size(); hh++) {
        Match_Node_t  match;

        match.Start  = hits[0].aBgn;   //  Begin position in a
        match.Offset = hits[0].bBgn;   //  Begin position in b
        match.Len    = merSize;        //  tLen can include mismatches
        match.Next   = 0;              //  Not used here

        int32      errors  = 0;
        Overlap_t  ovltype = Extend_Alignment(&match,       //  Initial exact match, relative to start of string
                                              aStr, aLen,
                                              bStr, bLen,
                                              aLo,  aHi,    //  Output: Regions which the match extends
                                              bLo,  bHi,
                                              errors,
                                              WA);

        int32  olapLen = 1 + min(aHi - aLo, bHi - bLo);
        double quality = (double)errors / olapLen;

        if (G.Doing_Partial_Overlaps == true)
          //  Not set for partial overlaps; just don't set it to 'none'.
          ovltype = DOVETAIL;

        if (olapLen < 40)
          ovltype = NONE;

#if 0
        fprintf(stderr, "thread %2u hit %2u  A %6u %5d-%5d %s B %6u %5d-%5d  %s -- ",
                WA->thread_id,
                hh,
                ovl->a_iid, ovl->a_bgn(WA->gkpStore), ovl->a_end(WA->gkpStore),
                ovl->flipped() ? "<-" : "->",
                ovl->b_iid, ovl->b_bgn(WA->gkpStore), ovl->b_end(WA->gkpStore));
        fprintf(stderr, "type %d A %5d-%5d  B %5d-%5d  errors %4d  quality %6.3f  llen %4d rlen %4d%s\n",
                ovltype,
                aLo, aHi, bLo, bHi,
                errors,
                quality,
                WA->editDist->Left_Delta_Len,
                WA->editDist->Right_Delta_Len,
                (ovltype == 0) ? "  FAILED" : "");
#endif

        if (ovltype != NONE) {
          ovl->dat.ovl.forOBT = (G.Doing_Partial_Overlaps == false);
          ovl->dat.ovl.forDUP = (G.Doing_Partial_Overlaps == true);
          ovl->dat.ovl.forUTG = (G.Doing_Partial_Overlaps == true);
          //  Need more...
          break;
        }
      }
    }
  }

  //  Do something.

  delete [] bRev;  //  Well, somthing more than that...

  fprintf(stderr, "Thread %u finished.\n", WA->thread_id);
}







#define MAX_LEN 40960

int
main(int argc, char **argv) {
  char  *gkpName      = NULL;
  char  *ovlName      = NULL;
  uint32  bgnID       = 0;
  uint32  endID       = UINT32_MAX;
  uint32  numThreads  = 1;

  uint64  memLimit    = 4;

  argc = AS_configure(argc, argv);

  G.initialize();

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {
      numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-erate") == 0) {
      G.maxErate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-g") == 0) {
      G.Doing_Partial_Overlaps = true;

    } else if (strcmp(argv[arg], "-memory") == 0) {
      memLimit = atoi(argv[++arg]);

    } else {
      err++;
    }

    arg++;
  }
  if (err) {
    exit(1);
  }

  fprintf(stderr, "Opening gkpStore '%s'\n", gkpName);
  gkStore  *gkpStore = new gkStore(gkpName);

  fprintf(stderr, "Opening ovlStore '%s'\n", ovlName);
  ovStore  *ovlStore = new ovStore(ovlName);

  gkRead           *read;
  gkReadData        data;

  Work_Area_t      *WA  = new Work_Area_t [numThreads];
  pthread_t        *tID = new pthread_t   [numThreads];
  pthread_attr_t    attr;

  pthread_attr_init(&attr);
  pthread_attr_setstacksize(&attr, THREAD_STACKSIZE);
  pthread_mutex_init(&balanceMutex, NULL);

  //  Initialize thread work areas.  Mirrored from overlapInCore.C

  for (uint32 tt=0; tt<numThreads; tt++) {
    fprintf(stderr, "Initialize thread %u\n", tt);

    WA[tt].String_Olap_Size  = INIT_STRING_OLAP_SIZE;
    WA[tt].String_Olap_Space = new String_Olap_t [WA[tt].String_Olap_Size];

    WA[tt].Match_Node_Size   = INIT_MATCH_NODE_SIZE;
    WA[tt].Match_Node_Space  = new Match_Node_t [WA[tt].Match_Node_Size];

    WA[tt].gkpStore    = gkpStore;

    WA[tt].status      = 0;
    WA[tt].thread_id   = tt;

    WA[tt].bgnID       = 0;
    WA[tt].endID       = 0;

    WA[tt].overlapsLen = 0;
    WA[tt].overlapsMax = 0;
    WA[tt].overlaps    = NULL;

    WA[tt].editDist    = new prefixEditDistance(G.Doing_Partial_Overlaps, G.maxErate);
  }


  //  Thread flow:
  //
  //  for reads bgn to end {
  //    Load N overlaps
  //    Load new reads - set touched reads to age=0
  //    Wait for threads to finish
  //    Increment age of reads in cache, delete reads that are too old
  //    Launch threads
  //  }
  //
  //  instead of fixed cutoff on age, use max memory usage and cull the oldest to remain below

  uint32       overlapsMax = BATCH_SIZE;

  uint32       overlapsALen = 0;
  uint32       overlapsBLen = 0;
  ovsOverlap  *overlapsA    = new ovsOverlap [overlapsMax];
  ovsOverlap  *overlapsB    = new ovsOverlap [overlapsMax];

  //  Set the globals

  uint32      *overlapsLen  = &overlapsALen;
  ovsOverlap  *overlaps     =  overlapsA;

  cache        =  new readCache(gkpStore, memLimit);

  //  Load the first batch of overlaps and reads.

  *overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax, false);
  fprintf(stderr, "Loaded %u overlaps.\n", *overlapsLen);

  cache->loadReads(overlaps, *overlapsLen);

  //  Loop over all the overlaps.

  while (*overlapsLen > 0) {

    //  Launch next batch of threads
    fprintf(stderr, "LAUNCH THREADS\n");

    //uint32  nPerThread = *overlapsLen / numThreads + 1;

    //  Globals, ugh.  These limit the threads to the range of overlaps we have loaded.  Each thread
    //  will pull out THREAD_SIZE overlaps at a time to compute, updating batchBgnID as it does so.
    //  Each thread will stop when batchBgnID > batchEndID.

    batchBgnID =  0;
    batchEndID = *overlapsLen;

    for (uint32 tt=0; tt<numThreads; tt++) {
      //WA[tt].bgnID = nPerThread * tt;
      //WA[tt].endID = nPerThread * tt + nPerThread;

      WA[tt].overlapsLen = *overlapsLen;
      WA[tt].overlaps    =  overlaps;

      if (WA[tt].endID > *overlapsLen)
        WA[tt].endID = *overlapsLen;

      int32 status = pthread_create(tID + tt, &attr, recomputeOverlaps, WA + tt);

      if (status != 0)
        fprintf(stderr, "pthread_create error:  %s\n", strerror(status)), exit(1);
    }

    //  Load next batch ov overlaps, after flipping to the unused space
    fprintf(stderr, "LOAD OVERLAPS\n");

    if (overlaps == overlapsA) {
      overlapsLen = &overlapsBLen;
      overlaps    =  overlapsB;

    } else {
      overlapsLen = &overlapsALen;
      overlaps    =  overlapsA;
    }

    *overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax, false);
    fprintf(stderr, "Loaded %u overlaps.\n", *overlapsLen);

    //  Load new reads
    fprintf(stderr, "LOAD READS\n");

    cache->loadReads(overlaps, *overlapsLen);

    //  Wait for threads to finish
    fprintf(stderr, "WAIT\n");

    for (uint32 tt=0; tt<numThreads; tt++) {
      int32 status = pthread_join(tID[tt], NULL);

      if (status != 0)
        fprintf(stderr, "pthread_join error: %s\n", strerror(status)), exit(1);
    }

    //  Expire old reads

    cache->purgeReads();
  }

  //  Goodbye.

  return(0);
}




const char *mainid = "$Id:  $";

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include "prefixEditDistance.H"

#include "AS_UTL_reverseComplement.H"

#include "Binomial_Bound.H"  //  liboverlap

#include "kMer.H"
#include "merStream.H"

#include "overlapPair-readCache.H"

#include <map>
#include <vector>
#include <algorithm>

using namespace std;


#define THREAD_SIZE  512
#define BATCH_SIZE   512 * 256


readCache        *cache         = NULL;
pthread_mutex_t   balanceMutex;
uint32            batchBgnID    = 0;
uint32            batchEndID    = 0;




class workSpace {
public:
  workSpace() {
    threadID        = 0;
    maxErate        = 0;
    partialOverlaps = false;
    gkpStore        = NULL;
    editDist        = NULL;
    overlapsLen     = 0;
    overlaps        = NULL;
  };
  ~workSpace() {
    delete editDist;
  };

public:
  uint32                 threadID;
  double                 maxErate;
  bool                   partialOverlaps;

  gkStore               *gkpStore;
  prefixEditDistance    *editDist;

  uint32                 overlapsLen;       //  Not used.
  ovOverlap            *overlaps;
};



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
  workSpace    *WA = (workSpace *)ptr;

  char         *bRev = new char [AS_MAX_READLEN];

  uint32        bgnID = 0;
  uint32        endID = 0;

  uint32        nPassed = 0;
  uint32        nFailed = 0;

  //  Lazy allocation of the prefixEditDistance structure; it's slow.

  if (WA->editDist == NULL) {
    fprintf(stderr, "initi thrad %u\n", WA->threadID);
    WA->editDist = new prefixEditDistance(WA->partialOverlaps, WA->maxErate);
    fprintf(stderr, "initi thrad %u -- DONE\n", WA->threadID);
  }


  while (getRange(bgnID, endID)) {
    //fprintf(stderr, "Thread %u computes range %u - %u\n", WA->threadID, bgnID, endID);

    for (uint32 oo=bgnID; oo<endID; oo++) {
      ovOverlap  *ovl = WA->overlaps + oo;

      //  Mark this overlap as bad.  If it recomputes, it will be marked good.

      ovl->dat.ovl.forOBT = false;
      ovl->dat.ovl.forDUP = false;
      ovl->dat.ovl.forUTG = false;

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

        bLo = bLen - ovl->b_bgn(WA->gkpStore);  //  Now correct for the reverse complemented sequence
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

        if (alignLen < 40) {
          fprintf(stderr, "A %6u %5d-%5d ->   B %6u %5d-%5d %s ALIGN LEN %d TOO SHORT.\n",
                  ovl->a_iid, ovl->a_bgn(WA->gkpStore), ovl->a_end(WA->gkpStore),
                  ovl->b_iid, ovl->b_bgn(WA->gkpStore), ovl->b_end(WA->gkpStore),
                  ovl->flipped() ? "<-" : "->", alignLen);
          continue;
        }

        int32  bgnDiag = aLo - bLo;
        int32  endDiag = aHi - bHi;

        if (bgnDiag < endDiag) {
          minDiag = bgnDiag - WA->maxErate * alignLen / 2;
          maxDiag = endDiag + WA->maxErate * alignLen / 2;

        } else {
          minDiag = endDiag - WA->maxErate * alignLen / 2;
          maxDiag = bgnDiag + WA->maxErate * alignLen / 2;
        }

        //  For very very short overlaps (mhap kindly reports 4 bp overlaps) reset the min/max
        //  to something a little more permissive.
        //if (minDiag > -5)  minDiag = -5;
        //if (maxDiag <  5)  maxDiag =  5;

        if (minDiag > maxDiag)
          fprintf(stderr, "ERROR: minDiag=%d >= maxDiag=%d\n", minDiag, maxDiag);
        assert(minDiag <= maxDiag);
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

      //  Recompute.

      for (uint32 hh=0; hh<hits.size(); hh++) {
        Match_Node_t  match;

        match.Start  = hits[0].aBgn;   //  Begin position in a
        match.Offset = hits[0].bBgn;   //  Begin position in b
        match.Len    = merSize;        //  tLen can include mismatches
        match.Next   = 0;              //  Not used here

        int32      errors  = 0;
        Overlap_t  ovltype = WA->editDist->Extend_Alignment(&match,       //  Initial exact match, relative to start of string
                                                            aStr, aLen,
                                                            bStr, bLen,
                                                            aLo,  aHi,    //  Output: Regions which the match extends
                                                            bLo,  bHi,
                                                            errors,
                                                            WA->partialOverlaps);

        int32  olapLen = 1 + min(aHi - aLo, bHi - bLo);
        double quality = (double)errors / olapLen;

        if (WA->partialOverlaps == true)
          //  ovltype isn't set for partial overlaps; just don't set it to 'none' so we can get
          //  through the 'good overlap' check below.
          ovltype = DOVETAIL;

        if (olapLen < 40)
          ovltype = NONE;

        if (ovl->dat.ovl.flipped == true) {
          bLo = bLen - bLo;  //  Now correct for the original forward sequence
          bHi = bLen - bHi;  //  Done early just for the print below
        }

#if 0
        fprintf(stderr, "thread %2u hit %2u  A %6u %5d-%5d %s B %6u %5d-%5d -- ",
                WA->threadID,
                hh,
                ovl->a_iid, ovl->a_bgn(WA->gkpStore), ovl->a_end(WA->gkpStore),
                ovl->flipped() ? "<-" : "->",
                ovl->b_iid, ovl->b_bgn(WA->gkpStore), ovl->b_end(WA->gkpStore));
        fprintf(stderr, "type %d A %5d-%5d (%5d)  B %5d-%5d (%5d)  errors %4d  quality %6.3f  llen %4d rlen %4d%s\n",
                ovltype,
                aLo, aHi, aLen, bLo, bHi, bLen,
                errors,
                quality,
                WA->editDist->Left_Delta_Len,
                WA->editDist->Right_Delta_Len,
                (ovltype == 0) ? "  FAILED" : "");
#endif

        if (ovltype == NONE)
          //  Not a good overlap, keep searching for one.
          continue;

        //  Got a good overlap.  Save it, and get out of this loop.

        if (ovl->dat.ovl.flipped == false) {
          ovl->dat.ovl.ahg5 =        (aLo);
          ovl->dat.ovl.ahg3 = aLen - (aHi);  //+1?
          ovl->dat.ovl.bhg5 =        (bLo);
          ovl->dat.ovl.bhg3 = bLen - (bHi);  //+1?

        } else {
          ovl->dat.ovl.ahg5 =        (aLo);
          ovl->dat.ovl.ahg3 = aLen - (aHi);  //+1?
          ovl->dat.ovl.bhg5 = bLen - (bLo);  //+1?
          ovl->dat.ovl.bhg3 =        (bHi);
        }

        ovl->erate(quality);

        //  Good for UTG if we aren't computing partial overlaps, and the overlap came out dovetail

        ovl->dat.ovl.forOBT = (WA->partialOverlaps == true);
        ovl->dat.ovl.forDUP = (WA->partialOverlaps == true);
        ovl->dat.ovl.forUTG = (WA->partialOverlaps == false) && (ovl->overlapIsDovetail() == true);

        //  Stop searching the hits for a good overlap.

        break;
      }

      if ((ovl->dat.ovl.forOBT == false) &&
          (ovl->dat.ovl.forDUP == false) &&
          (ovl->dat.ovl.forUTG == false))
        nFailed++;
      else
        nPassed++;
    }
  }

  //  All done.

  delete [] bRev;

  fprintf(stderr, "Thread %u finished -- %u failed %u passed.\n", WA->threadID, nFailed, nPassed);
}







#define MAX_LEN 40960

int
main(int argc, char **argv) {
  char    *gkpName         = NULL;
  char    *ovlName         = NULL;
  char    *ovlNameOut      = NULL;

  uint32   bgnID           = 0;
  uint32   endID           = UINT32_MAX;

  uint32   numThreads      = 1;

  double   maxErate        = 0.12;
  bool     partialOverlaps = false;
  uint64   memLimit        = 4;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      ovlNameOut = argv[++arg];

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {
      numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-erate") == 0) {
      maxErate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-partial") == 0) {
      partialOverlaps = true;

    } else if (strcmp(argv[arg], "-memory") == 0) {
      memLimit = atoi(argv[++arg]);

    } else {
      err++;
    }

    arg++;
  }

  if (gkpName == NULL)
    err++;
  if (ovlName == NULL)
    err++;
  if (ovlNameOut == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -G gkpStore     Mandatory, path to gkpStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Inputs can come from either a store or a file.\n");
    fprintf(stderr, "  -O ovlStore     \n");
    fprintf(stderr, "  -O ovlFile      \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "If from an ovlStore, the range of reads processed can be restricted.\n");
    fprintf(stderr, "  -b bgnID        \n");
    fprintf(stderr, "  -e endID        \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Outputs will be written to a store or file, depending on the input type\n");
    fprintf(stderr, "  -o ovlStore     \n");
    fprintf(stderr, "  -o ovlFile      \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -erate e        Overlaps are computed at 'e' fraction error\n");
    fprintf(stderr, "  -partial        Overlaps are 'overlapInCore -G' partial overlaps\n");
    fprintf(stderr, "  -memory m       Use up to 'm' GB of memory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t n            Use up to 'n' cores\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  gkStore          *gkpStore = new gkStore(gkpName);

  ovStore          *ovlStore = NULL,  *ovlStoreOut = NULL;
  ovFile           *ovlFile  = NULL,  *ovlFileOut  = NULL;

  if (AS_UTL_fileExists(ovlName, true)) {
    fprintf(stderr, "Reading overlaps from store '%s' and writing to '%s'\n",
            ovlName, ovlNameOut);
    ovlStore    = new ovStore(ovlName);
    ovlStoreOut = new ovStore(ovlNameOut, ovStoreWrite);

    if (bgnID < 1)
      bgnID = 1;
    if (endID > gkpStore->gkStore_getNumReads())
      endID = gkpStore->gkStore_getNumReads();

    ovlStore->setRange(bgnID, endID);

  } else {
    fprintf(stderr, "Reading overlaps from file '%s' and writing to '%s'\n",
            ovlName, ovlNameOut);
    ovlFile     = new ovFile(ovlName,    ovFileFull);
    ovlFileOut  = new ovFile(ovlNameOut, ovFileFullWrite);
  }

  workSpace        *WA  = new workSpace [numThreads];
  pthread_t        *tID = new pthread_t [numThreads];
  pthread_attr_t    attr;

  pthread_attr_init(&attr);
  //pthread_attr_setstacksize(&attr,  THREAD_STACKSIZE);
  pthread_mutex_init(&balanceMutex, NULL);

  //  Initialize thread work areas.  Mirrored from overlapInCore.C

  for (uint32 tt=0; tt<numThreads; tt++) {
    fprintf(stderr, "Initialize thread %u\n", tt);


    WA[tt].threadID         = tt;
    WA[tt].maxErate         = maxErate;
    WA[tt].partialOverlaps  = partialOverlaps;

    WA[tt].gkpStore         = gkpStore;
    WA[tt].editDist         = NULL;
    WA[tt].overlaps         = NULL;
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
  ovOverlap  *overlapsA    = new ovOverlap [overlapsMax];
  ovOverlap  *overlapsB    = new ovOverlap [overlapsMax];

  //  Set the globals

  uint32      *overlapsLen  = &overlapsALen;
  ovOverlap  *overlaps     =  overlapsA;

  cache        =  new readCache(gkpStore, memLimit);

  //  Load the first batch of overlaps and reads.

  if (ovlStore)
    *overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax, false);
  if (ovlFile)
    *overlapsLen = ovlFile->readOverlaps(overlaps, overlapsMax);

  fprintf(stderr, "Loaded %u overlaps.\n", *overlapsLen);


  cache->loadReads(overlaps, *overlapsLen);

  //  Loop over all the overlaps.

  while (overlapsALen + overlapsBLen > 0) {

    //  Launch next batch of threads
    fprintf(stderr, "LAUNCH THREADS\n");

    //  Globals, ugh.  These limit the threads to the range of overlaps we have loaded.  Each thread
    //  will pull out THREAD_SIZE overlaps at a time to compute, updating batchBgnID as it does so.
    //  Each thread will stop when batchBgnID > batchEndID.

    batchBgnID =  0;
    batchEndID = *overlapsLen;

    for (uint32 tt=0; tt<numThreads; tt++) {
      WA[tt].overlapsLen = *overlapsLen;
      WA[tt].overlaps    =  overlaps;

      int32 status = pthread_create(tID + tt, &attr, recomputeOverlaps, WA + tt);

      if (status != 0)
        fprintf(stderr, "pthread_create error:  %s\n", strerror(status)), exit(1);
    }

    //  Flip back to the now computed overlaps

    if (overlaps == overlapsA) {
      overlapsLen = &overlapsBLen;
      overlaps    =  overlapsB;

    } else {
      overlapsLen = &overlapsALen;
      overlaps    =  overlapsA;
    }

    //  Write recomputed overlaps - if this is the first pass through the loop,
    //  then overlapsLen will be zero
    fprintf(stderr, "WRITE OVERLAPS\n");

#warning Should we output overlaps that failed to recompute?

    if (ovlStore)
      for (uint64 oo=0; oo<*overlapsLen; oo++)
        ovlStoreOut->writeOverlap(overlaps + oo);
    if (ovlFile)
      ovlFileOut->writeOverlaps(overlaps, *overlapsLen);

    //  Load more overlaps
    fprintf(stderr, "LOAD OVERLAPS\n");

    if (ovlStore)
      *overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax, false);
    if (ovlFile)
      *overlapsLen = ovlFile->readOverlaps(overlaps, overlapsMax);

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

  delete    cache;

  delete    gkpStore;

  delete    ovlStore;
  delete    ovlStoreOut;

  delete    ovlFile;
  delete    ovlFileOut;

  delete [] overlapsA;
  delete [] overlapsB;

  return(0);
}



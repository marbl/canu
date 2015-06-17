
const char *mainid = "$Id:  $";

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include "overlapAlign.H"

#include "AS_UTL_reverseComplement.H"

#include "overlapPair-readCache.H"



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
    align           = NULL;
    overlapsLen     = 0;
    overlaps        = NULL;
  };
  ~workSpace() {
    delete align;
  };

public:
  uint32                 threadID;
  double                 maxErate;
  bool                   partialOverlaps;

  gkStore               *gkpStore;
  overlapAlign           *align;

  uint32                 overlapsLen;       //  Not used.
  ovOverlap            *overlaps;
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

  if (WA->align == NULL)
    WA->align = new overlapAlign(WA->partialOverlaps, WA->maxErate, 18);

  while (getRange(bgnID, endID)) {
    //fprintf(stderr, "Thread %u computes range %u - %u\n", WA->threadID, bgnID, endID);

    for (uint32 oo=bgnID; oo<endID; oo++) {
      ovOverlap  *ovl = WA->overlaps + oo;

#if 0
#warning SKIPPING SOME
      if ((ovl->a_iid != 1) || (ovl->b_iid != 34025))
        continue;
#endif

      //  Load A.

      char   *aStr = cache->getRead  (ovl->a_iid);
      uint32  aLen = cache->getLength(ovl->a_iid);

      int32   aLo = ovl->a_bgn();
      int32   aHi = ovl->a_end();

      assert(aLo < aHi);

      //  Load B.

      char   *bStr = cache->getRead  (ovl->b_iid);
      uint32  bLen = cache->getLength(ovl->b_iid);

      int32   bLo = ovl->b_bgn();
      int32   bHi = ovl->b_end();

      //  Make B the correct orientation, and adjust coordinates.

      if (ovl->flipped() == true) {
        memcpy(bRev, bStr, sizeof(char) * (bLen + 1));

        reverseComplementSequence(bRev, bLen);

        bStr = bRev;

        bLo = bLen - ovl->b_bgn();  //  Now correct for the reverse complemented sequence
        bHi = bLen - ovl->b_end();
      }

      assert(bLo < bHi);

      //  Compute the overlap

      //fprintf(stderr, "START %d vs %d\n", ovl->a_iid, ovl->b_iid);

      WA->align->initialize(aStr, aLen, aLo, aHi,
                            bStr, bLen, bLo, bHi);

      if (WA->align->findMinMaxDiagonal(40) == false) {
        fprintf(stderr, "A %6u %5d-%5d ->   B %6u %5d-%5d %s ALIGN LENGTH TOO SHORT.\n",
                ovl->a_iid, ovl->a_bgn(), ovl->a_end(),
                ovl->b_iid, ovl->b_bgn(), ovl->b_end(),
                ovl->flipped() ? "<-" : "->");
        continue;
      }

      if (WA->align->findSeeds(false) == false) {
        fprintf(stderr, "A %6u %5d-%5d ->   B %6u %5d-%5d %s NO SEEDS.\n",
                ovl->a_iid, ovl->a_bgn(), ovl->a_end(),
                ovl->b_iid, ovl->b_bgn(), ovl->b_end(),
                ovl->flipped() ? "<-" : "->");
        continue;
      }

      WA->align->findHits();
      WA->align->chainHits();

      if (WA->align->processHits(ovl) == true)
        nPassed++;
      else
        nFailed++;
    }
  }

  //  All done.

  delete [] bRev;

  //  Report.  The last batch has no work to do.

  if (nFailed + nPassed > 0)
    fprintf(stderr, "Thread %u finished -- %u failed %u passed.\n", WA->threadID, nFailed, nPassed);
}







#define MAX_LEN 40960

int
main(int argc, char **argv) {
  char    *gkpName         = NULL;
  char    *ovlName         = NULL;
  char    *outName         = NULL;

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
      outName = argv[++arg];

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
  if (outName == NULL)
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
    fprintf(stderr, "  -erate e        Overlaps are computed at 'e' fraction error; must be larger than the original erate\n");
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
            ovlName, outName);
    ovlStore    = new ovStore(ovlName, gkpStore);
    ovlStoreOut = new ovStore(outName, gkpStore, ovStoreWrite);

    if (bgnID < 1)
      bgnID = 1;
    if (endID > gkpStore->gkStore_getNumReads())
      endID = gkpStore->gkStore_getNumReads();

    ovlStore->setRange(bgnID, endID);

  } else {
    fprintf(stderr, "Reading overlaps from file '%s' and writing to '%s'\n",
            ovlName, outName);
    ovlFile     = new ovFile(ovlName, ovFileFull);
    ovlFileOut  = new ovFile(outName, ovFileFullWrite);
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
    WA[tt].align            = NULL;
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
    //fprintf(stderr, "LAUNCH THREADS\n");

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
    //
    //  Should we output overlaps that failed to recompute?

    if (ovlStore)
      for (uint64 oo=0; oo<*overlapsLen; oo++)
        ovlStoreOut->writeOverlap(overlaps + oo);
    if (ovlFile)
      ovlFileOut->writeOverlaps(overlaps, *overlapsLen);

    //  Load more overlaps

    if (ovlStore)
      *overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax, false);
    if (ovlFile)
      *overlapsLen = ovlFile->readOverlaps(overlaps, overlapsMax);

    fprintf(stderr, "Loaded %u overlaps.\n", *overlapsLen);

    cache->loadReads(overlaps, *overlapsLen);

    //  Wait for threads to finish

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



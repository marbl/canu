
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
#include "system.H"

#include <pthread.h>

#include "sqStore.H"
#include "ovStore.H"

#include "edlib.H"

#include "overlapReadCache.H"

#include "sequence.H"

//  The process will load BATCH_SIZE overlaps into memory, then load all the reads referenced by
//  those overlaps.  Once all data is loaded, compute threads are spawned.  Each thread will reserve
//  THREAD_SIZE overlaps to compute.  A small THREAD_SIZE relative to BATCH_SIZE will result in
//  better load balancing, but too small and the overhead of reserving overlaps will dominate (too
//  small is on the order of 1).  While threads are computing, the next batch of overlaps and reads
//  is loaded.
//
//  A large BATCH_SIZE will make startup cost large - no computes are started until the initial load
//  is finished.  To alleivate this (a little bit), the initial load is only 1/8 of the full
//  BATCH_SIZE.

#define BATCH_SIZE   1024 * 1024
#define THREAD_SIZE  128

//  Does slightly better with 2550 than 500.  Speed takes a slight hit.
#define MHAP_SLOP       500



class alignStats {
public:
  alignStats() {
    startTime       = getTime();
    reportThreshold = 0;

    clear();
  };
  ~alignStats() {
  };

  void        clear(void) {
    nSkipped = 0;
    nPassed  = 0;
    nFailed  = 0;

    nFailExtA = 0;
    nFailExtB = 0;
    nFailExt  = 0;

    nExtendedA = 0;
    nExtendedB = 0;

    nPartial  = 0;
    nDovetail = 0;

    nExt5a = 0;
    nExt3a = 0;
    nExt5b = 0;
    nExt3b = 0;
  };


  alignStats &operator+=(alignStats &that) {
    nSkipped += that.nSkipped;
    nPassed  += that.nPassed;
    nFailed  += that.nFailed;

    nFailExtA += that.nFailExtA;
    nFailExtB += that.nFailExtB;
    nFailExt  += that.nFailExt;

    nExtendedA += that.nExtendedA;
    nExtendedB += that.nExtendedB;

    nPartial  += that.nPartial;
    nDovetail += that.nDovetail;

    nExt5a += that.nExt5a;
    nExt3a += that.nExt3a;
    nExt5b += that.nExt5b;
    nExt3b += that.nExt3b;

    return(*this);
  };

  void   reportStatus(void) {

    if (nPassed + nFailed < reportThreshold)
      return;

    reportThreshold += 10000;

    fprintf(stderr, "Tested %9" F_U64P " olaps -- Skipped %8.4f%% -- Passed %8.4f%% -- %8.2f olaps/sec\n",
            nPassed + nFailed,
            100.0 * nSkipped / (nPassed + nFailed),
            100.0 * nPassed  / (nPassed + nFailed),
            (nPassed + nFailed) / (getTime() - startTime));
  };

  void    reportFinal(void) {
    fprintf(stderr, "\n");
    fprintf(stderr, " -- %" F_U64P " overlaps processed.\n", nPassed + nFailed);
    fprintf(stderr, " -- %" F_U64P " skipped.\n", nSkipped);
    fprintf(stderr, " -- %" F_U64P " failed %" F_U64P " passed (%.4f%%).\n", nFailed, nPassed, 100.0 * nPassed / (nPassed + nFailed));
    fprintf(stderr, " --\n");
    fprintf(stderr, " -- %" F_U64P " failed initial alignment, allowing A to extend\n", nFailExtA);
    fprintf(stderr, " -- %" F_U64P " failed initial alignment, allowing B to extend\n", nFailExtB);
    fprintf(stderr, " -- %" F_U64P " failed initial alignment\n", nFailExt);
    fprintf(stderr, " --\n");
    fprintf(stderr, " -- %" F_U64P " partial overlaps (of any quality)\n", nPartial);
    fprintf(stderr, " -- %" F_U64P " dovetail overlaps (before extensions, of any quality)\n", nDovetail);
    fprintf(stderr, " --\n");
    fprintf(stderr, " -- %" F_U64P "/%" F_U64P " A read dovetail extensions\n", nExt5a, nExt3a);
    fprintf(stderr, " -- %" F_U64P "/%" F_U64P " B read dovetail extensions\n", nExt5b, nExt3b);
  };

  double        startTime;
  uint64        reportThreshold;

  uint64        nSkipped;
  uint64        nPassed;
  uint64        nFailed;

  uint64        nFailExtA;
  uint64        nFailExtB;
  uint64        nFailExt;

  uint64        nExtendedA;
  uint64        nExtendedB;

  uint64        nPartial;
  uint64        nDovetail;

  uint64        nExt5a;
  uint64        nExt3a;
  uint64        nExt5b;
  uint64        nExt3b;
};



class workSpace {
public:
  workSpace() {
    threadID        = 0;

    maxErate        = 0;
    partialOverlaps = false;
    invertOverlaps  = false;

    seqStore        = NULL;
    overlapsLen     = 0;
    overlaps        = NULL;
    readSeq         = NULL;
  };
  ~workSpace() {
    delete[] readSeq;
  };

public:
  uint32                 threadID;
  double                 maxErate;
  bool                   partialOverlaps;
  bool                   invertOverlaps;
  char*                  readSeq;

  sqStore               *seqStore;

  uint32                 overlapsLen;       //  Not used.
  ovOverlap             *overlaps;
};





overlapReadCache  *rcache        = NULL;  //  Used to be just 'cache', but that conflicted with -pg: /usr/lib/libc_p.a(msgcat.po):(.bss+0x0): multiple definition of `cache'
uint32             batchPrtID    = 0;  //  When to report progress
uint32             batchPosID    = 0;  //  The current position of the batch
uint32             batchEndID    = 0;  //  The end of the batch
pthread_mutex_t    balanceMutex;

uint32             minOverlapLength = 0;

alignStats         globalStats;

bool               debug         = false;




bool
getRange(uint32 &bgnID, uint32 &endID) {

  pthread_mutex_lock(&balanceMutex);

  bgnID       = batchPosID;
  batchPosID += THREAD_SIZE;         //  Supposed to overflow.
  endID       = batchPosID;

  if (endID > batchEndID)
    endID = batchEndID;

  pthread_mutex_unlock(&balanceMutex);

  //  If we're out of overlaps, batchPosID is more than batchEndID (from the last call to this
  //  function), which makes bgnID > endID (in this call).

  return(bgnID < endID);
}




//  Try to extend the overlap on the B read.  If successful, returns new bbgn,bend and editDist and alignLen.
//
bool
extendAlignment(char  *aRead,  int32   abgn,  int32   aend,  int32  UNUSED(alen),  char const *Alabel,  uint32 Aid,
                char  *bRead,  int32  &bbgn,  int32  &bend,  int32         blen,   char const *Blabel,  uint32 Bid,
                double  maxErate,
                int32   slop,
                int32  &editDist,
                int32  &alignLen) {
  alignStats        threadStats;
  EdlibAlignResult  result  = { 0, NULL, NULL, 0, NULL, 0, 0 };
  bool              success = false;

  //  Find an alignment, allowing extensions on the B read.

  int32   bbgnExt   = max(0,    bbgn - slop);
  int32   bendExt   = min(blen, bend + slop);

  //  This probably isn't exactly correct, but close enough.
  int32   maxEdit  = (int32)ceil(max(aend - abgn, bendExt - bbgnExt) * maxErate * 1.1);

  if (debug)
    fprintf(stderr, "  align %s %6u %6d-%-6d to %s %6u %6d-%-6d", Alabel, Aid, abgn, aend, Blabel, Bid, bbgnExt, bendExt);

  result = edlibAlign(aRead + abgn,    aend    - abgn,
                      bRead + bbgnExt, bendExt - bbgnExt,
                      edlibNewAlignConfig(maxEdit, EDLIB_MODE_HW, EDLIB_TASK_LOC));

  //  Change the overlap for any extension found.

  if (result.numLocations > 0) {
    bbgn = bbgnExt + result.startLocations[0];
    bend = bbgnExt + result.endLocations[0] + 1;    //  Edlib returns 0-based positions, add one to end to get space-based.

    editDist = result.editDistance;
    alignLen = result.alignmentLength;              //  Edlib 'alignmentLength' isn't populated for TASK_LOC, so we approximate it.
    alignLen = ((aend - abgn) + (bend - bbgn) + (editDist)) / 2;

    if (debug)
      fprintf(stderr, "    aligned to %s at %6d-%-6d editDist %5d alignLen %6d qual %6.4f\n",
              Blabel, bbgn, bend, editDist, alignLen, 1.0 - editDist / (double)alignLen);

    success = true;
  } else {
    if (debug)
      fprintf(stderr, "\n");
  }

  edlibFreeAlignResult(result);

  return(success);
}



bool
finalAlignment(char *aRead, int32 alen,// char *Alabel, uint32 Aid,
               char *bRead, int32 blen,// char *Blabel, uint32 Bid,
               ovOverlap *ovl,
               double  maxErate,
               int32  &editDist,
               int32  &alignLen) {
  EdlibAlignResult  result  = { 0, NULL, NULL, 0, NULL, 0, 0 };
  bool              success = false;

  int32   abgn      = (int32)       ovl->dat.ovl.ahg5;
  int32   aend      = (int32)alen - ovl->dat.ovl.ahg3;
  int32   bbgn      = (int32)       ovl->dat.ovl.bhg5;
  int32   bend      = (int32)blen - ovl->dat.ovl.bhg3;

  int32   maxEdit  = (int32)ceil(max(aend - abgn, bend - bbgn) * maxErate * 1.1);

  result = edlibAlign(aRead + abgn, aend - abgn,
                      bRead + bbgn, bend - bbgn,
                      edlibNewAlignConfig(maxEdit, EDLIB_MODE_NW, EDLIB_TASK_LOC));  //  NOTE!  Global alignment.

  if (result.numLocations > 0) {
    editDist = result.editDistance;
    alignLen = result.alignmentLength;              //  Edlib 'alignmentLength' isn't populated for TASK_LOC, so we approximate it.
    alignLen = ((aend - abgn) + (bend - bbgn) + (editDist)) / 2;

    success = true;
  } else {
  }

  edlibFreeAlignResult(result);

  return(success);
}



void *
recomputeOverlaps(void *ptr) {
  workSpace    *WA = (workSpace *)ptr;

  uint32        bgnID = 0;
  uint32        endID = 0;

  while (getRange(bgnID, endID)) {
    alignStats  localStats;

    for (uint32 oo=bgnID; oo<endID; oo++) {
      ovOverlap  *ovl = WA->overlaps + oo;

      //  Swap IDs if requested (why would anyone want to do this?)

      if (WA->invertOverlaps) {
        ovOverlap  swapped = WA->overlaps[oo];

        WA->overlaps[oo].swapIDs(swapped);  //  Needs to be from a temporary!
      }

      //  Initialize early, just so we can use goto.

      uint32  aID       = ovl->a_iid;
      char   *aRead     = rcache->getRead(aID);
      int32   alen      = (int32)rcache->getLength(aID);
      int32   abgn      = (int32)       ovl->dat.ovl.ahg5;
      int32   aend      = (int32)alen - ovl->dat.ovl.ahg3;

      uint32  bID       = ovl->b_iid;
      char   *bRead     = WA->readSeq;
      int32   blen      = (int32)rcache->getLength(bID);
      int32   bbgn      = (int32)       ovl->dat.ovl.bhg5;
      int32   bend      = (int32)blen - ovl->dat.ovl.bhg3;

      int32   alignLen  = 1;
      int32   editDist  = INT32_MAX;

      EdlibAlignResult  result = { 0, NULL, NULL, 0, NULL, 0, 0 };

      if (debug) {
        fprintf(stderr, "--------\n");
        fprintf(stderr, "OLAP A %7" F_U32P " %6d-%-6d\n",    aID, abgn, aend);
        fprintf(stderr, "     B %7" F_U32P " %6d-%-6d %s\n", bID, bbgn, bend, (ovl->flipped() == false) ? "" : " flipped");
        fprintf(stderr, "\n");
      }

      //  Invalidate the overlap.

      ovl->evalue(AS_MAX_EVALUE);
      ovl->dat.ovl.forOBT = false;
      ovl->dat.ovl.forDUP = false;
      ovl->dat.ovl.forUTG = false;

      //  Make some bad changes, for testing
#if 0
      abgn += 100;
      aend -= 100;
      bbgn += 100;
      bend -= 100;
#endif

      //  Too short?  Don't bother doing anything.
      //
      //  Warning!  Edlib failed on a 10bp to 10bp (extended to 5kbp) alignment.

      if ((aend - abgn < minOverlapLength) ||
          (bend - bbgn < minOverlapLength)) {
        localStats.nSkipped++;
        goto finished;
      }

      //  Grab the B read sequence.

      strcpy(bRead, rcache->getRead(bID));

      //  If flipped, reverse complement the B read.

      if (ovl->flipped() == true)
        reverseComplementSequence(bRead, blen);

      //
      //  Find initial alignments, allowing one, then the other, sequence to be extended as needed.
      //

      if (extendAlignment(bRead, bbgn, bend, blen, "B", bID,
                          aRead, abgn, aend, alen, "A", aID,
                          WA->maxErate, MHAP_SLOP,
                          editDist,
                          alignLen) == false) {
        localStats.nFailExtA++;
      }

      if (extendAlignment(aRead, abgn, aend, alen, "A", aID,
                          bRead, bbgn, bend, blen, "B", bID,
                          WA->maxErate, MHAP_SLOP,
                          editDist,
                          alignLen) == false) {
        localStats.nFailExtB++;
      }

      //  If no alignments were found, fail.

      if (alignLen == 1) {
        localStats.nFailExt++;
        goto finished;
      }

      //  Update the overlap.

      ovl->dat.ovl.ahg5 = abgn;
      ovl->dat.ovl.ahg3 = alen - aend;

      ovl->dat.ovl.bhg5 = bbgn;
      ovl->dat.ovl.bhg3 = blen - bend;

      if (debug) {
        fprintf(stderr, "\n");
        fprintf(stderr, "init A %7" F_U32P " %6d-%-6d\n", aID, abgn, aend);
        fprintf(stderr, "     B %7" F_U32P " %6d-%-6d\n", bID, bbgn, bend);
        fprintf(stderr, "\n");
      }

      //  If we're just doing partial alignments or if we've found a dovetail, we're all done.

      if (WA->partialOverlaps == true) {
        localStats.nPartial++;
        goto finished;
      }

      if (ovl->overlapIsDovetail() == true) {
        localStats.nDovetail++;
        goto finished;
      }

#warning do we need to check for contained too?



      //  Otherwise, try to extend the alignment to make a dovetail overlap.

      {
        int32  ahg5 = ovl->dat.ovl.ahg5;
        int32  ahg3 = ovl->dat.ovl.ahg3;

        int32  bhg5 = ovl->dat.ovl.bhg5;
        int32  bhg3 = ovl->dat.ovl.bhg3;

        int32  slop = 0;

        if ((ahg5 >= bhg5) && (bhg5 > 0)) {
          //fprintf(stderr, "extend 5' by B=%d\n", bhg5);
          ahg5 -= bhg5;
          bhg5 -= bhg5;   //  Now zero.
          slop  = bhg5 * WA->maxErate + 100;

          abgn = (int32)       ahg5;
          aend = (int32)alen - ahg3;

          bbgn = (int32)       bhg5;
          bend = (int32)blen - bhg3;

          if (extendAlignment(bRead, bbgn, bend, blen, "Bb5", bID,
                              aRead, abgn, aend, alen, "Ab5", aID,
                              WA->maxErate, slop,
                              editDist,
                              alignLen) == true) {
            ahg5 = abgn;
            //ahg3 = alen - aend;
          } else {
            ahg5 = ovl->dat.ovl.ahg5;
            bhg5 = ovl->dat.ovl.bhg5;
          }
          localStats.nExt5b++;
        }

        if ((bhg5 >= ahg5) && (ahg5 > 0)) {
          //fprintf(stderr, "extend 5' by A=%d\n", ahg5);
          bhg5 -= ahg5;
          ahg5 -= ahg5;   //  Now zero.
          slop  = ahg5 * WA->maxErate + 100;

          abgn = (int32)       ahg5;
          aend = (int32)alen - ahg3;

          bbgn = (int32)       bhg5;
          bend = (int32)blen - bhg3;

          if (extendAlignment(aRead, abgn, aend, alen, "Aa5", aID,
                              bRead, bbgn, bend, blen, "Ba5", bID,
                              WA->maxErate, slop,
                              editDist,
                              alignLen) == true) {
            bhg5 = bbgn;
            //bhg3 = blen - bend;
          } else {
            bhg5 = ovl->dat.ovl.bhg5;
            ahg5 = ovl->dat.ovl.ahg5;
          }
          localStats.nExt5a++;
        }



        if ((bhg3 >= ahg3) && (ahg3 > 0)) {
          //fprintf(stderr, "extend 3' by A=%d\n", ahg3);
          bhg3 -= ahg3;
          ahg3 -= ahg3;   //  Now zero.
          slop  = ahg3 * WA->maxErate + 100;

          abgn = (int32)       ahg5;
          aend = (int32)alen - ahg3;

          bbgn = (int32)       bhg5;
          bend = (int32)blen - bhg3;

          if (extendAlignment(aRead, abgn, aend, alen, "Aa3", aID,
                              bRead, bbgn, bend, blen, "Ba3", bID,
                              WA->maxErate, slop,
                              editDist,
                              alignLen) == true) {
            //bhg5 = bbgn;
            bhg3 = blen - bend;
          } else {
            bhg3 = ovl->dat.ovl.bhg3;
            ahg3 = ovl->dat.ovl.ahg3;
          }
          localStats.nExt3a++;
        }

        if ((ahg3 >= bhg3) && (bhg3 > 0)) {
          //fprintf(stderr, "extend 3' by B=%d\n", bhg3);
          ahg3 -= bhg3;
          bhg3 -= bhg3;   //  Now zero.
          slop  = bhg3 * WA->maxErate + 100;

          abgn = (int32)       ahg5;
          aend = (int32)alen - ahg3;

          bbgn = (int32)       bhg5;
          bend = (int32)blen - bhg3;

          if (extendAlignment(bRead, bbgn, bend, blen, "Bb3", bID,
                              aRead, abgn, aend, alen, "Ab3", aID,
                              WA->maxErate, slop,
                              editDist,
                              alignLen) == true) {
            //ahg5 = abgn;
            ahg3 = alen - aend;
          } else {
            ahg3 = ovl->dat.ovl.ahg3;
            bhg3 = ovl->dat.ovl.bhg3;
          }
          localStats.nExt3b++;
        }

        //  Now reset the overlap.

        ovl->dat.ovl.ahg5 = ahg5;
        ovl->dat.ovl.ahg3 = ahg3;

        ovl->dat.ovl.bhg5 = bhg5;
        ovl->dat.ovl.bhg3 = bhg3;
      }  //  If not a contained overlap



      //  If we're still not dovetail, nothing more we want to do.  Let the overlap be trashed.


      if (debug) {
        fprintf(stderr, "\n");
        fprintf(stderr, "fini A %7" F_U32P " %6d-%-6d %d %d\n",    aID, abgn, aend, ovl->a_bgn(), ovl->a_end());
        fprintf(stderr, "     B %7" F_U32P " %6d-%-6d %d %d %s\n", bID, bbgn, bend, ovl->b_bgn(), ovl->b_end(), (ovl->flipped() == false) ? "" : " flipped");
        fprintf(stderr, "\n");
      }

      finalAlignment(aRead, alen,// "A", aID,
                     bRead, blen,// "B", bID,
                     ovl, WA->maxErate, editDist, alignLen);


    finished:

      //  Trash the overlap if it's junky quality.

      double  eRate = editDist / (double)alignLen;

      if ((alignLen < minOverlapLength) ||
          (eRate    > WA->maxErate)) {
        localStats.nFailed++;
        ovl->evalue(AS_MAX_EVALUE);
        ovl->dat.ovl.forOBT = false;
        ovl->dat.ovl.forDUP = false;
        ovl->dat.ovl.forUTG = false;

      } else {
        localStats.nPassed++;
        ovl->erate(eRate);
        ovl->dat.ovl.forOBT = (WA->partialOverlaps == true);
        ovl->dat.ovl.forDUP = (WA->partialOverlaps == true);
        ovl->dat.ovl.forUTG = (WA->partialOverlaps == false) && (ovl->overlapIsDovetail() == true);
      }

    }  //  Over all overlaps in this range


    //  Log that we've done stuff

    pthread_mutex_lock(&balanceMutex);
    globalStats += localStats;
    globalStats.reportStatus();
    localStats.clear();
    pthread_mutex_unlock(&balanceMutex);
  }  //  Over all ranges

  return(NULL);
}






int
main(int argc, char **argv) {
  char    *seqName         = NULL;
  char    *ovlName         = NULL;
  char    *outName         = NULL;

  uint32   bgnID           = 0;
  uint32   endID           = UINT32_MAX;

  uint32   numThreads      = 1;

  double   maxErate        = 0.12;
  bool     partialOverlaps = false;
  bool     invertOverlaps  = false;

  uint64   memLimit        = 4;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

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

    } else if (strcmp(argv[arg], "-invert") == 0) {
      invertOverlaps = true;

    } else if (strcmp(argv[arg], "-memory") == 0) {
      memLimit = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-len") == 0) {
      minOverlapLength = atoi(argv[++arg]);

    } else {
      err++;
    }

    arg++;
  }

  if (seqName == NULL)
    err++;
  if (ovlName == NULL)
    err++;
  if (outName == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -S seqStore     Mandatory, path to seqStore\n");
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
    fprintf(stderr, "  -partial        Overlaps are 'overlapInCore -S' partial overlaps\n");
    fprintf(stderr, "  -memory m       Use up to 'm' GB of memory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t n            Use up to 'n' cores\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Advanced options:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -invert         Invert the overlap A <-> B before aligning (they are not re-inverted before output)\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  sqStore          *seqStore = new sqStore(seqName);

  ovStore          *ovlStore = NULL;
  ovStoreWriter    *outStore = NULL;
  ovFile           *ovlFile  = NULL;
  ovFile           *outFile  = NULL;

  if (directoryExists(ovlName)) {
    fprintf(stderr, "Reading overlaps from store '%s' and writing to '%s'\n",
            ovlName, outName);
    ovlStore = new ovStore(ovlName, seqStore);
    outStore = new ovStoreWriter(outName, seqStore);

    if (bgnID < 1)
      bgnID = 1;
    if (endID > seqStore->sqStore_lastReadID())
      endID = seqStore->sqStore_lastReadID();

    ovlStore->setRange(bgnID, endID);

  } else {
    fprintf(stderr, "Reading overlaps from file '%s' and writing to '%s'\n",
            ovlName, outName);
    ovlFile = new ovFile(seqStore, ovlName, ovFileFull);
    outFile = new ovFile(seqStore, outName, ovFileFullWrite);
  }

  workSpace        *WA  = new workSpace [numThreads];
  pthread_t        *tID = new pthread_t [numThreads];
  pthread_attr_t    attr;

  pthread_attr_init(&attr);
  pthread_attr_setstacksize(&attr,  12 * 131072);
  pthread_mutex_init(&balanceMutex, NULL);

  //  Initialize thread work areas.  Mirrored from overlapInCore.C

  for (uint32 tt=0; tt<numThreads; tt++) {
    fprintf(stderr, "Initialize thread %u\n", tt);

    WA[tt].threadID         = tt;
    WA[tt].maxErate         = maxErate;
    WA[tt].partialOverlaps  = partialOverlaps;
    WA[tt].invertOverlaps   = invertOverlaps;

    WA[tt].seqStore         = seqStore;
    WA[tt].overlaps         = NULL;

    // preallocate some work thread memory for common tasks to avoid allocation
    WA[tt].readSeq = new char[AS_MAX_READLEN+1];
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

  class overlapBlock {
  public:
    overlapBlock() {
      _len = 0;
      _max = BATCH_SIZE;
      _ovl = new ovOverlap[BATCH_SIZE];
    }
    ~overlapBlock() {
      delete [] _ovl;
    };
    uint32      _len;
    uint32      _max;
    ovOverlap  *_ovl;
  };

  overlapBlock  overlapsA;
  overlapBlock  overlapsB;

  overlapBlock *overlaps  = &overlapsA;

  rcache = new overlapReadCache(seqStore, memLimit);

  //  Load the first batch of overlaps and reads.

  if (ovlStore)
    overlaps->_len = ovlStore->loadBlockOfOverlaps(overlaps->_ovl, overlaps->_max);

  if (ovlFile)
    overlaps->_len = ovlFile->readOverlaps(overlaps->_ovl, overlaps->_max);

  fprintf(stderr, "Loaded %u overlaps.\n", overlaps->_len);

  rcache->loadReads(overlaps->_ovl, overlaps->_len);

  //  Loop over all the overlaps.

  while (overlapsA._len + overlapsB._len > 0) {

    //  Launch next batch of threads
    //fprintf(stderr, "LAUNCH THREADS\n");

    //  Globals, ugh.  These limit the threads to the range of overlaps we have loaded.  Each thread
    //  will pull out THREAD_SIZE overlaps at a time to compute, updating batchPosID as it does so.
    //  Each thread will stop when batchPosID > batchEndID.

    batchPrtID =  0;
    batchPosID =  0;
    batchEndID = overlaps->_len;

    for (uint32 tt=0; tt<numThreads; tt++) {
      WA[tt].overlapsLen = overlaps->_len;
      WA[tt].overlaps    = overlaps->_ovl;

      int32 status = pthread_create(tID + tt, &attr, recomputeOverlaps, WA + tt);

      if (status != 0)
        fprintf(stderr, "pthread_create error:  %s\n", strerror(status)), exit(1);
    }

    //  Flip back to the now computed overlaps

    if (overlaps == &overlapsA)
      overlaps = &overlapsB;
    else
      overlaps = &overlapsA;

    //  Write recomputed overlaps - if this is the first pass through the loop,
    //  then overlapsLen will be zero
    //
    //  Should we output overlaps that failed to recompute?

    if (ovlStore)
      for (uint64 oo=0; oo<overlaps->_len; oo++)
        outStore->writeOverlap(overlaps->_ovl + oo);
    if (ovlFile)
      outFile->writeOverlaps(overlaps->_ovl, overlaps->_len);

    //  Load more overlaps

    if (ovlStore)
      overlaps->_len = ovlStore->loadBlockOfOverlaps(overlaps->_ovl, overlaps->_max);
    if (ovlFile)
      overlaps->_len = ovlFile->readOverlaps(overlaps->_ovl, overlaps->_max);

    fprintf(stderr, "Loaded %u overlaps.\n", overlaps->_len);

    rcache->loadReads(overlaps->_ovl, overlaps->_len);

    //  Wait for threads to finish

    for (uint32 tt=0; tt<numThreads; tt++) {
      //fprintf(stderr, "wait for thread %u\n", tt);
      int32 status = pthread_join(tID[tt], NULL);
      //fprintf(stderr, "joined thread %u\n", tt);

      if (status != 0)
        fprintf(stderr, "pthread_join error: %s\n", strerror(status)), exit(1);
    }

    //  Expire old reads

    rcache->purgeReads();
  }

  //  Report.  The last batch has no work to do.

  globalStats.reportFinal();

  //  Goodbye.

  delete    rcache;

  delete seqStore;

  delete    ovlStore;
  delete    outStore;

  delete    ovlFile;
  delete    outFile;

  delete [] WA;
  delete [] tID;

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  return(0);
}



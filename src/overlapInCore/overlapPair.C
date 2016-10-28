
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
 *    Brian P. Walenz from 2015-MAR-27 to 2015-JUL-20
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include <pthread.h>

#include "gkStore.H"
#include "ovStore.H"

#include "edlib.H"

#include "overlapReadCache.H"

#include "AS_UTL_reverseComplement.H"

#include "timeAndSize.H" //  getTime();

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

#define MHAP_SLOP    500
//#define DEBUG 1


overlapReadCache  *rcache        = NULL;  //  Used to be just 'cache', but that conflicted with -pg: /usr/lib/libc_p.a(msgcat.po):(.bss+0x0): multiple definition of `cache'
uint32             batchPrtID    = 0;  //  When to report progress
uint32             batchPosID    = 0;  //  The current position of the batch
uint32             batchEndID    = 0;  //  The end of the batch
pthread_mutex_t    balanceMutex;

uint32 minOverlapLength          = 0;


class workSpace {
public:
  workSpace() {
    threadID        = 0;

    maxErate        = 0;
    partialOverlaps = false;
    invertOverlaps  = false;

    gkpStore        = NULL;
    //analyze         = NULL;
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

  gkStore               *gkpStore;

  uint32                 overlapsLen;       //  Not used.
  ovOverlap             *overlaps;
};




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



void *
recomputeOverlaps(void *ptr) {
  workSpace    *WA = (workSpace *)ptr;

  uint32        bgnID = 0;
  uint32        endID = 0;

  uint32        nPassed = 0;
  uint32        nFailed = 0;
  uint32	nTested = 0;

  //  Lazy allocation of the prefixEditDistance structure; it's slow.

  //if (WA->analyze == NULL)
  //  WA->analyze = new analyzeAlignment();

  while (getRange(bgnID, endID)) {
    double  startTime = getTime();

    for (uint32 oo=bgnID; oo<endID; oo++) {
      ovOverlap  *ovl = WA->overlaps + oo;

      if (WA->invertOverlaps) {
        ovOverlap  swapped = WA->overlaps[oo];

        WA->overlaps[oo].swapIDs(swapped);  //  Needs to be from a temporary!
      }

      //  Invalidate the overlap.

      ovl->evalue(AS_MAX_EVALUE);
      ovl->dat.ovl.forOBT = false;
      ovl->dat.ovl.forDUP = false;
      ovl->dat.ovl.forUTG = false;


      uint32  aID  = ovl->a_iid;
      uint32  bID  = ovl->b_iid;

      //  Compute the overlap
      if (ovl->a_end() - ovl->a_bgn() + 1 < minOverlapLength && ovl->b_end() - ovl->b_bgn() + 1 < minOverlapLength) { 
         continue; 
      }

#ifdef DEBUG
nTested++;
if (nTested % 1000 == 0) {
   double  deltaTime = getTime() - startTime;
   fprintf(stderr, "*******Thread %2u computed overlaps %7u - %7u in %7.3f seconds - %6.2f olaps per second (%8u fail %8u pass)\n",
              WA->threadID, bgnID, endID, deltaTime, (endID - bgnID) / deltaTime, nFailed, nPassed);
}
#endif
  char *bRead          = WA->readSeq;
  int32 astart         = (int32)ovl->a_bgn();
  int32 aend           = (int32)ovl->a_end();
  int32 astartExtended = max((int32)0, (int32)ovl->a_bgn() - MHAP_SLOP);
  int32 aendExtended   = min((int32)rcache->getLength(aID), (int32)ovl->a_end() + MHAP_SLOP);
  int32 bstart         = (int32)ovl->b_bgn();
  int32 bend           = (int32)ovl->b_end();
  int32 bstartExtended = max((int32)0, (int32)ovl->b_bgn() - MHAP_SLOP);
  int32 bendExtended   = min((int32)rcache->getLength(bID), (int32)ovl->b_end() + MHAP_SLOP);
  strcpy(bRead, rcache->getRead(bID));
  if (ovl->flipped()) {
     reverseComplementSequence(bRead, rcache->getLength(bID));
     bstart         = (int32)rcache->getLength(bID) - (int32)ovl->b_bgn();
     bend           = (int32)rcache->getLength(bID) - (int32)ovl->b_end();
     bstartExtended = max((int32)0, (int32)rcache->getLength(bID) - (int32)ovl->b_bgn() - MHAP_SLOP);
     bendExtended   = min((int32)rcache->getLength(bID), (int32)rcache->getLength(bID) - (int32)ovl->b_end() + MHAP_SLOP);
  }

  int tolerance =  (int)ceil((double)max(aendExtended-astartExtended, bendExtended-bstartExtended)*WA->maxErate*1.1);
  EdlibAlignResult bQuery = edlibAlign(rcache->getRead(aID)+astart, aend-astart, bRead+bstartExtended, bendExtended-bstartExtended, edlibNewAlignConfig(tolerance, EDLIB_MODE_HW, EDLIB_TASK_LOC));
  EdlibAlignResult aQuery = edlibAlign(bRead+bstart, bend-bstart, rcache->getRead(aID)+astartExtended, aendExtended-astartExtended, edlibNewAlignConfig(tolerance, EDLIB_MODE_HW, EDLIB_TASK_LOC));

  uint32 alignmentLength = 0;
  double dist = 0;

#ifdef DEBUG
fprintf(stderr, "Overlap between %d and %d at %d found %d %d hits\n", aID, bID, tolerance, bQuery.numLocations, aQuery.numLocations);
#endif
  if (aQuery.numLocations >= 1 || bQuery.numLocations >= 1) {
     // if we couldn't find one of the options, try trimming and re-computing 
     if (bQuery.numLocations == 0) {
        ovl->dat.ovl.ahg5 = aQuery.startLocations[0] + astartExtended;
        ovl->dat.ovl.ahg3 = rcache->getLength(aID) - (aQuery.endLocations[0] + astartExtended + 1);
        alignmentLength = max(alignmentLength, (uint32)(aQuery.endLocations[0] - aQuery.startLocations[0]));
        dist = min(aQuery.editDistance, (int)dist);
        edlibFreeAlignResult(bQuery);
        bQuery = edlibAlign(rcache->getRead(aID)+astart, aend-astart, bRead+bstartExtended, bendExtended-bstartExtended, edlibNewAlignConfig(tolerance, EDLIB_MODE_HW, EDLIB_TASK_LOC));
     }
     if (aQuery.numLocations == 0) {
        ovl->dat.ovl.bhg5 = bQuery.startLocations[0] + bstartExtended;
        ovl->dat.ovl.bhg3 = rcache->getLength(bID) - (bQuery.endLocations[0] + bstartExtended + 1);
        alignmentLength = bQuery.endLocations[0] - bQuery.startLocations[0];
        dist = bQuery.editDistance;
        edlibFreeAlignResult(aQuery);
        aQuery = edlibAlign(bRead+bstart, bend-bstart, rcache->getRead(aID)+astartExtended, aendExtended-astartExtended, edlibNewAlignConfig(tolerance, EDLIB_MODE_HW, EDLIB_TASK_LOC));
     }

     // now update the trim points based on where the overlapping broke
     // the aligner computes 0-based end positions so for matching ACGTA to ACTGTA positiosn are 0-4 so we need to adjust for that
     if (bQuery.numLocations >= 1) {
        ovl->dat.ovl.bhg5 = bQuery.startLocations[0] + bstartExtended;
        ovl->dat.ovl.bhg3 = rcache->getLength(bID) - (bQuery.endLocations[0] + bstartExtended + 1);
        alignmentLength = bQuery.endLocations[0] - bQuery.startLocations[0];
        dist = bQuery.editDistance;
     }
     if (aQuery.numLocations >= 1) {
        ovl->dat.ovl.ahg5 = aQuery.startLocations[0] + astartExtended;
        ovl->dat.ovl.ahg3 = rcache->getLength(aID) - (aQuery.endLocations[0] + astartExtended + 1);
        alignmentLength = max(alignmentLength, (uint32)(aQuery.endLocations[0] - aQuery.startLocations[0]));
        dist = min(aQuery.editDistance, (int)dist);
     }
     edlibFreeAlignResult(aQuery);
     edlibFreeAlignResult(bQuery);

#ifdef DEBUG
fprintf(stderr, "Expected overlap between %d and %d from %d - %d and %d - %d found overlap from %d - %d and %d - %d length %d dist %f\n", aID, bID, astart, aend, bstart, bend, ovl->a_bgn(), ovl->a_end(), ovl->b_bgn(), ovl->b_end(), alignmentLength, dist);
#endif

     bool changed = true;
     tolerance = (int)(alignmentLength * WA->maxErate) + 1;
     EdlibAlignResult result;

     // extend to the ends if we are not looking for partial and we can, don't extend contains
     if (changed && WA->partialOverlaps == false && !ovl->overlapIsDovetail()) {
        bstart = ovl->flipped() ? rcache->getLength(bID) - ovl->b_bgn() : ovl->b_bgn();
        bend = ovl->flipped() ? rcache->getLength(bID) - ovl->b_end() : ovl->b_end();
#ifdef DEBUG
fprintf(stderr, "Overlap %d %d (%d %d) invert %d is %f error and not dovetail with %d %d and %d %d\n", aID, bID, rcache->getLength(aID) , rcache->getLength(bID), ovl->flipped(), dist/*result.editDistance*/,  ovl->dat.ovl.ahg5,  ovl->dat.ovl.ahg3,  ovl->dat.ovl.bhg5,  ovl->dat.ovl.bhg3);
#endif
        // dist = result.editDistance;
        // check these cases one by one and extend both concordantly with each other
        // first is a contained in b
        if (rcache->getLength(aID) <= rcache->getLength(bID) && ovl->dat.ovl.ahg5 >= 0 && ovl->dat.ovl.ahg3 >= 0 && ovl->dat.ovl.bhg5 >= ovl->dat.ovl.ahg5 && ovl->dat.ovl.bhg3 >= ovl->dat.ovl.ahg3 && ((double)(ovl->dat.ovl.ahg5 + ovl->dat.ovl.ahg3 + dist) / ((double)(alignmentLength + ovl->dat.ovl.ahg5 + ovl->dat.ovl.ahg3))) <= WA->maxErate) {
           ovl->dat.ovl.bhg5 = max(0, ovl->dat.ovl.bhg5 - ovl->dat.ovl.ahg5); ovl->dat.ovl.ahg5 = 0;
           ovl->dat.ovl.bhg3 = max(0, ovl->dat.ovl.bhg3 - ovl->dat.ovl.ahg3); ovl->dat.ovl.ahg3 = 0;
           changed = true;
#ifdef DEBUG
fprintf(stderr, "Overlap %d %d case 1 acontained\n", aID, bID);
#endif
        }
        // second is b contained (both b hangs can be extended)
        //
        else if (rcache->getLength(aID) >= rcache->getLength(bID) && ovl->dat.ovl.bhg5 >= 0 && ovl->dat.ovl.bhg3 >= 0 && ovl->dat.ovl.ahg5 >= ovl->dat.ovl.bhg5 && ovl->dat.ovl.ahg3 >= ovl->dat.ovl.bhg3 && ((double)(ovl->dat.ovl.bhg5 + ovl->dat.ovl.bhg3 + dist) / ((double)(alignmentLength + ovl->dat.ovl.bhg5 + ovl->dat.ovl.bhg3))) <= WA->maxErate) {
           ovl->dat.ovl.ahg5 = max(0, ovl->dat.ovl.ahg5 - ovl->dat.ovl.bhg5); ovl->dat.ovl.bhg5 = 0;
           ovl->dat.ovl.ahg3 = max(0, ovl->dat.ovl.ahg3 - ovl->dat.ovl.bhg3); ovl->dat.ovl.bhg3 = 0;
           changed = true;
#ifdef DEBUG
fprintf(stderr, "Overlap %d %d case 2 bconatined\n", aID, bID);
#endif
        }
        // third is 5' dovetal  ---------->
        //                          ---------->
        //                          or
        //                          <---------
        //                         bhg5 here is always first overhang on b read
        //
        else if (ovl->dat.ovl.ahg3 <= ovl->dat.ovl.bhg3 && (ovl->dat.ovl.ahg3 >= 0 && ((double)(ovl->dat.ovl.ahg3 + dist) / ((double)(alignmentLength + ovl->dat.ovl.ahg3))) <= WA->maxErate) &&
                (ovl->dat.ovl.bhg5 >= 0 && ((double)(ovl->dat.ovl.bhg5 + dist) / ((double)(alignmentLength + ovl->dat.ovl.bhg5))) <= WA->maxErate)) {
           ovl->dat.ovl.ahg5 = max(0, ovl->dat.ovl.ahg5 - ovl->dat.ovl.bhg5); ovl->dat.ovl.bhg5 = 0;
           ovl->dat.ovl.bhg3 = max(0, ovl->dat.ovl.bhg3 - ovl->dat.ovl.ahg3); ovl->dat.ovl.ahg3 = 0;
           changed = true;
#ifdef DEBUG
fprintf(stderr, "Overlap %d %d case 3 5' dovetail \n", aID, bID);
#endif
        }
        //
        // fourth is 3' dovetail    ---------->
        //                     ---------->
        //                     or
        //                     <----------
        //                     bhg5 is always first overhang on b read
        else if (ovl->dat.ovl.ahg5 <= ovl->dat.ovl.bhg5 && (ovl->dat.ovl.ahg5 >= 0 && ((double)(ovl->dat.ovl.ahg5 + dist) / ((double)(alignmentLength + ovl->dat.ovl.ahg5))) <= WA->maxErate) &&
                (ovl->dat.ovl.bhg3 >= 0 && ((double)(ovl->dat.ovl.bhg3 + dist) / ((double)(alignmentLength + ovl->dat.ovl.bhg3))) <= WA->maxErate)) {
           ovl->dat.ovl.bhg5 = max(0, ovl->dat.ovl.bhg5 - ovl->dat.ovl.ahg5); ovl->dat.ovl.ahg5 = 0;
           ovl->dat.ovl.ahg3 = max(0, ovl->dat.ovl.ahg3 - ovl->dat.ovl.bhg3); ovl->dat.ovl.bhg3 = 0;
           changed = true;
#ifdef DEBUG
fprintf(stderr, "Overlap %d %d case 4 3' dovetail \n", aID, bID);
#endif
        }
     }

#ifdef DEBUG
fprintf(stderr, "Recomputed overlap from %d to %d is %d - %d and %d - %d\n", aID, bID, ovl->a_bgn(), ovl->a_end(), (ovl->flipped() ? rcache->getLength(bID) - ovl->b_bgn() : ovl->b_bgn()), (ovl->flipped() ? rcache->getLength(bID) - ovl->b_end() : ovl->b_end()));
#endif
     // now compute the final
     if (changed) { 
        bstart = ovl->flipped() ? rcache->getLength(bID) - ovl->b_bgn() : ovl->b_bgn();
        bend = ovl->flipped() ? rcache->getLength(bID) - ovl->b_end() : ovl->b_end();
        result = edlibAlign(rcache->getRead(aID)+ovl->a_bgn(), ovl->a_end()-ovl->a_bgn(), bRead+bstart, bend-bstart, edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, EDLIB_TASK_LOC));
        if (result.numLocations >= 1) {
           dist = result.editDistance;
           alignmentLength = ovl->a_end() - ovl->a_bgn(); 
        } else { 
           dist = ovl->a_end() - ovl->a_bgn();
           alignmentLength = 0;
        }
        edlibFreeAlignResult(result);
     }
#ifdef DEBUG
fprintf(stderr, "Done and error rate for this overlap between %d and %d is %d bp and %f errors is dovetail %d\n", aID, bID, alignmentLength, dist, ovl->overlapIsDovetail());
#endif
  } else {
     edlibFreeAlignResult(aQuery);
     edlibFreeAlignResult(bQuery);
  }

  if (alignmentLength >= minOverlapLength && (dist / (double) (alignmentLength)) <= WA->maxErate) {
        nPassed++;

        ovl->erate((double)dist/(alignmentLength));
        ovl->dat.ovl.forOBT = (WA->partialOverlaps == true);
        ovl->dat.ovl.forDUP = (WA->partialOverlaps == true);
        ovl->dat.ovl.forUTG = (WA->partialOverlaps == false) && (ovl->overlapIsDovetail() == true);

      } else {
        nFailed++;
        ovl->evalue(AS_MAX_EVALUE);

        ovl->dat.ovl.forOBT = false;
        ovl->dat.ovl.forDUP = false;
        ovl->dat.ovl.forUTG = false;
      }
    }

#ifdef DEBUG
    double  deltaTime = getTime() - startTime;
    fprintf(stderr, "Thread %2u computed overlaps %7u - %7u in %7.3f seconds - %6.2f olaps per second (%8u fail %8u pass)\n",
            WA->threadID, bgnID, endID, deltaTime, (endID - bgnID) / deltaTime, nFailed, nPassed);
#endif
  }

  //  Report.  The last batch has no work to do.

  if (nFailed + nPassed > 0)
    fprintf(stderr, "Thread %u finished -- %u failed %u passed.\n", WA->threadID, nFailed, nPassed);

  return(NULL);
}






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
  bool     invertOverlaps  = false;

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
    fprintf(stderr, "Advanced options:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -invert         Invert the overlap A <-> B before aligning (they are not re-inverted before output)\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  gkStore          *gkpStore = gkStore::gkStore_open(gkpName);

  ovStore          *ovlStore = NULL;
  ovStoreWriter    *outStore = NULL;
  ovFile           *ovlFile  = NULL;
  ovFile           *outFile  = NULL;

  if (AS_UTL_fileExists(ovlName, true)) {
    fprintf(stderr, "Reading overlaps from store '%s' and writing to '%s'\n",
            ovlName, outName);
    ovlStore = new ovStore(ovlName, gkpStore);
    outStore = new ovStoreWriter(outName, gkpStore);

    if (bgnID < 1)
      bgnID = 1;
    if (endID > gkpStore->gkStore_getNumReads())
      endID = gkpStore->gkStore_getNumReads();

    ovlStore->setRange(bgnID, endID);

  } else {
    fprintf(stderr, "Reading overlaps from file '%s' and writing to '%s'\n",
            ovlName, outName);
    ovlFile = new ovFile(gkpStore, ovlName, ovFileFull);
    outFile = new ovFile(gkpStore, outName, ovFileFullWrite);
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

      WA[tt].gkpStore         = gkpStore;
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

  uint32       overlapsMax = BATCH_SIZE;

  uint32       overlapsALen = 0;
  uint32       overlapsBLen = 0;
  ovOverlap  *overlapsA    = ovOverlap::allocateOverlaps(gkpStore, overlapsMax);
  ovOverlap  *overlapsB    = ovOverlap::allocateOverlaps(gkpStore, overlapsMax);

  //  Set the globals

  uint32      *overlapsLen  = &overlapsALen;
  ovOverlap  *overlaps      =  overlapsA;

  rcache = new overlapReadCache(gkpStore, memLimit);

  //  Load the first batch of overlaps and reads.  Purposely loading only 1/8th the normal batch size, to
  //  get computes computing while the next full batch is loaded.

  overlapsMax /= 8;

  if (ovlStore)
    *overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax, false);
  if (ovlFile)
    *overlapsLen = ovlFile->readOverlaps(overlaps, overlapsMax);

  overlapsMax *= 8;  //  Back to the normal batch size.

  fprintf(stderr, "Loaded %u overlaps.\n", *overlapsLen);

  rcache->loadReads(overlaps, *overlapsLen);

  //  Loop over all the overlaps.

  while (overlapsALen + overlapsBLen > 0) {

    //  Launch next batch of threads
    //fprintf(stderr, "LAUNCH THREADS\n");

    //  Globals, ugh.  These limit the threads to the range of overlaps we have loaded.  Each thread
    //  will pull out THREAD_SIZE overlaps at a time to compute, updating batchPosID as it does so.
    //  Each thread will stop when batchPosID > batchEndID.

    batchPrtID =  0;
    batchPosID =  0;
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
        outStore->writeOverlap(overlaps + oo);
    if (ovlFile)
      outFile->writeOverlaps(overlaps, *overlapsLen);

    //  Load more overlaps

    if (ovlStore)
      *overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax, false);
    if (ovlFile)
      *overlapsLen = ovlFile->readOverlaps(overlaps, overlapsMax);

    fprintf(stderr, "Loaded %u overlaps.\n", *overlapsLen);

    rcache->loadReads(overlaps, *overlapsLen);

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

  //  Goodbye.

  delete    rcache;

  gkpStore->gkStore_close();

  delete    ovlStore;
  delete    outStore;

  delete    ovlFile;
  delete    outFile;

  delete [] overlapsA;
  delete [] overlapsB;

  delete [] WA;
  delete [] tID;

  return(0);
}



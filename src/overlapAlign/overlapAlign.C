
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
 *    Brian P. Walenz beginning on 2019-APR-22
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "sweatShop.H"

#include "system.H"
#include "sequence.H"

#include <pthread.h>

#include "sqStore.H"
#include "sqCache.H"
#include "ovStore.H"

#include "edlib.H"




#include "alignStats.H"
#include "overlapAlign-globalData.H"
#include "overlapAlign-threadData.H"
#include "overlapAlign-computation.H"


void
computeOverlapAlignment(ovOverlap   *ovl,
                        char        *aseq, int32 alen,
                        char        *bseq, int32 blen,
                        uint32       minOverlapLength,
                        double       maxErate,
                        bool         partialOverlaps,
                        alignStats  &localStats);




void *
overlapReader(void *G) {
  trGlobalData     *g = (trGlobalData  *)G;
  maComputation    *s = NULL;

  if (g->ovlStore) {
    while ((g->curID <= g->endID) &&
           (g->ovlStore->numOverlaps(g->curID) == 0))
      g->curID++;

    if (g->curID <= g->endID)
      s = new maComputation(g->curID++, g->seqCache, g->ovlStore);
  }

  //if (g->ovlFile) {
  //  s = new maComputation(g->ovlFile);
  //}

  return(s);
}



void
overlapWriter(void *G, void *S) {
  trGlobalData     *g = (trGlobalData  *)G;
  maComputation    *s = (maComputation *)S;

#if 1
  if (g->outStore)
    for (uint64 oo=0; oo<s->_overlapsLen; oo++)
      g->outStore->writeOverlap(s->_overlaps + oo);

  if (g->outFile)
    g->outFile->writeOverlaps(s->_overlaps, s->_overlapsLen);
#endif

  delete s;
}




void
maComputation::computeAlignments(uint32  minOverlapLength,
                                 double  maxErate) {
  alignStats  localStats;

  if (_overlapsLen == 0)
    return;

  //  Step 1:  Decide which overlaps to compute.  We don't care about overlaps that are:
  //               well contained in the read
  //               thin

  uint32  aID       = _aID;
  int32   alen      = _seqCache->sqCache_getLength(aID);

  int32   minLength  = 500;         //  Overlap seeds below this length are discarded.
  int32   maxEdge    = 250;         //  Overlap seeds contained by this length on both sides are discarded.

  uint32  nShort     = 0;
  uint32  nContained = 0;
  uint32  nUnaligned = 0;
  uint32  nThin      = 0;

  //  Decide if I am well contained in some other read.  If so, filter all overlaps.

  bool    wellContained = false;

  for (uint32 ii=0; ii<_overlapsLen; ii++) {
    ovOverlap *ov = _overlaps + ii;

    int32   a5        = (int32)ov->dat.ovl.ahg5;
    int32   a3        = (int32)ov->dat.ovl.ahg3;

    int32   b5        = (int32)ov->dat.ovl.bhg5;
    int32   b3        = (int32)ov->dat.ovl.bhg3;

    //  Discard if the implied overlap is well contained.
    //
    //                a5          a3
    //               ----|------|----
    //      -------------|------|---------
    //            b5                b3
    //
    if ((b5 - a5 > maxEdge) && (b3 - a3 > maxEdge)) {
      fprintf(stderr, "contained in b %u\n", ov->b_iid);
      wellContained = true;
    }
  }

#if 0
  if (wellContained) {
    fprintf(stderr, "read %8u well contained\n", aID);

    for (uint32 ii=0; ii<_overlapsLen; ii++)
      _overlaps[ii].evalue(AS_MAX_EVALUE);

    return;
  }
#endif

  //  Decide on a plausible thickest overlap on each end.  Any overlap smaller than this
  //  will be filtered.
  //
  //  The threshold length is picked to allow X * coverage overlaps per end.

  int32    thick5len = 0;
  int32   *thick5arr = new int32 [_overlapsLen];

  int32    thick3len = 0;
  int32   *thick3arr = new int32 [_overlapsLen];

  for (uint32 ii=0; ii<_overlapsLen; ii++) {
    ovOverlap *ov = _overlaps + ii;

    int32   a5        = (int32)ov->dat.ovl.ahg5;
    int32   a3        = (int32)ov->dat.ovl.ahg3;
    int32   abgn      =        a5;
    int32   aend      = alen - a3;

    int32   b5        = (int32)ov->dat.ovl.bhg5;
    int32   b3        = (int32)ov->dat.ovl.bhg3;
    uint32  bID       = ov->b_iid;
    int32   blen      = _seqCache->sqCache_getLength(bID);
    int32   bbgn      =        b5;
    int32   bend      = blen - b3;

    //    --------------|-------------------|-----
    //              ----|-------------------|---------------
    //
    if ((a5 > b5) &&
        (a3 < b3)) {
      thick3arr[thick3len++] = alen - a5 + b5;
      fprintf(stderr, "overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d thick3 %d\n", ii, aID, abgn, aend, bID, bbgn, bend, alen - a5 + b5);
    }

    //              ----|-------------------|---------------
    //    --------------|-------------------|-----
    //
    if ((a5 < b5) &&
        (a3 > b3)) {
      thick5arr[thick5len++] = alen - a3 + b3;
      fprintf(stderr, "overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d thick5 %d\n", ii, aID, abgn, aend, bID, bbgn, bend, alen - a3 + b3);
    }
  }

#warning hardcoded coverage
  double   coverage = 30;
  int32    novl     = (int32)ceil(0.50 * coverage);

  sort(thick5arr, thick5arr + thick5len);
  sort(thick3arr, thick3arr + thick3len);

  int32    thick3 = 0;
  int32    thick5 = 0;

  if (novl < thick5len)
    thick5 = thick5arr[thick5len - 1 - novl];

  if (novl < thick3len)
    thick3 = thick3arr[thick3len - 1 - novl];

  delete [] thick5arr;
  delete [] thick3arr;

  fprintf(stderr, "read %8u thick5 %d (out of %d)\n", aID, thick5, thick5len);
  fprintf(stderr, "              thick3 %d (out of %d) novl %d\n", thick3, thick3len, novl);





  //  Filter individual overlaps.

  for (uint32 ii=0; ii<_overlapsLen; ii++) {
    ovOverlap *ov = _overlaps + ii;

    int32   a5        = (int32)ov->dat.ovl.ahg5;
    int32   a3        = (int32)ov->dat.ovl.ahg3;
    int32   abgn      =        a5;
    int32   aend      = alen - a3;

    int32   b5        = (int32)ov->dat.ovl.bhg5;
    int32   b3        = (int32)ov->dat.ovl.bhg3;
    uint32  bID       = ov->b_iid;
    int32   blen      = _seqCache->sqCache_getLength(bID);
    int32   bbgn      =        b5;
    int32   bend      = blen - b3;

    //  Discard if the overlap seed is too small, regardless of the implied overlap length.
    //
    if ((aend - abgn < minLength) ||
        (bend - bbgn < minLength)) {
      nShort++;
      fprintf(stderr, "overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d short.\n", ii, aID, abgn, aend, bID, bbgn, bend);
      ov->evalue(AS_MAX_EVALUE);
      continue;
    }

    //  Discard if the implied overlap is well contained.
    //
    //            a5                a3
    //      -------------|------|---------
    //               ----|------|----
    //                b5          b3
    //
    if ((a5 - b5 > maxEdge) && (a3 - b3 > maxEdge)) {
      nContained++;
      fprintf(stderr, "overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d contained.\n", ii, aID, abgn, aend, bID, bbgn, bend);
      ov->evalue(AS_MAX_EVALUE);
      continue;
    }


    //
    //  Discard if the overlap seed is a thin dovetail.
    //

    //    --------------|-------------------|-----
    //              ----|-------------------|---------------
    //
    if ((a5 > b5) &&
        (a3 < b3) &&
        (alen - a5 + b5 < thick3)) {
      nThin++;
      fprintf(stderr, "overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d thick3 %d < %d\n", ii, aID, abgn, aend, bID, bbgn, bend, alen - a5 + b5, thick3);
      ov->evalue(AS_MAX_EVALUE);
      continue;
    }

    //              ----|-------------------|---------------
    //    --------------|-------------------|-----
    //
    if ((a5 < b5) &&
        (a3 > b3) &&
        (alen - a3 + b3 < thick5)) {
      nThin++;
      fprintf(stderr, "overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d thick5 %d < %d\n", ii, aID, abgn, aend, bID, bbgn, bend, alen - a3 + b3, thick5);
      ov->evalue(AS_MAX_EVALUE);
      continue;
    }


    fprintf(stderr, "overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d\n", ii, aID, abgn, aend, bID, bbgn, bend);
  }

  fprintf(stderr, "read %8u nShort %5u nContained %5u nUnaligned %5u nThin %5u filtered %.2f%%\n",
          _aID, nShort, nContained, nUnaligned, nThin, 100.0 * (nShort + nContained + nUnaligned + nThin) / _overlapsLen);

  //  Step 1:  Find, precisely, the optimal overlap for each read.
  //  No alignments, just end points.  The ovOverlap is updated
  //  with correct end points.

  _seqCache->sqCache_getSequence(_aID, _aRead, _aLen, _aMax);

  for (uint32 ii=0; ii<_overlapsLen; ii++) {
    ovOverlap *ov = _overlaps + ii;

    //  Skip the overlap if it was filtered already.
    if (ov->evalue() == AS_MAX_EVALUE)
      continue;

    //  Load the b read.
    fprintf(stderr, "Processing overlap #%u to read %u with minOlap %u maxErate %f\n", ii, ov->b_iid, minOverlapLength, maxErate);
    _seqCache->sqCache_getSequence(ov->b_iid, _bRead, _bLen, _bMax);

    //  Reverse complement, if needed (should be part of sqCache, really).
    if (ov->flipped() == true)
      reverseComplementSequence(_bRead, _bLen);

    // 
    computeOverlapAlignment(ov,
                            _aRead, _aLen,
                            _bRead, _bLen,
                            minOverlapLength,
                            maxErate,
                            true,
                            localStats);
  }

  //  Step 2:  Transfer to a tig.  generateCorrectionLayouts.C generateLayout().

  //  Step 3:  Align all reads in tig to tig sequence.  Save alignments
  //  in conveninent multialign structure.  This should be part of tgTig.

  //  Step 4:  

}


void
overlapRecompute(void *G, void *T, void *S) {
  trGlobalData     *g = (trGlobalData  *)G;
  maThreadData     *t = (maThreadData  *)T;
  maComputation    *s = (maComputation *)S;

  //fprintf(stderr, "Processing read %u with %u overlaps.\n", s->_aID, s->_overlapsLen);

  s->computeAlignments(g->minOverlapLength, g->maxErate);
};



int
main(int argc, char **argv) {
  trGlobalData   *g = new trGlobalData;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0)
      g->seqStoreName = argv[++arg];

    else if (strcmp(argv[arg], "-O") == 0)
      g->ovlStoreName = argv[++arg];

    else if (strcmp(argv[arg], "-o") == 0)
      g->outStoreName = argv[++arg];

    else if (strcmp(argv[arg], "-r") == 0)
      decodeRange(argv[++arg], g->bgnID, g->endID);

    else if (strcmp(argv[arg], "-t") == 0)
      g->numThreads = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-erate") == 0)
      g->maxErate = atof(argv[++arg]);

    else if (strcmp(argv[arg], "-memory") == 0)
      g->memLimit = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-len") == 0)
      g->minOverlapLength = atoi(argv[++arg]);

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (g->seqStoreName == NULL)   err.push_back("No sequence store (-S option) supplied.\n");
  if (g->ovlStoreName == NULL)   err.push_back("No overlap store (-O option) supplied.\n");
  if (g->outStoreName == NULL)   err.push_back("No output store (-o option) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -S seqStore       Mandatory, path to seqStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Inputs can come from either a store or a file.\n");
    fprintf(stderr, "  -O ovlStore       \n");
    fprintf(stderr, "  -O ovlFile        \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "If from an ovlStore, the range of reads processed can be restricted.\n");
    fprintf(stderr, "  -r bgnID[-endID]  process reads bgnID to endID, inclusive\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Outputs will be written to a store or file, depending on the input type\n");
    fprintf(stderr, "  -o ovlStore       \n");
    fprintf(stderr, "  -o ovlFile        \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -erate e          Overlaps are computed at 'e' fraction error; must be larger than the original erate\n");
    fprintf(stderr, "  -partial          Overlaps are 'overlapInCore -S' partial overlaps\n");
    fprintf(stderr, "  -memory m         Use up to 'm' GB of memory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t n              Use up to 'n' cores\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Advanced options:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  g->initialize();



#if 0

  maThreadData  *t = new maThreadData(g, 0);

  while (1) {
    maComputation *c = (maComputation *)overlapReader(g);

    if (c == NULL)
      break;

    overlapRecompute(g, t, c);
    overlapWriter(g, c);
  }

  delete t;

#else

  maThreadData **td = new maThreadData * [g->numThreads];
  sweatShop     *ss = new sweatShop(overlapReader, overlapRecompute, overlapWriter);

  ss->setLoaderQueueSize(128);
  ss->setWriterQueueSize(128);

  ss->setNumberOfWorkers(g->numThreads);

  for (uint32 w=0; w<g->numThreads; w++)
    ss->setThreadData(w, td[w] = new maThreadData(g, w));  //  these leak

  ss->run(g, false);

  delete ss;

  for (uint32 w=0; w<g->numThreads; w++)
    delete td[w];

#endif

  delete g;

  fprintf(stderr, "\nSuccess!  Bye.\n");

  return(0);
}


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

#include "system.H"
#include "sequence.H"

#include "sqStore.H"
#include "sqCache.H"
#include "ovStore.H"

#include "alignStats.H"
#include "overlapAlign-globalData.H"
#include "overlapAlign-threadData.H"
#include "overlapAlign-computation.H"

#include <tuple>



void
computeOverlapAlignment(ovOverlap   *ovl,
                        char        *aseq, int32 alen,
                        char        *bseq, int32 blen,
                        uint32       minOverlapLength,
                        double       maxErate,
                        bool         partialOverlaps,
                        alignStats  &localStats);



//  Returns true if the read is 'well' contained in ANY implied overlap.
//
//                a5          a3
//               ----|------|----
//      -------------|------|---------
//            b5                b3
//
bool
maComputation::isWellContained(int32 maxEdge) {
  bool   wellContained = false;

  for (uint32 ii=0; ii<_overlapsLen; ii++) {
    ovOverlap *ov = _overlaps + ii;

    int32   a5        = (int32)ov->dat.ovl.ahg5;
    int32   a3        = (int32)ov->dat.ovl.ahg3;

    int32   b5        = (int32)ov->dat.ovl.bhg5;
    int32   b3        = (int32)ov->dat.ovl.bhg3;

    if ((b5 - a5 > maxEdge) && (b3 - a3 > maxEdge)) {
      fprintf(stderr, "contained in b %u\n", ov->b_iid);
      wellContained = true;
    }
  }

  return(wellContained);
}



//  Return the length of the shortest overlap we care about
//  of the 5' and 3' ends.
//
tuple<int32, int32>
maComputation::findMinThickestEdge(double coverage) {
  int32    thick5len = 0,  *thick5arr = new int32 [_overlapsLen];
  int32    thick3len = 0,  *thick3arr = new int32 [_overlapsLen];

  for (uint32 ii=0; ii<_overlapsLen; ii++) {
    ovOverlap *ov = _overlaps + ii;

    int32   alen  = _seqCache->sqCache_getLength(_aID);

    int32   a5    = (int32)ov->dat.ovl.ahg5;
    int32   a3    = (int32)ov->dat.ovl.ahg3;
    int32   b5    = (int32)ov->dat.ovl.bhg5;
    int32   b3    = (int32)ov->dat.ovl.bhg3;

    //    --------------|-------------------|-----
    //              ----|-------------------|---------------
    //
    if ((a5 > b5) &&
        (a3 < b3))
      thick3arr[thick3len++] = alen - a5 + b5;

    //              ----|-------------------|---------------
    //    --------------|-------------------|-----
    //
    if ((a5 < b5) &&
        (a3 > b3))
      thick5arr[thick5len++] = alen - a3 + b3;
  }

  sort(thick5arr, thick5arr + thick5len);
  sort(thick3arr, thick3arr + thick3len);

  int32     novl   = (int32)ceil(0.50 * coverage);
  uint32    thick3 = 0;
  uint32    thick5 = 0;

  if (novl < thick5len)
    thick5 = thick5arr[thick5len - 1 - novl];

  if (novl < thick3len)
    thick3 = thick3arr[thick3len - 1 - novl];

  delete [] thick5arr;
  delete [] thick3arr;

  //fprintf(stderr, "read %8u thick5 %d (out of %d)\n", aID, thick5, thick5len);
  //fprintf(stderr, "              thick3 %d (out of %d) novl %d\n", thick3, thick3len, novl);

  return(make_tuple(thick5, thick3));
}


void
maComputation::computeAlignments(uint32  minOverlapLength,
                                 double  maxErate) {
  alignStats  localStats;

  uint32  aID        = _aID;
  int32   alen       = _seqCache->sqCache_getLength(aID);

  //  Set parameters.

#warning hardcoded coverage
  double   coverage  = 30;

  int32   minLength  = 500;         //  Overlap seeds below this length are discarded.
  int32   maxEdge    = 2500;        //  Overlap seeds contained by this length on both sides are discarded.

  uint32  nShort     = 0;
  uint32  nContained = 0;
  uint32  nUnaligned = 0;
  uint32  nThin      = 0;

  //  If there are no overlaps, there are no overlaps to align and output.

  if (_overlapsLen == 0)
    return;

  //  If this read is 'well' contained in some other read, assume it's
  //  useless for assembly and throw out all of its overlaps.

  if (isWellContained(maxEdge) == true) {
    fprintf(stderr, "read %8u well contained\n", _aID);

    for (uint32 ii=0; ii<_overlapsLen; ii++)
      _overlaps[ii].evalue(AS_MAX_EVALUE);

    return;
  }

  //  Decide which overlaps to compute alignments for:  only the thickest on each side.

  uint32   thick5;
  uint32   thick3;

  tie(thick5, thick3) = findMinThickestEdge(coverage);

  //  Examine each overlap.  Either deicde it isn't useful for assembly, or compute, precisely,
  //  the overlap end points and save the alignment.

  _seqCache->sqCache_getSequence(_aID, _aRead, _aLen, _aMax);

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

    //  Discard if the implied overlapping read is well contained in us.
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

    //  Discard if the overlap seed is a thin dovetail on the 3' end.
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

    //  Discard if the overlap seed is a thin dovetail on the 5' end.
    //
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

    //  A good overlap to process!
    //
    //  Find, precisely, the optimal overlap for each read.  No alignments,
    //  just end points.  The ovOverlap is updated with correct end points.

    fprintf(stderr, "overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d\n", ii, aID, abgn, aend, bID, bbgn, bend);

    //  Load the b read.
    _seqCache->sqCache_getSequence(ov->b_iid, _bRead, _bLen, _bMax);

    //  Reverse complement, if needed (should be part of sqCache, really).
    if (ov->flipped() == true)
      reverseComplementSequence(_bRead, _bLen);

    //  Compute the alignment.
    computeOverlapAlignment(ov,
                            _aRead, _aLen,
                            _bRead, _bLen,
                            minOverlapLength,
                            maxErate,
                            true,
                            localStats);
  }

  fprintf(stderr, "read %8u nShort %5u nContained %5u nUnaligned %5u nThin %5u filtered %.2f%%\n",
          _aID, nShort, nContained, nUnaligned, nThin, 100.0 * (nShort + nContained + nUnaligned + nThin) / _overlapsLen);

  //  Step 2:  Transfer to a tig.  generateCorrectionLayouts.C generateLayout().

  //  Step 3:  Align all reads in tig to tig sequence.  Save alignments
  //  in conveninent multialign structure.  This should be part of tgTig.

  //  Step 4:  

}

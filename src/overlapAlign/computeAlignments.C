
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



void
maComputation::fetchUntrimmedRead(uint32 id_,
                                  bool   isA_,
                                  bool   revComp_) {
  uint32  &id   = ((isA_ == true) ? _aID    : _bID);
  char   *&read = ((isA_ == true) ? _aRead  : _bRead);
  uint32  &max  = ((isA_ == true) ? _aMax   : _bMax);
  uint32   len  = 0;

  id = id_;

  //  Fetch the read from the store.
  _seqCache->sqCache_getSequence(id, read, len, max);

  //  Make sure lengths agree.
  assert(len == _readData[id].rawLength);
  assert(len == _seqCache->sqCache_getLength(id));

  //  Reverse complement the trimmed read.
  if (revComp_)
    reverseComplementSequence(read, _readData[id].rawLength);
}



void
maComputation::fetchTrimmedRead(uint32 id_,
                                bool   isA_,
                                bool   revComp_) {
  uint32  &id   = ((isA_ == true) ? _aID    : _bID);
  char   *&read = ((isA_ == true) ? _aRead  : _bRead);
  uint32  &max  = ((isA_ == true) ? _aMax   : _bMax);
  uint32   len  = 0;

  id = id_;

  //  Fetch the read from the store.
  _seqCache->sqCache_getSequence(id, read, len, max);

  //  Make sure lengths agree.
  assert(len == _readData[id].rawLength);
  assert(len == _seqCache->sqCache_getLength(id));

  //  Trim the read.  If the clrBgn is non-zero, shift all the bases
  //  to the left.

  if (_readData[id].clrBgn > 0)
    memmove(read, read + _readData[id].clrBgn, _readData[id].trimmedLength);

  read[_readData[id].trimmedLength] = 0;  //  maComputation allocates one extra byte for each read.

  //  Reverse complement the trimmed read.
  if (revComp_)
    reverseComplementSequence(read, _readData[id].trimmedLength);
}



bool
maComputation::trimOverlap_Normal(ovOverlap *ovl) {

  int32  aclrbgn = _readData[_aID].clrBgn;
  int32  aclrend = _readData[_aID].clrEnd;
  int32  aovlbgn =                             ovl->dat.ovl.ahg5;
  int32  aovlend = _readData[_aID].rawLength - ovl->dat.ovl.ahg3;
  int32  alen    = aovlend - aovlbgn;

  int32  bclrbgn = _readData[_bID].clrBgn;
  int32  bclrend = _readData[_bID].clrEnd;
  int32  bovlbgn =                             ovl->dat.ovl.bhg5;
  int32  bovlend = _readData[_bID].rawLength - ovl->dat.ovl.bhg3;
  int32  blen    = bovlend - bovlbgn;

  assert(aovlbgn < aovlend);
  assert(bovlbgn < bovlend);

#if 1
  fprintf(stderr, "trimOverlap_Normal()-- ORIG A %8u clear %6u-%-6u length %6u overlap %6u-%-6u\n", ovl->a_iid, aclrbgn, aclrend, _readData[_aID].rawLength, aovlbgn, aovlend);
  fprintf(stderr, "trimOverlap_Normal()--      B %8u clear %6u-%-6u length %6u overlap %6u-%-6u\n", ovl->b_iid, bclrbgn, bclrend, _readData[_bID].rawLength, bovlbgn, bovlend);
#endif

  //  If the overlap doesn't intersect the clear range, the overlap disappears.
  if ((aclrend <= aovlbgn) || (aovlend <= aclrbgn) ||
      (bclrend <= bovlbgn) || (bovlend <= bclrbgn)) {
    return(false);
  }

  double  afracbgn = (double)((aclrbgn < aovlbgn) ? (0) : (aclrbgn - aovlbgn)) / alen;
  double  afracend = (double)((aclrend > aovlend) ? (0) : (aovlend - aclrend)) / alen;

  double  bfracbgn = (double)((bclrbgn < bovlbgn) ? (0) : (bclrbgn - bovlbgn)) / blen;
  double  bfracend = (double)((bclrend > bovlend) ? (0) : (bovlend - bclrend)) / blen;

  double  maxbgn   = max(afracbgn, bfracbgn);
  double  maxend   = max(afracend, bfracend);

  assert(maxbgn < 1.0);
  assert(maxend < 1.0);

  ovl->dat.ovl.ahg5 += (uint32)round(maxbgn * alen);
  ovl->dat.ovl.ahg3 += (uint32)round(maxend * alen);

  ovl->dat.ovl.bhg5 += (uint32)round(maxbgn * blen);
  ovl->dat.ovl.bhg3 += (uint32)round(maxend * blen);

#if 1
  aovlbgn =                                 ovl->dat.ovl.ahg5;
  aovlend = _readData[_aID].trimmedLength - ovl->dat.ovl.ahg3;

  bovlbgn =                                 ovl->dat.ovl.bhg5;
  bovlend = _readData[_bID].trimmedLength - ovl->dat.ovl.bhg3;

  fprintf(stderr, "trimOverlap_Normal()-- TRIM A %8u                     length %6u overlap %6u-%-6u\n", ovl->a_iid, _readData[_aID].trimmedLength, aovlbgn, aovlend);
  fprintf(stderr, "trimOverlap_Normal()--      B %8u                     length %6u overlap %6u-%-6u\n", ovl->b_iid, _readData[_bID].trimmedLength, bovlbgn, bovlend);
#endif

  return(true);
}



bool
maComputation::trimOverlap_Flipped(ovOverlap *ovl) {

  int32  aclrbgn = _readData[_aID].clrBgn;
  int32  aclrend = _readData[_aID].clrEnd;
  int32  aovlbgn =                             ovl->dat.ovl.ahg5;
  int32  aovlend = _readData[_aID].rawLength - ovl->dat.ovl.ahg3;
  int32  alen    = aovlend - aovlbgn;

  int32  bclrbgn = _readData[_bID].rawLength - _readData[_bID].clrEnd;
  int32  bclrend = _readData[_bID].rawLength - _readData[_bID].clrBgn;
  int32  bovlbgn =                             ovl->dat.ovl.bhg5;
  int32  bovlend = _readData[_bID].rawLength - ovl->dat.ovl.bhg3;
  int32  blen    = bovlend - bovlbgn;

  assert(aovlbgn < aovlend);
  assert(bovlbgn < bovlend);

#if 1
  fprintf(stderr, "trimOverlap_Flipped()-- ORIG A %8u clear %6u-%-6u length %6u overlap %6u-%-6u\n", ovl->a_iid, aclrbgn, aclrend, _readData[_aID].rawLength, aovlbgn, aovlend);
  fprintf(stderr, "trimOverlap_Flipped()--      B %8u clear %6u-%-6u length %6u overlap %6u-%-6u\n", ovl->b_iid, bclrbgn, bclrend, _readData[_bID].rawLength, bovlbgn, bovlend);
#endif

  //  If the overlap doesn't intersect the clear range, the overlap disappears.
  if ((aclrend <= aovlbgn) || (aovlend <= aclrbgn) ||
      (bclrend <= bovlbgn) || (bovlend <= bclrbgn)) {
    return(false);
  }

  double  afracbgn = (double)((aclrbgn < aovlbgn) ? (0) : (aclrbgn - aovlbgn)) / alen;
  double  afracend = (double)((aclrend > aovlend) ? (0) : (aovlend - aclrend)) / alen;

  double  bfracbgn = (double)((bclrbgn < bovlbgn) ? (0) : (bclrbgn - bovlbgn)) / blen;
  double  bfracend = (double)((bclrend > bovlend) ? (0) : (bovlend - bclrend)) / blen;

  double  maxbgn   = max(afracbgn, bfracbgn);
  double  maxend   = max(afracend, bfracend);

  assert(maxbgn < 1.0);
  assert(maxend < 1.0);

  ovl->dat.ovl.ahg5 += (uint32)round(maxbgn * alen);
  ovl->dat.ovl.ahg3 += (uint32)round(maxend * alen);

  ovl->dat.ovl.bhg5 += (uint32)round(maxbgn * blen);
  ovl->dat.ovl.bhg3 += (uint32)round(maxend * blen);

#if 1
  aovlbgn =                                 ovl->dat.ovl.ahg5;
  aovlend = _readData[_aID].trimmedLength - ovl->dat.ovl.ahg3;

  bovlbgn =                                 ovl->dat.ovl.bhg5;
  bovlend = _readData[_bID].trimmedLength - ovl->dat.ovl.bhg3;

  fprintf(stderr, "trimOverlap_Flipped()-- TRIM A %8u                     length %6u overlap %6u-%-6u\n", ovl->a_iid, _readData[_aID].trimmedLength, aovlbgn, aovlend);
  fprintf(stderr, "trimOverlap_Flipped()--      B %8u                     length %6u overlap %6u-%-6u\n", ovl->b_iid, _readData[_bID].trimmedLength, bovlbgn, bovlend);
#endif

  return(true);
}



bool
maComputation::trimOverlap(ovOverlap *ovl) {
  bool  success;

  assert(_aID == ovl->a_iid);
  assert(_bID == ovl->b_iid);

  if (ovl->flipped() == true)
    success = trimOverlap_Normal(ovl);
  else
    success = trimOverlap_Flipped(ovl);

  return(success);
}





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

    int32   alen  = _readData[_aID].trimmedLength;

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

  fprintf(stderr, "read %8u thick5 %d (out of %d)\n", _aID, thick5, thick5len);
  fprintf(stderr, "              thick3 %d (out of %d) novl %d\n", thick3, thick3len, novl);

  return(make_tuple(thick5, thick3));
}


void
maComputation::computeAlignments(uint32  minOverlapLength,
                                 double  maxErate) {
  alignStats  localStats;

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

  fprintf(stderr, "computeAlignments()-- \n");
  fprintf(stderr, "computeAlignments()-- ========================================BEGIN");

  //  If this read is 'well' contained in some other read, assume it's
  //  useless for assembly and throw out all of its overlaps.

  if (isWellContained(maxEdge) == true) {
    fprintf(stderr, "computeAlignments()-- read %8u well contained\n", _aID);

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

  fetchTrimmedRead(_aID);

  for (uint32 ii=0; ii<_overlapsLen; ii++) {
    ovOverlap *ov = _overlaps + ii;

    fprintf(stderr, "computeAlignments()-- \n");
    fprintf(stderr, "computeAlignments()-- ----------------------------------------OVERLAP\n");

    assert(_aID == ov->a_iid);
    _bID = ov->b_iid;

    //  Adjust the overlap for the trimming we have already done.

    if (trimOverlap(ov) == false) {
      fprintf(stderr, "computeAlignments()-- Overlap trimmed out.\n");
      ov->evalue(AS_MAX_EVALUE);
      continue;
    }

    //  Harvest overlap parameters.  We do this OUTSIDE of ovOverlap so we
    //  can use the TRIMMED length of each read, which ovOverlap (and
    //  sqStore) don't know about yet.

    int32   a5        = (int32)ov->dat.ovl.ahg5;
    int32   a3        = (int32)ov->dat.ovl.ahg3;
    int32   abgn      =                                 a5;
    int32   aend      = _readData[_aID].trimmedLength - a3;

    int32   b5        = (int32)ov->dat.ovl.bhg5;
    int32   b3        = (int32)ov->dat.ovl.bhg3;
    int32   bbgn      =                                 b5;
    int32   bend      = _readData[_bID].trimmedLength - b3;

    //  Discard if the overlap seed is too small, regardless of the implied overlap length.
    //
    if ((aend - abgn < minLength) ||
        (bend - bbgn < minLength)) {
      nShort++;
      fprintf(stderr, "computeAlignments()-- overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d short.\n", ii, _aID, abgn, aend, _bID, bbgn, bend);
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
      fprintf(stderr, "computeAlignments()-- overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d contained.\n", ii, _aID, abgn, aend, _bID, bbgn, bend);
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
        (_readData[_aID].trimmedLength - a5 + b5 < thick3)) {
      nThin++;
      fprintf(stderr, "computeAlignments()-- overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d thick3 %d < %d\n", ii, _aID, abgn, aend, _bID, bbgn, bend, _readData[_aID].trimmedLength - a5 + b5, thick3);
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
        (_readData[_aID].trimmedLength - a3 + b3 < thick5)) {
      nThin++;
      fprintf(stderr, "computeAlignments()-- overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d thick5 %d < %d\n", ii, _aID, abgn, aend, _bID, bbgn, bend, _readData[_aID].trimmedLength - a3 + b3, thick5);
      ov->evalue(AS_MAX_EVALUE);
      continue;
    }

    //  A good overlap to process!
    //
    //  Find, precisely, the optimal overlap for each read.  No alignments,
    //  just end points.  The ovOverlap is updated with correct end points.

    //  Reported during computeOverlapAlignment()
    //fprintf(stderr, "computeAlignments()-- overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d valid\n", ii, _aID, abgn, aend, _bID, bbgn, bend);

    //  Load the b read and compute the alignment.
    fetchTrimmedRead(_bID, false, ov->flipped());

    computeOverlapAlignment(ov,
                            _aRead, _readData[_aID].trimmedLength,
                            _bRead, _readData[_bID].trimmedLength,
                            minOverlapLength,
                            maxErate,
                            true,
                            localStats);
  }

  fprintf(stderr, "computeAlignments()-- read %8u nShort %5u nContained %5u nUnaligned %5u nThin %5u filtered %.2f%%\n",
          _aID, nShort, nContained, nUnaligned, nThin, 100.0 * (nShort + nContained + nUnaligned + nThin) / _overlapsLen);

  //  Step 2:  Transfer to a tig.  generateCorrectionLayouts.C generateLayout().

  //  Step 3:  Align all reads in tig to tig sequence.  Save alignments
  //  in conveninent multialign structure.  This should be part of tgTig.

  //  Step 4:  

}

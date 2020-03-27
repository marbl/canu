
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
#include "sequence.H"

#include "sqStore.H"
#include "sqCache.H"
#include "ovStore.H"

#include "overlapAlign-globalData.H"
#include "overlapAlign-threadData.H"
#include "overlapAlign-computation.H"



void
maComputation::fetchUntrimmedRead(uint32 id_,
                                  bool   isA_,
                                  bool   revComp_) {
  uint32  &id   = ((isA_ == true) ? _aID    : _bID);
  char   *&read = ((isA_ == true) ? _aRead  : _bRead);
  uint32  &max  = ((isA_ == true) ? _aMax   : _bMax);
  uint32   len  = 0;

  id = id_;

  //fprintf(stderr, "Fetch untrimmed %c read %c %u\n",
  //        (revComp_ == true) ? 'R' : 'F',
  //        (isA_ == true) ? 'A' : 'B',
  //        id_);

  //  Fetch the read from the store.
  _seqCache->sqCache_getSequence(id, read, len, max);

  //  Make sure lengths agree.
  assert(len == _readData[id].rawLength);
  assert(len == _seqCache->sqCache_getLength(id));

  //  Reverse complement the read.
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

  //fprintf(stderr, "Fetch trimmed %c read %c %u coords %d-%d out of untrimmed length %d\n",
  //        (revComp_ == true) ? 'R' : 'F',
  //        (isA_ == true) ? 'A' : 'B',
  //        id_,
  //        _readData[id].clrBgn, _readData[id].clrEnd,
  //        _readData[id].rawLength);

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

  //  Reverse complement the read.
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

  if (_verboseAlign > 0) {
    fprintf(stderr, "\n");
    fprintf(stderr, "trimOverlap_Normal()-- ORIG A %8u clear %6d-%-6d length %6u overlap %6d-%-6d\n", ovl->a_iid, aclrbgn, aclrend, _readData[_aID].rawLength, aovlbgn, aovlend);
    fprintf(stderr, "trimOverlap_Normal()--      B %8u clear %6d-%-6d length %6u overlap %6d-%-6d\n", ovl->b_iid, bclrbgn, bclrend, _readData[_bID].rawLength, bovlbgn, bovlend);
  }

  //  If the overlap doesn't intersect the clear range, the overlap disappears.
  if ((aclrend <= aovlbgn) || (aovlend <= aclrbgn) ||
      (bclrend <= bovlbgn) || (bovlend <= bclrbgn)) {
    if (_verboseAlign > 0)
      fprintf(stderr, "trimOverlap_Normal()--      No intersection with clear range.\n");
    return(false);
  }

  //  Compute how much we need to shift the overlap to make it fall within both clear ranges.

  int32  ashift5 = (double)((aclrbgn < aovlbgn) ? (0) : (aclrbgn - aovlbgn));
  int32  ashift3 = (double)((aclrend > aovlend) ? (0) : (aovlend - aclrend));

  int32  bshift5 = (double)((bclrbgn < bovlbgn) ? (0) : (bclrbgn - bovlbgn));
  int32  bshift3 = (double)((bclrend > bovlend) ? (0) : (bovlend - bclrend));

  assert(ashift5 >= 0);
  assert(ashift3 >= 0);
  assert(bshift5 >= 0);
  assert(bshift3 >= 0);

  //  But the reads can contain different gappiness, so we need to scale the shift differently for each read.

  double abscale   = (double)(aovlend - aovlbgn) / (bovlend - bovlbgn);   //  One base in B is this many bases in A.
  double bascale   = (double)(bovlend - bovlbgn) / (aovlend - aovlbgn);   //  One base in B is this many bases in A.

  //  Shift the overlap coordinates, in raw read space.

  if (ashift5 < bshift5) {
    aovlbgn += bshift5 * abscale;
    bovlbgn += bshift5;
  } else {
    aovlbgn += ashift5;
    bovlbgn += ashift5 * bascale;
  }


  if (ashift3 < bshift3) {
    aovlend -= bshift3 * abscale;
    bovlend -= bshift3;
  } else {
    aovlend -= ashift3;
    bovlend -= ashift3 * bascale;
  }

  //  Adjust coordinates for the begin clear range.

  aovlbgn -= aclrbgn;
  aovlend -= aclrbgn;

  bovlbgn -= bclrbgn;
  bovlend -= bclrbgn;

  //  Report the overlap.

  if (_verboseAlign > 0) {
    fprintf(stderr, "trimOverlap_Normal()-- TRIM A %8u                     length %6u overlap %6d-%-6d\n", ovl->a_iid, _readData[_aID].trimmedLength, aovlbgn, aovlend);
    fprintf(stderr, "trimOverlap_Normal()--      B %8u                     length %6u overlap %6d-%-6d\n", ovl->b_iid, _readData[_bID].trimmedLength, bovlbgn, bovlend);
  }

  //  If the resulting overlap is too small (or negative!), the overlap cannot be formed from the trimmed reads.

  if ((aovlend - aovlbgn < 1000) ||
      (bovlend - bovlbgn < 1000)) {
    if (_verboseAlign > 0)
      fprintf(stderr, "trimOverlap_Normal()--      Too short after adjusting.\n");
    return(false);
  }

  //  Update the overlap with the new hangs.

  ovl->dat.ovl.ahg5 =                                 aovlbgn;
  ovl->dat.ovl.ahg3 = _readData[_aID].trimmedLength - aovlend;

  ovl->dat.ovl.bhg5 =                                 bovlbgn;
  ovl->dat.ovl.bhg3 = _readData[_bID].trimmedLength - bovlend;

  return(true);
}



bool
maComputation::trimOverlap_Flipped(ovOverlap *ovl) {

  int32  aclrbgn = _readData[_aID].clrBgn;
  int32  aclrend = _readData[_aID].clrEnd;
  int32  aovl5   =                             ovl->dat.ovl.ahg5;
  int32  aovl3   = _readData[_aID].rawLength - ovl->dat.ovl.ahg3;
  int32  alen    = aovl3 - aovl5;

  int32  bclrbgn = _readData[_bID].clrBgn;
  int32  bclrend = _readData[_bID].clrEnd;
  int32  bovl5   = _readData[_bID].rawLength - ovl->dat.ovl.bhg5;   //  This is the high coordinate, on the left  end of the overlap.
  int32  bovl3   =                             ovl->dat.ovl.bhg3;   //  This is the low coordinate,  on the right end of the overlap.
  int32  blen    = bovl5 - bovl3;

  assert(aovl5 < aovl3);
  assert(bovl3 < bovl5);

  if (_verboseAlign > 0) {
    fprintf(stderr, "\n");
    fprintf(stderr, "trimOverlap_Flipped()-- ORIG A %8u clear %6d-%-6d length %6u overlap %6d-%-6d\n", ovl->a_iid, aclrbgn, aclrend, _readData[_aID].rawLength, aovl5, aovl3);
    fprintf(stderr, "trimOverlap_Flipped()--      B %8u clear %6d-%-6d length %6u overlap %6d-%-6d\n", ovl->b_iid, bclrbgn, bclrend, _readData[_bID].rawLength, bovl5, bovl3);
  }

  //  If the overlap doesn't intersect the clear range, the overlap disappears.
  if ((aclrend <= aovl5) || (aovl3 <= aclrbgn) ||
      (bclrend <= bovl3) || (bovl5 <= bclrbgn)) {
    if (_verboseAlign > 0)
      fprintf(stderr, "trimOverlap_Flipped()--      No intersection with clear range.\n");
    return(false);
  }

  //  Compute how much we need to shift the overlap to make it fall within both clear ranges.
  //
  //  Coordinates on the B read are confusing.
  //
  //            | overlap |
  //  ----[-----|---------|----------]-->
  //          <-|---[-----|------------------------------------]-------
  //            ^   ^     ^                                    ^
  //            |   |     +bovl3                             +bclrbgn (smallest number)
  //            |   +bclrend
  //            +bovl5 (biggest number)
  //
  //  'bgn' and 'end' refer to coordinates in the read.
  //  '5'   and '3'   refer to 'left' and 'right' end of the overlap.

  int32  ashift5 = (double)((aclrbgn < aovl5) ? (0) : (aclrbgn - aovl5));
  int32  ashift3 = (double)((aclrend > aovl3) ? (0) : (aovl3 - aclrend));

  int32  bshift5 = (double)((bovl5 < bclrend) ? (0) : (bovl5 - bclrend));
  int32  bshift3 = (double)((bovl3 > bclrbgn) ? (0) : (bclrbgn - bovl3));

  //fprintf(stderr, "trimOverlap_Flipped()--      ashift %d %d  bshift %d %d\n", ashift5, ashift3, bshift5, bshift3);

  assert(ashift5 >= 0);
  assert(ashift3 >= 0);
  assert(bshift5 >= 0);
  assert(bshift3 >= 0);

  //  But the reads can contain different gappiness, so we need to scale the shift differently for each read.

  double abscale   = (double)(aovl3 - aovl5) / (bovl5 - bovl3);   //  One base in B is this many bases in A.
  double bascale   = (double)(bovl5 - bovl3) / (aovl3 - aovl5);   //  One base in B is this many bases in A.

  //  Shift the overlap coordinates, using the largest shift on each end.
  //  Scale the shift on the other read in an attempt to adjust for differing
  //  gappiness.

  if (ashift5 < bshift5) {
    aovl5 += bshift5 * abscale;
    bovl5 -= bshift5;
  } else {
    aovl5 += ashift5;
    bovl5 -= ashift5 * bascale;
  }


  if (ashift3 < bshift3) {
    aovl3 -= bshift3 * abscale;
    bovl3 += bshift3;
  } else {
    aovl3 -= ashift3;
    bovl3 += ashift3 * bascale;
  }

  //  Adjust coordinates for the begin clear range.

  aovl5 -= aclrbgn;
  aovl3 -= aclrbgn;

  bovl3 -= bclrbgn;
  bovl5 -= bclrbgn;

  //  Report the overlap.

  if (_verboseAlign > 0) {
    fprintf(stderr, "trimOverlap_Flipped()-- TRIM A %8u                     length %6u overlap %6d-%-6d\n", ovl->a_iid, _readData[_aID].trimmedLength, aovl5, aovl3);
    fprintf(stderr, "trimOverlap_Flipped()--      B %8u                     length %6u overlap %6d-%-6d\n", ovl->b_iid, _readData[_bID].trimmedLength, bovl3, bovl5);
  }

  //  If the resulting overlap is too small (or negative!), the overlap cannot be formed from the trimmed reads.

  if ((aovl3 - aovl5 < 1000) ||
      (bovl5 - bovl3 < 1000)) {
    if (_verboseAlign > 0)
      fprintf(stderr, "trimOverlap_Flipped()--      Too short after adjusting.\n");
    return(false);
  }

  //  Update the overlap with the new hangs.

  ovl->dat.ovl.ahg5 =                                 aovl5;
  ovl->dat.ovl.ahg3 = _readData[_aID].trimmedLength - aovl3;

  ovl->dat.ovl.bhg5 = _readData[_bID].trimmedLength - bovl5;
  ovl->dat.ovl.bhg3 =                                 bovl3;

  return(true);
}



bool
maComputation::trimOverlap(ovOverlap *ovl) {
  bool  success;

  _bID = ovl->b_iid;

  assert(_aID == ovl->a_iid);
  assert(_bID == ovl->b_iid);    //  Well, duh.

  if (ovl->flipped() == true)
    success = trimOverlap_Flipped(ovl);
  else
    success = trimOverlap_Normal(ovl);

  if ((_verboseAlign > 0) && (success == false))
    fprintf(stderr, "computeAlignments()-- Overlap trimmed out.\n");

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
maComputation::isWellContained(void) {
  uint32   wellContained          = 0;

  for (uint32 ii=0; ii<_overlapsLen; ii++) {
    ovOverlap *ov = _overlaps + ii;

    int32   a5        = (int32)ov->dat.ovl.ahg5;
    int32   a3        = (int32)ov->dat.ovl.ahg3;

    int32   b5        = (int32)ov->dat.ovl.bhg5;
    int32   b3        = (int32)ov->dat.ovl.bhg3;

    if ((b5 - a5 > _maxEdge) &&
        (b3 - a3 > _maxEdge)) {
      if (_verboseAlign > 0)
        fprintf(stderr, "isWellContained()-- read %8u well contained in read %8u\n", ov->a_iid, ov->b_iid);

      wellContained++;
    }
  }

  if ((_verboseAlign > 0) && (wellContained > _wellContainedThreshold))
    fprintf(stderr, "isWellContained()-- read %8u well contained\n", _aID);

  return(wellContained > _wellContainedThreshold);
}



//  Return the length of the shortest overlap we care about
//  of the 5' and 3' ends.
//
void
maComputation::findMinThickestEdge(int32   &thick5,
                                   int32   &thick3) {
  int32    thick5len = 0,  *thick5arr = new int32 [_overlapsLen];
  int32    thick3len = 0,  *thick3arr = new int32 [_overlapsLen];

  for (uint32 ii=0; ii<_overlapsLen; ii++) {
    ovOverlap *ov = _overlaps + ii;

    if (ov->evalue() == AS_MAX_EVALUE)    //  Skip overlaps that are flagged as trimmed out
      continue;                           //  or otherwise not useful.

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

#if 0
  for (int32 ii=0; ii<max(thick5len, thick3len); ii++) {
    if      ((ii < thick5len) && (ii < thick3len))
      fprintf(stderr, "%5u  %8d  %8d\n", ii, thick5arr[ii], thick3arr[ii]);
    else if (ii < thick5len)
      fprintf(stderr, "%5u  %8d  %8d\n", ii, thick5arr[ii], 0);
    else if (ii < thick3len)
      fprintf(stderr, "%5u  %8d  %8d\n", ii, 0,             thick3arr[ii]);
  }
#endif

  int32     novl   = (int32)ceil(_dovetailFraction * _coverage);

  thick5 = (novl < thick5len) ? thick5arr[thick5len - 1 - novl] : 0;
  thick3 = (novl < thick3len) ? thick3arr[thick3len - 1 - novl] : 0;

  delete [] thick5arr;
  delete [] thick3arr;

  if (_verboseAlign > 0) {
    fprintf(stderr, "findMinThickestEdge()-- read %8u thick5 %d (out of %d)\n", _aID, thick5, thick5len);
    fprintf(stderr, "findMinThickestEdge()--               thick3 %d (out of %d) novl %d\n", thick3, thick3len, novl);
  }
}


void
maComputation::computeAlignments(uint32  minOverlapLength,
                                 double  maxErate) {

  //  If there are no overlaps, there are no overlaps to align and output.

  if (_overlapsLen == 0)
    return;

  if (_verboseAlign > 0) {
    fprintf(stderr, "computeAlignments()-- \n");
    fprintf(stderr, "computeAlignments()-- ========================================BEGIN read %u\n", _overlaps[0].a_iid);
  }

  //  Adjust all the overlaps to trimmed reads.  We need to do this before we
  //  can start filtering weak overlaps.  Any overlap that is trimmed out
  //  needs to be flagged as bad.

  for (uint32 ii=0; ii<_overlapsLen; ii++)
    if (trimOverlap(_overlaps + ii) == false)
      _overlaps[ii].evalue(AS_MAX_EVALUE);

  //  If this read is 'well' contained in some other read, assume it's
  //  useless for assembly and throw out all of its overlaps.

  if (isWellContained() == true) {
    for (uint32 ii=0; ii<_overlapsLen; ii++)
      _overlaps[ii].evalue(AS_MAX_EVALUE);

    return;
  }

  //  Decide which overlaps to compute alignments for:  only the thickest on each side.

  int32   thick5;
  int32   thick3;

  findMinThickestEdge(thick5, thick3);

  //  Examine each overlap.  Either deicde it isn't useful for assembly, or compute, precisely,
  //  the overlap end points and save the alignment.

  fetchTrimmedRead(_aID);

  //
  //  Process each overlap.
  //

  for (uint32 ovlid=0; ovlid<_overlapsLen; ovlid++) {
    ovOverlap *ov = _overlaps + ovlid;

    if (ov->evalue() == AS_MAX_EVALUE)    //  Skip overlaps that are flagged as trimmed out
      continue;                           //  or otherwise not useful.

    if (_verboseAlign > 0) {
      fprintf(stderr, "\n");
      fprintf(stderr, "computeAlignments()-- ----------------------------------------OVERLAP %u/%u\n", ovlid, _overlapsLen);
    }

    assert(_aID == ov->a_iid);
    _bID = ov->b_iid;

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
    if ((aend - abgn < _minLength) ||
        (bend - bbgn < _minLength)) {
      if (_verboseAlign > 0)
        fprintf(stderr, "computeAlignments()-- overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d short.\n", ovlid, _aID, abgn, aend, _bID, bbgn, bend);
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
    if ((a5 - b5 > _maxEdge) &&
        (a3 - b3 > _maxEdge)) {
      if (_verboseAlign > 0)
        fprintf(stderr, "computeAlignments()-- overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d contained.\n", ovlid, _aID, abgn, aend, _bID, bbgn, bend);
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
      if (_verboseAlign > 0)
        fprintf(stderr, "computeAlignments()-- overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d thick3 %d < %d\n", ovlid, _aID, abgn, aend, _bID, bbgn, bend, _readData[_aID].trimmedLength - a5 + b5, thick3);
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
      if (_verboseAlign > 0)
        fprintf(stderr, "computeAlignments()-- overlap %4u a=%8u %6d-%-6d b=%8u %6d-%-6d thick5 %d < %d\n", ovlid, _aID, abgn, aend, _bID, bbgn, bend, _readData[_aID].trimmedLength - a3 + b3, thick5);
      ov->evalue(AS_MAX_EVALUE);
      continue;
    }

    //
    //  A good overlap to process!
    //


    //  Find, precisely, the optimal overlap for each read.  No alignments,
    //  just end points.  The ovOverlap is updated with correct end points.

    //  Load the b read and compute the alignment.
    fetchTrimmedRead(_bID, false, ov->flipped());

    computeOverlapAlignment(ovlid,
                            minOverlapLength,
                            maxErate,
                            _overlapSlop,
                            _maxRepeatLength);
  }

  //if (_verboseAlign > 0)
  //  fprintf(stderr, "computeAlignments()-- read %8u nShort %5u nContained %5u nUnaligned %5u nThin %5u filtered %.2f%%\n",
  //          _aID, nShort, nContained, nUnaligned, nThin, 100.0 * (nShort + nContained + nUnaligned + nThin) / _overlapsLen);



  //  Step 2:  Transfer to a tig.  generateCorrectionLayouts.C generateLayout().

  //  Step 3:  Align all reads in tig to tig sequence.  Save alignments
  //  in conveninent multialign structure.  This should be part of tgTig.

  //  Step 4:

}

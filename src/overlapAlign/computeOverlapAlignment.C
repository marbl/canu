
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

#include "edlib.H"

#include "overlapAlign-globalData.H"
#include "overlapAlign-threadData.H"
#include "overlapAlign-computation.H"



double   localMatch     =  1;
double   localMismatch  = -1;
double   localGap       = -1;

double
getScore(EdlibAlignResult result, int32 pp) {

  if ((pp < 0) ||
      (pp >= result.alignmentLength))
    return(localMatch);

  switch (result.alignment[pp]) {
    case EDLIB_EDOP_MATCH:     return(localMatch);     break;
    case EDLIB_EDOP_INSERT:    return(localGap);       break;
    case EDLIB_EDOP_DELETE:    return(localGap);       break;
    case EDLIB_EDOP_MISMATCH:  return(localMismatch);  break;
    default:                   return(0);              break;
  }
}


bool
testAlignment(char   *aRead,  int32  &abgn,  int32  &aend,  int32         alen,   uint32 Aid,
              char   *bRead,  int32  &bbgn,  int32  &bend,  int32         blen,   uint32 Bid,
              double  maxAlignErate,
              double  maxAcceptErate,
              double &erate) {

  assert(abgn >= 0);
  assert(aend <= alen);
  assert(bbgn >= 0);
  assert(bend <= blen);

  erate = 1.0;   //  Set the default return value, 100% error.

  //  Compute an alignment with the alignment path.

  EdlibAlignResult result = edlibAlign(aRead + abgn, aend - abgn,
                                       bRead + bbgn, bend - bbgn,
                                       edlibNewAlignConfig((int32)ceil(1.1 * maxAlignErate * ((aend - abgn) + (bend - bbgn)) / 2.0),
                                                           EDLIB_MODE_HW,
                                                           EDLIB_TASK_PATH));

  //  If no result, free it and return that the alignment failed.

  if (result.numLocations == 0) {
    edlibFreeAlignResult(result);
    return(false);
  }

  //  Compute the (approximate) length of the alignment, for EDLIB_TASK_LOC
  //result.alignmentLength = ((aend - abgn) + (result.endLocations[0] + 1 - result.startLocations[0]) + (result.editDistance)) / 2;

  //  Compute the global error rate.  If it's terrible, return that the alignment failed.

  erate = (double)result.editDistance / result.alignmentLength;

  if (erate > maxAcceptErate) {
    edlibFreeAlignResult(result);
    return(false);
  }

  //  Good quality.  Save the positions.
  //  The A positions don't change (yet).

  bend = bbgn + result.endLocations[0] + 1;    //  Edlib returns 0-based positions, add one to end to get space-based.
  bbgn = bbgn + result.startLocations[0];

  edlibFreeAlignResult(result);

  return(true);
}







int32
extend5(char   *aRead,  int32  &abgn,  int32  &aend,  int32         alen,   uint32 Aid,
        char   *bRead,  int32  &bbgn,  int32  &bend,  int32         blen,   uint32 Bid,
        double  maxAlignErate,
        double  maxAcceptErate,
        double &erate) {

  assert(abgn >= 0);
  assert(aend <= alen);
  assert(bbgn >= 0);
  assert(bend <= blen);

  erate = 1.0;   //  Set the default return value, 100% error.

  //  Compute an alignment with the alignment path.

  EdlibAlignResult result = edlibAlign(aRead + abgn, aend - abgn,   //  Query
                                       bRead + bbgn, bend - bbgn,   //  Target
                                       edlibNewAlignConfig((int32)ceil(1.1 * maxAlignErate * ((aend - abgn) + (bend - bbgn)) / 2.0),
                                                           EDLIB_MODE_HW,
                                                           EDLIB_TASK_PATH));

  //  If no result, free it and return that the alignment failed.

  if (result.numLocations == 0) {
    edlibFreeAlignResult(result);
    return(-1);
  }

#if 0
  char *aAln = new char [result.alignmentLength + 1];
  char *bAln = new char [result.alignmentLength + 1];

  edlibAlignmentToStrings(result,
                          aRead + abgn, aend - abgn,   //  Query
                          bRead + bbgn, bend - bbgn,   //  Target
                          aAln,
                          bAln);

  fprintf(stderr, "A %s\n", aAln);
  fprintf(stderr, "B %s\n", bAln);

  delete [] aAln;
  delete [] bAln;
#endif


  //  Compute the (approximate) length of the alignment, for EDLIB_TASK_LOC
  //result.alignmentLength = ((aend - abgn) + (result.endLocations[0] + 1 - result.startLocations[0]) + (result.editDistance)) / 2;

  //  Compute the global error rate.  If it's terrible, return that the alignment failed.

  erate = (double)result.editDistance / result.alignmentLength;

  //if (erate > maxAcceptErate) {
  //  edlibFreeAlignResult(result);
  //  return(-1);
  //}

  //  Good quality.  Save the positions.
  //  The A positions don't change (yet).

  bend = bbgn + result.endLocations[0] + 1;    //  Edlib returns 0-based positions, add one to end to get space-based.
  bbgn = bbgn + result.startLocations[0];

  //  Analyze the alignment in detail.  Trim back bad ends.
  //
  //  We're trying to extend the 5' end of the A sequence.
  //  Scan the alignment backwards to check if we enter a region of bad alignment.

  int32    localWindow    = 50;
  double   localScore     = localMatch * localWindow;

  double   localThreshold = 0.50;

  int32  abgnN = aend;
  int32  bbgnN = bend;
  int32  ii    = result.alignmentLength - 1;   //  Window start (left most position of window)

  //  Scan towards the 5' end until we hit a region of local garbage.

  for (; ii >= 0; ii--) {
    localScore += getScore(result, ii);
    localScore -= getScore(result, ii + localWindow);

    abgnN -= (result.alignment[ii] != EDLIB_EDOP_DELETE);   //  Move left if not gap in B (delete in B)?
    bbgnN -= (result.alignment[ii] != EDLIB_EDOP_INSERT);   //  Move left if not gap in A (insert in B)?

    assert(abgnN >= 0);
    assert(bbgnN >= 0);
    //fprintf(stderr, "NEW %d %d after align %d\n", abgnN, bbgnN, result.alignment[ii]);

    //fprintf(stderr, "extend5      at ii %d = %f   A %d-%d  B %d-%d\n",
    //        ii, localScore / localWindow, abgnN, aend, bbgnN, bend);

    if (localScore / localWindow < localThreshold) {
      //fprintf(stderr, "extend5 fail at ii %d = %f   A %d-%d  B %d-%d\n",
      //        ii, localScore / localWindow, abgnN, aend, bbgnN, bend);
      break;
    }
  }

  //  If the score is good, the whole extension is good.

  if (localScore / localWindow >= localThreshold) {
    edlibFreeAlignResult(result);
    return(0);
  }

  //  Otherwise, scan the window to find the maximal score and trim to there.

  double  maxWindowScore = localScore;
  int32   ww             = ii;

  for (int32 kk=0; kk<localWindow; kk++) {
    localScore -= getScore(result, ii + kk);
    localScore += getScore(result, ii + kk + localWindow);

    if (maxWindowScore < localScore) {
      maxWindowScore = localScore;
      ww             = ii + kk;
    }
  }

  //  Now that we know where we should terminate the alignment, back up to that spot,
  //  adjusting the coordinates.

  for (int32 kk=ii; kk<=ww; kk++) {
    abgnN += (result.alignment[kk] != EDLIB_EDOP_DELETE);   //  Delete in B?
    bbgnN += (result.alignment[kk] != EDLIB_EDOP_INSERT);   //  Insert in B?
  }

  //fprintf(stderr, "extend5 fail at ii %d = %f   abgn %d  bbgn %d\n",
  //        ii, localScore / localWindow, abgnN, bbgnN);

  abgn = abgnN;
  bbgn = bbgnN;

  edlibFreeAlignResult(result);
  return(1);
}











int32
extend3(char   *aRead,  int32  &abgn,  int32  &aend,  int32         alen,   uint32 Aid,
        char   *bRead,  int32  &bbgn,  int32  &bend,  int32         blen,   uint32 Bid,
        double  maxAlignErate,
        double  maxAcceptErate,
        double &erate) {

  assert(abgn >= 0);
  assert(aend <= alen);
  assert(bbgn >= 0);
  assert(bend <= blen);

  erate = 1.0;   //  Set the default return value, 100% error.

  //  Compute an alignment with the alignment path.

  EdlibAlignResult result = edlibAlign(aRead + abgn, aend - abgn,
                                       bRead + bbgn, bend - bbgn,
                                       edlibNewAlignConfig((int32)ceil(1.1 * maxAlignErate * ((aend - abgn) + (bend - bbgn)) / 2.0),
                                                           EDLIB_MODE_HW,
                                                           EDLIB_TASK_PATH));

  //  If no result, free it and return that the alignment failed.

  if (result.numLocations == 0) {
    edlibFreeAlignResult(result);
    return(-1);
  }

#if 0
  char *aAln = new char [result.alignmentLength + 1];
  char *bAln = new char [result.alignmentLength + 1];

  edlibAlignmentToStrings(result,
                          aRead + abgn, aend - abgn,   //  Query
                          bRead + bbgn, bend - bbgn,   //  Target
                          aAln,
                          bAln);

  fprintf(stderr, "A %s\n", aAln);
  fprintf(stderr, "B %s\n", bAln);

  delete [] aAln;
  delete [] bAln;
#endif

  //  Compute the (approximate) length of the alignment, for EDLIB_TASK_LOC
  //result.alignmentLength = ((aend - abgn) + (result.endLocations[0] + 1 - result.startLocations[0]) + (result.editDistance)) / 2;

  //  Compute the global error rate.  If it's terrible, return that the alignment failed.

  erate = (double)result.editDistance / result.alignmentLength;

  //if (erate > maxAcceptErate) {
  //  edlibFreeAlignResult(result);
  //  return(-1);
  //}

  //  Good quality.  Save the positions.
  //  The A positions don't change (yet).

  bend = bbgn + result.endLocations[0] + 1;    //  Edlib returns 0-based positions, add one to end to get space-based.
  bbgn = bbgn + result.startLocations[0];

  //  Analyze the alignment in detail.  Trim back bad ends.
  //
  //  We're trying to extend the 5' end of the A sequence.
  //  Scan the alignment backwards to check if we enter a region of bad alignment.

  int32    localWindow    = 50;
  double   localScore     = localMatch * localWindow;

  double   localThreshold = 0.50;

  int32  aendN = abgn;
  int32  bendN = bbgn;
  int32  ii    = 0;        //  Window end (right most position of window)

  //  Scan towards the 3' end until we hit a region of local garbage.

  for (; ii < result.alignmentLength; ii++) {
    localScore += getScore(result, ii);
    localScore -= getScore(result, ii - localWindow);

    aendN += (result.alignment[ii] != EDLIB_EDOP_DELETE);   //  Delete in B?
    bendN += (result.alignment[ii] != EDLIB_EDOP_INSERT);   //  Insert in B?

    assert(aendN <= alen);
    assert(bendN <= blen);

    if (localScore / localWindow < localThreshold) {
      //fprintf(stderr, "extend3 fail at ii %d = %f\n", ii, localScore / localWindow);
      break;
    }
  }

  //  If the score is good, the whole extension is good.

  if (localScore / localWindow >= localThreshold) {
    edlibFreeAlignResult(result);
    //fprintf(stderr, "extend3 all good!\n");
    return(0);
  }

  //  Otherwise, scan the window to find the maximal score and trim to there.

  double  maxWindowScore = localScore;
  int32   ww             = ii;

  for (int32 kk=0; kk<localWindow; kk++) {
    localScore -= getScore(result, ii - kk);
    localScore += getScore(result, ii - kk - localWindow);

    if (maxWindowScore < localScore) {
      maxWindowScore = localScore;
      ww             = ii - kk;
    }
  }

  //  Now that we know where we should terminate the alignment, back up to that spot,
  //  adjusting the coordinates.

  for (int32 kk=ii; kk > ww; kk--) {
    aendN -= (result.alignment[kk] != EDLIB_EDOP_DELETE);   //  Delete in B?
    bendN -= (result.alignment[kk] != EDLIB_EDOP_INSERT);   //  Insert in B?
  }

  aend = aendN;
  bend = bendN;

  edlibFreeAlignResult(result);
  return(1);
}















//  Compute an alignment between A and B,
//    requring all of A to be in the alignment,
//    returning new end points on B, allowing free gaps on the ends.
//  The editDist and alignLen of the alignment are returned.
//
bool
computeAlignment(char  *aRead,  int32   abgn,  int32   aend,  int32  UNUSED(alen),  char const *Alabel,  uint32 Aid,
                 char  *bRead,  int32  &bbgn,  int32  &bend,  int32         blen,   char const *Blabel,  uint32 Bid,
                 double  maxErate,
                 int32  &editDist,
                 int32  &alignLen,
                 uint32  verbose) {
  bool  success = false;

  editDist = 0;
  alignLen = 0;

  if (verbose > 0) {
    fprintf(stderr, "computeAlignment()--         align %s %6u %6d-%-6d  %.2f%% error\n", Alabel, Aid, abgn, aend, 100.0 * maxErate);
    fprintf(stderr, "computeAlignment()--            vs %s %6u %6d-%-6d\n",               Blabel, Bid, bbgn, bend);
  }

  EdlibAlignResult result = edlibAlign(aRead + abgn, aend - abgn,
                                       bRead + bbgn, bend - bbgn,
                                       edlibNewAlignConfig((int32)ceil(1.1 * maxErate * ((aend - abgn) + (bend - bbgn)) / 2.0),
                                                           EDLIB_MODE_HW,
                                                           EDLIB_TASK_LOC));

  //  If there is a result, compute the (approximate) length of the alignment.
  //  Edlib mode TASK_LOC doesn't populate this field.

  if (result.numLocations > 0)
    result.alignmentLength = ((aend - abgn) + (result.endLocations[0] + 1 - result.startLocations[0]) + (result.editDistance)) / 2;

  //  Save the result if it is of acceptable quality.

  //fprintf(stderr, "XX numLocations %d editDist %d maxErate %f alignLength %d\n",
  //        result.numLocations, result.editDistance, maxErate, result.alignmentLength);

  if ((result.numLocations > 0) &&
      (result.editDistance <= maxErate * result.alignmentLength)) {
    success  = true;

    bend     = bbgn + result.endLocations[0] + 1;    //  Edlib returns 0-based positions, add one to end to get space-based.
    bbgn     = bbgn + result.startLocations[0];

    editDist = result.editDistance;
    alignLen = result.alignmentLength;

    if (verbose > 0) {
      fprintf(stderr, "computeAlignment()--         PASS  %s %6u %6d-%-6d editDist %5d alignLen %6d error %6.4f\n",
              Blabel, Bid, bbgn, bend,
              result.editDistance,
              result.alignmentLength,
              result.editDistance / (double)result.alignmentLength);
    }
  }

  else {
    if ((verbose > 0) && (result.numLocations > 0)) {
      fprintf(stderr, "computeAlignment()--         FAIL  %s %6u %6d-%-6d editDist %5d alignLen %6d error %6.4f\n",
              Blabel, Bid, bbgn + result.endLocations[0] + 1, bbgn + result.startLocations[0],
              result.editDistance,
              result.alignmentLength,
              result.editDistance / (double)result.alignmentLength);
    }

    if ((verbose > 0) && (result.numLocations == 0)) {
      fprintf(stderr, "computeAlignment()--         FAIL  %s %6u %6d-%-6d no alignment\n",
              Blabel, Bid, bbgn, bend);
    }
  }

  edlibFreeAlignResult(result);

  return(success);
}




//  Our goal is to convert an alignment-free overlap with imprecise edges into
//  an alignment-based overlap with precise edges.
//
//  The strategy is to, for each end of the overlap, extend the alignment to the end
//  of the sequence with the shorter hang, then compute the alignment of
//  that to the other read allowing free gaps.
//
//                { alignment-free overlap }
//      ---(a5----{-----a3)----------------}------
//           (b5--{---b3)------------------}-------------------
//
//  Since the 5' hang on the B read is smaller, we use b5 to compute the other values.
//    b5 = known from alignment-free overlap
//    b3 = max(4*b5, 2*maxRepeatLength)
//    a5 = b5 * (1 + maxErate) + SLOP
//    a3 = b3 * (1 + maxErate) + SLOP
//
//  For b3: this is our anchor in assumed good sequence.  it needs to be larger than
//          a repeat, and should be larger than the unknown sequence
//
//  For a5: we need to allow space for all of b5 to align, plus any errors.  The
//          SLOP adjustment allows for shifting of the imprecise overlap edge.
//
//  Similar for a3.
//
//  Output of the alignment will be updated coordinates for a5 and a3,
//  which we use to update this end of the overlap.
//
void
maComputation::computeOverlapAlignment(uint32       ovlid,
                                       uint32       minOverlapLength,
                                       double       maxErate,
                                       uint32       overlapSlop,
                                       uint32       maxRepeat) {
  ovOverlap *ovl    = &_overlaps[ovlid];   //  Convenience pointer to the overlap.
  ovOverlap  ori    =  _overlaps[ovlid];   //  Copy of the original overlap.

  int32   alen      = _readData[_aID].trimmedLength;
  int32   abgn      = (int32)       ovl->dat.ovl.ahg5;
  int32   aend      = (int32)alen - ovl->dat.ovl.ahg3;

  int32   blen      = _readData[_bID].trimmedLength;
  int32   bbgn      = (int32)       ovl->dat.ovl.bhg5;
  int32   bend      = (int32)blen - ovl->dat.ovl.bhg3;

  int32   alignLen  = 0;
  int32   editDist  = INT32_MAX;

  EdlibAlignResult  result = { 0, NULL, NULL, 0, NULL, 0, 0 };

  if (_verboseAlign > 0) {
    fprintf(stderr, "computeOverlapAlignment()-- A %8u %6d-%-6d %d\n", _aID, abgn, aend, alen);
    fprintf(stderr, "computeOverlapAlignment()-- B %8u %6d-%-6d %d\n", _bID, bbgn, bend, blen);
  }

  //  A    ------------{------...
  //  B          ------{------...
  if (ovl->dat.ovl.bhg5 < ovl->dat.ovl.ahg5) {
    int32   ahg5 = ovl->dat.ovl.ahg5;
    int32   ahg3 = ovl->dat.ovl.ahg3;
    int32   bhg5 = ovl->dat.ovl.bhg5;
    int32   bhg3 = ovl->dat.ovl.bhg3;

    int32   b5   = bhg5;
    int32   b3   = max(4 * bhg5, (int32)(2 * maxRepeat));

    int32   a5   = b5 * (1 + maxErate) + overlapSlop;
    int32   a3   = b3 * (1 + maxErate) + overlapSlop;

    int32   bbgn = max(0,    bhg5 - b5);    //  Now zero.
    int32   bend = min(blen, bhg5 + b3);

    int32   abgn = max(0,    ahg5 - a5);
    int32   aend = min(alen, ahg5 + a3);

    if (_verboseAlign > 0)
      fprintf(stderr, "computeOverlapAlignment()-- bhg5:  B %d-%d onto A %d-%d\n", bbgn, bend, abgn, aend);

    if (computeAlignment(_bRead, bbgn, bend, blen, "B", _bID,    //  Align all of sequence B into
                         _aRead, abgn, aend, alen, "A", _aID,    //  sequence A with free ends.
                         maxErate,
                         editDist,
                         alignLen, _verboseAlign) == true) {
      if (_verboseAlign > 0)
        fprintf(stderr, "computeOverlapAlignment()--        B %d-%d onto A %d-%d\n", bbgn, bend, abgn, aend);

      ovl->dat.ovl.ahg5 = abgn;
      ovl->dat.ovl.bhg5 = bbgn;

      if ((bend == blen) ||                      //  Update the other end if the B read is contained, or
          (alen - ovl->dat.ovl.ahg3 < aend)) {   //  if the alignment extends past our existing end.
        ovl->dat.ovl.ahg3 = alen - aend;
        ovl->dat.ovl.bhg3 = blen - bend;
      }
    }

    else {
    }
  }

  //  A          ------{------...
  //  B    ------------{------...
  else {
    int32   ahg5 = ovl->dat.ovl.ahg5;
    int32   ahg3 = ovl->dat.ovl.ahg3;
    int32   bhg5 = ovl->dat.ovl.bhg5;
    int32   bhg3 = ovl->dat.ovl.bhg3;

    int32   a5   = ahg5;
    int32   a3   = max(4 * ahg5, (int32)(2 * maxRepeat));

    int32   b5   = a5 * (1 + maxErate) + overlapSlop;
    int32   b3   = a3 * (1 + maxErate) + overlapSlop;

    int32   bbgn = max(0,    bhg5 - b5);
    int32   bend = min(blen, bhg5 + b3);

    int32   abgn = max(0,    ahg5 - a5);    //  Now zero.
    int32   aend = min(alen, ahg5 + a3);

    if (_verboseAlign > 0)
      fprintf(stderr, "computeOverlapAlignment()-- ahg5:  A %d-%d onto B %d-%d\n", abgn, aend, bbgn, bend);

    if (computeAlignment(_aRead, abgn, aend, alen, "A", _aID,    //  Align all of sequence A into
                         _bRead, bbgn, bend, blen, "B", _bID,    //  sequence B with free ends.
                         maxErate,
                         editDist,
                         alignLen, _verboseAlign) == true) {
      if (_verboseAlign > 0)
        fprintf(stderr, "computeOverlapAlignment()--        A %d-%d onto B %d-%d\n", abgn, aend, bbgn, bend);

      ovl->dat.ovl.ahg5 = abgn;
      ovl->dat.ovl.bhg5 = bbgn;

      if ((aend == alen) ||                      //  Update the other end if the A read is contained, or
          (blen - ovl->dat.ovl.bhg3 < bend)) {   //  if the alignment extends past our existing end.
        ovl->dat.ovl.ahg3 = alen - aend;
        ovl->dat.ovl.bhg3 = blen - bend;
      }
    }

    else {
    }
  }




  //  A    ...------}------------
  //  B    ...------}------
  if (ovl->dat.ovl.bhg3 < ovl->dat.ovl.ahg3) {
    int32   ahg5 = ovl->dat.ovl.ahg5;
    int32   ahg3 = ovl->dat.ovl.ahg3;
    int32   bhg5 = ovl->dat.ovl.bhg5;
    int32   bhg3 = ovl->dat.ovl.bhg3;

    int32   b5   = max(4 * bhg3, (int32)(2 * maxRepeat));
    int32   b3   = bhg3;

    int32   a5   = b5 * (1 + maxErate) + overlapSlop;
    int32   a3   = b3 * (1 + maxErate) + overlapSlop;

    int32   bbgn = max(0,    blen - bhg3 - b5);
    int32   bend = min(blen, blen - bhg3 + b3);    //  Now blen.

    int32   abgn = max(0,    alen - ahg3 - a5);
    int32   aend = min(alen, alen - ahg3 + a3);

    if (_verboseAlign > 0)
      fprintf(stderr, "computeOverlapAlignment()-- bhg3:  B %d-%d onto A %d-%d\n", bbgn, bend, abgn, aend);

    if (computeAlignment(_bRead, bbgn, bend, blen, "B", _bID,    //  Align all of sequence B into
                         _aRead, abgn, aend, alen, "A", _aID,    //  sequence A with free ends.
                         maxErate,
                         editDist,
                         alignLen, _verboseAlign) == true) {
      if (_verboseAlign > 0)
        fprintf(stderr, "computeOverlapAlignment()--        B %d-%d onto A %d-%d\n", bbgn, bend, abgn, aend);

      if ((bbgn == 0) ||                  //  Update the other end if the B read is contained, or
          (abgn < ovl->dat.ovl.ahg5)) {   //  if the alignment extends past our existing end.
        ovl->dat.ovl.ahg5 = abgn;
        ovl->dat.ovl.bhg5 = bbgn;
      }

      ovl->dat.ovl.ahg3 = alen - aend;
      ovl->dat.ovl.bhg3 = blen - bend;
    }

    else {
    }
  }

  //  A    ...------}------
  //  B    ...------}------------
  else {
    int32   ahg5 = ovl->dat.ovl.ahg5;
    int32   ahg3 = ovl->dat.ovl.ahg3;
    int32   bhg5 = ovl->dat.ovl.bhg5;
    int32   bhg3 = ovl->dat.ovl.bhg3;

    int32   a5   = max(4 * ahg3, (int32)(2 * maxRepeat));
    int32   a3   = ahg3;

    int32   b5   = a5 * (1 + maxErate) + overlapSlop;
    int32   b3   = a3 * (1 + maxErate) + overlapSlop;

    int32   bbgn = max(0,    blen - bhg3 - b5);
    int32   bend = min(blen, blen - bhg3 + b3);

    int32   abgn = max(0,    alen - ahg3 - a5);    //  Now alen.
    int32   aend = min(alen, alen - ahg3 + a3);

    if (_verboseAlign > 0)
      fprintf(stderr, "computeOverlapAlignment()-- ahg3:  A %d-%d onto B %d-%d\n", abgn, aend, bbgn, bend);

    if (computeAlignment(_aRead, abgn, aend, alen, "A", _aID,    //  Align all of sequence A into
                         _bRead, bbgn, bend, blen, "B", _bID,    //  sequence B with free ends.
                         maxErate,
                         editDist,
                         alignLen, _verboseAlign) == true) {
      if (_verboseAlign > 0)
        fprintf(stderr, "computeOverlapAlignment()--        A %d-%d onto B %d-%d\n", abgn, aend, bbgn, bend);

      if ((abgn == 0) ||                  //  Update the other end if the A read is contained, or
          (bbgn < ovl->dat.ovl.bhg5)) {   //  if the alignment extends past our existing end.
        ovl->dat.ovl.ahg5 = abgn;
        ovl->dat.ovl.bhg5 = bbgn;
      }

      ovl->dat.ovl.ahg3 = alen - aend;
      ovl->dat.ovl.bhg3 = blen - bend;
    }

    else {
    }
  }

  //
  //  Compute the final alignment, then copy it to our output buffer if it is
  //  of good quality.
  //

  if (1) {
    int32   abgn      = (int32)       ovl->dat.ovl.ahg5;
    int32   aend      = (int32)alen - ovl->dat.ovl.ahg3;

    int32   bbgn      = (int32)       ovl->dat.ovl.bhg5;
    int32   bend      = (int32)blen - ovl->dat.ovl.bhg3;

    if (_verboseAlign > 0)
      fprintf(stderr, "computeOverlapAlignment()-- final:   A %d-%d vs B %d-%d\n", abgn, aend, bbgn, bend);

    EdlibAlignResult result = edlibAlign(_aRead + abgn, aend - abgn,
                                         _bRead + bbgn, bend - bbgn,
                                         edlibNewAlignConfig((int32)ceil(1.1 * maxErate * ((aend - abgn) + (bend - bbgn)) / 2.0),
                                                             EDLIB_MODE_NW,
                                                             EDLIB_TASK_PATH));

    //  Decide, based on the edit distance and alignment length, if we should
    //  retain or discard the overlap.

    if ((result.alignmentLength < minOverlapLength) ||                  //  Alignment is too short, or
        (result.editDistance > maxErate * result.alignmentLength)) {    //               too noisy.
      ovl->dat.ovl.forOBT = false;
      ovl->dat.ovl.forDUP = false;
      ovl->dat.ovl.forUTG = false;

      ovl->evalue(AS_MAX_EVALUE);
    }

    else {
      ovl->dat.ovl.forOBT = false;                                //  Flag as useful for contigging
      ovl->dat.ovl.forDUP = false;                                //  if it is a dovetail overlap.
      ovl->dat.ovl.forUTG = (ovl->overlapIsDovetail() == true);

      ovl->erate(editDist / (double)alignLen);

      _alignsA  [ovlid] = new char [alen + 1];                    //  Allocate space for the alignment output.
      _alignsB  [ovlid] = new char [alen + 1];

      for (uint32 ii=0; ii<alen; ii++) {
        _alignsA[ovlid][ii] = '-';
        _alignsB[ovlid][ii] = '-';
      }

      _alignsA[ovlid][alen] = 0;
      _alignsB[ovlid][alen] = 0;

      char *aaln = new char [result.alignmentLength + 1];         //  Convert the alignment to a string.
      char *baln = new char [result.alignmentLength + 1];

      memset(aaln, 0, sizeof(char) * (result.alignmentLength + 1));
      memset(baln, 0, sizeof(char) * (result.alignmentLength + 1));

      edlibAlignmentToStrings(result.alignment,                   //  Alignment
                              result.alignmentLength,             //    and length
                              result.startLocations[0],           //  tgtStart (_bRead)
                              result.endLocations[0]+1,           //  tgtEnd
                              0,                                  //  qryStart (_aRead)
                              alen,                               //  qryEnd
                              _bRead + bbgn,                      //  tgt sequence (_bRead)
                              _aRead + abgn,                      //  qry sequence (_aRead)
                              baln,                               //  output tgt alignment string
                              aaln);                              //  output qry alignment string

      //  Dump

      //fprintf(stdout, "B %u    abgn %d aend %d  ahg5 %d ahg5 %d alen %d\n", _bID, abgn, aend, (int32)ovl->dat.ovl.ahg5, (int32)ovl->dat.ovl.ahg3, alen);
      //fprintf(stdout, "aaln %s\n", aaln);
      //fprintf(stdout, "baln %s\n", baln);

      //  Copy it into the output space.

      uint32 pp = 0;   //  Position in the emitted alignment string.
      uint32 aa = 0;   //  Position in sequence A.

      for (uint32 ii=0; ii<ovl->dat.ovl.ahg5; ii++, pp++) {
        _alignsA[ovlid][pp] = _aRead[aa];
        _alignsB[ovlid][pp] = '-';
        aa++;
      }

      for (uint32 ii=0; ii<result.alignmentLength; ii++) {
        if (aaln[ii] != '-') {
          _alignsA[ovlid][pp] = aaln[ii];
          _alignsB[ovlid][pp] = baln[ii];
          aa++;
          pp++;
        }
      }

      for (uint32 ii=0; ii<ovl->dat.ovl.ahg3; ii++, pp++) {
        _alignsA[ovlid][pp] = _aRead[aa];
        _alignsB[ovlid][pp] = '-';
        aa++;
      }

      assert(pp == alen);     //  If correct, we should have walked over the
      assert(aa == alen);     //  full trimmed read with both indices.

      delete [] aaln;
      delete [] baln;

      _alignsA[ovlid][alen] = 0;
      _alignsB[ovlid][alen] = 0;
    }

    //  Trash the alignment results.

    edlibFreeAlignResult(result);
  }

  //  More logging.

  if (_verboseAlign > 0) {
    fprintf(stderr, "computeOverlapAlignment()-- A %7u %6d-%-6d -> %6d-%-6d\n",       _aID, ori.a_bgn(), ori.a_end(), ovl->a_bgn(), ovl->a_end());
    fprintf(stderr, "computeOverlapAlignment()-- B %7u %6d-%-6d -> %6d-%-6d %5.4f\n", _bID, ori.b_bgn(), ori.b_end(), ovl->b_bgn(), ovl->b_end(), ovl->erate());
  }
}



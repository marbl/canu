
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

//#include <pthread.h>

#include "sqStore.H"
#include "sqCache.H"
#include "ovStore.H"

#include "edlib.H"
#include "alignStats.H"

#define SLOP       500
#define MAX_REPEAT 2000




bool
testAlignment(char   *aRead,  int32   abgn,  int32   aend,  int32         alen,   uint32 Aid,
              char   *bRead,  int32  &bbgn,  int32  &bend,  int32         blen,   uint32 Bid,
              double  maxAlignErate,
              double  maxAcceptErate,
              double &erate) {

  assert(abgn >= 0);
  assert(aend <= alen);
  assert(bbgn >= 0);
  assert(bend <= blen);

  erate = 1.0;   //  Set the default return value, 100% error.

  EdlibAlignResult result = edlibAlign(aRead + abgn, aend - abgn,
                                       bRead + bbgn, bend - bbgn,
                                       edlibNewAlignConfig((int32)ceil(1.1 * maxAlignErate * ((aend - abgn) + (bend - bbgn)) / 2.0),
                                                           EDLIB_MODE_HW,
                                                           EDLIB_TASK_LOC));

  //  If there is a result, compute the (approximate) length of the alignment and the error rate.
  //  Edlib mode TASK_LOC doesn't populate this field.

  if (result.numLocations > 0) {
    result.alignmentLength = ((aend - abgn) + (result.endLocations[0] + 1 - result.startLocations[0]) + (result.editDistance)) / 2;

    erate = (double)result.editDistance / result.alignmentLength;

    //  Save the result if it is of acceptable quality.
    if (erate < maxAcceptErate) {
      bend     = bbgn + result.endLocations[0] + 1;    //  Edlib returns 0-based positions, add one to end to get space-based.
      bbgn     = bbgn + result.startLocations[0];
    }
  }

  edlibFreeAlignResult(result);

  return(erate < maxAcceptErate);
}




//  Compute an alignment between A and B,
//    requring all of A to be in the alignment,
//    returning new end points on B, allowing free gaps on the ends.
//  The editDist and alignLen of the alignment are returned.
//
bool
computeAlignment(char  *aRead,  int32   abgn,  int32   aend,  int32  UNUSED(alen),  char *Alabel,  uint32 Aid,
                 char  *bRead,  int32  &bbgn,  int32  &bend,  int32         blen,   char *Blabel,  uint32 Bid,
                 double  maxErate,
                 int32  &editDist,
                 int32  &alignLen) {
  bool  success = false;
  bool  debug   = false;

  if (debug)
    fprintf(stderr, "  align %s %6u %6d-%-6d to %s %6u %6d-%-6d",
            Alabel, Aid, abgn, aend,
            Blabel, Bid, bbgn, bend);

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

  if ((result.numLocations > 0) &&
      (result.editDistance <= maxErate * result.alignmentLength)) {
    success  = true;

    bend     = bbgn + result.endLocations[0] + 1;    //  Edlib returns 0-based positions, add one to end to get space-based.
    bbgn     = bbgn + result.startLocations[0];

    editDist = result.editDistance;
    alignLen = result.alignmentLength;

    if (debug)
      fprintf(stderr, "    aligned to %s at %6d-%-6d editDist %5d alignLen %6d qual %6.4f\n",
              Blabel, bbgn, bend, editDist, alignLen, 1.0 - editDist / (double)alignLen);
  }

  else {
    if (debug)
      fprintf(stderr, "    FAILED  to %s at %6d-%-6d editDist %5d alignLen %6d qual %6.4f\n",
              Blabel, bbgn, bend, editDist, alignLen, 1.0 - editDist / (double)alignLen);
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
//      ---(a5----{-----a3)-----------------}------
//           (b5--{---b3)-------------------}-------------------
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
computeOverlapAlignment(ovOverlap   *ovl,
                        char        *aseq, int32 alen,
                        char        *bseq, int32 blen,
                        uint32       minOverlapLength,
                        double       maxErate,
                        bool         partialOverlaps,
                        alignStats  &localStats) {
  bool    debug     = false;

  ovOverlap ori     = *ovl;

  uint32  aID       = ovl->a_iid;
  int32   abgn      = (int32)       ovl->dat.ovl.ahg5;
  int32   aend      = (int32)alen - ovl->dat.ovl.ahg3;

  uint32  bID       = ovl->b_iid;
  int32   bbgn      = (int32)       ovl->dat.ovl.bhg5;
  int32   bend      = (int32)blen - ovl->dat.ovl.bhg3;

  int32   alignLen  = 0;
  int32   editDist  = INT32_MAX;

  EdlibAlignResult  result = { 0, NULL, NULL, 0, NULL, 0, 0 };

  if (debug) {
    fprintf(stderr, "\n");
    fprintf(stderr, "A %8u %6d-%-6d %d\n", aID, abgn, aend, alen);
    fprintf(stderr, "B %8u %6d-%-6d %d\n", bID, bbgn, bend, blen);
  }

  //  A    ------------{------...
  //  B          ------{------...
  if (ovl->dat.ovl.bhg5 < ovl->dat.ovl.ahg5) {
    int32   ahg5 = ovl->dat.ovl.ahg5;
    int32   ahg3 = ovl->dat.ovl.ahg3;
    int32   bhg5 = ovl->dat.ovl.bhg5;
    int32   bhg3 = ovl->dat.ovl.bhg3;

    int32   b5   = bhg5;
    int32   b3   = max(4 * bhg5, 2 * MAX_REPEAT);

    int32   a5   = b5 * (1 + maxErate) + SLOP;
    int32   a3   = b3 * (1 + maxErate) + SLOP;

    int32   bbgn = max(0,    bhg5 - b5);    //  Now zero.
    int32   bend = min(blen, bhg5 + b3);

    int32   abgn = max(0,    ahg5 - a5);
    int32   aend = min(alen, ahg5 + a3);

    if (debug)
      fprintf(stderr, "bhg5:  B %d-%d onto A %d-%d\n", bbgn, bend, abgn, aend);

    if (computeAlignment(bseq, bbgn, bend, blen, "B", bID,    //  Align all of sequence B into
                         aseq, abgn, aend, alen, "A", aID,    //  sequence A with free ends.
                         maxErate,
                         editDist,
                         alignLen) == true) {
      if (debug)
        fprintf(stderr, "       B %d-%d onto A %d-%d\n", bbgn, bend, abgn, aend);

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
    int32   a3   = max(4 * ahg5, 2 * MAX_REPEAT);

    int32   b5   = a5 * (1 + maxErate) + SLOP;
    int32   b3   = a3 * (1 + maxErate) + SLOP;

    int32   bbgn = max(0,    bhg5 - b5);
    int32   bend = min(blen, bhg5 + b3);

    int32   abgn = max(0,    ahg5 - a5);    //  Now zero.
    int32   aend = min(alen, ahg5 + a3);

    if (debug)
      fprintf(stderr, "ahg5:  A %d-%d onto B %d-%d\n", abgn, aend, bbgn, bend);

    if (computeAlignment(aseq, abgn, aend, alen, "A", aID,    //  Align all of sequence A into
                         bseq, bbgn, bend, blen, "B", bID,    //  sequence B with free ends.
                         maxErate,
                         editDist,
                         alignLen) == true) {
      if (debug)
        fprintf(stderr, "       A %d-%d onto B %d-%d\n", abgn, aend, bbgn, bend);

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

    int32   b5   = max(4 * bhg3, 2 * MAX_REPEAT);
    int32   b3   = bhg3;

    int32   a5   = b5 * (1 + maxErate) + SLOP;
    int32   a3   = b3 * (1 + maxErate) + SLOP;

    int32   bbgn = max(0,    blen - bhg3 - b5);
    int32   bend = min(blen, blen - bhg3 + b3);    //  Now blen.

    int32   abgn = max(0,    alen - ahg3 - a5);
    int32   aend = min(alen, alen - ahg3 + a3);

    if (debug)
      fprintf(stderr, "bhg3:  B %d-%d onto A %d-%d\n", bbgn, bend, abgn, aend);

    if (computeAlignment(bseq, bbgn, bend, blen, "B", bID,    //  Align all of sequence B into
                         aseq, abgn, aend, alen, "A", aID,    //  sequence A with free ends.
                         maxErate,
                         editDist,
                         alignLen) == true) {
      if (debug)
        fprintf(stderr, "       B %d-%d onto A %d-%d\n", bbgn, bend, abgn, aend);

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

    int32   a5   = max(4 * ahg3, 2 * MAX_REPEAT);
    int32   a3   = ahg3;

    int32   b5   = a5 * (1 + maxErate) + SLOP;
    int32   b3   = a3 * (1 + maxErate) + SLOP;

    int32   bbgn = max(0,    blen - bhg3 - b5);
    int32   bend = min(blen, blen - bhg3 + b3);

    int32   abgn = max(0,    alen - ahg3 - a5);    //  Now alen.
    int32   aend = min(alen, alen - ahg3 + a3);

    if (debug)
      fprintf(stderr, "ahg3:  A %d-%d onto B %d-%d\n", abgn, aend, bbgn, bend);

    if (computeAlignment(aseq, abgn, aend, alen, "A", aID,    //  Align all of sequence A into
                         bseq, bbgn, bend, blen, "B", bID,    //  sequence B with free ends.
                         maxErate,
                         editDist,
                         alignLen) == true) {
      if (debug)
        fprintf(stderr, "       A %d-%d onto B %d-%d\n", abgn, aend, bbgn, bend);

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




  if (1) {
    int32   abgn      = (int32)       ovl->dat.ovl.ahg5;
    int32   aend      = (int32)alen - ovl->dat.ovl.ahg3;

    int32   bbgn      = (int32)       ovl->dat.ovl.bhg5;
    int32   bend      = (int32)blen - ovl->dat.ovl.bhg3;

    if (debug)
      fprintf(stderr, "final:   A %d-%d vs B %d-%d\n", abgn, aend, bbgn, bend);

    EdlibAlignResult result = edlibAlign(aseq + abgn, aend - abgn,
                                         bseq + bbgn, bend - bbgn,
                                         edlibNewAlignConfig((int32)ceil(1.1 * maxErate * ((aend - abgn) + (bend - bbgn)) / 2.0),
                                                             EDLIB_MODE_NW,
                                                             EDLIB_TASK_PATH));

    //  Decide, based on the edit distance and alignment length, if we should retain
    //  or discard the overlap.
    //
    //  If it's saved, flag it as 'not useful for unitigging' if the overlap isn't dovetail.

    if ((result.alignmentLength < minOverlapLength) ||                  //  Alignment is too short, or
        (result.editDistance > maxErate * result.alignmentLength)) {    //               too noisy.
      localStats.nFailed++;

      ovl->evalue(AS_MAX_EVALUE);

      ovl->dat.ovl.forOBT = false;
      ovl->dat.ovl.forDUP = false;
      ovl->dat.ovl.forUTG = false;
    }

    else {
      localStats.nPassed++;

      ovl->erate(editDist / (double)alignLen);

      ovl->dat.ovl.forOBT = (partialOverlaps == true);
      ovl->dat.ovl.forDUP = (partialOverlaps == true);
      ovl->dat.ovl.forUTG = (partialOverlaps == false) && (ovl->overlapIsDovetail() == true);
    }

    edlibFreeAlignResult(result);
  }

  //  More logging.

  if (debug) {
    fprintf(stderr, "A %7u %6d-%-6d -> %6d-%-6d\n",       aID, ori.a_bgn(), ori.a_end(), ovl->a_bgn(), ovl->a_end());
    fprintf(stderr, "B %7u %6d-%-6d -> %6d-%-6d %5.4f\n", bID, ori.b_bgn(), ori.b_end(), ovl->b_bgn(), ovl->b_end(), ovl->erate());
  }
}



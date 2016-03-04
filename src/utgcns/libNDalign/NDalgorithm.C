
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
 *  This file is derived from:
 *
 *    src/utgcns/libNDalign/prefixEditDistance.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-JUL-20 to 2015-AUG-05
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-FEB-25
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "NDalgorithm.H"

#include "Binomial_Bound.H"


const char *
toString(pedAlignType at) {
  const char *ret = NULL;

  switch (at) {
    case pedOverlap:  ret = "full-overlap";      break;
    case pedLocal:    ret = "local-overlap";     break;
    case pedGlobal:   ret = "global-alignment";  break;
    default:
      assert(0);
      break;
  }

  return(ret);
}


const char *
toString(pedOverlapType ot) {
  const char *ret = NULL;

  switch (ot) {
    case pedBothBranch:   ret = "both-branch";       break;
    case pedLeftBranch:   ret = "left-branch";       break;
    case pedRightBranch:  ret = "right-branch";      break;
    case pedDovetail:     ret = "dovetail-overlap";  break;
    default:
      assert(0);
      break;
  }

  return(ret);
}


NDalgorithm::NDalgorithm(pedAlignType alignType_, double maxErate_) {
  alignType            = alignType_;
  maxErate             = maxErate_;

  ERRORS_FOR_FREE        = 1;
  MIN_BRANCH_END_DIST    = 20;
  MIN_BRANCH_TAIL_SLOPE  = ((maxErate > 0.06) ? 1.0 : 0.20);

  Left_Delta  = new int32  [AS_MAX_READLEN];
  Right_Delta = new int32  [AS_MAX_READLEN];
  Delta_Stack = new int32  [AS_MAX_READLEN];

  allocated = 3 * AS_MAX_READLEN * sizeof(int32);

  Edit_Space_Max  = AS_MAX_READLEN;  //(alignType == pedGlobal) ? (AS_MAX_READLEN) : (1 + (int32)ceil(maxErate * AS_MAX_READLEN));
  Edit_Space_Lazy = new pedEdit *  [Edit_Space_Max];
  Edit_Array_Lazy = new pedEdit *  [Edit_Space_Max];

  memset(Edit_Space_Lazy, 0, sizeof(pedEdit *) * Edit_Space_Max);
  memset(Edit_Array_Lazy, 0, sizeof(pedEdit *) * Edit_Space_Max);

  allocated += Edit_Space_Max * sizeof (pedEdit *);
  allocated += Edit_Space_Max * sizeof (pedEdit *);

  int32   dataIndex = (int)ceil(maxErate * 100) - 1;

  if ((dataIndex < 0) || (50 <= dataIndex))
    fprintf(stderr, "NDalgorithm()-- Invalid maxErate=%f -> dataIndex=%d\n",
            maxErate, dataIndex);
  assert(0 <= dataIndex);
  assert(dataIndex < 50);

#if 0

  //  Use the precomputed values.
  {
    Edit_Match_Limit_Allocation = NULL;
    Edit_Match_Limit            = Edit_Match_Limit_Data[dataIndex];

    fprintf(stderr, "NDalgorithm()-- Set Edit_Match_Limit to %p; dataIndex=%d 6 = %p\n",
            Edit_Match_Limit, dataIndex, Edit_Match_Limit_0600);
  }

#else

  //  Compute values on the fly.

  {
    int32 MAX_ERRORS = (1 + (int32)ceil(maxErate * AS_MAX_READLEN));

    Edit_Match_Limit_Allocation = new int32 [MAX_ERRORS + 1];

    for (int32 e=0;  e<= ERRORS_FOR_FREE; e++)
      Edit_Match_Limit_Allocation[e] = 0;

    int Start = 1;

    for (int32 e=ERRORS_FOR_FREE + 1; e<MAX_ERRORS; e++) {
      Start = Binomial_Bound(e - ERRORS_FOR_FREE,
                             maxErate,
                             Start);
      Edit_Match_Limit_Allocation[e] = Start - 1;

      assert(Edit_Match_Limit_Allocation[e] >= Edit_Match_Limit_Allocation[e-1]);
    }

    Edit_Match_Limit = Edit_Match_Limit_Allocation;
  }

#endif



  for (int32 i=0; i <= AS_MAX_READLEN; i++) {
    //Error_Bound[i] = (int32) (i * maxErate + 0.0000000000001);
    Error_Bound[i] = (int32)ceil(i * maxErate);
  }


  //  Value to add for a match in finding branch points.
  //
  //  ALH: Note that maxErate also affects what overlaps get found
  //
  //  ALH: Scoring seems to be unusual: given an alignment of length l with k mismatches, the score
  //  seems to be computed as l + k * error value and NOT (l-k)*match+k*error
  //
  //  I.e. letting x := DEFAULT_BRANCH_MATCH_VAL, the max mismatch fraction p to give a non-negative
  //  score would be p = x/(1-x); conversely, to compute x for a goal p, we have x = p/(1+p).  E.g.
  //
  //  for p=0.06, x = .06 / (1.06) = .0566038
  //  for p=0.35, x = .35 / (1.35) = .259259
  //  for p=0.2,  x = .20 / (1.20) = .166667
  //  for p=0.15, x = .15 / (1.15) = .130435
  //
  //  Value was for 6% vs 35% error discrimination.
  //
  //  Converting to integers didn't make it faster.
  //
  //  Corresponding error value is this value minus 1.0

  Branch_Match_Value = maxErate / (1 + maxErate);

  for (uint32 ii=0; ii<256; ii++)
    tolower[ii] = ::tolower(ii);
};



NDalgorithm::~NDalgorithm() {
  delete [] Left_Delta;
  delete [] Right_Delta;

  delete [] Delta_Stack;

  for (uint32 i=0; i<Edit_Space_Max; i++)
    if (Edit_Space_Lazy[i])
      delete [] Edit_Space_Lazy[i];

  delete [] Edit_Space_Lazy;
  delete [] Edit_Array_Lazy;

  delete [] Edit_Match_Limit_Allocation;
};


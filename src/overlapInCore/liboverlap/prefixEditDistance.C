
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

#include "prefixEditDistance.H"

#include "Binomial_Bound.H"


prefixEditDistance::prefixEditDistance(bool doingPartialOverlaps_, double maxErate_) {
  maxErate             = maxErate_;
  doingPartialOverlaps = doingPartialOverlaps_;

  MAX_ERRORS             = (1 + (int)ceil(maxErate * AS_MAX_READLEN));
  MIN_BRANCH_END_DIST    = 20;
  MIN_BRANCH_TAIL_SLOPE  = ((maxErate > 0.06) ? 1.0 : 0.20);

  Left_Delta_Len  = 0;
  Left_Delta      = new int  [MAX_ERRORS];

  Right_Delta_Len = 0;
  Right_Delta     = new int  [MAX_ERRORS];

  allocated  = 3 * MAX_ERRORS * sizeof(int);

  Delta_Stack = new int  [MAX_ERRORS];

  Edit_Space_Lazy = new int *  [MAX_ERRORS];
  Edit_Array_Lazy = new int *  [MAX_ERRORS];

  memset(Edit_Space_Lazy, 0, sizeof(int *) * MAX_ERRORS);
  memset(Edit_Array_Lazy, 0, sizeof(int *) * MAX_ERRORS);

  allocated += MAX_ERRORS * sizeof (int);
  allocated += MAX_ERRORS * sizeof (int);

  //

  Edit_Match_Limit_Allocation = new int32 [MAX_ERRORS + 1];
  Edit_Match_Limit = Edit_Match_Limit_Allocation;

  Initialize_Match_Limit(Edit_Match_Limit_Allocation, maxErate, MAX_ERRORS);


  for (int32 i=0; i <= AS_MAX_READLEN; i++) {
    //Error_Bound[i] = (int32) (i * maxErate + 0.0000000000001);
    Error_Bound[i] = (int32)ceil(i * maxErate);
  }


  //  Value to add for a match in finding branch points.
  //
  //  ALH: Note that maxErate also affects what overlaps get found
  //
  //  ALH: Scoring seems to be unusual: given an alignment
  //  of length l with k mismatches, the score seems to be
  //  computed as l + k * error value and NOT (l-k)*match+k*error
  //
  //  I.e. letting x := DEFAULT_BRANCH_MATCH_VAL,
  //  the max mismatch fraction p to give a non-negative score
  //  would be p = x/(1-x); conversely, to compute x for a
  //  goal p, we have x = p/(1+p).  E.g.
  //
  //  for p=0.06, x = .06 / (1.06) = .0566038
  //  for p=0.35, x = .35 / (1.35) = .259259
  //  for p=0.2,  x = .20 / (1.20) = .166667
  //  for p=0.15, x = .15 / (1.15) = .130435
  //
  //  Value was for 6% vs 35% error discrimination.
  //  Converting to integers didn't make it faster.
  //  Corresponding error value is this value minus 1.0

  Branch_Match_Value = maxErate / (1 + maxErate);
  Branch_Error_Value = Branch_Match_Value - 1.0;
};



prefixEditDistance::~prefixEditDistance() {
  delete [] Left_Delta;
  delete [] Right_Delta;

  delete [] Delta_Stack;

  for (uint32 i=0; i<MAX_ERRORS; i++)
    if (Edit_Space_Lazy[i])
      delete [] Edit_Space_Lazy[i];

  delete [] Edit_Space_Lazy;
  delete [] Edit_Array_Lazy;

  delete [] Edit_Match_Limit_Allocation;
};


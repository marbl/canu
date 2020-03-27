
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

#include "overlapInCore.H"
#include "prefixEditDistance.H"



//  See how far the exact match in  Match  extends.  The match
//  refers to strings  S  and  T  with lengths  S_Len  and  T_Len ,
//  respectively.  Set  S_Lo ,  S_Hi ,  T_Lo  and  T_Hi  to the
//  regions within  S  and  T  to which the match extends.
//  Return the type of overlap:  NONE if doesn't extend to
//  the end of either fragment;  LEFT_BRANCH_PT
//  or  RIGHT_BRANCH_PT  if match extends to the end of just one fragment;
//  or  DOVETAIL  if it extends to the end of both fragments, i.e.,
//  it is a complete overlap.
//  Set  Errors  to the number of errors in the alignment if it is
//  a  DOVETAIL  overlap.

Overlap_t
prefixEditDistance::Extend_Alignment(Match_Node_t *Match,
                                     char         *S,      uint32   S_ID,   int32   S_Len,
                                     char         *T,      uint32   T_ID,   int32   T_Len,
                                     int32        &S_Lo,   int32   &S_Hi,
                                     int32        &T_Lo,   int32   &T_Hi,
                                     int32        &Errors) {
  int32  Right_Errors = 0;
  int32  Left_Errors  = 0;
  int32  Leftover     = 0;

  bool   rMatchToEnd = true;
  bool   lMatchToEnd = true;

  int32  S_Left_Begin  = Match->Start - 1;
  int32  S_Right_Begin = Match->Start + Match->Len;
  int32  S_Right_Len   = S_Len - S_Right_Begin;

  int32  T_Left_Begin  = Match->Offset - 1;
  int32  T_Right_Begin = Match->Offset + Match->Len;
  int32  T_Right_Len   = T_Len - T_Right_Begin;

  int32  Total_Olap = (min(Match->Start, Match->Offset) +
                       Match->Len +
                       min(S_Right_Len, T_Right_Len));

  int32  Error_Limit = Error_Bound[Total_Olap];

#ifdef SHOW_EXTEND_ALIGN
  fprintf(stdout, "prefixEditDistance::Extend_Alignment()--  limit olap of %u bases to %u errors - %f%%\n",
          Total_Olap, Error_Limit, 100.0 * Error_Limit / Total_Olap);
  fprintf(stdout, "prefixEditDistance::Extend_Alignment()--  S: %d-%d and %d-%d  T: %d-%d and %d-%d\n",
          0, S_Left_Begin, S_Right_Begin, S_Right_Begin + S_Right_Len,
          0, T_Left_Begin, T_Right_Begin, T_Right_Begin + T_Right_Len);
#endif

  Left_Delta_Len  = 0;
  Right_Delta_Len = 0;

  if ((S_Right_Len == 0) ||
      (T_Right_Len == 0)) {
    S_Hi        = 0;
    T_Hi        = 0;
    rMatchToEnd = true;
  }

  else if (S_Right_Len <= T_Right_Len) {
    Right_Errors = forward(S + S_Right_Begin, S_Right_Len,
                           T + T_Right_Begin, T_Right_Len,
                           Error_Limit,
                           S_Hi, T_Hi,
                           rMatchToEnd);
    for (int32 i=0; i<Right_Delta_Len; i++)
      Right_Delta[i] *= -1;
  }

  else {
    Right_Errors = forward(T + T_Right_Begin, T_Right_Len,
                           S + S_Right_Begin, S_Right_Len,
                           Error_Limit,
                           T_Hi, S_Hi,
                           rMatchToEnd);
    //for (int32 i=0; i<Right_Delta_Len; i++)
    //  Right_Delta[i] *= -1;
  }

  S_Hi += S_Right_Begin - 1;
  T_Hi += T_Right_Begin - 1;


  if ((S_Left_Begin < 0) ||
      (T_Left_Begin < 0)) {
    S_Lo        = 0;
    T_Lo        = 0;
    lMatchToEnd = true;
  }

  else if (S_Right_Begin <= T_Right_Begin) {
    Left_Errors = reverse(S + S_Left_Begin, S_Left_Begin + 1,
                          T + T_Left_Begin, T_Left_Begin + 1,
                          Error_Limit - Right_Errors,
                          S_Lo, T_Lo,
                          Leftover,
                          lMatchToEnd);
    //for (int32 i=0; i<Left_Delta_Len; i++)
    //  Left_Delta[i] *= -1;
  }

  else {
    Left_Errors = reverse(T + T_Left_Begin,  T_Left_Begin + 1,
                          S + S_Left_Begin,  S_Left_Begin + 1,
                          Error_Limit - Right_Errors,
                          T_Lo, S_Lo,
                          Leftover,
                          lMatchToEnd);
    for (int32 i=0; i<Left_Delta_Len; i++)
      Left_Delta[i] *= -1;
  }

  S_Lo += S_Left_Begin + 1;
  T_Lo += T_Left_Begin + 1;

  Errors = Left_Errors + Right_Errors;

  if ((Right_Errors > Error_Limit) ||
      (Left_Errors  > Error_Limit) ||
      (Errors       > Error_Limit))
    fprintf(stderr, "Extend_Alignment()- S_ID %u T_ID %u Error_Limit %d Left_Errors %d Right_Errors %d\n",
            S_ID, T_ID, Error_Limit, Left_Errors, Right_Errors);

  assert(Right_Errors <= Error_Limit);
  assert(Left_Errors  <= Error_Limit);
  assert(Errors       <= Error_Limit);

#ifdef SHOW_EXTEND_ALIGN
  fprintf(stdout, "WorkArea %2d OVERLAP %6d %6d - %5d + %5d = %5d errors out of %d possible\n", omp_get_thread_num(), S_ID, T_ID, Left_Errors, Right_Errors, Errors, Edit_Space_Lazy_Max);
#endif

  //  No overlap if both right and left don't match to end, otherwise a branch point if only one.
  //  If both match to end, a dovetail overlap.  Indenting is all screwed up here.

  Overlap_t return_type = (rMatchToEnd == false) ? ((lMatchToEnd == false) ? NONE : RIGHT_BRANCH_PT) :
                                                   ((lMatchToEnd == false) ? LEFT_BRANCH_PT : DOVETAIL);

  //  Merge the deltas.  Previously, this wouldn't be done if the overlap wasn't dovetail (and not partial).


  if (Right_Delta_Len > 0) {
    if (Right_Delta[0] > 0)
      Left_Delta[Left_Delta_Len++] = -(Right_Delta[0] + Leftover + Match->Len);
    else
      Left_Delta[Left_Delta_Len++] = -(Right_Delta[0] - Leftover - Match->Len);
  }

  //  WHY?!  Does this mean the invesion on the forward() calls is backwards?
  for (int32 i=1; i<Right_Delta_Len; i++)
    Left_Delta[Left_Delta_Len++] = -Right_Delta[i];

  Right_Delta_Len = 0;  //  Copied into left_delta!

  //  Return.

  return(return_type);
}

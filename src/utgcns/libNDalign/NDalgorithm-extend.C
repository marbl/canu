
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
 *    src/utgcns/libNDalign/prefixEditDistance-extend.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-JUL-20 to 2015-JUL-28
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-13
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "NDalgorithm.H"

#undef  DEBUG


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

pedOverlapType
NDalgorithm::Extend_Alignment(Match_Node_t *Match,
                              char         *S,     int32   S_Len,
                              char         *T,     int32   T_Len,
                              int32        &S_Lo,  int32   &S_Hi,
                              int32        &T_Lo,  int32   &T_Hi) {
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

#ifdef DEBUG
  fprintf(stderr, "NDalgorithm::Extend_Alignment()--  S: %d-%d and %d-%d  T: %d-%d and %d-%d\n",
          0, S_Left_Begin, S_Right_Begin, S_Right_Begin + S_Right_Len,
          0, T_Left_Begin, T_Right_Begin, T_Right_Begin + T_Right_Len);
#endif

  Right_Score       = 0;
  Right_Delta_Len   = 0;

  Left_Score        = 0;
  Left_Delta_Len    = 0;

  bool   invertLeftDeltas  = false;
  bool   invertRightDeltas = false;


  if ((S_Right_Len == 0) ||
      (T_Right_Len == 0)) {
    S_Hi        = 0;
    T_Hi        = 0;
    rMatchToEnd = true;
  }

  else if (S_Right_Len <= T_Right_Len) {
#ifdef DEBUG
    fprintf(stderr, "NDalgorithm::Extend_Alignment()--  FORWARD S T\n");
#endif
    forward(S + S_Right_Begin, S_Right_Len,
            T + T_Right_Begin, T_Right_Len,
            S_Hi, T_Hi,
            rMatchToEnd);
    for (int32 i=0; i<Right_Delta_Len; i++)
      Right_Delta[i] *= -1;
  }

  else {
#ifdef DEBUG
    fprintf(stderr, "NDalgorithm::Extend_Alignment()--  FORWARD T S\n");
#endif
    forward(T + T_Right_Begin, T_Right_Len,
            S + S_Right_Begin, S_Right_Len,
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
#ifdef DEBUG
    fprintf(stderr, "NDalgorithm::Extend_Alignment()--  REVERSE S T\n");
#endif
    reverse(S + S_Left_Begin, S_Left_Begin + 1,
            T + T_Left_Begin, T_Left_Begin + 1,
            S_Lo, T_Lo,
            Leftover,
            lMatchToEnd);
    //for (int32 i=0; i<Left_Delta_Len; i++)
    //  Left_Delta[i] *= -1;
  }

  else {
#ifdef DEBUG
    fprintf(stderr, "NDalgorithm::Extend_Alignment()--  REVERSE T S\n");
#endif
    reverse(T + T_Left_Begin,  T_Left_Begin + 1,
            S + S_Left_Begin,  S_Left_Begin + 1,
            T_Lo, S_Lo,
            Leftover,
            lMatchToEnd);
    for (int32 i=0; i<Left_Delta_Len; i++)
      Left_Delta[i] *= -1;
  }

  S_Lo += S_Left_Begin + 1;
  T_Lo += T_Left_Begin + 1;

  //  Report.

#ifdef DEBUG
  fprintf(stderr, "NDalgorithm::Extend_Alignment()--  LEFT deltaLen %d matchToEnd %s leftover %d -- RIGHT deltaLen %d matchToEnd %s\n",
          Left_Delta_Len,  lMatchToEnd ? "true" : "false", Leftover,
          Right_Delta_Len, rMatchToEnd ? "true" : "false");
#endif

  //  Merge the deltas.

  if (Right_Delta_Len > 0) {
    if (Right_Delta[0] > 0)
      Left_Delta[Left_Delta_Len++] = -(Right_Delta[0] + Leftover + Match->Len);
    else
      Left_Delta[Left_Delta_Len++] = -(Right_Delta[0] - Leftover - Match->Len);
  }

  //  WHY?!  Does this mean the invesion on the forward() calls is backwards?
  //  But note interaction with the if test just above here!
  for (int32 i=1; i<Right_Delta_Len; i++)
    Left_Delta[Left_Delta_Len++] = -Right_Delta[i];

  Right_Delta_Len = 0;  //  Copied into left_delta!

  //  Return.

  return((rMatchToEnd == false) ? ((lMatchToEnd == false) ? pedBothBranch : pedRightBranch) :
                                  ((lMatchToEnd == false) ? pedLeftBranch : pedDovetail));
}


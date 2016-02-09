
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
 *    src/overlapInCore/prefixEditDistance.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-FEB-11 to 2015-JUL-20
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef PREFIX_EDIT_DISTANCE_H
#define PREFIX_EDIT_DISTANCE_H


#include "AS_global.H"
#include "gkStore.H"  //  For AS_MAX_READLEN


//  Used in -forward and -reverse
#define Sign(a) ( ((a) > 0) - ((a) < 0) )



enum Overlap_t {
  NONE,
  LEFT_BRANCH_PT,
  RIGHT_BRANCH_PT,
  DOVETAIL
};


//  the input to Extend_Alignment.
struct Match_Node_t {
  int32  Offset;              // To start of exact match in  hash-table frag
  int32  Len;                 // Of exact match
  int32  Start;               // Of exact match in current (new) frag
  int32  Next;                // Subscript of next match in list
};




class prefixEditDistance {
public:
  prefixEditDistance(bool doingPartialOverlaps_, double maxErate_);
  ~prefixEditDistance();

  void   Allocate_More_Edit_Space(void);

  void   Set_Right_Delta(int32 e, int32 d);
  int32  forward(char    *A,   int32 m,
                 char    *T,   int32 n,
                 int32    Error_Limit,
                 int32   &A_End,
                 int32   &T_End,
                 bool    &Match_To_End);


  void   Set_Left_Delta(int32 e, int32 d, int32 &leftover, int32 &t_end, int32 t_len);
  int32  reverse(char    *A,   int32 m,
                 char    *T,   int32 n,
                 int32    Error_Limit,
                 int32   &A_End,
                 int32   &T_End,
                 int32   &Leftover,      //  <- novel
                 bool    &Match_To_End);


  Overlap_t  Extend_Alignment(Match_Node_t  *Match,
                              char          *S,     int32    S_Len,
                              char          *T,     int32    T_Len,
                              int32         &S_Lo,  int32   &S_Hi,
                              int32         &T_Lo,  int32   &T_Hi,
                              int32         &Errors);

public:
  //  The four below were global #defines, two depended on the error rate which is now local.

  //  Most errors in any edit distance computation
  uint32   MAX_ERRORS;

  //  Branch points must be at least this many bases from the end of the fragment to be reported
  uint32   MIN_BRANCH_END_DIST;


  //  Branch point tails must fall off from the max by at least this rate
  double   MIN_BRANCH_TAIL_SLOPE;

  double   maxErate;
  bool     doingPartialOverlaps;

  uint64   allocated;

  int32    Left_Delta_Len;
  int32   *Left_Delta;

  int32    Right_Delta_Len;
  int32   *Right_Delta;

  int32   *Delta_Stack;

  int32  **Edit_Space_Lazy;        //  Array of pointers, if set, it was a new'd allocation
  int32  **Edit_Array_Lazy;        //  Array of pointers, some are not new'd allocations

  //  This array [e] is the minimum value of  Edit_Array[e][d]
  //  to be worth pursuing in edit-distance computations between reads
  const
  int32   *Edit_Match_Limit;
  int32   *Edit_Match_Limit_Allocation;

  //  The maximum number of errors allowed in a match between reads of length i,
  //  which is i * AS_OVL_ERROR_RATE.
  int32    Error_Bound[AS_MAX_READLEN + 1];

  //  Scores of matches and mismatches in alignments.  Alignment ends at maximum score.
  double   Branch_Match_Value;
  double   Branch_Error_Value;

};


#endif

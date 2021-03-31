
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

#ifndef PREFIX_EDIT_DISTANCE_H
#define PREFIX_EDIT_DISTANCE_H


#include "types.H"
#include "sqStore.H"  //  For AS_MAX_READLEN


#undef  DEBUG_EDIT_SPACE_ALLOC
#undef  SHOW_EXTEND_ALIGN
#undef  SHOW_BRI

//  Used in -forward and -reverse
#define Sign(a) ( ((a) > 0) - ((a) < 0) )



enum Overlap_t {
  NONE,
  LEFT_BRANCH_PT,
  RIGHT_BRANCH_PT,
  DOVETAIL
};

static
const
char *
toString(Overlap_t k) {
  switch (k) {
    case  NONE:             return("NONE");             break;
    case  LEFT_BRANCH_PT:   return("LEFT_BRANCH_PT");   break;
    case  RIGHT_BRANCH_PT:  return("RIGHT_BRANCH_PT");  break;
    case  DOVETAIL:         return("DOVETAIL");         break;
    default:                return("UNKNOWN");          break;
  }
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
  prefixEditDistance(bool doingPartialOverlaps_, double maxErate_, double maxAlignErate_);
  ~prefixEditDistance();

  void   Allocate_More_Edit_Space(int e);

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


  Overlap_t  Extend_Alignment(Match_Node_t *Match,
                              char         *S,      uint32   S_ID,   int32   S_Len,
                              char         *T,      uint32   T_ID,   int32   T_Len,
                              int32        &S_Lo,   int32   &S_Hi,
                              int32        &T_Lo,   int32   &T_Hi,
                              int32        &Errors);


  //  pruneAlign() returns true if this DP path should be pruned, false if
  //  it should be kept.  It doesn't quite work; the original method is
  //  faster for PacBio and sufficient for HiFi.
private:

#if 0
  //  HiFi - this works pretty good, in that it finds the one test overlap
  //  that the original didn't find.  Not tested for performance otherwise.
  double    alignFraction1    = 0.10;   //  current align len  0 .. f1
  double    alignFraction2    = 0.75;   //  current align len f1 .. f2
  double    alignFraction3    = 1.00;   //  current align len f2 .. f3  (assumed to be 1.0)

  double    errorFraction1    = 1.0;    //  prune if currentErrorRateabove maximum of
  double    errorFraction2    = 0.0;    //   f1 fraction error
  double    errorFraction3    = 0.0;    //   m1 * maxErate

  double    errorMultiplier1  = 5.0;    //
  double    errorMultiplier2  = 1.5;    //
  double    errorMultiplier3  = 1.1;    //
#endif

#if 0
  //  PacBio - fast but misses overlaps.
  //  3938.950u 1.810s 4:54.92 1336.2%        365+172k 2+42io 4pf+0w
  double    alignFraction1    = 0.01;   //  current align len  0 .. f1
  double    alignFraction2    = 0.10;   //  current align len f1 .. f2
  double    alignFraction3    = 1.00;   //  current align len f2 .. f3  (assumed to be 1.0)

  double    errorFraction1    = 0.0;    //  prune if currentErrorRateabove maximum of
  double    errorFraction2    = 0.0;    //   f1 fraction error
  double    errorFraction3    = 0.0;    //   m1 * maxErate

  double    errorMultiplier1  = 3.0;    //
  double    errorMultiplier2  = 1.5;    //
  double    errorMultiplier3  = 1.0;    //
#endif

#if 0
  //  PacBio - fast but misses overlaps.
  //  4051.327u 1.897s 5:00.32 1349.6%        365+172k 2+36io 4pf+0w
  double    alignFraction1    = 0.02;   //  current align len  0 .. f1
  double    alignFraction2    = 0.06;   //  current align len f1 .. f2
  double    alignFraction3    = 1.00;   //  current align len f2 .. f3  (assumed to be 1.0)

  double    errorFraction1    = 0.0;    //  prune if currentErrorRateabove maximum of
  double    errorFraction2    = 0.0;    //   f1 fraction error
  double    errorFraction3    = 0.0;    //   m1 * maxErate

  double    errorMultiplier1  = 3.0;    //
  double    errorMultiplier2  = 1.2;    //
  double    errorMultiplier3  = 1.0;    //
#endif

#if 1
  //  PacBio - no pruning, finds lots of overlaps.
  //  266216.973u 4.504s 4:47:37.12 1542.6%   365+172k 2+39io 4pf+0w
  double    alignFraction1    = 0.50;   //  current align len  0 .. f1
  double    alignFraction2    = 0.90;   //  current align len f1 .. f2
  double    alignFraction3    = 1.00;   //  current align len f2 .. f3  (assumed to be 1.0)

  double    errorFraction1    = 1.0;    //  prune if currentErrorRateabove maximum of
  double    errorFraction2    = 0.0;    //   f1 fraction error
  double    errorFraction3    = 0.0;    //   m1 * maxErate

  double    errorMultiplier1  = 3.0;    //
  double    errorMultiplier2  = 4.0;    //
  double    errorMultiplier3  = 2.0;    //
#endif


  //  Returns a percent error that we'd consider good for a partial alignment
  //  of length EAL[e][i] given the full alignment is expected to be
  //  totalAlignLen.
  //
  //    |
  //    |-------+-
  //    |         --
  //    |           --
  //    |             --
  //    |               -+---
  //    |                    ---
  //    |  R1   |   R2   |   R3 ---
  //    +-------+---------+-----------+--------
  //           aF1       aF2         aF3
  //
  double  pruneAlign_pos(double l, double m, double r)   { return( (m-l) / (r-l) ); }

  double
  pruneAlign_pe(int32 e, int32 i, int32 extra, int32 totalAlignLen) {
    int32  currentAlignLen  = Edit_Array_Lazy[e][i] + extra;
    double alignFraction    = (double)currentAlignLen / totalAlignLen;
    double allowedErrorRate = maxErate;

    //  Region 1 is expected to be constant with 'l' and 'r' both set to 1.0
    if      (alignFraction <= alignFraction1) {
      double  l = 1.00;
      double  r = std::max(errorFraction1, errorMultiplier1 * maxErate);
      double  p = pruneAlign_pos(0, alignFraction, alignFraction1);

      allowedErrorRate = p * (r - l) + l;
    }

    //  Region 2 should be decreasing from 'l' (at aF1) to 'r' (at aF2).
    else if (alignFraction <= alignFraction2) {
      double  l = std::max(errorFraction1, errorMultiplier1 * maxErate);
      double  r = std::max(errorFraction2, errorMultiplier2 * maxErate);
      double  p = pruneAlign_pos(alignFraction1, alignFraction, alignFraction2);

      allowedErrorRate = p * (r - l) + l;
    }

    //  Region 3 is similar.
    else if (alignFraction <= alignFraction3) {
      double  l = std::max(errorFraction2, errorMultiplier2 * maxErate);
      double  r = std::max(errorFraction3, errorMultiplier3 * maxErate);
      double  p = pruneAlign_pos(alignFraction2, alignFraction, alignFraction3);

      allowedErrorRate = p * (r - l) + l;
    }

    //  Alignments longer than expected must be below maxErate.
    else {
      allowedErrorRate = maxErate;
    }

    return(allowedErrorRate);
  }

  //  Return the length of alignment needed to get 'e' errors at 'maxErate'
  //  percent error.  Just rearranging:  errorFraction = errors / alignLen
  //  to get:                            minLen        = errors / maxErate
  //
  int32
  pruneAlign_ml(int32 e, int32 i, int32 extra, int32 totalAlignLen) {
    return(e / maxErate);
  }

  //  Return the previous minimum alignment length needed for 'e' errors.
  //
  int32
  pruneAlign_ML(int32 e, int32 i, int32 extra, int32 totalAlignLen) {
    return(Edit_Match_Limit[e]);
  }

  //  Returns true if the current alignment length is below what we'd
  //  consider to be a good alignment.
  //
  bool
  pruneAlign(int32 e,                 //  Index into EAL[e][i]  (number of errors)
             int32 i,                 //  Index into EAL[e][i]  (diagonal)
             int32 extra,             //
             int32 totalAlignLen) {   //  Expected total alignment length

#if 1

    //  Original method.
    return(Edit_Array_Lazy[e][i] + extra < Edit_Match_Limit[e]);

#else

    int32  currentAlignLen      = Edit_Array_Lazy[e][i] + extra;
    double currentErrorFraction = (double)e / currentAlignLen;

    //  Don't prune if low error.
    if (currentErrorFraction < pruneAlign_pe(e, i, extra, totalAlignLen))
      return(false);

    //  Don't prune if long.
    if (currentAlignLen > pruneAlign_ml(e, i, extra, totalAlignLen))
      return(false);

    //  Short and low quality.  Bye.
    return(true);

#endif
  }


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

#ifdef DEBUG_EDIT_SPACE_ALLOC
  int32    Edit_Space_Lazy_Max;  //  Last allocated block, DEBUG ONLY
#endif

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

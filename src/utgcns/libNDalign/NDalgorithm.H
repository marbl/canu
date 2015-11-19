
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
 *    src/utgcns/libNDalign/prefixEditDistance.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-JUL-20 to 2015-AUG-05
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-19
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef NDALGORITHM_H
#define NDALGORITHM_H


#include "AS_global.H"
#include "gkStore.H"  //  For AS_MAX_READLEN


//  Used in -forward and -reverse
#define Sign(a) ( ((a) > 0) - ((a) < 0) )

//  Scores used in the dynamic programming.
//
#define PEDMATCH        2    //  Match from UPPERCASE to UPPERCASE
#define PEDGAPMATCH     1    //  Match from LOWERCASE to UPPERCASE
#define PEDMISMATCH    -3    //  Mismatch
#define PEDFREEGAP      1    //  Gap from LOWERCASE
#define PEDGAP         -5    //  Gap

//  PEDMINSCORE wants to be INT32_MIN, except that it will underflow, and wrap around to the highest
//  score possible, on the first insert.  I think it suffices to set this just slightly higher than
//  the minimum, since it should never be picked as a maximum, and thus never see more than one
//  insert.  Just to be safe, we give it lots of space.
//
#define PEDMINSCORE   (INT32_MIN - PEDGAP * 32768)

enum pedAlignType {
  pedOverlap     = 0,  //  Former normal overlap
  pedLocal       = 1,  //  Former partial overlap
  pedGlobal      = 2   //  New forced global overlap
};

enum pedOverlapType {
  pedBothBranch  = 0,
  pedLeftBranch  = 1,
  pedRightBranch = 2,
  pedDovetail    = 3
};

const char *toString(pedAlignType at);
const char *toString(pedOverlapType ot);


//  the input to Extend_Alignment.
class Match_Node_t {
public:
  int32  Offset;              // To start of exact match in  hash-table frag
  int32  Len;                 // Of exact match
  int32  Start;               // Of exact match in current (new) frag
  int32  Next;                // Subscript of next match in list
};



class pedEdit {
public:
  int32  row;     //  Position in the A sequence; was previously the score too
  int32  dist;    //  The number of non-free bases in the alignment, used to decide band
  int32  errs;    //  The number of non-free errors in the alignment, used to decide band
  int32  score;   //  The dynamic programming score, used to decide best align
  int32  fromd;   //  For backtracking, where we came from

  void   init(void) {
    row   = -2;
    dist  = -2;
    errs  =  0;
    score =  PEDMINSCORE;
    fromd =  INT32_MAX;
  };

  void   display(int32 e, int32 d) {
    fprintf(stderr, "e=%d d=%d - ptr %p - ", e, d, this);
    fprintf(stderr, "row=%d ", row);
    fprintf(stderr, "dist=%d ", dist);
    fprintf(stderr, "errs=%d ", errs);
    fprintf(stderr, "score=%d ", score);
    fprintf(stderr, "fromd=%d\n", fromd);
  };
};



class NDalgorithm {
public:
  NDalgorithm(pedAlignType alignType_, double maxErate_);
  ~NDalgorithm();

private:
  bool   allocateMoreEditSpace(void);

  void   Set_Right_Delta(int32  e, int32  d);

  void   forward(char    *A,   int32 m,
                 char    *T,   int32 n,
                 int32   &A_End,
                 int32   &T_End,
                 bool    &Match_To_End);


  void   Set_Left_Delta(int32  e, int32 d,
                        int32 &leftover);

  void   reverse(char    *A,   int32 m,
                 char    *T,   int32 n,
                 int32   &A_End,
                 int32   &T_End,
                 int32   &Leftover,      //  <- novel
                 bool    &Match_To_End);


  //  Returns true if letter 'a' from sequence A matches letter 't' from sequence T.
  //  Wanted to allow lowercase as free matches, but the O(ND) algorithm doesn't support that.
  //
  bool   isMatch(char a, char t) {
    return(a == t);
  };

  //  Returns the score of a match.  Pretty basic, was written to support 'a' == 'A' matches, but
  //  the O(ND) algorithm cannot support that.
  //
  int32  matchScore(char a, char t) {

    if (isMatch(a, t))
      return(PEDMATCH);

    fprintf(stderr, "ERROR: a=%c/%d does not match t=%c/%d\n", a, a, t, t);
    assert(isMatch(a, t) == true);
    return(0);
  };

  //  Returns the score of a mismatch.  As with gaps, we can get free mismatches to lower case letters,
  //  which should prevent alignments such as:
  //
  //  A: GtTTgaTGACCTGG  (the scoring should be 2 mismatches then 3 gaps)
  //  B: GTTTGA---CCTGG
  //
  //  Favoring instead
  //
  //  A: GtTTgaTGACCTGG  (the scoring should be 2 free gaps, 1 gap and 2 matches)
  //  B: GTTT---GACCTGG
  //
  int32  mismatchScore(char a, char t) {
    if ((a == tolower[a]) || (t == tolower[t]))
      return(PEDGAPMATCH);

    return(PEDMISMATCH);
  }


  //  Returns true if the gap we'd be adding to T is zero cost:
  //    lower case 'a' is a don't case position in sequence A
  //
  bool   isFreeGap(char a) {
    if (a == 0)
      return(false);

    return(a == tolower[a]);
  }

public:
  pedOverlapType  Extend_Alignment(Match_Node_t *Match,
                                   char         *S,     int32   S_Len,
                                   char          *T,     int32   T_Len,
                                   int32         &S_Lo,  int32   &S_Hi,
                                   int32         &T_Lo,  int32   &T_Hi);

  int32           score(void) {  return(Left_Score + Right_Score);  };

public:
  //  The three below were global #defines, two depended on the error rate which is now local.

  //  The number of errors that are ignored in setting probability bound for terminating alignment
  //  extensions in edit distance calculations
  uint32                  ERRORS_FOR_FREE;

  //  Branch points must be at least this many bases from the end of the fragment to be reported
  uint32                  MIN_BRANCH_END_DIST;

  //  Branch point tails must fall off from the max by at least this rate
  double                  MIN_BRANCH_TAIL_SLOPE;


  pedAlignType            alignType;
  double                  maxErate;

  uint64                  allocated;

  int32                   Left_Score;
  int32                   Left_Delta_Len;
  int32                  *Left_Delta;

  int32                   Right_Score;
  int32                   Right_Delta_Len;
  int32                  *Right_Delta;

  int32                  *Delta_Stack;

  int32                   Edit_Space_Max;
  pedEdit               **Edit_Space_Lazy;        //  Array of pointers, if set, it was a new'd allocation
  pedEdit               **Edit_Array_Lazy;        //  Array of pointers, some are not new'd allocations

  //  This array [e] is the minimum value of  Edit_Array[e][d]
  //  to be worth pursuing in edit-distance computations between reads
  const
  int32                  *Edit_Match_Limit;
  int32                  *Edit_Match_Limit_Allocation;

  //  The maximum number of errors allowed in a match between reads of length i,
  //  which is i * AS_OVL_ERROR_RATE.
  int32                   Error_Bound[AS_MAX_READLEN + 1];

  //  Scores of matches and mismatches in alignments.  Alignment ends at maximum score.
  double                  Branch_Match_Value;

  char                    tolower[256];
};


#endif  //  NDALGORITHM_H

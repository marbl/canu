

#ifndef PREFIX_EDIT_DISTANCE_H
#define PREFIX_EDIT_DISTANCE_H


#include "AS_global.H"
#include "gkStore.H"
#include "Binomial_Bound.H"


#define Sign(a) ( ((a) > 0) - ((a) < 0) ) 


class prefixEditDistance {
public:
  prefixEditDistance(bool doingPartialOverlaps_, double maxErate_) {
    maxErate             = maxErate_;
    doingPartialOverlaps = doingPartialOverlaps_;

    MAX_ERRORS             = (1 + (int) (maxErate * AS_MAX_READLEN));
    ERRORS_FOR_FREE        = 1;
    MIN_BRANCH_END_DIST    = 20;
    MIN_BRANCH_TAIL_SLOPE  = ((maxErate > 0.06) ? 1.0 : 0.20);

    Left_Delta  = new int  [MAX_ERRORS];
    Right_Delta = new int  [MAX_ERRORS];

    allocated += 3 * MAX_ERRORS * sizeof(int);

    Delta_Stack = new int  [MAX_ERRORS];

    Edit_Space_Lazy = new int *  [MAX_ERRORS];
    Edit_Array_Lazy = new int *  [MAX_ERRORS];

    memset(Edit_Space_Lazy, 0, sizeof(int *) * MAX_ERRORS);
    memset(Edit_Array_Lazy, 0, sizeof(int *) * MAX_ERRORS);

    allocated += MAX_ERRORS * sizeof (int);
    allocated += MAX_ERRORS * sizeof (int);


    for (int32 e=0;  e<= ERRORS_FOR_FREE; e++)
      Edit_Match_Limit[e] = 0;

    int Start = 1;

    for (int32 e=ERRORS_FOR_FREE + 1; e<MAX_ERRORS; e++) {
      Start = Binomial_Bound(e - ERRORS_FOR_FREE,
                             maxErate,
                             Start);
      Edit_Match_Limit[e] = Start - 1;

      assert(Edit_Match_Limit[e] >= Edit_Match_Limit[e-1]);
    }

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

  ~prefixEditDistance() {
    delete [] Left_Delta;
    delete [] Right_Delta;

    delete [] Delta_Stack;

    for (uint32 i=0; i<MAX_ERRORS; i++)
      if (Edit_Space_Lazy[i])
        delete [] Edit_Space_Lazy[i];

    delete [] Edit_Space_Lazy;
    delete [] Edit_Array_Lazy;
  };


  void  Allocate_More_Edit_Space(void);


  void  Set_Right_Delta(int32 e, int32 d);

  int32  forward(char    *A,   int32 m,
                 char    *T,   int32 n,
                 int32    Error_Limit,
                 int32   &A_End,
                 int32   &T_End,
                 bool    &Match_To_End);


  void  Set_Left_Delta(int32 e, int32 d, int32 &leftover, int32 &t_end, int32 t_len);

  int32  reverse(char    *A,   int32 m,
                 char    *T,   int32 n,
                 int32    Error_Limit,
                 int32   &A_End,
                 int32   &T_End,
                 int32   &Leftover,      //  <- novel
                 bool    &Match_To_End);


public:
  //  The four below were global #defines, two depended on the error rate which is now local.

  //  Most errors in any edit distance computation
  uint32  MAX_ERRORS;

  //  The number of errors that are ignored in setting probability bound for terminating alignment
  //  extensions in edit distance calculations
  uint32  ERRORS_FOR_FREE;

  //  Branch points must be at least this many bases from the end of the fragment to be reported
  uint32  MIN_BRANCH_END_DIST;


  //  Branch point tails must fall off from the max by at least this rate
  double  MIN_BRANCH_TAIL_SLOPE;

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
  //  (only MAX_ERRORS needed)
  int32  Edit_Match_Limit[AS_MAX_READLEN];

  //  The maximum number of errors allowed in a match between reads of length i,
  //  which is i * AS_OVL_ERROR_RATE.
  int32  Error_Bound[AS_MAX_READLEN + 1];

  //  Scores of matches and mismatches in alignments.  Alignment ends at maximum score.
  double  Branch_Match_Value;
  double  Branch_Error_Value;
  
};


#endif

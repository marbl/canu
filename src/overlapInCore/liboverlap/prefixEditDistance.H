

#ifndef PREFIX_EDIT_DISTANCE_H
#define PREFIX_EDIT_DISTANCE_H


#include "AS_global.H"
#include "gkStore.H"
#include "Binomial_Bound.H"


#define Sign(a) ( ((a) > 0) - ((a) < 0) ) 


extern const int32 *Edit_Match_Limit_0100;
extern const int32 *Edit_Match_Limit_0200;
extern const int32 *Edit_Match_Limit_0300;
extern const int32 *Edit_Match_Limit_0400;
extern const int32 *Edit_Match_Limit_0500;
extern const int32 *Edit_Match_Limit_0600;
extern const int32 *Edit_Match_Limit_0700;
extern const int32 *Edit_Match_Limit_0800;
extern const int32 *Edit_Match_Limit_0900;
extern const int32 *Edit_Match_Limit_1000;
extern const int32 *Edit_Match_Limit_1100;
extern const int32 *Edit_Match_Limit_1200;
extern const int32 *Edit_Match_Limit_1300;
extern const int32 *Edit_Match_Limit_1400;
extern const int32 *Edit_Match_Limit_1500;
extern const int32 *Edit_Match_Limit_1600;
extern const int32 *Edit_Match_Limit_1700;
extern const int32 *Edit_Match_Limit_1800;
extern const int32 *Edit_Match_Limit_1900;
extern const int32 *Edit_Match_Limit_2000;
extern const int32 *Edit_Match_Limit_2100;
extern const int32 *Edit_Match_Limit_2200;
extern const int32 *Edit_Match_Limit_2300;
extern const int32 *Edit_Match_Limit_2400;
extern const int32 *Edit_Match_Limit_2500;
extern const int32 *Edit_Match_Limit_2600;
extern const int32 *Edit_Match_Limit_2700;
extern const int32 *Edit_Match_Limit_2800;
extern const int32 *Edit_Match_Limit_2900;
extern const int32 *Edit_Match_Limit_3000;
extern const int32 *Edit_Match_Limit_3100;
extern const int32 *Edit_Match_Limit_3200;
extern const int32 *Edit_Match_Limit_3300;
extern const int32 *Edit_Match_Limit_3400;
extern const int32 *Edit_Match_Limit_3500;
extern const int32 *Edit_Match_Limit_3600;
extern const int32 *Edit_Match_Limit_3700;
extern const int32 *Edit_Match_Limit_3800;
extern const int32 *Edit_Match_Limit_3900;
extern const int32 *Edit_Match_Limit_4000;
extern const int32 *Edit_Match_Limit_4100;
extern const int32 *Edit_Match_Limit_4200;
extern const int32 *Edit_Match_Limit_4300;
extern const int32 *Edit_Match_Limit_4400;
extern const int32 *Edit_Match_Limit_4500;
extern const int32 *Edit_Match_Limit_4600;
extern const int32 *Edit_Match_Limit_4700;
extern const int32 *Edit_Match_Limit_4800;
extern const int32 *Edit_Match_Limit_4900;
extern const int32 *Edit_Match_Limit_5000;


static
const
int32 *
Edit_Match_Limit_Data[50] = {
  Edit_Match_Limit_0100,
  Edit_Match_Limit_0200,
  Edit_Match_Limit_0300,
  Edit_Match_Limit_0400,
  Edit_Match_Limit_0500,
  Edit_Match_Limit_0600,
  Edit_Match_Limit_0700,
  Edit_Match_Limit_0800,
  Edit_Match_Limit_0900,
  Edit_Match_Limit_1000,
  Edit_Match_Limit_1100,
  Edit_Match_Limit_1200,
  Edit_Match_Limit_1300,
  Edit_Match_Limit_1400,
  Edit_Match_Limit_1500,
  Edit_Match_Limit_1600,
  Edit_Match_Limit_1700,
  Edit_Match_Limit_1800,
  Edit_Match_Limit_1900,
  Edit_Match_Limit_2000,
  Edit_Match_Limit_2100,
  Edit_Match_Limit_2200,
  Edit_Match_Limit_2300,
  Edit_Match_Limit_2400,
  Edit_Match_Limit_2500,
  Edit_Match_Limit_2600,
  Edit_Match_Limit_2700,
  Edit_Match_Limit_2800,
  Edit_Match_Limit_2900,
  Edit_Match_Limit_3000,
  Edit_Match_Limit_3100,
  Edit_Match_Limit_3200,
  Edit_Match_Limit_3300,
  Edit_Match_Limit_3400,
  Edit_Match_Limit_3500,
  Edit_Match_Limit_3600,
  Edit_Match_Limit_3700,
  Edit_Match_Limit_3800,
  Edit_Match_Limit_3900,
  Edit_Match_Limit_4000,
  Edit_Match_Limit_4100,
  Edit_Match_Limit_4200,
  Edit_Match_Limit_4300,
  Edit_Match_Limit_4400,
  Edit_Match_Limit_4500,
  Edit_Match_Limit_4600,
  Edit_Match_Limit_4700,
  Edit_Match_Limit_4800,
  Edit_Match_Limit_4900,
  Edit_Match_Limit_5000
};


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
  prefixEditDistance(bool doingPartialOverlaps_, double maxErate_) {
    maxErate             = ceil(maxErate_);
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

    fprintf(stderr, "maxErate %f -- %d\n", maxErate, (int)(maxErate * 100));

    Edit_Match_Limit = Edit_Match_Limit_Data[(int)(maxErate * 100)];

#if 0
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
#endif

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


  Overlap_t
  Extend_Alignment(Match_Node_t *Match,
                   char         *S,     int32   S_Len,
                   char          *T,     int32   T_Len,
                   int32         &S_Lo,  int32   &S_Hi,
                   int32         &T_Lo,  int32   &T_Hi,
                   int32         &Errors,
                   bool           partialOverlaps);

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
  const
  int32   *Edit_Match_Limit;

  //  The maximum number of errors allowed in a match between reads of length i,
  //  which is i * AS_OVL_ERROR_RATE.
  int32  Error_Bound[AS_MAX_READLEN + 1];

  //  Scores of matches and mismatches in alignments.  Alignment ends at maximum score.
  double  Branch_Match_Value;
  double  Branch_Error_Value;
  
};


#endif

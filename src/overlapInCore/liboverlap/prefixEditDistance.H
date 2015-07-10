

#ifndef PREFIX_EDIT_DISTANCE_H
#define PREFIX_EDIT_DISTANCE_H


#include "AS_global.H"
#include "gkStore.H"  //  For AS_MAX_READLEN


//  Used in -forward and -reverse
#define Sign(a) ( ((a) > 0) - ((a) < 0) ) 

//  Scores used in the dynamic programming.  PEDMINSCORE wants to be INT32_MIN, except that it will
//  underflow, and wrap around to the highest score possible, on the first insert.  I think it
//  suffices to set this just slightly higher than the minimum, since it should never be picked as a
//  maximum, and thus never see more than one insert.  Just to be safe, we give it lots of space.
//
#define PEDMATCH       1
#define PEDMISMATCH   -2
#define PEDGAP        -3

#define PEDMINSCORE   (INT32_MIN - PEDGAP * 32768)

#define USE_SCORE

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
struct Match_Node_t {
  int32  Offset;              // To start of exact match in  hash-table frag
  int32  Len;                 // Of exact match
  int32  Start;               // Of exact match in current (new) frag
  int32  Next;                // Subscript of next match in list
};



struct pedEdit {
  int32  row;
  int32  score;
};



class prefixEditDistance {
public:
  prefixEditDistance(pedAlignType alignType_, double maxErate_);
  ~prefixEditDistance();


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


  pedOverlapType
  Extend_Alignment(Match_Node_t *Match,
                   char         *S,     int32   S_Len,
                   char          *T,     int32   T_Len,
                   int32         &S_Lo,  int32   &S_Hi,
                   int32         &T_Lo,  int32   &T_Hi,
                   int32         &Errors);

public:
  //  The four below were global #defines, two depended on the error rate which is now local.

  //  Most errors in any edit distance computation
  uint32                  MAX_ERRORS;

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

  int32                   Left_Delta_Len;
  int32                  *Left_Delta;

  int32                   Right_Delta_Len;
  int32                  *Right_Delta;

  int32                  *Delta_Stack;

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
};


#endif

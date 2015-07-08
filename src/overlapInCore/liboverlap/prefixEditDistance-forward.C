#include "prefixEditDistance.H"

#undef DEBUG



//  Put the delta encoding of the alignment represented in Edit_Array
//  starting at row e (which is the number of errors) and column d
//  (which is the diagonal) and working back to the start, into
//  Right_Delta. Set Right_Delta_Len to the number of
//  delta entries.
//
//  Used in Prefix_Edit_Distance
//
void
prefixEditDistance::Set_Right_Delta (int e, int d) {
  int  from, last, max;
  int  i, j, k;

  assert(Edit_Array_Lazy[e] != NULL);

  last = Edit_Array_Lazy[e][d];
  Right_Delta_Len = 0;

  for (k = e;  k > 0;  k--) {
    assert(Edit_Array_Lazy[k] != NULL);

    from = d;
    max = 1 + Edit_Array_Lazy[k - 1][d];
    if ((j = Edit_Array_Lazy[k - 1][d - 1]) > max) {
      from = d - 1;
      max = j;
    }
    if ((j = 1 + Edit_Array_Lazy[k - 1][d + 1]) > max) {
      from = d + 1;
      max = j;
    }
    if (from == d - 1) {
      Delta_Stack[Right_Delta_Len++] = max - last - 1;
      d--;
      last = Edit_Array_Lazy[k - 1][from];
    } else if (from == d + 1) {
      Delta_Stack[Right_Delta_Len++] = last - (max - 1);
      d++;
      last = Edit_Array_Lazy[k - 1][from];
    }
  }
  Delta_Stack[Right_Delta_Len++] = last + 1;

  k = 0;
  for (i = Right_Delta_Len - 1;  i > 0;  i--)
    Right_Delta[k++]
      = abs (Delta_Stack[i]) * Sign (Delta_Stack[i - 1]);
  Right_Delta_Len--;
}







//  Return the minimum number of changes (inserts, deletes, replacements)
//  needed to match string  A[0 .. (m-1)]  with a prefix of string
//   T[0 .. (n-1)]  if it's not more than  Error_Limit .
//  If no match, return the number of errors for the best match
//  up to a branch point.
//  Put delta description of alignment in  Right_Delta  and set
//  Right_Delta_Len  to the number of entries there if it's a complete
//  match.
//  Set  A_End  and  T_End  to the rightmost positions where the
//  alignment ended in  A  and  T , respectively.
//  Set  Match_To_End  true if the match extended to the end
//  of at least one string; otherwise, set it false to indicate
//  a branch point.

int32
prefixEditDistance::forward(char    *A,   int32 m,
                            char    *T,   int32 n,
                            int32    Error_Limit,
                            int32   &A_End,
                            int32   &T_End,
                            bool    &Match_To_End) {
  double  Score;
  int  Max_Score_Len = 0, Max_Score_Best_d = 0, Max_Score_Best_e = 0;
  int  Best_d, Best_e, From, Last, Longest, Max, Row;
  int  d, e, i, j, k;

  assert (m <= n);
  Best_d = Best_e = Longest = 0;
  Right_Delta_Len = 0;

  for (Row = 0;  Row < m
          && (A[Row] == T[Row]
              || A[Row] == 'N'
              || T[Row] == 'N');  Row++)
    ;

  if (Edit_Array_Lazy[0] == NULL)
    Allocate_More_Edit_Space();

  Edit_Array_Lazy[0][0] = Row;

  if (Row == m) {
    // Exact match
    A_End = T_End = m;
    Match_To_End = TRUE;
#ifdef DEBUG
    fprintf(stderr, "forward()- exact match\n");
#endif
    return  0;
  }

  int32 Left  = 0;
  int32 Right = 0;

  double Max_Score = 0.0;

  for (e = 1;  e <= Error_Limit;  e++) {
    if (Edit_Array_Lazy[e] == NULL)
      Allocate_More_Edit_Space();

    Left  = MAX (Left  - 1, -e);
    Right = MIN (Right + 1,  e);

    Edit_Array_Lazy[e - 1][Left     ] = -2;
    Edit_Array_Lazy[e - 1][Left  - 1] = -2;
    Edit_Array_Lazy[e - 1][Right    ] = -2;
    Edit_Array_Lazy[e - 1][Right + 1] = -2;

    for (d = Left;  d <= Right;  d++) {
      Row = 1 + Edit_Array_Lazy[e - 1][d];

      if ((j = Edit_Array_Lazy[e - 1][d - 1]) > Row)
        Row = j;

      if ((j = 1 + Edit_Array_Lazy[e - 1][d + 1]) > Row)
        Row = j;

      while  (Row < m && Row + d < n && (A[Row] == T[Row + d] || A[Row] == 'N' || T[Row + d] == 'N'))
        Row++;

      Edit_Array_Lazy[e][d] = Row;

      if (Row == m || Row + d == n) {
        //  Check for branch point here caused by uneven distribution of errors
        Score = Row * Branch_Match_Value - e;  //  Assumes Branch_Match_Value - Branch_Error_Value == 1.0

        int32  Tail_Len = Row - Max_Score_Len;
        bool   abort    = false;

        double slope    = (double)(Max_Score - Score) / Tail_Len;

#ifdef DEBUG 
        fprintf(stderr, "prefixEditDistance()-- e=%d MIN=%d Tail_Len=%d Max_Score=%d Score=%d slope=%f SLOPE=%f\n",
                e, MIN_BRANCH_END_DIST, Tail_Len, Max_Score, Score, slope, MIN_BRANCH_TAIL_SLOPE);
#endif

        //  If we're looking for local (former partial) overlaps, stop as soon as the score decreases.

        if ((alignType == pedLocal) &&
            (Score < Max_Score))
          abort = true;

        //  If we're looking for overlaps, use a more complicated rule that....does something.

        if ((alignType == pedOverlap) &&
            (e        >  MIN_BRANCH_END_DIST / 2) &&  //  e == number of errors in the current alignment?
            (Tail_Len >= MIN_BRANCH_END_DIST) &&      //  tail_len == amount extended since last best alignment
            (slope    >= MIN_BRANCH_TAIL_SLOPE))      //  slope == rate of score drop in the extenstion
          abort = true;

        //  If we're looking for global overlaps, don't stop.

        if ((alignType == pedGlobal))
          abort = false;


        if (abort) {
          A_End = Max_Score_Len;
          T_End = Max_Score_Len + Max_Score_Best_d;

          Set_Right_Delta (Max_Score_Best_e, Max_Score_Best_d);

          Match_To_End = FALSE;

#ifdef DEBUG
          fprintf(stderr, "forward()- ABORT alignment\n");
#endif
          return(Max_Score_Best_e);
        }

        // Force last error to be mismatch rather than insertion
        if ((Row == m) &&
            (1 + Edit_Array_Lazy[e - 1][d + 1] == Edit_Array_Lazy[e][d]) &&
            (d < Right)) {
          d++;
          Edit_Array_Lazy[e][d] = Edit_Array_Lazy[e][d - 1];
        }

        A_End = Row;           // One past last align position
        T_End = Row + d;

        Set_Right_Delta (e, d);

        Match_To_End = TRUE;

#ifdef DEBUG
        fprintf(stderr, "forward()- END alignment\n");
#endif
        return(e);
      }
    }

    while  ((Left <= Right) && (Left < 0) && (Edit_Array_Lazy[e][Left] < Edit_Match_Limit[e]))
      Left++;

    if (Left >= 0)
      while  ((Left <= Right) && (Edit_Array_Lazy[e][Left] + Left < Edit_Match_Limit[e]))
        Left++;

    if (Left > Right) {
#ifdef DEBUG
      fprintf(stderr, "reverse()- Left=%d Right=%d BREAK\n", Left, Right);
#endif
      break;
    }

    while  ((Right > 0) && (Edit_Array_Lazy[e][Right] + Right < Edit_Match_Limit[e]))
      Right--;

    if (Right <= 0)
      while  (Edit_Array_Lazy[e][Right] < Edit_Match_Limit[e])
        Right--;

    assert (Left <= Right);

    for (d = Left;  d <= Right;  d++)
      if (Edit_Array_Lazy[e][d] > Longest) {
        Best_d = d;
        Best_e = e;
        Longest = Edit_Array_Lazy[e][d];
      }

    Score = Longest * Branch_Match_Value - e;

    // Assumes  Branch_Match_Value - Branch_Error_Value == 1.0
    if (Score > Max_Score) {
      Max_Score = Score;
      Max_Score_Len = Longest;
      Max_Score_Best_d = Best_d;
      Max_Score_Best_e = Best_e;
    }
  }

#ifdef DEBUG
  fprintf(stderr, "forward()- return e=%d Error_Limit=%d\n",
          e, Error_Limit);
#endif

  A_End = Max_Score_Len;
  T_End = Max_Score_Len + Max_Score_Best_d;
  Set_Right_Delta (Max_Score_Best_e, Max_Score_Best_d);
  Match_To_End = FALSE;
  return  Max_Score_Best_e;
}


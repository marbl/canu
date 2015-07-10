#include "prefixEditDistance.H"

#undef DEBUG



//  Put the delta encoding of the alignment represented in Edit_Array
//  starting at row e (which is the number of errors) and column d
//  (which is the diagonal) and working back to the start, into
//  Left_Delta. Set Left_Delta_Len to the number of
//  delta entries and set leftover to the number of
//  characters that match after the last Left_Delta entry.
//  Don't allow the first delta to be an indel if it can be
//  converted to a substitution by adjusting t_end which
//  is where the alignment ends in the T string, which has length
//   t_len .
//
void
prefixEditDistance::Set_Left_Delta(int32 e, int32 d, int32 &leftover, int32 &t_end, int32 t_len) {

  Left_Delta_Len = 0;

  int32  lastr = Edit_Array_Lazy[e][d].row;
  int32  lasts = Edit_Array_Lazy[e][d].score;

  for (int32 k=e; k>0; k--) {
    assert(Edit_Array_Lazy[k] != NULL);

    int32  from = d;
    int32  maxs =     Edit_Array_Lazy[k-1][d].score + PEDMISMATCH;
    int32  maxr = 1 + Edit_Array_Lazy[k-1][d].row;

    //  Figure out which sequence has the insertion, if any.

#ifndef USE_SCORE
    if (0 + Edit_Array_Lazy[k-1][d-1].row > maxr) {
      from = d - 1;
      maxs = 0 + Edit_Array_Lazy[k-1][d-1].score;
      maxr = 0 + Edit_Array_Lazy[k-1][d-1].row;
    }

    if (1 + Edit_Array_Lazy[k-1][d+1].row > maxr) {
      from = d + 1;
      maxs = 1 + Edit_Array_Lazy[k-1][d+1].score;
      maxr = 1 + Edit_Array_Lazy[k-1][d+1].row;
    }
#else
    if (    Edit_Array_Lazy[k-1][d-1].score + PEDGAP > maxs) {
      from = d - 1;
      maxs =     Edit_Array_Lazy[k-1][d-1].score + PEDGAP;
      maxr = 0 + Edit_Array_Lazy[k-1][d-1].row;
    }

    if (    Edit_Array_Lazy[k-1][d+1].score + PEDGAP > maxs) {
      from = d + 1;
      maxs =     Edit_Array_Lazy[k-1][d+1].score + PEDGAP;
      maxr = 1 + Edit_Array_Lazy[k-1][d+1].row;
    }
#endif

    //  And make an insertion.  

    if (from == d - 1) {
      //fprintf(stderr, "LeftDelta:  insert gap in T at %d max=%d last=%d\n", maxr - lastr - 1, maxr, lastr);
      Left_Delta[Left_Delta_Len++] = maxr - lastr - 1;
      d--;
      lastr = Edit_Array_Lazy[k-1][from].row;
      lasts = Edit_Array_Lazy[k-1][from].score;
    }

    else if (from == d + 1) {
      //fprintf(stderr, "LeftDelta:  insert gap in A at %d max=%d last=%d\n", lastr - (maxr - 1), maxr, lastr);
      Left_Delta[Left_Delta_Len++] = lastr - (maxr - 1);
      d++;
      lastr = Edit_Array_Lazy[k-1][from].row;
      lasts = Edit_Array_Lazy[k-1][from].score;
    }

    else {
      //fprintf(stderr, "LeftDelta:  mismatch        at %d max=%d last=%d\n", maxr - lastr - 1, maxr, lastr);
    }
  }

  leftover = lastr;

  // Don't allow first delta to be +1 or -1
  assert((Left_Delta_Len == 0) || (Left_Delta[0] != -1));

  if ((Left_Delta_Len > 1) && (Left_Delta[0] == 1) && (t_end + t_len > 0)) {
    if (Left_Delta[1] > 0)
      Left_Delta[0] = Left_Delta[1] + 1;
    else
      Left_Delta[0] = Left_Delta[1] - 1;

    for (int32 i=2; i<Left_Delta_Len; i++)
      Left_Delta[i-1] = Left_Delta[i];

    Left_Delta_Len--;

    t_end--;

    if (Left_Delta_Len == 0)
      leftover++;
  }
}



//  Return the minimum number of changes (inserts, deletes, replacements)
//  needed to match string  A[0 .. (1-m)]  right-to-left with a prefix of string
//   T[0 .. (1-n)]  right-to-left if it's not more than  Error_Limit .
//  If no match, return the number of errors for the best match
//  up to a branch point.
//  Put delta description of alignment in  Left_Delta  and set
//  Left_Delta_Len  to the number of entries there.
//  Set  A_End  and  T_End  to the leftmost positions where the
//  alignment ended in  A  and  T , respectively.
//  If the alignment succeeds set  Leftover  to the number of
//  characters that match after the last  Left_Delta  entry;
//  otherwise, set  Leftover  to zero.
//  Set  Match_To_End  true if the match extended to the end
//  of at least one string; otherwise, set it false to indicate
//  a branch point.

int32
prefixEditDistance::reverse(char    *A,   int32 Alen,   //  first sequence and length
                            char    *T,   int32 Tlen,   //  second sequence and length
                            int32    Error_Limit,
                            int32   &A_End,
                            int32   &T_End,
                            int32   &Leftover,      //  <- novel
                            bool    &Match_To_End) {

  assert (Alen <= Tlen);

  int32  Best_d     = 0;
  int32  Best_e     = 0;
  int32  Best_row   = 0;
  int32  Best_score = 0;

  Left_Delta_Len = 0;

  int32  Row = 0;
  int32  Sco = 0;
  int32  Max = 0;

  //  Skip ahead over matches.  The original used to also skip if either sequence was N.
  while ((Row < Alen) && (A[-Row] == T[-Row])) {
    Row++;
    Sco += PEDMATCH;
  }

  if (Edit_Array_Lazy[0] == NULL)
    Allocate_More_Edit_Space();

  Edit_Array_Lazy[0][0].score = Sco;
  Edit_Array_Lazy[0][0].row   = Row;

  //  Exact match?

  if (Row == Alen) {
    A_End        = -Alen;
    T_End        = -Alen;
    Leftover     =  Alen;
    Match_To_End = true;
#ifdef DEBUG
    fprintf(stderr, "reverse()- exact match\n");
#endif
    return(0);
  }

  int32  Left  = 0;
  int32  Right = 0;

#ifndef USE_SCORE
  double Max_Score         = 0.0;
#else
  int32  Max_Score         = PEDMINSCORE;
#endif
  int32  Max_Score_Len     = 0;
  int32  Max_Score_Best_d  = 0;
  int32  Max_Score_Best_e  = 0;

  for (int32 e = 1;  e <= Error_Limit;  e++) {
    if (Edit_Array_Lazy[e] == NULL)
      Allocate_More_Edit_Space();

    Left  = MAX (Left  - 1, -e);
    Right = MIN (Right + 1,  e);

    Edit_Array_Lazy[e-1][Left  - 1].row = -2;    Edit_Array_Lazy[e-1][Left  - 1].score = PEDMINSCORE;
    Edit_Array_Lazy[e-1][Left     ].row = -2;    Edit_Array_Lazy[e-1][Left     ].score = PEDMINSCORE;
    //  Of note, [0][0] on the first iteration is not reset here.
    Edit_Array_Lazy[e-1][Right    ].row = -2;    Edit_Array_Lazy[e-1][Right    ].score = PEDMINSCORE;
    Edit_Array_Lazy[e-1][Right + 1].row = -2;    Edit_Array_Lazy[e-1][Right + 1].score = PEDMINSCORE;

    for (int32 d = Left;  d <= Right;  d++) {

#ifndef USE_SCORE
      //  A mismatch.
      Row = 1 + Edit_Array_Lazy[e-1][d].row;
      Sco =     Edit_Array_Lazy[e-1][d].score + PEDMISMATCH;

      //  Insert a gap in T.  (opposite of the 'forward')
      if (0 + Edit_Array_Lazy[e-1][d-1].row > Row) {
        Row = 0 + Edit_Array_Lazy[e-1][d-1].row;
        Sco =     Edit_Array_Lazy[e-1][d-1].score + PEDGAP;
      }

      //  Insert a gap in A.  (opposite of the 'forward')
      if (1 + Edit_Array_Lazy[e-1][d+1].row > Row) {
        Row = 1 + Edit_Array_Lazy[e-1][d+1].row;
        Sco =     Edit_Array_Lazy[e-1][d+1].score + PEDGAP;
      }
#else
      //  A mismatch.
      Row = 1 + Edit_Array_Lazy[e-1][d].row;
      Sco =     Edit_Array_Lazy[e-1][d].score + PEDMISMATCH;

      //  Insert a gap in T.  (opposite of the 'forward')
      if (    Edit_Array_Lazy[e-1][d-1].score + PEDGAP > Sco) {
        Row = 0 + Edit_Array_Lazy[e-1][d-1].row;
        Sco =     Edit_Array_Lazy[e-1][d-1].score + PEDGAP;
      }

      //  Insert a gap in A.  (opposite of the 'forward')
      if (    Edit_Array_Lazy[e-1][d+1].score + PEDGAP > Sco) {
        Row = 1 + Edit_Array_Lazy[e-1][d+1].row;
        Sco =     Edit_Array_Lazy[e-1][d+1].score + PEDGAP;
      }
#endif

      //  If A or B is N, that isn't a mismatch.
      //  If A is lowercase and T is uppercase, it's a match.
      //  If A is lowercase and T doesn't match, ignore the cost of the gap in B

      while ((Row < Alen) && (Row + d < Tlen) && (A[-Row] == T[-Row - d])) {
        Row++;
        Sco += PEDMATCH;
      }

      Edit_Array_Lazy[e][d].row   = Row;
      Edit_Array_Lazy[e][d].score = Sco;


      if (Row == Alen || Row + d == Tlen) {

        //  Check for branch point here caused by uneven distribution of errors
        //double Score = Row * Branch_Match_Value - e;
        int32  Score = Sco;

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
          A_End = - Max_Score_Len;
          T_End = - Max_Score_Len - Max_Score_Best_d;

          Set_Left_Delta (Max_Score_Best_e, Max_Score_Best_d, Leftover, T_End, Tlen);

          Match_To_End = FALSE;

#ifdef DEBUG
          fprintf(stderr, "reverse()- ABORT alignment\n");
#endif
          return(Max_Score_Best_e);
        }

        A_End = - Row;           // One past last align position
        T_End = - Row - d;

        Set_Left_Delta (e, d, Leftover, T_End, Tlen);

        Match_To_End = TRUE;

#ifdef DEBUG
        fprintf(stderr, "reverse()- END alignment\n");
#endif
        return(e);
      }
    }

    //  Reset the band

    while  ((Left <= Right) && (Left < 0) && (Edit_Array_Lazy[e][Left].row < Edit_Match_Limit[e]))
      Left++;

    if (Left >= 0)
      while  ((Left <= Right) && (Edit_Array_Lazy[e][Left].row + Left < Edit_Match_Limit[e]))
        Left++;

    if (Left > Right) {
#ifdef DEBUG
      fprintf(stderr, "reverse()- Left=%d Right=%d BREAK\n", Left, Right);
#endif
      break;
    }

    while  ((Right > 0) && (Edit_Array_Lazy[e][Right].row + Right < Edit_Match_Limit[e]))
      Right--;

    if (Right <= 0)
      while  (Edit_Array_Lazy[e][Right].row < Edit_Match_Limit[e])
        Right--;

    assert (Left <= Right);

#ifndef USE_SCORE
    for (int32 d = Left;  d <= Right;  d++)
      if (Edit_Array_Lazy[e][d].row > Best_row) {
        Best_d     = d;
        Best_e     = e;
        Best_row   = Edit_Array_Lazy[e][d].row;
        Best_score = Edit_Array_Lazy[e][d].score;
      }

    if (Best_row * Branch_Match_Value - e > Max_Score) {
      Max_Score_Best_d = Best_d;
      Max_Score_Best_e = Best_e;
      Max_Score        = Best_row * Branch_Match_Value - e;
      Max_Score_Len    = Best_row;
    }
#else
    for (int32 d = Left;  d <= Right;  d++)
      if (Edit_Array_Lazy[e][d].score > Best_score) {
        Best_d     = d;
        Best_e     = e;
        Best_row   = Edit_Array_Lazy[e][d].row;
        Best_score = Edit_Array_Lazy[e][d].score;
      }

    if (Best_score > Max_Score) {
      Max_Score_Best_d = Best_d;
      Max_Score_Best_e = Best_e;
      Max_Score        = Best_score;
      Max_Score_Len    = Best_row;
    }
#endif
  }

#ifdef DEBUG  
  fprintf(stderr, "reverse()- return e=%d Error_Limit=%d\n",
          e, Error_Limit);
#endif

  A_End = - Max_Score_Len;
  T_End = - Max_Score_Len - Max_Score_Best_d;

  Set_Left_Delta (Max_Score_Best_e, Max_Score_Best_d, Leftover, T_End, Tlen);

  Match_To_End = FALSE;

  return  Max_Score_Best_e;
}

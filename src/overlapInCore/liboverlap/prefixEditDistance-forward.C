#include "prefixEditDistance.H"

#undef DEBUG



//  Put the delta encoding of the alignment represented in Edit_Array
//  starting at row e (which is the number of errors) and column d
//  (which is the diagonal) and working back to the start, into
//  Right_Delta. Set Right_Delta_Len to the number of
//  delta entries.
//
void
prefixEditDistance::Set_Right_Delta(char  *A, char  *T,
                                    int32  e, int32  d) {

  Right_Delta_Len = 0;

  int32  lastr = Edit_Array_Lazy[e][d].row;
  int32  fromd = Edit_Array_Lazy[e][d].fromd;  //   The cell that we're moving to in [e-1].

#if 0
  fprintf(stderr, "prefixEditDistance::Set_Right_Delta()-- e=%11d d=%11d lastr=%11d - r=%11d s=%11d d=%11d - r=%11d s=%11d d=%11d - r=%11d s=%11d d=%11d\n",
          e, d, lastr,
          Edit_Array_Lazy[e][d-1].row, Edit_Array_Lazy[e][d-1].score, Edit_Array_Lazy[e][d-1].fromd,
          Edit_Array_Lazy[e][d+0].row, Edit_Array_Lazy[e][d+0].score, Edit_Array_Lazy[e][d+0].fromd,
          Edit_Array_Lazy[e][d+1].row, Edit_Array_Lazy[e][d+1].score, Edit_Array_Lazy[e][d+1].fromd);
#endif

  for (int32 k=e; k>0; k--) {
    assert(Edit_Array_Lazy[k] != NULL);

    //  Analyze cells at errors = k-1 for the maximum -- no analysis needed, since we stored this cell as fromd.

#if 0
    fprintf(stderr, "prefixEditDistance::Set_Right_Delta()-- e=%11d d=%11d lastr=%11d - r=%11d s=%11d d=%11d - r=%11d s=%11d d=%11d - r=%11d s=%11d d=%11d\n",
            k-1, d, lastr,
            Edit_Array_Lazy[k-1][d-1].row, Edit_Array_Lazy[k-1][d-1].score, Edit_Array_Lazy[k-1][d-1].fromd,
            Edit_Array_Lazy[k-1][d+0].row, Edit_Array_Lazy[k-1][d+0].score, Edit_Array_Lazy[k-1][d+0].fromd,
            Edit_Array_Lazy[k-1][d+1].row, Edit_Array_Lazy[k-1][d+1].score, Edit_Array_Lazy[k-1][d+1].fromd);
#endif


#if 1
    int32   from = Edit_Array_Lazy[k][d].fromd;

    if (from == d - 1) {
      Delta_Stack[Right_Delta_Len++] = Edit_Array_Lazy[k-1][d-1].row - lastr - 1;
      d--;
      lastr = Edit_Array_Lazy[k-1][from].row;
    }

    else if (from == d + 1) {
      Delta_Stack[Right_Delta_Len++] = lastr - Edit_Array_Lazy[k-1][d+1].row;
      d++;
      lastr = Edit_Array_Lazy[k-1][from].row;
    }

    else {
      //fprintf(stderr, "LeftDelta:  mismatch        at %d max=%d last=%d\n", maxr - lastr - 1, maxr, lastr);
    }
#endif


#if 0
    //  Original
    int32  from = d;
    int32  maxs =     Edit_Array_Lazy[k-1][d].score;// + PEDMISMATCH;
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
    if (    Edit_Array_Lazy[k-1][d-1].score /*+ PEDGAP*/ > maxs) {
      from = d - 1;
      maxs =     Edit_Array_Lazy[k-1][d-1].score;// + PEDGAP;
      maxr = 0 + Edit_Array_Lazy[k-1][d-1].row;
    }

    if (    Edit_Array_Lazy[k-1][d+1].score /*+ PEDGAP*/ > maxs) {
      from = d + 1;
      maxs =     Edit_Array_Lazy[k-1][d+1].score;// + PEDGAP;
      maxr = 1 + Edit_Array_Lazy[k-1][d+1].row;
    }
#endif

    //  And make an insertion.  

    if (from == d - 1) {
      //fprintf(stderr, "RightDelta: insert gap in A at %d max=%d last=%d\n", maxr - lastr - 1, maxr, lastr);
      Delta_Stack[Right_Delta_Len++] = maxr - lastr - 1;
      d--;
      lastr = Edit_Array_Lazy[k-1][from].row;
    }

    else if (from == d + 1) {
      //fprintf(stderr, "RightDelta: insert gap in T at %d max=%d last=%d\n", lastr - (maxr - 1), maxr, lastr);
      Delta_Stack[Right_Delta_Len++] = lastr - (maxr - 1);
      d++;
      lastr = Edit_Array_Lazy[k-1][from].row;
    }

    else {
      //fprintf(stderr, "RightDelta: mismatch        at %d max=%d last=%d\n", maxr - lastr - 1, maxr, lastr);
    }
#endif

  }

  Delta_Stack[Right_Delta_Len++] = lastr + 1;

  for (int32 k=0, i=Right_Delta_Len-1; i>0; i--)
    Right_Delta[k++] = abs(Delta_Stack[i]) * Sign(Delta_Stack[i-1]);

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
prefixEditDistance::forward(char    *A,   int32 Alen,
                            char    *T,   int32 Tlen,
                            int32    Error_Limit,
                            int32   &A_End,
                            int32   &T_End,
                            bool    &Match_To_End) {

  assert (Alen <= Tlen);

  int32  Best_d      = 0;
  int32  Best_e      = 0;
  int32  Best_row    = 0;
  int32  Best_score  = 0;

  Right_Delta_Len = 0;

  int32  Row = 0;
  int32  Dst = 0;
  int32  Err = 0;
  int32  Sco = 0;

  int32  fromd = 0;

  //  Skip ahead over matches.  The original used to also skip if either sequence was N.
  while ((Row < Alen) && (isMatch(A[Row], T[Row]))) {
    Sco += matchScore(A[Row], T[Row]);  //PEDMATCH;
    Row++;
  }

  if (Edit_Array_Lazy[0] == NULL)
    Allocate_More_Edit_Space();

  Edit_Array_Lazy[0][0].score  = Sco;
  Edit_Array_Lazy[0][0].dist   = Dst;
  Edit_Array_Lazy[0][0].errs   = 0;
  Edit_Array_Lazy[0][0].row    = Row;

  // Exact match?

  if (Row == Alen) {
    A_End        = Alen;
    T_End        = Alen;
    Match_To_End = true;
#ifdef DEBUG
    fprintf(stderr, "prefixEditDistance::forward()- exact match\n");
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

  for (int32 ei=1; ei <= Edit_Space_Max; ei++) {
    if (Edit_Array_Lazy[ei] == NULL)
      Allocate_More_Edit_Space();

    Left  = MAX (Left  - 1, -ei);
    Right = MIN (Right + 1,  ei);

    Edit_Array_Lazy[ei-1][Left  - 1].init();
    Edit_Array_Lazy[ei-1][Left     ].init();
    //  Of note, [0][0] on the first iteration is not reset here.
    Edit_Array_Lazy[ei-1][Right    ].init();
    Edit_Array_Lazy[ei-1][Right + 1].init();

    for (int32 d = Left;  d <= Right;  d++) {

#ifndef USE_SCORE
      //  A mismatch.
      Row   = 1 + Edit_Array_Lazy[ei-1][d].row;
      Dst   =     Edit_Array_Lazy[ei-1][d].dist  + 1;
      Err   =     Edit_Array_Lazy[ei-1][d].errs  + 1;
      Sco   =     Edit_Array_Lazy[ei-1][d].score + PEDMISMATCH;
      fromd =     d;

      //  Insert a gap in A.
      if (0 + Edit_Array_Lazy[ei-1][d-1].row > Row) {
        Row   = 0 + Edit_Array_Lazy[ei-1][d-1].row;  //  +0 because row is the index into A, and the A base doesn't change.
        Dst   =     Edit_Array_Lazy[ei-1][d-1].dist  + 0;
        Err   =     Edit_Array_Lazy[ei-1][d-1].errs  + 0;
        Sco   =     Edit_Array_Lazy[ei-1][d-1].score + PEDGAP;
        fromd =     d-1;
      }

      //  Insert a gap in T.
      if (1 + Edit_Array_Lazy[ei-1][d+1].row > Row) {
        Row   = 1 + Edit_Array_Lazy[ei-1][d+1].row;  //  +1 because we ate up a base in A.
        Dst   =     Edit_Array_Lazy[ei-1][d+1].dist  + 1;
        Err   =     Edit_Array_Lazy[ei-1][d+1].errs  + 1;
        Sco   =     Edit_Array_Lazy[ei-1][d+1].score + PEDGAP;
        fromd =     d+1;
      }
#else
      //  A mismatch.
      Row   = 1 + Edit_Array_Lazy[ei-1][d].row;
      Dst   =     Edit_Array_Lazy[ei-1][d].dist  + 1;
      Err   =     Edit_Array_Lazy[ei-1][d].errs  + 1;
      Sco   =     Edit_Array_Lazy[ei-1][d].score + PEDMISMATCH;
      fromd =     d;

      //  Insert a gap in A.  Check the other sequence to see if this is a zero-cost gap.  Note
      //  agreement with future value of Row and what is used in isMatch() below.
      {
        int32  tPos     = 0 + Edit_Array_Lazy[ei-1][d-1].row + d;
        assert(tPos <= Tlen);

        int32  gapCost = isFreeGap( T[tPos] ) ? 0 : PEDGAP;

        if (Edit_Array_Lazy[ei-1][d-1].score + gapCost > Sco) {
          Row   =     Edit_Array_Lazy[ei-1][d-1].row;
          Dst   =     Edit_Array_Lazy[ei-1][d-1].dist  + (gapCost == 0) ? 0 : 0;
          Err   =     Edit_Array_Lazy[ei-1][d-1].errs  + (gapCost == 0) ? 0 : 0;
          Sco   =     Edit_Array_Lazy[ei-1][d-1].score +  gapCost;
          fromd =     d-1;
        }
      }

      //  Insert a gap in T.
      {
        int32  aPos = 1 + Edit_Array_Lazy[ei-1][d+1].row;
        assert(aPos <= Alen);

        int32  gapCost = isFreeGap( A[aPos] ) ? 0 : PEDGAP;

        if (Edit_Array_Lazy[ei-1][d+1].score + gapCost > Sco) {
          Row   = 1 + Edit_Array_Lazy[ei-1][d+1].row;
          Dst   =     Edit_Array_Lazy[ei-1][d+1].dist  + (gapCost == 0) ? 0 : 1;
          Err   =     Edit_Array_Lazy[ei-1][d+1].errs  + (gapCost == 0) ? 0 : 1;
          Sco   =     Edit_Array_Lazy[ei-1][d+1].score +  gapCost;
          fromd =     d+1;
        }
      }
#endif

      //  If A or B is N, that isn't a mismatch.
      //  If A is lowercase and T is uppercase, it's a match.
      //  If A is lowercase and T doesn't match, ignore the cost of the gap in B

      while ((Row < Alen) && (Row + d < Tlen) && (isMatch(A[Row], T[Row + d]))) {
        Sco += matchScore(A[Row], T[Row + d]);  //PEDMATCH;
        Row += 1;
        Dst += 1;
        Err += 0;
      }

      Edit_Array_Lazy[ei][d].row   = Row;
      Edit_Array_Lazy[ei][d].dist  = Dst;
      Edit_Array_Lazy[ei][d].errs  = Err;
      Edit_Array_Lazy[ei][d].score = Sco;
      Edit_Array_Lazy[ei][d].fromd = fromd;


      if (Row == Alen || Row + d == Tlen) {

        //  Check for branch point here caused by uneven distribution of errors
        //double Score = Row * Branch_Match_Value - e;
        int32  Score = Sco;

        int32  Tail_Len = Row - Max_Score_Len;
        bool   abort    = false;

        double slope    = (double)(Max_Score - Score) / Tail_Len;

#ifdef DEBUG 
        fprintf(stderr, "prefixEditDistance::forward()-- e=%d MIN=%d Tail_Len=%d Max_Score=%d Score=%d slope=%f SLOPE=%f\n",
                ei, MIN_BRANCH_END_DIST, Tail_Len, Max_Score, Score, slope, MIN_BRANCH_TAIL_SLOPE);
#endif

        //  If we're looking for local (former partial) overlaps, stop as soon as the score decreases.

        if ((alignType == pedLocal) &&
            (Score < Max_Score))
          abort = true;

        //  If we're looking for overlaps, use a more complicated rule that....does something.

        if ((alignType == pedOverlap) &&
            (ei       >  MIN_BRANCH_END_DIST / 2) &&  //  e == number of errors in the current alignment?
            (Tail_Len >= MIN_BRANCH_END_DIST) &&      //  tail_len == amount extended since last best alignment
            (slope    >= MIN_BRANCH_TAIL_SLOPE))      //  slope == rate of score drop in the extenstion
          abort = true;

        //  If we're looking for global overlaps, don't stop.

        if ((alignType == pedGlobal))
          abort = false;


        if (abort) {
          A_End = Max_Score_Len;
          T_End = Max_Score_Len + Max_Score_Best_d;

          Set_Right_Delta(A, T, Max_Score_Best_e, Max_Score_Best_d);

          Match_To_End = FALSE;

#ifdef DEBUG
          fprintf(stderr, "prefixEditDistance::forward()- ABORT alignment\n");
#endif
          return(Max_Score_Best_e);
        }

        // Force last error to be mismatch rather than insertion
        if ((Row == Alen) &&
            (1 + Edit_Array_Lazy[ei-1][d+1].row == Edit_Array_Lazy[ei][d].row) &&
            (d < Right)) {
          d++;
          Edit_Array_Lazy[ei][d].score  = Edit_Array_Lazy[ei][d-1].score;   //  ??
          Edit_Array_Lazy[ei][d].dist   = Edit_Array_Lazy[ei][d-1].dist;    //  ??
          Edit_Array_Lazy[ei][d].row    = Edit_Array_Lazy[ei][d-1].row;
        }

        A_End = Row;           // One past last align position
        T_End = Row + d;

        Set_Right_Delta(A, T, ei, d);

        Match_To_End = TRUE;

#ifdef DEBUG
        fprintf(stderr, "prefixEditDistance::forward()- END alignment\n");
#endif
        return(ei);
      }
    }  //  Over all diagonals.

    //  Reset the band
    //
    //  The .dist used to be .row.

#ifdef DEBUG
    //fprintf(stderr, "prefixEditDistance::forward()- Left=%d Right=%d e=%d row=%d limit=%d\n",
    //        Left, Right, Edit_Array_Lazy[ei][Left].errs, Edit_Array_Lazy[ei][Left].dist, Edit_Match_Limit[ Edit_Array_Lazy[ei][Left].errs ]);
#endif

    while  ((Left <= Right) && (Left < 0) && (Edit_Array_Lazy[ei][Left].dist < Edit_Match_Limit[ Edit_Array_Lazy[ei][Left].errs ]))
      Left++;

    if (Left >= 0)
      while  ((Left <= Right) && (Edit_Array_Lazy[ei][Left].dist + Left < Edit_Match_Limit[ Edit_Array_Lazy[ei][Left].errs ]))
        Left++;

    if (Left > Right)
      break;

#ifdef DEBUG
    //fprintf(stderr, "prefixEditDistance::forward()- Left=%d Right=%d e=%d row=%d limit=%d\n",
    //        Left, Right, Edit_Array_Lazy[ei][Left].errs, Edit_Array_Lazy[ei][Left].dist, Edit_Match_Limit[ Edit_Array_Lazy[ei][Left].errs ]);
#endif

    while  ((Right > 0) && (Edit_Array_Lazy[ei][Right].dist + Right < Edit_Match_Limit[ Edit_Array_Lazy[ei][Right].errs ]))
      Right--;

    if (Right <= 0)
      while  (Edit_Array_Lazy[ei][Right].dist < Edit_Match_Limit[ Edit_Array_Lazy[ei][Right].errs ])
        Right--;

    assert (Left <= Right);

#ifndef USE_SCORE
    for (int32 d = Left;  d <= Right;  d++)
      if (Edit_Array_Lazy[ei][d].row > Best_row) {
        Best_d      = d;
        Best_e      = ei;
        Best_row    = Edit_Array_Lazy[ei][d].row;
        Best_score  = Edit_Array_Lazy[ei][d].score;
      }

    if (Best_row * Branch_Match_Value - e > Max_Score) {
      Max_Score_Best_d = Best_d;
      Max_Score_Best_e = Best_e;
      Max_Score        = Best_row * Branch_Match_Value - e;
      Max_Score_Len    = Best_row;
    }
#else
    for (int32 d = Left;  d <= Right;  d++)
      if (Edit_Array_Lazy[ei][d].score > Best_score) {
        Best_d      = d;
        Best_e      = ei;
        Best_row    = Edit_Array_Lazy[ei][d].row;
        Best_score  = Edit_Array_Lazy[ei][d].score;
      }

    if (Best_score > Max_Score) {
      Max_Score_Best_d = Best_d;
      Max_Score_Best_e = Best_e;
      Max_Score        = Best_score;
      Max_Score_Len    = Best_row;
    }
#endif

  }  //  Over all possible number of errors

#ifdef DEBUG
  fprintf(stderr, "prefixEditDistance::forward()- iterated over all errors, return best found\n");
#endif

  A_End = Max_Score_Len;
  T_End = Max_Score_Len + Max_Score_Best_d;

  Set_Right_Delta(A, T, Max_Score_Best_e, Max_Score_Best_d);

  Match_To_End = false;

  return(Max_Score_Best_e);
}

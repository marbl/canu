
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

#include "NDalgorithm.H"



//  Put the delta encoding of the alignment represented in Edit_Array
//  starting at row e (which is the number of errors) and column d
//  (which is the diagonal) and working back to the start, into
//  Right_Delta. Set Right_Delta_Len to the number of
//  delta entries.
//
void
NDalgorithm::Set_Right_Delta(int32  e, int32  d) {

  Right_Score       = Edit_Array_Lazy[e][d].score;
  Right_Delta_Len   = 0;

  int32  lastr = Edit_Array_Lazy[e][d].row;

  //fprintf(stderr, "NDalgorithm::Set_Right_Delta()-- e  =%5d d=%5d lastr=%5d\n",
  //        e, d, lastr);

  for (int32 k=e; k>0; k--) {
    assert(Edit_Array_Lazy[k] != NULL);

    //  Analyze cells at errors = k-1 for the maximum -- no analysis needed, since we stored this cell as fromd.

    int32   from  = Edit_Array_Lazy[k][d].fromd;
    int32   lasts = Edit_Array_Lazy[k][d].score;

    //Edit_Array_Lazy[k-1][from].display(k-1, from);

    //fprintf(stderr, "NDalgorithm::Set_Right_Delta()-- k-1=%5d d=%5d from=%5d lastr=%5d - r=%5d s=%5d d=%5d\n",
    //        k-1, d, from, lastr,
    //        Edit_Array_Lazy[k-1][from].row, Edit_Array_Lazy[k-1][from].score, Edit_Array_Lazy[k-1][from].fromd);

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

void
NDalgorithm::forward(char    *A,   int32 Alen,
                     char    *T,   int32 Tlen,
                     int32   &A_End,
                     int32   &T_End,
                     bool    &Match_To_End) {

  assert (Alen <= Tlen);

  int32  Best_d      = 0;
  int32  Best_e      = 0;
  int32  Best_row    = 0;
  int32  Best_score  = 0;


  int32  Row = 0;
  int32  Dst = 0;
  int32  Err = 0;
  int32  Sco = 0;

  int32  fromd = 0;

  //  Skip ahead over matches.  The original used to also skip if either sequence was N.
  while ((Row < Alen) && (isMatch(A[Row], T[Row]))) {
    Sco += matchScore(A[Row], T[Row]);
    Row++;
  }

  if (Edit_Array_Lazy[0] == NULL)
    allocateMoreEditSpace();

  Edit_Array_Lazy[0][0].row    = Row;
  Edit_Array_Lazy[0][0].dist   = Dst;
  Edit_Array_Lazy[0][0].errs   = 0;
  Edit_Array_Lazy[0][0].score  = Sco;
  Edit_Array_Lazy[0][0].fromd  = INT32_MAX;

  // Exact match?

  if (Row == Alen) {
    A_End        = Alen;
    T_End        = Alen;
    Match_To_End = true;

    Right_Score       = Sco;
    Right_Delta_Len   = 0;

    return;
  }

  int32  Left  = 0;
  int32  Right = 0;

  int32  Max_Score         = PEDMINSCORE;
  int32  Max_Score_Len     = 0;
  int32  Max_Score_Best_d  = 0;
  int32  Max_Score_Best_e  = 0;

  for (int32 ei=1; ei <= Edit_Space_Max; ei++) {
    if (Edit_Array_Lazy[ei] == NULL)
      if (allocateMoreEditSpace() == false) {
        //  FAIL
        return;
      }

    Left  = std::max(Left  - 1, -ei);
    Right = std::min(Right + 1,  ei);

    //fprintf(stderr, "FORWARD ei=%d Left=%d Right=%d\n", ei, Left, Right);

    Edit_Array_Lazy[ei-1][Left  - 1].init();
    Edit_Array_Lazy[ei-1][Left     ].init();
    //  Of note, [0][0] on the first iteration is not reset here.
    Edit_Array_Lazy[ei-1][Right    ].init();
    Edit_Array_Lazy[ei-1][Right + 1].init();

    for (int32 d = Left;  d <= Right;  d++) {

      //  A mismatch.
      {
        int32  aPos         =  (1 + Edit_Array_Lazy[ei-1][d].row)     - 1;  //  -1 because we need to compare the base we are at,
        int32  tPos         =  (1 + Edit_Array_Lazy[ei-1][d].row) + d - 1;  //  not the base we will be at after the mismatch

        Row   = 1 + Edit_Array_Lazy[ei-1][d].row;
        Dst   =     Edit_Array_Lazy[ei-1][d].dist  + 1;
        Err   =     Edit_Array_Lazy[ei-1][d].errs  + 1;
        fromd =     d;

        //  If positive, we have a pointer into valid sequence.  If not, this mismatch
        //  doesn't make sense, and the row/score are set to bogus values.

        if ((aPos >= 0) && (tPos >= 0)) {
          assert (aPos <= Alen);
          assert( tPos <= Tlen);

          assert(A[aPos] != T[tPos]);

          Sco = Edit_Array_Lazy[ei-1][d].score + mismatchScore(A[aPos], T[tPos]);

        } else {
          Sco = PEDMINSCORE;
        }
      }

      //  Insert a gap in A.  Check the other sequence to see if this is a zero-cost gap.  Note
      //  agreement with future value of Row and what is used in isMatch() below.

      {
        int32  tPos    = 0 + Edit_Array_Lazy[ei-1][d-1].row + d;

        //assert(tPos >= 0);
        //assert(tPos < Tlen);

        if ((tPos >= 0) && (tPos <= Tlen)) {
          int32  gapCost = isFreeGap( T[tPos] ) ? PEDFREEGAP : PEDGAP;

          //if (gapCost == 0)
          //  fprintf(stderr, "NDalgorithm::forward()--  free A gap for aPos=%d tPos=%d t=%c/%d\n", tPos - d, tPos, T[tPos], T[tPos]);

          if (Edit_Array_Lazy[ei-1][d-1].score + gapCost > Sco) {
            Row   =     Edit_Array_Lazy[ei-1][d-1].row;
            Dst   =     Edit_Array_Lazy[ei-1][d-1].dist  + (gapCost == PEDFREEGAP) ? 0 : 0;
            Err   =     Edit_Array_Lazy[ei-1][d-1].errs  + (gapCost == PEDFREEGAP) ? 0 : 0;
            Sco   =     Edit_Array_Lazy[ei-1][d-1].score +  gapCost;
            fromd =     d-1;
          }
        }
      }

      //  Insert a gap in T.
      //  Testcase test-st-ts shows this works.

      {
        int32  aPos    = 1 + Edit_Array_Lazy[ei-1][d+1].row;

        //assert(aPos >= 0);
        //assert(aPos < Tlen);

        if ((aPos >= 0) && (aPos <= Alen)) {
          int32  gapCost = isFreeGap( A[aPos] ) ? 0 : PEDGAP;

          //if (gapCost == 0)
          //  fprintf(stderr, "NDalgorithm::forward()--  free T gap for aPos=%d tPos=%d a=%c/%d\n", aPos, aPos + d, A[aPos], A[aPos]);

          if (Edit_Array_Lazy[ei-1][d+1].score + gapCost > Sco) {
            Row   = 1 + Edit_Array_Lazy[ei-1][d+1].row;
            Dst   =     Edit_Array_Lazy[ei-1][d+1].dist  + (gapCost == PEDFREEGAP) ? 0 : 1;
            Err   =     Edit_Array_Lazy[ei-1][d+1].errs  + (gapCost == PEDFREEGAP) ? 0 : 1;
            Sco   =     Edit_Array_Lazy[ei-1][d+1].score +  gapCost;
            fromd =     d+1;
          }
        }
      }

      //  If A or B is N, that isn't a mismatch.
      //  If A is lowercase and T is uppercase, it's a match.
      //  If A is lowercase and T doesn't match, ignore the cost of the gap in B

      while ((Row < Alen) && (Row + d < Tlen) && (isMatch(A[Row], T[Row + d]))) {
        Sco += matchScore(A[Row], T[Row + d]);
        Row += 1;
        Dst += 1;
        Err += 0;
      }

      Edit_Array_Lazy[ei][d].row   = Row;
      Edit_Array_Lazy[ei][d].dist  = Dst;
      Edit_Array_Lazy[ei][d].errs  = Err;
      Edit_Array_Lazy[ei][d].score = Sco;
      Edit_Array_Lazy[ei][d].fromd = fromd;

      //fprintf(stderr, "SET ei=%d d=%d -- row=%d dist=%d errs=%d score=%d fromd=%d\n", ei, d, Row, Dst, Err, Sco, fromd);

      if (Row == Alen || Row + d == Tlen) {
        A_End = Row;           // One past last align position
        T_End = Row + d;

        Set_Right_Delta(ei, d);

        Match_To_End = true;

        return;  //return(ei);
      }
    }  //  Over all diagonals.

    //  Reset the band
    //
    //  The .dist used to be .row.

    while  ((Left <= Right) && (Left < 0) && (Edit_Array_Lazy[ei][Left].dist < Edit_Match_Limit[ Edit_Array_Lazy[ei][Left].errs ]))
      Left++;

    if (Left >= 0)
      while  ((Left <= Right) && (Edit_Array_Lazy[ei][Left].dist + Left < Edit_Match_Limit[ Edit_Array_Lazy[ei][Left].errs ]))
        Left++;

    if (Left > Right)
      break;

    while  ((Right > 0) && (Edit_Array_Lazy[ei][Right].dist + Right < Edit_Match_Limit[ Edit_Array_Lazy[ei][Right].errs ]))
      Right--;

    if (Right <= 0)
      while  (Edit_Array_Lazy[ei][Right].dist < Edit_Match_Limit[ Edit_Array_Lazy[ei][Right].errs ])
        Right--;

    assert (Left <= Right);

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
  }  //  Over all possible number of errors

  //fprintf(stderr, "NDalgorithm::forward()- iterated over all errors, return best found\n");

  A_End = Max_Score_Len;
  T_End = Max_Score_Len + Max_Score_Best_d;

  Set_Right_Delta(Max_Score_Best_e, Max_Score_Best_d);

  Match_To_End = false;

  return;  //return(Max_Score_Best_e);
}

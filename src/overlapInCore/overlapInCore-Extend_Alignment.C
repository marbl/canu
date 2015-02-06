
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static const char *rcsid = "$Id$";

#include "overlapInCore.H"




//  Put the delta encoding of the alignment represented in WA->Edit_Array
//  starting at row e (which is the number of errors) and column d
//  (which is the diagonal) and working back to the start, into
//  WA->Left_Delta. Set WA->Left_Delta_Len to the number of
//  delta entries and set (* leftover) to the number of
//  characters that match after the last WA->Left_Delta entry.
//  Don't allow the first delta to be an indel if it can be
//  converted to a substitution by adjusting (* t_end) which
//  is where the alignment ends in the T string, which has length
//   t_len .

static
void
Set_Left_Delta (int e, int d, int * leftover, int * t_end, int t_len, Work_Area_t * WA) {
  int  from, last, max;
  int  j, k;

  assert(WA->Edit_Array_Lazy[e] != NULL);

  last = WA->Edit_Array_Lazy[e][d];
  WA->Left_Delta_Len = 0;

  for  (k = e;  k > 0;  k--) {
    assert(WA->Edit_Array_Lazy[k] != NULL);

    from = d;
    max = 1 + WA->Edit_Array_Lazy[k - 1][d];
    if  ((j = WA->Edit_Array_Lazy[k - 1][d - 1]) > max) {
      from = d - 1;
      max = j;
    }
    if  ((j = 1 + WA->Edit_Array_Lazy[k - 1][d + 1]) > max) {
      from = d + 1;
      max = j;
    }
    if  (from == d - 1) {
      WA->Left_Delta[WA->Left_Delta_Len++] = max - last - 1;
      d--;
      last = WA->Edit_Array_Lazy[k - 1][from];
    } else if  (from == d + 1) {
      WA->Left_Delta[WA->Left_Delta_Len++] = last - (max - 1);
      d++;
      last = WA->Edit_Array_Lazy[k - 1][from];
    }
  }
  (* leftover) = last;

  // Don't allow first delta to be +1 or -1
  assert (WA->Left_Delta_Len == 0 || WA->Left_Delta[0] != -1);

  //  BPW - The original test was Left_Delta_Len>0, but that hit uninitialized data when Len == 1.

  if  (WA->Left_Delta_Len > 1 && WA->Left_Delta[0] == 1 && (* t_end) + t_len > 0) {
    int  i;

    if  (WA->Left_Delta[1] > 0)  //  <- uninitialized here
      WA->Left_Delta[0] = WA->Left_Delta[1] + 1;
    else
      WA->Left_Delta[0] = WA->Left_Delta[1] - 1;

    for  (i = 2;  i < WA->Left_Delta_Len;  i++)
      WA->Left_Delta[i - 1] = WA->Left_Delta[i];

    WA->Left_Delta_Len--;

    (* t_end)--;
    if  (WA->Left_Delta_Len == 0)
      (* leftover)++;
  }
}



//  Put the delta encoding of the alignment represented in WA->Edit_Array
//  starting at row e (which is the number of errors) and column d
//  (which is the diagonal) and working back to the start, into
//  WA->Right_Delta. Set WA->Right_Delta_Len to the number of
//  delta entries.

static
void
Set_Right_Delta (int e, int d, Work_Area_t * WA) {
  int  from, last, max;
  int  i, j, k;

  assert(WA->Edit_Array_Lazy[e] != NULL);

  last = WA->Edit_Array_Lazy[e][d];
  WA->Right_Delta_Len = 0;

  for  (k = e;  k > 0;  k--) {
    assert(WA->Edit_Array_Lazy[k] != NULL);

    from = d;
    max = 1 + WA->Edit_Array_Lazy[k - 1][d];
    if  ((j = WA->Edit_Array_Lazy[k - 1][d - 1]) > max) {
      from = d - 1;
      max = j;
    }
    if  ((j = 1 + WA->Edit_Array_Lazy[k - 1][d + 1]) > max) {
      from = d + 1;
      max = j;
    }
    if  (from == d - 1) {
      WA->Delta_Stack[WA->Right_Delta_Len++] = max - last - 1;
      d--;
      last = WA->Edit_Array_Lazy[k - 1][from];
    } else if  (from == d + 1) {
      WA->Delta_Stack[WA->Right_Delta_Len++] = last - (max - 1);
      d++;
      last = WA->Edit_Array_Lazy[k - 1][from];
    }
  }
  WA->Delta_Stack[WA->Right_Delta_Len++] = last + 1;

  k = 0;
  for  (i = WA->Right_Delta_Len - 1;  i > 0;  i--)
    WA->Right_Delta[k++]
      = abs (WA->Delta_Stack[i]) * Sign (WA->Delta_Stack[i - 1]);
  WA->Right_Delta_Len--;
}







//  Allocate another block of 64mb for edits

//  Needs to be at least:
//       52,432 to handle 40% error at  64k overlap
//      104,860 to handle 80% error at  64k overlap
//      209,718 to handle 40% error at 256k overlap
//      419,434 to handle 80% error at 256k overlap
//    3,355,446 to handle 40% error at   4m overlap
//    6,710,890 to handle 80% error at   4m overlap
//  Bigger means we can assign more than one Edit_Array[] in one allocation.

uint32  EDIT_SPACE_SIZE  = 1 * 1024 * 1024;

static
int *
Allocate_More_Edit_Space(Work_Area_t *WA) {

  //  Determine the last allocated block, and the last assigned block

  int32  b = 0;  //  Last edit array assigned
  int32  e = 0;  //  Last edit array assigned more space
  int32  a = 0;  //  Last allocated block

  while (WA->Edit_Array_Lazy[b] != NULL)
    b++;

  while (WA->Edit_Space_Lazy[a] != NULL)
    a++;

  //  Fill in the edit space array.  Well, not quite yet.  First, decide the minimum size.
  //
  //  Element [0] can access from [-2] to [2] = 5 elements.
  //  Element [1] can access from [-3] to [3] = 7 elements.
  //
  //  Element [e] can access from [-2-e] to [2+e] = 5 + e * 2 elements
  //
  //  So, our offset for this new block needs to put [e][0] at offset...

  int32 Offset = 2 + b;
  int32 Del    = 6 + b * 2;
  int32 Size   = EDIT_SPACE_SIZE;

  while (Size < Offset + Del)
    Size *= 2;

  //  Allocate another block

  WA->Edit_Space_Lazy[a] = new int [Size];

  //  And, now, fill in the edit space array.

  e = b;

  while (Offset + Del < Size) {
    WA->Edit_Array_Lazy[e++] = WA->Edit_Space_Lazy[a] + Offset;

    Offset += Del;
    Del    += 2;
  }

  if (e == b)
    fprintf(stderr, "Allocate_More_Edit_Space()-- ERROR: couldn't allocate enough space for even one more entry!  e=%d\n", e);
  assert(e != b);

  fprintf(stderr, "WorkArea %d allocates space %d of size %d for array %d through %d\n", WA->thread_id, a, Size, b, e-1);
}




//  Return the minimum number of changes (inserts, deletes, replacements)
//  needed to match string  A[0 .. (m-1)]  with a prefix of string
//   T[0 .. (n-1)]  if it's not more than  Error_Limit .
//  If no match, return the number of errors for the best match
//  up to a branch point.
//  Put delta description of alignment in  WA->Right_Delta  and set
//  WA->Right_Delta_Len  to the number of entries there if it's a complete
//  match.
//  Set  A_End  and  T_End  to the rightmost positions where the
//  alignment ended in  A  and  T , respectively.
//  Set  Match_To_End  true if the match extended to the end
//  of at least one string; otherwise, set it false to indicate
//  a branch point.

static
int
Prefix_Edit_Dist(char A[], int m, char T[], int n, int Error_Limit,
                 int * A_End, int * T_End, int * Match_To_End, Work_Area_t * WA) {
  //int  Delta_Stack[MAX_ERRORS];
  double  Score, Max_Score;
  int  Max_Score_Len = 0, Max_Score_Best_d = 0, Max_Score_Best_e = 0;
  int  Tail_Len;
  int  Best_d, Best_e, From, Last, Longest, Max, Row;
  int  Left, Right;
  int  d, e, i, j, k;

  assert (m <= n);
  Best_d = Best_e = Longest = 0;
  WA->Right_Delta_Len = 0;

  for  (Row = 0;  Row < m
          && (A[Row] == T[Row]
              || A[Row] == DONT_KNOW_CHAR
              || T[Row] == DONT_KNOW_CHAR);  Row++)
    ;

  if (WA->Edit_Array_Lazy[0] == NULL)
    Allocate_More_Edit_Space(WA);

  WA->Edit_Array_Lazy[0][0] = Row;

  if  (Row == m) {
    // Exact match
    (* A_End) = (* T_End) = m;
    (* Match_To_End) = TRUE;
    return  0;
  }

  Left = Right = 0;
  Max_Score = 0.0;
  for  (e = 1;  e <= Error_Limit;  e++) {
    if (WA->Edit_Array_Lazy[e] == NULL)
      Allocate_More_Edit_Space(WA);

    Left = MAX (Left - 1, -e);
    Right = MIN (Right + 1, e);
    WA->Edit_Array_Lazy[e - 1][Left] = -2;
    WA->Edit_Array_Lazy[e - 1][Left - 1] = -2;
    WA->Edit_Array_Lazy[e - 1][Right] = -2;
    WA->Edit_Array_Lazy[e - 1][Right + 1] = -2;

    for  (d = Left;  d <= Right;  d++) {
      Row = 1 + WA->Edit_Array_Lazy[e - 1][d];
      if  ((j = WA->Edit_Array_Lazy[e - 1][d - 1]) > Row)
        Row = j;
      if  ((j = 1 + WA->Edit_Array_Lazy[e - 1][d + 1]) > Row)
        Row = j;
      while  (Row < m && Row + d < n
              && (A[Row] == T[Row + d]
                  || A[Row] == DONT_KNOW_CHAR
                  || T[Row + d] == DONT_KNOW_CHAR))
        Row++;

      WA->Edit_Array_Lazy[e][d] = Row;

      if  (Row == m || Row + d == n) {
        //  Check for branch point here caused by uneven
        //  distribution of errors
        Score = Row * Branch_Match_Value - e;
        // Assumes  Branch_Match_Value
        //             - Branch_Error_Value == 1.0
        Tail_Len = Row - Max_Score_Len;
        if  ((G.Doing_Partial_Overlaps && Score < Max_Score)
             ||  (e > MIN_BRANCH_END_DIST / 2
                  && Tail_Len >= MIN_BRANCH_END_DIST
                  && (Max_Score - Score) / Tail_Len >= MIN_BRANCH_TAIL_SLOPE)) {
          (* A_End) = Max_Score_Len;
          (* T_End) = Max_Score_Len + Max_Score_Best_d;
          Set_Right_Delta (Max_Score_Best_e, Max_Score_Best_d, WA);
          (* Match_To_End) = FALSE;
          return  Max_Score_Best_e;
        }

        // Force last error to be mismatch rather than insertion
        if  (Row == m && 1 + WA->Edit_Array_Lazy[e - 1][d + 1] == WA->Edit_Array_Lazy[e][d] && d < Right) {
          d++;
          WA->Edit_Array_Lazy[e][d] = WA->Edit_Array_Lazy[e][d - 1];
        }

        (* A_End) = Row;           // One past last align position
        (* T_End) = Row + d;
        Set_Right_Delta (e, d, WA);
        (* Match_To_End) = TRUE;
        return  e;
      }
    }

    while  (Left <= Right && Left < 0
            && WA->Edit_Array_Lazy[e][Left] < WA->Edit_Match_Limit[e])
      Left++;

    if  (Left >= 0)
      while  (Left <= Right
              && WA->Edit_Array_Lazy[e][Left] + Left < WA->Edit_Match_Limit[e])
        Left++;

    if  (Left > Right)
      break;

    while  (Right > 0
            && WA->Edit_Array_Lazy[e][Right] + Right < WA->Edit_Match_Limit[e])
      Right--;

    if  (Right <= 0)
      while  (WA->Edit_Array_Lazy[e][Right] < WA->Edit_Match_Limit[e])
        Right--;

    assert (Left <= Right);

    for  (d = Left;  d <= Right;  d++)
      if  (WA->Edit_Array_Lazy[e][d] > Longest) {
        Best_d = d;
        Best_e = e;
        Longest = WA->Edit_Array_Lazy[e][d];
      }

    Score = Longest * Branch_Match_Value - e;

    // Assumes  Branch_Match_Value - Branch_Error_Value == 1.0
    if  (Score > Max_Score) {
      Max_Score = Score;
      Max_Score_Len = Longest;
      Max_Score_Best_d = Best_d;
      Max_Score_Best_e = Best_e;
    }
  }

  (* A_End) = Max_Score_Len;
  (* T_End) = Max_Score_Len + Max_Score_Best_d;
  Set_Right_Delta (Max_Score_Best_e, Max_Score_Best_d, WA);
  (* Match_To_End) = FALSE;
  return  Max_Score_Best_e;
}







//  Return the minimum number of changes (inserts, deletes, replacements)
//  needed to match string  A[0 .. (1-m)]  right-to-left with a prefix of string
//   T[0 .. (1-n)]  right-to-left if it's not more than  Error_Limit .
//  If no match, return the number of errors for the best match
//  up to a branch point.
//  Put delta description of alignment in  WA->Left_Delta  and set
//  WA->Left_Delta_Len  to the number of entries there.
//  Set  A_End  and  T_End  to the leftmost positions where the
//  alignment ended in  A  and  T , respectively.
//  If the alignment succeeds set  Leftover  to the number of
//  characters that match after the last  WA->Left_Delta  entry;
//  otherwise, set  Leftover  to zero.
//  Set  Match_To_End  true if the match extended to the end
//  of at least one string; otherwise, set it false to indicate
//  a branch point.

static
int
Rev_Prefix_Edit_Dist (char A[], int m, char T[], int n, int Error_Limit,
                      int * A_End, int * T_End, int * Leftover, int * Match_To_End,
                      Work_Area_t * WA) {

  double  Score, Max_Score;
  int  Max_Score_Len = 0, Max_Score_Best_d = 0, Max_Score_Best_e = 0;
  int  Tail_Len;
  int  Best_d, Best_e, From, Last, Longest, Max, Row;
  int  Left, Right;
  int  d, e, j, k;


  assert (m <= n);
  Best_d = Best_e = Longest = 0;
  WA->Left_Delta_Len = 0;

  for  (Row = 0;  Row < m
          && (A[- Row] == T[- Row]
              || A[- Row] == DONT_KNOW_CHAR
              || T[- Row] == DONT_KNOW_CHAR);  Row++)
    ;

  if (WA->Edit_Array_Lazy[0] == NULL)
    Allocate_More_Edit_Space(WA);

  WA->Edit_Array_Lazy[0][0] = Row;

  if  (Row == m) {
    (* A_End) = (* T_End) = - m;
    (* Leftover) = m;
    (* Match_To_End) = TRUE;
    return  0;
  }

  Left = Right = 0;
  Max_Score = 0.0;
  for  (e = 1;  e <= Error_Limit;  e++) {
    if (WA->Edit_Array_Lazy[e] == NULL)
      Allocate_More_Edit_Space(WA);

    Left = MAX (Left - 1, -e);
    Right = MIN (Right + 1, e);
    WA->Edit_Array_Lazy[e - 1][Left] = -2;
    WA->Edit_Array_Lazy[e - 1][Left - 1] = -2;
    WA->Edit_Array_Lazy[e - 1][Right] = -2;
    WA->Edit_Array_Lazy[e - 1][Right + 1] = -2;

    for  (d = Left;  d <= Right;  d++) {
      Row = 1 + WA->Edit_Array_Lazy[e - 1][d];
      if  ((j = WA->Edit_Array_Lazy[e - 1][d - 1]) > Row)
        Row = j;
      if  ((j = 1 + WA->Edit_Array_Lazy[e - 1][d + 1]) > Row)
        Row = j;
      while  (Row < m && Row + d < n
              && (A[- Row] == T[- Row - d]
                  || A[- Row] == DONT_KNOW_CHAR
                  || T[- Row - d] == DONT_KNOW_CHAR))
        Row++;

      WA->Edit_Array_Lazy[e][d] = Row;

      if  (Row == m || Row + d == n) {

        //  Check for branch point here caused by uneven
        //  distribution of errors

        Score = Row * Branch_Match_Value - e;
        // Assumes  Branch_Match_Value
        //             - Branch_Error_Value == 1.0
        Tail_Len = Row - Max_Score_Len;
        if  ((G.Doing_Partial_Overlaps && Score < Max_Score)
             || (e > MIN_BRANCH_END_DIST / 2
                 && Tail_Len >= MIN_BRANCH_END_DIST
                 && (Max_Score - Score) / Tail_Len >= MIN_BRANCH_TAIL_SLOPE)) {
          (* A_End) = - Max_Score_Len;
          (* T_End) = - Max_Score_Len - Max_Score_Best_d;
          Set_Left_Delta (Max_Score_Best_e, Max_Score_Best_d,
                          Leftover, T_End, n, WA);
          (* Match_To_End) = FALSE;
          return  Max_Score_Best_e;
        }

        (* A_End) = - Row;           // One past last align position
        (* T_End) = - Row - d;
        Set_Left_Delta (e, d, Leftover, T_End, n, WA);
        (* Match_To_End) = TRUE;
        return  e;
      }
    }

    while  (Left <= Right && Left < 0
            && WA->Edit_Array_Lazy[e][Left] < WA->Edit_Match_Limit[e])
      Left++;

    if  (Left >= 0)
      while  (Left <= Right
              && WA->Edit_Array_Lazy[e][Left] + Left < WA->Edit_Match_Limit[e])
        Left++;

    if  (Left > Right)
      break;

    while  (Right > 0
            && WA->Edit_Array_Lazy[e][Right] + Right < WA->Edit_Match_Limit[e])
      Right--;

    if  (Right <= 0)
      while  (WA->Edit_Array_Lazy[e][Right] < WA->Edit_Match_Limit[e])
        Right--;

    assert (Left <= Right);

    for  (d = Left;  d <= Right;  d++)
      if  (WA->Edit_Array_Lazy[e][d] > Longest) {
        Best_d = d;
        Best_e = e;
        Longest = WA->Edit_Array_Lazy[e][d];
      }

    Score = Longest * Branch_Match_Value - e;

    // Assumes  Branch_Match_Value - Branch_Error_Value == 1.0
    if  (Score > Max_Score) {
      Max_Score = Score;
      Max_Score_Len = Longest;
      Max_Score_Best_d = Best_d;
      Max_Score_Best_e = Best_e;
    }
  }

  (* A_End) = - Max_Score_Len;
  (* T_End) = - Max_Score_Len - Max_Score_Best_d;
  Set_Left_Delta (Max_Score_Best_e, Max_Score_Best_d, Leftover, T_End, n, WA);
  (* Match_To_End) = FALSE;
  return  Max_Score_Best_e;
}



//  See how far the exact match in  Match  extends.  The match
//  refers to strings  S  and  T  with lengths  S_Len  and  T_Len ,
//  respectively.  Set  S_Lo ,  S_Hi ,  T_Lo  and  T_Hi  to the
//  regions within  S  and  T  to which the match extends.
//  Return the type of overlap:  NONE if doesn't extend to
//  the end of either fragment;  LEFT_BRANCH_PT
//  or  RIGHT_BRANCH_PT  if match extends to the end of just one fragment;
//  or  DOVETAIL  if it extends to the end of both fragments, i.e.,
//  it is a complete overlap.
//  Set  Errors  to the number of errors in the alignment if it is
//  a  DOVETAIL  overlap.

Overlap_t
Extend_Alignment(Match_Node_t * Match, char * S, int S_Len, char * T, int T_Len,
                 int * S_Lo, int * S_Hi, int * T_Lo, int * T_Hi, int * Errors,
                 Work_Area_t * WA) {
  Overlap_t  return_type;
  int  S_Left_Begin, S_Right_Begin, S_Right_Len;
  int  T_Left_Begin, T_Right_Begin, T_Right_Len;
  int  Error_Limit, Left_Errors, Right_Errors, Total_Olap;
  int  i, Leftover, Right_Match_To_End, Left_Match_To_End;

  S_Left_Begin  = Match->Start - 1;
  S_Right_Begin = Match->Start + Match->Len;
  S_Right_Len   = S_Len - S_Right_Begin;

  T_Left_Begin  = Match->Offset - 1;
  T_Right_Begin = Match->Offset + Match->Len;
  T_Right_Len   = T_Len - T_Right_Begin;

  Total_Olap = (MIN(Match->Start, Match->Offset) +
                Match->Len +
                MIN(S_Right_Len, T_Right_Len));

  Error_Limit = WA->Error_Bound[Total_Olap];

  if  (S_Right_Len == 0 || T_Right_Len == 0) {
    Right_Errors = 0;
    WA->Right_Delta_Len = 0;
    (* S_Hi) = (* T_Hi) = 0;
    Right_Match_To_End = TRUE;

  } else if  (S_Right_Len <= T_Right_Len) {
    Right_Errors = Prefix_Edit_Dist (S + S_Right_Begin, S_Right_Len,
                                     T + T_Right_Begin, T_Right_Len, Error_Limit,
                                     S_Hi, T_Hi, & Right_Match_To_End, WA);
  } else {
    Right_Errors = Prefix_Edit_Dist (T + T_Right_Begin, T_Right_Len,
                                     S + S_Right_Begin, S_Right_Len, Error_Limit,
                                     T_Hi, S_Hi, & Right_Match_To_End, WA);
  }

  for  (i = 0;  i < WA->Right_Delta_Len;  i++)
    WA->Right_Delta[i] *= -1;

  (* S_Hi) += S_Right_Begin - 1;
  (* T_Hi) += T_Right_Begin - 1;

  assert (Right_Errors <= Error_Limit);

  if  (S_Left_Begin < 0 || T_Left_Begin < 0) {
    Left_Errors = 0;
    WA->Left_Delta_Len = 0;
    (* S_Lo) = (* T_Lo) = 0;
    Leftover = 0;
    Left_Match_To_End = TRUE;
  } else if  (S_Right_Begin <= T_Right_Begin) {
    Left_Errors = Rev_Prefix_Edit_Dist (S + S_Left_Begin,
                                        S_Left_Begin + 1, T + T_Left_Begin,
                                        T_Left_Begin + 1,
                                        Error_Limit - Right_Errors,
                                        S_Lo, T_Lo, & Leftover, & Left_Match_To_End,
                                        WA);
  } else {
    Left_Errors = Rev_Prefix_Edit_Dist (T + T_Left_Begin,
                                        T_Left_Begin + 1, S + S_Left_Begin,
                                        S_Left_Begin + 1,
                                        Error_Limit - Right_Errors,
                                        T_Lo, S_Lo, & Leftover, & Left_Match_To_End,
                                        WA);
  }

  for  (i = 0;  i < WA->Left_Delta_Len;  i++)
    WA->Left_Delta[i] *= -1;

  (* S_Lo) += S_Left_Begin + 1;        // Check later for branch points
  (* T_Lo) += T_Left_Begin + 1;        // Check later for branch points

  if  (! Right_Match_To_End) {
    if  (! G.Doing_Partial_Overlaps)
      WA->Left_Delta_Len = 0;
    if  (! Left_Match_To_End)
      return_type = NONE;
    else
      return_type = RIGHT_BRANCH_PT;
  } else {
    if  (! Left_Match_To_End)
      return_type = LEFT_BRANCH_PT;
    else
      return_type = DOVETAIL;
  }

  if  (return_type == DOVETAIL || G.Doing_Partial_Overlaps) {
    (* Errors) = Left_Errors + Right_Errors;
    assert ((* Errors) <= Error_Limit);

    if  (WA->Right_Delta_Len > 0) {
      if  (WA->Right_Delta[0] > 0)
        WA->Left_Delta[WA->Left_Delta_Len++] = WA->Right_Delta[0] + Leftover + Match->Len;
      else
        WA->Left_Delta[WA->Left_Delta_Len++] = WA->Right_Delta[0] - Leftover - Match->Len;
    }
    for  (i = 1;  i < WA->Right_Delta_Len;  i++)
      WA->Left_Delta[WA->Left_Delta_Len++] = WA->Right_Delta[i];
  }

  return  return_type;
}


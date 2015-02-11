
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
/*************************************************
 * Module:  CorrectOlapsOVL.c
 * Description:
 *   Based on overlaps between DNA fragment sequences, make corrections
 *   to single bases in the sequences.
 *
 *    Programmer:  A. Delcher
 *       Started:  11 Dec 2000
 *
 * Assumptions:
 *
 * Notes:
 *
 *************************************************/

const char *mainid = "$Id: CorrectOlapsOVL.C 6709 2015-02-06 09:39:26Z bri $";

#include  "AS_global.H"

#include  "AS_PER_gkpStore.H"
#include  "AS_MSG_pmesg.H"
#include  "AS_UTL_reverseComplement.H"
#include  "FragCorrectOVL.H"
#include  "AS_OVS_overlapStore.H"





int Verbose_Level = 0;



static
void
Compute_Delta(Thread_Work_Area_t   *WA,
              int32                 e,
              int32                 d,
              int32                 row) {
  int32  last = row;

  WA->deltaLen = 0;

  for (int32 k=e; k>0; k--) {
    int32  from = d;
    int32  max  = 1 + WA->Edit_Array_Lazy[k-1][d];
    if ((j = WA->Edit_Array_Lazy[k-1][d - 1]) > max) {
      from = d - 1;
      max = j;
    }
    if ((j = 1 + WA->Edit_Array_Lazy[k-1][d + 1]) > max) {
      from = d + 1;
      max = j;
    }
    if (from == d - 1) {
      WA->deltaStack[WA->deltaLen++] = max - last - 1;
      d--;
      last = WA->Edit_Array_Lazy[k-1][from];
    } else if (from == d + 1) {
      WA->deltaStack[WA->deltaLen++] = last - (max - 1);
      d++;
      last = WA->Edit_Array_Lazy[k-1][from];
    }
  }

  WA->deltaStack[WA->deltaLen++] = last + 1;
  
  k = 0;

  for (i = WA->deltaLen - 1;  i > 0;  i--)
    Delta[k++]
      = abs (WA->deltaStack[i]) * Sign (WA->deltaStack[i - 1]);

  WA->deltaLen--;
}






//  Allocate another block of 64mb for edits

//  Needs to be at least:
//       52,432 to handle 40% error at  64k overlap
//      104,860 to handle 80% error at  64k overlap
//      209,718 to handle 40% error at 256k overlap
//      419,434 to handle 80% error at 256k overlap
//    3,355,446 to handle 40% error at   4m overlap
//    6,710,890 to handle 80% error at   4m overlap
//  Bigger means we can assign more than one WA->Edit_Array[] in one allocation.

uint32  EDIT_SPACE_SIZE  = 16 * 1024 * 1024;

static
void
Allocate_More_Edit_Space(void) {

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
  //  Element[0] can access from[-2] to[2] = 5 elements.
  //  Element[1] can access from[-3] to[3] = 7 elements.
  //
  //  Element[e] can access from[-2-e] to[2+e] = 5 + e * 2 elements
  //
  //  So, our offset for this new block needs to put[e][0] at offset...

  int32 Offset = 2 + b;
  int32 Del    = 6 + b * 2;
  int32 Size   = EDIT_SPACE_SIZE;

  while (Size < Offset + Del)
    Size *= 2;

  //  Allocate another block

  WA->Edit_Space_Lazy[a] = new int[Size];

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
}





//  Return the minimum number of changes (inserts, deletes, replacements)
//  needed to match string  A[0 .. (m-1)]  with a prefix of string
//   T[0 .. (n-1)]  if it's not more than  Error_Limit .
//  Put delta description of alignment in  Delta  and set
//  Delta_Len  to the number of entries there if it's a complete
//  match.
//  Set  A_End  and  T_End  to the rightmost positions where the
//  alignment ended in  A  and  T , respectively.
//  Set  Match_To_End  true if the match extended to the end
//  of at least one string; otherwise, set it false to indicate
//  a branch point.

static
int
Prefix_Edit_Dist(char A[], int m, char T[], int n, int Error_Limit,
                     int * A_End, int * T_End, int * Match_To_End,
                     int * Delta, int * Delta_Len) {

  int  Delta_Stack[AS_READ_MAX_NORMAL_LEN+1];  //  only MAX_ERRORS needed
  double  Score, Max_Score;
  int  Max_Score_Len, Max_Score_Best_d, Max_Score_Best_e;
  int  Best_d, Best_e, From, Last, Longest, Max, Row;
  int  Left, Right;
  int  d, e, i, j, k, shorter;

  //   assert (m <= n);
  Best_d = Best_e = Longest = 0;
  Delta_Len = 0;

  shorter = OVL_Min_int (m, n);
  for (Row = 0;  Row < shorter && A[Row] == T[Row];  Row++)
    ;

  if (WA->Edit_Array_Lazy[0] == NULL)
    Allocate_More_Edit_Space(WA);

  WA->Edit_Array_Lazy[0][0] = Row;

  // Exact match?
  if (Row == shorter) {
    A_End        = Row;
    T_End        = Row;
    Match_To_End = true;
    return(0);
  }

  int32  Left = 0;
  int32  Right = 0;
  double Max_Score = 0.0;
  int32  Max_Score_Len = 0;
  int32  Max_Score_Best_d = 0;
  int32  Max_Score_Best_e = 0;

  for (int32 e=1;  e<=Error_Limit;  e++) {
    if (WA->Edit_Array_Lazy[e] == NULL)
      Allocate_More_Edit_Space(WA);

    Left  = max(Left  - 1, -e);
    Right = min(Right + 1,  e);

    WA->Edit_Array_Lazy[e-1][Left]    = -2;
    WA->Edit_Array_Lazy[e-1][Left-1]  = -2;
    WA->Edit_Array_Lazy[e-1][Right]   = -2;
    WA->Edit_Array_Lazy[e-1][Right+1] = -2;

    for (int32 d=Left; d<=Right; d++) {
      Row = 1 + WA->Edit_Array_Lazy[e-1][d];
      Row = max(Row, WA->Edit_Array_Lazy[e-1][d-1]);
      Row = max(Row, WA->Edit_Array_Lazy[e-1][d+1] + 1);

      while ((Row < m) && (Row + d < n) && (A[Row] == T[Row + d]))
        Row++;

      assert(e < MAX_ERRORS);
      //assert(d < ??);

      WA->Edit_Array_Lazy[e][d] = Row;

      if (Row == m || Row + d == n) {

        // Force last error to be mismatch rather than insertion
        if (Row == m
            && 1 + WA->Edit_Array_Lazy[e-1][d + 1]
            == WA->Edit_Array_Lazy[e][d]
            && d < Right) {
          d++;
          WA->Edit_Array_Lazy[e][d] = WA->Edit_Array_Lazy[e][d - 1];
        }

        A_End        = Row;           // One past last align position
        T_End        = Row + d;
        Match_To_End = true;

        Compute_Delta(WA, e, d, Row);

        return(e);
      }
    }

    while (Left <= Right && Left < 0
           && WA->Edit_Array_Lazy[e][Left] < WA->G->Edit_Match_Limit[e])
      Left++;

    if (Left >= 0)
      while (Left <= Right
             && WA->Edit_Array_Lazy[e][Left] + Left < WA->G->Edit_Match_Limit[e])
        Left++;

    if (Left > Right)
      break;

    while (Right > 0
           && WA->Edit_Array_Lazy[e][Right] + Right < WA->G->Edit_Match_Limit[e])
      Right--;

    if (Right <= 0)
      while (WA->Edit_Array_Lazy[e][Right] < WA->G->Edit_Match_Limit[e])
        Right--;

    assert (Left <= Right);

    for (int32 d=Left;  d <= Right;  d++)
      if (WA->Edit_Array_Lazy[e][d] > Longest) {
        Best_d = d;
        Best_e = e;
        Longest = WA->Edit_Array_Lazy[e][d];
      }

#if  1
    int32 Score = Longest * BRANCH_PT_MATCH_VALUE - e;

    // Assumes  BRANCH_PT_MATCH_VALUE - BRANCH_PT_ERROR_VALUE == 1.0
    if (Score > Max_Score) {
      Max_Score = Score;
      Max_Score_Len = Longest;
      Max_Score_Best_d = Best_d;
      Max_Score_Best_e = Best_e;
    }
#endif
  }

  A_End        = Max_Score_Len;
  T_End        = Max_Score_Len + Max_Score_Best_d;
  Match_To_End = false;

  return(e);
}


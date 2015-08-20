
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
 *  Modifications by:
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "findErrors.H"

//  Set  delta  to the entries indicating the insertions/deletions
//  in the alignment encoded in  edit_array  ending at position
//  edit_array[e][d].  row  is the position in the first
//  string where the alignment ended.  Set  (*delta_len)  to
//  the number of entries in  delta .

static
void
Compute_Delta(pedWorkArea_t   *WA,
              int32            e,
              int32            d,
              int32            row) {
  int32  last     = row;
  int32  stackLen = 0;

  for (int32 k=e; k>0; k--) {
    int32  from = d;
    int32  max  = 1 + WA->Edit_Array_Lazy[k-1][d];

    int32  j = WA->Edit_Array_Lazy[k-1][d-1];

    if (j > max) {
      from = d-1;
      max = j;
    }

    j = 1 + WA->Edit_Array_Lazy[k-1][d+1];

    if (j > max) {
      from = d+1;
      max = j;
    }

    if (from == d-1) {
      WA->deltaStack[stackLen++] = max - last - 1;
      d--;
      last = WA->Edit_Array_Lazy[k-1][from];
    }

    else if (from == d+1) {
      WA->deltaStack[stackLen++] = last - (max - 1);
      d++;
      last = WA->Edit_Array_Lazy[k-1][from];
    }
  }

  WA->deltaStack[stackLen++] = last + 1;

  for (int32 k=0, i=stackLen-1; i>0; i--)
    WA->delta[k++] = abs(WA->deltaStack[i]) * Sign(WA->deltaStack[i-1]);

  WA->deltaLen = stackLen - 1;
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

uint32  EDIT_SPACE_SIZE  = 16 * 1024 * 1024;

static
void
Allocate_More_Edit_Space(pedWorkArea_t *WA) {

  //  Determine the last allocated block, and the last assigned block

  int32  b = 0;                 //  Last edit array assigned
  int32  e = 0;                 //  Last edit array assigned more space
  int32  a = WA->alloc.size();  //  First unallocated block

  while (WA->Edit_Array_Lazy[b] != NULL)
    b++;

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

  int32  *alloc = new int32 [Size];

  WA->alloc.push_back(alloc);

  //  And, now, fill in the edit space array.

  e = b;

  while ((Offset + Del < Size) &&
         (e < WA->Edit_Array_Max)) {
    WA->Edit_Array_Lazy[e++] = alloc + Offset;

    Offset += Del;
    Del    += 2;
  }

  if (e == b)
    fprintf(stderr, "Allocate_More_Edit_Space()-- ERROR: couldn't allocate enough space for even one more entry!  e=%d\n", e);
  assert(e != b);

  //fprintf(stderr, "WorkArea %d allocates space %d of size %d for array %d through %d\n", WA->thread_id, a, Size, b, e-1);
}





//  Return the minimum number of changes (inserts, deletes, replacements)
//  needed to match string  A[0 .. (m-1)]  with a prefix of string
//   T[0 .. (n-1)]  if it's not more than  Error_Limit .
//
//  Put delta description of alignment in  WA->delta  and set
//  WA->deltaLen to the number of entries there if it's a complete
//  match.
//  Set  A_End  and  T_End  to the rightmost positions where the
//  alignment ended in  A  and  T , respectively.
//  Set  Match_To_End  true if the match extended to the end
//  of at least one string; otherwise, set it false to indicate
//  a branch point.

int32
Prefix_Edit_Dist(char   *A, int32 m,
                 char   *T, int32 n,
                 int32   Error_Limit,
                 int32  &A_End,
                 int32  &T_End,
                 bool   &Match_To_End,
                 pedWorkArea_t *WA) {

  //assert (m <= n);

  int32 Best_d  = 0;
  int32 Best_e  = 0;
  int32 Longest = 0;

  WA->deltaLen = 0;

  int32 shorter = min(m, n);

  int32 Row = 0;
  while ((Row < shorter) && (A[Row] == T[Row]))
    Row++;

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

  for (int32 e=1; e<=Error_Limit; e++) {
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

      assert(e < WA->Edit_Array_Max);

      WA->Edit_Array_Lazy[e][d] = Row;

      if (Row == m || Row + d == n) {

        // Force last error to be mismatch rather than insertion
        if ((Row == m) &&
            (1 + WA->Edit_Array_Lazy[e-1][d+1] == WA->Edit_Array_Lazy[e][d]) &&
            (d < Right)) {
          d++;
          WA->Edit_Array_Lazy[e][d] = WA->Edit_Array_Lazy[e][d-1];
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
        Best_d  = d;
        Best_e  = e;
        Longest = WA->Edit_Array_Lazy[e][d];
      }

    int32  Score = Longest * BRANCH_PT_MATCH_VALUE - e;

    // Assumes  BRANCH_PT_MATCH_VALUE - BRANCH_PT_ERROR_VALUE == 1.0

    //  CorrectOverlaps didn't have the second clause.
    //  Neither did overlapper.

    if ((Score > Max_Score) &&
        (Best_e <= WA->G->Error_Bound[min(Longest, Longest + Best_d)])) {
      Max_Score        = Score;
      Max_Score_Len    = Longest;
      Max_Score_Best_d = Best_d;
      Max_Score_Best_e = Best_e;
    }
  }

  //  CorrectOverlaps doesn't have this call
  Compute_Delta(WA, Max_Score_Best_e, Max_Score_Best_d, Max_Score_Len);

  A_End        = Max_Score_Len;
  T_End        = Max_Score_Len + Max_Score_Best_d;
  Match_To_End = false;

  //  CorrectOverlaps was returning just 'e', not the best as below.
  //  It used this return value to compute the new error rate.

  return(Max_Score_Best_e);
}

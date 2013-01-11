
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

static const char *rcsid = "$Id: overlapInCore-Process_String_Overlaps.C,v 1.2 2013-01-11 11:11:04 brianwalenz Exp $";

#include "overlapInCore.H"


static void  Combine_Into_One_Olap
(Olap_Info_t olap [], int ct, int deleted [])

//  Choose the best overlap in  olap [0 .. (ct - 1)] .
//  Mark all others as deleted (by setting  deleted []  true for them)
//  and combine their information in the min/max entries in the
//  best one.

{
  int  min_diag, max_diag;
  int  s_left_boundary, s_right_boundary;
  int  t_left_boundary, t_right_boundary;
  int  i, best;

  best = 0;
  min_diag = olap [0] . min_diag;
  max_diag = olap [0] . max_diag;
  s_left_boundary = olap [0] . s_left_boundary;
  s_right_boundary = olap [0] . s_right_boundary;
  t_left_boundary = olap [0] . t_left_boundary;
  t_right_boundary = olap [0] . t_right_boundary;

  for  (i = 1;  i < ct;  i ++)
    {
      if  (olap [i] . quality < olap [best] . quality)
        best = i;
      if  (olap [i] . min_diag < min_diag)
        min_diag = olap [i] . min_diag;
      if  (olap [i] . max_diag > max_diag)
        max_diag = olap [i] . max_diag;
      if  (olap [i] . s_left_boundary < s_left_boundary)
        s_left_boundary = olap [i] . s_left_boundary;
      if  (olap [i] . s_right_boundary > s_right_boundary)
        s_right_boundary = olap [i] . s_right_boundary;
      if  (olap [i] . t_left_boundary < t_left_boundary)
        t_left_boundary = olap [i] . t_left_boundary;
      if  (olap [i] . t_right_boundary > t_right_boundary)
        t_right_boundary = olap [i] . t_right_boundary;
    }

  olap [best] . min_diag = min_diag;
  olap [best] . max_diag = max_diag;
  olap [best] . s_left_boundary = s_left_boundary;
  olap [best] . s_right_boundary = s_right_boundary;
  olap [best] . t_left_boundary = t_left_boundary;
  olap [best] . t_right_boundary = t_right_boundary;

  for  (i = 0;  i < ct;  i ++)
    deleted [i] = (i != best);

  return;
}











static void  Merge_Intersecting_Olaps
(Olap_Info_t p [], int ct, int deleted [])

//  Combine overlaps whose overlap regions intersect sufficiently
//  in  p [0 .. (ct - 1)]  by marking the poorer quality one
//  deleted (by setting  deleted []  true for it) and combining
//  its min/max info in the other.  Assume all entries in
//  deleted are 0 initially.

{
  int  i, j, lo_diag, hi_diag;

  for  (i = 0;  i < ct - 1;  i ++)
    for  (j = i + 1;  j < ct;  j ++)
      {
        if  (deleted [i] || deleted [j])
          continue;
        lo_diag = p [i] . min_diag;
        hi_diag = p [i] . max_diag;
        if  ((lo_diag <= 0 && p [j] . min_diag > 0)
             || (lo_diag > 0 && p [j] . min_diag <= 0))
          continue;
        if  ((lo_diag >= 0
              && p [j] .  t_right_boundary - lo_diag
              - p [j] . s_left_boundary
              >= MIN_INTERSECTION)
             || (lo_diag <= 0
                 && p [j] .  s_right_boundary + lo_diag
                 - p [j] .  t_left_boundary
                 >= MIN_INTERSECTION)
             || (hi_diag >= 0
                 && p [j] .  t_right_boundary - hi_diag
                 - p [j] . s_left_boundary
                 >= MIN_INTERSECTION)
             || (hi_diag <= 0
                 && p [j] .  s_right_boundary + hi_diag
                 - p [j] .  t_left_boundary
                 >= MIN_INTERSECTION))
          {
            Olap_Info_t  * discard, * keep;

            if  (p [i] . quality < p [j] . quality)
              {
                keep = p + i;
                discard = p + j;
                deleted [j] = TRUE;
              }
            else
              {
                keep = p + j;
                discard = p + i;
                deleted [i] = TRUE;
              }
            if  (discard -> min_diag < keep -> min_diag)
              keep -> min_diag = discard -> min_diag;
            if  (discard -> max_diag > keep -> max_diag)
              keep -> max_diag = discard -> max_diag;
            if  (discard -> s_left_boundary < keep -> s_left_boundary)
              keep -> s_left_boundary = discard -> s_left_boundary;
            if  (discard -> s_right_boundary > keep -> s_right_boundary)
              keep -> s_right_boundary = discard -> s_right_boundary;
            if  (discard -> t_left_boundary < keep -> t_left_boundary)
              keep -> t_left_boundary = discard -> t_left_boundary;
            if  (discard -> t_right_boundary > keep -> t_right_boundary)
              keep -> t_right_boundary = discard -> t_right_boundary;
          }
      }

  return;
}





static void  Add_Overlap
(int s_lo, int s_hi, int t_lo, int t_hi, double qual,
 Olap_Info_t * olap, int * ct, Work_Area_t * WA)

//  Add information for the overlap between strings  S  and  T
//  at positions  s_lo .. s_hi  and  t_lo .. t_hi , resp., and
//  with quality  qual  to the array  olap []  which
//  currently has  (* ct)  entries.  Increment  (* ct)  if this
//  is a new, distinct overlap; otherwise, modify an existing
//  entry if this is just a "slide" of an existing overlap.

{
  int  i, new_diag;

#ifdef COMPARE
  fprintf(stderr, "Add_Overlap()  %d %d %d %d qual %f\n", s_lo, s_hi, t_lo, t_hi, qual);
#endif

  if  (Doing_Partial_Overlaps)
    olap += (* ct);
  // Don't combine overlaps when doing partials
  else
    {

      new_diag = t_lo - s_lo;
      for  (i = 0;  i < (* ct);  i ++)
        {
          int  old_diag = olap -> t_lo - olap -> s_lo;

          if  ((new_diag > 0
                && old_diag > 0
                && olap -> t_right_boundary - new_diag - olap -> s_left_boundary
                >= MIN_INTERSECTION)
               || (new_diag <= 0
                   && old_diag <= 0
                   && olap -> s_right_boundary + new_diag - olap -> t_left_boundary
                   >= MIN_INTERSECTION))
            {
              if  (new_diag < olap -> min_diag)
                olap -> min_diag = new_diag;
              if  (new_diag > olap -> max_diag)
                olap -> max_diag = new_diag;
              if  (s_lo < olap -> s_left_boundary)
                olap -> s_left_boundary = s_lo;
              if  (s_hi > olap -> s_right_boundary)
                olap -> s_right_boundary = s_hi;
              if  (t_lo < olap -> t_left_boundary)
                olap -> t_left_boundary = t_lo;
              if  (t_hi > olap -> t_right_boundary)
                olap -> t_right_boundary = t_hi;
              if  (qual < olap -> quality)      // lower value is better
                {
                  olap -> s_lo = s_lo;
                  olap -> s_hi = s_hi;
                  olap -> t_lo = t_lo;
                  olap -> t_hi = t_hi;
                  olap -> quality = qual;
                  memcpy (& (olap ->  delta), WA -> Left_Delta,
                          WA -> Left_Delta_Len * sizeof (int));
                  olap -> delta_ct = WA -> Left_Delta_Len;
                }

              //  check for intersections before outputting
              return;
            }

          olap ++;
        }
    }

  if  ((* ct) >= MAX_DISTINCT_OLAPS)
    return;   // no room for a new entry; this shouldn't happen

  olap -> s_lo = olap -> s_left_boundary = s_lo;
  olap -> s_hi = olap -> s_right_boundary = s_hi;
  olap -> t_lo = olap -> t_left_boundary = t_lo;
  olap -> t_hi = olap -> t_right_boundary = t_hi;
  olap -> quality = qual;
  memcpy (& (olap ->  delta), WA -> Left_Delta,
          WA -> Left_Delta_Len * sizeof (int));
  olap -> delta_ct = WA -> Left_Delta_Len;
  olap -> min_diag = olap -> max_diag = t_lo - s_lo;
  (* ct) ++;

  return;
}







static int  Lies_On_Alignment
(int start, int offset, int s_lo, int t_lo, Work_Area_t * WA)

//  Return  TRUE  iff the exact match region beginning at
//  position  start  in the first string and  offset  in
//  the second string lies along the alignment from
//   lo .. hi  on the first string where the delta-encoding
//  of the alignment is given by
//   WA -> Left_Delta [0 .. (WA -> Left_Delta_Len-1)] .

{
  int  i, diag, new_diag;

  diag = t_lo - s_lo;
  new_diag = offset - start;

  for  (i = 0;  i < WA -> Left_Delta_Len;  i ++)
    {
      s_lo += abs (WA -> Left_Delta [i]);
      if  (start < s_lo)
        return  (abs (new_diag - diag) <= SHIFT_SLACK);
      if  (WA -> Left_Delta [i] < 0)
        diag ++;
      else
        {
          s_lo ++;
          diag --;
        }
    }

  return  (abs (new_diag - diag) <= SHIFT_SLACK);
}



static void  Choose_Best_Partial
(Olap_Info_t * olap, int ct, int deleted [])

//  Choose the best partial overlap in  olap [0 .. (ct - 1)] .
//  Set the corresponding  deleted  entry to false for the others.
//  Best is the greatest number of matching bases.

{
  double  matching_bases;
  int  i, best;

  best = 0;
  matching_bases = (1.0 - olap [0] . quality) * (2 + olap [0] . s_hi - olap [0] . s_lo
                                                 + olap [0] . t_hi - olap [0] . t_lo);
  // actually twice the number of matching bases but the max will be
  // the same overlap

  for  (i = 1;  i < ct;  i ++)
    {
      double  mb;

      mb = (1.0 - olap [i] . quality) * (2 + olap [i] . s_hi - olap [i] . s_lo
                                         + olap [i] . t_hi - olap [i] . t_lo);
      if  (matching_bases < mb || (matching_bases == mb
                                   && olap [i] . quality < olap [best] . quality))
        best = i;
    }

  for  (i = 0;  i < ct;  i ++)
    deleted [i] = (i != best);

  return;
}




//  Return  TRUE  iff there is any length-( window_len ) subarray
//  of  a [0 .. (n-1)]  that sums to  threshold  or higher.
static
int
Has_Bad_Window(char *a, int n, int window_len, int threshold) {

  if (n < window_len)
    return(FALSE);

  int32 sum = 0;
  int32 i=0;
  int32 j=0;

  for (i=0; i<window_len; i++)
    sum += a[i];

  if (sum >= threshold)
    return(TRUE);

  while (i < n) {
    sum -= a[j++];
    sum += a[i++];
    if (sum >= threshold)
      return(TRUE);
  }

  return(FALSE);
}



//  Find and report all overlaps and branch points between string  S
//  (with length  S_Len  and id  S_ID ) and string  T  (with
//  length & screen info in  t_info  and id  T_ID ) using the exact
//  matches in the list beginning at subscript  (* Start) .   Dir  is
//  the orientation of  S .

static
void
Process_Matches (int * Start,
                 char * S,
                 int S_Len,
                 char * S_quality,
                 AS_IID S_ID,
                 Direction_t Dir,
                 char * T,
                 Hash_Frag_Info_t t_info,
                 char * T_quality,
                 AS_IID T_ID,
                 Work_Area_t * WA,
                 int consistent) {
  int  P, * Ref;
  Olap_Info_t  distinct_olap [MAX_DISTINCT_OLAPS];
  Match_Node_t  * Longest_Match, * Ptr;
  Overlap_t  Kind_Of_Olap = NONE;
  double  Quality;
  int  Olap_Len;
  int  overlaps_output = 0;
  int  distinct_olap_ct;
  int  Max_Len, S_Lo, S_Hi, T_Lo, T_Hi;
  int  t_len;
  int  Done_S_Left, Done_S_Right;
  int  Errors;

  Done_S_Left = Done_S_Right = FALSE;
  t_len = t_info . length;

  assert ((* Start) != 0);

  // If a singleton match is hopeless on either side
  // it needn't be processed

  if  (Use_Hopeless_Check
       && WA -> Match_Node_Space [(* Start)] . Next == 0
       && ! Doing_Partial_Overlaps)
    {
      int  s_head, t_head, s_tail, t_tail;
      int  is_hopeless = FALSE;

      s_head = WA -> Match_Node_Space [(* Start)] . Start;
      t_head = WA -> Match_Node_Space [(* Start)] . Offset;
      if  (s_head <= t_head)
        {
          if  (s_head > HOPELESS_MATCH && ! WA -> left_end_screened)
            is_hopeless = TRUE;
        }
      else
        {
          if  (t_head > HOPELESS_MATCH  && ! t_info.lfrag_end_screened)
            is_hopeless = TRUE;
        }

      s_tail = S_Len - s_head - WA -> Match_Node_Space [(* Start)] . Len + 1;
      t_tail = t_len - t_head - WA -> Match_Node_Space [(* Start)] . Len + 1;
      if  (s_tail <= t_tail)
        {
          if  (s_tail > HOPELESS_MATCH && ! WA -> right_end_screened)
            is_hopeless = TRUE;
        }
      else
        {
          if  (t_tail > HOPELESS_MATCH && ! t_info.rfrag_end_screened)
            is_hopeless = TRUE;
        }

#ifdef COMPARE
      fprintf(stderr, "IS HOPELESS %d %d %d %d\n",
              WA->left_end_screened, WA->right_end_screened,
              t_info.lfrag_end_screened,
              t_info.rfrag_end_screened);
#endif

      if  (is_hopeless)
        {
          (* Start) = 0;
          WA->Kmer_Hits_Without_Olap_Ct ++;
          return;
        }
    }

  distinct_olap_ct = 0;

  while  ((* Start) != 0)
    {
      int  a_hang, b_hang;
      int  hit_limit = FALSE;

      Max_Len = WA -> Match_Node_Space [(* Start)] . Len;
      Longest_Match = WA -> Match_Node_Space + (* Start);
      for  (P = WA -> Match_Node_Space [(* Start)] . Next;  P != 0;
            P = WA -> Match_Node_Space [P] . Next)
        if  (WA -> Match_Node_Space [P] . Len > Max_Len)
          {
            Max_Len = WA -> Match_Node_Space [P] . Len;
            Longest_Match = WA -> Match_Node_Space + P;
          }

      a_hang = Longest_Match -> Start - Longest_Match -> Offset;
      b_hang = a_hang + S_Len - t_len;
      hit_limit =  ((WA -> A_Olaps_For_Frag >= Frag_Olap_Limit
                     && a_hang <= 0)
                    ||
                    (WA -> B_Olaps_For_Frag >= Frag_Olap_Limit
                     && b_hang <= 0));

      if  (! hit_limit)
        {
          Kind_Of_Olap = Extend_Alignment (Longest_Match, S, S_Len, T, t_len,
                                           & S_Lo, & S_Hi, & T_Lo, & T_Hi, & Errors, WA);

#ifdef COMPARE
          fprintf(stderr, "Extend_Alignment S: %d %d %d T %d %d %d kind %d\n",
                  S_Len, S_Lo, S_Hi, t_len, T_Lo, T_Hi, Kind_Of_Olap);
#endif

          if  (Kind_Of_Olap == DOVETAIL || Doing_Partial_Overlaps)
            {
              if  (1 + S_Hi - S_Lo >= Min_Olap_Len
                   && 1 + T_Hi - T_Lo >= Min_Olap_Len)
                {
                  Olap_Len = 1 + MIN (S_Hi - S_Lo, T_Hi - T_Lo);
                  Quality = (double) Errors / Olap_Len;

#ifdef COMPARE
                  fprintf(stderr, "Add_Overlap: scr_sub 0 erorrs %d bound %d passes 1\n",
                          Errors, WA->Error_Bound [Olap_Len]);
#endif

                  if  (Errors <= WA -> Error_Bound [Olap_Len])
                    {
                      Add_Overlap (S_Lo, S_Hi, T_Lo, T_Hi, Quality,
                                   distinct_olap, & distinct_olap_ct, WA);
                    }
                }
            }
        }

      if  (consistent)
        (* Start) = 0;

      for  (Ref = Start;  (* Ref) != 0;  )
        {
          Ptr = WA -> Match_Node_Space + (* Ref);
          if  (Ptr == Longest_Match
               || ((Kind_Of_Olap == DOVETAIL || Doing_Partial_Overlaps)
                   && S_Lo - SHIFT_SLACK <= Ptr -> Start
                   && Ptr -> Start + Ptr -> Len
                   <= (S_Hi + 1) + SHIFT_SLACK - 1
                   && Lies_On_Alignment
                   (Ptr -> Start, Ptr -> Offset,
                    S_Lo, T_Lo, WA)
                   ))
            (* Ref) = Ptr -> Next;         // Remove this node, it matches
          //   the alignment
          else
            Ref = & (Ptr -> Next);
        }
    }

#ifdef COMPARE
  fprintf(stderr, "distinct_olap_ct = %d\n", distinct_olap_ct);
#endif

  if  (distinct_olap_ct > 0)
    {
      int  deleted [MAX_DISTINCT_OLAPS] = {0};
      Olap_Info_t  * p;
      int  i;

      //  Check if any previously distinct overlaps should be merged because
      //  of other merges.

      if  (Doing_Partial_Overlaps)
        {
          if  (Unique_Olap_Per_Pair)
            Choose_Best_Partial (distinct_olap, distinct_olap_ct, deleted);
          //else
          //  Do nothing, output them all
        }
      else
        {
          if  (Unique_Olap_Per_Pair)
            Combine_Into_One_Olap (distinct_olap, distinct_olap_ct,
                                   deleted);
          else
            Merge_Intersecting_Olaps (distinct_olap, distinct_olap_ct,
                                      deleted);
        }

      p = distinct_olap;

      for  (i = 0;  i < distinct_olap_ct;  i ++) {
        if  (! deleted [i]) {
          bool rejected = FALSE;

          if  (Use_Window_Filter) {
            int32  d;
            int32  i = p -> s_lo;;
            int32  j = p -> t_lo;
            int32  q_len = 0;

#if AS_READ_MAX_NORMAL_LEN > 5000
            char  *q_diff = new char [AS_READ_MAX_NORMAL_LEN];
#else
            char   q_diff[AS_READ_MAX_NORMAL_LEN];
#endif
            for  (int32 k=0;  k<p->delta_ct;  k++) {
              int32 len = abs(p->delta[k]);

              for (int32 n = 1;  n < len;  n++) {
                if  (S[i] == T[j] || S[i] == DONT_KNOW_CHAR || T[j] == DONT_KNOW_CHAR) {
                  d = 0;
                } else {
                  d = MIN (S_quality[i], T_quality[j]);
                  d = MIN (d, QUALITY_CUTOFF);
                }

                q_diff[q_len++] = d;

                i++;
                j++;
              }

              if  (p -> delta[k] > 0) {
                d = S_quality[i];
                i++;
              } else {
                d = T_quality[j];
                j++;
              }
              q_diff[q_len++] = MIN (d, QUALITY_CUTOFF);
            }

            while  (i <= p -> s_hi) {
              if  (S[i] == T[j] || S[i] == DONT_KNOW_CHAR || T[j] == DONT_KNOW_CHAR) {
                d = 0;
              } else {
                d = MIN(S_quality[i], T_quality[j]);
                d = MIN(d, QUALITY_CUTOFF);
              }

              q_diff[q_len++] = d;

              i++;
              j++;
            }
            
            if (Has_Bad_Window(q_diff, q_len, BAD_WINDOW_LEN, BAD_WINDOW_VALUE)) {
              rejected = TRUE;
              Bad_Short_Window_Ct++;
            }

            else if (Has_Bad_Window(q_diff, q_len, 100, 240)) {
              rejected = TRUE;
              Bad_Long_Window_Ct++;
            }
            
#if AS_READ_MAX_NORMAL_LEN > 5000
            delete [] q_diff;
#endif
          }

          if  (! rejected) {
            if  (Doing_Partial_Overlaps)
              Output_Partial_Overlap(S_ID, T_ID, Dir, p, S_Len, t_len, WA);
            else
              Output_Overlap(S_ID, S_Len, Dir, T_ID, t_len, p, WA);

            overlaps_output++;

            if  (p->s_lo == 0)
              WA->A_Olaps_For_Frag++;

            if  (p->s_hi >= S_Len - 1)
              WA->B_Olaps_For_Frag++;
          }
        }
        p++;
      }
    }


  if  (overlaps_output == 0)
    WA->Kmer_Hits_Without_Olap_Ct++;
  else
    {
      WA->Kmer_Hits_With_Olap_Ct++;
      if  (overlaps_output > 1)
        WA->Multi_Overlap_Ct++;
    }


  return;
}





//  Compare the  diag_sum  fields  in  a  and  b  as  (String_Olap_t *) 's and
//  return  -1  if  a < b ,  0  if  a == b , and  1  if  a > b .
//  Used for  qsort .
static
int
By_Diag_Sum (const void * a, const void * b) {
  String_Olap_t  * x, * y;

  x = (String_Olap_t *) a;
  y = (String_Olap_t *) b;

  if  (x -> diag_sum < y -> diag_sum)
    return  -1;
  else if  (x -> diag_sum > y -> diag_sum)
    return  1;

  return  0;
}




//  Find and report all overlaps and branch points between string  S
//  and all strings in global  String_Olap_Space .
//  Return the number of entries processed.
//  Len  is the length of  S ,  ID  is its fragment ID  and
//  Dir  indicates if  S  is forward, or reverse-complemented.
int
Process_String_Olaps (char * S,
                      int Len,
                      char * S_quality,
                      AS_IID ID,
                      Direction_t Dir,
                      Work_Area_t * WA) {
  int32  i, ct, root_num, start, processed_ct;

  //  Move all full entries to front of String_Olap_Space and set
  //  diag_sum to average diagonal.  if enough entries to bother,
  //  sort by average diagonal.  Then process positive & negative diagonals
  //  separately in order from longest to shortest overlap.  Stop
  //  processing when output limit has been reached.

  for  (i = ct = 0;  i < WA -> Next_Avail_String_Olap;  i ++)
    if  (WA -> String_Olap_Space [i] . Full) {
      root_num = WA -> String_Olap_Space [i] . String_Num;
      if  (root_num + Hash_String_Num_Offset > ID) {
        if  (WA -> String_Olap_Space [i] . Match_List == 0) {
          fprintf (stderr, " Curr_String_Num = %d  root_num  %d have no matches\n",
                   ID, root_num);
          exit (-2);
        }
        if  (i != ct)
          WA -> String_Olap_Space [ct] = WA -> String_Olap_Space [i];
        assert (WA -> String_Olap_Space [ct] . diag_ct > 0);
        WA -> String_Olap_Space [ct] . diag_sum
          /= WA -> String_Olap_Space [ct] . diag_ct;
        ct ++;
      }
    }

  if  (ct == 0)
    return  ct;

  if  (ct <= Frag_Olap_Limit) {
    for  (i = 0;  i < ct;  i ++) {
      root_num = WA -> String_Olap_Space [i] . String_Num;

      WA -> Edit_Match_Limit = Read_Edit_Match_Limit;
      WA -> Error_Bound = Read_Error_Bound;

      Process_Matches
        (& (WA -> String_Olap_Space [i] . Match_List),
         S, Len, S_quality, ID, Dir,
         Data + String_Start [root_num],
         String_Info [root_num],
         Quality_Data + String_Start [root_num],
         root_num + Hash_String_Num_Offset, WA,
         WA -> String_Olap_Space [i] . consistent);
      assert (WA -> String_Olap_Space [i] . Match_List == 0);
    }

    return  ct;
  }

#if  1
  qsort (WA -> String_Olap_Space, ct, sizeof (String_Olap_t), By_Diag_Sum);
#else
  for  (i = 0;  i < ct - 1;  i ++)
    for  (j = i + 1;  j < ct;  j ++)
      if  (WA -> String_Olap_Space [i] . diag_sum
           > WA -> String_Olap_Space [j] . diag_sum)
        {
          String_Olap_t  save = WA -> String_Olap_Space [i];

          WA -> String_Olap_Space [i] = WA -> String_Olap_Space [j];
          WA -> String_Olap_Space [j] = save;
        }
#endif

  for  (start = 0;
        start < ct && WA -> String_Olap_Space [start] . diag_sum < 0;
        start ++)
    ;

  processed_ct = 0;

  for  (i = start;  i < ct && WA -> A_Olaps_For_Frag < Frag_Olap_Limit ;  i ++) {
    root_num = WA -> String_Olap_Space [i] . String_Num;

    WA -> Edit_Match_Limit = Read_Edit_Match_Limit;
    WA -> Error_Bound = Read_Error_Bound;

    Process_Matches
      (& (WA -> String_Olap_Space [i] . Match_List),
       S, Len, S_quality, ID, Dir,
       Data + String_Start [root_num],
       String_Info [root_num],
       Quality_Data + String_Start [root_num],
       root_num + Hash_String_Num_Offset, WA,
       WA -> String_Olap_Space [i] . consistent);
    assert (WA -> String_Olap_Space [i] . Match_List == 0);

    processed_ct ++;
  }
  for  (i = start - 1;  i >= 0 && WA -> B_Olaps_For_Frag < Frag_Olap_Limit ;  i --) {
    root_num = WA -> String_Olap_Space [i] . String_Num;

    WA -> Edit_Match_Limit = Read_Edit_Match_Limit;
    WA -> Error_Bound = Read_Error_Bound;

    Process_Matches
      (& (WA -> String_Olap_Space [i] . Match_List),
       S, Len, S_quality, ID, Dir,
       Data + String_Start [root_num],
       String_Info [root_num],
       Quality_Data + String_Start [root_num],
       root_num + Hash_String_Num_Offset, WA,
       WA -> String_Olap_Space [i] . consistent);
    assert (WA -> String_Olap_Space [i] . Match_List == 0);

    processed_ct ++;
  }

  return  processed_ct;
}


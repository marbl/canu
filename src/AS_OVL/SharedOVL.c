
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
* Module:  SharedOVL.c
* Description:
*   Functions shared by overlap programs
* 
*    Programmer:  A. Delcher
*       Started:   15 February 2007
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: SharedOVL.c,v 1.3 2007-08-02 21:23:52 adelcher Exp $
 * $Revision: 1.3 $
*/


#include  "SharedOVL.h"


int  Fwd_Homopoly_Prefix_Match
  (const char * A, int m, const char * T, int n, int error_limit,
   int * a_end, int * t_end, int * match_to_end, int * delta,
   int * delta_len, Homopoly_Match_Entry_t * * edit_array)

// Return the minimum number of non-homopoly insertions and substitutions
// needed to match a prefix of string  A  of length  m  with a
// prefix of string  T  of length  n , extending to the end of
// one of those strings if possible without exceeding  error_limit
// differences.  Homopolymer-type insertions/deletions do not count.
// Set  a_end  and  b_end  to where the match ends (each is the number
// of characters aligned in the respective string).  Set
//  match_to_end  true iff the alignment reached the end.
// Return the delta-coded alignment in  delta  with  delta_len  set
// to the number of entries in  delta .  Use  edit_array  to store
// intermediate values in computing the alignment (in the Vishkin-Schieber
// style of a pyrimidal array).

  {
   int  score, shorter;
   int  row, left, right;
   int  d, e, j, k, s;

   (* delta_len) = 0;

   shorter = OVL_Min_int (m, n);
   for (row = 0; row < shorter && A [row] == T [row]; row ++)
     ;
   edit_array [0] [0] . len = row;
   edit_array [0] [0] . score = 0;
   edit_array [0] [0] . at_end = (row == shorter);

   if (row == shorter)     // Exact match to end
     {
      (* a_end) = (* t_end) = row;
      (* match_to_end) = TRUE;
      return 0;
     }

   left = right = 0;
   for (e = 1; e <= error_limit; e ++)
     {
      left = OVL_Max_int (left - 1, -e);
      right = OVL_Min_int (right + 1, e);

      for (d = left; d <= right; d ++)
        {
         // first find the best score in the preceding row
         row = 0;
         score = INT_MAX;

         // check cell directly above
         if (left < row && row < right && ! edit_array [e - 1] [d] . at_end)
           {
            row = 1 + edit_array [e - 1] [d] . len;
            score = edit_array [e - 1] [d] . score + HP_SUBST_SCORE;
           }

         // check upper left cell
         if (left + 1 < row && ! edit_array [e - 1] [d - 1] . at_end)
           {
            k = edit_array [e - 1] [d - 1] . len + d - 1;
            if (k > 0 && T [k - 1] == T [k])
              s = edit_array [e - 1] [d - 1] . score + HP_INDEL_SCORE;
                  // homopoly insert in T
            else
              s = edit_array [e - 1] [d - 1] . score + NON_HP_INDEL_SCORE;
                  // non-homopoly insert in T
            if (s < score)   // low scores are better
              {
               row = edit_array [e - 1] [d - 1] . len;
               score = s;
              }
           }

         // check upper right cell
         if (row < right - 1 && ! edit_array [e - 1] [d + 1] . at_end)
           {
            k = edit_array [e - 1] [d + 1] . len;
            if (k > 0 && A [k - 1] == A [k])
              s = edit_array [e - 1] [d + 1] . score + HP_INDEL_SCORE;
                   // homopoly insert in A
            else
              s = edit_array [e - 1] [d + 1] . score + NON_HP_INDEL_SCORE;
                   // non-homopoly insert in A
            if (s < score)   // low scores are better
              {
               row = edit_array [e - 1] [d - 1] . len + 1;
               score = s;
              }
           }

         // now extend as far as matches allow
         while (row < m && row + d < n && A [row] == T [row + d])
            row ++;

         edit_array [e] [d] . len = row;
         edit_array [e] [d] . score = score;
         edit_array [e] [d] . at_end = (row == m || row + d == n);
        }

      //**ALD  left off here
      // stop and return best at_end score if all other scores are
      // worse than it.

      // shrink left..right range here based on hopelessly bad scores

      assert (left <= right);
     }

   // For failed alignments just return the exact match at the beginning
   (* a_end) = edit_array [0] [0] . len;
   (* t_end) = edit_array [0] [0] . len;
//**ALD  need new version of this
//   Set_Fwd_Delta (delta, delta_len, edit_array, 0, 0);
   (* match_to_end) = FALSE;

   return 0;
  }



int  Fwd_Prefix_Edit_Dist
  (char a_string [], int m, char t_string [], int n, int error_limit,
   int * a_end, int * t_end, int * match_to_end,
   double match_value, int * delta, int * delta_len, int * * edit_array,
   int edit_match_limit [], int error_bound [], int doing_partial)

// Return the minimum number of changes (inserts, deletes, replacements)
// needed to match prefixes of strings  a_string [0 .. (m-1)]  and
//  t_string [0 .. (n-1)]  if it's not more than  error_limit , where
// the match must extend to the end of one of those strings.
// Put delta description of alignment in  delta  and set
// (* delta_len)  to the number of entries there if it's a complete
// match.   match_value  is the score for each matching character
// and is normalized so that  1.0 - match_value  is the score for
// a mismatch (including indels).
// Set  a_end  and  t_end  to the rightmost positions where the
// alignment ended in  a_string  and  t_string , respectively.
// Set  match_to_end  true if the match extended to the end
// of at least one string; otherwise, set it false to indicate
// a branch point.
//  edit_array  has storage preallocated storage that can be used
// for this computation (including in a threaded environment).
//  edit_match_limit [e]  is the minimum match length worth attempting to
// extend containing  e  errors.   error_bound [i]  has the most errors
// that can be tolerated in a match of length  i .  If  doing_partial
// is true return the best match, whether it extends to the end of either
// string or not.

  {
   double  score, max_score;
   int  max_score_len, max_score_best_d, max_score_best_e;
   int  best_d, best_e, longest, row, tail_len;
   int  left, right;
   int  d, e, j, shorter;

   best_d = best_e = longest = 0;
   (* delta_len) = 0;

   shorter = OVL_Min_int (m, n);
   for (row = 0; row < shorter && a_string [row] == t_string [row]; row ++)
      ;

   edit_array [0] [0] = row;

   if (row == shorter)                              // Exact match
     {
      (* a_end) = (* t_end) = row;
      (* match_to_end) = TRUE;
      return 0;
     }

   left = right = 0;
   max_score = 0.0;
   max_score_len = max_score_best_d = max_score_best_e = 0;
   for (e = 1; e <= error_limit; e ++)
     {
      left = OVL_Max_int (left - 1, -e);
      right = OVL_Min_int (right + 1, e);
      edit_array [e - 1] [left] = -2;
      edit_array [e - 1] [left - 1] = -2;
      edit_array [e - 1] [right] = -2;
      edit_array [e - 1] [right + 1] = -2;

      for (d = left; d <= right; d ++)
        {
         row = 1 + edit_array [e - 1] [d];
         if ((j = edit_array [e - 1] [d - 1]) > row)
            row = j;
         if ((j = 1 + edit_array [e - 1] [d + 1]) > row)
            row = j;
         while (row < m && row + d < n
              && a_string [row] == t_string [row + d])
            row ++;

         edit_array [e] [d] = row;

         if (row == m || row + d == n)
           {
            // Force last error to be mismatch rather than insertion
            if (row == m && 1 + edit_array [e - 1] [d + 1] == edit_array [e] [d]
                 && d < right)
              {
               d ++;
               edit_array [e] [d] = edit_array [e] [d - 1];
              }
            (* a_end) = row;           // One past last align position
            (* t_end) = row + d;

#if 1
            //  Check for branch point here caused by uneven
            //  distribution of errors
            score = row * match_value - e;
                 // Assumes  match_value - mismatch_value == 1.0
            tail_len = row - max_score_len;
            if ((doing_partial && score < max_score)
                 || (e > MIN_BRANCH_END_DIST / 2
                    && tail_len >= MIN_BRANCH_END_DIST
                    && (max_score - score) / tail_len >= MIN_BRANCH_TAIL_SLOPE))
              {
               (* a_end) = max_score_len;
               (* t_end) = max_score_len + max_score_best_d;
               Set_Fwd_Delta (delta, delta_len, edit_array,
                    max_score_best_e, max_score_best_d);
               (* match_to_end) = FALSE;
               return max_score_best_e;
              }
#endif
            Set_Fwd_Delta (delta, delta_len, edit_array, e, d);
            (* match_to_end) = TRUE;
            return e;
           }
        }

      while (left <= right && left < 0
           && edit_array [e] [left] < edit_match_limit [e])
        left ++;
      if (left >= 0)
         while (left <= right
              && edit_array [e] [left] + left < edit_match_limit [e])
           left ++;
      if (left > right)
          break;
      while (right > 0
           && edit_array [e] [right] + right < edit_match_limit [e])
        right --;
      if (right <= 0)
          while (edit_array [e] [right] < edit_match_limit [e])
            right --;
      assert (left <= right);

      for (d = left; d <= right; d ++)
         if (edit_array [e] [d] > longest)
           {
            best_d = d;
            best_e = e;
            longest = edit_array [e] [d];
           }

      score = longest * match_value - e;
               // Assumes  match_value - mismatch_value == 1.0
      if (score > max_score
               && best_e <= error_bound [OVL_Min_int (longest, longest + best_d)])
        {
         max_score = score;
         max_score_len = longest;
         max_score_best_d = best_d;
         max_score_best_e = best_e;
        }
     }

   (* a_end) = max_score_len;
   (* t_end) = max_score_len + max_score_best_d;
   Set_Fwd_Delta (delta, delta_len, edit_array, max_score_best_e,
        max_score_best_d);
   (* match_to_end) = FALSE;

   return max_score_best_e;
  }


int  OVL_Max_int
  (int a, int b)

// Return the larger of  a  and  b .

  {
   if (a < b)
      return b;

   return a;
  }


int  OVL_Min_int
  (int a, int b)

// Return the smaller of  a  and  b .

  {
   if (a < b)
      return a;

   return b;
  }


int  Rev_Prefix_Edit_Dist
  (char a_string [], int m, char t_string [], int n, int error_limit,
   int * a_end, int * t_end, int * leftover, int * match_to_end,
   double match_value, int * delta, int * delta_len, int * * edit_array,
   int edit_match_limit [], int error_bound [], int doing_partial)

// Return the minimum number of changes (inserts, deletes, replacements)
// needed to match a prefix of string  a_string [0 .. (1-m)]  right-to-left
// with a prefix of string  t_string [0 .. (1-n)] , also right-to-left,
// if it's not more than  error_limit .  Note the subscripts are negative.
// The match must extend to the beginning of one of the two strings.
// Put the delta description of the (forward) alignment in  delta  and set
// (* delta_len)  to the number of entries there if it's a complete
// match.  Set  leftover  to the number of characters that match after
// the last  delta  entry.  match_value  is the score for each matching
// character  and is normalized so that  1.0 - match_value  is the score
// for a mismatch (including indels).
// Set  a_end  and  t_end  to the leftmost positions where the
// alignment ended in  a_string  and  t_string , respectively.
// Set  match_to_end  true if the match extended to the end
// of at least one string; otherwise, set it false to indicate
// a branch point.
//  edit_array  has storage preallocated storage that can be used
// for this computation (including in a threaded environment).
//  edit_match_limit [e]  is the minimum match length worth attempting to
// extend containing  e  errors.   error_bound [i]  has the most errors
// that can be tolerated in a match of length  i .  If  doing_partial
// is true return the best match, whether it extends to the end of either
// string or not.

  {
   double  score, max_score;
   int  max_score_len, max_score_best_d, max_score_best_e;
   int  best_d, best_e, longest, row, tail_len;
   int  left, right;
   int  d, e, j, shorter;

   best_d = best_e = longest = 0;
   (* delta_len) = 0;

   shorter = OVL_Min_int (m, n);
   for (row = 0; row < shorter && a_string [- row] == t_string [- row]; row ++)
      ;    //**ALD maybe should allow dont_know characters here?

   edit_array [0] [0] = row;

   if (row == shorter)                              // Exact match
     {
      (* a_end) = (* t_end) = - row;
      (* leftover) = shorter;
      (* match_to_end) = TRUE;
      return 0;
     }

   left = right = 0;
   max_score = 0.0;
   max_score_len = max_score_best_d = max_score_best_e = 0;
   for (e = 1; e <= error_limit; e ++)
     {
      left = OVL_Max_int (left - 1, -e);
      right = OVL_Min_int (right + 1, e);
      edit_array [e - 1] [left] = -2;
      edit_array [e - 1] [left - 1] = -2;
      edit_array [e - 1] [right] = -2;
      edit_array [e - 1] [right + 1] = -2;

      for (d = left; d <= right; d ++)
        {
         row = 1 + edit_array [e - 1] [d];
         if ((j = edit_array [e - 1] [d - 1]) > row)
            row = j;
         if ((j = 1 + edit_array [e - 1] [d + 1]) > row)
            row = j;
         while (row < m && row + d < n
              && a_string [- row] == t_string [- row - d])
            row ++;     //**ALD  dont_know characters??

         edit_array [e] [d] = row;

         if (row == m || row + d == n)
           {
            // Force last error to be mismatch rather than insertion
            if (row == m && 1 + edit_array [e - 1] [d + 1] == edit_array [e] [d]
                 && d < right)
              {
               d ++;
               edit_array [e] [d] = edit_array [e] [d - 1];
              }
            (* a_end) = - row;           // One past last align position
            (* t_end) = - row - d;

            //  Check for branch point here caused by uneven
            //  distribution of errors
            score = row * match_value - e;
                      // Assumes  match_value - mismatch_value == 1.0
            tail_len = row - max_score_len;

            if ((doing_partial && score < max_score) ||
                 (e > MIN_BRANCH_END_DIST / 2
                   && tail_len >= MIN_BRANCH_END_DIST
                   && (max_score - score) / tail_len >= MIN_BRANCH_TAIL_SLOPE))
              {
               (* a_end) = - max_score_len;
               (* t_end) = - max_score_len - max_score_best_d;
               Set_Rev_Delta (delta, delta_len, edit_array, max_score_best_e,
                    max_score_best_d, leftover, t_end, n);
               (* match_to_end) = FALSE;
               return max_score_best_e;
              }

            Set_Rev_Delta (delta, delta_len, edit_array, e, d,
                 leftover, t_end, n);
            (* match_to_end) = TRUE;
            return e;
           }
        }

      while (left <= right && left < 0
           && edit_array [e] [left] < edit_match_limit [e])
        left ++;
      if (left >= 0)
         while (left <= right
              && edit_array [e] [left] + left < edit_match_limit [e])
           left ++;
      if (left > right)
          break;
      while (right > 0
           && edit_array [e] [right] + right < edit_match_limit [e])
        right --;
      if (right <= 0)
         while (edit_array [e] [right] < edit_match_limit [e])
            right --;
      assert (left <= right);

      for (d = left; d <= right; d ++)
         if (edit_array [e] [d] > longest)
           {
            best_d = d;
            best_e = e;
            longest = edit_array [e] [d];
           }

      score = longest * match_value - e;
               // Assumes  match_value - mismatch_value == 1.0
      if (score > max_score
               && best_e <= error_bound [OVL_Min_int (longest, longest + best_d)])
        {
         max_score = score;
         max_score_len = longest;
         max_score_best_d = best_d;
         max_score_best_e = best_e;
        }
     }

   (* a_end) = - max_score_len;
   (* t_end) = - max_score_len - max_score_best_d;
   Set_Rev_Delta (delta, delta_len, edit_array, max_score_best_e,
        max_score_best_d, leftover, t_end, n);
   (* match_to_end) = FALSE;

   return max_score_best_e;
  }


void  Set_Fwd_Delta
  (int delta [], int * delta_len, int * * edit_array,
   int e, int d)

// Set  delta  to the entries indicating the insertions/deletions
// in the alignment encoded in  edit_array  ending at position
//  edit_array [e] [d] .  This is the position in the first
// string where the alignment ended.  Set  (* delta_len)  to
// the number of entries in  delta .

  {
   int  d_len, from, last, max;
   int  i, j, k;

   last = edit_array [e] [d];
   d_len = 0;

   // temporarily put values in  delta  starting at subscript  e
   // and working down.

   for (k = e; k > 0; k --)
     {
      // Find best entry on prior row  k - 1
      from = d;
      max = 1 + edit_array [k - 1] [d];
      if ((j = edit_array [k - 1] [d - 1]) > max)
        {
         from = d - 1;
         max = j;
        }
      if ((j = 1 + edit_array [k - 1] [d + 1]) > max)
        {
         from = d + 1;
         max = j;
        }

      // Save the  appropriate delta value
      if (from == d - 1)
        {
         delta [e - d_len] = max - last - 1;
         d_len ++;
         d --;
         last = edit_array [k - 1] [from];
        }
      else if (from == d + 1)
        {
         delta [e - d_len] = last - (max - 1);
         d_len ++;
         d ++;
         last = edit_array [k - 1] [from];
        }
     }
   delta [e - d_len] = last + 1;

   // Now shift values to the front of  delta
   k = 0;
   for (i = e - d_len; i < e; i ++)
     delta [k ++] = abs (delta [i]) * Sign (delta [i + 1]);

   (* delta_len) = d_len;

   return;
  }


void  Set_Rev_Delta
  (int delta [], int * delta_len, int * * edit_array,
   int e, int d, int * leftover, int * t_end, int t_len)

// Set  delta  to the entries indicating the insertions/deletions
// in the alignment encoded in  edit_array  ending at position
//  edit_array [e] [d] .  This is the position in the first
// string where the alignment ended.  This alignment was in the
// reverse direction of the original strings, so the delta-encoding
// returned will be in the forward direction of those strings.
// Set  (* delta_len)  to  the number of entries in  delta .
// Set  (* leftover)  to the number of characters that match after
// the last  delta  entry.  Don't allow the first delta entry to be
// an insertion in the T string if it can be converted to a substitution
// by adjusting  (* t_end)  which is the negative (because the alignment
// was in the reverse direction) value where the alignment ended in the
// T string, which has length  t_len .

  {
   int  d_len, from, last, max;
   int  j, k;

   last = edit_array [e] [d];
   d_len = 0;

   // put values in  delta  tracing the alignment backwards

   for (k = e; k > 0; k --)
     {
      // Find best entry on prior row  k - 1
      from = d;
      max = 1 + edit_array [k - 1] [d];
      if ((j = edit_array [k - 1] [d - 1]) > max)
        {
         from = d - 1;
         max = j;
        }
      if ((j = 1 + edit_array [k - 1] [d + 1]) > max)
        {
         from = d + 1;
         max = j;
        }

      // Save the  appropriate delta value
      if (from == d - 1)
        {
         delta [d_len ++] = max - last - 1;
         d --;
         last = edit_array [k - 1] [from];
        }
      else if (from == d + 1)
        {
         delta [d_len ++] = last - (max - 1);
         d ++;
         last = edit_array [k - 1] [from];
        }
     }
   (* leftover) = last;

   // Don't allow first delta to be +1 or -1
   assert (d_len == 0 || delta [0] != -1);
   if (d_len > 0 && delta [0] == 1 && (* t_end) + t_len > 0)
     {  //**ALD I think this case is prevented in the calling routine
      int  i;

      if (delta [1] > 0)
         delta [0] = delta [1] + 1;
      else
         delta [0] = delta [1] - 1;
      for (i = 2;  i < d_len;  i ++)
        delta [i - 1] = delta [i];
      d_len --;
      (* t_end) --;
      if (d_len == 0)
         (* leftover) ++;
     }

   (* delta_len) = d_len;

   return;
  }


int  Sign
  (int a)

// Return the algebraic sign of  a .

  {
   if (a > 0)
     return  1;
   else if (a < 0)
     return  -1;

   return  0;
  }




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
 * $Id: SharedOVL.c,v 1.4 2007-08-23 15:09:06 adelcher Exp $
 * $Revision: 1.4 $
*/


#include  "SharedOVL.h"


//**ALD
#define  MAX_SEQUENCE_LEN  2048

int  Fwd_Banded_Homopoly_Prefix_Match
  (const char * AA, int m, const char * TT, int n, int score_limit,
   int * return_score, int * a_end, int * t_end, int * match_to_end, int * delta,
   int * delta_len, Alignment_Cell_t edit_space [])

// Return the number of homopoly-type changes (indels and substitutions)
// in the best match a prefix of string  AA  of length  m  with a
// prefix of string  TT  of length  n , extending to the end of
// one of those strings if possible without exceeding  score_limit
// as the score of the alignment.  The return value is the number of
// changes; set  return_score  to the score of the best alignment.
// Set  a_end  and  b_end  to where the match ends (each is the number
// of characters aligned in the respective string).  Set
//  match_to_end  true iff the alignment reached the end.
// Return the delta-coded alignment in  delta  with  delta_len  set
// to the number of entries in  delta .  Use array  edit_space  to store
// intermediate values in computing the alignment (as a simple banded
// alignment).

  {
   Alignment_Cell_t  * curr, * prev;
   const char  * A, * T;
   int  indent [MAX_SEQUENCE_LEN], width [MAX_SEQUENCE_LEN];
   int  best_row, best_col, best_score, best_errors;
   int  cutoff_score, global_best, last_row = 0;
   int  te, ts;  // temporary values for errors and score
   int  i, j;

   // offset subscripts to start at 1
   A = AA - 1; 
   T = TT - 1;

   global_best = INT_MAX;

   // fill in the first row/band
   indent [0] = 0;
   curr = edit_space;
   ts = best_score = curr [0] . score = curr [0] . errors = 0;
   cutoff_score = best_score + BAND_SCORE_DELTA;
   for (j = 1; ts < cutoff_score && j <= n; j ++)
     {
      curr [j] . from = -1;   // left
      ts = curr [j - 1] . score;
      if (1 < j && T [j] == T [j - 1])
         ts += HP_INDEL_SCORE;   // homopoly insert in T
      else
        ts += NON_HP_INDEL_SCORE;   // non-homopoly insert in T
      curr [j] . score = ts;
      curr [j] . errors = j;
      // best score is automatically the first score on the first row
      if (j == n)
        {
         global_best = ts;
         best_row = 0;
         best_col = j;
         best_errors = j;
        }
     }
   width [0] = j;  // width of this row of the alignment

   for (i = 1; i <= m; i ++)
     {
      int  from, offset, shift, score, errors, min_width;

      prev = curr;
      curr += width [i - 1];

      // determine how many cells to shift the left end of this row
      // from the left end of the prior row
      // can't shift further than previous best
      for (shift = 0; cutoff_score < prev [shift] . score; shift ++)
        ;
      offset = indent [i] = indent [i - 1] + shift;

      // set values in left-most cell
      from = 1;   // top--is always a cell directly above the left-most cell
      score = prev [shift] . score;
      if (1 < i && A [i] == A [i - 1] && (0 < offset && A [i] == T [offset]
            || prev [shift] . from == 1))
        score += HP_INDEL_SCORE;
      else
        score += NON_HP_INDEL_SCORE;
      errors = prev [shift] . errors + 1;
      if (0 < shift)
        {  // ok to check top left
         ts = prev [0] . score;
         te = prev [0] . errors;
         if (A [i] != T [offset])
           {
            ts += HP_SUBST_SCORE;
            te ++;
           }
         if (ts < score)
           {
            score = ts;
            errors = te;
            from = 0;
           }
        }
      curr [0] . score = score;
      curr [0] . errors = errors;
      curr [0] . from = from;

      // just in case this is the best;  j=0 here
      if ((i == m || indent [i] == n) && score < global_best)
        {
         best_row = i;
         best_col = indent [i];
         best_errors = errors;
         global_best = score;
        }

      min_width = width [i - 1] - shift;   // at least this many cells on current row
      cutoff_score = best_score + BAND_SCORE_DELTA;   // determined from prior row
      best_score = score;

      // set remaining cells in this row/band
      for (j = 1; indent [i] + j <= n && (j < min_width || score < cutoff_score); j ++)
        {
         // start with score from left
         from = -1;
         score = curr [j - 1] . score;
         if (1 < offset + j && T [offset + j] == T [offset + j - 1]
               && (T [offset + j] == A [i] || curr [j - 1] . from == -1))
           score += HP_INDEL_SCORE;
         else
           score += NON_HP_INDEL_SCORE;
         errors = curr [j - 1] . errors + 1;

         // see if top score is better
         if (j < min_width)
           {
            ts = prev [j + shift] . score;
            if (1 < i && A [i] == A [i - 1]
                  && (A [i] == T [j + offset] || prev [j + shift] . from == 1))
              ts += HP_INDEL_SCORE;
            else
              ts += NON_HP_INDEL_SCORE;
            te = prev [j + shift] . errors + 1;
            if (ts < score)
              {
               score = ts;
               errors = te;
               from = 1;
              }
           }

         // see if top-left score is better (if it's available)
         if (j <= min_width)
           {
            ts = prev [j + shift - 1] . score;
            te = prev [j + shift - 1] . errors;
            if (A [i] != T [j + offset])
              {
               ts += HP_SUBST_SCORE;
               te ++;
              }
            if (ts < score)
              {
               score = ts;
               errors = te;
               from = 0;
              }
           }

         curr [j] . score = score;
         curr [j] . errors = errors;
         curr [j] . from = from;
         if (score < best_score)
           best_score = score;
         if ((i == m || indent [i] + j == n) && score < global_best)
           {
            best_row = i;
            best_col = indent [i] + j;
            best_errors = errors;
            global_best = score;
           }
        }
      width [i] = j;
      last_row = i;

      // stop if best score on the row is no better than the best end
      // alignment score
      if (global_best <= best_score)
        break;
     }

   if (Verbose_Level > 2)
     {
      curr = edit_space;
      for (i = 0; i <= last_row; i ++)
        {
         printf ("%3d:  %c  %3d ", i, (0 < i ? A [i] : '-'), indent [i]);
         for (j = 0; j < width [i]; j ++)
           printf (" %3d%+2d%c", curr [j] . score, curr [j] . from,
                (0 < j + indent [i] && j + indent [i] <= n ? T [j + indent [i]] : '-'));
         putchar ('\n');
         curr += width [i];
        }
     }

   if (score_limit < global_best)
     {
      (* a_end) = (* t_end) = 0;
      (* match_to_end) = FALSE;
      (* return_score) = global_best;
      (* delta_len) = 0;
      return  0;
     }

   Set_Fwd_Banded_Delta (delta, delta_len, edit_space, best_row, best_col,
        width, indent);

   (* a_end) = best_row;
   (* t_end) = best_col;
   (* match_to_end) = TRUE;
   (* return_score) = global_best;

   return  best_errors;
  }


int  Fwd_Homopoly_Prefix_Match
  (const char * A, int m, const char * T, int n, int score_limit,
   int * return_score, int * a_end, int * t_end, int * match_to_end, int * delta,
   int * delta_len, Homopoly_Match_Entry_t ** edit_array)

// Return the number of homopoly-type changes (indels and substitutions)
// in the best match a prefix of string  A  of length  m  with a
// prefix of string  T  of length  n , extending to the end of
// one of those strings if possible without exceeding  score_limit
// as the score of the alignment.  The return value is the number of
// changes; set  return_score  to the score of the best alignment.
// Set  a_end  and  b_end  to where the match ends (each is the number
// of characters aligned in the respective string).  Set
//  match_to_end  true iff the alignment reached the end.
// Return the delta-coded alignment in  delta  with  delta_len  set
// to the number of entries in  delta .  Use  edit_array  to store
// intermediate values in computing the alignment (in the Vishkin-Schieber
// style of a pyramidal array).

  {
   Homopoly_Match_Entry_t  * ep;
   int  best_end_d, best_end_e, best_end_score;
   double  ratio, best_ratio, best_end_ratio;
   int  from, lowest_score, score, shorter;
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
      (* return_score) = 0;
      return 0;
     }

   lowest_score = 0;
   best_end_score = MAX_HOMOPOLY_SCORE;
   best_end_ratio = MAX_HOMOPOLY_SCORE;
   left = right = 0;
   for (e = 1; lowest_score < score_limit && e < MAX_ERRORS; e ++)
     {
      left = OVL_Max_int (left - 1, -e);
      right = OVL_Min_int (right + 1, e);
      lowest_score = MAX_HOMOPOLY_SCORE;
      best_ratio = DBL_MAX;

      for (d = left; d <= right; d ++)
        {
         // first find the best score in the preceding row
         row = 0;
         score = MAX_HOMOPOLY_SCORE;
         from = 2;

//printf ("e=%d  d=%d\n", e, d);
         // check cell directly above
         ep = edit_array [e - 1] + d;
         if (left < d && d < right && ! ep -> at_end)
           {
//printf ("check above\n");
            row = 1 + ep -> len;
            score = ep -> score + HP_SUBST_SCORE;
            from = 0;
           }

         // check upper left cell
         ep = edit_array [e - 1] + d - 1;
         if (left + 1 < d && ! ep -> at_end)
           {
//printf ("check upper left\n");
            j = ep -> len - 1;;
            k = j + d;
            if (k > 0 && T [k - 1] == T [k] && T [k] == A [j])
              s = ep -> score + HP_INDEL_SCORE;
                  // homopoly insert in T
            else
              s = ep -> score + NON_HP_INDEL_SCORE;
                  // non-homopoly insert in T
            if (s < score || (s == score && row < ep -> len))   // low scores are better
              {
//printf ("is better\n");
               row = ep -> len;
               score = s;
               from = -1;
              }
           }

         // check upper right cell
         ep = edit_array [e - 1] + d + 1;
         if (d < right - 1 && ! ep -> at_end)
           {
//printf ("check upper right\n");
            j = ep -> len;
            k = j + d;
            if (j > 0 && A [j - 1] == A [j] && A [j] == T [k])
              s = ep -> score + HP_INDEL_SCORE;
                   // homopoly insert in A
            else
              s = ep -> score + NON_HP_INDEL_SCORE;
                   // non-homopoly insert in A
            if (s < score || (s == score && row < ep -> len))   // low scores are better
              {
//printf ("is better\n");
               row = ep -> len + 1;
               score = s;
               from = +1;
              }
           }

         // now extend as far as matches allow
         while (row < m && row + d < n && A [row] == T [row + d])
            row ++;

         edit_array [e] [d] . len = row;
         edit_array [e] [d] . score = score;
         edit_array [e] [d] . at_end = (row == m || row + d == n);
         edit_array [e] [d] . from = from;

         if (score == MAX_HOMOPOLY_SCORE)
           {  // there was no legal value on the prior row so
              // set the at_end of this cell true so it won't be
              // used by any cells below it
            edit_array [e] [d] . at_end = TRUE;
            continue;
           }

         if (MIN_RATIO_E <= e)
           {
            ratio = (double) score / (row + EPSILON);  // add EPSION in case row=0
            if (ratio < best_ratio)
              best_ratio = ratio;
           }

         if (edit_array [e] [d] . at_end && score < best_end_score)
           {
            best_end_score = score;
            best_end_e = e;
            best_end_d = d;
            best_end_ratio = (double) score / (edit_array [e] [d] . len + EPSILON);
            score_limit = best_end_score;
           }
         if (score < lowest_score)
           lowest_score = score;
        }

//**ALD
if (0 && best_end_score < MAX_HOMOPOLY_SCORE)
  {
   printf ("best_end_e=%d  best_end_d=%d  best_end_score=%d\n",
        best_end_e, best_end_d, best_end_score);
   Show_Homopoly_Match_Array (stdout, edit_array, e);
  }

      // if lowest score on row is worse (higher) than the best
      // match_to_end score, then have found the best match
      if (best_end_score <= lowest_score)
        {
         Homopoly_Match_Entry_t  * bp = edit_array [best_end_e] + best_end_d;

         (* a_end) = bp -> len;
         (* t_end) = bp -> len + best_end_d;
         Set_Fwd_Homopoly_Delta (delta, delta_len, edit_array, best_end_e, best_end_d);
         (* match_to_end) = TRUE;
         (* return_score) = best_end_score;
        
         return  best_end_e;
        }

      //**ALD  left off here
      // shrink left..right range here based on hopelessly bad scores
      if (MIN_RATIO_E < e)
        {
         // trim cells whose score/len ratio is worse than 1.5 times the min
         best_ratio = 1.5 * best_ratio + EPSILON;  
         if (best_ratio < best_end_ratio)
           best_ratio = best_end_ratio;

         while (left < right)
           {
            ratio = (double) edit_array [e] [left] . score
                 / (edit_array [e] [left] . len + EPSILON);
            if (! edit_array [e] [left] . at_end && ratio <= best_ratio)
              break;
            edit_array [e] [left] . at_end = TRUE;
            left ++;
           }

         while (left < right)
           {
            ratio = (double) edit_array [e] [right] . score
                 / (edit_array [e] [right] . len + EPSILON);
            if (! edit_array [e] [right] . at_end && ratio <= best_ratio)
              break;
            edit_array [e] [right] . at_end = TRUE;
            right --;
           }
        }

      assert (left <= right);
     }

   // if there was a match to the end at least as good as the score limit
   // return it
   if (best_end_score <= score_limit)
     {
      Homopoly_Match_Entry_t  * bp = edit_array [best_end_e] + best_end_d;

      (* a_end) = bp -> len;
      (* t_end) = bp -> len + best_end_d;
      Set_Fwd_Homopoly_Delta (delta, delta_len, edit_array, best_end_e, best_end_d);
      (* match_to_end) = TRUE;
      (* return_score) = best_end_score;

      return  best_end_e;
     }

   // For failed alignments just return the exact match at the beginning
   (* a_end) = edit_array [0] [0] . len;
   (* t_end) = edit_array [0] [0] . len;
   Set_Fwd_Homopoly_Delta (delta, delta_len, edit_array, 0, 0);
   (* match_to_end) = FALSE;
   (* return_score) = 0;

   return 0;
  }


int  Fwd_Prefix_Edit_Dist
  (char a_string [], int m, char t_string [], int n, int error_limit,
   int * a_end, int * t_end, int * match_to_end,
   double match_value, int * delta, int * delta_len, int ** edit_array,
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


int  Rev_Homopoly_Match_Start
  (const char * A, int m, const char * T, int n, int score_limit,
   int * return_score, int * a_end, int * t_end, int * match_to_end,
   Homopoly_Match_Entry_t ** edit_array)

// Return the minimum number of homopolymer-type changes (indels & substitutions)
// for the best match a prefix of string  A [0 .. (1-m)]  right-to-left
// with a prefix of string  T [0 .. (1-n)] , also right-to-left,
// with alignment score  <= score_limit .  Note the subscripts are negative.
// The return value is the number of changes and  return_score  is set to the
// score of the best alignment.
// The match must extend to the beginning of one of the two strings.
// Set  a_end  and  t_end  to the leftmost positions where the
// alignment ended in  A  and  T , respectively, as negative values.
// Set  match_to_end  true iff the match extended to the end
// of at least one string.
//  edit_array  has preallocated storage that can be used
// for this computation (including in a threaded environment).

  {
   Homopoly_Match_Entry_t  * ep;
   int  best_end_d, best_end_e, best_end_score;
   double  ratio, best_ratio, best_end_ratio;
   int  from, lowest_score, score, shorter;
   int  row, left, right;
   int  d, e, j, k, s;

   shorter = OVL_Min_int (m, n);
   for (row = 0; row < shorter && A [- row] == T [- row]; row ++)
     ;
   edit_array [0] [0] . len = row;
   edit_array [0] [0] . score = 0;
   edit_array [0] [0] . at_end = (row == shorter);

   if (row == shorter)     // Exact match to end
     {
      (* a_end) = (* t_end) = - row;
      (* match_to_end) = TRUE;
      (* return_score) = 0;
      return 0;
     }

   lowest_score = 0;
   best_end_score = MAX_HOMOPOLY_SCORE;
   best_end_ratio = MAX_HOMOPOLY_SCORE;
   left = right = 0;
   for (e = 1; lowest_score < score_limit && e < MAX_ERRORS; e ++)
     {
      left = OVL_Max_int (left - 1, -e);
      right = OVL_Min_int (right + 1, e);
      lowest_score = MAX_HOMOPOLY_SCORE;
      best_ratio = DBL_MAX;

      for (d = left; d <= right; d ++)
        {
         // first find the best score in the preceding row
         row = 0;
         score = MAX_HOMOPOLY_SCORE;
         from = 2;

//printf ("e=%d  d=%d  left=%d  right=%d\n", e, d, left, right);
         // check cell directly above
         ep = edit_array [e - 1] + d;
         if (left < d && d < right && ! ep -> at_end)
           {
//printf ("check above\n");
            row = 1 + ep -> len;
            score = ep -> score + HP_SUBST_SCORE;
            from = 0;
           }

         // check upper left cell
         ep = edit_array [e - 1] + d - 1;
         if (left + 1 < d && ! ep -> at_end)
           {
//printf ("check upper left  e/d=%d/%d  left=%d\n", e, d, left);
            j = ep -> len - 1;;
            k = j + d;
            if (k > 0 && T [- k + 1] == T [- k]
                  && (T [- k] == A [- j] || ep -> from == -1))
                      // if ep -> from is -1 then either T [-k] itself
                      // was inserted or else it matched A [-j]
              s = ep -> score + HP_INDEL_SCORE;
                  // homopoly insert in T
            else
              s = ep -> score + NON_HP_INDEL_SCORE;
                  // non-homopoly insert in T
            if (s < score || (s == score && row < ep -> len))   // low scores are better
              {
//printf ("is better\n");
               row = ep -> len;
               score = s;
               from = -1;
              }
           }

         // check upper right cell
         ep = edit_array [e - 1] + d + 1;
         if (d < right - 1 && ! ep -> at_end)
           {
//printf ("check upper right\n");
            j = ep -> len;
            k = j + d;
            if (j > 0 && A [- j + 1] == A [- j]
                  && (A [- j] == T [- k] || ep -> from == +1))
              s = ep -> score + HP_INDEL_SCORE;
                   // homopoly insert in A
            else
              s = ep -> score + NON_HP_INDEL_SCORE;
                   // non-homopoly insert in A
            if (s < score || (s == score && row < ep -> len))   // low scores are better
              {
//printf ("is better\n");
               row = ep -> len + 1;
               score = s;
               from = +1;
              }
           }

         // now extend as far as matches allow
         while (row < m && row + d < n && A [- row] == T [- row - d])
            row ++;

         edit_array [e] [d] . len = row;
         edit_array [e] [d] . score = score;
         edit_array [e] [d] . at_end = (row == m || row + d == n);
         edit_array [e] [d] . from = from;

         if (score == MAX_HOMOPOLY_SCORE)
           {  // there was no legal value on the prior row so
              // set the at_end of this cell true so it won't be
              // used by any cells below it
//printf ("is hopeless\n");
            edit_array [e] [d] . at_end = TRUE;
            continue;
           }

         if (MIN_RATIO_E <= e)
           {
            ratio = (double) score / (row + EPSILON);  // add EPSION in case row=0
            if (ratio < best_ratio)
              best_ratio = ratio;
           }

         if (edit_array [e] [d] . at_end && score < best_end_score)
           {
            best_end_score = score;
            best_end_e = e;
            best_end_d = d;
            best_end_ratio = (double) score / (edit_array [e] [d] . len + EPSILON);
            score_limit = best_end_score;
           }
         if (score < lowest_score)
           lowest_score = score;
        }

//**ALD
if (0 && best_end_score < MAX_HOMOPOLY_SCORE)
  {
   printf ("Global_Debug_Flag=%d\n", Global_Debug_Flag);
   printf ("Rev_Homopoly_Match_Start:  best_end_e=%d  best_end_d=%d"
        "  best_end_score=%d\n", best_end_e, best_end_d, best_end_score);
   Show_Homopoly_Match_Array (stdout, edit_array, e);
  }

      // if lowest score on row is worse (higher) than the best
      // match_to_end score, then have found the best match
      if (best_end_score <= lowest_score)
        {
         Homopoly_Match_Entry_t  * bp = edit_array [best_end_e] + best_end_d;

         (* a_end) = - (bp -> len);
         (* t_end) = - (bp -> len) - best_end_d;
         (* match_to_end) = TRUE;
         (* return_score) = best_end_score;
        
         return  best_end_e;
        }

      // shrink left..right range here based on hopelessly bad scores
      if (MIN_RATIO_E < e)
        {
         // trim cells whose score/len ratio is worse than 1.5 times the min
         best_ratio = 1.5 * best_ratio + EPSILON;  
         if (best_end_ratio < best_ratio)
           best_ratio = best_end_ratio;

         while (left < right)
           {
            ratio = (double) edit_array [e] [left] . score
                 / (edit_array [e] [left] . len + EPSILON);
            if (! edit_array [e] [left] . at_end && ratio <= best_ratio)
              break;
            edit_array [e] [left] . at_end = TRUE;
//**ALD
//printf ("shift left  e=%d  left=%d  ratio=%.3f  best_ratio=%.3f  best_end_ratio=%.3f\n",
//     e, left, ratio, best_ratio, best_end_ratio);
            left ++;
           }

         while (left < right)
           {
            ratio = (double) edit_array [e] [right] . score
                 / (edit_array [e] [right] . len + EPSILON);
            if (! edit_array [e] [right] . at_end && ratio <= best_ratio)
              break;
            edit_array [e] [right] . at_end = TRUE;
//**ALD
//printf ("shift right  e=%d  left=%d  ratio=%.3f  best_ratio=%.3f  best_end_ratio=%.3f\n",
//     e, left, ratio, best_ratio, best_end_ratio);
            right --;
           }
        }

      assert (left <= right);
     }

   // if there was a match to the end at least as good as the score limit
   // return it
   if (best_end_score <= score_limit)
     {
      Homopoly_Match_Entry_t  * bp = edit_array [best_end_e] + best_end_d;

      (* a_end) = - (bp -> len);
      (* t_end) = - (bp -> len) - best_end_d;
      (* match_to_end) = TRUE;
      (* return_score) = best_end_score;

      return  best_end_e;
     }

   // For failed alignments just return the exact match at the beginning
   (* a_end) = - edit_array [0] [0] . len;
   (* t_end) = - edit_array [0] [0] . len;
   (* match_to_end) = FALSE;
   (* return_score) = 0;

   return 0;
  }


int  Rev_Prefix_Edit_Dist
  (char a_string [], int m, char t_string [], int n, int error_limit,
   int * a_end, int * t_end, int * leftover, int * match_to_end,
   double match_value, int * delta, int * delta_len, int ** edit_array,
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
//  edit_array  has preallocated storage that can be used
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


void  Set_Fwd_Banded_Delta
  (int delta [], int * delta_len, Alignment_Cell_t edit_space [],
   int row, int col, int width [], int indent [])

// Set  delta  to the entries indicating the insertions/deletions
// in the alignment encoded in  edit_space  ending at logical
// cell  row,col  and going back to cell  0,0 .   width  has the
// actual number of cells on each row and  indent  has the offset
// to the first actual cell on each row.
// Set  (* delta_len)  to the number of entries in  delta .

  {
   Alignment_Cell_t  * curr;
   int  first = 1;
   int  ct, d_len, offset, sign;
   int  i, j, r, c;

   // get offset to start
   offset = 0;
   for (i = 0; i < row; i ++)
     offset += width [i];

   // will put values on delta backwards starting at delta [0]
   // and then reverse it at the end
   ct = d_len = 0;

//**ALD
   if (Verbose_Level > 2)
     printf ("offset=%d  indent[%d]=%d\n", offset, row, indent [row]);
   r = row;
   c = col;
   curr = edit_space + offset + c - indent [r];

   while (0 < r || 0 < c)
     {
//**ALD
      if (Verbose_Level > 2)
        printf ("r=%d  c=%d  score=%d  from=%d  offset=%d  first=%d  d_len=%d\n", r, c,
             curr -> score, curr -> from, offset, first, d_len);
      switch (curr -> from)
        {
         case -1 :  // go left
           assert (0 < c);
           if (first)
             first = 0;
           else
             delta [d_len ++] = sign * ct;
           sign = -1;
           ct = 1;
           c --;
           curr --;
           break;

         case 0 :  // go top left
           assert (0 < r && 0 < c);
           ct ++;
           r --;
           c --;
           offset -= width [r];
           curr = edit_space + offset + c - indent [r];
           break;
           
         case +1 :  // go top
           assert (0 < r);
           if (first)
             first = 0;
           else
             delta [d_len ++] = sign * ct;
           sign = +1;
           ct = 1;
           r --;
           offset -= width [r];
           curr = edit_space + offset + c - indent [r];
           break;

         default :
           fprintf (stderr, "ERROR:  line %d  file %s\n", __LINE__, __FILE__);
           fprintf (stderr, "Bad from value = %d  r = %d  c = %d\n", curr -> from,
                r, c);
           exit (EXIT_FAILURE);
        }
     }

   if (! first)
     delta [d_len ++] = sign * ct;

   // now reverse the contents of delta
   for (i = 0, j = d_len - 1; i < j; i ++, j --)
     {
      int  save = delta [i];

      delta [i] = delta [j];
      delta [j] = save;
     }
   (* delta_len) = d_len;

   return;
  }


void  Set_Fwd_Delta
  (int delta [], int * delta_len, int ** edit_array,
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


void  Set_Fwd_Homopoly_Delta
  (int delta [], int * delta_len, Homopoly_Match_Entry_t ** edit_array,
   int e, int d)

// Set  delta  to the entries indicating the insertions/deletions
// in the alignment encoded in  edit_array  ending at position
//  edit_array [e] [d] .  This is the position in the first
// string where the alignment ended.  Set  (* delta_len)  to
// the number of entries in  delta .

  {
   Homopoly_Match_Entry_t  * last, * curr, * next;
   int  d_len, from, max;
   int  i, j, k;

   last = edit_array [e] + d;
   d_len = 0;

   // temporarily put values in  delta  starting at subscript  e
   // and working down.

   for (k = e; k > 0; k --)
     {
      curr = edit_array [k] + d;
      next = edit_array [k - 1] + d + curr -> from;

//**ALD
//printf ("e=%d  d=%d  len=%u  from=%d\n", k, d, curr -> len, curr -> from);
      switch (curr -> from)
        {
         case -1 :
           delta [e - d_len] = next -> len - last -> len - 1;
           d_len ++;
           d --;
           last = next;
           break;
         case +1 :
           delta [e - d_len] = last -> len - next -> len;
           d_len ++;
           d ++;
           last = next;
           break;
         case 0 :
           // do nothing
           break;
        }
     }

   delta [e - d_len] = last -> len + 1;

   // Now shift values to the front of  delta
   k = 0;
   for (i = e - d_len; i < e; i ++)
     delta [k ++] = abs (delta [i]) * Sign (delta [i + 1]);

   (* delta_len) = d_len;

//**ALD
if (0)
{
 printf ("delta_len = %d\n", * delta_len);
 for (i = 0; i < d_len; i ++)
   printf ("delta [%d] = %d\n", i, delta [i]);
}

   return;
  }


void  Set_Rev_Delta
  (int delta [], int * delta_len, int ** edit_array,
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


void  Show_Homopoly_Match_Array
  (FILE * fp, Homopoly_Match_Entry_t ** hp, int e)

// Display the values in pyramidal array  hp  in rows  0 .. e  to
// file  fp .

  {
   int  i, j;

fprintf (fp, "e=%d\n", e);
   for (i = 0; i <= e; i ++)
     {
      fprintf (fp, "%2d: ", i);
      for (j = - e; j < - i; j ++)
        fprintf (fp, " %8s", "");
      for (j = - i; j <= i; j ++)
        fprintf (fp, " %3d%+2d%c%-2d", hp [i] [j] . len, hp [i] [j] . from,
             (hp [i] [j] . at_end ? 'X' : ':'),
             hp [i] [j] . score);
      fputc ('\n', fp);
     }

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



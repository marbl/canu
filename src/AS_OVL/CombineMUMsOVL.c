
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
* Module:  CombineMUMsOVL.c
* Description:
*   Take clusters of MUMs (from  mgaps  program) and do alignments
*   to combine them into matches.  Use overlapper extension
*   parameters.  Also extend matches off ends.
*   Created by modifying  AutoScreenOVL.c
* 
*    Programmer:  A. Delcher
*       Started:   26 Apr 2001
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: CombineMUMsOVL.c,v 1.6 2007-02-12 22:16:57 brianwalenz Exp $
 * $Revision: 1.6 $
*/

static char CM_ID[] = "$Id: CombineMUMsOVL.c,v 1.6 2007-02-12 22:16:57 brianwalenz Exp $";


//  System include files

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <fcntl.h>
#include  <sys/types.h>
#include  <string.h>
#include  <dirent.h>
#include  <sys/stat.h>
#include  <unistd.h>


//  Local include files

#include  "AutoScreenOVL.h"
#include  "AS_PER_ReadStruct.h"
#include  "AS_PER_genericStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_UTL_version.h"
#include  "AS_CGW_dataTypes.h"
#include  "Globals_CGW.h"
#include  "ScaffoldGraph_CGW.h"



//  Constants

#define  BRANCH_PT_MATCH_VALUE    0.272
    //  Value to add for a match in finding branch points
    //  1.20 was the calculated value for 6% vs 35% error discrimination
    //  Converting to integers didn't make it faster
static float branch_pt_match_value=BRANCH_PT_MATCH_VALUE;

#define  BRANCH_PT_ERROR_VALUE    -0.728
    //  Value to add for a mismatch in finding branch points
    //   -2.19 was the calculated value for 6% vs 35% error discrimination
    //  Converting to integers didn't make it faster
static float  branch_pt_error_value=BRANCH_PT_ERROR_VALUE;

static double pid;

#define  EDIT_DIST_PROB_BOUND        1e-4
    //  Probability limit to "band" edit-distance calculation
    //  Determines  NORMAL_DISTRIB_THOLD
#define  ERRORS_FOR_FREE             1
    //  The number of errors that are ignored in setting probability
    //  bound for terminating alignment extensions in edit distance
    //  calculations
#define  EXPANSION_FACTOR            1.4
    //  Factor by which to grow memory when realloc'ing
#define  GIVE_UP_LEN                 200
    //  Stop alignment when go this many bases past max-score value
#define  INITIAL_BUFFER_SIZE         10000
    //  Initial number of bytes to use for sequence strings
#define  MATCH_SCORE                 1.0
    //  Score to add for match in finding length of alignment
#define  MAX_ERROR_RATE              AS_GUIDE_ERROR_RATE
    //  The largest error allowed in overlaps
#define  MAX_FRAG_LEN                2048
    //  The longest fragment allowed
#define  MAX_ERRORS                  500
    //  Most errors in any edit distance computation
#define  MAX_EXTENSION               10000
    //  Maximum extension match out from a match
#define  MAX_HDR_LEN                 10000
    //  Longest allowable fasta header line
#define  MAX_LINE                    10000
    //  Longest possible fasta data line
#define  MAX_MEMORY_STORE            50000
    //  The most fragments allowed in memory store
#define  MIN_BRANCH_END_DIST         20
    //  Branch points must be at least this many bases from the
    //  end of the fragment to be reported
#define  MIN_BRANCH_TAIL_SLOPE       0.20
    //  Branch point tails must fall off from the max by at least
    //  this rate
#define  MIN_MATCH_LEN               40
    //  Number of bases in match region in order to count it
#define  MISMATCH_SCORE              -3.0
    //  Score to add for non-match in finding length of alignment
#define  NORMAL_DISTRIB_THOLD        3.62
    //  Determined by  EDIT_DIST_PROB_BOUND
#define  REALLY_VERBOSE              0
    //  If  1  prints tons of stuff
#define  VERBOSE                     0
    //  If  1  prints lots of stuff



//  Type definitions

typedef  struct s_Cover_t
  {
   long int  lo, hi;
   struct s_Cover_t  * next;
  }  Cover_t;



//  Static Globals

static int  Consec_Non_ACGT = 0;
    // Stop alignments when encounter at least this many non-acgt characters.
static int  * Edit_Array [MAX_ERRORS];
    // Use for alignment calculation.  Points into  Edit_Space .
static int  Edit_Match_Limit [MAX_ERRORS] = {0};
    // This array [e] is the minimum value of  Edit_Array [e] [d]
    // to be worth pursuing in edit-distance computations between guides
static int  Edit_Space [(MAX_ERRORS + 4) * MAX_ERRORS];
    // Memory used by alignment calculation
static int  Error_Bound [MAX_FRAG_LEN + 1];
    //  This array [i]  is the maximum number of errors allowed
    //  in a match between sequences of length  i , which is
    //  i * MAXERROR_RATE .
static int  Fill_Ct = 0;
    // Number of non-acgt bases in ref sequence
static char  * Gaps_File_Path = NULL;
    // Name of file produced by  mgaps  program
static char  * Match_File_Path = NULL;
    // Name of multifasta file of sequences to compare against the reference
    // sequence
static int UserScoring = FALSE;
    // If TRUE, then user specified a percent ID cutoff and scoring
    // is adjusted to enforce this in extensions
static int  Nucleotides_Only = FALSE;
    // If  TRUE , then only acgt's can match
static int  Output_Cover_Files = TRUE;
    // If  TRUE , output files showing coverage of each genome.
static char  * Query = NULL;
    // The query sequence
static long int  Query_Len;
    // The length of the query sequence
static char  * Query_Suffix = "Query";
    // Suffix for query tag
static char  * Ref = NULL;
    // The reference sequence
static char  * Ref_File_Path = NULL;
    // Name of (single) fasta file of reference sequence
static long int  Ref_Len;
    // The length of the reference sequence
static long int  Ref_Size;
    // The size of the reference sequence buffer
static char  * Ref_Suffix = "Ref";
    // Suffix for reference tag
static int  Show_Differences = FALSE;
    // If  TRUE  then show differences in all alignments
static int  Tag_From_Fasta_Line = FALSE;
    // If  TRUE  then use fasta tag from ref & query sequences as
    // 1st & 2nd column, resp., on match line to identify matches
static int  Verbose = 0;
    // Controls printing of extra debuggin information



//  Static Functions

static void  Add_Coverage
    (Cover_t * * list, long int lo, long int hi);
static int  Binomial_Bound
    (int, double, int, double);
static char  Complement
    (char);
static void  Display_Alignment
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct);
static int  Extend_Backward
    (long int * ref_lo, long int * query_lo);
static int  Extend_Forward
    (long int * ref_hi, long int * query_hi);
static void  Fasta_Print
    (FILE * fp, char * s, char * hdr);
static char  Filter
    (char ch);
static Overlap *  Find_Overlap
    (char * seq1, char * seq2, ChunkOrientationType orientation, 
     int min_ahang, int max_ahang,
     double erate, double thresh, int minlen, CompareOptions what);
static void  Initialize_Globals
    (void);
extern FILE *  Local_File_Open
    (const char * filename, const char * mode);
static int  OVL_Max_int
    (int a, int b);
static int  OVL_Min_int
    (int a, int b);
static void  Parse_Command_Line
    (int argc, char * argv []);
static int  Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Match_To_End,
     int Delta [MAX_ERRORS], int * Delta_Len, int extending);
static int  Read_String
    (FILE * fp, char * * T, long int * Size, char header []);
static void  Rev_Complement
    (char * s);
static void  Rev_Display_Alignment
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct);
static int  Rev_Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Match_To_End,
     int Delta [MAX_ERRORS], int * Delta_Len, int extending);
static void  Rev_Show_Diffs
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct,
     long int a_ref, long int b_ref);
static void  Set_Deltas
    (int delta [], int * delta_len, int row, int d, int e);
static void  Show_Coverage
    (Cover_t * list, char * filename, char * tag, char * seq);
static void  Show_Diffs
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct,
     long int a_ref, long int b_ref);
static int  Sign
    (int a);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   FILE  * fp, * match_fp;
   char  ref_header [MAX_HDR_LEN];
   char  header [MAX_HDR_LEN], ref_tag [MAX_HDR_LEN], query_tag [MAX_HDR_LEN];
   char  line [MAX_LINE], filename [MAX_LINE];
   long int  ref_pos, ref_lo, ref_hi;
   long int  query_pos, query_lo, query_hi, query_size;
   long int  match_len;
   double  avg_len, error = 0.0;
   long int  total_errors = 0;
   int  match_ct = 0;
   double  match_total = 0.0;
   int  is_forward = FALSE;
   int  first_match = TRUE;
   Cover_t  * ref_cover_list = NULL;
   Cover_t  * query_cover_list = NULL;
   char  * p;
   int  i;

   Parse_Command_Line  (argc, argv);

   if(UserScoring){
     // percentID = mismatch penalty / (match bonus + mismatch penalty)
     // or, more compactly, p = mm / (m+mm)
     // Since we require that m+mm = 1, m=1-mm
     // and so we get the trivial mm = p
     //
     // This means that branch_pt_error_value = -p
     // and then branch_pt_match_value = 1-p
     // 
     // Make it so:
     branch_pt_error_value = -pid;
     branch_pt_match_value = 1.-pid;
   }

   Initialize_Globals ();

   p = strrchr (Ref_File_Path, '/');
   if  (p == NULL)
       strcpy (ref_tag, Ref_File_Path);
     else
       strcpy (ref_tag, p + 1);
   p = strrchr (ref_tag, '.');
   if  (p != NULL)
       (* p) = '\0';
   strcat (ref_tag, Ref_Suffix);

   p = strrchr (Match_File_Path, '/');
   if  (p == NULL)
       strcpy (query_tag, Match_File_Path);
     else
       strcpy (query_tag, p + 1);
   p = strrchr (query_tag, '.');
   if  (p != NULL)
       (* p) = '\0';
   strcat (query_tag, Query_Suffix);

   fp = Local_File_Open (Ref_File_Path, "r");
   Ref_Size = 100000;
   Ref = (char *) safe_malloc (Ref_Size);

   if(Read_String (fp, & Ref, & Ref_Size, ref_header) == 0)
     assert(0);
   if  (Tag_From_Fasta_Line)
       strcpy (ref_tag, strtok (ref_header, " \t\n>"));
   fclose (fp);

   Ref_Len = strlen (Ref + 1);
   Ref = (char *) safe_realloc (Ref, 1 + Ref_Len);
   for  (i = 1;  i <= Ref_Len;  i ++)
     switch  (Ref [i])
       {
        case  'a' :
        case  'c' :
        case  'g' :
        case  't' :
          break;
        default :
          if  (Nucleotides_Only)
              Ref [i] = '2';
          Fill_Ct ++;
       }

   fp = Local_File_Open (Gaps_File_Path, "r");
   match_fp = Local_File_Open (Match_File_Path, "r");
   query_size = 100000;
   Query = (char *) safe_malloc (query_size);

   while  (fgets (line, MAX_LINE, fp) != NULL)
     {
      if  (line [0] == '>')
          {
           if  (! first_match)
               {
                total_errors += Extend_Forward (& ref_hi, & query_hi);
                avg_len = (1 + ref_hi - ref_lo + 1 + abs (query_hi - query_lo)) / 2.0;
                error = avg_len == 0.0 ? 0.0 : total_errors / avg_len;
                printf ("%-10s %-10s  %9ld %9ld  %9ld %9ld  %9ld %9ld  %5.2f\n",
                        ref_tag, query_tag, Ref_Len, Query_Len,
                        ref_lo, ref_hi,
                        is_forward ? query_lo : 1 + Query_Len - query_lo,
                        is_forward ? query_hi : 1 + Query_Len - query_hi,
                        100.0 * error);
                match_ct ++;
                match_total += 1 + ref_hi - ref_lo;
                total_errors = 0;
                Add_Coverage (& ref_cover_list, ref_lo, ref_hi);
                if  (is_forward)
                    Add_Coverage (& query_cover_list, query_lo, query_hi);
                  else
                    Add_Coverage (& query_cover_list, 1 + Query_Len - query_hi,
                                  1 + Query_Len - query_lo);
               }
           if  (Tag_From_Fasta_Line)
               {
                char  * p;

                p = strtok (line + 1, " >\t\n");
                if  (strlen (p) > 0)
                    strcpy (query_tag, p);
                  else
                    {
                     fprintf (stderr, "No tag on line %s", line);
                     strcpy (query_tag, "???");
                    }
               }
           is_forward = ! is_forward;
           if  (is_forward)
               {
                 if(Read_String (match_fp, & Query, & query_size, header) == 0)
                   assert(0);
                Query_Len = strlen (Query + 1);
                if  (Nucleotides_Only)
                    {
                     for  (i = 1;  i <= Query_Len;  i ++)
                       switch  (Query [i])
                         {
                          case  'a' :
                          case  'c' :
                          case  'g' :
                          case  't' :
                            break;
                          default :
                            Query [i] = '1';
                         }
                    }
               }
             else
               Rev_Complement (Query + 1);
           first_match = TRUE;
          }
      else if  (line [0] == '#')
          {
           if  (! first_match)
               {
                total_errors += Extend_Forward (& ref_hi, & query_hi);
                avg_len = (1 + ref_hi - ref_lo + 1 + abs (query_hi - query_lo)) / 2.0;
                error = avg_len == 0.0 ? 0.0 : total_errors / avg_len;
                printf ("%-10s %-10s  %9ld %9ld  %9ld %9ld  %9ld %9ld  %5.2f\n",
                        ref_tag, query_tag, Ref_Len, Query_Len,
                        ref_lo, ref_hi,
                        is_forward ? query_lo : 1 + Query_Len - query_lo,
                        is_forward ? query_hi : 1 + Query_Len - query_hi,
                        100.0 * error);
                match_ct ++;
                match_total += 1 + ref_hi - ref_lo;
                total_errors = 0;
                Add_Coverage (& ref_cover_list, ref_lo, ref_hi);
                if  (is_forward)
                    Add_Coverage (& query_cover_list, query_lo, query_hi);
                  else
                    Add_Coverage (& query_cover_list, 1 + Query_Len - query_hi,
                                  1 + Query_Len - query_lo);
               }
           first_match = TRUE;
          }
        else
          {
            if(sscanf (line, "%ld %ld %ld",
                       & ref_pos, & query_pos, & match_len) != 3)
              assert(0);
           if  (first_match)
               {
                ref_lo = ref_pos;
                query_lo = query_pos;
                total_errors += Extend_Backward (& ref_lo, & query_lo);
               }
             else
               {
                int  errors;
                int  a_len, b_len;
                int  a_end, b_end, match_to_end;
                int  delta [MAX_ERRORS], delta_len;

                a_len = ref_pos - ref_hi - 1;
                b_len = query_pos - query_hi - 1;
                
                errors = Prefix_Edit_Dist
                           (Ref + ref_hi + 1, a_len,
                            Query + query_hi + 1, b_len,
                            MAX_ERRORS - 1, & a_end, & b_end,
                            & match_to_end, delta, & delta_len, FALSE);
                if  (Show_Differences)
                    Show_Diffs (Ref + ref_hi + 1, a_end,
                                Query + query_hi + 1, b_end,
                                delta, delta_len, ref_hi + 1, query_hi + 1);
                if  (Verbose > 0)
                    {
                     printf ("%9ld %9ld  %4d %4d  %3d  %c  %4d %4d\n",
                             ref_pos, query_pos, a_len, b_len, errors,
                             match_to_end ? 'T' : 'F', a_end, b_end);
                     Display_Alignment (Ref + ref_hi + 1, a_end,
                                        Query + query_hi + 1, b_end,
                                        delta, delta_len);
                    }
                total_errors += errors;

                if  (! match_to_end)
                    {
                     ref_hi += a_end;
                     query_hi += b_end;
                     avg_len = (1 + ref_hi - ref_lo + 1
                                 + abs (query_hi - query_lo)) / 2.0;
                     error = avg_len == 0.0 ? 0.0 : total_errors / avg_len;
                     printf ("%-10s %-10s  %9ld %9ld  %9ld %9ld  %9ld %9ld  %5.2f\n",
                             ref_tag, query_tag, Ref_Len, Query_Len,
                             ref_lo, ref_hi,
                             is_forward ? query_lo : 1 + Query_Len - query_lo,
                             is_forward ? query_hi : 1 + Query_Len - query_hi,
                             100.0 * error);
                     match_ct ++;
                     match_total += 1 + ref_hi - ref_lo;
                     total_errors = 0;
                     Add_Coverage (& ref_cover_list, ref_lo, ref_hi);
                     if  (is_forward)
                         Add_Coverage (& query_cover_list, query_lo, query_hi);
                       else
                         Add_Coverage (& query_cover_list, 1 + Query_Len - query_hi,
                                       1 + Query_Len - query_lo);
                     ref_lo = ref_pos;
                     query_lo = query_pos;
                     total_errors += Extend_Backward (& ref_lo, & query_lo);
                    }
               }

           ref_hi = ref_pos + match_len - 1;
           query_hi = query_pos + match_len - 1;
           first_match = FALSE;
          }
     }

   if  (! first_match)
       {
        total_errors += Extend_Forward (& ref_hi, & query_hi);
        avg_len = (1 + ref_hi - ref_lo + 1
                    + abs (query_hi - query_lo)) / 2.0;
        error = avg_len == 0.0 ? 0.0 : total_errors / avg_len;
        printf ("%-10s %-10s  %9ld %9ld  %9ld %9ld  %9ld %9ld  %5.2f\n",
                ref_tag, query_tag, Ref_Len, Query_Len,
                ref_lo, ref_hi,
                is_forward ? query_lo : 1 + Query_Len - query_lo,
                is_forward ? query_hi : 1 + Query_Len - query_hi,
                100.0 * error);
        match_ct ++;
        match_total += 1 + ref_hi - ref_lo;
        Add_Coverage (& ref_cover_list, ref_lo, ref_hi);
        if  (is_forward)
            Add_Coverage (& query_cover_list, query_lo, query_hi);
          else
            Add_Coverage (& query_cover_list, 1 + Query_Len - query_hi,
                          1 + Query_Len - query_lo);
       }

   fprintf (stderr, "           Ref len = %ld\n", Ref_Len);
   fprintf (stderr, "            acgt's = %ld\n", Ref_Len - Fill_Ct);
   fprintf (stderr, "        Non acgt's = %d\n", Fill_Ct);
   fprintf (stderr, " Number of matches = %d\n", match_ct);
   fprintf (stderr, "Sum of match bases = %.0f\n", match_total);
   fprintf (stderr, "   Avg match bases = %.0f\n",
            match_ct == 0 ? 0.0 : match_total / match_ct);

   if  (Output_Cover_Files)
       {
        strcpy (filename, ref_tag);
        strcat (filename, ".cover");
        Show_Coverage (ref_cover_list, filename, ref_tag, Ref);

        strcpy (filename, query_tag);
        strcat (filename, ".cover");
        Rev_Complement (Query + 1);
        Show_Coverage (query_cover_list, filename, query_tag, Query);
       }

   return  0;
  }



static void  Add_Coverage
    (Cover_t * * list, long int lo, long int hi)

//  Add  lo .. hi  to list of regions covered in  (* list) .
//  Combine nodes when appropriate.

  {
   Cover_t  * new_node, * p, * prev = NULL;

   if  ((* list) == NULL || hi + 1 < (* list) -> lo)
       {
        new_node = (Cover_t *) safe_malloc (sizeof (Cover_t));
        new_node -> lo = lo;
        new_node -> hi = hi;
        new_node -> next = (* list);
        (* list) = new_node;
        return;
       }

   for  (p = (* list);  p != NULL && lo - 1 > p -> hi;  p = p -> next)
     prev = p;

   if  (p == NULL || hi + 1 < p -> lo)
       {  // insert between or on end
        assert (prev != NULL);
        new_node = (Cover_t *) safe_malloc (sizeof (Cover_t));
        new_node -> lo = lo;
        new_node -> hi = hi;
        new_node -> next = p;
        prev -> next = new_node;
        return;
       }

   if  (lo < p -> lo)
       p -> lo = lo;
   while  (p -> next != NULL && hi + 1 >= p -> next -> lo)
     {
      Cover_t  * save;

      p -> hi = p -> next -> hi;
      save = p -> next;
      p -> next = save -> next;
      free (save);
     }
   if  (hi > p -> hi)
       p -> hi = hi;

   return;
  }



static int  Binomial_Bound
    (int e, double p, int Start, double Limit)

//  Return the smallest  n >= Start  s.t.
//    prob [>= e  errors in  n  binomial trials (p = error prob)]
//          > Limit

  {
   double  Normal_Z, Mu_Power, Factorial, Poisson_Coeff;
   double  q, Sum, P_Power, Q_Power, X;
   int  k, n, Bin_Coeff, Ct;

   q = 1.0 - p;
   if  (Start < e)
       Start = e;

   for  (n = Start;  n < MAX_FRAG_LEN;  n ++)
     {
      if  (n <= 35)
          {
           Sum = 0.0;
           Bin_Coeff = 1;
           Ct = 0;
           P_Power = 1.0;
           Q_Power = pow (q, n);

           for  (k = 0;  k < e && 1.0 - Sum > Limit;  k ++)
             {
              X = Bin_Coeff * P_Power * Q_Power;
              Sum += X;
              Bin_Coeff *= n - Ct;
              Bin_Coeff /= ++ Ct;
              P_Power *= p;
              Q_Power /= q;
             }
           if  (1.0 - Sum > Limit)
               return  n;
          }
        else
          {
           Normal_Z = (e - 0.5 - n * p) / sqrt (n * p * q);
           if  (Normal_Z <= NORMAL_DISTRIB_THOLD)
               return  n;
           Sum = 0.0;
           Mu_Power = 1.0;
           Factorial = 1.0;
           Poisson_Coeff = exp (- n * p);
           for  (k = 0;  k < e;  k ++)
             {
              Sum += Mu_Power * Poisson_Coeff / Factorial;
              Mu_Power *= n * p;
              Factorial *= k + 1;
             }
           if  (1.0 - Sum > Limit)
               return  n;
          }
     }

   return  MAX_FRAG_LEN;
  }



static char  Complement
    (char ch)

/*  Return the DNA complement of  ch . */

  {
   switch  (tolower ((int) ch))
     {
      case  'a' :
        return  't';
      case  'c' :
        return  'g';
      case  'g' :
        return  'c';
      case  't' :
        return  'a';
      default :
        if  (! isalpha (ch) && (! Nucleotides_Only || (ch != '1' && ch != '2')))
            {
             fprintf (stderr, "ERROR(complement):  Unexpected character `%c\'\n",
                      ch);
             exit (-1);
            }
     }

   return  tolower (ch);
  }



#define  DISPLAY_WIDTH   60

static void  Display_Alignment
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct)

//  Show (to  stdout ) the alignment encoded in  delta [0 .. (delta_ct - 1)]
//  between strings  a [0 .. (a_len - 1)]  and  b [0 .. (b_len - 1)] .

  {
   int  i, j, k, m, top_len, bottom_len;
   char  top [MAX_EXTENSION + 2 * MAX_ERRORS + 1];
   char  bottom [MAX_EXTENSION + 2 * MAX_ERRORS + 1];

   i = j = top_len = bottom_len = 0;
   for  (k = 0;  k < delta_ct;  k ++)
     {
      for  (m = 1;  m < abs (delta [k]);  m ++)
        {
         top [top_len ++] = a [i ++];
         j ++;
        }
      if  (delta [k] < 0)
          {
           top [top_len ++] = '-';
           j ++;
          }
        else
          {
           top [top_len ++] = a [i ++];
          }
     }
   while  (i < a_len && j < b_len)
     {
      top [top_len ++] = a [i ++];
      j ++;
     }
   top [top_len] = '\0';
     

   i = j = 0;
   for  (k = 0;  k < delta_ct;  k ++)
     {
      for  (m = 1;  m < abs (delta [k]);  m ++)
        {
         bottom [bottom_len ++] = b [j ++];
         i ++;
        }
      if  (delta [k] > 0)
          {
           bottom [bottom_len ++] = '-';
           i ++;
          }
        else
          {
           bottom [bottom_len ++] = b [j ++];
          }
     }
   while  (j < b_len && i < a_len)
     {
      bottom [bottom_len ++] = b [j ++];
      i ++;
     }
   bottom [bottom_len] = '\0';


   for  (i = 0;  i < top_len || i < bottom_len;  i += DISPLAY_WIDTH)
     {
      putchar ('\n');
      printf ("A: ");
      for  (j = 0;  j < DISPLAY_WIDTH && i + j < top_len;  j ++)
        putchar (top [i + j]);
      putchar ('\n');
      printf ("B: ");
      for  (j = 0;  j < DISPLAY_WIDTH && i + j < bottom_len;  j ++)
        putchar (bottom [i + j]);
      putchar ('\n');
      printf ("   ");
      for  (j = 0;  j < DISPLAY_WIDTH && i + j < bottom_len && i + j < top_len;
                j ++)
        if  (top [i + j] != ' ' && bottom [i + j] != ' '
                 && top [i + j] != bottom [i + j])
            putchar ('^');
          else
            putchar (' ');
      putchar ('\n');
     }

   return;
  }



static int  Extend_Backward
    (long int * ref_lo, long int * query_lo)

//  Do edit distance off end of match beginning at  (* ref_lo)  and
//  (* query_lo)  in the reference and query sequences, resp.,
//  in the reverse direction.
//  Return the number of errors in the extension and set  (* ref_hi)
//  and  (* query_hi)  to the end of the extension.

  {
   int  ct = 0, errors, sum = 0;
   int  a_len, b_len;
   int  a_end, b_end, match_to_end;
   int  delta [MAX_ERRORS], delta_len;

   do
     {
      a_len = (* ref_lo) - 1;
      if  (a_len > MAX_EXTENSION)
          a_len = MAX_EXTENSION;
      b_len = (* query_lo) - 1;
      if  (b_len > MAX_EXTENSION)
          b_len = MAX_EXTENSION;

      errors = Rev_Prefix_Edit_Dist
                 (Ref + (* ref_lo) - 1, a_len,
                  Query + (* query_lo) - 1, b_len,
                  MAX_ERRORS - 1, & a_end, & b_end,
                  & match_to_end, delta, & delta_len, TRUE);
      if  (Show_Differences)
          Rev_Show_Diffs
              (Ref + (* ref_lo - 1), a_end,
               Query + (* query_lo - 1), b_end,
               delta, delta_len, (* ref_lo) - 1, (* query_lo) - 1);
      if  (Verbose > 0)
          {
           printf ("revextend#%d  %9ld %9ld  %5d %5d  %3d  %c  %4d %4d\n",
                   ++ ct, (* ref_lo) - 1, (* query_lo) - 1, a_len, b_len,
                   errors, match_to_end ? 'T' : 'F', a_end, b_end);
           Rev_Display_Alignment
               (Ref + (* ref_lo) - 1, a_end, Query + (* query_lo) - 1, b_end,
                delta, delta_len);
          }

      (* ref_lo) -= a_end;
      (* query_lo) -= b_end;
      sum += errors;
     }  while  (a_end > 0.9 * MAX_EXTENSION || b_end > MAX_EXTENSION);

   return  sum;
  }



static int  Extend_Forward
    (long int * ref_hi, long int * query_hi)

//  Do edit distance off end of match ending at  (* ref_hi)  and
//  (* query_hi)  in the reference and query sequences, resp.
//  Return the number of errors in the extension and set  (* ref_hi)
//  and  (* query_hi)  to the end of the extension.

  {
   int  ct = 0, errors, sum = 0;
   int  a_end, b_end, match_to_end;
   int  a_len, b_len;
   int  delta [MAX_ERRORS], delta_len;

   do
     {
      a_len = Ref_Len - (* ref_hi);
      if  (a_len > MAX_EXTENSION)
          a_len = MAX_EXTENSION;
      b_len = Query_Len - (* query_hi);
      if  (b_len > MAX_EXTENSION)
          b_len = MAX_EXTENSION;
      
      errors = Prefix_Edit_Dist
                 (Ref + (* ref_hi) + 1, a_len,
                  Query + (* query_hi) + 1, b_len,
                  MAX_ERRORS - 1, & a_end, & b_end,
                  & match_to_end, delta, & delta_len, TRUE);
      if  (Show_Differences)
          Show_Diffs (Ref + (* ref_hi) + 1, a_end,
                      Query + (* query_hi) + 1, b_end,
                      delta, delta_len, (* ref_hi) + 1, (* query_hi) + 1);
      if  (Verbose > 0)
          {
           printf ("extend#%d  %9ld %9ld  %5d %5d  %3d  %c  %4d %4d\n",
                   ++ ct, (* ref_hi) + 1, (* query_hi) + 1, a_len, b_len,
                   errors, match_to_end ? 'T' : 'F', a_end, b_end);
           Display_Alignment (Ref + (* ref_hi) + 1, a_end,
                              Query + (* query_hi) + 1, b_end,
                              delta, delta_len);
          }

      (* ref_hi) += a_end;
      (* query_hi) += b_end;
      sum += errors;
     }  while  (a_end > 0.9 * MAX_EXTENSION || b_end > MAX_EXTENSION);

   return  sum;
  }



static void  Fasta_Print
    (FILE * fp, char * s, char * hdr)

//  Print string  s  in fasta format to  fp .  Put string  hdr
//  on header line.

  {
   int  ct = 0;

   fprintf (fp, ">%s\n", hdr);

   while  (* s != '\0')
     {
      if  (ct == 60)
          {
           fputc ('\n', fp);
           ct = 0;
          }
      fputc (* s, fp);
      s ++;
      ct ++;
     }

   fputc ('\n', fp);

   return;
  }



static char  Filter
    (char ch)

//  Convert  ch  to lowercase if necessary and if not 'a', 'c', 'g' or 't'
//  make it an 'a'.

  {
   ch = tolower (ch);

   switch  (ch)
     {
      case  'a' :
      case  'c' :
      case  'g' :
      case  't' :
        return  ch;
     }

   return  'a';
  }



static Overlap *  Find_Overlap
    (char * seq1, char * seq2, ChunkOrientationType orientation, 
     int min_ahang, int max_ahang,
     double erate, double thresh, int minlen, CompareOptions what)

//  Same as  OverlapSequences  in CGW

{
  Overlap *omesg;
  int flip = 0;
  
  // if the orientation is BA_AB or BA_BA, we need to reverse complement the first contig
  if (orientation == BA_AB || orientation == BA_BA)
	Complement_Seq( seq1 );
  
  // if the orientation is AB_BA or BA_BA, we need to set the flip variable for the second contig
  if (orientation == AB_BA || orientation == BA_BA)
	flip = 1;

  // min_ahang and end are essentially bounds on the a-hang
  omesg = DP_Compare(seq1, seq2,
					 min_ahang, max_ahang, (int) flip,
					 erate, thresh, minlen,
					 what);

  // return seq1 to its original state
  if (orientation == BA_AB || orientation == BA_BA)
	Complement_Seq( seq1 );
  
  // omesg->begpos is the a-hang, omesg->endpos is the b-hang
  return omesg;
}



static void  Initialize_Globals
    (void)

//  Initialize global variables used in this program

  {
   int  i, offset, del;
   int  e, start;

   offset = 2;
   del = 6;
   for  (i = 0;  i < MAX_ERRORS;  i ++)
     {
       Edit_Array [i] = Edit_Space + offset;
       offset += del;
       del += 2;
     }


   assert (MAX_ERROR_RATE >= AS_READ_ERROR_RATE
             && MAX_ERROR_RATE >= AS_GUIDE_ERROR_RATE);

   for  (i = 0;  i <= ERRORS_FOR_FREE;  i ++)
     Edit_Match_Limit [i] = 0;

   start = 1;
   for  (e = ERRORS_FOR_FREE + 1;  e < MAX_ERRORS;  e ++)
     {
      start = Binomial_Bound (e - ERRORS_FOR_FREE, MAX_ERROR_RATE,
                  start, EDIT_DIST_PROB_BOUND);
      Edit_Match_Limit [e] = start - 1;
      assert (Edit_Match_Limit [e] >= Edit_Match_Limit [e - 1]);
     }

   for  (i = 0;  i <= MAX_FRAG_LEN;  i ++)
     Error_Bound [i] = OVL_Max_int (10, 1 + (int) (2 * i * MAX_ERROR_RATE));

   return;
  }





static int  OVL_Max_int
    (int a, int b)

//  Return the larger of  a  and  b .

  {
   if  (a < b)
       return  b;

   return  a;
  }



static int  OVL_Min_int
    (int a, int b)

//  Return the smaller of  a  and  b .

  {
   if  (a < b)
       return  a;

   return  b;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
    int  ch, errflg = FALSE;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "nN:q:r:tv:xe:S")) != EOF))
     switch  (ch)
       {
        case 'e' : 
	  pid = 1.-atof(optarg);
	  UserScoring = TRUE;
	  break;

        case  'n' :
          Nucleotides_Only = TRUE;
          break;

        case  'N' :
          Consec_Non_ACGT = strtol (optarg, NULL, 10);
          break;

        case  'q' :
          Query_Suffix = optarg;
          break;

        case  'r' :
          Ref_Suffix = optarg;
          break;

        case  'S' :
          Show_Differences = TRUE;
          break;

        case  't' :
          Tag_From_Fasta_Line = TRUE;
          break;

        case  'v' :
          Verbose = strtol (optarg, NULL, 10);
          break;

        case  'x' :
          Output_Cover_Files = FALSE;
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc - 3)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   Ref_File_Path = argv [optind ++];

   Match_File_Path = argv [optind ++];

   Gaps_File_Path = argv [optind ++];

   return;
  }



static int  Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Match_To_End,
     int Delta [MAX_ERRORS], int * Delta_Len, int extending)

//  Return the minimum number of changes (inserts, deletes, replacements)
//  needed to match string  A [0 .. (m-1)]  with a prefix of string
//   T [0 .. (n-1)]  if it's not more than  Error_Limit .
//  Put delta description of alignment in  Delta  and set
//  (* Delta_Len)  to the number of entries there if it's a complete
//  match.
//  Set  A_End  and  T_End  to the rightmost positions where the
//  alignment ended in  A  and  T , respectively.
//  Set  Match_To_End  true if the match extended to the end
//  of at least one string; otherwise, set it false to indicate
//  a branch point.
//  extending  indicates whether the match is for a fixed region (FALSE)
//  or is extending off the end of a match as far as possible (TRUE)

  {
   double  Score, Max_Score;
   int  Max_Score_Len = 0, Max_Score_Best_d = 0, Max_Score_Best_e = 0;
   int  Best_d, Best_e, Longest, Row;
   int  Left, Right;
   int  d, e, i, j, shorter;

//   assert (m <= n);
   Best_d = Best_e = Longest = 0;
   (* Delta_Len) = 0;

   if  (Consec_Non_ACGT > 0)
       {
        int  ct;

        for  (i = ct = 0;  i < m && ct < Consec_Non_ACGT;  i ++)
          {
           switch  (A [i])
             {
              case  'a' :
              case  'c' :
              case  'g' :
              case  't' :
                ct = 0;
                break;
              default :
                ct ++;
             }
          }
        if  (ct >= Consec_Non_ACGT)
            {
             m = i - ct;
             extending = TRUE;
            }
        
        for  (i = ct = 0;  i < n && ct < Consec_Non_ACGT;  i ++)
          {
           switch  (T [i])
             {
              case  'a' :
              case  'c' :
              case  'g' :
              case  't' :
                ct = 0;
                break;
              default :
                ct ++;
             }
          }
        if  (ct >= Consec_Non_ACGT)
            {
             n = i - ct;
             extending = TRUE;
            }
       }

   shorter = OVL_Min_int (m, n);
   for  (Row = 0;  Row < shorter && A [Row] == T [Row];  Row ++)
     ;

   Edit_Array [0] [0] = Row;

   if  (Row == m && Row == n)      // Exact match
       {
        (* Delta_Len) = 0;
        (* A_End) = (* T_End) = Row;
        (* Match_To_End) = ! extending;
        return  0;
       }

   Left = Right = 0;
   Max_Score = Row * branch_pt_match_value;
   Max_Score_Len = Row;
   Max_Score_Best_d = 0;
   Max_Score_Best_e = 0;

   for  (e = 1;  e <= Error_Limit;  e ++)
     {
      int  cutoff;

      Left = OVL_Max_int (Left - 1, -e);
      Right = OVL_Min_int (Right + 1, e);
      Edit_Array [e - 1] [Left] = -2;
      Edit_Array [e - 1] [Left - 1] = -2;
      Edit_Array [e - 1] [Right] = -2;
      Edit_Array [e - 1] [Right + 1] = -2;

      for  (d = Left;  d <= Right;  d ++)
        {
         Row = 1 + Edit_Array [e - 1] [d];
         if  ((j = Edit_Array [e - 1] [d - 1]) > Row)
             Row = j;
         if  ((j = 1 + Edit_Array [e - 1] [d + 1]) > Row)
             Row = j;
         while  (Row < m && Row + d < n
                  && A [Row] == T [Row + d])
           Row ++;

         Edit_Array [e] [d] = Row;

         if  (Row == m && Row + d == n)
             {
              if  (extending)
                  {
                   for  (j = Left;  j <= d;  j ++)
                     if  (Edit_Array [e] [j] > Longest)
                         {
                          Best_d = j;
                          Longest = Edit_Array [e] [j];
                         }
                   Score = Longest * branch_pt_match_value - e;
                   if  (Score < Max_Score)
                       {
                        (* A_End) = Max_Score_Len;
                        (* T_End) = Max_Score_Len + Max_Score_Best_d;
                        Set_Deltas (Delta, Delta_Len, Max_Score_Len,
                                    Max_Score_Best_d, Max_Score_Best_e);
                        (* Match_To_End) = FALSE;
                        return  Max_Score_Best_e;
                       }
                  }

              (* A_End) = Row;           // One past last align position
              (* T_End) = Row + d;
              Set_Deltas (Delta, Delta_Len, Row, d, e);
              (* Match_To_End) = ! extending;
              return  e;
             }
        }

#if  1
      cutoff = Longest - GIVE_UP_LEN;
      while  (Left <= Right && Edit_Array [e] [Left] < cutoff)
        Left ++;
      if  (Left > Right)
          break;
      while  (Right > Left && Edit_Array [e] [Right] < cutoff)
        Right --;
#else
      while  (Left <= Right && Left < 0
                  && Edit_Array [e] [Left] < Edit_Match_Limit [e])
        Left ++;
      if  (Left >= 0)
          while  (Left <= Right
                    && Edit_Array [e] [Left] + Left < Edit_Match_Limit [e])
            Left ++;
      if  (Left > Right)
          break;
      while  (Right > 0
                  && Edit_Array [e] [Right] + Right < Edit_Match_Limit [e])
        Right --;
      if  (Right <= 0)
          while  (Edit_Array [e] [Right] < Edit_Match_Limit [e])
            Right --;
#endif
      assert (Left <= Right);

      for  (d = Left;  d <= Right;  d ++)
        if  (Edit_Array [e] [d] > Longest)
            {
             Best_d = d;
             Best_e = e;
             Longest = Edit_Array [e] [d];
            }
#if  1
      Score = Longest * branch_pt_match_value - e;
               // Assumes  branch_pt_match_value - branch_pt_error_value == 1.0
      if  (Score > Max_Score)
          {
           Max_Score = Score;
           Max_Score_Len = Longest;
           Max_Score_Best_d = Best_d;
           Max_Score_Best_e = Best_e;
          }
#endif

      if  (Longest - Max_Score_Len >= GIVE_UP_LEN)
          break;
     }

   (* A_End) = Max_Score_Len;
   (* T_End) = Max_Score_Len + Max_Score_Best_d;
   Set_Deltas (Delta, Delta_Len, Max_Score_Len, Max_Score_Best_d, Max_Score_Best_e);
   (* Match_To_End) = FALSE;

   return  Max_Score_Best_e;
  }



static int  Read_String
    (FILE * fp, char * * T, long int * Size, char header [])

/* Read next string from  fp  (assuming FASTA format) into  (* T) [1 ..]
*  which has  Size  characters.  Allocate extra memory if needed
*  and adjust  Size  accordingly.  Return  TRUE  if successful,  FALSE
*  otherwise (e.g., EOF).  Set  header  to the contents of the FASTA
*  header line. */

  {
   char  Line [MAX_LINE];
   long int  Len;
   int  Ch, Ct;

   if  (feof (fp))
       return  FALSE;

   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     ;

   if  (Ch != '>')
       return  FALSE;

   fgets (Line, MAX_LINE, fp);
   Len = strlen (Line);
   assert (Len > 0 && Line [Len - 1] == '\n');
   Line [Len - 1] = '\0';
   strcpy (header, Line);


   Ct = 0;
   Len = 0;
   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     {
      if  (isspace (Ch))
          continue;

      Ct ++;

      if  (Ct >= (* Size) - 1)
          {
           (* Size) *= EXPANSION_FACTOR;
           (* T) = (char *) safe_realloc ((* T), (* Size));
          }
      Ch = tolower (Ch);
      if  (! isalpha (Ch))
          {
           fprintf (stderr, "Unexpected character `%c\' in string %s\n",
                                 Ch, header);
           Ch = 'x';
          }
      (* T) [++ Len] = Ch;
     }

   (* T) [Len + 1] = '\0';
   if  (Ch == '>')
       ungetc (Ch, fp);

   return  TRUE;
  }



static void  Rev_Complement
    (char * s)

/* Set string  s  to its DNA reverse complement. */

  {
   char  ch;
   int  i, j, len;

   len = strlen (s);

   for  (i = 0, j = len - 1;  i < j;  i ++, j --)
     {
      ch = Complement (s [i]);
      s [i] = Complement (s [j]);
      s [j] = ch;
     }

   if  (i == j)
       s [i] = Complement (s [i]);

   return;
  }



static void  Rev_Display_Alignment
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct)

//  Show (to  stdout ) the alignment encoded in  delta [0 .. (delta_ct - 1)]
//  between strings  a [0 .. (a_len - 1)]  and  b [0 .. (b_len - 1)]
//  in the reverse direction.

  {
   int  i, j, k, m, top_len, bottom_len;
   char  top [MAX_EXTENSION + 2 * MAX_ERRORS + 1];
   char  bottom [MAX_EXTENSION + 2 * MAX_ERRORS + 1];

   i = j = top_len = bottom_len = 0;
   for  (k = 0;  k < delta_ct;  k ++)
     {
      for  (m = 1;  m < abs (delta [k]);  m ++)
        {
         top [top_len ++] = a [- i ++];
         j ++;
        }
      if  (delta [k] < 0)
          {
           top [top_len ++] = '-';
           j ++;
          }
        else
          {
           top [top_len ++] = a [- i ++];
          }
     }
   while  (i < a_len && j < b_len)
     {
      top [top_len ++] = a [- i ++];
      j ++;
     }
   top [top_len] = '\0';
     

   i = j = 0;
   for  (k = 0;  k < delta_ct;  k ++)
     {
      for  (m = 1;  m < abs (delta [k]);  m ++)
        {
         bottom [bottom_len ++] = b [- j ++];
         i ++;
        }
      if  (delta [k] > 0)
          {
           bottom [bottom_len ++] = '-';
           i ++;
          }
        else
          {
           bottom [bottom_len ++] = b [- j ++];
          }
     }
   while  (j < b_len && i < a_len)
     {
      bottom [bottom_len ++] = b [- j ++];
      i ++;
     }
   bottom [bottom_len] = '\0';


   for  (i = 0;  i < top_len || i < bottom_len;  i += DISPLAY_WIDTH)
     {
      putchar ('\n');
      printf ("A: ");
      for  (j = 0;  j < DISPLAY_WIDTH && i + j < top_len;  j ++)
        putchar (top [i + j]);
      putchar ('\n');
      printf ("B: ");
      for  (j = 0;  j < DISPLAY_WIDTH && i + j < bottom_len;  j ++)
        putchar (bottom [i + j]);
      putchar ('\n');
      printf ("   ");
      for  (j = 0;  j < DISPLAY_WIDTH && i + j < bottom_len && i + j < top_len;
                j ++)
        if  (top [i + j] != ' ' && bottom [i + j] != ' '
                 && top [i + j] != bottom [i + j])
            putchar ('^');
          else
            putchar (' ');
      putchar ('\n');
     }

   return;
  }



static int  Rev_Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Match_To_End,
     int Delta [MAX_ERRORS], int * Delta_Len, int extending)

//  Return the minimum number of changes (inserts, deletes, replacements)
//  needed to match string  A [0 .. -(m-1)]  with a prefix of string
//   T [0 .. -(n-1)]  if it's not more than  Error_Limit .
//  Note:  Match is reverse direction, right to left.
//  Put delta description of alignment in  Delta  and set
//  (* Delta_Len)  to the number of entries there if it's a complete
//  match.
//  Set  A_End  and  T_End  to the rightmost positions where the
//  alignment ended in  A  and  T , respectively.
//  Set  Match_To_End  true if the match extended to the end
//  of at least one string; otherwise, set it false to indicate
//  a branch point.
//  extending  indicates whether the match is for a fixed region (FALSE)
//  or is extending off the end of a match as far as possible (TRUE)

  {
   double  Score, Max_Score;
   int  Max_Score_Len = 0, Max_Score_Best_d = 0, Max_Score_Best_e = 0;
   int  Best_d, Best_e, Longest, Row;
   int  Left, Right;
   int  d, e, i, j, shorter;

//   assert (m <= n);
   Best_d = Best_e = Longest = 0;
   (* Delta_Len) = 0;

   if  (Consec_Non_ACGT > 0)
       {
        int  ct;

        for  (i = ct = 0;  i < m && ct < Consec_Non_ACGT;  i ++)
          {
           switch  (A [- i])
             {
              case  'a' :
              case  'c' :
              case  'g' :
              case  't' :
                ct = 0;
                break;
              default :
                ct ++;
             }
          }
        if  (ct >= Consec_Non_ACGT)
            {
             m = i - ct;
             extending = TRUE;
            }
        
        for  (i = ct = 0;  i < n && ct < Consec_Non_ACGT;  i ++)
          {
           switch  (T [- i])
             {
              case  'a' :
              case  'c' :
              case  'g' :
              case  't' :
                ct = 0;
                break;
              default :
                ct ++;
             }
          }
        if  (ct >= Consec_Non_ACGT)
            {
             n = i - ct;
             extending = TRUE;
            }
       }

   shorter = OVL_Min_int (m, n);
   for  (Row = 0;  Row < shorter && A [- Row] == T [- Row];  Row ++)
     ;

   Edit_Array [0] [0] = Row;

   if  (Row == m && Row == n)      // Exact match
       {
        (* Delta_Len) = 0;
        (* A_End) = (* T_End) = Row;
        (* Match_To_End) = ! extending;
        return  0;
       }

   Left = Right = 0;
   Max_Score = Row * branch_pt_match_value;
   Max_Score_Len = Row;
   Max_Score_Best_d = 0;
   Max_Score_Best_e = 0;

   for  (e = 1;  e <= Error_Limit;  e ++)
     {
      int  cutoff;

      Left = OVL_Max_int (Left - 1, -e);
      Right = OVL_Min_int (Right + 1, e);
      Edit_Array [e - 1] [Left] = -2;
      Edit_Array [e - 1] [Left - 1] = -2;
      Edit_Array [e - 1] [Right] = -2;
      Edit_Array [e - 1] [Right + 1] = -2;

      for  (d = Left;  d <= Right;  d ++)
        {
         Row = 1 + Edit_Array [e - 1] [d];
         if  ((j = Edit_Array [e - 1] [d - 1]) > Row)
             Row = j;
         if  ((j = 1 + Edit_Array [e - 1] [d + 1]) > Row)
             Row = j;
         while  (Row < m && Row + d < n
                  && A [- Row] == T [- Row - d])
           Row ++;

         Edit_Array [e] [d] = Row;

         if  (Row == m && Row + d == n)
             {
              if  (extending)
                  {
                   for  (j = Left;  j <= d;  j ++)
                     if  (Edit_Array [e] [j] > Longest)
                         {
                          Best_d = j;
                          Longest = Edit_Array [e] [j];
                         }
                   Score = Longest * branch_pt_match_value - e;
                   if  (Score < Max_Score)
                       {
                        (* A_End) = Max_Score_Len;
                        (* T_End) = Max_Score_Len + Max_Score_Best_d;
                        Set_Deltas (Delta, Delta_Len, Max_Score_Len,
                                    Max_Score_Best_d, Max_Score_Best_e);
                        (* Match_To_End) = FALSE;
                        return  Max_Score_Best_e;
                       }
                  }

              (* A_End) = Row;           // One past last align position
              (* T_End) = Row + d;
              Set_Deltas (Delta, Delta_Len, Row, d, e);
              (* Match_To_End) = ! extending;
              return  e;
             }
        }

      cutoff = Longest - GIVE_UP_LEN;
      while  (Left <= Right && Edit_Array [e] [Left] < cutoff)
        Left ++;
      if  (Left > Right)
          break;
      while  (Right > Left && Edit_Array [e] [Right] < cutoff)
        Right --;
      assert (Left <= Right);

      for  (d = Left;  d <= Right;  d ++)
        if  (Edit_Array [e] [d] > Longest)
            {
             Best_d = d;
             Best_e = e;
             Longest = Edit_Array [e] [d];
            }
#if  1
      Score = Longest * branch_pt_match_value - e;
               // Assumes  branch_pt_match_value - branch_pt_error_value == 1.0
      if  (Score > Max_Score)
          {
           Max_Score = Score;
           Max_Score_Len = Longest;
           Max_Score_Best_d = Best_d;
           Max_Score_Best_e = Best_e;
          }
#endif

      if  (Longest - Max_Score_Len >= GIVE_UP_LEN)
          break;
     }

   (* A_End) = Max_Score_Len;
   (* T_End) = Max_Score_Len + Max_Score_Best_d;
   Set_Deltas (Delta, Delta_Len, Max_Score_Len, Max_Score_Best_d, Max_Score_Best_e);
   (* Match_To_End) = FALSE;

   return  Max_Score_Best_e;
  }



static void  Rev_Show_Diffs
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct,
     long int a_ref, long int b_ref)

//  Show (to  stdout ) the differences
//  between strings  a [0 .. (a_len - 1)]  and  b [0 .. (b_len - 1)]
//  encoded in  delta [0 .. (delta_ct - 1)]
//  in the reverse direction.   a_ref  is the position of the beginning
//  of string  a .   b_ref  is the offset to the start of
//  string  b .

  {
   int  i, j, k, m;

   i = j = 0;
   for  (k = 0;  k < delta_ct;  k ++)
     {
      for  (m = 1;  m < abs (delta [k]);  m ++)
        {
         if  (a [- i] != b [- j])
             printf ("%8ld %8ld:  %c  %c\n",
                     a_ref - i, b_ref - j, a [- i], b [- j]);
         i ++;
         j ++;
        }
      if  (delta [k] < 0)
          {
           printf ("%8ld %8ld: ins %c\n", a_ref - i, b_ref - j, b [- j]);
           j ++;
          }
        else
          {
           printf ("%8ld %8ld: del %c\n", a_ref - i, b_ref - j, a [- i]);
           i ++;
          }
     }
   while  (i < a_len && j < b_len)
     {
      if  (a [- i] != b [- j])
          printf ("%8ld %8ld:  %c  %c\n",
                  a_ref - i, b_ref - j, a [- i], b [- j]);
      i ++;
      j ++;
     }

   return;
  }



static void  Set_Deltas
    (int delta [], int * delta_len, int row, int d, int e)

//  Set  delta  to values that show the alignment recorded in
//  global  Edit_Array .  Set  (* delta_len)  to the number of
//  entries set.   row  is the length of the match in the first
//  string.   d  is the offset at the end of the match between
//  the first and second string.   e  is the number of errors.

  {
   int  delta_stack [MAX_ERRORS];
   int  from, last, max;
   int  i, j, k;

   last = row;
   (* delta_len) = 0;

   for  (k = e;  k > 0;  k --)
     {
      from = d;
      max = 1 + Edit_Array [k - 1] [d];
      if  ((j = Edit_Array [k - 1] [d - 1]) > max)
          {
           from = d - 1;
           max = j;
          }
      if  ((j = 1 + Edit_Array [k - 1] [d + 1]) > max)
          {
           from = d + 1;
           max = j;
          }
      if  (from == d - 1)
          {
           delta_stack [(* delta_len) ++] = max - last - 1;
           d --;
           last = Edit_Array [k - 1] [from];
          }
      else if  (from == d + 1)
          {
           delta_stack [(* delta_len) ++] = last - (max - 1);
           d ++;
           last = Edit_Array [k - 1] [from];
          }
     }
   delta_stack [(* delta_len) ++] = last + 1;

   k = 0;
   for  (i = (* delta_len) - 1;  i > 0;  i --)
     delta [k ++]
         = abs (delta_stack [i]) * Sign (delta_stack [i - 1]);
   (* delta_len) --;

   return;
  }



static void  Show_Coverage
    (Cover_t * list, char * filename, char * tag, char * seq)

//  Print to stderr summary stats of regions in  list .
//  Send to  filename  a list of region info suitable for  celagram .
//  Use  tag  to print the headings.  Regions refer to  seq .

  {
   FILE  * fp;
   Cover_t  * p;
   int  ct = 0;
   long int  i, prev_hi = 0, n_ct, len, seq_len;
   long int  cov_ns = 0, gap_ns = 0, total_ns;
   double  cov_total = 0.0, gap_total = 0.0;

   fp = Local_File_Open (filename, "w");
   fprintf (fp, "%-9s %9s %9s %8s %8s   %6s\n",
            "Region", "Start", "End", "Len", "N's", "%N's");

   for  (p = list;  p != NULL;  p = p -> next)
     {
      n_ct = 0;
      for  (i = prev_hi + 1;  i < p -> lo;  i ++)
        if  (strchr ("acgt", seq [i]) == NULL)
            n_ct ++;
      len = p -> lo - prev_hi - 1;
      fprintf (fp, "gap%-6d %9ld %9ld %8ld %8ld  %5.1f%%\n",
               ct, prev_hi + 1, p -> lo - 1, len, n_ct,
               len == 0 ? 0.0 : 100.0 * n_ct / len);
      gap_total += len;
      gap_ns += n_ct;

      ct ++;

      n_ct = 0;
      for  (i = p -> lo;  i <= p -> hi;  i ++)
        if  (strchr ("acgt", seq [i]) == NULL)
            n_ct ++;
      len = 1 + p -> hi - p -> lo;
      fprintf (fp, "cov%-6d %9ld %9ld %8ld %8ld  %5.1f%%\n",
               ct, p -> lo, p -> hi, len, n_ct,
               len == 0 ? 0.0 : 100.0 * n_ct / len);
      cov_total += len;
      cov_ns += n_ct;

      prev_hi = p -> hi;
     }
     
   n_ct = 0;
   seq_len = strlen (seq + 1);
   for  (i = prev_hi + 1;  i <= seq_len;  i ++)
     if  (strchr ("acgt", seq [i]) == NULL)
         n_ct ++;
   len = seq_len - prev_hi;
   fprintf (fp, "gap%-6d %9ld %9ld %8ld %8ld  %5.1f%%\n",
            ct, prev_hi + 1, seq_len, len, n_ct,
            len == 0 ? 0.0 : 100.0 * n_ct / len);
   gap_total += len;
   gap_ns += n_ct;

   total_ns = cov_ns + gap_ns;

   fclose (fp);

   fprintf (stderr, "\n%s Sequence Coverage:\n", tag);
   fprintf (stderr, "   Sequence length = %ld\n", seq_len);
   fprintf (stderr, "            acgt's = %ld\n", seq_len - total_ns);
   fprintf (stderr, "        Non acgt's = %ld\n", total_ns);
   fprintf (stderr, " Number of regions = %d\n", ct);
   fprintf (stderr, "     Matched bases = %.0f  (%.1f%%, %.1f%% of acgt's)\n",
            cov_total,
            seq_len == 0.0 ? 0.0 : 100.0 * cov_total / seq_len,
            seq_len - total_ns == 0.0 ? 0.0 :
                100.0 * (cov_total - cov_ns) / (seq_len - total_ns));
   fprintf (stderr, "     Avg match len = %.0f\n",
            ct == 0 ? 0.0 : cov_total / ct);
   fprintf (stderr, "    N's in matches = %ld  (%.1f%%)\n",
            cov_ns,
            cov_total == 0.0 ? 0.0 : 100.0 * cov_ns / cov_total);
   fprintf (stderr, "   Unmatched bases = %.0f  (%.1f%%, %.1f%% of acgt's)\n",
            gap_total,
            seq_len == 0.0 ? 0.0 : 100.0 * gap_total / seq_len,
            seq_len - total_ns == 0.0 ? 0.0 :
                100.0 * (gap_total - gap_ns) / (seq_len - total_ns));
   fprintf (stderr, "       Avg gap len = %.0f\n",
            gap_total / (1.0 + ct));
   fprintf (stderr, "       N's in gaps = %ld  (%.1f%%)\n",
            gap_ns,
            gap_total == 0.0 ? 0.0 : 100.0 * gap_ns / gap_total);

   return;
  }



static void  Show_Diffs
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct,
     long int a_ref, long int b_ref)

//  Show (to  stdout ) the differences
//  between strings  a [0 .. (a_len - 1)]  and  b [0 .. (b_len - 1)]
//  encoded in  delta [0 .. (delta_ct - 1)] .  a_ref  is the offset to
//  the start of string  a .   b_ref  is the offset to the start of
//  string  b .

  {
   int  i, j, k, m;

   i = j = 0;
   for  (k = 0;  k < delta_ct;  k ++)
     {
      for  (m = 1;  m < abs (delta [k]);  m ++)
        {
         if  (a [i] != b [j])
             printf ("%8ld %8ld:  %c  %c\n",
                     a_ref + i, b_ref + j, a [i], b [j]);
         i ++;
         j ++;
        }
      if  (delta [k] < 0)
          {
           printf ("%8ld %8ld: ins %c\n", a_ref + i - 1, b_ref + j, b [j]);
           j ++;
          }
        else
          {
           printf ("%8ld %8ld: del %c\n", a_ref + i, b_ref + j - 1, a [i]);
           i ++;
          }
     }
   while  (i < a_len && j < b_len)
     {
      if  (a [i] != b [j])
          printf ("%8ld %8ld:  %c  %c\n", a_ref + i, b_ref + j, a [i], b [j]);
      i ++;
      j ++;
     }

   return;
  }



static int  Sign
    (int a)

//  Return the algebraic sign of  a .

  {
   if  (a > 0)
       return  1;
   else if  (a < 0)
       return  -1;

   return  0;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s <RefSequence> <MatchSequences> <GapsFile>\n"
           "\n"
           "Combines MUMs in <GapsFile> by extending matches off\n"
           "ends and between MUMs.  <RefSequence> is a fasta file\n"
           "of the reference sequence.  <MatchSequences> is a\n"
           "multi-fasta file of the sequences matched against the\n"
           "reference\n"
           "\n"
           "Options:\n"
           "-n      Allow matches only between nucleotides, i.e., ACGTs\n"
           "-N num  Break matches at <num> or more consecutive non-ACGTs \n"
           "-q tag  Used to label query match\n"
           "-r tag  Used to label reference match\n"
           "-S      Output all differences in strings\n"
           "-t      Label query matches with query fasta header\n"
           "-v num  Set verbose level for extra output\n"
           "-x      Don't output .cover files\n"
	   "-e      Set error-rate cutoff to e (e.g. 0.02 is two percent)\n",
	    
           command);

   return;
  }


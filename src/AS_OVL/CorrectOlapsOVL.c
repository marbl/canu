
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

const char *mainid = "$Id: CorrectOlapsOVL.c,v 1.36 2008-10-28 15:42:10 brianwalenz Exp $";

//  System include files

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <fcntl.h>
#include  <string.h>
#include  <unistd.h>


//  Local include files

#include  "AS_OVL_delcher.h"
#include  "AS_PER_gkpStore.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_UTL_reverseComplement.h"
#include  "AS_UTL_qsort_mt.h"
#include  "FragCorrectOVL.h"
#include  "AS_OVS_overlapStore.h"


//  Constants

#define  BRANCH_PT_MATCH_VALUE    0.272
    //  Value to add for a match in finding branch points
    //  1.20 was the calculated value for 6% vs 35% error discrimination
    //  Converting to integers didn't make it faster
#define  BRANCH_PT_ERROR_VALUE    -0.728
    //  Value to add for a mismatch in finding branch points
    //   -2.19 was the calculated value for 6% vs 35% error discrimination
    //  Converting to integers didn't make it faster
#define  DEFAULT_END_EXCLUDE_LEN     3
    //  Default value for  End_Exclude_Len
#define  DEFAULT_HALF_LEN            4
    //  Default value for bases on each side of SNP to vote for change
#define  DEFAULT_KMER_LEN            9
    //  Default value for  Kmer_Len
#define  DEFAULT_QUALITY_THRESHOLD   0.015
    //  Default value for  Quality_Threshold
#define  EDIT_DIST_PROB_BOUND        1e-4
    //  Probability limit to "band" edit-distance calculation
    //  Determines  NORMAL_DISTRIB_THOLD
#define  ERRORS_FOR_FREE             1
    //  The number of errors that are ignored in setting probability
    //  bound for terminating alignment extensions in edit distance
    //  calculations
#define  MAX_ERROR_RATE              AS_GUIDE_ERROR_RATE
    //  The largest error allowed in overlaps
#define  MAX_FASTA_LINE              2048
    //  Most bytes allowed in line of fasta file
#define  MAX_FILENAME_LEN            1000
    //  Longest name allowed for a file in the overlap store
#define  MAX_FRAG_LEN                2048
    //  The longest fragment allowed
#define  MAX_ERRORS                  (1 + (int) (AS_OVL_ERROR_RATE * MAX_FRAG_LEN))
    //  Most errors in any edit distance computation // 0.40
    //  KNOWN ONLY AT RUN TIME
#define  EXPANSION_FACTOR            1.4
    // Factor by which to grow memory in olap array when reading it
#define  MIN_BRANCH_END_DIST     20
    //  Branch points must be at least this many bases from the
    //  end of the fragment to be reported
#define  MIN_BRANCH_TAIL_SLOPE   0.20
    //  Branch point tails must fall off from the max by at least
    //  this rate
#define  NORMAL_DISTRIB_THOLD    3.62
    //  Determined by  EDIT_DIST_PROB_BOUND
#define  VERBOSE                 0
    //  If  1  will print lots of extra output



//  Type definitions

typedef  int32  Int_Frag_ID_t;

typedef  struct
  {
   int  pos;
   int  adjust;
  }  Adjust_t;

typedef  struct
  {
   char  * sequence;
   Adjust_t  * adjust;
   int  unitig1, lo1, hi1;
   int  unitig2, lo2, hi2;
   unsigned  keep_right : 1;    // set true if right overlap degree is low
   unsigned  keep_left : 1;     // set true if left overlap degree is low
   int16  adjust_ct;
  }  Frag_Info_t;

typedef  struct
  {
   int32   a_iid, b_iid;
   int16   a_hang, b_hang;
   int32   place;                // position in array before sort
   uint16  corr_erate;
   char    orient;
  }  Olap_Info_t;

typedef  struct
  {
   unsigned  lo_iid : 31;
   unsigned  hi_iid : 31;
   unsigned  confirmed : 1;
  }  Frag_Pair_t;

typedef  struct
  {
   int  * list;
   int  len;
  }  Int_List_t;



//  Static Globals

static char  * Correct_File_Path;
    // Name of file containing fragment corrections
static FILE  * Delete_fp = NULL;
    // File to which list of overlaps to delete is written if  -x  option is specified
static int  * Edit_Array [AS_FRAG_MAX_LEN];
    // Use for alignment calculation.  Points into  Edit_Space .
    // (only MAX_ERRORS needed)
static int  Edit_Match_Limit [AS_FRAG_MAX_LEN] = {0};
    // This array [e] is the minimum value of  Edit_Array [e] [d]
    // to be worth pursuing in edit-distance computations between guides
    // (only MAX_ERRORS needed)
static int  Edit_Space [(AS_FRAG_MAX_LEN + 4) * AS_FRAG_MAX_LEN];
    // Memory used by alignment calculation
    // (only (MAX_ERRORS + 4) * MAX_ERRORS needed)
static int  End_Exclude_Len = DEFAULT_END_EXCLUDE_LEN;
    // Length of ends of exact-match regions not used in preventing
    // sequence correction
static char  * Erate_Path = NULL;
    // Name of binary file to which to dump revised error rates
    // This allows the program to run in parallel (e.g., under LSF)
    // Presumably the error-rate files will be added to the overlap
    // store later by the  update-erates  program.
static int  Error_Bound [MAX_FRAG_LEN + 1];
    // This array [i]  is the maximum number of errors allowed
    // in a match between sequences of length  i , which is
    //  i * MAXERROR_RATE .
static int  Failed_Alignments_Ct = 0;
    // Count the number of alignments that failed
static Frag_Info_t  * Frag;
    // Sequence and vote information for current range of fragments
    // being corrected
static int32  Frags_Per_File;
    // Maximum number of fragments in each data file of fragment store.
    // This is read from the store
static GateKeeperStore  *gkpStore;
    // Internal fragment store where fragments are loaded
static FragStream  *Frag_Stream;
    // Stream to extract fragments from internal store
static char  * gkpStore_Path;
    // Name of directory containing fragment store from which to get fragments
static int  Half_Len = DEFAULT_HALF_LEN;
    // Number of bases on each side of SNP to vote for change
static int32  Hi_Frag_IID;
    // Internal ID of last fragment in frag store to process
static int  Hi_Unitig = -1;
    // The highest numbered unitig.
static unsigned  Highest_Frag = 0;
    // The highest numbered fragment in the unitig messages.
static int  * IUM = NULL;
    // Has unitig ID for each fragment.
static int  IUM_Size = 0;
    // Number of entries in  IUM .
static Int_List_t  * Keep_Pair = NULL;
    // Array holding pairs of unitigs that corrected overlaps imply
    // should overlap.
static int  Kmer_Len = DEFAULT_KMER_LEN;
    // Length of minimum exact match in overlap to confirm base pairs
static int32  Lo_Frag_IID;
    // Internal ID of first fragment in frag store to process
static int  Num_Frags = 0;
    // Number of fragments being corrected
static uint64  Num_Olaps;
    // Number of overlaps being used
static Olap_Info_t  * Olap = NULL;
    // Array of overlaps being used
static uint32  * Olap_Offset = NULL;
    // Indicates the first overlap of each fragment
static char  * Olap_Path;
    // Name of file containing a sorted list of overlaps
static int  Olaps_From_Store = FALSE;
    // Indicates if overlap info comes from  get-olaps  or from
    // a binary overlap store
static FILE  * OVL_fp = NULL;
    // File to which OVL messages are written if  -o  option is specified
static double  Quality_Threshold = DEFAULT_QUALITY_THRESHOLD;
    // Overlaps better than this error rate will be output
static int  Total_Alignments_Ct = 0;
    // Count the number of alignments attempted
static int  * UF = NULL;
    // For Union-Find data structure for sets of fragments.


//  Static Functions

static void  Add_To_List
    (Int_List_t * a, int b);
static void  Apply_Seq_Corrects
    (char * * seq, Adjust_t * * adjust, int16 * adjust_ct,
     Correction_t correct [], int n, int fixed);
static int  Binomial_Bound
    (int, double, int, double);
static int  By_B_IID
    (const void * a, const void * b);
static int  By_Lo_IID
    (const void * a, const void * b);
static int  By_Place
    (const void * a, const void * b);
static int  Compare_Frags
    (char a [], char b []);
static void  Correct_Frags
    (void);
static void  Display_Alignment
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct);
static void  Display_Frags
    (void);
static void  Dump_Erate_File
    (char * path, int32 lo_id, int32 hi_id, Olap_Info_t * olap, uint64 num);
static void  Fasta_Print
    (FILE * fp, char * s, char * hdr);
static char  Filter
    (char ch);
static int  Find
    (int w, int a []);
static void  Get_Canonical_Olap_Region
    (Olap_Info_t * olap, int sub, char * a_seq, char * b_seq,
     Adjust_t forw_adj [], int adj_ct,
     int frag_len, char * * a_part, char * * b_part);
static void  Get_Olaps_From_Store
    (char * path, int32 lo_id, int32 hi_id, Olap_Info_t * * olap, uint64 * num);
static int  Hang_Adjust
    (int hang, Adjust_t adjust [], int adjust_ct);
static void  Initialize_Globals
    (void);
static int  Intersect_Len
    (int a, int b, int c, int d);
static void  Keep_Olap
    (Olap_Info_t * olap);
static void  Make_Rev_Adjust
    (Adjust_t rev_adj [], Adjust_t forw_adj [], int adj_ct, int frag_len);
static int  OVL_Max_int
    (int a, int b);
static int  OVL_Min_int
    (int a, int b);
static int  Olap_In_Unitig
    (Olap_Info_t * olap);
static void  Output_Delete_OVLs
    (void);
static int  Output_OVL
    (Olap_Info_t * olap, double quality);
static void  Parse_Command_Line
    (int argc, char * argv []);
static int  Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Match_To_End,
     int * Delta, int * Delta_Len);
static void  Process_Olap
    (Olap_Info_t * olap, char * b_seq, Adjust_t forw_adj [], int adj_ct,
     int frag_len);
static char *  Read_Fasta
    (FILE * fp);
static void  Read_Frags
    (void);
static void  Read_Olaps
    (void);
static void  Redo_Olaps
    (void);
static int  Rev_Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Match_To_End,
     int * Delta, int * Delta_Len);
static int  Sign
    (int a);
static int  Union
    (int i, int j, int a []);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   int  i;

   Parse_Command_Line  (argc, argv);

   Initialize_Globals ();

   fprintf (stderr, "Starting Read_Frags ()\n");
   Read_Frags ();

//**ALD
//printf ("Before corrections applied:\n");
//   Display_Frags ();

   fprintf (stderr, "Starting Correct_Frags ()\n");
   Correct_Frags ();

//**ALD
//printf ("After corrections applied:\n");
//   Display_Frags ();

   fprintf (stderr, "Starting Read_Olaps ()\n");
   Read_Olaps ();

   if (Num_Olaps != 0)
   {
     fprintf (stderr, "Starting qsort ()\n");
     qsort (Olap, Num_Olaps, sizeof (Olap_Info_t), By_B_IID);

     fprintf (stderr, "Starting Redo_Olaps ()\n");
     Redo_Olaps ();
   }
   else
   {
     fprintf(stderr, "No overlaps in fragment range, skipping Redo_Olaps()\n");
   }

   if  (OVL_fp != NULL)
       fclose (OVL_fp);

   if  (Olaps_From_Store)
       {
        fprintf (stderr, "Starting re-sort\n");
        qsort (Olap, Num_Olaps, sizeof (Olap_Info_t), By_Place);
        if  (Erate_Path != NULL)
            {
             fprintf (stderr, "Saving corrected error rates to file %s\n",
                      Erate_Path);
             Dump_Erate_File (Erate_Path, Lo_Frag_IID, Hi_Frag_IID, Olap, Num_Olaps);
            }
       }

   fprintf (stderr, "%d/%d failed/total alignments (%.1f%%)\n",
            Failed_Alignments_Ct, Total_Alignments_Ct,
            Total_Alignments_Ct == 0 ? 0.0
                : (100.0 * Failed_Alignments_Ct) / Total_Alignments_Ct);
   fprintf (stderr, "Finished\n");

   return  0;
  }



static void  Add_To_List
    (Int_List_t * a, int b)

//  Add  b  to list  (* a)  if it's not there already.

  {
   int  i;

   for  (i = 0;  i < a -> len;  i ++)
     if  (b == a -> list [i])
         return;

   a -> len ++;
   a -> list = (int *) safe_realloc (a -> list, a -> len * sizeof (int));
   a -> list [a -> len - 1] = b;

   return;
  }



static void  Apply_Seq_Corrects
    (char * * seq, Adjust_t * * adjust, int16 * adjust_ct,
     Correction_t correct [], int n, int fixed)

//  Apply the corrections in  correct [0 .. (n - 1)]  to
//   sequence (* seq) .  Set  (* adjust)  to values
//  needed to adjust offsets for overlaps because of indels
//  in corrections.  If  fixed  is true, then apply corrections
//  in existing sequence space (which is assumed to be large enough);
//  otherwise, realloc a new string for the revised sequence.

  {
   char  buff [2 * MAX_FRAG_LEN];
   Adjust_t  adj_buff [2 * MAX_FRAG_LEN];
   int  adj_val, len, new_len;
   int  i, j, k, ct;

   if  ((* seq) == NULL)
       { // Skip this sequence
        (* adjust_ct) = 0;
        return;
       }
   len = strlen (* seq);
   ct = adj_val = 0;

   i = j = k = 0;
   while  (i < len)
     {
      if  (k >= n || i < correct [k] . pos)
          buff [j ++] = (* seq) [i ++];
      else
        {
         assert (i == correct [k] . pos || i == correct [k] . pos + 1);
         switch ((Vote_Value_t) correct [k] . type)
           {
            case  DELETE :
              // skip  (* seq) [i]
              adj_buff [ct] . pos = i + 1;
              adj_buff [ct ++] . adjust = -- adj_val;
              i ++;
              break;
            case  A_SUBST :
              buff [j ++] = 'a';
              i ++;
              break;
            case  C_SUBST :
              buff [j ++] = 'c';
              i ++;
              break;
            case  G_SUBST :
              buff [j ++] = 'g';
              i ++;
              break;
            case  T_SUBST :
              buff [j ++] = 't';
              i ++;
              break;
            case  A_INSERT :
              if  (i != correct [k] . pos + 1)   // Insert not immediately after subst
                  buff [j ++] = (* seq) [i ++];
              adj_buff [ct] . pos = i + 1;
              buff [j ++] = 'a';
              adj_buff [ct ++] . adjust = ++ adj_val;
              break;
            case  C_INSERT :
              if  (i != correct [k] . pos + 1)   // Insert not immediately after subst
                  buff [j ++] = (* seq) [i ++];
              adj_buff [ct] . pos = i + 1;
              buff [j ++] = 'c';
              adj_buff [ct ++] . adjust = ++ adj_val;
              break;
            case  G_INSERT :
              if  (i != correct [k] . pos + 1)   // Insert not immediately after subst
                  buff [j ++] = (* seq) [i ++];
              adj_buff [ct] . pos = i + 1;
              buff [j ++] = 'g';
              adj_buff [ct ++] . adjust = ++ adj_val;
              break;
            case  T_INSERT :
              if  (i != correct [k] . pos + 1)   // Insert not immediately after subst
                  buff [j ++] = (* seq) [i ++];
              adj_buff [ct] . pos = i + 1;
              buff [j ++] = 't';
              adj_buff [ct ++] . adjust = ++ adj_val;
              break;
            default :
              fprintf (stderr, "ERROR:  Illegal vote type\n");
           }
         k ++;
        }
     }

   buff [j] = '\0';
   new_len = j;

   if  (! fixed)
       {
        int  size;

        (* seq) = (char *) safe_realloc ((* seq), 1 + new_len);
        size = OVL_Max_int (ct, 1);   // Prevent realloc failures
        (* adjust) = (Adjust_t *) safe_realloc ((* adjust), size * sizeof (Adjust_t));
       }
   strcpy ((* seq), buff);
   for  (i = 0;  i < ct;  i ++)
     (* adjust) [i] = adj_buff [i];
   (* adjust_ct) = ct;

#if  0
fprintf (stderr,
         "Apply_Seq_Corrects:  len = %d  new_len = %d  n = %d\n",
         len, new_len, n);
#endif

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






static int  By_B_IID
    (const void * a, const void * b)

//  Compare the values in  a  and  b  as  (* Olap_Info_t) 's,
//  first by  b_iid , then by  a_iid.
//  Return  -1  if  a < b ,  0  if  a == b , and  1  if  a > b .
//  Used for  qsort .

  {
   Olap_Info_t  * x, * y;

   x = (Olap_Info_t *) a;
   y = (Olap_Info_t *) b;

   if  (x -> b_iid < y -> b_iid)
       return  -1;
   else if  (x -> b_iid > y -> b_iid)
       return  1;
   else if  (x -> a_iid < y -> a_iid)
       return  -1;
   else if  (x -> a_iid > y -> a_iid)
       return  1;

   return  0;
  }



static int  By_Lo_IID
    (const void * a, const void * b)

//  Compare the values in  a  and  b  as  (* Frag_Pair_t) 's,
//  first by  lo_iid , then by  hi_iid.
//  Return  -1  if  a < b ,  0  if  a == b , and  1  if  a > b .
//  Used for  qsort .

  {
   Frag_Pair_t  * x, * y;

   x = (Frag_Pair_t *) a;
   y = (Frag_Pair_t *) b;

   if  (x -> lo_iid < y -> lo_iid)
       return  -1;
   else if  (x -> lo_iid > y -> lo_iid)
       return  1;
   else if  (x -> hi_iid < y -> hi_iid)
       return  -1;
   else if  (x -> hi_iid > y -> hi_iid)
       return  1;

   return  0;
  }



static int  By_Place
    (const void * a, const void * b)

//  Compare the values in  a  and  b  as  (* Olap_Info_t) 's,
//  by  place  field.
//  Return  -1  if  a < b ,  0  if  a == b , and  1  if  a > b .
//  Used for  qsort .

  {
   Olap_Info_t  * x, * y;

   x = (Olap_Info_t *) a;
   y = (Olap_Info_t *) b;

   if  (x -> place < y -> place)
       return  -1;
   else if  (x -> place > y -> place)
       return  1;

   return  0;
  }



static int  Compare_Frags
    (char a [], char b [])

//  Display alignment between strings  a  and  b .

  {
   int delta [AS_FRAG_MAX_LEN];  //  only MAX_ERRORS needed
   int  a_end, b_end, delta_len, olap_len, match_to_end;
   int  a_len, b_len, errors;

   a_len = strlen (a);
   b_len = strlen (b);
   olap_len = OVL_Min_int (a_len, b_len);

   errors = Prefix_Edit_Dist
              (a, a_len, b, b_len,
               Error_Bound [olap_len], & a_end, & b_end,
               & match_to_end, delta, & delta_len);

   if  (Verbose_Level > 0)
       {
        printf ("  errors = %d  delta_len = %d\n", errors, delta_len);
        Display_Alignment (a, a_len, b, b_len, delta, delta_len);
       }

   return  errors;
  }



static void  Correct_Frags
    (void)

//  Open and read corrections from  Correct_File_Path  and
//  apply them to sequences in  Frag .

  {
   FILE  * fp;
   Correction_Output_t  msg;
   Correction_t  correct [MAX_FRAG_LEN];
   int  before_errors = 0, after_errors;
   int  num_corrects = 0;
   int  correcting = FALSE;
   uint32  iid = 0;

   fp = File_Open (Correct_File_Path, "rb");

   while  (fread (& msg, sizeof (Correction_Output_t), 1, fp) == 1)
     {
      if  (msg . frag . is_ID)
          {
           if  (correcting)
               {
                Apply_Seq_Corrects (& (Frag [iid - Lo_Frag_IID] . sequence),
                                    & (Frag [iid - Lo_Frag_IID] . adjust),
                                    & (Frag [iid - Lo_Frag_IID] . adjust_ct),
                                    correct, num_corrects, FALSE);
               }

           iid = msg . frag . iid;
           if  (iid < Lo_Frag_IID)
               continue;
           else if  (iid > Hi_Frag_IID)
               {
                correcting = FALSE;
                break;
               }
             else
               {
                correcting = TRUE;
                num_corrects = 0;
                Frag [iid - Lo_Frag_IID] . keep_left = msg . frag . keep_left;
                Frag [iid - Lo_Frag_IID] . keep_right = msg . frag . keep_right;
               }
          }
      else if  (correcting)
          correct [num_corrects ++] = msg . corr;
     }

   if  (correcting)
       {
        Apply_Seq_Corrects (& (Frag [iid - Lo_Frag_IID] . sequence),
                            & (Frag [iid - Lo_Frag_IID] . adjust),
                            & (Frag [iid - Lo_Frag_IID] . adjust_ct),
                            correct, num_corrects, FALSE);
       }

   fclose (fp);
   return;
  }



#define  DISPLAY_WIDTH   60

static void  Display_Alignment
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct)

//  Show (to  stdout ) the alignment encoded in  delta [0 .. (delta_ct - 1)]
//  between strings  a [0 .. (a_len - 1)]  and  b [0 .. (b_len - 1)] .

  {
   int  i, j, k, m, top_len, bottom_len;
   char  top [2000], bottom [2000];

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



static void  Display_Frags
    (void)

//  List selected fragments in fasta format to stdout

  {
   int  i;

   for  (i = 0;  i < Num_Frags;  i ++)
     {
      if (Frag [i] . sequence != NULL)
        {
          int  j, ct;

          printf (">%d\n", Lo_Frag_IID + i);
          ct = 0;
          for  (j = 0;  Frag [i] . sequence [j] != '\0';  j ++)
            {
              if  (ct == 60)
                {
                  putchar ('\n');
                  ct = 0;
                }
              putchar (Frag [i] . sequence [j]);
              ct ++;
            }
          putchar ('\n');
        }
     }

   return;
  }



static void  Dump_Erate_File
    (char * path, int32 lo_id, int32 hi_id, Olap_Info_t * olap, uint64 num)

//  Create a binary file of new error rates in  path .  The format
//  is  lo_id , then  hi_id , then  num , followed by an array of  num
//  error rates.

  {
   FILE  * fp;
   int32   header [2];
   uint16  * erate = NULL;
   int  i;

   fp = File_Open (path, "wb");

   header [0] = lo_id;
   header [1] = hi_id;
   Safe_fwrite (header, sizeof (int32), 2, fp);
	Safe_fwrite (&num,   sizeof (uint64), 1, fp);

   erate = (uint16 *) safe_malloc (num * sizeof(uint16));
   for  (i = 0;  i < num;  i ++)
     erate [i] = olap [i] . corr_erate;

   Safe_fwrite (erate, sizeof (uint16), num, fp);

   safe_free(erate);

   fclose (fp);

   return;
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



static int  Find
    (int w, int a [])

//  Return the number of the set containing  w  in the Union-Find
//  array  a [] .

  {
   int  i, j, k;

   for  (i = w;  a [i] >= 0;  i = a [i])
     ;

   for  (j = w;  a [j] > 0;  j = k)
     {
      k = a [j];
      a [j] = i;
     }

   return  i;
  }



static void  Get_Canonical_Olap_Region
    (Olap_Info_t * olap, int sub, char * a_seq, char * b_seq,
     Adjust_t forw_adj [], int adj_ct,
     int frag_len, char * * a_part, char * * b_part)

//  Set  (* a_part)  and  (* b_part)  to the start of the region
//  to be aligned for the overlap in  (* olap) .   a_seq  is the
//  a-fragment sequence.   b_seq  is the b-fragment sequence.
//  forw_adj [0 .. (adj_ct - 1)]  has adjustment values caused by
//  corrections in the B sequence in the forward orientation.
//  frag_len  is the length of the B sequence.   sub  is the subscript
//  of the a-fragment in the global  Frag  array.

  {
   static int32  b_rev_id = -1;
   static char  b_rev_seq [AS_READ_MAX_LEN + 1];
   static Adjust_t  b_rev_adj [MAX_FRAG_LEN];
#if 0
   static char  a_rev_seq [AS_READ_MAX_LEN + 1];
   static Adjust_t  a_rev_adj [MAX_FRAG_LEN];
#endif

//   if  (olap -> a_iid < olap -> b_iid)
       {
        (* a_part) = a_seq;
        if  (olap -> a_hang > 0)
            (* a_part) += Hang_Adjust (olap -> a_hang, Frag [sub] . adjust,
                                       Frag [sub] . adjust_ct);

        if  (olap -> orient == 'N')
            (* b_part) = b_seq;
          else
            {
             if  (b_rev_id != olap -> b_iid)
                 {
                  strcpy (b_rev_seq, b_seq);
                  reverseComplementSequence (b_rev_seq, 0);
                  b_rev_id = olap -> b_iid;
                  Make_Rev_Adjust (b_rev_adj, forw_adj, adj_ct, frag_len);
                 }
             (* b_part) = b_rev_seq;
            }

        if  (olap -> a_hang < 0)
            {
             int  adjustment;

             if  (olap -> orient == 'N')
                 adjustment = Hang_Adjust (- olap -> a_hang,
                                           forw_adj, adj_ct);
               else
                 adjustment = Hang_Adjust (- olap -> a_hang,
                                           b_rev_adj, adj_ct);
             (* b_part) += adjustment;
            }
       }
#if  0   // new part that I'm adding
     else
       {
        (* a_part) = b_seq;

        if  (olap -> orient == 'N')
            {
             (* b_part) = a_seq;
             if  (olap -> a_hang < 0)
                 (* a_part) += Hang_Adjust (- olap -> a_hang, forw_adj, adj_ct);
             else if  (olap -> a_hang > 0)
                 (* b_part) += Hang_Adjust (olap -> a_hang, Frag [sub] . adjust,
                                            Frag [sub] . adjust_ct);
            }
          else
            {
             strcpy (a_rev_seq, a_seq);
             reverseComplement (a_rev_seq);
             (* b_part) = a_rev_seq;
             if  (olap -> b_hang < 0)
                 {
                  Make_Rev_Adjust (a_rev_adj, Frag [sub] . adjust,
                                   Frag [sub] . adjust_ct, strlen (a_seq));
// stopped here
                 }
            }

       }
#endif

   if  (Verbose_Level > 0)
       {
        Fasta_Print (stdout, (* a_part), "a_part");
        Fasta_Print (stdout, (* b_part), "b_part");
       }

   return;
  }



static void  Get_Olaps_From_Store
    (char * path, int32 lo_id, int32 hi_id, Olap_Info_t * * olap, uint64 * num)

//  Open overlap store  path  and read from it the overlaps for fragments
//   lo_id .. hi_id , putting them in  (* olap)  for which space
//  is dynamically allocated.  Set  (* num)  to the number of entries
//  in  (* olap) .

  {
    OverlapStore  *ovs = NULL;
    OVSoverlap     ovl;
    uint64         numolaps = 0;
    uint64         numread  = 0;

    assert (1 <= lo_id && lo_id <= hi_id);

    ovs = AS_OVS_openOverlapStore(path);

    AS_OVS_setRangeOverlapStore(ovs, lo_id, hi_id);

    numolaps = AS_OVS_numOverlapsInRange(ovs);

    *olap = (Olap_Info_t *)safe_malloc(numolaps * sizeof(Olap_Info_t));
    *num  = 0;

    while (AS_OVS_readOverlapFromStore(ovs, &ovl, AS_OVS_TYPE_OVL)) {
      (*olap)[numread].a_iid  = ovl.a_iid;
      (*olap)[numread].b_iid  = ovl.b_iid;
      (*olap)[numread].a_hang = ovl.dat.ovl.a_hang;
      (*olap)[numread].b_hang = ovl.dat.ovl.b_hang;

      (*olap)[numread].place = numread;

      (*olap)[numread].corr_erate = ovl.dat.ovl.corr_erate;

      (*olap)[numread].orient = 'N';
      if  (ovl.dat.ovl.flipped)
        (*olap)[numread].orient = 'I';

      numread++;
    }

		// numread is defined as uint64, we were assigning it to int before, which would truncate
    (*num) = numread;
  }



static int  Hang_Adjust
    (int hang, Adjust_t adjust [], int adjust_ct)

//  Return the adjusted value of  hang  based on
//   adjust [0 .. (adjust_ct - 1)] .

  {
   int  i, delta = 0;

   for  (i = 0;  i < adjust_ct && hang >= adjust [i] . pos;  i ++)
     delta = adjust [i] . adjust;

   return  hang + delta;
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
     Error_Bound [i] = (int) (i * MAX_ERROR_RATE);

   return;
  }



static int  Intersect_Len
    (int a, int b, int c, int d)

//  Return the number of bases in the intersection between
//  ranges  a .. b  and  c .. d .  Note that maybe  b < a  and
//  d < c  and that coordinates are "gap-based".

  {
   int  save;

   if  (b < a)
       {
        save = a;
        a = b;
        b = save;
       }
   if  (d < c)
       {
        save = c;
        c = d;
        d = save;
       }

   if  (b <= c || d <= a)
       return  0;

   return  OVL_Min_int (b, d) - OVL_Max_int (a, c);
  }



static void  Keep_Olap
    (Olap_Info_t * olap)

//  Check if the pair of unitigs whose overlap is implied by  olap
//  is already in the global  Keep_Pair  list.

  {
   int  a_uni, b_uni, error = FALSE;

   a_uni = IUM [olap -> a_iid];
   b_uni = IUM [olap -> b_iid];

   if  (a_uni < 0)
       {
        fprintf (stderr, "No unitig for frag %d\n", olap -> a_iid);
        error = TRUE;
       }
   if  (b_uni < 0)
       {
        fprintf (stderr, "No unitig for frag %d\n", olap -> b_iid);
        error = TRUE;
       }

   if  (a_uni == b_uni || error)
       return;

   if  (a_uni < b_uni)
       Add_To_List (Keep_Pair + a_uni, b_uni);
     else
       Add_To_List (Keep_Pair + b_uni, a_uni);

   return;
  }



static void  Make_Rev_Adjust
    (Adjust_t rev_adj [], Adjust_t forw_adj [], int adj_ct, int frag_len)

//  Set hanging offset values for reversed fragment in
//   rev_adj [0 .. (adj_ct - 1)]  based on corresponding forward
//  values in  forw_adj [0 .. (adj_ct - 1)] .   frag_len  is the length
//  of the fragment.

  {
   int  i, j, prev;

   if  (adj_ct == 0)
       return;

   j = prev = 0;
   for  (i = adj_ct - 1;  i > 0;  i --)
     {
      if  (forw_adj [i] . adjust == forw_adj [i - 1] . adjust + 1)
          {
           rev_adj [j] . pos = 2 + frag_len - forw_adj [i] . pos;
           rev_adj [j] . adjust = prev + 1;
           prev = rev_adj [j] . adjust;
          }
      else if  (forw_adj [i] . adjust == forw_adj [i - 1] . adjust - 1)
          {
           rev_adj [j] . pos = 3 + frag_len - forw_adj [i] . pos;
           rev_adj [j] . adjust = prev - 1;
           prev = rev_adj [j] . adjust;
          }
        else
          {
           fprintf (stderr, "ERROR:  Bad adjustment value\n");
           fprintf (stderr,
                    "  i = %d  adj_ct = %d  adjust [i] = %d  adjust [i - 1] = %d\n",
                    i, adj_ct, forw_adj [i] . adjust,
                    forw_adj [i - 1] . adjust);
           exit (EXIT_FAILURE);
          }
      j ++;
     }

   if  (forw_adj [i] . adjust == 1)
       {
        rev_adj [j] . pos = 2 + frag_len - forw_adj [i] . pos;
        rev_adj [j] . adjust = prev + 1;
       }
   else if  (forw_adj [i] . adjust == -1)
       {
        rev_adj [j] . pos = 3 + frag_len - forw_adj [i] . pos;
        rev_adj [j] . adjust = prev - 1;
       }
     else
       {
        fprintf (stderr, "ERROR:  Bad adjustment value\n");
        fprintf (stderr,
                 "  i = %d  adj_ct = %d  adjust [i] = %d\n",
                 i, adj_ct, forw_adj [i] . adjust);
        exit (EXIT_FAILURE);
       }

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



static int  Olap_In_Unitig
    (Olap_Info_t * olap)

//  Return  TRUE  iff olap is consistent with either of the two
//  unitigs stored in  Frag .

  {
   int  a_sub, b_sub;
   int  a, b, c, d, p, q;

   a_sub = olap -> a_iid - Lo_Frag_IID;
   b_sub = olap -> b_iid - Lo_Frag_IID;

   a = Frag [a_sub] . lo1;
   b = Frag [a_sub] . hi1;
   p = Frag [a_sub] . unitig1;
   c = Frag [b_sub] . lo1;
   d = Frag [b_sub] . hi1;
   q = Frag [b_sub] . unitig1;

   if  (p == q && c <= b && a <= d)
       return  TRUE;

   a = Frag [a_sub] . lo2;
   b = Frag [a_sub] . hi2;
   p = Frag [a_sub] . unitig2;
   c = Frag [b_sub] . lo2;
   d = Frag [b_sub] . hi2;
   q = Frag [b_sub] . unitig2;

   if  (p == q && c <= b && a <= d)
       return  TRUE;

   return  FALSE;
  }




static int  Output_OVL
    (Olap_Info_t * olap, double quality)

//  Output an OVL message for  olap  with error rate  quality
//  if it's the canonical order.  Return  TRUE  if output is done;
//  otherwise, return  FALSE .

  {
   OverlapMesg  ovMesg;
   GenericMesg  outputMesg;
   Frag_Pair_t  pair, * found;

   switch  (olap -> orient)
     {
      case  'N' :
        if  (olap -> a_hang < 0
             || (olap -> a_hang == 0 && olap -> b_hang > 0))
            return  FALSE;     // not canonical

        ovMesg . orientation = AS_NORMAL;
        ovMesg . ahg = olap -> a_hang;
        ovMesg . bhg = olap -> b_hang;
        break;

      case  'I' :
        if  (olap -> a_hang >= 0)
            {
              if  ((olap -> b_hang >= 0 && olap -> a_iid < olap -> b_iid)
                   || (olap -> a_hang == 0 && olap -> b_hang > 0))
                 return  FALSE;     // not canonical

             ovMesg . orientation = AS_INNIE;
             ovMesg . ahg = olap -> a_hang;
             ovMesg . bhg = olap -> b_hang;
            }
          else
            {
             if  (olap -> b_hang >= 0)
                 return  FALSE;     // not canonical
             if  (olap -> a_iid > olap -> b_iid)
                 return  FALSE;     // not canonical

             ovMesg . ahg = - olap -> b_hang;
             ovMesg . bhg = - olap -> a_hang;
             ovMesg . orientation = AS_OUTTIE;
            }
        break;

      default :
        fprintf (stderr, "ERROR:  Bad orientation = \'%d\'\n",
                 olap -> orient);
        assert (FALSE);
     }

   if  (OVL_fp != NULL)
       {
        outputMesg . m = & ovMesg;
        outputMesg . t = MESG_OVL;
        ovMesg . alignment_trace = NULL;
        ovMesg . aifrag = (Int_Frag_ID_t) olap -> a_iid;
        ovMesg . bifrag = (Int_Frag_ID_t) olap -> b_iid;
        if  (ovMesg . bhg <= 0)
            ovMesg . overlap_type = AS_CONTAINMENT;
        else
            ovMesg . overlap_type = AS_DOVETAIL;
        ovMesg . quality = quality;
        ovMesg . min_offset = ovMesg . max_offset = ovMesg . ahg;   // kludge for now
        ovMesg . polymorph_ct = 0;

        WriteProtoMesg_AS (OVL_fp, & outputMesg);
       }

   return  TRUE;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   FILE  * fp;
   int  ch, errflg = FALSE;
   char  * p;

   argc = AS_configure(argc, argv);

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "e:F:o:Pq:S:v:X:")) != EOF))
     switch  (ch)
       {
        case  'e' :
          Erate_Path = optarg;
          break;

        case  'F' :
          Olap_Path = optarg;
          break;

        case  'o' :
          OVL_fp = File_Open (optarg, "w");
          break;

        case  'q' :
          Quality_Threshold = strtod (optarg, NULL);
          break;

        case  'S' :
          Olap_Path = optarg;
          Olaps_From_Store = TRUE;
          break;

        case  'v' :
          Verbose_Level = (int) strtol (optarg, & p, 10);
          fprintf (stderr, "Verbose level set to %d\n", Verbose_Level);
          break;

        case  'X' :
          Delete_fp = File_Open (optarg, "w");
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc - 4)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   if  (Olap_Path == NULL)
       {
        fprintf (stderr, "ERROR:  Must specify overlaps with -F or -S\n");
        exit (EXIT_FAILURE);
       }

   gkpStore_Path = argv [optind ++];

   Correct_File_Path = argv [optind ++];

   Lo_Frag_IID = (int) strtol (argv [optind], & p, 10);
   if  (p == optarg || Lo_Frag_IID < 1)
       {
        fprintf (stderr, "ERROR:  Illegal low fragment IID \"%s\"\n",
                 argv [optind]);
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }
   optind ++;

   if  (strcmp (argv [optind], "end") == 0)
       {
        Hi_Frag_IID = INT_MAX;
        p = NULL;
       }
     else
       Hi_Frag_IID = (int) strtol (argv [optind], & p, 10);
   if  (p == argv [optind] || Hi_Frag_IID < Lo_Frag_IID)
       {
        fprintf (stderr, "ERROR:  Illegal high fragment IID \"%s\"\n",
                 argv [optind]);
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   fprintf (stderr, "Quality Threshold = %.2f%%\n", 100.0 * Quality_Threshold);

   return;
  }



static int  Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Match_To_End,
     int * Delta, int * Delta_Len)

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

  {
   int  Delta_Stack [AS_FRAG_MAX_LEN];  //  only MAX_ERRORS needed
   double  Score, Max_Score;
   int  Max_Score_Len, Max_Score_Best_d, Max_Score_Best_e;
   int  Best_d, Best_e, From, Last, Longest, Max, Row;
   int  Left, Right;
   int  d, e, i, j, k, shorter;

//   assert (m <= n);
   Best_d = Best_e = Longest = 0;
   (* Delta_Len) = 0;

   shorter = OVL_Min_int (m, n);
   for  (Row = 0;  Row < shorter && A [Row] == T [Row];  Row ++)
     ;

   Edit_Array [0] [0] = Row;

   if  (Row == shorter)                              // Exact match
       {
        (* A_End) = (* T_End) = Row;
        (* Match_To_End) = TRUE;
        return  0;
       }

   Left = Right = 0;
   Max_Score = 0.0;
   Max_Score_Len = Max_Score_Best_d = Max_Score_Best_e = 0;
   for  (e = 1;  e <= Error_Limit;  e ++)
     {
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

         if  (Row == m || Row + d == n)
             {
#if  1
              // Force last error to be mismatch rather than insertion
              if  (Row == m
                     && 1 + Edit_Array [e - 1] [d + 1]
                          == Edit_Array [e] [d]
                     && d < Right)
                  {
                   d ++;
                   Edit_Array [e] [d] = Edit_Array [e] [d - 1];
                  }
#endif
              (* A_End) = Row;           // One past last align position
              (* T_End) = Row + d;
              // Compute Delta
              Last = Row;
              (* Delta_Len) = 0;
              for  (k = e;  k > 0;  k --)
                {
                 From = d;
                 Max = 1 + Edit_Array [k - 1] [d];
                 if  ((j = Edit_Array [k - 1] [d - 1]) > Max)
                     {
                      From = d - 1;
                      Max = j;
                     }
                 if  ((j = 1 + Edit_Array [k - 1] [d + 1]) > Max)
                     {
                      From = d + 1;
                      Max = j;
                     }
                 if  (From == d - 1)
                     {
                      Delta_Stack [(* Delta_Len) ++] = Max - Last - 1;
                      d --;
                      Last = Edit_Array [k - 1] [From];
                     }
                 else if  (From == d + 1)
                     {
                      Delta_Stack [(* Delta_Len) ++] = Last - (Max - 1);
                      d ++;
                      Last = Edit_Array [k - 1] [From];
                     }
                }
              Delta_Stack [(* Delta_Len) ++] = Last + 1;

              k = 0;
              for  (i = (* Delta_Len) - 1;  i > 0;  i --)
                Delta [k ++]
                    = abs (Delta_Stack [i]) * Sign (Delta_Stack [i - 1]);
              (* Delta_Len) --;

#if  0
              //  Check for branch point here caused by uneven
              //  distribution of errors

              Score = Row * BRANCH_PT_MATCH_VALUE - e;
                        // Assumes  BRANCH_PT_MATCH_VALUE
                        //             - BRANCH_PT_ERROR_VALUE == 1.0
              Tail_Len = Row - Max_Score_Len;
              if  (e > MIN_BRANCH_END_DIST / 2
                       && Tail_Len >= MIN_BRANCH_END_DIST
                       && (Max_Score - Score) / Tail_Len >= MIN_BRANCH_TAIL_SLOPE)
                  {
                   (* A_End) = Max_Score_Len;
                   (* T_End) = Max_Score_Len + Max_Score_Best_d;
                   (* Match_To_End) = FALSE;
                   return  Max_Score_Best_e;
                  }
#endif

              (* Match_To_End) = TRUE;
              return  e;
             }
        }

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
      assert (Left <= Right);

      for  (d = Left;  d <= Right;  d ++)
        if  (Edit_Array [e] [d] > Longest)
            {
             Best_d = d;
             Best_e = e;
             Longest = Edit_Array [e] [d];
            }
#if  1
      Score = Longest * BRANCH_PT_MATCH_VALUE - e;
               // Assumes  BRANCH_PT_MATCH_VALUE - BRANCH_PT_ERROR_VALUE == 1.0
      if  (Score > Max_Score)
          {
           Max_Score = Score;
           Max_Score_Len = Longest;
           Max_Score_Best_d = Best_d;
           Max_Score_Best_e = Best_e;
          }
#endif
     }

   (* A_End) = Max_Score_Len;
   (* T_End) = Max_Score_Len + Max_Score_Best_d;
   (* Match_To_End) = FALSE;

   return  e;
  }



static void  Process_Olap
    (Olap_Info_t * olap, char * b_seq, Adjust_t forw_adj [], int adj_ct,
     int frag_len)

//  Find the alignment referred to in  olap , where the  a_iid
//  fragment is in  Frag  and the  b_iid  sequence is in  b_seq .
//  forw_adj [0 .. (adj_ct - 1)]  has
//  adjustment values caused by corrections in the B sequence in
//  the forward orientation.   frag_len  is the length of the B sequence.

  {
   char  * a_part, * b_part, * a_seq;
   double  denom, quality;
   int  a_part_len, b_part_len, a_end, b_end, olap_len;
   int  match_to_end, delta [AS_FRAG_MAX_LEN], delta_len, errors;  //  only MAX_ERRORS needed
   int  sub;

   if  (Verbose_Level > 0)
       printf ("Process_Olap:  %8d %8d %5d %5d  %c\n",
               olap -> a_iid, olap -> b_iid,
               olap -> a_hang, olap -> b_hang,
               olap -> orient);

   sub = olap -> a_iid - Lo_Frag_IID;
   a_seq = Frag [sub] . sequence;

   if  (a_seq == NULL)
       {
        // Deleted fragment, nothing to do
        return;
       }

   Get_Canonical_Olap_Region
       (olap, sub, a_seq, b_seq, forw_adj, adj_ct, frag_len, & a_part, & b_part);

   // Get the alignment

   a_part_len = strlen (a_part);
   b_part_len = strlen (b_part);
   olap_len = OVL_Min_int (a_part_len, b_part_len);

//**ALD
// may want to increase the Error_Bound value for homopolymer type reads

   errors = Prefix_Edit_Dist
              (a_part, a_part_len, b_part, b_part_len,
               Error_Bound [olap_len], & a_end, & b_end,
               & match_to_end, delta, & delta_len);

#if  0
{
 int  rev_delta [MAX_ERRORS], rev_delta_len;
 int  rev_errors, rev_a_end, rev_b_end, rev_match_to_end;
 int  i;

 rev_errors = Rev_Prefix_Edit_Dist
                  (b_part + b_end - 1, b_end + b_adjustment,
                   a_part + a_end - 1, a_end + a_adjustment,
                   1 + errors, & rev_b_end, & rev_a_end,
                   & rev_match_to_end, rev_delta, & rev_delta_len);

 printf (">>> %2d %2d   %4d %4d   %4d %4d   %2d: ",
         errors, rev_errors, a_end, rev_a_end, b_end, rev_b_end, rev_delta_len);
 for  (i = 0;  i < rev_delta_len;  i ++)
   printf (" %3d", rev_delta [i]);
 putchar ('\n');

 rev_errors = Rev_Prefix_Edit_Dist
                  (a_part + a_end - 1, a_end + a_adjustment,
                   b_part + b_end - 1, b_end + b_adjustment,
                   1 + errors, & rev_a_end, & rev_b_end,
                   & rev_match_to_end, rev_delta, & rev_delta_len);

 printf (">>> %2d %2d   %4d %4d   %4d %4d   %2d: ",
         errors, rev_errors, a_end, rev_a_end, b_end, rev_b_end, rev_delta_len);
 for  (i = 0;  i < rev_delta_len;  i ++)
   printf (" %3d", rev_delta [i]);
 putchar ('\n');
}
#endif

   if  (delta_len > 0 && delta [0] == 1 && 0 < olap -> a_hang)
       {
        int  i, stop;

        stop = OVL_Min_int (delta_len, olap -> a_hang);
        for  (i = 0;  i < stop && delta [i] == 1;  i ++)
          ;
        assert (i == stop || delta [i] != -1);
        delta_len -= i;
        memmove (delta, delta + i, delta_len * sizeof (int));
        a_part += i;
        a_end -= i;
        a_part_len -= i;
        errors -= i;
       }
   else if  (delta_len > 0 && delta [0] == -1 && olap -> a_hang < 0)
       {
        int  i, stop;

        stop = OVL_Min_int (delta_len, - olap -> a_hang);
        for  (i = 0;  i < stop && delta [i] == -1;  i ++)
          ;
        assert (i == stop || delta [i] != 1);
        delta_len -= i;
        memmove (delta, delta + i, delta_len * sizeof (int));
        b_part += i;
        b_end -= i;
        b_part_len -= i;
        errors -= i;
       }

   if  (Verbose_Level > 0)
       printf ("  errors = %d  delta_len = %d\n", errors, delta_len);

   Total_Alignments_Ct ++;
   if  (! match_to_end)
       {
        Failed_Alignments_Ct ++;
        if  (Verbose_Level > 0)
            printf ("    alignment failed\n");
        return;
       }

   if  (Verbose_Level > 0)
       Display_Alignment (a_part, a_part_len, b_part, b_part_len,
                          delta, delta_len);

   denom = OVL_Min_int (a_end, b_end);
   if  (denom <= 0.0)
       {
        fprintf (stderr, "ERROR:  Bad alignment ends  a_end = %d  b_end = %d\n",
                 a_end, b_end);
        fprintf (stderr, "   a_iid = %d  b_iid = %d  errors = %d\n",
                 olap -> a_iid, olap -> b_iid, errors);
        return;
       }

   if  (0)
     printf ("errors = %d  denom = %.1f\n", errors, denom);

   quality = errors / denom;
   olap -> corr_erate = AS_OVS_encodeQuality (quality);

#if  1
   if  (quality <= Quality_Threshold
        || (olap -> a_hang <= 0 && Frag [sub] . keep_left)
        || (olap -> b_hang >= 0 && Frag [sub] . keep_right))
#else
//   if  (Find (olap -> a_iid, UF) == Find (olap -> b_iid, UF))
     if  (Olap_In_Unitig (olap))
#endif
       {
        Output_OVL (olap, quality);
       }

   return;
  }



static char *  Read_Fasta
    (FILE * fp)

//  Read string in  fp  (assuming fasta format with  #  comments)
//  and return a pointer to it.

  {
   int  buff_size = 50000;
   char  * buff = (char *) safe_malloc (buff_size);
   int  len, buff_len = 0;
   char  line [MAX_FASTA_LINE];

   while  (fgets (line, MAX_FASTA_LINE, fp) != NULL)
     {
      int  i;

      switch  (line [0])
        {
         case  '\0' :
         case  '>' :
         case  '#' :
           break;

         default :
           len = strlen (line);
           if  (line [len - 1] == '\n')
               {
                len --;
                line [len] = '\0';
               }
           while  (buff_size <= buff_len + len)
             {
              buff_size *= EXPANSION_FACTOR;
              buff = (char *) safe_realloc (buff, buff_size);
             }
           for  (i = 0;  i < len;  i ++)
             line [i] = tolower (line [i]);
           strcpy (buff + buff_len, line);
           buff_len += len;
           assert (buff [buff_len] == '\0');
        }
     }
   buff = (char *) safe_realloc (buff, buff_len + 1);
   fprintf (stderr, "buff_len = %d\n", buff_len);

   return  buff;
  }



static void  Read_Frags
    (void)

//  Open and read fragments with IIDs from  Lo_Frag_IID  to
//  Hi_Frag_IID  from  gkpStore_Path  and store them in
//  global  Frag .

  {
   fragRecord frag_read;
   unsigned  clear_start, clear_end;
   int  i, j;

   if  (Hi_Frag_IID == INT_MAX)
       {
        gkpStore = openGateKeeperStore (gkpStore_Path, FALSE);
        assert (gkpStore != NULL);
        Hi_Frag_IID = getLastElemFragStore (gkpStore);
        closeGateKeeperStore (gkpStore);
       }

   Num_Frags = 1 + Hi_Frag_IID - Lo_Frag_IID;
   Frag = (Frag_Info_t *) safe_calloc (Num_Frags, sizeof (Frag_Info_t));

   gkpStore = openGateKeeperStore(gkpStore_Path, FALSE);
   loadGateKeeperStorePartial(gkpStore, Lo_Frag_IID, Hi_Frag_IID, FRAG_S_SEQ);

   Frag_Stream = openFragStream (gkpStore, FRAG_S_SEQ);
   resetFragStream (Frag_Stream, Lo_Frag_IID, Hi_Frag_IID);

   for  (i = 0;  nextFragStream (Frag_Stream, &frag_read);
           i ++)
     {
      char  seq_buff [AS_READ_MAX_LEN + 1];
      char *seqptr;
      unsigned  deleted;
      int  result;

      deleted = getFragRecordIsDeleted (&frag_read);
      if  (deleted)
          {
           Frag [i] . sequence = NULL;
           Frag [i] . adjust = NULL;
           Frag [i] . adjust_ct = 0;
           continue;
          }

      clear_start = getFragRecordClearRegionBegin(&frag_read, AS_READ_CLEAR_OBT);
      clear_end   = getFragRecordClearRegionEnd  (&frag_read, AS_READ_CLEAR_OBT);

      seqptr = getFragRecordSequence(&frag_read);

      // Make sure that we have a legal lowercase sequence string

      for  (j = clear_start;  j < clear_end;  j ++)
         seq_buff [j] = Filter (seqptr [j]);

      seq_buff [clear_end] = '\0';

      Frag [i] . sequence = strdup (seq_buff + clear_start);
     }

   closeFragStream (Frag_Stream);
   closeGateKeeperStore (gkpStore);

   return;
  }



static void  Read_Olaps
    (void)

//  Open and read those overlaps with first IIDs from  Lo_Frag_IID  to
//  Hi_Frag_IID  from  Olap_Path  and store them in
//  global  Olaps .  If  Olap_From_Store  is true, then the overlaps
//  are read from a binary overlap store; otherwise, they are from
//  a text file in the format produced by
//  get-olaps and each overlap must appear twice, once in each order.

  {
   FILE  * fp;
   int32  a_iid, b_iid;
   int  a_hang, b_hang;
   char  orient [10];
   double  error_rate;
   long int  olap_size;
   long int  ct = 0;

   if  (Olaps_From_Store)
       Get_Olaps_From_Store (Olap_Path, Lo_Frag_IID, Hi_Frag_IID,
                             & Olap, & Num_Olaps);
     else
       {
        fp = File_Open (Olap_Path, "r");

        olap_size = 1000;
        Olap = (Olap_Info_t *) safe_malloc
                   (olap_size * sizeof (Olap_Info_t));

        while  (fscanf (fp, "%d %d %d %d %s %lf",
                        & a_iid, & b_iid, & a_hang, & b_hang,
                        orient, & error_rate)
                  == 6)
          {
           if  (Lo_Frag_IID <= a_iid && a_iid <= Hi_Frag_IID)
               {
                if  (ct >= olap_size)
                    {
                     olap_size *= EXPANSION_FACTOR;
                     Olap = (Olap_Info_t*) safe_realloc
		       (Olap, olap_size * sizeof (Olap_Info_t));
                    }
                Olap [ct] . a_iid = a_iid;
                Olap [ct] . b_iid = b_iid;
                if  (orient [0] == 'O')
                    {
                     Olap [ct] . a_hang = - b_hang;
                     Olap [ct] . b_hang = - a_hang;
                     Olap [ct] . orient = 'I';
                    }
                  else
                    {
                     Olap [ct] . a_hang = a_hang;
                     Olap [ct] . b_hang = b_hang;
                     Olap [ct] . orient = orient [0];
                    }
                Olap [ct] . corr_erate = AS_OVS_encodeQuality(error_rate);
                ct ++;
               }

           if  (a_iid > Hi_Frag_IID)   // Speed up if file is sorted
               break;
          }

        Num_Olaps = ct;
        Olap = (Olap_Info_t*) safe_realloc
	  (Olap, Num_Olaps * sizeof (Olap_Info_t));

        fclose (fp);
       }

   return;
  }



static void  Redo_Olaps
    (void)

//  Read old fragments in  gkpStore  and choose the ones that
//  have overlaps with fragments in  Frag .  Recompute the
//  overlaps, using fragment corrections and output the revised error.

  {
   FILE  * fp;
   fragRecord frag_read;
   unsigned  clear_start, clear_end;
   int  lo_frag, hi_frag;
   uint64  next_olap;
   Correction_Output_t  msg;
   Correction_t  correct [MAX_FRAG_LEN];
   Adjust_t  adjust [MAX_FRAG_LEN];
   int16  adjust_ct;
   int  num_corrects;
   uint32  correct_iid = 0, next_iid;
   int  i, j;

   gkpStore = openGateKeeperStore (gkpStore_Path, FALSE);
   Frag_Stream = openFragStream (gkpStore, FRAG_S_SEQ);

   lo_frag = Olap [0] . b_iid;
   hi_frag = Olap [Num_Olaps - 1] . b_iid;

   resetFragStream (Frag_Stream, lo_frag, hi_frag);

   fp = File_Open (Correct_File_Path, "rb");

   next_olap = 0;
   for  (i = 0;  nextFragStream (Frag_Stream, &frag_read)
                   && next_olap < Num_Olaps;
           i ++)
     {
      char  seq_buff [AS_READ_MAX_LEN + 1];
      char *seqptr;
      char  * seq_ptr = seq_buff;
      Adjust_t  * adjust_ptr = adjust;
      uint32  frag_iid;
      unsigned  deleted;
      int  frag_len, result;

      frag_iid = getFragRecordIID (&frag_read);
      if  (frag_iid < Olap [next_olap] . b_iid)
          continue;

      deleted = getFragRecordIsDeleted (&frag_read);
      if  (deleted)
          continue;

      clear_start = getFragRecordClearRegionBegin(&frag_read, AS_READ_CLEAR_OBT);
      clear_end   = getFragRecordClearRegionEnd  (&frag_read, AS_READ_CLEAR_OBT);

      seqptr = getFragRecordSequence(&frag_read);

      // Make sure that we have a legal lowercase sequence string

      frag_len = 0;
      for  (j = clear_start;  j < clear_end;  j ++)
         seq_buff [frag_len ++] = Filter (seqptr [j]);

      seq_buff [frag_len] = '\0';

      num_corrects = 0;
      next_iid = correct_iid;
      while  (next_iid <= frag_iid)
        {
         if  (fread (& msg, sizeof (Correction_Output_t), 1, fp) != 1)
             {
              next_iid = INT_MAX;
              break;
             }
         if  (msg . frag . is_ID)
             {
              next_iid = msg . frag . iid;
              if  (next_iid <= frag_iid)
                  correct_iid = next_iid;
             }
         else if  (correct_iid == frag_iid)
             correct [num_corrects ++] = msg . corr;
        }
      if  (correct_iid == frag_iid && num_corrects > 0)
          {
//**ALD
//printf ("frag_iid = %d  before corrections\n", frag_iid);
//Fasta_Print (stdout, seq_ptr, "Before");

           Apply_Seq_Corrects (& seq_ptr, & adjust_ptr, & adjust_ct,
                               correct, num_corrects, TRUE);

//**ALD
//printf ("frag_iid = %d  after corrections\n", frag_iid);
//Fasta_Print (stdout, seq_ptr, "After");
          }
        else
          adjust_ct = 0;
      correct_iid = next_iid;

      while  (next_olap < Num_Olaps
                && Olap [next_olap] . b_iid == frag_iid)
        {
         Process_Olap (Olap + next_olap, seq_buff, adjust, adjust_ct,
                       frag_len);
         next_olap ++;
        }
     }

   closeFragStream (Frag_Stream);
   closeGateKeeperStore (gkpStore);

   return;
  }



static int  Rev_Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Match_To_End,
     int * Delta, int * Delta_Len)

//  Return the minimum number of changes (inserts, deletes, replacements)
//  needed to match string  A [0 .. -(m-1)]  with a prefix of string
//   T [0 .. -(n-1)]  if it's not more than  Error_Limit .
//  Note that the match is done in the reverse direction.
//  Put delta description of alignment in  Delta  and set
//  (* Delta_Len)  to the number of entries there if it's a complete
//  match.
//  Set  A_End  and  T_End  to the rightmost positions where the
//  alignment ended in  A  and  T , respectively.
//  Set  Match_To_End  true if the match extended to the end
//  of at least one string; otherwise, set it false to indicate
//  a branch point.

  {
   int  Delta_Stack [AS_FRAG_MAX_LEN];  //  only MAX_ERRORS needed
   double  Score, Max_Score;
   int  Max_Score_Len, Max_Score_Best_d, Max_Score_Best_e;
   int  Best_d, Best_e, From, Last, Longest, Max, Row;
   int  Left, Right;
   int  d, e, i, j, k, shorter;

//   assert (m <= n);
   Best_d = Best_e = Longest = 0;
   (* Delta_Len) = 0;

   shorter = OVL_Min_int (m, n);
   for  (Row = 0;  Row < shorter && A [- Row] == T [- Row];  Row ++)
     ;

   Edit_Array [0] [0] = Row;

   if  (Row == shorter)                              // Exact match
       {
        (* A_End) = (* T_End) = - Row;
        (* Match_To_End) = TRUE;
        return  0;
       }

   Left = Right = 0;
   Max_Score = 0.0;
   Max_Score_Len = Max_Score_Best_d = Max_Score_Best_e = 0;
   for  (e = 1;  e <= Error_Limit;  e ++)
     {
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

         if  (Row == m || Row + d == n)
             {
              // Force last error to be mismatch rather than insertion
              if  (Row == m
                     && 1 + Edit_Array [e - 1] [d + 1]
                          == Edit_Array [e] [d]
                     && d < Right)
                  {
                   d ++;
                   Edit_Array [e] [d] = Edit_Array [e] [d - 1];
                  }

              (* A_End) = - Row;           // One past last align position
              (* T_End) = - Row - d;
              // Compute Delta
              Last = Row;
              (* Delta_Len) = 0;
              for  (k = e;  k > 0;  k --)
                {
                 From = d;
                 Max = 1 + Edit_Array [k - 1] [d];
                 if  ((j = Edit_Array [k - 1] [d - 1]) > Max)
                     {
                      From = d - 1;
                      Max = j;
                     }
                 if  ((j = 1 + Edit_Array [k - 1] [d + 1]) > Max)
                     {
                      From = d + 1;
                      Max = j;
                     }
                 if  (From == d - 1)
                     {
                      Delta_Stack [(* Delta_Len) ++] = Max - Last - 1;
                      d --;
                      Last = Edit_Array [k - 1] [From];
                     }
                 else if  (From == d + 1)
                     {
                      Delta_Stack [(* Delta_Len) ++] = Last - (Max - 1);
                      d ++;
                      Last = Edit_Array [k - 1] [From];
                     }
                }
              Delta_Stack [(* Delta_Len) ++] = Last + 1;

              k = 0;
              for  (i = (* Delta_Len) - 1;  i > 0;  i --)
                Delta [k ++]
                    = abs (Delta_Stack [i]) * Sign (Delta_Stack [i - 1]);
              (* Delta_Len) --;

              //  Check for branch point here caused by uneven
              //  distribution of errors

#if  0
              Score = Row * BRANCH_PT_MATCH_VALUE - e;
                        // Assumes  BRANCH_PT_MATCH_VALUE
                        //             - BRANCH_PT_ERROR_VALUE == 1.0
              Tail_Len = Row - Max_Score_Len;
              if  (e > MIN_BRANCH_END_DIST / 2
                       && Tail_Len >= MIN_BRANCH_END_DIST
                       && (Max_Score - Score) / Tail_Len >= MIN_BRANCH_TAIL_SLOPE)
                  {
                   (* A_End) = Max_Score_Len;
                   (* T_End) = Max_Score_Len + Max_Score_Best_d;
                   (* Match_To_End) = FALSE;
                   return  Max_Score_Best_e;
                  }
#endif

              (* Match_To_End) = TRUE;
              return  e;
             }
        }

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
      assert (Left <= Right);

      for  (d = Left;  d <= Right;  d ++)
        if  (Edit_Array [e] [d] > Longest)
            {
             Best_d = d;
             Best_e = e;
             Longest = Edit_Array [e] [d];
            }
#if  1
      Score = Longest * BRANCH_PT_MATCH_VALUE - e;
               // Assumes  BRANCH_PT_MATCH_VALUE - BRANCH_PT_ERROR_VALUE == 1.0
      if  (Score > Max_Score)
          {
           Max_Score = Score;
           Max_Score_Len = Longest;
           Max_Score_Best_d = Best_d;
           Max_Score_Best_e = Best_e;
          }
#endif
     }

   (* A_End) = - Max_Score_Len;
   (* T_End) = - Max_Score_Len - Max_Score_Best_d;
   (* Match_To_End) = FALSE;

   return  e;
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



static int  Union
    (int i, int j, int a [])

//  Union sets  i  and  j  in Union-Find array  a [] .
//  Return the ID of the resulting set.

  {
   assert (a [i] < 0 && a [j] < 0);

   if  (i == j)
       return  i;

   if  (a [i] <= a [j])
       {
        a [i] += a [j];
        a [j] = i;
        return  i;
       }
     else
       {
        a [j] += a [i];
        a [i] = j;
        return  j;
       }
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
       "USAGE:  %s [-d <dna-file>] [-o <ovl_file>] [-q <quality>]\n"
       "            [-x <del_file>] [-F OlapFile] [-S OlapStore]\n"
       "            [-c <cgb_file>] [-e <erate_file>\n"
       "           <gkpStore> <CorrectFile> <lo> <hi>\n"
       "\n"
       "Recalculates overlaps for frags  <lo> .. <hi>  in\n"
       " <gkpStore>  using corrections in  <CorrectFile> \n"
       "\n"
       "Options:\n"
       "-e <erate-file>  specifies binary file to dump corrected erates to\n"
       "                 for later updating of olap store by  update-erates \n"
       "-F             specify file of sorted overlaps to use (in the format\n"
       "               produced by  get-olaps\n"
       "-o <ovl_file>  specifies name of file to which OVL messages go\n"
       "-q <quality>   overlaps less than this error rate are\n"
       "               automatically output\n"
       "-S             specify the binary overlap store containing overlaps to use\n"
       "-v <num>       specify level of verbose outputs, higher is more\n"
       "-X <del_file>  specifies name of file where list of ovl's to delete goes\n",
       command);

   return;
  }


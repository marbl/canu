
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
* Module:  AutoScreenOVL.c
* Description:
*   Based on overlaps between DNA fragment sequences, create a set
*   of sequences that cover the high-overlap regions.
* 
*    Programmer:  A. Delcher
*       Started:   15 Jan 2001
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: AutoScreenOVL.c,v 1.9 2007-03-06 01:02:44 brianwalenz Exp $
 * $Revision: 1.9 $
*/

static char CM_ID[] = "$Id: AutoScreenOVL.c,v 1.9 2007-03-06 01:02:44 brianwalenz Exp $";


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
#include  "OlapStoreOVL.h"
#include  "AS_PER_gkpStore.h"
#include  "AS_PER_genericStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_UTL_version.h"
#include  "AS_CGW_dataTypes.h"
#include  "Globals_CGW.h"
#include  "ScaffoldGraph_CGW.h"



//  Constants

#define  AGREEMENT_LEVEL          0.80
    //  Fraction of votes that must agree to count as a match
#define  BRANCH_PT_MATCH_VALUE    0.272
    //  Value to add for a match in finding branch points
    //  1.20 was the calculated value for 6% vs 35% error discrimination
    //  Converting to integers didn't make it faster
#define  BRANCH_PT_ERROR_VALUE    -0.728
    //  Value to add for a mismatch in finding branch points
    //   -2.19 was the calculated value for 6% vs 35% error discrimination
    //  Converting to integers didn't make it faster
#define  DEFAULT_DEGREE_THRESHOLD    100
    //  Default value for  Degree_Threshold
#define  DEFAULT_END_EXCLUDE_LEN     3
    //  Default value for  End_Exclude_Len
#define  DEFAULT_HALF_LEN            4
    //  Default value for bases on each side of SNP to vote for change
#define  DEFAULT_WINDOW_LEN            9
    //  Default value for  Window_Len
#define  DEFAULT_VOTE_QUALIFY_LEN    9
    //  Default value for bases surrounding SNP to vote for change
#define  EDIT_DIST_PROB_BOUND        1e-4
    //  Probability limit to "band" edit-distance calculation
    //  Determines  NORMAL_DISTRIB_THOLD
#define  ERRORS_FOR_FREE             1
    //  The number of errors that are ignored in setting probability
    //  bound for terminating alignment extensions in edit distance
    //  calculations
#define  EXPANSION_FACTOR            1.4
    //  Factor by which to grow memory in olap array when reading it
#define  INITIAL_BUFFER_SIZE         10000
    //  Initial number of bytes to use for sequence strings
#define  MATCH_SCORE                 1.0
    //  Score to add for match in finding length of alignment
#define  MAX_EDGE_ERATE              0.04
    //  Largest error overlap used to build graph edges
#define  MAX_ERROR_RATE              AS_GUIDE_ERROR_RATE
    //  The largest error allowed in overlaps
#define  MAX_FRAG_LEN                2048
    //  The longest fragment allowed
#define  MAX_ERRORS                  (1 + (int) (MAX_ERROR_RATE * MAX_FRAG_LEN))
    //  Most errors in any edit distance computation
#define  MAX_HDR_LEN                 256
    //  Longest allowable fasta header line
#define  MAX_MEMORY_STORE            50000
    //  The most fragments allowed in memory store
#define  MAX_PROGRESS_DIFF           3
    //  Maximum allowable difference between expected progress (i.e., hang)
    //  from an overlap, and the progress found by DP_Compare
#define  MIN_BRANCH_END_DIST         20
    //  Branch points must be at least this many bases from the
    //  end of the fragment to be reported
#define  MIN_BRANCH_TAIL_SLOPE       0.20
    //  Branch point tails must fall off from the max by at least
    //  this rate
#define  MIN_MATCH_LEN               40
    //  Number of bases in match region in order to count it
#define  TARGET_PROGRESS                20
    //  Target for choosing overlaps to extend the current sequence
    //  Want overlaps that are close to this.
#define  MISMATCH_SCORE              -3.0
    //  Score to add for non-match in finding length of alignment
#define  NORMAL_DISTRIB_THOLD        3.62
    //  Determined by  EDIT_DIST_PROB_BOUND



//  Type definitions

typedef  struct
  {
   uint16  confirmed;
   uint16  deletes;
   uint16  a_subst;
   uint16  c_subst;
   uint16  g_subst;
   uint16  t_subst;
   uint16  a_insert;
   uint16  c_insert;
   uint16  g_insert;
   uint16  t_insert;
  }  Vote_Tally_t;

typedef  struct
  {
   unsigned  iid : 31;
   unsigned  over_left : 1;   // If true, overlap involves left (A or 5') end
                              //   of frag  iid
   int16  progress;           // Bases overlap extends sequence.  Uses only
                              //   high-degree portion.
   int16  olap_len;           // Number of bases in overlap region
  }  Graph_Edge_t;

typedef  struct
  {
   char  * sequence;
   Vote_Tally_t  * vote;
   Graph_Edge_t  * left_list, * right_list;
   uint16  left_match_len, right_match_len;
   uint16  left_list_len, right_list_len;
   uint16  left_degree, right_degree;
   uint16  len;
   unsigned  contained : 1;
   unsigned  need_sequence : 1;
   unsigned  coupled_ends : 1;
   unsigned  left_covered : 1;
   unsigned  right_covered : 1;
   unsigned  high_left_degree : 1;
   unsigned  high_right_degree : 1;
   unsigned  used : 1;
  }  Frag_Info_t;

typedef  struct
  {
   int32  iid;
   unsigned  on_left : 1;
   unsigned  degree : 31;
  }  End_Ref_t;

typedef  struct                 
  {
   int32  a_iid, b_iid;
   int16  a_hang, b_hang;
   ERATE_TYPE  e_rate : ERATE_BITS;   // original error rate (aka quality)
   char  orient;
  }  Olap_Info_t;



//  Static Globals

static char  * Degree_Path;
    // Name of file containing fragment overlap degrees
static int  Degree_Threshold = DEFAULT_DEGREE_THRESHOLD;
    // Minimum overlap degree on fragment end to process
static int  * Edit_Array [MAX_ERRORS];
    // Use for alignment calculation.  Points into  Edit_Space .
static int  Edit_Match_Limit [MAX_ERRORS] = {0};
    // This array [e] is the minimum value of  Edit_Array [e] [d]
    // to be worth pursuing in edit-distance computations between guides
static int  Edit_Space [(MAX_ERRORS + 4) * MAX_ERRORS];
    // Memory used by alignment calculation
static int  End_Exclude_Len = DEFAULT_END_EXCLUDE_LEN;
    // Length of ends of exact-match regions not used in preventing
    // sequence correction
static int  Error_Bound [MAX_FRAG_LEN + 1];
    //  This array [i]  is the maximum number of errors allowed
    //  in a match between sequences of length  i , which is
    //  i * MAXERROR_RATE .
static int  Failed_Olaps = 0;
    // Counts overlaps that didn't make the error bound
static Frag_Info_t  * Frag;
    // Sequence and vote information for current range of fragments
    // being corrected
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
static int32  Lo_Frag_IID;
    // Internal ID of first fragment in frag store to process
static int  Num_Frags;
    // Number of fragments being corrected
static int  Num_Graph_Nodes;
    // Number of non-contained fragments with at least one high-degree end
static char  * Olap_File_Path;
    // Name of file containing a sorted list of overlaps
static int  Olaps_From_Store = FALSE;
    // If set true by  -S  option, then overlaps will be read from
    // the specified overlap store.
static FILE  * Screen_File;
    // File for outputing screen sequences
static int  Verbose = 0;
    // Determines level of debugging messages
static int  Vote_Qualify_Len = DEFAULT_VOTE_QUALIFY_LEN;
    // Number of bases surrounding a SNP to vote for change
static int  Window_Len = DEFAULT_WINDOW_LEN;
    // Length of window to scan consensus regions



//  Static Functions

static void  Add_Olap_To_List
    (Graph_Edge_t * * list, uint16 * list_len, int iid,
     int over_left, int progress, int olap_len);
static void  Analyze_Alignment
    (int delta [], int delta_len, char * a_part, char * b_part,
     int a_len, int a_offset, int sub);
static int  Binomial_Bound
    (int, double, int, double);
static int  By_B_IID
    (const void * a, const void * b);
static int  By_Degree
    (const void * a, const void * b);
static char  Complement
    (char);
static void  Display_Alignment
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct);
static void  Display_Frags
    (void);
static void  Dump_Frag
    (int sub, int lo, int hi, char * tag);
static void  Fasta_Print
    (FILE * fp, char * s, char * hdr);
static char  Filter
    (char ch);
static Overlap *  Find_Overlap
    (char * seq1, char * seq2, ChunkOrientationType orientation, 
     int min_ahang, int max_ahang,
     double erate, double thresh, int minlen, CompareOptions what);
static void  Get_Consensus
    (int iid, int use_left, int complement, char * * buff, int * buff_size);
static int  Get_Olap
    (Olap_Info_t * olap, FILE * fp, OVL_Stream_t * stream);
static void  Get_Start
    (int iid, int over_left, int * start_iid, int * start_left);
static void  Initialize_Globals
    (void);
static int  OVL_Max_int
    (int a, int b);
static int  OVL_Min_int
    (int a, int b);
static void  Make_Graph_Edges
    (void);
static void  Make_Graph_Nodes
    (void);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Process_Graph
    (void);
static void  Process_Olap
    (Olap_Info_t * olap);
static void  Read_Degrees
    (void);
static void  Read_Frags
    (void);
static void  Read_Olaps
    (void);
static void  Rev_Complement
    (char * s);
static int  Sign
    (int a);
static void  Traverse
    (int iid, int use_left);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   Parse_Command_Line  (argc, argv);

   Initialize_Globals ();

   Read_Degrees ();

   Read_Frags ();

   Read_Olaps ();

   fprintf (stderr, "Failed overlaps = %d\n", Failed_Olaps);

   Make_Graph_Nodes ();

   Make_Graph_Edges ();

   Screen_File = Local_File_Open ("screen.fasta", "w");

   Process_Graph ();

   fclose (Screen_File);

   return  0;
  }



static void  Add_Olap_To_List
    (Graph_Edge_t * * list, uint16 * list_len, int iid,
     int over_left, int progress, int olap_len)

//  Add overlap to frag  iid  with orientation indicated by  over_left
//  to list  (* list) and increment  (* list_len) , the number of
//  entries on that list.   progress  indicates the number of bases
//  the overlap extends beyond the end of the first fragment.
//  olap_len  is the number of bases in the overlap region

  {
   int  new_len = (* list_len) + 1;

   * list = (Graph_Edge_t *) safe_realloc
                (* list, new_len * sizeof (Graph_Edge_t));
   (* list) [* list_len] . iid = iid;
   (* list) [* list_len] . over_left = over_left;
   (* list) [* list_len] . progress = progress;
   (* list) [* list_len] . olap_len = olap_len;
   (* list_len) = new_len;

   return;
  }



static void  Analyze_Alignment
    (int delta [], int delta_len, char * a_part, char * b_part,
     int  a_len, int a_offset, int sub)

//  Analyze the delta-encoded alignment in  delta [0 .. (delta_len - 1)]
//  between  a_part  and  b_part  and store the resulting votes
//  about the a sequence in  Frag [sub] .  The alignment starts
//  a_offset  bytes in from the start of the a sequence in  Frag [sub] .

  {
   int  i, j, k, m;

   i = j = 0;

   for  (k = 0;  k < delta_len;  k ++)
     {
      assert (delta [k] != 0);

      for  (m = 1;  m < abs (delta [k]);  m ++)
        {
         switch  (b_part [j])
           {
            case  'a' :
              Frag [sub] . vote [a_offset + i] . a_subst ++;
              break;
            case  'c' :
              Frag [sub] . vote [a_offset + i] . c_subst ++;
              break;
            case  'g' :
              Frag [sub] . vote [a_offset + i] . g_subst ++;
              break;
            case  't' :
              Frag [sub] . vote [a_offset + i] . t_subst ++;
              break;
            default :
              fprintf (stderr, "ERROR:  Bad sequence char \'%c\'\n",
                       b_part [j]);
           }
         i ++;
         j ++;
        }
      if  (delta [k] < 0)
          {
           // Only one insert allowed at a position
           if  (delta [k] != -1 || k == 0 || delta [k - 1] > 0)
               {
                switch  (b_part [j])
                  {
                   case  'a' :
                     Frag [sub] . vote [a_offset + i - 1] . a_insert ++;
                     break;
                   case  'c' :
                     Frag [sub] . vote [a_offset + i - 1] . c_insert ++;
                     break;
                   case  'g' :
                     Frag [sub] . vote [a_offset + i - 1] . g_insert ++;
                     break;
                   case  't' :
                     Frag [sub] . vote [a_offset + i - 1] . t_insert ++;
                     break;
                   default :
                     fprintf (stderr, "ERROR:  Bad sequence char \'%c\'\n",
                              b_part [j]);
                  }
               }
           j ++;
          }
        else
          {
           Frag [sub] . vote [a_offset + i] . deletes ++;
           i ++;
          }
     }

   while  (i < a_len)
     {
      switch  (b_part [j])
        {
         case  'a' :
           Frag [sub] . vote [a_offset + i] . a_subst ++;
           break;
         case  'c' :
           Frag [sub] . vote [a_offset + i] . c_subst ++;
           break;
         case  'g' :
           Frag [sub] . vote [a_offset + i] . g_subst ++;
           break;
         case  't' :
           Frag [sub] . vote [a_offset + i] . t_subst ++;
           break;
         default :
           fprintf (stderr, "ERROR:  Bad sequence char \'%c\'\n",
                    b_part [j]);
        }
      i ++;
      j ++;
     }

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



static int  By_Degree
    (const void * a, const void * b)

//  Compare the values in  a  and  b  as  (* End_Ref_t) 's,
//  by  degree .
//  Return  -1  if  a < b ,  0  if  a == b , and  1  if  a > b .
//  Used for  qsort .

  {
   End_Ref_t  * x, * y;

   x = (End_Ref_t *) a;
   y = (End_Ref_t *) b;

   if  (x -> degree < y -> degree)
       return  -1;
   else if  (x -> degree > y -> degree)
       return  1;

   return  0;
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
        fprintf (stderr, "ERROR(complement):  Unexpected character `%c\'\n", ch);
        exit (-1);
     }

   return  'x';    // Just to make the compiler happy
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

   return;
  }



static void  Dump_Frag
    (int sub, int lo, int hi, char * tag)

//  Print a fasta version of consensus votes for  Frag [sub]
//  at positions  lo .. hi .  Print  tag  on the header line.

  {
   int  i, ct;

   printf (">%d  %d..%d  %s  len = %d\n", sub, lo, hi, tag, Frag [sub] . len);

   if  (Verbose > 1)
       {
        printf ("pos: ch %4s %4s %4s %4s %4s %4s %4s %4s %4s\n",
                "del", "asub", "csub", "gsub", "tsub",
                "ains", "cins", "gins", "tins");
//      for  (i = lo;  i <= hi;  i ++)
        for  (i = 0;  i < Frag [sub] . len;  i ++)
          printf ("%3d:  %c %4d %4d %4d %4d %4d %4d %4d %4d %4d\n",
                  i, Frag [sub] . sequence [i],
                  Frag [sub] . vote [i] . deletes,
                  Frag [sub] . vote [i] . a_subst,
                  Frag [sub] . vote [i] . c_subst,
                  Frag [sub] . vote [i] . g_subst,
                  Frag [sub] . vote [i] . t_subst,
                  Frag [sub] . vote [i] . a_insert,
                  Frag [sub] . vote [i] . c_insert,
                  Frag [sub] . vote [i] . g_insert,
                  Frag [sub] . vote [i] . t_insert);
       }

   ct = 0;
   for  (i = lo;  i <= hi;  i ++)
     {
      int  max, subst_total, insert_total;
      char  ch;

      ch = ' ';
      max = Frag [sub] . vote [i] . deletes;
      if  (Frag [sub] . vote [i] . a_subst > max)
          {
           max = Frag [sub] . vote [i] . a_subst;
           ch = 'a';
          }
      if  (Frag [sub] . vote [i] . c_subst > max)
          {
           max = Frag [sub] . vote [i] . c_subst;
           ch = 'c';
          }
      if  (Frag [sub] . vote [i] . g_subst > max)
          {
           max = Frag [sub] . vote [i] . g_subst;
           ch = 'g';
          }
      if  (Frag [sub] . vote [i] . t_subst > max)
          {
           max = Frag [sub] . vote [i] . t_subst;
           ch = 't';
          }
      if  (ch != ' ')
          {
           if  (ct == 60)
               {
                putchar ('\n');
                ct = 0;
               }
           ct ++;
           putchar (ch);
          }

      if  (i == hi)
          break;

      subst_total = Frag [sub] . vote [i] . deletes
                      + Frag [sub] . vote [i] . a_subst
                      + Frag [sub] . vote [i] . c_subst
                      + Frag [sub] . vote [i] . g_subst
                      + Frag [sub] . vote [i] . t_subst;
      insert_total = Frag [sub] . vote [i] . a_insert
                       + Frag [sub] . vote [i] . c_insert
                       + Frag [sub] . vote [i] . g_insert
                       + Frag [sub] . vote [i] . t_insert;
      if  (subst_total > 0
             && (double) insert_total / subst_total >= AGREEMENT_LEVEL)
          {
           ch = 'a';
           max = Frag [sub] . vote [i] . a_insert;
           if  (Frag [sub] . vote [i] . c_insert > max)
               {
                max = Frag [sub] . vote [i] . c_insert;
                ch = 'c';
               }
           if  (Frag [sub] . vote [i] . g_insert > max)
               {
                max = Frag [sub] . vote [i] . g_insert;
                ch = 'g';
               }
           if  (Frag [sub] . vote [i] . t_insert > max)
               {
                max = Frag [sub] . vote [i] . t_insert;
                ch = 't';
               }
           if  (ct == 60)
               {
                putchar ('\n');
                ct = 0;
               }
           ct ++;
           putchar (ch);
          }
     }
   putchar ('\n');

   return;
  }



static void  Fasta_Print
    (FILE * fp, char * s, char * hdr)

//  Print string  s  in fasta format to  fp .  Put string  hdr
//  on header line.

  {
   int  ct = 0;

   if  (hdr != NULL)
       fprintf (fp, "> %s\n", hdr);

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

  if  (Verbose > 0)
      {
       Print_Overlap (stdout, seq1, seq2, omesg);
      }

  // return seq1 to its original state
  if (orientation == BA_AB || orientation == BA_BA)
	Complement_Seq( seq1 );
  
  // omesg->begpos is the a-hang, omesg->endpos is the b-hang
  return omesg;
}



static void  Get_Consensus
    (int iid, int use_left, int complement, char * * buff, int * buff_size)

//  Extract consensus sequence from votes in  Frag [iid]  for end
//  indicated by  use_left .   Put it in  (* buff) .  If necessary
//  Expand  (* buff)  and put the new size in  (* buff_size) .
//  If  complement , then reverse-complement the sequence.

  {
   Frag_Info_t  * frag = Frag + iid;
   int  i, ct, lo, hi;

   if  (use_left)
       {
        lo = 0;
        hi = frag -> left_match_len;
       }
     else
       {
        lo = frag -> len - frag -> right_match_len;
        hi = frag -> len - 1;
       }

   if  (2 * (hi - lo) > (* buff_size))
       {
        int  new_size = (* buff_size) * EXPANSION_FACTOR;

        if  (new_size < 2 * (hi - lo))
            new_size = 2 * (hi - lo);

        (* buff) = (char *) safe_realloc ((* buff), new_size);
        (* buff_size) = new_size;
       }

   ct = 0;
   for  (i = lo;  i <= hi;  i ++)
     {
      int  max, subst_total, insert_total;
      char  ch;

      ch = ' ';
      max = frag -> vote [i] . deletes;
      if  (frag -> vote [i] . a_subst > max)
          {
           max = frag -> vote [i] . a_subst;
           ch = 'a';
          }
      if  (frag -> vote [i] . c_subst > max)
          {
           max = frag -> vote [i] . c_subst;
           ch = 'c';
          }
      if  (frag -> vote [i] . g_subst > max)
          {
           max = frag -> vote [i] . g_subst;
           ch = 'g';
          }
      if  (frag -> vote [i] . t_subst > max)
          {
           max = frag -> vote [i] . t_subst;
           ch = 't';
          }
      if  (ch != ' ')
          (* buff) [ct ++] = ch;

      if  (i == hi)
          break;

      subst_total = frag -> vote [i] . deletes
                      + frag -> vote [i] . a_subst
                      + frag -> vote [i] . c_subst
                      + frag -> vote [i] . g_subst
                      + frag -> vote [i] . t_subst;
      insert_total = frag -> vote [i] . a_insert
                       + frag -> vote [i] . c_insert
                       + frag -> vote [i] . g_insert
                       + frag -> vote [i] . t_insert;
      if  (subst_total > 0
             && (double) insert_total / subst_total >= AGREEMENT_LEVEL)
          {
           ch = 'a';
           max = frag -> vote [i] . a_insert;
           if  (frag -> vote [i] . c_insert > max)
               {
                max = frag -> vote [i] . c_insert;
                ch = 'c';
               }
           if  (frag -> vote [i] . g_insert > max)
               {
                max = frag -> vote [i] . g_insert;
                ch = 'g';
               }
           if  (frag -> vote [i] . t_insert > max)
               {
                max = frag -> vote [i] . t_insert;
                ch = 't';
               }
           (* buff) [ct ++] = ch;
          }
     }

   (* buff) [ct] = '\0';
 
   if  (complement)
       Rev_Complement ((* buff));

   return;
  }



static int  Get_Olap
    (Olap_Info_t * olap, FILE * fp, OVL_Stream_t * stream)

//  Read the next overlap into  olap .  It comes from file  stream
//  if  Olaps_From_Store  is true; otherwise, it comes from  fp .
//  Return  TRUE  is an overlap is successfully read; otherwise,
//  return  FALSE .

  {
   Long_Olap_Data_t  stream_olap;
   char  orient [10];
   double  error_rate;

   if  (Olaps_From_Store)
       {
        if  (! Next_From_OVL_Stream (& stream_olap, stream))
            return  FALSE;
        olap -> a_iid = stream_olap . a_iid;
        olap -> b_iid = stream_olap . b_iid;
        olap -> a_hang = stream_olap . a_hang;
        olap -> b_hang = stream_olap . b_hang;
        olap -> orient = stream_olap . flipped ? 'I' : 'N';
        olap -> e_rate = stream_olap . orig_erate;
       }
     else
       {
        int  a_hang, b_hang;

        if  (fscanf (fp, "%u %u %d %d %s %lf",
                     & (olap -> a_iid), & (olap -> b_iid),
                     & a_hang, & b_hang,
                     orient, & error_rate)
                != 6)
            return  FALSE;
        olap -> a_hang = a_hang;
        olap -> b_hang = b_hang;
        olap -> orient = orient [0];
        olap -> e_rate = Shrink_Quality (error_rate);
       }

   return  TRUE;
  }



static void  Get_Start
    (int iid, int over_left, int * start_iid, int * start_left)

//  Follow the set of overlaps from  iid  over the end indicated by
//  over_left  until it ends.  Then set  (* start_iid)  to the
//  fragment where the path ended and  (* start_left)  to which end
//  of that fragment the overlap passed over.  The path maximizes
//  progress in a greedy fashion.

  {
   Graph_Edge_t  * edge, * best_edge = NULL;
   int  stack [Num_Graph_Nodes];
   int  top = 0;
   int  i, j, using_left, max;

   using_left = over_left;
   for  (i = iid;  i > 0; )
     {
      Frag [i] . used = TRUE;
      stack [top ++] = i;

      max = 0;
      if  (using_left)
          {
           for  (j = 0;  j < Frag [i] . left_list_len;  j ++)
             {
              edge = Frag [i] . left_list + j;

              if  (Frag [edge -> iid] . used)
                  continue;
              if  ((edge -> over_left && Frag [edge -> iid] . left_covered)
                     || (! edge -> over_left
                           && Frag [edge -> iid] . right_covered))
                  continue;
              if  (edge -> progress > max)
                  {
                   max = edge -> progress;
                   best_edge = edge;
                  }
             }
          }
        else
          {
           for  (j = 0;  j < Frag [i] . right_list_len;  j ++)
             {
              edge = Frag [i] . right_list + j;

              if  (Frag [edge -> iid] . used)
                  continue;
              if  ((edge -> over_left && Frag [edge -> iid] . left_covered)
                     || (! edge -> over_left
                           && Frag [edge -> iid] . right_covered))
                  continue;
              if  (edge -> progress > max)
                  {
                   max = edge -> progress;
                   best_edge = edge;
                  }
             }
          }
      if  (max == 0)
          {
           if  (Frag [i] . coupled_ends)
               using_left = ! using_left;            // turn around
           break;
          }

      i = best_edge -> iid;
      if  (Frag [best_edge -> iid] . coupled_ends)
          using_left = ! best_edge -> over_left;     // keep going
        else
          {
           using_left = best_edge -> over_left;     // stop & turn around
           break;
          }
     }

   (* start_iid) = i;
   (* start_left) = using_left;

   while  (top > 0)
     Frag [stack [-- top]] . used = FALSE;

   return;
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



static void  Make_Graph_Edges
    (void)

//  Create overlap edges between fragments with  high_left_degree
//  or high_right_degree flags set.

  {
   FILE  * fp = NULL;
   OVL_Store_t  * ovl_store = NULL;
   OVL_Stream_t  * ovl_stream = NULL;
   Olap_Info_t  olap;
   int  progress, olap_len;
   int  ct = 0;
   ERATE_TYPE  max_erate;

   max_erate = Shrink_Quality (MAX_EDGE_ERATE);

   if  (Olaps_From_Store)
       {
        ovl_store = New_OVL_Store ();
        if  (Open_OVL_Store (ovl_store, Olap_File_Path) != 0)
            {
             fprintf (stderr, "ERROR:  Failed to open overlap store \"%s\"\n",
                      Olap_File_Path);
             exit (EXIT_FAILURE);
            }
        ovl_stream = New_OVL_Stream ();
        Init_OVL_Stream (ovl_stream, 1, Last_Frag_In_OVL_Store (ovl_store),
                         ovl_store);
       }
     else
       fp = Local_File_Open (Olap_File_Path, "r");

   while  (Get_Olap (& olap, fp, ovl_stream))
     {
      if  (olap . e_rate > max_erate)
          continue;

      if  (olap . a_iid < Lo_Frag_IID || Hi_Frag_IID < olap . a_iid
             || olap . b_iid < Lo_Frag_IID || Hi_Frag_IID < olap . b_iid)
          continue;
             
      if  (! Frag [olap . a_iid] . high_left_degree
             && ! Frag [olap . a_iid] . high_right_degree)
          continue;

      if  (olap . orient == 'O')
          {
           int  save = olap . a_hang;
           olap . a_hang = - olap . b_hang;
           olap . b_hang = - save;
           olap . orient = 'I';
          }

      if  ((olap . a_hang <= 0 && olap . b_hang >= 0)
             || (olap . a_hang >= 0 && olap . b_hang <= 0))
          continue;           // contained

      switch (olap . orient)
        {
         case  'I' :
           if  ((olap . a_hang >= 0
                      && ! (Frag [olap . b_iid] . high_right_degree
                              && Frag [olap . a_iid] . high_right_degree))
                 || (olap . a_hang <= 0
                      && ! (Frag [olap . a_iid] . high_left_degree
                              && Frag [olap . b_iid] . high_left_degree)))
               continue;
           if  (olap . a_hang >= 0)
               {
                progress = olap . b_hang - Frag [olap . b_iid] . len
                             + Frag [olap . b_iid] . right_match_len;
                olap_len = Frag [olap . b_iid] . right_match_len;
                if  (progress > 0)
                    olap_len -= progress;
                Add_Olap_To_List (& Frag [olap . a_iid] . right_list,
                                  & Frag [olap . a_iid] . right_list_len,
                                  olap . b_iid, FALSE, progress, olap_len);
               }
             else
               {
                progress = - olap . a_hang - Frag [olap . b_iid] . len
                             + Frag [olap . b_iid] . left_match_len;
                olap_len = Frag [olap . b_iid] . left_match_len;
                if  (progress > 0)
                    olap_len -= progress;
                if  (olap . b_hang > 0)
                    olap_len -= olap . b_hang;
                Add_Olap_To_List (& Frag [olap . a_iid] . left_list,
                                  & Frag [olap . a_iid] . left_list_len,
                                  olap . b_iid, TRUE, progress, olap_len);
               }
           break;
         case  'N' :
           if  ((olap . a_hang >= 0
                      && ! (Frag [olap . b_iid] . high_left_degree
                              && Frag [olap . a_iid] . high_right_degree))
                 || (olap . a_hang <= 0
                      && ! (Frag [olap . a_iid] . high_left_degree
                              && Frag [olap . b_iid] . high_right_degree)))
               continue;
           if  (olap . a_hang >= 0)
               {
                progress = olap . b_hang - Frag [olap . b_iid] . len
                             + Frag [olap . b_iid] . left_match_len;
                olap_len = Frag [olap . b_iid] . left_match_len;
                if  (progress > 0)
                    olap_len -= progress;
                Add_Olap_To_List (& Frag [olap . a_iid] . right_list,
                                  & Frag [olap . a_iid] . right_list_len,
                                  olap . b_iid, TRUE, progress, olap_len);
               }
             else
               {
                progress = - olap . a_hang - Frag [olap . b_iid] . len
                             + Frag [olap . b_iid] . right_match_len;
                olap_len = Frag [olap . b_iid] . right_match_len;
                if  (progress > 0)
                    olap_len -= progress;
                if  (olap . b_hang > 0)
                    olap_len -= olap . b_hang;
                Add_Olap_To_List (& Frag [olap . a_iid] . left_list,
                                  & Frag [olap . a_iid] . left_list_len,
                                  olap . b_iid, FALSE, progress, olap_len);
               }
           break;
         default :
           fprintf (stderr, "ERROR:  Unexpected orientation %c\n",
                    olap . orient);
           exit (EXIT_FAILURE);
        }

      ct ++;
     }

   if  (Olaps_From_Store)
       {
        Free_OVL_Stream (ovl_stream);
        Free_OVL_Store (ovl_store);
       }
     else
       fclose (fp);

   fprintf (stderr, "Created %d edges\n", ct);

   if  (Verbose > 0)
       {
        int  i, j;

        for  (i = Lo_Frag_IID;  i <= Hi_Frag_IID;  i ++)
          {
           if  (Frag [i] . high_left_degree)
               {
                printf ("> %d left edges:\n", i);
                for  (j = 0;  j < Frag [i] . left_list_len;  j ++)
                  printf (" %6d  %c  %4d\n",
                          Frag [i] . left_list [j] . iid,
                          Frag [i] . left_list [j] . over_left ? 'L' : 'R',
                          Frag [i] . left_list [j] . progress);
               }
           if  (Frag [i] . high_right_degree)
               {
                printf ("> %d right edges:\n", i);
                for  (j = 0;  j < Frag [i] . right_list_len;  j ++)
                  printf (" %6d  %c  %4d\n",
                          Frag [i] . right_list [j] . iid,
                          Frag [i] . right_list [j] . over_left ? 'L' : 'R',
                          Frag [i] . right_list [j] . progress);
               }
          }
       }

   return;
  }



static void  Make_Graph_Nodes
    (void)

//  Create the non-contained, high-degree conserved sequences to
//  process as nodes in overlap graph.

  {
   double  score, max_score;
   int  left_max, right_min;
   int  max_sub, consec_bad;
   int  len;
   int  i, j;

   for  (i = Lo_Frag_IID;  i <= Hi_Frag_IID;  i ++)
     {
      if  (Frag [i] . contained)
          {
           Frag [i] . high_left_degree
               = Frag [i] . high_right_degree = FALSE;
           continue;
          }

      left_max = -1;
      right_min = AS_READ_MAX_LEN;

      len = Frag [i] . len;
      if  (Frag [i] . high_left_degree)
          {
           score = max_score = 0.0;
           max_sub = -1;
           consec_bad = 0;
           for  (j = 0;  Frag [i] . sequence [j] != '\0'
                           && consec_bad < Degree_Threshold;  j ++)
             {
              int  total_votes;
              int  max_votes;

              total_votes = Frag [i] . vote [j] . deletes
                              + Frag [i] . vote [j] . a_subst
                              + Frag [i] . vote [j] . c_subst
                              + Frag [i] . vote [j] . g_subst
                              + Frag [i] . vote [j] . t_subst;
              max_votes = OVL_Max_int (Frag [i] . vote [j] . deletes,
                                   Frag [i] . vote [j] . a_subst);
              max_votes = OVL_Max_int (max_votes, Frag [i] . vote [j] . c_subst);
              max_votes = OVL_Max_int (max_votes, Frag [i] . vote [j] . g_subst);
              max_votes = OVL_Max_int (max_votes, Frag [i] . vote [j] . t_subst);

              if  (total_votes >= Degree_Threshold
                     && (double) max_votes / total_votes >= AGREEMENT_LEVEL)
                  {
                   score += MATCH_SCORE;
                   if  (score > max_score)
                       {
                        max_score = score;
                        max_sub = j;
                       }
                   consec_bad = 0;
                  }
                else
                  {
                   score += MISMATCH_SCORE;
                   consec_bad ++;
                  }
             }

           if  (max_sub >= MIN_MATCH_LEN - 1)
               left_max = max_sub;
          }

      if  (Frag [i] . high_right_degree)
          {
           score = max_score = 0.0;
           max_sub = len;
           consec_bad = 0;
           for  (j = len - 1;  j >= 0 && consec_bad < Degree_Threshold;  j --)
             {
              int  total_votes;
              int  max_votes;

              total_votes = Frag [i] . vote [j] . deletes
                              + Frag [i] . vote [j] . a_subst
                              + Frag [i] . vote [j] . c_subst
                              + Frag [i] . vote [j] . g_subst
                              + Frag [i] . vote [j] . t_subst;
              max_votes = OVL_Max_int (Frag [i] . vote [j] . deletes,
                                   Frag [i] . vote [j] . a_subst);
              max_votes = OVL_Max_int (max_votes, Frag [i] . vote [j] . c_subst);
              max_votes = OVL_Max_int (max_votes, Frag [i] . vote [j] . g_subst);
              max_votes = OVL_Max_int (max_votes, Frag [i] . vote [j] . t_subst);

              if  (total_votes >= Degree_Threshold
                     && (double) max_votes / total_votes >= AGREEMENT_LEVEL)
                  {
                   score += MATCH_SCORE;
                   if  (score > max_score)
                       {
                        max_score = score;
                        max_sub = j;
                       }
                   consec_bad = 0;
                  }
                else
                  {
                   score += MISMATCH_SCORE;
                   consec_bad ++;
                  }
             }

           if  (len - max_sub >= MIN_MATCH_LEN)
               right_min = max_sub;
          }

      Frag [i] . high_left_degree
          = Frag [i] . high_right_degree
          = Frag [i] . left_covered
          = Frag [i] . right_covered
          = Frag [i] . coupled_ends
          = Frag [i] . used = FALSE;

      if  (Verbose > 0)
          printf ("right_min = %d  left_max = %d\n", right_min, left_max);

      if  (right_min <= left_max)
          {
           if  (Verbose > 0)
               Dump_Frag (i, 0, len - 1, "AB");

           Frag [i] . left_match_len = Frag [i] . right_match_len = len;
           Frag [i] . high_left_degree
               = Frag [i] . high_right_degree = TRUE;
           Frag [i] . coupled_ends = TRUE;
           Frag [i] . left_list
               = Frag [i] . right_list = NULL;
           Frag [i] . left_list_len
               = Frag [ i] . right_list_len = 0;
          }
        else
          {
           if  (left_max >= 0)
               {
                if  (Verbose > 0)
                    Dump_Frag (i, 0, left_max, "A");

                Frag [i] . left_match_len = left_max + 1;
                Frag [i] . high_left_degree = TRUE;
                Frag [i] . left_list = NULL;
                Frag [i] . left_list_len = 0;
               }
           if  (right_min < AS_READ_MAX_LEN)
               {
                if  (Verbose > 0)
                    Dump_Frag (i, right_min, len - 1, "B");

                Frag [i] . right_match_len = len - right_min;
                Frag [i] . high_right_degree = TRUE;
                Frag [i] . right_list = NULL;
                Frag [i] . right_list_len = 0;
               }
          }
     }

   return;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   int  ch, errflg = FALSE;
   char  * p;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "d:Sv:w:")) != EOF))
     switch  (ch)
       {
        case  'd' :
          Degree_Threshold = (int) strtol (optarg, & p, 10);
          if  (p == optarg)
              {
               fprintf (stderr, "ERROR:  Illegal degree threshold \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
          break;

        case  'S' :
          Olaps_From_Store = TRUE;
          break;

        case  'v' :
          Verbose = (int) strtol (optarg, & p, 10);
          break;

        case  'w' :
          Window_Len = (int) strtol (optarg, & p, 10);
          if  (p == optarg || Window_Len <= 1)
              {
               fprintf (stderr, "ERROR:  Illegal window length \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
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

   Degree_Path = argv [optind ++];

   gkpStore_Path = argv [optind ++];

   Olap_File_Path = argv [optind ++];

   return;
  }



static int  Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Match_To_End,
     int Delta [MAX_ERRORS], int * Delta_Len)

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
   int  Delta_Stack [MAX_ERRORS];
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
   Max_Score_Len = Max_Score_Best_d = 0;
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

   (* A_End) = Max_Score_Len;
   (* T_End) = Max_Score_Len + Max_Score_Best_d;
   (* Match_To_End) = FALSE;

   return  e;
  }



static void  Process_Graph
    (void)

//  Traverse the graph of high-degree segments and create
//  supersequence chains containing them.

  {
   End_Ref_t  * ref;
   int  ref_ct = 0;
   int  start_iid, start_left;
   int  i, j;

   Num_Graph_Nodes = 0;
   for  (i = Lo_Frag_IID;  i <= Hi_Frag_IID;  i ++)
     if  (Frag [i] . high_left_degree
            || Frag [i] . high_right_degree)
         Num_Graph_Nodes ++;

   for  (i = Lo_Frag_IID;  i <= Hi_Frag_IID;  i ++)
     {
      if  (Frag [i] . high_left_degree)
          ref_ct ++;
      if  (Frag [i] . high_right_degree)
          ref_ct ++;
     }
   ref = (End_Ref_t *) safe_malloc (ref_ct * sizeof (End_Ref_t));
   ref_ct = 0;
   for  (i = Lo_Frag_IID;  i <= Hi_Frag_IID;  i ++)
     {
      if  (Frag [i] . high_left_degree)
          {
           ref [ref_ct] . iid = i;
           ref [ref_ct] . on_left = TRUE;
           ref [ref_ct] . degree = Frag [i] . left_degree;
           ref_ct ++;
          }
      if  (Frag [i] . high_right_degree)
          {
           ref [ref_ct] . iid = i;
           ref [ref_ct] . on_left = FALSE;
           ref [ref_ct] . degree = Frag [i] . right_degree;
           ref_ct ++;
          }
     }
   qsort (ref, ref_ct, sizeof (End_Ref_t), By_Degree);

   for  (i = 0;  i < ref_ct;  i ++)
     if  (ref [i] . on_left)
         {
          j = ref [i] . iid;
          if  (! Frag [j] . left_covered)
              {
               Get_Start (j, TRUE, & start_iid, & start_left);
               Traverse (start_iid, start_left);
              }
         }
       else
         {
          j = ref [i] . iid;
          if  (! Frag [j] . right_covered)
              {
               Get_Start (j, FALSE, & start_iid, & start_left);
               Traverse (start_iid, start_left);
              }
         }

   free (ref);

   return;
  }



static void  Process_Olap
    (Olap_Info_t * olap)

//  Find the alignment referred to in  olap , using the sequence
//  information in  Frag .
//  Use the alignment to increment the appropriate vote fields
//  for the a fragment.

  {
   char  * a_part, * b_part;
   int  a_part_len, b_part_len, a_end, b_end, olap_len;
   int  match_to_end, delta [MAX_ERRORS], delta_len, errors;
   int  a_offset;
   static int32  rev_id = -1;
   static char  rev_seq [AS_READ_MAX_LEN + 1] = "acgt";

   if  (Verbose > 1)
       printf ("Process_Olap:  %8d %8d %5d %5d  %c\n",
               olap -> a_iid, olap -> b_iid,
               olap -> a_hang, olap -> b_hang,
               olap -> orient);

   a_part = Frag [olap -> a_iid] . sequence;
   a_offset = 0;

   if  (olap -> a_hang > 0)
       {
        a_part += olap -> a_hang;
        a_offset += olap -> a_hang;
       }

   if  (olap -> orient == 'N')
       b_part = Frag [olap -> b_iid] . sequence;
     else
       {
        if  (rev_id != olap -> b_iid)
            {
             strcpy (rev_seq, Frag [olap -> b_iid] . sequence);
             Rev_Complement (rev_seq);
             rev_id = olap -> b_iid;
            }
        b_part = rev_seq;
       }
   if  (olap -> a_hang < 0)
       b_part -= olap -> a_hang;

   if  (Verbose > 1)
       printf ("b_part = %p  is ascii %d  rev_seq is %d\n",
               b_part, (int) (* b_part), (int) (* rev_seq));

   if  (! isalpha (* b_part) || ! isalpha (* rev_seq))
       exit (-1);

   if  (Verbose > 1)
       {
        int  j, ct;

        printf (">a_part\n");
        ct = 0;
        for  (j = 0;  a_part [j] != '\0';  j ++)
          {
           if  (ct == 60)
               {
                putchar ('\n');
                ct = 0;
               }
           putchar (a_part [j]);
           ct ++;
          }
        putchar ('\n');

        printf (">b_part\n");
        ct = 0;
        for  (j = 0;  b_part [j] != '\0';  j ++)
          {
           if  (ct == 60)
               {
                putchar ('\n');
                ct = 0;
               }
           putchar (b_part [j]);
           ct ++;
          }
        putchar ('\n');

       }

   // Get the alignment

   a_part_len = strlen (a_part);
   b_part_len = strlen (b_part);
   olap_len = OVL_Min_int (a_part_len, b_part_len);

   errors = Prefix_Edit_Dist
              (a_part, a_part_len, b_part, b_part_len,
               Error_Bound [olap_len], & a_end, & b_end,
               & match_to_end, delta, & delta_len);


   if  (Verbose > 1)
       {
        printf ("  errors = %d  delta_len = %d\n", errors, delta_len);
        Display_Alignment (a_part, a_part_len, b_part, b_part_len, delta, delta_len);
       }

   if  (errors <= Error_Bound [olap_len] && match_to_end)
       Analyze_Alignment (delta, delta_len, a_part, b_part,
                          a_end, a_offset, olap -> a_iid);
     else
       Failed_Olaps ++;
                      
   if  (Verbose > 1)
       {
        Frag_Info_t  * frag = Frag + olap -> a_iid;
        int  i;

        printf ("Votes:\n");
        for  (i = 0;  i < frag -> len;  i ++)
          printf ("%3d: %3d  %3d %3d %3d %3d  %3d %3d %3d %3d\n",
                  i,
                  frag -> vote [i] . deletes,
                  frag -> vote [i] . a_subst,
                  frag -> vote [i] . c_subst,
                  frag -> vote [i] . g_subst,
                  frag -> vote [i] . t_subst,
                  frag -> vote [i] . a_insert,
                  frag -> vote [i] . c_insert,
                  frag -> vote [i] . g_insert,
                  frag -> vote [i] . t_insert);
       }

   return;
  }



static void  Read_Degrees
    (void)

//  Read fragment overlap degrees from  Degree_Path  and store them
//  in global  Frag .

  {
   FILE  * fp;
   int32  iid;
   int  left_degree, right_degree;
   int  i, frag_size;
   int  left_ct = 0, right_ct = 0, both_ct = 0, either_ct = 0;

   fp = Local_File_Open (Degree_Path, "r");

   Lo_Frag_IID = 1;

   frag_size = 1000;
   Frag = (Frag_Info_t *) safe_malloc (frag_size * sizeof (Frag_Info_t));
   for  (i = 0;  i < frag_size;  i ++)
     {
      Frag [i] . sequence = NULL;
      Frag [i] . left_degree = Frag [i] . right_degree = 0;
      Frag [i] . need_sequence = Frag [i] . contained = FALSE;
      Frag [i] . high_left_degree = Frag [i] . high_right_degree = FALSE;
     }

   while  (fscanf (fp, "%d %d %d", & iid, & left_degree, & right_degree)
             == 3)
     {
      if  (iid >= frag_size)
          {
           int  new_size = frag_size * EXPANSION_FACTOR;

           if  (iid >= new_size)
               new_size = iid + 1;

           Frag = (Frag_Info_t *) safe_realloc (Frag,
                     new_size * sizeof (Frag_Info_t));

           for (i = frag_size;  i < new_size;  i ++)
             {
              Frag [i] . sequence = NULL;
              Frag [i] . left_degree = Frag [i] . right_degree = 0;
              Frag [i] . need_sequence = Frag [i] . contained = FALSE;
              Frag [i] . high_left_degree
                  = Frag [i] . high_right_degree = FALSE;
             }
           frag_size = new_size;
          }

      Frag [iid] . high_left_degree = (left_degree >= Degree_Threshold);
      Frag [iid] . high_right_degree = (right_degree >= Degree_Threshold);
      Frag [iid] . need_sequence = Frag [iid] . high_left_degree
                                     || Frag [iid] . high_right_degree;

if  (Frag [iid] . high_left_degree)
    {
     left_ct ++;
     either_ct ++;
     if  (Frag [iid] . high_right_degree)
         {
          both_ct ++;
          right_ct ++;
         }
    }
else if  (Frag [iid] . high_right_degree)
    {
     right_ct ++;
     either_ct ++;
    }

      if  (iid > Hi_Frag_IID)
          Hi_Frag_IID = iid;
     }

   Frag = (Frag_Info_t *) safe_realloc (Frag, (1 + Hi_Frag_IID) * sizeof (Frag_Info_t));

   fclose (fp);

printf ("left_ct = %d  right_ct = %d  both_ct = %d  either_ct = %d\n",
        left_ct, right_ct, both_ct, either_ct);

   return;
  }




static void  Read_Frags
    (void)

//  Open and read fragments with  Frag [] . need_sequence == TRUE
//  from  gkpStore_Path  and store their sequence in
//  Frag .

  {
   ReadStructp  frag_read;
   unsigned  clear_start, clear_end;
   int32  first_iid, last_iid, hi;
   int  iid, j, check_iid;

   gkpStore = openGateKeeperStore (gkpStore_Path, FALSE);
   assert (gkpStore != NULL);
   hi = getLastElemFragStore (gkpStore);
   if  (hi < Hi_Frag_IID)
       {
        fprintf (stderr,
                 "ERROR:  High overlap frag = %d  Last frag in store = %d\n",
                 Hi_Frag_IID, hi);
        exit (EXIT_FAILURE);
       }
   closeGateKeeperStore (gkpStore);

   frag_read = new_ReadStruct ();

   last_iid = 0;

   for  (first_iid = last_iid + 1;  first_iid <= Hi_Frag_IID;  first_iid ++)
     if  (Frag [first_iid] . need_sequence)
         break;

   while  (first_iid <= Hi_Frag_IID)
     {
      for  (check_iid = first_iid + 1;
              check_iid <= Hi_Frag_IID && check_iid - first_iid < MAX_MEMORY_STORE;
              check_iid ++)
        if  (Frag [check_iid] . need_sequence)
            last_iid = check_iid;

      gkpStore = openGateKeeperStore(gkpStore_Path, FALSE);
      loadGateKeeperStorePartial(gkpStore, first_iid, last_iid, FRAG_S_ALL);
   
      for  (iid = first_iid;  iid <= last_iid;  iid ++)
        {
         char  seq_buff [AS_READ_MAX_LEN + 1];
         char  qual_buff [AS_READ_MAX_LEN + 1];
         unsigned  deleted;
         int  result;

         if  (! Frag [iid] . need_sequence)
             continue;    // skip it, overlap degrees too low

         getFrag (gkpStore, iid, frag_read, FRAG_S_ALL);

         getIsDeleted_ReadStruct (frag_read, & deleted);
         if  (deleted)
             {
              Frag [iid] . sequence = NULL;
              continue;
             }

         result = getSequence_ReadStruct
                      (frag_read, seq_buff, qual_buff, AS_READ_MAX_LEN + 1);
         if  (result != 0)
             {
              fprintf (stderr,
                       "Error reading frag store tried to fit %d in a buffer of %d\n",
                       result, AS_READ_MAX_LEN + 1);
              exit (EXIT_FAILURE);
             }

         result = getClearRegion_ReadStruct
                      (frag_read, & clear_start, & clear_end, READSTRUCT_LATEST);
         if  (result != 0)
             {
              fprintf (stderr, "Error reading clear range from frag store\n");
              exit (EXIT_FAILURE);
             }

         // Make sure that we have a legal lowercase sequence string

         for  (j = clear_start;  j < clear_end;  j ++)
            seq_buff [j] = Filter (seq_buff [j]);

         seq_buff [clear_end] = '\0';

         Frag [iid] . sequence = strdup (seq_buff + clear_start);
         Frag [iid] . vote = (Vote_Tally_t *) safe_calloc (clear_end - clear_start,
                                                         sizeof (Vote_Tally_t));
         Frag [iid] . len = clear_end - clear_start;
        }

      closeGateKeeperStore (gkpStore);

      for  (first_iid = check_iid;  first_iid <= Hi_Frag_IID;  first_iid ++)
        if  (Frag [first_iid] . need_sequence)
            break;
     }

   delete_ReadStruct (frag_read);

   return;
  }



static void  Read_Olaps
    (void)

//  Open and read overlaps from  Olap_File_Path  and if they
//  join high-degree fragment ends, record the alignment votes
//  for the first fragment.  Overlaps are in the format produced by
//  get-olaps and each overlap must appear twice, once in each order.

  {
   FILE  * fp = NULL;
   OVL_Store_t  * ovl_store = NULL;
   OVL_Stream_t  * ovl_stream = NULL;
   Olap_Info_t  olap;
   long int  ct = 0;

   if  (Olaps_From_Store)
       {
        ovl_store = New_OVL_Store ();
        if  (Open_OVL_Store (ovl_store, Olap_File_Path) != 0)
            {
             fprintf (stderr, "ERROR:  Failed to open overlap store \"%s\"\n",
                      Olap_File_Path);
             exit (EXIT_FAILURE);
            }
        ovl_stream = New_OVL_Stream ();
        Init_OVL_Stream (ovl_stream, 1, Last_Frag_In_OVL_Store (ovl_store),
                         ovl_store);
       }
     else
       fp = Local_File_Open (Olap_File_Path, "r");

   while  (Get_Olap (& olap, fp, ovl_stream))
     {
      if  (olap . a_iid < Lo_Frag_IID || Hi_Frag_IID < olap . a_iid
             || olap . b_iid < Lo_Frag_IID || Hi_Frag_IID < olap . b_iid)
          continue;
             
      if  (Frag [olap . a_iid] . contained)
          continue;
      if  (Frag [olap . a_iid] . sequence == NULL
             || Frag [olap . b_iid] . sequence == NULL)
          continue;

      if  (olap . orient == 'O')
          {
           int  save;

           save = olap . a_hang;
           olap . a_hang = - olap . b_hang;
           olap . b_hang = - save;
           olap . orient = 'I';
          }

      if  (olap . a_hang <= 0 && olap . b_hang >= 0)
          {
           Frag [olap . a_iid] . contained = TRUE;
           continue;
          }

      switch (olap . orient)
        {
         case  'I' :
           if  ((olap . a_hang >= 0
                      && ! Frag [olap . b_iid] . high_right_degree)
                 || (olap . a_hang <= 0
                      && ! Frag [olap . a_iid] . high_left_degree)
                 || (olap . b_hang >= 0
                      && ! Frag [olap . a_iid] . high_right_degree)
                 || (olap . b_hang <= 0
                      && ! Frag [olap . b_iid] . high_left_degree))
               continue;
           break;
         case  'N' :
           if  ((olap . a_hang >= 0
                      && ! Frag [olap . b_iid] . high_left_degree)
                 || (olap . a_hang <= 0
                      && ! Frag [olap . a_iid] . high_left_degree)
                 || (olap . b_hang >= 0
                      && ! Frag [olap . a_iid] . high_right_degree)
                 || (olap . b_hang <= 0
                      && ! Frag [olap . b_iid] . high_right_degree))
               continue;
           break;
         default :
           fprintf (stderr, "ERROR:  Unexpected orientation %c\n",
                    olap . orient);
           exit (EXIT_FAILURE);
        }

      ct ++;

      Process_Olap (& olap);
     }

   if  (Olaps_From_Store)
       {
        Free_OVL_Stream (ovl_stream);
        Free_OVL_Store (ovl_store);
       }
     else
       fclose (fp);

   return;
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



static void  Traverse
    (int iid, int use_left)

//  Traverse graph starting at node  iid  over end indicated by  use_left .
//  Create and output sequence represented by overlap path and
//  eliminate nodes on it and with sequence contained in this sequence.

  {
   Graph_Edge_t  * edge = NULL, * best_edge = NULL;
   static char  * cum_buff = NULL, * seq_buff = NULL;
   static int  cum_size, seq_size;
   char  header [MAX_HDR_LEN];
   int  start_iid = iid;
   int  cum_len, seq_len;
   int  j, best_progress, best_olap, extension_ct = 0;

   if  (cum_buff == NULL)
       {
        cum_size = INITIAL_BUFFER_SIZE;
        cum_buff = (char *) safe_malloc (cum_size);
        seq_size = INITIAL_BUFFER_SIZE;
        seq_buff = (char *) safe_malloc (seq_size);
       }

   Get_Consensus (iid, use_left, use_left, & seq_buff, & seq_size);

   cum_len = strlen (seq_buff);
   if  (cum_len >= cum_size)
       {
        cum_size = 1 + cum_len;
        cum_buff = (char *) safe_realloc (cum_buff, cum_size);
       }
   strcpy (cum_buff, seq_buff);

   if  (Verbose > 0)
       {
        printf ("Traverse from %d over %s end  len = %d\n",
                iid, use_left ? "left" : "right", cum_len);
        Fasta_Print (stdout, seq_buff, NULL);
       }

   do
     {
      Overlap  * result;
      char  * cat_seq;
      int  best_iid;
      int  expected_olap, min_ahang, max_ahang;

      if  (use_left)
          Frag [iid] . left_covered = TRUE;
        else
          Frag [iid] . right_covered = TRUE;

      best_progress = best_olap = 0;
      if  (use_left)
          {
           for  (j = 0;  j < Frag [iid] . left_list_len;  j ++)
             {
              edge = Frag [iid] . left_list + j;

              if  ((edge -> over_left && Frag [edge -> iid] . left_covered)
                     || (! edge -> over_left
                           && Frag [edge -> iid] . right_covered))
                  continue;
              if  ((edge -> progress >= TARGET_PROGRESS
                          && edge -> olap_len > best_olap)
                     || (edge -> progress < TARGET_PROGRESS
                               && edge -> progress > best_progress))
                  {
                   best_progress = edge -> progress;
                   best_olap = edge -> olap_len;
                   best_edge = edge;
                  }
             }
//  Eliminate fragment ends contained within this overlap.
//  Temporary kludge using progress--would be better to use
//  actual overlaps.
           for  (j = 0;  j < Frag [iid] . left_list_len;  j ++)
             if  (edge -> progress <= best_progress)
                 {
                  if  (Frag [edge -> iid] . coupled_ends)
                      Frag [edge -> iid] . left_covered
                          = Frag [edge -> iid] . right_covered = TRUE;
                  else if  (edge -> over_left)
                      Frag [edge -> iid] . left_covered = TRUE;
                    else
                      Frag [edge -> iid] . right_covered = TRUE;
                 }
          }
        else
          {
           for  (j = 0;  j < Frag [iid] . right_list_len;  j ++)
             {
              edge = Frag [iid] . right_list + j;

              if  ((edge -> over_left && Frag [edge -> iid] . left_covered)
                     || (! edge -> over_left
                           && Frag [edge -> iid] . right_covered))
                  continue;
              if  ((edge -> progress >= TARGET_PROGRESS
                          && edge -> olap_len > best_olap)
                     || (edge -> progress < TARGET_PROGRESS
                               && edge -> progress > best_progress))
                  {
                   best_progress = edge -> progress;
                   best_olap = edge -> olap_len;
                   best_edge = edge;
                  }
             }
//  Eliminate fragment ends contained within this overlap.
//  Temporary kludge using progress--would be better to use
//  actual overlaps.
           for  (j = 0;  j < Frag [iid] . right_list_len;  j ++)
             if  (edge -> progress <= best_progress)
                 {
                  if  (Frag [edge -> iid] . coupled_ends)
                      Frag [edge -> iid] . left_covered
                          = Frag [edge -> iid] . right_covered = TRUE;
                  else if  (edge -> over_left)
                      Frag [edge -> iid] . left_covered = TRUE;
                    else
                      Frag [edge -> iid] . right_covered = TRUE;
                 }
          }
      if  (best_progress == 0)
          break;

      best_iid = best_edge -> iid;

      Get_Consensus (best_iid, best_edge -> over_left,
                     ! best_edge -> over_left,
                     & seq_buff, & seq_size);
      if  (best_edge -> over_left)
          expected_olap = Frag [best_iid] . left_match_len
                            - best_edge -> progress;
        else
          expected_olap = Frag [best_iid] . right_match_len
                            - best_edge -> progress;
      max_ahang = cum_len - 0.94 * expected_olap;
      min_ahang = cum_len - 1.06 * expected_olap;

      if  (Verbose > 0)
          printf ("Find_Overlap %d -> %d\n", iid, best_iid);
      result = Find_Overlap
                  (cum_buff, seq_buff, AB_AB, min_ahang, max_ahang,
                   CGW_DP_ERATE, CGW_DP_THRESH, CGW_DP_MINLEN,
                   AS_FIND_ALIGN);
      if  (result == NULL)
          {
           int  offset;
           
           printf ("BUMMER:  No overlap found %d -> %d  expected = %d  progress = %d\n",
                   iid, best_iid, expected_olap, best_edge -> progress);
           if  (Verbose > 0)
               {
                offset = min_ahang;
                if  (offset < 0)
                    offset = 0;
                Fasta_Print (stdout, cum_buff + offset , "From");
                Fasta_Print (stdout, seq_buff , "To");
               }

           // Stop extending here  If some other fragment overlaps
           // we'll get it later and  remove-dup-screen  will combine
           // the two regions

           break;
          }

      if  (Verbose > 0)
          printf ("  endpos = %d\n", result -> endpos);

      if  (result -> endpos <= 0)
          {
           printf ("  Olap %d -> %d  progress shrank from %d to %d\n",
                   iid, best_iid, best_edge -> progress, result -> endpos);
           // Eliminate this fragment; it's contained
           Frag [best_iid] . left_covered
               = Frag [best_iid] . right_covered = TRUE;
           continue;
          }
      if  (abs (best_edge -> progress - result -> endpos) > MAX_PROGRESS_DIFF)
          {
           printf ("  Olap %d -> %d  progress est/actual = %d/%d  eliminate %d\n",
                   iid, best_iid, best_edge -> progress, result -> endpos,
                   best_iid);
           // Eliminate  best_iid ; it's a problem
           Frag [best_iid] . left_covered
               = Frag [best_iid] . right_covered = TRUE;
           continue;
          }
          
      // Eliminate other nodes contained in the overlap between these two
      // Temporary kludge version
      if  (best_edge -> over_left)
          {
           for  (j = 0;  j < Frag [best_iid] . left_list_len;  j ++)
             {
              edge = Frag [best_iid] . left_list + j;
              if  (Frag [edge -> iid] . coupled_ends)
                  {
                   Frag [edge -> iid] . left_covered
                       = Frag [edge -> iid] . right_covered = TRUE;
                  }
              else if  (edge -> over_left)
                  Frag [edge -> iid] . left_covered = TRUE;
                else
                  Frag [edge -> iid] . right_covered = TRUE;
             }
          }
        else
          {
           for  (j = 0;  j < Frag [best_iid] . right_list_len;  j ++)
             {
              edge = Frag [best_iid] . right_list + j;
              if  (Frag [edge -> iid] . coupled_ends)
                  {
                   Frag [edge -> iid] . left_covered
                       = Frag [edge -> iid] . right_covered = TRUE;
                  }
              else if  (edge -> over_left)
                  Frag [edge -> iid] . left_covered = TRUE;
                else
                  Frag [edge -> iid] . right_covered = TRUE;
             }
          }
      
      seq_len = strlen (seq_buff);
      cat_seq = seq_buff + seq_len - result -> endpos;
      seq_len = strlen (cat_seq);
      if  (cum_len + seq_len + 1 > cum_size)
          {
           cum_size = OVL_Max_int ((int) (cum_size * EXPANSION_FACTOR),
                               cum_len + seq_len + 1);
           cum_buff = (char *) safe_realloc (cum_buff, cum_size);
          }
      strcpy (cum_buff + cum_len, cat_seq);
      cum_len += seq_len;

      iid = best_iid;
      extension_ct ++;
      if  (Verbose > 0)
          {
           printf ("  to %d over %s end  est/actual prog = %d/%d"
                   "  cum_len = %d  len = %d\n",
                   iid, best_edge -> over_left ? "left" : "right",
                   best_edge -> progress, result -> endpos, cum_len, seq_len);
           Fasta_Print (stdout, cat_seq, NULL);
          }

      if  (best_edge -> over_left)
          Frag [iid] . left_covered = TRUE;
        else
          Frag [iid] . right_covered = TRUE;

      if  (Frag [iid] . coupled_ends)
          use_left = ! best_edge -> over_left;     // keep going
        else
          break;                                   // stop
     }  while  (iid > 0);

   if  (extension_ct == 0)
       sprintf (header, "iid %d  len = %d", start_iid, cum_len);
     else
       sprintf (header, "iids %d..%d  len = %d", start_iid, iid, cum_len);

   // Skip ones that are too short (presumably caused by deletions in getting
   // consensus
   if  (cum_len >= MIN_MATCH_LEN)
       Fasta_Print (Screen_File, cum_buff, header);

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s [-d DegrThresh] [-w WindowLen] [-v VerboseLevel]\n"
           "           <DegreeFile> <FragStore> <OlapFile>\n"
           "\n"
           "Finds high-overlap consensus sequences and creates\n"
           "a screen library from them.  <DegreeFile> has the\n"
           "overlap degree at both ends of each fragment.\n"
           "Fragments come from <FragStore>.  Overlaps are in\n"
           "sorted order in <OlapFile>.\n"
           "\n"
           "Options:\n"
           "-d   set minimum overlap degree to be screened\n"
           "-v   set the level of verbose debugging messages to print\n"
           "-w   set length of window to scan consensus regions\n",
           command);

   return;
  }


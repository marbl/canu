
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
* Module:  FragCorrectOVL.c
* Description:
*   Based on overlaps between DNA fragment sequences, make corrections
*   to single bases in the sequences.
* 
*    Programmer:  A. Delcher
*       Started:   4 Dec 2000
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: FragCorrectOVL.c,v 1.14 2007-02-24 15:42:33 brianwalenz Exp $
 * $Revision: 1.14 $
*/

static char CM_ID[] = "$Id: FragCorrectOVL.c,v 1.14 2007-02-24 15:42:33 brianwalenz Exp $";


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
#include  <pthread.h>


//  Local include files

#include  "AS_OVL_delcher.h"
#include  "AS_PER_gkpStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_UTL_version.h"
#include  "FragCorrectOVL.h"



//  Constants

#define  BRANCH_PT_MATCH_VALUE    0.272
    //  Value to add for a match in finding branch points
    //  1.20 was the calculated value for 6% vs 35% error discrimination
    //  Converting to integers didn't make it faster
#define  BRANCH_PT_ERROR_VALUE    -0.728
    //  Value to add for a mismatch in finding branch points
    //   -2.19 was the calculated value for 6% vs 35% error discrimination
    //  Converting to integers didn't make it faster
#define  DEFAULT_CORRECTION_FILENAME  "frag.cor"
    //  Default name of file where corrections are sent
#define  DEFAULT_DEGREE_THRESHOLD    2
    //  Default value for  Degree_Threshold
#define  DEFAULT_END_EXCLUDE_LEN     3
    //  Default value for  End_Exclude_Len
#define  DEFAULT_KMER_LEN            9
    //  Default value for  Kmer_Len
#define  DEFAULT_NUM_PTHREADS        2
    //  Default number of pthreads to use
#define  DEFAULT_VOTE_QUALIFY_LEN    9
    //  Default value for bases surrounding SNP to vote for change
#define  EDIT_DIST_PROB_BOUND        1e-4
    //  Probability limit to "band" edit-distance calculation
    //  Determines  NORMAL_DISTRIB_THOLD
#define  ERATE_BITS                  16
    //  Number of bits to store integer versions of error rates
#define  ERRORS_FOR_FREE             1
    //  The number of errors that are ignored in setting probability
    //  bound for terminating alignment extensions in edit distance
    //  calculations
#define  EXPANSION_FACTOR            1.4
    //  Factor by which to grow memory in olap array when reading it
#define  FRAGS_PER_BATCH             100000
    //  Number of old fragments to read into memory-based fragment
    //  store at a time for processing
#define  MAX_ERROR_RATE              AS_GUIDE_ERROR_RATE
    //  The largest error allowed in overlaps
#define  MAX_FILENAME_LEN            1000
    //  Longest name allowed for a file in the overlap store
#define  MAX_FRAG_LEN                2048
    //  The longest fragment allowed
#define  MAX_ERRORS                  (1 + (int) (MAX_ERROR_RATE * MAX_FRAG_LEN))
    //  Most errors in any edit distance computation
#define  MAX_DEGREE                  32767
    //  Highest number of votes before overflow
#define  MAX_VOTE                    255
    //  Highest number of votes before overflow
#define  MIN_BRANCH_END_DIST         20
    //  Branch points must be at least this many bases from the
    //  end of the fragment to be reported
#define  MIN_BRANCH_TAIL_SLOPE       0.20
    //  Branch point tails must fall off from the max by at least
    //  this rate
#define  MIN_HAPLO_OCCURS            3
    //  This many or more votes at the same base indicate
    //  a separate haplotype
#define  NORMAL_DISTRIB_THOLD        3.62
    //  Determined by  EDIT_DIST_PROB_BOUND
#define  THREAD_STACKSIZE        (16 * 512 * 512)
    //  The amount of memory to allocate for the stack of each thread

#define USE_STORE_DIRECTLY_READ
    //  Use the store directly during the initial load -- we aren't doing
    //  random access, just streaming through and loading.  This lets us
    //  load about 2x the frags.

//#define USE_STORE_DIRECTLY_STREAM
    //  Use the store directly during the stream -- good if you don't have
    //  lots of frags loaded.

//#define USE_STREAM_FOR_EXTRACT
    //  When loading frags during the stream, use a fragStream instead
    //  of random access to the store.  Useless unless USE_STORE_DIRECTLY_STREAM
    //  is enabled.

//  Type definitions

typedef  struct
  {
   unsigned  int confirmed : 8;
   unsigned  int deletes : 8;
   unsigned  int a_subst : 8;
   unsigned  int c_subst : 8;
   unsigned  int g_subst : 8;
   unsigned  int t_subst : 8;
   unsigned  int no_insert : 8;
   unsigned  int a_insert : 8;
   unsigned  int c_insert : 8;
   unsigned  int g_insert : 8;
   unsigned  int t_insert : 8;
  }  Vote_Tally_t;

typedef  struct
  {
   int  frag_sub;
   int  align_sub;
   Vote_Value_t  vote_val;
  }  Vote_t;

typedef  struct
  {
   char  * sequence;
   Vote_Tally_t  * vote;
   unsigned  clear_len : 15;
   unsigned  left_degree : 15;
   unsigned  right_degree : 15;
   unsigned  shredded : 1;    // True if shredded read
   unsigned  unused : 1;
  }  Frag_Info_t;

const int  INNIE = 0;
const int  NORMAL = 1;

typedef  struct                 
  {
   int32  a_iid, b_iid;
   signed int  a_hang : 15;
   signed int  b_hang : 15;
   signed int  orient : 2;
  }  Olap_Info_t;

typedef  struct
  {
   unsigned  b_iid : 31;
   unsigned  flipped : 1;
   signed int  a_hang : 16;
   signed int  b_hang : 16;
   unsigned  orig_erate : ERATE_BITS;   // original error rate (aka quality)
   unsigned  corr_erate : ERATE_BITS;   // error rate after fragment correction
  }  Short_Olap_Data_t;

typedef  struct
  {
   unsigned  id : 31;
   unsigned  shredded : 1;
   int  start;              // position of beginning of sequence in  buffer
  }  Frag_List_Entry_t;

typedef  struct
  {
   Frag_List_Entry_t  * entry;
   char  * buffer;
   int  size, ct, buffer_size;
  }  Frag_List_t;

typedef  struct
  {
   int  thread_id;
   int32  lo_frag, hi_frag;
   int  next_olap;
   FragStream  *frag_stream;
   fragRecord *frag_read;
   Frag_List_t  * frag_list;
   char  rev_seq [AS_READ_MAX_LEN + 1];
   int  rev_id;
   int  * edit_array [MAX_ERRORS];
   int  edit_space [(MAX_ERRORS + 4) * MAX_ERRORS];
  }  Thread_Work_Area_t;



//  Static Globals

static char  * Correction_Filename = DEFAULT_CORRECTION_FILENAME;
    // Name of file to which correction information is sent
static int  Degree_Threshold = DEFAULT_DEGREE_THRESHOLD;
    // Set keep flag on end of fragment if number of olaps < this value
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
static int  Extend_Fragments = FALSE;
    // If true, try to extend clear range of fragments.
    // Set by  -e  option
static int  Failed_Olaps = 0;
    // Counts overlaps that didn't make the error bound
static Frag_Info_t  * Frag;
    // Sequence and vote information for current range of fragments
    // being corrected
static Frag_List_t  Frag_List;
    // List of ids and sequences of fragments with overlaps to fragments
    // in  Frag .  Allows simultaneous access by threads.
static GateKeeperStore  *gkpStore;
    // Fragment store from which fragments are loaded
static FragStream  *Frag_Stream;
    // Stream to extract fragments from internal store
static char  * gkpStore_Path;
    // Name of directory containing fragment store from which to get fragments
static int32  Hi_Frag_IID;
    // Internal ID of last fragment in frag store to process
static GateKeeperStore  *Internal_gkpStore;
    // Holds partial frag store to be processed simultanously by
    // multiple threads
static int  Kmer_Len = DEFAULT_KMER_LEN;
    // Length of minimum exact match in overlap to confirm base pairs
static int32  Lo_Frag_IID;
    // Internal ID of first fragment in frag store to process
static time_t  Now;
    // Used to get current time
static int  Num_Frags;
    // Number of fragments being corrected
static int  Num_Olaps;
    // Number of overlaps being used
static int  Num_PThreads = DEFAULT_NUM_PTHREADS;
    // Number of pthreads to process overlaps/corrections;
static Olap_Info_t  * Olap = NULL;
    // Array of overlaps being used
static int  Olaps_From_Store = FALSE;
    // Indicates if overlap info comes from  get-olaps  or from
    // a binary overlap store
static char  * Olap_Path;
    // Name of file containing a sorted list of overlaps
static pthread_mutex_t  Print_Mutex;
    // To make debugging printout come out together
static int  Use_Haplo_Ct = TRUE;
    // Set false by  -h  option to ignore haplotype counts
    // when correcting
static int  Verbose_Level = 0;
    // Determines number of extra outputs
static int  Vote_Qualify_Len = DEFAULT_VOTE_QUALIFY_LEN;
    // Number of bases surrounding a SNP to vote for change



//  Static Functions

static void  Analyze_Alignment
    (int delta [], int delta_len, char * a_part, char * b_part,
     int a_len, int b_len, int a_offset, int sub);
static int  Binomial_Bound
    (int, double, int, double);
static int  By_B_IID
    (const void * a, const void * b);
static void  Cast_Vote
    (Vote_Value_t val, int p, int sub);
static char  Complement
    (char);
static void  Compute_Delta
    (int delta [], int * delta_len, int * edit_array [MAX_ERRORS],
     int e, int d, int row);
static void  Display_Alignment
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct,
     int capitalize_start);
static void  Display_Frags
    (void);
static void  Extract_Needed_Frags
    (GateKeeperStore *store, int32 lo_frag, int32 hi_frag,
     Frag_List_t * list, int * next_olap);
static char  Filter
    (char ch);
static void  Get_Olaps_From_Store
    (char * path, int32 lo_id, int32 hi_id, Olap_Info_t * * olap, int * num);
static void  Init_Frag_List
    (Frag_List_t * list);
static void  Initialize_Globals
    (void);
static void  Init_Thread_Work_Area
    (Thread_Work_Area_t * wa, int id);
static Vote_Value_t  Matching_Vote
    (char ch);
static int  OVL_Max_int
    (int a, int b);
static int  OVL_Min_int
    (int a, int b);
static void  Output_Corrections
    (FILE * fp);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Process_Olap
    (Olap_Info_t * olap, char * b_seq, char * rev_seq, int * rev_id,
     int shredded, Thread_Work_Area_t * wa);
static void  Read_Frags
    (void);
static void  Read_Olaps
    (void);
static void  Rev_Complement
    (char * s);
static void  Stream_Old_Frags
    (void);
static int  Sign
    (int a);
void *  Threaded_Process_Stream
    (void * ptr);
static void  Threaded_Stream_Old_Frags
    (void);
static void  Usage
    (char * command);


int  main
    (int argc, char * argv [])

  {
   FILE  * fp;

   Parse_Command_Line  (argc, argv);

   Now = time (NULL);
   fprintf (stderr, "### Starting at  %s", ctime (& Now));

   Initialize_Globals ();

   gkpStore = openGateKeeperStore(gkpStore_Path, FALSE);
   assert (gkpStore != NULL);

   fprintf (stderr, "Starting Read_Frags ()\n");
   Read_Frags ();

   fprintf (stderr, "Starting Read_Olaps ()\n");
   Read_Olaps ();

   fprintf (stderr, "Before sort\n");
   qsort (Olap, Num_Olaps, sizeof (Olap_Info_t), By_B_IID);

   if  (Verbose_Level > 2)
       {
        int  i;

        for  (i = 0;  i < Num_Olaps;  i ++)
          printf ("%8d %8d %5d %5d  %c\n",
                  Olap [i] . a_iid, Olap [i] . b_iid,
                  Olap [i] . a_hang, Olap [i] . b_hang,
                  Olap [i] . orient == INNIE ? 'I' : 'N');
       }

   if  (Num_Olaps > 0)
       {
        fprintf (stderr, "Before Stream_Old_Frags  Num_Olaps = %d\n", Num_Olaps);
        if  (Num_PThreads > 0)
            Threaded_Stream_Old_Frags ();
          else
            Stream_Old_Frags ();
        fprintf (stderr, "                   Failed overlaps = %d\n", Failed_Olaps);
       }

   closeGateKeeperStore (gkpStore);

   if  (Verbose_Level > 1)
       {
        int  i, j;

        for  (i = 0;  i < Num_Frags;  i ++)
          {
           printf (">%d\n", Lo_Frag_IID + i);
           for  (j = 0;  Frag [i] . sequence [j] != '\0';  j ++)
             printf ("%3d: %c  %3d  %3d | %3d %3d %3d %3d | %3d %3d %3d %3d %3d\n",
                     j,
                     j >= Frag [i] . clear_len ?
                         toupper (Frag [i] . sequence [j]) : Frag [i] . sequence [j],
                     Frag [i] . vote [j] . confirmed,
                     Frag [i] . vote [j] . deletes,
                     Frag [i] . vote [j] . a_subst,
                     Frag [i] . vote [j] . c_subst,
                     Frag [i] . vote [j] . g_subst,
                     Frag [i] . vote [j] . t_subst,
                     Frag [i] . vote [j] . no_insert,
                     Frag [i] . vote [j] . a_insert,
                     Frag [i] . vote [j] . c_insert,
                     Frag [i] . vote [j] . g_insert,
                     Frag [i] . vote [j] . t_insert);
          }
       }

   fprintf (stderr, "Before Output_Corrections  Num_Frags = %d\n", Num_Frags);
   fp = File_Open (Correction_Filename, "wb");
   Output_Corrections (fp);
   fclose (fp);

   Now = time (NULL);
   fprintf (stderr, "### Finished at  %s", ctime (& Now));

   return  0;
  }



static void  Analyze_Alignment
    (int delta [], int delta_len, char * a_part, char * b_part,
     int  a_len, int b_len, int a_offset, int sub)

//  Analyze the delta-encoded alignment in  delta [0 .. (delta_len - 1)]
//  between  a_part  and  b_part  and store the resulting votes
//  about the a sequence in  Frag [sub] .  The alignment starts
//   a_offset  bytes in from the start of the a sequence in  Frag [sub] .
//   a_len  and  b_len  are the lengths of the prefixes of  a_part  and
//   b_part , resp., that align.

  {
   int  prev_match, next_match;
   Vote_t  vote [MAX_FRAG_LEN];
   int  ct;
   int  i, j, k, m, p;


   if  (a_len < 0 || b_len < 0)
       {
        fprintf (stderr, "ERROR:  a_len = %d  b_len = %d  sub = %d\n",
                 a_len, b_len, sub);
        exit (-3);
       }

   vote [0] . frag_sub = -1;
   vote [0] . align_sub = -1;
   vote [0] . vote_val = A_SUBST;   // Dummy value
   ct = 1;
   i = j = p = 0;

   for  (k = 0;  k < delta_len;  k ++)
     {
      for  (m = 1;  m < abs (delta [k]);  m ++)
        {
         if  (a_part [i] != b_part [j])
             {
              vote [ct] . frag_sub = i;
              vote [ct] . align_sub = p;
              switch  (b_part [j])
                {
                 case  'a' :
                   vote [ct] . vote_val = A_SUBST;
                   break;
                 case  'c' :
                   vote [ct] . vote_val = C_SUBST;
                   break;
                 case  'g' :
                   vote [ct] . vote_val = G_SUBST;
                   break;
                 case  't' :
                   vote [ct] . vote_val = T_SUBST;
                   break;
                 default :
                   fprintf (stderr, "ERROR:  [1] Bad sequence char \'%c\' (ASCII %d)\n",
                            b_part [j], (int) b_part [j]);
                   exit (1);
                }
              ct ++;           
             }
         i ++;
         j ++;
         p ++;
        }
      if  (delta [k] < 0)
          {
           vote [ct] . frag_sub = i - 1;
           vote [ct] . align_sub = p;
           switch  (b_part [j])
              {
               case  'a' :
                 vote [ct] . vote_val = A_INSERT;
                 break;
               case  'c' :
                 vote [ct] . vote_val = C_INSERT;
                 break;
               case  'g' :
                 vote [ct] . vote_val = G_INSERT;
                 break;
               case  't' :
                 vote [ct] . vote_val = T_INSERT;
                 break;
               default :
                 fprintf (stderr, "ERROR:  [2] Bad sequence char \'%c\' (ASCII %d)\n",
                          b_part [j], (int) b_part [j]);
                   exit (2);
              }
           ct ++;
           j ++;
           p ++;
          }
        else
          {
           vote [ct] . frag_sub = i;
           vote [ct] . align_sub = p;
           vote [ct] . vote_val = DELETE;
           ct ++;
           i ++;
           p ++;
          }
     }

   while  (i < a_len)
     {
      if  (a_part [i] != b_part [j])
          {
           vote [ct] . frag_sub = i;
           vote [ct] . align_sub = p;
           switch  (b_part [j])
             {
              case  'a' :
                vote [ct] . vote_val = A_SUBST;
                break;
              case  'c' :
                vote [ct] . vote_val = C_SUBST;
                break;
              case  'g' :
                vote [ct] . vote_val = G_SUBST;
                break;
              case  't' :
                vote [ct] . vote_val = T_SUBST;
                break;
              default :
                fprintf (stderr, "ERROR:  [3] Bad sequence char \'%c\' (ASCII %d)\n",
                         b_part [j], (int) b_part [j]);
                fprintf (stderr, "i = %d  a_len = %d  j = %d  b_len = %d\n",
                         i, a_len, j, b_len);
                exit (3);
             }
           ct ++;           
          }
      i ++;
      j ++;
      p ++;
     }

   vote [ct] . frag_sub = i;
   vote [ct] . align_sub = p;

   for  (i = 1;  i <= ct;  i ++)
     {
      int  k, p_lo, p_hi;

      prev_match = vote [i] . align_sub - vote [i - 1] . align_sub - 1;
      p_lo = (i == 1 ? 0 : End_Exclude_Len);
      p_hi = (i == ct ? prev_match : prev_match - End_Exclude_Len);
      if  (prev_match >= Kmer_Len)
          {
           for  (p = 0;  p < p_lo;  p ++)
             Cast_Vote
                 (Matching_Vote (a_part [vote [i - 1] . frag_sub + p + 1]),
                  a_offset + vote [i - 1] . frag_sub + p + 1, sub);

           for  (p = p_lo;  p < p_hi;  p ++)
             {
              k = a_offset + vote [i - 1] . frag_sub + p + 1;
              if  (Frag [sub] . vote [k] . confirmed < MAX_VOTE)
                  Frag [sub] . vote [k] . confirmed ++;
              if  (p < p_hi - 1
                       && Frag [sub] . vote [k] . no_insert < MAX_VOTE)
                  Frag [sub] . vote [k] . no_insert ++;
             }

           for  (p = p_hi;  p < prev_match;  p ++)
             Cast_Vote
                 (Matching_Vote (a_part [vote [i - 1] . frag_sub + p + 1]),
                  a_offset + vote [i - 1] . frag_sub + p + 1, sub);
          }
      if  (i < ct
            && (prev_match > 0
                  || vote [i - 1] . vote_val <= T_SUBST
                  || vote [i] . vote_val <= T_SUBST))
               // Don't allow consecutive inserts
          {
           next_match = vote [i + 1] . align_sub - vote [i] . align_sub - 1;
           if  (prev_match + next_match >= Vote_Qualify_Len)
               Cast_Vote (vote [i] . vote_val, a_offset + vote [i] . frag_sub, sub);
          }
     }

   if  (Verbose_Level > 0)
       {
        int  ct = 0;

        printf (">a_part\n");
        for  (j = 0;  a_part [j] != '\0';  j ++)
          {
           if  (ct == 60)
               {
                putchar ('\n');
                ct = 0;
               }
           if  (ct == 0)
               printf ("   ");
           putchar (Frag [sub] . vote [a_offset + j] . confirmed ? '*' : ' ');
           ct ++;
          }
        putchar ('\n');
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



static void  Cast_Vote
    (Vote_Value_t val, int p, int sub)

//  Add vote  val  to  Frag [sub]  at sequence position  p

  {
   switch  (val)
     {
      case  DELETE :
        if  (Frag [sub] . vote [p] . deletes < MAX_VOTE)
            Frag [sub] . vote [p] . deletes ++;
        break;
      case  A_SUBST :
        if  (Frag [sub] . vote [p] . a_subst < MAX_VOTE)
            Frag [sub] . vote [p] . a_subst ++;
        break;
      case  C_SUBST :
        if  (Frag [sub] . vote [p] . c_subst < MAX_VOTE)
            Frag [sub] . vote [p] . c_subst ++;
        break;
      case  G_SUBST :
        if  (Frag [sub] . vote [p] . g_subst < MAX_VOTE)
            Frag [sub] . vote [p] . g_subst ++;
        break;
      case  T_SUBST :
        if  (Frag [sub] . vote [p] . t_subst < MAX_VOTE)
            Frag [sub] . vote [p] . t_subst ++;
        break;
      case  A_INSERT :
        if  (Frag [sub] . vote [p] . a_insert < MAX_VOTE)
            Frag [sub] . vote [p] . a_insert ++;
        break;
      case  C_INSERT :
        if  (Frag [sub] . vote [p] . c_insert < MAX_VOTE)
            Frag [sub] . vote [p] . c_insert ++;
        break;
      case  G_INSERT :
        if  (Frag [sub] . vote [p] . g_insert < MAX_VOTE)
            Frag [sub] . vote [p] . g_insert ++;
        break;
      case  T_INSERT :
        if  (Frag [sub] . vote [p] . t_insert < MAX_VOTE)
            Frag [sub] . vote [p] . t_insert ++;
        break;
      case  NO_VOTE :
        // do nothing
        break;
      default :
        fprintf (stderr, "ERROR:  Illegal vote type\n");
     }

   return;
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
      case  'n' :
        return  'n';
      default :
        fprintf (stderr, "ERROR(complement):  Unexpected character `%c\'\n", ch);
        exit (-1);
     }

   return  'x';    // Just to make the compiler happy
  }



static void  Compute_Delta
    (int delta [], int * delta_len, int * edit_array [MAX_ERRORS],
     int e, int d, int row)

//  Set  delta  to the entries indicating the insertions/deletions
//  in the alignment encoded in  edit_array  ending at position
//  edit_array [e] [d] .   row  is the position in the first
//  string where the alignment ended.  Set  (* delta_len)  to
//  the number of entries in  delta .

  {
   int  delta_stack [MAX_ERRORS];
   int  from, last, max;
   int  i, j, k;

   last = row;
   (* delta_len) = 0;

   for  (k = e;  k > 0;  k --)
     {
      from = d;
      max = 1 + edit_array [k - 1] [d];
      if  ((j = edit_array [k - 1] [d - 1]) > max)
          {
           from = d - 1;
           max = j;
          }
      if  ((j = 1 + edit_array [k - 1] [d + 1]) > max)
          {
           from = d + 1;
           max = j;
          }
      if  (from == d - 1)
          {
           delta_stack [(* delta_len) ++] = max - last - 1;
           d --;
           last = edit_array [k - 1] [from];
          }
      else if  (from == d + 1)
          {
           delta_stack [(* delta_len) ++] = last - (max - 1);
           d ++;
           last = edit_array [k - 1] [from];
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



#define  DISPLAY_WIDTH   60

static void  Display_Alignment
    (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct,
     int capitalize_start)

//  Show (to  stdout ) the alignment encoded in  delta [0 .. (delta_ct - 1)]
//  between strings  a [0 .. (a_len - 1)]  and  b [0 .. (b_len - 1)] .
//  Capitialize  a  characters for positions at and after  capitalize_start .

  {
   int  i, j, k, m, top_len, bottom_len;
   char  top [2000], bottom [2000];

   i = j = top_len = bottom_len = 0;
   for  (k = 0;  k < delta_ct;  k ++)
     {
      for  (m = 1;  m < abs (delta [k]);  m ++)
        {
         if  (i >= capitalize_start)
             top [top_len ++] = toupper (a [i ++]);
           else
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
           if  (i >= capitalize_start)
               top [top_len ++] = toupper (a [i ++]);
             else
               top [top_len ++] = a [i ++];
          }
     }
   while  (i < a_len && j < b_len)
     {
      if  (i >= capitalize_start)
          top [top_len ++] = toupper (a [i ++]);
        else
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
                 && tolower (top [i + j]) != tolower (bottom [i + j]))
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



static void  Extract_Needed_Frags
    (GateKeeperStore *store, int32 lo_frag, int32 hi_frag,
     Frag_List_t * list, int * next_olap)

//  Read fragments  lo_frag .. hi_frag  from  store  and save
//  the ids and sequences of those with overlaps to fragments in
//  global  Frag .

  {

#ifdef USE_STREAM_FOR_EXTRACT
   FragStream  *frag_stream;
   int i;
#endif
   static fragRecord  *frag_read = NULL;
   uint32  frag_iid;
   int  bytes_used, total_len, new_total;
   int  extract_ct, stream_ct;
   int  j;

   if  (frag_read == NULL)
       frag_read = new_fragRecord ();

#ifdef USE_STREAM_FOR_EXTRACT
   frag_stream = openFragStream (store, FRAG_S_SEQ);
   resetFragStream (frag_stream, lo_frag, hi_frag);
#endif

   list -> ct = 0;
   total_len = 0;
   extract_ct = stream_ct = 0;

#ifdef USE_STREAM_FOR_EXTRACT
   for  (i = 0;
           nextFragStream (frag_stream, frag_read) && (* next_olap) < Num_Olaps;
           i ++)
#else
   frag_iid = Olap [(* next_olap)] . b_iid;
   while  ((* next_olap) < Num_Olaps
               &&  frag_iid <= hi_frag)
#endif
     {
      char *seqptr = NULL;
      char  seq_buff[AS_READ_MAX_LEN];

      FragType  read_type;
      unsigned  deleted, clear_start, clear_end;
      int  result, shredded;

      stream_ct ++;

#ifdef USE_STREAM_FOR_EXTRACT
      frag_iid = getFragRecordIID(frag_read);
      if  (frag_iid < Olap [(* next_olap)] . b_iid)
          continue;
#else
      getFrag (store, frag_iid, frag_read, FRAG_S_SEQ);
#endif

      deleted = getFragRecordIsDeleted(frag_read);
      if  (deleted)
          goto  Advance_Next_Olap;

      //getReadType_ReadStruct (frag_read, & read_type);
      read_type = AS_READ;
      shredded = (AS_FA_SHREDDED(read_type))? TRUE : FALSE; 

      clear_start = getFragRecordClearRegionBegin(frag_read, AS_READ_CLEAR_OBT);
      clear_end   = getFragRecordClearRegionEnd  (frag_read, AS_READ_CLEAR_OBT);

      // Make sure that we have a valid lowercase sequence string

      seqptr = getFragRecordSequence(frag_read);

      for  (j = clear_start;  j < clear_end;  j ++)
         seq_buff [j] = Filter (seqptr [j]);

      seq_buff [clear_end] = '\0';

      if  (list -> ct >= list -> size)
          {
           list -> size *= 2;
           assert (list -> size > list -> ct);
           list -> entry = (Frag_List_Entry_t *) safe_realloc
                             (list -> entry, list -> size * sizeof (Frag_List_Entry_t));
          }

      list -> entry [list -> ct] . id = frag_iid;
      list -> entry [list -> ct] . shredded = shredded;
      bytes_used = 1 + strlen (seq_buff + clear_start);
      new_total = total_len + bytes_used;
      if  (new_total > list -> buffer_size)
          {
           list -> buffer_size *= 2;
           assert (list -> buffer_size >= new_total);
           list -> buffer = (char *) safe_realloc
                               (list -> buffer, list -> buffer_size);
          }
      list -> entry [list -> ct] . start = total_len;
      strcpy (list -> buffer + total_len, seq_buff + clear_start);
      list -> ct ++;
      total_len = new_total;

      extract_ct ++;

   Advance_Next_Olap:
      while  ((* next_olap) < Num_Olaps
                && Olap [(* next_olap)] . b_iid == frag_iid)
        (* next_olap) ++;
      frag_iid = Olap [(* next_olap)] . b_iid;
     }

#ifdef USE_STREAM_FOR_EXTRACT
   closeFragStream (frag_stream);
#endif

   fprintf (stderr, "Extracted %d of %d fragments in iid range %d .. %d\n",
            extract_ct, stream_ct, lo_frag, hi_frag);

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



static void  Get_Olaps_From_Store
    (char * path, int32 lo_id, int32 hi_id, Olap_Info_t * * olap, int * num)

//  Open overlap store  path  and read from it the overlaps for fragments
//   lo_id .. hi_id , putting them in  (* olap)  for which space
//  is dynamically allocated.  Set  (* num)  to the number of entries
//  in  (* olap) .

  {
   FILE  * fp;
   Short_Olap_Data_t  buff;
   uint32  * offset = NULL;
   uint32  header [3], max_frag, frags_per_file, last_frag_in_file;
   int  file_index, num_frags, olap_size;
   char  filename [MAX_FILENAME_LEN];
   int  first = TRUE;
   int  ct;
   int  i, j;

   sprintf (filename, "%s/offset.olap", path);
   fp = File_Open (filename, "rb");

   Safe_fread (header, sizeof (uint32), 3, fp);
   max_frag = header [0];
   assert (header [1] == sizeof (Short_Olap_Data_t));
   frags_per_file = header [2];
   fprintf (stderr, "sizeof (Short_Olap_Data_t) = %u\n", header [1]);
   fprintf (stderr, "Frags per File = %u\n", frags_per_file);

   assert (1 <= lo_id && lo_id <= hi_id);

   if  (max_frag < lo_id)
       {
        fprintf (stderr, "No overlaps for frags in this range--nothing to do\n");
        (* num) = 0;
        return;
       }
   if  (max_frag < hi_id)
       {
        fprintf (stderr, "Hi frag %d past last ovlStore frag %d\n",
                 hi_id, max_frag);

         hi_id = max_frag;
       }

   num_frags = 2 + hi_id - lo_id;   // go 1 past the end
   offset = (uint32 *) safe_malloc (num_frags * sizeof (uint32));
   CDS_FSEEK (fp, (off_t) (lo_id * sizeof (uint32)), SEEK_CUR);
   Safe_fread (offset, sizeof (uint32), num_frags, fp);

   fclose (fp);

   file_index = (int) ceil ((double) lo_id / frags_per_file);

    // Don't know how many overlaps the last frag in file has
    // without reading the file to see where it ends
    // 1000 is a reasonable upper bound, I hope.
    // Will check below to make sure don't overflow
   last_frag_in_file = (file_index - 1) * frags_per_file;
   olap_size = 0;
   while  (hi_id > last_frag_in_file)
     {
      uint32  first_frag_in_file, lo_offset;

      first_frag_in_file = last_frag_in_file + 1;
      if  (lo_id <= first_frag_in_file)
          lo_offset = 0;
        else
          lo_offset = offset [0];
      last_frag_in_file += frags_per_file;
      if  (max_frag < last_frag_in_file)
          last_frag_in_file = max_frag;
      if  (hi_id < last_frag_in_file
             && last_frag_in_file != first_frag_in_file)
          olap_size += offset [hi_id + 1 - lo_id] - lo_offset;
        else
          olap_size += offset [last_frag_in_file - lo_id] + 1000 - lo_offset;
     }

   (* olap) = (Olap_Info_t*) safe_malloc (olap_size * sizeof (Olap_Info_t));

   ct = 0;
   for  (i = lo_id;  i <= hi_id;  i ++)
     {
      if  (first || i % frags_per_file == 1)
          {
           sprintf (filename, "%s/data%02d.olap", path, file_index);
           fp = File_Open (filename, "rb");
           if  (first)
               CDS_FSEEK (fp, (off_t) (offset [i - lo_id] * sizeof (Short_Olap_Data_t)), SEEK_SET);
           first = FALSE;
          }

      if  (i % frags_per_file == 0)
          while  (fread (& buff, sizeof (Short_Olap_Data_t), 1, fp) == 1)
            {
             if  (ct >= olap_size)
                 {
                  olap_size *= EXPANSION_FACTOR;
                  (* olap) = (Olap_Info_t *) safe_realloc ((* olap),
                                           olap_size * sizeof (Olap_Info_t));
                 }
             (* olap) [ct] . a_iid = i;
             (* olap) [ct] . b_iid = buff . b_iid;
             (* olap) [ct] . a_hang = buff . a_hang;
             (* olap) [ct] . b_hang = buff . b_hang;
             if  (buff . flipped)
                 (* olap) [ct] . orient = INNIE;
               else
                 (* olap) [ct] . orient = NORMAL;
             ct ++;
            }
        else
          for  (j = offset [i - lo_id];  j < offset [i + 1 - lo_id];  j ++)
            {
             Safe_fread (& buff, sizeof (Short_Olap_Data_t), 1, fp);
             if  (ct >= olap_size)
                 {
                  olap_size *= EXPANSION_FACTOR;
                  (* olap) = (Olap_Info_t *) safe_realloc ((* olap),
                                           olap_size * sizeof (Olap_Info_t));
                 }
             (* olap) [ct] . a_iid = i;
             (* olap) [ct] . b_iid = buff . b_iid;
             (* olap) [ct] . a_hang = buff . a_hang;
             (* olap) [ct] . b_hang = buff . b_hang;
             if  (buff . flipped)
                 (* olap) [ct] . orient = INNIE;
               else
                 (* olap) [ct] . orient = NORMAL;
             ct ++;
            }

      if  (i % frags_per_file == 0 || i == hi_id)
          {
           fclose (fp);
           file_index ++;
          }
     }

   (* num) = ct;
   (* olap) = (Olap_Info_t *) safe_realloc ((* olap), ct * sizeof (Olap_Info_t));

   return;
  }



static void  Init_Frag_List
    (Frag_List_t * list)

//  Initilize the entries in fragment list  (* list)

 {
  list -> ct = 0;
  list -> size = 1000;
  list -> entry = (Frag_List_Entry_t *) safe_malloc
                        (Frag_List . size * sizeof (Frag_List_Entry_t));
  list -> buffer_size = Frag_List . size * 550;
  list -> buffer = (char *) safe_malloc (Frag_List . buffer_size);

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

   Frag_List . ct = 0;
   Frag_List . size = 1000;
   Frag_List . entry = (Frag_List_Entry_t *) safe_malloc
                         (Frag_List . size * sizeof (Frag_List_Entry_t));
   Frag_List . buffer_size = Frag_List . size * 550;
   Frag_List . buffer = (char *) safe_malloc (Frag_List . buffer_size);

   return;
  }



static void  Init_Thread_Work_Area
    (Thread_Work_Area_t * wa, int id)

//  Initialize variables in work area  (* wa)  used by thread
//  number  i .

  {
   int  del, offset;
   int  i;

   wa -> thread_id = id;
   wa -> frag_read = new_fragRecord ();
   strcpy (wa -> rev_seq, "acgt");

   offset = 2;
   del = 6;
   for  (i = 0;  i < MAX_ERRORS;  i ++)
     {
      wa -> edit_array [i] = wa -> edit_space + offset;
      offset += del;
      del += 2;
     }

   return;
  }



static Vote_Value_t  Matching_Vote
    (char ch)

//  Return the substitution vote corresponding to  Ch .

  {
   switch  (tolower (ch))
     {
      case  'a' :
        return  A_SUBST;
      case  'c' :
        return  C_SUBST;
      case  'g' :
        return  G_SUBST;
      case  't' :
        return  T_SUBST;
      default :
        return  NO_VOTE;
     }
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



static void  Output_Corrections
    (FILE  * fp)

//  Output the corrections in  Frag  to  fp .

  {
   Correction_Output_t  out;
   double  extension_sum = 0.0;
   int  extension_ct = 0;
   int  i, j;

   for  (i = 0;  i < Num_Frags;  i ++)
     {
      int  clear_extension, last_conf;

      out . frag . is_ID = TRUE;
      out . frag . keep_left = (Frag [i] . left_degree < Degree_Threshold);
      out . frag . keep_right = (Frag [i] . right_degree < Degree_Threshold);
      out . frag . iid = Lo_Frag_IID + i;
      fwrite (& out, sizeof (Correction_Output_t), 1, fp);
      if  (Frag [i] . sequence == NULL)
          continue;   // Deleted fragment

      last_conf = Frag [i] . clear_len - 1;
      if  (Extend_Fragments)
          {
           for  (j = Frag [i] . clear_len;  Frag [i] . sequence [j] != 0;  j ++)
             if  (Frag [i] . vote [j] . confirmed > 0)
                 last_conf = j;
             else if  (j - last_conf > 2 * End_Exclude_Len + 1)
                 break;
           clear_extension = 1 + last_conf - Frag [i] . clear_len;
           extension_sum += clear_extension;
           extension_ct ++;
           out . corr . is_ID = FALSE;
           out . corr . pos = clear_extension;
           out . corr . type = (int) EXTENSION;
           fwrite (& out, sizeof (Correction_Output_t), 1, fp);
          }

      for  (j = 0;  j <= last_conf;  j ++)
        {
         Vote_Value_t  vote, ins_vote;
         int  haplo_ct, ins_haplo_ct;
         int  max, total, tmp;
         int  ins_max, ins_total;
         int  is_change = TRUE;

         if  (Frag [i] . vote [j] . confirmed < 2)
             {
              haplo_ct = 0;
              vote = DELETE;
              total = max = Frag [i] . vote [j] . deletes;
              if  (max >= MIN_HAPLO_OCCURS)
                  haplo_ct ++;

              tmp = Frag [i] . vote [j] . a_subst;
              total += tmp;
              if  (tmp > max)
                  {
                   max = tmp;
                   vote = A_SUBST;
                   is_change = (Frag [i] . sequence [j] != 'a');
                  }
              if  (tmp >= MIN_HAPLO_OCCURS)
                  haplo_ct ++;
              
              tmp = Frag [i] . vote [j] . c_subst;
              total += tmp;
              if  (tmp > max)
                  {
                   max = tmp;
                   vote = C_SUBST;
                   is_change = (Frag [i] . sequence [j] != 'c');
                  }
              if  (tmp >= MIN_HAPLO_OCCURS)
                  haplo_ct ++;
              
              tmp = Frag [i] . vote [j] . g_subst;
              total += tmp;
              if  (tmp > max)
                  {
                   max = tmp;
                   vote = G_SUBST;
                   is_change = (Frag [i] . sequence [j] != 'g');
                  }
              if  (tmp >= MIN_HAPLO_OCCURS)
                  haplo_ct ++;
              
              tmp = Frag [i] . vote [j] . t_subst;
              total += tmp;
              if  (tmp > max)
                  {
                   max = tmp;
                   vote = T_SUBST;
                   is_change = (Frag [i] . sequence [j] != 't');
                  }
              if  (tmp >= MIN_HAPLO_OCCURS)
                  haplo_ct ++;

              if  (2 * max > total
                     && total > 1
                     && is_change
                     && (haplo_ct < 2 || ! Use_Haplo_Ct)
                     && (Frag [i] . vote [j] . confirmed == 0
                           || (Frag [i] . vote [j] . confirmed == 1
                               && max > 6)))
                  {
                   out . corr . is_ID = FALSE;
                   out . corr . pos = j;
                   out . corr . type = (int) vote;
                   fwrite (& out, sizeof (Correction_Output_t), 1, fp);
                  }
             }
         if  (Frag [i] . vote [j] . no_insert < 2)
             {
              ins_haplo_ct = 0;
              ins_vote = A_INSERT;
              ins_total = ins_max = Frag [i] . vote [j] . a_insert;
              if  (ins_max >= MIN_HAPLO_OCCURS)
                  ins_haplo_ct ++;

              tmp = Frag [i] . vote [j] . c_insert;
              ins_total += tmp;
              if  (tmp > ins_max)
                  {
                   ins_max = tmp;
                   ins_vote = C_INSERT;
                  }
              if  (tmp >= MIN_HAPLO_OCCURS)
                  ins_haplo_ct ++;
              
              tmp = Frag [i] . vote [j] . g_insert;
              ins_total += tmp;
              if  (tmp > ins_max)
                  {
                   ins_max = tmp;
                   ins_vote = G_INSERT;
                  }
              if  (tmp >= MIN_HAPLO_OCCURS)
                  ins_haplo_ct ++;
              
              tmp = Frag [i] . vote [j] . t_insert;
              ins_total += tmp;
              if  (tmp > ins_max)
                  {
                   ins_max = tmp;
                   ins_vote = T_INSERT;
                  }
              if  (tmp >= MIN_HAPLO_OCCURS)
                  ins_haplo_ct ++;

              if  (2 * ins_max > ins_total
                     && ins_total > 1
                     && (ins_haplo_ct < 2 || ! Use_Haplo_Ct)
                     && (Frag [i] . vote [j] . no_insert == 0
                           || (Frag [i] . vote [j] . no_insert == 1
                               && ins_max > 6)))
                  {
                   out . corr . is_ID = FALSE;
                   out . corr . pos = j;
                   out . corr . type = (int) ins_vote;
                   fwrite (& out, sizeof (Correction_Output_t), 1, fp);
                  }
             }
        }
     }

   fprintf (stderr, "Fragments processed = %d\n", extension_ct);
   if  (Extend_Fragments)
       fprintf (stderr, "   Avg 3' extension = %.1f bases\n",
                extension_ct == 0 ? 0.0 : extension_sum / extension_ct);

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
             && ((ch = getopt (argc, argv, "d:eF:hk:o:pS:t:v:V:x:")) != EOF))
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

        case  'e' :
          Extend_Fragments = TRUE;
          break;

        case  'F' :
          Olap_Path = optarg;
          break;

        case  'h' :
          errflg = TRUE;
          break;

        case  'k' :
          Kmer_Len = (int) strtol (optarg, & p, 10);
          if  (p == optarg || Kmer_Len <= 1)
              {
               fprintf (stderr, "ERROR:  Illegal k-mer length \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
          break;

        case  'o' :
          Correction_Filename = optarg;
          break;

        case  'p' :
          Use_Haplo_Ct = FALSE;
          break;

        case  'S' :
          Olap_Path = optarg;
          Olaps_From_Store = TRUE;
          break;

        case  't' :
          Num_PThreads = (int) strtol (optarg, & p, 10);
          fprintf (stderr, "Number of pthreads set to %d\n", Num_PThreads);
          break;

        case  'v' :
          Verbose_Level = (int) strtol (optarg, & p, 10);
          fprintf (stderr, "Verbose level set to %d\n", Verbose_Level);
          break;

        case  'V' :
          Vote_Qualify_Len = (int) strtol (optarg, & p, 10);
          fprintf (stderr, "Correction min-match len set to %d\n", Vote_Qualify_Len);
          break;

        case  'x' :
          End_Exclude_Len = (int) strtol (optarg, & p, 10);
          if  (p == optarg || End_Exclude_Len < 0)
              {
               fprintf (stderr, "ERROR:  Illegal end-exclude length \"%s\"\n",
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

   if  (Olap_Path == NULL)
       {
        fprintf (stderr, "ERROR:  Must specify overlaps with -F or -S\n");
        exit (EXIT_FAILURE);
       }

   gkpStore_Path = argv [optind ++];

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

   return;
  }



static int  Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Match_To_End,
     int Delta [MAX_ERRORS], int * Delta_Len, Thread_Work_Area_t * wa)

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
//  (* wa) has storage used by this thread

  {
   double  Score, Max_Score;
   int  Max_Score_Len, Max_Score_Best_d, Max_Score_Best_e;
#if 0
   int Tail_Len;
#endif
   int  Best_d, Best_e, Longest, Row;
   int  Left, Right;
   int  d, e, j, shorter;

//   assert (m <= n);
   Best_d = Best_e = Longest = 0;
   (* Delta_Len) = 0;

   shorter = OVL_Min_int (m, n);
   for  (Row = 0;  Row < shorter && A [Row] == T [Row];  Row ++)
     ;

   wa -> edit_array [0] [0] = Row;

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
      wa -> edit_array [e - 1] [Left] = -2;
      wa -> edit_array [e - 1] [Left - 1] = -2;
      wa -> edit_array [e - 1] [Right] = -2;
      wa -> edit_array [e - 1] [Right + 1] = -2;

      for  (d = Left;  d <= Right;  d ++)
        {
         Row = 1 + wa -> edit_array [e - 1] [d];
         if  ((j = wa -> edit_array [e - 1] [d - 1]) > Row)
             Row = j;
         if  ((j = 1 + wa -> edit_array [e - 1] [d + 1]) > Row)
             Row = j;
         while  (Row < m && Row + d < n
                  && A [Row] == T [Row + d])
           Row ++;

         wa -> edit_array [e] [d] = Row;

         if  (Row == m || Row + d == n)
             {
#if  1
              // Force last error to be mismatch rather than insertion
              if  (Row == m
                     && 1 + wa -> edit_array [e - 1] [d + 1]
                          == wa -> edit_array [e] [d]
                     && d < Right)
                  {
                   d ++;
                   wa -> edit_array [e] [d] = wa -> edit_array [e] [d - 1];
                  }
#endif
              (* A_End) = Row;           // One past last align position
              (* T_End) = Row + d;

              Compute_Delta
                  (Delta, Delta_Len, wa -> edit_array, e, d, Row);

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
                  && wa -> edit_array [e] [Left] < Edit_Match_Limit [e])
        Left ++;
      if  (Left >= 0)
          while  (Left <= Right
                    && wa -> edit_array [e] [Left] + Left < Edit_Match_Limit [e])
            Left ++;
      if  (Left > Right)
          break;
      while  (Right > 0
                  && wa -> edit_array [e] [Right] + Right < Edit_Match_Limit [e])
        Right --;
      if  (Right <= 0)
          while  (wa -> edit_array [e] [Right] < Edit_Match_Limit [e])
            Right --;
      assert (Left <= Right);

      for  (d = Left;  d <= Right;  d ++)
        if  (wa -> edit_array [e] [d] > Longest)
            {
             Best_d = d;
             Best_e = e;
             Longest = wa -> edit_array [e] [d];
            }
#if  1
      Score = Longest * BRANCH_PT_MATCH_VALUE - e;
               // Assumes  BRANCH_PT_MATCH_VALUE - BRANCH_PT_ERROR_VALUE == 1.0
      if  (Score > Max_Score
               && Best_e <= Error_Bound [OVL_Min_int (Longest, Longest + Best_d)])
          {
           Max_Score = Score;
           Max_Score_Len = Longest;
           Max_Score_Best_d = Best_d;
           Max_Score_Best_e = Best_e;
          }
#endif
     }

   Compute_Delta
       (Delta, Delta_Len, wa -> edit_array, Max_Score_Best_e,
        Max_Score_Best_d, Max_Score_Len);

   (* A_End) = Max_Score_Len;
   (* T_End) = Max_Score_Len + Max_Score_Best_d;
   (* Match_To_End) = FALSE;

   return  Max_Score_Best_e;
  }



static void  Process_Olap
    (Olap_Info_t * olap, char * b_seq, char * rev_seq, int * rev_id,
     int shredded, Thread_Work_Area_t * wa)

//  Find the alignment referred to in  olap , where the  a_iid
//  fragment is in  Frag  and the  b_iid  sequence is in  b_seq .
//  Use the alignment to increment the appropriate vote fields
//  for the a fragment.   shredded  is true iff the b fragment
//  is from shredded data, in which case the overlap will be
//  ignored if the a fragment is also shredded.
//  rev_seq  is a buffer to hold the reverse complement of  b_seq
//  if needed.  (* rev_id) is used to keep track of whether
//  rev_seq  has been created yet.  (* wa) is the work-area
//  containing space for the process to use in case of multi-threading.

  {
   char  * a_part, * b_part;
   int  a_part_len, b_part_len, a_end, b_end, olap_len;
   int  match_to_end, delta [MAX_ERRORS], delta_len, errors;
   int  a_offset, sub;

   if  (Verbose_Level > 0)
       printf ("Process_Olap:  %8d %8d %5d %5d  %c\n",
               olap -> a_iid, olap -> b_iid,
               olap -> a_hang, olap -> b_hang,
               olap -> orient == INNIE ? 'I' : 'N');

   sub = olap -> a_iid - Lo_Frag_IID;

   if  (shredded && Frag [sub] . shredded)
       return;

   a_part = Frag [sub] . sequence;
   a_offset = 0;

   if  (olap -> a_hang <= 0)
       {
        if  (Frag [sub] . left_degree < MAX_DEGREE)
            Frag [sub] . left_degree ++;
       }
   if  (olap -> b_hang >= 0)
       {
        if  (Frag [sub] . right_degree < MAX_DEGREE)
            Frag [sub] . right_degree ++;
       }

   if  (olap -> a_hang > 0)
       {
        a_part += olap -> a_hang;
        a_offset += olap -> a_hang;
       }

   if  (olap -> orient == NORMAL)
       b_part = b_seq;
     else
       {
        if  ((* rev_id) != olap -> b_iid)
            {
             strcpy (rev_seq, b_seq);
             Rev_Complement (rev_seq);
             (* rev_id) = olap -> b_iid;
            }
        b_part = rev_seq;
       }
   if  (olap -> a_hang < 0)
       b_part -= olap -> a_hang;

   if  (Verbose_Level > 0)
       printf ("b_part = %p  is ascii %d  rev_seq is %d\n",
               b_part, (int) (* b_part), (int) (* rev_seq));
   if  (! isalpha (* b_part) || ! isalpha (* rev_seq))
       exit (-1);

   if  (Verbose_Level > 0)
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
           if  (j + a_offset >= Frag [sub] . clear_len )
               putchar (toupper (a_part [j]));
             else
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
               & match_to_end, delta, & delta_len, wa);
if  (a_end < 0 || a_end > a_part_len || b_end < 0 || b_end > b_part_len)
    {
     fprintf (stderr, "ERROR:  Bad edit distance\n");
     fprintf (stderr, "errors = %d  a_end = %d  b_end = %d\n",
              errors, a_end, b_end);
     fprintf (stderr, "a_part_len = %d  b_part_len = %d\n",
              a_part_len, b_part_len);
     fprintf (stderr, "a_iid = %d  b_iid = %d  match_to_end = %c\n",
              olap -> a_iid, olap -> b_iid, match_to_end ? 'T' : 'F');
     exit (-3);
    }


   if  (Verbose_Level > 0)
       {
        printf ("  errors = %d  delta_len = %d\n", errors, delta_len);
        printf ("  a_align = %d/%d  b_align = %d/%d\n",
                a_end, a_part_len, b_end, b_part_len);
        Display_Alignment
            (a_part, a_end, b_part, b_end, delta, delta_len,
             Frag [sub] . clear_len - a_offset);
       }

   if  (! match_to_end && a_end + a_offset >= Frag [sub] . clear_len - 1)
       {
        olap_len = OVL_Min_int (a_end, b_end);
        match_to_end = TRUE;
       }
       
   if  (errors <= Error_Bound [olap_len] && match_to_end)
       {
        Analyze_Alignment (delta, delta_len, a_part, b_part,
                           a_end, b_end, a_offset, sub);
       }
     else
       Failed_Olaps ++;
                      
   return;
  }



static void  Read_Frags
    (void)

//  Open and read fragments with IIDs from  Lo_Frag_IID  to
//  Hi_Frag_IID  from  gkpStore_Path  and store them in
//  global  Frag .

  {
   char  seq_buff [AS_READ_MAX_LEN + 1];
   char  qual_buff [AS_READ_MAX_LEN + 1];
   fragRecord  *frag_read;
   unsigned  clear_start, clear_end;
   int32  high_store_frag;
   int  i, j;

   high_store_frag = getLastElemFragStore (gkpStore);
   if  (Hi_Frag_IID == INT_MAX)
       Hi_Frag_IID = high_store_frag;
   if  (Hi_Frag_IID > high_store_frag)
       {
        fprintf (stderr, "ERROR:  Hi frag %d is past last store frag %d\n",
                 Hi_Frag_IID, high_store_frag);
        exit (EXIT_FAILURE);
       }

   Num_Frags = 1 + Hi_Frag_IID - Lo_Frag_IID;
   Frag = (Frag_Info_t *) safe_calloc (Num_Frags, sizeof (Frag_Info_t));

   frag_read = new_fragRecord ();

#ifdef USE_STORE_DIRECTLY_READ
  Internal_gkpStore = openGateKeeperStore (gkpStore_Path, FALSE);
  assert (Internal_gkpStore != NULL);
#else
   Internal_gkpStore
       = loadFragStorePartial (gkpStore_Path,
                               Lo_Frag_IID, Hi_Frag_IID, FRAG_S_SEQ);
#endif
   
   Frag_Stream = openFragStream (Internal_gkpStore, FRAG_S_SEQ);
   resetFragStream (Frag_Stream, Lo_Frag_IID, Hi_Frag_IID);

   for  (i = 0;  nextFragStream (Frag_Stream, frag_read);
           i ++)
     {
      FragType  read_type;
      unsigned  deleted;
      int  result, frag_len;

      deleted = getFragRecordIsDeleted(frag_read);
      if  (deleted)
          {
           Frag [i] . sequence = NULL;
           Frag [i] . vote = NULL;
           continue;
          }

      //getReadType_ReadStruct (frag_read, & read_type);
      read_type = AS_READ;
      Frag [i] . shredded = (AS_FA_SHREDDED(read_type))? TRUE : FALSE; 

      strcpy(seq_buff, getFragRecordSequence(frag_read));

      clear_start = getFragRecordClearRegionBegin(frag_read, AS_READ_CLEAR_OBT);
      clear_end   = getFragRecordClearRegionEnd  (frag_read, AS_READ_CLEAR_OBT);

      // Make sure that we have a legal lowercase sequence string

      Frag [i] . clear_len = clear_end - clear_start;
      if  (Extend_Fragments)
          frag_len = strlen (seq_buff);
        else
          {
           frag_len = clear_end;
           seq_buff [clear_end] = '\0';
          }
      for  (j = clear_start;  j < frag_len;  j ++)
         seq_buff [j] = Filter (seq_buff [j]);

      Frag [i] . sequence = strdup (seq_buff + clear_start);
      Frag [i] . vote = (Vote_Tally_t *) safe_calloc (frag_len - clear_start,
                                                      sizeof (Vote_Tally_t));
      Frag [i] . left_degree = Frag [i] . right_degree = 0;
     }

   del_fragRecord (frag_read);
   closeFragStream (Frag_Stream);
   closeGateKeeperStore (Internal_gkpStore);

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
        olap_size = 1000;
        Olap = (Olap_Info_t*) safe_malloc (olap_size * sizeof (Olap_Info_t));

        fp = File_Open (Olap_Path, "r");

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
                     Olap = (Olap_Info_t *) safe_realloc (Olap,
                                olap_size * sizeof (Olap_Info_t));
                    }
                Olap [ct] . a_iid = a_iid;
                Olap [ct] . b_iid = b_iid;
                if  (orient [0] == 'O')
                    {
                     Olap [ct] . a_hang = - b_hang;
                     Olap [ct] . b_hang = - a_hang;
                     Olap [ct] . orient = INNIE;
                    }
                  else
                    {
                     Olap [ct] . a_hang = a_hang;
                     Olap [ct] . b_hang = b_hang;
                     Olap [ct] . orient = NORMAL;
                    }
                ct ++;
               }

           if  (a_iid > Hi_Frag_IID)   // Speed up if file is sorted
               break;
          }

        Num_Olaps = ct;
        fclose (fp);

        if  (Num_Olaps == 0)
            {
             fprintf (stderr, "No overlaps read, nothing to do\n");
             exit (EXIT_FAILURE);
            }

        Olap = (Olap_Info_t *) safe_realloc (Olap, Num_Olaps * sizeof (Olap_Info_t));
       }

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



static void  Stream_Old_Frags
    (void)

//  Read old fragments in  gkpStore  and choose the ones that
//  have overlaps with fragments in  Frag .  Recompute the
//  overlaps and record the vote information about changes to
//  make (or not) to fragments in  Frag .

  {
   Thread_Work_Area_t  wa;
   fragRecord *frag_read;
   char  seq_buff [AS_READ_MAX_LEN + 1];
   char *seqptr;
   char  rev_seq [AS_READ_MAX_LEN + 1] = "acgt";
   unsigned  clear_start, clear_end;
   int32  lo_frag, hi_frag;
   int  next_olap;
   int  i, j;

   Init_Thread_Work_Area (& wa, 0);
   frag_read = new_fragRecord ();

   Frag_Stream = openFragStream (gkpStore, FRAG_S_SEQ);

   lo_frag = Olap [0] . b_iid;
   hi_frag = Olap [Num_Olaps - 1] . b_iid;

   resetFragStream (Frag_Stream, lo_frag, hi_frag);
   
   next_olap = 0;
   for  (i = 0;  nextFragStream (Frag_Stream, frag_read)
                   && next_olap < Num_Olaps;
           i ++)
     {
      FragType  read_type;
      int32  rev_id;
      uint32  frag_iid;
      unsigned  deleted;
      int  result, shredded;

      frag_iid = getFragRecordIID(frag_read);
      if  (frag_iid < Olap [next_olap] . b_iid)
          continue;

      deleted = getFragRecordIsDeleted (frag_read);
      if  (deleted)
          continue;

      //getReadType_ReadStruct (frag_read, & read_type);
      read_type = AS_READ;
      shredded = (AS_FA_SHREDDED(read_type))? TRUE : FALSE; 

      clear_start = getFragRecordClearRegionBegin(frag_read, AS_READ_CLEAR_OBT);
      clear_end   = getFragRecordClearRegionEnd  (frag_read, AS_READ_CLEAR_OBT);

      seqptr = getFragRecordSequence(frag_read);

      // Make sure that we have a valid lowercase sequence string

      for  (j = clear_start;  j < clear_end;  j ++)
         seq_buff [j] = Filter (seqptr [j]);

      seq_buff [clear_end] = '\0';

      rev_id = -1;
      while  (next_olap < Num_Olaps
                && Olap [next_olap] . b_iid == frag_iid)
        {
         Process_Olap (Olap + next_olap, seq_buff + clear_start,
                       rev_seq, & rev_id, shredded, & wa);
         next_olap ++;
        }
     }

   del_fragRecord (frag_read);
   closeFragStream (Frag_Stream);

   return;
  }



void *  Threaded_Process_Stream
    (void * ptr)

//  Process all old fragments in  Internal_gkpStore .  Only
//  do overlaps/corrections with fragments where
//    frag_iid % Num_PThreads == thread_id

  {
   Thread_Work_Area_t  * wa = (Thread_Work_Area_t *) ptr;
   int  olap_ct;
   int  i;

   olap_ct = 0;

   for  (i = 0;  i < wa -> frag_list -> ct;  i ++)
     {
      int32  skip_id = -1;

      while  (wa -> frag_list -> entry [i] . id > Olap [wa -> next_olap] . b_iid)
        {
         if  (Olap [wa -> next_olap] . b_iid != skip_id)
             {
              fprintf (stderr, "SKIP:  b_iid = %d\n",
                       Olap [wa -> next_olap] . b_iid);
              skip_id = Olap [wa -> next_olap] . b_iid;
             }
         wa -> next_olap ++;
        }
      if  (wa -> frag_list -> entry [i] . id != Olap [wa -> next_olap] . b_iid)
          {
           fprintf (stderr, "ERROR:  Lists don't match\n");
           fprintf (stderr, "frag_list iid = %d  next_olap = %d  i = %d\n",
                    wa -> frag_list -> entry [i] . id,
                    Olap [wa -> next_olap] . b_iid, i);
           exit (EXIT_FAILURE);
          }

      wa -> rev_id = -1;
      while  (wa -> next_olap < Num_Olaps
                && Olap [wa -> next_olap] . b_iid == wa -> frag_list -> entry [i] . id)
        {
         if  (Olap [wa -> next_olap] . a_iid % Num_PThreads == wa -> thread_id)
             {
              Process_Olap
                  (Olap + wa -> next_olap,
                   wa -> frag_list -> buffer + wa -> frag_list -> entry [i] . start,
                   wa -> rev_seq, & (wa -> rev_id),
                   wa -> frag_list -> entry [i] . shredded, wa);
              olap_ct ++;
             }
         wa -> next_olap ++;
        }
     }

pthread_mutex_lock (& Print_Mutex);
Now = time (NULL);
fprintf (stderr, "Thread %d processed %d olaps at %s",
         wa -> thread_id, olap_ct, ctime (& Now));
pthread_mutex_unlock (& Print_Mutex);

//   if  (wa -> thread_id > 0)
       pthread_exit (ptr);

   return  ptr;
  }



static void  Threaded_Stream_Old_Frags
    (void)

//  Read old fragments in  gkpStore  that have overlaps with
//  fragments in  Frag .  Read a batch at a time and process them
//  with multiple pthreads.  Each thread processes all the old fragments
//  but only changes entries in  Frag  that correspond to its thread
//  ID.  Recomputes the overlaps and records the vote information about
//  changes to make (or not) to fragments in  Frag .

  {
   pthread_attr_t  attr;
   pthread_t  * thread_id;
   Frag_List_t  frag_list_1, frag_list_2;
   Frag_List_t  * curr_frag_list, * next_frag_list, * save_frag_list;
   Thread_Work_Area_t  * thread_wa;
   int  next_olap, save_olap, status;
   int32  first_frag, last_frag, lo_frag, hi_frag;
   int  i;

   fprintf (stderr, "### Using %d pthreads (new version)\n", Num_PThreads);

   pthread_mutex_init (& Print_Mutex, NULL);
   pthread_attr_init (& attr);
   pthread_attr_setstacksize (& attr, THREAD_STACKSIZE);
   thread_id = (pthread_t *) safe_calloc
                   (Num_PThreads, sizeof (pthread_t));
   thread_wa = (Thread_Work_Area_t *) safe_malloc
                   (Num_PThreads * sizeof (Thread_Work_Area_t));

   for  (i = 0;  i < Num_PThreads;  i ++)
     Init_Thread_Work_Area (thread_wa + i, i);
   Init_Frag_List (& frag_list_1);
   Init_Frag_List (& frag_list_2);

   first_frag = Olap [0] . b_iid;
   last_frag = Olap [Num_Olaps - 1] . b_iid;
   assert (first_frag <= last_frag);

   lo_frag = first_frag;
   hi_frag = lo_frag + FRAGS_PER_BATCH - 1;
   if  (hi_frag > last_frag)
       hi_frag = last_frag;
   next_olap = 0;

#ifdef USE_STORE_DIRECTLY_STREAM
  Internal_gkpStore = openGateKeeperStore (gkpStore_Path, FALSE);
  assert (Internal_gkpStore != NULL);
#else
   Internal_gkpStore
       = loadFragStorePartial (gkpStore_Path, lo_frag, hi_frag, FRAG_S_SEQ);
#endif

   curr_frag_list = & frag_list_1;
   next_frag_list = & frag_list_2;
   save_olap = next_olap;

   Extract_Needed_Frags (Internal_gkpStore, lo_frag, hi_frag,
                         curr_frag_list, & next_olap);

#ifndef USE_STORE_DIRECTLY_STREAM
   closeGateKeeperStore (Internal_gkpStore);
#endif

   while  (lo_frag <= last_frag)
     {
      // Process fragments in  curr_frag_list  in background
      for  (i = 0;  i < Num_PThreads;  i ++)
        {
         thread_wa [i] . lo_frag = lo_frag;
         thread_wa [i] . hi_frag = hi_frag;
         thread_wa [i] . next_olap = save_olap;
         thread_wa [i] . frag_list = curr_frag_list;
         status = pthread_create
                      (thread_id + i, & attr, Threaded_Process_Stream,
                       thread_wa + i);
         if  (status != 0)
             {
              fprintf (stderr, "pthread_create error at line %d:  %s\n",
                       __LINE__, strerror (status));
              exit (-3);
             }
        }

      // Read next batch of fragments
      lo_frag = hi_frag + 1;
      if  (lo_frag <= last_frag)
          {
           hi_frag = lo_frag + FRAGS_PER_BATCH - 1;
           if  (hi_frag > last_frag)
               hi_frag = last_frag;

#ifndef USE_STORE_DIRECTLY_STREAM
           Internal_gkpStore
               = loadFragStorePartial (gkpStore_Path, lo_frag, hi_frag, FRAG_S_SEQ);
#endif

           save_olap = next_olap;

           Extract_Needed_Frags (Internal_gkpStore, lo_frag, hi_frag,
                                 next_frag_list, & next_olap);

#ifndef USE_STORE_DIRECTLY_STREAM
           closeGateKeeperStore (Internal_gkpStore);
#endif
          }

      // Wait for background processing to finish
      for  (i = 0;  i < Num_PThreads;  i ++)
        {
         void  * ptr;

         status = pthread_join (thread_id [i], & ptr);
         if  (status != 0)
             {
              fprintf (stderr, "pthread_join error at line %d:  %s\n",
                       __LINE__, strerror (status));
              exit (-3);
             }
        }

      save_frag_list = curr_frag_list;
      curr_frag_list = next_frag_list;
      next_frag_list = save_frag_list;
     }
   
#ifdef USE_STORE_DIRECTLY_STREAM
   closeGateKeeperStore (Internal_gkpStore);
#endif

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s [-ehp] [-d DegrThresh] [-k KmerLen] [-x ExcludeLen]\n"
           "           [-F OlapFile] [-S OlapStore] [-o CorrectFile]\n"
           "           [-t NumPThreads] [-v VerboseLevel]\n"
           "           [-V Vote_Qualify_Len]\n"
           "           <FragStore> <lo> <hi>\n"
           "\n"
           "Makes corrections to fragment sequence based on overlaps\n"
           "and recomputes overlaps on corrected fragments\n"
           "Fragments come from <FragStore> <lo> and <hi> specify\n"
           "the range of fragments to modify\n"
           "\n"
           "Options:\n"
           "-d   set keep flag on end of frags with less than this many olaps\n"
           "-F   specify file of sorted overlaps to use (in the format produced\n"
           "     by  get-olaps\n"
           "-h   print this message\n"
           "-k   minimum exact-match region to prevent change\n"
           "-o   specify output file to hold correction info\n"
           "-p   don't use haplotype counts to correct\n"
           "-S   specify the binary overlap store containing overlaps to use\n"
           "-t   set number of p-threads to use\n"
           "-v   specify level of verbose outputs, higher is more\n"
           "-V   specify number of exact match bases around an error to vote to change\n"
           "-x   length of end of exact match to exclude in preventing change\n",
           command);

   return;
  }



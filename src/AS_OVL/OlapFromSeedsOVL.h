
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
* Module:  OlapFromSeeds.h
* Description:
*   Declarations for olap-from-seeds program
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
 * $Id: OlapFromSeedsOVL.h,v 1.1 2007-05-01 22:53:09 adelcher Exp $
 * $Revision: 1.1 $
*/


#ifndef  __OLAPFROMSEEDS_H_INCLUDED
#define  __OLAPFROMSEEDS_H_INCLUDED

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
#include  "AS_OVS_overlapStore.h"
#include  "SharedOVL.h"

#if 0
#include  "AS_PER_ReadStruct.h"
#include  "AS_PER_genericStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_UTL_version.h"
#include  "OlapStoreOVL.h"
#endif


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
#define  DEFAULT_MIN_OLAP_LEN        40
  //  Default value for  Min_Olap_Len
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
#define  FRAG_LEN_BITS               15
  //  Number of bits to store lengths and positions on fragments
#define  FRAGS_PER_BATCH             100000
  //  Number of old fragments to read into memory-based fragment
  //  store at a time for processing
#define  MAX_ERROR_RATE              AS_GUIDE_ERROR_RATE
  //  The largest error allowed in overlaps
#define  MAX_FRAG_LEN                2048
  //  The longest fragment allowed
#define  MAX_DEGREE                  32767
  //  Highest number of votes before overflow
#define  MAX_LINE                    1000
  //  Longest possible input line
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

#define USE_STORE_DIRECTLY_STREAM
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
   unsigned  clear_len : FRAG_LEN_BITS;
   unsigned  trim_5p : FRAG_LEN_BITS;
   unsigned  trim_3p : FRAG_LEN_BITS;
   unsigned  left_degree : FRAG_LEN_BITS;
   unsigned  right_degree : FRAG_LEN_BITS;
   unsigned  shredded : 1;    // True if shredded read
   unsigned  unused : 1;
  }  Frag_Info_t;

const int  INNIE = 0;
const int  NORMAL = 1;

typedef  struct                 
  {
   int32  a_iid, b_iid;
   signed int  a_hang : FRAG_LEN_BITS;
   signed int  b_hang : FRAG_LEN_BITS;
   signed int  orient : 2;
  }  Olap_Info_t;

typedef  struct
  {
   unsigned  id : 31;
   unsigned  shredded : 1;
   unsigned  trim_5p : FRAG_LEN_BITS;
   unsigned  trim_3p : FRAG_LEN_BITS;
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
   int  failed_olaps;
   FragStream  * frag_stream;
   fragRecord  * frag_read;
   Frag_List_t  * frag_list;
   char  rev_seq [AS_READ_MAX_LEN + 1];
   int  rev_id;
   int  * edit_array [MAX_ERRORS];
   int  edit_space [(MAX_ERRORS + 4) * MAX_ERRORS];
  }  Thread_Work_Area_t;

typedef enum
  {TEXT_FILE, BINARY_FILE, OVL_STORE}  OVL_Output_t;



//  Static Globals

static BinaryOverlapFile  * Binary_OVL_Output_fp = NULL;
  // Pointer for binary overlap outputs
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
static int  Min_Olap_Len = DEFAULT_MIN_OLAP_LEN;
  // The minimum number of bp in each read to report an overlap
static time_t  Now;
  // Used to get current time
static int  Num_Frags;
  // Number of fragments being corrected
static int  Num_Olaps;
  // Number of overlaps being used
static int  Num_PThreads = DEFAULT_NUM_PTHREADS;
  // Number of pthreads to process overlaps/corrections;
static int  Offsets_WRT_Raw = TRUE;
  // Indicates if offset information of seeds is relative to the
  // entire raw read, or just the clear range
static Olap_Info_t  * Olap = NULL;
  // Array of overlaps being used
static FILE  * OVL_Output_fp = NULL;
  // File where output overlaps should be sent
static char  * OVL_Output_Path = NULL;
  // Name of file or store to which overlaps should be sent
static OVL_Output_t  OVL_Output_Type = TEXT_FILE;
  // Type of output for overlaps
static int  Seeds_From_Store = FALSE;
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
static void  Get_Seeds_From_Store
  (char * path, int32 lo_id, int32 hi_id, Olap_Info_t * * olap, int * num);
static void  Init_Frag_List
  (Frag_List_t * list);
static void  Initialize_Globals
  (void);
static void  Init_Thread_Work_Area
  (Thread_Work_Area_t * wa, int id);
static Vote_Value_t  Matching_Vote
  (char ch);
static void  Output_Corrections
  (FILE * fp);
static void  Output_Olap
  (FILE * fp, Olap_Info_t * olap, int a_lo, int a_hi, int a_len,
   int b_lo, int b_hi, int b_len, int errors);
static void  Parse_Command_Line
  (int argc, char * argv []);
static void  Process_Olap
  (Olap_Info_t * olap, char * b_seq, unsigned b_len, char * rev_seq, int * rev_id,
   int shredded, Thread_Work_Area_t * wa);
static void  Read_Frags
  (void);
static void  Read_Seeds
  (void);
static void  Rev_Complement
  (char * s);
static void  Show_Edit_Array
  (int * * ea, int errs);
static void  Stream_Old_Frags
  (void);
void *  Threaded_Process_Stream
  (void * ptr);
static void  Threaded_Stream_Old_Frags
  (void);
static void  Tidy_Up
  (void);
static void  Usage
  (char * command);


#endif

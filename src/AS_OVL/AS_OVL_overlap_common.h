
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
* Module:  AS_OVL_overlap.c
* Description:
*   Find matches between pairs of DNA strings and output overlaps
*   (matches that go to an end of both strings involved) and branch points
*   (matches that do not go to the end of either string).
*   Use hashing to select candidates and simple edit distance
*   to confirm matches.
* 
*    Programmer:  A. Delcher
*       Written:  29 May 98
*  Last Revised:  30 Nov 99
* 
* 
* Assumptions:
*  argv[1] is <filename>[.<ext>]
*  Input is from file  argv [1] 
*  Output is to file  <filename> . OUTPUT_FILENAME_EXTENSION
* 
* Notes:
*   - Removed  ceil  function everywhere.  Overlaps are all <= ERROR_RATE
*     except where indels change length of olap region.
*   - Use hash table for index of strings that match
*   - Do memory allocation by self.
*   - "Band" edit distance calculations
*   - Avoid "hopeless" edit distance calculations
*   - Find branch points in overlaps.
*
*************************************************/

/* RCS info
 * $Id: AS_OVL_overlap_common.h,v 1.3 2005-03-22 19:06:43 jason_miller Exp $
 * $Revision: 1.3 $
*/


/*************************************************************************/
/* System include files */
/*************************************************************************/

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <fcntl.h>
#include  <sys/types.h>
#include  <string.h>
#include  <dirent.h>
#include  <sys/stat.h>
#include  <unistd.h>
#include  <float.h>

/*************************************************************************/
/* Local include files */
/*************************************************************************/

#include  "AS_OVL_delcher.h"
#include  "AS_PER_ReadStruct.h"
#include  "AS_PER_genericStore.h"
#include  "AS_PER_fragStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_OVL_overlap.h"
#include  "AS_UTL_version.h"



/*************************************************************************/
/* Type definitions */
/*************************************************************************/

typedef  uint32  Check_Vector_t;
    // Bit vector to see if hash bucket could possibly contain a match

typedef  enum Overlap_Type
  {NONE, LEFT_BRANCH_PT, RIGHT_BRANCH_PT, DOVETAIL}  Overlap_t;

typedef  struct Olap_Info
  {
   int  s_lo, s_hi, t_lo, t_hi;
   double  quality;
   int  delta [MAX_ERRORS];
   int  delta_ct;
   int  s_left_boundary, s_right_boundary;
   int  t_left_boundary, t_right_boundary;
   int  min_diag, max_diag;
  }  Olap_Info_t;

typedef  struct String_Ref
  {
   unsigned int  String_Num : STRING_NUM_BITS;
   unsigned int  Offset : OFFSET_BITS;
   unsigned int  Empty : 1;
   unsigned int  Last : 1;
  }  String_Ref_t;

typedef  struct Hash_Bucket
  {
   String_Ref_t  Entry [ENTRIES_PER_BUCKET];
   unsigned char  Check [ENTRIES_PER_BUCKET];
#if  ANALYZE_HITS
   unsigned  Hits [ENTRIES_PER_BUCKET];
#else
   unsigned char  Hits [ENTRIES_PER_BUCKET];
#endif
   int16  Entry_Ct;
  }  Hash_Bucket_t;

typedef  struct Hash_Frag_Info
  {
   unsigned int  length : 30;
   unsigned int  left_end_screened : 1;
   unsigned int  right_end_screened : 1;
  }  Hash_Frag_Info_t;

#if  DO_OLAP_ALIGN_PROFILE
typedef  int  Align_Entry_t [6];
    //  Holds counts of what matched at a base position in an
    //  overlap alignment.  Entries  [0 .. 3]  represent
    //  a, c, g, t, resp.;  entry  [4]  is for insertions;
    //  entry  [5]  is for deletions.
typedef  Align_Entry_t  * Align_Ct_t;
#endif


/*************************************************************************/
/* Static Globals */
/*************************************************************************/

#if  SHOW_STATS
FILE  * Stat_File = NULL;

#include  "AS_OVL_olapstats.h"

static Distrib_t  Olap_Len_Dist, Num_Olaps_Dist, Kmer_Freq_Dist, Kmer_Hits_Dist;
static Distrib_t  String_Olap_Dist, Hits_Per_Olap_Dist, Diag_Dist, Gap_Dist;
static Distrib_t  Exacts_Per_Olap_Dist, Edit_Depth_Dist, Hash_Hits_Dist;
#endif

#if  LIST_TOTALLY_SCREENED
static FILE  * Total_Screen_File = NULL;
#endif

static int64  Bad_Short_Window_Ct = 0;
    //  The number of overlaps rejected because of too many errors in
    //  a small window
static int64  Bad_Long_Window_Ct = 0;
    //  The number of overlaps rejected because of too many errors in
    //  a long window
static double  Branch_Match_Value = DEFAULT_BRANCH_MATCH_VAL;
static double  Branch_Error_Value = DEFAULT_BRANCH_MATCH_VAL - 1.0;
    //  Scores of matches and mismatches in alignments.  Alignment
    //  ends at maximum score.
static char  * Data = NULL;
    //  Stores sequence data of fragments in hash table
static char  * Quality_Data = NULL;
    //  Stores quality data of fragments in hash table
static size_t  Data_Len = 0;
static int  Doing_Partial_Overlaps = FALSE;
    //  If set true by the G option (G for Granger)
    //  then allow overlaps that do not extend to the end
    //  of either read.
static int64  Extra_Ref_Ct = 0;
static String_Ref_t  * Extra_Ref_Space = NULL;
static int  Frag_Olap_Limit = FRAG_OLAP_LIMIT;
    //  Maximum number of overlaps for end of an old fragment against
    //  a single hash table of frags, in each orientation
static double  Genome_Len = -1.0;
    //  Estimated genome length read by  -K  option.  Used to compute
    //  k-mer hit limit.
static Check_Vector_t  * Hash_Check_Array = NULL;
    //  Bit vector to eliminate impossible hash matches
static int  Hash_String_Num_Offset = 1;
static Hash_Bucket_t  * Hash_Table;
static int  Hi_Hit_Limit = DEFAULT_HI_HIT_LIMIT;
    //  Any kmer in the hash table with this many or more references
    //  is not used to initiate overlaps.  Set by  -K  option.
static int  Ignore_Screen_Info = FALSE;
    //  If true will ignore screen messages with fragments
static double  Kmer_Freq_Bound = -1.0;
    //  Target number of kmers in genome to compute  Hi_Hit_Limit
    //  from.  Gotten from  -K  option.
static int64  Kmer_Hits_With_Olap_Ct = 0;
static int64  Kmer_Hits_Without_Olap_Ct = 0;
static double  Kmer_Prob_Limit = -1.0;
    //  Probability threshold for setting  Hi_Hit_Limit .
    //  Gotten from  -K  option.
static uint64  * Loc_ID = NULL;
    //  Locale ID field of each frag in hash table if in  Contig_Mode .
static int64  Multi_Overlap_Ct = 0;
static String_Ref_t  * Next_Ref = NULL;
static int  Single_Line_Output = FALSE;
    //  If set true by -q option, output a single ASCII line
    //  for each overlap in the same format as used by
    //  partial overlaps
static int  String_Ct;
    //  Number of fragments in the hash table
static Hash_Frag_Info_t  * String_Info = NULL;
static int64  * String_Start = NULL;
static int  Unique_Olap_Per_Pair = TRUE;
    //  If true will allow at most
    //  one overlap output message per oriented fragment pair
    //  Set true by  -u  command-line option; set false by  -m
static int  Use_Window_Filter = FALSE;

static int  Read_Edit_Match_Limit [MAX_ERRORS] = {0};
    //  This array [e] is the minimum value of  Edit_Array [e] [d]
    //  to be worth pursuing in edit-distance computations between reads
static int  Guide_Edit_Match_Limit [MAX_ERRORS] = {0};
    //  This array [e] is the minimum value of  Edit_Array [e] [d]
    //  to be worth pursuing in edit-distance computations between guides
static int  Read_Error_Bound [AS_READ_MAX_LEN + 1];
    //  This array [i]  is the maximum number of errors allowed
    //  in a match between reads of length  i , which is
    //  i * AS_READ_ERROR_RATE .
static int  Guide_Error_Bound [AS_READ_MAX_LEN + 1];
    //  This array [i]  is the maximum number of errors allowed
    //  in a match between guides of length  i , which is
    //  i * AS_GUIDE_ERROR_RATE .
static double  Branch_Cost [AS_READ_MAX_LEN + 1];
    //  Branch_Cost [i]  is the "goodness" of matching i characters
    //  after a single error in determining branch points.
static int  Bit_Equivalent [256] = {0};
    //  Table to convert characters to 2-bit integer code
static int  Char_Is_Bad [256] = {0};
    //  Table to check if character is not a, c, g or t.
static FragType  * Kind_Of_Frag = NULL;
    //  Type of fragment in hash table, read or guide

static int64  Hash_Entries = 0;
#if  SHOW_STATS
static int64  Collision_Ct = 0, Match_Ct = 0;
static int64  Hash_Find_Ct = 0, Edit_Dist_Ct = 0;
static int32  Overlap_Ct;
#endif

static int64  Total_Overlaps = 0;
static int64  Contained_Overlap_Ct = 0;
static int64  Dovetail_Overlap_Ct = 0;

#if  USE_SOURCE_FIELD
#define  MAX_SOURCE_LENGTH    2048
static char  Global_Left_Annotation, Global_Right_Annotation;
static char  Global_Other_Annotation;
static char  * Left_Repeat_Tag = NULL;
static char  * Right_Repeat_Tag = NULL;
static char  * Other_Repeat_Tag = NULL;
FILE  * Source_Log_File = NULL;
#endif

#if  SHOW_SNPS
static int  Global_SNP_A_Id, Global_SNP_B_Id;
static int  Global_Match_SNP_Bin [20] = {0};
static int  Global_Indel_SNP_Bin [20] = {0};
static int  Global_All_Match_Bin [20] = {0};
static int  Global_Hi_Qual_Bin [20] = {0};
static int  Global_Olap_Ct = 0;
static int  Global_Unscreened_Ct = 0;
#endif

static int  Global_Debug_Flag = FALSE;

#if  DO_OLAP_ALIGN_PROFILE
Align_Ct_t  * Align_Ct = NULL;
int  Align_Ct_Size = 0;
Align_Entry_t  * Align_P;
#endif


/*************************************************************************/
/* External Global Definitions */
/*************************************************************************/

int  Hash_Mask_Bits;
int  Max_Hash_Strings;
int  Max_Hash_Data_Len;
int  Max_Frags_In_Memory_Store;
  // The number of fragments to read in a batch when streaming
  // the old reads against the hash table.

int  Contig_Mode = FALSE;
uint32  Last_Hash_Frag_Read;
int  LSF_Mode = FALSE;
int  Lo_Hash_Frag = 0;
int  Hi_Hash_Frag = INT_MAX;
int  Lo_Old_Frag = 0;
int  Hi_Old_Frag = INT_MAX;
int  Num_PThreads = DEFAULT_NUM_PTHREADS;
int64  Olap_Ct = 0;
int  Table_Ct = 0;
clock_t  Start_Time = 0, Stop_Time = 0;


Screen_Range_t  * Screen_Space;

int  Screen_Space_Size;
int  * Screen_Sub = NULL;

char  Sequence_Buffer [2 * AS_READ_MAX_LEN];
char  Quality_Buffer [2 * AS_READ_MAX_LEN];

FILE  * BOL_File = NULL;
FILE  * Kmer_Skip_File = NULL;
    // File of kmers to be ignored in the hash table
    // Specified by the  -k  option
Input_Stream  In_Stream = NULL;
Output_Stream  Out_Stream = NULL;
FragStore  OldFragStore;
FragStore  BACtigStore;
DistStore  OldDistStore;
char  * Frag_Store_Path;
char  * BACtig_Store_Path = NULL;
uint32  * IID_List = NULL;
int  IID_List_Len = 0;

MesgReader  Read_Msg_Fn;
MesgWriter  Write_Msg_Fn, Error_Write;

#if  USE_THREADS
pthread_mutex_t  FragStore_Mutex;
pthread_mutex_t  Write_Proto_Mutex;
pthread_mutex_t  Log_Msg_Mutex;
#endif

#if  ANALYZE_HITS
FILE  * High_Hits_File = NULL;
#endif


/*************************************************************************/
/* Function prototypes for internal static functions */
/*************************************************************************/

static void  Add_Match
    (String_Ref_t, int *, int, int *, Work_Area_t *);
static void  Add_Overlap
    (int, int, int, int, double, Olap_Info_t *, int *, Work_Area_t *);
static void  Add_Ref
    (String_Ref_t, int, Work_Area_t *);
static int  Binomial_Bound
    (int, double, int, double);
static int  Binomial_Hit_Limit
    (double p, double trials, double prob_bound);
static int  By_Diag_Sum
    (const void * a, const void * b);
static int  By_uint32
    (const void * a, const void * b);
static void  Capped_Increment
    (unsigned char *, unsigned char, unsigned char);
static char  char_Min
    (char a, char b);
static void  Combine_Into_One_Olap
    (Olap_Info_t * olap, int ct, int deleted []);
static char  Complement
    (char);
DistStore  Dist_Store_Open
    (char * storename, char * mode);
static void  Dump_Screen_Info
    (int frag_id, Screen_Info_t * screen, char dir);
static Overlap_t  Extend_Alignment
    (Match_Node_t *, char *, int, char *, int, int *, int *,
     int *, int *, int *, Work_Area_t *);
#if  USE_SOURCE_FIELD
static void  Find_End_Annotations
    (char * s, int frag_len, char * L, char * R);
#endif
static void  Find_Overlaps
    (char [], int, char [], Int_Frag_ID_t, Direction_t, Work_Area_t *);
static void  Flip_Screen_Range
    (Screen_Info_t * screen, int len);
FragStore  Frag_Store_Open
    (char * storename, char * mode);
static void  Get_Range
    (char * s, char ch, int * lo, int * hi);
static void  Get_Relevant_Screen_Matches
    (Screen_Info_t * screen);
static int  Has_Bad_Window
    (char a [], int n, int window_len, int threshold);
static String_Ref_t  Hash_Find
    (uint64 Key, int64 Sub, char * S, int64 * Where, int * hi_hits);
static void  Hash_Insert
    (String_Ref_t, uint64, char *);
static void  Hash_Mark_Empty
    (uint64 key, char * s);
static void  Initialize_Globals
    (void);
static int  Interval_Intersection
    (int, int, int, int);
static int  Lies_On_Alignment
    (int, int, int, int, Work_Area_t *);
static void  Mail_Error_Message
    (char * person, int num, char * store);
static void  Mark_Screened_Ends_Chain
    (String_Ref_t ref);
static void  Mark_Screened_Ends_Single
    (String_Ref_t ref);
static void  Mark_Skip_Kmers
    (void);
static void  Merge_Intersecting_Olaps
    (Olap_Info_t * p, int ct, int deleted []);
static int64  Next_Odd_Prime
    (int64);
static void  Output_Overlap
    (Int_Frag_ID_t, int, Direction_t, Int_Frag_ID_t,
     int, Olap_Info_t *);
static void  Output_Partial_Overlap
    (Int_Frag_ID_t s_id, Int_Frag_ID_t t_id, Direction_t dir,
     const Olap_Info_t * p, int s_len, int t_len);
static int  Passes_Screen
    (int, int, int);
static int  Prefix_Edit_Dist
    (char *, int, char *, int, int, int *, int *, int *, Work_Area_t *);
static void  Process_Matches
    (int *, char *, int, char * S_quality, Int_Frag_ID_t, Direction_t,
     char *, Hash_Frag_Info_t t_info, char * T_quality, Int_Frag_ID_t,
     Work_Area_t *, int);
static int  Process_String_Olaps
    (char *, int, char * quality, Int_Frag_ID_t, Direction_t,
     Work_Area_t *);
static void  Put_String_In_Hash
    (int i);
static int  Read_Next_Frag
    (char frag [AS_READ_MAX_LEN + 1], char quality [AS_READ_MAX_LEN + 1],
     FragStreamHandle stream, ReadStructp, Screen_Info_t *,
     uint32 * last_frag_read);
static void  Read_uint32_List
    (char * file_name, uint32 * * list, int * n);
static void  Rev_Complement
    (char *, int);
static int  Rev_Prefix_Edit_Dist
    (char *, int, char *, int, int, int *, int *, int *, int *,
     Work_Area_t *);
static void  Reverse_String
    (char *, int);
static void  Set_Left_Delta
    (int e, int d, int * leftover, int * t_end, int t_len,
     Work_Area_t * WA);
static void  Set_Right_Delta
    (int e, int d, Work_Area_t * WA);
static void  Set_Screened_Ends
    (void);
static void  Show_Alignment
    (char *, char *, Olap_Info_t *);
static void  Show_Match
    (Match_Node_t *, char *, char *);
static void  Show_Overlap
    (char *, int, char * a_quality,
     char *, int, char * b_quality, Olap_Info_t *);
#if  SHOW_SNPS
static void  Show_SNPs
    (char * a, int a_len, char * a_quality,
     char * b, int b_len, char * b_quality, Olap_Info_t * p,
     Work_Area_t * WA, int * mismatch_ct, int * indel_ct,
     int * olap_bases, int * unscreened_olap_bases,
     int * all_match_ct, int * hi_qual_ct,
     int local_msnp_bin [], int local_isnp_bin [],
     int local_match_bin [], int local_other_bin []);
#endif



int  main  (int argc, char * argv [])

  {
   char  File_Name [MAX_NAME_LEN];
   char  * bolfile_name = NULL, * Outfile_Name = NULL;
   char  * iidlist_file_name = NULL;
   int  illegal, create, append, force, exists;
   char  * p;
   int  noOverlaps;  /* If 1, run but don't compute/generate overlaps */

   assert (8 * sizeof (uint64) > 2 * WINDOW_SIZE);
//   assert (sizeof (Hash_Bucket_t) <= 64);

   // Set default hash-table parameters;  can be changed by -M option
   Hash_Mask_Bits = DEF_HASH_MASK_BITS;
   Max_Hash_Strings = DEF_MAX_HASH_STRINGS;
   Max_Hash_Data_Len = DEF_MAX_HASH_DATA_LEN;
   Max_Frags_In_Memory_Store
        = OVL_Min_int (Max_Hash_Strings, MAX_OLD_BATCH_SIZE);

#ifdef CONTIG_OVERLAPPER_VERSION
fprintf (stderr, "Running Contig version, AS_READ_MAX_LEN = %d\n",
         AS_READ_MAX_LEN);
#else
fprintf (stderr, "Running Fragment version, AS_READ_MAX_LEN = %d\n",
         AS_READ_MAX_LEN);
#endif

{
 time_t  now = time (NULL);
 fprintf (stderr, "### Starting at  %s\n", ctime (& now));
}

fprintf (stderr, "### Bucket size = " F_SIZE_T " bytes\n", sizeof (Hash_Bucket_t));
fprintf (stderr, "### Read error rate = %.2f%%\n", 100.0 * AS_READ_ERROR_RATE);
fprintf (stderr, "### Guide error rate = %.2f%%\n", 100.0 * AS_GUIDE_ERROR_RATE);
   noOverlaps = 0; /* If 1, don't compute/generate overlaps */
   Write_Msg_Fn = OutputFileType_AS (AS_BINARY_OUTPUT);

   illegal = 0;
   create = 1;
   append = 0;
   force = 0;
   exists = 0;

   { /* Parse the argument list using "man 3 getopt". */ 
     int ch, errflg = 0;
     optarg = NULL;
     while  (! errflg
               && ((ch = getopt (argc, argv, "ab:cfGh:I:k:K:l:mM:no:Pqr:st:uw")) != EOF))
       switch  (ch)
         {
          case  'a' :
            append = 1;
            create = 0;
            force  = 0;
            break;
          case  'b' :
            bolfile_name = strdup (optarg);
            assert (bolfile_name != NULL);
            break;
          case  'c':
            Contig_Mode = TRUE;
            break;
          case  'f' :
            force  = 1;
            append = 0;
            break;
          case  'G' :
            Doing_Partial_Overlaps = TRUE;
            Branch_Match_Value = PARTIAL_BRANCH_MATCH_VAL;
            Branch_Error_Value = Branch_Match_Value - 1.0;
            break;
          case  'h' :
            Get_Range (optarg, ch, & Lo_Hash_Frag, & Hi_Hash_Frag);
            LSF_Mode = TRUE;
            force = create = append = 0;
            break;
          case  'I' :
            iidlist_file_name = strdup (optarg);
            assert (iidlist_file_name != NULL);
            break;
          case  'k' :
            Kmer_Skip_File = File_Open (optarg, "r");
            break;
          case  'K' :
            if  (optarg [0] == '[')
                {
                 if  (sscanf (optarg, "[%lf,%lf,%lf]", & Genome_Len,
                              & Kmer_Freq_Bound, & Kmer_Prob_Limit)
                        != 3)
                     {
                      fprintf (stderr,
                               "ERROR:  Bad format for K option.  \n"
                               "Should be \"[%%lf,%%d,%%lf]\"\n");
                      exit (EXIT_FAILURE);
                     }

                 // Hi_Hit_Limit  needs to be calculated for each hash table
                 // since it depends on the amount of data in the table.

                 fprintf (stderr,
                          "Genome_Len = %.0lf  Kmer_Freq_Bound = %.0lf"
                          "  Kmer_Prob_Limit = %.4e\n",
                          Genome_Len, Kmer_Freq_Bound, Kmer_Prob_Limit);
                }
              else
                {
                 Hi_Hit_Limit = (int) strtol (optarg, & p, 10);
                 if  (p == optarg || Hi_Hit_Limit <= 1)
                     {
                      fprintf (stderr, "ERROR:  Illegal number kmer hit limit \"%s\"\n",
                               optarg);
                      exit (-1);
                     }
                 if  (Hi_Hit_Limit > HIGHEST_KMER_LIMIT)
                     fprintf (stderr,
                              "WARNING:  kmer hit limit = %d is too large, ignored\n",
                              Hi_Hit_Limit);
                }
            break;
          case  'l' :
            Frag_Olap_Limit = (int) strtol (optarg, & p, 10);
            if  (p == optarg)
                {
                 fprintf (stderr,
                          "ERROR:  Illegal number for frag olap limit \"%s\"\n",
                          optarg);
                 exit (-1);
                }
            if  (Frag_Olap_Limit < 1)
                Frag_Olap_Limit = INT_MAX;
            break;
          case  'm' :
            Unique_Olap_Per_Pair = FALSE;
            break;
          case  'M' :
#ifdef CONTIG_OVERLAPPER_VERSION
            fprintf (stderr, "M option not allowed for Contig Version--ignored\n");
#else
            if  (strcmp (optarg, "8GB") == 0)
                { // Set parameters for 8GB memory machine
                 Hash_Mask_Bits = 24;
                 Max_Hash_Strings = 500000;
                 Max_Hash_Data_Len = 360000000;
                 Max_Frags_In_Memory_Store
                      = OVL_Min_int (Max_Hash_Strings, MAX_OLD_BATCH_SIZE);
                }
            else if  (strcmp (optarg, "4GB") == 0)
                { // Set parameters for 4GB memory machine
                 Hash_Mask_Bits = 23;
                 Max_Hash_Strings = 250000;
                 Max_Hash_Data_Len = 180000000;
                 Max_Frags_In_Memory_Store
                      = OVL_Min_int (Max_Hash_Strings, MAX_OLD_BATCH_SIZE);
                }
            else if  (strcmp (optarg, "2GB") == 0)
                { // Set parameters for 2GB memory machine
                 Hash_Mask_Bits = 22;
                 Max_Hash_Strings = 100000;
                 Max_Hash_Data_Len = 75000000;
                 Max_Frags_In_Memory_Store
                      = OVL_Min_int (Max_Hash_Strings, MAX_OLD_BATCH_SIZE);
                }
            else if  (strcmp (optarg, "1GB") == 0)
                { // Set parameters for 1GB memory machine
                 Hash_Mask_Bits = 21;
                 Max_Hash_Strings = 40000;
                 Max_Hash_Data_Len = 30000000;
                 Max_Frags_In_Memory_Store
                      = OVL_Min_int (Max_Hash_Strings, MAX_OLD_BATCH_SIZE);
                }
            else if  (strcmp (optarg, "256MB") == 0)
                { // Set parameters for 256MB memory machine
                 Hash_Mask_Bits = 19;
                 Max_Hash_Strings = 10000;
                 Max_Hash_Data_Len = 6000000;
                 Max_Frags_In_Memory_Store
                      = OVL_Min_int (Max_Hash_Strings, MAX_OLD_BATCH_SIZE);
                }
              else
                {
                 fprintf (stderr, "ERROR:  Unknown memory size \"%s\"\n",
                      optarg);
                 errflg ++;
                }
#endif
            break;
          case  'n' :
            noOverlaps = 1;
            fprintf(stderr,"* No Overlaps will be generated!\n");
            break;
          case  'o' :
            Outfile_Name = strdup (optarg);
            assert (Outfile_Name != NULL);
            break;
          case  'P' :
            Write_Msg_Fn = OutputFileType_AS(AS_PROTO_OUTPUT);
            break;
          case  'q' :
            Single_Line_Output = TRUE;
            break;
          case  'r' :
            Get_Range (optarg, ch, & Lo_Old_Frag, & Hi_Old_Frag);
            break;
          case  's' :
            Ignore_Screen_Info = TRUE;
            break;
          case  't' :
            Num_PThreads = (int) strtol (optarg, & p, 10);
            if  (p == optarg)
                {
                 fprintf (stderr, "ERROR:  Illegal number of threads \"%s\"\n",
                          optarg);
                 exit (-1);
                }
            if  (Num_PThreads < 1)
                Num_PThreads = 1;
            break;
          case  'u' :
            Unique_Olap_Per_Pair = TRUE;
            break;
          case  'w' :
            Use_Window_Filter = TRUE;
            break;
          case  '?' :
            fprintf (stderr, "Unrecognized option -%c\n", optopt);
          default :
            errflg++;
         }

#if  ! USE_THREADS
     Num_PThreads = 1;
#endif

     if(force == 1 && append  == 1)
       {
         fprintf (stderr,
                  "* Illegal combination of command line flags"
                  "-- they are mutually exclusive\n");
         illegal = 1;
       }
     
     if  ((illegal == 1)
          || (argc - optind != 2 && (! LSF_Mode || Contig_Mode))
          || (argc - optind != 1 && (LSF_Mode && ! Contig_Mode)))
       {
         fprintf (stderr,
                  "USAGE:  %s [options] "
                  "<FragStorePath> <InputFilename>[.<ext>]\n"
                  "Opens <InputFilename>.<ext> to read .inp/.urc input\n"
                  "Creates/updates FragStore in <FragStorePath>\n"
                  "Writes .ovl output to <InputFilename>.ovl\n"
                  "Use -a to append to frag store\n"
                  "Use -f to force a new frag store\n"
                  "Use -G to do partial overlaps\n"
                  "Use -h <range> to specify fragments to put in hash table\n"
                  "    Implies LSF mode (no changes to frag store)\n"
                  "Use -I to designate a file of frag iids to limit olaps to\n"
                  "    (Contig mode only)\n"
                  "Use -k to specify a filename containing a list of kmers\n"
                  "    to ignore in the hash table\n"
                  "Use -K to designate limit on repetitive kmer hits to\n"
                  "    ignore in hash table\n"
                  "    [g,c,p] is special format to set K value automatically\n"
                  "    so that the probability of K or more hits in a genome\n"
                  "    of length g containing c copies of the kmer is\n"
                  "    less than p  (Put in quotes to get past the shell.)\n"
                  "Use -l to specify the maximum number of overlaps per\n"
                  "    fragment-end per batch of fragments.\n"
                  "Use -m to allow multiple overlaps per oriented fragment pair\n"
                  "Use -M to specify memory size.  Valid values are '8GB', '4GB',\n"
                  "    '2GB', '1GB', '256MB'.  (Not for Contig mode)\n"
                  "Use -n to skip overlaps (only update frag store)\n"
                  "Use -o to specify output file name\n"
                  "Use -P to force ASCII output\n"
                  "Use -q to make single-line output (same format as -G)\n"
                  "Use -r <range> to specify old fragments to overlap\n"
                  "Use -s to ignore screen information with fragments\n"
                  "Use -t to designate number of parallel threads\n"
                  "Use -u to allow only 1 overlap per oriented fragment pair\n"
                  "Use -w to filter out overlaps with too many errors in a window\n",
                  argv [0]);
         exit (EXIT_FAILURE);
       }

     Frag_Store_Path = argv [optind ++];

     if  (optind < argc)
         {
          if  (LSF_Mode)
              {
               if  (! Contig_Mode)
                   {
                    fprintf (stderr,
                             "ERROR:  Extraneous command-line argument \"%s\"\n",
                             argv [optind]);
                    exit (EXIT_FAILURE);
                   }

               if  (bolfile_name == NULL)
                   {
                    fprintf (stderr,
                             "ERROR:  No .BOL file name specified\n");
                    exit (EXIT_FAILURE);
                   }
               BOL_File = File_Open (bolfile_name, "w");

               if  (Outfile_Name != NULL)
                   {
                    fprintf (stderr,
                             "WARNING:  Output file name \"%s\" ignored\n",
                             Outfile_Name);
                    // Only .bol file is produced
                   }

               BACtig_Store_Path = argv [optind];
              }
            else
              {
               char  * suffix;

               suffix = strrchr (argv [optind], (int) '.');
               strcpy (File_Name, argv [optind]);
	       assert(NULL == In_Stream);
               In_Stream = File_Open (File_Name, "r");     // inp file
	       assert(NULL != In_Stream);
               Read_Msg_Fn = InputFileType_AS (In_Stream);

               if  (Outfile_Name == NULL)
                   {
                    if  (suffix)
                        * suffix = '\0';
                    if  (Contig_Mode)
                        {
                         sprintf (File_Name, "%s.bol", argv [optind]);
			 assert(NULL == BOL_File);
                         BOL_File = File_Open (File_Name, "w");
                        }
                    sprintf (File_Name, "%s.ovl", argv [optind]);
		    assert(NULL == Out_Stream);
                    Out_Stream = File_Open (File_Name, "w");     // ovl file
                   }
                 else
		   assert(NULL == Out_Stream);
                   Out_Stream = File_Open (Outfile_Name, "w");     // ovl file
              }

          optind ++;
         }
     else if  (! LSF_Mode)
         {
          fprintf (stderr, "ERROR:  No input file name specified\n");
          exit (EXIT_FAILURE);
         }
     else if  (Outfile_Name == NULL)
         {
          fprintf (stderr, "ERROR:  No output file name specified\n");
          exit (EXIT_FAILURE);
         }
       else
         {
	   assert(NULL == Out_Stream);
	   Out_Stream = File_Open (Outfile_Name, "w");     // ovl file
	   free (Outfile_Name);
         }

     /* End of command line parsing */
   }
   
#if  SHOW_STATS
 assert(NULL == Stat_File);
 Stat_File = File_Open (STAT_FILE_NAME, "w");
 {
   int  i;
   fprintf (Stat_File, ">");
   for  (i = 0;  i < argc;  i ++)
     fprintf (Stat_File, " %s", argv [i]);
   fprintf (Stat_File, "\n");
 }
#endif

#if  LIST_TOTALLY_SCREENED
   Total_Screen_File = File_Open (TOTAL_SCREEN_FILE_NAME, "w");
#endif

   fprintf (stderr, "    Hash_Mask_Bits = %d\n", Hash_Mask_Bits);
   fprintf (stderr, "  Max_Hash_Strings = %d\n", Max_Hash_Strings);
   fprintf (stderr, " Max_Hash_Data_Len = %d\n", Max_Hash_Data_Len);
   fprintf (stderr, "       Kmer Length = %d\n", WINDOW_SIZE);

   Initialize_Globals ();

   if  (LSF_Mode && Contig_Mode)
       {
        if  (BACtig_Store_Path == NULL)
            {
             fprintf (stderr, "ERROR:  No BACtig store specified\n");
             exit (EXIT_FAILURE);
            }
        if  (! existsFragStore (BACtig_Store_Path))
            {
             fprintf (stderr, "ERROR:  BACtig store \"%s\" does NOT exist\n",
                      BACtig_Store_Path);
             exit (EXIT_FAILURE);
            }
        BACtigStore = Frag_Store_Open (BACtig_Store_Path, "r");

        if  (iidlist_file_name != NULL)
            {
             fprintf (stderr, "### Using IID list = %s\n",
                      iidlist_file_name);
             Read_uint32_List (iidlist_file_name, & IID_List, & IID_List_Len);
             if  (IID_List_Len <= 0)
                 {
                  fprintf (stderr, "ERROR:  Empty IID_List\n");
                  exit (EXIT_FAILURE);
                 }
             if  (Lo_Old_Frag != 0 || Hi_Old_Frag != INT_MAX)
                 {
                  fprintf (stderr,
                           "WARNING:  IID list specified, range from -r ignored\n");
                  Lo_Old_Frag = 0;
                  Hi_Old_Frag = INT_MAX;
                 }
            }
       }

   exists = existsFragStore (Frag_Store_Path);
   if  (exists == 1)
       {
        char  buffer [MAX_NAME_LEN];

        sprintf (buffer, "%s/db.dst", Frag_Store_Path);
        if  (! File_Exists (buffer))
            exists = -1;
       }

   if  (force && exists)
       {
        char cmd[2048];

        fprintf (stderr, "* Fragment Store %s exists ...nuking\n", Frag_Store_Path);
        removeFragStore (Frag_Store_Path);
        sprintf (cmd,"rm %s/db.dst", Frag_Store_Path);
        if(system(cmd) != 0) assert(0);
        create = 1;
        append = 0;
       }
   else if  (create && exists)
       {
        fprintf (stderr, "\a\aFragment Store %s exists, use  -f  to overwrite\n",
                 Frag_Store_Path);
        exit (-1);
       }

   if  (append)
       {
        if  (exists == 0)
            {
             fprintf (stderr, "* Directory %s DOES NOT exist ...creating before append\n",
                      Frag_Store_Path);
             create = 1;
             append = 0;
            }
        else if  (exists == -1)
            {
             fprintf (stderr, "* Directory %s contains some, but not all files, exiting\n",
                      Frag_Store_Path);
             exit(1);
            }
          else
            {
             fprintf (stderr, "* Directory %s DOES  exist ...\n",
                      Frag_Store_Path);
             append = 1;
             create = 0;
            }
       }

   if  (append)
       {
        char buffer[FILENAME_MAX];

        // fprintf(stderr,"* Appending to fragstore %s\n", Frag_Store_Path);
        sprintf(buffer,"%s/db.dst", Frag_Store_Path);
        OldDistStore = Dist_Store_Open (buffer, "r+"); 
        OldFragStore = Frag_Store_Open (Frag_Store_Path, "r+");
        // fprintf(stderr,"** Opened store %d for OldFragStore \n", OldFragStore);
       }
   else if  (create)
       {
        char buffer[FILENAME_MAX];

        // fprintf (stderr, "Frag_Store_Path = %p\n", Frag_Store_Path);
        //fprintf (stderr, "Frag_Store_Path [0] = %d\n", (int) (Frag_Store_Path [0]));
        fprintf(stderr,"* Creating fragstore %s\n", Frag_Store_Path);
        sprintf(buffer,"%s/db.dst", Frag_Store_Path);
        OldFragStore = createFragStore(Frag_Store_Path, "Just testing", 1);
        // fprintf(stderr,"** Allocated store %d for OldFragStore \n", OldFragStore);
        OldDistStore = createDistStore(buffer,1); 
       }
   else if  (LSF_Mode)
       {
        char  buffer [FILENAME_MAX];

        fprintf (stderr, "* Using fragstore %s\n", Frag_Store_Path);
        sprintf (buffer, "%s/db.dst", Frag_Store_Path);
//        OldDistStore = Dist_Store_Open (buffer, "r"); 
        OldFragStore = Frag_Store_Open (Frag_Store_Path, "r");
       }
     else
       {
        fprintf (stderr,"** Serious error...bye\n");
        exit (-1);
       }
   
   if  (LSF_Mode && ! Contig_Mode)
       {
        GenericMesg   * pmesg = (GenericMesg *) Safe_malloc (sizeof (GenericMesg));
        AuditMesg  * adtmesg = (AuditMesg *) Safe_malloc (sizeof (AuditMesg));

        pmesg->t = MESG_ADT;
        adtmesg -> list = NULL;
        VersionStampADT(adtmesg, argc, argv);

	pmesg -> m = adtmesg;
#if  ! (FOR_CARL_FOSLER || SHOW_SNPS)
        if  (! Doing_Partial_Overlaps && ! Single_Line_Output)
            Write_Msg_Fn (Out_Stream, pmesg);
#endif

        free (adtmesg);
        free (pmesg);
       }

   /****************************************/
   OverlapDriver(noOverlaps, argc, argv);
   /****************************************/

#ifdef  CHECK_SCREENING
Dump_Screen_Info (0, NULL, 'x');
#endif

   if(NULL != In_Stream) { fclose (In_Stream); In_Stream = NULL;}
   assert(NULL != Out_Stream);
   fclose (Out_Stream);
   Out_Stream = NULL;
   if  (Contig_Mode) {
     assert(NULL != BOL_File);
     fclose (BOL_File);
     BOL_File = NULL;
   }

#if  SHOW_STATS
fprintf (Stat_File, "Regular Overlaps:\n %10llu\n",
           Regular_Olap_Ct);
fprintf (Stat_File, "Duplicate Olaps:\n %10llu\n",
           Duplicate_Olap_Ct);
fprintf (Stat_File, "Too Short:\n %10llu\n",
           Too_Short_Ct);
Print_Distrib (Olap_Len_Dist, "Overlap Lengths");
Print_Distrib (Num_Olaps_Dist, "Numbers of Overlaps:");
Print_Distrib (Kmer_Freq_Dist, "Kmer Frequency in Hash Table:");
Print_Distrib (Kmer_Hits_Dist, "Kmer Hits per Fragment:");
Print_Distrib (String_Olap_Dist, "String Olaps per Fragment (each direction):");
Print_Distrib (Hits_Per_Olap_Dist, "Kmer Hits\n  per String Olap:");
Print_Distrib (Diag_Dist, "Diagonals per String Olap:");
Print_Distrib (Gap_Dist, "Num Gaps in\n  Nice Olaps:");
Print_Distrib (Exacts_Per_Olap_Dist, "Exact Matches per\n  String Olap:");
Print_Distrib (Edit_Depth_Dist, "Edit Distance Depth:");
fprintf (Stat_File, "Edit_Dist_Ct = %lld\n", Edit_Dist_Ct);
fprintf (Stat_File, "Matches = %lld\n", Match_Ct);
fprintf (Stat_File, "%lld collisions in %lld finds\n", Collision_Ct, Hash_Find_Ct);
// fprintf (Stat_File, "String_Olap_Size = %ld\n", WA -> String_Olap_Size);
fprintf (stderr, "Stats written to file \"%s\"\n", STAT_FILE_NAME);
fclose (Stat_File);
#endif

#if  LIST_TOTALLY_SCREENED
   fclose (Total_Screen_File);
#endif

   fprintf (stderr, " Kmer hits without olaps = %lld\n", Kmer_Hits_Without_Olap_Ct);
   fprintf (stderr, "    Kmer hits with olaps = %lld\n", Kmer_Hits_With_Olap_Ct);
   fprintf (stderr, "  Multiple overlaps/pair = %lld\n", Multi_Overlap_Ct);
   fprintf (stderr, " Total overlaps produced = %lld\n", Total_Overlaps);
   fprintf (stderr, "      Contained overlaps = %lld\n", Contained_Overlap_Ct);
   fprintf (stderr, "       Dovetail overlaps = %lld\n", Dovetail_Overlap_Ct);
   fprintf (stderr, "Rejected by short window = %lld\n", Bad_Short_Window_Ct);
   fprintf (stderr, " Rejected by long window = %lld\n", Bad_Long_Window_Ct);

   closeFragStore (OldFragStore);

#if  SHOW_SNPS
{
 int  total [20] = {0};
 int  a, b, c, d, e;
 int  i;

 fprintf (stderr, "%5s  %9s  %9s  %9s  %9s  %9s\n",
          "Qual", "MatchSNPs", "IndelSNPs", "AllMatch", "Other", "Total");
 for  (i = 3;  i < 13;  i ++)
   {
    total [i] = Global_Match_SNP_Bin [i]
                  + Global_Indel_SNP_Bin [i]
                  + Global_All_Match_Bin [i]
                  + Global_Hi_Qual_Bin [i];
    fprintf (stderr, "%2d-%2d  %9d  %9d  %9d  %9d  %9d\n",
             5 * i, 5 * i + 4,
             Global_Match_SNP_Bin [i],
             Global_Indel_SNP_Bin [i],
             Global_All_Match_Bin [i],
             Global_Hi_Qual_Bin [i],
             total [i]);
   }

 a = b = c = d = e = 0;
 fprintf (stderr, "\n%5s  %9s  %9s  %9s  %9s  %9s\n",
          "Qual", "MatchSNPs", "IndelSNPs", "AllMatch", "Other", "Total");
 for  (i = 12;  i >= 3;  i --)
   {
    a += Global_Match_SNP_Bin [i];
    b += Global_Indel_SNP_Bin [i];
    c += Global_All_Match_Bin [i];
    d += Global_Hi_Qual_Bin [i];
    e += total [i];
    fprintf (stderr, " >=%2d  %9d  %9d  %9d  %9d  %9d\n",
             5 * i, a, b, c, d, e);
   }

 fprintf (stderr, "\nTotal overlap bases = %d\n", Global_Olap_Ct);
 fprintf (stderr, "Total unscreened bases = %d\n", Global_Unscreened_Ct);
}
#endif

{
 time_t  now = time (NULL);
 fprintf (stderr, "### Finished at  %s\n", ctime (& now));
}

#if  DO_OLAP_ALIGN_PROFILE
//  Print the alignment counts for the fragments in the hash table
{
 int  i, k;
 int  len;

 for  (k = Lo_Hash_Frag;  k <= Hi_Hash_Frag;  k ++)
   {
    char  * Align_String = Data + String_Start [k - Hash_String_Num_Offset];

    printf ("\n> Frag %d:\n", k);

    len = strlen (Align_String);
    Align_P = Align_Ct [k - Hash_String_Num_Offset];

    for  (i = 0;  i < len;  i ++)
      {
       int  j, sum = 0, max = -1;

       printf ("%c", Align_String [i]);
       for  (j = 0;  j < 6;  j ++)
         {
          printf (" %6d", Align_P [i] [j]);
          if  (Align_P [i] [j] > max)
              max = Align_P [i] [j];
          sum += Align_P [i] [j];
         }
       if  (Align_P [i] [Bit_Equivalent [Align_String [i]]] != max)
           printf ("  **** not majority");
       else if  (Align_P [i] [Bit_Equivalent [Align_String [i]]] < 0.9 * sum)
           printf ("  <<<< %4.1f%%",
                   100.0 * Align_P [i] [Bit_Equivalent [Align_String [i]]] / sum);
       putchar ('\n');
      }
   }
}
#endif

fprintf (stderr, "### Return from main\n");
   return  0;
  }



static void  Add_Match
    (String_Ref_t ref, int * start, int offset, int * consistent,
     Work_Area_t * wa)

//  Add information for the match in  ref  to the list
//  starting at subscript  (* start) .  The matching window begins
//  offset  bytes from the beginning of this string.

  {
   int  * p, save;
   int  diag = 0, new_diag, expected_start = 0, num_checked = 0;
   int  move_to_front = FALSE;

#if  SHOW_STATS
Kmer_Hits_Ct ++;
#endif

   new_diag = ref . Offset - offset;

   for  (p = start;  (* p) != 0;  p = & (wa -> Match_Node_Space [(* p)] . Next))
     {
      expected_start = wa -> Match_Node_Space [(* p)] . Start
                         + wa -> Match_Node_Space [(* p)] . Len
                         - WINDOW_SIZE + 1 + HASH_KMER_SKIP;
//  Added  HASH_KMER_SKIP here --------------^^^^^^^^^^^^^^
      diag = wa -> Match_Node_Space [(* p)] . Offset
               - wa -> Match_Node_Space [(* p)] . Start;

      if  (expected_start < offset)
          break;
      if  (expected_start == offset)
          {
           if  (new_diag == diag)
               {
                wa -> Match_Node_Space [(* p)] . Len += 1 + HASH_KMER_SKIP;
                if  (move_to_front)
                    {
                     save = (* p);
                     (* p) = wa -> Match_Node_Space [(* p)] . Next;
                     wa -> Match_Node_Space [save] . Next = (* start);
                     (* start) = save;
                    }
                return;
               }
             else
               move_to_front = TRUE;
          }
      num_checked ++;
     }

   if  (wa -> Next_Avail_Match_Node == wa -> Match_Node_Size)
       {
        wa -> Match_Node_Size = (int32) (wa -> Match_Node_Size *
               MEMORY_EXPANSION_FACTOR);
        fprintf (stderr, "### reallocing  Match_Node_Space  Size = %d\n",
                 wa -> Match_Node_Size);
        wa -> Match_Node_Space = (Match_Node_t *) Safe_realloc
               (wa -> Match_Node_Space,
                wa -> Match_Node_Size * sizeof (Match_Node_t));
       }

   if  ((* start) != 0
          && (num_checked > 0
                || abs (diag - new_diag) > 3
                || offset < expected_start + WINDOW_SIZE - 2))
       (* consistent) = FALSE;

   save = (* start);
   (* start) = wa -> Next_Avail_Match_Node;
   wa -> Next_Avail_Match_Node ++;

   wa -> Match_Node_Space [(* start)] . Offset = ref . Offset;
   wa -> Match_Node_Space [(* start)] . Len = WINDOW_SIZE;
   wa -> Match_Node_Space [(* start)] . Start = offset;
   wa -> Match_Node_Space [(* start)] . Next = save;

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


static void  Add_Ref
    (String_Ref_t Ref, int Offset, Work_Area_t * WA)

//  Add information for  Ref  and all its matches to the global
//  hash table in  String_Olap_Space .  Grow the space if necessary
//  by  MEMORY_EXPANSION_FACTOR .  The matching window begins
//  Offset  bytes from the beginning of this string.

  {
   int32  Prev, Sub;
   int  consistent;

   Sub = (Ref . String_Num ^ (Ref . String_Num >> STRING_OLAP_SHIFT))
                  & STRING_OLAP_MASK;
   
   while  (WA -> String_Olap_Space [Sub] . Full
              && WA -> String_Olap_Space [Sub] . String_Num != Ref . String_Num)
     {
      Prev = Sub;
      Sub = WA -> String_Olap_Space [Sub] . Next;
      if  (Sub == 0)
          {
           if  (WA -> Next_Avail_String_Olap == WA -> String_Olap_Size)
               {
                WA -> String_Olap_Size = (int32) (WA -> String_Olap_Size *
                       MEMORY_EXPANSION_FACTOR);
                fprintf (stderr, "### reallocing  String_Olap_Space  Size = %d\n",
                         WA -> String_Olap_Size);
                WA -> String_Olap_Space = (String_Olap_t *) Safe_realloc
                       (WA -> String_Olap_Space,
                        WA -> String_Olap_Size * sizeof (String_Olap_t));
               }
           Sub = WA -> Next_Avail_String_Olap ++;
           WA -> String_Olap_Space [Prev] . Next = Sub;
           WA -> String_Olap_Space [Sub] . Full = FALSE;
           break;
          }
     }

   if  (! WA -> String_Olap_Space [Sub] . Full)
       {
        WA -> String_Olap_Space [Sub] . String_Num = Ref . String_Num;
        WA -> String_Olap_Space [Sub] . Match_List = 0;
        WA -> String_Olap_Space [Sub] . diag_sum = 0.0;
        WA -> String_Olap_Space [Sub] . diag_ct = 0;
        WA -> String_Olap_Space [Sub] . Next = 0;
        WA -> String_Olap_Space [Sub] . Full = TRUE;
        WA -> String_Olap_Space [Sub] . consistent = TRUE;
#if  SHOW_STATS
WA -> String_Olap_Space [Sub] . Kmer_Hits = 0;
#endif
       }

#if  SHOW_STATS
WA -> String_Olap_Space [Sub] . Kmer_Hits ++;
#endif

   consistent = WA -> String_Olap_Space [Sub] . consistent;

   WA -> String_Olap_Space [Sub] . diag_sum += Ref . Offset - Offset;
   WA -> String_Olap_Space [Sub] . diag_ct ++;
   Add_Match (Ref, & (WA -> String_Olap_Space [Sub] . Match_List),
              Offset, & consistent, WA);

   WA -> String_Olap_Space [Sub] . consistent = consistent;

   return;
  }



int  Build_Hash_Index
    (FragStreamHandle stream, int32 first_frag_id, ReadStructp myRead)

/* Read the next batch of strings from  stream  and create a hash
*  table index of their  WINDOW_SIZE -mers.  Return  1  if successful;
*  0 otherwise.  The batch ends when either end-of-file is encountered
*  or  Max_Hash_Strings  have been read in.   first_frag_id  is the
*  internal ID of the first fragment in the hash table. */

  {
   String_Ref_t  ref;
   Screen_Info_t  screen;
   int64  total_len;
   static int64  max_extra_ref_ct = 0;
   static int64  old_ref_len, new_ref_len;
   int  frag_status;
   int64  i;
   int  screen_blocks_used = 0;
   int  hash_entry_limit;
   int  j;

   Hash_String_Num_Offset = first_frag_id;
   String_Ct = 0;
   total_len = 0;
   if  (Data == NULL)
       {
        Data_Len = Max_Hash_Data_Len + AS_READ_MAX_LEN;
        Data = (char *) Safe_realloc (Data, Data_Len);
        Quality_Data = (char *) Safe_realloc (Quality_Data, Data_Len);
        old_ref_len = Data_Len / (HASH_KMER_SKIP + 1);
        Next_Ref = (String_Ref_t *) Safe_realloc
             (Next_Ref, old_ref_len * sizeof (String_Ref_t));
       }

   memset (Next_Ref, '\377', old_ref_len * sizeof (String_Ref_t));
   memset (Hash_Table, 0, HASH_TABLE_SIZE * sizeof (Hash_Bucket_t));
   memset (Hash_Check_Array, 0, HASH_TABLE_SIZE * sizeof (Check_Vector_t));

   screen . match = (IntScreenMatch *) Safe_malloc
                         (INIT_SCREEN_MATCHES * sizeof (IntScreenMatch));
   screen . range = (Screen_Range_t *) Safe_malloc
                         (INIT_SCREEN_MATCHES * sizeof (Screen_Range_t));
   screen . match_len = INIT_SCREEN_MATCHES;
   screen . num_matches = 0;

   fprintf (stderr, "### Build_Hash:  first_frag_id = %d  Max_Hash_Strings = %d\n",
        first_frag_id, Max_Hash_Strings);

   if  (LSF_Mode)
       screen_blocks_used = 1;

   Extra_Ref_Ct = 0;
   Hash_Entries = 0;
   hash_entry_limit = MAX_HASH_LOAD * HASH_TABLE_SIZE * ENTRIES_PER_BUCKET;

   while  (String_Ct < Max_Hash_Strings
             && total_len < Max_Hash_Data_Len
             && Hash_Entries < hash_entry_limit
             && (frag_status
                     = Read_Next_Frag (Sequence_Buffer, Quality_Buffer, stream,
                                       myRead, & screen, & Last_Hash_Frag_Read)))
     {
      int  extra, len;
      size_t  new_len;

#if  USE_SOURCE_FIELD
      Left_Repeat_Tag [String_Ct] = Global_Left_Annotation;
      Right_Repeat_Tag [String_Ct] = Global_Right_Annotation;
      Other_Repeat_Tag [String_Ct] = Global_Other_Annotation;
#endif

      if  (frag_status == DELETED_FRAG)
          {
           Sequence_Buffer [0] = '\0';
           Quality_Buffer [0] = '\0';
          }
        else
          {
           if  (LSF_Mode)
               {
                if  (screen . num_matches == 0)
                    Screen_Sub [String_Ct] = 0;
                  else
                    {
                     int  new_size = screen_blocks_used + screen . num_matches;

                     while  (new_size >= Screen_Space_Size)
                         {
                          Screen_Space_Size *= MEMORY_EXPANSION_FACTOR;
                          fprintf (stderr, "### reallocing  Screen_Space  Size = %d\n",
                                   Screen_Space_Size);
                          Screen_Space
                              = (Screen_Range_t *) Safe_realloc
                                    (Screen_Space,
                                     Screen_Space_Size * sizeof (Screen_Range_t));
                         }

                     memmove (Screen_Space + screen_blocks_used, screen . range,
                              screen . num_matches * sizeof (Screen_Range_t));

                     Screen_Sub [String_Ct] = screen_blocks_used;
                     screen_blocks_used = new_size;
                    }
               }

           getReadType_ReadStruct (myRead, Kind_Of_Frag + String_Ct);
           if  (Contig_Mode)
               {
                getLocID_ReadStruct (myRead, Loc_ID + String_Ct);
               }
          }

      String_Start [String_Ct] = total_len;
      len = strlen (Sequence_Buffer);
      String_Info [String_Ct] . length = len;
      String_Info [String_Ct] . left_end_screened
          = screen . left_end_screened;
      String_Info [String_Ct] . right_end_screened
          = screen . right_end_screened;
      new_len = total_len + len + 1;
      extra = new_len % (HASH_KMER_SKIP + 1);
      if  (extra > 0)
          new_len += 1 + HASH_KMER_SKIP - extra;

#if  DO_OLAP_ALIGN_PROFILE
if  (Align_Ct_Size == 0)
    {
     Align_Ct_Size = 200;
     Align_Ct = (Align_Ct_t *)
                    Safe_malloc (Align_Ct_Size * sizeof (Align_Ct_t));
    }
if  (String_Ct >= Align_Ct_Size - 1)
    {
     Align_Ct_Size *= 2;
     Align_Ct = (Align_Ct_t *)
                    Safe_realloc (Align_Ct, Align_Ct_Size * sizeof (Align_Ct_t));
    }
Align_Ct [String_Ct] = (Align_Entry_t *)
                               Safe_calloc (1 + len, sizeof (Align_Entry_t));
#endif

      if  (new_len > Data_Len)
          {
           Data_Len = (size_t) (Data_Len * MEMORY_EXPANSION_FACTOR);
           if  (new_len > Data_Len)
               Data_Len = new_len;
           fprintf (stderr, "### reallocing  Data and Quality_Data  Data_Len = " F_SIZE_T "\n",
                    Data_Len);
           Data = (char *) Safe_realloc (Data, Data_Len);
           Quality_Data = (char *) Safe_realloc (Quality_Data, Data_Len);
           new_ref_len = Data_Len / (HASH_KMER_SKIP + 1);
           Next_Ref = (String_Ref_t *) Safe_realloc
                (Next_Ref, new_ref_len * sizeof (String_Ref_t));
           memset (Next_Ref + old_ref_len, '\377',
                (new_ref_len - old_ref_len) * sizeof (String_Ref_t));
           old_ref_len = new_ref_len;
          }

      strcpy (Data + total_len, Sequence_Buffer);
      memcpy (Quality_Data + total_len, Quality_Buffer, len + 1);
      total_len = new_len;

      Put_String_In_Hash (String_Ct);

      String_Ct ++;
     }

   if  (String_Ct == 0)
       {
        free (screen . match);
        free (screen . range);
        return  0;
       }

   fprintf (stderr, "strings read = %d  total_len = " F_S64 "\n",
            String_Ct, total_len);

   if  (Genome_Len > 0.0)
       {
        double  avg_frag_kmers;
        double  genome_kmers;

        avg_frag_kmers
            = ((double) total_len / String_Ct - WINDOW_SIZE)
                 / (1 + HASH_KMER_SKIP);
        genome_kmers = 2 * (Genome_Len - (WINDOW_SIZE - 1));
        Hi_Hit_Limit
            = 1 + Binomial_Hit_Limit (Kmer_Freq_Bound * avg_frag_kmers/genome_kmers,
                                      String_Ct - 1, Kmer_Prob_Limit);
        // "1 +" because one kmer is automatically present from first fragment
        // probability is for other fragments to hit it too

        if  (Hi_Hit_Limit < MIN_CALC_KMER)
            Hi_Hit_Limit = MIN_CALC_KMER;
       }
   fprintf (stderr, "### Kmer hash hit limit = %d\n", Hi_Hit_Limit);
   if  (Hi_Hit_Limit > HIGHEST_KMER_LIMIT)
       fprintf (stderr, "  is too large and has no effect\n");
     else
       Set_Screened_Ends ();

#if  SHOW_STATS
fprintf (Stat_File, "Num_Strings:\n %10ld\n", String_Ct);
fprintf (Stat_File, "total_len:\n %10ld\n", total_len);
#endif

   fprintf (stderr, "Hash_Entries = " F_S64 "  Load = %.1f%%\n",
        Hash_Entries, (100.0 * Hash_Entries) / (HASH_TABLE_SIZE * ENTRIES_PER_BUCKET));

   if  (Extra_Ref_Ct > max_extra_ref_ct)
       {
        max_extra_ref_ct *= MEMORY_EXPANSION_FACTOR;
        if  (Extra_Ref_Ct > max_extra_ref_ct)
            max_extra_ref_ct = Extra_Ref_Ct;
        fprintf (stderr,
                 "### realloc  Extra_Ref_Space  max_extra_ref_ct = " F_S64 "\n",
                 max_extra_ref_ct);
        Extra_Ref_Space = (String_Ref_t *) Safe_realloc (Extra_Ref_Space,
                               max_extra_ref_ct * sizeof (String_Ref_t));
       }

#if  SHOW_STATS
fprintf (Stat_File, "Hash Table Size:\n %10ld\n",
            HASH_TABLE_SIZE * ENTRIES_PER_BUCKET);
fprintf (Stat_File, "Collisions:\n %10ld\n", Collision_Ct);
fprintf (Stat_File, "Hash entries:\n %10ld\n", Hash_Entries);
Collision_Ct = 0;
#endif

#if  SHOW_STATS
{
 int64  Ct, sub;

 Hash_Entries = 0;
 for  (i = 0;  i < HASH_TABLE_SIZE;  i ++)
   for  (j = 0;  j < Hash_Table [i] . Entry_Ct;  j ++)
     {
      Ct = 0;
      for  (ref = Hash_Table [i] . Entry [j];  ! ref . Empty;
                        ref = Next_Ref [sub])
        {
         Ct ++;
         sub = (String_Start [ref . String_Num] + ref . Offset)
                 / (HASH_KMER_SKIP + 1);
        }
      if  (Ct > 0)
          Incr_Distrib (& Kmer_Freq_Dist, Ct);
     }
}
#endif

   if  (Kmer_Skip_File != NULL)
       Mark_Skip_Kmers ();

   // Coalesce reference chain into adjacent entries in  Extra_Ref_Space
   Extra_Ref_Ct = 0;
   for  (i = 0;  i < HASH_TABLE_SIZE;  i ++)
     for  (j = 0;  j < Hash_Table [i] . Entry_Ct;  j ++)
       {
        ref = Hash_Table [i] . Entry [j];
        if  (! ref . Last && ! ref . Empty)
            {
             Extra_Ref_Space [Extra_Ref_Ct] = ref;
             Hash_Table [i] . Entry [j] . String_Num = (Extra_Ref_Ct >> OFFSET_BITS);
             Hash_Table [i] . Entry [j] . Offset = (Extra_Ref_Ct & OFFSET_MASK);
             Extra_Ref_Ct ++;
             do
               {
                ref = Next_Ref [(String_Start [ref . String_Num] + ref . Offset)
                                  / (HASH_KMER_SKIP + 1)];
                Extra_Ref_Space [Extra_Ref_Ct ++] = ref;
               }  while  (! ref . Last);
            }
       }

#if  SHOW_PROGRESS
Stop_Time = clock ();
fprintf (stderr, "Time to build hash table for %d fragments = %.1f sec\n",
            String_Ct, (double) (Stop_Time - Start_Time) / CLOCKS_PER_SEC);
Start_Time = clock ();
#endif

   free (screen . match);
   free (screen . range);

   return  1;
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



static int  Binomial_Hit_Limit
    (double p, double n, double prob_bound)

//  Return  k  such that the probability of  >= k  successes in
//   n  trials is  < prob_bound , where  p  is the probability
//  of success of each trial.

  {
   double  lambda, target, sum, term, q;
   int  i;

   lambda = n * p;
   q = 1.0 - p;

   if  (lambda <= 5.0 && n >= 50.0)
       {  // use Poisson approximation
        target = (1.0 - prob_bound) * exp (lambda);
        sum = term = 1.0;
        for  (i = 1;  sum <= target && i < 50;  i ++)
          {
           term = term * lambda / (double) i;
           sum += term;
          }
        if  (sum > target)
            return  i;
       }
   if  (n >= 30.0)
       {  // use Normal approximation
        double  t, z;
        double  c [3] = {2.515517, 0.802853, 0.010328};
        double  d [4] = {1.0, 1.432788, 0.189269, 0.001308};

        if  (prob_bound <= 0.5)
            t = sqrt (-2.0 * log (prob_bound));
          else
            t = sqrt (-2.0 * log (1.0 - prob_bound));

        z = t - ((c [2] * t + c [1]) * t + c [0])
                  / (((d [3] * t + d [2]) * t + d [1]) * t + d [0]);

        if  (prob_bound <= 0.5)
            target = z;
          else
            target = -z;

        return  (int) ceil (lambda + sqrt (lambda * q) * target);
       }
     else
       {  // brute force
        target = 1.0 - prob_bound;
        sum = term = pow (q, n);
        for  (i = 1;  sum <= target && i < n;  i ++)
          {
           term *= (n + 1 - i) / i;
           term *= p / q;
           sum += term;
          }
        return  i;
       }
  }



static int  By_Diag_Sum
    (const void * a, const void * b)

//  Compare the  diag_sum  fields  in  a  and  b  as  (String_Olap_t *) 's and
//  return  -1  if  a < b ,  0  if  a == b , and  1  if  a > b .
//  Used for  qsort .

  {
   String_Olap_t  * x, * y;

   x = (String_Olap_t *) a;
   y = (String_Olap_t *) b;

   if  (x -> diag_sum < y -> diag_sum)
       return  -1;
   else if  (x -> diag_sum > y -> diag_sum)
       return  1;

   return  0;
  }



static int  By_uint32
    (const void * a, const void * b)

//  Compare the values in  a  and  b  as  (uint32 *) 's and
//  return  -1  if  a < b ,  0  if  a == b , and  1  if  a > b .
//  Used for  qsort .

  {
   uint32  * x, * y;

   x = (uint32 *) a;
   y = (uint32 *) b;

   if  ((* x) < (* y))
       return  -1;
   else if  ((* x) > (* y))
       return  1;

   return  0;
  }



static void  Capped_Increment
    (unsigned char * a, unsigned char b, unsigned char cap)

//  Set  (* a)  to the min of  a + b  and  cap .

  {
   int  sum;

   sum = (* a) + b;
   if  (sum < cap)
       (* a) = sum;
     else
       (* a) = cap;

   return;
  }



static char  char_Min
    (char a, char b)

//  Return the smaller of  a  and  b .

  {
   if  (a < b)
       return  a;

   return  b;
  }



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



static char  Complement
    (char Ch)

/*  Return the DNA complement of  Ch . */

  {
   switch  (tolower ((int) Ch))
     {
      case  'a' :
        return  't';
      case  'c' :
        return  'g';
      case  'g' :
        return  'c';
      case  't' :
        return  'a';
      case  DONT_KNOW_CHAR :
        return  DONT_KNOW_CHAR;
      default :
        fprintf (stderr, "ERROR(complement):  Unexpected character `%c\'\n", Ch);
        exit (-1);
     }

   return  'x';    // Just to make the compiler happy
  }



DistStore  Dist_Store_Open
    (char * storename, char * mode)

/* Open  storename  in  mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   DistStore  fp;
   int  retry;

   fp = openDistStore (storename, mode);
   for  (retry = 0;  fp == NULLSTOREHANDLE && retry < 3;  retry ++)
     {
      sleep (10);
      fp = openDistStore (storename, mode);
     }
   if  (fp == NULLSTOREHANDLE)
       {
        fprintf (stderr, "ERROR %d:  Could not open diststore  %s \n",
                 errno, storename);
        perror (strerror (errno));

        exit (FILE_OPEN_FAILURE);
       }

   return  fp;
  }



static void  Dump_Screen_Info
    (int frag_id, Screen_Info_t * screen, char dir)

//  Show screen information (if any) in  screen  for frag number
//   frag_id .   dir  indicates whether the fragment screen info
//  has been reversed (because the fragment was reverse-complemented).


  {
   static FILE  * dump_file = NULL;
   int  i;

   if  (dir == 'x')
       {
        fclose (dump_file);
        return;
       }

   if  (dump_file == NULL)
       {
        dump_file = File_Open ("screendump.out", "w");
       }

//   if  (screen -> num_matches == 0)
//       return;

   fprintf (dump_file, "%7d%c: ", frag_id, dir);
   for  (i = 0;  i < screen -> num_matches;  i ++)
     {
      fprintf (dump_file, "  %3d,%-3d%c",
               screen -> range [i] . bgn,
               screen -> range [i] . end,
               screen -> range [i] . last ? '*' : ' ');
     }
   fprintf (dump_file, "\n");

   return;
  }


static Overlap_t  Extend_Alignment
    (Match_Node_t * Match, char * S, int S_Len, char * T, int T_Len,
     int * S_Lo, int * S_Hi, int * T_Lo, int * T_Hi, int * Errors,
     Work_Area_t * WA)

//  See how far the exact match in  Match  extends.  The match
//  refers to strings  S  and  T  with lengths  S_Len  and  T_Len ,
//  respectively.  Set  S_Lo ,  S_Hi ,  T_Lo  and  T_Hi  to the
//  regions within  S  and  T  to which the match extends.
//  Return the type of overlap:  NONE if doesn't extend to
//  the end of either fragment;  LEFT_BRANCH_PT
//  or  RIGHT_BRANCH_PT  if match extends to the end of just one fragment;
//  or  DOVETAIL  if it extends to the end of both fragments, i.e.,
//  it is a complete overlap.
//  Set  Errors  to the number of errors in the alignment if it is
//  a  DOVETAIL  overlap.

  {
   Overlap_t  return_type;
   int  S_Left_Begin, S_Right_Begin, S_Right_Len;
   int  T_Left_Begin, T_Right_Begin, T_Right_Len;
   int  Error_Limit, Left_Errors, Right_Errors, Total_Olap;
   int  i, Leftover, Right_Match_To_End, Left_Match_To_End;

   S_Left_Begin = Match -> Start - 1;
   S_Right_Begin = Match -> Start + Match -> Len;
   S_Right_Len = S_Len - S_Right_Begin;
   T_Left_Begin = Match -> Offset - 1;
   T_Right_Begin = Match -> Offset + Match -> Len;
   T_Right_Len = T_Len - T_Right_Begin;
   Total_Olap = OVL_Min_int (Match -> Start, Match -> Offset)
                     + OVL_Min_int (S_Right_Len, T_Right_Len)
                     + Match -> Len;
   Error_Limit = WA -> Error_Bound [Total_Olap];

   if  (S_Right_Len == 0 || T_Right_Len == 0)
       {
        Right_Errors = 0;
        WA -> Right_Delta_Len = 0;
        (* S_Hi) = (* T_Hi) = 0;
        Right_Match_To_End = TRUE;
       }
   else if  (S_Right_Len <= T_Right_Len)
       Right_Errors = Prefix_Edit_Dist (S + S_Right_Begin, S_Right_Len,
                          T + T_Right_Begin, T_Right_Len, Error_Limit,
                          S_Hi, T_Hi, & Right_Match_To_End, WA);
     else
       {
        Right_Errors = Prefix_Edit_Dist (T + T_Right_Begin, T_Right_Len,
                          S + S_Right_Begin, S_Right_Len, Error_Limit,
                          T_Hi, S_Hi, & Right_Match_To_End, WA);
        for  (i = 0;  i < WA -> Right_Delta_Len;  i ++)
          WA -> Right_Delta [i] *= -1;
       }

   (* S_Hi) += S_Right_Begin - 1;
   (* T_Hi) += T_Right_Begin - 1;

   assert (Right_Errors <= Error_Limit);

   if  (S_Left_Begin < 0 || T_Left_Begin < 0)
       {
        Left_Errors = 0;
        WA -> Left_Delta_Len = 0;
        (* S_Lo) = (* T_Lo) = 0;
        Leftover = 0;
        Left_Match_To_End = TRUE;
       }
   else if  (S_Right_Begin <= T_Right_Begin)
       Left_Errors = Rev_Prefix_Edit_Dist (S + S_Left_Begin,
                        S_Left_Begin + 1, T + T_Left_Begin,
                        T_Left_Begin + 1,
                        Error_Limit - Right_Errors,
                        S_Lo, T_Lo, & Leftover, & Left_Match_To_End,
                        WA);
     else
       {
        Left_Errors = Rev_Prefix_Edit_Dist (T + T_Left_Begin,
                        T_Left_Begin + 1, S + S_Left_Begin,
                        S_Left_Begin + 1,
                        Error_Limit - Right_Errors,
                        T_Lo, S_Lo, & Leftover, & Left_Match_To_End,
                        WA);
        for  (i = 0;  i < WA -> Left_Delta_Len;  i ++)
          WA -> Left_Delta [i] *= -1;
       }

   (* S_Lo) += S_Left_Begin + 1;        // Check later for branch points
   (* T_Lo) += T_Left_Begin + 1;        // Check later for branch points

   if  (! Right_Match_To_End)
       {
        if  (! Doing_Partial_Overlaps)
            WA -> Left_Delta_Len = 0;
        if  (! Left_Match_To_End)
            return_type = NONE;
          else
            return_type = RIGHT_BRANCH_PT;
       }
     else
       {
        if  (! Left_Match_To_End)
            return_type = LEFT_BRANCH_PT;
          else
            return_type = DOVETAIL;
       }

   if  (return_type == DOVETAIL || Doing_Partial_Overlaps)
       {
        (* Errors) = Left_Errors + Right_Errors;
        assert ((* Errors) <= Error_Limit);

        if  (WA -> Right_Delta_Len > 0)
            {
             if  (WA -> Right_Delta [0] > 0)
                 WA -> Left_Delta [WA -> Left_Delta_Len ++]
                      = WA -> Right_Delta [0] + Leftover + Match -> Len;
               else
                 WA -> Left_Delta [WA -> Left_Delta_Len ++]
                      = WA -> Right_Delta [0] - Leftover - Match -> Len;
            }
        for  (i = 1;  i < WA -> Right_Delta_Len;  i ++)
          WA -> Left_Delta [WA -> Left_Delta_Len ++] = WA -> Right_Delta [i];
       }

   return  return_type;
  }



#if  USE_SOURCE_FIELD
static void  Find_End_Annotations
    (char * s, int frag_len, char * left, char * right, char * other)

//  Find annotations in celsim source field  s  and put character code
//  of one on the left end of fragment in  (* left)  and one on the right
//  end of the fragment in  (* right) .  If there is a non-end annotation
//  put it into  (* Other) .  If any is missing, set it to a blank
//  character.  Modifies string  s .   frag_len  is the length of the
//  fragment.

  {
   char  * p;

   fprintf (stderr, "source = \"%s\"\n", s);

   (* left) = (* right) = (* other) = ' ';
   p = strtok (s, "\n\r");        // frag id tag
   if  (p == NULL)
       return;
   p = strtok (NULL, "\n\r");     // frag coordinates in genome
   if  (p == NULL)
       return;
   while  ((p = strtok (NULL, "\n\r")) != NULL)
     {
      int  n, repeat_start, repeat_stop, frag_start, frag_stop;
      char  tag [MAX_SOURCE_LENGTH];

      n = sscanf (p, "%s [%d,%d] [%d,%d]", tag, & repeat_start, & repeat_stop,
                  & frag_start, & frag_stop);
      if  (n != 5)
          break;
      if  (frag_start == 0 || frag_stop == 0)
          (* left) = tag [0];
      if  (frag_start == frag_len || frag_stop == frag_len)
          (* right) = tag [0];
      if  (frag_start != 0 && frag_stop != 0
               && frag_start != frag_len && frag_stop != frag_len)
          (* other) = tag [0];
     }

   return;
  }
#endif



static void  Find_Overlaps
    (char Frag [], int Frag_Len, char quality [],
     Int_Frag_ID_t Frag_Num, Direction_t Dir,
     Work_Area_t * WA)

//  Find and output all overlaps and branch points between string
//   Frag  and any fragment currently in the global hash table.
//   Frag_Len  is the length of  Frag  and  Frag_Num  is its ID number.
//   Dir  is the orientation of  Frag .

  {
   String_Ref_t  Ref;
   char  * P, * Window;
   uint64  Key, Next_Key;
   int64  Sub, Next_Sub, Where;
   Check_Vector_t  This_Check, Next_Check;
   int  Offset, Shift, Next_Shift;
   int  hi_hits;
   int  screen_sub, screen_lo, screen_hi;
   int  j;

   memset (WA -> String_Olap_Space, 0, STRING_OLAP_MODULUS * sizeof (String_Olap_t));
   WA -> Next_Avail_String_Olap = STRING_OLAP_MODULUS;
   WA -> Next_Avail_Match_Node = 1;

   assert (Frag_Len >= WINDOW_SIZE);

   Offset = 0;
   P = Window = Frag;
   screen_sub = 0;
   if  (WA -> screen_info . num_matches == 0)
       screen_lo = screen_hi = INT_MAX;
     else
       {
        screen_lo = WA -> screen_info . range [0] . bgn;
        screen_hi = WA -> screen_info . range [0] . end;
       }

   Key = 0;
   for  (j = 0;  j < WINDOW_SIZE;  j ++)
     Key |= (uint64) (Bit_Equivalent [(int) * (P ++)]) << (2 * j);

   Sub = HASH_FUNCTION (Key);
   Shift = HASH_CHECK_FUNCTION (Key);
   Next_Key = (Key >> 2);
   Next_Key |= ((uint64)
                 (Bit_Equivalent [(int) * P])) << (2 * (WINDOW_SIZE - 1));
   Next_Sub = HASH_FUNCTION (Next_Key);
   Next_Shift = HASH_CHECK_FUNCTION (Next_Key);
   Next_Check = Hash_Check_Array [Next_Sub];

   if  ((Hash_Check_Array [Sub] & (((Check_Vector_t) 1) << Shift)) != 0
          && Offset <= screen_lo - WINDOW_SIZE + WINDOW_SCREEN_OLAP)
       {
        Ref = Hash_Find (Key, Sub, Window, & Where, & hi_hits);
        if  (hi_hits)
            WA -> screen_info . left_end_screened = TRUE;
        if  (! Ref . Empty)
          {
           while  (TRUE)
             {
#if  SHOW_STATS
Match_Ct ++;
#endif
              if  (Contig_Mode
                     || Frag_Num < (int) Ref . String_Num + Hash_String_Num_Offset)
                  Add_Ref  (Ref, Offset, WA);

              if  (Ref . Last)
                  break;
                else
                  {
                   Ref = Extra_Ref_Space [++ Where];
                   assert (! Ref . Empty);
                  }
             }
          }
       }

   while  ((* P) != '\0')
     {
      if  (Offset == screen_hi - 1 - WINDOW_SCREEN_OLAP)
          {
           if  (WA -> screen_info . range [screen_sub] . last)
               screen_lo = screen_hi = INT_MAX;
             else
               {
                screen_sub ++;
                screen_lo = WA -> screen_info . range [screen_sub] . bgn;
                screen_hi = WA -> screen_info . range [screen_sub] . end;
               }
          }

      Window ++;
      Offset ++;

      Key = Next_Key;
      Shift = Next_Shift;
      Sub = Next_Sub;
      This_Check = Next_Check;
      P ++;
      Next_Key = (Key >> 2);
      Next_Key |= ((uint64)
                    (Bit_Equivalent [(int) * P])) << (2 * (WINDOW_SIZE - 1));
      Next_Sub = HASH_FUNCTION (Next_Key);
      Next_Shift = HASH_CHECK_FUNCTION (Next_Key);
      Next_Check = Hash_Check_Array [Next_Sub];

      if  ((This_Check & (((Check_Vector_t) 1) << Shift)) != 0
             && Offset <= screen_lo - WINDOW_SIZE + WINDOW_SCREEN_OLAP)
          {
           Ref = Hash_Find (Key, Sub, Window, & Where, & hi_hits);
           if  (hi_hits)
               {
                if  (Offset < HOPELESS_MATCH)
                    WA -> screen_info . left_end_screened = TRUE;
                if  (Frag_Len - Offset - WINDOW_SIZE + 1 < HOPELESS_MATCH)
                    WA -> screen_info . right_end_screened = TRUE;
               }
           if  (! Ref . Empty)
             {
              while  (TRUE)
                {
#if  SHOW_STATS
Match_Ct ++;
#endif
                 if  (Contig_Mode
                        || Frag_Num < (int) Ref . String_Num + Hash_String_Num_Offset)
                     Add_Ref  (Ref, Offset, WA);

                 if  (Ref . Last)
                     break;
                   else
                     {
                      Ref = Extra_Ref_Space [++ Where];
                      assert (! Ref . Empty);
                     }
                }
             }
          }
     }

#if  SHOW_STATS
{
 int32  i, Ct = 0;

 for  (i = 0;  i < WA -> Next_Avail_String_Olap;  i ++)
   if  (WA -> String_Olap_Space [i] . Full)
       Ct ++;

 Incr_Distrib (& String_Olap_Dist, Ct);
}
#endif

#if  0
fprintf (stderr, "Frag %6d %c has %6d matches\n",
         Frag_Num, Dir == FORWARD ? 'f' : 'r', WA -> Next_Avail_Match_Node - 1);
#endif

   Process_String_Olaps  (Frag, Frag_Len, quality, Frag_Num, Dir, WA);

   return;
  }



static void  Flip_Screen_Range
    (Screen_Info_t * screen, int len)

//  Change the  range  entries in  (* screen)  to reflect that
//  the current fragment has been reverse-complemented.
//  len  is the length of the fragment.  Also reverse the
//  order of the entries.

  {
   int  i, j, k;
   Screen_Range_t  save;

   screen -> left_end_screened = FALSE;
   screen -> right_end_screened = FALSE;

   if  (screen -> num_matches < 1)
       return;

   for  (i = 0, j = screen -> num_matches - 1;  i < j;  i ++, j --)
     {
      save = screen -> range [j];
      screen -> range [j] . bgn = len - screen -> range [i] . end;
      screen -> range [j] . end = len - screen -> range [i] . bgn;
      screen -> range [i] . bgn = len - save . end;
      screen -> range [i] . end = len - save . bgn;
     }

   if  (i == j)
       {
        k = screen -> range [i] . end;
        screen -> range [i] . end = len - screen -> range [i] . bgn;
        screen -> range [i] . bgn = len - k;
       }

   screen -> left_end_screened
       = screen -> range [0] . bgn < HOPELESS_MATCH;
   screen -> right_end_screened
       = screen -> range [screen -> num_matches - 1] . end
                > len - HOPELESS_MATCH;

   return;
  }



FragStore  Frag_Store_Open
    (char * storename, char * mode)

/* Open  storename  in  mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   FragStore  fp;
   int  retry, max_sleep = 10;

   fp = openFragStore (storename, mode);
   for  (retry = 0;  fp == NULLSTOREHANDLE && retry < 5;  retry ++)
     {
      int  sleep_time;

      sleep_time = 1 + rand () % max_sleep;
      sleep (sleep_time);
      max_sleep *= 2;

      fprintf (stderr, "ERROR %d:  Could not open fragstore  %s  Retrying ... \n",
                 errno, storename);
      fp = openFragStore (storename, mode);
     }
   if  (fp == NULLSTOREHANDLE)
       {
        fprintf (stderr, "ERROR %d:  Could not open fragstore  %s \n",
                 errno, storename);
        perror (strerror (errno));

        exit (FILE_OPEN_FAILURE);
       }

   return  fp;
  }



static void  Get_Range
    (char * s, char ch, int * lo, int * hi)

//  Extract the string range specified in string  s  into
//  (* lo)  and  (* hi) .   ch  is the option letter associated
//  with this range.

  {
   char  * p;

   p = strchr (s, '-');
   if  (p == NULL)
       {
        fprintf (stderr, "ERROR:  No hyphen in  -%c  range \"%s\"\n",
                 ch, s);
        exit (-1);
       }
   if  (p == s)
       (* lo) = 0;
   else if  (sscanf (optarg, "%d", lo) != 1)
       {
        fprintf (stderr, "ERROR:  Bad first number in  -%c  range \"%s\"\n",
                 ch, s);
        exit (-1);
       }
   p ++;
   if  ((* p) == '\0')
       (* hi) = INT_MAX;
   else if  (sscanf (p, "%d", hi) != 1)
       {
        fprintf (stderr, "ERROR:  Bad first number in  -%c  range \"%s\"\n",
                 ch, s);
        exit (-1);
       }

   if  ((* lo) > (* hi))
       {
        fprintf (stderr, "ERROR:  Numbers reversed in  -%c  range \"%s\"\n",
                 ch, s);
        exit (-1);
       }

   return;
  }


static void  Get_Relevant_Screen_Matches
    (Screen_Info_t * screen)

//  Move relevant matches in  (* screen)  from  match  array
//  to  range  array and coalesce the overlapping ones.

  {
   int  i, j;

#if  0
fprintf (stderr, "Num matches = %d\n", screen -> num_matches);
for  (i = 0;  i < screen -> num_matches;  i ++)
  fprintf (stderr, "  %4d %4d %s\n",
           screen -> match [i] . where . bgn,
           screen -> match [i] . where . end,
           screen -> match [i] . relevance & AS_OVL_HEED_RPT ?
              "hard" : "soft");
#endif

   for  (i = j = 0;  i < screen -> num_matches;  i ++)
     if  (screen -> match [i] . relevance & AS_OVL_HEED_RPT)
         {
          screen -> range [j] . bgn = screen -> match [i] . where . bgn;
          screen -> range [j] . end = screen -> match [i] . where . end;
          screen -> range [j] . last = FALSE;
          j ++;
         }

   Coalesce_Screen_Info (screen -> range, 0, & j);

   screen -> num_matches = j;
   if  (j > 0)
       screen -> range [j - 1] . last = TRUE;

   return;
  }




static String_Ref_t  Hash_Find
    (uint64 Key, int64 Sub, char * S, int64 * Where, int * hi_hits)

//  Search for string  S  with hash key  Key  in the global
//  Hash_Table  starting at subscript  Sub .  Return the matching
//  reference in the hash table if there is one, or else a reference
//  with the  Empty bit set true.  Set  (* Where)  to the subscript in
//  Extra_Ref_Space  where the reference was found if it was found there.
//  Set  (* hi_hits)  to  TRUE  if hash table entry is found and has
//  at least  Hi_Hit_Limit  references; otherwise, set to  FALSE .

  {
   String_Ref_t  H_Ref;
   char  * T;
   unsigned char  Key_Check;
   int64  Ct, Probe;
   int  i;

#if  SHOW_STATS
Hash_Find_Ct ++;
#endif
   Key_Check = KEY_CHECK_FUNCTION (Key);
   Probe = PROBE_FUNCTION (Key);

   (* hi_hits) = FALSE;
   Ct = 0;
   do
     {
      for  (i = 0;  i < Hash_Table [Sub] . Entry_Ct;  i ++)
        if  (Hash_Table [Sub] . Check [i] == Key_Check)
            {
             int  is_empty;

             H_Ref = Hash_Table [Sub] . Entry [i];
             is_empty = H_Ref . Empty;
             if  (! H_Ref . Last && ! is_empty)
                 {
                  (* Where) = (H_Ref . String_Num << OFFSET_BITS) + H_Ref . Offset;
                  H_Ref = Extra_Ref_Space [(* Where)];
                 }
             T = Data + String_Start [H_Ref . String_Num] + H_Ref . Offset;
             if  (strncmp (S, T, WINDOW_SIZE) == 0)
                 {
#if  ANALYZE_HITS && ! DO_KMER_HITS_PROFILE
                  if  (Hash_Table [Sub] . Hits [i] < 255)
                      Hash_Table [Sub] . Hits [i] ++;
#endif
                  if  (Hash_Table [Sub] . Hits [i] >= Hi_Hit_Limit
                         || is_empty)
                      {
                       H_Ref . Empty = TRUE;
                       (* hi_hits) = TRUE;
                      }
                  return  H_Ref;
                 }
            }
      if  (Hash_Table [Sub] . Entry_Ct < ENTRIES_PER_BUCKET)
          {
           H_Ref . Empty = TRUE;
           return  H_Ref;
          }
#if  SHOW_STATS
Collision_Ct ++;
#endif
      Sub = (Sub + Probe) % HASH_TABLE_SIZE;
     }  while  (++ Ct < HASH_TABLE_SIZE);

   H_Ref . Empty = TRUE;
   return  H_Ref;
  }



static int  Hash_Hits
    (uint64 key, char * s)

//  Search for string  s  with hash key  key  in the global
//  Hash_Table  and, if found, return the number of hits
//  recorded for it.  If not found, return  0 .

  {
#if  ANALYZE_HITS
   unsigned char  key_check;
   int64  ct, probe, sub;

   sub = HASH_FUNCTION (key);
   key_check = KEY_CHECK_FUNCTION (key);
   probe = PROBE_FUNCTION (key);

   ct = 0;
   do
     {
      int  i;

      for  (i = 0;  i < Hash_Table [sub] . Entry_Ct;  i ++)
        if  (Hash_Table [sub] . Check [i] == key_check)
            {
             String_Ref_t  h_ref;
             char  * t;

             h_ref = Hash_Table [sub] . Entry [i];

             if  (! h_ref . Last && ! h_ref . Empty)
                 h_ref = Extra_Ref_Space
                             [(h_ref . String_Num << OFFSET_BITS) + h_ref . Offset];

             t = Data + String_Start [h_ref . String_Num] + h_ref . Offset;
             if  (strncmp (s, t, WINDOW_SIZE) == 0)
                 return  Hash_Table [sub] . Hits [i];
            }

      if  (Hash_Table [sub] . Entry_Ct < ENTRIES_PER_BUCKET)
          return  0;

      sub = (sub + probe) % HASH_TABLE_SIZE;
     }  while  (++ ct < HASH_TABLE_SIZE);
#endif

   return  0;
  }



static int  Has_Bad_Window
    (char a [], int n, int window_len, int threshold)

//  Return  TRUE  iff there is any length-( window_len ) subarray
//  of  a [0 .. (n-1)]  that sums to  threshold  or higher.

  {
   int  i, j, sum;

   if  (n < window_len)
       return  FALSE;

   sum = 0;
   for  (i = 0;  i < window_len;  i ++)
     sum += a [i];

   j = 0;
   if  (sum >= threshold)
       {
#if  SHOW_BAD_WINDOWS
        printf (">>> Window [%d .. %d]  sum = %d\n",
                j, i - 1, sum);
#endif
        return  TRUE;
       }

   while  (i < n)
     {
      sum -= a [j ++];
      sum += a [i ++];
      if  (sum >= threshold)
          {
#if  SHOW_BAD_WINDOWS
           printf (">>> Window [%d .. %d]  sum = %d\n",
                   j, i - 1, sum);
#endif
           return  TRUE;
          }
     }

   return  FALSE;
  }



static void  Hash_Insert
    (String_Ref_t Ref, uint64 Key, char * S)

//  Insert  Ref  with hash key  Key  into global  Hash_Table .
//  Ref  represents string  S .

  {
   String_Ref_t  H_Ref;
   char  * T;
   int  Shift;
   unsigned char  Key_Check;
   int64  Ct, Probe, Sub;
   int  i;

   Sub = HASH_FUNCTION (Key);
   Shift = HASH_CHECK_FUNCTION (Key);
   Hash_Check_Array [Sub] |= (((Check_Vector_t) 1) << Shift);
   Key_Check = KEY_CHECK_FUNCTION (Key);
   Probe = PROBE_FUNCTION (Key);

   Ct = 0;
   do
     {
      for  (i = 0;  i < Hash_Table [Sub] . Entry_Ct;  i ++)
        if  (Hash_Table [Sub] . Check [i] == Key_Check)
            {
             H_Ref = Hash_Table [Sub] . Entry [i];
             T = Data + String_Start [H_Ref . String_Num] + H_Ref . Offset;
             if  (strncmp (S, T, WINDOW_SIZE) == 0)
                 {
                  if  (H_Ref . Last)
                      Extra_Ref_Ct ++;
                  Next_Ref [(String_Start [Ref . String_Num] + Ref . Offset)
                              / (HASH_KMER_SKIP + 1)]
                                = H_Ref;
                  Extra_Ref_Ct ++;
                  Ref . Last = FALSE;
                  Hash_Table [Sub] . Entry [i] = Ref;

#if  ANALYZE_HITS && DO_KMER_HITS_PROFILE
                  Hash_Table [Sub] . Hits [i] ++;
#else
                  if  (Hash_Table [Sub] . Hits [i] < HIGHEST_KMER_LIMIT)
                      Hash_Table [Sub] . Hits [i] ++;
#endif

                  return;
                 }
            }
if  (i != Hash_Table [Sub] . Entry_Ct)
    {
     fprintf (stderr, "i = %d  Sub = " F_S64 "  Entry_Ct = %d\n",
                  i, Sub, Hash_Table [Sub] . Entry_Ct);
    }
      assert (i == Hash_Table [Sub] . Entry_Ct);
      if  (Hash_Table [Sub] . Entry_Ct < ENTRIES_PER_BUCKET)
          {
           Ref . Last = TRUE;
           Hash_Table [Sub] . Entry [i] = Ref;
           Hash_Table [Sub] . Check [i] = Key_Check;
           Hash_Table [Sub] . Entry_Ct ++;
           Hash_Entries ++;
           Hash_Table [Sub] . Hits [i] = 1;
           return;
          }
#if  SHOW_STATS
Collision_Ct ++;
#endif
      Sub = (Sub + Probe) % HASH_TABLE_SIZE;
     }  while  (++ Ct < HASH_TABLE_SIZE);

   fprintf (stderr, "ERROR:  Hash table full\n");
   assert (FALSE);
   
   return;
  }



static void  Hash_Mark_Empty
    (uint64 key, char * s)

//  Set the  empty  bit to true for the hash table entry
//  corresponding to string  s  whose hash key is  key .
//  Also set global  String_Info . left/right_end_screened
//  true if the entry occurs near the left/right end, resp.,
//  of the string in the hash table.

  {
   String_Ref_t  h_ref;
   char  * t;
   unsigned char  key_check;
   int64  ct, probe;
   int64  sub;
   int  i;

   sub = HASH_FUNCTION (key);
   key_check = KEY_CHECK_FUNCTION (key);
   probe = PROBE_FUNCTION (key);

   ct = 0;
   do
     {
      for  (i = 0;  i < Hash_Table [sub] . Entry_Ct;  i ++)
        if  (Hash_Table [sub] . Check [i] == key_check)
            {
             h_ref = Hash_Table [sub] . Entry [i];
             t = Data + String_Start [h_ref . String_Num] + h_ref . Offset;
             if  (strncmp (s, t, WINDOW_SIZE) == 0)
                 {
                  Mark_Screened_Ends_Chain (Hash_Table [sub] . Entry [i]);
                  Hash_Table [sub] . Entry [i] . Empty = TRUE;
                  return;
                 }
            }
      assert (i == Hash_Table [sub] . Entry_Ct);
      if  (Hash_Table [sub] . Entry_Ct < ENTRIES_PER_BUCKET)
          return;     // Not found
      sub = (sub + probe) % HASH_TABLE_SIZE;
     }  while  (++ ct < HASH_TABLE_SIZE);

   fprintf (stderr, "ERROR:  Hash table full\n");
   assert (FALSE);
   
   return;
  }



static void  Initialize_Globals
    (void)

//  Set initial values of global variables.

  {
   int  e, i, Start;

   assert (MAX_ERROR_RATE >= AS_READ_ERROR_RATE
             && MAX_ERROR_RATE >= AS_GUIDE_ERROR_RATE);

   for  (i = 0;  i <= ERRORS_FOR_FREE;  i ++)
     Read_Edit_Match_Limit [i] = 0;

   Start = 1;
   for  (e = ERRORS_FOR_FREE + 1;  e < MAX_ERRORS;  e ++)
     {
      Start = Binomial_Bound (e - ERRORS_FOR_FREE, AS_READ_ERROR_RATE,
                  Start, EDIT_DIST_PROB_BOUND);
      Read_Edit_Match_Limit [e] = Start - 1;
      assert (Read_Edit_Match_Limit [e] >= Read_Edit_Match_Limit [e - 1]);
     }

   for  (i = 0;  i <= ERRORS_FOR_FREE;  i ++)
     Guide_Edit_Match_Limit [i] = 0;

   Start = 1;
   for  (e = ERRORS_FOR_FREE + 1;  e < MAX_ERRORS;  e ++)
     {
      Start = Binomial_Bound (e - ERRORS_FOR_FREE, AS_GUIDE_ERROR_RATE,
                  Start, EDIT_DIST_PROB_BOUND);
      Guide_Edit_Match_Limit [e] = Start - 1;
      assert (Guide_Edit_Match_Limit [e] >= Guide_Edit_Match_Limit [e - 1]);
     }

   for  (i = 0;  i <= AS_READ_MAX_LEN;  i ++)
     Read_Error_Bound [i] = (int) (i * AS_READ_ERROR_RATE + 0.0000000000001);

   for  (i = 0;  i <= AS_READ_MAX_LEN;  i ++)
     Guide_Error_Bound [i] = (int) (i * AS_GUIDE_ERROR_RATE + 0.0000000000001);

   for  (i = 0;  i <= AS_READ_MAX_LEN;  i ++)
     Branch_Cost [i] = i * Branch_Match_Value + Branch_Error_Value;

   Bit_Equivalent ['a'] = Bit_Equivalent ['A'] = 0;
   Bit_Equivalent ['c'] = Bit_Equivalent ['C'] = 1;
   Bit_Equivalent ['g'] = Bit_Equivalent ['G'] = 2;
   Bit_Equivalent ['t'] = Bit_Equivalent ['T'] = 3;

   for  (i = 0;  i < 256;  i ++)
     {
      char  ch = tolower ((char) i);

      if  (ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
          Char_Is_Bad [i] = 0;
        else
          Char_Is_Bad [i] = 1;
     }

   Screen_Space_Size = Max_Hash_Strings;
   Screen_Space
       = (Screen_Range_t *) Safe_malloc
             (Screen_Space_Size * sizeof (Screen_Range_t));

   Hash_Table = (Hash_Bucket_t *) Safe_malloc
                    (HASH_TABLE_SIZE * sizeof (Hash_Bucket_t));
   fprintf (stderr, "### Bytes in hash table = " F_SIZE_T "\n",
            HASH_TABLE_SIZE * sizeof (Hash_Bucket_t));
   
   Hash_Check_Array = (Check_Vector_t *) Safe_malloc
                          (HASH_TABLE_SIZE * sizeof (Check_Vector_t));
   Loc_ID = (uint64 *) Safe_calloc (Max_Hash_Strings, sizeof (uint64));
   String_Info = (Hash_Frag_Info_t *) Safe_calloc (Max_Hash_Strings,
                     sizeof (Hash_Frag_Info_t));
   String_Start = (int64 *) Safe_calloc (Max_Hash_Strings, sizeof (int64));
   Kind_Of_Frag = (FragType *) Safe_calloc (Max_Hash_Strings, sizeof (FragType));
   Screen_Sub = (int *) Safe_calloc (Max_Hash_Strings, sizeof (int));

#if  USE_SOURCE_FIELD
   Left_Repeat_Tag = (char *) Safe_malloc (Max_Hash_Strings);
   Right_Repeat_Tag = (char *) Safe_malloc (Max_Hash_Strings);
   Other_Repeat_Tag = (char *) Safe_malloc (Max_Hash_Strings);
#endif

#if  SHOW_STATS
 Init_Distrib (& Olap_Len_Dist, 14);
 {
   int  i;
   for  (i = 0;  i < 14;  i ++)
     Olap_Len_Dist . Thold [i] = 50.0 * (i + 1);
 }
 Init_Distrib (& Num_Olaps_Dist, 25);
 {
   int  i;
   for  (i = 0;  i <= 10;  i ++)
     Num_Olaps_Dist . Thold [i] = 1.0 * i;
   for  (i = 11;  i <= 18;  i ++)
     Num_Olaps_Dist . Thold [i] = 5.0 * (i - 8);
   for  (i = 19;  i <= 23;  i ++)
     Num_Olaps_Dist . Thold [i] = 100.0 * (i - 18);
   Num_Olaps_Dist . Thold [24] = 1000.0;
 }
 Init_Distrib (& Kmer_Freq_Dist, 25);
 {
   int  i;
   for  (i = 0;  i <= 10;  i ++)
     Kmer_Freq_Dist . Thold [i] = 1.0 * i;
   for  (i = 11;  i <= 18;  i ++)
     Kmer_Freq_Dist . Thold [i] = 5.0 * (i - 8);
   for  (i = 19;  i <= 23;  i ++)
     Kmer_Freq_Dist . Thold [i] = 100.0 * (i - 18);
   Kmer_Freq_Dist . Thold [24] = 1000.0;
 }
 Init_Distrib (& Kmer_Hits_Dist, 21);
 {
   int  i;
   for  (i = 0;  i <= 20;  i ++)
     Kmer_Hits_Dist . Thold [i] = 100.0 * i;
 }
#if  ANALYZE_HITS
 Init_Distrib (& Hash_Hits_Dist, 32);
 {
   int  i;
   for  (i = 0;  i <= 10;  i ++)
     Hash_Hits_Dist . Thold [i] = 10.0 * i;
   for  (i = 11;  i <= 20;  i ++)
     Hash_Hits_Dist . Thold [i] = 100.0 * (i - 9);
   for  (i = 21;  i <= 31;  i ++)
     Hash_Hits_Dist . Thold [i] = 1000.0 * (i - 19);
 }
#endif
 Init_Distrib (& String_Olap_Dist, 28);
 {
   int  i;
   for  (i = 0;  i <= 10;  i ++)
     String_Olap_Dist . Thold [i] = 1.0 * i;
   for  (i = 11;  i <= 18;  i ++)
     String_Olap_Dist . Thold [i] = 5.0 * (i - 8);
   for  (i = 19;  i <= 27;  i ++)
     String_Olap_Dist . Thold [i] = 50.0 * (i - 17);
 }
 Init_Distrib (& Hits_Per_Olap_Dist, 28);
 {
   int  i;
   for  (i = 0;  i <= 10;  i ++)
     Hits_Per_Olap_Dist . Thold [i] = 1.0 * i;
   for  (i = 11;  i <= 18;  i ++)
     Hits_Per_Olap_Dist . Thold [i] = 5.0 * (i - 8);
   for  (i = 19;  i <= 27;  i ++)
     Hits_Per_Olap_Dist . Thold [i] = 50.0 * (i - 17);
 }
Init_Distrib (& Diag_Dist, 21);
{
 int  i;
 for  (i = 0;  i <= 20;  i ++)
   Diag_Dist . Thold [i] = 1.0 * i;
}
Init_Distrib (& Gap_Dist, 21);
{
 int  i;
 for  (i = 0;  i <= 20;  i ++)
   Gap_Dist . Thold [i] = 1.0 * i;
}
Init_Distrib (& Exacts_Per_Olap_Dist, 21);
{
 int  i;
 for  (i = 0;  i <= 20;  i ++)
   Exacts_Per_Olap_Dist . Thold [i] = 1.0 * i;
}
Init_Distrib (& Edit_Depth_Dist, 26);
{
 int  i;
 for  (i = 0;  i <= 25;  i ++)
   Edit_Depth_Dist . Thold [i] = 1.0 * i;
}
#endif

#if  SHOW_PROGRESS
Table_Ct = 1;
Olap_Ct = 0;
Start_Time = clock ();
#endif
   return;
  }



void  Initialize_Work_Area
    (Work_Area_t * WA, int id)

//  Allocate memory for  (* WA)  and set initial values.
//  Set  thread_id  field to  id .

  {
   int  i, Offset, Del;

   WA -> String_Olap_Size = INIT_STRING_OLAP_SIZE;
   WA -> String_Olap_Space = (String_Olap_t *) Safe_malloc
           (WA -> String_Olap_Size * sizeof (String_Olap_t));
   WA -> Match_Node_Size = INIT_MATCH_NODE_SIZE;
   WA -> Match_Node_Space = (Match_Node_t *) Safe_malloc
           (WA -> Match_Node_Size * sizeof (Match_Node_t));

   Offset = 2;
   Del = 6;
   for  (i = 0;  i < MAX_ERRORS;  i ++)
     {
       WA -> Edit_Array [i] = WA -> Edit_Space + Offset;
       Offset += Del;
       Del += 2;
     }

   WA -> status = 0;
   WA -> myRead = new_ReadStruct ();

   WA -> screen_info . match = (IntScreenMatch *) Safe_malloc
                                   (INIT_SCREEN_MATCHES * sizeof (IntScreenMatch));
   WA -> screen_info . range = (Screen_Range_t *) Safe_malloc
                                   (INIT_SCREEN_MATCHES * sizeof (Screen_Range_t));
   WA -> screen_info . match_len = INIT_SCREEN_MATCHES;
   WA -> screen_info . num_matches = 0;

   WA -> thread_id = id;

   return;
  }



static int  Interval_Intersection
    (int a, int b, int c, int d)

//  Return the number of bases by which the closed interval  [a, b]
//  intersects the closed interval  [c, d] .

  {
   if  (d < a || b < c)
       return  0;

   return  1 + OVL_Min_int (b, d) - OVL_Max_int (a, c);
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



int  OVL_Max_int
    (int a, int b)

//  Return the larger of  a  and  b .

  {
   if  (a < b)
       return  b;

   return  a;
  }



int  OVL_Min_int
    (int a, int b)

//  Return the smaller of  a  and  b .

  {
   if  (a < b)
       return  a;

   return  b;
  }



static void  Mail_Error_Message
    (char * person, int num, char * store)

//  E-mail an error message to  person  about error number
//  num  on store  store .

  {
   FILE  * message_file = NULL;
   char  filename [1000];

   sprintf (filename, "/tmp/lsfovlerr.%d.msg", getpid ());
   message_file = fopen (filename, "w");
   if  (message_file != NULL)
       {
        char  buff [1000];

        fprintf (message_file, 
                 "Overlap error %d  store %s\n"
                 "%s\n",
                 num, store, strerror (num));
        fclose (message_file);
        sprintf (buff, "mail %s < %s",
                 person, filename);
        system (buff);
        sprintf (buff, "rm -f %s", filename);
        system (buff);
       }

   return;
  }



static void  Mark_Screened_Ends_Chain
    (String_Ref_t ref)

//  Mark  left/right_end_screened in global  String_Info  for
//   ref  and everything in its list, if they occur near
//  enough to the end of the string.

  {
   Mark_Screened_Ends_Single (ref);

   while  (! ref . Last)
       {
        ref = Next_Ref [(String_Start [ref . String_Num] + ref . Offset)
                          / (HASH_KMER_SKIP + 1)];
        Mark_Screened_Ends_Single (ref);
       }

   return;
  }



static void  Mark_Screened_Ends_Single
    (String_Ref_t ref)

//  Mark  left/right_end_screened in global  String_Info  for
//  single reference  ref , if it occurs near
//  enough to the end of its string.

  {
   int  s_num, len;

   s_num = ref . String_Num;
   len = String_Info [s_num] . length;

   if  (ref . Offset < HOPELESS_MATCH)
       String_Info [s_num] . left_end_screened = TRUE;
   if  (len - ref . Offset - WINDOW_SIZE + 1 < HOPELESS_MATCH)
       String_Info [s_num] . right_end_screened = TRUE;

   return;
  }



static void  Mark_Skip_Kmers
    (void)

//  Set  Empty  bit true for all entries in global  Hash_Table
//  that match a kmer in file  Kmer_Skip_File .

  {
   uint64  key;
   char  line [MAX_LINE_LEN];
   int  ct = 0;

   rewind (Kmer_Skip_File);

   while  (fgets (line, MAX_LINE_LEN, Kmer_Skip_File) != NULL)
     {
      int  i, len;

      ct ++;
      len = strlen (line) - 1;
      if  (line [0] != '>' || line [len] != '\n')
          {
           fprintf (stderr, "ERROR:  Bad line %d in kmer skip file\n", ct);
           fputs (line, stderr);
           exit (EXIT_FAILURE);
          }

      if  (fgets (line, MAX_LINE_LEN, Kmer_Skip_File) == NULL)
          {
           fprintf (stderr, "ERROR:  Bad line after %d in kmer skip file\n", ct);
           exit (EXIT_FAILURE);
          }
      ct ++;
      len = strlen (line) - 1;
      if  (len != WINDOW_SIZE || line [len] != '\n')
          {
           fprintf (stderr, "ERROR:  Bad line %d in kmer skip file\n", ct);
           fputs (line, stderr);
           exit (EXIT_FAILURE);
          }

      key = 0;
      for  (i = 0;  i < len;  i ++)
        {
         line [i] = tolower (line [i]);
         key |= (uint64) (Bit_Equivalent [(int) line [i]]) << (2 * i);
        }
      Hash_Mark_Empty (key, line);

      Rev_Complement (line, len);
      key = 0;
      for  (i = 0;  i < len;  i ++)
        key |= (uint64) (Bit_Equivalent [(int) line [i]]) << (2 * i);
      Hash_Mark_Empty (key, line);
     }

   fprintf (stderr, "Read %d kmers to mark to skip\n", ct / 2);

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



static int64  Next_Odd_Prime
    (int64 N)

//  Return the first odd prime  >= N .

  {
   int64  Div, Last;

   if  (N % 2 == 0)
       N ++;
   while  (TRUE)
     {
      Last = (int64) (sqrt ((double) N));
      for  (Div = 3;  Div <= Last;  Div += 2)
        if  (N % Div == 0)
            break;
      if  (Div > Last)
          return  N;
      N += 2;
     }

   return  0;     // Just to make the compiler happy
  }



#if  ANALYZE_HITS
#define  MIN_HI_HIT_KMERS     200
#define  HIT_THRESHOLD        100

void  Output_High_Hit_Frags
    (void)

//  From the fragments used to build the hash table print those that
//  contain at least  MIN_HI_HIT_KMERS  kmers that were hit at least
//   HIT_THRESHOLD  times.

  {
   int  i, j;

   if  (High_Hits_File == NULL)
       {
        High_Hits_File = File_Open ("highhits.out", "w");
        fprintf (High_Hits_File, "Fragments with >= %d kmers hit >= %d times\n",
                 MIN_HI_HIT_KMERS, HIT_THRESHOLD);
       }
   for  (i = 0;  i < String_Ct;  i ++)
     {
      char  * p, * window;
      int  curr_sub, offset, screen_lo, screen_hi;
      int  high_hits = 0;
      uint64  key, key_is_bad;
      int  frag_hits [AS_READ_MAX_LEN + 1];

      if  (String_Info [i] . length < WINDOW_SIZE)
          continue;

      curr_sub = Screen_Sub [i];
      if  (curr_sub == 0)
          screen_lo = screen_hi = INT_MAX;
        else
          {
           screen_lo = Screen_Space [curr_sub] . bgn;
           screen_hi = Screen_Space [curr_sub] . end;
          }

      p = window = Data + String_Start [i];
      key = key_is_bad = 0;
      offset = 0;
      for  (j = 0;  j < WINDOW_SIZE;  j ++)
        {
         key_is_bad |= (uint64) (Char_Is_Bad [* p]) << j;
         key |= (uint64) (Bit_Equivalent [* (p ++)]) << (2 * j);
        }

      while  ((* p) != '\0')
        {
         if  (offset <= screen_lo - WINDOW_SIZE + WINDOW_SCREEN_OLAP
                  && ! key_is_bad)
             frag_hits [offset] = Hash_Hits (key, window);
           else
             frag_hits [offset] = -1;

         if  (frag_hits [offset] >= HIT_THRESHOLD)
             high_hits ++;

         if  (offset == screen_hi - 1 - WINDOW_SCREEN_OLAP)
             {
              if  (Screen_Space [curr_sub] . last)
                  screen_lo = screen_hi = INT_MAX;
                else
                  {
                   curr_sub ++;
                   screen_lo = Screen_Space [curr_sub] . bgn;
                   screen_hi = Screen_Space [curr_sub] . end;
                  }
             }

         window ++;
         offset ++;

         key_is_bad >>= 1;
         key_is_bad |= (uint64) (Char_Is_Bad [* p]) << (WINDOW_SIZE - 1);
         key >>= 2;
         key |= (uint64)
                    (Bit_Equivalent [* (p ++)]) << (2 * (WINDOW_SIZE - 1));
        }

      if  (high_hits >= MIN_HI_HIT_KMERS)
          {
           int  ct = 0, lo = -1, hi;
           long int  sum;
           double  avg;

           for  (j = 0;  j < String_Info [i] . length;  j ++)
             if  (frag_hits [j] >= HIT_THRESHOLD)
                 {
                  if  (lo == -1)
                      lo = j;
                  hi = j + WINDOW_SIZE - 1;
                 }
           sum = 0;
           for  (j = lo;  j < hi - WINDOW_SIZE;  j ++)
             sum += frag_hits [j];
           avg = (double) sum / (1 + hi - lo);
           fprintf (High_Hits_File,
                    ">Frag #%d  %d .. %d  hi_ct_positions = %d  avg_hits = %.1f\n",
                    i + Hash_String_Num_Offset, lo, hi, high_hits, avg);
           ct = 0;
           for  (j = lo;  j <= hi;  j ++)
             {
              fputc (Data [String_Start [i] + j], High_Hits_File);
              ct ++;
              if  (ct % DISPLAY_WIDTH == 0 || j == hi)
                  fputc ('\n', High_Hits_File);
             }
          }
     }

   return;
  }
#endif



static void  Output_Overlap
    (Int_Frag_ID_t S_ID, int S_Len, Direction_t S_Dir,
     Int_Frag_ID_t T_ID, int T_Len, Olap_Info_t * olap)

//  Output the overlap between strings  S_ID  and  T_ID  which
//  have lengths  S_Len  and  T_Len , respectively.
//  The overlap information is in  (* olap) .
//  S_Dir  indicates the orientation of  S .

  {
   int  S_Right_Hang, T_Right_Hang;
   int  this_diag;
   OverlapMesg ovMesg;
   GenericMesg outputMesg;
   signed char deltas[2 * AS_READ_MAX_LEN];
   signed char *deltaCursor = deltas;
   outputMesg.m = &ovMesg; 
   outputMesg.t = MESG_OVL;

#if  SHOW_PROGRESS
Olap_Ct ++;
#endif
#if  SHOW_STATS
Overlap_Ct ++;
#endif

if  (Contig_Mode)
    {
     int  start, stop;

     if  (olap -> s_lo > 0)
         {
          olap -> t_lo -= olap -> s_lo;
          Dovetail_Overlap_Ct ++;
         }
     else if  (olap -> s_hi < S_Len - 1)
         {
          olap -> t_hi += S_Len - 1 - olap -> s_hi;
          Dovetail_Overlap_Ct ++;
         }
       else
         Contained_Overlap_Ct ++;

     if  (S_Dir == FORWARD)
         {
          start = olap -> t_lo;
          stop = olap -> t_hi + 1;
         }
       else
         {
          start = olap -> t_hi + 1;
          stop = olap -> t_lo;
         }

#if  1
     fprintf (BOL_File, "%11" F_U64P " %8d %8d %7d %7d\n",
             Loc_ID [T_ID - Hash_String_Num_Offset], T_ID, S_ID, start, stop);
#else
//  Show lengths of fragments, too.  Used to determine which overlaps
//  are A-end dovetails, B-end dovetails, contained or spanning.
     fprintf (BOL_File, "%11" F_U64P " %8d %8d %7d %7d %7d %7d\n",
             Loc_ID [T_ID - Hash_String_Num_Offset], T_ID, S_ID,
             start, stop, T_Len, S_Len);
#endif

     Total_Overlaps ++;
     return;
    }

   ovMesg.delta = deltas;
   *deltas = '\0';
   assert (S_ID < T_ID);

   S_Right_Hang = S_Len - olap -> s_hi - 1;
   T_Right_Hang = T_Len - olap -> t_hi - 1;
   this_diag = olap -> t_lo - olap -> s_lo;

   if  (olap -> s_lo > olap -> t_lo
            || (olap -> s_lo == olap -> t_lo
                   && S_Right_Hang > T_Right_Hang))
       {       // S is on the left
        ovMesg.aifrag = (Int_Frag_ID_t) S_ID;
        ovMesg.bifrag = (Int_Frag_ID_t) T_ID;

        if  (S_Dir == FORWARD)
          ovMesg.orientation = AS_NORMAL;
        else
          ovMesg.orientation = AS_OUTTIE;

        if  (S_Right_Hang >= T_Right_Hang)
            ovMesg.overlap_type = AS_CONTAINMENT;
        else
            ovMesg.overlap_type =  AS_DOVETAIL;

        ovMesg.ahg = olap -> s_lo;
        ovMesg.bhg = T_Right_Hang - S_Right_Hang;
        ovMesg.quality = olap -> quality;
        ovMesg.min_offset = olap -> s_lo - (olap -> max_diag - this_diag);
        ovMesg.max_offset = olap -> s_lo + (this_diag - olap -> min_diag);
#if  USE_SOURCE_FIELD
if  (ovMesg.min_offset + 3 < ovMesg.ahg)
    fprintf (Source_Log_File, "%4d %4d %4d %c\n", ovMesg.min_offset,
             ovMesg.ahg, ovMesg.max_offset,
             Left_Repeat_Tag [T_ID - Hash_String_Num_Offset]);
#endif
#if  FOR_CARL_FOSLER
        {
         double  a_percen, b_percen;
         int  a_start, a_end, b_start, b_end;

         a_percen = S_Len;
         if  (ovMesg . ahg > 0)
             {
              a_percen -= ovMesg . ahg;
              a_start = ovMesg . ahg;
              b_start = 0;
             }
           else
             {
              a_start = 0;
              b_start = - ovMesg . ahg;
             }
         if  (ovMesg . bhg < 0)
             {
              a_percen += ovMesg . bhg;
              a_end = S_Len + ovMesg . bhg;
              b_end = T_Len;
             }
           else
             {
              a_end = S_Len;
              b_end = T_Len - ovMesg . bhg;
             }
         a_percen *= 100.0 / S_Len;

         if  (S_Dir != FORWARD)
             {
              int  save = a_start;

              a_start = S_Len - a_end;
              a_end = S_Len - save;
             }
         
         b_percen = T_Len;
         if  (ovMesg . ahg < 0)
             b_percen += ovMesg . ahg;
         if  (ovMesg . bhg > 0)
             b_percen -= ovMesg . bhg;
         b_percen *= 100.0 / T_Len;
         
         fprintf (Out_Stream,
                  "%9d  %9d  %c  %5d  %5d  %5.1f  %5d  %5d  %5.1f  %5.3f\n",
                  S_ID, T_ID, (char) ovMesg . orientation,
                  a_start, a_end, a_percen,
                  b_start, b_end, b_percen,
                  100.0 * ovMesg . quality);
         return;
        }
#endif
       }
     else
       {
#if  OUTPUT_OVERLAP_DELTAS
        for  (i = 0;  i < olap -> delta_ct;  i ++)
          olap -> delta [i] *= -1;
#endif

         ovMesg.bifrag = (Int_Frag_ID_t) S_ID;
         ovMesg.aifrag = (Int_Frag_ID_t) T_ID;
        if  (S_Dir == FORWARD)
          ovMesg.orientation = AS_NORMAL;
          else
            ovMesg.orientation = AS_INNIE;

        if  (T_Right_Hang >= S_Right_Hang)
            ovMesg.overlap_type = AS_CONTAINMENT;
        else
            ovMesg.overlap_type =  AS_DOVETAIL;
        ovMesg.ahg = olap -> t_lo;
        ovMesg.bhg = S_Right_Hang - T_Right_Hang;
        ovMesg.quality =  olap -> quality;
        ovMesg.min_offset = olap -> t_lo - (this_diag - olap -> min_diag);
        ovMesg.max_offset = olap -> t_lo + (olap -> max_diag - this_diag);
#if  USE_SOURCE_FIELD
if  (ovMesg.min_offset + 3 < ovMesg.ahg)
    fprintf (Source_Log_File, "%4d %4d %4d %c\n", ovMesg.min_offset,
             ovMesg.ahg, ovMesg.max_offset,
             Right_Repeat_Tag [T_ID - Hash_String_Num_Offset]);
#endif
#if  FOR_CARL_FOSLER
        {
         double  a_percen, b_percen;
         int  a_start, a_end, b_start, b_end;

         a_percen = T_Len;
         if  (ovMesg . ahg > 0)
             {
              a_percen -= ovMesg . ahg;
              a_start = ovMesg . ahg;
              b_start = 0;
             }
           else
             {
              a_start = 0;
              b_start = - ovMesg . ahg;
             }
         if  (ovMesg . bhg < 0)
             {
              a_percen += ovMesg . bhg;
              a_end = T_Len + ovMesg . bhg;
              b_end = S_Len;
             }
           else
             {
              a_end = T_Len;
              b_end = S_Len - ovMesg . bhg;
             }
         a_percen *= 100.0 / T_Len;
         
         if  (S_Dir != FORWARD)
             {
              int  save = b_start;

              b_start = S_Len - b_end;
              b_end = S_Len - save;
             }
         
         b_percen = S_Len;
         if  (ovMesg . ahg < 0)
             b_percen += ovMesg . ahg;
         if  (ovMesg . bhg > 0)
             b_percen -= ovMesg . bhg;
         b_percen *= 100.0 / S_Len;
         
#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_lock (& Write_Proto_Mutex);
#endif
         fprintf (Out_Stream,
                  "%9d  %9d  %c  %5d  %5d  %5.1f  %5d  %5d  %5.1f  %5.3f\n",
                  T_ID, S_ID, (char) ovMesg . orientation,
                  a_start, a_end, a_percen,
                  b_start, b_end, b_percen,
                  100.0 * ovMesg . quality);
#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_unlock (& Write_Proto_Mutex);
#endif
         return;
        }
#endif
       }

   assert (ovMesg.min_offset <= ovMesg.ahg && ovMesg.ahg <= ovMesg.max_offset);
   ovMesg.polymorph_ct = 0;

#if  OUTPUT_OVERLAP_DELTAS
   for  (i = 0;  i < olap -> delta_ct;  i ++)
     {
      int  j;

      for  (j = abs (olap -> delta [i]);  j > 0;  j -= AS_LONGEST_DELTA)
        {
         if  (j > AS_LONGEST_DELTA)
           *deltaCursor++ = AS_LONG_DELTA_CODE;
         else{
           *deltaCursor++ = j * Sign (olap -> delta [i]);
         }
        }
     }
#endif
   *deltaCursor = AS_ENDOF_DELTA_CODE;


#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_lock (& Write_Proto_Mutex);
#endif

   if((ovMesg.overlap_type == AS_CONTAINMENT) &&
      (ovMesg.orientation == AS_OUTTIE)) {
     // CMM: Regularize the reverse orientated containment overlaps to
     // a common orientation.
     const int ahg = ovMesg.ahg;
     const int bhg = ovMesg.bhg;
     const int min_delta = ovMesg.min_offset - ahg;
     const int max_delta = ovMesg.max_offset - ahg;

     ovMesg.orientation = AS_INNIE;
     ovMesg.ahg = -bhg;
     ovMesg.bhg = -ahg;
     ovMesg.min_offset = ovMesg.ahg - max_delta;
     ovMesg.max_offset = ovMesg.ahg - min_delta;
   }

   Total_Overlaps ++;
   if  (ovMesg . bhg <= 0)
       Contained_Overlap_Ct ++;
     else
       Dovetail_Overlap_Ct ++;
   Write_Msg_Fn (Out_Stream, & outputMesg);

#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_unlock (& Write_Proto_Mutex);
#endif

   return;
  }



static void  Output_Partial_Overlap
    (Int_Frag_ID_t s_id, Int_Frag_ID_t t_id, Direction_t dir,
     const Olap_Info_t * p, int s_len, int t_len)

//  Output the overlap between strings  s_id  and  t_id  which
//  have lengths  s_len  and  t_len , respectively.
//  The overlap information is in  p .
//  dir  indicates the orientation of the S string.

  {
   int  a, b, c, d;
   char  dir_ch;

   Total_Overlaps ++;

   // Convert to canonical form with s forward and use space-based
   // coordinates
   if  (dir == FORWARD)
       {
        a = p -> s_lo;
        b = p -> s_hi + 1;
        c = p -> t_lo;
        d = p -> t_hi + 1;
        dir_ch = 'f';
       }
     else
       {
        a = s_len - p -> s_hi - 1;
        b = s_len - p -> s_lo;
        c = p -> t_hi + 1;
        d = p -> t_lo;
        dir_ch = 'r';
       }
#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_lock (& Write_Proto_Mutex);
#endif
   fprintf (Out_Stream, "%7d %7d  %c %4d %4d %4d  %4d %4d %4d  %5.2f\n",
        s_id, t_id, dir_ch, a, b, s_len, c, d, t_len,
        100.0 * p -> quality);
#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_unlock (& Write_Proto_Mutex);
#endif

   return;
  }



static int  Passes_Screen
    (int lo, int hi, int sub)

//  Return TRUE iff the range  lo .. hi  contains enough bases
//  outside of the screen ranges in  Screen_Space [sub .. ]
//  to be a valid overlap.

  {
   int  ct;

   ct = Interval_Intersection (lo, hi, Screen_Space [sub] . bgn,
                               Screen_Space [sub] . end - 1);

   while  (! Screen_Space [sub] . last)
     {
      sub ++;
      ct += Interval_Intersection
                (lo, hi, Screen_Space [sub] . bgn,
                 Screen_Space [sub] . end - 1);
     }

   return (1 + hi - lo - ct >= MIN_OLAP_OUTSIDE_SCREEN);
  }



static int  Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Match_To_End, Work_Area_t * WA)

//  Return the minimum number of changes (inserts, deletes, replacements)
//  needed to match string  A [0 .. (m-1)]  with a prefix of string
//   T [0 .. (n-1)]  if it's not more than  Error_Limit .
//  If no match, return the number of errors for the best match
//  up to a branch point.
//  Put delta description of alignment in  WA -> Right_Delta  and set
//  WA -> Right_Delta_Len  to the number of entries there if it's a complete
//  match.
//  Set  A_End  and  T_End  to the rightmost positions where the
//  alignment ended in  A  and  T , respectively.
//  Set  Match_To_End  true if the match extended to the end
//  of at least one string; otherwise, set it false to indicate
//  a branch point.

  {
   int  Delta_Stack [MAX_ERRORS];
#if 0
   double  Cost_Sum, Min_Cost_Sum;
   int Adjustment, Cost_Sub;
#endif
   double  Score, Max_Score;
   int  Max_Score_Len = 0, Max_Score_Best_d = 0, Max_Score_Best_e = 0;
   int  Tail_Len;
   int  Best_d, Best_e, From, Last, Longest, Max, Row;
   int  Left, Right;
   int  d, e, i, j, k;

#if  SHOW_STATS
Edit_Dist_Ct ++;
#endif
   assert (m <= n);
   Best_d = Best_e = Longest = 0;
   WA -> Right_Delta_Len = 0;

   for  (Row = 0;  Row < m
                     && (A [Row] == T [Row]
                           || A [Row] == DONT_KNOW_CHAR
                           || T [Row] == DONT_KNOW_CHAR);  Row ++)
     ;

   WA -> Edit_Array [0] [0] = Row;

   if  (Row == m)                              // Exact match
       {
        (* A_End) = (* T_End) = m;
#if  SHOW_STATS
Incr_Distrib (& Edit_Depth_Dist, 0);
#endif
        (* Match_To_End) = TRUE;
        return  0;
       }

   Left = Right = 0;
   Max_Score = 0.0;
   for  (e = 1;  e <= Error_Limit;  e ++)
     {
      Left = OVL_Max_int (Left - 1, -e);
      Right = OVL_Min_int (Right + 1, e);
      WA -> Edit_Array [e - 1] [Left] = -2;
      WA -> Edit_Array [e - 1] [Left - 1] = -2;
      WA -> Edit_Array [e - 1] [Right] = -2;
      WA -> Edit_Array [e - 1] [Right + 1] = -2;

      for  (d = Left;  d <= Right;  d ++)
        {
         Row = 1 + WA -> Edit_Array [e - 1] [d];
         if  ((j = WA -> Edit_Array [e - 1] [d - 1]) > Row)
             Row = j;
         if  ((j = 1 + WA -> Edit_Array [e - 1] [d + 1]) > Row)
             Row = j;
         while  (Row < m && Row + d < n
                  && (A [Row] == T [Row + d]
                        || A [Row] == DONT_KNOW_CHAR
                        || T [Row + d] == DONT_KNOW_CHAR))
           Row ++;

         WA -> Edit_Array [e] [d] = Row;

         if  (Row == m || Row + d == n)
             {
#if  SHOW_STATS
Incr_Distrib (& Edit_Depth_Dist, e);
#endif
              //  Check for branch point here caused by uneven
              //  distribution of errors
              Score = Row * Branch_Match_Value - e;
                        // Assumes  Branch_Match_Value
                        //             - Branch_Error_Value == 1.0
              Tail_Len = Row - Max_Score_Len;
              if  ((Doing_Partial_Overlaps && Score < Max_Score)
                     ||  (e > MIN_BRANCH_END_DIST / 2
                           && Tail_Len >= MIN_BRANCH_END_DIST
                           && (Max_Score - Score) / Tail_Len >= MIN_BRANCH_TAIL_SLOPE))
                  {
                   (* A_End) = Max_Score_Len;
                   (* T_End) = Max_Score_Len + Max_Score_Best_d;
                   Set_Right_Delta (Max_Score_Best_e, Max_Score_Best_d, WA);
                   (* Match_To_End) = FALSE;
                   return  Max_Score_Best_e;
                  }

              // Force last error to be mismatch rather than insertion
              if  (Row == m
                     && 1 + WA -> Edit_Array [e - 1] [d + 1]
                          == WA -> Edit_Array [e] [d]
                     && d < Right)
                  {
                   d ++;
                   WA -> Edit_Array [e] [d] = WA -> Edit_Array [e] [d - 1];
                  }

              (* A_End) = Row;           // One past last align position
              (* T_End) = Row + d;
              Set_Right_Delta (e, d, WA);
              (* Match_To_End) = TRUE;
              return  e;
             }
        }

      while  (Left <= Right && Left < 0
                  && WA -> Edit_Array [e] [Left] < WA -> Edit_Match_Limit [e])
        Left ++;
      if  (Left >= 0)
          while  (Left <= Right
                    && WA -> Edit_Array [e] [Left] + Left < WA -> Edit_Match_Limit [e])
            Left ++;
      if  (Left > Right)
          break;
      while  (Right > 0
                  && WA -> Edit_Array [e] [Right] + Right < WA -> Edit_Match_Limit [e])
        Right --;
      if  (Right <= 0)
          while  (WA -> Edit_Array [e] [Right] < WA -> Edit_Match_Limit [e])
            Right --;
      assert (Left <= Right);

      for  (d = Left;  d <= Right;  d ++)
        if  (WA -> Edit_Array [e] [d] > Longest)
            {
             Best_d = d;
             Best_e = e;
             Longest = WA -> Edit_Array [e] [d];
            }

      Score = Longest * Branch_Match_Value - e;
               // Assumes  Branch_Match_Value - Branch_Error_Value == 1.0
      if  (Score > Max_Score)
          {
           Max_Score = Score;
           Max_Score_Len = Longest;
           Max_Score_Best_d = Best_d;
           Max_Score_Best_e = Best_e;
          }
     }

#if  SHOW_STATS
Incr_Distrib (& Edit_Depth_Dist, 1 + Error_Limit);
#endif

   (* A_End) = Max_Score_Len;
   (* T_End) = Max_Score_Len + Max_Score_Best_d;
   Set_Right_Delta (Max_Score_Best_e, Max_Score_Best_d, WA);
   (* Match_To_End) = FALSE;
   return  Max_Score_Best_e;
  }



static void  Process_Matches
    (int * Start, char * S, int S_Len, char * S_quality,
     Int_Frag_ID_t S_ID, Direction_t Dir, char * T, Hash_Frag_Info_t t_info,
     char * T_quality, Int_Frag_ID_t T_ID,
     Work_Area_t * WA, int consistent)

//  Find and report all overlaps and branch points between string  S
//  (with length  S_Len  and id  S_ID ) and string  T  (with
//  length & screen info in  t_info  and id  T_ID ) using the exact
//  matches in the list beginning at subscript  (* Start) .   Dir  is
//  the orientation of  S .

  {
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
   int  Errors, rejected;

   Done_S_Left = Done_S_Right = FALSE;
   t_len = t_info . length;

#if  SHOW_STATS
Is_Duplicate_Olap = FALSE;
{
 int  Space [2 * AS_READ_MAX_LEN + 3] = {0};
 int  * Ct = Space + AS_READ_MAX_LEN + 1;
 int  i, Diag_Ct, Gap_Ct, In_Order, Prev;
 
 In_Order = TRUE;
 Prev = INT_MAX;
 Gap_Ct = 0;
 for  (i = (* Start);  i != 0;  i = WA -> Match_Node_Space [i] . Next)
   {
    Gap_Ct ++;
    if  (WA -> Match_Node_Space [i] . Start >= Prev)
        In_Order = FALSE;
      else
        Prev = WA -> Match_Node_Space [i] . Start;
   }
 Incr_Distrib (& Exacts_Per_Olap_Dist, Gap_Ct);
 if  (In_Order)
     Incr_Distrib (& Gap_Dist, Gap_Ct - 1);
 Diag_Ct = 0;
 for  (i = - AS_READ_MAX_LEN;  i <= AS_READ_MAX_LEN;  i ++)
   if  (Ct [i] > 0)
       Diag_Ct ++;
 Incr_Distrib (& Diag_Dist, Diag_Ct);
}
#endif

   assert ((* Start) != 0);

       // If a singleton match is hopeless on either side
       // it needn't be processed
   if  (WA -> Match_Node_Space [(* Start)] . Next == 0
          && ! Doing_Partial_Overlaps)
       {
        int  s_head, t_head, s_tail, t_tail;
        int  is_hopeless = FALSE;

        s_head = WA -> Match_Node_Space [(* Start)] . Start;
        t_head = WA -> Match_Node_Space [(* Start)] . Offset;
        if  (s_head <= t_head)
            {
             if  (s_head > HOPELESS_MATCH
                    && ! WA -> screen_info . left_end_screened)
                 is_hopeless = TRUE;
            }
          else
            {
             if  (t_head > HOPELESS_MATCH
                    && ! t_info . left_end_screened)
                 is_hopeless = TRUE;
            }

        s_tail = S_Len - s_head - WA -> Match_Node_Space [(* Start)] . Len + 1;
        t_tail = t_len - t_head - WA -> Match_Node_Space [(* Start)] . Len + 1;
        if  (s_tail <= t_tail)
            {
             if  (s_tail > HOPELESS_MATCH
                    && ! WA -> screen_info . right_end_screened)
                 is_hopeless = TRUE;
            }
          else
            {
             if  (t_tail > HOPELESS_MATCH
                    && ! t_info . right_end_screened)
                 is_hopeless = TRUE;
            }

        if  (is_hopeless)
            {
             (* Start) = 0;
             Kmer_Hits_Without_Olap_Ct ++;
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

           if  (Kind_Of_Olap == DOVETAIL || Doing_Partial_Overlaps)
               {
                if  (1 + S_Hi - S_Lo >= MIN_OLAP_LEN
                       && 1 + T_Hi - T_Lo >= MIN_OLAP_LEN)
                    {
                     int  scr_sub;

                     Olap_Len = 1 + OVL_Min_int (S_Hi - S_Lo, T_Hi - T_Lo);
                     Quality = (double) Errors / Olap_Len;
                     scr_sub = Screen_Sub [T_ID - Hash_String_Num_Offset];
                     if  (Errors <= WA -> Error_Bound [Olap_Len]
                            && (scr_sub == 0
                                  || Passes_Screen (T_Lo, T_Hi, scr_sub)))
                         {
                          Add_Overlap (S_Lo, S_Hi, T_Lo, T_Hi, Quality,
                                       distinct_olap, & distinct_olap_ct, WA);
#if  SHOW_STATS
Incr_Distrib (& Olap_Len_Dist, 1 + S_Hi - S_Lo);
if  (Is_Duplicate_Olap)
    Duplicate_Olap_Ct ++;
  else
    {
     Regular_Olap_Ct ++;
     Is_Duplicate_Olap = TRUE;
    }
#endif
                         }
                    }
#if  SHOW_STATS
                  else
                    Too_Short_Ct ++;
#endif
               }
          }

#if  1
      if  (consistent)      
          (* Start) = 0;
#endif

      for  (Ref = Start;  (* Ref) != 0;  )
        {
         Ptr = WA -> Match_Node_Space + (* Ref);
         if  (Ptr == Longest_Match
                 || (Kind_Of_Olap == DOVETAIL
                     && S_Lo - SHIFT_SLACK <= Ptr -> Start
                     && Ptr -> Start + Ptr -> Len
                                     <= S_Hi + SHIFT_SLACK - 1
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
   
   if  (distinct_olap_ct > 0)
       {
        int  deleted [MAX_DISTINCT_OLAPS] = {0};
        Olap_Info_t  * p;
        int  i;

//  Check if any previously distinct overlaps should be merged because
//  of other merges.

        if  (Unique_Olap_Per_Pair)
            Combine_Into_One_Olap (distinct_olap, distinct_olap_ct, deleted);
          else
            Merge_Intersecting_Olaps (distinct_olap, distinct_olap_ct, deleted);

        p = distinct_olap;
        for  (i = 0;  i < distinct_olap_ct;  i ++)
          {
           if  (! deleted [i])
               {
#if  SHOW_OVERLAPS
    {
#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_lock (& Write_Proto_Mutex);
#endif

#if  DO_OLAP_ALIGN_PROFILE
Align_P = Align_Ct [T_ID - Hash_String_Num_Offset];
#else
     printf ("\nA = %d  B = %d  quality = %.4f  A_Len = %d  B_Len = %d\n",
             S_ID, T_ID, p -> quality, S_Len, t_len);
#endif

     Show_Overlap (S, S_Len, S_quality, T, t_len, T_quality, p);

#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_unlock (& Write_Proto_Mutex);
#endif
    }
#endif

#define  SHOW_V_ALIGN  0
                rejected = FALSE;
                if  (Use_Window_Filter)
                    {
                     int  d, i, j, k, q_len;
                     char  q_diff [MAX_FRAG_LEN];

                     i = p -> s_lo;
                     j = p -> t_lo;
                     q_len = 0;

                     for  (k = 0;  k < p -> delta_ct;  k ++)
                       {
                        int  n, len;

                        len = abs (p -> delta [k]);
                        for  (n = 1;  n < len;  n ++)
                          {
                           if  (S [i] == T [j]
                                  || S [i] == DONT_KNOW_CHAR
                                  || T [j] == DONT_KNOW_CHAR)
                               d = 0;
                             else
                               {
                                d = char_Min (S_quality [i], T_quality [j]);
                                d = char_Min (d, QUALITY_CUTOFF);
                               }
                    #if  SHOW_V_ALIGN
                           fprintf (stderr, "%3d: %c %2d  %3d: %c %2d  %3d\n",
                                    i, S [i], S_quality [i],
                                    j, T [j], T_quality [j], d);
                    #endif
                           q_diff [q_len ++] = d;
                           i ++;
                           j ++;
                          }
                        if  (p -> delta [k] > 0)
                            {
                             d = S_quality [i];
                    #if  SHOW_V_ALIGN
                             fprintf (stderr, "%3d: %c %2d  %3s: %c %2c  %3d\n",
                                      i, S [i], S_quality [i],
                                      " . ", '-', '-', d);
                    #endif
                             i ++;
                            }
                          else
                            {
                             d = T_quality [j];
                    #if  SHOW_V_ALIGN
                             fprintf (stderr, "%3s: %c %2c  %3d: %c %2d  %3d\n",
                                      " . ", '-', '-',
                                      j, T [j], T_quality [j], d);
                    #endif
                             j ++;
                            }
                        q_diff [q_len ++] = char_Min (d, QUALITY_CUTOFF);
                       }
                     while  (i <= p -> s_hi)
                       {
                        if  (S [i] == T [j]
                               || S [i] == DONT_KNOW_CHAR
                               || T [j] == DONT_KNOW_CHAR)
                            d = 0;
                          else
                            {
                             d = char_Min (S_quality [i], T_quality [j]);
                             d = char_Min (d, QUALITY_CUTOFF);
                            }
                    #if  SHOW_V_ALIGN
                        fprintf (stderr, "%3d: %c %2d  %3d: %c %2d  %3d\n",
                                 i, S [i], S_quality [i],
                                 j, T [j], T_quality [j], d);
                    #endif
                        q_diff [q_len ++] = d;
                        i ++;
                        j ++;
                       }

                     if  (Has_Bad_Window (q_diff, q_len,
                                          BAD_WINDOW_LEN, BAD_WINDOW_VALUE))
                         {
#if  SHOW_BAD_WINDOWS
                          printf ("\n>>>At least %d errors in window of length %d\n",
                                  BAD_WINDOW_VALUE / QUALITY_CUTOFF, BAD_WINDOW_LEN);
                          printf ("\nA = %d  B = %d  quality = %.4f\n",
                                  S_ID, T_ID, p -> quality);
                          Show_Overlap (S, S_Len, S_quality, T, t_len, T_quality, p);
#endif
                          rejected = TRUE;
                          Bad_Short_Window_Ct ++;
                         }
                     else if  (Has_Bad_Window (q_diff, q_len, 100, 240))
                         {
#if  SHOW_BAD_WINDOWS
                          printf ("\n###At least %d errors in window of length %d\n",
                                  240 / QUALITY_CUTOFF, 100);
                          printf ("\nA = %d  B = %d  quality = %.4f\n",
                                  S_ID, T_ID, p -> quality);
                          Show_Overlap (S, S_Len, S_quality, T, t_len, T_quality, p);
#endif
                          rejected = TRUE;
                          Bad_Long_Window_Ct ++;
                         }
                    }

#if  SHOW_SNPS
  {
   int  mismatch_ct, indel_ct;
   int  olap_bases, unscreened_olap_bases, all_match_ct, hi_qual_ct;
   int  local_msnp_bin [20] = {0};
   int  local_isnp_bin [20] = {0};
   int  local_match_bin [20] = {0};
   int  local_other_bin [20] = {0};
   int  i;

#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_lock (& Write_Proto_Mutex);
#endif

   Show_SNPs (S, S_Len, S_quality, T, t_len, T_quality, p, WA,
              & mismatch_ct, & indel_ct, & olap_bases,
              & unscreened_olap_bases, & all_match_ct, & hi_qual_ct,
              local_msnp_bin, local_isnp_bin,
              local_match_bin, local_other_bin);
   fprintf (Out_Stream, "%8d %c %8d %3d %3d %4d %4d %4d %4d",
            S_ID, Dir == FORWARD ? 'f' : 'r', T_ID,
            mismatch_ct, indel_ct, olap_bases, unscreened_olap_bases,
            all_match_ct, hi_qual_ct);
   fprintf (Out_Stream, "  ");
   for  (i = 3;  i <= 9;  i ++)
     fprintf (Out_Stream, " %3d", local_msnp_bin [i]);
   fprintf (Out_Stream, "  ");
   for  (i = 3;  i <= 9;  i ++)
     fprintf (Out_Stream, " %3d", local_isnp_bin [i]);
   fprintf (Out_Stream, "  ");
   for  (i = 3;  i <= 9;  i ++)
     fprintf (Out_Stream, " %3d", local_match_bin [i]);
   fprintf (Out_Stream, "  ");
   for  (i = 3;  i <= 9;  i ++)
     fprintf (Out_Stream, " %3d", local_other_bin [i]);
   fprintf (Out_Stream, "\n");

#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_unlock (& Write_Proto_Mutex);
#endif
   rejected = TRUE;
  }
#endif

                if  (! rejected)
                    {
                     if  (Doing_Partial_Overlaps || Single_Line_Output)
                         Output_Partial_Overlap (S_ID, T_ID, Dir, p,
                              S_Len, t_len);
                       else
                         Output_Overlap (S_ID, S_Len, Dir, T_ID, t_len, p);
                     overlaps_output ++;
                     if  (p -> s_lo == 0)
                         WA -> A_Olaps_For_Frag ++;
                     if  (p -> s_hi >= S_Len - 1)
                         WA -> B_Olaps_For_Frag ++;
                    }
               }
           p ++;
          }
       }
         
#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_lock (& Write_Proto_Mutex);
#endif

   if  (overlaps_output == 0)
       Kmer_Hits_Without_Olap_Ct ++;
     else
       {
        Kmer_Hits_With_Olap_Ct ++;
        if  (overlaps_output > 1)
            Multi_Overlap_Ct ++;
       }

#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_unlock (& Write_Proto_Mutex);
#endif

   return;
  }



void  Process_Overlaps
    (FragStreamHandle stream, Work_Area_t * WA)

//  Find and output all overlaps between strings in  stream 
//  and those in the global hash table.   (* WA)  has the
//  data structures that used to be global.

  {
   char  Frag [AS_READ_MAX_LEN + 1];
   char  quality [AS_READ_MAX_LEN + 1];
   Int_Frag_ID_t  Curr_String_Num, start_string_num;
   uint32  last_old_frag_read;
   int  frag_status;
   int  Len;

#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_lock (& FragStore_Mutex);
#endif

   Curr_String_Num = getStartIndexFragStream(stream);
   start_string_num = Curr_String_Num;

#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_unlock (& FragStore_Mutex);
#endif

   while  ((frag_status
              = Read_Next_Frag (Frag, quality, stream, WA -> myRead,
                                & (WA -> screen_info), & last_old_frag_read)))
     {

#if  SHOW_PROGRESS
if  (Curr_String_Num % 10000 == 0)
    {
     Stop_Time = clock ();
     fprintf (stderr, "Matching string %6d\n", Curr_String_Num);
     fprintf (stderr, 
              "%7.1f sec %7" F_S64P " olaps\n",
              (double) (Stop_Time - Start_Time) / CLOCKS_PER_SEC,
              Olap_Ct);
     Olap_Ct = 0;
     Start_Time = clock ();
    }
#endif

      if  (frag_status == DELETED_FRAG)
          {
           Curr_String_Num ++;
           continue;
          }

      getReadType_ReadStruct (WA -> myRead, & (WA -> curr_frag_type));
      Len = strlen (Frag);
      if  (Len < MIN_OLAP_LEN)
          fprintf (stderr, "Frag %d is too short.  Len = %d\n",
                      Curr_String_Num, Len);
        else
          {
#if DEBUGSTREAM
fprintf(stderr,"* Calling Find_Overlaps for string FORWARD%d olaps = %ld\n",
        Curr_String_Num, Olap_Ct);
#endif      

#if  SHOW_STATS
Overlap_Ct = 0;
Kmer_Hits_Ct = 0;
#endif


#ifdef  CHECK_SCREENING
Dump_Screen_Info (Curr_String_Num, & (WA -> screen_info), 'f');
#endif

           WA -> A_Olaps_For_Frag = WA -> B_Olaps_For_Frag = 0;

#if  0  // Special Mouse
{
 int  j;
 char  * p, * window;
 uint64  key;

 printf (">> Kmer frequencies in hash table for frag %d\n", Curr_String_Num);
 p = window = Frag;
 key = 0;

 for  (j = 0;  j < WINDOW_SIZE - 1;  j ++)
   {
    key |= (uint64) (Bit_Equivalent [* (p ++)]) << (2 * j);
    printf ("%c %6s\n", * (p - 1), "-");
   }

 while  ((* p) != '\0')
   {
    key |= (uint64) (Bit_Equivalent [* (p ++)]) << (2 * (WINDOW_SIZE - 1));
    printf ("%c %6d\n", * (p - 1), Hash_Hits (key, window ++));
    key >>= 2;
   }
}
#endif
           Find_Overlaps (Frag, Len, quality, Curr_String_Num, FORWARD, WA);

#if  FOR_CARL_FOSLER
           if  (WA -> A_Olaps_For_Frag >= Frag_Olap_Limit
                 || WA -> B_Olaps_For_Frag >= Frag_Olap_Limit)
               {
#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_lock (& Log_Msg_Mutex);
#endif
                if  (WA -> A_Olaps_For_Frag >= Frag_Olap_Limit)
                    fprintf (stderr,
                             "### Hit forward A-end olap limit for frag %d\n",
                         Curr_String_Num);
                if  (WA -> B_Olaps_For_Frag >= Frag_Olap_Limit)
                    fprintf (stderr,
                             "### Hit forward B-end olap limit for frag %d\n",
                         Curr_String_Num);
#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_unlock (& Log_Msg_Mutex);
#endif
               }
#endif


           Rev_Complement (Frag, Len);
           Reverse_String (quality, Len);
           Flip_Screen_Range (& (WA -> screen_info), Len);
#ifdef  CHECK_SCREENING
Dump_Screen_Info (Curr_String_Num, & (WA -> screen_info), 'r');
#endif

           WA -> A_Olaps_For_Frag = WA -> B_Olaps_For_Frag = 0;

           Find_Overlaps (Frag, Len, quality, Curr_String_Num, REVERSE, WA);

#if  FOR_CARL_FOSLER
           if  (WA -> A_Olaps_For_Frag >= Frag_Olap_Limit
                 || WA -> B_Olaps_For_Frag >= Frag_Olap_Limit)
               {
#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_lock (& Log_Msg_Mutex);
#endif
                if  (WA -> A_Olaps_For_Frag >= Frag_Olap_Limit)
                    fprintf (stderr,
                             "### Hit reverse A-end olap limit for frag %d\n",
                         Curr_String_Num);
                if  (WA -> B_Olaps_For_Frag >= Frag_Olap_Limit)
                    fprintf (stderr,
                             "### Hit reverse B-end olap limit for frag %d\n",
                         Curr_String_Num);
#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_unlock (& Log_Msg_Mutex);
#endif
               }
#endif

#if  SHOW_STATS
Incr_Distrib (& Num_Olaps_Dist, Overlap_Ct);
Incr_Distrib (& Kmer_Hits_Dist, Kmer_Hits_Ct);
#endif

          }

      Curr_String_Num ++;
     }

#if  SHOW_THREAD_PROGRESS
fprintf (stderr, "### Thread #%d processed overlaps for frags %ld .. %ld\n",
         WA -> thread_id, start_string_num, Curr_String_Num - 1);
#endif

   return;
  }



static int  Process_String_Olaps
    (char * S, int Len, char * S_quality, Int_Frag_ID_t ID, 
     Direction_t Dir, Work_Area_t * WA)

//  Find and report all overlaps and branch points between string  S
//  and all strings in global  String_Olap_Space .
//  Return the number of entries processed.
//  Len  is the length of  S ,  ID  is its fragment ID  and
//  Dir  indicates if  S  is forward, or reverse-complemented.

  {
   int32  i, ct, root_num, start, processed_ct;

//  Move all full entries to front of String_Olap_Space and set
//  diag_sum to average diagonal.  if enough entries to bother,
//  sort by average diagonal.  Then process positive & negative diagonals
//  separately in order from longest to shortest overlap.  Stop
//  processing when output limit has been reached.

   for  (i = ct = 0;  i < WA -> Next_Avail_String_Olap;  i ++)
     if  (WA -> String_Olap_Space [i] . Full)
         {
          root_num = WA -> String_Olap_Space [i] . String_Num;
          if  (Contig_Mode || root_num + Hash_String_Num_Offset > ID)
              {
#if  SHOW_STATS
Incr_Distrib (& Hits_Per_Olap_Dist, WA -> String_Olap_Space [i] . Kmer_Hits);
#endif
if  (WA -> String_Olap_Space [i] . Match_List == 0)
    {
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

   if  (ct <= Frag_Olap_Limit)
       {
        for  (i = 0;  i < ct;  i ++)
          {
           root_num = WA -> String_Olap_Space [i] . String_Num;
           if  ((WA -> curr_frag_type == AS_READ
                       || WA -> curr_frag_type == AS_EXTR
                       || WA -> curr_frag_type == AS_TRNR
                       || WA -> curr_frag_type == AS_FBAC)
                    && (Kind_Of_Frag [root_num] == AS_READ
                            || Kind_Of_Frag [root_num] == AS_EXTR
                            || Kind_Of_Frag [root_num] == AS_TRNR
                            || Kind_Of_Frag [root_num] == AS_FBAC))
               {
                WA -> Edit_Match_Limit = Read_Edit_Match_Limit;
                WA -> Error_Bound = Read_Error_Bound;
               }
             else
               {
                WA -> Edit_Match_Limit = Guide_Edit_Match_Limit;
                WA -> Error_Bound = Guide_Error_Bound;
               }
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

   for  (i = start;  i < ct && WA -> A_Olaps_For_Frag < Frag_Olap_Limit ;  i ++)
     {
      root_num = WA -> String_Olap_Space [i] . String_Num;
      if  ((WA -> curr_frag_type == AS_READ
                  || WA -> curr_frag_type == AS_EXTR
                  || WA -> curr_frag_type == AS_TRNR
                  || WA -> curr_frag_type == AS_FBAC)
               && (Kind_Of_Frag [root_num] == AS_READ
                       || Kind_Of_Frag [root_num] == AS_EXTR
                       || Kind_Of_Frag [root_num] == AS_TRNR
                       || Kind_Of_Frag [root_num] == AS_FBAC))
          {
           WA -> Edit_Match_Limit = Read_Edit_Match_Limit;
           WA -> Error_Bound = Read_Error_Bound;
          }
        else
          {
           WA -> Edit_Match_Limit = Guide_Edit_Match_Limit;
           WA -> Error_Bound = Guide_Error_Bound;
          }
      Process_Matches
             (& (WA -> String_Olap_Space [i] . Match_List),
              S, Len, S_quality, ID, Dir,
              Data + String_Start [root_num],
              String_Info [root_num],
              Quality_Data + String_Start [root_num],
              root_num + Hash_String_Num_Offset, WA,
              WA -> String_Olap_Space [i] . consistent);
      assert (WA -> String_Olap_Space [i] . Match_List == 0);

#if  FOR_CARL_FOSLER
{
 if  (ID == 517)
     fprintf (stderr, "### 517A  A_Olaps = %d  B_Olaps = %d\n",
              WA -> A_Olaps_For_Frag, WA -> B_Olaps_For_Frag);
}
#endif
      processed_ct ++;
     }
   for  (i = start - 1;  i >= 0 && WA -> B_Olaps_For_Frag < Frag_Olap_Limit ;  i --)
     {
      root_num = WA -> String_Olap_Space [i] . String_Num;
      if  ((WA -> curr_frag_type == AS_READ
                  || WA -> curr_frag_type == AS_EXTR
                  || WA -> curr_frag_type == AS_TRNR
                  || WA -> curr_frag_type == AS_FBAC)
               && (Kind_Of_Frag [root_num] == AS_READ
                       || Kind_Of_Frag [root_num] == AS_EXTR
                       || Kind_Of_Frag [root_num] == AS_TRNR
                       || Kind_Of_Frag [root_num] == AS_FBAC))
          {
           WA -> Edit_Match_Limit = Read_Edit_Match_Limit;
           WA -> Error_Bound = Read_Error_Bound;
          }
        else
          {
           WA -> Edit_Match_Limit = Guide_Edit_Match_Limit;
           WA -> Error_Bound = Guide_Error_Bound;
          }
      Process_Matches
             (& (WA -> String_Olap_Space [i] . Match_List),
              S, Len, S_quality, ID, Dir,
              Data + String_Start [root_num],
              String_Info [root_num],
              Quality_Data + String_Start [root_num],
              root_num + Hash_String_Num_Offset, WA,
              WA -> String_Olap_Space [i] . consistent);
      assert (WA -> String_Olap_Space [i] . Match_List == 0);

#if  FOR_CARL_FOSLER
{
 if  (ID == 517)
     fprintf (stderr, "### 517B  A_Olaps = %d  B_Olaps = %d\n",
              WA -> A_Olaps_For_Frag, WA -> B_Olaps_For_Frag);
}
#endif
      processed_ct ++;
     }

   return  processed_ct;
  }



#define  PROFILE_INFILE_NAME  dan5.inp
    //  Name of file of fragments to profile kmer hits on

void  Profile_Hits
    (void)

//  Show for each fragment in  PROFILE_INFILE_NAME  at each kmer
//  position the number of kmer occurrences in the hash table.
//  Does both forward and reverse complement sense of each fragment.
//  After that, separately print all kmers in the hash table
//  that occur at least  MIN_HITS  times.
//  Output goes to standard output.

  {
#if  ANALYZE_HITS && ! SHOW_HI_HIT_KMERS
// #if  0  // Special Mouse
   Input_Stream  profile_stream;
   MessageType  imesgtype;
   GenericMesg  * pmesg;
   InternalFragMesg  * ifg_mesg;
   char  frag [2 * AS_READ_MAX_LEN];
   char  rev_frag [2 * AS_READ_MAX_LEN];
   int  ct [2 * AS_READ_MAX_LEN];
   int  rev_ct [2 * AS_READ_MAX_LEN];
   Screen_Info_t  screen;
   const int  MIN_HITS = 100;
   int  status;
   int  kmer_ct, i, j;

   screen . match = (IntScreenMatch *) Safe_malloc
                         (INIT_SCREEN_MATCHES * sizeof (IntScreenMatch));
   screen . range = (Screen_Range_t *) Safe_malloc
                         (INIT_SCREEN_MATCHES * sizeof (Screen_Range_t));
   screen . match_len = INIT_SCREEN_MATCHES;
   screen . num_matches = 0;

   printf ("\nStart Profile\n");

   profile_stream = File_Open ("dan5.inp", "r");     // inp file
   Read_Msg_Fn = InputFileType_AS (profile_stream);

   while  (EOF != Read_Msg_Fn (profile_stream, & pmesg))
     {
      imesgtype = pmesg -> t;
      switch  (imesgtype)
        {
         case MESG_IFG:
           {
            int  j, len;
            char  * p, * window;
            uint64  key;

            ifg_mesg = pmesg -> m;
            stripWhiteSpace (frag, ifg_mesg -> sequence, AS_READ_MAX_LEN * 2);
            len = strlen (frag);
            
            for  (j = 0;  j < len;  j ++)
              frag [j] = tolower (frag [j]);

            p = window = frag;
            key = 0;

            for  (j = 0;  j < WINDOW_SIZE - 1;  j ++)
              {
               ct [j] = -1;
               key |= (uint64) (Bit_Equivalent [* (p ++)]) << (2 * j);
              }

            while  ((* p) != '\0')
              {
               key |= (uint64) (Bit_Equivalent [* (p ++)]) << (2 * (WINDOW_SIZE - 1));
               ct [j ++] = Hash_Hits (key, window ++);
               key >>= 2;
              }


            strcpy (rev_frag, frag);
            Rev_Complement (rev_frag, len);

            p = window = rev_frag;
            key = 0;

            for  (j = 0;  j < WINDOW_SIZE - 1;  j ++)
              {
               rev_ct [j] = -1;
               key |= (uint64) (Bit_Equivalent [* (p ++)]) << (2 * j);
              }

            while  ((* p) != '\0')
              {
               key |= (uint64) (Bit_Equivalent [* (p ++)]) << (2 * (WINDOW_SIZE - 1));
               rev_ct [j ++] = Hash_Hits (key, window ++);
               key >>= 2;
              }


            printf ("\n> Frag %d:  %s\n", ifg_mesg -> iaccession, ifg_mesg -> source);
            printf ("%8s   %8s\n", "Forward", "Reverse");
            for  (j = 0;  j < len;  j ++)
              {
               if  (j >= WINDOW_SIZE - 1)
                   printf ("%c %6d", frag [j], ct [j]);
                 else
                   printf ("%c %6s", frag [j], "-");
               if  (j >= WINDOW_SIZE - 1)
                   printf ("     %c %6d\n",
                           rev_frag [len - 1 - j], rev_ct [len - 1 - j + WINDOW_SIZE - 1]);
                 else
                   printf ("     %c %6s\n", rev_frag [len - 1 - j], "-");
              }
           }
        }
     }

   fclose (profile_stream);

   printf ("\nKmers with at least %d occurrences\n", MIN_HITS);
   kmer_ct = 0;
   for  (i = 0;  i < HASH_TABLE_SIZE;  i ++)
     for  (j = 0;  j < Hash_Table [i] . Entry_Ct;  j ++)
       {
        String_Ref_t  h_ref;
        char  * p;

        kmer_ct ++;
        if  (Hash_Table [i] . Hits [j] >= MIN_HITS)
            {
             h_ref = Hash_Table [i] . Entry [j];
             if  (! h_ref . Last)
                 h_ref = Extra_Ref_Space
                             [(h_ref . String_Num << 10) + h_ref . Offset];

             p = Data + String_Start [h_ref . String_Num] + h_ref . Offset;
             printf ("%.20s %7d\n", p, Hash_Table [i] . Hits [j]);
            }
       }
   printf ("\nTotal kmers = %d\n", kmer_ct);
#endif

#if  SHOW_HI_HIT_KMERS
// #if  0   // Special Mouse
{
 Distrib_t  Hash_Hits_Dist;
 FILE  * fp = NULL;
 const int  MIN_HITS = 100;
 int  i, j;

 Init_Distrib (& Hash_Hits_Dist, 381);
 {
   int  i;
   for  (i = 0;  i <= 100;  i ++)
     Hash_Hits_Dist . Thold [i] = i;
   for  (i = 101;  i <= 190;  i ++)
     Hash_Hits_Dist . Thold [i] = 10.0 * (i - 90);
   for  (i = 191;  i <= 380;  i ++)
     Hash_Hits_Dist . Thold [i] = 100.0 * (i - 180);
 }

 fp = File_Open ("hihits.kmer", "w");
 for  (i = 0;  i < HASH_TABLE_SIZE;  i ++)
   for  (j = 0;  j < Hash_Table [i] . Entry_Ct;  j ++)
     {
      String_Ref_t  h_ref;
      char  * p;

      Incr_Distrib (& Hash_Hits_Dist, Hash_Table [i] . Hits [j]);

      if  (Hash_Table [i] . Hits [j] >= MIN_HITS)
          {
           String_Ref_t  h_ref;
           char  * p;

           h_ref = Hash_Table [i] . Entry [j];
           if  (! h_ref . Last)
               h_ref = Extra_Ref_Space
                           [(h_ref . String_Num << 10) + h_ref . Offset];

           p = Data + String_Start [h_ref . String_Num] + h_ref . Offset;
           fprintf (fp, "%.20s %7d\n", p, Hash_Table [i] . Hits [j]);
          }
     }
 fclose (fp);

 Print_Distrib (Hash_Hits_Dist, "Hash Table Kmer Hits");

 fclose (Stat_File);
}
#endif

   return;
  }



static void  Put_String_In_Hash
    (int i)

//  Insert string subscript  i  into the global hash table.
//  Sequence and information about the string are in
//  global variables  Data, String_Start, String_Info, ....

  {
   String_Ref_t  ref;
   char  * p, * window;
   int  kmers_inserted = 0;
   int  skip_ct;
   int  screen_sub, screen_lo, screen_hi;
#if  SCREEN_CHECK_ONLY
   int  kmers_screened = 0;
#endif
   uint64  key, key_is_bad;
   int  j;

   if  (String_Info [i] . length < WINDOW_SIZE)
       return;

   screen_sub = Screen_Sub [i];
   if  (screen_sub == 0)
       screen_lo = screen_hi = INT_MAX;
     else
       {
        screen_lo = Screen_Space [screen_sub] . bgn;
        screen_hi = Screen_Space [screen_sub] . end;
       }

   p = window = Data + String_Start [i];
   key = key_is_bad = 0;
   for  (j = 0;  j < WINDOW_SIZE;  j ++)
     {
      key_is_bad |= (uint64) (Char_Is_Bad [(int) * p]) << j;
      key |= (uint64) (Bit_Equivalent [(int) * (p ++)]) << (2 * j);
     }

   ref . String_Num = i;
   ref . Offset = 0;
   skip_ct = 0;
   ref . Empty = FALSE;

   if  ((int) (ref . Offset) <= screen_lo - WINDOW_SIZE + WINDOW_SCREEN_OLAP
            && ! key_is_bad)
       {
#if  ! SCREEN_CHECK_ONLY
        Hash_Insert (ref, key, window);
#endif
        kmers_inserted ++;
       }
#if  SCREEN_CHECK_ONLY
     else
       kmers_screened ++;
#endif

   while  ((* p) != '\0')
     {
      if  ((int) (ref . Offset) == screen_hi - 1 - WINDOW_SCREEN_OLAP)
          {
           if  (Screen_Space [screen_sub] . last)
               screen_lo = screen_hi = INT_MAX;
             else
               {
                screen_sub ++;
                screen_lo = Screen_Space [screen_sub] . bgn;
                screen_hi = Screen_Space [screen_sub] . end;
               }
          }

      window ++;
      ref . Offset ++;
      if  (++ skip_ct > HASH_KMER_SKIP)
          skip_ct = 0;

      key_is_bad >>= 1;
      key_is_bad |= (uint64) (Char_Is_Bad [(int) * p]) << (WINDOW_SIZE - 1);
      key >>= 2;
      key |= (uint64)
                 (Bit_Equivalent [(int) * (p ++)]) << (2 * (WINDOW_SIZE - 1));

      if  (skip_ct == 0
               && (int) (ref . Offset)
                     <= screen_lo - WINDOW_SIZE + WINDOW_SCREEN_OLAP
               && ! key_is_bad)
          {
#if  ! SCREEN_CHECK_ONLY
           Hash_Insert (ref, key, window);
#endif
           kmers_inserted ++;
          }
#if  SCREEN_CHECK_ONLY
        else
          kmers_screened ++;
#endif
     }
#if  LIST_TOTALLY_SCREENED
   if  (kmers_inserted == 0)
       fprintf (Total_Screen_File, "%ld\n", first_frag_id + i);
#endif

   return;
  }



static int  Read_Next_Frag
    (char frag [AS_READ_MAX_LEN + 1], 
     char quality [AS_READ_MAX_LEN + 1], 
     FragStreamHandle stream,
     ReadStructp myRead,
     Screen_Info_t * screen,
     uint32 * last_frag_read)

/* Read the next fragment from fragment stream  stream  and store it in  frag ,
*  with its quality values in  quality .  Put the read itself in  myRead
*  and the screen information in  screen .  Set  (* last_frag_read)
*  to the iid of the fragment read.
*  If successful, return  TRUE ; if there is no fragment, return  FALSE . */

  {
   int  i, return_val, success, match_ct;
   uint32  deleted;
   size_t  frag_len;
#if  USE_SOURCE_FIELD
   char  frag_source [MAX_SOURCE_LENGTH + 1];
#endif
   
#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_lock (& FragStore_Mutex);
#endif

#if  USE_SOURCE_FIELD
   success = nextFragStream (stream, myRead, FRAG_S_SEQUENCE | FRAG_S_SOURCE);
#else
   success = nextFragStream (stream, myRead, FRAG_S_SEQUENCE);
#endif
   getIsDeleted_ReadStruct (myRead, & deleted);
   if  (! success)
       return_val = 0;
   else if  (deleted)
       {
        getReadIndex_ReadStruct (myRead, last_frag_read);
        return_val = DELETED_FRAG;
       }
     else
       {
        int  result;
        uint  clear_start, clear_end;

        getReadIndex_ReadStruct (myRead, last_frag_read);
        result = getSequence_ReadStruct
                     (myRead, frag, quality, AS_READ_MAX_LEN + 1);
        if  (result != 0)
            {
             fprintf (stderr,
                      "Error reading frag store tried to fit %d in a buffer of %d\n",
                      result, AS_READ_MAX_LEN + 1);
             exit (EXIT_FAILURE);
            }

        // Make sure that we have a lowercase sequence string
        {
         int i;
         int len = strlen (frag);

         for  (i = 0;  i < len;  i ++)
           frag [i] = tolower (frag [i]);
        }

        if  (Doing_Partial_Overlaps)
            result = getClearRegion_ReadStruct
                         (myRead, & clear_start, & clear_end, READSTRUCT_ORIGINAL);
          else
            result = getClearRegion_ReadStruct
                         (myRead, & clear_start, & clear_end, READSTRUCT_OVL);
        if  (result != 0)
            {
             fprintf (stderr, "Error reading frag store\n");
             exit (EXIT_FAILURE);
            }

        frag_len = clear_end - clear_start;
        if  (clear_start > 0)
            {
             memmove (frag, frag + clear_start, frag_len);
             memmove (quality, quality + clear_start, frag_len);
            }
        frag [frag_len] = '\0';
        quality [frag_len] = '\0';
        for  (i = 0;  i < frag_len;  i ++)
          quality [i] -= QUALITY_BASE_CHAR;

        match_ct = getNumScreenMatches_ReadStruct (myRead);
        if  (match_ct > screen -> match_len)
            {
             screen -> match_len = match_ct;
             fprintf (stderr,
                      "### realloc  screen match & range  match_ct = %d\n",
                      match_ct);
             screen -> match = (IntScreenMatch *) Safe_realloc
                                   (screen -> match,
                                    match_ct * sizeof (IntScreenMatch));
             screen -> range = (Screen_Range_t *) Safe_realloc
                                   (screen -> range,
                                    match_ct * sizeof (Screen_Range_t));
            }
        getScreenMatches_ReadStruct (myRead, screen -> match,
                                     match_ct);

        if  (Ignore_Screen_Info)
            screen -> num_matches = 0;
          else
            screen -> num_matches = match_ct;

        Get_Relevant_Screen_Matches (screen);
        screen -> left_end_screened
            = (screen -> num_matches > 0
                && screen -> range [0] . bgn < HOPELESS_MATCH);
        screen -> right_end_screened
            = (screen -> num_matches > 0
                && screen -> range [screen -> num_matches - 1] . end
                     > frag_len - HOPELESS_MATCH);

#if  USE_SOURCE_FIELD
        getSource_ReadStruct (myRead, frag_source, MAX_SOURCE_LENGTH);
        Find_End_Annotations (frag_source, frag_len,
                              & Global_Left_Annotation,
                              & Global_Right_Annotation);
#endif

        return_val = VALID_FRAG;
       }

#if  USE_THREADS
   if  (Num_PThreads > 1)
       pthread_mutex_unlock (& FragStore_Mutex);
#endif

   return  return_val;
  }



static void  Read_uint32_List
    (char * file_name, uint32 * * list, int * n)

//  Open  file_name  and use the list of integers in it to build
//  a sorted (ascending) list of integers in  (* list) .
//  Set  (* n)  to the number of entries in that list.

  {
   FILE  * fp = File_Open (file_name, "r");
   int  num;


   (* n) = 0;
   while  (fscanf (fp, "%u", & num) == 1)
     {
      (* n) ++;
      (* list) = (uint32 *) Safe_realloc ((* list), (* n) * sizeof (uint32));
      (* list) [(* n) - 1] = num;
     }

   if  ((* n) > 0)
       qsort ((* list), (* n), sizeof (uint32), By_uint32);

   fclose (fp);

   return;
  }




static void  Rev_Complement
    (char * S, int Len)

/* Set  S [0 .. Len - 1]  to its DNA reverse complement. */

  {
   char  Ch;
   int  i, j;

   for  (i = 0, j = Len - 1;  i < j;  i ++, j --)
     {
      Ch = Complement (S [i]);
      S [i] = Complement (S [j]);
      S [j] = Ch;
     }

   if  (i == j)
       S [i] = Complement (S [i]);

   return;
  }



static int  Rev_Prefix_Edit_Dist
    (char A [], int m, char T [], int n, int Error_Limit,
     int * A_End, int * T_End, int * Leftover, int * Match_To_End,
     Work_Area_t * WA)

//  Return the minimum number of changes (inserts, deletes, replacements)
//  needed to match string  A [0 .. (1-m)]  right-to-left with a prefix of string
//   T [0 .. (1-n)]  right-to-left if it's not more than  Error_Limit .
//  If no match, return the number of errors for the best match
//  up to a branch point.
//  Put delta description of alignment in  WA -> Left_Delta  and set
//  WA -> Left_Delta_Len  to the number of entries there.
//  Set  A_End  and  T_End  to the leftmost positions where the
//  alignment ended in  A  and  T , respectively.
//  If the alignment succeeds set  Leftover  to the number of
//  characters that match after the last  WA -> Left_Delta  entry;
//  otherwise, set  Leftover  to zero.
//  Set  Match_To_End  true if the match extended to the end
//  of at least one string; otherwise, set it false to indicate
//  a branch point.

  {
#if 0
   double  Cost_Sum, Min_Cost_Sum;
   int Adjustment, Cost_Sub;
#endif
   double  Score, Max_Score;
   int  Max_Score_Len = 0, Max_Score_Best_d = 0, Max_Score_Best_e = 0;
   int  Tail_Len;
   int  Best_d, Best_e, From, Last, Longest, Max, Row;
   int  Left, Right;
   int  d, e, j, k;


#if  SHOW_STATS
Edit_Dist_Ct ++;
#endif
   assert (m <= n);
   Best_d = Best_e = Longest = 0;
   WA -> Left_Delta_Len = 0;

   for  (Row = 0;  Row < m
                     && (A [- Row] == T [- Row]
                          || A [- Row] == DONT_KNOW_CHAR
                          || T [- Row] == DONT_KNOW_CHAR);  Row ++)
     ;

   WA -> Edit_Array [0] [0] = Row;

   if  (Row == m)
       {
        (* A_End) = (* T_End) = - m;
        (* Leftover) = m;
#if  SHOW_STATS
Incr_Distrib (& Edit_Depth_Dist, 0);
#endif
        (* Match_To_End) = TRUE;
        return  0;
       }

   Left = Right = 0;
   Max_Score = 0.0;
   for  (e = 1;  e <= Error_Limit;  e ++)
     {
      Left = OVL_Max_int (Left - 1, -e);
      Right = OVL_Min_int (Right + 1, e);
      WA -> Edit_Array [e - 1] [Left] = -2;
      WA -> Edit_Array [e - 1] [Left - 1] = -2;
      WA -> Edit_Array [e - 1] [Right] = -2;
      WA -> Edit_Array [e - 1] [Right + 1] = -2;

      for  (d = Left;  d <= Right;  d ++)
        {
         Row = 1 + WA -> Edit_Array [e - 1] [d];
         if  ((j = WA -> Edit_Array [e - 1] [d - 1]) > Row)
             Row = j;
         if  ((j = 1 + WA -> Edit_Array [e - 1] [d + 1]) > Row)
             Row = j;
         while  (Row < m && Row + d < n
                  && (A [- Row] == T [- Row - d]
                        || A [- Row] == DONT_KNOW_CHAR
                        || T [- Row - d] == DONT_KNOW_CHAR))
           Row ++;

         WA -> Edit_Array [e] [d] = Row;

         if  (Row == m || Row + d == n)
             {
#if  SHOW_STATS
Incr_Distrib (& Edit_Depth_Dist, e);
#endif

              //  Check for branch point here caused by uneven
              //  distribution of errors

              Score = Row * Branch_Match_Value - e;
                        // Assumes  Branch_Match_Value
                        //             - Branch_Error_Value == 1.0
              Tail_Len = Row - Max_Score_Len;
              if  ((Doing_Partial_Overlaps && Score < Max_Score)
                     || (e > MIN_BRANCH_END_DIST / 2
                           && Tail_Len >= MIN_BRANCH_END_DIST
                           && (Max_Score - Score) / Tail_Len >= MIN_BRANCH_TAIL_SLOPE))
                  {
                   (* A_End) = - Max_Score_Len;
                   (* T_End) = - Max_Score_Len - Max_Score_Best_d;
                   Set_Left_Delta (Max_Score_Best_e, Max_Score_Best_d,
                        Leftover, T_End, n, WA);
                   (* Match_To_End) = FALSE;
                   return  Max_Score_Best_e;
                  }

              (* A_End) = - Row;           // One past last align position
              (* T_End) = - Row - d;
              Set_Left_Delta (e, d, Leftover, T_End, n, WA);
              (* Match_To_End) = TRUE;
              return  e;
             }
        }

      while  (Left <= Right && Left < 0
                  && WA -> Edit_Array [e] [Left] < WA -> Edit_Match_Limit [e])
        Left ++;
      if  (Left >= 0)
          while  (Left <= Right
                    && WA -> Edit_Array [e] [Left] + Left < WA -> Edit_Match_Limit [e])
            Left ++;
      if  (Left > Right)
          break;
      while  (Right > 0
                  && WA -> Edit_Array [e] [Right] + Right < WA -> Edit_Match_Limit [e])
        Right --;
      if  (Right <= 0)
          while  (WA -> Edit_Array [e] [Right] < WA -> Edit_Match_Limit [e])
            Right --;
      assert (Left <= Right);

      for  (d = Left;  d <= Right;  d ++)
        if  (WA -> Edit_Array [e] [d] > Longest)
            {
             Best_d = d;
             Best_e = e;
             Longest = WA -> Edit_Array [e] [d];
            }

      Score = Longest * Branch_Match_Value - e;
               // Assumes  Branch_Match_Value - Branch_Error_Value == 1.0
      if  (Score > Max_Score)
          {
           Max_Score = Score;
           Max_Score_Len = Longest;
           Max_Score_Best_d = Best_d;
           Max_Score_Best_e = Best_e;
          }
     }

#if  SHOW_STATS
Incr_Distrib (& Edit_Depth_Dist, 1 + Error_Limit);
#endif

   (* A_End) = - Max_Score_Len;
   (* T_End) = - Max_Score_Len - Max_Score_Best_d;
   Set_Left_Delta (Max_Score_Best_e, Max_Score_Best_d,
        Leftover, T_End, n, WA);
   (* Match_To_End) = FALSE;
   return  Max_Score_Best_e;
  }



static void  Reverse_String
    (char * s, int len)

//  Reverse the order of characters in  s [0 .. len - 1] .

  {
   int  i, j;

   for  (i = 0, j = len - 1;  i < j;  i ++, j --)
     {
      char  ch;

      ch = s [i];
      s [i] = s [j];
      s [j] = ch;
     }

   return;
  }



static void  Set_Left_Delta
    (int e, int d, int * leftover, int * t_end, int t_len,
     Work_Area_t * WA)

//  Put the delta encoding of the alignment represented in  WA -> Edit_Array
//  starting at row  e  (which is the number of errors) and column  d
//  (which is the diagonal) and working back to the start, into
//  WA -> Left_Delta .  Set  WA -> Left_Delta_Len  to the number of
//  delta entries and set  (* leftover)  to the number of
//  characters that match after the last  WA -> Left_Delta  entry.
//  Don't allow the first delta to be an indel if it can be
//  converted to a substitution by adjusting  (* t_end)  which
//  is where the alignment ends in the T string, which has length
//   t_len .

  {
   int  from, last, max;
   int  j, k;

   last = WA -> Edit_Array [e] [d];
   WA -> Left_Delta_Len = 0;
   for  (k = e;  k > 0;  k --)
     {
      from = d;
      max = 1 + WA -> Edit_Array [k - 1] [d];
      if  ((j = WA -> Edit_Array [k - 1] [d - 1]) > max)
          {
           from = d - 1;
           max = j;
          }
      if  ((j = 1 + WA -> Edit_Array [k - 1] [d + 1]) > max)
          {
           from = d + 1;
           max = j;
          }
      if  (from == d - 1)
          {
           WA -> Left_Delta [WA -> Left_Delta_Len ++] = max - last - 1;
           d --;
           last = WA -> Edit_Array [k - 1] [from];
          }
      else if  (from == d + 1)
          {
           WA -> Left_Delta [WA -> Left_Delta_Len ++] = last - (max - 1);
           d ++;
           last = WA -> Edit_Array [k - 1] [from];
          }
     }
   (* leftover) = last;

   // Don't allow first delta to be +1 or -1
   assert (WA -> Left_Delta_Len == 0 || WA -> Left_Delta [0] != -1);
   if  (WA -> Left_Delta_Len > 0 && WA -> Left_Delta [0] == 1
          && (* t_end) + t_len > 0)
       {
        int  i;

        if  (WA -> Left_Delta [1] > 0)
            WA -> Left_Delta [0] = WA -> Left_Delta [1] + 1;
          else
            WA -> Left_Delta [0] = WA -> Left_Delta [1] - 1;
        for  (i = 2;  i < WA -> Left_Delta_Len;  i ++)
          WA -> Left_Delta [i - 1] = WA -> Left_Delta [i];
        WA -> Left_Delta_Len --;
        (* t_end) --;
        if  (WA -> Left_Delta_Len == 0)
            (* leftover) ++;
       }

   return;
  }



static void  Set_Right_Delta
    (int e, int d, Work_Area_t * WA)

//  Put the delta encoding of the alignment represented in  WA -> Edit_Array
//  starting at row  e  (which is the number of errors) and column  d
//  (which is the diagonal) and working back to the start, into
//  WA -> Right_Delta .  Set  WA -> Right_Delta_Len  to the number of
//  delta entries.

  {
   int  delta_stack [MAX_ERRORS];
   int  from, last, max;
   int  i, j, k;

   last = WA -> Edit_Array [e] [d];
   WA -> Right_Delta_Len = 0;
   for  (k = e;  k > 0;  k --)
     {
      from = d;
      max = 1 + WA -> Edit_Array [k - 1] [d];
      if  ((j = WA -> Edit_Array [k - 1] [d - 1]) > max)
          {
           from = d - 1;
           max = j;
          }
      if  ((j = 1 + WA -> Edit_Array [k - 1] [d + 1]) > max)
          {
           from = d + 1;
           max = j;
          }
      if  (from == d - 1)
          {
           delta_stack [WA -> Right_Delta_Len ++] = max - last - 1;
           d --;
           last = WA -> Edit_Array [k - 1] [from];
          }
      else if  (from == d + 1)
          {
           delta_stack [WA -> Right_Delta_Len ++] = last - (max - 1);
           d ++;
           last = WA -> Edit_Array [k - 1] [from];
          }
     }
   delta_stack [WA -> Right_Delta_Len ++] = last + 1;

   k = 0;
   for  (i = WA -> Right_Delta_Len - 1;  i > 0;  i --)
     WA -> Right_Delta [k ++]
         = abs (delta_stack [i]) * Sign (delta_stack [i - 1]);
   WA -> Right_Delta_Len --;

   return;
  }



static void  Set_Screened_Ends
    (void)

//  Set  left/right_end_screened  fields true in global  String_Info
//  for reads with a kmer count higher than  Hi_Hit_Limit
//  sufficiently near the respective end of the the string.
//  Also mark the hash table entry as empty so it won't be
//  coalesced.

  {
   int  i, j;

   for  (i = 0;  i < HASH_TABLE_SIZE;  i ++)
     for  (j = 0;  j < Hash_Table [i] . Entry_Ct;  j ++)
       if  (Hash_Table [i] . Hits [j] >= Hi_Hit_Limit)
           {
            Mark_Screened_Ends_Chain (Hash_Table [i] . Entry [j]);
            Hash_Table [i] . Entry [j] . Empty = TRUE;
           }

   return;
  }



void  Show_Alignment
    (char * S, char * T, Olap_Info_t * p)

//  Display the alignment between strings  S  and  T  indicated
//  in  (* p) .

  {
   char  S_Line [AS_READ_MAX_LEN + 1], T_Line [AS_READ_MAX_LEN + 1];
   char  X_Line [AS_READ_MAX_LEN + 1];
   int  i, j, ks, kt, ns, nt, ct;

   i = p -> s_lo;
   j = p -> t_lo;
   ks = kt = 0;
   ns = nt = 1;
   while  (i <= p -> s_hi)
     {
      ct = 0;
      while  (i <= p -> s_hi && ct < DISPLAY_WIDTH)
        {
         if  (ks < p -> delta_ct && ns == abs (p -> delta [ks]))
             {
              if  (p -> delta [ks] < 0)
                  S_Line [ct] = '-';
                else
                  S_Line [ct] = S [i ++];
              ks ++;
              ns = 1;
             }
           else
             {
              S_Line [ct] = S [i ++];
              ns ++;
             }
         ct ++;
        }
      S_Line [ct] = '\0';
      printf ("%s\n", S_Line);

      ct = 0;
      while  (j <= p -> t_hi && ct < DISPLAY_WIDTH)
        {
         if  (kt < p -> delta_ct && nt == abs (p -> delta [kt]))
             {
              if  (p -> delta [kt] > 0)
                  T_Line [ct] = '-';
                else
                  T_Line [ct] = T [j ++];
              kt ++;
              nt = 1;
             }
           else
             {
              T_Line [ct] = T [j ++];
              nt ++;
             }
         if  (S_Line [ct] == T_Line [ct])
             X_Line [ct] =  ' ';
         else if  (S_Line [ct] == DONT_KNOW_CHAR
                     || T_Line [ct] == DONT_KNOW_CHAR)
             X_Line [ct] =  '?';
           else
             X_Line [ct] =  '^';
         ct ++;
        }
      T_Line [ct] = X_Line [ct] = '\0';
      printf ("%s\n", T_Line);
      printf ("%s\n", X_Line);

      putchar ('\n');
     }
  }



static void  Show_Match
    (Match_Node_t * P, char * S, char * T)

//  Show the match information in  P 's match node for strings
//  S  and  T .

  {
   printf (
   "  Start = %3d  Offset = %3d  Len = %3d\n",
                P -> Start, P -> Offset, P -> Len);
   printf ("    %-s\n", S + P -> Start);
   printf ("    %-s\n", T + P -> Offset);
  }



static void  Show_Overlap
    (char * a, int a_len, char * a_quality,
     char * b, int b_len, char * b_quality, Olap_Info_t * p)

//  Show the overlap recorded in  (* p)  between strings  A  and  B ,
//  with lengths  A_Len  and  B_Len, respectively.

  {
   int  i, j, k, m, diag, top_len, bottom_len;
   char  top [2000], bottom [2000];
   char  q_top [2000], q_bottom [2000];

#if  DO_OLAP_ALIGN_PROFILE
i = p -> s_lo;
j = p -> t_lo;
for  (k = 0;  k < p -> delta_ct;  k ++)
  {
   for  (m = 1;  m < abs (p -> delta [k]);  m ++)
     {
      Align_P [j] [Bit_Equivalent [a [i]]] ++;
      i ++;
      j ++;
     }
   if  (p -> delta [k] > 0)
       {
        Align_P [j] [4] ++;          // insert
        i ++;
       }
     else
       {
        Align_P [j] [5] ++;          // delete
        j ++;
       }
  }
while  (i < a_len && j < b_len)
  {
   Align_P [j] [Bit_Equivalent [a [i]]] ++;
   i ++;
   j ++;
  }

return;
#endif


   top_len = bottom_len = 0;
   diag = p -> t_lo - p -> s_lo;
   if  (diag >= 0)
       for  (i = j = 0;  j < diag;  j ++)
         {
          q_top [top_len] = ' ';
          q_bottom [bottom_len] = b_quality [j] + QUALITY_BASE_CHAR;
          top [top_len ++] = ' ';
          bottom [bottom_len ++] = b [j];
         }
     else
       for  (i = j = 0;  i < - diag;  i ++)
         {
          q_top [top_len] = a_quality [i] + QUALITY_BASE_CHAR;
          q_bottom [bottom_len] = ' ';
          top [top_len ++] = a [i];
          bottom [bottom_len ++] = ' ';
         }

   while  (i < p -> s_lo && j < p -> t_lo)
     {
      q_top [top_len] = a_quality [i] + QUALITY_BASE_CHAR;
      q_bottom [bottom_len] = b_quality [j] + QUALITY_BASE_CHAR;
      top [top_len ++] = a [i ++];
      bottom [bottom_len ++] = b [j ++];
     }
   if  (i != p -> s_lo || j != p -> t_lo)
       {
        printf ("OOOPS\n");
        exit (-4);
       }


   for  (k = 0;  k < p -> delta_ct;  k ++)
     {
      for  (m = 1;  m < abs (p -> delta [k]);  m ++)
        {
         q_top [top_len] = a_quality [i] + QUALITY_BASE_CHAR;
         top [top_len ++] = a [i ++];
        }
      if  (p -> delta [k] < 0)
          {
           q_top [top_len] = ' ';
           top [top_len ++] = '-';
          }
        else
          {
           q_top [top_len] = a_quality [i] + QUALITY_BASE_CHAR;
           top [top_len ++] = a [i ++];
          }
     }
   while  (i < a_len)
     {
      q_top [top_len] = a_quality [i] + QUALITY_BASE_CHAR;
      top [top_len ++] = a [i ++];
     }
   q_top [top_len] = '\0';
   top [top_len] = '\0';
     

   for  (k = 0;  k < p -> delta_ct;  k ++)
     {
      for  (m = 1;  m < abs (p -> delta [k]);  m ++)
        {
         q_bottom [bottom_len] = b_quality [j] + QUALITY_BASE_CHAR;
         bottom [bottom_len ++] = b [j ++];
        }
      if  (p -> delta [k] > 0)
          {
           q_bottom [bottom_len] = ' ';
           bottom [bottom_len ++] = '-';
          }
        else
          {
           q_bottom [bottom_len] = b_quality [j] + QUALITY_BASE_CHAR;
           bottom [bottom_len ++] = b [j ++];
          }
     }
   while  (j < b_len)
     {
      q_bottom [bottom_len] = b_quality [j] + QUALITY_BASE_CHAR;
      bottom [bottom_len ++] = b [j ++];
     }
   q_bottom [bottom_len] = '\0';
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
        if  (top [i + j] == ' ' || bottom [i + j] == ' '
               ||  top [i + j] == bottom [i + j])
            putchar (' ');
        else if  (top [i + j] == DONT_KNOW_CHAR
                    || bottom [i + j] == DONT_KNOW_CHAR)
            putchar ('?');
          else
            putchar ('^');
      putchar ('\n');
#if  0
      printf ("   ");
      for  (j = 0;  j < DISPLAY_WIDTH && i + j < top_len;  j ++)
        putchar (q_top [i + j]);
      putchar ('\n');
      printf ("   ");
      for  (j = 0;  j < DISPLAY_WIDTH && i + j < bottom_len;  j ++)
        putchar (q_bottom [i + j]);
      putchar ('\n');
#endif
     }

   return;
  }



#if  SHOW_SNPS
#define  SHOW_SNP_DETAILS     0
#define  SIDE_BASES           4
#define  MIN_SIDE_QUALITY    20
#define  INIT_MASK            (1 << (2 * SIDE_BASES)) - 1
#define  QUAL_MATCH_MASK      ((1 << (2 * SIDE_BASES + 1)) - 1)
#define  SNP_MATCH_MASK       (QUAL_MATCH_MASK ^ (1 << SIDE_BASES))
#define  END_AVOID_LEN       30

static void  Show_SNPs
    (char * a, int a_len, char * a_quality,
     char * b, int b_len, char * b_quality, Olap_Info_t * p,
     Work_Area_t * WA, int * mismatch_ct, int * indel_ct,
     int * olap_bases, int * unscreened_olap_bases,
     int * all_match_ct, int * hi_qual_ct,
     int local_msnp_bin [], int local_isnp_bin [],
     int local_match_bin [], int local_other_bin [])

//  Show the overlap recorded in  (* p)  between strings  A  and  B ,
//  with lengths  A_Len  and  B_Len, respectively.  WA has screen
//  information about  A .  Set  (* mismatch_ct)  to the number
//  of mismatch SNPs and  (* indel_ct)  to the number of insert/delete
//  SNPs.  Set  (* olap_bases)  to the number of bases in  A  involved
//  in the overlap.  Set  (* unscreened_olap_bases)  to the number
//  of bases in the overlap that were not in screened regions.

  {
   int  i, j, k, m, n;
   int  indel_ring [SIDE_BASES + 1] = {0};
   int  quality_ring [SIDE_BASES + 1] = {0};
   int  screen_sub, screen_lo, screen_hi;
   unsigned int  match_mask = 0, quality_mask = 0, bit;

   (* mismatch_ct) = 0;
   (* indel_ct) = 0;
   (* olap_bases) = 0;
   (* unscreened_olap_bases) = 0;
   (* all_match_ct) = 0;
   (* hi_qual_ct) = 0;

   for  (k = 0;  k <= SIDE_BASES;  k ++)
     indel_ring [k] = 3;

   screen_sub = 0;
   if  (WA -> screen_info . num_matches == 0)
       screen_lo = screen_hi = INT_MAX;
     else
       {
        screen_lo = WA -> screen_info . range [0] . bgn;
        screen_hi = WA -> screen_info . range [0] . end;
       }
#if  SHOW_SNP_DETAILS
   fprintf (Out_Stream, "Screen range: %10d %10d\n",
            screen_lo, screen_hi);
#endif

   i = p -> s_lo;
   j = p -> t_lo;

   while  (i > screen_hi)
     {
      if  (WA -> screen_info . range [screen_sub] . last)
          screen_lo = screen_hi = INT_MAX;
        else
          {
           screen_sub ++;
           screen_lo = WA -> screen_info . range [screen_sub] . bgn;
           screen_hi = WA -> screen_info . range [screen_sub] . end;
          }
#if  SHOW_SNP_DETAILS
      fprintf (Out_Stream, "Screen range: %10d %10d\n",
               screen_lo, screen_hi);
#endif
     }

   n = 0;
   for  (k = 0;  k < p -> delta_ct;  k ++)
     {
      for  (m = 1;  m < abs (p -> delta [k]);  m ++)
        {
         match_mask &= INIT_MASK;
         bit = ((a [i] == b [j]) ? 1 : 0);
         match_mask = (match_mask << 1) | bit;

         quality_mask &= INIT_MASK;
         bit = ((a_quality [i] >= MIN_SIDE_QUALITY
                   && b_quality [j] >= MIN_SIDE_QUALITY) ? 1 : 0);
         quality_mask = (quality_mask << 1) | bit;

         if  (i < END_AVOID_LEN
                || i >= a_len - END_AVOID_LEN
                || j < END_AVOID_LEN
                || j >= b_len - END_AVOID_LEN)
             indel_ring [n % (SIDE_BASES + 1)] = 3;  // too near end
         else if  (i < screen_lo)
             indel_ring [n % (SIDE_BASES + 1)] = 0;  // mismatch
           else
             indel_ring [n % (SIDE_BASES + 1)] = 2;  // screened
         quality_ring [n % (SIDE_BASES + 1)]
             = OVL_Min_int (a_quality [i], b_quality [j]);

         if  (quality_mask == QUAL_MATCH_MASK
                && indel_ring [(n + 1) % (SIDE_BASES + 1)] < 2)
             {
              if  (match_mask == SNP_MATCH_MASK)
                  {
                   switch  (indel_ring [(n + 1) % (SIDE_BASES + 1)])
                     {
                      case  0 :    // mismatch
#if  SHOW_SNP_DETAILS
                        fprintf (Out_Stream, "SNP  %4d  %4d  %.*s/%.*s\n",
                                 i - SIDE_BASES, j - SIDE_BASES,
                                 2 * SIDE_BASES + 1, a + i - 2 * SIDE_BASES,
                                 2 * SIDE_BASES + 1, b + j - 2 * SIDE_BASES);
                        fprintf (Out_Stream,
                                 "     i = %4d  j = %4d  qual_mask = %04o"
                                 "  match_mask = %04o\n",
                                 i, j, quality_mask, match_mask);
#endif
                        (* mismatch_ct) ++;
                        Global_Match_SNP_Bin
                            [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                        local_msnp_bin
                            [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                        break;
                      case  1 :    // insert
#if  SHOW_SNP_DETAILS
                        fprintf (Out_Stream, "SNP  %4d  %4d  %.*s-%.*s/%.*s\n",
                                 i - SIDE_BASES, j - SIDE_BASES,
                                 SIDE_BASES, a + i - 2 * SIDE_BASES + 1,
                                 SIDE_BASES, a + i - SIDE_BASES + 1,
                                 2 * SIDE_BASES + 1, b + j - 2 * SIDE_BASES);
                        fprintf (Out_Stream,
                                 "     i = %4d  j = %4d  qual_mask = %04o"
                                 "  match_mask = %04o\n",
                                 i, j, quality_mask, match_mask);
#endif
                        (* indel_ct) ++;
                        Global_Indel_SNP_Bin
                            [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                        local_isnp_bin
                            [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                        break;
                      case  -1 :    // delete
#if  SHOW_SNP_DETAILS
                        fprintf (Out_Stream, "SNP  %4d  %4d  %.*s/%.*s-%.*s\n",
                                 i - SIDE_BASES, j - SIDE_BASES,
                                 2 * SIDE_BASES + 1, a + i - 2 * SIDE_BASES,
                                 SIDE_BASES, b + j - 2 * SIDE_BASES + 1,
                                 SIDE_BASES, b + j - SIDE_BASES + 1);
                        fprintf (Out_Stream,
                                 "     i = %4d  j = %4d  qual_mask = %04o"
                                 "  match_mask = %04o\n",
                                 i, j, quality_mask, match_mask);
#endif
                        (* indel_ct) ++;
                        Global_Indel_SNP_Bin
                            [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                        local_isnp_bin
                            [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                        break;
                      case  2 :    // screened
#if  SHOW_SNP_DETAILS
                        fprintf (Out_Stream, "SNP  %4d  %4d  screened\n",
                                 i - SIDE_BASES, j - SIDE_BASES);
#endif
                        break;
                      case  3 :   // too near end
#if  SHOW_SNP_DETAILS
                        fprintf (Out_Stream, "SNP  %4d  %4d  too near end\n",
                                 i - SIDE_BASES, j - SIDE_BASES);
#endif
                        break;
                      default :
                        assert (FALSE);
                     }
                  }
              else if  (match_mask == QUAL_MATCH_MASK)    // perfect match
                  {
                   (* all_match_ct) ++;
                   Global_All_Match_Bin
                       [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                   local_match_bin
                       [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                  }
                else
                  {
                   (* hi_qual_ct) ++;
                   Global_Hi_Qual_Bin
                       [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                   local_other_bin
                       [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                  }
             }

         if  (indel_ring [(n + 1) % (SIDE_BASES + 1)] != 3)
             {
              (* olap_bases) ++;
              if  (indel_ring [(n + 1) % (SIDE_BASES + 1)] != 2)
                  (* unscreened_olap_bases) ++;
             }

         i ++;
         j ++;
         n ++;

         while  (i > screen_hi)
           {
            if  (WA -> screen_info . range [screen_sub] . last)
                screen_lo = screen_hi = INT_MAX;
              else
                {
                 screen_sub ++;
                 screen_lo = WA -> screen_info . range [screen_sub] . bgn;
                 screen_hi = WA -> screen_info . range [screen_sub] . end;
                }
#if  SHOW_SNP_DETAILS
            fprintf (Out_Stream, "Screen range: %10d %10d\n",
                     screen_lo, screen_hi);
#endif
           }
        }
      if  (p -> delta [k] > 0)
          {
           match_mask &= INIT_MASK;
           bit = 0;
           match_mask = (match_mask << 1) | bit;

           quality_mask &= INIT_MASK;
           bit = ((a_quality [i] >= MIN_SIDE_QUALITY) ? 1 : 0);
           quality_mask = (quality_mask << 1) | bit;

           if  (i < END_AVOID_LEN
                  || i >= a_len - END_AVOID_LEN
                  || j < END_AVOID_LEN
                  || j >= b_len - END_AVOID_LEN)
               indel_ring [n % (SIDE_BASES + 1)] = 3;  // too near end
           else if  (i < screen_lo)
               indel_ring [n % (SIDE_BASES + 1)] = -1;
             else
               indel_ring [n % (SIDE_BASES + 1)] = 2;  // screened

           quality_ring [n % (SIDE_BASES + 1)]
               = a_quality [i];
           if  (quality_mask == QUAL_MATCH_MASK
                  && indel_ring [(n + 1) % (SIDE_BASES + 1)] < 2)
               {
                (* hi_qual_ct) ++;
                Global_Hi_Qual_Bin
                    [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                local_other_bin
                    [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
               }

           if  (indel_ring [(n + 1) % (SIDE_BASES + 1)] != 3)
               {
                (* olap_bases) ++;
                if  (indel_ring [(n + 1) % (SIDE_BASES + 1)] != 2)
                    (* unscreened_olap_bases) ++;
               }

           i ++;
           n ++;

           while  (i > screen_hi)
             {
              if  (WA -> screen_info . range [screen_sub] . last)
                  screen_lo = screen_hi = INT_MAX;
                else
                  {
                   screen_sub ++;
                   screen_lo = WA -> screen_info . range [screen_sub] . bgn;
                   screen_hi = WA -> screen_info . range [screen_sub] . end;
                  }
#if  SHOW_SNP_DETAILS
              fprintf (Out_Stream, "Screen range: %10d %10d\n",
                       screen_lo, screen_hi);
#endif
             }
          }
        else
          {
           match_mask &= INIT_MASK;
           bit = 0;
           match_mask = (match_mask << 1) | bit;

           quality_mask &= INIT_MASK;
           bit = ((b_quality [j] >= MIN_SIDE_QUALITY) ? 1 : 0);
           quality_mask = (quality_mask << 1) | bit;

           if  (i < END_AVOID_LEN
                  || i >= a_len - END_AVOID_LEN
                  || j < END_AVOID_LEN
                  || j >= b_len - END_AVOID_LEN)
               indel_ring [n % (SIDE_BASES + 1)] = 3;  // too near end
           else if  (i < screen_lo)
               indel_ring [n % (SIDE_BASES + 1)] = 1;
             else
               indel_ring [n % (SIDE_BASES + 1)] = 2;  // screened

           quality_ring [n % (SIDE_BASES + 1)]
               = b_quality [j];

           j ++;
           n ++;
          }
     }
   while  (i <= a_len - END_AVOID_LEN && j <= b_len - END_AVOID_LEN)
     {
      match_mask &= INIT_MASK;
      bit = ((a [i] == b [j]) ? 1 : 0);
      match_mask = (match_mask << 1) | bit;

      quality_mask &= INIT_MASK;
      bit = ((a_quality [i] >= MIN_SIDE_QUALITY
                && b_quality [j] >= MIN_SIDE_QUALITY) ? 1 : 0);
      quality_mask = (quality_mask << 1) | bit;

      if  (i < END_AVOID_LEN
             || i >= a_len - END_AVOID_LEN
             || j < END_AVOID_LEN
             || j >= b_len - END_AVOID_LEN)
          indel_ring [n % (SIDE_BASES + 1)] = 3;  // too near end
      else if  (i < screen_lo)
          indel_ring [n % (SIDE_BASES + 1)] = 0;
        else
          indel_ring [n % (SIDE_BASES + 1)] = 2;  // screened
      quality_ring [n % (SIDE_BASES + 1)]
          = OVL_Min_int (a_quality [i], b_quality [j]);

      if  (quality_mask == QUAL_MATCH_MASK
             && indel_ring [(n + 1) % (SIDE_BASES + 1)] < 2)
          {
           if  (match_mask == SNP_MATCH_MASK)
               {
                switch  (indel_ring [(n + 1) % (SIDE_BASES + 1)])
                  {
                   case  0 :    // mismatch
#if  SHOW_SNP_DETAILS
                     fprintf (Out_Stream, "SNP  %4d  %4d  %.*s/%.*s\n",
                              i - SIDE_BASES, j - SIDE_BASES,
                              2 * SIDE_BASES + 1, a + i - 2 * SIDE_BASES,
                              2 * SIDE_BASES + 1, b + j - 2 * SIDE_BASES);
                     fprintf (Out_Stream,
                              "     i = %4d  j = %4d  qual_mask = %04o"
                              "  match_mask = %04o\n",
                              i, j, quality_mask, match_mask);
#endif
                     (* mismatch_ct) ++;
                     Global_Match_SNP_Bin
                         [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                     local_msnp_bin
                         [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                     break;
                   case  1 :    // insert
#if  SHOW_SNP_DETAILS
                     fprintf (Out_Stream, "SNP  %4d  %4d  %.*s-%.*s/%.*s\n",
                              i - SIDE_BASES, j - SIDE_BASES,
                              SIDE_BASES, a + i - 2 * SIDE_BASES + 1,
                              SIDE_BASES, a + i - SIDE_BASES + 1,
                              2 * SIDE_BASES + 1, b + j - 2 * SIDE_BASES);
                     fprintf (Out_Stream,
                              "     i = %4d  j = %4d  qual_mask = %04o"
                              "  match_mask = %04o\n",
                              i, j, quality_mask, match_mask);
#endif
                     (* indel_ct) ++;
                     Global_Indel_SNP_Bin
                         [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                     local_isnp_bin
                         [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                     break;
                   case  -1 :    // delete
#if  SHOW_SNP_DETAILS
                     fprintf (Out_Stream, "SNP  %4d  %4d  %.*s/%.*s-%.*s\n",
                              i - SIDE_BASES, j - SIDE_BASES,
                              2 * SIDE_BASES + 1, a + i - 2 * SIDE_BASES,
                              SIDE_BASES, b + j - 2 * SIDE_BASES + 1,
                              SIDE_BASES, b + j - SIDE_BASES + 1);
                     fprintf (Out_Stream,
                              "     i = %4d  j = %4d  qual_mask = %04o"
                              "  match_mask = %04o\n",
                              i, j, quality_mask, match_mask);
#endif
                     (* indel_ct) ++;
                     Global_Indel_SNP_Bin
                         [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                     local_isnp_bin
                         [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                     break;
                   case  2 :    // screened
#if  SHOW_SNP_DETAILS
                     fprintf (Out_Stream, "SNP  %4d  %4d  screened\n",
                              i - SIDE_BASES, j - SIDE_BASES);
#endif
                     break;
                   case  3 :   // too near end
#if  SHOW_SNP_DETAILS
                     fprintf (Out_Stream, "SNP  %4d  %4d  too near end\n",
                              i - SIDE_BASES, j - SIDE_BASES);
#endif
                     break;
                   default :
                     assert (FALSE);
                  }
               }
           else if  (match_mask == QUAL_MATCH_MASK)    // perfect match
               {
                (* all_match_ct) ++;
                Global_All_Match_Bin
                    [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                local_match_bin
                    [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
               }
             else
               {
                (* hi_qual_ct) ++;
                Global_Hi_Qual_Bin
                    [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
                local_other_bin
                    [quality_ring [(n + 1) % (SIDE_BASES + 1)] / 5] ++;
               }
          }

      if  (indel_ring [(n + 1) % (SIDE_BASES + 1)] != 3)
          {
           (* olap_bases) ++;
           if  (indel_ring [(n + 1) % (SIDE_BASES + 1)] != 2)
               (* unscreened_olap_bases) ++;
          }

      i ++;
      j ++;
      n ++;

      while  (i > screen_hi)
        {
         if  (WA -> screen_info . range [screen_sub] . last)
             screen_lo = screen_hi = INT_MAX;
           else
             {
              screen_sub ++;
              screen_lo = WA -> screen_info . range [screen_sub] . bgn;
              screen_hi = WA -> screen_info . range [screen_sub] . end;
             }
#if  SHOW_SNP_DETAILS
         fprintf (Out_Stream, "Screen range: %10d %10d\n",
                  screen_lo, screen_hi);
#endif
        }
     }

   Global_Olap_Ct += (* olap_bases);
   Global_Unscreened_Ct += (* unscreened_olap_bases);

   return;
  }
#endif



int  Sign
    (int a)

//  Return the algebraic sign of  a .

  {
   if  (a > 0)
       return  1;
   else if  (a < 0)
       return  -1;

   return  0;
  }


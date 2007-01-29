
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
/*********************************************************************
   Module:  AS_OVL
   Description:  Assembly Overlap Module.  Computes overlaps between
      pairs of DNA strings.
   Assumptions:  Input meets specifications in the ProtoIO documents
 *********************************************************************/

/* RCS info
 * $Id: AS_OVL_overlap.h,v 1.12 2007-01-29 20:41:17 brianwalenz Exp $
 * $Revision: 1.12 $
*/


#ifndef AS_OVL_OVERLAP_H
#define AS_OVL_OVERLAP_H

//
// Component:
//   AS_OVL_overlap
//
//   Art. Delcher w/ much assistance form Saul Kravitz
//   Last revised:  7 Jan 99
//
//   Reference:  "Overlap Module.rtf"
// 
//  Description:
// 
//   This component reads a stream of fragments, computes which of
//   them overlap, i.e., match (within error tolerance) to the
//   end of a fragment in both directins.  The output is a stream
//   of overlap messages.  The input fragments are also saved in
//   a fragment store.
//      
// 
// Memory Usage:
// 
//   This component allocates its own memory.
// 
// Interface:
// 
//   Uses the persistent fragment store created and maintained by
//   AS_OVL_driver.c .  Outputs overlap messages
//   to proto-I/O stream.
// 
// Design:
// 
//   Uses 2 streams of fragments.  From the New stream, a
//   hash-table index is built for each k-mer occuring in those
//   fragments.  The hash-table entries point to the start position
//   of that k-mer in a fragment.  When a k-mer occurs in multiple
//   fragments, there is a reference in  Next_Ref  to the next
//   fragment position with the same k-mer.   Next_Ref  is currently
//   a parallel array alongside the actual fragment data.
//   References are stored as  String_Ref  structures containing the
//   fragment number and offset within that fragment of the
//   referenced position.  The start position of each fragment is
//   contained in the global array  String_Start .
// 
//   Hash entries are in buckets.  Each entry also has a  Check  byte.
//   The  Check  byte is a small secondary hash key to minimize
//   unnecessary accesses to the main fragment array.  When a
//   k-mer is hashed, the primary hash key determines the bucket
//   its in.  The  Check  byte is stored with it.  When a lookup
//   is done on a k-mer, if the  Check  byte does not match, the
//   entry does not match.  When the  Check  byte does match, the
//   k-mer must still be verified by following the reference to
//   the actual fragment data.
// 
//   Collisions are handled by secondary hashing.  A second hash
//   value is computed from the k-mer, and this value is used as
//   the distance to go to the next bucket (modulo the hash table
//   size).
// 
//   The second stream of fragments is the Old stream.  For each
//   fragment in it, all the k-mers are looked up in the hash
//   table.  Each match is stored for later processing.  Matches
//   are grouped together first by fragment number (in a small hash
//   table), and then coallesced into contiguous exact-match
//   regions of length longer than  k .  These exact-match regions
//   are kept as a linked list associated with each fragment.
// 
//   For each fragment the longest exact-match region is found.
//   Then an edit-distance calculation (aka Smith-Waterman)
//   is done in the forward direction to see how far that match
//   can be extended within the error-bound.  If the match
//   extends to the end of either fragment, then a similar
//   edit-distance calculation is done starting at the exact-match
//   region, going in the reverse direction, to see if the match
//   extends to the end of either fragment.  If the match reaches
//   the end in both directions, then the alignment is computed
//   (as a sequence of Delta values, see the  .rtf  document) and
//   output.
//
//   To speed overlap calculations, this edit distance calculation
//   is now "banded" in the sense that it only pursues alignments
//   that are likely to reach the end of a fragment.
//   Specifically, suppose a particular alignment "thread"
//   includes  e  errors in an alignment of length  L .  If the
//   binomial probability of  e  or more errors in  L  independent
//   trials (with ERROR_RATE as probability of each error) is less
//   than EDIT_DIST_PROB_BOUND, then this alignment thread is
//   abandoned.  For each  e  value, the minimum  L  to pursue the
//   alignment thread is kept in the global array  Edit_Match_Limit .
//   If these values are all set to zero, then banding is
//   eliminated.
//
//   Each fragment in the Old stream is processed twice.  Once in
//   the input orientation, and then again after being
//   reverse-complemented.
// 
// Limitations:
// 
//   The number of fragments in the New stream stored in the hash
//   table is assumed to be at most  Max_Hash_Strings .
// 
//   MAX_ERRORS  must be at least  ERROR_RATE * AS_READ_MAX_LEN .
//
//   Max_Hash_Strings  must be  < 2^21 .
// 
// Status:
// 
//   Full overlaps are working.
//   Branch points are removed.
//   No quality values are being used.
//   UR detector screened regions are ignored in constructing the
//     hash table.
//   SuperRepeats are not filtered.
//   Nothing special is done to detect microsat/polymorphic overlaps.
//   No special handling is given to periodic k-mers/overlap regions.
// 
// Architecture and Dependencies:
// 
//   AS_OVL_delcher.h   Generic  #includes  and simple function prototypes
//   AS_OVL_delcher.c   Corresponding function definitions
//   AS_OVL_driver.c    Routines to create and update fragment store
//                      Depends on  AS_PER , AS_UTL  and  AS_MSG
//   AS_OVL_overlap.h   This file
//   AS_OVL_overlap.c   main  and functions to compute overlaps
//



/***************************************************************/
/* Constant definitions; Macro definitions; type definitions */
/***************************************************************/

#define  FOR_CARL_FOSLER          0
    //  If  1  outputs a single line per overlap in ASCII form
    //  Also restricts overlaps more
#define  HASH_KMER_SKIP           0
    //  Skip this many kmers between the kmers put into the hash
    //  table.  Setting this to 0 will make every kmer go
    //  into the hash table.


#define  BAD_WINDOW_LEN           50
    //  Length of window in which to look for clustered errors
    //  to invalidate an overlap
#define  BAD_WINDOW_VALUE         (8 * QUALITY_CUTOFF)
    //  This many or more errors in a window of  BAD_WINDOW_LEN
    //  invalidates an overlap
#define  CHECK_MASK              0xff
    //  To set Check field in hash bucket

#define  DEFAULT_HI_HIT_LIMIT    INT_MAX
    //  Any kmer in hash table with at least this many hits
    //  cannot initiate an overlap.  Can be changed on the command
    //  line with  -K  option


#define  DEFAULT_BRANCH_MATCH_VAL    ( ERR_FRACTION_IN_AS_GLOBAL_H / (1.+ERR_FRACTION_IN_AS_GLOBAL_H) )
#define  PARTIAL_BRANCH_MATCH_VAL    DEFAULT_BRANCH_MATCH_VAL
    //  Value to add for a match in finding branch points.
    //  ALH: Note that AS_READ_ERROR_RATE also affects what overlaps get found
    //  ALH: Scoring seems to be unusual: given an alignment
    //  of length l with k mismatches, the score seems to be
    //  computed as l + k * error value and NOT (l-k)*match+k*error
    // 
    //  I.e. letting x := DEFAULT_BRANCH_MATCH_VAL,
    //  the max mismatch fraction p to give a non-negative score
    //  would be p = x/(1-x); conversely, to compute x for a
    //  goal p, we have x = p/(1+p).  E.g. 
    //  for p=0.06, x = .06/(1.06) = .0566038; 
    //  for p=0.35, x = .35/(1.35) = .259259
    //  for p=0.2, x = .2/(1.2) = .166667
    //  for p=0.15, x = .15/(1.15) = .130435
    //
    //  Value was for 6% vs 35% error discrimination.
    //  Converting to integers didn't make it faster.
    //  Corresponding error value is this value minus 1.0
#define  DEFAULT_NUM_PTHREADS    4
    //  Number of threads to use if not overridden by  -t
    //  command-line option
#define  DELETED_FRAG            2
    //  Indicates fragment was marked as deleted in the store
#define  DISPLAY_WIDTH           60
    //  Number of characters per line when displaying sequences
#define  DONT_KNOW_CHAR          'n'
    //  Unknown input character
#define  EDIT_DIST_PROB_BOUND    1e-4
    //  Probability limit to "band" edit-distance calculation
    //  Determines  NORMAL_DISTRIB_THOLD
#define  ENTRIES_PER_BUCKET      21
    //  In main hash table.  Recommended values are 21, 31 or 42
    //  depending on cache line size.
#define  ERRORS_FOR_FREE         1
    //  The number of errors that are ignored in setting probability
    //  bound for terminating alignment extensions in edit distance
    //  calculations
#define  FRAG_OLAP_LIMIT         INT_MAX
    //  Most overlaps any single old fragment can have with a single
    //  set of new fragments in the hash table.
#define  HASH_CHECK_MASK         0x1f
    //  Used to set and check bit in Hash_Check_Array
    //  Change if change  Check_Vector_t
#define  HASH_EXPANSION_FACTOR   1.4
    //  Hash table size is >= this times  MAX_HASH_STRINGS
#define  HASH_MASK               ((1 << Hash_Mask_Bits) - 1)
    //  Extract right Hash_Mask_Bits bits of hash key
#define  HASH_TABLE_SIZE         (1 + HASH_MASK)
    //  Number of buckets in hash table
#define  HIGHEST_KMER_LIMIT      255
    //  If  Hi_Hit_Limit  is more than this, it's ignored
#define  HOPELESS_MATCH          90
    //  A string this long or longer without an exact kmer
    //  match is assumed to be hopeless to find a match
    //  within the error threshold
#define  IID_GAP_LIMIT           100
    //  When using a list of fragment IID's, gaps between
    //  IID's this long or longer force a new load partial
    //  store to be done
#define  INIT_MATCH_NODE_SIZE    10000
    //  Initial number of nodes to hold exact matches
#define  INIT_SCREEN_MATCHES     50
    //  Initial number of screen-match entries per fragment
#define  INIT_STRING_OLAP_SIZE   5000
    //  Initial number of different New fragments that
    //  overlap a single Old fragment
#define  K_MER_STEP          1
    //  1 = every k-mer in search
    //  2 = every other k-mer
    //  3 = every third, ...
    //  Used to skip some k-mers in finding matches
#define  MAX_BRANCH_COUNT        UCHAR_MAX
    //  The largest branch point count before an overflow
#define  MAX_DISTINCT_OLAPS      3
    //  Most possible really different overlaps (i.e., not
    //  just shifts from periodic regions) between 2 fragments
    //  in a given orientation.  For fragments of approximately
    //  same size, should never be more than 2.
#define  MAX_ERROR_RATE          AS_GUIDE_ERROR_RATE
    //  The largest error allowed in overlaps
#define  MAX_FRAG_LEN            2048
    //  The longest fragment allowed
#define  MAX_ERRORS              (1 + (int) (MAX_ERROR_RATE * MAX_FRAG_LEN))
    //  Most errors in any edit distance computation
#define  MAX_FRAGS_PER_THREAD    500
    //  The number of fragments each parallel thread tries to
    //  process in a "round"

#ifdef  CONTIG_OVERLAPPER_VERSION
#define  STRING_NUM_BITS         11
#else
#define  STRING_NUM_BITS         19
#endif
    //  Number of bits used to store the string number in the
    //  hash table
#define  MAX_STRING_NUM          ((1 << STRING_NUM_BITS) - 1)
   //   Largest string number that can fit in the hash table

#define  OFFSET_BITS             (30 - STRING_NUM_BITS)
    //  Number of bits used to store lengths of strings stored
    //  in the hash table
#define  OFFSET_MASK             ((1 << OFFSET_BITS) - 1)
    //  Mask used to extract bits to put in  Offset  field

//  xxxx constants were here

#ifdef  CONTIG_OVERLAPPER_VERSION
#define  EXPECTED_STRING_LEN     (AS_READ_MAX_LEN / 2)
#else
#define  EXPECTED_STRING_LEN     (MAX_FRAG_LEN / 2)
#endif
#define  INITIAL_DATA_LEN        (EXPECTED_STRING_LEN * Max_Hash_Strings)
    //  The number of bytes to allocate initially for hash-table sequence

#define  MAX_OLD_BATCH_SIZE      100000
    //  Most old fragments to read at a time from the fragment store
#define  MAX_HASH_LOAD           0.60
    //  Maximum portion of hash table allowed to be filled
#define  MAX_LINE_LEN            1000
    //  Maximum input line when reading FASTA file
#define  MAX_NAME_LEN            500
    //  Longest file name allowed
#define  MEMORY_EXPANSION_FACTOR  1.2
    //  How much to grow malloc'ed space if necessary
#define  MIN_BRANCH_END_DIST     20
    //  Branch points must be at least this many bases from the
    //  end of the fragment to be reported
#if ERR_MODEL_IN_AS_GLOBAL_H > 6
  #define  MIN_BRANCH_TAIL_SLOPE   1.0
#else
  #define  MIN_BRANCH_TAIL_SLOPE   0.20
#endif
    //  Branch point tails must fall off from the max by at least
    //  this rate
#define  MIN_CALC_KMER           4
    //  When calculating the  Hi_Hit_Limit  based on genome length, etc,
    //  don't set it below this
#define  MIN_INTERSECTION        10
    //  Minimum number of bases periodic overlaps must have in common
#define  MIN_OLAP_LEN            40
    //  Minimum length of match region to be worth reporting
#define  MIN_OLAP_OUTSIDE_SCREEN 30
    //  Minimum number of bases outside of screened regions
    //  to be a reportable overlap.  Entire overlap (including screened
    //  portion) must still be  MIN_OLAP_LEN .
#define  NORMAL_DISTRIB_THOLD    3.62
    //  Determined by  EDIT_DIST_PROB_BOUND
#define  OUTPUT_OVERLAP_DELTAS   0
    //  If true include delta-encoding of overlap alignment
    //  in overlap messages.  Otherwise, omit them.
#define  PROBE_MASK              0x3e
    //  Used to determine probe step to resolve collisions
#define  QUALITY_BASE_CHAR       '0'
    //  Quality values were added to this to create printable
    //  characters
#define  QUALITY_CUTOFF          20
    //  Regard quality values higher than this as equal to this
    //  for purposes of finding bad windows
#define  SCRIPT_NAME             "lsf-ovl"
    //  Default name of script produced by  make-ovl-script
#define  SHIFT_SLACK  1
    // Allow to be off by this many bases in combining/comparing alignments
#define  STAT_FILE_NAME          "overlap.stats"
    //  Where statistics info is written if  SHOW_STATS != 0
    //  This version is closest to standalone version olapv109.c

#ifdef  CONTIG_OVERLAPPER_VERSION
#define  STRING_OLAP_SHIFT       8
#else
#define  STRING_OLAP_SHIFT       10
#endif
    //  To compute hash function into the String_Olap hash table.

#define  STRING_OLAP_MODULUS     (1 << STRING_OLAP_SHIFT)
    //  The size of the  String_Olap  hash table.  The rest of
    //  the space is used for chaining.  This number should be
    //  relatively small to reflect the number of fragments a
    //  given fragment has exact matches with.
#define  STRING_OLAP_MASK        (STRING_OLAP_MODULUS - 1)
    //  To compute hash function into the String_Olap hash table.
#define  THREAD_STACKSIZE        (16 * 512 * 512)
//#define  THREAD_STACKSIZE        (2 * 512 * 512)
    //  The amount of stack space to allocate to each thread.
#define  TOTAL_SCREEN_FILE_NAME          "total-screen.lst"
    //  Name of file where list of totally screened fragments put
    //  if  LIST_TOTALLY_SCREENED != 0
#define  VALID_FRAG              1
    //  Indicates fragment was valid in the fragment store
#define  WINDOW_SCREEN_OLAP      10
    //  Amount by which k-mers can overlap a screen region and still
    //  be added to the hash table.

#if  FOR_CARL_FOSLER
  #define  WINDOW_SIZE             28
#else
  #if ERR_MODEL_IN_AS_GLOBAL_H > 6
    #define  WINDOW_SIZE             14
  #else
    #define  WINDOW_SIZE             22
  #endif
#endif
    //  Length of segments hashed, i.e., the  k  value in k-mer.
    //  There must be an exact match of this length or more to
    //  find a match.
#define  KMER_LEN  WINDOW_SIZE
    //  A more descriptive alias
#define  MAX_EXTRA_SUBCOUNT        (AS_FRAG_MAX_LEN / KMER_LEN)


#define  HSF1                    (WINDOW_SIZE - (Hash_Mask_Bits / 2))
    //  First shift value for HASH_FUNCTION
#define  HSF2                    (2 * WINDOW_SIZE - Hash_Mask_Bits)
    //  Second shift value for HASH_FUNCTION
#define  SV1                     (HSF1 + 2)
#define  SV2                     ((HSF1 + HSF2) / 2)
#define  SV3                     (HSF2 - 2)
#define  HASH_FUNCTION(k)        (((k) ^ ((k) >> HSF1) ^ ((k) >> HSF2)) & HASH_MASK)
    //  Gives subscript in hash table for key  k
#define  HASH_CHECK_FUNCTION(k)  (((k) ^ ((k) >> SV1) ^ ((k) >> SV2)) & HASH_CHECK_MASK)
    //  Gives bit position to see if key could be in bucket
#define  KEY_CHECK_FUNCTION(k)   (((k) ^ ((k) >> SV1) ^ ((k) >> SV3)) & CHECK_MASK)
    //  Gives bit pattern to see if key could match
#define  PROBE_FUNCTION(k)       ((((k) ^ ((k) >> SV2) ^ ((k) >> SV3)) & PROBE_MASK) | 1)
    //  Gives secondary hash function.  Force to be odd so that will be relatively
    //  prime wrt the hash table size, which is a power of 2.

#define  INPUT_FILENAME_EXTENSION   ".urc"
#define  OUTPUT_FILENAME_EXTENSION  ".ovl"
    //  For proto-I/O files

#define  ANALYZE_HITS             0
#define  DEBUG                    0
#define  DEBUG7                   0
#define  DEBUG8                   0
#define  DO_OLAP_ALIGN_PROFILE    0
    //  If  1  take overlap alignments and count frequency of occurrence
    //  of bases at each base of hash-table fragments
    //  ** Need to set  SHOW_OVERLAPS  to  1  also **
#define  DO_KMER_HITS_PROFILE     0
    //  If  1  show hash-table kmer frequency at each base of fragments
    //  read in separately.
    //  ** Need to set  ANALYZE_HITS  to  1  also **
#define  LIST_TOTALLY_SCREENED    0
    //  If  1  creates a file named  TOTAL_SCREEN_FILE_NAME  with IID's
    //  of frag's that were totally screened when being put in the
    //  hash table
#define  SCREEN_CHECK_ONLY        0
    //  If  1  only check screening on hash table frags.  Don't
    //  really put them in the hash table and don't overlap anything
    //  against them
#define  SHOW_BAD_WINDOWS         0
#define  SHOW_HI_HIT_KMERS        0
    //  If  1  dump frequent kmers from hash table after building it
    //  and then exit immediately.  Need to set  ANALYZE_HITS  to  1
    //  and  DO_KMER_HITS_PROFILE  also
#define  SHOW_OVERLAPS            0
#define  SHOW_PROGRESS            0
#define  SHOW_SNPS                0
#define  SHOW_STATS               0
#define  SHOW_THREAD_PROGRESS     0
#define  DEBUGSTREAM              0
#define  TEST_DP_BPT_CODE         0
#define  UPDATE_FRAG_STORE        1     // Temporarily set to 0 for contigs
#define  USE_SOURCE_FIELD         0

#include <pthread.h>

typedef  FragStreamHandle  Frag_Stream;
typedef  FragStoreHandle  FragStore;
typedef  FILE *  Output_Stream;
typedef  FILE *  Input_Stream;


typedef  int64  Ext_Frag_ID_t;
typedef  int32  Int_Frag_ID_t;

typedef  enum Direction_Type
  {FORWARD, REVERSE}  Direction_t;

typedef  struct Match_Node
  {
   int32  Offset;              // To start of exact match in  hash-table frag
   int32  Len;                 // Of exact match
   int32  Start;               // Of exact match in current (new) frag
   int32  Next;                // Subscript of next match in list
  }  Match_Node_t;

typedef  struct String_Olap_Node
  {
   uint32  String_Num;                // Of hash-table frag that
                                      //   have exact match with
   int32  Match_List;                 // Subscript of start of list
                                      //   of exact matches
   double  diag_sum;                  // Sum of diagonals of all k-mer
                                      //   matches to this frag
   int32  diag_ct;                    // Count of all k-mer matches to this
                                      //   frag
#if  SHOW_STATS                          
   int32  Kmer_Hits;
#endif
   signed int  Next : 29;             // Next match if this is a collision
   unsigned  Full : 1;
   unsigned  consistent : 1;
  }  String_Olap_t;

typedef  struct
  {
   unsigned  bgn : 15;
   unsigned  end : 15;
   unsigned  last : 1;
  }  Screen_Range_t;

typedef  struct
  {
   Screen_Range_t  * range;
   int  match_len;                   // Number of entries allocated in  match
   int  num_matches;                 // Number of entries occupied in match
   int  left_end_screened;
   int  right_end_screened;
  }  Screen_Info_t;

//  The following structure holds what used to be global information, but
//  is now encapsulated so that multiple copies can be made for multiple
//  parallel threads.

typedef  struct Work_Area
  {
   int  Left_Delta_Len;
   int  Left_Delta [MAX_ERRORS];
   int  Right_Delta_Len;
   int  Right_Delta [MAX_ERRORS];
   int  Edit_Space [(MAX_ERRORS + 4) * MAX_ERRORS];
   int  * Edit_Array [MAX_ERRORS];
   int  * Edit_Match_Limit;
   int  * Error_Bound;
   String_Olap_t  * String_Olap_Space;
   int32  String_Olap_Size;
   int32  Next_Avail_String_Olap;
   Match_Node_t  * Match_Node_Space;
   int32  Match_Node_Size;
   int32  Next_Avail_Match_Node;
   int32  A_Olaps_For_Frag, B_Olaps_For_Frag;
            //  Counts the number of overlaps for each fragment.  Cuts off
            //  overlaps above a limit.
   Frag_Stream  stream_segment;
   Screen_Info_t  screen_info;
   int  status;
   ReadStructp  myRead;
   FragType  curr_frag_type;
   int  thread_id;
  }  Work_Area_t;


extern int  Hash_Mask_Bits;
extern int  Max_Hash_Data_Len;
extern int  Max_Hash_Strings;
extern int  Max_Frags_In_Memory_Store;

extern int  Contig_Mode;
extern uint32  Last_Hash_Frag_Read;
extern int  LSF_Mode;
extern int  Lo_Hash_Frag, Hi_Hash_Frag;
extern int  Lo_Old_Frag, Hi_Old_Frag;
extern int  Num_PThreads;
extern int64  Olap_Ct;
extern int  Table_Ct;
extern clock_t  Start_Time, Stop_Time;
    //  Used to output info to  stderr  to track progress

extern Screen_Range_t  * Screen_Space;
extern int  Screen_Space_Size;
extern int  * Screen_Sub;
    //  To keep track of screened sections of fragments in the
    //  hash table.  K-mers from screened sections are *NOT*
    //  put into the hash table.   For fragment  i  in the hash
    //  table, if  j = Screen_Sub [i] , then continguous entries
    //  Screen_Space [j .. ]  hold the regions screened from
    //  that fragment.

extern char  Sequence_Buffer [];
extern char  Quality_Buffer [];
    //  Used to read fragments and to build hash table which
    //  are done single threaded.

extern FragStore  OldFragStore;
extern FragStore  BACtigStore;
extern DistStore  OldDistStore;
extern char  * BACtig_Store_Path;
extern char  * Frag_Store_Path;
extern Input_Stream  In_Stream;
extern Output_Stream  Out_Stream;
    //  To handle I/O
extern uint32  * IID_List;
extern int  IID_List_Len;

extern pthread_mutex_t  FragStore_Mutex;
extern pthread_mutex_t  Write_Proto_Mutex;
extern pthread_mutex_t  Log_Msg_Mutex;

//
//  Prototypes of functions used by both  AS_OVL_driver.c  and
//  AS_OVL_overlap.c .
//

int  Build_Hash_Index
    (FragStreamHandle stream, int32 First_Frag_ID, ReadStructp myRead);
void  Coalesce_Screen_Info
    (Screen_Range_t space [], int lo, int * hi);
void  Copy_Branches_To_Store
    (Frag_Stream stream, FragStore store, ReadStructp myRead);
void  Dump_Stored_Branch_Info
    (Frag_Stream stream, ReadStructp myRead);
void  Initialize_Work_Area
    (Work_Area_t *, int);
int  OVL_Max_int
    (int, int);
int  OVL_Min_int
    (int, int);
void  Output_Hash_Table_Branch_Pts
    (void);
void  Output_High_Hit_Frags
    (void);
int  OverlapDriver
    (int noOverlaps, int argc, char **argv);
void  Process_Overlaps
    (FragStreamHandle stream, Work_Area_t *);
void  Profile_Hits
    (void);
int  Sign
    (int);
void  stripWhiteSpace
    (char * target, char * source, int maxlen);

#if  USE_SOURCE_FIELD
extern FILE  * Source_Log_File;
#endif
#if  ANALYZE_HITS
extern FILE  * High_Hits_File;
#endif


#endif


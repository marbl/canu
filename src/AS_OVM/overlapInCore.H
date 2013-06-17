
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

//static const char *rcsid = "$Id: overlapInCore.H,v 1.5 2013-01-11 11:11:04 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <pthread.h>

#include "AS_UTL_fileIO.h"
#include "AS_UTL_reverseComplement.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlapStore.h"


#ifndef OVERLAPINCORE_H
#define OVERLAPINCORE_H


////////////////////////////////////////
//
//  AS_OVL_overlap.h BEGINS
//
////////////////////////////////////////

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

#define  DEFAULT_BRANCH_MATCH_VAL    ( AS_OVL_ERROR_RATE / (1.+AS_OVL_ERROR_RATE) )
#define  PARTIAL_BRANCH_MATCH_VAL    DEFAULT_BRANCH_MATCH_VAL
//  Value to add for a match in finding branch points.
//  ALH: Note that AS_OVL_ERROR_RATE also affects what overlaps get found
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

#define  MAX_ERRORS              (1 + (int) (AS_OVL_ERROR_RATE * AS_READ_MAX_NORMAL_LEN))
//  Most errors in any edit distance computation
//  THIS VALUE IS KNOWN ONLY AT RUN TIME!

#define  EXPECTED_STRING_LEN     (AS_READ_MAX_NORMAL_LEN / 2)
#define  INITIAL_DATA_LEN        (EXPECTED_STRING_LEN * Max_Hash_Strings)
//  The number of bytes to allocate initially for hash-table sequence

#define  MAX_LINE_LEN            1000
//  Maximum input line when reading FASTA file

#define  MAX_NAME_LEN            500
//  Longest file name allowed

#define  MEMORY_EXPANSION_FACTOR  1.2
//  How much to grow malloc'ed space if necessary

#define  MIN_BRANCH_END_DIST     20
//  Branch points must be at least this many bases from the
//  end of the fragment to be reported

#define  MIN_BRANCH_TAIL_SLOPE   ((AS_OVL_ERROR_RATE > 0.06) ? 1.0 : 0.20)
//  Branch point tails must fall off from the max by at least
//  this rate
//  THIS VALUE IS KNOWN ONLY AT RUN TIME!

#define  MIN_CALC_KMER           4
//  When calculating the  Hi_Hit_Limit  based on genome length, etc,
//  don't set it below this

#define  MIN_INTERSECTION        10
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
//  As of 6 Oct 2008, support for overlap deltas has been removed.
//  However, there are enough remnants in AS_MSG to output them.
//  Just enabling OUTPUT_OVERLAP_DELTAS will not compile; see
//  AS_MSG_USE_OVL_DELTA in AS_MSG.

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

#define  STRING_OLAP_SHIFT       8
//  To compute hash function into the String_Olap hash table.

#define  STRING_OLAP_MODULUS     (1 << STRING_OLAP_SHIFT)
//  The size of the  String_Olap  hash table.  The rest of
//  the space is used for chaining.  This number should be
//  relatively small to reflect the number of fragments a
//  given fragment has exact matches with.

#define  STRING_OLAP_MASK        (STRING_OLAP_MODULUS - 1)
//  To compute hash function into the String_Olap hash table.

#define  THREAD_STACKSIZE        (16 * 512 * 512)
//  The amount of stack space to allocate to each thread.

#define  VALID_FRAG              1
//  Indicates fragment was valid in the fragment store

//#define  WINDOW_SCREEN_OLAP      10
//  Amount by which k-mers can overlap a screen region and still
//  be added to the hash table.

#define  MAX_EXTRA_SUBCOUNT        (AS_READ_MAX_NORMAL_LEN / Kmer_Len)


#define  HASH_FUNCTION(k)        (((k) ^ ((k) >> HSF1) ^ ((k) >> HSF2)) & HASH_MASK)
//  Gives subscript in hash table for key  k

#define  HASH_CHECK_FUNCTION(k)  (((k) ^ ((k) >> SV1) ^ ((k) >> SV2)) & HASH_CHECK_MASK)
//  Gives bit position to see if key could be in bucket

#define  KEY_CHECK_FUNCTION(k)   (((k) ^ ((k) >> SV1) ^ ((k) >> SV3)) & CHECK_MASK)
//  Gives bit pattern to see if key could match

#define  PROBE_FUNCTION(k)       ((((k) ^ ((k) >> SV2) ^ ((k) >> SV3)) & PROBE_MASK) | 1)
//  Gives secondary hash function.  Force to be odd so that will be relatively
//  prime wrt the hash table size, which is a power of 2.



typedef  enum Direction_Type {
  FORWARD,
  REVERSE
} Direction_t;

typedef  struct Match_Node {
  int32  Offset;              // To start of exact match in  hash-table frag
  int32  Len;                 // Of exact match
  int32  Start;               // Of exact match in current (new) frag
  int32  Next;                // Subscript of next match in list
}  Match_Node_t;

typedef  struct String_Olap_Node {
  uint32  String_Num;                // Of hash-table frag that have exact match with
  int32  Match_List;                 // Subscript of start of list of exact matches
  double  diag_sum;                  // Sum of diagonals of all k-mer matches to this frag
  int32  diag_ct;                    // Count of all k-mer matches to this frag
  signed int  Next : 29;             // Next match if this is a collision
  unsigned  Full : 1;
  unsigned  consistent : 1;
}  String_Olap_t;

//  The following structure holds what used to be global information, but
//  is now encapsulated so that multiple copies can be made for multiple
//  parallel threads.

typedef  struct Work_Area {
  int    Left_Delta_Len;
  int  * Left_Delta;
  int    Right_Delta_Len;
  int  * Right_Delta;
  int  * Delta_Stack;
  int  * Edit_Space;
  int  ** Edit_Array;
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
  gkStream  *stream_segment;

  int  left_end_screened;
  int  right_end_screened;

  int  status;
  gkFragment  myRead;
  int  thread_id;

  //  Instead of outputting each overlap as we create it, we
  //  buffer them and output blocks of overlaps.
  uint64         overlapsLen;
  uint64         overlapsMax;
  OVSoverlap   *overlaps;

  //  Various stats that used to be global and updated whenever we
  //  output an overlap or finished processing a set of hits.
  //  Needed a mutex to update.
  uint64         Total_Overlaps;
  uint64         Contained_Overlap_Ct;
  uint64         Dovetail_Overlap_Ct;

  uint64         Kmer_Hits_Without_Olap_Ct;
  uint64         Kmer_Hits_With_Olap_Ct;
  uint64         Multi_Overlap_Ct;
}  Work_Area_t;




typedef  uint32  Check_Vector_t;
// Bit vector to see if hash bucket could possibly contain a match

typedef  enum Overlap_Type
  {NONE, LEFT_BRANCH_PT, RIGHT_BRANCH_PT, DOVETAIL}  Overlap_t;

typedef  struct Olap_Info {
  int  s_lo, s_hi;
  int  t_lo, t_hi;
  double  quality;
  int  delta [AS_READ_MAX_NORMAL_LEN+1];  //  needs only MAX_ERRORS
  int  delta_ct;
  int  s_left_boundary, s_right_boundary;
  int  t_left_boundary, t_right_boundary;
  int  min_diag, max_diag;
}  Olap_Info_t;




typedef  uint64   String_Ref_t;

#define  BIT_EMPT  62
#define  BIT_LAST  63

extern uint32 STRING_NUM_BITS;
extern uint32 OFFSET_BITS;

extern uint64 TRUELY_ZERO;
extern uint64 TRUELY_ONE;

extern uint64 STRING_NUM_MASK;
extern uint64 OFFSET_MASK;

extern uint64 MAX_STRING_NUM;



//
//  [ Last (1) ][ Empty (1) ][ Offset (11) ][ StringNum (19) ]
//

#define getStringRefStringNum(X)      (((X)                   ) & STRING_NUM_MASK)
#define getStringRefOffset(X)         (((X) >> STRING_NUM_BITS) & OFFSET_MASK)
#define getStringRefEmpty(X)          (((X) >> BIT_EMPT       ) & TRUELY_ONE)
#define getStringRefLast(X)           (((X) >> BIT_LAST       ) & TRUELY_ONE)

#define setStringRefStringNum(X, Y)   ((X) = (((X) & ~(STRING_NUM_MASK                   )) | ((Y))))
#define setStringRefOffset(X, Y)      ((X) = (((X) & ~(OFFSET_MASK     << STRING_NUM_BITS)) | ((Y) << STRING_NUM_BITS)))
#define setStringRefEmpty(X, Y)       ((X) = (((X) & ~(TRUELY_ONE      << BIT_EMPT       )) | ((Y) << BIT_EMPT)))
#define setStringRefLast(X, Y)        ((X) = (((X) & ~(TRUELY_ONE      << BIT_LAST       )) | ((Y) << BIT_LAST)))


typedef  struct Hash_Bucket {
  String_Ref_t  Entry [ENTRIES_PER_BUCKET];
  unsigned char  Check [ENTRIES_PER_BUCKET];
  unsigned char  Hits [ENTRIES_PER_BUCKET];
  int16  Entry_Ct;
}  Hash_Bucket_t;

typedef  struct Hash_Frag_Info {
  uint32  length             : 30;
  uint32  lfrag_end_screened : 1;
  uint32  rfrag_end_screened : 1;
}  Hash_Frag_Info_t;


extern int64   Bad_Short_Window_Ct;
extern int64   Bad_Long_Window_Ct;
extern double  Branch_Match_Value;
extern double  Branch_Error_Value;
extern char  * Data;
extern char  * Quality_Data;
extern size_t  Data_Len;
extern bool  Doing_Partial_Overlaps;
extern size_t  Extra_Data_Len;
extern uint64  Extra_Ref_Ct;
extern String_Ref_t  * Extra_Ref_Space;
extern uint64  Extra_String_Ct;
extern uint64  Extra_String_Subcount;
extern uint64  Frag_Olap_Limit;
extern Check_Vector_t  * Hash_Check_Array;
extern uint64  Hash_String_Num_Offset;
extern Hash_Bucket_t  * Hash_Table;
extern bool  Ignore_Clear_Range;
extern uint64  Kmer_Hits_With_Olap_Ct;
extern uint64  Kmer_Hits_Without_Olap_Ct;
extern int32  Min_Olap_Len;
extern uint64  Multi_Overlap_Ct;
extern String_Ref_t  * Next_Ref;
extern uint64  String_Ct;
extern Hash_Frag_Info_t  * String_Info;
extern int64  * String_Start;
extern uint32  String_Start_Size;
extern bool  Unique_Olap_Per_Pair;
extern size_t  Used_Data_Len;
extern bool  Use_Hopeless_Check;
extern bool  Use_Window_Filter;
extern int32  Read_Edit_Match_Limit [AS_READ_MAX_NORMAL_LEN];
extern int32  Guide_Edit_Match_Limit [AS_READ_MAX_NORMAL_LEN];
extern int32  Read_Error_Bound [AS_READ_MAX_NORMAL_LEN + 1];
extern int32  Guide_Error_Bound [AS_READ_MAX_NORMAL_LEN + 1];
extern double  Branch_Cost [AS_READ_MAX_NORMAL_LEN + 1];
extern int32  Bit_Equivalent [256];
extern int32  Char_Is_Bad [256];
extern uint64  Hash_Entries;
extern uint64  Total_Overlaps;
extern uint64  Contained_Overlap_Ct;
extern uint64  Dovetail_Overlap_Ct;
extern uint32 minLibToHash;
extern uint32 maxLibToHash;
extern uint32 minLibToRef;
extern uint32 maxLibToRef;
extern uint64  Kmer_Len;
extern uint64  HSF1;
extern uint64  HSF2;
extern uint64  SV1;
extern uint64  SV2;
extern uint64  SV3;
extern uint32  Hash_Mask_Bits;
extern double  Max_Hash_Load;
extern uint32  Max_Hash_Strings;
extern uint64  Max_Hash_Data_Len;
extern AS_IID  Max_Frags_In_Memory_Store;
extern AS_IID  Last_Hash_Frag_Read;
extern AS_IID  Lo_Hash_Frag;
extern AS_IID  Hi_Hash_Frag;
extern AS_IID  Lo_Old_Frag;
extern AS_IID  Hi_Old_Frag;
extern uint32  Num_PThreads;
extern char  Sequence_Buffer [2 * AS_READ_MAX_NORMAL_LEN];
extern char  Quality_Buffer [2 * AS_READ_MAX_NORMAL_LEN];
extern FILE  * Kmer_Skip_File;
extern BinaryOverlapFile  *Out_BOF;
extern gkStore  *OldFragStore;
extern char  * Frag_Store_Path;
extern pthread_mutex_t  FragStore_Mutex;
extern pthread_mutex_t  Write_Proto_Mutex;
extern AS_IID  First_Hash_Frag;
extern AS_IID  Last_Hash_Frag;
extern gkFragment  myRead;
extern AS_IID  Frag_Segment_Lo;
extern AS_IID  Frag_Segment_Hi;




#define Sign(a) ( ((a) > 0) - ((a) < 0) ) 


int
Read_Next_Frag(char frag [AS_READ_MAX_NORMAL_LEN + 1],
               char quality [AS_READ_MAX_NORMAL_LEN + 1],
               gkStream *stream,
               gkFragment *myRead,
               AS_IID * last_frag_read,
               uint32 minLibToRead, uint32 maxLibToRead);


void
Output_Overlap(AS_IID S_ID, int S_Len, Direction_t S_Dir,
               AS_IID T_ID, int T_Len, Olap_Info_t * olap,
               Work_Area_t *WA);

void
Output_Partial_Overlap(AS_IID s_id, AS_IID t_id, Direction_t dir,
                       const Olap_Info_t * p, int s_len, int t_len,
                       Work_Area_t  *WA);



Overlap_t  Extend_Alignment
(Match_Node_t * Match, char * S, int S_Len, char * T, int T_Len,
 int * S_Lo, int * S_Hi, int * T_Lo, int * T_Hi, int * Errors,
 Work_Area_t * WA);


int
Process_String_Olaps (char * S,
                      int Len,
                      char * S_quality,
                      AS_IID ID,
                      Direction_t Dir,
                      Work_Area_t * WA);

void
Find_Overlaps (char Frag [], int Frag_Len, char quality [], AS_IID Frag_Num, Direction_t Dir, Work_Area_t * WA);

void
Process_Overlaps (gkStream *stream, Work_Area_t * WA);

int
Build_Hash_Index(gkStream *stream, int32 first_frag_id, gkFragment *myRead);

#endif  //  OVERLAPINCORE_H


/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_OVL/AS_OVL_overlap.h
 *    src/AS_OVM/overlapInCore.H
 *
 *  Modifications by:
 *
 *    Michael Schatz on 2004-SEP-23
 *      are Copyright 2004 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2005-MAR-22
 *      are Copyright 2005 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUN-16 to 2013-AUG-01
 *      are Copyright 2005,2007-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Aaron Halpern from 2005-SEP-29 to 2006-AUG-23
 *      are Copyright 2005-2006 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Art Delcher on 2007-FEB-13
 *      are Copyright 2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2011-MAR-08
 *      are Copyright 2011 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz from 2014-DEC-15 to 2015-AUG-25
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-27
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-NOV-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include "prefixEditDistance.H"

#include <pthread.h>


#ifndef OVERLAPINCORE_H
#define OVERLAPINCORE_H

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

#define  DELETED_FRAG            2
//  Indicates fragment was marked as deleted in the store

#define  DISPLAY_WIDTH           60
//  Number of characters per line when displaying sequences

#define  ENTRIES_PER_BUCKET      21
//  In main hash table.  Recommended values are 21, 31 or 42
//  depending on cache line size.

#define  HASH_CHECK_MASK         0x1f
//  Used to set and check bit in Hash_Check_Array
//  Change if change  Check_Vector_t

#define  HASH_EXPANSION_FACTOR   1.4
//  Hash table size is >= this times  MAX_HASH_STRINGS

#define  HASH_MASK               ((1 << G.Hash_Mask_Bits) - 1)
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

#define  EXPECTED_STRING_LEN     (AS_MAX_READLEN / 2)
#define  INITIAL_DATA_LEN        (EXPECTED_STRING_LEN * Max_Hash_Strings)
//  The number of bytes to allocate initially for hash-table sequence

#define  MAX_LINE_LEN            1000
//  Maximum input line when reading FASTA file

#define  MAX_NAME_LEN            500
//  Longest file name allowed

#define  MIN_CALC_KMER           4
//  When calculating the  Hi_Hit_Limit  based on genome length, etc,
//  don't set it below this

#define  MIN_INTERSECTION        10
//  Minimum length of match region to be worth reporting

#define  MIN_OLAP_OUTSIDE_SCREEN 30
//  Minimum number of bases outside of screened regions
//  to be a reportable overlap.  Entire overlap (including screened
//  portion) must still be  MIN_OLAP_LEN .

#define  OUTPUT_OVERLAP_DELTAS   0
//  If true include delta-encoding of overlap alignment
//  in overlap messages.  Otherwise, omit them.
//  As of 6 Oct 2008, support for overlap deltas has been removed.
//  However, there are enough remnants in AS_MSG to output them.
//  Just enabling OUTPUT_OVERLAP_DELTAS will not compile; see
//  AS_MSG_USE_OVL_DELTA in AS_MSG.

#define  PROBE_MASK              0x3e
//  Used to determine probe step to resolve collisions

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

#define  MAX_EXTRA_SUBCOUNT        (AS_MAX_READLEN / G.Kmer_Len)


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

typedef  struct String_Olap_Node {
  uint32  String_Num;                // Of hash-table frag that have exact match with
  int32  Match_List;                 // Subscript of start of list of exact matches
  double  diag_sum;                  // Sum of diagonals of all k-mer matches to this frag
  int32  diag_ct;                    // Count of all k-mer matches to this frag
  signed int  Next : 29;             // Next match if this is a collision
  unsigned  Full : 1;
  unsigned  consistent : 1;
}  String_Olap_t;


typedef  struct Olap_Info {
  int  s_lo, s_hi;
  int  t_lo, t_hi;
  double  quality;
  int  delta [AS_MAX_READLEN+1];  //  needs only MAX_ERRORS
  int  delta_ct;
  int  s_left_boundary, s_right_boundary;
  int  t_left_boundary, t_right_boundary;
  int  min_diag, max_diag;
}  Olap_Info_t;

//  The following structure holds what used to be global information, but
//  is now encapsulated so that multiple copies can be made for multiple
//  parallel threads.

typedef  struct Work_Area {

  //int   *Error_Bound;

  String_Olap_t  * String_Olap_Space;
  int32  String_Olap_Size;
  int32  Next_Avail_String_Olap;

  Match_Node_t  * Match_Node_Space;
  int32  Match_Node_Size;
  int32  Next_Avail_Match_Node;

  //  Counts the number of overlaps for each fragment.  Cuts off
  //  overlaps above a limit.
  int32  A_Olaps_For_Frag;
  int32  B_Olaps_For_Frag;

  gkStore  *gkpStore;

  int    left_end_screened;
  int    right_end_screened;

  int    status;
  int    thread_id;

  uint32         bgnID;  //  Range of reads we are processing
  uint32         endID;  //  was frag_segment_lo and frag_segment_hi (all lowercase)

  //  Instead of outputting each overlap as we create it, we
  //  buffer them and output blocks of overlaps.
  uint64         overlapsLen;
  uint64         overlapsMax;
  ovOverlap     *overlaps;

  //  Various stats that used to be global and updated whenever we
  //  output an overlap or finished processing a set of hits.
  //  Needed a mutex to update.
  uint64         Total_Overlaps;
  uint64         Contained_Overlap_Ct;
  uint64         Dovetail_Overlap_Ct;

  uint64         Kmer_Hits_Without_Olap_Ct;
  uint64         Kmer_Hits_With_Olap_Ct;
  uint64         Multi_Overlap_Ct;

  prefixEditDistance  *editDist;


   char * q_diff;
   Olap_Info_t  *distinct_olap;
}  Work_Area_t;




typedef  uint32  Check_Vector_t;
// Bit vector to see if hash bucket could possibly contain a match


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


extern char           *basesData;
extern char           *qualsData;
extern String_Ref_t   *nextRef;
extern size_t          Data_Len;

extern int64   Bad_Short_Window_Ct;
extern int64   Bad_Long_Window_Ct;

extern size_t  Extra_Data_Len;

extern uint64  Max_Extra_Ref_Space;
extern uint64  Extra_Ref_Ct;
extern String_Ref_t  * Extra_Ref_Space;
extern uint64  Extra_String_Ct;
extern uint64  Extra_String_Subcount;

extern Check_Vector_t  * Hash_Check_Array;
extern uint64  Hash_String_Num_Offset;
extern Hash_Bucket_t  * Hash_Table;
extern uint64  Kmer_Hits_With_Olap_Ct;
extern uint64  Kmer_Hits_Without_Olap_Ct;
extern uint64  Multi_Overlap_Ct;
extern uint64  String_Ct;
extern Hash_Frag_Info_t  * String_Info;

extern int64  * String_Start;
extern uint32  String_Start_Size;

extern size_t  Used_Data_Len;

extern int32  Bit_Equivalent [256];
extern int32  Char_Is_Bad [256];
extern uint64  Hash_Entries;
extern uint64  Total_Overlaps;
extern uint64  Contained_Overlap_Ct;
extern uint64  Dovetail_Overlap_Ct;

class oicParameters {
public:
  oicParameters() {};
  ~oicParameters() {};

  void   initialize(void) {
    maxErate               = 0.06;
    Doing_Partial_Overlaps = false;

    bgnHashID = 1;
    endHashID = UINT32_MAX;
    minLibToHash = 0;
    maxLibToHash = UINT32_MAX;

    bgnRefID = 1;
    endRefID = UINT32_MAX;
    minLibToRef  = 0;
    maxLibToRef  = UINT32_MAX;

    Kmer_Len = 0;
    Kmer_Skip_File = NULL;

    Frag_Olap_Limit = UINT64_MAX;

    Unique_Olap_Per_Pair = true;

    Hash_Mask_Bits       = 22;
    Max_Hash_Load        = 0.6;
    Max_Hash_Strings     = 10000;
    Max_Hash_Data_Len    = 100000000;

    Outfile_Name = NULL;
    Outstat_Name = NULL;

    Num_PThreads = 1;

    Min_Olap_Len = 0;

    Use_Window_Filter  = false;

    Use_Hopeless_Check = true;

    Frag_Store_Path = NULL;
  };

  double maxErate;

  //  If set true by the G option (G for Granger)
  //  then allow overlaps that do not extend to the end
  //  of either read.
  bool  Doing_Partial_Overlaps;  //  -G

  uint32  bgnHashID;     //  -h
  uint32  endHashID;
  uint32  minLibToHash;  //  -H
  uint32  maxLibToHash;


  //  The range of reads we need to process

  uint32         frag_segment_lo;
  uint32         frag_segment_hi;

  uint32  bgnRefID;      //  -r
  uint32  curRefID;     //  When processing, this is where we are at, bgn < cur <= end.
  uint32  endRefID;
  uint32  minLibToRef;   //  -R
  uint32  maxLibToRef;

  uint32  perThread;        //  When processing, how many to do per block

  uint64  Kmer_Len;         //  -k
  FILE   *Kmer_Skip_File;   //  -k

  //  Maximum number of overlaps for end of an old fragment against
  //  a single hash table of frags, in each orientation
  uint64  Frag_Olap_Limit;  //  -l

  //  If true will allow at most
  //  one overlap output message per oriented fragment pair
  //  Set true by  -u  command-line option; set false by  -m
  bool  Unique_Olap_Per_Pair;  //  -m and -u

  uint32  Hash_Mask_Bits;  //  --hashbits

  uint32  Max_Hash_Strings;  //  --hashstrings
  uint64  Max_Hash_Data_Len;  //  --hashdatalen
  double  Max_Hash_Load;  //  --hashload

  //  --maxreadlen sets OFFSET_BITS, STRING_NUM_BITS, STRING_NUM_MASK and MAX_STRING_NUM.

  char  *Outfile_Name;  //  -o
  char  *Outstat_Name;  //  -s

  uint32  Num_PThreads;  //  -t

  int32  Min_Olap_Len;  //  --minlength, former -v

  //  Determines whether check for a window containing too many
  //  errors is used to disqualify overlaps.
  bool  Use_Window_Filter;  //  -w

  //  Determines whether check for absence of kmer matches
  //  at the end of a read is used to abort the overlap before
  //  the extension from a single kmer match is attempted.
  bool  Use_Hopeless_Check;  //  -z

  char *Frag_Store_Path;
};

extern oicParameters G;


extern uint64  HSF1;
extern uint64  HSF2;
extern uint64  SV1;
extern uint64  SV2;
extern uint64  SV3;

extern ovFile  *Out_BOF;
extern pthread_mutex_t  Frag_Store_Mutex;
extern pthread_mutex_t  Write_Proto_Mutex;




void
Output_Overlap(uint32 S_ID, int S_Len, Direction_t S_Dir,
               uint32 T_ID, int T_Len, Olap_Info_t * olap,
               Work_Area_t *WA);

void
Output_Partial_Overlap(uint32 s_id, uint32 t_id, Direction_t dir,
                       const Olap_Info_t * p, int s_len, int t_len,
                       Work_Area_t  *WA);


int
Process_String_Olaps (char * S,
                      int Len,
                      char * S_quality,
                      uint32 ID,
                      Direction_t Dir,
                      Work_Area_t * WA);

void
Find_Overlaps (char Frag [], int Frag_Len, char quality [], uint32 Frag_Num, Direction_t Dir, Work_Area_t * WA);

void *
Process_Overlaps (void *);

int
Build_Hash_Index(gkStore *store, uint32 bgnID, uint32 endID);

#endif  //  OVERLAPINCORE_H

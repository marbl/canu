
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
 *  Modifications by:
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include <pthread.h>

#include "gkStore.H"
#include "ovStore.H"

#include "correctionOutput.H"

#include <algorithm>

using namespace std;






//  Private stuff



#define Sign(a) ( ((a) > 0) - ((a) < 0) )





//  Value to add for a match in finding branch points
//  1.20 was the calculated value for 6% vs 35% error discrimination
//  Converting to integers didn't make it faster
#define  BRANCH_PT_MATCH_VALUE    0.272

//  Value to add for a mismatch in finding branch points
//   -2.19 was the calculated value for 6% vs 35% error discrimination
//  Converting to integers didn't make it faster
#define  BRANCH_PT_ERROR_VALUE    -0.728

//  Number of bits to store integer versions of error rates
#define  ERATE_BITS                  16

//  The number of errors that are ignored in setting probability
//  bound for terminating alignment extensions in edit distance
//  calculations
#define  ERRORS_FOR_FREE             1

//  Factor by which to grow memory in olap array when reading it
#define  EXPANSION_FACTOR            1.4

//  Number of old fragments to read into memory-based fragment
//  store at a time for processing
#define  FRAGS_PER_BATCH             100000

//  Longest name allowed for a file in the overlap store
#define  MAX_FILENAME_LEN            1000

//  Most errors in any edit distance computation
//  KNOWN ONLY AT RUN TIME
//#define  MAX_ERRORS                  (1 + (int) (AS_OVL_ERROR_RATE * AS_MAX_READLEN))

//  Highest number of votes before overflow
#define  MAX_DEGREE                  32767

//  Highest number of votes before overflow
#define  MAX_VOTE                    255

//  Branch points must be at least this many bases from the
//  end of the fragment to be reported
#define  MIN_BRANCH_END_DIST         20

//  Branch point tails must fall off from the max by at least
//  this rate
#define  MIN_BRANCH_TAIL_SLOPE       0.20

//  This many or more votes at the same base indicate
//  a separate haplotype
#define  MIN_HAPLO_OCCURS            3

//  The amount of memory to allocate for the stack of each thread
#define  THREAD_STACKSIZE        (128 * 512 * 512)





struct Vote_Tally_t {
  uint32  confirmed : 8;
  uint32  deletes   : 8;
  uint32  a_subst   : 8;
  uint32  c_subst   : 8;

  uint32  g_subst   : 8;
  uint32  t_subst   : 8;
  uint32  no_insert : 8;
  uint32  a_insert  : 8;

  uint32  c_insert  : 8;
  uint32  g_insert  : 8;
  uint32  t_insert  : 8;
};


struct Vote_t {
  int32         frag_sub;
  int32         align_sub;
  Vote_Value_t  vote_val;
};



class Frag_Info_t {
public:
  Frag_Info_t() {
    sequence     = NULL;
    vote         = NULL;
    clear_len    = 0;
    left_degree  = 0;
    right_degree = 0;
    shredded     = false;
    unused       = false;
  };
  ~Frag_Info_t() {
  };

  char          *sequence;
  Vote_Tally_t  *vote;
  uint64         clear_len     : 31;
  uint64         left_degree   : 31;
  uint64         right_degree  : 31;
  uint64         shredded      : 1;    // True if shredded read
  uint64         unused        : 1;
};

class Olap_Info_t {
public:
  Olap_Info_t() {};
  ~Olap_Info_t() {};

  uint32      a_iid;
  uint32      b_iid;
  int64       a_hang : 31;
  int64       b_hang : 31;
  uint64      innie  : 1;  //  was 'orient' with choice INNIE=0 or NORMAL=1
  uint64      normal : 1;  //  so 'normal' always != 'innie'

  //  Sort by increasing b_iid, then increasing a_iid.
  bool  operator<(Olap_Info_t const &that) const {
    if (b_iid < that.b_iid)      return(true);
    if (b_iid > that.b_iid)      return(false);

    if (a_iid < that.a_iid)      return(true);
    if (a_iid > that.a_iid)      return(false);

    //  It is possible, but unlikely, to have two overlaps to the same pair of reads,
    //  if we overlap a5'-b3' and a3'-b5'.  I think.

    return(innie != that.innie);
  };
};



class Frag_List_t {
public:
  Frag_List_t() {
    readsMax    = 0;
    readsLen    = 0;
    readIDs     = NULL;
    readBases   = NULL;
    basesMax    = 0;
    basesLen    = 0;
    bases       = NULL;
  };

  ~Frag_List_t() {
    delete [] readIDs;
    delete [] readBases;
    delete [] bases;
  };

  uint32             readsMax;
  uint32             readsLen;
  uint32            *readIDs;
  char             **readBases;

  uint64             basesMax;
  uint64             basesLen;
  char              *bases;        //  Read sequences, 0 terminated
};



class feParameters;


class pedWorkArea_t {
public:
  pedWorkArea_t() {
    G        = NULL;

    memset(delta,      0, sizeof(int32) * AS_MAX_READLEN);
    memset(deltaStack, 0, sizeof(int32) * AS_MAX_READLEN);

    deltaLen = 0;

    Edit_Array_Lazy = NULL;
  };

  ~pedWorkArea_t() {

    for (uint32 xx=0; xx < alloc.size(); xx++)
      delete [] alloc[xx];

    delete [] Edit_Array_Lazy;
  };

  void          initialize(feParameters *G_, double errorRate) {
    G = G_;

    Edit_Array_Max  = 1 + (uint32)(errorRate * AS_MAX_READLEN);
    Edit_Array_Lazy = new int32 * [Edit_Array_Max];

    memset(Edit_Array_Lazy, 0, sizeof(int32 *) * Edit_Array_Max);
  };

public:
  feParameters *G;

  int32              delta[AS_MAX_READLEN];  //  Only need ERATE * READLEN
  int32              deltaStack[AS_MAX_READLEN];
  int32              deltaLen;

  vector<int32 *>    alloc;            //  Allocated blocks, don't use directly.

  int32            **Edit_Array_Lazy;  //  Doled out space.
  int32              Edit_Array_Max;   //  Former MAX_ERRORS
};



struct Thread_Work_Area_t {
  int32         thread_id;
  uint32        loID;
  uint32        hiID;
  uint64        nextOlap;

  feParameters *G;

  Frag_List_t  *frag_list;

  char          rev_seq[AS_MAX_READLEN + 1];  //  Used in Process_Olap to hold RC of the B read
  uint32        rev_id;                       //  Ident of the rev_seq read.

  Vote_t        globalvote[AS_MAX_READLEN];

  uint32        failedOlaps;

  pedWorkArea_t ped;
};










class feParameters {
public:
  feParameters() {
    gkpStorePath   = NULL;
    ovlStorePath   = NULL;

    bgnID          = 0;
    endID          = UINT32_MAX;

    readBases      = NULL;
    readVotes      = NULL;
    reads          = NULL;
    readsLen       = 0;

    olaps          = NULL;
    olapsLen       = 0;

    outputFileName = NULL;

    numThreads     = 4;
    errorRate      = 0.06;
    minOverlap     = 0;

    //  Parameters

    //  Output
    Degree_Threshold = 2;  //, DEFAULT_DEGREE_THRESHOLD;
    Use_Haplo_Ct     = TRUE;

    //  Analyze_Alignment
    End_Exclude_Len   = 3;  //DEFAULT_END_EXCLUDE_LEN;
    Kmer_Len          = 9;  //DEFAULT_KMER_LEN;
    Vote_Qualify_Len  = 9; //DEFAULT_VOTE_QUALIFY_LEN;
  };
  ~feParameters() {
    delete [] readBases;
    delete [] readVotes;
    delete [] reads;
    delete [] olaps;
  };


  //  Paths to stores
  char         *gkpStorePath;
  char         *ovlStorePath;

  // Range of IDs to process
  uint32        bgnID;
  uint32        endID;

  char         *readBases;
  Vote_Tally_t *readVotes;
  Frag_Info_t  *reads;
  uint32        readsLen;  // Number of fragments being corrected

  Olap_Info_t  *olaps;
  uint64        olapsLen;  // Number of overlaps being used

  char         *outputFileName;

  uint32        numThreads;

  double        errorRate;
  uint32        minOverlap;

  // This array [e] is the minimum value of Edit_Array [e] [d] to be worth pursuing in edit-distance
  // computations between guides (only MAX_ERRORS needed)
  int  Edit_Match_Limit [AS_MAX_READLEN + 1];


  // Set keep flag on end of fragment if number of olaps < this value
  int  Degree_Threshold;
  // Set false by -h option to ignore haplotype counts when correcting
  int  Use_Haplo_Ct;


  // Length of ends of exact-match regions not used in preventing sequence correction
  int  End_Exclude_Len;
  // Length of minimum exact match in overlap to confirm base pairs
  int  Kmer_Len;
  // Number of bases surrounding a SNP to vote for change
  int  Vote_Qualify_Len;


  //  This array [i] is the maximum number of errors allowed in a match between sequences of length
  //  i , which is i * MAXERROR_RATE .
  int  Error_Bound [AS_MAX_READLEN + 1];

  // To make debugging printout come out together
  pthread_mutex_t  Print_Mutex;
};


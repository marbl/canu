
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "system.H"

#include <pthread.h>

#include "sqStore.H"
#include "ovStore.H"

#include "correctionOutput.H"

#include <algorithm>
#include <string>
#include <vector>







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

//  Longest name allowed for a file in the overlap store
#define  MAX_FILENAME_LEN            1000

//  Most errors in any edit distance computation
//  KNOWN ONLY AT RUN TIME
//#define  MAX_ERRORS                  (1 + (int) (AS_OVL_ERROR_RATE * AS_MAX_READLEN))

//  Highest number of votes before overflow
#define  MAX_DEGREE                  32767

//  Highest number of votes before overflow
#define  MAX_VOTE                    65535

//  Branch points must be at least this many bases from the
//  end of the fragment to be reported
#define  MIN_BRANCH_END_DIST         20

//  Branch point tails must fall off from the max by at least
//  this rate
#define  MIN_BRANCH_TAIL_SLOPE       0.20

//  Votes at this distance from the extremities of aligned region are ignored
//  Leads to start/end of reads not being corrected,
//  which should be accounted by overlap scoring (see DEFAULT_FLANK_IGNORE)
#define  IGNORE_FLANK_VOTE           5

//  The amount of memory to allocate for the stack of each thread
#define  THREAD_STACKSIZE        (128 * 512 * 512)

struct Vote_Tally_t {
  Vote_Tally_t() {
     confirmed = 0;
     conf_no_insert = 0;
     deletes   = 0;
     a_subst   = 0;
     c_subst   = 0;
     g_subst   = 0;
     t_subst   = 0;
     no_insert = 0;

     insertion_cnt = 0;
  };

  static const char INSERTIONS_DELIM = '$';

  uint32  confirmed      : 16;
  uint32  conf_no_insert : 16;
  uint32  deletes        : 16;
  uint32  a_subst        : 16;

  uint32  c_subst        : 16;
  uint32  g_subst        : 16;
  uint32  t_subst        : 16;
  uint32  no_insert      : 16;

  uint32 insertion_cnt;
  std::string insertions;

  std::vector<std::string> insertions_list() const {
    std::vector<std::string> answer;
    std::size_t start;
    std::size_t end = 0;

    while ((start = insertions.find_first_not_of(INSERTIONS_DELIM, end)) != std::string::npos) {
      end = insertions.find(INSERTIONS_DELIM, start);
      answer.push_back(insertions.substr(start, end - start));
    }
    return answer;
  }

  //NB: total does not consider insertions
  uint32 total() const {
    return deletes + a_subst + c_subst + g_subst + t_subst;
  }

  uint32 subst(char bp) const {
    switch (bp) {
      case 'a':
        return a_subst;
      case 'c':
        return c_subst;
      case 'g':
        return g_subst;
      case 't':
        return t_subst;
      default:
        assert(false);
    }
    return 0;
  }

  uint32 all_but(char bp) const {
    return total() + insertion_cnt - subst(bp);
  }
};


struct Vote_t {
  int32         frag_sub;
  int32         align_sub;
  Vote_Value_t  vote_val;
};



struct Frag_Info_t {
  Frag_Info_t() {
    sequence     = NULL;
    vote         = NULL;
    clear_len    = 0;
    left_degree  = 0;
    right_degree = 0;
    shredded     = false;
    unused       = false;
  };

  char          *sequence;
  Vote_Tally_t  *vote;
  uint64         clear_len     : 31;
  uint64         left_degree   : 31;
  uint64         right_degree  : 31;
  uint64         shredded      : 1;    // True if shredded read
  uint64         unused        : 1;
};

struct Olap_Info_t {
  Olap_Info_t() {
    a_iid   = 0;
    b_iid   = 0;
    a_hang  = 0;
    b_hang  = 0;
    innie   = false;
    normal  = false;
  };

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

    return(innie < that.innie);
  };
};



struct Frag_List_t {
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



struct feParameters;


struct pedWorkArea_t {
  pedWorkArea_t() {
    G        = NULL;

    memset(delta,      0, sizeof(int32) * AS_MAX_READLEN);
    memset(deltaStack, 0, sizeof(int32) * AS_MAX_READLEN);

    deltaLen = 0;

    Edit_Array_Lazy = NULL;
    Edit_Array_Max  = 0;
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

  feParameters *G;

  int32                delta[AS_MAX_READLEN];  //  Only need ERATE * READLEN
  int32                deltaStack[AS_MAX_READLEN];
  int32                deltaLen;

  std::vector<int32 *> alloc;            //  Allocated blocks, don't use directly.

  int32              **Edit_Array_Lazy;  //  Doled out space.
  int32                Edit_Array_Max;   //  Former MAX_ERRORS
};



struct Thread_Work_Area_t {
  int32         thread_id;
  uint64        nextOlap;

  feParameters *G;

  Frag_List_t  *frag_list;

  char          rev_seq[AS_MAX_READLEN + 1];  //  Used in Process_Olap to hold RC of the B read
  uint32        rev_id;                       //  Ident of the rev_seq read.

  Vote_t        globalvote[AS_MAX_READLEN];

  uint64        passedOlaps;
  uint64        failedOlaps;

  pedWorkArea_t ped;
};










struct feParameters {
  feParameters() {
    seqStorePath   = NULL;
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

    numThreads     = getMaxThreadsAllowed();
    errorRate      = 0.06;
    minOverlap     = 0;

    //  Parameters

    //  Output
    Degree_Threshold = 2;  //, DEFAULT_DEGREE_THRESHOLD;

    Haplo_Expected   = 2;
    Haplo_Strong     = 2;
    Haplo_Confirm    = 5;
    Haplo_Freeze     = 1;

    checkTrivialDNA = false;
    maskedErrorRate = 0.06;

    //  Analyze_Alignment

    End_Exclude_Len   = 3;  //DEFAULT_END_EXCLUDE_LEN;
    //Kmer_Len          = 9;  //DEFAULT_KMER_LEN;
    Vote_Qualify_Len  = 9; //DEFAULT_VOTE_QUALIFY_LEN;
  };

  ~feParameters() {
    delete [] readBases;
    delete [] readVotes;
    delete [] reads;
    delete [] olaps;
  };


  //  Paths to stores
  char         *seqStorePath;
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

  //  Require at least Haplo_Confirm reads to confirm a haplotype, but if
  //  there are more than Haplo_Expected haplotypes detected, treat the
  //  column as a true error.
  //
  int32   Haplo_Expected;  //  (2) Number of haplotypes expected; previously hard-coded inline '2'.
  int32   Haplo_Strong;    //  (2) Number of reads to confirm insert/delete; previously STRONG_CONFIRMATION_READ_CNT.
  int32   Haplo_Confirm;   //  (5) Number of reads to confirm a haplotype; previously MIN_HAPLO_OCCURS.
  uint32  Haplo_Freeze;    //  (1) Radius around heterozygous position that is considered too risky to correct.

  bool    checkTrivialDNA;
  double  maskedErrorRate;

  // Length of ends of exact-match regions not used in preventing sequence correction
  int  End_Exclude_Len;
  // Length of minimum exact match in overlap to confirm base pairs
  //int  Kmer_Len;
  // Number of bases surrounding a SNP to vote for change
  int  Vote_Qualify_Len;


  //  This array [i] is the maximum number of errors allowed in a match between sequences of length
  //  i , which is i * MAXERROR_RATE .
  int  Error_Bound [AS_MAX_READLEN + 1];
};


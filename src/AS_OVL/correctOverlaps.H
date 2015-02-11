
#include "AS_global.H"

#include <pthread.h>

#include "gkStore.H"
#include "ovStore.H"

#include <algorithm>

using namespace std;






#define Sign(a) ( ((a) > 0) - ((a) < 0) ) 





//  Value to add for a match in finding branch points
//  1.20 was the calculated value for 6% vs 35% error discrimination
//  Converting to integers didn't make it faster
#define  BRANCH_PT_MATCH_VALUE    0.272

//  Value to add for a mismatch in finding branch points
//   -2.19 was the calculated value for 6% vs 35% error discrimination
//  Converting to integers didn't make it faster
#define  BRANCH_PT_ERROR_VALUE    -0.728

//  Default value for  End_Exclude_Len
#define  DEFAULT_END_EXCLUDE_LEN     3

//  Default value for bases on each side of SNP to vote for change
#define  DEFAULT_HALF_LEN            4

//  Default value for  Kmer_Len
#define  DEFAULT_KMER_LEN            9

//  Default value for  Quality_Threshold
#define  DEFAULT_QUALITY_THRESHOLD   0.015

//  Probability limit to "band" edit-distance calculation
//  Determines  NORMAL_DISTRIB_THOLD
#define  EDIT_DIST_PROB_BOUND        1e-4

//  The number of errors that are ignored in setting probability
//  bound for terminating alignment extensions in edit distance
//  calculations
#define  ERRORS_FOR_FREE             1

//  Most bytes allowed in line of fasta file
#define  MAX_FASTA_LINE              2048

//  Longest name allowed for a file in the overlap store
#define  MAX_FILENAME_LEN            1000

//  Most errors in any edit distance computation // 0.40
//  KNOWN ONLY AT RUN TIME
#define  MAX_ERRORS                  (1 + (int) (AS_OVL_ERROR_RATE * AS_MAX_READLEN))

// Factor by which to grow memory in olap array when reading it
#define  EXPANSION_FACTOR            1.4

//  Branch points must be at least this many bases from the
//  end of the fragment to be reported
#define  MIN_BRANCH_END_DIST     20

//  Branch point tails must fall off from the max by at least
//  this rate
#define  MIN_BRANCH_TAIL_SLOPE   0.20

//  Determined by  EDIT_DIST_PROB_BOUND
#define  NORMAL_DISTRIB_THOLD    3.62





struct Adjust_t {
  int32  pos;
  int32  adjust;
};


class Frag_Info_t {
public:
  Frag_Info_t() {
    bases        = NULL;
    basesLen     = 0;

    adjusts      = NULL;
    adjustsLen   = 0;

    keep_left    = false;
    keep_right   = false;
  };
  ~Frag_Info_t() {
  };

  char          *bases;
  Adjust_t      *adjusts;

  uint32         basesLen;
  uint32         adjustsLen;

  uint32         keep_right     : 1;   //  I think these are unused.
  uint32         keep_left      : 1;   //  If so, we get back 8 bytes.  If not, redo Correct_Frags to use a temporary, and use 31 bits for the two lengths.
};


class Olap_Info_t {
public:
  Olap_Info_t() {};
  ~Olap_Info_t() {};

  uint32      a_iid;
  uint32      b_iid;
  int64       a_hang : 31;
  int64       b_hang : 31;
  int64       innie  : 1;  //  was 'orient' with choice INNIE=0 or NORMAL=1
  int64       normal : 1;  //  so 'normal' always != 'innie'

  uint64      order;
  uint32      evalue;
};

//  Sort by increasing b_iid.
//
//  It is possible, but unlikely, to have two overlaps to the same pair of reads,
//  if we overlap a5'-b3' and a3'-b5'.  I think.
//
class Olap_Info_t_by_bID {
public:
  inline bool  operator()(const Olap_Info_t &a, const Olap_Info_t &b) {
    if (a.b_iid < b.b_iid)      return(true);
    if (a.b_iid > b.b_iid)      return(false);

    if (a.a_iid < b.a_iid)      return(true);
    if (a.a_iid > b.a_iid)      return(false);

    return(a.innie != b.innie);
  };
};


class Olap_Info_t_by_Order {
public:
  inline bool  operator()(const Olap_Info_t &a, const Olap_Info_t &b) {
    return(a.order < b.order);
  };
};








class coParameters {
public:
  coParameters() {
  };
  ~coParameters() {
    delete [] bases;
    delete [] adjusts;
    delete [] reads;
    delete [] olaps;
  };


  //  Paths to stores
  char         *gkpStorePath;
  char         *ovlStorePath;

  //  Input read corrections, output overlap corrections
  char         *correctionsName;
  char         *eratesName;

  // Range of IDs to process
  uint32        bgnID;
  uint32        endID;

  char         *bases;
  uint64        basesLen;

  Adjust_t     *adjusts;
  uint64        adjustsLen;

  Frag_Info_t  *reads;
  uint32        readsLen;  // Number of fragments being corrected

  Olap_Info_t  *olaps;
  uint64        olapsLen;  // Number of overlaps being used

  uint32        numThreads;

 
//  Globals

int  Edit_Match_Limit[AS_MAX_READLEN+1];
    // This array[e] is the minimum value of  Edit_Array[e][d]
    // to be worth pursuing in edit-distance computations between guides
    // (only MAX_ERRORS needed)

int  Error_Bound[AS_MAX_READLEN + 1];
    // This array[i]  is the maximum number of errors allowed
    // in a match between sequences of length  i , which is
    //  i * MAXERROR_RATE .




#if 0

char  * Correct_File_Path;
    // Name of file containing fragment corrections

FILE  * Delete_fp = NULL;
    // File to which list of overlaps to delete is written if  -x  option is specified

  int  * Edit_Array_Lazy[AS_MAX_READLEN+1];
    // Use for alignment calculation.  Points into  Edit_Space .
    // (only MAX_ERRORS needed)


int * Edit_Space_Lazy[AS_MAX_READLEN+1];
    // Memory used by alignment calculation
    // (only (MAX_ERRORS + 4) * MAX_ERRORS needed)

int  End_Exclude_Len = DEFAULT_END_EXCLUDE_LEN;
    // Length of ends of exact-match regions not used in preventing
    // sequence correction

char  * Erate_Path = NULL;
    // Name of binary file to which to dump revised error rates
    // This allows the program to run in parallel (e.g., under LSF)
    // Presumably the error-rate files will be added to the overlap
    // store later by the  update-erates  program.

int  Failed_Alignments_Ct = 0;
    // Count the number of alignments that failed

Frag_Info_t  * Frag;
    // Sequence and vote information for current range of fragments
    // being corrected

int32  Frags_Per_File;
    // Maximum number of fragments in each data file of fragment store.
    // This is read from the store

int  Half_Len = DEFAULT_HALF_LEN;
    // Number of bases on each side of SNP to vote for change

AS_IID  Hi_Frag_IID;
    // Internal ID of last fragment in frag store to process

int  Hi_Unitig = -1;
    // The highest numbered unitig.

AS_IID  Highest_Frag = 0;
    // The highest numbered fragment in the unitig messages.

AS_IID* IUM = NULL;
    // Has unitig ID for each fragment.

uint32 IUM_Size = 0;
    // Number of entries in  IUM .

Int_List_t  * Keep_Pair = NULL;
    // Array holding pairs of unitigs that corrected overlaps imply
    // should overlap.

int  Kmer_Len = DEFAULT_KMER_LEN;
    // Length of minimum exact match in overlap to confirm base pairs

AS_IID Lo_Frag_IID;
    // Internal ID of first fragment in frag store to process

uint32 Num_Frags = 0;
    // Number of fragments being corrected

uint64  Num_Olaps;
    // Number of overlaps being used

Olap_Info_t  * Olap = NULL;
    // Array of overlaps being used

uint32  * Olap_Offset = NULL;
    // Indicates the first overlap of each fragment

char  * Olap_Path;
    // Name of file containing a sorted list of overlaps

bool Olaps_From_Store = FALSE;
    // Indicates if overlap info comes from  get-olaps  or from
    // a binary overlap store

FILE  * OVL_fp = NULL;
    // File to which OVL messages are written if  -o  option is specified


int  Total_Alignments_Ct = 0;
    // Count the number of alignments attempted

#endif

};

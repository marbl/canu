
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
 *    Brian P. Walenz beginning on 2016-MAY-02
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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
//#define  MAX_ERRORS                  (1 + (int) (AS_OVL_ERROR_RATE * AS_MAX_READLEN))

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
  int32  adjpos;
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
  uint64      innie  : 1;  //  was 'orient' with choice INNIE=0 or NORMAL=1
  uint64      normal : 1;  //  so 'normal' always != 'innie'

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




class coParameters;


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

  void          initialize(coParameters *G_, double errorRate) {
    G = G_;

    fprintf(stderr, "-- Allocate "F_U64" MB for Edit_Array pointers.\n", (sizeof(int32 *) * Edit_Array_Max) >> 20);

    Edit_Array_Max  = 1 + (uint32)(errorRate * AS_MAX_READLEN);
    Edit_Array_Lazy = new int32 * [Edit_Array_Max];

    memset(Edit_Array_Lazy, 0, sizeof(int32 *) * Edit_Array_Max);
  };

public:
  coParameters *G;

  int32             delta[AS_MAX_READLEN];  //  Only need ERATE * READLEN
  int32             deltaStack[AS_MAX_READLEN];
  int32             deltaLen;

  vector<int32 *>   alloc;            //  Allocated blocks, don't use directly.

  int32           **Edit_Array_Lazy;  //  Doled out space.
  int32             Edit_Array_Max;   //  Former MAX_ERRORS
};




class coParameters {
public:
  coParameters() {
    gkpStorePath = NULL;
    ovlStorePath = NULL;

    //  Input read corrections, output overlap corrections
    correctionsName = NULL;
    eratesName      = NULL;

    // Range of IDs to process
    bgnID = 0;
    endID = UINT32_MAX;

    bases    = NULL;
    basesLen = 0;

    adjusts    = NULL;
    adjustsLen = 0;

    reads    = NULL;
    readsLen = 0;

    olaps    = NULL;
    olapsLen = 0;

    numThreads = 1;
    errorRate  = 0.06;
    minOverlap = 0;
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

  //  Range of IDs to process
  uint32        bgnID;
  uint32        endID;

  char         *bases;
  uint64        basesLen;

  Adjust_t     *adjusts;
  uint64        adjustsLen;

  Frag_Info_t  *reads;     //  These are relative to bgnID!
  uint32        readsLen;  //  Number of fragments being corrected

  Olap_Info_t  *olaps;
  uint64        olapsLen;  //  Number of overlaps being used

  uint32        numThreads;  //  Only one thread supported.

  double        errorRate;
  uint32        minOverlap;

  pedWorkArea_t ped;

  //  Globals

  // This array[e] is the minimum value of  Edit_Array[e][d]
  // to be worth pursuing in edit-distance computations between guides
  // (only MAX_ERRORS needed)
  int  Edit_Match_Limit[AS_MAX_READLEN+1];

  // This array[i]  is the maximum number of errors allowed
  // in a match between sequences of length  i , which is
  //  i * MAXERROR_RATE .
  int  Error_Bound[AS_MAX_READLEN + 1];
};

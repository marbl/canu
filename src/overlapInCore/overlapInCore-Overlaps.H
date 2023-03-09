


#if 0
typedef  enum Direction_Type {
  FORWARD,
  REVERSE
} Direction_t;
#endif


typedef  struct String_Olap_Node {
  uint32  String_Num;                // Of hash-table frag that have exact match with
  int32  Match_List;                 // Subscript of start of list of exact matches
  double  diag_sum;                  // Sum of diagonals of all k-mer matches to this frag
  int32   diag_ct;                   // Count of all k-mer matches to this frag
  int     diag_bgn;
  int     diag_end;
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



//  Various stats that used to be global and updated whenever we output an
//  overlap or finished processing a set of hits.  Needed a mutex to update.
//
class overlapHitStats {
public:
  overlapHitStats &operator+=(overlapHitStats const A) {
    totalOverlaps      += A.totalOverlaps;
    containedOverlaps  += A.containedOverlaps;
    dovetailOverlaps   += A.dovetailOverlaps;

    hitsWithoutOverlap += A.hitsWithoutOverlap;
    hitsWithOverlap    += A.hitsWithOverlap;
    hitsSkipped        += A.hitsSkipped;

    multipleOverlap    += A.multipleOverlap;

    return(*this);
  };

  overlapHitStats &operator=(uint32 A) {
    totalOverlaps       = 0;
    containedOverlaps   = 0;
    dovetailOverlaps    = 0;

    hitsWithoutOverlap  = 0;
    hitsWithOverlap     = 0;
    hitsSkipped         = 0;

    multipleOverlap     = 0;

    return(*this);
  };


  uint64   totalOverlaps;
  uint64   containedOverlaps;
  uint64   dovetailOverlaps;

  uint64   hitsWithoutOverlap;
  uint64   hitsWithOverlap;
  uint64   hitsSkipped;

  uint64   multipleOverlap;
};








class workArea {
public:
  workArea(sqStore *readStore,
           sqCache *readCache,
           bool     partialOverlaps,
           double   maxErate,
           double   maxAlignErate) {
    String_Olap_Size  = 5000;
    String_Olap_Space = new String_Olap_t [String_Olap_Size];

    Match_Node_Size  = 10000;
    Match_Node_Space = new Match_Node_t [Match_Node_Size];

    readStore     = readStore;
    readCache     = readCache;

    overlapsLen   = 0;
    overlapsMax   = 1024 * 1024 / sizeof(ovOverlap);
    overlaps      = new ovOverlap [overlapsMax];

    editDist      = new prefixEditDistance(partialOverlaps, maxErate, maxAlignErate);

    q_diff        = new char        [AS_MAX_READLEN];
    distinct_olap = new Olap_Info_t [MAX_DISTINCT_OLAPS];
  };

  ~workArea() {
    delete    editDist;
    delete [] String_Olap_Space;
    delete [] Match_Node_Space;
    delete [] overlaps;

    delete [] distinct_olap;
    delete [] q_diff;
  };


public:
  void
  findOverlaps(char        *readSeq,
               uint32       readLen,
               uint32       readID,
               bool         doForward,
               hashTable   *HT);

  void
  processOverlaps(char   *readSeq,
                  uint32  readLen,
                  uint32  readID,
                  bool    doForward);

  const
  overlapHitStats &
  outputOverlaps(uint32 batch,
                 uint32 readBgn, uint32 readEnd,
                 ovFile *OF) {

    fprintf(stderr, "Batch %03u writes    reads " F_U32 "-" F_U32 " (" F_U64 " overlaps " F_U64 "/" F_U64 "/" F_U64 " kmer hits with/without overlap/skipped)\n",
            batch, readBgn, readEnd, overlapsLen,
            statistics.hitsWithOverlap,
            statistics.hitsWithoutOverlap,
            statistics.hitsSkipped);

    for (int zz=0; zz<overlapsLen; zz++)
      OF->writeOverlap(overlaps + zz);

    overlapsLen = 0;

    //globalStats += statistics;

    return(statistics);
  };

  //  Add a reference ... ???
private:
  void
   addMatch(kmerRef kmer);

  void
  addReference(kmerRef  kmer,
               uint64   offset);

  //map<uint32, String_Olap_t>  matches;

  String_Olap_t  *String_Olap_Space;        //  One entry per stream read with hits
  int32           String_Olap_Size;
  int32           Next_Avail_String_Olap;


  //  Defined in prefixEditDistance.
  //struct Match_Node_t {
  //  int32  Offset;              // To start of exact match in  hash-table frag
  //  int32  Len;                 // Of exact match
  //  int32  Start;               // Of exact match in current (new) frag
  //  int32  Next;                // Subscript of next match in list
  //};

  Match_Node_t   *Match_Node_Space;
  int32           Match_Node_Size;
  int32           Next_Avail_Match_Node;

  //  Counts the number of overlaps for each fragment.  Cuts off
  //  overlaps above a limit.
  int32  A_Olaps_For_Frag;
  int32  B_Olaps_For_Frag;

  sqStore  *readStore;
  sqCache  *readCache;

  int    left_end_screened;
  int    right_end_screened;

  //uint32         bgnID;  //  Range of reads we are processing
  //uint32         endID;  //  was frag_segment_lo and frag_segment_hi (all lowercase)

  //  Instead of outputting each overlap as we create it, we
  //  buffer them and output blocks of overlaps.
  uint64         overlapsLen;
  uint64         overlapsMax;
  ovOverlap     *overlaps;

  prefixEditDistance  *editDist;

  char * q_diff;
  Olap_Info_t  *distinct_olap;

  overlapHitStats   statistics;
};

#include <pthread.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <fcntl.h>
#include <assert.h>
#include <math.h>

#include <new>

#include "bio++.H"

#include "positionDB.H"
#include "existDB.H"
#include "sim4.H"

#include "logMsg.H"


#define MAX_THREADS  64


//  A meta-option.  Enable all the good stuff.  Or not, I guess.
//
#if 0

//  Define this to print a message for each search, showing times and
//  stats, unless it finished quicker than some minimum time.
//
#define VERBOSE_SEARCH
#define VERBOSE_SEARCH_MINIMUM_TIME  1.0

//  Define this to print the number of hits (raw and filtered) for each
//  sequence, unless it has fewer than some minimum number of raw hits.
//
#define VERBOSE_FILTER
#define VERBOSE_FILTER_MINIMUM  10

//  Define this to show the hitMatrix
#define SHOW_HITMATRIX

//  Define this to print a message whenever a polish starts.
//
#define SHOW_POLISHING

//  Define these to show polishes that take a long time -- individual
//  polishes, not all polishes for a single sequence.  The time is in
//  seconds.
//
#define SHOW_POLISHING_EXPENSIVE  0.5

//  Define this to show the exon discarding and match splitting.
//
#define SHOW_MATCH_SPLITTING

//  Define this to show the "hit discarding" results.  Any hits
//  that look like they are repeats are re-searched using nearly
//  unique mers.
//
#define SHOW_HIT_DISCARDING

//  Define this to show the kmers being added to the sim4command from
//  thr-polish.C.  This generates a lot of output!
//
#define SHOW_HITS_ADDED
//#define SHOW_HITS_ADDED_AFTER_QUERY   0

//  END OF GOOD STUFF!
#endif



class configuration {
public:
  configuration();
  ~configuration();

  void  read(int argc, char **argv);

  void  setTime(struct timespec *ts, double t) {
    ts->tv_sec  = (time_t)floor(t);
    ts->tv_nsec = (long)((t - ts->tv_sec) * 1e9);
  };

public:
  bool             _beVerbose;

  //  These are private to the kMerBuilder -- in particular, we cannot
  //  use _merSize as a surrogate for the length of the seed,
  //  compressed seeds can be much longer.
  //
  uint32           _KBmerSize;
  uint32           _KBcompression;
  char            *_KBspacingTemplate;

  uint32           _merSkip;

  uint32           _numSearchThreads;

  bool             _doReverse;
  bool             _doForward;
  bool             _doValidation;
  char            *_doValidationFileName;

  bool             _doAlignments;

  double           _Lo;
  double           _Hi;
  double           _Va;

  uint32           _maxDiagonal;

  double           _extendWeight;
  uint32           _extendMinimum;
  uint32           _extendMaximum;

  uint32           _repeatThreshold;

  //  Minimums for hits
  double           _minHitCoverage;
  uint32           _minHitLength;

  //  Minimums for matches
  uint32           _minMatchIdentity;
  uint32           _minMatchCoverage;

  //  Filtering of hits
  bool             _afEnabled;
  double           _afThreshold;
  uint32           _afLength;
  uint32           _afInit;

  //  Filtering and/or cleanup of matches
  uint32           _discardExonLength;
  uint32           _discardExonQuality;
  bool             _splitMatches;
  bool             _polishOptimally;

  char            *_dbFileName;
  char            *_psFileName;
  char            *_qsFileName;

  char            *_maskFileName;
  char            *_onlyFileName;

  bool             _buildOnly;

  uint32           _ignoreThreshold;
  char            *_maskPrefix;
  uint32           _maskThreshold;
  char            *_onlyPrefix;
  uint32           _onlyThreshold;

  char            *_outputFileName;
  char            *_logmsgFileName;
  char            *_statsFileName;
};






//  Shared data
//
extern configuration          config;
extern sim4parameters         sim4params;

extern seqCache              *genome;
extern seqStream             *genomeMap;

extern seqCache              *qsFASTA;  //  Used exclusively by thr-loader.C

extern existDB               *maskDB;  //  thr-search.C
extern existDB               *onlyDB;  //  thr-search.C

extern positionDB            *positions;

extern volatile uint32        numberOfQueries;



//  aHit -- storing the internal hits
//
//  _status
//  & 0x00000001 -- direction, forward if set, otherwise reverse
//    0x00000002 -- not filtered, if set, polish the hit
//    0x00000004 -- not filtered, if set, polish the hit because it has something unique-ish
//    0x00000008 -- match verified via polishng
//    0x00000010 -- match discarded via hit refinement
//    0x0000fff8 -- unused
//    0x00ff0000 -- percent idendity of match
//    0xff000000 -- percent coverage of match
//
#define AHIT_DIRECTION_MASK  0x00000001
#define AHIT_POLISHABLE      0x00000002
#define AHIT_HAS_UNIQUE      0x00000004
#define AHIT_VERIFIED        0x00000008
#define AHIT_DISCARDED       0x00000010


struct aHit {
  uint32     _status;
  uint32     _qsIdx;
  uint32     _dsIdx;
  uint32     _dsLo;
  uint32     _dsHi;
  uint32     _covered;
  uint32     _matched;
  uint32     _numMers;
  merList   *_ML;
};



class query {
public:
  query() {
    seq          = 0L;

    theHitsLen   = 0;
    theHitsMax   = 4;
    theHits      = new aHit [theHitsMax];

    theLog       = 0L;
    if (config._logmsgFileName)
      theLog     = new logMsg(true);

    theOutputLen = 0;
    theOutputMax = 0;
    theOutput    = 0L;
  };

  ~query() {
    delete    seq;
    delete [] theHits;
    delete [] theLog;
    delete [] theOutput;
  };

  bool         loadSequence(seqCache *qs) {
    seq = qs->getSequenceInCore();
    return(seq != 0L);
  };

  seqInCore   *seq;
  uint32       theHitsLen;
  uint32       theHitsMax;
  aHit        *theHits;

  logMsg      *theLog;

  uint32       theOutputLen;
  uint32       theOutputMax;
  char        *theOutput;
};



class searcherState {
public:
  uint64         threadID;

  uint64         posnMax;
  uint64         posnLen;
  uint64        *posn;

  kMerBuilder   *KB;
  void          *DP;

  searcherState(uint64 U) {
    threadID    = U;

    posnMax     = 0;
    posnLen     = 0;
    posn        = 0L;

    KB          = 0L;
    DP          = 0L;
  };

  ~searcherState() {
    delete [] posn;
    delete    KB;
  };
};




struct diagonalLine {
  union {
    uint64  all;
    struct {
      uint64     uniq : 10;  //  uniqueness score for the mer here
      uint64     qPos : 22;  //  position in the query, 0 to 4M
      uint64     dPos : 32;  //  position in the genome, 0 to 4G
    } val;
  };
};



class trapMatch {
public:
  uint32         _unique;
  uint32         _qsLo;
  uint32         _qsHi;
  uint32         _dsLo;
  uint32         _dsHi;

  merCovering   *_IL;
  merList       *_ML;

  trapMatch     *_next;

  trapMatch(uint32         isunique,
            uint32         qsLo,
            uint32         qsHi,
            uint32         dsLo,
            uint32         dsHi,
            merCovering   *IL,
            merList       *ML) {
    _unique     = isunique;
    _qsLo       = qsLo;
    _qsHi       = qsHi;
    _dsLo       = dsLo;
    _dsHi       = dsHi;
    _IL         = IL;
    _ML         = ML;
    _next       = 0L;
  };
};




class hitMatrix {
public:
  hitMatrix(uint32 qsLen, uint32 q, uint32 qsIdx, logMsg *theLog);
  ~hitMatrix();

  void    addHits(uint32  qi,
                  uint64 *ps,
                  uint64  cn,
                  uint64  ad=0);

  void    sort_diagonal(void);
  void    sort_dsPos(void);

  void    filter(char      direction,
                 double    minHitCoverage,
                 uint32    minHitLength,
                 aHit    *&theOutput,
                 uint32   &theOutputPos,
                 uint32   &theOutputMax);
private:
  uint32                        _qsLen;   //  Seq Len of Q
  uint32                        _qsMers;  //  Valid mers in Q
  uint32                        _qsIdx;   //  Index of Q in the FastA

  //  Instead of building the lines during add(), we store
  //  the information used to build lines, and then build them
  //  in chain().  This was done to reduce simultaneous memory
  //  usage, as the lineArrayMap and etc take up considerable space.
  //
  uint32                        _hitsLen;
  uint32                        _hitsMax;
  diagonalLine                 *_hits;

  logMsg                       *_theLog;

  //  Making sense of the raw output from the search is not a trivial
  //  task for perl.  SMALL searches (dbEST vs 0.5MB sequence) used more
  //  than 4GB of memory in perl.
  //
  //  So, we bite the bullet and do it here.
  //
  //  _matches is a sorted linked list of the regions we have found.
  //  The list is kept in REVERSE order, as we usually add regions
  //  in the correct order (correct reverse order), occasionally
  //  we need to swap the last two.
  //
  //  The list is deleted in filter()
  //
  trapMatch                    *_matches;

  void addMatch(uint32         isunique,
                uint32         qsLo,
                uint32         qsHi,
                uint32         dsLo,
                uint32         dsHi,
                merCovering   *IL,
                merList       *ML);

};


inline
void
hitMatrix::addHits(uint32  qi,
                   uint64 *ps,
                   uint64  cn,
                   uint64  ad) {

  if ((_hitsLen + cn) >= _hitsMax) {
    _hitsMax = _hitsMax + _hitsMax + (uint32)cn;

    diagonalLine *h;
    try {
      h = new diagonalLine [_hitsMax];
    } catch (std::bad_alloc) {
      fprintf(stderr, "hitMatrix::addHits()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
      fprintf(stderr, "hitMatrix::addHits()-- have "uint32FMT" hits, tried to add "uint64FMT" more\n", _hitsLen, cn);
      exit(1);
    }

    memcpy(h, _hits, sizeof(uint64) * _hitsLen);
    delete [] _hits;
    _hits = h;
  }

  uint64  uniq = cn;
  if (ad > 0)
    uniq = ad;
  if (uniq > 0x000003ff)
    uniq = 0x000003ff;

  for (uint64 i=0; i<cn; i++) {
    _hits[_hitsLen].val.uniq = uniq;
    _hits[_hitsLen].val.qPos = qi;
    _hits[_hitsLen].val.dPos = ps[i];
    _hitsLen++;
  }
}



void
doSearch(searcherState       *state,
         query               *qry,
         bool                 rc);

void
doFilter(searcherState       *state,
         query               *qry);

void
doPolishS4(searcherState       *state,
           query               *qry);

void
doPolishDP(searcherState       *state,
           query               *qry);


uint32    configureFilter(double L,
                          double H,
                          double V,
                          aHit  *theHits,
                          uint32 theHitsLen);

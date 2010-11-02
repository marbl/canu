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
  u32bit           _KBmerSize;
  u32bit           _KBcompression;
  char            *_KBspacingTemplate;

  u32bit           _merSkip;

  u32bit           _numSearchThreads;

  bool             _doReverse;
  bool             _doForward;
  bool             _doValidation;
  char            *_doValidationFileName;

  bool             _doAlignments;

  double           _Lo;
  double           _Hi;
  double           _Va;

  u32bit           _maxDiagonal;

  double           _extendWeight;
  u32bit           _extendMinimum;
  u32bit           _extendMaximum;

  u32bit           _repeatThreshold;

  //  Minimums for hits
  double           _minHitCoverage;
  u32bit           _minHitLength;

  //  Minimums for matches
  u32bit           _minMatchIdentity;
  u32bit           _minMatchCoverage;

  //  Filtering of hits
  bool             _afEnabled;
  double           _afThreshold;
  u32bit           _afLength;
  u32bit           _afInit;

  //  Filtering and/or cleanup of matches
  u32bit           _discardExonLength;
  u32bit           _discardExonQuality;
  bool             _splitMatches;
  bool             _polishOptimally;

  char            *_dbFileName;
  char            *_psFileName;
  char            *_qsFileName;

  char            *_maskFileName;
  char            *_onlyFileName;

  bool             _buildOnly;

  u32bit           _ignoreThreshold;
  char            *_maskPrefix;
  u32bit           _maskThreshold;
  char            *_onlyPrefix;
  u32bit           _onlyThreshold;

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

extern volatile u32bit        numberOfQueries;



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
  u32bit     _status;
  u32bit     _qsIdx;
  u32bit     _dsIdx;
  u32bit     _dsLo;
  u32bit     _dsHi;
  u32bit     _covered;
  u32bit     _matched;
  u32bit     _numMers;
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
  u32bit       theHitsLen;
  u32bit       theHitsMax;
  aHit        *theHits;

  logMsg      *theLog;

  u32bit       theOutputLen;
  u32bit       theOutputMax;
  char        *theOutput;
};



class searcherState {
public:
  u64bit         threadID;

  u64bit         posnMax;
  u64bit         posnLen;
  u64bit        *posn;

  kMerBuilder   *KB;
  void          *DP;

  searcherState(u64bit U) {
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
    u64bit  all;
    struct {
      u64bit     uniq : 10;  //  uniqueness score for the mer here
      u64bit     qPos : 22;  //  position in the query, 0 to 4M
      u64bit     dPos : 32;  //  position in the genome, 0 to 4G
    } val;
  };
};



class trapMatch {
public:
  u32bit         _unique;
  u32bit         _qsLo;
  u32bit         _qsHi;
  u32bit         _dsLo;
  u32bit         _dsHi;

  merCovering   *_IL;
  merList       *_ML;

  trapMatch     *_next;

  trapMatch(u32bit         isunique,
            u32bit         qsLo,
            u32bit         qsHi,
            u32bit         dsLo,
            u32bit         dsHi,
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
  hitMatrix(u32bit qsLen, u32bit q, u32bit qsIdx, logMsg *theLog);
  ~hitMatrix();

  void    addHits(u32bit  qi,
                  u64bit *ps,
                  u64bit  cn,
                  u64bit  ad=0);

  void    sort_diagonal(void);
  void    sort_dsPos(void);

  void    filter(char      direction,
                 double    minHitCoverage,
                 u32bit    minHitLength,
                 aHit    *&theOutput,
                 u32bit   &theOutputPos,
                 u32bit   &theOutputMax);
private:
  u32bit                        _qsLen;   //  Seq Len of Q
  u32bit                        _qsMers;  //  Valid mers in Q
  u32bit                        _qsIdx;   //  Index of Q in the FastA

  //  Instead of building the lines during add(), we store
  //  the information used to build lines, and then build them
  //  in chain().  This was done to reduce simultaneous memory
  //  usage, as the lineArrayMap and etc take up considerable space.
  //
  u32bit                        _hitsLen;
  u32bit                        _hitsMax;
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

  void addMatch(u32bit         isunique,
                u32bit         qsLo,
                u32bit         qsHi,
                u32bit         dsLo,
                u32bit         dsHi,
                merCovering   *IL,
                merList       *ML);

};


inline
void
hitMatrix::addHits(u32bit  qi,
                   u64bit *ps,
                   u64bit  cn,
                   u64bit  ad) {

  if ((_hitsLen + cn) >= _hitsMax) {
    _hitsMax = _hitsMax + _hitsMax + (u32bit)cn;

    diagonalLine *h;
    try {
      h = new diagonalLine [_hitsMax];
    } catch (std::bad_alloc) {
      fprintf(stderr, "hitMatrix::addHits()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
      fprintf(stderr, "hitMatrix::addHits()-- have "u32bitFMT" hits, tried to add "u64bitFMT" more\n", _hitsLen, cn);
      exit(1);
    }

    memcpy(h, _hits, sizeof(u64bit) * _hitsLen);
    delete [] _hits;
    _hits = h;
  }

  u64bit  uniq = cn;
  if (ad > 0)
    uniq = ad;
  if (uniq > 0x000003ff)
    uniq = 0x000003ff;

  for (u64bit i=0; i<cn; i++) {
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


u32bit    configureFilter(double L,
                          double H,
                          double V,
                          aHit  *theHits,
                          u32bit theHitsLen);

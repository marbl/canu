#ifndef THR_H
#define THR_H


#include "posix.H"
#include "snapper2.H"
#include "bio++.H"
#include "sim4.H"

//  If you really, really, really want to know the exact number
//  of bases left in the query, use the interval list.  Otherwise,
//  it's faster to guess.
//
//#define USEEXACTSIZE

//  Define this to print a message whenever a search starts.
//
//#define VERBOSE_SEARCH

//  Define this to print a message whenever a polish starts.
//
//#define SHOW_POLISHING

//  Define these to show polishes that take a long time -- individual
//  polishes, not all polishes for a single sequence.  The time is in
//  seconds.
//
//#define SHOW_POLISHING_EXPENSIVE  0.5

//  Define this to create a bazillion files that contain information
//  about the signal for each match.  One file per query sequence.
//
//#define SAVE_HITS_TO_FILES

//  Define this to show the "hit discarding" results.  Any hits
//  that look like they are repeats are re-searched using nearly
//  unique mers.
//
//#define SHOW_HIT_DISCARDING

//  Define this to show the hits being added to the sim4command from
//  thr-polish.C.  This generates a lot of output!
//
//#define SHOW_HITS_ADDED
//#define SHOW_HITS_ADDED_AFTER_QUERY   4600


class searcherState {
public:
  u64bit         threadID;

  u64bit         posnMax;
  u64bit         posnLen;
  u64bit        *posn;

  double         maskTime;
  double         searchTime;
  double         chainTime;
  double         polishTime;

  u64bit         searched;
  u64bit         polished;

#ifdef SAVE_HITS_TO_FILES
  FILE          *hitsDump;
#endif

  //  This guy is owned by the main thread.  It's here for convenience.
  //

  searcherState(u64bit U) {
    threadID    = U;

    posnMax     = 16384;
    posnLen     = 0;
    posn        = new u64bit [ posnMax ];

    maskTime    = 0.0;
    searchTime  = 0.0;
    chainTime   = 0.0;
    polishTime  = 0.0;

    searched    = 0;
    polished    = 0;

#ifdef SAVE_HITS_TO_FILES
    hitsDump    = 0L;
#endif
  };

  ~searcherState() {
    delete [] posn;
  };
};




void
doSearch(searcherState       *state,
         FastASequenceInCore *seq,
         u32bit               idx,
         bool                 rc,
         aHit               *&theHits,
         u32bit              &theHitsLen,
         u32bit              &theHitsMax,
         logMsg              *theLog);



u32bit
doFilter(searcherState       *state,
         aHit               *&theHits,
         u32bit              &theHitsLen,
         logMsg              *theLog);



char*
doPolish(searcherState       *state,
         FastASequenceInCore *seq,
         aHit               *&theHits,
         u32bit              &theHitsLen,
         logMsg              *theLog);




#endif  //  THR_H

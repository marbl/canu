#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <fcntl.h>
#include <assert.h>
#include <math.h>

//  At one time, this was needed for pthread.h or semaphore.h
typedef unsigned short ushort;

#include <pthread.h>
#include <semaphore.h>

#include "bri++.H"

#include "fasta.H"
#include "fasta-cache.H"

#include "positionDB.H"
#include "existDB.H"
#include "sim4.H"
#include "hitMatrix.H"



class encodedQuery {
private:
  u32bit    _mersTotal;  //  Number of mers possible in the query
  u32bit    _mersAvail;  //  Number of mers not masked out
  u64bit   *_mers;       //  List of mers
  u32bit   *_skip;
public:
  encodedQuery(char const           *seq,
               u32bit                seqLen,
               u32bit                k,
               bool                  rc);
  ~encodedQuery();

  u32bit           numberOfMers(void)        { return(_mersTotal); };
  u32bit           numberOfValidMers(void)   { return(_mersAvail); };
  u64bit           getMer(u32bit i)          { return(_mers[i]); };

  u32bit           getSkip(u32bit i)         { return(_skip[i]); };
  void             setSkip(u32bit i)         {
    //  If skip[i] is already set, don't adjust; otherwise, subtract
    //  one from the available.
    _mersAvail -= 1 - _skip[i];
    _skip[i]    = 1;
  };
};














//
//  A singleton for working with the command line parameters.
//
#define MAX_THREADS  64

//  Stuff for -use
//
struct use_s {
  u32bit   seq;
  u32bit   size;
  u32bit   start;
};

class configuration {
public:
  bool             _beVerbose;

  u32bit           _merSize;
  u32bit           _merSkip;

  u32bit           _numSearchThreads;

  bool             _doReverse;
  bool             _doForward;
  bool             _doValidation;
  char            *_doValidationFileName;

  double           _Lo;
  double           _Hi;
  double           _Va;

  u32bit           _maxDiagonal;

  double           _extendWeight;
  u32bit           _extendMinimum;

  //  Minimums for hits
  u32bit           _minHitLength;
  double           _minHitCoverage;

  //  Minimums for matches
  u32bit           _minMatchIdentity;
  u32bit           _minMatchCoverage;

  char            *_dbFileName;
  char            *_psFileName;
  char            *_qsFileName;

  char            *_maskFileName;
  char            *_onlyFileName;

  char            *_maskPrefix;
  u32bit           _maskThreshold;
  char            *_onlyPrefix;
  u32bit           _onlyThreshold;

  char            *_outputFileName;
  char            *_statsFileName;

  u32bit           _useListLen;
  u32bit           _useListMax;
  use_s           *_useList;

  //  Wall clock times
  //
  double           _startTime;
  double           _initTime;
  double           _buildTime;
  double           _searchTime;
  double           _totalTime;

  //  Loader parameters
  //
  u32bit           _loaderHighWaterMark;
  struct timespec  _loaderSleep;
  bool             _loaderWarnings;

  //  Search parameters
  //
  struct timespec  _searchSleep;

  //  Output parameters
  //
  u32bit           _writerHighWaterMark;
  struct timespec  _writerSleep;
  bool             _writerWarnings;

  configuration();
  ~configuration();

  void  usage(char *name);
  void  read(int argc, char **argv);
  void  display(FILE *out=stdout);

  void  completeUseList(u32bit numSeqs);

  void  setTime(struct timespec *ts, double t) {
    ts->tv_sec  = (time_t)floor(t);
    ts->tv_nsec = (long)((t - ts->tv_sec) * 1e9);
  };
private:
  void  addToUse(u32bit v);
  void  parseUseLine(char *line);
};







//  Shared data
//
extern configuration          config;
extern sim4parameters         sim4params;
extern FastACache            *cache;

extern FastAWrapper          *qsFASTA;  //  Used exclusively by thr-loader.C

extern existDB               *maskDB;  //  thr-search.C
extern existDB               *onlyDB;  //  thr-search.C

extern positionDB            *positions;


extern volatile u32bit        numberOfQueries;

extern aHit                 **answer;
extern u32bit                *answerLen;
extern char                 **output;
extern u32bit                *outputLen;

extern pthread_mutex_t        inputTailMutex;
extern FastASequenceInCore  **input;
extern volatile u32bit        inputHead;
extern volatile u32bit        inputTail;

extern volatile u32bit        outputPos;

extern char                  *threadStats[MAX_THREADS];

void     *deadlockDetector(void *U);
void     *deadlockChecker(void *U);

void     *loaderThread(void *U);
void     *searchThread(void *U);

u32bit    configureFilter(double L,
                          double H,
                          double V,
                          aHit  *theHits,
                          u32bit theHitsLen);

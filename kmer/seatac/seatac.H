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

#include "libbri.H"
#include "fasta.H"

#include "positionDB.H"
#include "existDB.H"

#include "hitMatrix.H"

#include "filterObj.H"
#include "statObj.H"

class encodedQuery {
private:
  char const  *_seq;
  u32bit       _seqLen;
  u32bit       _merSize;
  bool         _rc;

  u32bit       _seqPos;

  u64bit       _substring;
  u64bit       _mermask;
  s32bit       _timeUntilValid;
public:
  encodedQuery(char const           *seq,
               u32bit                seqLen,
               u32bit                k,
               bool                  rc);
  ~encodedQuery();

  bool   getMer(u64bit &mer, u32bit &pos);
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
  u32bit           _numSearchThreads;

  u32bit           _merSkip;

  bool             _doReverse;
  bool             _doForward;

  u32bit           _maxDiagonal;
  u32bit           _maxGap;
  u32bit           _qsOverlap;
  u32bit           _dsOverlap;

  u32bit           _minLength;

  char            *_dbFileName;
  char            *_qsFileName;
  char            *_maskFileName;
  char            *_onlyFileName;
  char            *_tmpFileName;
  char            *_outputFileName;
  char            *_statsFileName;

  char            *_tableFileName;

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

  //  Filter parameters
  //
  char            *_filterpath;
  sharedObj       *_filtername;
  char            *_filteropts;

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

  void  completeUseList(FastAWrapper *db);

  void  setTime(struct timespec *ts, double t) {
    ts->tv_sec = floor(t);
    t -= ts->tv_sec;
    t *= 1e9;
    ts->tv_nsec = t;
  };
private:
  void  addToUse(u32bit v);
  void  parseUseLine(char *line);
};







//  Shared data
//
extern configuration          config;

extern FastAWrapper          *qsFASTA;  //  Used exclusively by thr-loader.C

extern positionDB            *positions;

extern volatile u32bit        numberOfQueries;

extern filterObj            **output;

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

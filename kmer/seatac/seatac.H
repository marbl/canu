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

#include "bio++.H"

#include "seqCache.H"

#include "existDB.H"
#include "positionDB.H"

#include "hitMatrix.H"

#include "filterObj.H"
#include "statObj.H"

class encodedQuery {
private:
  char const  *_seq;
  uint32       _seqLen;
  uint32       _merSize;
  bool         _rc;

  uint32       _seqPos;

  int32       _timeUntilValid;

  uint64       _substring;
  uint64       _mermask;
public:
  encodedQuery(char const           *seq,
               uint32                seqLen,
               uint32                k,
               bool                  rc);
  ~encodedQuery();

  bool   getMer(uint64 &mer, uint32 &pos);
};




//
//  A singleton for working with the command line parameters.
//
#define MAX_THREADS  64


class configuration {
public:
  bool             _beVerbose;

  uint32           _merSize;
  uint32           _merSkip;
  uint32           _numSearchThreads;

  bool             _doReverse;
  bool             _doForward;

  uint32           _maxDiagonal;
  uint32           _maxGap;
  uint32           _qsOverlap;
  uint32           _dsOverlap;

  uint32           _minLength;

  char            *_dbFileName;
  char            *_qsFileName;
  char            *_maskFileName;
  char            *_onlyFileName;
  char            *_outputFileName;
  char            *_statsFileName;

  char            *_tableFileName;
  bool             _tableBuildOnly;

  seqStream       *_genome;

  //  Filter parameters
  //
  char            *_filtername;
  char            *_filteropts;
  sharedObj       *_filterObj;

  //  Wall clock times
  //
  double           _startTime;
  double           _initTime;
  double           _buildTime;
  double           _searchTime;
  double           _totalTime;

  //  Loader parameters
  //
  uint32           _loaderHighWaterMark;
  struct timespec  _loaderSleep;
  bool             _loaderWarnings;

  //  Search parameters
  //
  struct timespec  _searchSleep;

  //  Output parameters
  //
  uint32           _writerHighWaterMark;
  struct timespec  _writerSleep;
  bool             _writerWarnings;

  configuration();
  ~configuration();

  void  usage(char *name);
  void  read(int argc, char **argv);
  void  writeATACheader(FILE *out);

  void  setTime(struct timespec *ts, double t) {
    ts->tv_sec  = (time_t)floor(t);
    ts->tv_nsec = (long)((t - ts->tv_sec) * 1e9);
  };
};







//  Shared data
//
extern configuration          config;

extern seqCache              *qsFASTA;  //  Used exclusively by thr-loader.C

extern positionDB            *positions;

extern volatile uint32        numberOfQueries;

extern filterObj            **output;

extern pthread_mutex_t        inputTailMutex;
extern seqInCore            **input;
extern volatile uint32        inputHead;
extern volatile uint32        inputTail;

extern volatile uint32        outputPos;

extern char                  *threadStats[MAX_THREADS];

void     *deadlockDetector(void *U);
void     *deadlockChecker(void *U);

void     *loaderThread(void *U);
void     *searchThread(void *U);

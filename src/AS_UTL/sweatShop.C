
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
 *  This file is derived from:
 *
 *    kmer/libutil/sweatShop.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2006-MAR-02 to 2014-APR-11
 *      are Copyright 2006,2008,2010-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2015-JUN-24
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sweatShop.H"
#include "timeAndSize.H"

#include <sched.h>  //  pthread scheduling stuff


class sweatShopWorker {
public:
  sweatShopWorker() {
    shop            = 0L;
    threadUserData  = 0L;
    numComputed     = 0;
    workerQueue     = 0L;
    workerQueueLen  = 0L;
  };

  sweatShop        *shop;
  void             *threadUserData;
  pthread_t         threadID;
  uint32            numComputed;
  sweatShopState  **workerQueue;
  uint32            workerQueueLen;
};


//  This gets created by the loader, passed to the worker, and printed
//  by the writer.  userData is controlled by the user.
//
class sweatShopState {
public:
  sweatShopState(void *userData) {
    _user     = userData;
    _computed = false;
    _next     = 0L;
  };
  ~sweatShopState() {
  };

  void             *_user;
  bool              _computed;
  sweatShopState   *_next;
};




//  Simply forwards control to the class
void*
_sweatshop_loaderThread(void *ss_) {
  sweatShop *ss = (sweatShop *)ss_;
  return(ss->loader());
}

void*
_sweatshop_workerThread(void *sw_) {
  sweatShopWorker *sw = (sweatShopWorker *)sw_;
  return(sw->shop->worker(sw));
}

void*
_sweatshop_writerThread(void *ss_) {
  sweatShop *ss = (sweatShop *)ss_;
  return(ss->writer());
}

void*
_sweatshop_statusThread(void *ss_) {
  sweatShop *ss = (sweatShop *)ss_;
  return(ss->status());
}



sweatShop::sweatShop(void*(*loaderfcn)(void *G),
                     void (*workerfcn)(void *G, void *T, void *S),
                     void (*writerfcn)(void *G, void *S)) {

  _userLoader       = loaderfcn;
  _userWorker       = workerfcn;
  _userWriter       = writerfcn;

  _globalUserData   = 0L;

  _writerP          = 0L;
  _workerP          = 0L;
  _loaderP          = 0L;

  _showStatus       = false;

  _loaderQueueSize  = 1024;
  _loaderQueueMax   = 10240;
  _loaderQueueMin   = 4;  //  _numberOfWorkers * 2, reset when that changes
  _loaderBatchSize  = 1;
  _workerBatchSize  = 1;
  _writerQueueSize  = 4096;
  _writerQueueMax   = 10240;

  _numberOfWorkers  = 2;

  _workerData       = 0L;

  _numberLoaded     = 0;
  _numberComputed   = 0;
  _numberOutput     = 0;
}


sweatShop::~sweatShop() {
  delete [] _workerData;
}



void
sweatShop::setThreadData(uint32 t, void *x) {
  if (_workerData == 0L)
    _workerData = new sweatShopWorker [_numberOfWorkers];

  if (t >= _numberOfWorkers)
    fprintf(stderr, "sweatShop::setThreadData()-- worker ID "F_U32" more than number of workers="F_U32"\n", t, _numberOfWorkers), exit(1);

  _workerData[t].threadUserData = x;
}



//  Build a list of states to add in one swoop
//
void
sweatShop::loaderSave(sweatShopState *&tail, sweatShopState *&head, sweatShopState *thisState) {

  thisState->_next  = 0L;

  if (tail) {
    head->_next = thisState;
    head        = thisState;
  } else {
    tail = head = thisState;
  }
  _numberLoaded++;
}


//  Add a bunch of new states to the queue.
//
void
sweatShop::loaderAppend(sweatShopState *&tail, sweatShopState *&head) {
  int err;

  if ((tail == 0L) || (head == 0L))
    return;

  err = pthread_mutex_lock(&_stateMutex);
  if (err != 0)
    fprintf(stderr, "sweatShop::loaderAppend()--  Failed to lock mutex (%d).  Fail.\n", err), exit(1);

  if (_loaderP == 0L) {
    _writerP      = tail;
    _workerP      = tail;
    _loaderP      = head;
  } else {
    _loaderP->_next = tail;
  }
  _loaderP        = head;

  err = pthread_mutex_unlock(&_stateMutex);
  if (err != 0)
    fprintf(stderr, "sweatShop::loaderAppend()--  Failed to unlock mutex (%d).  Fail.\n", err), exit(1);

  tail = 0L;
  head = 0L;
}



void*
sweatShop::loader(void) {

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 166666666ULL;  //  1/6 second

  //  We can batch several loads together before we push them onto the
  //  queue, this should reduce the number of times the loader needs to
  //  lock the queue.
  //
  //  But it also increases the latency, so it's disabled by default.
  //
  sweatShopState        *tail       = 0L;  //  The first thing loaded
  sweatShopState        *head       = 0L;  //  The last thing loaded
  uint32                 numLoaded  = 0;

  bool  moreToLoad = true;

  while (moreToLoad) {

    //  Zzzzzzz....
    while (_numberLoaded > _numberComputed + _loaderQueueSize)
      nanosleep(&naptime, 0L);

    sweatShopState  *thisState = new sweatShopState((*_userLoader)(_globalUserData));

    //  If we actually loaded a new state, add it
    //
    if (thisState->_user) {
      loaderSave(tail, head, thisState);
      numLoaded++;
      if (numLoaded >= _loaderBatchSize)
        loaderAppend(tail, head);
    } else {
      //  Didn't read, must be all done!  Push on the end-of-input marker state.
      //
      loaderSave(tail, head, new sweatShopState(0L));
      loaderAppend(tail, head);

      moreToLoad = false;
      delete thisState;
    }
  }

  //fprintf(stderr, "sweatShop::reader exits.\n");
  return(0L);
}



void*
sweatShop::worker(sweatShopWorker *workerData) {

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 50000000ULL;

  bool    moreToCompute = true;
  int     err;

  while (moreToCompute) {

    //  Usually beacuse some worker is taking a long time, and the
    //  output queue isn't big enough.
    //
    while (_numberOutput + _writerQueueSize < _numberComputed)
      nanosleep(&naptime, 0L);

    //  Grab the next state.  We don't grab it if it's the last in the
    //  queue (else we would fall off the end) UNLESS it really is the
    //  last one.
    //
    err = pthread_mutex_lock(&_stateMutex);
    if (err != 0)
      fprintf(stderr, "sweatShop::worker()--  Failed to lock mutex (%d).  Fail.\n", err), exit(1);

    for (workerData->workerQueueLen = 0; ((workerData->workerQueueLen < _workerBatchSize) &&
                                          (_workerP) &&
                                          ((_workerP->_next != 0L) || (_workerP->_user == 0L))); workerData->workerQueueLen++) {
      workerData->workerQueue[workerData->workerQueueLen] = _workerP;
      _workerP = _workerP->_next;
    }

    if (_workerP == 0L)
      moreToCompute = false;

    err = pthread_mutex_unlock(&_stateMutex);
    if (err != 0)
      fprintf(stderr, "sweatShop::worler()--  Failed to lock mutex (%d).  Fail.\n", err), exit(1);


    if (workerData->workerQueueLen == 0) {
      //  No work, sleep a bit to prevent thrashing the mutex and resume.
      nanosleep(&naptime, 0L);
      continue;
    }

    //  Execute
    //
    for (uint32 x=0; x<workerData->workerQueueLen; x++) {
      sweatShopState *ts = workerData->workerQueue[x];

      if (ts && ts->_user) {
        (*_userWorker)(_globalUserData, workerData->threadUserData, ts->_user);
        ts->_computed = true;
        workerData->numComputed++;
      } else {
        //  When we really do run out of stuff to do, we'll end up here
        //  (only one thread will end up in the other case, with
        //  something to do and moreToCompute=false).  If it's actually
        //  the end, skip the sleep and just get outta here.
        //
        if (moreToCompute == true) {
          fprintf(stderr, "WARNING!  Worker is sleeping because the reader is slow!\n");
          nanosleep(&naptime, 0L);
        }
      }
    }
  }

  //fprintf(stderr, "sweatShop::worker exits.\n");
  return(0L);
}


void*
sweatShop::writer(void) {
  sweatShopState  *deleteState = 0L;

  //  Wait for output to appear, then write.
  //
  while (_writerP && _writerP->_user) {

    if        (_writerP->_computed == false) {
      //  Wait for a slow computation.
      struct timespec   naptime;
      naptime.tv_sec      = 0;
      naptime.tv_nsec     = 5000000ULL;

      //fprintf(stderr, "Writer waits for slow thread at "F_U64".\n", _numberOutput);
      nanosleep(&naptime, 0L);
    } else if (_writerP->_next == 0L) {
      //  Wait for the input.
      struct timespec   naptime;
      naptime.tv_sec      = 0;
      naptime.tv_nsec     = 5000000ULL;

      //fprintf(stderr, "Writer waits for all threads at "F_U64".\n", _numberOutput);
      nanosleep(&naptime, 0L);
    } else {
      (*_userWriter)(_globalUserData, _writerP->_user);
      _numberOutput++;

      deleteState = _writerP;
      _writerP    = _writerP->_next;
      delete deleteState;
    }
  }

  //  Tell status to stop.
  _writerP = 0L;

  //fprintf(stderr, "sweatShop::writer exits.\n");
  return(0L);
}


//  This thread not only shows a status message, but it also updates the critical shared variable
//  _numberComputed.  Worker threads use this to throttle themselves.  Thus, even if _showStatus is
//  not set, and this thread doesn't _appear_ to be doing anything useful....it is.
//
void*
sweatShop::status(void) {

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 250000000ULL;

  double  startTime = getTime() - 0.001;
  double  thisTime  = 0;

  uint64  deltaOut = 0;
  uint64  deltaCPU = 0;

  double  cpuPerSec = 0;

  uint64  readjustAt = 16384;

  while (_writerP) {
    uint32 nc = 0;
    for (uint32 i=0; i<_numberOfWorkers; i++)
      nc += _workerData[i].numComputed;
    _numberComputed = nc;

    deltaOut = deltaCPU = 0;

    thisTime = getTime();

    if (_numberComputed > _numberOutput)
      deltaOut = _numberComputed - _numberOutput;
    if (_numberLoaded > _numberComputed)
      deltaCPU = _numberLoaded - _numberComputed;

    cpuPerSec = _numberComputed / (thisTime - startTime);

    if (_showStatus) {
      fprintf(stderr, " %6.1f/s - %8"F_U64P" loaded; %8"F_U64P" queued for compute; %08"F_U64P" finished; %8"F_U64P" written; %8"F_U64P" queued for output)\r",
              cpuPerSec, _numberLoaded, deltaCPU, _numberComputed, _numberOutput, deltaOut);
      fflush(stderr);
    }

    //  Readjust queue sizes based on current performance, but don't let it get too big or small.
    //  In particular, don't let it get below 2*numberOfWorkers.
    //
     if (_numberComputed > readjustAt) {
       readjustAt       += (uint64)(2 * cpuPerSec);
       _loaderQueueSize  = (uint32)(5 * cpuPerSec);
     }

    if (_loaderQueueSize < _loaderQueueMin)
      _loaderQueueSize = _loaderQueueMin;

    if (_loaderQueueSize < 2 * _numberOfWorkers)
      _loaderQueueSize = 2 * _numberOfWorkers;

    if (_loaderQueueSize > _loaderQueueMax)
      _loaderQueueSize = _loaderQueueMax;

    nanosleep(&naptime, 0L);
  }

  if (_showStatus) {
    thisTime = getTime();

    if (_numberComputed > _numberOutput)
      deltaOut = _numberComputed - _numberOutput;
    if (_numberLoaded > _numberComputed)
      deltaCPU = _numberLoaded - _numberComputed;

    cpuPerSec = _numberComputed / (thisTime - startTime);

    fprintf(stderr, " %6.1f/s - %08"F_U64P" queued for compute; %08"F_U64P" finished; %08"F_U64P" queued for output)\n",
            cpuPerSec, deltaCPU, _numberComputed, deltaOut);
  }

  //fprintf(stderr, "sweatShop::status exits.\n");
  return(0L);
}





void
sweatShop::run(void *user, bool beVerbose) {
  pthread_attr_t      threadAttr;
  pthread_t           threadIDloader;
  pthread_t           threadIDwriter;
  pthread_t           threadIDstats;
#if 0
  int                 threadSchedPolicy = 0;
  struct sched_param  threadSchedParamDef;
  struct sched_param  threadSchedParamMax;
#endif
  int                 err = 0;

  _globalUserData = user;
  _showStatus     = beVerbose;

  //  Configure everything ahead of time.

  if (_workerBatchSize < 1)
    _workerBatchSize = 1;

  if (_workerData == 0L)
    _workerData = new sweatShopWorker [_numberOfWorkers];

  for (uint32 i=0; i<_numberOfWorkers; i++) {
    _workerData[i].shop        = this;
    _workerData[i].workerQueue = new sweatShopState * [_workerBatchSize];
  }

  //  Open the doors.

  errno = 0;

  err = pthread_mutex_init(&_stateMutex, NULL);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to configure pthreads (state mutex): %s.\n", strerror(err)), exit(1);

  err = pthread_attr_init(&threadAttr);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to configure pthreads (attr init): %s.\n", strerror(err)), exit(1);

  err = pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to configure pthreads (set scope): %s.\n", strerror(err)), exit(1);

  err = pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to configure pthreads (joinable): %s.\n", strerror(err)), exit(1);

#if 0
  err = pthread_attr_getschedparam(&threadAttr, &threadSchedParamDef);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to configure pthreads (get default param): %s.\n", strerror(err)), exit(1);

  err = pthread_attr_getschedparam(&threadAttr, &threadSchedParamMax);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to configure pthreads (get max param): %s.\n", strerror(err)), exit(1);
#endif

  //  SCHED_RR needs root privs to run on FreeBSD.
  //
  //err = pthread_attr_setschedpolicy(&threadAttr, SCHED_RR);
  //if (err)
  //  fprintf(stderr, "sweatShop::run()--  Failed to configure pthreads (sched policy): %s.\n", strerror(err)), exit(1);

#if 0
  err = pthread_attr_getschedpolicy(&threadAttr, &threadSchedPolicy);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to configure pthreads (sched policy): %s.\n", strerror(err)), exit(1);

  errno = 0;
  threadSchedParamMax.sched_priority = sched_get_priority_max(threadSchedPolicy);
  if (errno)
    fprintf(stderr, "sweatShop::run()--  WARNING: Failed to configure pthreads (set max param priority): %s.\n", strerror(errno));

  //  Fire off the loader

  err = pthread_attr_setschedparam(&threadAttr, &threadSchedParamMax);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to set loader priority: %s.\n", strerror(err)), exit(1);
#endif

  err = pthread_create(&threadIDloader, &threadAttr, _sweatshop_loaderThread, this);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to launch loader thread: %s.\n", strerror(err)), exit(1);

  //  Wait for it to actually load something (otherwise all the
  //  workers immediately go home)

  while (!_writerP && !_workerP && !_loaderP) {
    struct timespec   naptime;
    naptime.tv_sec      = 0;
    naptime.tv_nsec     = 250000ULL;
    nanosleep(&naptime, 0L);
  }

  //  Start the statistics and writer

#if 0
  err = pthread_attr_setschedparam(&threadAttr, &threadSchedParamMax);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to set status and writer priority: %s.\n", strerror(err)), exit(1);
#endif

  err = pthread_create(&threadIDstats,  &threadAttr, _sweatshop_statusThread, this);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to launch status thread: %s.\n", strerror(err)), exit(1);

  err = pthread_create(&threadIDwriter, &threadAttr, _sweatshop_writerThread, this);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to launch writer thread: %s.\n", strerror(err)), exit(1);

  //  And some labor

#if 0
  err = pthread_attr_setschedparam(&threadAttr, &threadSchedParamDef);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to set worker priority: %s.\n", strerror(err)), exit(1);
#endif

  for (uint32 i=0; i<_numberOfWorkers; i++) {
    err = pthread_create(&_workerData[i].threadID, &threadAttr, _sweatshop_workerThread, _workerData + i);
    if (err)
      fprintf(stderr, "sweatShop::run()--  Failed to launch worker thread "F_U32": %s.\n", i, strerror(err)), exit(1);
  }

  //  Now sit back and relax.

  err = pthread_join(threadIDloader, 0L);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to join loader thread: %s.\n", strerror(err)), exit(1);

  err = pthread_join(threadIDwriter, 0L);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to join writer thread: %s.\n", strerror(err)), exit(1);

  err = pthread_join(threadIDstats,  0L);
  if (err)
    fprintf(stderr, "sweatShop::run()--  Failed to join status thread: %s.\n", strerror(err)), exit(1);

  for (uint32 i=0; i<_numberOfWorkers; i++) {
    err = pthread_join(_workerData[i].threadID, 0L);
    if (err)
      fprintf(stderr, "sweatShop::run()--  Failed to join worker thread "F_U32": %s.\n", i, strerror(err)), exit(1);
  }

  //  Cleanup.

  delete _loaderP;
  _loaderP = _workerP = _writerP = 0L;
}

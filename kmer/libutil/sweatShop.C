#include "sweatShop.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

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
  u32bit            numComputed;
  sweatShopState  **workerQueue;
  u32bit            workerQueueLen;
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
_sweatshop_loaderThread(void *ss) {
  sweatShop *SS = (sweatShop *)ss;
  return(SS->loader());
}

void*
_sweatshop_workerThread(void *x) {
  sweatShopWorker *SW = (sweatShopWorker *)x;
  return(SW->shop->worker(SW));
}

void*
_sweatshop_writerThread(void *x) {
  sweatShop *SS = (sweatShop *)x;
  return(SS->writer());
}

void*
_sweatshop_statusThread(void *x) {
  sweatShop *SS = (sweatShop *)x;
  return(SS->status());
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
  _loaderBatchSize  = 1;
  _workerBatchSize  = 1;
  _writerQueueSize  = 128;

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
sweatShop::setThreadData(u32bit t, void *x) {
  if (_workerData == 0L)
    _workerData = new sweatShopWorker [_numberOfWorkers];

  if (t >= _numberOfWorkers)
    fprintf(stderr, "sweatShop::setThreadData()-- worker ID "u32bitFMT" more than number of workers="u32bitFMT"\n", t, _numberOfWorkers), exit(1);

  _workerData[t].threadUserData = x;
}



//  Add a new state to the head of the queue.  This is just gross.
//
#if 0
void
sweatShop::loaderAdd(sweatShopState *thisState) {
  pthread_mutex_lock(&_stateMutex);
  thisState->_next  = 0L;
  if (_loaderP == 0L) {
    _writerP      = thisState;
    _workerP      = thisState;
    _loaderP      = thisState;
  } else {
    _loaderP->_next = thisState;
  }
  _loaderP        = thisState;
  pthread_mutex_unlock(&_stateMutex);
  _numberLoaded++;
}
#endif

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

  if ((tail == 0L) || (head == 0L))
    return;

  pthread_mutex_lock(&_stateMutex);

  if (_loaderP == 0L) {
    _writerP      = tail;
    _workerP      = tail;
    _loaderP      = head;
  } else {
    _loaderP->_next = tail;
  }
  _loaderP        = head;

  pthread_mutex_unlock(&_stateMutex);

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
  u32bit                 numLoaded  = 0;

  bool  moreToLoad = true;

  while (moreToLoad) {

    //  Zzzzzzz....
    while (_numberLoaded - _numberComputed > _loaderQueueSize)
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
  naptime.tv_sec      = 2;
  naptime.tv_nsec     = 50000000ULL;

  bool    moreToCompute = true;

  while (moreToCompute) {

    //  Usually beacuse some worker is taking a long time, and the
    //  output queue isn't big enough.
    //
    while (_numberComputed - _numberOutput > _writerQueueSize)
      nanosleep(&naptime, 0L);

    //fprintf(stderr, "Worker %0x16p starting.\n", workerData);

    //  Grab the next state.  We don't grab it if it's the last in the
    //  queue (else we would fall off the end) UNLESS it really is the
    //  last one.
    //
    pthread_mutex_lock(&_stateMutex);

    for (workerData->workerQueueLen = 0; ((workerData->workerQueueLen < _workerBatchSize) &&
                                          (_workerP) &&
                                          ((_workerP->_next != 0L) || (_workerP->_user == 0L))); workerData->workerQueueLen++) {
      workerData->workerQueue[workerData->workerQueueLen] = _workerP;
      _workerP = _workerP->_next;
    }

    if (_workerP == 0L)
      moreToCompute = false;

    pthread_mutex_unlock(&_stateMutex);

    //  Execute
    //
    for (u32bit x=0; x<workerData->workerQueueLen; x++) {
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
      naptime.tv_nsec     = 50000000ULL;

      //fprintf(stderr, "Writer waits for slow thread at "u64bitFMT".\n", _numberOutput);
      nanosleep(&naptime, 0L);
    } else if (_writerP->_next == 0L) {
      //  Wait for the input.
      struct timespec   naptime;
      naptime.tv_sec      = 50000000;
      naptime.tv_nsec     = 0ULL;

      //fprintf(stderr, "Writer waits for all threads at "u64bitFMT".\n", _numberOutput);
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


void*
sweatShop::status(void) {

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 250000000ULL;

  double  startTime = getTime() - 0.001;
  double  thisTime  = 0;

  u64bit  deltaOut = 0;
  u64bit  deltaCPU = 0;

  double  cpuPerSec = 0;

  u64bit  readjustAt = 16384;

  while (_writerP) {
    u32bit nc = 0;
    for (u32bit i=0; i<_numberOfWorkers; i++)
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
      fprintf(stderr, " %6.1f/s - "u64bitFMTW(8)" loaded; "u64bitFMTW(8)" queued for compute; "u64bitFMTW(8)" finished; "u64bitFMTW(8)" written; "u64bitFMTW(8)" queued for output)\r",
              cpuPerSec, _numberLoaded, deltaCPU, _numberComputed, _numberOutput, deltaOut);
      fflush(stderr);
    }

    //  Readjust queue sizes based on current performance.
    //
    if (_numberComputed > readjustAt) {
      readjustAt       += (u64bit)(2 * cpuPerSec);
      _loaderQueueSize  = (u32bit)(5 * cpuPerSec);
    }

    nanosleep(&naptime, 0L);
  }

  if (_showStatus) {
    thisTime = getTime();

    if (_numberComputed > _numberOutput)
      deltaOut = _numberComputed - _numberOutput;
    if (_numberLoaded > _numberComputed)
      deltaCPU = _numberLoaded - _numberComputed;

    cpuPerSec = _numberComputed / (thisTime - startTime);

    fprintf(stderr, " %6.1f/s - "u64bitFMTW(8)" queued for compute; "u64bitFMTW(8)" finished; "u64bitFMTW(8)" queued for output)\n",
            cpuPerSec, deltaCPU, _numberComputed, deltaOut);
  }

  fprintf(stderr, "sweatShop::status exits.\n");
  return(0L);
}





void
sweatShop::run(void *user, bool beVerbose) {
  pthread_attr_t      threadAttr;
  pthread_t           threadIDloader;
  pthread_t           threadIDwriter;
  pthread_t           threadIDstats;
  struct sched_param  threadSchedParamDef;
  struct sched_param  threadSchedParamMax;

  _globalUserData = user;
  _showStatus     = beVerbose;

  pthread_mutex_init(&_stateMutex, NULL);

  pthread_attr_init(&threadAttr);
  pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);
  pthread_attr_setschedpolicy(&threadAttr, SCHED_RR);

  pthread_attr_getschedparam(&threadAttr, &threadSchedParamDef);
  pthread_attr_getschedparam(&threadAttr, &threadSchedParamMax);

  threadSchedParamMax.sched_priority = sched_get_priority_max(SCHED_RR);

  //  Fire off the loader

  pthread_attr_setschedparam(&threadAttr, &threadSchedParamMax);

  pthread_create(&threadIDloader, &threadAttr, _sweatshop_loaderThread, this);

  //  Wait for it to actually load something (otherwise all the
  //  workers immediately go home)

  while (!_writerP && !_workerP && !_loaderP) {
    struct timespec   naptime;
    naptime.tv_sec      = 0;
    naptime.tv_nsec     = 250000ULL;
    nanosleep(&naptime, 0L);
  }

  //  Start the statistics and writer

  pthread_attr_setschedparam(&threadAttr, &threadSchedParamMax);

  pthread_create(&threadIDstats,  &threadAttr, _sweatshop_statusThread, this);
  pthread_create(&threadIDwriter, &threadAttr, _sweatshop_writerThread, this);

  //  And some labor

  if (_workerBatchSize < 1)
    _workerBatchSize = 1;

  if (_workerData == 0L)
    _workerData = new sweatShopWorker [_numberOfWorkers];

  pthread_attr_setschedparam(&threadAttr, &threadSchedParamDef);

  for (u32bit i=0; i<_numberOfWorkers; i++) {
    _workerData[i].shop        = this;
    _workerData[i].workerQueue = new sweatShopState * [_workerBatchSize];

    pthread_create(&_workerData[i].threadID, &threadAttr, _sweatshop_workerThread, _workerData + i);
  }

  //  Now sit back and relax.

  pthread_join(threadIDloader, 0L);
  pthread_join(threadIDwriter, 0L);
  pthread_join(threadIDstats,  0L);

  for (u32bit i=0; i<_numberOfWorkers; i++)
    pthread_join(_workerData[i].threadID, 0L);

  //  Cleanup.

  delete _loaderP;
  _loaderP = _workerP = _writerP = 0L;
}

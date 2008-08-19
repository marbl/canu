#include "sweatShop.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>


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

#if 0
//  unused, we do this in 'main()'
void*
_sweatshop_writerThread(void *x) {
  sweatShop *SS = (sweatShop *)x;
  return(SS->writer());
}
#endif

void*
_sweatshop_statusThread(void *x) {
  sweatShop *SS = (sweatShop *)x;
  return(SS->status());
}



sweatShop::sweatShop(void*(*loader)(void *G),
                     void (*worker)(void *G, void *T, void *S),
                     void (*writer)(void *G, void *S)) {
  _userLoader = loader;
  _userWorker = worker;
  _userWriter = writer;

  _globalUserData = 0L;

  _writerP = 0L;  //  Where output takes stuff from, the tail
  _workerP = 0L;  //  Where computes happen, the middle
  _loaderP = 0L;  //  Where input is put, the head

  _loadBatches = false;

  _loaderQueueSize  =  128 * 1024;
  _writerQueueSize  =   64 * 1024;

  _numberOfWorkers = 3;

  _workerData     = new sweatShopWorker [_numberOfWorkers];
  _workerDataSet  = false;

  _numberLoaded   = 0;
  _numberComputed = 0;
  _numberOutput   = 0;
}


sweatShop::~sweatShop() {
}


u32bit
sweatShop::loaderQueueSize(u32bit x) {
  if (x > 0)
    _loaderQueueSize = x;
  return(_loaderQueueSize);
}

u32bit
sweatShop::writerQueueSize(u32bit x) {
  if (x > 0)
    _writerQueueSize = x;
  return(_writerQueueSize);
}

u32bit
sweatShop::numberOfWorkers(u32bit x) {
  if (x > 0) {
    if (_workerDataSet) {
      fprintf(stderr, "sweatShop::numberOfWorkers()-- ERROR!  This cannot be called after setThreadData()\n");
      exit(1);
    }
    _numberOfWorkers = x;
    delete [] _workerData;
    _workerData    = new sweatShopWorker [_numberOfWorkers];
  }
  return(_numberOfWorkers);
}

void
sweatShop::setThreadData(u32bit t, void *x) {
  if (t >= _numberOfWorkers) {
    fprintf(stderr, "sweatShop::setThreadData()-- worker ID "u32bitFMT" more than number of workers="u32bitFMT"\n",
            t, _numberOfWorkers);
    exit(1);
  }
  _workerData[t].threadUserData = x;
  _workerDataSet = true;
}


//  Add a new state to the head of the queue.  This is just gross.
//
void
sweatShop::loaderAdd(state_s *thisState) {
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


//  Build a list of states to add in one swoop
//
void
sweatShop::loaderSave(state_s *&tail, state_s *&head, state_s *thisState) {

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
sweatShop::loaderAppend(state_s *&tail, state_s *&head) {

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
  state_s               *tail = 0L;  //  The first thing loaded
  state_s               *head = 0L;  //  The last thing loaded
  u32bit                 batchSize  = 0;
  u32bit                 batchLimit = 128;

  bool  moreToLoad = true;

  while (moreToLoad) {

    //  Zzzzzzz....
    while (_numberLoaded - _numberComputed > _loaderQueueSize)
      nanosleep(&naptime, 0L);

    state_s  *thisState = new state_s( (*_userLoader)(_globalUserData) );

    //  If we actually loaded a new state, add it
    //
    if (thisState->_user) {
      if (_loadBatches) {
        loaderSave(tail, head, thisState);
        batchSize++;
        if (batchSize >= batchLimit)
          loaderAppend(tail, head);
      } else {
        loaderAdd(thisState);
      }
    } else {
      //  Didn't read, must be all done!  Push on the end-of-input marker state.
      //
      if (_loadBatches) {
        loaderSave(tail, head, new state_s(0L));
        loaderAppend(tail, head);
      } else {
        loaderAdd(new state_s(0L));
      }
      moreToLoad = false;
      delete thisState;
    }
  }

  //fprintf(stderr, "sweatShop::reader exits.\n");
  return(0L);
}



void*
sweatShop::worker(sweatShopWorker *workerData) {

  //fprintf(stderr, "Worker %d! (userData=%p)\n", workerData->threadID, workerData->threadUserData);

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 500000000ULL;  //  1/2 second

  bool    moreToCompute = true;

  while (moreToCompute) {

    while (_numberComputed - _numberOutput > _writerQueueSize)
      nanosleep(&naptime, 0L);

    //  Grab the next state.  We don't grab it if it's the last in the
    //  queue (else we would fall off the end) UNLESS it really is the
    //  last one.
    //
    state_s  *thisState = 0L;

    pthread_mutex_lock(&_stateMutex);

    if ((_workerP) && ((_workerP->_next != 0L) || (_workerP->_user == 0L))) {
      thisState = _workerP;
      _workerP = _workerP->_next;
    }

    if (_workerP == 0L)
      moreToCompute = false;

    pthread_mutex_unlock(&_stateMutex);

    //  Execute
    //
    if (thisState && thisState->_user) {
      (*_userWorker)(_globalUserData, workerData->threadUserData, thisState->_user);
      thisState->_computed = true;
      workerData->numComputed++;
    } else {
      //  When we really do run out of stuff to do, we'll end up here
      //  (only one thread will end up in the other case, with
      //  something to do and moreToCompute=false).  If it's actually
      //  the end, skip the sleep and just get outta here.
      //
      if (moreToCompute == true)
        nanosleep(&naptime, 0L);
    }
  }

  //fprintf(stderr, "sweatShop::worker exits.\n");
  return(0L);
}


void*
sweatShop::writer(void) {
  state_s  *deleteState = 0L;

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 333333333ULL;  //  1/3 second 10000000ULL;

  //  Wait for output to appear.
  //
  while (_writerP && _writerP->_user) {

    //  SOMEONE needs to update the number of things computed,
    //  otherwise everyone stalls.  This used to be in stats(), but if
    //  we disable those, then obviously we don't update this.
    //
    u32bit nc = 0;
    for (u32bit i=0; i<_numberOfWorkers; i++)
      nc += _workerData[i].numComputed;
    _numberComputed = nc;

    //  Wait a bit if there is no output on this node (we caught up to
    //  the computes).  We now don't need to explicitly check if we catch
    //  the input (yikes!) because it still won't have been computed.
    //
    if ((_writerP->_computed == false) || (_writerP->_next == 0L)) {
      nanosleep(&naptime, 0L);
    } else {
      (*_userWriter)(_globalUserData, _writerP->_user);
      _numberOutput++;

      deleteState = _writerP;
      _writerP    = _writerP->_next;
      delete deleteState;
    }
  }

  //fprintf(stderr, "sweatShop::writer exits.\n");
  return(0L);
}


void*
sweatShop::status(void) {

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 250000000ULL;

  double  startTime = getTime() - 0.001;

  u64bit  deltaOut = 0;
  u64bit  deltaCPU = 0;
  double  perSec   = 0;

  while (_writerP && _writerP->_user) {
    deltaOut = deltaCPU = 0;

    if (_numberComputed > _numberOutput)
      deltaOut = _numberComputed - _numberOutput;
    if (_numberLoaded > _numberComputed)
      deltaCPU = _numberLoaded - _numberComputed;

    perSec = _numberOutput / (getTime() - startTime);

    fprintf(stderr, "%6.1f/s (out="u64bitFMTW(8)") + "u64bitFMTW(8)" = (cpu = "u64bitFMTW(8)") + "u64bitFMTW(8)" = (in = "u64bitFMTW(8)")\r",
            perSec,
            _numberOutput,
            deltaOut,
            _numberComputed,
            deltaCPU,
            _numberLoaded);
    fflush(stderr);
    nanosleep(&naptime, 0L);

    //  Too big?  Shrink.
    if (10 * perSec < _loaderQueueSize)
      _loaderQueueSize = (u32bit)(5.0 * perSec);

    //  Too small?  Grow.
    if (_loaderQueueSize < 4.0 * perSec)
      _loaderQueueSize = (u32bit)(5.0 * perSec);
  }

  if (_numberComputed > _numberOutput)
    deltaOut = _numberComputed - _numberOutput;
  if (_numberLoaded > _numberComputed)
    deltaCPU = _numberLoaded - _numberComputed;

  fprintf(stderr, "%6.1f/s (out="u64bitFMTW(8)") + "u64bitFMTW(8)" = (cpu = "u64bitFMTW(8)") + "u64bitFMTW(8)" = (in = "u64bitFMTW(8)")\r",
          _numberOutput / (getTime() - startTime),
          _numberOutput,
          _numberComputed - _numberOutput,
          _numberComputed,
          _numberLoaded - _numberComputed,
          _numberLoaded);

  //fprintf(stderr, "sweatShop::status exits.\n");
  return(0L);
}





void
sweatShop::run(void *user, bool beVerbose) {
  pthread_attr_t   threadAttr;
  pthread_t        threadID;

  _globalUserData = user;

  pthread_mutex_init(&_stateMutex, NULL);

  pthread_attr_init(&threadAttr);
  pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_DETACHED);
  pthread_attr_setschedpolicy(&threadAttr, SCHED_OTHER);

  //  Fire off the loader, or the all-vs-all loader
  //
  pthread_create(&threadID, &threadAttr, _sweatshop_loaderThread, this);

  //  Wait for it to actually load something (otherwise all the
  //  workers immediately go home)
  //
  while (!_writerP && !_workerP && !_loaderP) {
    struct timespec   naptime;
    naptime.tv_sec      = 0;
    naptime.tv_nsec     = 250000ULL;
    nanosleep(&naptime, 0L);
  }


  //  Fire off some workers
  //
  for (u32bit i=0; i<_numberOfWorkers; i++) {
    _workerData[i].shop           = this;
    _workerData[i].threadID       = i;

    pthread_create(&threadID, &threadAttr, _sweatshop_workerThread, _workerData + i);
  }

  //  And the stats
  //
  if (beVerbose)
    pthread_create(&threadID, &threadAttr, _sweatshop_statusThread, this);

  //  We run the writer in the main, right here.  We could probably
  //  run stats here, but we certainly aren't done until the worker
  //  finishes, so just do that one.
  //
  writer();

  //  If the writer exits, then someone signalled that we're done.
  //  Clean up -- make sure all the queues are free'd.
  //
  delete _loaderP;
  _loaderP = _workerP = _writerP = 0L;
}

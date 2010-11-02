#ifndef SWEATSHOP_H
#define SWEATSHOP_H

#include <pthread.h>
#include <semaphore.h>

#include "util++.H"

class sweatShopWorker;
class sweatShopState;

class sweatShop {
public:
  sweatShop(void*(*loaderfcn)(void *G),
            void (*workerfcn)(void *G, void *T, void *S),
            void (*writerfcn)(void *G, void *S));
  ~sweatShop();

  void        setNumberOfWorkers(u32bit x) { _numberOfWorkers = x; };

  void        setThreadData(u32bit t, void *x);

  void        setLoaderBatchSize(u32bit batchSize) { _loaderBatchSize = batchSize; };
  void        setLoaderQueueSize(u32bit queueSize) { _loaderQueueSize = queueSize;  _loaderQueueMax = queueSize; };

  void        setWorkerBatchSize(u32bit batchSize) { _workerBatchSize = batchSize; };

  void        setWriterQueueSize(u32bit queueSize) { _writerQueueSize = queueSize;  _writerQueueMax = queueSize; };

  void        run(void *user=0L, bool beVerbose=false);
private:

  //  Stubs that forward control from the c-based pthread to this class
  friend void  *_sweatshop_loaderThread(void *ss);
  friend void  *_sweatshop_workerThread(void *ss);
  friend void  *_sweatshop_writerThread(void *ss);
  friend void  *_sweatshop_statusThread(void *ss);

  //  The threaded routines
  void   *loader(void);
  void   *worker(sweatShopWorker *workerData);
  void   *writer(void);
  void   *status(void);

  //  Utilities for the loader thread
  //void    loaderAdd(sweatShopState *thisState);
  void    loaderSave(sweatShopState *&tail, sweatShopState *&head, sweatShopState *thisState);
  void    loaderAppend(sweatShopState *&tail, sweatShopState *&head);

  pthread_mutex_t        _stateMutex;

  void                *(*_userLoader)(void *global);
  void                 (*_userWorker)(void *global, void *thread, void *thing);
  void                 (*_userWriter)(void *global, void *thing);

  void                  *_globalUserData;

  sweatShopState        *_writerP;  //  Where output takes stuff from, the tail
  sweatShopState        *_workerP;  //  Where computes happen, the middle
  sweatShopState        *_loaderP;  //  Where input is put, the head

  bool                   _showStatus;

  u32bit                 _loaderQueueSize, _loaderQueueMax;
  u32bit                 _loaderBatchSize;
  u32bit                 _workerBatchSize;
  u32bit                 _writerQueueSize, _writerQueueMax;

  u32bit                 _numberOfWorkers;

  sweatShopWorker       *_workerData;

  u64bit                 _numberLoaded;
  u64bit                 _numberComputed;
  u64bit                 _numberOutput;
};

#endif  //  SWEATSHOP_H

#ifndef SWEATSHOP_H
#define SWEATSHOP_H

#if defined(__bsd) || defined(__alpha)

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199506L
#endif

#endif

#include <pthread.h>
#include <semaphore.h>

//  For u32bit, getTime() -- really only used by the status thread
#include "util++.H"

class sweatShop;

struct sweatShopWorker {
  sweatShop   *shop;
  void        *threadUserData;
  u32bit       threadID;
  u32bit       numComputed;
};

class sweatShop {
public:
  sweatShop(void*(*loader)(void *G),
            void (*worker)(void *G, void *T, void *S),
            void (*writer)(void *G, void *S));
  ~sweatShop();

  void        setLoader(void*(*f)(void *G))                   { _userLoader = f; };
  void        setWorker(void (*f)(void *G, void *T, void *S)) { _userWorker = f; };
  void        setWriter(void (*f)(void *G, void *S))          { _userWriter = f; };

  u32bit      numberOfWorkers(u32bit x=0);

  void        setThreadData(u32bit t, void *x);

  u32bit      loaderQueueSize(u32bit x=0);
  u32bit      writerQueueSize(u32bit x=0);

  void        run(void *user=0L, bool beVerbose=false);
private:

  //  This gets created by the loader, passed to the worker, and printed
  //  by the writer.  userData is controlled by the user.
  class state_s {
  public:
    state_s(void *userData) {
      _user     = userData;
      _computed = false;
      _next     = 0L;
    };
    ~state_s() {
    };

    void        *_user;
    bool         _computed;
    state_s     *_next;
  };

  //  Stubs that forward control from the c-based pthread to this class
  friend void  *loaderThread(void *ss);
  friend void  *workerThread(void *ss);
  friend void  *writerThread(void *ss);
  friend void  *statusThread(void *ss);

  //  The threaded routines
  void   *loader(void);
  void   *worker(sweatShopWorker *workerData);
  void   *writer(void);
  void   *status(void);

  //  Utilities for the loader thread
  void    loaderAdd(state_s *thisState);
  void    loaderSave(state_s *&tail, state_s *&head, state_s *thisState);
  void    loaderAppend(state_s *&tail, state_s *&head);

  pthread_mutex_t        _stateMutex;

  void                *(*_userLoader)(void *global);
  void                 (*_userWorker)(void *global, void *thread, void *thing);
  void                 (*_userWriter)(void *global, void *thing);

  void                  *_globalUserData;

  state_s               *_writerP;  //  Where output takes stuff from, the tail
  state_s               *_workerP;  //  Where computes happen, the middle
  state_s               *_loaderP;  //  Where input is put, the head

  bool                   _loadBatches;

  u32bit                 _loaderQueueSize;
  u32bit                 _writerQueueSize;

  u32bit                 _numberOfWorkers;

  sweatShopWorker       *_workerData;

  u64bit                 _numberLoaded;
  u64bit                 _numberComputed;
  u64bit                 _numberOutput;
};

#endif  //  SWEATSHOP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "meryl.H"
#include "merstreamfile.H"
#include "libmeryl.H"
#include "britime.H"





#ifndef ENABLE_THREADS

void
runThreaded(merylArgs *args) {
}

#else

//  Rest of this file is for threading support.
#include <pthread.h>



//  launch t threads
//    get the segmentNext (protected by a mutex)
//    compute
//    touch segmentFinished[s]
//    repeat or exit if segmentNext >= segmentMax
//  main waits for all segmentFinished[] to be touched
//  close threads, cleanup, return


void
runSegment(merylArgs *args, u64bit segment);



//  Shared data
//
pthread_mutex_t        segmentMutex;
u64bit                 segmentNext;
u64bit                 segmentMax;
u32bit                *segmentDone;




void*
buildThread(void *U) {
  u64bit      segment = u32bitZERO;
  merylArgs  *args = (merylArgs *)U;

  while (segment < segmentMax) {
    pthread_mutex_lock(&segmentMutex);
    segment = segmentNext++;
    pthread_mutex_unlock(&segmentMutex);

    if (segment < segmentMax) {
      runSegment(args, segment);
      segmentDone[segment]++;
    }
  }

  fprintf(stderr, "Thread exits.\n");

  return(0L);
}





void
runThreaded(merylArgs *args) {

  //  Clear stuff
  //
  segmentNext = u64bitZERO;
  segmentMax  = args->segmentLimit;
  segmentDone = new u32bit [segmentMax];
  for (u64bit s=0; s<segmentMax; s++)
    segmentDone[s] = u32bitZERO;


  //  Initialize threads
  //
  pthread_attr_t   threadAttr;
  pthread_t        threadID;

  pthread_mutex_init(&segmentMutex, NULL);

  pthread_attr_init(&threadAttr);
  pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_DETACHED);
  pthread_attr_setschedpolicy(&threadAttr, SCHED_OTHER);


  //  Start the threads
  //
  for (u64bit i=0; i<args->numThreads; i++)
    pthread_create(&threadID, &threadAttr, buildThread, (void *)args);


  //  Wait for the threads to complete
  //
  struct timespec  shortNap;
  shortNap.tv_sec  = 1;
  shortNap.tv_nsec = 0;

  u64bit s=0;
  while (s < segmentMax) {
    if (segmentDone[s] == 0)
      nanosleep(&shortNap, 0L);
    else
      s++;
  }


  fprintf(stderr, "Threads all done, cleaning up.\n");


  //  Cleanup
  //
  pthread_attr_destroy(&threadAttr);
  pthread_mutex_destroy(&segmentMutex);

  delete [] segmentDone;
}



#endif  //  ENABLE_THREADS

#include <stdio.h>
#include <stdlib.h>
#include <new.h>

#include "posix.H"
#include "snapper2.H"

//  Define this to print a message whenever a sequence is loaded.
//  Useful for testing the loader with large sequences (scaffolds,
//  chromosomes).
//
//#define VERBOSE_LOADER

void*
loaderThread(void *) {
  u32bit               waterLevel = 0;
  FastASequenceInCore *B          = 0L;
  bool                 slept      = false;

  while (inputHead < numberOfQueries) {

    //  We fill the input as fast as we can, up to the high water
    //  mark, then we take a little snooze to let the workers catch up.
    //
    pthread_mutex_lock(&inputTailMutex);
    waterLevel = inputHead - inputTail;
    pthread_mutex_unlock(&inputTailMutex);

    //  Warn if we're too small.
    //
    if ((slept) && (waterLevel <= 1)) {

      u32bit i = 0.1 * config._loaderHighWaterMark;
      if (i == 0)
        i = 1;
      config._loaderHighWaterMark += i;

      config.setTime(&config._loaderSleep,
                     0.9 * ((double)config._loaderSleep.tv_sec + (double)config._loaderSleep.tv_nsec * 1e-9));

      if (config._loaderWarnings)
        fprintf(stderr,
#ifdef TRUE64BIT
                "WARNING: Loader ran dry.  Increasing limit to %u sequences, decreasing sleep to %f.\n",
#else
                "WARNING: Loader ran dry.  Increasing limit to %lu sequences, decreasing sleep to %f.\n",
#endif
                config._loaderHighWaterMark,
                ((double)config._loaderSleep.tv_sec + (double)config._loaderSleep.tv_nsec * 1e-9));
    }

    //  Sleep, if we need to, otherwise, get the next sequence and
    //  push it onto the input list at inputHead.  This alloc is
    //  deleted by the output thread.
    //
    if (waterLevel >= config._loaderHighWaterMark) {
      slept = true;
      nanosleep(&config._loaderSleep, 0L);
    } else {
      slept = false;

#ifdef VERBOSE_LOADER
      fprintf(stderr, "Loading sequence %u (tail = %u)\n", inputHead, inputTail);
#endif

      try {
        B = qsFASTA->getSequence();
      } catch (std::bad_alloc) {
        fprintf(stderr, "loaderThread()-- Failed to load next query sequence\ncaught bad_alloc in %s at line %d\n", __FILE__, __LINE__);
        exit(1);
      }

      pthread_mutex_lock(&inputTailMutex);
      input[inputHead] = B;
      inputHead++;
      pthread_mutex_unlock(&inputTailMutex);
    }
  }

  return(0L);
}

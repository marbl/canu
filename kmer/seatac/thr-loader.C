#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "posix.H"
#include "seatac.H"

//  Define this to print a message whenever a sequence is loaded.
//  Useful for testing the loader with large sequences (scaffolds,
//  chromosomes).
//
//#define VERBOSE_LOADER

#ifdef TRUE64BIT
char const *loadDesc = "WARNING: Loader ran dry.  Increasing limit to %u sequences, decreasing sleep to %f.\n";
#else
char const *loadDesc = "WARNING: Loader ran dry.  Increasing limit to %lu sequences, decreasing sleep to %f.\n";
#endif

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

      u32bit i = (u32bit) (0.1 * config._loaderHighWaterMark);
      if (i == 0)
        i = 1;
      config._loaderHighWaterMark += i;

      config.setTime(&config._loaderSleep,
                     0.9 * ((double)config._loaderSleep.tv_sec + (double)config._loaderSleep.tv_nsec * 1e-9));

      if (config._loaderWarnings)
        fprintf(stderr, loadDesc,
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

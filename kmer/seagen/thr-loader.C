#include <stdio.h>
#include <stdlib.h>
#include <new.h>

#include "posix.H"
#include "searchGENOME.H"

//  Define this to get low water level warnings
//
#define WATERLEVELWARNINGS

#ifdef WATERLEVELWARNINGS
#ifdef TRUE64BIT
char const *loadDesc = "WARNING: Loader is at %3.0f%% capacity. (cDNA %7u out of %7u).  Load/Search might be unbalanced.\n";
#else
char const *loadDesc = "WARNING: Loader is at %3.0f%% capacity. (cDNA %7lu out of %7lu).  Load/Search might be unbalanced.\n";
#endif
#endif

//  XXX:  We probably should make this produce FastABuffer references,
//  instead of copying the sequence.  Would want to do it so that the
//  buffers are reused, which makes it a little hard.

void*
loaderThread(void *) {

  try {
    u32bit               LOADER_HIGH_WATER_MARK = 16384;
    u32bit               waterLevel             = 0;
    FastASequenceInCore *B;

    struct timespec loaderSleep = { 0, 25000000 };

    while (inputHead < numberOfQueries) {

      //  We fill the input as fast as we can, up to the high water
      //  mark, then we take a little snooze to let the workers catch up.
      //
      pthread_mutex_lock(&inputTailMutex);
      waterLevel = inputHead - inputTail;
      pthread_mutex_unlock(&inputTailMutex);

      //  Now, sleep, if we need to.
      //
      if (waterLevel >= LOADER_HIGH_WATER_MARK) {
        nanosleep(&loaderSleep, 0L);

#ifdef WATERLEVELWARNINGS
        //  Warn if we're too small.
        //
        if (waterLevel < (LOADER_HIGH_WATER_MARK >> 1))
          fprintf(stderr, loadDesc,
                  100.0 * waterLevel / LOADER_HIGH_WATER_MARK,
                  inputHead, numberOfQueries);
#endif
      } else {

        //  Get the next cDNA and push it onto the input list at
        //  inputHead.  This alloc is deleted by the output thread.
        //
        B = qsFASTA->getSequence();

        pthread_mutex_lock(&inputTailMutex);
        input[inputHead] = B;
        inputHead++;
        pthread_mutex_unlock(&inputTailMutex);
      }
    }
  } catch (std::bad_alloc) {
    fprintf(stderr, "loaderThread()-- caught bad_alloc in %s at line %d\n", __FILE__, __LINE__);
  }

  return(0L);
}

#include "thr.H"





void*
searchThread(void *U) {
  u32bit               idx      = 0;
  FastASequenceInCore *seq      = 0L;
  u32bit               blockedI = 0;
  u32bit               blockedO = 0;

  searcherState       *state    = new searcherState((u64bit)U);

  fprintf(stderr, "Hello!  I'm a searchThread (number "u64bitFMT")!\n", state->threadID);

  //  Allocate and fill out the thread stats -- this ensures that we
  //  always have stats (even if they're bogus).
  //
  threadStats[state->threadID] = new char [1025];
  sprintf(threadStats[state->threadID], "No status.\n");

  while (inputTail < numberOfQueries) {

    //  Grab the next sequence.
    //
    pthread_mutex_lock(&inputTailMutex);
    idx = inputTail;
    if (idx < numberOfQueries) {
      seq = input[idx];
      input[idx] = 0L;
      if (seq)
        inputTail++;
    }
    pthread_mutex_unlock(&inputTailMutex);

    //  Still need to check that the index is valid.  Another thread
    //  could (and does) steal execution between the while and the
    //  mutex lock.
    //
    if (idx < numberOfQueries) {

      //  If there is no sequence, oh boy, we are in bad shape.  Sleep a
      //  little bit to let the loader catch up, then try again.
      //
      if (seq == 0L) {
        //if (config._loaderWarnings)
        //  fprintf(stderr, "Thread "u64bitFMT" Blocked by input.\n", state->threadID);
        blockedI++;
        nanosleep(&config._searchSleep, 0L);
      } else {

        //  If our idx is too far away from the output thread, sleep
        //  a little bit.  We keep the idx and seq that we have obtained,
        //  though.
        //
        while (idx > (outputPos + config._writerHighWaterMark)) {
          if (config._writerWarnings)
            fprintf(stderr, "Thread "u64bitFMT" blocked by output on idx = "u32bitFMT", outputPos = "u32bitFMT").\n", state->threadID, idx, outputPos);
          blockedO++;
          nanosleep(&config._searchSleep, 0L);
        }

        u32bit               theHitsLen = 0;
        u32bit               theHitsMax = 4;
        aHit                *theHits    = new aHit [theHitsMax];

        logMsg              *theLog     = new logMsg;

#ifdef VERBOSE_SEARCH
        double startTime = getTime();
#endif

        ///////////////////////////////////////
        //
        //  Do searches.
        //
        if (config._doForward)
          doSearch(state, seq, idx, false, theHits, theHitsLen, theHitsMax, theLog);
        if (config._doReverse)
          doSearch(state, seq, idx, true,  theHits, theHitsLen, theHitsMax, theLog);

#ifdef VERBOSE_SEARCH
        double searchTime = getTime();
#endif

        ///////////////////////////////////////
        //
        //  Filter the hits
        //
        u32bit numF = doFilter(state, theHits, theHitsLen, theLog);

#ifdef VERBOSE_SEARCH
        double filterTime = getTime();
#endif

        ///////////////////////////////////////
        //
        //  Polish the filtered hits
        //
        char *o = doPolish(state, seq, theHits, theHitsLen, theLog);


        //  Clean up
        //
        for (u32bit h=0; h<theHitsLen; h++) {
          delete theHits[h]._ML;
          theHits[h]._ML = 0L;
        }


        ///////////////////////////////////////
        //
        //  Signal that we are done.
        //
        answerLen[idx] = theHitsLen;
        answer[idx]    = theHits;
        logmsg[idx]    = theLog;
        outputLen[idx] = (u32bit)strlen(o);
        output[idx]    = o;

        state->searched++;

        //  Clean up stuff
        //
        delete    seq;

#ifdef VERBOSE_SEARCH
        double endTime = getTime();
        if ((endTime - startTime) > 5.0) {
          fprintf(stdout, u64bitFMT" done with idx = "u32bitFMT" ("u32bitFMT" bp, "u32bitFMT" hits) search %f filter %f polish %f  search %f chain %f polish %f\n",
                  state->threadID, idx, seq->sequenceLength(), theHitsLen,
                  searchTime - startTime,
                  filterTime - searchTime,
                  endTime - filterTime,
                  state->searchTime,
                  state->chainTime,
                  state->polishTime);
          fflush(stdout);
        }
#endif




#ifdef MEMORY_DEBUG
        //fprintf(stderr, "Memory dump in searchThread\n");
        //_dump_allocated_delta(fileno(stderr));
#endif

      } // end of seq != 0L
    } // end of idx < numberOfQueries
  } // end of inputTail < numberOfQueries


  //  OK, now fill out the read thread stats
  //
  sprintf(threadStats[state->threadID],
          "["u64bitFMT"] searched: "u64bitFMTW(8)" polished: "u64bitFMTW(8)"  blockedI/O: "u32bitFMTW(4)"/"u32bitFMTW(4)" Time: search:%8.2f chain:%8.2f polish:%8.2f\n",
          state->threadID,
          state->searched, state->polished,
          blockedI, blockedO,
          state->searchTime,
          state->chainTime,
          state->polishTime);

  delete state;

  return(0L);
}

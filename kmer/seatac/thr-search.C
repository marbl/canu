#include "posix.H"
#include "seatac.H"
#include "time.H"


#ifdef TRUE64BIT
char const *srchGbye = "[%lu] computed: %8lu  blocked: %4lu/%4lu  encodeTime: %7.2f   searchTime: %7.2f   processTime: %7.2f\n";
#else
char const *srchGbye = "[%llu] computed: %8lu  blocked: %4lu/%4lu  encodeTime: %7.2f   searchTime: %7.2f   processTime: %7.2f\n";
#endif


class searcherState {
public:
  u64bit         posnMax;
  u64bit         posnLen;
  u64bit        *posn;

  double         encodeTime;
  double         maskTime;
  double         searchTime;
  double         processTime;

  searcherState() {
    posnMax = 16384;
    posnLen = 0;
    posn    = new u64bit [ posnMax ];

    encodeTime = 0.0;
    maskTime   = 0.0;
    searchTime = 0.0;
    processTime = 0.0;
  };

  ~searcherState() {
    delete [] posn;
  };
};


void
doSearch(searcherState *state,
         FastASequenceInCore *seq,
         u32bit idx,
         bool rc,
         filterObj *FO) {
  encodedQuery  *query  = 0L;
  hitMatrix     *matrix = 0L;
  double         startTime  = 0.0;

  //  Build and mask the query
  //
  startTime = getTime();
  query = new encodedQuery(seq->sequence(), seq->sequenceLength(), config._merSize, rc);
  state->encodeTime += getTime() - startTime;

  startTime = getTime();

  state->maskTime += getTime() - startTime;

  //  Get the hits
  //
  startTime = getTime();
  matrix = new hitMatrix(seq->sequenceLength(), idx);

  u64bit mer = u64bitZERO;
  u32bit pos = u32bitZERO;

  fprintf(stderr, "Retrieving hits for rc=%d %s\n", rc, seq->header());

  while (query->getMer(mer, pos) == true)
    if (positions->get(mer, state->posn, state->posnMax, state->posnLen))
      matrix->addHits(pos, state->posn, state->posnLen);

  state->searchTime += getTime() - startTime;

  fprintf(stderr, "Filtering hits for  rc=%d %s\n", rc, seq->header());

  startTime = getTime();
  matrix->processMatrix(rc ? 'r' : 'f', FO);
  state->processTime += getTime() - startTime;

  fprintf(stderr, "Finished with       rc=%d %s\n", rc, seq->header());

  delete matrix;
  delete query;
}



void*
searchThread(void *U) {
  u32bit               idx      = 0;
  FastASequenceInCore *seq      = 0L;
  u32bit               blockedI = 0;
  u32bit               blockedO = 0;
  u32bit               computed = 0;

  searcherState       *state    = new searcherState;

  //  Allocate and fill out the thread stats -- this ensures that we
  //  always have stats (even if they're bogus).
  //
  threadStats[(u64bit)U] = new char [1025];
  sprintf(threadStats[(u64bit)U], srchGbye,
          (u64bit)U,
          (u32bit)0, (u32bit)0, (u32bit)0,
          0.0, 0.0, 0.0);

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
        //  fprintf(stderr, "%lu Blocked by input.\n", (u64bit)U);
        blockedI++;
        nanosleep(&config._searchSleep, 0L);
      } else {

        //  If our idx is too far away from the output thread, sleep
        //  a little bit.  We keep the idx and seq that we have obtained,
        //  though.
        //
        while (idx > (outputPos + config._writerHighWaterMark)) {
          if (config._writerWarnings)
            fprintf(stderr, "%lu Blocked by output (idx = %d, outputPos = %d).\n", (u64bit)U, idx, outputPos);
          blockedO++;
          nanosleep(&config._searchSleep, 0L);
        }

        //  Construct a filter object
        //
        filterObj *FO = new filterObj(config._filtername, config._filteropts);
        fprintf(stderr, "Created filterObj 0x%016lx\n", FO);

        //  Do searches.
        //
        if (config._doForward)
          doSearch(state, seq, idx, false, FO);
        if (config._doReverse)
          doSearch(state, seq, idx, true, FO);

        // CMM
        FO->filter();
        //  Signal that we are done.
        //
        output[idx] = FO;
        computed++;

        delete seq;

#if 0
        fprintf(stderr, "I sleep\n");
        for (int i=0; i<1000; i++) {
          nanosleep(&config._searchSleep, 0L);
        }
        fprintf(stderr, "I awake\n");
#endif

      } // end of seq != 0L
    } // end of idx < numberOfQueries
  } // end of inputTail < numberOfQueries


  //  OK, now fill out the read thread stats
  //
  sprintf(threadStats[(u64bit)U], srchGbye, (u64bit)U,
          computed, blockedI, blockedO,
          state->encodeTime,
          state->searchTime,
          state->processTime);

  delete state;

  return(0L);
}

#include "posix.H"
#include "searchGENOME.H"
#include "time.H"

//  If you really, really, really want to know the exact number
//  of bases left in the query, use the interval list.  Otherwise,
//  it's faster to guess.
//
//#define USEEXACTSIZE


#ifdef TRUE64BIT
char const *srchGbye = "[%lu] computed: %8lu  blocked: %4lu/%4lu  encodeTime: %7.2f   searchTime: %7.2f   filterTime: %7.2f\n";
#else
char const *srchGbye = "[%llu] computed: %8lu  blocked: %4lu/%4lu  encodeTime: %7.2f   searchTime: %7.2f   filterTime: %7.2f\n";
#endif


class searcherState {
public:
  u64bit         posnMax;
  u64bit         posnLen;
  u64bit        *posn;

  double         encodeTime;
  double         maskTime;
  double         searchTime;
  double         filterTime;

  searcherState() {
    posnMax = 16384;
    posnLen = 0;
    posn    = new u64bit [ posnMax ];

    encodeTime = 0.0;
    maskTime   = 0.0;
    searchTime = 0.0;
    filterTime = 0.0;
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
         char *&theOutput, u32bit &theOutputPos, u32bit &theOutputMax) {
  encodedQuery  *query  = 0L;
  hitMatrix     *matrix = 0L;
  u32bit         qMers  = 0;
  double         startTime  = 0.0;

  //  Build and mask the query
  //
  startTime = getTime();
  query = new encodedQuery(seq->sequence(), seq->sequenceLength(), config._merSize, rc);
  state->encodeTime += getTime() - startTime;


  startTime = getTime();
  if (maskDB)
    for (u32bit qi=0; qi<query->numberOfMers(); qi++)
      if ((query->getSkip(qi) == false) &&
          (maskDB->exists(query->getMer(qi))))
        query->setSkip(qi);

  if (onlyDB)
    for (u32bit qi=0; qi<query->numberOfMers(); qi++)
      if ((query->getSkip(qi) == false) &&
          (!onlyDB->exists(query->getMer(qi))))
        query->setSkip(qi);

#ifdef USEEXACTSIZE
  intervalList   *IL = new intervalList(config._merSize);

  for (u32bit qi=0; qi<query->numberOfMers(); qi++) {
    if (query->getSkip(qi) == false)
      IL->addInterval(qi);
  }

  qMers = IL->sumIntervalLengths();
  delete IL;
#else
  qMers = query->numberOfValidMers();
#endif
  state->maskTime += getTime() - startTime;



  //  Get the hits
  //
  startTime = getTime();
  matrix = new hitMatrix(seq->sequenceLength(), qMers, idx);
  for (u32bit qi=0; qi<query->numberOfMers(); qi++)
    if ((query->getSkip(qi) == false) &&
        (positions->get(query->getMer(qi), state->posn, state->posnMax, state->posnLen)))
      matrix->addHits(qi, state->posn, state->posnLen);
  state->searchTime += getTime() - startTime;

  //  Filter, storing the resutls into theOutput
  //
  startTime = getTime();
  matrix->filter(rc ? 'r' : 'f', theOutput, theOutputPos, theOutputMax);
  state->filterTime += getTime() - startTime;

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

  u32bit               theOutputPos = 0;
  u32bit               theOutputMax = 0;
  char                *theOutput    = 0L;

  struct timespec      searchSleepLong  = { 0, 500000000 };
  struct timespec      searchSleepShort = { 0,  10000000 };

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
        //fprintf(stderr, "%lu Blocked by input.\n", (u64bit)U);
        blockedI++;
        nanosleep(&searchSleepLong, 0L);
      } else {

        //  If our idx is too far away from the output thread, sleep
        //  a little bit.  We keep the idx and seq that we have obtained,
        //  though.
        //
        while (idx > (outputPos + 32 * 1024)) {
          //fprintf(stderr, "%lu Blocked by output (idx = %d, outputPos = %d).\n", (u64bit)U, idx, outputPos);
          blockedO++;
          nanosleep(&searchSleepShort, 0L);
        }

        //  Allocate space for the output -- 1MB should be enough for
        //  about 29000 signals.  Make it 32K -> 900 signals.
        //
        theOutputPos = 0;
        theOutputMax = 32 * 1024;
        theOutput    = new char [theOutputMax];

        //  Do searches.
        //
        if (config._doForward)
          doSearch(state, seq, idx, false, theOutput, theOutputPos, theOutputMax);
        if (config._doReverse)
          doSearch(state, seq, idx, true,  theOutput, theOutputPos, theOutputMax);

        //  Signal that we are done.
        //
        outputLen[idx] = theOutputPos;
        output[idx]    = theOutput;
        computed++;

        delete seq;
      } // end of seq != 0L
    } // end of idx < numberOfQueries
  } // end of inputTail < numberOfQueries


  //  OK, now fill out the read thread stats
  //
  sprintf(threadStats[(u64bit)U], srchGbye, (u64bit)U,
          computed, blockedI, blockedO,
          state->encodeTime,
          state->searchTime,
          state->filterTime);

  delete state;

  return(0L);
}

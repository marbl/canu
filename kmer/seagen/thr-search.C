#include "posix.H"
#include "searchGENOME.H"

//  $Id$

#define VERBOSE
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
         unsigned char *seq, u32bit seqLen, u32bit idx,
         bool rc,
         char *&theOutput, u32bit &theOutputPos, u32bit &theOutputMax) {
  encodedQuery  *query  = 0L;
  hitMatrix     *matrix = 0L;
  u32bit         qMers  = 0;
  double         startTime  = 0.0;

  //  Build and mask the query
  //
  startTime = getTime();
  query = new encodedQuery(seq,
                           seqLen,
                           config._merSize,
                           rc);
  state->encodeTime += getTime() - startTime;

  qMers = query->numberOfMers();

  if (maskDB) {
    //fprintf(stderr, "Begin masking qMers=%u\n", qMers);

    //  If you really, really, really want to know the exact number
    //  of bases left in the query, use the interval list.  Otherwise,
    //  it's faster to guess.
    //
    //intervalList   *IL = new intervalList(config._merSize);

    startTime = getTime();
    for (u32bit qi=0; qi<query->numberOfMers(); qi++) {
      if ((query->getSkip(qi) == false) &&
          (maskDB->exists(query->getMer(qi)))) {
        qMers--;
        query->setSkip(qi);
      }

      //if (query->getSkip(qi) == false)
      //  IL->addInterval(qi);
    }

    //qMers = IL->sumIntervalLengths();
    //delete IL;

    state->maskTime += getTime() - startTime;
  }

  if (onlyDB) {
    startTime = getTime();
    for (u32bit qi=0; qi<query->numberOfMers(); qi++) {
      if ((query->getSkip(qi) == false) &&
          (!onlyDB->exists(query->getMer(qi)))) {
        qMers--;
        query->setSkip(qi);
      }
    }

    state->maskTime += getTime() - startTime;
  }

  //  Construct a new hitMatrix.  There is only one of these for all sequences in the
  //  positionDB.
  //
  matrix = new hitMatrix(seqLen, qMers, idx);

  //  Get the hits
  //
  startTime = getTime();

  for (u32bit qi=0; qi<query->numberOfMers(); qi++) {
    if ((query->getSkip(qi) == false) &&
        (positions->get(query->getMer(qi), state->posn, state->posnMax, state->posnLen))) {
#if 0
      if (state->posnLen != 1) {
        fprintf(stderr, "ERROR: Got %lu hits for mer qi=%u\n", state->posnLen, qi);
        for (int x=0; x<state->posnLen; x++)
          fprintf(stderr, "  %lu\n", state->posn[x]);
      }
#endif

      matrix->addHits(qi, state->posn, state->posnLen);
    }
  }

  state->searchTime += getTime() - startTime;

  //  Filter, storing the resutls
  //
  startTime = getTime();
  matrix->filter(rc ? 'r' : 'f', theOutput, theOutputPos, theOutputMax);
  state->filterTime += getTime() - startTime;

  delete matrix;
  delete query;
}



void*
searchThread(void *U) {
  u32bit         idx      = 0;
  unsigned char *seq      = 0L;
  u32bit         seqLen   = 0;
  u32bit         blockedI = 0;
  u32bit         blockedO = 0;
  u32bit         computed = 0;

  //fprintf(stderr, "Hello!  I'm a searcher!\n");

  searcherState *state    = new searcherState;

  u32bit          theOutputPos = 0;
  u32bit          theOutputMax = 0;
  char           *theOutput    = 0L;

  struct timespec searchSleepLong  = { 0, 500000000 };
  struct timespec searchSleepShort = { 0,  10000000 };

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
      seq    = input[idx];
      seqLen = inputLen[idx];
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

#if 0
        if (config._beVerbose)
          fprintf(stderr, "%lu Working on %u of length %u.\n", (u64bit)U, idx, seqLen);
#endif

        //  Allocate space for the output -- 1MB should be enough for
        //  about 29000 signals.  Make it 32K -> 900 signals.
        //
        theOutputPos = 0;
        theOutputMax = 32 * 1024;
        theOutput    = new char [theOutputMax];
      
        //  Do searches.
        //
        if (config._doForward)
          doSearch(state, seq, seqLen, idx, false, theOutput, theOutputPos, theOutputMax);
        if (config._doReverse)
          doSearch(state, seq, seqLen, idx, true,  theOutput, theOutputPos, theOutputMax);

        //  Signal that we are done.
        //
        outputLen[idx] = theOutputPos;
        output[idx]    = theOutput;
        computed++;
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

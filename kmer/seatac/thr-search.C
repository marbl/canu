#include "seatac.H"

char const *srchGbye = "[%ld] computed: "uint64FMTW(8)"  blocked: "uint64FMTW(4)"/"uint64FMTW(4)"  encodeTime: %7.2f   searchTime: %7.2f   processTime: %7.2f\n";

class searcherState {
public:
  uint64         posnMax;
  uint64         posnLen;
  uint64        *posn;

  double         encodeTime;
  double         maskTime;
  double         searchTime;
  double         processTime;

  searcherState() {
    posnMax = 16384;
    posnLen = 0;
    posn    = new uint64 [ posnMax ];

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
         seqInCore *seq,
         uint32 idx,
         bool rc,
         filterObj *FO) {
  encodedQuery  *query      = 0L;
  hitMatrix     *matrix     = 0L;
  double         startTime  = 0.0;
  uint64         mer        = uint64ZERO;
  uint32         pos        = uint32ZERO;
  uint64         count      = 0;

  //  Build and mask the query
  //
  startTime = getTime();
  query = new encodedQuery(seq->sequence(), seq->sequenceLength(), config._merSize, rc);
  state->encodeTime += getTime() - startTime;

  //  Get the hits
  //
  startTime = getTime();
  matrix = new hitMatrix(seq->sequenceLength(), idx);

  while (query->getMer(mer, pos) == true)
    if (positions->getExact(mer, state->posn, state->posnMax, state->posnLen, count))
      matrix->addHits(pos, state->posn, state->posnLen);

  state->searchTime += getTime() - startTime;

  //  Begin processing
  //
  startTime = getTime();
  matrix->processMatrix(rc ? 'r' : 'f', FO);
  state->processTime += getTime() - startTime;

  delete matrix;
  delete query;
}



void*
searchThread(void *U) {
  uint32               idx      = 0;
  seqInCore           *seq      = 0L;
  uint32               blockedI = 0;
  uint32               blockedO = 0;
  uint32               computed = 0;

  searcherState       *state    = new searcherState;

  //  Allocate and fill out the thread stats -- this ensures that we
  //  always have stats (even if they're bogus).
  //
  threadStats[(long)U] = new char [1025];
  sprintf(threadStats[(long)U], srchGbye,
          (long)U,
          (uint32)0, (uint32)0, (uint32)0,
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
        //  fprintf(stderr, "%lu Blocked by input.\n", (uint64)U);
        blockedI++;
        nanosleep(&config._searchSleep, 0L);
      } else {

        //  If our idx is too far away from the output thread, sleep
        //  a little bit.  We keep the idx and seq that we have obtained,
        //  though.
        //
        while (idx > (outputPos + config._writerHighWaterMark)) {
          if (config._writerWarnings)
            fprintf(stderr, uint64FMT" Blocked by output (idx = "uint32FMT", outputPos = "uint32FMT").\n", (long)U, idx, outputPos);
          blockedO++;
          nanosleep(&config._searchSleep, 0L);
        }

        //  Construct a filter object
        //
        filterObj *FO = new filterObj(config._filterObj, config._filteropts);

        //  Do searches.
        //
        if (config._doForward)
          doSearch(state, seq, idx, false, FO);
        if (config._doReverse)
          doSearch(state, seq, idx, true, FO);

        //  Do filtering.
        //
        FO->filter();

        //  Signal that we are done.
        //
        output[idx] = FO;
        computed++;

        delete seq;

      } // end of seq != 0L
    } // end of idx < numberOfQueries
  } // end of inputTail < numberOfQueries


  //  OK, now fill out the read thread stats
  //
  sprintf(threadStats[(long)U], srchGbye, (long)U,
          computed, blockedI, blockedO,
          state->encodeTime,
          state->searchTime,
          state->processTime);

  delete state;

  return(0L);
}

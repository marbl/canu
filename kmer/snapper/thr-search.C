#include "posix.H"
#include "snapper2.H"
#include "time.H"
#include "sim4.H"

//  If you really, really, really want to know the exact number
//  of bases left in the query, use the interval list.  Otherwise,
//  it's faster to guess.
//
//#define USEEXACTSIZE


#ifdef TRUE64BIT
char const *srchGbye = "[%lu] processed: %8lu@%8lu/%8lu blocked: %4lu/%4lu Time: encode:%8.2f search:%8.2f chain:%8.2f filter:%8.2f polish:%8.2f\n";
#else
char const *srchGbye = "[%llu] processed: %8llu@%8llu/%8lu blocked: %4lu/%4lu Time: encode:%8.2f search:%8.2f chain:%8.2f filter:%8.2f polish:%8.2f\n";
#endif


class searcherState {
public:
  u64bit         posnMax;
  u64bit         posnLen;
  u64bit        *posn;

  double         encodeTime;
  double         maskTime;
  double         searchTime;
  double         chainTime;
  double         filterTime;
  double         polishTime;

  u64bit         searched;
  u64bit         polished;
  u64bit         discarded;

  searcherState() {
    posnMax    = 16384;
    posnLen    = 0;
    posn       = new u64bit [ posnMax ];

    encodeTime = 0.0;
    maskTime   = 0.0;
    searchTime = 0.0;
    chainTime  = 0.0;
    filterTime = 0.0;
    polishTime = 0.0;

    searched   = 0;
    polished   = 0;
    discarded  = 0;
  };

  ~searcherState() {
    delete [] posn;
  };
};



void
doSearch(searcherState       *state,
         FastASequenceInCore *seq,
         u32bit               idx,
         bool                 rc,
         aHit               *&theHits,
         u32bit              &theHitsLen,
         u32bit              &theHitsMax) {
  encodedQuery  *query  = 0L;
  hitMatrix     *matrix = 0L;
  u32bit         qMers  = 0;
  double         startTime  = 0.0;

  ///////////////////////////////////////
  //
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

  ///////////////////////////////////////
  //
  //  Get the hits
  //
  startTime = getTime();
  matrix = new hitMatrix(seq->sequenceLength(), qMers, idx);
  for (u32bit qi=0; qi<query->numberOfMers(); qi++)
    if ((query->getSkip(qi) == false) &&
        (positions->get(query->getMer(qi), state->posn, state->posnMax, state->posnLen)))
      matrix->addHits(qi, state->posn, state->posnLen);
  state->searchTime += getTime() - startTime;

  ///////////////////////////////////////
  //
  //  Chain the hits
  //
  startTime = getTime();
  matrix->filter(rc ? 'r' : 'f', theHits, theHitsLen, theHitsMax);
  state->chainTime += getTime() - startTime;

  delete matrix;
  delete query;
}






void
doFilter(searcherState       *state,
         FastASequenceInCore *seq,
         u32bit               idx,
         aHit               *&theHits,
         u32bit              &theHitsLen,
         u32bit              &theHitsMax) {
  double   startTime = getTime();

  if (theHitsLen == 0)
    return;

  //fprintf(stderr, "Found %u hits for %u\n", theHitsLen, idx);

  double L  = 0.5;
  double H  = 1.0;
  double V  = 0.6;

  //  Find the highest and lowest quality hit
  //
  u32bit  hiQ = theHits[0]._covered;
  u32bit  loQ = hiQ;

  for (u32bit i=0; i<theHitsLen; i++) {
    if (hiQ < theHits[i]._covered)
      hiQ = theHits[i]._covered;
    if (loQ > theHits[i]._covered)
      loQ = theHits[i]._covered;
  }

  double h = (double)(hiQ - loQ) / (double)theHits[0]._numMers;
  double p = 0.0;

  //  _numMers is not the same as the number covered.
  if (h > 1.0)
    h = 1.0;

  if (h <= L)    p = 1.0;
  if (h >= H)    p = V;
  if (p == 0.0)  p = 1.0 - (1.0 - V) * (h - L) / (H - L);

  //  check p; it should be between V and 1.0
  if (p > 1.0) {
    fprintf(stderr, "error in p; p=%f h=%f (%f %f %f)\n", p, h, L, H, V);
    p = 1.0;
  }

  if (p < V) {
    fprintf(stderr, "error in p; p=%f h=%f (%f %f %f)\n", p, h, L, H, V);
    p = V;
  }

  //  Filter the hits.  Any thing above cutL is good, and we should
  //  polish it.  Anything below is junk, and we should ignore it.
  //
  u32bit cutL = (u32bit)floor(hiQ - p * h);

  for (u32bit i=0; i < theHitsLen; i++)
    if (theHits[i]._covered < cutL)
      theHits[i]._direction |= 0x00000002;

  state->filterTime += getTime() - startTime;
}






char*
doPolish(searcherState       *state,
         FastASequenceInCore *seq,
         u32bit               idx,
         aHit               *&theHits,
         u32bit              &theHitsLen,
         u32bit              &theHitsMax) {
  double   startTime = getTime();
  u32bit   outputLen = 0;
  u32bit   outputMax = 1024 * 1024;
  char    *output    = new char [outputMax];

  output[0] = 0;

  for (u32bit i=0; i<theHitsLen; i++) {

    


    if ((config._doValidation) ||
        (theHits[i]._direction & 0x00000002)) {
      FastASequenceInCore  *ESTseq = seq;
      FastASequenceInCore  *GENseq = cache->getSequence(theHits[i]._dsIdx);
      u32bit                GENlo  = theHits[i]._dsLo;
      u32bit                GENhi  = theHits[i]._dsHi;

      //fprintf(stderr, "Polish %u vs %u[%u-%u]\n", idx, theHits[i]._dsIdx, GENlo, GENhi);
      //fprintf(stderr, "edef=%s\nddef=%s\n", ESTseq->header(), GENseq->header());

      bool    doForward =  theHits[i]._direction & 0x00000001;
      bool    doReverse = !doForward;

      doForward = doReverse = true;

      sim4command     *P4 = new sim4command(ESTseq,
                                            GENseq,
                                            GENlo,
                                            GENhi,
                                            doForward,
                                            doReverse);
      Sim4            *S4 = new Sim4(&sim4params);
      char            *O4 = S4->run(P4);

      if (O4 && O4[0]) {
        //fputs(O4, stderr);

        u32bit l = strlen(O4);
        if (outputLen + l + 1 > outputMax) {
          outputMax <<= 1;
          char *o = new char [outputMax];
          memcpy(o, output, sizeof(char) * outputLen);
          delete [] output;
          output = o;
        }
        memcpy(output + outputLen, O4, sizeof(char) * l);
        outputLen += l;
        output[outputLen] = 0;

        //  XXX:  Need to set the quality values, and dump the answer

        theHits[i]._direction |= 0x00000004;  //  Set the 'match' flag.
      } else {
        theHits[i]._direction &= 0x00000003;  //  Clear the 'match' flag, set everything else to zero.
      }

      delete [] O4;
      delete    S4;
      delete    P4;

      state->polished++;
    } else {
      state->discarded++;
    }
  }

  state->polishTime += getTime() - startTime;

  return(output);
}






void*
searchThread(void *U) {
  u32bit               idx      = 0;
  FastASequenceInCore *seq      = 0L;
  u32bit               blockedI = 0;
  u32bit               blockedO = 0;

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


        u32bit               theHitsLen = 0;
        u32bit               theHitsMax = 4;
        aHit                *theHits    = new aHit [theHitsMax];


        ///////////////////////////////////////
        //
        //  Do searches.
        //
        if (config._doForward)
          doSearch(state, seq, idx, false, theHits, theHitsLen, theHitsMax);
        if (config._doReverse)
          doSearch(state, seq, idx, true,  theHits, theHitsLen, theHitsMax);


        ///////////////////////////////////////
        //
        //  Filter the hits
        //
        doFilter(state, seq, idx, theHits, theHitsLen, theHitsMax);


        ///////////////////////////////////////
        //
        //  Polish the filtered hits
        //
        char *o = doPolish(state, seq, idx, theHits, theHitsLen, theHitsMax);


        ///////////////////////////////////////
        //
        //  Signal that we are done.
        //
        outputLen[idx] = strlen(o);
        output[idx]    = o;

        state->searched++;

        //  Clean up stuff
        //
        delete [] theHits;
        delete    seq;
      } // end of seq != 0L
    } // end of idx < numberOfQueries
  } // end of inputTail < numberOfQueries


  //  OK, now fill out the read thread stats
  //
  sprintf(threadStats[(u64bit)U], srchGbye, (u64bit)U,
          state->searched, state->polished, state->discarded,
          blockedI, blockedO,
          state->encodeTime,
          state->searchTime,
          state->chainTime,
          state->filterTime,
          state->polishTime);

  delete state;

  return(0L);
}

#include "posix.H"
#include "snapper2.H"
#include "bri++.H"
#include "sim4.H"

//  If you really, really, really want to know the exact number
//  of bases left in the query, use the interval list.  Otherwise,
//  it's faster to guess.
//
//#define USEEXACTSIZE

//  Define this to print a message whenever a search starts.
//
//#define VERBOSE_SEARCH

//  Define this to print a message whenever a polish starts.
//
//#define SHOW_POLISHING

//  Define these to show polishes that take a long time -- individual
//  polishes, not all polishes for a single sequence.  The time is in
//  seconds.
//
//#define SHOW_POLISHING_EXPENSIVE  0.5


char const *srchGbye = "["u64bitFMT"] processed: "u64bitFMTW(8)"@"u64bitFMTW(8)"/"u64bitFMTW(8)" blocked: "u64bitFMTW(4)"/"u64bitFMTW(4)" Time: encode:%8.2f search:%8.2f chain:%8.2f polish:%8.2f\n";


class searcherState {
public:
  u32bit         threadID;
  u64bit         posnMax;
  u64bit         posnLen;
  u64bit        *posn;

  double         encodeTime;
  double         maskTime;
  double         searchTime;
  double         chainTime;
  double         polishTime;

  u64bit         searched;
  u64bit         polished;
  u64bit         discarded;

  searcherState(u32bit U) {
    threadID   = U;

    posnMax    = 16384;
    posnLen    = 0;
    posn       = new u64bit [ posnMax ];

    encodeTime = 0.0;
    maskTime   = 0.0;
    searchTime = 0.0;
    chainTime  = 0.0;
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
  merCovering   *IL = new merCovering(config._merSize);

  for (u32bit qi=0; qi<query->numberOfMers(); qi++) {
    if (query->getSkip(qi) == false)
      IL->addMer(qi);
  }

  qMers = IL->sumLengths();
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





u32bit
configureFilter(double L,
                double H,
                double V,
                aHit  *theHits,
                u32bit theHitsLen) {

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

  double p = 0.0;

  //  _numMers is not the same as the number covered, so we should
  //  ensure that h is in range.
  //
  //  Note: _numMers is constant for all hits, so we can use any of them
  //
  double h = (double)(hiQ - loQ) / (double)theHits[0]._numMers;
  if (h > 1.0)
    h = 1.0;

  if (h <= L)    p = 1.0;
  if (h >= H)    p = V;
  if (p == 0.0)  p = 1.0 - (1.0 - V) * (h - L) / (H - L);

  if (p > 1.0) {
    fprintf(stderr, "error in p; p=%f h=%f (%f %f %f)\n", p, h, L, H, V);
    p = 1.0;
  }

  if (p < V) {
    fprintf(stderr, "error in p; p=%f h=%f (%f %f %f)\n", p, h, L, H, V);
    p = V;
  }

  //  Any thing at or above cutL is good, and we should polish it.
  //  Anything below is junk, and we should ignore it.
  //
  return((u32bit)floor(hiQ - p * h * theHits[0]._numMers));
}


u32bit
doFilter(searcherState       *state,
         aHit               *&theHits,
         u32bit              &theHitsLen) {

  if (theHitsLen == 0)
    return(0);

  u32bit numF = 0;
  u32bit cutL = configureFilter(config._Lo,
                                config._Hi,
                                config._Va, theHits, theHitsLen);

  for (u32bit i=0; i < theHitsLen; i++) {
    //  XXX:  This should be useless, but we also should guarantee that the flag is clear
    theHits[i]._status &= 0x00000001;

    //  If the coverage of the hit is more than the minimum, mark the
    //  hit as polishable.
    //
    if (theHits[i]._covered >= cutL) {
      numF++;
      theHits[i]._status |= 0x00000002;
    }
  }

  return(numF);
}






char*
doPolish(searcherState       *state,
         FastASequenceInCore *seq,
         aHit               *&theHits,
         u32bit              &theHitsLen) {
  double   startTime = getTime();
  u32bit   outputLen = 0;
  u32bit   outputMax = 2 * 1024 * theHitsLen;
  char    *output    = 0L;

  try {
    output = new char [outputMax];
  } catch (...) {
    fprintf(stderr, "doPolish()-- Can't allocate space for the output string (%u bytes) in thread %lu \n", outputMax, state->threadID);
    abort();
  }

  output[0] = 0;

#ifdef SHOW_POLISHING
  fprintf(stderr, "Need to polish "u32bitFMT" things.\n", theHitsLen);
#endif

  for (u32bit h=0; h<theHitsLen; h++) {
    if ((config._doValidation) ||
        (theHits[h]._status & 0x00000002)) {
      FastASequenceInCore  *ESTseq = seq;
      FastASequenceInCore  *GENseq = cache->getSequence(theHits[h]._dsIdx);
      u32bit                GENlo  = theHits[h]._dsLo;
      u32bit                GENhi  = theHits[h]._dsHi;

      bool    doForward =  theHits[h]._status & 0x00000001;
      bool    doReverse = !doForward;

#ifdef SHOW_POLISHING
      fprintf(stdout, "Hit %u out of %u (%u -> %u[%u-%u]) cov=%u matched=%u numMers=%u\n",
              h, theHitsLen,
              ESTseq->getIID(),
              theHits[h]._dsIdx,
              theHits[h]._dsLo,
              theHits[h]._dsHi,
              theHits[h]._covered,
              theHits[h]._matched,
              theHits[h]._numMers);
      fflush(stdout);
#endif


#ifdef SHOW_POLISHING_EXPENSIVE
      double startTime = getTime();
#endif

      sim4command     *P4 = new sim4command(ESTseq,
                                            GENseq,
                                            GENlo,
                                            GENhi,
                                            doForward,
                                            doReverse);
      Sim4            *S4 = new Sim4(&sim4params);
      sim4polishList  *l4 = S4->run(P4);
      sim4polishList  &L4 = *l4;

      //  Even though we don't expect multiple polishes, we still have to deal with
      //  them.  :-(

      //  Clear the 'match' flag and set qualities to zero.  XXX:
      //  Again, this should be already done, but we need to guarantee
      //  it.
      //
      theHits[h]._status &= 0x00000003;

      u32bit  pi = 0;
      u32bit  pc = 0;

      for (u32bit i=0; L4[i]; i++) {

        //  We need to remember the best pair of percent
        //  identity/coverage.  These wil be stored in the hit after
        //  we process all matches.
        //
        if ((L4[i]->percentIdentity >= pi) && (L4[i]->querySeqIdentity >= pc)) {
          pi = L4[i]->percentIdentity;
          pc = L4[i]->querySeqIdentity;
        }

#ifdef SHOW_POLISHING
        fprintf(stdout, "  match["u32bitFMT"] : id=%u cov=%u\n", i, L4[i]->percentIdentity, L4[i]->querySeqIdentity);
        fflush(stdout);
#endif

        //  If we have a real hit, set the flag and save the output
        //
        if ((L4[i]->percentIdentity  >= config._minMatchIdentity) &&
            (L4[i]->querySeqIdentity >= config._minMatchCoverage)) {

#ifdef SHOW_POLISHING
          fprintf(stderr, "    saving hit 1\n");
#endif

          theHits[h]._status |= 0x00000004;

          char *pstr = s4p_polishToString(L4[i]);

#ifdef SHOW_POLISHING
          fprintf(stderr, "    saving hit 2\n");
#endif

          u32bit l = (u32bit)strlen(pstr);
          if (outputLen + l + 1 >= outputMax) {
            outputMax = outputMax + outputMax + l;
            char *o = 0L;
            try {
              o = new char [outputMax];
            } catch (...) {
              fprintf(stderr, "doPolish()-- Can't reallocate space for the output string ("u32bitFMT" bytes) in thread "u64bitFMT"\n", outputMax, state->threadID);
              abort();
            }
            memcpy(o, output, sizeof(char) * outputLen);
            delete [] output;
            output = o;
          }

#ifdef SHOW_POLISHING
          fprintf(stderr, "    saving hit 3 len="u32bitFMT", max="u32bitFMT", l="u32bitFMT"\n", outputLen, outputMax, l);
#endif

          memcpy(output + outputLen, pstr, sizeof(char) * l);
          outputLen += l;
          output[outputLen] = 0;

#ifdef SHOW_POLISHING
          fprintf(stderr, "    saving hit 4\n");
#endif

          free(pstr);

#ifdef SHOW_POLISHING
          fprintf(stderr, "    saving hit 5\n");
#endif
        }
      }

#ifdef SHOW_POLISHING
      fprintf(stderr, "    saving hit 6\n");
#endif

      //  Save the best scores
      //
      theHits[h]._status |= pi << 16;
      theHits[h]._status |= pc << 24;

      delete l4;
      delete S4;
      delete P4;

#ifdef SHOW_POLISHING_EXPENSIVE
      double elapsedTime = getTime() - startTime;
      if (elapsedTime >= SHOW_POLISHING_EXPENSIVE) {
        fprintf(stdout, "Hit %u out of %u (%u -> %u[%u-%u]) took %f seconds ().\n",
                h, theHitsLen,
                ESTseq->getIID(), GENseq->getIID(), theHits[h]._dsLo, theHits[h]._dsHi,
                elapsedTime);
        fflush(stdout);
      }
#endif

      state->polished++;
    } else {
      state->discarded++;
    }
  }

#ifdef SHOW_POLISHING
  fprintf(stderr, "Done polishing.\n");
#endif

  state->polishTime += getTime() - startTime;

  return(output);
}






void*
searchThread(void *U) {
  u32bit               idx      = 0;
  FastASequenceInCore *seq      = 0L;
  u32bit               blockedI = 0;
  u32bit               blockedO = 0;

  searcherState       *state    = new searcherState((u64bit)U);

  fprintf(stderr, "Hello!  I'm a searchThread (number %u)!\n", (u64bit)U);

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

#ifdef VERBOSE_SEARCH
        double startTime = getTime();
#endif

        ///////////////////////////////////////
        //
        //  Do searches.
        //
        if (config._doForward)
          doSearch(state, seq, idx, false, theHits, theHitsLen, theHitsMax);
        if (config._doReverse)
          doSearch(state, seq, idx, true,  theHits, theHitsLen, theHitsMax);

#ifdef VERBOSE_SEARCH
        double searchTime = getTime();
#endif

        ///////////////////////////////////////
        //
        //  Filter the hits
        //
        u32bit numF = doFilter(state, theHits, theHitsLen);

#ifdef VERBOSE_SEARCH
        fprintf(stdout, "seq %u has %u decent hits out of %u total\n", idx, numF, theHitsLen);
        fflush(stdout);

        double filterTime = getTime();
#endif

        ///////////////////////////////////////
        //
        //  Polish the filtered hits
        //
        char *o = doPolish(state, seq, theHits, theHitsLen);

        ///////////////////////////////////////
        //
        //  Signal that we are done.
        //
        answerLen[idx] = theHitsLen;
        answer[idx]    = theHits;
        outputLen[idx] = (u32bit)strlen(o);
        output[idx]    = o;

        state->searched++;

        //  Clean up stuff
        //
        delete    seq;

#ifdef VERBOSE_SEARCH
        double endTime = getTime();
        if ((endTime - startTime) > 5.0) {
          fprintf(stdout, "%lu done with idx = %d (%u bp, %u hits) search %f filter %f polish %f  encode %f search %f chain %f polish %f\n",
                  (u64bit)U, idx, seq->sequenceLength(), theHitsLen,
                  searchTime - startTime,
                  filterTime - searchTime,
                  endTime - filterTime,
                  state->encodeTime,
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
  sprintf(threadStats[(u64bit)U], srchGbye, (u64bit)U,
          state->searched, state->polished, state->discarded,
          blockedI, blockedO,
          state->encodeTime,
          state->searchTime,
          state->chainTime,
          state->polishTime);

  delete state;

  return(0L);
}

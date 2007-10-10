#include "snapper2.H"


void
doSearch(searcherState       *state,
         seqInCore           *seq,
         u32bit               idx,
         bool                 rc,
         aHit               *&theHits,
         u32bit              &theHitsLen,
         u32bit              &theHitsMax,
         logMsg              *theLog) {
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
      IL->add(qi, config._merSize);
  }

  IL->merge();

  qMers = IL->sumOfLengths();
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
  for (u32bit qi=0; qi<query->numberOfMers(); qi++) {
    if ((query->getSkip(qi) == false) &&
        (positions->get(query->getMer(qi), state->posn, state->posnMax, state->posnLen))) {
      matrix->addHits(qi, state->posn, state->posnLen);
    }
  }
  state->searchTime += getTime() - startTime;




  ///////////////////////////////////////
  //
  //  Chain the hits
  //
  startTime = getTime();
  matrix->filter(rc ? 'r' : 'f', config._minHitCoverage, config._minHitLength,
                 theHits, theHitsLen, theHitsMax);
  state->chainTime += getTime() - startTime;



  ////////////////////////////////////////
  //
  //  Refine the hits -- if any hit looks like it contains a repeat,
  //  rebuild it using an adaptive mask threshold.
  //
  //  We work backwards because we add on new hits to the end of our
  //  list.
  //
  for (u32bit h=theHitsLen; h--; ) {

    //  The first test eliminates hits that were not generated for the
    //  complementarity used in this search (e.g., the first search
    //  does rc=forward, adds some hits, the second search does
    //  rc=reverse, and we should skip all the rc=forward hits.
    //  
    if (((theHits[h]._status & AHIT_DIRECTION_MASK) == !rc) && 
        (theHits[h]._matched > 2 * theHits[h]._numMers)) {

#ifdef SHOW_HIT_DISCARDING
      theLog->add("Seq "u32bitFMT" Hit "u32bitFMT" (%c) has "u32bitFMT" matched, but only "u32bitFMT" mers.\n",
                 seq->getIID(), h, rc ? 'r' : 'f', theHits[h]._matched, theHits[h]._numMers);
#endif

      //  Grab the genomic sequence.
      //  Construct a merstream for the region.
      //  Build a positionDB of the region (both positions and counts).
      //  Fill out another hitMatrix using about 2*length mers.
      //
      seqInCore            *GENseq = cache->getSequenceInCore(theHits[h]._dsIdx);
      u32bit                GENlo  = theHits[h]._dsLo;
      u32bit                GENhi  = theHits[h]._dsHi;

#warning this kMerBuilder should be in config
      kMerBuilder           KB(config._merSize);

      merStream            *MS     = new merStream(&KB, GENseq, GENlo, GENhi - GENlo);
      positionDB           *PS     = new positionDB(MS, config._merSize, 0, 0L, 0L, 0, false);
      hitMatrix            *HM     = new hitMatrix(seq->sequenceLength(), qMers, idx);

      //  We find the number of hits we would get if we use a
      //  countLimit of i.
      //
#define COUNT_MAX   256

      u32bit numHitsAtCount[COUNT_MAX] = { 0 };
      u32bit countLimit                = 0;

      u32bit numMers = 0;
#ifdef SHOW_HIT_DISCARDING
      u32bit numHits = 0;
      u32bit minNum  = ~u32bitZERO;
      u32bit maxNum  = 0;
#endif

      for (u32bit qi=0; qi<query->numberOfMers(); qi++) {
        if ((query->getSkip(qi) == false) &&
            (PS->get(query->getMer(qi), state->posn, state->posnMax, state->posnLen))) {
          numMers++;

          if (state->posnLen < COUNT_MAX)
            numHitsAtCount[state->posnLen] += state->posnLen;

#ifdef SHOW_HIT_DISCARDING
          numHits += state->posnLen;
          if (minNum > state->posnLen)  minNum = state->posnLen;
          if (maxNum < state->posnLen)  maxNum = state->posnLen;
#endif
        }
      }

      //  Scan the number of hits at count, pick the first highest
      //  count such that the number of hits is below our threshold.
      //
      for (u32bit qi=1; qi<COUNT_MAX; qi++) {
        numHitsAtCount[qi] = numHitsAtCount[qi-1] + numHitsAtCount[qi];

        if (numHitsAtCount[qi] <= numMers * config._repeatThreshold)
          countLimit = qi;
      }

#ifdef SHOW_HIT_DISCARDING
      theLog->add(" -- found "u32bitFMT" hits in "u32bitFMT" mers, min="u32bitFMT" max="u32bitFMT" avg=%.5f hits/mer.\n",
                 numHits, numMers, minNum, maxNum, (double)numHits / (double)numMers);
      theLog->add(" -- using a countLimit of "u32bitFMT" which gets us "u32bitFMT" mers\n",
                 countLimit, numHitsAtCount[countLimit]);
#endif

      for (u32bit qi=0; qi<query->numberOfMers(); qi++) {
        if ((query->getSkip(qi) == false) &&
            (PS->get(query->getMer(qi), state->posn, state->posnMax, state->posnLen))) {
          if (state->posnLen <= countLimit) {
            for (u32bit x=0; x<state->posnLen; x++)
              state->posn[x] += GENlo + config._useList.startOf(theHits[h]._dsIdx);

            //  The kmer counts for these mers are relative to the
            //  sub-regions, not the global, so we want to disable any
            //  filtering by kmer counts.  We could add a flag to the filter
            //  to stop this, or we can reset the counts here to large
            //  values.  Or we could simply reset the counts to the global
            //  value.
            //
            HM->addHits(qi, state->posn, state->posnLen, positions->count(query->getMer(qi)));
          }
        }
      }

      //  Chain the hits
      //
      HM->filter(rc ? 'r' : 'f', 0.01, 0, theHits, theHitsLen, theHitsMax);

      //  Mark this hit as dead
      //
      theHits[h]._status |= AHIT_DISCARDED;

      delete HM;
      delete PS;
      delete MS;
    }
  }

  delete matrix;
  delete query;
}





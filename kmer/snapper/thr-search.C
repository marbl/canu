#include "snapper2.H"





class encodedQuery {
private:
  u64bit   *_mers;
  u32bit   *_posn;
  u32bit   *_span;
  u32bit    _mersActive;
  u32bit    _mersInQuery;
public:

  encodedQuery(seqInCore           *seq,
               kMerBuilder         *KB,
               bool                 rc) {
    _mers        = new u64bit [seq->sequenceLength()];
    _posn        = new u32bit [seq->sequenceLength()];
    _span        = new u32bit [seq->sequenceLength()];
    _mersActive  = 0;
    _mersInQuery = 0;

    //  Unfortunately, we need to use the slightly heavyweight merStream
    //  and kMerBuilder to get mers.  We used to build mers in a tight
    //  loop, but with the inclusion of spacing and compression, we
    //  cannot do that anymore.

    merStream  *MS = new merStream(KB, seq);
    u64bit      mer;
    u32bit      val;

    //  The rc flag tells us if we should build for the forward or
    //  reverse strand.  If forward (rc == false) the mers are in the
    //  same order.  If reverse, the mers are both reverse-complemented,
    //  and appear in our mers[] and skip[] lists reversed.

    if (rc == false) {
      while (MS->nextMer()) {
        mer = MS->theFMer();

        if ((maskDB && (maskDB->exists(mer) == true)) ||
            (onlyDB && (onlyDB->exists(mer) == false)))
          ;  //  Don't use it.
        else {
          _mers[_mersActive] = mer;
          _posn[_mersActive] = MS->thePositionInSequence();
          _span[_mersActive] = MS->theFMer().getMerSpan();
          _mersActive++;
        }

        _mersInQuery++;
      }
    } else {
      while (MS->nextMer()) {
        mer = MS->theRMer();

        if ((maskDB && (maskDB->exists(mer) == true)) ||
            (onlyDB && (onlyDB->exists(mer) == false)))
          ;  //  Don't use it.
        else {
          //  We die horribly unless we do the goofy math to get the
          //  _posn.  I'm sure that could be cleaned up, but it'd take
          //  more effort than I want now (being we'd have to figure
          //  out what the search/hitMatrix stuff is doing).
          _mers[_mersActive] = mer;
          _posn[_mersActive] = seq->sequenceLength() - MS->thePositionInSequence() - MS->theRMer().getMerSpan();
          _span[_mersActive] = MS->theRMer().getMerSpan();
          _mersActive++;
        }

        _mersInQuery++;
      }

      //  Reverse the array -- this appears to be optional.
#if 1
      if (_mersActive > 0)
        for (u32bit i=0, j=_mersActive-1; i<j; i++, j--) {
          mer      = _mers[i];
          _mers[i] = _mers[j];
          _mers[j] = mer;

          val      = _posn[i];
          _posn[i] = _posn[j];
          _posn[j] = val;

          val      = _span[i];
          _span[i] = _span[j];
          _span[j] = val;
        }
#endif
    }

    delete MS;
  };


  ~encodedQuery() {
    delete [] _mers;
    delete [] _posn;
    delete [] _span;
  };

  u32bit           numberOfMersActive(void)    { return(_mersActive);  };
  u32bit           numberOfMersInQuery(void)   { return(_mersInQuery); };

  u64bit           getMer(u32bit i)            { return(_mers[i]);     };
  u32bit           getPosn(u32bit i)           { return(_posn[i]);     };
  u32bit           getSpan(u32bit i)           { return(_span[i]);     };
};











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
  double         startTime  = 0.0;

  if (state->KB == 0L)
    state->KB = new kMerBuilder(config._KBmerSize,
                                config._KBcompression,
                                config._KBspacingTemplate);

  startTime = getTime();

  ///////////////////////////////////////
  //
  //  Build and mask the query
  //
  query = new encodedQuery(seq, state->KB, rc);


  ///////////////////////////////////////
  //
  //  Get the hits
  //
  startTime = getTime();
  matrix = new hitMatrix(seq->sequenceLength(), query->numberOfMersInQuery(), idx, theLog);
  for (u32bit qidx=0; qidx<query->numberOfMersActive(); qidx++) {
    if (positions->get(query->getMer(qidx), state->posn, state->posnMax, state->posnLen)) {
#if 0
      fprintf(stderr, "rc=%d qidx="u32bitFMT" pos="u32bitFMT" mer="u64bitHEX" hits="u64bitFMT"\n",
              rc, qidx, query->getPosn(qidx), query->getMer(qidx), state->posnLen);
#endif
      matrix->addHits(query->getPosn(qidx), state->posn, state->posnLen);
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

      merStream            *MS     = new merStream(state->KB, GENseq, GENlo, GENhi - GENlo);
      positionDB           *PS     = new positionDB(MS, config._KBmerSize, 0, 0L, 0L, 0, false);
      hitMatrix            *HM     = new hitMatrix(seq->sequenceLength(), query->numberOfMersInQuery(), idx, theLog);

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

      for (u32bit qidx=0; qidx<query->numberOfMersActive(); qidx++) {
        if (PS->get(query->getMer(qidx), state->posn, state->posnMax, state->posnLen)) {
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
      for (u32bit qidx=1; qidx<COUNT_MAX; qidx++) {
        numHitsAtCount[qidx] = numHitsAtCount[qidx-1] + numHitsAtCount[qidx];

        if (numHitsAtCount[qidx] <= numMers * config._repeatThreshold)
          countLimit = qidx;
      }

#ifdef SHOW_HIT_DISCARDING
      theLog->add(" -- found "u32bitFMT" hits in "u32bitFMT" mers, min="u32bitFMT" max="u32bitFMT" avg=%.5f hits/mer.\n",
                 numHits, numMers, minNum, maxNum, (double)numHits / (double)numMers);
      theLog->add(" -- using a countLimit of "u32bitFMT" which gets us "u32bitFMT" mers\n",
                 countLimit, numHitsAtCount[countLimit]);
#endif

      for (u32bit qidx=0; qidx<query->numberOfMersActive(); qidx++) {
        if (PS->get(query->getMer(qidx), state->posn, state->posnMax, state->posnLen)) {
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
            HM->addHits(query->getPosn(qidx), state->posn, state->posnLen, positions->count(query->getMer(qidx)));
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





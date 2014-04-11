#include "snapper2.H"

#if defined (__SVR4) && defined (__sun) 
// Solaris defines SS in sys/regset.h
#undef SS
#endif

class encodedQuery {
private:
  uint64   *_mers;
  uint32   *_posn;
  uint32   *_span;
  uint32    _mersActive;
  uint32    _mersInQuery;

public:
  encodedQuery(seqInCore           *seq,
               kMerBuilder         *KB,
               bool                 rc) {
    _mers        = new uint64 [seq->sequenceLength()];
    _posn        = new uint32 [seq->sequenceLength()];
    _span        = new uint32 [seq->sequenceLength()];
    _mersActive  = 0;
    _mersInQuery = 0;

    //  Unfortunately, we need to use the slightly heavyweight merStream
    //  and kMerBuilder to get mers.  We used to build mers in a tight
    //  loop, but with the inclusion of spacing and compression, we
    //  cannot do that anymore.

    seqStream  *SS = new seqStream(seq->sequence(), seq->sequenceLength());
    merStream  *MS = new merStream(KB, SS);
    uint64      mer;
    uint32      val;

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
        for (uint32 i=0, j=_mersActive-1; i<j; i++, j--) {
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
    delete SS;
  };


  ~encodedQuery() {
    delete [] _mers;
    delete [] _posn;
    delete [] _span;
  };

  uint32           numberOfMersActive(void)    { return(_mersActive);  };
  uint32           numberOfMersInQuery(void)   { return(_mersInQuery); };

  uint64           getMer(uint32 i)            { return(_mers[i]);     };
  uint32           getPosn(uint32 i)           { return(_posn[i]);     };
  uint32           getSpan(uint32 i)           { return(_span[i]);     };
};





void
doSearch(searcherState       *state,
         query               *qry,
         bool                 rc) {

  if (state->KB == 0L)
    state->KB = new kMerBuilder(config._KBmerSize,
                                config._KBcompression,
                                config._KBspacingTemplate);

  encodedQuery *encqry = new encodedQuery(qry->seq, state->KB, rc);

  hitMatrix    *matrix = new hitMatrix(qry->seq->sequenceLength(),
                                       encqry->numberOfMersInQuery(),
                                       qry->seq->getIID(),
                                       qry->theLog);

  for (uint32 qidx=0; qidx<encqry->numberOfMersActive(); qidx++) {
    uint64  count = 0;

    if (positions->getExact(encqry->getMer(qidx), state->posn, state->posnMax, state->posnLen, count))
      matrix->addHits(encqry->getPosn(qidx), state->posn, state->posnLen);
  }

  //  Chain the hits
  //
  matrix->filter(rc ? 'r' : 'f', config._minHitCoverage, config._minHitLength, qry->theHits, qry->theHitsLen, qry->theHitsMax);


  ////////////////////////////////////////
  //
  //  Refine the hits -- if any hit looks like it contains a repeat,
  //  rebuild it using an adaptive mask threshold.
  //
  //  We work backwards because we add on new hits to the end of our
  //  list.
  //
  for (uint32 h=qry->theHitsLen; h--; ) {

    //  The first test eliminates hits that were not generated for the
    //  complementarity used in this search (e.g., the first search
    //  does rc=forward, adds some hits, the second search does
    //  rc=reverse, and we should skip all the rc=forward hits.
    //  
    if (((qry->theHits[h]._status & AHIT_DIRECTION_MASK) == !rc) && 
        (qry->theHits[h]._matched > 2 * qry->theHits[h]._numMers)) {

#ifdef SHOW_HIT_DISCARDING
      qry->theLog->add("Seq "uint32FMT" Hit "uint32FMT" (%c) has "uint32FMT" matched, but only "uint32FMT" mers.\n",
                 seq->getIID(), h, rc ? 'r' : 'f', qry->theHits[h]._matched, qry->theHits[h]._numMers);
#endif

      //  Grab the genomic sequence.
      //  Construct a merstream for the region.
      //  Build a positionDB of the region (both positions and counts).
      //  Fill out another hitMatrix using about 2*length mers.
      //
      seqInCore            *GENseq = genome->getSequenceInCore(qry->theHits[h]._dsIdx);
      uint32                GENlo  = qry->theHits[h]._dsLo;
      uint32                GENhi  = qry->theHits[h]._dsHi;

      merStream            *MS     = new merStream(state->KB,
                                                   new seqStream(GENseq->sequence(), GENseq->sequenceLength()),
                                                   false, true);

      MS->setBaseRange(GENlo, GENhi);

      positionDB           *PS     = new positionDB(MS, config._KBmerSize, 0, 0L, 0L, 0L, 0, 0, 0, 0, false);
      hitMatrix            *HM     = new hitMatrix(qry->seq->sequenceLength(),
                                                   encqry->numberOfMersInQuery(),
                                                   qry->seq->getIID(),
                                                   qry->theLog);

      //  We find the number of hits we would get if we use a
      //  countLimit of i.
      //
#define COUNT_MAX   256

      uint32 numHitsAtCount[COUNT_MAX] = { 0 };
      uint32 countLimit                = 0;
      uint64 count                     = 0;

      uint32 numMers = 0;
#ifdef SHOW_HIT_DISCARDING
      uint32 numHits = 0;
      uint32 minNum  = ~uint32ZERO;
      uint32 maxNum  = 0;
#endif

      for (uint32 qidx=0; qidx<encqry->numberOfMersActive(); qidx++) {
        if (PS->getExact(encqry->getMer(qidx), state->posn, state->posnMax, state->posnLen, count)) {
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
      for (uint32 qidx=1; qidx<COUNT_MAX; qidx++) {
        numHitsAtCount[qidx] = numHitsAtCount[qidx-1] + numHitsAtCount[qidx];

        if (numHitsAtCount[qidx] <= numMers * config._repeatThreshold)
          countLimit = qidx;
      }

#ifdef SHOW_HIT_DISCARDING
      qry->theLog->add(" -- found "uint32FMT" hits in "uint32FMT" mers, min="uint32FMT" max="uint32FMT" avg=%.5f hits/mer.\n",
                 numHits, numMers, minNum, maxNum, (double)numHits / (double)numMers);
      qry->theLog->add(" -- using a countLimit of "uint32FMT" which gets us "uint32FMT" mers\n",
                 countLimit, numHitsAtCount[countLimit]);
#endif

      for (uint32 qidx=0; qidx<encqry->numberOfMersActive(); qidx++) {
        if (PS->getExact(encqry->getMer(qidx), state->posn, state->posnMax, state->posnLen, count)) {
          if (state->posnLen <= countLimit) {
            for (uint32 x=0; x<state->posnLen; x++)
              state->posn[x] += genomeMap->startOf(qry->theHits[h]._dsIdx);

            //  The kmer counts for these mers are relative to the
            //  sub-regions, not the global, so we want to disable any
            //  filtering by kmer counts.  We could add a flag to the filter
            //  to stop this, or we can reset the counts here to large
            //  values.  Or we could simply reset the counts to the global
            //  value.
            //
            HM->addHits(encqry->getPosn(qidx), state->posn, state->posnLen, positions->countExact(encqry->getMer(qidx)));
          }
        }
      }

      //  Chain the hits
      //
      HM->filter(rc ? 'r' : 'f', 0.01, 0, qry->theHits, qry->theHitsLen, qry->theHitsMax);

      //  Mark this hit as dead
      //
      qry->theHits[h]._status |= AHIT_DISCARDED;

      delete HM;
      delete PS;
      delete MS;

      delete GENseq;
    }
  }

  delete matrix;
  delete encqry;
}





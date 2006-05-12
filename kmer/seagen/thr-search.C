#include "searchGENOME.H"
#include "encodedQuery.H"

//  If you really, really, really want to know the exact number
//  of bases left in the query, use the interval list.  Otherwise,
//  it's faster to guess.
//
//#define USEEXACTSIZE

void
doSearch(searcherState *state,
         encodedQuery  *query,
         bool           isReverse) {

  //  Get the hits
  double startTime = getTime();

  hitMatrix *matrix = new hitMatrix(query->bpTotal(),
                                    query->bpCovered(false),
                                    query->IID());

  for (u32bit qi=0; qi<query->numberOfMers(); qi++)
    if ((query->getSkip(qi, isReverse) == false) &&
        (config._positions->get(query->getMer(qi, isReverse),
                                state->posn,
                                state->posnMax,
                                state->posnLen)))
      matrix->addHits(qi, state->posn, state->posnLen);

  state->searchTime += getTime() - startTime;


  //  Filter, storing the resutls into theOutput
  startTime = getTime();

  matrix->filter(query, isReverse);
  delete matrix;

  state->filterTime += getTime() - startTime;
}



void
searchThread(void *U, void *T, void *Q) {
  searcherState *state = (searcherState *)T;
  encodedQuery  *query = (encodedQuery *)Q;

  //  Finish building the query -- mask out repetitive junk
  //
  double startTime = getTime();

  if (config._maskDB)
    for (u32bit qi=0; qi<query->numberOfMers(); qi++)
      if ((query->getSkip(qi, false) == false) &&
          (config._maskDB->exists(query->getMer(qi, false))))
        query->setSkip(qi, false);

  if (config._onlyDB)
    for (u32bit qi=0; qi<query->numberOfMers(); qi++)
      if ((query->getSkip(qi, false) == false) &&
          (!config._onlyDB->exists(query->getMer(qi, false))))
        query->setSkip(qi, false);

  state->maskTime += getTime() - startTime;


  //  Do searches.
  //
  if (config._doForward)
    doSearch(state, query, false);
  if (config._doReverse)
    doSearch(state, query, true);
}

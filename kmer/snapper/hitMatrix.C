#include "snapper2.H"

#define MINCOUNT  3

hitMatrix::hitMatrix(uint32 qsLen, uint32 qsMers, uint32 qsIdx, logMsg *theLog) {
  _qsLen    = qsLen;
  _qsMers   = qsMers;
  _qsIdx    = qsIdx;

  _hitsLen  = 0;
  _hitsMax  = 8;
  _hits     = new diagonalLine [_hitsMax];

  _matches  = 0L;

  _theLog   = theLog;
}

hitMatrix::~hitMatrix() {
  delete [] _hits;
}


void
hitMatrix::addMatch(uint32         isunique,
                    uint32         qsLo,
                    uint32         qsHi,
                    uint32         dsLo,
                    uint32         dsHi,
                    merCovering   *IL,
                    merList       *ML) {
  uint32 offset = 0;

  offset = (uint32)(config._extendWeight * qsLo);
  if (offset < config._extendMinimum)
    offset = config._extendMinimum;
  if (offset > config._extendMaximum)
    offset = config._extendMaximum;
  if (dsLo < offset)
    dsLo = 0;
  else
    dsLo -= offset;

  offset = (uint32)(config._extendWeight * (_qsLen - qsHi));
  if (offset < config._extendMinimum)
    offset = config._extendMinimum;
  if (offset > config._extendMaximum)
    offset = config._extendMaximum;
  dsHi += offset;


  //  Create a new match
  //
  //  n = new match
  //  m = current match
  //  l = last match
  //
  trapMatch *n = new trapMatch(isunique, qsLo, qsHi, dsLo, dsHi, IL, ML);

#ifdef SHOW_HITMATRIX
  _theLog->add("chained:  Q::"uint32FMT"-"uint32FMT"("uint32FMT") G::"uint32FMT"-"uint32FMT"("uint32FMT")\n",
               qsLo, qsHi, qsHi - qsLo,
               dsLo, dsHi, dsHi - dsLo);
#endif

  //  And find a home for it in the list.  No merging of matches is done here.  It's
  //  too hard.
  //
  if ((_matches == 0L) || (n->_dsHi > _matches->_dsHi)) {
    n->_next = _matches;
    _matches = n;
  } else {
    trapMatch *l = _matches;
    trapMatch *m = _matches->_next;

    while ((m) && (n->_dsHi < m->_dsHi)) {
      l = m;
      m = m->_next;
    }

    n->_next = m;
    l->_next = n;
  }
}


//  Utility for sorting the diagonal lines in the hitMatrix
//
//  The two comparison functions return true if the first line
//  is less than the second line.

inline
int
compareLines(diagonalLine *A, diagonalLine *B, uint32 qsLen) {
  uint32 a = qsLen - A->val.qPos - 1 + A->val.dPos;
  uint32 b = qsLen - B->val.qPos - 1 + B->val.dPos;

  return(((a  < b)) ||
         ((a == b) && (A->val.qPos < B->val.qPos)));
}

inline
int
compareLines(uint32 l, uint32 q, diagonalLine *B, uint32 qsLen) {
  uint32 b = qsLen - B->val.qPos - 1 + B->val.dPos;

  return(((l  < b)) ||
         ((l == b) && (q < B->val.qPos)));
}

inline
void
adjustHeap(diagonalLine *L, int32 p, int32 n, uint32 qsLen) {
  uint64  v = L[p].all;
  uint32  q = L[p].val.qPos;
  uint32  l = qsLen - q - 1 + L[p].val.dPos;
  int32  c = (p << 1) + 1;  //  let c be the left child of p

  while (c < n) {

    //  Find the larger of the two children
    //
    if ((c+1 < n) && compareLines(L+c, L+c+1, qsLen))
      c++;

    //  Does the node in question fit here?
    //
    if (compareLines(l, q, L+c, qsLen) == false)
      break;

    //  Else, swap the parent and the child
    //
    L[p].all = L[c].all;

    //  Move down the tree
    //
    p = c;
    c = (p << 1) + 1;
  }

  L[p].all = v;
}




void
hitMatrix::filter(char      direction,
                  double    minHitCoverage,
                  uint32    minHitLength,
                  aHit    *&theHits,
                  uint32   &theHitsPos,
                  uint32   &theHitsMax) {

  if (_hitsLen == 0)
    return;

  //  Decide on the minimum quality values; we pick the larger of
  //  the fixed lengths, and the sequence length * coverage.
  //
  uint32   minLength = (uint32)(minHitCoverage * _qsLen);
  if (minLength < minHitLength)
    minLength = minHitLength;

  //  First, sort by the dsPos.  This is done so that we can find all the hits for
  //  a specific scaffold.
  //
  sort_dsPos();

  //  Now, while there are hits left....
  //
  uint32  firstHit   = 0;
  uint32  lastHit    = 0;
  uint32  currentSeq = 0;

  //
  //  Step 1:  Sort the mer-hits, chain, promote decent ones to matches
  //

  while (firstHit < _hitsLen) {

    //  Move the currentSeq until the firstHit is below it.  After
    //  this loop, currentSeq is the sequence AFTER the one that we
    //  want hits in.
    //
    while ((currentSeq < genomeMap->numberOfSequences()) &&
           (genomeMap->startOf(currentSeq) <= _hits[firstHit].val.dPos))
      currentSeq++;

    //  Find the first hit that is in currentSeq.  If this is the last sequence,
    //  then, of course, all remaining hits are in it.
    //
    if (currentSeq < genomeMap->numberOfSequences()) {
      lastHit = firstHit + 1;
      while ((lastHit < _hitsLen) &&
             (_hits[lastHit].val.dPos < genomeMap->startOf(currentSeq)))
        lastHit++;
    } else {
      lastHit = _hitsLen;
    }

    //  Drop back one sequence; this is the sequence the hits are in.
    //
    currentSeq--;

    //  Adjust the hits to be relative to the start of this sequence
    //
    for (uint32 i=firstHit; i<lastHit; i++)
      _hits[i].val.dPos -= genomeMap->startOf(currentSeq);

    //  Sort them, if needed.
    //
    if (lastHit - firstHit > 1) {

      //  We cheat; heapsort isn't too friendly to sorting the middle of
      //  an array, so we make a new array in the middle!
      //
      diagonalLine  *hitsToSort = _hits + firstHit;

      //  Build the heap.  I initially thought this could be done at the
      //  same time as the scan for the last hit, but it can't (easily)
      //
      for (int32 i=(lastHit - firstHit)/2 - 1; i>=0; i--)
        adjustHeap(hitsToSort, i, lastHit - firstHit, _qsLen);

      //  Sort the hits by diagonal.  This is the second part of
      //  heap sort -- Interchange the new maximum with the element
      //  at the end of the tree
      //
      for (uint32 i=lastHit - firstHit - 1; i>0; i--) {
        uint64 v            = hitsToSort[i].all;
        hitsToSort[i].all   = hitsToSort[0].all;
        hitsToSort[0].all   = v;
      
        adjustHeap(hitsToSort, 0, i, _qsLen);
      }
    }

    //  Filter them
    //
    uint32  frstDiagonal = _qsLen - _hits[firstHit].val.qPos - 1 + _hits[firstHit].val.dPos;
    uint32  lastDiagonal = frstDiagonal;
    uint32  unique       = uint32ZERO;
    uint32  qsLow        = _hits[firstHit].val.qPos;
    uint32  qsHigh       = _hits[firstHit].val.qPos;
    uint32  dsLow        = _hits[firstHit].val.dPos;
    uint32  dsHigh       = _hits[firstHit].val.dPos;
    uint32  minCount     = ~uint32ZERO;

    merCovering   *IL = new merCovering(config._KBmerSize);
    merList       *ML = new merList();

    for (uint32 i=firstHit; i<lastHit; i++) {
      uint32 thisDiagonalID = _qsLen - _hits[i].val.qPos - 1 + _hits[i].val.dPos;

      //  Unconditionally extend if the diagonal difference is small.
      //
      if (lastDiagonal + config._maxDiagonal >= thisDiagonalID) {
        lastDiagonal = thisDiagonalID;
        if (qsLow    > _hits[i].val.qPos)   qsLow    = _hits[i].val.qPos;
        if (qsHigh   < _hits[i].val.qPos)   qsHigh   = _hits[i].val.qPos;
        if (dsLow    > _hits[i].val.dPos)   dsLow    = _hits[i].val.dPos;
        if (dsHigh   < _hits[i].val.dPos)   dsHigh   = _hits[i].val.dPos;
        if (minCount > _hits[i].val.uniq)   minCount = _hits[i].val.uniq;
        IL->addMer(_hits[i].val.qPos);
        ML->addMer(_hits[i].val.qPos, _hits[i].val.dPos);
        continue;
      }

      //  Doesn't look like these hits belong together.  Promote the hit
      //  to a match if it's decent.

      IL->merge();

      if ((minCount <= MINCOUNT) || (minLength <= IL->sumOfLengths())) {
        addMatch(minCount <= MINCOUNT,
                 qsLow,
                 qsHigh + config._KBmerSize,
                 dsLow,
                 dsHigh + config._KBmerSize,
                 IL,
                 ML);
        IL = new merCovering(config._KBmerSize);
        ML = new merList();
      } else {
        IL->clear();
        ML->clear();
      }

      frstDiagonal = thisDiagonalID;
      lastDiagonal = thisDiagonalID;
      qsLow        = _hits[i].val.qPos;
      qsHigh       = _hits[i].val.qPos;
      dsLow        = _hits[i].val.dPos;
      dsHigh       = _hits[i].val.dPos;
      minCount     = _hits[i].val.uniq;

      IL->addMer(_hits[i].val.qPos);
      ML->addMer(_hits[i].val.qPos, _hits[i].val.dPos);
    }

    //  Save the final cluster?

    IL->merge();

    if ((minCount <= MINCOUNT) || (minLength <= IL->sumOfLengths())) {
      addMatch(minCount <= MINCOUNT,
               qsLow,
               qsHigh + config._KBmerSize,
               dsLow,
               dsHigh + config._KBmerSize,
               IL,
               ML);
    } else {
      delete IL;
      delete ML;
    }


    //
    //  Step 2:  Merge matches into, sigh, hits, stuff them into the output
    //


    while (_matches) {

      //  Save the current match, then delete it.
      //
      unique     = _matches->_unique;
      dsLow      = _matches->_dsLo;
      dsHigh     = _matches->_dsHi;
      IL         = _matches->_IL;
      ML         = _matches->_ML;

      {
        trapMatch *n = _matches;
        _matches = _matches->_next;
        delete n;
      }

      //  Assimilate as many of the remaining matches as possible.
      //
      //  Think of this as first reversing the list, then merging as
      //  long as (dsHigh + 1000 > _matches->_dsLo).  But since we
      //  don't reverse the list, we can map:
      //    dsHigh            --> _matches->dsHi
      //    _matches->_dsLo   --> dsLow
      //  where dsHigh and dsLow are the values for the extended match.
      //
      while (_matches && (dsLow < _matches->_dsHi + 5000)) {

        //  Combine the two merCoverings
        //
        IL->merge(_matches->_IL);
        ML->merge(_matches->_ML);

        delete _matches->_IL;
        delete _matches->_ML;

        //  The start of the new match might be after the start of the
        //  merged region.  (Only rarely is it before)
        //
        //  The end of current match is always greater than the end of the
        //  new match!
        //
        if (dsLow > _matches->_dsLo)
          dsLow = _matches->_dsLo;

        unique |= _matches->_unique;

        {
          trapMatch *n = _matches;
          _matches = _matches->_next;
          delete n;
        }
      }

      if (theHitsPos >= theHitsMax) {
        theHitsMax <<= 1;
        aHit *o = 0L;
        try {
          o = new aHit [theHitsMax];
        } catch (std::bad_alloc) {
          fprintf(stderr, "hitMatrix::filter()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
          fprintf(stderr, "hitMatrix::filter()-- tried to extend output string from "uint32FMT" to "uint32FMT".\n", theHitsPos, theHitsMax);
          exit(1);
        }
        memcpy(o, theHits, theHitsPos * sizeof(aHit));
        delete [] theHits;
        theHits = o;
      }

      IL->merge();

      aHit *a = theHits + theHitsPos++;

      a->_status    = (direction == 'f');
      a->_status   |= (unique ? AHIT_HAS_UNIQUE : 0);
      a->_qsIdx     = _qsIdx;
      a->_dsIdx     = genomeMap->IIDOf(currentSeq);
      a->_dsLo      = dsLow;
      a->_dsHi      = dsHigh;
      a->_covered   = IL->sumOfLengths();
      a->_matched   = IL->numberOfPieces();  //numberOfIntervals();
      a->_numMers   = _qsMers;
      a->_ML        = ML;

      assert(a->_dsLo < a->_dsHi);

#ifdef SHOW_HITMATRIX
      _theLog->add("merged:   G::"uint32FMT"-"uint32FMT"("uint32FMT")  q:"uint32FMT" g:"uint32FMT" cov:"uint32FMT" mat:"uint32FMT" mer:"uint32FMT"\n",
              a->_dsLo, a->_dsHi, a->_dsHi - a->_dsLo,
              a->_qsIdx,
              a->_dsIdx,
              a->_covered, a->_matched, a->_numMers);
#endif

      delete IL;
    }

    //  All done with these hits.  Move to the next set.
    //
    firstHit = lastHit;
  }
}


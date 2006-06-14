#include "posix.H"
#include "snapper2.H"


//  Reports debugging information on decoding the hits
//
//#define DEBUG_HIT_DECODE

//  Reports when any chained hit is saved
//
//#define REPORT_CHAINED_HITS


hitMatrix::hitMatrix(u32bit qsLen, u32bit qsMers, u32bit qsIdx) {
  _qsLen    = qsLen;
  _qsMers   = qsMers;
  _qsIdx    = qsIdx;

  _hitsLen  = 0;
  _hitsMax  = 128;
  _hits     = new diagonalLine [_hitsMax];

  _matches  = 0L;
}


hitMatrix::~hitMatrix() {
  delete [] _hits;
}


void
hitMatrix::addMatch(u32bit         qsLo,
                    u32bit         qsHi,
                    u32bit         dsLo,
                    u32bit         dsHi,
                    merCovering   *IL,
                    merList       *ML) {
  u32bit offset = 0;

  offset = (u32bit)(config._extendWeight * qsLo);
  if (offset < config._extendMinimum)
    offset = config._extendMinimum;
  if (dsLo < offset)
    dsLo = 0;
  else
    dsLo -= offset;

  offset = (u32bit)(config._extendWeight * (_qsLen - qsHi));
  if (offset < config._extendMinimum)
    offset = config._extendMinimum;
  dsHi += offset;


  //  Create a new match
  //
  //  n = new match
  //  m = current match
  //  l = last match
  //
  trapMatch *n = new trapMatch(qsLo, qsHi, dsLo, dsHi, IL, ML);

#ifdef REPORT_CHAINED_HITS
  fprintf(stderr, "chained:  Q::"u32bitFMT"-"u32bitFMT"("u32bitFMT") G::"u32bitFMT"-"u32bitFMT"("u32bitFMT")\n",
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

#ifdef WITHOUT_DIAGONALID

inline
int
compareLines(diagonalLine *A, diagonalLine *B, u32bit qsLen) {
  u32bit a = qsLen - A->_qsPos - 1 + A->_dsPos;
  u32bit b = qsLen - B->_qsPos - 1 + B->_dsPos;

  return(((a  < b)) ||
         ((a == b) && (A->_qsPos < B->_qsPos)));
}

inline
int
compareLines(u32bit l, u32bit q, diagonalLine *B, u32bit qsLen) {
  u32bit b = qsLen - B->_qsPos - 1 + B->_dsPos;

  return(((l  < b)) ||
         ((l == b) && (q < B->_qsPos)));
}

inline
void
adjustHeap(diagonalLine *L, s32bit p, s32bit n, u32bit qsLen) {
  u32bit  q = L[p]._qsPos;
  u32bit  d = L[p]._dsPos;
  u32bit  l = qsLen - q - 1 + d;
  s32bit  c = (p << 1) + 1;  //  let c be the left child of p

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
    L[p]._qsPos      = L[c]._qsPos;
    L[p]._dsPos      = L[c]._dsPos;

    //  Move down the tree
    //
    p = c;
    c = (p << 1) + 1;
  }

  L[p]._qsPos      = q;
  L[p]._dsPos      = d;
}


#else // WITH_DIAGONALID


inline
int
compareLines(diagonalLine *A, diagonalLine *B) {
  return(((A->_diagonalID  < B->_diagonalID)) ||
         ((A->_diagonalID == B->_diagonalID) && (A->_qsPos < B->_qsPos)));
}

inline
int
compareLines(u32bit l, u32bit q, diagonalLine *B) {
  return(((l  < B->_diagonalID)) ||
         ((l == B->_diagonalID) && (q < B->_qsPos)));
}

inline
void
adjustHeap(diagonalLine *L, s32bit p, s32bit n) {
  u32bit  q = L[p]._qsPos;
  u32bit  d = L[p]._dsPos;
  u32bit  l = L[p]._diagonalID;
  s32bit  c = (p << 1) + 1;  //  let c be the left child of p

  while (c < n) {

    //  Find the larger of the two children
    //
    if ((c+1 < n) && compareLines(L+c, L+c+1))
      c++;

    //  Does the node in question fit here?
    //
    if (compareLines(l, q, L+c) == false)
      break;

    //  Else, swap the parent and the child
    //
    L[p]._qsPos      = L[c]._qsPos;
    L[p]._dsPos      = L[c]._dsPos;
    L[p]._diagonalID = L[c]._diagonalID;

    //  Move down the tree
    //
    p = c;
    c = (p << 1) + 1;
  }

  L[p]._qsPos      = q;
  L[p]._dsPos      = d;
  L[p]._diagonalID = l;
}


#endif


void
hitMatrix::filter(char      direction,
                  double    minHitCoverage,
                  u32bit    minHitLength,
                  aHit    *&theHits,
                  u32bit   &theHitsPos,
                  u32bit   &theHitsMax) {

  if (_hitsLen == 0)
    return;

  //  Decide on the minimum quality values; we pick the larger of
  //  the fixed lengths, and the sequence length * coverage.
  //
  u32bit   minLength = (u32bit)(minHitCoverage * _qsLen);
  if (minLength < minHitLength)
    minLength = minHitLength;

  //  First, sort by the dsPos.  This is done so that we can find all the hits for
  //  a specific scaffold.
  //
  sort_dsPos();

  //  Now, while there are hits left....
  //
  u32bit  firstHit   = 0;
  u32bit  lastHit    = 0;
  u32bit  currentSeq = 0;

#ifdef DEBUG_HIT_DECODE
  fprintf(stderr, "filter: got "u32bitFMT" hits\n", _hitsLen);
#endif

  while (firstHit < _hitsLen) {

    //  Move the currentSeq until the firstHit is below it.
    //
    while ((currentSeq < config._useList.numberOfSequences()) &&
           (config._useList.startOf(currentSeq) <= _hits[firstHit]._dsPos)) {
#ifdef DEBUG_HIT_DECODE_DETAILED
      fprintf(stderr, "currentSeq: "u32bitFMT" length "u64bitFMT" hit "u32bitFMT"\n",
              currentSeq,
              config._useList.startOf(currentSeq),
              _hits[firstHit]._dsPos);
#endif
      currentSeq++;
    }

    //
    //  currentSeq is now the sequence AFTER the one that we want hits in.
    //

    //  Find the first hit that is in currentSeq.  If this is the last sequence,
    //  then, of course, all remaining hits are in it.
    //
    if (currentSeq < config._useList.numberOfSequences()) {
      lastHit = firstHit + 1;
      while ((lastHit < _hitsLen) &&
             (_hits[lastHit]._dsPos < config._useList.startOf(currentSeq)))
        lastHit++;
    } else {
      lastHit = _hitsLen;
    }

    //  Drop back one sequence; this is the sequence the hits are in.
    //
    currentSeq--;

#ifdef DEBUG_HIT_DECODE
    fprintf(stderr, "Found sequence "u32bitFMT" for hits "u32bitFMT" to "u32bitFMT"\n", currentSeq, firstHit, lastHit);
#endif

    //  Adjust the hits to be relative to the start of this sequence
    //
    for (u32bit i=firstHit; i<lastHit; i++)
      _hits[i]._dsPos -= config._useList.startOf(currentSeq);

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
      for (s32bit i=(lastHit - firstHit)/2 - 1; i>=0; i--)
#ifdef WITHOUT_DIAGONALID
        adjustHeap(hitsToSort, i, lastHit - firstHit, _qsLen);
#else
        adjustHeap(hitsToSort, i, lastHit - firstHit);
#endif

      //  Sort the hits by diagonal.  This is the second part of
      //  heap sort -- Interchange the new maximum with the element
      //  at the end of the tree
      //
      for (u32bit i=lastHit - firstHit - 1; i>0; i--) {
        u32bit  q  = hitsToSort[i]._qsPos;
        u32bit  d  = hitsToSort[i]._dsPos;
#ifndef WITHOUT_DIAGONALID
        u32bit  l  = hitsToSort[i]._diagonalID;
#endif

        hitsToSort[i]._qsPos      = hitsToSort[0]._qsPos;
        hitsToSort[i]._dsPos      = hitsToSort[0]._dsPos;
#ifndef WITHOUT_DIAGONALID
        hitsToSort[i]._diagonalID = hitsToSort[0]._diagonalID;
#endif

        hitsToSort[0]._qsPos      = q;
        hitsToSort[0]._dsPos      = d;
#ifndef WITHOUT_DIAGONALID
        hitsToSort[0]._diagonalID = l;
#endif
      
#ifdef WITHOUT_DIAGONALID
        adjustHeap(hitsToSort, 0, i, _qsLen);
#else
        adjustHeap(hitsToSort, 0, i);
#endif
      }
    }

    //  Filter them
    //
#ifdef WITHOUT_DIAGONALID
    u32bit  frstDiagonal = _qsLen - _hits[firstHit]._qsPos - 1 + _hits[firstHit]._dsPos;
    u32bit  lastDiagonal = frstDiagonal;
#else
    u32bit  frstDiagonal = _hits[firstHit]._diagonalID;
    u32bit  lastDiagonal = _hits[firstHit]._diagonalID;
#endif
    u32bit  qsLow        = _hits[firstHit]._qsPos;
    u32bit  qsHigh       = _hits[firstHit]._qsPos;
    u32bit  dsLow        = _hits[firstHit]._dsPos;
    u32bit  dsHigh       = _hits[firstHit]._dsPos;

    merCovering   *IL = new merCovering(config._merSize);
    merList       *ML = new merList();

    for (u32bit i=firstHit; i<lastHit; i++) {
#ifdef WITHOUT_DIAGONALID
      u32bit thisDiagonalID = _qsLen - _hits[i]._qsPos - 1 + _hits[i]._dsPos;
#else
      u32bit thisDiagonalID = _hits[i]._diagonalID;
#endif

      //  Unconditionally extend if the diagonal difference is small.
      //
      if (lastDiagonal + config._maxDiagonal >= thisDiagonalID) {
        lastDiagonal = thisDiagonalID;
        if (qsLow  > _hits[i]._qsPos)   qsLow  = _hits[i]._qsPos;
        if (qsHigh < _hits[i]._qsPos)   qsHigh = _hits[i]._qsPos;
        if (dsLow  > _hits[i]._dsPos)   dsLow  = _hits[i]._dsPos;
        if (dsHigh < _hits[i]._dsPos)   dsHigh = _hits[i]._dsPos;
        IL->addMer(_hits[i]._qsPos);
        ML->addMer(_hits[i]._qsPos, _hits[i]._dsPos);
        continue;
      }

      //  Save the current cluster and start a new one?
      //
      u32bit qCov = IL->sumLengths();
      if (qCov >= minLength) {
        addMatch(qsLow,
                 qsHigh + config._merSize,
                 dsLow,
                 dsHigh + config._merSize,
                 IL,
                 ML);
        IL = new merCovering(config._merSize);
        ML = new merList();
      }

      if (IL)
        IL->clear();
      if (ML)
        ML->clear();

      frstDiagonal = thisDiagonalID;
      lastDiagonal = thisDiagonalID;
      qsLow        = _hits[i]._qsPos;
      qsHigh       = _hits[i]._qsPos;
      dsLow        = _hits[i]._dsPos;
      dsHigh       = _hits[i]._dsPos;

      IL->addMer(_hits[i]._qsPos);
      ML->addMer(_hits[i]._qsPos, _hits[i]._dsPos);
    }

    //  Save the final cluster?
    //
    u32bit qCov = IL->sumLengths();
    if (qCov >= minLength) {
      addMatch(qsLow,
               qsHigh + config._merSize,
               dsLow,
               dsHigh + config._merSize,
               IL,
               ML);
        IL = 0L;
        ML = 0L;
    }

    //  Delete any remaining IL
    //
    delete IL;
    delete ML;


    //  Merge and print the matches
    //
    trapMatch     *n        = 0L;

    while (_matches) {

      //  Save the current match, then delete it.
      //
      dsLow      = _matches->_dsLo;
      dsHigh     = _matches->_dsHi;
      IL         = _matches->_IL;
      ML         = _matches->_ML;

      n = _matches;
      _matches = _matches->_next;
      delete n;

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

        //  The start of the new match might be after the start of the
        //  merged region.  (Only rarely is it before)
        //
        if (dsLow > _matches->_dsLo)
          dsLow = _matches->_dsLo;

        //  The end of current match is always greater than the end of the
        //  new match!
        //
        n = _matches;
        _matches = _matches->_next;
        delete n->_IL;
        delete n->_ML;
        delete n;
      }

      if (theHitsPos >= theHitsMax) {
        theHitsMax <<= 1;
        aHit *o = 0L;
        try {
          o = new aHit [theHitsMax];
        } catch (std::bad_alloc) {
          fprintf(stderr, "hitMatrix::filter()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
          fprintf(stderr, "hitMatrix::filter()-- tried to extend output string from "u32bitFMT" to "u32bitFMT".\n", theHitsPos, theHitsMax);
          exit(1);
        }
        memcpy(o, theHits, theHitsPos * sizeof(aHit));
        delete [] theHits;
        theHits = o;
      }

      aHit *a = theHits + theHitsPos++;

      a->_status    = (direction == 'f');
      a->_qsIdx     = _qsIdx;
      a->_dsIdx     = config._useList.IIDOf(currentSeq);
      a->_dsLo      = dsLow;
      a->_dsHi      = dsHigh;
      a->_covered   = IL->sumLengths();
      a->_matched   = IL->numberOfPieces();
      a->_numMers   = _qsMers;
      a->_ML        = ML;

#ifdef REPORT_CHAINED_HITS
      fprintf(stderr, "merged:   G::"u32bitFMT"-"u32bitFMT"("u32bitFMT")  q:"u32bitFMT" g:"u32bitFMT" cov:"u32bitFMT" mat:"u32bitFMT" mer:"u32bitFMT"\n",
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


#include "posix.H"
#include "snapper2.H"
#include "intervalList.H"

#ifdef TRUE64BIT
#define HITOUTPUTLINE "-%c -e %u -D %u %u %u -M %u %u %u\n"
#else
#define HITOUTPUTLINE "-%c -e %lu -D %lu %lu %lu -M %lu %lu %lu\n"
#endif

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
                    intervalList  *IL) {
  u32bit offset = 0;

  offset = config._extendWeight * qsLo;
  if (offset < config._extendMinimum)
    offset = config._extendMinimum;
  if (dsLo < offset)
    dsLo = 0;
  else
    dsLo -= offset;

  offset = config._extendWeight * (_qsLen - qsHi);
  if (offset < config._extendMinimum)
    offset = config._extendMinimum;
  dsHi += offset;

  //  Create a new match
  //
  //  n = new match
  //  m = current match
  //  l = last match
  //
  trapMatch *n = new trapMatch(qsLo, qsHi, dsLo, dsHi, IL);


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
//
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

void
hitMatrix::filter(char      direction,
                  aHit    *&theOutput,
                  u32bit   &theOutputPos,
                  u32bit   &theOutputMax) {

  if (_hitsLen == 0)
    return;

  //  Decide on the minimum quality values; we pick the larger of
  //  the fixed lengths, and the sequence length * coverage.
  //
  u32bit   minLength = (u32bit)(config._minHitCoverage * _qsLen);
  if (minLength < config._minHitLength)
    minLength = config._minHitLength;

  //  First, sort by the dsPos.  This is done so that we can find all the hits for
  //  a specific scaffold.
  //
  sort_dsPos();

  //  Now, while there are hits left....
  //
  u32bit  firstHit   = 0;
  u32bit  lastHit    = 0;
  u32bit  currentSeq = 0;

  while (firstHit < _hitsLen) {

    //  Move the currentSeq until the firstHit is below it.
    //
    while ((currentSeq < config._useListLen) &&
           (config._useList[currentSeq].start <= _hits[firstHit]._dsPos))
      currentSeq++;

    //
    //  currentSeq is now the sequence AFTER the one that we want hits in.
    //

    //  Find the first hit that is in currentSeq.  If this is the last sequence,
    //  then, of course, all remaining hits are in it.
    //
    if (currentSeq < config._useListLen) {
      lastHit = firstHit + 1;
      while ((lastHit < _hitsLen) &&
             (_hits[lastHit]._dsPos < config._useList[currentSeq].start))
        lastHit++;
    } else {
      lastHit = _hitsLen;
    }

    //  Drop back one sequence; this is the sequence the hits are in.
    //
    currentSeq--;


    //  Adjust the hits to be relative to the start of this sequence
    //
    for (u32bit i=firstHit; i<lastHit; i++)
      _hits[i]._dsPos -= config._useList[currentSeq].start;

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
        adjustHeap(hitsToSort, i, lastHit - firstHit);

      //  Sort the hits be diagonal.  This is the second part of
      //  heap sort -- Interchange the new maximum with the element
      //  at the end of the tree
      //
      for (u32bit i=lastHit - firstHit - 1; i>0; i--) {
        u32bit  q  = hitsToSort[i]._qsPos;
        u32bit  d  = hitsToSort[i]._dsPos;
        u32bit  l  = hitsToSort[i]._diagonalID;
        
        hitsToSort[i]._qsPos      = hitsToSort[0]._qsPos;
        hitsToSort[i]._dsPos      = hitsToSort[0]._dsPos;
        hitsToSort[i]._diagonalID = hitsToSort[0]._diagonalID;

        hitsToSort[0]._qsPos      = q;
        hitsToSort[0]._dsPos      = d;
        hitsToSort[0]._diagonalID = l;
      
        adjustHeap(hitsToSort, 0, i);
      }
    }

    //  Filter them
    //
    u32bit  frstDiagonal = _hits[firstHit]._diagonalID;
    u32bit  lastDiagonal = _hits[firstHit]._diagonalID;
    u32bit  qsLow        = _hits[firstHit]._qsPos;
    u32bit  qsHigh       = _hits[firstHit]._qsPos;
    u32bit  dsLow        = _hits[firstHit]._dsPos;
    u32bit  dsHigh       = _hits[firstHit]._dsPos;

    //  Create a new interval list, and space to count the number of mers in a match
    //
    intervalList   *IL = new intervalList(config._merSize);

    for (u32bit i=firstHit; i<lastHit; i++) {

      //  Unconditionally extend if the diagonal difference is small.
      //
      if (lastDiagonal + config._maxDiagonal >= _hits[i]._diagonalID) {
        lastDiagonal = _hits[i]._diagonalID;
        if (qsLow  > _hits[i]._qsPos)   qsLow  = _hits[i]._qsPos;
        if (qsHigh < _hits[i]._qsPos)   qsHigh = _hits[i]._qsPos;
        if (dsLow  > _hits[i]._dsPos)   dsLow  = _hits[i]._dsPos;
        if (dsHigh < _hits[i]._dsPos)   dsHigh = _hits[i]._dsPos;
        IL->addInterval(_hits[i]._qsPos);
        continue;
      }

      //  Save the current cluster and start a new one?
      //
      u32bit qCov = IL->sumIntervalLengths();
      if (qCov >= minLength) {
        addMatch(qsLow,
                 qsHigh + config._merSize,
                 dsLow,
                 dsHigh + config._merSize,
                 IL);
        IL = new intervalList(config._merSize);
      }

      if (IL)
        IL->clear();

      frstDiagonal = _hits[i]._diagonalID;
      lastDiagonal = _hits[i]._diagonalID;
      qsLow        = _hits[i]._qsPos;
      qsHigh       = _hits[i]._qsPos;
      dsLow        = _hits[i]._dsPos;
      dsHigh       = _hits[i]._dsPos;

      IL->addInterval(_hits[i]._qsPos);
    }

    //  Save the final cluster?
    //
    u32bit qCov = IL->sumIntervalLengths();
    if (qCov >= minLength) {
      addMatch(qsLow,
               qsHigh + config._merSize,
               dsLow,
               dsHigh + config._merSize,
               IL);
        IL = 0;
    }

    //  Delete any remaining IL
    //
    delete IL;


    //  Merge and print the matches
    //
    trapMatch     *n  = 0L;
    u32bit         ML = 0;

    while (_matches) {

      //  Save the current match, then delete it.
      //
      dsLow      = _matches->_dsLo;
      dsHigh     = _matches->_dsHi;
      IL         = _matches->_IL;
      ML         = IL->sumIntervalLengths();

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

        //  Combine the two interval lists
        //
        IL->merge(_matches->_IL);
        ML += _matches->_IL->sumIntervalLengths();

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
        delete n;
      }

      if (theOutputPos >= theOutputMax) {
        theOutputMax <<= 1;
        aHit *o = 0L;
        try {
          o = new aHit [theOutputMax];
        } catch (std::bad_alloc) {
          fprintf(stderr, "hitMatrix::filter()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
          fprintf(stderr, "hitMatrix::filter()-- tried to extend output string from %lu to %lu.\n", theOutputPos, theOutputMax);
          exit(1);
        }
        memcpy(o, theOutput, theOutputPos * sizeof(aHit));
        delete [] theOutput;
        theOutput = o;
      }

      aHit *a = theOutput + theOutputPos++;

      a->_status    = (direction == 'f');
      a->_qsIdx     = _qsIdx;
      a->_dsIdx     = config._useList[currentSeq].seq;
      a->_dsLo      = dsLow;
      a->_dsHi      = dsHigh;
      a->_covered   = IL->sumIntervalLengths();
      a->_matched   = ML;
      a->_numMers   = _qsMers;

      delete IL;
    }

    //  All done with these hits.  Move to the next set.
    //
    firstHit = lastHit;
  }
}


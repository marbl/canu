#include "posix.H"
#include "searchGENOME.H"
#include "aHit.H"

#define TRACE    0

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
                    merCovering   *IL) {
  u32bit offset = 0;

  //  Extend the match
  //
  //  Two methods: the first uses hardcoded parameters, has two
  //  plateau's, and is the one to use for ESTs and mRNA in ESTmapper.
  //  The second is paramterized, and has a single plateau.
  //

  if (config._extendAlternate) {
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
  } else {
    //  If the start of the match is near the start of the EST, we do
    //  not need to search very far in the genome.
    //
    offset = 0;
    if (qsLo < 50)
      offset = 2000;
    else
      if (qsLo < 100)
        offset = 5000;
      else
        offset = 50 * qsLo;

    if (dsLo < offset)
      dsLo = 0;
    else
      dsLo -= offset;

    //  Likewise, if the match is near the end of the EST, extend.  We
    //  don't know the length of the genomic sequence, so we can't check
    //  for "overflow".
    //
    offset = _qsLen - qsHi;
    if (offset < 50)
      dsHi += 2000;
    else
      if (offset < 100)
        dsHi += 5000;
      else
        dsHi += 50 * offset;
  }



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
hitMatrix::filter(encodedQuery *query, bool isReverse) {

  if (_hitsLen == 0)
    return;

  //  Decide on the minimum quality values; we pick the larger of
  //  the fixed lengths, and the sequence length * coverage.
  //
  u32bit   minLengthSingle   = (u32bit)(config._minCoverageSingle   * _qsLen);
  u32bit   minLengthMultiple = (u32bit)(config._minCoverageMultiple * _qsLen);

  if (minLengthSingle < config._minLengthSingle)
    minLengthSingle = config._minLengthSingle;

  if (minLengthMultiple < config._minLengthMultiple)
    minLengthMultiple = config._minLengthMultiple;



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
    while ((currentSeq < config._useList.numberOfSequences()) &&
           (config._useList.startOf(currentSeq) <= _hits[firstHit]._dsPos))
      currentSeq++;

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

#if TRACE
    fprintf(stdout, "Hits are in sequence %d\n", config._useList.IIDOf(currentSeq));
    fprintf(stdout, "filtering %u hits -- first = %u last = %u.\n", _hitsLen, firstHit, lastHit);

#if 0
    fprintf(stdout, "UNSORTED\n");
    for (u32bit i=firstHit; i<lastHit; i++)
      fprintf(stdout, "hit at qs=%4u ds=%6u diag=%6u\n",
              _hits[i]._qsPos,
              _hits[i]._dsPos,
              _hits[i]._diagonalID);
#endif
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

      //  Sort the hits be diagonal.  This is the second part of
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


    //  Check the sorting
    //
#if 0
#if 0
    fprintf(stderr, "sort by diagonal:\n");
    for (u32bit i=firstHit; i<lastHit; i++)
      fprintf(stderr, "%8u %8u %8u\n", _hits[i]._diagonalID, _hits[i]._qsPos, _hits[i]._dsPos);
#endif
    for (u32bit i=firstHit; i<lastHit-1; i++) {
      if (_hits[i]._diagonalID > _hits[i+1]._diagonalID) {
        fprintf(stderr, "sort by diagonal failed.\n");
        exit(1);
      }
    }
#endif



#if TRACE
#if 0
    fprintf(stdout, "SORTED\n");
    for (u32bit i=firstHit; i<lastHit; i++)
      fprintf(stdout, "hit at qs=%4u ds=%6u diag=%6u\n",
              _hits[i]._qsPos,
              _hits[i]._dsPos,
              _hits[i]._diagonalID);
#endif

    fprintf(stdout, "FILTERED\n");
#endif

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

    //  Create a new merCovering, and space to count the number of mers in a match
    //
    merCovering *IL = new merCovering(config._merSize);

    for (u32bit i=firstHit; i<lastHit; i++) {
#ifdef WITHOUT_DIAGONALID
      u32bit thisDiagonalID = _qsLen - _hits[i]._qsPos - 1 + _hits[i]._dsPos;
#else
      u32bit thisDiagonalID = _hits[i]._diagonalID;
#endif



#if TRACE
      fprintf(stdout, "hit[qs=%6u ds=%7u d=%7u]  box[qs=%6u-%6u ds=%7u-%7u d=%7u-%7u]  ",
              _hits[i]._qsPos,
              _hits[i]._dsPos,
              thisDiagonalID,
              qsLow, qsHigh, dsLow, dsHigh, frstDiagonal, lastDiagonal);
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
#if TRACE
        fprintf(stdout, "extend qs=%9u-%9u ds=%9u-%9u  diag=%9u-%9u (diagonal)\n",
                qsLow, qsHigh, dsLow, dsHigh, frstDiagonal, lastDiagonal);
#endif
        continue;
      }


      //  XXX:  Prototype for extending only if the next hit is near
      //  the last hit.
      //
      if (((dsHigh <= _hits[i]._dsPos) && (_hits[i]._dsPos - dsHigh <= config._maxIntronLength)) ||
          ((dsHigh >= _hits[i]._dsPos) && (dsHigh - _hits[i]._dsPos <= config._maxIntronLength))) {

        //  Extend into multiple-exon like things only if the input
        //  sequence is long.
        //
        if (_qsLen > config._smallSequenceCutoff) {

          //  Extend if the qsOverlap is small (or nonexistant)
          //
          if ((qsHigh + config._merSize) < (_hits[i]._qsPos + config._qsOverlap)) {
            lastDiagonal = thisDiagonalID;
            if (qsLow  > _hits[i]._qsPos)   qsLow  = _hits[i]._qsPos;
            if (qsHigh < _hits[i]._qsPos)   qsHigh = _hits[i]._qsPos;
            if (dsLow  > _hits[i]._dsPos)   dsLow  = _hits[i]._dsPos;
            if (dsHigh < _hits[i]._dsPos)   dsHigh = _hits[i]._dsPos;
            IL->addMer(_hits[i]._qsPos);
#if TRACE
            fprintf(stdout, "extend qs=%9u-%9u ds=%9u-%9u  diag=%9u-%9u (qsOverlap)\n",
                    qsLow, qsHigh, dsLow, dsHigh, frstDiagonal, lastDiagonal);
#endif
            continue;
          }

          //  Extend if the dsOverlap is small (or nonexistant)
          //
          if (_hits[i]._dsPos < (dsLow + config._dsOverlap)) {
            lastDiagonal = thisDiagonalID;
            if (qsLow  > _hits[i]._qsPos)   qsLow  = _hits[i]._qsPos;
            if (qsHigh < _hits[i]._qsPos)   qsHigh = _hits[i]._qsPos;
            if (dsLow  > _hits[i]._dsPos)   dsLow  = _hits[i]._dsPos;
            if (dsHigh < _hits[i]._dsPos)   dsHigh = _hits[i]._dsPos;
            IL->addMer(_hits[i]._qsPos);
#if TRACE
            fprintf(stdout, "extend qs=%9u-%9u ds=%9u-%9u  diag=%9u-%9u (dsOverlap)\n",
                    qsLow, qsHigh, dsLow, dsHigh, frstDiagonal, lastDiagonal);
#endif
            continue;
          }
        }
      }  //  XXX:  End prototype

#if TRACE
      fprintf(stdout, "close current cluster.\nGOOD?  qsCov=%u; >= %u or %u?  diag: %u < 25?\n",
              qsHigh - qsLow,
              minLengthSingle,
              minLengthMultiple,
              lastDiagonal - frstDiagonal);
#endif


      //  Save the current cluster and start a new one?
      //
      u32bit qCov = IL->sumOfLengths();
      if ((qCov >= minLengthMultiple) ||
          ((lastDiagonal - frstDiagonal < 25) && (qCov >= minLengthSingle))) {
#if TRACE
        fprintf(stdout, "add match!\n");
#endif
        addMatch(qsLow,
                 qsHigh + config._merSize,
                 dsLow,
                 dsHigh + config._merSize,
                 IL);
        IL = new merCovering(config._merSize);
      }

      if (IL)
        IL->clear();

#if TRACE
      fprintf(stdout, "reset!\n");
#endif

      frstDiagonal = thisDiagonalID;
      lastDiagonal = thisDiagonalID;
      qsLow        = _hits[i]._qsPos;
      qsHigh       = _hits[i]._qsPos;
      dsLow        = _hits[i]._dsPos;
      dsHigh       = _hits[i]._dsPos;

#if TRACE
      fprintf(stdout, "hit[qs=%6u ds=%7u d=%7u]  box[qs=%6u-%6u ds=%7u-%7u d=%7u-%7u]  (initial hit)\n",
              _hits[i]._qsPos,
              _hits[i]._dsPos,
              _qsLen - _hits[i]._qsPos - 1 + _hits[i]._dsPos,
              qsLow, qsHigh, dsLow, dsHigh, frstDiagonal, lastDiagonal);
#endif

      IL->addMer(_hits[i]._qsPos);
    }

    //  Save the final cluster?
    //
    u32bit qCov = IL->sumOfLengths();
    if ((qCov >= minLengthMultiple) ||
        ((lastDiagonal - frstDiagonal < 21) && (qCov >= minLengthSingle))) {
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
      ML         = IL->sumOfLengths();

      n = _matches;
      _matches = _matches->_next;
      delete n;

#if TRACE
      fprintf(stdout, "Merge: %8u %8u\n", dsLow, dsHigh);
#endif

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
        ML += _matches->_IL->sumOfLengths();

        //  The start of the new match might be after the start of the
        //  merged region.  (Only rarely is it before)
        //
        if (dsLow > _matches->_dsLo)
          dsLow = _matches->_dsLo;

        //  The end of current match is always greater than the end of the
        //  new match!
        //
        //dsHigh = _matches->_dsHi;

#if TRACE
        fprintf(stdout, "Merge: %8u %8u -> %8u %8u\n", _matches->_dsLo, _matches->_dsHi, dsLow, dsHigh);
#endif

        n = _matches;
        _matches = _matches->_next;
        delete n->_IL;
        delete n;
      }


      if (config._binaryOutput) {
        aHit a;
      
        a._forward   = !isReverse;
        a._merged    = false;
        a._qsIdx     = _qsIdx;
        a._dsIdx     = config._useList.IIDOf(currentSeq);
        a._dsLo      = dsLow;
        a._dsHi      = dsHigh;
        a._covered   = IL->sumOfLengths();
        a._matched   = ML;
        a._numMers   = _qsMers;

        query->addOutput(&a, sizeof(aHit));
      } else {
        char  line[128];

        sprintf(line, "-%c -e "u32bitFMT" -D "u64bitFMT" "u32bitFMT" "u32bitFMT" -M "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
                isReverse ? 'r' : 'f', _qsIdx,
                config._useList.IIDOf(currentSeq),
                dsLow, dsHigh, IL->sumOfLengths(), ML, _qsMers);

        query->addOutput(line, 0);
      }

      delete IL;
    }

    //  All done with these hits.  Move to the next set.
    //
    firstHit = lastHit;
  }
}


#include "posix.H"
#include "searchGENOME.H"
#include "aHit.H"

#define TRACE    0

hitMatrix::hitMatrix(u32bit qsLen, u32bit qsMers, u32bit qsIdx, bool reversed) {
  _qsLen    = qsLen;
  _qsMers   = qsMers;
  _qsIdx    = qsIdx;

  _hitsLen  = 0;
  _hitsMax  = 128;
  _hits     = new diagonalLine [_hitsMax];

  //_reversed = reversed;

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
hitMatrix::filter(char direction, char *&theOutput, u32bit &theOutputPos, u32bit &theOutputMax) {

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
    u32bit  frstDiagonal = _hits[firstHit]._diagonalID;
    u32bit  lastDiagonal = _hits[firstHit]._diagonalID;
    u32bit  qsLow        = _hits[firstHit]._qsPos;
    u32bit  qsHigh       = _hits[firstHit]._qsPos;
    u32bit  dsLow        = _hits[firstHit]._dsPos;
    u32bit  dsHigh       = _hits[firstHit]._dsPos;

    //  Create a new merCovering, and space to count the number of mers in a match
    //
    merCovering *IL = new merCovering(config._merSize);

    for (u32bit i=firstHit; i<lastHit; i++) {

#if TRACE
      fprintf(stdout, "hit[qs=%6u ds=%7u d=%7u]  box[qs=%6u-%6u ds=%7u-%7u d=%7u-%7u]  ",
              _hits[i]._qsPos,
              _hits[i]._dsPos,
              _hits[i]._diagonalID,
              qsLow, qsHigh, dsLow, dsHigh, frstDiagonal, lastDiagonal);
#endif

      //  Unconditionally extend if the diagonal difference is small.
      //
      if (lastDiagonal + config._maxDiagonal >= _hits[i]._diagonalID) {
        lastDiagonal = _hits[i]._diagonalID;
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
          lastDiagonal = _hits[i]._diagonalID;
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
          lastDiagonal = _hits[i]._diagonalID;
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

      //  XXX:  End prototype
      }

#if TRACE
      fprintf(stdout, "close current cluster.\nGOOD?  qsCov=%u; >= %u or %u?  diag: %u < 25?\n",
              qsHigh - qsLow,
              minLengthSingle,
              minLengthMultiple,
              lastDiagonal - frstDiagonal);
#endif


      //  Save the current cluster and start a new one?
      //
      u32bit qCov = IL->sumLengths();
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

      frstDiagonal = _hits[i]._diagonalID;
      lastDiagonal = _hits[i]._diagonalID;
      qsLow        = _hits[i]._qsPos;
      qsHigh       = _hits[i]._qsPos;
      dsLow        = _hits[i]._dsPos;
      dsHigh       = _hits[i]._dsPos;

#if TRACE
      fprintf(stdout, "hit[qs=%6u ds=%7u d=%7u]  box[qs=%6u-%6u ds=%7u-%7u d=%7u-%7u]  (initial hit)\n",
              _hits[i]._qsPos,
              _hits[i]._dsPos,
              _hits[i]._diagonalID,
              qsLow, qsHigh, dsLow, dsHigh, frstDiagonal, lastDiagonal);
#endif

      IL->addMer(_hits[i]._qsPos);
    }

    //  Save the final cluster?
    //
    u32bit qCov = IL->sumLengths();
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
      ML         = IL->sumLengths();

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
        ML += _matches->_IL->sumLengths();

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

      if (theOutputPos + 128 >= theOutputMax) {
        theOutputMax <<= 1;
        char *o = 0L;
        try {
          o = new char [theOutputMax];
        } catch (std::bad_alloc) {
          fprintf(stderr, "hitMatrix::filter()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
          fprintf(stderr, "hitMatrix::filter()-- tried to extend output string from %lu to %lu bytes.\n", theOutputPos, theOutputMax);
          exit(1);
        }
        memcpy(o, theOutput, theOutputPos);
        delete [] theOutput;
        theOutput = o;
      }


      if (config._binaryOutput) {
        //  This little bit of nastiness is to make the compiler stop warning about
        //  alignment when we just cast a char pointer + offset to u32bit:
        //
        //  u32bit *O = (u32bit *)(theOutput + theOutputPos);
        //
        //  hitMatrix.C:534: warning: cast from `char*' to `u32bit*' increases required 
        //  alignment of target type
        //
        //  Plus, it lets us define a type for hits.
        //
        aHit a;
      
        a._forward   = (direction == 'f');
        a._merged    = false;
        a._qsIdx     = _qsIdx;
        a._dsIdx     = config._useList.IIDOf(currentSeq);
        a._dsLo      = dsLow;
        a._dsHi      = dsHigh;
        a._covered   = IL->sumLengths();
        a._matched   = ML;
        a._numMers   = _qsMers;

        memcpy(theOutput + theOutputPos, &a, sizeof(aHit));

        theOutputPos += (u32bit)sizeof(aHit);
      } else {

#ifdef TRUE64BIT
#define HITOUTPUTLINE "-%c -e %u -D %u %u %u -M %u %u %u\n"
#else
#define HITOUTPUTLINE 
#endif



        sprintf(theOutput + theOutputPos, "-%c -e "u32bitFMT" -D "u64bitFMT" "u32bitFMT" "u32bitFMT" -M "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
                direction, _qsIdx,
                config._useList.IIDOf(currentSeq),
                dsLow, dsHigh, IL->sumLengths(), ML, _qsMers);

        while (theOutput[theOutputPos])
          theOutputPos++;
      }

      delete IL;

      //pthread_mutex_lock(&queryMatchMutex);
      queryMatchCounts[_qsIdx]++;
      //pthread_mutex_unlock(&queryMatchMutex);
    }

    //  All done with these hits.  Move to the next set.
    //
    firstHit = lastHit;
  }
}


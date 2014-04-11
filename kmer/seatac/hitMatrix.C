#include "seatac.H"


hitMatrix::hitMatrix(uint32 qsLen, uint32 qsIdx) {
  _qsLen    = qsLen;
  _qsIdx    = qsIdx;

  //  Because this is doing scaffolds or chromosomes against more than
  //  1/4 a genome, we expect a LOT of hits.  Start off with a good
  //  amount of memory.
  //
  //  At 8 bytes per diagonalLine, 128M of these is 1GB.  Which works
  //  great for aligning mamalian chromosomes and stinks for microbes.
  //
  _hitsLen  = 0;
  _hitsMax  = 32 * 1024 * 1024;
  _hits     = new diagonalLine [_hitsMax];
}


hitMatrix::~hitMatrix() {
  delete [] _hits;
}


//  Utility for sorting the diagonal lines in the hitMatrix
//
//  The two comparison functions return true if the first line
//  is less than the second line.

#ifdef WITHOUT_DIAGONALID

inline
int
compareLines(diagonalLine *A, diagonalLine *B, uint32 qsLen) {
  uint32 a = qsLen - A->_qsPos - 1 + A->_dsPos;
  uint32 b = qsLen - B->_qsPos - 1 + B->_dsPos;

  return(((a  < b)) ||
         ((a == b) && (A->_qsPos < B->_qsPos)));
}

inline
int
compareLines(uint32 l, uint32 q, diagonalLine *B, uint32 qsLen) {
  uint32 b = qsLen - B->_qsPos - 1 + B->_dsPos;

  return(((l  < b)) ||
         ((l == b) && (q < B->_qsPos)));
}

inline
void
adjustHeap(diagonalLine *L, int32 p, int32 n, uint32 qsLen) {
  uint32  q = L[p]._qsPos;
  uint32  d = L[p]._dsPos;
  uint32  l = qsLen - q - 1 + d;
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
compareLines(uint32 l, uint32 q, diagonalLine *B) {
  return(((l  < B->_diagonalID)) ||
         ((l == B->_diagonalID) && (q < B->_qsPos)));
}

inline
void
adjustHeap(diagonalLine *L, int32 p, int32 n) {
  uint32  q = L[p]._qsPos;
  uint32  d = L[p]._dsPos;
  uint32  l = L[p]._diagonalID;
  int32  c = (p << 1) + 1;  //  let c be the left child of p

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
hitMatrix::processMatrix(char direction, filterObj *FO) {

  if (_hitsLen == 0)
    return;

  //  First, sort by the dsPos.  This is done so that we can find all the hits for
  //  a specific scaffold.
  //
  sort_dsPos();


  merCovering    IL(config._merSize);
  uint32         ILlength = 0;

  //  Now, while there are hits left....
  //
  uint32  firstHit   = 0;
  uint32  lastHit    = 0;
  uint32  currentSeq = 0;

  while (firstHit < _hitsLen) {

    //  Move the currentSeq until the firstHit is below it.
    //
    while ((currentSeq < config._genome->numberOfSequences()) &&
           (config._genome->startOf(currentSeq) <= _hits[firstHit]._dsPos))
      currentSeq++;

    //
    //  currentSeq is now the sequence AFTER the one that we want hits in.
    //

    //  Find the first hit that is in currentSeq.  If this is the last sequence,
    //  then, of course, all remaining hits are in it.
    //
    if (currentSeq < config._genome->numberOfSequences()) {
      lastHit = firstHit + 1;
      while ((lastHit < _hitsLen) &&
             (_hits[lastHit]._dsPos < config._genome->startOf(currentSeq)))
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
      _hits[i]._dsPos -= config._genome->startOf(currentSeq);

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
#ifdef WITHOUT_DIAGONALID
        adjustHeap(hitsToSort, i, lastHit - firstHit, _qsLen);
#else
        adjustHeap(hitsToSort, i, lastHit - firstHit);
#endif

      //  Sort the hits be diagonal.  This is the second part of
      //  heap sort -- Interchange the new maximum with the element
      //  at the end of the tree
      //
      for (uint32 i=lastHit - firstHit - 1; i>0; i--) {
        uint32  q  = hitsToSort[i]._qsPos;
        uint32  d  = hitsToSort[i]._dsPos;
#ifndef WITHOUT_DIAGONALID
        uint32  l  = hitsToSort[i]._diagonalID;
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
    uint32  lastDiagonal = _qsLen - _hits[firstHit]._qsPos - 1 + _hits[firstHit]._dsPos;
#else
    uint32  lastDiagonal = _hits[firstHit]._diagonalID;
#endif
    uint32  qsLow        = _hits[firstHit]._qsPos;
    uint32  qsHigh       = _hits[firstHit]._qsPos;
    uint32  dsLow        = _hits[firstHit]._dsPos;
    uint32  dsHigh       = _hits[firstHit]._dsPos;

    IL.clear();

    for (uint32 i=firstHit; i<lastHit; i++) {
      //fprintf(stdout, "hit[%6u] seq=%8u qs=%5u ds=%5u\n", i, currentSeq, _hits[i]._qsPos, _hits[i]._dsPos);

      //
      //  Extend if on the same diagonal, and consecutive sequence.
      //
      if ((lastDiagonal ==
#ifdef WITHOUT_DIAGONALID
           (_qsLen - _hits[i]._qsPos - 1 + _hits[i]._dsPos)
#else
           _hits[i]._diagonalID
#endif
           ) &&
          (qsLow <= _hits[i]._qsPos) &&
          (_hits[i]._qsPos <= qsHigh + config._merSize + config._maxGap)) {
        if (qsLow  > _hits[i]._qsPos)   qsLow  = _hits[i]._qsPos;
        if (qsHigh < _hits[i]._qsPos)   qsHigh = _hits[i]._qsPos;
        if (dsLow  > _hits[i]._dsPos)   dsLow  = _hits[i]._dsPos;
        if (dsHigh < _hits[i]._dsPos)   dsHigh = _hits[i]._dsPos;
        IL.addMer(_hits[i]._qsPos);
      } else {

        //
        //  Save the match.  cut-n-paste with below.
        //

        ILlength = IL.sumOfLengths();
        IL.clear();

        if (ILlength >= config._minLength) {
          if (direction == 'r') {
            FO->addHit(direction,
                       config._genome->IIDOf(currentSeq),
                       dsLow,
                       dsHigh - dsLow + config._merSize,
                       _qsIdx,
                       _qsLen - qsHigh - config._merSize,
                       qsHigh - qsLow + config._merSize,
                       ILlength);
          } else {
            FO->addHit(direction,
                       config._genome->IIDOf(currentSeq),
                       dsLow,
                       dsHigh - dsLow + config._merSize,
                       _qsIdx,
                       qsLow,
                       qsHigh - qsLow + config._merSize,
                       ILlength);
          }
        }

#ifdef WITHOUT_DIAGONALID
        lastDiagonal = _qsLen - _hits[i]._qsPos - 1 + _hits[i]._dsPos;
#else
        lastDiagonal = _hits[i]._diagonalID;
#endif
        qsLow        = _hits[i]._qsPos;
        qsHigh       = _hits[i]._qsPos;
        dsLow        = _hits[i]._dsPos;
        dsHigh       = _hits[i]._dsPos;
        IL.addMer(_hits[i]._qsPos);
      }
    }

    //  Save the final cluster?  (cut-n-paste from above)
    //
    ILlength = IL.sumOfLengths();
    IL.clear();

    if (ILlength >= config._minLength) {
      if (direction == 'r') {
        FO->addHit(direction,
                   config._genome->IIDOf(currentSeq),
                   dsLow,
                   dsHigh - dsLow + config._merSize,
                   _qsIdx,
                   _qsLen - qsHigh - config._merSize,
                   qsHigh - qsLow + config._merSize,
                   ILlength);
      } else {
        FO->addHit(direction,
                   config._genome->IIDOf(currentSeq),
                   dsLow,
                   dsHigh - dsLow + config._merSize,
                   _qsIdx,
                   qsLow,
                   qsHigh - qsLow + config._merSize,
                   ILlength);
      }
    }

    //  All done with these hits.  Move to the next set.
    //
    firstHit = lastHit;
  }
}

  

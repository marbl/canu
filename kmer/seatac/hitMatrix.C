#include "posix.H"
#include "seatac.H"


hitMatrix::hitMatrix(u32bit qsLen, u32bit qsIdx, bool reversed) {
  _qsLen    = qsLen;
  _qsIdx    = qsIdx;

  //  Because this is doing scaffolds or chromosomes against more than
  //  1/4 a genome, we expect a LOT of hits.  Start off with a good
  //  amount of memory.
  //
  _hitsLen  = 0;
  _hitsMax  = 16 * 1024 * 1024;
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
hitMatrix::processMatrix(char direction, filterObj *FO) {

  if (_hitsLen == 0)
    return;

  //  First, sort by the dsPos.  This is done so that we can find all the hits for
  //  a specific scaffold.
  //
  sort_dsPos();


  merCovering    IL(config._merSize);
  u32bit         ILlength = 0;

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


    //  Filter them
    //
#ifdef WITHOUT_DIAGONALID
    u32bit  lastDiagonal = _qsLen - _hits[firstHit]._qsPos - 1 + _hits[firstHit]._dsPos;
#else
    u32bit  lastDiagonal = _hits[firstHit]._diagonalID;
#endif
    u32bit  qsLow        = _hits[firstHit]._qsPos;
    u32bit  qsHigh       = _hits[firstHit]._qsPos;
    u32bit  dsLow        = _hits[firstHit]._dsPos;
    u32bit  dsHigh       = _hits[firstHit]._dsPos;

    IL.clear();

    for (u32bit i=firstHit; i<lastHit; i++) {

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

        ILlength = IL.sumLengths();
        IL.clear();

        if (ILlength >= config._minLength) {
          if (direction == 'r') {
            FO->addHit(direction,
                       _qsIdx,
                       _qsLen - qsHigh - config._merSize,
                       qsHigh - qsLow + config._merSize,
                       config._useList[currentSeq].seq,
                       dsLow,
                       dsHigh - dsLow + config._merSize,
                       ILlength);
          } else {
            FO->addHit(direction,
                       _qsIdx,
                       qsLow,
                       qsHigh - qsLow + config._merSize,
                       config._useList[currentSeq].seq,
                       dsLow,
                       dsHigh - dsLow + config._merSize,
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
    ILlength = IL.sumLengths();
    IL.clear();

    if (ILlength >= config._minLength) {
      if (direction == 'r') {
        FO->addHit(direction,
                   _qsIdx,
                   _qsLen - qsHigh - config._merSize,
                   qsHigh - qsLow + config._merSize,
                   config._useList[currentSeq].seq,
                   dsLow,
                   dsHigh - dsLow + config._merSize,
                   ILlength);
      } else {
        FO->addHit(direction,
                   _qsIdx,
                   qsLow,
                   qsHigh - qsLow + config._merSize,
                   config._useList[currentSeq].seq,
                   dsLow,
                   dsHigh - dsLow + config._merSize,
                   ILlength);
      }
    }

    //  All done with these hits.  Move to the next set.
    //
    firstHit = lastHit;
  }
}

  

#include "posix.H"
#include "searchGENOME.H"
#include "intervalList.H"
#include "aHit.H"

//  $Id$

#ifdef TRUE64BIT
#define HITOUTPUTLINE  "-%c -e %u %u %u -D %u %u %u\n"
#else
#define HITOUTPUTLINE  "-%c -e %lu %lu %lu -D %lu %lu %lu\n"
#endif

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
                    intervalList  *IL,
                    char           direction) {
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
    u32bit  lastDiagonal = _hits[firstHit]._diagonalID;
    u32bit  qsLow        = _hits[firstHit]._qsPos;
    u32bit  qsHigh       = _hits[firstHit]._qsPos;
    u32bit  dsLow        = _hits[firstHit]._dsPos;
    u32bit  dsHigh       = _hits[firstHit]._dsPos;

    for (u32bit i=firstHit; i<lastHit; i++) {


      //  Error check
      //
#if 0
      if (lastDiagonal == _hits[i]._diagonalID) {
        if (_hits[i]._dsPos < dsLow) {
          fprintf(stderr, "ERROR:  UNSORTED HIT!\n");
          exit(1);
        }
      }
#endif


      //
      //  Extend if on the same diagonal, and consecutive sequence.
      //
      if ((lastDiagonal == _hits[i]._diagonalID) &&
          (qsLow <= _hits[i]._qsPos) &&
          (_hits[i]._qsPos <= qsHigh + config._merSize)) {
        if (qsLow  > _hits[i]._qsPos)   qsLow  = _hits[i]._qsPos;
        if (qsHigh < _hits[i]._qsPos)   qsHigh = _hits[i]._qsPos;
        if (dsLow  > _hits[i]._dsPos)   dsLow  = _hits[i]._dsPos;
        if (dsHigh < _hits[i]._dsPos)   dsHigh = _hits[i]._dsPos;
      } else {

        //
        //  Save the match.  cut-n-paste with below.
        //

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

        if (direction == 'r') {
          sprintf(theOutput + theOutputPos, HITOUTPUTLINE,
                  direction,
                  _qsIdx,
                  _qsLen - qsHigh - config._merSize,
                  qsHigh - qsLow + config._merSize,
                  config._useList[currentSeq].seq,
                  dsLow,
                  dsHigh - dsLow + config._merSize);
        } else {
          sprintf(theOutput + theOutputPos, HITOUTPUTLINE,
                  direction,
                  _qsIdx,
                  qsLow,
                  qsHigh - qsLow + config._merSize,
                  config._useList[currentSeq].seq,
                  dsLow,
                  dsHigh - dsLow + config._merSize);
        }
        while (theOutput[theOutputPos])
          theOutputPos++;

        pthread_mutex_lock(&queryMatchMutex);
        queryMatchCounts[_qsIdx]++;
        pthread_mutex_unlock(&queryMatchMutex);


        lastDiagonal = _hits[i]._diagonalID;
        qsLow        = _hits[i]._qsPos;
        qsHigh       = _hits[i]._qsPos;
        dsLow        = _hits[i]._dsPos;
        dsHigh       = _hits[i]._dsPos;
      }
    }

    //  Save the final cluster?  (cut-n-paste from above)
    //
    if (theOutputPos + 128 >= theOutputMax) {
      theOutputMax <<= 1;
      char *o = new char [theOutputMax];
      memcpy(o, theOutput, theOutputPos);
      delete [] theOutput;
      theOutput = o;
    }

    if (direction == 'r') {
      sprintf(theOutput + theOutputPos, HITOUTPUTLINE,
              direction,
              _qsIdx,
              _qsLen - qsHigh - config._merSize,
              qsHigh - qsLow + config._merSize,
              config._useList[currentSeq].seq,
              dsLow,
              dsHigh - dsLow + config._merSize);
    } else {
      sprintf(theOutput + theOutputPos, HITOUTPUTLINE,
              direction,
              _qsIdx,
              qsLow,
              qsHigh - qsLow + config._merSize,
              config._useList[currentSeq].seq,
              dsLow,
              dsHigh - dsLow + config._merSize);
    }
    while (theOutput[theOutputPos])
      theOutputPos++;

    pthread_mutex_lock(&queryMatchMutex);
    queryMatchCounts[_qsIdx]++;
    pthread_mutex_unlock(&queryMatchMutex);

    //  All done with these hits.  Move to the next set.
    //
    firstHit = lastHit;
  }
}


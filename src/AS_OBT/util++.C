#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>     //  gettimeofday()
#include <sys/utsname.h>  //  uname()
#include <sys/resource.h> //  getrusage()

#include "util++.H"


//  Define this to print some debugging information
//
//#define DEBUG_LIST


intervalList::intervalList() {
  _isSorted = true;
  _listLen  = 0;
  _listMax  = 16;
  _list     = new _intervalPair [_listMax];
}


intervalList::~intervalList() {
  delete [] _list;
}


intervalList &
intervalList::operator=(intervalList &src) {
  _isSorted = src._isSorted;
  _listLen = src._listLen;

  if (_listMax < src._listMax) {
    delete [] _list;
    _listMax = src._listMax;
    _list    = new _intervalPair [_listMax];
  }

  memcpy(_list, src._list, _listLen * sizeof(_intervalPair));

  return(*this);
}

void
intervalList::add(intervalNumber position, intervalNumber length) {

#ifdef DEBUG_LIST
  fprintf(stderr, "Adding %u - %u\n", position, position+length);
#endif

  if (_listLen >= _listMax) {
    _listMax *= 2;
    _intervalPair *l = new _intervalPair [_listMax];
    memcpy(l, _list, sizeof(_intervalPair) * _listLen);
    delete [] _list;
    _list = l;
  }

  _list[_listLen].lo = position;
  _list[_listLen].hi = position + length;

  if ((_listLen > 0) &&
      (_list[_listLen-1].lo > _list[_listLen].lo)) {
#ifdef DEBUG_LIST
    fprintf(stderr, "list isn't sorted\n");
#endif
    _isSorted = false;
  }

  _listLen++;
}


static
int
intervalList_sort_helper(const void *a, const void *b) {
  _intervalPair *A = (_intervalPair *)a;
  _intervalPair *B = (_intervalPair *)b;

  if (A->lo < B->lo) return(-1);
  if (A->lo > B->lo) return(1);
  if (A->hi < B->hi) return(-1);
  if (A->hi > B->hi) return(1);
  return(0);
}


void
intervalList::sort(void) {

  if (_isSorted) {
#ifdef DEBUG_LIST
    fprintf(stderr, "List is already sorted!\n");
#endif
    return;
  }

  if (_listLen > 1)
    qsort(_list, _listLen, sizeof(_intervalPair), intervalList_sort_helper);

  _isSorted = true;
}


void
intervalList::merge(void) {
  uint32  thisInterval  = 0;
  uint32  nextInterval = 1;

  if (_listLen < 2)
    return;

  sort();

  while (nextInterval < _listLen) {

    if ((_list[thisInterval].lo == 0) &&
        (_list[thisInterval].hi == 0)) {

      //  Our interval is empty.  Copy in the interval we are
      //  examining and move to the next.

      //  XXX This is probably useless, thisInterval should always be
      //  valid.

      _list[thisInterval].lo = _list[nextInterval].lo;
      _list[thisInterval].hi = _list[nextInterval].hi;
      nextInterval++;
    } else {

      //  This interval is valid.  See if it overlaps with the next
      //  interval.

      if (_list[thisInterval].hi >= _list[nextInterval].lo) {

        //  Got an intersection.

        //  Merge nextInterval into thisInterval -- the hi range
        //  is extended if the nextInterval range is larger.
        //
        if (_list[thisInterval].hi < _list[nextInterval].hi)
          _list[thisInterval].hi = _list[nextInterval].hi;
        
        //  Clear the just merged nextInterval and move to the next one.
        //
        _list[nextInterval].lo = 0;
        _list[nextInterval].hi = 0;
        nextInterval++;
      } else {

        //  No intersection.  Move along.  Nothing to see here.

        //  If there is a gap between the target and the examine (we
        //  must have merged sometime in the past), copy examine to
        //  the next target.

        thisInterval++;

        if (thisInterval != nextInterval) {
          _list[thisInterval].lo = _list[nextInterval].lo;
          _list[thisInterval].hi = _list[nextInterval].hi;
        }

        nextInterval++;
      }
    }
  }

  if (thisInterval+1 < _listLen)
    _listLen = thisInterval + 1;
}




uint32
intervalList::overlapping(intervalNumber    rangelo,
                          intervalNumber    rangehi,
                          uint32          *&intervals,
                          uint32           &intervalsLen,
                          uint32           &intervalsMax) {


  //  XXX: Naive implementation that is easy to verify (and that works
  //  on an unsorted list).

  if (intervals == 0L) {
    intervalsMax = 256;
    intervals    = new uint32 [intervalsMax];
  }

  intervalsLen = 0;

  for (uint32 i=0; i<_listLen; i++) {
    if ((rangelo <= _list[i].hi) &&
        (rangehi >= _list[i].lo)) {
      if (intervalsLen >= intervalsMax) {
        intervalsMax *= 2;
        uint32 *X = new uint32 [intervalsMax];
        memcpy(X, intervals, sizeof(uint32) * intervalsLen);
        delete [] intervals;
        intervals = X;
      }

      intervals[intervalsLen++] = i;
    }
  }

  return(intervalsLen);
}

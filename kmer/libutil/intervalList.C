#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util++.H"


intervalList::intervalList() {
  _isSorted = true;
  _isMerged = true;
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
  _isMerged = src._isMerged;
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

  if (_listLen >= _listMax) {
    _listMax *= 2;
    _intervalPair *l = new _intervalPair [_listMax];
    memcpy(l, _list, sizeof(_intervalPair) * _listLen);
    delete [] _list;
    _list = l;
  }

  _list[_listLen].lo   = position;
  _list[_listLen].hi   = position + length;

#if 0
  //  Aborted attempt to add a data field here.  Got stuck
  //  deciding how to handle merges lightweight

  _list[_listLen].data = 0L;

  if (data != ~u64bitZERO) {
    _list[_listLen].dataLen = 1;
    _list[_listLen].dataMax = 4;
    _list[_listLen].data    = new u64bit [_list[_listLen].dataMax];
    _list[_listLen].data[0] = data;
  }
#endif
    
  if ((_listLen > 0) &&
      (_list[_listLen-1].lo > _list[_listLen].lo)) {
    _isSorted = false;
  }

  //  XXX  not sure if we need >= here or not
  //
  if ((_listLen > 0) &&
      (_list[_listLen-1].lo >= _list[_listLen].lo)) {
    _isMerged = false;
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

  if (_isSorted)
    return;

  if (_listLen > 1)
    qsort(_list, _listLen, sizeof(_intervalPair), intervalList_sort_helper);

  _isSorted = true;
}


void
intervalList::merge(void) {
  u32bit  thisInterval  = 0;
  u32bit  nextInterval = 1;

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

  _isMerged = true;
}


void
intervalList::invert(intervalNumber lo, intervalNumber hi) {

  if (!_isSorted || !_isMerged) {
    fprintf(stderr, "intervalList::invert()--  ERROR!  List is not sorted or not merged!\n");
    exit(1);
  }

  //  Create a new list to store the inversion
  //
  u32bit             invLen = 0;
  u32bit             invMax = _listLen + 2;
  _intervalPair     *inv    = new _intervalPair [invMax];

  //  Add the first
  //
  if (lo < _list[0].lo) {
    inv[invLen].lo = lo;
    inv[invLen].hi = _list[0].lo;
    invLen++;
  }

  //  Add the pieces
  for (u32bit i=1; i<_listLen; i++) {
    if (_list[i-1].hi < _list[i].lo) {
      inv[invLen].lo = _list[i-1].hi;
      inv[invLen].hi = _list[i].lo;
      invLen++;
    }
  }

  //  Add the last
  if (_list[_listLen-1].hi < hi) {
    inv[invLen].lo = _list[_listLen-1].hi;
    inv[invLen].hi = hi;
    invLen++;
  }

  //  Nuke the old list, swap in the new one
  delete [] _list;

  _list = inv;
  _listLen = invLen;
  _listMax = invMax;
}



u32bit
intervalList::overlapping(intervalNumber    rangelo,
                          intervalNumber    rangehi,
                          u32bit          *&intervals,
                          u32bit           &intervalsLen,
                          u32bit           &intervalsMax) {


  //  XXX: Naive implementation that is easy to verify (and that works
  //  on an unsorted list).

  if (intervals == 0L) {
    intervalsMax = 256;
    intervals    = new u32bit [intervalsMax];
  }

  intervalsLen = 0;

  for (u32bit i=0; i<_listLen; i++) {
    if ((rangelo <= _list[i].hi) &&
        (rangehi >= _list[i].lo)) {
      if (intervalsLen >= intervalsMax) {
        intervalsMax *= 2;
        u32bit *X = new u32bit [intervalsMax];
        memcpy(X, intervals, sizeof(u32bit) * intervalsLen);
        delete [] intervals;
        intervals = X;
      }

      intervals[intervalsLen++] = i;
    }
  }

  return(intervalsLen);
}



void
intervalList::intersect(intervalList &A,
                        intervalList &B) {
  A.merge();
  B.merge();

  u32bit  ai = 0;
  u32bit  bi = 0;

  while ((ai < A.numberOfIntervals()) &&
         (bi < B.numberOfIntervals())) {
    u32bit   al = A.lo(ai);
    u32bit   ah = A.hi(ai);
    u32bit   bl = B.lo(bi);
    u32bit   bh = B.hi(bi);
    u32bit   nl = 0;
    u32bit   nh = 0;

    //  If they intersect, make a new region
    //
    if ((al <= bl) && (bl < ah)) {
      nl = bl;
      nh = (ah < bh) ? ah : bh;
    }

    if ((bl <= al) && (al < bh)) {
      nl = al;
      nh = (ah < bh) ? ah : bh;
    }

    if (nh - nl > 0)
      add(nl, nh - nl);

    //  Advance the list with the earlier region.
    //
    if        (al < bl) {
      //  A begins before B
      ai++;
    } else if (al > bl) {
      //  B begins before A
      bi++;
    } else if (ah < bh) {
      //  A and B begin at the same point, A ends first
      ai++;
    } else if (ah > bh) {
      //  A and B begin at the same point, B ends first
      bi++;
    } else {
      //  Exactly the same!
      ai++;
      bi++;
    }
  }
}








static
int
intervalDepth_sort_helper(const void *a, const void *b) {
  _intervalDepth *A = (_intervalDepth *)a;
  _intervalDepth *B = (_intervalDepth *)b;

  if (A->lo < B->lo) return(-1);
  if (A->lo > B->lo) return(1);
  if (A->hi < B->hi) return(-1);
  if (A->hi > B->hi) return(1);
  return(0);
}


intervalDepth::intervalDepth(intervalList &IL) {

  u32bit           idlen = IL.numberOfIntervals() * 2;
  _intervalDepth  *id    = new _intervalDepth [idlen];

  for (u32bit i=0; i<IL.numberOfIntervals(); i++) {
    id[2*i  ].lo = IL.lo(i);
    id[2*i  ].hi = 0;
    id[2*i  ].de = 1;
    id[2*i+1].lo = IL.hi(i);
    id[2*i+1].hi = 0;
    id[2*i+1].de = 0;
  }

  qsort(id, idlen, sizeof(_intervalDepth), intervalDepth_sort_helper);

  //  Scan the list, counting how many times we change depth.
  //
  _listMax = 1;
  for (u32bit i=1; i<idlen; i++) {
    if (id[i-1].lo != id[i].lo)
      _listMax++;
  }

  //  Allocate the real depth of coverage intervals
  //
  _listLen = 0;
  _list    = new _intervalDepth [_listMax];

  //  Build new intervals
  //
  //  Initialize the first interval
  //
  _list[_listLen].lo = id[0].lo;
  _list[_listLen].hi = id[0].lo;
  _list[_listLen].de = 1;

  for (u32bit i=1; i<idlen; i++) {

    if (_list[_listLen].de == 0) {
      //  Update the start position if the current interval is at zero
      //  depth.
      //
      _list[_listLen].lo = id[i].lo;
    } else {

      //  If we are at a position different from the start, we need to
      //  close out the current interval and make a new one.
      //
      if (id[i-1].lo != id[i].lo) {
        _list[_listLen].hi = id[i].lo;

        _listLen++;

        _list[_listLen].lo = id[i].lo;
        _list[_listLen].hi = id[i].lo;
        _list[_listLen].de = _list[_listLen-1].de;
      }
    }

    //  Finally, update the depth of the current interval
    //
    if (id[i].de)
      _list[_listLen].de++;
    else
      _list[_listLen].de--;
  }

  delete [] id;
}

intervalDepth::~intervalDepth() {
  delete [] _list;
}

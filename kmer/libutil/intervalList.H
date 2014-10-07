#ifndef INTERVALLIST_H
#define INTERVALLIST_H

#include <algorithm>

//using namespace std;

template <class iNum>
class  _intervalPair {
public:
  iNum      lo;
  iNum      hi;
  uint32    ct;

  bool  operator<(const _intervalPair &that) const {
    if (lo != that.lo)
      return(lo < that.lo);
    return(hi < that.hi);

#if 0
    if (lo < that.lo) return(-1);
    if (lo > that.lo) return(1);
    if (hi < that.hi) return(-1);
    if (hi > that.hi) return(1);
    return(0);
#endif
  };
};



template <class iNum>
class  _intervalDepth {
public:
  iNum      lo;
  iNum      hi;
  uint32    de;
};



template <class iNum>
class  intervalDepthRegions {
public:
  iNum      pos;  //  Position of the change in depth
  int32     cha;  //  The change in depth, e.g, +1, -1

  bool  operator<(const intervalDepthRegions &that) const {
    return(pos < that.pos);
  };
};



template <class iNum>
class intervalList {
public:
  intervalList(uint32 initialSize=32) {
    _isSorted = true;
    _isMerged = true;
    _listLen  = 0;
    _listMax  = initialSize;
    _list     = new _intervalPair<iNum> [_listMax];
  };

  ~intervalList() {
    delete [] _list;
  };

  intervalList<iNum> &operator=(intervalList<iNum> &src);

  void        clear(void) {
    _isSorted = true;
    _isMerged = true;
    _listLen  = 0;
  }

  void        add(iNum position, iNum length);
  void        sort(void);
  void        merge(uint32 minOverlap=0);         //  Merge overlapping regions
  void        merge(intervalList<iNum> *IL);      //  Insert IL into this list

  void        intersect(intervalList<iNum> &A,
                        intervalList<iNum> &B);

  uint32      overlapping(iNum      lo,
                          iNum      hi,
                          uint32  *&intervals,
                          uint32   &intervalsLen,
                          uint32   &intervalsMax);

  //  Populates this intervalList with regions in A that are completely
  //  contained in a region in B.
  //
  //  Both A and B call merge().
  //
  void        contained(intervalList<iNum> &A,
                        intervalList<iNum> &B);

  void        invert(iNum lo, iNum hi);

  uint32      numberOfIntervals(void)   { return(_listLen); };

  iNum      sumOfLengths(void) {
    iNum len = 0;
    uint32         i   = numberOfIntervals();

    if (i > 0)
      while (i--)
        len += _list[i].hi - _list[i].lo;

    return(len);
  };

  iNum               &lo(uint32 i) { return(_list[i].lo); };
  iNum               &hi(uint32 i) { return(_list[i].hi); };
  uint32             &ct(uint32 i) { return(_list[i].ct); };

private:
  bool                      _isSorted;
  bool                      _isMerged;
  uint32                    _listLen;
  uint32                    _listMax;
  _intervalPair<iNum>      *_list;
};



//  Takes as input an intervalList, computes the number of intervals
//  covering every position in there, stores this as a new set of
//  intervals, annotated with the depth.
//
//  This is a static object, initialized once by the intervalList.
//
template <class iNum>
class intervalDepth {
public:
  intervalDepth(intervalList<iNum> &IL);
  intervalDepth(intervalDepthRegions<iNum> *id, uint32 idlen);
  ~intervalDepth();

  uint32                     numberOfIntervals(void)   { return(_listLen); };

  iNum                       lo(uint32 i) { return(_list[i].lo); };
  iNum                       hi(uint32 i) { return(_list[i].hi); };
  uint32                     de(uint32 i) { return(_list[i].de); };

private:
  void                       computeIntervals(intervalDepthRegions<iNum> *id, uint32 idlen);

  uint32                    _listLen;
  uint32                    _listMax;
  _intervalDepth<iNum>     *_list;
};







template <class iNum>
intervalList<iNum> &
intervalList<iNum>::operator=(intervalList &src) {
  _isSorted = src._isSorted;
  _isMerged = src._isMerged;
  _listLen = src._listLen;

  if (_listMax < src._listMax) {
    delete [] _list;
    _listMax = src._listMax;
    _list    = new _intervalPair<iNum> [_listMax];
  }

  memcpy(_list, src._list, _listLen * sizeof(_intervalPair<iNum>));

  return(*this);
}


template <class iNum>
void
intervalList<iNum>::add(iNum position, iNum length) {

  if (_listLen >= _listMax) {
    _listMax *= 2;
    _intervalPair<iNum> *l = new _intervalPair<iNum> [_listMax];
    memcpy(l, _list, sizeof(_intervalPair<iNum>) * _listLen);
    delete [] _list;
    _list = l;
  }

  _list[_listLen].lo   = position;
  _list[_listLen].hi   = position + length;
  _list[_listLen].ct   = 1;

#if 0
  //  Aborted attempt to add a data field here.  Got stuck
  //  deciding how to handle merges lightweight

  _list[_listLen].data = 0L;

  if (data != ~uint64ZERO) {
    _list[_listLen].dataLen = 1;
    _list[_listLen].dataMax = 4;
    _list[_listLen].data    = new uint64 [_list[_listLen].dataMax];
    _list[_listLen].data[0] = data;
  }
#endif
    
  if ((_listLen > 0) &&
      (_list[_listLen-1].lo > _list[_listLen].lo)) {
    _isSorted = false;
    _isMerged = false;
  }

  _listLen++;
}


#if 0
template <class iNum>
static
int
intervalList_sort_helper(const void *a, const void *b) {
  _intervalPair<iNum> *A = (_intervalPair<iNum> *)a;
  _intervalPair<iNum> *B = (_intervalPair<iNum> *)b;

  if (A->lo < B->lo) return(-1);
  if (A->lo > B->lo) return(1);
  if (A->hi < B->hi) return(-1);
  if (A->hi > B->hi) return(1);
  return(0);
}
#endif


template <class iNum>
void
intervalList<iNum>::sort(void) {

  if (_isSorted)
    return;

  if (_listLen > 1)
    std::sort(_list, _list + _listLen);

  _isSorted = true;
}


template <class iNum>
void
intervalList<iNum>::merge(uint32 minOverlap) {
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
      _list[thisInterval].ct = _list[nextInterval].ct;

      _list[nextInterval].lo = 0;
      _list[nextInterval].hi = 0;
      nextInterval++;
    } else {

      //  This interval is valid.  See if it overlaps with the next
      //  interval.

      bool  intersects = false;

      if ((_list[thisInterval].lo <= _list[nextInterval].lo) &&
          (_list[nextInterval].hi <= _list[thisInterval].hi))
        //  next is contained in this
        intersects = true;

      if (_list[thisInterval].hi - minOverlap >= _list[nextInterval].lo)
        //  next has thick overlap to this
        intersects = true;


      if (intersects) {

        //  Got an intersection.

        //  Merge nextInterval into thisInterval -- the hi range
        //  is extended if the nextInterval range is larger.
        //
        if (_list[thisInterval].hi < _list[nextInterval].hi)
          _list[thisInterval].hi = _list[nextInterval].hi;

        _list[thisInterval].ct += _list[nextInterval].ct;
        
        //  Clear the just merged nextInterval and move to the next one.
        //
        _list[nextInterval].lo = 0;
        _list[nextInterval].hi = 0;
        _list[nextInterval].ct = 0;
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
          _list[thisInterval].ct = _list[nextInterval].ct;
        }

        nextInterval++;
      }
    }
  }

  if (thisInterval+1 < _listLen)
    _listLen = thisInterval + 1;

  _isMerged = true;
}



template <class iNum>
void
intervalList<iNum>::merge(intervalList<iNum> *IL) {
  for (uint32 i=0; i<IL->_listLen; i++)
    add(IL->_list[i].lo, IL->_list[i].hi - IL->_list[i].lo);
}



template <class iNum>
void
intervalList<iNum>::invert(iNum lo, iNum hi) {

  if (!_isSorted || !_isMerged) {
    fprintf(stderr, "intervalList<iNum>::invert()--  ERROR!  List is not sorted or not merged!\n");
    exit(1);
  }

  //  Create a new list to store the inversion
  //
  uint32                   invLen = 0;
  uint32                   invMax = _listLen + 2;
  _intervalPair<iNum>     *inv    = new _intervalPair<iNum> [invMax];

  //  Add the zeroth and only?
  if (_listLen == 0) {
    inv[invLen].lo = lo;
    inv[invLen].hi = hi;
    invLen++;
  }

  //  Add the first, then the pieces, then the last
  //
  else {
    if (lo < _list[0].lo) {
      inv[invLen].lo = lo;
      inv[invLen].hi = _list[0].lo;
      invLen++;
    }

    for (uint32 i=1; i<_listLen; i++) {
      if (_list[i-1].hi < _list[i].lo) {
        inv[invLen].lo = _list[i-1].hi;
        inv[invLen].hi = _list[i].lo;
        invLen++;
      }
    }

    if (_list[_listLen-1].hi < hi) {
      inv[invLen].lo = _list[_listLen-1].hi;
      inv[invLen].hi = hi;
      invLen++;
    }
  }

  //  Nuke the old list, swap in the new one
  delete [] _list;

  assert(invLen <= invMax);

  _list = inv;
  _listLen = invLen;
  _listMax = invMax;
}



template <class iNum>
void
intervalList<iNum>::intersect(intervalList<iNum> &A,
                              intervalList<iNum> &B) {
  A.merge();
  B.merge();

  uint32  ai = 0;
  uint32  bi = 0;

  while ((ai < A.numberOfIntervals()) &&
         (bi < B.numberOfIntervals())) {
    uint32   al = A.lo(ai);
    uint32   ah = A.hi(ai);
    uint32   bl = B.lo(bi);
    uint32   bh = B.hi(bi);
    uint32   nl = 0;
    uint32   nh = 0;

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

    if (nl < nh)
      add(nl, nh - nl);

    //  Advance the list with the earlier region.
    //
    if        (ah < bh) {
      //  A ends before B
      ai++;
    } else if (ah > bh) {
      //  B ends before A
      bi++;
    } else {
      //  Exactly the same ending!
      ai++;
      bi++;
    }
  }
}



//  Populates an array with the intervals that are within the supplied interval.
//
//  Naive implementation that is easy to verify (and that works on an unsorted list).
//
template <class iNum>
uint32
intervalList<iNum>::overlapping(iNum      rangelo,
                                iNum      rangehi,
                                uint32  *&intervals,
                                uint32   &intervalsLen,
                                uint32   &intervalsMax) {

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



template <class iNum>
void
intervalList<iNum>::contained(intervalList<iNum> &A,
                              intervalList<iNum> &B) {
  A.merge();
  B.merge();

  uint32  ai = 0;
  uint32  bi = 0;

  while ((ai < A.numberOfIntervals()) &&
         (bi < B.numberOfIntervals())) {
    uint32   al = A.lo(ai);
    uint32   ah = A.hi(ai);
    uint32   bl = B.lo(bi);
    uint32   bh = B.hi(bi);

    //  If A is contained in B, make a new region.
    //
    if ((bl <= al) && (ah <= bh))
      add(bl, bh - bl);

#if 0
    if ((al <= bl) && (bh <= ah))
      add(al, ah - al);
#endif

    //  Advance the list with the earlier region.
    //
    if        (ah < bh) {
      //  A ends before B
      ai++;
    } else if (ah > bh) {
      //  B ends before A
      bi++;
    } else {
      //  Exactly the same ending!
      ai++;
      bi++;
    }
  }
}


#if 0
template <class iNum>
static
int
intervalDepth_sort_helper(const void *a, const void *b) {
  intervalDepthRegions<iNum> *A = (intervalDepthRegions<iNum> *)a;
  intervalDepthRegions<iNum> *B = (intervalDepthRegions<iNum> *)b;

  if (A->pos < B->pos) return(-1);
  if (A->pos > B->pos) return(1);

  return(0);
}
#endif

template <class iNum>
intervalDepth<iNum>::intervalDepth(intervalList<iNum> &IL) {
  uint32                       idlen = IL.numberOfIntervals() * 2;
  intervalDepthRegions<iNum>  *id    = new intervalDepthRegions<iNum> [idlen];

  for (uint32 i=0; i<IL.numberOfIntervals(); i++) {
    id[2*i  ].pos = IL.lo(i);
    id[2*i  ].cha = 1;
    id[2*i+1].pos = IL.hi(i);
    id[2*i+1].cha = -1;
  }

  std::sort(id, id + idlen);

  computeIntervals(id, idlen);

  delete [] id;
}


template <class iNum>
intervalDepth<iNum>::intervalDepth(intervalDepthRegions<iNum> *id, uint32 idlen) {
  std::sort(id, id + idlen);  //, sizeof(intervalDepthRegions<iNum>), intervalDepth_sort_helper);
  computeIntervals(id, idlen);
}


template <class iNum>
void
intervalDepth<iNum>::computeIntervals(intervalDepthRegions<iNum> *id, uint32 idlen) {

  //  No intervals input?  No intervals output.

  if (idlen == 0) {
    _listLen = 0;
    _listMax = 0;
    _list    = NULL;

    return;
  }

  //  Scan the list, counting how many times we change depth.

  _listMax = 1;
  for (uint32 i=1; i<idlen; i++) {
    if (id[i-1].pos != id[i].pos)
      _listMax++;
  }

  //  Allocate the real depth of coverage intervals

  _listLen = 0;
  _list    = new _intervalDepth<iNum> [_listMax];

  //  Build new intervals

  //  Init first interval.
  _list[_listLen].lo = id[0].pos;
  _list[_listLen].hi = id[0].pos;
  _list[_listLen].de = id[0].cha;

  for (uint32 i=1; i<idlen; i++) {
    //  Update the end of the current interval.
    _list[_listLen].hi = id[i].pos;

    //fprintf(stderr, "ID %d %d - cur %d %d-%d %d\n",
    //        id[i].pos, id[i].cha,
    //        _listLen, _list[_listLen].lo, _list[_listLen].hi, _list[_listLen].de);

    //  If the position is different than the last, make a new interval
    if (id[i-1].pos != id[i].pos) {
      _listLen++;

      _list[_listLen].lo = id[i].pos;
      _list[_listLen].de = _list[_listLen-1].de;
    }

    //  Finally, update the depth and end of the current interval
    _list[_listLen].hi  = id[i].pos;
    _list[_listLen].de += id[i].cha;

    //fprintf(stderr, "ID %d %d - cur %d %d-%d %d POST\n",
    //        id[i].pos, id[i].cha,
    //        _listLen, _list[_listLen].lo, _list[_listLen].hi, _list[_listLen].de);
  }

  assert(_listLen > 0);
  assert(_listLen < _listMax);
}

template <class iNum>
intervalDepth<iNum>::~intervalDepth() {
  delete [] _list;
}




#endif  //  INTERVALLIST_H

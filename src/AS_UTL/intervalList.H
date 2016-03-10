
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    kmer/libutil/intervalList.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-AUG-11 to 2003-AUG-13
 *      are Copyright 2003 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2004-APR-21
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUL-12 to 2014-APR-11
 *      are Copyright 2005-2007,2010,2013-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-07 to 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-30
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INTERVALLIST_H
#define INTERVALLIST_H

#include <algorithm>

//  iNum   - lo, hi - coordinates of the interval
//  iVal   - va     - data stored at each interval
//  uint32 - ct     - number of elements in this interval
//                  - when merged, needs function that converts multiple iVal and a uint32 into a single iVal

template <class iNum, class iVal=int32>
class  _intervalPair {
public:
  iNum      lo;
  iNum      hi;
  uint32    ct;  //  Number of source intervals
  iVal      va;  //  Value at this interval; default is 1

  bool  operator<(const _intervalPair &that) const {
    if (lo != that.lo)
      return(lo < that.lo);
    return(hi < that.hi);
  };
};


template <class iNum, class iVal=int32>
class  intervalDepthRegions {
public:
  iNum      pos;     //  Position of the change in depth
  iVal      change;  //  The value associated with this object; added or subtracted from 'va'.
  bool      open;    //  If true, the start of a new interval

  bool  operator<(const intervalDepthRegions &that) const {
    if (pos != that.pos)
      return(pos < that.pos);
    return(open > that.open);
  };
};





template <class iNum, class iVal=int32>
class intervalList {
public:
  intervalList(uint32 initialSize=32) {
    _isSorted = true;
    _isMerged = true;
    _listLen  = 0;
    _listMax  = initialSize;
    _list     = new _intervalPair<iNum, iVal> [_listMax];
  };

  //  Takes as input an unmerged intervalList, returns to a new set of intervals, one
  //  for each 'depth'.  Two intervals, (1,4) and (2,6) would return 'depths':
  //    1,2,1  bgn=1, end=2, depth=1
  //    2,4,2
  //    4,6,1
  //
  intervalList(intervalList<iNum, iVal> &IL) {
    _isSorted = false;
    _isMerged = false;
    _listLen  = 0;
    _listMax  = 0;
    _list     = 0L;

    depth(IL);
  };

  intervalList(intervalDepthRegions<iNum, iVal> *id, uint32 idlen) {
    _isSorted = false;
    _isMerged = false;
    _listLen  = 0;
    _listMax  = 0;
    _list     = 0L;

#ifdef _GLIBCXX_PARALLEL
    //  Don't use the parallel sort, not with the expense of starting threads.
    __gnu_sequential::sort(id, id + idlen);
#else
    std::sort(id, id + idlen);
#endif

    computeDepth(id, idlen);
  };

  ~intervalList() {
    delete [] _list;
  };

  intervalList<iNum, iVal> &operator=(intervalList<iNum, iVal> &src);

  void      clear(void) {
    _isSorted = true;
    _isMerged = true;
    _listLen  = 0;
  }

  void      add(iNum position, iNum length, iVal value=0);
  void      sort(void);
  void      merge(iNum minOverlap=0);                 //  Merge overlapping regions
  void      merge(intervalList<iNum, iVal> *IL);      //  Insert IL into this list

  void      filterShort(iNum minLength);

  void      intersect(intervalList<iNum, iVal> &A,
                      intervalList<iNum, iVal> &B);

  uint32    overlapping(iNum      lo,
                        iNum      hi,
                        uint32  *&intervals,
                        uint32   &intervalsLen,
                        uint32   &intervalsMax);

  //  Populates this intervalList with regions in A that are completely
  //  contained in a region in B.
  //
  //  Both A and B call merge().
  //
  void      contained(intervalList<iNum, iVal> &A,
                      intervalList<iNum, iVal> &B);

  void      invert(iNum lo, iNum hi);

  void      depth(intervalList<iNum, iVal> &A);

  uint32    numberOfIntervals(void)   { return(_listLen); };

  iNum      sumOfLengths(void) {
    iNum len = 0;
    uint32         i   = numberOfIntervals();

    if (i > 0)
      while (i--)
        len += _list[i].hi - _list[i].lo;

    return(len);
  };

  iNum     &lo(uint32 i)    { return(_list[i].lo); };
  iNum     &hi(uint32 i)    { return(_list[i].hi); };

  uint32   &count(uint32 i) { return(_list[i].ct); };  //  Number of source intervals.
  uint32   &depth(uint32 i) { return(_list[i].ct); };  //  Depth, if converted.
  iVal     &value(uint32 i) { return(_list[i].va); };  //  Value or sum of values.

private:
  void     computeDepth(intervalDepthRegions<iNum, iVal> *id, uint32 idlen);


  bool                         _isSorted;
  bool                         _isMerged;

  uint32                       _listMax;
  uint32                       _listLen;
  _intervalPair<iNum, iVal>   *_list;
};






template <class iNum, class iVal>
intervalList<iNum, iVal> &
intervalList<iNum, iVal>::operator=(intervalList &src) {
  _isSorted = src._isSorted;
  _isMerged = src._isMerged;


  if (_listMax < src._listMax) {
    delete [] _list;
    _listMax = src._listMax;
    _list    = new _intervalPair<iNum, iVal> [_listMax];
  }

  _listLen  = src._listLen;

  memcpy(_list, src._list, _listLen * sizeof(_intervalPair<iNum, iVal>));

  return(*this);
}


template <class iNum, class iVal>
void
intervalList<iNum, iVal>::add(iNum position, iNum length, iVal val) {

  if (_listLen >= _listMax) {
    _listMax *= 2;
    _intervalPair<iNum, iVal> *l = new _intervalPair<iNum, iVal> [_listMax];
    memcpy(l, _list, sizeof(_intervalPair<iNum, iVal>) * _listLen);
    delete [] _list;
    _list = l;
  }

  _list[_listLen].lo   = position;
  _list[_listLen].hi   = position + length;
  _list[_listLen].ct   = 1;
  _list[_listLen].va   = val;

  //  Could optimize, and search the list to see if these are false,
  //  but that's rather expensive.
  _isSorted = false;
  _isMerged = false;

  _listLen++;
}



template <class iNum, class iVal>
void
intervalList<iNum, iVal>::sort(void) {

  if (_isSorted)
    return;

  if (_listLen > 1)
#ifdef _GLIBCXX_PARALLEL
    //  Don't use the parallel sort, not with the expense of starting threads.
    __gnu_sequential::sort(_list, _list + _listLen);
#else
    std::sort(_list, _list + _listLen);
#endif

  _isSorted = true;
}


template <class iNum, class iVal>
void
intervalList<iNum, iVal>::merge(iNum minOverlap) {
  uint32  thisInterval = 0;
  uint32  nextInterval = 1;

  if (_isMerged)
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
      _list[thisInterval].va = _list[nextInterval].va;

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
        _list[thisInterval].va += _list[nextInterval].va;

        //  Clear the just merged nextInterval and move to the next one.
        //
        _list[nextInterval].lo = 0;
        _list[nextInterval].hi = 0;
        _list[nextInterval].ct = 0;
        _list[nextInterval].va = 0;

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
          _list[thisInterval].va = _list[nextInterval].va;
        }

        nextInterval++;
      }
    }
  }

  if (thisInterval+1 < _listLen)
    _listLen = thisInterval + 1;

  _isMerged = true;
}



template <class iNum, class iVal>
void
intervalList<iNum, iVal>::merge(intervalList<iNum, iVal> *IL) {
  for (uint32 i=0; i<IL->_listLen; i++)
    add(IL->_list[i].lo, IL->_list[i].hi - IL->_list[i].lo);
}



template <class iNum, class iVal>
void
intervalList<iNum, iVal>::filterShort(iNum minLength) {
  uint32   fr = 0;
  uint32   to = 0;

  for (; fr < _listLen; fr++) {
    if (hi(fr) - lo(fr) <= minLength)
      continue;

    if (to != fr)
      _list[to] = _list[fr];

    to++;
  }

  _listLen = to;
}



template <class iNum, class iVal>
void
intervalList<iNum, iVal>::invert(iNum invlo, iNum invhi) {

  merge();

  //  Create a new list to store the inversion
  //
  uint32                         invLen = 0;
  uint32                         invMax = _listLen + 2;
  _intervalPair<iNum, iVal>     *inv    = new _intervalPair<iNum, iVal> [invMax];

  //  Add the zeroth and only?
  if (_listLen == 0) {
    inv[invLen].lo = invlo;
    inv[invLen].hi = invhi;
    inv[invLen].ct = 1;
    inv[invLen].va = 0;
    invLen++;
  }

  //  Add the first, then the pieces, then the last
  //
  else {
    if (invlo < _list[0].lo) {
      inv[invLen].lo = invlo;
      inv[invLen].hi = _list[0].lo;
      inv[invLen].ct = 1;
      inv[invLen].va = 0;
      invLen++;
    }

    for (uint32 i=1; i<_listLen; i++) {
      if (_list[i-1].hi < _list[i].lo) {
        inv[invLen].lo = _list[i-1].hi;
        inv[invLen].hi = _list[i].lo;
        inv[invLen].ct = 1;
        inv[invLen].va = 0;
        invLen++;
      }
    }

    if (_list[_listLen-1].hi < invhi) {
      inv[invLen].lo = _list[_listLen-1].hi;
      inv[invLen].hi = invhi;
      inv[invLen].ct = 1;
      inv[invLen].va = 0;
      invLen++;
    }
  }

  assert(invLen <= invMax);

  //  Nuke the old list, swap in the new one
  delete [] _list;

  _list    = inv;
  _listLen = invLen;
  _listMax = invMax;
}



template <class iNum, class iVal>
void
intervalList<iNum, iVal>::intersect(intervalList<iNum, iVal> &A,
                                    intervalList<iNum, iVal> &B) {
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
template <class iNum, class iVal>
uint32
intervalList<iNum, iVal>::overlapping(iNum      rangelo,
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



template <class iNum, class iVal>
void
intervalList<iNum, iVal>::contained(intervalList<iNum, iVal> &A,
                                    intervalList<iNum, iVal> &B) {
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







template <class iNum, class iVal>
void
intervalList<iNum, iVal>::depth(intervalList<iNum, iVal> &IL) {
  uint32                             idlen = IL.numberOfIntervals() * 2;
  intervalDepthRegions<iNum, iVal>  *id    = new intervalDepthRegions<iNum, iVal> [idlen];

  for (uint32 i=0; i<IL.numberOfIntervals(); i++) {
    id[2*i  ].pos    = IL.lo(i);
    id[2*i  ].change = IL.value(i);
    id[2*i  ].open   = true;

    id[2*i+1].pos    = IL.hi(i);
    id[2*i+1].change = IL.value(i);
    id[2*i+1].open   = false;
  }

  computeDepth(id, idlen);

  delete [] id;
}



template <class iNum, class iVal>
void
intervalList<iNum, iVal>::computeDepth(intervalDepthRegions<iNum, iVal> *id, uint32 idlen) {

  //  No intervals input?  No intervals output.

  _listLen = 0;

  if (idlen == 0)
    return;

  //  Sort by coordinate.

#ifdef _GLIBCXX_PARALLEL
  //  Don't use the parallel sort, not with the expense of starting threads.
  __gnu_sequential::sort(id, id + idlen);
#else
  std::sort(id, id + idlen);
#endif

  //  Allocate the (maximum possible) depth of coverage intervals

  if (_listMax < idlen) {
    delete [] _list;

    _listMax = idlen;
    _list    = new _intervalPair<iNum, iVal> [_listMax];
  }

  //  Init first interval.

  assert(id[0].open == true);

  _list[_listLen].lo = id[0].pos;
  _list[_listLen].hi = id[0].pos;
  _list[_listLen].ct = 1;
  _list[_listLen].va = id[0].change;

  uint32  nct;
  iVal    nva;

  for (uint32 i=1; i<idlen; i++) {
    //  Update the end of the current interval.
    _list[_listLen].hi = id[i].pos;

    //  Compute the count and value of the next interval.
    if (id[i].open == true) {
      nct = _list[_listLen].ct + 1;
      nva = _list[_listLen].va + id[i].change;
    } else {
      nct = _list[_listLen].ct - 1;
      nva = _list[_listLen].va - id[i].change;
    }

    //  If the position or value is different, make a new interval,
    //  But only if this interval is not null length.
    if (((id[i-1].pos        != id[i].pos) ||
         (_list[_listLen].va != nva)) &&
        (_list[_listLen].lo  != _list[_listLen].hi)) {
      _listLen++;

      _list[_listLen].lo = id[i].pos;
      _list[_listLen].ct = _list[_listLen-1].ct;
      _list[_listLen].va = _list[_listLen-1].va;
    }

    //  Finally, update whatver interval is current.
    _list[_listLen].hi  = id[i].pos;
    _list[_listLen].ct  = nct;
    _list[_listLen].va  = nva;

    //  Now, if this interval's begin is the same as the last interval's end,
    //  we need to merge.
    if ((_listLen > 1) &&
        (_list[_listLen-1].hi == _list[_listLen].lo) &&
        (_list[_listLen-1].ct == _list[_listLen].ct) &&
        (_list[_listLen-1].va == _list[_listLen].va)) {
      _list[_listLen-1].hi = _list[_listLen].hi;
      _listLen--;
    }

#if 0
    fprintf(stderr, "id[%2d] - list[%u] = lo=%u hi=%u ct=%u va=%d\n",
            i,
            _listLen,
            _list[_listLen].lo,
            _list[_listLen].hi,
            _list[_listLen].ct,
            _list[_listLen].va);
#endif
  }

  assert(_listLen >  0);
  assert(_listLen <= _listMax);
}



#endif  //  INTERVALLIST_H

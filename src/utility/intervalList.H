
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INTERVALLIST_H
#define INTERVALLIST_H

#include "AS_global.H"

#include <algorithm>


template <class iNum>
class intervalList {
public:
  intervalList(uint32 initialSize=32) {
    _isSorted = true;
    _isMerged = true;
    _listLen  = 0;
    _listMax  = initialSize;
    _list     = new _ip [_listMax];
  };

  ~intervalList() {
    delete [] _list;
  };

  //  Used in OBT
  intervalList<iNum> &operator=(intervalList<iNum> &src);

  void      clear(void) {
    _isSorted = true;
    _isMerged = true;
    _listLen  = 0;
  }

  void      add(iNum position, iNum length);
  void      sort(void);
  void      merge(iNum minOverlap=0);                 //  Merge overlapping regions

  void      filterShort(iNum minLength);

  void      invert(iNum lo, iNum hi);

  uint32    numberOfIntervals(void)   { return(_listLen); };

  iNum     &lo(uint32 i)    { return(_list[i]._bgn); };
  iNum     &hi(uint32 i)    { return(_list[i]._end); };

  //  Used in OBT
  uint32   &count(uint32 i) { return(_list[i]._cnt); };  //  Number of source intervals.

private:
  struct _ip {
    iNum      _bgn;
    iNum      _end;
    uint32    _cnt;  //  Number of source intervals
  };

  bool        _isSorted;
  bool        _isMerged;

  uint32      _listMax;
  uint32      _listLen;
  _ip        *_list;
};




template <class iNum>
intervalList<iNum> &
intervalList<iNum>::operator=(intervalList &src) {
  _isSorted = src._isSorted;
  _isMerged = src._isMerged;

  duplicateArray(_list, _listLen, _listMax, src._list, src._listLen, src._listMax);

  return(*this);
}



template <class iNum>
void
intervalList<iNum>::add(iNum position, iNum length) {

  increaseArray(_list, _listLen, _listMax, _listMax / 4);

  _list[_listLen++] = { position, position + length, 1 };

  _isSorted = false;
  _isMerged = false;

}



template <class iNum>
void
intervalList<iNum>::sort(void) {

  if ((_isSorted == false) && (_listLen > 1))
    std::sort(_list, _list + _listLen, [](_ip const &a, _ip const &b) {
                                         return((a._bgn < b._bgn) || ((a._bgn == b._bgn) && (a._end < b._end))); } );
  _isSorted = true;
}


template <class iNum>
void
intervalList<iNum>::merge(iNum minOverlap) {

  if (_isMerged)
    return;

  sort();

  uint32  thisI = 0;  //  Interval we're merging into.
  uint32  nextI = 1;  //  Interval we're merging from.

  while (nextI < _listLen) {
    assert(_list[thisI]._bgn <  _list[thisI]._end);   //  Basic checks.
    assert(_list[nextI]._bgn <  _list[nextI]._end);   //  Both intervals cannot be empty,
    assert(_list[thisI]._bgn <= _list[nextI]._bgn);   //  and thisI must be before nextI.

    //  If the nextI intersects with thisI -- either contained in thisI, or
    //  has a thick overlap to thisI -- merge it in.  We're guaranteed that
    //  this.bgn is before next.bgn, so all we need to do is extend this.end
    //  to cover the next interval.

    if ((_list[thisI]._end >= _list[nextI]._end) ||
        (_list[thisI]._end >= _list[nextI]._bgn + minOverlap)) {
      _list[thisI]._end  = max(_list[nextI]._end, _list[thisI]._end);
      _list[thisI]._cnt +=     _list[nextI]._cnt;
      nextI++;
    }

    //  Otherwise, move to the next thisI, copy the current nextI to it, and
    //  then move to the next nextI.  We should, to be pedantic, check that
    //  thisI != nextI before the copy, but no harm if we don't.

    else {
      _list[++thisI] = _list[nextI++];
    }
  }

  _listLen  = thisI + 1;   //  Update the length of the list,
  _isMerged = true;        //  and note that it's now merged.
}




//  Filter intervals shorter than 'minLength'.

template <class iNum>
void
intervalList<iNum>::filterShort(iNum minLength) {
  uint32  to = 0;

  for (uint32 fr=0; fr < _listLen; fr++)                //  Over every interval,
    if (_list[fr]._bgn + minLength <= _list[fr]._end)   //  If long enough, copy it
      _list[to++] = _list[fr];                          //  onto the 'new' list.

  _listLen = to;
}





template <class iNum>
void
intervalList<iNum>::invert(iNum invlo, iNum invhi) {

  merge();

  uint32   invLen = 0;                  //  Create a new list to store the
  uint32   invMax = _listLen + 1;       //  inversion.  We need at most one
  _ip     *inv    = new _ip [invMax];   //  more interval than the original.

  //  If no existing list, just add a single interval covering the universe.
  //  If the inversion range falls entirely inside a gap in the original list
  //  (which would also result in the inverted list having one interval
  //  covering the whole range) we'll catch it in the loop below.

  if (_listLen == 0) {
    inv[invLen++] = { invlo, invhi, 1 };
    return;
  }

  //  Add an interval for the first gap, if it's inside the inersion range.

  if (invlo < _list[0]._bgn)
    inv[invLen++] = { invlo, _list[0]._bgn, 1 };

  //  Add intervals covering the middle gaps.  Threshold each endpoint by the
  //  inversion range, and only add a new interval if it is of positive
  //  length.

  for (uint32 i=1; i<_listLen; i++) {
    iNum  bgn = max(invlo,         _list[i-1]._end);
    iNum  end = min(_list[i]._bgn, invhi);

    if (bgn < end)
      inv[invLen++] = { bgn, end, 1 };
  }

  //  Add an interval for the last gap, if it's inside the inversion range.

  if (_list[_listLen-1]._end < invhi)
    inv[invLen++] = { _list[_listLen-1]._end, invhi, 1 };

  //  If there is no inverted list here, then the inversion range was
  //  entirely inside an original interval, i.e., no gaps in the orignal
  //  list.

  delete [] _list;

  _list    = inv;
  _listLen = invLen;
  _listMax = invMax;

  assert(invLen <= invMax);
}







//
//  Takes as input an unmerged intervalList, returns to a new set of intervals, one
//  for each 'depth'.  Two intervals, (1,4) and (2,6) would return 'depths':
//    1,2,1  bgn=1, end=2, depth=1
//    2,4,2
//    4,6,1
//


template <class iNum>
class intervalDepth {
public:
  intervalDepth(intervalList<iNum> &IL) {
    uint32    idlen = IL.numberOfIntervals() * 2;
    _idp     *id    = new _idp [idlen];

    for (uint32 i=0; i<IL.numberOfIntervals(); i++) {
      id[2*i  ].pos = IL.lo(i);     //  Enter into an inteval, change
      id[2*i  ].dlt = 1;            //  depth by +1.

      id[2*i+1].pos = IL.hi(i);     //  Leave an interval, change
      id[2*i+1].dlt = -1;           //  depth by -1.
    }

    _listLen  = 0;
    _listMax  = 0;
    _list     = NULL;

    if (idlen > 0)
      computeDepth(id, idlen);

    delete [] id;
  };

  ~intervalDepth() {
    delete [] _list;
  };

  uint32    numberOfIntervals(void)   { return(_listLen); };

  iNum     &lo(uint32 i)    { return(_list[i]._bgn); };
  iNum     &hi(uint32 i)    { return(_list[i]._end); };
  uint32   &depth(uint32 i) { return(_list[i]._dpt); };

private:
  struct _idp {        //  An intervalDepthPosition stores the position
    iNum      pos;     //  of a change in depth, and the delta of that
    int32     dlt;     //  change (which is either +1 or -1).
  };

  struct _idr {
    iNum      _bgn;
    iNum      _end;
    uint32    _dpt;
  };

private:
  void     computeDepth(_idp *id, uint32 idlen);

  uint32     _listMax;
  uint32     _listLen;
  _idr      *_list;
};



template <class iNum>
void
intervalDepth<iNum>::computeDepth(_idp *id, uint32 idlen) {

  //  Sort regions so that earlier positions are first, and so that depth
  //  increases (+1) are before decreases (-1).

  std::sort(id, id + idlen, [](_idp const &a, _idp const &b) {
                              return((a.pos < b.pos) || ((a.pos == b.pos) && (a.dlt > b.dlt))); } );

  //  The first thing must be an 'open' event.  If not, someone supplied a
  //  negative length to the original intervalList.  Or, possibly, two
  //  zero-length intervals.

  if (id[0].dlt == -1)
    for (uint32 ii=0; ii<idlen; ii++)
      fprintf(stderr, "id[%u] pos %d dlt %d\n", ii, id[ii].pos, id[ii].dlt);

  assert(id[0].dlt == 1);

  //  Init first interval.

  _listLen  = 0;
  _listMax  = idlen + 1;
  _list     = new _idr [_listMax];

  _list[_listLen]._bgn = id[0].pos;
  _list[_listLen]._end = id[0].pos;
  _list[_listLen]._dpt = 1;

  for (uint32 i=1; i<idlen; i++) {

    //  Update the end of the current interval to this position.

    _list[_listLen]._end = id[i].pos;

    //  If this position is different than the last position, make
    //  a new depth interval.

    if (id[i-1].pos != id[i].pos) {
      _listLen++;

      _list[_listLen]._bgn = id[i].pos;
      _list[_listLen]._end = id[i].pos;
      _list[_listLen]._dpt = _list[_listLen-1]._dpt;
    }

    //  Process any depth change associated with this position.

    _list[_listLen]._dpt += id[i].dlt;

    //  Is it safe to blindly change the depth of this region?  Yes.  If the
    //  position is different than the last, we've already made a new depth
    //  region.  And if it wasn't different, one of the previous positions
    //  must have made a new depth region.  In any case, the depth region
    //  we're currently at must always be length 0 -- we're never able to
    //  change the depth of a region when we're at the end coordinate.

    assert(_list[_listLen]._bgn == _list[_listLen]._end);
  }

  assert(_listLen < _listMax);
  assert(_list[_listLen]._bgn == _list[_listLen]._end);
  assert(_list[_listLen]._dpt == 0);
}


#endif  //  INTERVALLIST_H

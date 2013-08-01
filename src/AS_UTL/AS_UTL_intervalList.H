
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

#ifndef AS_UTL_INTERVALLIST_H
#define AS_UTL_INTERVALLIST_H

#include "AS_global.H"

//  This is kmer/libutil/intervalList.[CH], replacing u32bit with uint32.


//  I don't want this to be public, but the sort function
//  in intervalList.C needs it so.
//
typedef int64  intervalNumber;

struct _intervalPair {
  intervalNumber    lo;
  intervalNumber    hi;
  uint32            ct;
};

struct _intervalDepth {
  intervalNumber    lo;
  intervalNumber    hi;
  uint32            de;
};

struct intervalDepthRegions {
  intervalNumber    pos;  //  Position of the change in depth
  int32             cha;  //  The change in depth, e.g, +1, -1
};


class intervalList {
public:
  intervalList(uint32 initialSize=32);
  ~intervalList();

  intervalList &operator=(intervalList &src);

  //  Clear a list
  void        clear(void) {
    _isSorted = true;
    _isMerged = true;
    _listLen  = 0;
  }

  //  Insert a new interval into the list
  void        add(intervalNumber position, intervalNumber length);

  //  Sort the set of intervals by the lo value
  void        sort(void);

  //  Merge overlapping or adjacent intervals together.
  void        merge(uint32 minOverlap=0);

  void        invert(intervalNumber lo, intervalNumber hi);

  //  Returns the number of intervals
  uint32      numberOfIntervals(void) {
    return(_listLen);
  };

  //  Returns the sum of the length of all intervals
  intervalNumber      sumOfLengths(void) {
    intervalNumber len = 0;
    uint32         i   = numberOfIntervals();

    if (i > 0)
      while (i--)
        len += _list[i].hi - _list[i].lo;

    return(len);
  };


  //  Populates an array with the intervals that are within the
  //  supplied interval.  Return
  //
  uint32      overlapping(intervalNumber    lo,
                          intervalNumber    hi,
                          uint32          *&intervals,
                          uint32           &intervalsLen,
                          uint32           &intervalsMax);

  //  Adds all the intervals in IL to this list
  //
  void        merge(intervalList *IL);

  //  Populates this intervalList with the intersection of A and B.
  //  This intervalList is not cleared prior to adding new intervals.
  //
  //  Both A and B call merge().
  //
  void                intersect(intervalList &A,
                                intervalList &B);

  //  Populates this intervalList with regions in A that are completely
  //  contained in a region in B.
  //
  //  Both A and B call merge().
  //
  void                contained(intervalList &A,
                                intervalList &B);


  intervalNumber      lo(uint32 i) { return(_list[i].lo); };
  intervalNumber      hi(uint32 i) { return(_list[i].hi); };
  uint32              ct(uint32 i) { return(_list[i].ct); };

private:
  bool                      _isSorted;
  bool                      _isMerged;
  uint32                    _listLen;
  uint32                    _listMax;
  _intervalPair            *_list;
};



//  Takes as input an intervalList, computes the number of intervals
//  covering every position in there, stores this as a new set of
//  intervals, annotated with the depth.
//
//  This is a static object, initialized once by the intervalList.
//
class intervalDepth {
public:
  intervalDepth(intervalList &IL);
  intervalDepth(intervalDepthRegions *id, uint32 idlen);
  ~intervalDepth();

  //  Returns the number of intervals
  uint32                     numberOfIntervals(void) {
    return(_listLen);
  };

  intervalNumber             lo(uint32 i) { return(_list[i].lo); };
  intervalNumber             hi(uint32 i) { return(_list[i].hi); };
  uint32                     de(uint32 i) { return(_list[i].de); };

private:
  void                       computeIntervals(intervalDepthRegions *id, uint32 idlen);

  uint32                    _listLen;
  uint32                    _listMax;
  _intervalDepth           *_list;
};





#endif  //  AS_UTL_INTERVALLIST_H

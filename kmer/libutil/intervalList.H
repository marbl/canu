#ifndef INTERVALLIST_H
#define INTERVALLIST_H


//  I don't want this to be public, but the sort function
//  in intervalList.C needs it so.
//
typedef u64bit  intervalNumber;

struct _intervalPair {
  intervalNumber    lo;
  intervalNumber    hi;
  u32bit            ct;
};

struct _intervalDepth {
  intervalNumber    lo;
  intervalNumber    hi;
  u32bit            de;
};

struct intervalDepthRegions {
  intervalNumber    pos;  //  Position of the change in depth
  s32bit            cha;  //  The change in depth, e.g, +1, -1
};


class intervalList {
public:
  intervalList(u32bit initialSize=32);
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
  void        merge(void);

  void        invert(intervalNumber lo, intervalNumber hi);

  //  Returns the number of intervals
  u32bit      numberOfIntervals(void) {
    return(_listLen);
  };

  //  Returns the sum of the length of all intervals
  intervalNumber      sumOfLengths(void) {
    intervalNumber len = 0;
    u32bit         i   = numberOfIntervals();

    if (i > 0)
      while (i--)
        len += _list[i].hi - _list[i].lo;

    return(len);
  };


  //  Populates an array with the intervals that are within the
  //  supplied interval.  Return
  //
  u32bit      overlapping(intervalNumber    lo,
                          intervalNumber    hi,
                          u32bit          *&intervals,
                          u32bit           &intervalsLen,
                          u32bit           &intervalsMax);

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


  intervalNumber      lo(u32bit i) { return(_list[i].lo); };
  intervalNumber      hi(u32bit i) { return(_list[i].hi); };
  u32bit              ct(u32bit i) { return(_list[i].ct); };

private:
  bool                      _isSorted;
  bool                      _isMerged;
  u32bit                    _listLen;
  u32bit                    _listMax;
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
  intervalDepth(intervalDepthRegions *id, u32bit idlen);
  ~intervalDepth();

  //  Returns the number of intervals
  u32bit                     numberOfIntervals(void) {
    return(_listLen);
  };

  intervalNumber             lo(u32bit i) { return(_list[i].lo); };
  intervalNumber             hi(u32bit i) { return(_list[i].hi); };
  u32bit                     de(u32bit i) { return(_list[i].de); };

private:
  void                       computeIntervals(intervalDepthRegions *id, u32bit idlen);

  u32bit                    _listLen;
  u32bit                    _listMax;
  _intervalDepth           *_list;
};




#endif  //  INTERVALLIST_H

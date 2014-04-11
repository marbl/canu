#ifndef INTERVALLIST_H
#define INTERVALLIST_H


//  I don't want this to be public, but the sort function
//  in intervalList.C needs it so.
//
typedef uint64  intervalNumber;

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
  int32            cha;  //  The change in depth, e.g, +1, -1
};


class intervalList {
public:
  intervalList(uint32 initialSize=32);
  ~intervalList();

  intervalList &operator=(intervalList &src);

  void        clear(void) {
    _isSorted = true;
    _isMerged = true;
    _listLen  = 0;
  }

  void        add(intervalNumber position, intervalNumber length);
  void        sort(void);
  void        merge(uint32 minOverlap=0);   //  Merge overlapping regions
  void        merge(intervalList *IL);      //  Insert IL into this list

  void        intersect(intervalList &A,
                        intervalList &B);

  uint32      overlapping(intervalNumber    lo,
                          intervalNumber    hi,
                          uint32          *&intervals,
                          uint32           &intervalsLen,
                          uint32           &intervalsMax);

  void        invert(intervalNumber lo, intervalNumber hi);

  uint32      numberOfIntervals(void)   { return(_listLen); };

  intervalNumber      sumOfLengths(void) {
    intervalNumber len = 0;
    uint32         i   = numberOfIntervals();

    if (i > 0)
      while (i--)
        len += _list[i].hi - _list[i].lo;

    return(len);
  };


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

  uint32                     numberOfIntervals(void)   { return(_listLen); };

  intervalNumber             lo(uint32 i) { return(_list[i].lo); };
  intervalNumber             hi(uint32 i) { return(_list[i].hi); };
  uint32                     de(uint32 i) { return(_list[i].de); };

private:
  void                       computeIntervals(intervalDepthRegions *id, uint32 idlen);

  uint32                    _listLen;
  uint32                    _listMax;
  _intervalDepth           *_list;
};


#endif  //  INTERVALLIST_H

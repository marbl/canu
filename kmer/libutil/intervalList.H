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

  void        clear(void) {
    _isSorted = true;
    _isMerged = true;
    _listLen  = 0;
  }

  void        add(intervalNumber position, intervalNumber length);
  void        sort(void);
  void        merge(u32bit minOverlap=0);   //  Merge overlapping regions
  void        merge(intervalList *IL);      //  Insert IL into this list

  void        intersect(intervalList &A,
                        intervalList &B);

  u32bit      overlapping(intervalNumber    lo,
                          intervalNumber    hi,
                          u32bit          *&intervals,
                          u32bit           &intervalsLen,
                          u32bit           &intervalsMax);

  void        invert(intervalNumber lo, intervalNumber hi);

  u32bit      numberOfIntervals(void)   { return(_listLen); };

  intervalNumber      sumOfLengths(void) {
    intervalNumber len = 0;
    u32bit         i   = numberOfIntervals();

    if (i > 0)
      while (i--)
        len += _list[i].hi - _list[i].lo;

    return(len);
  };


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

  u32bit                     numberOfIntervals(void)   { return(_listLen); };

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

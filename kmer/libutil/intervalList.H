#ifndef INTERVALLIST_H
#define INTERVALLIST_H


//  I don't want this to be public, but the sort function
//  in intervalList.C needs it so.
//
typedef u64bit  intervalNumber;

struct _intervalPair {
  intervalNumber    lo;
  intervalNumber    hi;
};


class intervalList {
public:
  intervalList();
  ~intervalList();

  intervalList &operator=(intervalList &src);

  //  Clear a list
  void        clear(void) {
    _isSorted = true;
    _listLen  = 0;
  }

  //  Insert a new interval into the list
  void        add(intervalNumber position, intervalNumber length);

  //  Sort the set of intervals by the lo value
  void        sort(void);

  //  Merge overlapping or adjacent intervals together.
  void        merge(void);

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

  intervalNumber      lo(u32bit i) { return(_list[i].lo); };
  intervalNumber      hi(u32bit i) { return(_list[i].hi); };

private:
  bool                      _isSorted;
  u32bit                    _listLen;
  u32bit                    _listMax;
  _intervalPair            *_list;
};


#endif  //  INTERVALLIST_H

#ifndef INTERVAL_LIST_H
#define INTERVAL_LIST_H

#include "libbritypes.h"
#include <stdlib.h>
#include <string.h>

//
//  An interval list -- 
//

//  I don't want this to be public, but the sort function
//  in intervalList.C needs it so.
//
struct _intervalPair {
  u32bit    lo;
  u32bit    hi;
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
  void        add(u32bit position, u32bit length);

  //  Sort the set of intervals by the lo value
  void        sort(void);

  //  Merge overlapping or adjacent intervals together.
  void        merge(void);

  //  Returns the number of intervals
  u32bit      numberOfIntervals(void) {
    return(_listLen);
  };

  //  Returns the sum of the length of all intervals
  u32bit      sumOfLengths(void) {
    u32bit len = 0;
    u32bit i   = numberOfIntervals();

    if (i > 0)
      while (i--)
        len += _list[i].hi - _list[i].lo;

    return(len);
  };

  u32bit      lo(u32bit i) {
    return(_list[i].lo);
  };
  u32bit      hi(u32bit i) {
    return(_list[i].hi);
  };

private:
  bool                      _isSorted;
  u32bit                    _listLen;
  u32bit                    _listMax;
  _intervalPair            *_list;
};


#endif  //  INTERVAL_LIST_H

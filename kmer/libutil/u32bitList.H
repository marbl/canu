#ifndef U32BITLIST_H
#define U32BITLIST_H

#include <string.h>

//  A very simple integer list.  Hopefully lighter weight than a
//  vector.

//  It might be useful to extend this to have 'undef' values,
//  and to allow shift(), pop().

class u32bitList {
public:
  u32bitList(u32bit max=16) {
    _len = 0;
    _max = max;
    _lst = new u32bit [_max];
  };
  ~u32bitList() {
    delete [] _lst;
  };

private:
  void     resize(u32bit idx) {
    if (idx >= _max) {
      _max *= 2;
      u32bit *L = new u32bit [_max];
      memcpy(L, _lst, sizeof(u32bit) * _len);
      delete [] _lst;
      _lst = L;
    }
    if (idx >= _len)
      _len = idx + 1;
  }

public:
  u32bit &operator[](u32bit idx) {
    resize(idx);
    return(_lst[idx]);
  }

  void    push(u32bit val) {
    resize(_len);
    _lst[_len++] = val;
  }

  u32bit  length(void) {
    return(_len);
  };

  void    clear(void) {
    _len = 0;
  }

private:
  u32bit  _len;
  u32bit  _max;
  u32bit *_lst;
};


#endif  //   U32BITLIST_H

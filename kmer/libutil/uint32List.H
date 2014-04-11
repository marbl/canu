#ifndef UINT32LIST_H
#define UINT32LIST_H

#include <string.h>

//  A very simple integer list.  Hopefully lighter weight than a
//  vector.

//  It might be useful to extend this to have 'undef' values,
//  and to allow shift(), pop().

class uint32List {
public:
  uint32List(uint32 max=16) {
    _len = 0;
    _max = max;
    _lst = new uint32 [_max];
  };
  ~uint32List() {
    delete [] _lst;
  };

private:
  void     resize(uint32 idx) {
    if (idx >= _max) {
      _max *= 2;
      uint32 *L = new uint32 [_max];
      memcpy(L, _lst, sizeof(uint32) * _len);
      delete [] _lst;
      _lst = L;
    }
    if (idx >= _len)
      _len = idx + 1;
  }

public:
  uint32 &operator[](uint32 idx) {
    resize(idx);
    return(_lst[idx]);
  }

  void    push(uint32 val) {
    resize(_len);
    _lst[_len++] = val;
  }

  uint32  length(void) {
    return(_len);
  };

  void    clear(void) {
    _len = 0;
  }

private:
  uint32  _len;
  uint32  _max;
  uint32 *_lst;
};


#endif  //   UINT32LIST_H

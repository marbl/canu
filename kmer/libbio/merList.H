#ifndef MER_LIST_H
#define MER_LIST_H

//  A simple list of mers.  Implemented as a list of lists.

class merList {
public:
  merList() {
    _ptrsMax = 8;
    _ptrsLen = 0;
    _ptrs    = new coord * [_ptrsMax];

    _mersWid = 12;
    _mersMax = u32bitMASK(_mersWid);
    _mersLen = 0;
    _ptrs[0] = new coord [_mersMax+1];
  };
  ~merList() {
    for (u32bit x=0; x<_ptrsLen+1; x++)
      delete [] _ptrs[x];
    delete [] _ptrs;
  };

  void  addMer(u32bit x, u32bit y) {
    if (_mersLen > _mersMax) {
      _ptrsLen++;

      if (_ptrsLen >= _ptrsMax) {
        _ptrsMax *= 2;
        coord **p = new coord * [_ptrsMax];
        memcpy(p, _ptrs, sizeof(coord*) * _ptrsLen);
        delete [] _ptrs;
        _ptrs = p;
      }

      _ptrs[_ptrsLen] = new coord [_mersMax+1];
      _mersLen = 0;
    }

    _ptrs[_ptrsLen][_mersLen]._qPos = x;
    _ptrs[_ptrsLen][_mersLen]._gPos = y;

    _mersLen++;
  };

  bool  getMer(u32bit i, u32bit &x, u32bit &y) {
    u32bit p = i >> _mersWid;
    u32bit a = i  & _mersMax;

    if ((p > _ptrsLen) || ((p == _ptrsLen) && (a >= _mersLen)))
      return(false);

    x = _ptrs[(i >> _mersWid)][i & _mersMax]._qPos;
    y = _ptrs[(i >> _mersWid)][i & _mersMax]._gPos;

    return(true);
  };

  void  clear(void) {
    // Don't delete the first guy!  We write into it blindly!
    for (u32bit x=1; x<_ptrsLen; x++)
      delete [] _ptrs[x];
    _ptrsLen = 0;
    _mersLen = 0;
  };

  void  merge(merList *ML) {
    u32bit i, x, y;
    for (i=0; ML->getMer(i, x, y); i++)
      addMer(x, y);
  };

private:
  struct coord {
    u32bit  _qPos;
    u32bit  _gPos;
  };

  //  The number of mer blocks we have space for, and the current mer
  //  block.
  //
  u32bit         _ptrsMax;
  u32bit         _ptrsLen;
  coord        **_ptrs;

  //  The number of mers available in each block, and the current mer
  //  we are at in the current block (for adding new mers).
  //
  u32bit         _mersWid;
  u32bit         _mersMax;
  u32bit         _mersLen;
};

#endif  //  MER_LIST_H

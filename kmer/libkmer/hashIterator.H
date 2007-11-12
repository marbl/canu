#ifndef HASHITERATOR_H
#define HASHITERATOR_H

#include "util++.H"


class hashIterator {
public:
  hashIterator(u32bit ts, u32bit ne) {
    _e1 = 0;  _s1 = 0;
    _e2 = 0;  _s2 = 2;
    _e3 = 0;  _s3 = 4;
    _e4 = 0;  _s4 = 6; 

    _tsiz = ts;
    _nerr = ne;

    _frst = 1;

    if (_nerr == 4) {
      _e4 = 0x00c0;
      _e3 = 0x0030;
      _e2 = 0x000c;
      _e1 = 0x0003;
      _stop = (_e1 << _tsiz) | _e3 | _e2 | _e1;
    }
    if (_nerr == 3) {
      _e3 = 0x0030;
      _e2 = 0x000c;
      _e1 = 0x0003;
      _stop = (_e1 << _tsiz) | _e2 | _e1;
    }
    if (_nerr == 2) {
      _e2 = 0x000c;
      _e1 = 0x0003;
      _stop = (_e1 << _tsiz) | _e1;
    }
    if (_nerr == 1) {
      _e1 = 0x0003;
      _stop = (_e1 << _tsiz);
    }
  };
  ~hashIterator() {};

  u64bit  error1(void) { return(_e1); };
  u64bit  error2(void) { return(_e2); };
  u64bit  error3(void) { return(_e3); };
  u64bit  error4(void) { return(_e4); };

  u32bit  shift1(void) { return(_s1); };
  u32bit  shift2(void) { return(_s2); };
  u32bit  shift3(void) { return(_s3); };
  u32bit  shift4(void) { return(_s4); };

  bool    next(void) {

    if (_frst) {
      _frst = 0;
      return(true);
    }

    _e1 <<= 2;
    _s1  += 2;
    if ((_nerr > 1) && (_e1 == _e2)) {
      _e1   = 0x0003;
      _s1   = 0;
      _e2 <<= 2;
      _s2  += 2;
      if ((_nerr > 2) && (_e2 == _e3)) {
        _e2   = 0x000c;
        _s2   = 2;
        _e3 <<= 2;
        _s3  += 2;
        if ((_nerr > 3) && (_e3 == _e4)) {
          _e3   = 0x0030;
          _s3   = 4;
          _e4 <<= 2;
          _s4  += 2;
        }
      }
    }

    return((_e1 | _e2 | _e3 | _e4) != _stop);
  };

private:
  u64bit  _e1;
  u64bit  _e2;
  u64bit  _e3;

  u64bit  _e4;

  u32bit  _s1;
  u32bit  _s2;
  u32bit  _s3;
  u32bit  _s4;

  u64bit  _stop;
  u64bit  _tsiz;
  u32bit  _nerr;
  u32bit  _frst;
};


#endif  //  HASHITERATOR_H

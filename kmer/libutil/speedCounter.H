#ifndef SPEEDCOUNTER_H
#define SPEEDCOUNTER_H

#include "libbri.H"

//
//  Format should be
//    "    %8f [unit]things (%8.5f [unit]things/sec)\r"
//

class speedCounter {
public:
  speedCounter(char const   *fmt,
               double        unit,
               u64bit        freq,
               bool          enabled=true) {
    _count     = 0;
    _unit      = unit;
    _freq      = freq;
    _startTime = getTime();
    _fmt       = fmt;
    _enabled   = enabled;

    //  Set all the bits below the hightest set in _freq --
    //  this allows us to do a super-fast test in tick().
    //
    _freq |= _freq >> 1;
    _freq |= _freq >> 2;
    _freq |= _freq >> 4;
    _freq |= _freq >> 8;
    _freq |= _freq >> 16;
    _freq |= _freq >> 32;
  };

  ~speedCounter() {
  };

  void   tick(void) {
    if (_enabled && ((++_count & _freq) == u64bitZERO)) {
      double  v = _count / _unit;
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fflush(stderr);
    }
  };

private:
  u64bit          _count;
  double          _unit;
  u64bit          _freq;
  double          _startTime;
  char const     *_fmt;
  bool            _enabled;
};


#endif  //  SPEEDCOUNTER_H

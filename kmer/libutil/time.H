#ifndef BRITIME_H
#define BRITIME_H

#include "libbritypes.h"

double  getTime(void);
void    write_rusage(FILE *F);


//  Format should be
//    "    %8f [unit]things (%8.5f [unit]things/sec)\r"
//
class speedCounter {
public:
  speedCounter(char const   *fmt,
               double        unit,
               u64bit        freq,
               bool          enabled=true);

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


#endif  //  BRITIME_H

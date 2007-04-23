#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <errno.h>

extern "C" {
#include "AS_global.h"
#include "AS_PER_gkpStore.h"
}

#include "util++.H"

void  doTrim(fragRecord *fr,
             double      minQuality,
             uint32     &left,
             uint32     &right);

void  findGoodQuality(double *qltD,
                      uint32  qltLen,
                      double  minQuality,
                      uint32 &qltL,
                      uint32 &qltR);

//  A simple initialized array -- performs a quality letter -> quality
//  value translation.
//
class qualityLookup {
public:
  qualityLookup() {
    for (uint32 i=0; i<='0'; i++)
      q[i] = 1.0;

    for (uint32 i='1'; i<255; i++)
      q[i] = 1 / pow(10, (i - '0') / 10.0);
  };
  ~qualityLookup() {
  };

  double lookupChar(char x) {
    return(q[(uint32)x]);
  };
  double lookupNumber(uint32 x) {
    return(q['0' + x]);
  };

private:
  double   q[255];
};


extern
qualityLookup  qual;


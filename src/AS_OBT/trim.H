#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <errno.h>

extern "C" {
#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
}

#include "util++.H"

void  doTrim(ReadStructp rd,
             double      minQuality,
             u32bit     &left,
             u32bit     &right);

void  findGoodQuality(double *qltD,
                      u32bit  qltLen,
                      double  minQuality,
                      u32bit &qltL,
                      u32bit &qltR);


//  A simple initialized array -- performs a quality letter -> quality
//  value translation.
//
class qualityLookup {
public:
  qualityLookup() {
    for (u32bit i=0; i<='0'; i++)
      q[i] = 1.0;
    for (u32bit i='1'; i<255; i++) {

      //  This is wrong, but consistent with what we had before now
      //q[i] = 1 / pow(10, log10((double)(i - '0')));

      //  This is right.
      q[i] = 1 / pow(10, (i - '0')/10.0 );
    }

    q[(u32bit)'0'] = 1.0;
  };
  ~qualityLookup() {
  };

  double lookupChar(char x) {
    return(q[(u32bit)x]);
  };
  double lookupNumber(u32bit x) {
    return(q['0' + x]);
  };

private:
  double   q[255];
};


extern
qualityLookup  qual;


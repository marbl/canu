#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#include <float.h>
#include <math.h>

#include <assert.h>
#include <errno.h>
#include <time.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/types.h>

typedef int8_t  int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t  uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

#include "AS_UTL_intervalList.H"
#include "AS_UTL_intervalList.C"

int
main(int argc, char **argv) {
  intervalList  L;

  L.add(0,10);
  L.add(1,10);
  L.add(2,10);
  L.add(3,10);

  L.add(20,20);
  L.add(20,20);
  L.add(20,20);

  L.add(60,0);

  L.add(60,20);

  for (uint32 dd=0; dd<L.numberOfIntervals(); dd++)
    fprintf(stderr, "inter[%d] %d-%d %d\n",
            dd, L.lo(dd), L.hi(dd), L.ct(dd));

  intervalDepth D(L);

  for (uint32 dd=0; dd<D.numberOfIntervals(); dd++)
    fprintf(stderr, "depth[%d] %d-%d %d\n",
            dd, D.lo(dd), D.hi(dd), D.de(dd));

  return(0);
}


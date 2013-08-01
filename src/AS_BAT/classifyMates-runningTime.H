
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

#ifndef CLASSIFYMATES_RUNNINGTIME_H
#define CLASSIFYMATES_RUNNINGTIME_H

static const char *rcsid_CLASSIFYMATES_RUNNINGTIME_H = "$Id$";

#include <stdio.h>
#include <math.h>

class onlineMeanStdDev {
public:
  onlineMeanStdDev() {
    hMin   = UINT32_MAX;
    hMax   = 0;
    hAlc   = 65536;
    hDat   = new uint32 [hAlc];

    memset(hDat, 0, sizeof(uint32) * hAlc);

    hCnt   = 0;

    computedMean   = 0.0;
    computedStdDev = 0.0;
  };
  ~onlineMeanStdDev() {
    delete [] hDat;
  };

  uint32    min(void)     { return(hMin);           };
  uint32    max(void)     { return(hMax);           };

  double    mean(void)    { return(computedMean);   };
  double    stddev(void)  { return(computedStdDev); };

  uint32    numData(void) { return(hCnt);           };

  void      addDataPoint(uint32 dp) {

    if (dp >= hAlc) {
      while (hAlc <= dp)
        hAlc *= 2;

      uint32  *hNew = new uint32 [hAlc];

      memset(hNew, 0,    sizeof(uint32) * (hAlc));
      memcpy(hNew, hDat, sizeof(uint32) * (hMax + 1));

      delete [] hDat;
      hDat = hNew;
    }

    if (hMin > dp)  hMin = dp;
    if (hMax < dp)  hMax = dp;

    hDat[dp]++;
    
    hCnt++;
  };

  void      recompute(void) {
    double m = 0;  //  Scratch copies; we update the real ones quickly at the end.
    double s = 0;

    uint32 c = 0;

    for (uint32 dp=hMin; dp<=hMax; dp++) {
      c += hDat[dp];
      m += hDat[dp] * dp;
    }

    if (c != hCnt)
      fprintf(stderr, "WARNING: c=%u hCnt=%u\n", c, hCnt);
    assert(c == hCnt);

    m /= hCnt;

    for (uint32 dp=hMin; dp<=hMax; dp++)
      s += hDat[dp] * (dp - m) * (dp - m);

    s /= hCnt - 1;
    s  = sqrt(s);

#if 0
    uint32    buckets[101] = {0};

    for (uint32 dp=0; dp<hMax; dp++) {
      uint32 b = dp * 100 / hMax;
      if (b > 100) b = 100;
      buckets[b] += hDat[dp];
    }

    errno = 0;
    FILE *F = fopen("distribution", "a");
    if (errno)
      fprintf(stderr, "FAILED to open 'distribution': %s\n", strerror(errno)), exit(1);
    for (uint32 bb=0; bb<100; bb++) {
      fprintf(F, "%u\t", buckets[bb]);
    }
    fprintf(F, "\n");
    fclose(F);
#endif

    computedMean   = m;
    computedStdDev = s;
  };

private:
  uint32    hMin;  //  Minimum value we've seen.
  uint32    hMax;  //  Maximum value we've seen.
  uint32    hAlc;  //  Maximum value we've allocated space for in hDat.
  uint32   *hDat;  //  Histogram of values.

  uint32    hCnt;  //  Number of datapoints - sum of hDat.

  double    computedMean;
  double    computedStdDev;
};

#endif  //  CLASSIFYMATES_RUNNINGTIME_H


/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

//  $Id: AS_UTL_histogram.h,v 1.4 2008-06-27 06:29:18 brianwalenz Exp $

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"


typedef struct {
  uint64   nSamples;
  uint64  *histogram;

  uint64   allocated;
  uint32   isMallocd;

  uint64   smallest;  //  The smallest/largest value inserted,
  uint64   largest;   //  not the smallest/largest count

  uint64   median;

  double   mean;
  double   stddev;

  uint64   mode;
  double   mad;
} AS_UTL_histogram;


typedef struct {
  uint64   nSamples;
  uint64 **histogram;
  uint64  *histogramData;

  uint64   allocated[2];  //  The X (aka [0]) dimension can be resized; the Y ([1]) cannot.
  uint32   isMallocd;

} AS_UTL_histogram3d;



AS_UTL_histogram *
AS_UTL_histogramAllocate(AS_UTL_histogram *h) {
  if (h == NULL) {
    h = safe_calloc(1, sizeof(AS_UTL_histogram));
    h->isMallocd = 1;
  }
  if (h->histogram == NULL) {
    h->allocated = 1048576;
    h->histogram = (uint64 *)safe_calloc(h->allocated, sizeof(uint64));
  }
  return(h);
}



void
AS_UTL_histogramFree(AS_UTL_histogram *h) {
  if (h == NULL)
    return;
  safe_free(h->histogram);
  if (h->isMallocd)
    safe_free(h);
}


void
AS_UTL_histogramAdd(AS_UTL_histogram *h, uint64 x) {

  h->nSamples++;

  if (x >= h->allocated) {
    h->allocated = x + x / 2;
    h->histogram = (uint64 *)safe_realloc(h->histogram, h->allocated * sizeof(uint64));
  }

  if (x < h->smallest)
    h->smallest = x;

  if (h->largest < x)
    h->largest = x;

  h->histogram[x]++;
}


void
AS_UTL_histogramCompute(AS_UTL_histogram *h) {
  uint32   i = 0;
  uint64   n = 0;
  uint64  *t = NULL;

  h->median  = 0;
  h->mean    = 0.0;
  h->stddev  = 0.0;
  h->mode    = 0;
  h->mad     = 0.0;

  if (h->nSamples != 0) {
     //  median - count up to nSamples/2
     //
     for (i=h->smallest, n=0; n < h->nSamples/2; i++)
       n += h->histogram[i];
     h->median = i - 1;

     //  mean - sum, divide by n
     //
     for (i=h->smallest; i <= h->largest; i++)
       h->mean += (double)h->histogram[i] * (double)i;
     h->mean /= h->nSamples;

     //  stddev -- we have h[i] values, and the value is i.
     //
     for (i=h->smallest; i <= h->largest; i++)
       h->stddev += h->histogram[i] * (h->mean - i) * (h->mean - i);
     h->stddev  = sqrt(h->stddev / h->nSamples);

     //  mode - just find the max
     //
     for (i=h->smallest, n=0; i <= h->largest; i++)
       if (h->histogram[n] < h->histogram[i])
         n = i;
     h->mode = n;

     //  mad - 1.4826*median(abs(median(v)-v))
     //
     //  really only need max(h->largest - h->median, h->median), I think
     //
     t = (uint64 *)safe_calloc(h->largest, sizeof(uint64));

     for (i=h->smallest; i <= h->largest; i++)
       t[(h->median < i) ? (i - h->median) : (h->median - i)] += h->histogram[i];
     for (i=0, n=0; n < h->nSamples/2; i++)
       n += t[i];
     h->mad = 1.4826 * (i - 1);
  }

  safe_free(t);
}


void
AS_UTL_histogramShow(AS_UTL_histogram *h, FILE *F, char *label) {
  int64   max, min;
  int64   rng, i;

  uint64  valFull[16][3] = {0};
  uint64  valZoom[16][3] = {0};
  uint64  valFullMax = 0;
  uint64  valZoomMax = 0;
  uint64  numBins = 16;  //  Note hardcoded in valFull and valZoom!

  fprintf(F, "\n");
  fprintf(F, "%s\n", label);
  fprintf(F, "nSamples      "F_U64"\n", h->nSamples);
  fprintf(F, "median        "F_U64"\n", h->median);
  fprintf(F, "mean/stddev   %.4f +- %.4f\n", h->mean, h->stddev);
  fprintf(F, "mode/mad      "F_U64" +- %.4f\n", h->mode, h->mad);

  //  A generic N-bucket histogram

  max = h->largest;
  min = h->smallest;
  rng = (max - min) / numBins + 1;
  i   = 0;
  while (min < max) {
    uint64  sum = 0;
    uint64  lim = min + rng;
    while (min < lim)
      sum += h->histogram[min++];
    valFull[i][0] = lim - rng;
    valFull[i][1] = lim - 1;
    valFull[i][2] = sum;
    i++;
  }
  valFullMax = i;

  //  Around the mode, we want to see some more detail

  max = h->mode + 3 * h->mad;
  min = h->mode - 3 * h->mad;

  if (h->mad < 1.0) {
    max = h->mode + 5;
    min = h->mode - 5;
  }

  if (max > (int64)h->largest)
    max = h->largest;

  if (min < (int64)h->smallest)
    min = h->smallest;

  rng = (max - min) / numBins + 1;
  i   = 0;

  //fprintf(stderr, "min %d max %d rng %d\n", min, max, rng);

  while (min < max) {
    uint64  sum = 0;
    uint64  lim = min + rng;
    while (min < lim)
      sum += h->histogram[min++];
    valZoom[i][0] = lim - rng;
    valZoom[i][1] = lim - 1;
    valZoom[i][2] = sum;
    i++;
  }
  valZoomMax = i;

  //  Finally, show it.

  for (i=0; i<valFullMax || i<valZoomMax; i++) {
    if (i < valFullMax) {
      fprintf(F, "%6llu - %6llu : %6llu    ",
              valFull[i][0], valFull[i][1], valFull[i][2]);
    } else {
      fprintf(F, "                            ");
    }

    if (i < valZoomMax) {
      fprintf(F, "%6llu - %6llu : %6llu\n",
              valZoom[i][0], valZoom[i][1], valZoom[i][2]);
    } else {
      fprintf(F, "\n");
    }
  }
}



void
AS_UTL_histogramDump(AS_UTL_histogram *h, char *filename, char *label) {
  FILE   *F;
  uint64  i;

  errno = 0;
  F = fopen(filename, "w");
  if (errno) {
    fprintf(stderr, "Failed to open '%s' for write: %s.\n", filename, strerror(errno));
    return;
  }

  //AS_UTL_histogramShow(h, F, label);

  for (i=h->smallest; i<=h->largest; i++)
    fprintf(F, F_U64"\t"F_U64"\n", i, h->histogram[i]);

  fclose(F);
}


////////////////////////////////////////////////////////////////////////////////

AS_UTL_histogram3d *
AS_UTL_histogram3dAllocate(AS_UTL_histogram3d *h, uint32 x, uint32 y) {
  if (h == NULL) {
    h = safe_calloc(1, sizeof(AS_UTL_histogram3d));
    h->isMallocd = 1;
  }
  if (h->histogramData == NULL) {
    int i;
    h->allocated[0] = x;
    h->allocated[1] = y;
    h->histogramData = (uint64  *)safe_calloc(h->allocated[0] * h->allocated[1], sizeof(uint64));
    h->histogram     = (uint64 **)safe_calloc(h->allocated[0], sizeof(uint64 *));
    for (i=0; i<h->allocated[0]; i++)
      h->histogram[i] = h->histogramData + i * h->allocated[1];
  }
  return(h);
}



void
AS_UTL_histogram3dFree(AS_UTL_histogram3d *h) {
  if (h == NULL)
    return;
  safe_free(h->histogramData);
  safe_free(h->histogram);
  if (h->isMallocd)
    safe_free(h);
}


void
AS_UTL_histogram3dAdd(AS_UTL_histogram3d *h, uint64 x, uint64 y) {

  h->nSamples++;

  if (x >= h->allocated[0]) {
    int i;
    h->allocated[0] = x * x / 2;
    h->histogramData = (uint64  *)safe_realloc(h->histogramData, h->allocated[0] * h->allocated[1] * sizeof(uint64));
    h->histogram     = (uint64 **)safe_realloc(h->histogram,     h->allocated[0] * sizeof(uint64 *));
    for (i=0; i<h->allocated[0]; i++)
      h->histogram[i] = h->histogramData + i * h->allocated[1];
  }
  if (y >= h->allocated[1]) {
    fprintf(stderr, "histogram3d: ("F_U64","F_U64") too big; only ("F_U64","F_U64") spots.\n",
            x, y, h->allocated[0], h->allocated[1]);
    return;
  }

#warning not remembering the smallest/largest

  h->histogram[x][y]++;
}


void
AS_UTL_histogram3dCompute(AS_UTL_histogram3d *h) {
}


void
AS_UTL_histogram3dShow(AS_UTL_histogram3d *h, char *label) {
}



void
AS_UTL_histogram3dDump(AS_UTL_histogram3d *h, char *filename, char *label) {
  FILE   *F;
  uint64  x, y;

  errno = 0;
  F = fopen(filename, "w");
  if (errno) {
    fprintf(stderr, "Failed to open '%s' for write: %s.\n", filename, strerror(errno));
    return;
  }

#if 1
  //  suitable for gnuplot 'splot'
  for (x=0; x<h->allocated[0]; x++) {
    for (y=0; y<h->allocated[1]; y++) {
      if (h->histogram[x][y] > 0)
        fprintf(F, F_U64" "F_U64" "F_U64"\n", x, y, h->histogram[x][y]);
    }
  }
#endif

#if 0
  //  Just a matrix dump.  gnuplot very slow!
  for (x=0; x<h->allocated[0]; x++) {
    for (y=0; y<h->allocated[1]; y++) {
      fprintf(F, F_U64" ", h->histogram[x][y]);
    }
    fprintf(F, "\n");
  }
#endif

#if 0
  //  Suitable for R
  for (y=0; y<h->allocated[1]; y++)
    fprintf(F, "\"%d\" ", y);
  fprintf(F, "\n");

  for (x=0; x<h->allocated[0]; x++) {
    fprintf(F, "\"%d\" ", x);
    for (y=0; y<h->allocated[1]; y++)
      fprintf(F, F_U64" ", h->histogram[x][y]);
    fprintf(F, "\n");
  }
#endif

  fclose(F);
}

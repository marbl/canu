
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
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

static const char *rcsid = "$Id: trim.C,v 1.9 2008-10-08 22:02:57 brianwalenz Exp $";

#include "trim.H"
#include "util++.H"

//  We're lazy; everyone needs to do quality letter -> quality value
//  translations.
//
qualityLookup   qual;


static
void
findGoodQuality(double  *qltD,
                uint32   qltLen,
                double   minQuality,
                uint32  &qltL,
                uint32  &qltR) {

  struct pair {
    uint32     start;
    uint32     end;
  };

  pair     f[AS_READ_MAX_LEN];  //  forward scan ranges
  pair     r[AS_READ_MAX_LEN];  //  reverse scan ranges

  uint32   fpos=0, flen=0;
  uint32   rpos=0, rlen=0;

  uint32   p = 0;
  double   q = 0;


  //  Scan forward, find first base with quality >= 20.
  //  Then, find the first base that makes cumulative expected #errors > 1/100
  //    #errors = 1 / 10 ^ log(q)
  //
  while (p < qltLen) {

    //  Find the next begin point
    //
    while ((p < qltLen) && (qltD[p] > minQuality))
      p++;

    //  Got a begin point!  Scan until the quality drops significantly.
    //
    f[fpos].start = p;
    f[fpos].end   = p;

    q = qltD[p];
    p++;

    while ((p < qltLen) && (q / (p - f[fpos].start) < minQuality)) {
      q += qltD[p];
      p++;
    }

    f[fpos].end = p;

    if (f[fpos].end - f[fpos].start > 10)
      fpos++;
  }

  //  Scan backward, just like the forward.
  //
  //  Stung by using uint32 for p; p is one more than it wants to be.
  //  Although, to be fair, there are just about as many cases of p+1
  //  as p below.
  //
  p = qltLen;
  q = 0;

  while (p > 0) {
    while ((p > 0) && (qltD[p-1] > minQuality))
      p--;

    r[rpos].start = p;
    r[rpos].end   = p;

    if (p > 0) {
      p--;
      q = qltD[p];

      while ((p > 0) && (q / (r[rpos].end - p) < minQuality)) {
        p--;
        q += qltD[p];
      }

      r[rpos].start = p;

      if (r[rpos].end - r[rpos].start > 10)
        rpos++;
    }
  }


  //  Now, just pick the largest overlap

  qltL = 0;
  qltR = 0;

  flen = fpos;
  rlen = rpos;

  //fprintf(stderr, "qltLen = "uint32FMT"  flen="uint32FMT"  rlen="uint32FMT"\n", qltLen, flen, rlen);

  uint32   winningFPos  = 0;
  uint32   winningRPos  = 0;
  uint32   winningStyle = 0;

  for (fpos=0; fpos<flen; fpos++) {
    for (rpos=0; rpos<rlen; rpos++) {

      //  Not all cases are needed.  Easier to take care of them all,
      //  than to figure out what the minimal tests are.

      //        fffffffffff
      //   rrrrrrrrrr
      //
      if ((r[rpos].start <= f[fpos].start) &&
          (f[fpos].start <= r[rpos].end) &&
          (r[rpos].end   <= f[fpos].end)) {
        if ((r[rpos].end - f[fpos].start) > (qltR - qltL)) {
          winningFPos  = fpos;
          winningRPos  = rpos;
          winningStyle = 0;
          qltL = f[fpos].start;
          qltR = r[rpos].end;
        }
      }

      //   fffffffffff
      //        rrrrrrrrrr
      //
      else if ((f[fpos].start <= r[rpos].start) &&
               (r[rpos].start <= f[fpos].end) &&
               (f[fpos].end   <= r[rpos].end)) {
        if ((f[fpos].end - r[rpos].start) > (qltR - qltL)) {
          winningFPos  = fpos;
          winningRPos  = rpos;
          winningStyle = 1;
          qltL = r[rpos].start;
          qltR = f[fpos].end;
        }
      }

      //   fffffffffffffffffff
      //        rrrrrrrrrr
      //
      else if ((f[fpos].start <= r[rpos].start) &&
               (r[rpos].end   <= f[fpos].end)) {
        if ((r[rpos].end - r[rpos].start) > (qltR - qltL)) {
          winningFPos  = fpos;
          winningRPos  = rpos;
          winningStyle = 2;
          qltL = r[rpos].start;
          qltR = r[rpos].end;
        }
      }

      //        fffffffffff
      //   rrrrrrrrrrrrrrrrrrrr
      //
      else if ((r[rpos].start <= f[fpos].start) &&
               (f[fpos].end   <= r[rpos].end)) {
        if ((f[fpos].end - f[fpos].start) > (qltR - qltL)) {
          winningFPos  = fpos;
          winningRPos  = rpos;
          winningStyle = 3;
          qltL = f[fpos].start;
          qltR = f[fpos].end;
        }
      }

      else if (f[fpos].end < r[rpos].start) {
        //  NOP, no overlap.
      }

      else if (r[rpos].end < f[fpos].start) {
        //  NOP, no overlap.
      }

      else {
        fprintf(stderr, "UNMATCHED OVERLAP\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\n",
                f[fpos].start, f[fpos].end, r[rpos].start, r[rpos].end);
      }
    }
  }
}



//  Takes a fragRecord, returns clear ranges for some quality
//  threshold.  Higher level than I wanted, but it obscures
//  everything, and is exactly the interface that this and
//  mergeTrimming.C want.
//
void
doTrim(fragRecord *fr, double minQuality, uint32 &left, uint32 &right) {
  static double  qltD[AS_READ_MAX_LEN];
  char          *qltC   = getFragRecordQuality(fr);
  uint32         qltLen = getFragRecordQualityLength(fr);

  for (uint32 i=0; i<qltLen; i++)
    qltD[i] = qual.lookupChar(qltC[i]);

  findGoodQuality(qltD, qltLen, minQuality, left, right);
}



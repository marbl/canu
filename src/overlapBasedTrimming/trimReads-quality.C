
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2015-NOV-23
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "trimReads.H"


//  This is mostly historical.
//
//  It needs to be updated to use sanger qv encoding, and to take that as input.


//  A simple initialized array -- performs a quality letter -> quality
//  value translation.
//
class qualityLookup {
public:
  qualityLookup() {
    for (uint32 i=0; i<255; i++)
      q[i] = 1 / pow(10, i / 10.0);
  };
  ~qualityLookup() {
  };

  double lookup(uint32 x)  { return(q[x]);          };

private:
  double   q[255];
};


qualityLookup   qual;




static
void
findGoodQuality(double  *qltD,
                uint32   qltLen,
                char     minQualityLetter,
                uint32  &qltL,
                uint32  &qltR) {

  struct pair {
    uint32     start;
    uint32     end;
  };

  pair    *f = new pair [qltLen + 1];
  pair    *r = new pair [qltLen + 1];

  uint32   fpos=0, flen=0;
  uint32   rpos=0, rlen=0;

  uint32   p = 0;
  double   q = 0;

  double   minQuality = qual.lookup(minQualityLetter - '!');


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

  //fprintf(stderr, "qltLen = "F_U32"  flen="F_U32"  rlen="F_U32"\n", qltLen, flen, rlen);

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

  delete [] f;
  delete [] r;
}



//  Takes a gkFragment, returns clear ranges for some quality
//  threshold.  Higher level than I wanted, but it obscures
//  everything, and is exactly the interface that this and
//  mergeTrimming.C want.
//
void
doTrim(gkRead      *read,
       gkReadData  *readData,
       double       minQuality,
       uint32      &left,
       uint32      &right) {
  uint32    qltLen = read->gkRead_sequenceLength();
  char     *qltC   = readData->gkReadData_getQualities();
  double   *qltD   = new double [qltLen];

  for (uint32 i=0; i<qltLen; i++)
    qltD[i] = qual.lookup(qltC[i]);

  findGoodQuality(qltD, qltLen, minQuality, left, right);

  delete [] qltD;
}



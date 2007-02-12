#include "trim.H"

//  We're lazy; everyone needs to do quality letter -> quality value
//  translations.
//
qualityLookup   qual;


//  Takes a ReadStructp, returns clear ranges for some quality
//  threshold.  Higher level than I wanted, but it obscures
//  everything, and is exactly the interface that this and
//  mergeTrimming.C want.
//
void
doTrim(ReadStructp rd, double minQuality, u32bit &left, u32bit &right) {
  u32bit  seqMax = 10240;
  u32bit  qltLen = 0;
  char   *seq    = new char   [seqMax];
  char   *qltC   = new char   [seqMax];
  double *qltD   = new double [seqMax];

  if (getSequence_ReadStruct(rd, seq, qltC, seqMax)) {
    fprintf(stderr, "getSequence_ReadStruct() failed.\n");
    exit(1);
  }

  //  XXX  Probably a better way to do this....
  qltLen = strlen(qltC);

  for (u32bit i=0; i<qltLen; i++)
    qltD[i] = qual.lookupChar(qltC[i]);

  findGoodQuality(qltD, qltLen, minQuality, left, right);

  delete [] seq;
  delete [] qltC;
  delete [] qltD;
}



void
findGoodQuality(double  *qltD,
                u32bit   qltLen,
                double   minQuality,
                u32bit  &qltL,
                u32bit  &qltR) {

  struct pair {
    u32bit     start;
    u32bit     end;
  };

  pair     f[2048];  //  forward scan ranges
  pair     r[2048];  //  reverse scan ranges

  u32bit   fpos=0, flen=0;
  u32bit   rpos=0, rlen=0;

  u32bit   p = 0;
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
  //  Stung by using u32bit for p; p is one more than it wants to be.
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

  //fprintf(stderr, "qltLen = "u32bitFMT"  flen="u32bitFMT"  rlen="u32bitFMT"\n", qltLen, flen, rlen);

  u32bit   winningFPos  = 0;
  u32bit   winningRPos  = 0;
  u32bit   winningStyle = 0;

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
        fprintf(stderr, "UNMATCHED OVERLAP "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
                f[fpos].start, f[fpos].end, r[rpos].start, r[rpos].end);
      }
    }
  }

#if 0
  //  Statistics - do we care?  Not really.  It's mostly for debugging stuff.
  //  (the variables should be globals...)
  //
  u32bit   overlapA = 0;  //  The winning style of overlap
  u32bit   overlapB = 0;
  u32bit   overlapC = 0;
  u32bit   overlapD = 0;

  switch (winningStyle) {
    case 0:
      overlapA++;
      break;
    case 1:
      fprintf(stderr, "overlapB:  "u32bitFMTW(4)"-"u32bitFMTW(4)"  "u32bitFMTW(4)"-"u32bitFMTW(4)"\n",
              f[winningFPos].start, f[winningFPos].end,
              r[winningRPos].start, r[winningRPos].end);
      overlapB++;
      break;
    case 2:
      fprintf(stderr, "overlapC:  "u32bitFMTW(4)"-"u32bitFMTW(4)"  "u32bitFMTW(4)"-"u32bitFMTW(4)"\n",
              f[winningFPos].start, f[winningFPos].end,
              r[winningRPos].start, r[winningRPos].end);
      overlapC++;
      break;
    case 3:
      fprintf(stderr, "overlapD:  "u32bitFMTW(4)"-"u32bitFMTW(4)"  "u32bitFMTW(4)"-"u32bitFMTW(4)"\n",
              f[winningFPos].start, f[winningFPos].end,
              r[winningRPos].start, r[winningRPos].end);
      overlapD++;
      break;
    default:
      fprintf(stderr, "OVERLAP HAD NO WINNING STYLE?\n");
      break;
  }
#endif
}

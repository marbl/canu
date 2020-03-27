
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "ovStore.H"
#include "sqStore.H"


sqStore *ovOverlap::g = NULL;

//  Even though the b_end_hi | b_end_lo is uint64 in the struct, the result
//  of combining them doesn't appear to be 64-bit.  The cast is necessary.

char *
ovOverlap::toString(char                  *str,
                    ovOverlapDisplayType   type,
                    bool                   newLine) {

  switch (type) {
    case ovOverlapAsHangs:
      sprintf(str, "%10" F_U32P " %10" F_U32P "  %c  %6" F_S32P " %6" F_U32P " %6" F_S32P "  %7.6f%s%s",
              a_iid, b_iid,
              flipped() ? 'I' : 'N',
              a_hang(), span(), b_hang(),
              erate(),
              (overlapIsDovetail()) ? "" : "  PARTIAL",
              (newLine) ? "\n" : "");
      break;

    case ovOverlapAsCoords:
      sprintf(str, "%10" F_U32P " %10" F_U32P "  %c  %6" F_U32P "  %6" F_U32P " %6" F_U32P "  %6" F_U32P " %6" F_U32P "  %7.6f%s",
              a_iid, b_iid,
              flipped() ? 'I' : 'N',
              span(),
              a_bgn(), a_end(),
              b_bgn(), b_end(),
              erate(),
              (newLine) ? "\n" : "");
      break;

    case ovOverlapAsUnaligned:
      sprintf(str, "%10" F_U32P " %10" F_U32P "  %c  %6" F_U32P "  %6" F_OVP " %6" F_OVP "  %6" F_OVP " %6" F_OVP "  %7.6f %s %s %s%s",
              a_iid, b_iid,
              flipped() ? 'I' : 'N',
              span(),
              dat.ovl.ahg5, dat.ovl.ahg3,
              dat.ovl.bhg5, dat.ovl.bhg3,
              erate(),
              dat.ovl.forOBT ? "OBT" : "   ",
              dat.ovl.forDUP ? "DUP" : "   ",
              dat.ovl.forUTG ? "UTG" : "   ",
              (newLine) ? "\n" : "");
      break;

    case ovOverlapAsPaf:
      // miniasm/map expects entries to be separated by tabs
      // no padding spaces on names we don't confuse read identifiers
      sprintf(str, "%" F_U32P "\t%6" F_U32P "\t%6" F_U32P "\t%6" F_U32P "\t%c\t%" F_U32P "\t%6" F_U32P "\t%6" F_U32P "\t%6" F_U32P "\t%6" F_U32P "\t%6" F_U32P "\t%6" F_U32P " \tdv:f:%6.4f%s",
              a_iid,
              g->sqStore_getReadLength(a_iid), a_bgn(), a_end(),
              flipped() ? '-' : '+',
              b_iid,
              g->sqStore_getReadLength(b_iid), flipped() ? b_end() : b_bgn(), flipped() ? b_bgn() : b_end(),
              (uint32)floor(span() == 0 ? ((1-erate()) * (a_end()-a_bgn())) : (1-erate()) * span()),
              span() == 0 ? a_end() - a_bgn() : span(),
              255, erate(),
              (newLine) ? "\n" : "");
      break;
  }

  return(str);
}



bool
ovOverlap::fromString(splitToWords          &W,
                      ovOverlapDisplayType   type) {


  switch (type) {
    case ovOverlapAsHangs:
      a_iid = W.touint32(0);
      b_iid = W.touint32(1);

      flipped(W[2][0] == 'I');

      a_hang(W.touint32(3));
      span(W.touint32(4));
      b_hang(W.touint32(5));

      erate(W.todouble(6));

      break;

    case ovOverlapAsCoords:
      a_iid = W.touint32(0);
      b_iid = W.touint32(1);

      flipped(W[2][0] == 'I');

      span(W.touint32(3));

      {
        uint32  alen = g->sqStore_getReadLength(a_iid);
        uint32  blen = g->sqStore_getReadLength(b_iid);

        uint32  abgn = W.touint32(4);
        uint32  aend = W.touint32(5);

        uint32  bbgn = W.touint32(6);
        uint32  bend = W.touint32(7);

        dat.ovl.ahg5 = abgn;
        dat.ovl.ahg3 = alen - aend;

        dat.ovl.bhg5 = (dat.ovl.flipped) ? blen - bbgn :        bbgn;
        dat.ovl.bhg3 = (dat.ovl.flipped) ?        bend : blen - bend;
      }

      erate(W.todouble(8));
      break;

    case ovOverlapAsUnaligned:
      a_iid = W.touint32(0);
      b_iid = W.touint32(1);

      flipped(W[2][0] == 'I');

      dat.ovl.span = W.touint32(3);

      dat.ovl.ahg5 = W.touint32(4);
      dat.ovl.ahg3 = W.touint32(5);

      dat.ovl.bhg5 = W.touint32(6);
      dat.ovl.bhg3 = W.touint32(7);

      erate(W.todouble(8));

      dat.ovl.forUTG = false;
      dat.ovl.forOBT = false;
      dat.ovl.forDUP = false;

      for (uint32 i = 9; i < W.numWords(); i++) {
        dat.ovl.forUTG |= ((W[i][0] == 'U') && (W[i][1] == 'T') && (W[i][2] == 'G'));  //  Fails if W[i] == "U".
        dat.ovl.forOBT |= ((W[i][0] == 'O') && (W[i][1] == 'B') && (W[i][2] == 'T'));
        dat.ovl.forDUP |= ((W[i][0] == 'D') && (W[i][1] == 'U') && (W[i][2] == 'P'));
      }
      break;

    case ovOverlapAsPaf:
      a_iid = W.touint32(0);
      b_iid = W.touint32(5);

      flipped(W[4][0] == '-');

      dat.ovl.span = W.touint32(3)-W.toint32(2);

      uint32  alen = W.toint32(1);
      uint32  blen = W.toint32(6);

      uint32  abgn = W.touint32(2);
      uint32  aend = W.touint32(3);

      // paf looks like our coord format but the start/end aren't decreasing like ours so flip them then we can use the same math below
      uint32  bbgn = (dat.ovl.flipped) ? W.touint32(8) : W.touint32(7);
      uint32  bend = (dat.ovl.flipped) ? W.touint32(7) : W.touint32(8);

      dat.ovl.ahg5 = abgn;
      dat.ovl.ahg3 = alen - aend;

      dat.ovl.bhg5 = (dat.ovl.flipped) ? blen-bbgn : bbgn;
      dat.ovl.bhg3 = (dat.ovl.flipped) ?      bend : blen - bend;

      for (int i = 0; i < W.numWords(); i++) {
         if (W[i][0] == 'd' && W[i][1] == 'v' && W[i][3] == 'f') {
            erate(strtodouble(W[i]+5));
            break;
          }
      }

      dat.ovl.forUTG = true;
      dat.ovl.forOBT = true;
      dat.ovl.forDUP = true;
      break;
  }

#warning NEEDS TO RETURN FALSE IF FAILED TO DECODE OVERLAP STRING

  return(true);
}



void
ovOverlap::swapIDs(ovOverlap const &orig) {

  a_iid = orig.b_iid;
  b_iid = orig.a_iid;

  //  Copy the overlap as is, then fix it for the ID swap.

  for (uint32 ii=0; ii<ovOverlapNWORDS; ii++)
    dat.dat[ii] = orig.dat.dat[ii];

  //  Swap the A and B hangs.  If the overlap is flipped, we also need to reverse 5' and 3' hangs to
  //  make the now-A read forward oriented.

  if (orig.dat.ovl.flipped == false) {
    dat.ovl.ahg5 = orig.dat.ovl.bhg5;
    dat.ovl.ahg3 = orig.dat.ovl.bhg3;
    dat.ovl.bhg5 = orig.dat.ovl.ahg5;
    dat.ovl.bhg3 = orig.dat.ovl.ahg3;
  } else {
    dat.ovl.ahg5 = orig.dat.ovl.bhg3;
    dat.ovl.ahg3 = orig.dat.ovl.bhg5;
    dat.ovl.bhg5 = orig.dat.ovl.ahg3;
    dat.ovl.bhg3 = orig.dat.ovl.ahg5;
  }

  //  Whatever alignment orientation was in the original, it is opposite now.

#ifndef DO_NOT_STORE_ALIGN_PTR
  dat.ovl.alignSwapped = ! orig.dat.ovl.alignSwapped;
#endif
}

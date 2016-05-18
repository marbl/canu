
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
 *    Brian P. Walenz from 2014-DEC-15 to 2015-JUL-01
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Sergey Koren beginning on 2016-MAR-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "ovStore.H"
#include "gkStore.H"

//  Even though the b_end_hi | b_end_lo is uint64 in the struct, the result
//  of combining them doesn't appear to be 64-bit.  The cast is necessary.

char *
ovOverlap::toString(char                  *str,
                    ovOverlapDisplayType   type,
                    bool                   newLine) {

  switch (type) {
    case ovOverlapAsHangs:
      sprintf(str, "%10"F_U32P" %10"F_U32P"  %c  %6"F_S32P" %6"F_U32P" %6"F_S32P"  %7.6f%s%s",
              a_iid, b_iid,
              flipped() ? 'I' : 'N',
              a_hang(), span(), b_hang(),
              erate(),
              (overlapIsDovetail()) ? "" : "  PARTIAL",
              (newLine) ? "\n" : "");
      break;

    case ovOverlapAsCoords:
      sprintf(str, "%10"F_U32P" %10"F_U32P"  %c  %6"F_U32P"  %6"F_U32P" %6"F_U32P"  %6"F_U32P" %6"F_U32P"  %7.6f%s",
              a_iid, b_iid,
              flipped() ? 'I' : 'N',
              span(),
              a_bgn(), a_end(),
              b_bgn(), b_end(),
              erate(),
              (newLine) ? "\n" : "");
      break;

    case ovOverlapAsRaw:
      sprintf(str, "%10"F_U32P" %10"F_U32P"  %c  %6"F_U32P"  %6"F_U64P" %6"F_U64P"  %6"F_U64P" %6"F_U64P"  %7.6f %s %s %s%s",
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

    case ovOverlapAsCompat:
      sprintf(str, "%8"F_U32P" %8"F_U32P"  %c  %6d  %6d  %5.2f  %5.2f%s",
              a_iid,
              b_iid,
              dat.ovl.flipped ? 'I' : 'N',
              a_hang(), b_hang(),
              erate() * 100.0,
              erate() * 100.0,
              (newLine) ? "\n" : "");
      break;
    case ovOverlapAsPaf:
      // miniasm/map expects entries to be separated by tabs
      // no padding spaces on names we don't confuse read identifiers
      sprintf(str, "%"F_U32P"\t%6"F_U32P"\t%6"F_U32P"\t%6"F_U32P"\t%c\t%"F_U32P"\t%6"F_U32P"\t%6"F_U32P"\t%6"F_U32P"\t%6"F_U32P"\t%6"F_U32P"\t%6"F_U32P" %s",
              a_iid,
              (g->gkStore_getRead(a_iid)->gkRead_sequenceLength()), a_bgn(), a_end(),
              flipped() ? '-' : '+',
              b_iid,
              (g->gkStore_getRead(b_iid)->gkRead_sequenceLength()), flipped() ? b_end() : b_bgn(), flipped() ? b_bgn() : b_end(),
              (uint32)floor(span() == 0 ? (1-erate() * (a_end()-a_bgn())) : (1-erate()) * span()),
              span() == 0 ? a_end() - a_bgn() : span(),
              255,
              (newLine) ? "\n" : "");
      break;

  }

  return(str);
}



void
ovOverlap::swapIDs(ovOverlap const &orig) {

  a_iid = orig.b_iid;
  b_iid = orig.a_iid;

  //  Copy the overlap as is, then fix it for the ID swap.

#if (ovOverlapNWORDS == 5)
  dat.dat[0] = orig.dat.dat[0];
  dat.dat[1] = orig.dat.dat[1];
  dat.dat[2] = orig.dat.dat[2];
  dat.dat[3] = orig.dat.dat[3];
  dat.dat[4] = orig.dat.dat[4];
#else
  dat.dat[0] = orig.dat.dat[0];
  dat.dat[1] = orig.dat.dat[1];
  dat.dat[2] = orig.dat.dat[2];
#endif

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

  dat.ovl.alignSwapped = ! orig.dat.ovl.alignSwapped;
}

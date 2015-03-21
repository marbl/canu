

#include "ovStore.H"
#include "gkStore.H"

//  Even though the b_end_hi | b_end_lo is uint64 in the struct, the result
//  of combining them doesn't appear to be 64-bit.  The cast is necessary.

char *
ovsOverlap::toString(char *str, gkStore *gkp, bool asCoords) {

#if 0

  //  Compatible with CA 8.2
  sprintf(str, "%8"F_U32P" %8"F_U32P"  %c  %6d  %6d  %4.2f  %4.2f",
          a_iid,
          b_iid,
          dat.ovl.flipped ? 'I' : 'N',
          a_hang(), b_hang(),
          AS_OVS_decodeQuality(dat.ovl.erate) * 100.0,
          AS_OVS_decodeQuality(dat.ovl.erate) * 100.0);

#else

  if (asCoords)
    sprintf(str, "%10"F_U32P" %10"F_U32P"  %c  %6"F_S32P"  %6"F_U32P" %6"F_U32P"  %6"F_U32P" %6"F_U32P"  %6.3f\n",
            a_iid, b_iid,
            flipped() ? 'I' : 'N',
            span(),
            a_bgn(gkp), a_end(gkp),
            b_bgn(gkp), b_end(gkp),
            erate());
  else
    sprintf(str, "%10"F_U32P" %10"F_U32P"  %c  %6"F_S32P" %6"F_S32P" %6"F_S32P"  %6.3f%s\n",
            a_iid, b_iid,
            flipped() ? 'I' : 'N',
            a_hang(), span(), b_hang(),
            erate(),
            (overlapIsDovetail()) ? "" : "  PARTIAL");

#endif

  return(str);
}



void
ovsOverlap::swapIDs(ovsOverlap const &orig) {

  a_iid = orig.b_iid;
  b_iid = orig.a_iid;

  //  Copy the overlap as is, then fix it for the ID swap.

#if (ovsOverlapNWORDS == 5)
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

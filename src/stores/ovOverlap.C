

#include "ovStore.H"


  //  Even though the b_end_hi | b_end_lo is uint64 in the struct, the result
  //  of combining them doesn't appear to be 64-bit.  The cast is necessary.

char *
ovsOverlap::toString(char *str) {

  sprintf(str, "%8"F_U32P" %8"F_U32P"  %c  %6"F_OVP"  %6"F_OVP" %6"F_OVP"  %6"F_OVP" %6"F_OVP"  %4.2f",
          a_iid,
          b_iid,
          dat.ovl.flipped ? 'I' : 'N',
          dat.ovl.span,
          dat.ovl.ahg5, dat.ovl.ahg3,
          dat.ovl.bhg5, dat.ovl.bhg3,
          AS_OVS_decodeQuality(dat.ovl.erate) * 100.0);

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

  //  If the overlap is flipped, we also need to reverse 5' and 3' hangs to make the now-A read
  //  forward oriented.  Otherwise, we just need to swap the A and B hangs.

  if (orig.dat.ovl.flipped) {
    dat.ovl.ahg5 = orig.dat.ovl.bhg3;
    dat.ovl.ahg3 = orig.dat.ovl.bhg5;
    dat.ovl.bhg5 = orig.dat.ovl.ahg3;
    dat.ovl.bhg3 = orig.dat.ovl.ahg5;
  } else {
    dat.ovl.ahg5 = orig.dat.ovl.bhg5;
    dat.ovl.ahg3 = orig.dat.ovl.bhg3;
    dat.ovl.bhg5 = orig.dat.ovl.ahg5;
    dat.ovl.bhg3 = orig.dat.ovl.ahg3;
  }

  //  Whatever alignment orientation was in the original, it is opposite now.

  dat.ovl.alignSwapped = ! orig.dat.ovl.alignSwapped;
}

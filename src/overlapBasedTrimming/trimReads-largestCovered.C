
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
 *    Brian P. Walenz from 2015-MAY-28 to 2015-JUN-16
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "trimReads.H"

#include "intervalList.H"


bool
largestCovered(ovOverlap  *ovl,
               uint32       ovlLen,
               gkRead      *read,
               uint32       ibgn,
               uint32       iend,
               uint32      &fbgn,
               uint32      &fend,
               char        *logMsg,
               uint32       errorValue,
               uint32       minOverlap,
               uint32       minCoverage,
               uint32       minReadLength) {

  logMsg[0] = 0;

  assert(read->gkRead_readID() == ovl[0].a_iid);
  assert(ovlLen > 0);

  intervalList<uint32>  IL;
  intervalList<uint32>  ID;
  int32                 iid = read->gkRead_readID();

  uint32                nSkip = 0;
  uint32                nUsed = 0;

  for (uint32 i=0; i<ovlLen; i++) {
    uint32 tbgn = ovl[i].a_bgn();
    uint32 tend = ovl[i].a_end();

    assert(tbgn < tend);
    assert(iid == ovl[i].a_iid);

    if (ovl[i].evalue() > errorValue) {
      //  Overlap is crappy.
      //fprintf(stderr, "skip %2u\n", i);
      nSkip++;
      continue;
    }

    //fprintf(stderr, "save %2u\n", i);
    nUsed++;

    IL.add(tbgn, tend - tbgn);
  }

#if 0
  for (uint32 it=0; it<IL.numberOfIntervals(); it++)
    fprintf(stderr, "IL - %d - "F_S64" "F_S64" "F_S64"\n", fr.gkFragment_getReadIID(), IL.lo(it), IL.hi(it), IL.ct(it));

  for (uint32 it=0; it<ID.numberOfIntervals(); it++)
    fprintf(stderr, "ID - %d - "F_S64" "F_S64" "F_S64"\n", fr.gkFragment_getReadIID(), ID.lo(it), ID.hi(it), ID.de(it));
#endif

  //  I thought I'd allow low coverage at the end of the read, but not internally, but that is hard,
  //  and largely unnecessary.  We'll just not be assembling at low coverage joins, which is
  //  acceptable.

  if (minCoverage > 0) {
    intervalList<uint32>  DE(IL);

    uint32  it = 0;
    uint32  ib = 0;
    uint32  ie = 0;

    while (it < DE.numberOfIntervals()) {
      //fprintf(stderr, "DE - %d - "F_S64" "F_S64" "F_U32"\n", fr.gkFragment_getReadIID(), DE.lo(it), DE.hi(it), DE.depth(it));

      if (DE.depth(it) < minCoverage) {
        //  Dropped below good coverage depth.  If we have an interval, save it.  Reset.
        if (ie > ib) {
          //fprintf(stderr, "AD1 %d-%d len %d\n", ib, ie, ie - ib);
          ID.add(ib, ie - ib);
        }
        ib = 0;
        ie = 0;

      } else if ((ib == 0) && (ie == 0)) {
        //  Depth is good.  If no current interval, make a new one.
        ib = DE.lo(it);
        ie = DE.hi(it);
        //fprintf(stderr, "NE1 %d-%d len %d\n", ib, ie, ie - ib);

      } else if (ie == DE.lo(it)) {
        //  Depth is good.  If this interval is adjacent to the current, extend.
        ie = DE.hi(it);
        //fprintf(stderr, "EXT %d-%d len %d\n", ib, ie, ie - ib);

      } else {
        //  Depth is good, but we just had a gap in coverage.  Save any current interval.  Reset.
        if (ie > ib) {
          //fprintf(stderr, "AD2 %d-%d len %d\n", ib, ie, ie - ib);
          ID.add(ib, ie - ib);
        }
        ib = DE.lo(it);
        ie = DE.hi(it);
        //fprintf(stderr, "NE2 %d-%d len %d\n", ib, ie, ie - ib);
      }

      it++;
    }

    if (ie > ib) {
      //fprintf(stderr, "AD3 %d-%d len %d\n", ib, ie, ie - ib);
      ID.add(ib, ie - ib);
    }
  }

  //  Now that we've created depth, merge the intervals.

  IL.merge(minOverlap);

  //  IL - covered interavls enforcing a minimum overlap size (these can overlap)
  //  ID - covered intervals enforcing a minimum depth (these cannot overlap)
  //
  //  Create new intervals from the intersection of IL and ID.
  //
  //  This catches one nasty case, where a thin overlap has more than minDepth coverage.
  //
  //         -------------               3x coverage
  //          -------------              all overlaps 1 or 2 dashes long
  //                     ---------
  //                      -----------

  if (minCoverage > 0) {
    intervalList<uint32> FI;

    uint32  li = 0;
    uint32  di = 0;

    while ((li < IL.numberOfIntervals()) &&
           (di < ID.numberOfIntervals())) {
      uint32   ll = IL.lo(li);
      uint32   lh = IL.hi(li);
      uint32   dl = ID.lo(di);
      uint32   dh = ID.hi(di);
      uint32   nl  = 0;
      uint32   nh  = 0;

      //  If they intersect, make a new region

      if ((ll <= dl) && (dl < lh)) {
        nl = dl;
        nh = (lh < dh) ? lh : dh;
      }

      if ((dl <= ll) && (ll < dh)) {
        nl = ll;
        nh = (lh < dh) ? lh : dh;
      }

      if (nl < nh)
        FI.add(nl, nh - nl);

      //  Advance the list with the earlier region.

      if (lh <= dh)
        //  IL ends at or before ID
        li++;

      if (dh <= lh) {
        //  ID ends at or before IL
        di++;
      }
    }

    //  Replace the intervals to use with the intersection.

    IL = FI;
  }

  ////////////////////////////////////////

  //  The IL.ct(it) is always 1 if we filter low coverage.  It is no longer reported.
#if 0
  if (IL.numberOfIntervals() > 1)
    for (uint32 it=0; it<IL.numberOfIntervals(); it++)
      fprintf(stderr, "IL[%02d] - iid %d - "F_S64" "F_S64"\n", it, read->gkRead_readID(), IL.lo(it), IL.hi(it));
#endif

  if (IL.numberOfIntervals() == 0) {
    strcpy(logMsg, "\tno high quality overlaps");
    return(false);
  }

  fbgn = IL.lo(0);
  fend = IL.hi(0);

  sprintf(logMsg, "\tskipped %u overlaps; used %u overlaps", nSkip, nUsed);

  for (uint32 it=0; it<IL.numberOfIntervals(); it++) {
    if (IL.hi(it) - IL.lo(it) > fend - fbgn) {
      fbgn = IL.lo(it);
      fend = IL.hi(it);
    }
  }

  return(true);
}


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

static const char *rcsid = "$Id: finalTrim-largestCovered.C,v 1.2 2011/06/23 15:25:23 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "AS_global.H"
#include "AS_PER_gkpStore.H"
#include "AS_OVS_overlapStore.H"
#include "AS_UTL_intervalList.H"


bool
largestCovered(OVSoverlap  *ovl,
               uint32       ovlLen,
               gkFragment  &fr,
               uint32       ibgn,
               uint32       iend,
               uint32      &fbgn,
               uint32      &fend,
               char        *logMsg,
               uint32       errorRate,
               double       errorLimit,
               bool         qvTrimAllowed,
               uint32       minOverlap,
               uint32       minCoverage) {

  fbgn      = ibgn;
  fend      = iend;
  logMsg[0] = 0;

  assert(fr.gkFragment_getReadIID() == ovl[0].a_iid);

  //  Guard against invalid initial clear ranges.  These should all be deleted already.
  if (ibgn + AS_READ_MIN_LEN > iend) {
    strcpy(logMsg, "\tinvalid initial clear range");
    return(false);
  }

  if (ovlLen == 0) {
    strcpy(logMsg, "\tno overlaps");
    return(false);
  }

  intervalList  IL;
  intervalList  ID;
  int32         iid = fr.gkFragment_getReadIID();

  for (uint32 i=0; i<ovlLen; i++) {
    uint32 tbgn = ibgn + ovl[i].dat.obt.a_beg;
    uint32 tend = ibgn + ovl[i].dat.obt.a_end;

    assert(tbgn < tend);
    assert(iid == ovl[i].a_iid);

    if ((ovl->dat.obt.erate > errorRate) &&
        (tend - tbgn) * AS_OVS_decodeQuality(ovl[i].dat.obt.erate) > errorLimit)
      //  Overlap is crappy.
      continue;

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
    intervalDepth  DE(IL);

    uint32  it = 0;
    uint32  ib = 0;
    uint32  ie = 0;

    while (it < DE.numberOfIntervals()) {
      //fprintf(stderr, "DE - %d - "F_S64" "F_S64" "F_S64"\n", fr.gkFragment_getReadIDE(), DE.lo(it), DE.hi(it), DE.de(it));

      if (DE.de(it) < minCoverage) {
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
    intervalList FI;

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

  //for (uint32 it=0; it<IL.numberOfIntervals(); it++)
  //  fprintf(stderr, "IL - %d - "F_S64" "F_S64" "F_S64"\n", fr.gkFragment_getReadIID(), IL.lo(it), IL.hi(it), IL.ct(it));

  if (IL.numberOfIntervals() == 0) {
    strcpy(logMsg, "\tno high quality overlaps");
    return(false);
  }

  fbgn = IL.lo(0);
  fend = IL.hi(0);

  for (uint32 it=0; it<IL.numberOfIntervals(); it++) {
    if (IL.hi(it) - IL.lo(it) > fend - fbgn) {
      fbgn = IL.lo(it);
      fend = IL.hi(it);
    }
  }

  return(true);
}

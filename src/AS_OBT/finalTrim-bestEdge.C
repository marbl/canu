
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2013, J. Craig Venter Institute.
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

static const char *rcsid = "$Id:  $";

#include "finalTrim.H"
#if 0
#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlapStore.h"
#include "AS_UTL_intervalList.H"
#endif

#include <vector>
#include <algorithm>

using namespace std;


bool
bestEdge(OVSoverlap  *ovl,
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

  //
  //  Trim once, using the largest covered rule.  This gets rid of any gaps in overlap coverage,
  //  etc.
  //

  uint32  lbgn = 0;
  uint32  lend = 0;

  if (largestCovered(ovl, ovlLen, fr, ibgn, iend, lbgn, lend, logMsg, errorRate, errorLimit, qvTrimAllowed, minOverlap, minCoverage) == false)
    return(false);

  //
  //  If the read is contained, multiply, return now.
  //

  uint32  nCont = 0;



  //
  //  Trim again, to maximize overlap length.
  //

  int32             iid = fr.gkFragment_getReadIID();
  uint32            len = fr.gkFragment_getSequenceLength();

  vector<uint32>    trim5;
  vector<uint32>    trim3;

  for (uint32 i=0; i<ovlLen; i++) {
    uint32 tbgn = ibgn + ovl[i].dat.obt.a_beg;
    uint32 tend = ibgn + ovl[i].dat.obt.a_end;

    if ((lend <= tbgn) ||
        (tend <= lbgn))
      //  Doesn't intersect the largest covered region.
      continue;

    assert(tbgn < tend);

    assert(ibgn <= tbgn);
    assert(tend <= iend);

    assert(iid == ovl[i].a_iid);

    if ((ovl->dat.obt.erate > errorRate) &&
        (tend - tbgn) * AS_OVS_decodeQuality(ovl[i].dat.obt.erate) > errorLimit)
      //  Overlap is crappy.
      continue;

    trim5.push_back(tbgn);
    trim3.push_back(tend);

    //  If overlap indicates this read is contained, we're done.  There isn't any
    //  trimming needed.

    if ((tbgn == ibgn) &&
        (tend == iend)) {
      trim5.clear();
      trim3.clear();
      break;
    }
  }

  trim5.push_back(ibgn);
  trim3.push_back(iend);

  sort(trim5.begin(), trim5.end(), std::less<uint32>());
  sort(trim3.begin(), trim3.end(), std::greater<uint32>());

  //  Remove duplicate points (easier here than in the loops below)

  {
    uint32 old = 0;
    uint32 sav = 0;

    for (old=0, sav=0; old<trim5.size(); old++)
      if (trim5[old] != trim5[sav])
        trim5[++sav] = trim5[old];
    trim5.resize(sav+1);

    for (old=0, sav=0; old<trim3.size(); old++)
      if (trim3[old] != trim3[sav])
        trim3[++sav] = trim3[old];
    trim3.resize(sav+1);
  }

  uint32   best5pt    = 0;
  uint32   best5score = 0;

  uint32   best3pt    = 0;
  uint32   best3score = 0;

  //  Find the best 5' point.

  for (uint32 pt=0; pt < trim5.size(); pt++) {
    uint32   triml    = trim5[pt];
    uint32   score    = 0;

    if (best5score >= len - triml)
      //  Not possible to get a higher score by trimming more.
      break;

    //fprintf(stderr, "trim5 pt %u out of %u\n", pt, trim5.size());

    for (uint32 i=0; i < ovlLen; i++) {
      uint32 tbgn = ibgn + ovl[i].dat.obt.a_beg;
      uint32 tend = ibgn + ovl[i].dat.obt.a_end;

      if ((triml <  tbgn) ||
          (tend  <= triml))
        //  Alignment starts after of the trim point; not a valid overlap.
        //  or, alignment ends before the trim point (trimmed out).
        continue;

      //fprintf(stderr, "trim5 iid %u %u pos %u-%u with trim %u\n",
      //        ovl[i].a_iid, ovl[i].b_iid, tbgn, tend, triml);

      uint32 tlen = tend - triml;

      assert(tlen <= len);

      if (score < tlen)
        score = tlen;
    }

    if (best5score < score) {
      best5pt    = pt;
      best5score = score;
    }
  }

  //fprintf(stderr, "BEST at %u position %u pt %u\n", best5score, trim5[best5pt], best5pt);

  //  Find the best 3' point.

  for (uint32 pt=0; pt<trim3.size(); pt++) {
    uint32   trimr    = trim3[pt];;
    uint32   score    = 0;

    if (best3score >= trimr - 0)
      //  Not possible to get a higher score by trimming more.
      break;

    //fprintf(stderr, "trim3 pt %u out of %u\n", pt, trim3.size());

    for (uint32 i=0; i < ovlLen; i++) {
      uint32 tbgn = ibgn + ovl[i].dat.obt.a_beg;
      uint32 tend = ibgn + ovl[i].dat.obt.a_end;

      if ((tend < trimr) ||
          (trimr <= tbgn))
        //  Alignment ends before the trim point; not a valid overlap,
        //  or, alignment starts after the trim point (trimmed out)
        continue;

      //fprintf(stderr, "trim3 iid %u %u pos %u-%u with trim %u\n",
      //        ovl[i].a_iid, ovl[i].b_iid, tbgn, tend, trimr);

      uint32 tlen = trimr - tbgn;

      assert(tlen <= len);

      if (score < tlen)
        score = tlen;
    }

    if (best3score < score) {
      best3pt    = pt;
      best3score = score;
    }
  }

  //fprintf(stderr, "BEST at %u position %u pt %u\n", best3score, trim3[best3pt], best3pt);

  //
  //  Set trimming.  Be just a little aggressive, and get rid of an extra base or two.
  //

  fbgn = trim5[best5pt] + 2;
  fend = trim3[best3pt] - 2;

  return(true);
}

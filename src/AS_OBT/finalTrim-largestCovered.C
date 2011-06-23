
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

static const char *rcsid = "$Id: finalTrim-largestCovered.C,v 1.2 2011-06-23 15:25:23 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlapStore.h"
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
               bool         qvTrimAllowed) {

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

  IL.merge();

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

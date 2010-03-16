
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

const char *mainid = "$Id: overlapMask.C,v 1.2 2010-03-16 05:27:57 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "util++.H"
#include "readOverlap.H"

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlapStore.h"

//  Implements a very simple trimming based on overlaps.  Trim to the largest region covered by
//  'good' overlap.  A 'good' overlap is the same as in BOG: less than 1.5% error or less than 3
//  errors.

#define MAX_OVERLAPS_PER_FRAG   (16 * 1024 * 1024)




int
main(int argc, char **argv) {
  uint32             errorRate    = AS_OVS_encodeQuality(0.015);
  double             errorLimit   = 2.5;
  gkStore           *gkp          = 0L;
  OverlapStore      *ovsprimary   = 0L;
  OverlapStore      *ovssecondary = 0L;

  bool               testing      = false;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-gkp", 2) == 0) {
      gkp = new gkStore(argv[++arg], FALSE, testing == false);
      gkp->gkStore_enableClearRange(AS_READ_CLEAR_OBTMERGE);

      //  The cache is not enabled, as we don't expect many changes to the store.
      gkp->gkStore_metadataCaching(false);

    } else if (strncmp(argv[arg], "-ovs", 2) == 0) {
      if (ovsprimary == NULL)
        ovsprimary = AS_OVS_openOverlapStore(argv[++arg]);
      else if (ovssecondary == NULL)
        ovssecondary = AS_OVS_openOverlapStore(argv[++arg]);
      else {
        fprintf(stderr, "Only two obtStores allowed.\n");
        err++;
      }

    } else if (strncmp(argv[arg], "-e", 2) == 0) {
      double erate = atof(argv[++arg]);
      if (erate >= AS_MAX_ERROR_RATE)
        fprintf(stderr, "Error rate %s too large; must be 'fraction error' and below %f\n", argv[arg], AS_MAX_ERROR_RATE), exit(1);
      errorRate = AS_OVS_encodeQuality(erate);

    } else if (strncmp(argv[arg], "-E", 2) == 0) {
      errorLimit = atof(argv[++arg]);

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((gkp == 0L) || (ovsprimary == 0L) || (err)) {
    fprintf(stderr, "usage: %s [-1] -gkp <gkpStore> -ovs <ovsStore> [opts]\n", argv[0]);
    fprintf(stderr, "  -e E        filter overlaps above this fraction error; default 0.015 (== 1.5%% error)\n");
    fprintf(stderr, "  -E E        filter overlaps above this number errors; default 2\n");
    exit(1);
  }


  bool            nothingToDo = true;

  for (uint32 i=1; i<=gkp->gkStore_getNumLibraries(); i++) {
    gkLibrary  *gkl = gkp->gkStore_getLibrary(i);

    //if (gkl->doOverlapMasking == true) {
    if (1) {
      fprintf(stderr, "Applying overlap masking to library %s.\n", AS_UID_toString(gkl->libraryUID));
      nothingToDo = false;
    } else {
      fprintf(stderr, "Ignoring library %s.\n", AS_UID_toString(gkl->libraryUID));
    }
  }



  uint32      ovlMax = MAX_OVERLAPS_PER_FRAG;
  uint32      ovlLen = 0;
  OVSoverlap *ovl    = (OVSoverlap *)safe_malloc(sizeof(OVSoverlap) * ovlMax);

  ovlLen += AS_OVS_readOverlapsFromStore(ovsprimary,   ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
  ovlLen += AS_OVS_readOverlapsFromStore(ovssecondary, ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);

  while (ovlLen > 0) {
    intervalList  IL;
    gkFragment    fr;

    int32         iid = ovl[0].a_iid;

    gkp->gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);

    int32 ibgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_OBTINITIAL);
    int32 iend = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_OBTINITIAL);

    if (ibgn >= iend)
      //  Reset invalid and/or unset clear ranges to null.  This shouldn't happen in a
      //  normal assembly (some flavor of initial trimming should have been run already).
      ibgn = iend = 0;

    for (uint32 i=0; i<ovlLen; i++) {
      assert(iid == ovl[i].a_iid);

      int32 tbgn = ibgn + ovl[i].dat.obt.a_beg;
      int32 tend = ibgn + ovl[i].dat.obt.a_end;

      assert(tbgn < tend);

      if (ovl->dat.obt.erate > errorRate)
        //  Overlap is crappy.
        continue;

      if ((tend - tbgn) * AS_OVS_decodeQuality(ovl[i].dat.obt.erate) > errorLimit)
        //  Overlap is crappy.
        continue;

      IL.add(tbgn, tend - tbgn);
    }

    IL.merge();

    if (IL.numberOfIntervals() == 0) {
      fprintf(stderr, "%d\tINIT\t%d\t%d\tTRIM\t%d\t%d (deleted)\n",
              iid, ibgn, iend, 0, 0);
      gkp->gkStore_delFragment(iid, false);
    } else {
      int32 tbgn = IL.lo(0);
      int32 tend = IL.hi(0);

      for (uint32 interval=1; interval<IL.numberOfIntervals(); interval++) {
        if (IL.hi(interval) - IL.lo(interval) > tend - tbgn) {
          tbgn = IL.lo(interval);
          tend = IL.hi(interval);
        }
      }

      IL.clear();

      assert(iid != -1);

      if ((ibgn != tbgn) || (iend != tend)) {
        fprintf(stderr, "%d\tINIT\t%d\t%d\tTRIM\t%d\t%d\n",
                iid, ibgn, iend, tbgn, tend);
        fr.gkFragment_setClearRegion(tbgn, tend, AS_READ_CLEAR_OBTMERGE);
        gkp->gkStore_setFragment(&fr);
      }
    }

    ovlLen  = 0;
    ovlLen += AS_OVS_readOverlapsFromStore(ovsprimary,   ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
    ovlLen += AS_OVS_readOverlapsFromStore(ovssecondary, ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
  }

  delete gkp;

  exit(0);
}

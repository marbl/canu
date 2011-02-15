
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id: classifyMatesSingle.C,v 1.1 2011-02-15 08:10:11 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"


void
fixOverlap(OVSoverlap &ovl) {

  //if (ovl.b_iid != aIID) {
    OVSoverlap  flipped = ovl;

    flipped.a_iid = ovl.b_iid;
    flipped.b_iid = ovl.a_iid;

    if (ovl.dat.ovl.flipped) {
      //  Different orientation.  Swap the hangs.
      flipped.dat.ovl.a_hang =  ovl.dat.ovl.b_hang;
      flipped.dat.ovl.b_hang =  ovl.dat.ovl.a_hang;

    } else {
      //  Same orientation.  Negate the hangs.
      flipped.dat.ovl.a_hang = -ovl.dat.ovl.a_hang;
      flipped.dat.ovl.b_hang = -ovl.dat.ovl.b_hang;
    }

    ovl = flipped;
    //}

#if 0
  if (AS_OVS_overlapBIsContained(ovl) == false) {
    //  B not contained in A
    ovl.a_iid = 0;
    ovl.b_iid = 0;
  }
#endif
}



int
main(int argc, char **argv) {
  char      *gkpStorePath = NULL;
  char      *ovlStorePath = NULL;

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlStorePath = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }
  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore\n", argv[0]);
    exit(1);
  }


  gkStore          *gkpStore  = new gkStore(gkpStorePath, FALSE, FALSE);
  gkStream         *gkpStream = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);

  OverlapStore     *ovlStore  = AS_OVS_openOverlapStore(ovlStorePath);

  gkFragment   aFrg;
  uint32       aOvlMax = 1048576;
  uint32       aOvlLen = 0;
  OVSoverlap  *aOvl    = new OVSoverlap [aOvlMax];

  gkFragment   bFrg;
  uint32       bOvlMax = 1048576;
  uint32       bOvlLen = 0;
  OVSoverlap  *bOvl    = new OVSoverlap [bOvlMax];

  uint32       verified = 0;
  uint32       tested   = 0;

  while (gkpStream->next(&aFrg)) {
    uint32 aIID = aFrg.gkFragment_getReadIID();
    uint32 bIID = aFrg.gkFragment_getMateIID();
    uint32 lib  = aFrg.gkFragment_getLibraryIID();

    if (gkpStore->gkStore_getLibrary(lib)->doMerBasedTrimming == 0)
      //  Not an illumina fragment
      //  TODO:  Make a real library feature for this
      continue;

    if (aIID + 1 != bIID)
      //  Illumina mated fragments are always adjacent.
      continue;

    //  And we can now get the bFrg - it's guaranteed to be the next in the stream.

    gkpStream->next(&bFrg);

    if (aFrg.gkFragment_getIsDeleted() || bFrg.gkFragment_getIsDeleted())
      continue;

    assert(aFrg.gkFragment_getMateIID() == bFrg.gkFragment_getReadIID());
    assert(bFrg.gkFragment_getMateIID() == aFrg.gkFragment_getReadIID());

    //  Load overlaps.

    AS_OVS_setRangeOverlapStore(ovlStore, aIID, bIID);

    aOvlLen = AS_OVS_readOverlapsFromStore(ovlStore, aOvl, aOvlMax, AS_OVS_TYPE_OVL);
    bOvlLen = AS_OVS_readOverlapsFromStore(ovlStore, bOvl, bOvlMax, AS_OVS_TYPE_OVL);

    if ((aOvlLen == 0) || (bOvlLen == 0))
      continue;

    //fprintf(stderr, "Read %u overlaps for %u -- %u overlaps for mate %u\n",
    //        aOvlLen, aIID,
    //        bOvlLen, bIID);

    //  Normalize overlaps.
    //  o Discard any that we are not contained in
    //  o Swap so that we are the b_iid fragment (makes things simpler later)

    for (uint32 i=0; i<aOvlLen; i++) {
      assert(aOvl[i].b_iid != aIID);
      fixOverlap(aOvl[i]);
    }

    for (uint32 i=0; i<bOvlLen; i++) {
      assert(bOvl[i].b_iid != bIID);
      fixOverlap(bOvl[i]);
    }


    //  Do they share a common fragment?

    bool    common    = false;
    uint32  commona   = 0;
    uint32  commonb   = 0;
    int32   commonLen = 0;

    uint32  aLen = aFrg.gkFragment_getClearRegionLength(AS_READ_CLEAR_OBTCHIMERA);
    uint32  bLen = bFrg.gkFragment_getClearRegionLength(AS_READ_CLEAR_OBTCHIMERA);

    for (uint32 a=0; a<aOvlLen; a++) {
      uint32  ao = aOvl[a].a_iid;
      if (ao == 0)
        continue;

      for (uint32 b=0; b<bOvlLen; b++) {
        uint32  bo = bOvl[b].a_iid;
        if (bo == 0)
          continue;

        if (ao != bo)
          //  Not the same overlapping fragment.
          continue;

        //  Compute placement on that fragment.  Since the larger fragment is
        //  A, and A is always forward, we can compute the start position of
        //  each of the two fragments simply -- it's the a_hang!

        int32  abgn  = aOvl[a].dat.ovl.a_hang;
        int32  aend  = abgn + aLen;
        int32  aflip = aOvl[a].dat.ovl.flipped;

        int32  bbgn  = bOvl[b].dat.ovl.a_hang;
        int32  bend  = aend + bLen;
        int32  bflip = bOvl[b].dat.ovl.flipped;

        //fprintf(stderr, "a: %d %d -- b %d %d\n",
        //       abgn, aflip, bbgn, bflip);

#if 0
        //  Check orientation.  These should be innie.
        //
        //       a----->     b<------
        //       b----->     a<------
        //
        if (!(((abgn < bbgn) && (aflip == false) && (bflip == true)) ||
              ((bbgn < abgn) && (bflip == false) && (aflip == true))))
          //  Nope, neither case was true.
          continue;
#endif

        //  Check orientation.  These should be OUTTIE.  fastqToCA & gatekeeper will flip
        //  the orientation of the reads, leaving the MPs innie and the PEs outtie.
        //
        //       a<-----     b------>
        //       b<-----     a------>
        //
        if (!(((abgn < bbgn) && (aflip == true) && (bflip == false)) ||
              ((bbgn < abgn) && (bflip == true) && (aflip == false))))
          //  Nope, neither case was true.
          continue;

        common    = true;
        commona   = a;
        commonb   = b;
        commonLen = MAX(aend, bend) - MIN(abgn, bbgn);
      }
    }

    tested++;

    if (common) {
      verified++;
      fprintf(stdout, "%u and %u VERIFIED with fragment %u insert size %d\n", aIID, bIID, aOvl[commona].a_iid, commonLen);
    } else {
      //fprintf(stdout, "%u and %u NOT verified\n", aIID, bIID);
    }

  }

  delete gkpStore;
  AS_OVS_closeOverlapStore(ovlStore);

  fprintf(stderr, "Found %d innie out of %d tested.\n", verified, tested);

  exit(0);
}


/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id: AS_BOG_MateBubble.cc,v 1.6 2010-10-27 04:31:20 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"


void
UnitigGraph::popMateBubbles(OverlapStore *ovlStoreUniq,
                            OverlapStore *ovlStoreRept) {

#define MAX_OVERLAPS_PER_FRAG   (16 * 1024 * 1024)

  uint32      ovlMax = MAX_OVERLAPS_PER_FRAG;
  uint32      ovlLen = 0;
  OVSoverlap *ovl    = (OVSoverlap *)safe_malloc(sizeof(OVSoverlap) * ovlMax);
  uint32     *ovlCnt = (uint32     *)safe_malloc(sizeof(uint32)     * AS_READ_MAX_NORMAL_LEN);

  uint32      nBubblePopped   = 0;
  uint32      nBubbleTooBig   = 0;
  uint32      nBubbleConflict = 0;

  fprintf(logFile, "==> SEARCHING FOR MATE BUBBLES\n");

  //  For each unitig, if all (or most) of the external mates are to a single other unitig (not
  //  counting singletons), then this is a potential bubble popping unitig.
  //
  //  At present, this is exploratory only.

  for (int  ti=0; ti<unitigs.size(); ti++) {
    Unitig        *tig = unitigs[ti];

    if ((tig == NULL) ||
        (tig->ufpath.size() == 0))
      //   No tig here.
      continue;

    if ((tig->getLength() > 1000) ||
        (tig->ufpath.size() >= 3000))
      //  Tig too big.
      continue;

    //if ((tig->getLength() < 150) ||
    //    (tig->ufpath.size() < 5))
    //  //  Tig too small.
    //  continue;

    uint32        *lkg    = new uint32 [tig->ufpath.size()];
    uint32         lkgLen = 0;
    uint32         lkgExt = 0;

    for (int fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode *frg = &tig->ufpath[fi];
      int32         frgID = frg->ident;
      int32         matID = FI->mateIID(frgID);

      int32         mtigID = 0;
      Unitig       *mtig   = 0L;

      if (matID == 0)
        //  No mate.
        continue;

      mtigID = tig->fragIn(matID);
      mtig   = unitigs[mtigID];

      if (mtigID == tig->id())
        //  Mate is not external.
        continue;

      lkgExt++;

      if (mtig->ufpath.size() < 2)
        //  Mate is in singleton.
        continue;

      lkg[lkgLen++] = mtigID;
    }

    if (lkgLen == 0)
      //  No external mates.
      continue;

    sort(lkg, lkg+lkgLen);

    int32  last = lkg[0];
    int32  lcnt = 1;

    for (uint32 i=1; i<lkgLen; i++) {
      if (last != lkg[i]) {
        if ((lcnt > 3))
          fprintf(logFile, "popMateBubble()-- tig %d len %d might pop bubble in tig %d (%d mates in there out of %d external mates)\n",
                  tig->id(), tig->getLength(), last, lcnt, lkgExt);
        last = lkg[i];
        lcnt = 0;
      }

      lcnt++;
    }

    if ((lcnt > 3))
      fprintf(logFile, "popMateBubble()-- tig %d len %d might pop bubble in tig %d (%d mates in there out of %d external mates)\n",
              tig->id(), tig->getLength(), last, lcnt, lkgExt);

    delete [] lkg;
  }
}    


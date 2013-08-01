
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

static const char *rcsid = "$Id$";

#include "AS_BOG_BestOverlapGraph.H"
#include "AS_BOG_UnitigGraph.H"
#include "AS_BOG_MateLocation.H"


void
UnitigGraph::evaluateMates(void) {

  //  [0] -- BOTH frag and mate are dovetail
  //  [1] -- ONE frag dovetail, ONE frag contained
  //  [2] -- BOTH frag and mate are contained

  uint64   unmated[3]         = { 0, 0, 0 };
  uint64   mated[3]           = { 0, 0, 0 };
  uint64   different[3][3]    = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
  uint64   happy[3]           = { 0, 0, 0 };
  uint64   grumpy[3]          = { 0, 0, 0 };

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *thisUtg = unitigs[ti];

    if ((thisUtg == NULL) ||
        (thisUtg->ufpath.size() < 2))
      continue;

    MateLocation          positions(thisUtg);

    for (uint32 fi=0; fi<thisUtg->ufpath.size(); fi++) {
      ufNode  *thisFrg = &thisUtg->ufpath[fi];

      uint32  thisFrgID = thisFrg->ident;
      uint32  mateFrgID = FI->mateIID(thisFrg->ident);

      BestContainment *thiscont = OG->getBestContainer(thisFrgID);
      BestContainment *matecont = OG->getBestContainer(mateFrgID);

      uint32  type = (thiscont != NULL) + (matecont != NULL);

      //  Trivial case, not a mated fragment.

      if (mateFrgID == 0) {
        unmated[type]++;
        continue;
      }

      uint32  thisUtgID = thisUtg->fragIn(thisFrgID);
      uint32  mateUtgID = thisUtg->fragIn(mateFrgID);

      MateLocationEntry  mloc     = positions.getById(thisFrg->ident);

      //  Skip this fragment, unless it is mleFrgID1.  Fragments with mates in other unitigs
      //  are always listed in ID1, and mates completely in this unitig should be counted once.
      if (mloc.mleFrgID1 != thisFrgID)
        continue;

      mated[type]++;

      //  Easy case, both fragments in the same unitig.  We can use the isGrumpy flag directly to
      //  decide if the mates are happy or not.

      if (thisUtgID == mateUtgID) {
        assert(mloc.mleUtgID1 == thisUtg->id());
        assert(mloc.mleUtgID2 == thisUtg->id());
        assert(mloc.mleFrgID1 == thisFrgID);
        assert(mloc.mleFrgID2 == mateFrgID);

        if (mloc.isGrumpy == false)
          happy[type]++;
        else
          grumpy[type]++;

        continue;
      }

      //  Hard case, fragments in different unitigs.  We want to distinguish between
      //  three cases:
      //    1) mates at the end that could potentially join unitigs across a gap
      //    2) mate at the end to an interior mate -- possibly a repeat
      //    3) both interior mates

      assert(mloc.mleUtgID1 == thisUtg->id());
      assert(mloc.mleUtgID2 == 0);
      assert(mloc.mleFrgID1 == thisFrgID);
      assert(mloc.mleFrgID2 == 0);

      //  Get the mate frag.

      Unitig        *mateUtg = unitigs[mateUtgID];
      ufNode        *mateFrg = &mateUtg->ufpath[mateUtg->pathPosition(mateFrgID)];

      //differentSum[type]++;

      bool  fragIsInterior = false;
      bool  mateIsInterior = false;

      uint32  lib           = FI->libraryIID(thisFrg->ident);
      int32   minInsertSize = IS->mean(lib) - BADMATE_INTRA_STDDEV * IS->stddev(lib);
      int32   maxInsertSize = IS->mean(lib) + BADMATE_INTRA_STDDEV * IS->stddev(lib);

      if (thisFrg->position.bgn < thisFrg->position.end) {
        //  Fragment is forward, so mate should be after it.
        if (thisUtg->getLength() - thisFrg->position.bgn > maxInsertSize)
          fragIsInterior = true;
      } else {
        //  Fragment is reverse, so mate should be before it.
        if (thisFrg->position.bgn > maxInsertSize)
          fragIsInterior = true;
      }

      if (mateFrg->position.bgn < mateFrg->position.end) {
        //  Fragment is forward, so mate should be after it.
        if (mateUtg->getLength() - mateFrg->position.bgn > maxInsertSize)
          mateIsInterior = true;
      } else {
        //  Fragment is reverse, so mate should be before it.
        if (mateFrg->position.bgn > maxInsertSize)
          mateIsInterior = true;
      }

      uint32  dtyp = (fragIsInterior == true) + (mateIsInterior == true);

      different[type][dtyp]++;
    }
  }

#warning THIS IS COMPLETELY BROKEN
  fprintf(logFile, "MATE HAPPINESS (dove/dove):  unmated %11"F_U64P"  mated %11"F_U64P"  sameTig: happy %11"F_U64P" grumpy %11"F_U64P"  diffTig: end-end %11"F_U64P" end-int %11"F_U64P" int-int %11"F_U64P"\n",
          unmated[0], mated[0], happy[0], grumpy[0], different[0][0], different[0][1], different[0][2]);
  fprintf(logFile, "MATE HAPPINESS (dove/cont):  unmated %11"F_U64P"  mated %11"F_U64P"  sameTig: happy %11"F_U64P" grumpy %11"F_U64P"  diffTig: end-end %11"F_U64P" end-int %11"F_U64P" int-int %11"F_U64P"\n",
          unmated[1], mated[1], happy[1], grumpy[1], different[1][0], different[1][1], different[1][2]);
  fprintf(logFile, "MATE HAPPINESS (cont/cont):  unmated %11"F_U64P"  mated %11"F_U64P"  sameTig: happy %11"F_U64P" grumpy %11"F_U64P"  diffTig: end-end %11"F_U64P" end-int %11"F_U64P" int-int %11"F_U64P"\n",
          unmated[2], mated[2], happy[2], grumpy[2], different[2][0], different[2][1], different[2][2]);
}


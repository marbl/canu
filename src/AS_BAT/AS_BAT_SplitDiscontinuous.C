
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

static const char *rcsid = "$Id: AS_BAT_SplitDiscontinuous.C,v 1.6 2012-03-13 21:43:43 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_MateLocation.H"



static
void
makeNewUnitig(UnitigVector &unitigs,
              uint32        splitFragsLen,
              ufNode       *splitFrags) {
  Unitig *dangler = new Unitig(false);

  if (logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS))
    fprintf(logFile, "splitDiscontinuous()--   new tig "F_U32" with "F_U32" fragments (starting at frag "F_U32").\n",
            dangler->id(), splitFragsLen, splitFrags[0].ident);

  int splitOffset = -MIN(splitFrags[0].position.bgn, splitFrags[0].position.end);

  //  This should already be true, but we force it still
  splitFrags[0].contained = 0;

  for (uint32 i=0; i<splitFragsLen; i++)
    dangler->addFrag(splitFrags[i], splitOffset, false);  //logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS));

  unitigs.push_back(dangler);
}




//  After splitting and ejecting some contains, check for discontinuous unitigs.
//
void splitDiscontinuousUnitigs(UnitigVector &unitigs) {

  fprintf(logFile, "==> SPLIT DISCONTINUOUS\n");

  uint32                numTested  = 0;
  uint32                numSplit   = 0;
  uint32                numCreated = 0;

  uint32                splitFragsLen = 0;
  uint32                splitFragsMax = 0;
  ufNode               *splitFrags    = NULL;

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) || (tig->ufpath.size() < 2))
      continue;

    //  Unitig must be sorted.  Someone upstream os screwing this up.
    tig->sort();

    //  We'll want to build an array of new fragments to split out.  This can be up
    //  to the size of the largest unitig.
    splitFragsMax = MAX(splitFragsMax, tig->ufpath.size());

    //  Check that the unitig starts at position zero.  Not critical for the next loop, but
    //  needs to be dome sometime.
    int32   minPos = MIN(tig->ufpath[0].position.bgn, tig->ufpath[0].position.end);

    if (minPos == 0)
      continue;

    fprintf(logFile, "splitDiscontinuous()-- tig "F_U32" offset messed up; reset by "F_S32".\n", tig->id(), minPos);

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode  *frg = &tig->ufpath[fi];

      frg->position.bgn -= minPos;
      frg->position.end -= minPos;
    }
  }

  splitFrags = new ufNode [splitFragsMax];

  //  Now, finally, we can check for gaps in unitigs.

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) || (tig->ufpath.size() < 2))
      continue;

    //  We don't expect many unitigs to be broken, so we'll do a first quick pass to just
    //  test if it is.

    int32  maxEnd   = MAX(tig->ufpath[0].position.bgn, tig->ufpath[0].position.end);
    bool   isBroken = false;

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode  *frg = &tig->ufpath[fi];

      int32    bgn = MIN(frg->position.bgn, frg->position.end);
      int32    end = MAX(frg->position.bgn, frg->position.end);

      if (bgn > maxEnd - AS_OVERLAP_MIN_LEN) {
        isBroken = true;
        break;
      }

      maxEnd = MAX(maxEnd, end);
    }

    numTested++;

    if (isBroken == false)
      continue;

    numSplit++;

    //  Dang, busted unitig.  Fix it up.

    splitFragsLen = 0;
    maxEnd        = 0;

    if (logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS))
      fprintf(logFile, "splitDiscontinuous()-- discontinuous tig "F_U32" with "F_SIZE_T" fragments broken into:\n",
              tig->id(), tig->ufpath.size());

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode  *frg = &tig->ufpath[fi];

      int32    bgn = MIN(frg->position.bgn, frg->position.end);
      int32    end = MAX(frg->position.bgn, frg->position.end);

      //  Good thick overlap exists to this fragment, save it.
      if (bgn <= maxEnd - AS_OVERLAP_MIN_LEN) {
        assert(splitFragsLen < splitFragsMax);
        splitFrags[splitFragsLen++] = *frg;
        maxEnd = MAX(maxEnd, end);
        continue;
      }

      //  No thick overlap found.  We need to break right here before the current fragment.

      //  If there is exactly one fragment, and it's contained, and it's not mated, move it to the
      //  container.  (This has a small positive benefit over just making every read a singleton).
      //
      if ((splitFragsLen == 1) &&
          (FI->mateIID(splitFrags[0].ident) == 0) &&
          (splitFrags[0].contained != 0)) {
        Unitig  *dangler  = unitigs[tig->fragIn(splitFrags[0].contained)];

        //  If the parent isn't in a unitig, we must have shattered the repeat unitig it was in.
        //  Do the same here.

        if (dangler == NULL) {
          if (logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS))
            fprintf(logFile, "splitDiscontinuous()--   singleton frag "F_U32" shattered.\n",
                    splitFrags[0].ident);
          Unitig::removeFrag(splitFrags[0].ident);

        } else {
          assert(dangler->id() == tig->fragIn(splitFrags[0].contained));

          if (logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS))
            fprintf(logFile, "splitDiscontinuous()--   old tig "F_U32" with "F_SIZE_T" fragments (contained frag "F_U32" moved here).\n",
                    dangler->id(), dangler->ufpath.size() + 1, splitFrags[0].ident);

          BestContainment  *bestcont = OG->getBestContainer(splitFrags[0].ident);

          assert(bestcont->isContained == true);

          dangler->addContainedFrag(splitFrags[0].ident, bestcont, false);
          dangler->bubbleSortLastFrag();

          assert(dangler->id() == Unitig::fragIn(splitFrags[0].ident));
        }
      }

      //  Otherwise, make an entirely new unitig for these fragments.
      else {
        numCreated++;
        makeNewUnitig(unitigs, splitFragsLen, splitFrags);
        tig = unitigs[ti];
      }

      //  Done with the split, save the current fragment.  This resets everything.

      splitFragsLen = 0;
      splitFrags[splitFragsLen++] = *frg;

      maxEnd = end;
    }


    //  If we did any splitting, then the length of the frags in splitFrags will be less than the length
    //  of the path in the current unitig.  Make a final new unitig for the remaining fragments.
    //
    if (splitFragsLen != tig->ufpath.size()) {
      numCreated++;
      makeNewUnitig(unitigs, splitFragsLen, splitFrags);

      delete unitigs[ti];
      unitigs[ti] = NULL;
    }
  }

  fprintf(logFile, "splitDiscontinuous()-- Tested "F_U32" unitigs, split "F_U32" into "F_U32" new unitigs.\n",
          numTested, numSplit, numCreated);

  delete [] splitFrags;
}


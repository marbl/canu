
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
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_SplitDiscontinuous.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2015-MAR-03
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"


static
void
makeNewUnitig(UnitigVector &unitigs,
              uint32        splitFragsLen,
              ufNode       *splitFrags) {
  Unitig *dangler = unitigs.newUnitig(false);

  if (logFileFlagSet(LOG_SPLIT_DISCONTINUOUS))
    writeLog("splitDiscontinuous()--   new tig "F_U32" with "F_U32" fragments (starting at frag "F_U32").\n",
            dangler->id(), splitFragsLen, splitFrags[0].ident);

  int splitOffset = -MIN(splitFrags[0].position.bgn, splitFrags[0].position.end);

  //  This should already be true, but we force it still
  splitFrags[0].contained = 0;

  for (uint32 i=0; i<splitFragsLen; i++)
    dangler->addFrag(splitFrags[i], splitOffset, false);  //logFileFlagSet(LOG_SPLIT_DISCONTINUOUS));
}




//  After splitting and ejecting some contains, check for discontinuous unitigs.
//
void splitDiscontinuousUnitigs(UnitigVector &unitigs, uint32 minOverlap) {

  writeLog("==> SPLIT DISCONTINUOUS\n");

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

    writeLog("splitDiscontinuous()-- tig "F_U32" offset messed up; reset by "F_S32".\n", tig->id(), minPos);

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

      if (bgn > maxEnd - minOverlap) {
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

    if (logFileFlagSet(LOG_SPLIT_DISCONTINUOUS))
      writeLog("splitDiscontinuous()-- discontinuous tig "F_U32" with "F_SIZE_T" fragments broken into:\n",
              tig->id(), tig->ufpath.size());

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode  *frg = &tig->ufpath[fi];

      int32    bgn = MIN(frg->position.bgn, frg->position.end);
      int32    end = MAX(frg->position.bgn, frg->position.end);

      //  Good thick overlap exists to this fragment, save it.
      if (bgn <= maxEnd - minOverlap) {
        assert(splitFragsLen < splitFragsMax);
        splitFrags[splitFragsLen++] = *frg;
        maxEnd = MAX(maxEnd, end);
        continue;
      }

      //  No thick overlap found.  We need to break right here before the current fragment.  We used
      //  to try to place contained reads with their container.  For simplicity, we instead just
      //  make a new unitig, letting the main() decide what to do with them (e.g., bubble pop or try
      //  to place all reads in singleton unitigs as contained reads again).

      numCreated++;
      makeNewUnitig(unitigs, splitFragsLen, splitFrags);
      tig = unitigs[ti];

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

  writeLog("splitDiscontinuous()-- Tested "F_U32" unitigs, split "F_U32" into "F_U32" new unitigs.\n",
          numTested, numSplit, numCreated);

  delete [] splitFrags;
}



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
 *    src/AS_BAT/AS_BAT_MateBubble.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz beginning on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static const char *rcsid = "$Id$";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"


void
popMateBubbles(UnitigVector &unitigs) {
  uint32      nBubblePopped   = 0;
  uint32      nBubbleTooBig   = 0;
  uint32      nBubbleConflict = 0;

  writeLog("==> SEARCHING FOR MATE BUBBLES\n");

  //  For each unitig, if all (or most) of the external mates are to a single other unitig (not
  //  counting singletons), then this is a potential bubble popping unitig.
  //
  //  At present, this is exploratory only.

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
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

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode *frg = &tig->ufpath[fi];
      int32         frgID = frg->ident;
      int32         matID = FI->mateIID(frgID);

      uint32        mtigID = 0;
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

    uint32  last = lkg[0];
    uint32  lcnt = 1;

    for (uint32 i=1; i<lkgLen; i++) {
      if (last != lkg[i]) {
        if ((lcnt > 3))
          writeLog("popMateBubble()-- tig %d len %d might pop bubble in tig %u (%u mates in there out of %d external mates)\n",
                  tig->id(), tig->getLength(), last, lcnt, lkgExt);
        last = lkg[i];
        lcnt = 0;
      }

      lcnt++;
    }

    if ((lcnt > 3))
      writeLog("popMateBubble()-- tig %d len %d might pop bubble in tig %u (%u mates in there out of %d external mates)\n",
              tig->id(), tig->getLength(), last, lcnt, lkgExt);

    delete [] lkg;
  }
}


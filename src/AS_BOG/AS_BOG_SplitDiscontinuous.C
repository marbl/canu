
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

static const char *rcsid = "$Id: AS_BOG_SplitDiscontinuous.cc,v 1.1 2010-10-01 13:43:49 brianwalenz Exp $";

#include "AS_BOG_BestOverlapGraph.H"
#include "AS_BOG_UnitigGraph.H"
#include "AS_BOG_MateLocation.H"


//  After splitting and ejecting some contains, check for discontinuous unitigs.
//
void UnitigGraph::splitDiscontinuousUnitigs(void) {

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) ||
        (tig->ufpath.size() < 2))
      continue;

    //  Check for discontinuities

    int32                 maxEnd   = 0;

    ufNode               *splitFrags    = new ufNode [tig->ufpath.size()];
    uint32                splitFragsLen = 0;

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode  *frg = &tig->ufpath[fi];

      //  If this is the first frag in this block (we are at
      //  the start of a unitig, or just split off a new
      //  unitig), remember the end location.
      //
      if (splitFragsLen == 0) {
        maxEnd =  MAX(frg->position.bgn, frg->position.end);
      }

      //  We require at least (currently 40bp, was 10bp hardcoded
      //  here) of overlap between fragments.  If we don't have that,
      //  split off the fragments we've seen.
      //
      //  10bp was a bad choice.  It caught most of the breaks, but
      //  missed one class; when a container fragment is moved out of
      //  the unitig, fragments contained in there are marked as
      //  uncontained.  That container fragment could have been the
      //  one holding the unitig together:
      //
      //  -----------------   <- container (removed)
      //    --------
      //      ---------
      //              -----------------
      //
      //  Because the two small guys are marked as uncontained, they
      //  are assumed to have a good dovetail overlap.
      //
      if (maxEnd - AS_OVERLAP_MIN_LEN < MIN(frg->position.bgn, frg->position.end)) {

        //  If there is exactly one fragment, and it's contained, and
        //  it's not mated, move it to the container.  (This has a
        //  small positive benefit over just making every read a
        //  singleton).
        //
        if ((splitFragsLen == 1) &&
            (FI->mateIID(splitFrags[0].ident) == 0) &&
            (splitFrags[0].contained != 0)) {

          Unitig           *dangler  = unitigs[tig->fragIn(splitFrags[0].contained)];
          BestContainment  *bestcont = OG->getBestContainer(splitFrags[0].ident);

          assert(dangler->id() == tig->fragIn(splitFrags[0].contained));

          if (logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS))
            fprintf(logFile, "Dangling contained fragment %d in unitig %d -> move them to container unitig %d\n",
                    splitFrags[0].ident, tig->id(), dangler->id());

          dangler->addContainedFrag(splitFrags[0].ident, bestcont, logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS));
          assert(dangler->id() == Unitig::fragIn(splitFrags[0].ident));

        } else {
          Unitig *dangler = new Unitig(logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS));

          if (logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS))
            fprintf(logFile, "Dangling fragments in unitig %d -> move them to unitig %d\n", tig->id(), dangler->id());

          int splitOffset = -MIN(splitFrags[0].position.bgn, splitFrags[0].position.end);

          //  This should already be true, but we force it still
          splitFrags[0].contained = 0;

          for (uint32 i=0; i<splitFragsLen; i++)
            dangler->addFrag(splitFrags[i], splitOffset, logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS));

          unitigs.push_back(dangler);
          tig = unitigs[ti];
        }

        //  We just split out these fragments.  Reset the list.
        splitFragsLen = 0;
      }  //  End break

      splitFrags[splitFragsLen++] = *frg;

      maxEnd = MAX(maxEnd, MAX(frg->position.bgn, frg->position.end));
    }  //  End of unitig fragment iteration

    //  If we split this unitig, the length of the
    //  frags in splitFrags will be less than the length of
    //  the path in this unitg.  If so, rebuild this unitig.
    //
    if (splitFragsLen != tig->ufpath.size()) {

      if (logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS))
        fprintf(logFile, "Rebuild unitig %d\n", tig->id());

      tig->ufpath.clear();

      int splitOffset = -MIN(splitFrags[0].position.bgn, splitFrags[0].position.end);

      //  This should already be true, but we force it still
      splitFrags[0].contained = 0;

      for (uint32 i=0; i<splitFragsLen; i++)
        tig->addFrag(splitFrags[i], splitOffset, logFileFlagSet(LOG_MATE_SPLIT_DISCONTINUOUS));
    }

    delete [] splitFrags;
    splitFrags    = NULL;
  }  //  End of discontinuity splitting
}


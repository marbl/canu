
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

static const char *rcsid = "$Id: AS_BAT_UnitigGraph.C,v 1.1 2010-11-24 01:03:31 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_ChunkGraph.H"
#include "AS_BAT_UnitigGraph.H"
#include "AS_BAT_InsertSizes.H"
#include "AS_BAT_MateLocation.H"

#include "MultiAlignStore.h"


UnitigGraph::UnitigGraph() {
}


UnitigGraph::~UnitigGraph() {
  for (uint32  ti=0; ti<unitigs.size(); ti++)
    delete unitigs[ti];
}


void
UnitigGraph::build(OverlapStore *ovlStoreUniq,
                   OverlapStore *ovlStoreRept,
                   bool enableIntersectionBreaking,
                   bool enableJoining,
                   bool enableBubblePopping,
                   int32 badMateBreakThreshold,
                   char *output_prefix) {

  // Initialize where we've been to nowhere
  Unitig::resetFragUnitigMap(FI->numFragments());

  //  There is no 0th unitig.
  unitigs.push_back(NULL);


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Build the initial unitig path from non-contained fragments.  The first pass is usually the
  //  only one needed, but occasionally (maybe) we miss fragments, so we make an explicit pass
  //  through all fragments and place whatever isn't already placed.
  //

  setLogFile(output_prefix, "buildUnitigs");
  fprintf(logFile, "==> BUILDING UNITIGS from %d fragments.\n", FI->numFragments());

  for (uint32 fi=CG->nextFragByChunkLength(); fi>0; fi=CG->nextFragByChunkLength())
    populateUnitig(fi);

  //setLogFile(output_prefix, "buildUnitigs-MissedFragments");
  fprintf(logFile, "==> BUILDING UNITIGS catching missed fragments.\n");

  for (uint32 fi=1; fi <= FI->numFragments(); fi++)
    populateUnitig(fi);

  reportOverlapsUsed(output_prefix, "buildUnitigs");
  reportUnitigs(output_prefix, "buildUnitigs");
  evaluateMates(output_prefix, "buildUnitigs");


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Place contained fragments.  There is a variety of evidence we can use to determine
  //  where to place these:
  //
  //    1) By using the best containment overlap (no way to measure if this is an uncontested
  //    placement)
  //
  //    2) By using all overlaps (we can measure if the placement is unique)
  //
  //    3) By using any mate relationship to a non-contained fragment
  //
  //    4) By using any mate relationship to an unambiguously placed (#2) contained fragment
  //
  //  For now, we're using the standard method #1.
  //
  //
  //  After placing contains, every fragment should be in a unitig.  Sometimes this
  //  is not true and we have zombie fragments (not dead, but not in a unitig).  We
  //  place any of those into their own unitig.
  //

  setLogFile(output_prefix, "placeContainsZombies");

  placeContainsUsingBestOverlaps();
  //placeContainsUsingAllOverlaps(ovlStoreUniq,
  //                              ovlStoreRept,
  //                              bool withMatesToNonContained,
  //                              bool withMatesToUnambiguousContain);

  placeZombies();

  checkUnitigMembership();
  reportOverlapsUsed(output_prefix, "placeContainsZombies");
  reportUnitigs(output_prefix, "placeContainsZombies");
  evaluateMates(output_prefix, "placeContainsZombies");


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Search for bubbles/staples.  Two methods so far.
  //
  //  1) By using best overlaps off the end of a small unitig that intersect a larger unitig to give
  //  an approximate placement, then using overlaps to place fragments into the larger unitig.
  //
  //  2) By using mates to give an approximate placement, then using overlaps to place fragments.
  //
  //  3) By blindly searching.  For each small unitig, attempt to place it using overlaps.  If the
  //  placement is consistent, the bubble is popped.  Works if there is no mate and no intersection.

  setLogFile(output_prefix, "bubblePopping");

  popIntersectionBubbles(ovlStoreUniq, ovlStoreRept);  //  Well supported as of Wed 13 Oct
  popMateBubbles(ovlStoreUniq, ovlStoreRept);          //  Exploratory as of Wed 13 Oct
  popOverlapBubbles(ovlStoreUniq, ovlStoreRept);       //  NOP as of Wed 13 Oct

  checkUnitigMembership();
  reportOverlapsUsed(output_prefix, "bubblePopping");
  reportUnitigs(output_prefix, "bubblePopping");
  evaluateMates(output_prefix, "bubblePopping");

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Search for intersection breaks and joins.  These are cases where we extended a BOG unitig into
  //  a spur, and fell off the best overlap path.  We hope to trim off the spur end, join the two
  //  unitigs into one, and maybe fold back the spur as a bubble.
  //
  //  We don't know how to handle contains when we split a unitig.  The splitting proceedure removes
  //  contains from unitigs that are split (it does NOT put them in singleton unitigs, it just
  //  removes them as if they were never placed), so after splitting we must do placeContains()
  //  again.

  setLogFile(output_prefix, "intersectionBreaking");

  breakUnitigs(output_prefix, enableIntersectionBreaking);
  placeContainsUsingBestOverlaps();

  checkUnitigMembership();
  reportOverlapsUsed(output_prefix, "intersectionBreaking");
  reportUnitigs(output_prefix, "intersectionBreaking");
  evaluateMates(output_prefix, "intersectionBreaking");

  setLogFile(output_prefix, "joinUnitigs");

  joinUnitigs(enableJoining);
  placeContainsUsingBestOverlaps();

  checkUnitigMembership();
  reportOverlapsUsed(output_prefix, "joinUnitigs");
  reportUnitigs(output_prefix, "joinUnitigs");
  evaluateMates(output_prefix, "joinUnitigs");

  //  Maybe another round of bubble popping as above.

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  This is the MateChecker
  //
  ////////////////////////////////////////////////////////////////////////////////

  //  Move unhappy contained fragments before attempting mate based
  //  splitting.  We know they're broken already, and if scaffolder
  //  can put them back, we'll let it.

  setLogFile(output_prefix, "moveContains1");
  fprintf(logFile, "==> MOVE CONTAINS #1\n");
  moveContains();

  reportOverlapsUsed(output_prefix, "moveContains1");
  checkUnitigMembership();
  evaluateMates(output_prefix, "moveContains1");

  //  This should do absolutely nothing.  If it does, something is
  //  broken.  By just ejecting unhappy contains, nothing should be
  //  disconnected.

  setLogFile(output_prefix, "splitDiscontinuous1");
  fprintf(logFile, "==> SPLIT DISCONTINUOUS #1\n");
  splitDiscontinuousUnitigs();

  reportOverlapsUsed(output_prefix, "splitDiscontinuous1");
  checkUnitigMembership();
  evaluateMates(output_prefix, "splitDiscontinuous1");

  setLogFile(output_prefix, "splitBadMates");
  fprintf(logFile, "==> SPLIT BAD MATES\n");

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) || (tig->getNumFrags() < 2))
      continue;

    MateLocation      *mates  = new MateLocation(tig);
    UnitigBreakPoints *breaks = computeMateCoverage(tig, mates, badMateBreakThreshold);
    UnitigVector      *newUs  = breakUnitigAt(tig, *breaks);

    if (logFileFlagSet(LOG_MATE_SPLIT_COVERAGE_PLOT))
      mates->dumpHappiness(output_prefix, "splitBadMates");

    if (newUs != NULL) {
      delete tig;
      unitigs[ti] = NULL;
      unitigs.insert(unitigs.end(), newUs->begin(), newUs->end());
    }

    delete newUs;
    delete breaks;
    delete mates;
  }

  placeContainsUsingBestOverlaps();

  reportOverlapsUsed(output_prefix, "splitBadMates");
  checkUnitigMembership();
  evaluateMates(output_prefix, "splitBadMates");

  //  The splitting code above is not smart enough to move contained
  //  fragments with containers.  This leaves unitigs disconnected.
  //  We break those unitigs here.

  setLogFile(output_prefix, "splitDiscontinuous2");
  fprintf(logFile, "==> SPLIT DISCONTINUOUS #2\n");
  splitDiscontinuousUnitigs();

  reportOverlapsUsed(output_prefix, "splitDiscontinuous2");
  checkUnitigMembership();
  evaluateMates(output_prefix, "splitDiscontinuous2");

  //  But now, all the splitting probably screwed up happiness of
  //  contained fragments, of left some unhappy fragments in a unitig
  //  that just lost the container.

  setLogFile(output_prefix, "moveContains2");
  fprintf(logFile, "==> MOVE CONTAINS #2\n");

  moveContains();

  reportOverlapsUsed(output_prefix, "moveContains2");
  checkUnitigMembership();
  evaluateMates(output_prefix, "moveContains2");

  //  Do one last check for disconnected unitigs.

  setLogFile(output_prefix, "splitDiscontinuous3");
  fprintf(logFile, "==> SPLIT DISCONTINUOUS #3\n");
  splitDiscontinuousUnitigs();

  reportOverlapsUsed(output_prefix, "splitDiscontinuous3");
  checkUnitigMembership();
  evaluateMates(output_prefix, "splitDiscontinuous3");

  setLogFile(output_prefix, "setParentAndHang");
  setParentAndHang();

  setLogFile(output_prefix, NULL);
}

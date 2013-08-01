
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

static const char *rcsid = "$Id: AS_BOG_UnitigGraph.cc,v 1.142 2012-02-02 09:11:21 brianwalenz Exp $";

#include "AS_BOG_Datatypes.H"
#include "AS_BOG_BestOverlapGraph.H"
#include "AS_BOG_ChunkGraph.H"
#include "AS_BOG_UnitigGraph.H"
#include "AS_BOG_InsertSizes.H"
#include "AS_BOG_MateLocation.H"

#include "MultiAlignStore.H"


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

  // Invert the containment map to key by container, instead of containee
  ContainerMap      cMap;

  for (uint32 c=0; c<=FI->numFragments(); c++)
    if (OG->isContained(c))
      cMap[OG->getBestContainer(c)->container].push_back(c);

  // Step through all the fragments

  setLogFile(output_prefix, "buildUnitigs");
  fprintf(logFile, "==> BUILDING UNITIGS from %d fragments.\n", FI->numFragments());

  //  There is no 0th unitig.
  unitigs.push_back(NULL);

  uint32 frag_idx;
  while ((frag_idx = CG->nextFragByChunkLength()) > 0) {
    if ((FI->fragmentLength(frag_idx) == 0) ||
        (Unitig::fragIn(frag_idx) != 0) ||
        (OG->isContained(frag_idx) == true))
      //  Skip deleted fragments, already placed fragments, and contained fragments.
      continue;

    populateUnitig(frag_idx);
  }

  //setLogFile(output_prefix, "buildUnitigs-MissedFragments");
  fprintf(logFile, "==> BUILDING UNITIGS catching missed fragments.\n");

  //  Pick up frags missed above, leftovers from possibly circular unitigs

  for (frag_idx=1; frag_idx <= FI->numFragments(); frag_idx++) {
    if ((FI->fragmentLength(frag_idx) == 0) ||
        (Unitig::fragIn(frag_idx) != 0) ||
        (OG->isContained(frag_idx) == true))
      continue;

    populateUnitig(frag_idx);
  }

  reportOverlapsUsed(output_prefix, "buildUnitigs");
  reportUnitigs(output_prefix, "buildUnitigs");

  //  DON'T POP BUBBLES
  //goto skipSteps;

  if (enableBubblePopping) {
    setLogFile(output_prefix, "bubblePopping");
    popIntersectionBubbles(ovlStoreUniq, ovlStoreRept);
    //popMateBubbles(ovlStoreUniq, ovlStoreRept);  NEED TO HAVE CONTAINS PLACED
    reportOverlapsUsed(output_prefix, "bubblePopping");
    reportUnitigs(output_prefix, "bubblePopping");
  }

  //  DON'T BREAK INTERSECTIONS
  //goto skipSteps;

  //  If enabled, break unitigs.  If not enabled, report on breaks.
  setLogFile(output_prefix, "intersectionBreaking");
  breakUnitigs(cMap, output_prefix, enableIntersectionBreaking);
  reportOverlapsUsed(output_prefix, "intersectionBreaking");
  reportUnitigs(output_prefix, "intersectionBreaking");

  //  DON'T REJOIN INTERSECTIONS
  //goto skipSteps;

  //  If enabled, join unitigs.  If not enabled, report on joins.
  //  (not supported, just crashes if called and not enabled)
  if (enableJoining) {
    setLogFile(output_prefix, "joinUnitigs");
    joinUnitigs(enableJoining);
    reportOverlapsUsed(output_prefix, "joinUnitigs");
    reportUnitigs(output_prefix, "joinUnitigs");
  }

  //  This will only analyze and report potential breaks.
  //setLogFile(output_prefix, "intersectionBreaking");
  //breakUnitigs(cMap, output_prefix, false);

 skipSteps:
  setLogFile(output_prefix, "placeContains");
  placeContainsUsingBestOverlaps();
  reportOverlapsUsed(output_prefix, "placeContains");
  reportUnitigs(output_prefix, "placeContains");

  setLogFile(output_prefix, "placeZombies");
  placeZombies();
  reportOverlapsUsed(output_prefix, "placeZombies");
  reportUnitigs(output_prefix, "placeZombies");

  //checkUnitigMembership();

  //  DONT DO BUBBLES OR SPLITTING
  //return;

  if (enableBubblePopping) {
    setLogFile(output_prefix, "bubblePopping");
    popIntersectionBubbles(ovlStoreUniq, ovlStoreRept);
    //popMateBubbles(ovlStoreUniq, ovlStoreRept);
    reportOverlapsUsed(output_prefix, "bubblePopping");
    reportUnitigs(output_prefix, "bubblePopping");
  }

  checkUnitigMembership();

  //  DON'T DO MATE BASED SPLITTING
  //return;

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  This is the MateChecker
  //
  ////////////////////////////////////////////////////////////////////////////////

  setLogFile(output_prefix, "libraryStats");
  IS = new InsertSizes();

  setLogFile(output_prefix, "evaluateMates");
  fprintf(logFile, "==> STARTING MATE BASED SPLITTING.\n");

  checkUnitigMembership();
  evaluateMates();

  //  Move unhappy contained fragments before attempting mate based
  //  splitting.  We know they're broken already, and if scaffolder
  //  can put them back, we'll let it.

  setLogFile(output_prefix, "moveContains1");
  fprintf(logFile, "==> MOVE CONTAINS #1\n");
  moveContains();

  reportOverlapsUsed(output_prefix, "moveContains1");
  checkUnitigMembership();
  evaluateMates();

  //  This should do absolutely nothing.  If it does, something is
  //  broken.  By just ejecting unhappy contains, nothing should be
  //  disconnected.

  setLogFile(output_prefix, "splitDiscontinuous1");
  fprintf(logFile, "==> SPLIT DISCONTINUOUS #1\n");
  splitDiscontinuousUnitigs();

  reportOverlapsUsed(output_prefix, "splitDiscontinuous2");
  checkUnitigMembership();
  evaluateMates();

  //  DON'T DO MATE BASED SPLITTING
  //return;

  setLogFile(output_prefix, "splitBadMates");
  fprintf(logFile, "==> SPLIT BAD MATES\n");

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) || (tig->getNumFrags() < 2))
      continue;

    UnitigBreakPoints* breaks = computeMateCoverage(tig, badMateBreakThreshold);
    UnitigVector*      newUs  = breakUnitigAt(tig, *breaks);

    if (newUs != NULL) {
      delete tig;
      unitigs[ti] = NULL;
      unitigs.insert(unitigs.end(), newUs->begin(), newUs->end());
    }

    delete newUs;
    delete breaks;
  }

  reportOverlapsUsed(output_prefix, "splitBadMates");
  checkUnitigMembership();
  evaluateMates();

  //  DO MATE BASED SPLITTING BUT NOTHING ELSE
  //return;

  //  The splitting code above is not smart enough to move contained
  //  fragments with containers.  This leaves unitigs disconnected.
  //  We break those unitigs here.

  setLogFile(output_prefix, "splitDiscontinuous2");
  fprintf(logFile, "==> SPLIT DISCONTINUOUS #2\n");
  splitDiscontinuousUnitigs();

  reportOverlapsUsed(output_prefix, "splitDiscontinuous2");
  checkUnitigMembership();
  evaluateMates();

  //  But now, all the splitting probably screwed up happiness of
  //  contained fragments, of left some unhappy fragments in a unitig
  //  that just lost the container.

  setLogFile(output_prefix, "moveContains2");
  fprintf(logFile, "==> MOVE CONTAINS #2\n");
  moveContains();

  reportOverlapsUsed(output_prefix, "moveContains2");
  checkUnitigMembership();
  evaluateMates();

  //  Do one last check for disconnected unitigs.

  setLogFile(output_prefix, "splitDiscontinuous3");
  fprintf(logFile, "==> SPLIT DISCONTINUOUS #3\n");
  splitDiscontinuousUnitigs();

  reportOverlapsUsed(output_prefix, "splitDiscontinuous3");
  checkUnitigMembership();
  evaluateMates();

  setLogFile(output_prefix, NULL);

  setParentAndHang();
}

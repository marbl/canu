
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

static const char *rcsid = "$Id: AS_BOG_UnitigGraph.cc,v 1.140 2010-10-27 04:15:06 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_InsertSizes.hh"
#include "AS_BOG_MateLocation.hh"

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

  // Invert the containment map to key by container, instead of containee
  ContainerMap      cMap;

  for (uint32 c=0; c<=FI->numFragments(); c++)
    if (OG->isContained(c))
      cMap[OG->getBestContainer(c)->container].push_back(c);

  // Step through all the fragments

  setLogFile("unitigger", "buildUnitigs");
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

  setLogFile("unitigger", "buildUnitigs-MissedFragments");
  fprintf(logFile, "==> BUILDING UNITIGS catching missed fragments.\n");

  //  Pick up frags missed above, leftovers from possibly circular unitigs

  for (frag_idx=1; frag_idx <= FI->numFragments(); frag_idx++) {
    if ((FI->fragmentLength(frag_idx) == 0) ||
        (Unitig::fragIn(frag_idx) != 0) ||
        (OG->isContained(frag_idx) == true))
      continue;

    populateUnitig(frag_idx);
  }

  reportOverlapsUsed("overlaps.afterbuild");
  reportUnitigs("unitigs.afterbuild");

  if (enableBubblePopping) {
    setLogFile("unitigger", "bubblePopping");
    popIntersectionBubbles(ovlStoreUniq, ovlStoreRept);
    //popMateBubbles(ovlStoreUniq, ovlStoreRept);  NEED TO HAVE CONTAINS PLACED
    reportOverlapsUsed("overlaps.afterbubbles1");
    reportUnitigs("unitigs.afterbubbles1");
  }

  //  If enabled, break unitigs.  If not enabled, report on breaks.
  setLogFile("unitigger", "intersectionBreaking");
  breakUnitigs(cMap, output_prefix, enableIntersectionBreaking);
  reportOverlapsUsed("overlaps.afterbreak");
  reportUnitigs("unitigs.afterbreak");

  //  If enabled, join unitigs.  If not enabled, report on joins.
  //  (not supported, just crashes if called and not enabled)
  if (enableJoining) {
    setLogFile("unitigger", "joinUnitigs");
    joinUnitigs(enableJoining);
    reportOverlapsUsed("overlaps.afterjoin");
    reportUnitigs("unitigs.afterjoin");
  }

  //  This will only analyze and report potential breaks.
  //setLogFile("unitigger", "intersectionBreaking");
  //breakUnitigs(cMap, output_prefix, false);

#if 0
  for (uint32 fid=1; fid<FI->numFragments()+1; fid++) {
    ufNode frag;

    if (Unitig::fragIn(fid) > 0)
      //  Fragment placed already.
      continue;

    frag.ident = fid;

    placeFragUsingOverlaps(frag, ovlStoreUniq, ovlStoreRept);
  }
#endif

  setLogFile("unitigger", "placeContains");
  placeContainsUsingBestOverlaps();
  reportOverlapsUsed("overlaps.aftercontains");
  reportUnitigs("unitigs.aftercontains");

#if 0
  for (uint32 fid=1; fid<FI->numFragments()+1; fid++) {
    ufNode frag;

    if (OG->getBestContainer(fid) == NULL)
      //  Dovetail node, don't try to re-place it
      continue;

    frag.ident = fid;

    placeFragUsingOverlaps(frag, ovlStoreUniq, ovlStoreRept);
  }
#endif

  setLogFile("unitigger", "placeZombies");
  placeZombies();
  reportOverlapsUsed("overlaps.afterzombies");
  reportUnitigs("unitigs.afterzombies");

  //checkUnitigMembership();

  if (enableBubblePopping) {
    setLogFile("unitigger", "bubblePopping");
    popIntersectionBubbles(ovlStoreUniq, ovlStoreRept);
    //popMateBubbles(ovlStoreUniq, ovlStoreRept);
    reportOverlapsUsed("overlaps.afterbubbles2");
    reportUnitigs("unitigs.afterbubbles2");
  }

  checkUnitigMembership();

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  This is the MateChecker
  //
  ////////////////////////////////////////////////////////////////////////////////

  setLogFile("unitigger", "libraryStats");
  IS = new InsertSizes();

  setLogFile("unitigger", "evaluateMates");
  fprintf(logFile, "==> STARTING MATE BASED SPLITTING.\n");

  checkUnitigMembership();
  evaluateMates();

  //  Move unhappy contained fragments before attempting mate based
  //  splitting.  We know they're broken already, and if scaffolder
  //  can put them back, we'll let it.

  setLogFile("unitigger", "moveContains1");
  fprintf(logFile, "==> MOVE CONTAINS #1\n");
  moveContains();

  reportOverlapsUsed("overlaps.aftermatecheck1");
  checkUnitigMembership();
  evaluateMates();

  //  This should do absolutely nothing.  If it does, something is
  //  broken.  By just ejecting unhappy contains, nothing should be
  //  disconnected.

  setLogFile("unitigger", "splitDiscontinuous1");
  fprintf(logFile, "==> SPLIT DISCONTINUOUS #1\n");
  splitDiscontinuousUnitigs();

  reportOverlapsUsed("overlaps.aftermatecheck2");
  checkUnitigMembership();
  evaluateMates();

  //  DON'T DO MATE BASED SPLITTING
  //return;

  setLogFile("unitigger", "splitBadMates");
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

  reportOverlapsUsed("overlaps.aftermatecheck3");
  checkUnitigMembership();
  evaluateMates();

  //  DO MATE BASED SPLITTING BUT NOTHING ELSE
  //return;

  //  The splitting code above is not smart enough to move contained
  //  fragments with containers.  This leaves unitigs disconnected.
  //  We break those unitigs here.

  setLogFile("unitigger", "splitDiscontinuous2");
  fprintf(logFile, "==> SPLIT DISCONTINUOUS #2\n");
  splitDiscontinuousUnitigs();

  reportOverlapsUsed("overlaps.aftermatecheck4");
  checkUnitigMembership();
  evaluateMates();

  //  But now, all the splitting probably screwed up happiness of
  //  contained fragments, of left some unhappy fragments in a unitig
  //  that just lost the container.

  setLogFile("unitigger", "moveContains2");
  fprintf(logFile, "==> MOVE CONTAINS #2\n");
  moveContains();

  reportOverlapsUsed("overlaps.aftermatecheck5");
  checkUnitigMembership();
  evaluateMates();

  //  Do one last check for disconnected unitigs.

  setLogFile("unitigger", "splitDiscontinuous3");
  fprintf(logFile, "==> SPLIT DISCONTINUOUS #3\n");
  splitDiscontinuousUnitigs();

  reportOverlapsUsed("overlaps.aftermatecheck6");
  checkUnitigMembership();
  evaluateMates();

  setLogFile("unitigger", NULL);

  setParentAndHang();
}

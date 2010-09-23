
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

static const char *rcsid = "$Id: AS_BOG_UnitigGraph.cc,v 1.135 2010-09-23 09:34:50 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"

#undef max

//  Logging
bool verboseBuild    = false;  //  Dovetail path construction
bool verboseMerge    = false;  //  Bubbles
bool verboseBreak    = false;  //  Intersection AND mate-based breaking
bool verboseJoin     = false;  //  Joining
bool verboseContains = false;  //  Containment placing


UnitigGraph::UnitigGraph(FragmentInfo *fi, BestOverlapGraph *bp) {
  unitigs = new UnitigVector;
  _fi     = fi;
  bog_ptr = bp;
}

UnitigGraph::~UnitigGraph() {
  for (int  ti=0; ti<unitigs->size(); ti++)
    delete (*unitigs)[ti];
  delete unitigs;
}


void UnitigGraph::build(ChunkGraph *cg_ptr,
                        OverlapStore *ovlStoreUniq,
                        OverlapStore *ovlStoreRept,
                        bool enableIntersectionBreaking,
                        bool enableJoining,
                        bool enableBubblePopping,
                        char *output_prefix) {

  // Initialize where we've been to nowhere
  Unitig::resetFragUnitigMap( _fi->numFragments() );

  // Invert the containment map to key by container, instead of containee
  ContainerMap      cMap;

  for (uint32 c=0; c<=_fi->numFragments(); c++)
    if (bog_ptr->isContained(c))
      cMap[bog_ptr->getBestContainer(c)->container].push_back(c);

  // Step through all the fragments

  fprintf(stderr, "==> BUILDING UNITIGS from %d fragments.\n", _fi->numFragments());

  //  There is no 0th unitig.
  unitigs->push_back(NULL);

  uint32 frag_idx;
  while ((frag_idx = cg_ptr->nextFragByChunkLength()) > 0) {
    if ((_fi->fragmentLength(frag_idx) == 0) ||
        (Unitig::fragIn(frag_idx) != 0) ||
        (bog_ptr->isContained(frag_idx) == true))
      //  Skip deleted fragments, already placed fragments, and contained fragments.
      continue;

    populateUnitig(frag_idx);
  }

  fprintf(stderr, "==> BUILDING UNITIGS catching missed fragments.\n");

  //  Pick up frags missed above, leftovers from possibly circular unitigs

  for (frag_idx=1; frag_idx <= _fi->numFragments(); frag_idx++) {
    if ((_fi->fragmentLength(frag_idx) == 0) ||
        (Unitig::fragIn(frag_idx) != 0) ||
        (bog_ptr->isContained(frag_idx) == true))
      continue;

    populateUnitig(frag_idx);
  }

  reportOverlapsUsed("overlaps.afterbuild");
  reportUnitigs("unitigs.afterbuild");

  if (enableBubblePopping) {
    popIntersectionBubbles(ovlStoreUniq, ovlStoreRept);
    //popMateBubbles(ovlStoreUniq, ovlStoreRept);  NEED TO HAVE CONTAINS PLACED
  }

  //  If enabled, break unitigs.  If not enabled, report on breaks.
  breakUnitigs(cMap, output_prefix, enableIntersectionBreaking);
  reportOverlapsUsed("overlaps.afterbreak");
  reportUnitigs("unitigs.afterbreak");

  //  If enabled, join unitigs.  If not enabled, report on joins.
  //  (not supported, just crashes if called and not enabled)
  if (enableJoining) {
    joinUnitigs(enableJoining);
    reportOverlapsUsed("overlaps.afterjoin");
    reportUnitigs("unitigs.afterjoin");
  }

  //  This will only analyze and report potential breaks.
  //breakUnitigs(cMap, output_prefix, false);

  placeContains();
  reportOverlapsUsed("overlaps.aftercontains");
  reportUnitigs("unitigs.aftercontains");

  placeZombies();

  checkUnitigMembership();

  if (enableBubblePopping) {
    popIntersectionBubbles(ovlStoreUniq, ovlStoreRept);
    //popMateBubbles(ovlStoreUniq, ovlStoreRept);
  }

  checkUnitigMembership();
}

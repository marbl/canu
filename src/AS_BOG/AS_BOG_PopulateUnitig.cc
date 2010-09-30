
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

static const char *rcsid = "$Id: AS_BOG_PopulateUnitig.cc,v 1.4 2010-09-30 05:50:17 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"



void
UnitigGraph::populateUnitig(Unitig           *unitig,
                            BestEdgeOverlap  *bestnext) {

  assert(unitig->getLength() > 0);

  if ((bestnext == NULL) || (bestnext->frag_b_id == 0))
    //  Nothing to add!
    return;

  DoveTailNode frag = unitig->dovetail_path_ptr->back();

  //  The ID of the last fragment in the unitig, and the end we should walk off of it.
  int32 lastID  = frag.ident;
  int32 lastEnd = (frag.position.bgn < frag.position.end) ? THREE_PRIME : FIVE_PRIME;

                   
  if (Unitig::fragIn(bestnext->frag_b_id) == unitig->id())
    //  Cicrular unitig.  Deal with later.
    return;

  if (Unitig::fragIn(bestnext->frag_b_id) != 0) {
    //  Intersection.  Remember.
    if (logFileFlagSet(LOG_INTERSECTIONS))
      fprintf(logFile,"unitigIntersect: unitig %d frag %d -> unitig %d frag %d (before construction)\n",
              unitig->id(), lastID, Unitig::fragIn(bestnext->frag_b_id), bestnext->frag_b_id);
    unitigIntersect[bestnext->frag_b_id].push_back(lastID);
    return;
  }

  //  While there are fragments to add, construct a reverse-edge, and add the fragment.

  while (bestnext->frag_b_id != 0) {
    BestEdgeOverlap  bestprev;

    //  Reverse nextedge (points from the unitig to the next fragment to add) so that it points from
    //  the next fragment to add back to something in the unitig.  If the fragments are
    //  innie/outtie, we need to reverse the overlap to maintain that the A fragment is forward.

    bestprev.frag_b_id = lastID;
    bestprev.bend      = lastEnd;
    bestprev.ahang     = -bestnext->ahang;
    bestprev.bhang     = -bestnext->bhang;

    if (bestprev.bend == bestnext->bend) {
      bestprev.ahang = bestnext->bhang;
      bestprev.bhang = bestnext->ahang;
    }

    //  Call the usual placement routine to place the next fragment relative to the last one.  This
    //  call depends on which end of the frag-to-be-added we are working with.

    frag.ident = bestnext->frag_b_id;

    int32  bidx5 = -1, bidx3 = -1;

    if (unitig->placeFrag(frag, bidx5, (bestnext->bend == THREE_PRIME) ? NULL : &bestprev,
                          frag, bidx3, (bestnext->bend == THREE_PRIME) ? &bestprev : NULL)) {
      unitig->addFrag(frag, 0, logFileFlagSet(LOG_POPULATE_UNITIG));

    } else {

      fprintf(logFile, "ERROR:  Failed to place frag %d into BOG path.\n", frag.ident);
      assert(0);
    }

    //  Set up for the next fragmnet

    lastID  = frag.ident;
    lastEnd = (frag.position.bgn < frag.position.end) ? THREE_PRIME : FIVE_PRIME;

    bestnext = OG->getBestEdgeOverlap(lastID, lastEnd);

    //  Abort if we intersect, or are circular.  Save the intersection, but not the circularity.

    if (Unitig::fragIn(bestnext->frag_b_id) == unitig->id()) {
      if (Unitig::pathPosition(bestnext->frag_b_id) == 0) {
        if (logFileFlagSet(LOG_INTERSECTIONS))
          fprintf(logFile,"unitigIntersect: unitig %d frag %d -> unitig %d frag %d (CIRCULAR)\n",
                  unitig->id(), lastID, Unitig::fragIn(bestnext->frag_b_id), bestnext->frag_b_id);
      } else {
        if (logFileFlagSet(LOG_INTERSECTIONS))
          fprintf(logFile,"unitigIntersect: unitig %d frag %d -> unitig %d frag %d (SELF)\n",
                  unitig->id(), lastID, Unitig::fragIn(bestnext->frag_b_id), bestnext->frag_b_id);
        unitigIntersect[bestnext->frag_b_id].push_back(lastID);
        selfIntersect[lastID] = true;
      }
      break;
    }

    if (Unitig::fragIn(bestnext->frag_b_id) != 0) {
      if (logFileFlagSet(LOG_INTERSECTIONS))
        fprintf(logFile,"unitigIntersect: unitig %d frag %d -> unitig %d frag %d (during construction)\n",
                unitig->id(), lastID, Unitig::fragIn(bestnext->frag_b_id), bestnext->frag_b_id);
      unitigIntersect[bestnext->frag_b_id].push_back(lastID);
      break;
    }
  }
}




void
UnitigGraph::populateUnitig(int32 frag_idx) {

  Unitig *utg = new Unitig(logFileFlagSet(LOG_POPULATE_UNITIG));

  unitigs.push_back(utg);

  //  Add a first fragment -- to be 'compatable' with the old code, the first fragment is added
  //  reversed, we walk off of its 5' end, flip it, and add the 3' walk.

  DoveTailNode  frag;

  frag.ident             = frag_idx;
  frag.contained         = 0;
  frag.parent            = 0;
  frag.ahang             = 0;
  frag.bhang             = 0;
  frag.position.bgn      = FI->fragmentLength(frag_idx);
  frag.position.end      = 0;
  frag.containment_depth = 0;

  utg->addFrag(frag, 0, logFileFlagSet(LOG_POPULATE_UNITIG));

  //  Add fragments as long as there is a path to follow...from the 3' end of the first fragment.

  BestEdgeOverlap  *bestedge5 = OG->getBestEdgeOverlap(frag_idx, FIVE_PRIME);
  BestEdgeOverlap  *bestedge3 = OG->getBestEdgeOverlap(frag_idx, THREE_PRIME);

  if (logFileFlagSet(LOG_POPULATE_UNITIG))
    fprintf(logFile, "Adding 5' edges off of frag %d in unitig %d\n",
            utg->dovetail_path_ptr->back().ident, utg->id());

  if (bestedge5->frag_b_id)
    populateUnitig(utg, bestedge5);

  utg->reverseComplement(false);

  if (logFileFlagSet(LOG_POPULATE_UNITIG))
    fprintf(logFile, "Adding 3' edges off of frag %d in unitig %d\n",
            utg->dovetail_path_ptr->back().ident, utg->id());

  if (bestedge3->frag_b_id)
    populateUnitig(utg, bestedge3);

  //  Enabling this reverse complement is known to degrade the assembly.  It is not known WHY it
  //  degrades the assembly.
  //
#warning Reverse complement degrades assemblies
  //utg->reverseComplement(false);
}


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

static const char *rcsid = "$Id: AS_BOG_PopulateUnitig.cc,v 1.6 2010-10-11 03:43:44 brianwalenz Exp $";

#include "AS_BOG_Datatypes.H"
#include "AS_BOG_UnitigGraph.H"
#include "AS_BOG_BestOverlapGraph.H"

#include "MultiAlignStore.H"



void
UnitigGraph::populateUnitig(Unitig           *unitig,
                            BestEdgeOverlap  *bestnext) {

  assert(unitig->getLength() > 0);

  if ((bestnext == NULL) || (bestnext->fragId() == 0))
    //  Nothing to add!
    return;

  ufNode frag = unitig->ufpath.back();

  //  The ID of the last fragment in the unitig, and the end we should walk off of it.
  int32 lastID  = frag.ident;
  bool  last3p  = (frag.position.bgn < frag.position.end);


  if (Unitig::fragIn(bestnext->fragId()) == unitig->id())
    //  Cicrular unitig.  Deal with later.
    return;

  if (Unitig::fragIn(bestnext->fragId()) != 0) {
    //  Intersection.  Remember.
    if (logFileFlagSet(LOG_INTERSECTIONS))
      fprintf(logFile,"unitigIntersect: unitig %d frag %d -> unitig %d frag %d (before construction)\n",
              unitig->id(), lastID, Unitig::fragIn(bestnext->fragId()), bestnext->fragId());
    unitigIntersect[bestnext->fragId()].push_back(lastID);
    return;
  }

  //  While there are fragments to add, construct a reverse-edge, and add the fragment.

  while (bestnext->fragId() != 0) {
    BestEdgeOverlap  bestprev;

    //  Reverse nextedge (points from the unitig to the next fragment to add) so that it points from
    //  the next fragment to add back to something in the unitig.  If the fragments are
    //  innie/outtie, we need to reverse the overlap to maintain that the A fragment is forward.

    if (last3p == bestnext->frag3p())
      bestprev.set(lastID, last3p, bestnext->bhang(), bestnext->ahang());
    else
      bestprev.set(lastID, last3p, -bestnext->ahang(), -bestnext->bhang());

    //  Call the usual placement routine to place the next fragment relative to the last one.  This
    //  call depends on which end of the frag-to-be-added we are working with.

    frag.ident = bestnext->fragId();

    int32  bidx5 = -1, bidx3 = -1;

    if (unitig->placeFrag(frag, bidx5, (bestnext->frag3p() ? NULL : &bestprev),
                          frag, bidx3, (bestnext->frag3p() ? &bestprev : NULL))) {
      unitig->addFrag(frag, 0, logFileFlagSet(LOG_POPULATE_UNITIG));

    } else {

      fprintf(logFile, "ERROR:  Failed to place frag %d into BOG path.\n", frag.ident);
      assert(0);
    }

    //  Set up for the next fragmnet

    lastID  = frag.ident;
    last3p  = (frag.position.bgn < frag.position.end);

    bestnext = OG->getBestEdgeOverlap(lastID, last3p);

    //  Abort if we intersect, or are circular.  Save the intersection, but not the circularity.

    if (Unitig::fragIn(bestnext->fragId()) == unitig->id()) {
      if (Unitig::pathPosition(bestnext->fragId()) == 0) {
        if (logFileFlagSet(LOG_INTERSECTIONS))
          fprintf(logFile,"unitigIntersect: unitig %d frag %d -> unitig %d frag %d (CIRCULAR)\n",
                  unitig->id(), lastID, Unitig::fragIn(bestnext->fragId()), bestnext->fragId());
      } else {
        if (logFileFlagSet(LOG_INTERSECTIONS))
          fprintf(logFile,"unitigIntersect: unitig %d frag %d -> unitig %d frag %d (SELF)\n",
                  unitig->id(), lastID, Unitig::fragIn(bestnext->fragId()), bestnext->fragId());
        unitigIntersect[bestnext->fragId()].push_back(lastID);
        selfIntersect[lastID] = true;
      }
      break;
    }

    if (Unitig::fragIn(bestnext->fragId()) != 0) {
      if (logFileFlagSet(LOG_INTERSECTIONS))
        fprintf(logFile,"unitigIntersect: unitig %d frag %d -> unitig %d frag %d (during construction)\n",
                unitig->id(), lastID, Unitig::fragIn(bestnext->fragId()), bestnext->fragId());
      unitigIntersect[bestnext->fragId()].push_back(lastID);
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

  ufNode  frag;

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

  BestEdgeOverlap  *bestedge5 = OG->getBestEdgeOverlap(frag_idx, false);
  BestEdgeOverlap  *bestedge3 = OG->getBestEdgeOverlap(frag_idx, true);

  if (logFileFlagSet(LOG_POPULATE_UNITIG))
    fprintf(logFile, "Adding 5' edges off of frag %d in unitig %d\n",
            utg->ufpath.back().ident, utg->id());

  if (bestedge5->fragId())
    populateUnitig(utg, bestedge5);

  utg->reverseComplement(false);

  if (logFileFlagSet(LOG_POPULATE_UNITIG))
    fprintf(logFile, "Adding 3' edges off of frag %d in unitig %d\n",
            utg->ufpath.back().ident, utg->id());

  if (bestedge3->fragId())
    populateUnitig(utg, bestedge3);

  //  Enabling this reverse complement is known to degrade the assembly.  It is not known WHY it
  //  degrades the assembly.
  //
#warning Reverse complement degrades assemblies
  //utg->reverseComplement(false);
}

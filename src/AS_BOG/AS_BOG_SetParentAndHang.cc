
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

static const char *rcsid = "$Id: AS_BOG_SetParentAndHang.cc,v 1.4 2010-09-30 05:50:17 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_UnitigGraph.hh"

#include "MultiAlignStore.h"

void
UnitigGraph::setParentAndHang(void) {

  for (int  ti=0; ti<unitigs.size(); ti++) {
    Unitig        *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    if (utg->dovetail_path_ptr->size() == 0)
      continue;

    //  Reset parent and hangs for everything.

    for (int fi=1; fi<utg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode *frg = &(*utg->dovetail_path_ptr)[fi];

      frg->parent       = 0;
      frg->ahang        = 0;
      frg->bhang        = 0;
    }

    //  For each fragment, set parent/hangs using the edges.

    for (int fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode *frg  = &(*utg->dovetail_path_ptr)[fi];

      //  If we're contained, gee, I sure hope the container is here!

      BestContainment *bestcont  = OG->getBestContainer(frg->ident);

      if ((bestcont) && (utg->fragIn(bestcont->container) == utg->id())) {
        int32         pi   = utg->pathPosition(bestcont->container);
        DoveTailNode *par  = &(*utg->dovetail_path_ptr)[pi];

        frg->parent = bestcont->container;

        //  The hangs assume the container is forward; adjust if not so.
        if (par->position.bgn < par->position.end) {
          frg->ahang  = bestcont->a_hang;
          frg->bhang  = bestcont->b_hang;
        } else {
          frg->ahang  = -bestcont->b_hang;
          frg->bhang  = -bestcont->a_hang;
        }

        continue;
      }

      //  Nope, not contained.  If we don't have a parent set, see if one of our best overlaps
      //  can set it.

      BestEdgeOverlap *bestedge5 = OG->getBestEdgeOverlap(frg->ident, FIVE_PRIME);
      BestEdgeOverlap *bestedge3 = OG->getBestEdgeOverlap(frg->ident, THREE_PRIME);

      if ((bestedge5->frag_b_id) && (utg->fragIn(bestedge5->frag_b_id) == utg->id())) {
        int32         pi5  = utg->pathPosition(bestedge5->frag_b_id);
        DoveTailNode *oth  = &(*utg->dovetail_path_ptr)[pi5];

        //  Consensus is expected parent/hangs to be relative to the parent fragment.  This is used
        //  ONLY to place the fragment, not to orient the fragment.  Orientation comes from the
        //  absolute positioning coordinates.
        //
        //  Interestingly, all four overlap transformations are used here.
        //
        //  The inner if tests (on fragment orientation) should be asserts, but due to imprecise
        //  layouts, they are sometimes violated:
        //    A fragment from       271-547 had a 5'overlap to something after it;
        //    the frag after was at 543-272, close enough to a tie to screw up placements
        //
        if (pi5 < fi) {
          //  We have an edge off our 5' end to something before us --> fragment MUST be forward.
          //  Flip the overlap so it is relative to the other fragment.
          if (frg->position.bgn < frg->position.end) {
            frg->parent = bestedge5->frag_b_id;
            frg->ahang  = -bestedge5->ahang;
            frg->bhang  = -bestedge5->bhang;
            assert(frg->ahang >= 0);
          }
        } else {
          //  We have an edge off our 5' end to something after us --> fragment MUST be reverse.
          //  Because our fragment is now reverse, we must reverse the overlap too.
          if (frg->position.end < frg->position.bgn) {
            oth->parent = frg->ident;
            oth->ahang  = -bestedge5->bhang;
            oth->bhang  = -bestedge5->ahang;
            assert(oth->ahang >= 0);
          }
        }
      }

      if ((bestedge3->frag_b_id) && (utg->fragIn(bestedge3->frag_b_id) == utg->id())) {
        int32         pi3  = utg->pathPosition(bestedge3->frag_b_id);
        DoveTailNode *oth  = &(*utg->dovetail_path_ptr)[pi3];

        if (pi3 < fi) {
          //  We have an edge off our 3' end to something before us --> fragment MUST be reverse.
          //  Flip the overlap so it is relative to the other fragment.
          //  Because our fragment is now reverse, we must reverse the overlap too.
          if (frg->position.end < frg->position.bgn) {
            frg->parent = bestedge3->frag_b_id;
            frg->ahang  = bestedge3->bhang;
            frg->bhang  = bestedge3->ahang;
            assert(frg->ahang >= 0);
          }
        } else {
          //  We have an edge off our 3' end to something after us --> fragment MUST be forward.
          //  This is the simplest case, the overlap is already correct.
          if (frg->position.bgn < frg->position.end) {
            oth->parent = frg->ident;
            oth->ahang  = bestedge3->ahang;
            oth->bhang  = bestedge3->bhang;
            assert(oth->ahang >= 0);
          }
        }
      }
    }
  }
}



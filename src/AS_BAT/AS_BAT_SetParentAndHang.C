
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

static const char *rcsid = "$Id: AS_BAT_SetParentAndHang.C,v 1.3 2010-12-17 09:55:56 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_SetParentAndHang.H"

void
setParentAndHang(UnitigVector &unitigs) {

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig        *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    if (utg->ufpath.size() == 0)
      continue;

    //  Reset parent and hangs for everything.

    for (uint32 fi=1; fi<utg->ufpath.size(); fi++) {
      ufNode *frg = &utg->ufpath[fi];

      frg->parent       = 0;
      frg->ahang        = 0;
      frg->bhang        = 0;
    }

    //  For each fragment, set parent/hangs using the edges.

    for (uint32 fi=0; fi<utg->ufpath.size(); fi++) {
      ufNode *frg  = &utg->ufpath[fi];

      //  If we're contained, gee, I sure hope the container is here!

      BestContainment *bestcont  = OG->getBestContainer(frg->ident);

      if (bestcont) {
        if (utg->fragIn(bestcont->container) == utg->id()) {
          int32   pi   = utg->pathPosition(bestcont->container);
          ufNode *par  = &utg->ufpath[pi];

          assert(par->ident == bestcont->container);

          frg->parent = par->ident;

          //  The hangs assume the container is forward; adjust if not so.
          if (par->position.bgn < par->position.end) {
            frg->ahang  = bestcont->a_hang;
            frg->bhang  = bestcont->b_hang;
          } else {
            frg->ahang  = -bestcont->b_hang;
            frg->bhang  = -bestcont->a_hang;
          }

          if (logFileFlags & LOG_SET_PARENT_AND_HANG)
            fprintf(logFile, "frag %d at %d,%d edge to cont frag %d at %d,%d -- hang %d,%d\n",
                    frg->ident, frg->position.bgn, frg->position.end,
                    par->ident, par->position.bgn, par->position.end,
                    frg->ahang, frg->bhang);
        } else {
          if (logFileFlags & LOG_SET_PARENT_AND_HANG)
            fprintf(logFile, "frag %d at %d,%d edge to cont frag %d in different unitig %d\n",
                    frg->ident, frg->position.bgn, frg->position.end,
                    bestcont->container, utg->fragIn(bestcont->container));
        }

        continue;
      }

      //  Nope, not contained.  If we don't have a parent set, see if one of our best overlaps
      //  can set it.

      BestEdgeOverlap *bestedge5 = OG->getBestEdgeOverlap(frg->ident, false);
      BestEdgeOverlap *bestedge3 = OG->getBestEdgeOverlap(frg->ident, true);

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

      if (bestedge5->fragId() > 0) {
        if (logFileFlags & LOG_SET_PARENT_AND_HANG)
          fprintf(logFile, "frag %d in unitig %d 5' to %d/%c' in unitig %d\n",
                  frg->ident, utg->id(),
                  bestedge5->fragId(), bestedge5->frag3p() ? '3' : '5',
                  utg->fragIn(bestedge5->fragId()));

        if (utg->fragIn(bestedge5->fragId()) == utg->id()) {
          uint32  pi5  = utg->pathPosition(bestedge5->fragId());
          ufNode *oth  = &utg->ufpath[pi5];

          assert(oth->ident == bestedge5->fragId());

          if ((pi5 < fi) && (isReverse(frg->position) == false)) {
            //  Edge is to a fragment before us, off our 5' end, and we are forward.
            frg->parent = bestedge5->fragId();
            frg->ahang  = -bestedge5->ahang();
            frg->bhang  = -bestedge5->bhang();
            assert(frg->ahang >= 0);

            if (logFileFlags & LOG_SET_PARENT_AND_HANG)
              fprintf(logFile, "->frag %d at %d,%d 5' edge to prev frag %d at %d,%d -- hang %d,%d\n",
                      frg->ident, frg->position.bgn, frg->position.end,
                      oth->ident, oth->position.bgn, oth->position.end,
                      frg->ahang, frg->bhang);

          } else if ((pi5 > fi) && (isReverse(frg->position) == true)) {
            //  Edge is to a fragment after us, off our 5' end, and we are reverse.
            //  Use this edge to set the other fragment parent and hang.
            //  That fragment must pass the same order/orient tests.
            //    Off the others 3' end, fragment must be reverse.
            //    Off the others 5' end, fragment must be forward.
            if (((bestedge5->frag3p() == true)  && (isReverse(oth->position) == true)) ||
                ((bestedge5->frag3p() == false) && (isReverse(oth->position) == false))) {
              oth->parent = frg->ident;
              oth->ahang  = -bestedge5->bhang();
              oth->bhang  = -bestedge5->ahang();
              assert(oth->ahang >= 0);

              if (logFileFlags & LOG_SET_PARENT_AND_HANG)
                fprintf(logFile, "<-frag %d at %d,%d %c' edge fr prev frag %d at %d,%d -- hang %d,%d\n",
                        oth->ident, oth->position.bgn, oth->position.end, bestedge5->frag3p() ? '3' : '5',
                        frg->ident, frg->position.bgn, frg->position.end,
                        frg->ahang, frg->bhang);
            } else {
              if (logFileFlags & LOG_SET_PARENT_AND_HANG)
                fprintf(logFile, "<-frag %d at %d,%d %c' edge fr prev frag %d at %d,%d -- NOT VALID\n",
                        oth->ident, oth->position.bgn, oth->position.end, bestedge5->frag3p() ? '3' : '5',
                        frg->ident, frg->position.bgn, frg->position.end);
            }

          } else {
            if (logFileFlags & LOG_SET_PARENT_AND_HANG)
              fprintf(logFile, "--frag %d at %d,%d 5' edge to prev frag %d at %d,%d -- NOT VALID\n",
                      frg->ident, frg->position.bgn, frg->position.end,
                      oth->ident, oth->position.bgn, oth->position.end);
          }
        }
      }

      if (bestedge3->fragId() > 0) {
        if (logFileFlags & LOG_SET_PARENT_AND_HANG)
          fprintf(logFile, "frag %d in unitig %d 3' to %d/%c' in unitig %d\n",
                  frg->ident, utg->id(),
                  bestedge3->fragId(), bestedge3->frag3p() ? '3' : '5',
                  utg->fragIn(bestedge5->fragId()));

        if (utg->fragIn(bestedge3->fragId()) == utg->id()) {
          uint32  pi3  = utg->pathPosition(bestedge3->fragId());
          ufNode *oth  = &utg->ufpath[pi3];

          assert(oth->ident == bestedge3->fragId());

          //  Edge is to a fragment before us, off our 3' end, and we are reverse.
          if        ((pi3 < fi) && (isReverse(frg->position) == true)) {
            frg->parent = oth->ident;
            frg->ahang  = bestedge3->bhang();
            frg->bhang  = bestedge3->ahang();
            assert(frg->ahang >= 0);

            if (logFileFlags & LOG_SET_PARENT_AND_HANG)
              fprintf(logFile, "->frag %d at %d,%d 3' edge to prev frag %d at %d,%d -- hang %d,%d\n",
                      frg->ident, frg->position.bgn, frg->position.end,
                      oth->ident, oth->position.bgn, oth->position.end,
                      frg->ahang, frg->bhang);

          } else if ((pi3 > fi) && (isReverse(frg->position) == false)) {
            if (((bestedge3->frag3p() == true)  && (isReverse(oth->position) == true)) ||
                ((bestedge3->frag3p() == false) && (isReverse(oth->position) == false))) {
              oth->parent = frg->ident;
              oth->ahang  = bestedge3->ahang();
              oth->bhang  = bestedge3->bhang();
              assert(oth->ahang >= 0);

              if (logFileFlags & LOG_SET_PARENT_AND_HANG)
                fprintf(logFile, "<-frag %d at %d,%d %c' edge fr prev frag %d at %d,%d -- hang %d,%d\n",
                        oth->ident, oth->position.bgn, oth->position.end, bestedge5->frag3p() ? '3' : '5',
                        frg->ident, frg->position.bgn, frg->position.end,
                        frg->ahang, frg->bhang);
            } else {
              if (logFileFlags & LOG_SET_PARENT_AND_HANG)
                fprintf(logFile, "<-frag %d at %d,%d %c' edge fr prev frag %d at %d,%d -- NOT VALID\n",
                        oth->ident, oth->position.bgn, oth->position.end, bestedge5->frag3p() ? '3' : '5',
                        frg->ident, frg->position.bgn, frg->position.end);
            }

          } else {
            if (logFileFlags & LOG_SET_PARENT_AND_HANG)
              fprintf(logFile, "--frag %d at %d,%d 3' edge to prev frag %d at %d,%d -- NOT VALID\n",
                      frg->ident, frg->position.bgn, frg->position.end,
                      oth->ident, oth->position.bgn, oth->position.end);
          }
        }
      }
    }  //  Over all fragment
  }  //  Over all unitigs
}

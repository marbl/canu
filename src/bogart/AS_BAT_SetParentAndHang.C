
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
 *    src/AS_BAT/AS_BAT_SetParentAndHang.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010,2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2015-AUG-05
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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

      if (bestcont->isContained == true) {
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
            writeLog("setParentAndHang()--  CONTAINED - frag %d at %d,%d edge to cont frag %d at %d,%d -- hang %d,%d\n",
                    frg->ident, frg->position.bgn, frg->position.end,
                    par->ident, par->position.bgn, par->position.end,
                    frg->ahang, frg->bhang);
        } else {
          if (logFileFlags & LOG_SET_PARENT_AND_HANG)
            writeLog("setParentAndHang()--  CONTAINED - frag %d at %d,%d edge to cont frag %d IN DIFFERENT UNITIG %d\n",
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
          writeLog("setParentAndHang()--  BEST5     - frag %d in unitig %d 5' to %d/%c' in unitig %d\n",
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
              writeLog("                             - -> frag %d at %d,%d 5' edge to prev frag %d at %d,%d -- hang %d,%d\n",
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
                writeLog("                                - <- frag %d at %d,%d %c' edge fr prev frag %d at %d,%d -- hang %d,%d\n",
                        oth->ident, oth->position.bgn, oth->position.end, bestedge5->frag3p() ? '3' : '5',
                        frg->ident, frg->position.bgn, frg->position.end,
                        frg->ahang, frg->bhang);
            } else {
              if (logFileFlags & LOG_SET_PARENT_AND_HANG)
                writeLog("                                - <- frag %d at %d,%d %c' edge fr prev frag %d at %d,%d -- NOT VALID\n",
                        oth->ident, oth->position.bgn, oth->position.end, bestedge5->frag3p() ? '3' : '5',
                        frg->ident, frg->position.bgn, frg->position.end);
            }

          } else {
            if (logFileFlags & LOG_SET_PARENT_AND_HANG)
              writeLog("                                - -- frag %d at %d,%d 5' edge to prev frag %d at %d,%d -- NOT VALID\n",
                      frg->ident, frg->position.bgn, frg->position.end,
                      oth->ident, oth->position.bgn, oth->position.end);
          }
        }
      }

      if (bestedge3->fragId() > 0) {
        if (logFileFlags & LOG_SET_PARENT_AND_HANG)
          writeLog("setParentAndHang()--  BEST3     - frag %d in unitig %d 3' to %d/%c' in unitig %d\n",
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
              writeLog("                                - -> frag %d at %d,%d 3' edge to prev frag %d at %d,%d -- hang %d,%d\n",
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
                writeLog("                                - <- frag %d at %d,%d %c' edge fr prev frag %d at %d,%d -- hang %d,%d\n",
                        oth->ident, oth->position.bgn, oth->position.end, bestedge5->frag3p() ? '3' : '5',
                        frg->ident, frg->position.bgn, frg->position.end,
                        frg->ahang, frg->bhang);
            } else {
              if (logFileFlags & LOG_SET_PARENT_AND_HANG)
                writeLog("                             - <- frag %d at %d,%d %c' edge fr prev frag %d at %d,%d -- NOT VALID\n",
                        oth->ident, oth->position.bgn, oth->position.end, bestedge5->frag3p() ? '3' : '5',
                        frg->ident, frg->position.bgn, frg->position.end);
            }

          } else {
            if (logFileFlags & LOG_SET_PARENT_AND_HANG)
              writeLog("                                - -- frag %d at %d,%d 3' edge to prev frag %d at %d,%d -- NOT VALID\n",
                      frg->ident, frg->position.bgn, frg->position.end,
                      oth->ident, oth->position.bgn, oth->position.end);
          }
        }
      }
    }  //  Over all fragment
  }  //  Over all unitigs
}

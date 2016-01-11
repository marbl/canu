
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
 *    src/AS_BAT/AS_BAT_IntersectSplit.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_Breaking.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_IntersectSplit.H"


//  The original version was filtering breakpoints.  It was accepting any break point with more than
//  MIN_BREAK_FRAGS fragments and longer than MIN_BREAK_LENGTH.  The shorter ones in between two
//  large break points were (I suspect) analyzed to see if many short break points were piling up in
//  one region.  If so, one was selected and accepted into the list of final break points.
//
//  This filtering was implemented as the first step in breakUnitigAt().  This turned
//  a (supposedly) general purpose unitig breaker into very special case.  And added a ton
//  of complexity to the UnitigBreakPoint structure -- it needed to keep all the info needed
//  for filtering.

static const int MIN_BREAK_FRAGS   = 1;
static const int MIN_BREAK_LENGTH  = 500;



intersectionList::intersectionList(UnitigVector &unitigs) {

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig             *tig = unitigs[ti];

    if (tig == NULL)
      continue;

    intersectionEvidence *evidence = new intersectionEvidence [tig->ufpath.size()];

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode  *frg = &tig->ufpath[fi];

      if (OG->isContained(frg->ident))
        continue;

      //  For my best overlap, the ID of the unitig that the overlapping fragment is in.

      evidence[fi].edge5 = *OG->getBestEdgeOverlap(frg->ident, false);
      evidence[fi].edge3 = *OG->getBestEdgeOverlap(frg->ident, true);

      evidence[fi].frag5tig = tig->fragIn(evidence[fi].edge5.fragId());
      evidence[fi].frag3tig = tig->fragIn(evidence[fi].edge3.fragId());

      //  Do NOT initialize these!  An earlier fragment could have already confirmed an end.
      //  Properly, only the 5' end of a forward fragment (or 3' end of a reverse fragment) can be
      //  confirmed already (otherwise the tig is nonsense), but we don't yet check that.
      //
      //evidence[fi].frag5confirmed = false;
      //evidence[fi].frag3confirmed = false;

      //  But, because the path could be promiscuous, not every overlap to a different tig is bad.
      //
      //  If my best overlap is to a different tig, but there is an overlapping fragment (in the
      //  unitig placement) with a best edge to me, I'm still good.  The BOG build this unitig using
      //  the edge from the other fragment to me.
      //
      //  If the fragments do not overlap in the layout (yet the best edge still exists) that is a
      //  self-intersection.
      //
      //  The two blocks are identical, except for 'edge3' and 'edge5'.

      if (evidence[fi].frag5tig == tig->id()) {
        uint32   ti  = tig->pathPosition(evidence[fi].edge5.fragId());
        ufNode  *trg = &tig->ufpath[ti];

        uint32  minf = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
        uint32  maxf = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

        uint32  mint = (trg->position.bgn < trg->position.end) ? trg->position.bgn : trg->position.end;
        uint32  maxt = (trg->position.bgn < trg->position.end) ? trg->position.end : trg->position.bgn;

        //  If they overlap, mark as confirmed, else remember an intersection.

        if (((minf < mint) && (mint < maxf)) ||  //  t begins inside f
            ((mint < minf) && (minf < maxt))) {  //  f begins inside t
          if (evidence[fi].edge5.frag3p())
            evidence[ti].frag3confirmed = true;
          else
            evidence[ti].frag5confirmed = true;

        } else {
          evidence[fi].frag5self = true;

          //  Not the correct place to report this.  Some of these get confirmed by later fragments.
          //writeLog("BUG1 F: %d,%d T %d,%d\n", minf, maxf, mint, maxt);
          //writeLog("INTERSECT from unitig %d frag %d end %d TO unitig %d frag %d end %d (SELF)\n",
          //        tig->id(), frg->ident, 5, evidence[fi].frag5tig, evidence[fi].edge5.fragId(), evidence[fi].edge5.frag3p() ? 3 : 5);
        }
      }



      if (evidence[fi].frag3tig == tig->id()) {
        uint32   ti  = tig->pathPosition(evidence[fi].edge3.fragId());
        ufNode  *trg = &tig->ufpath[ti];

        uint32  minf = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
        uint32  maxf = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

        uint32  mint = (trg->position.bgn < trg->position.end) ? trg->position.bgn : trg->position.end;
        uint32  maxt = (trg->position.bgn < trg->position.end) ? trg->position.end : trg->position.bgn;

        if (((minf < mint) && (mint < maxf)) ||  //  t begins inside f
            ((mint < minf) && (minf < maxt))) {  //  f begins inside t
          if (evidence[fi].edge3.frag3p())
            evidence[ti].frag3confirmed = true;
          else
            evidence[ti].frag5confirmed = true;

        } else {
          evidence[fi].frag3self = true;

          //  Not the correct place to report this.  Some of these get confirmed by later fragments.
          //writeLog("BUG2 F: %d,%d T %d,%d\n", minf, maxf, mint, maxt);
          //writeLog("INTERSECT from unitig %d frag %d end %d TO unitig %d frag %d end %d (SELF)\n",
          //        tig->id(), frg->ident, 3, evidence[fi].frag3tig, evidence[fi].edge3.fragId(), evidence[fi].edge3.frag3p() ? 3 : 5);
        }
      }
    }

    //
    //  Build the list.
    //

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode             *frg = &tig->ufpath[fi];

      if ((evidence[fi].frag5tig != 0) &&
          (evidence[fi].frag5tig != tig->id()) &&
          (evidence[fi].frag5confirmed == false))
        isects.push_back(intersectionPoint(evidence[fi].edge5, frg->ident, false, false));

      if ((evidence[fi].frag5tig == tig->id()) &&
          (evidence[fi].frag5self == true) &&
          (evidence[fi].frag5confirmed == false))
        isects.push_back(intersectionPoint(evidence[fi].edge5, frg->ident, false, true));

      if ((evidence[fi].frag3tig != 0) &&
          (evidence[fi].frag3tig != tig->id()) &&
          (evidence[fi].frag3confirmed == false))
        isects.push_back(intersectionPoint(evidence[fi].edge3, frg->ident, true, false));

      if ((evidence[fi].frag3tig == tig->id()) &&
          (evidence[fi].frag3self == true) &&
          (evidence[fi].frag3confirmed == false))
        isects.push_back(intersectionPoint(evidence[fi].edge3, frg->ident, true, true));
    }

    delete [] evidence;
  }


  //  Sort the intersections by the ID of the intersected fragment, then build an index into the array.

  std::sort(isects.begin(), isects.end());

  //  Terminate the intersection list with a sentinal intersection.  This is CRITICAL
  //  to the way we iterate over intersections.

  isects.push_back(intersectionPoint(BestEdgeOverlap(), 0, true, true));

  //  Build a map from fragment id to the first intersection in the list.

  for (uint32 i=0; i<isects.size(); i++) {
    isectsNum[isects[i].isectFrg]++;

    if (isectsMap.find(isects[i].isectFrg) == isectsMap.end())
      isectsMap[isects[i].isectFrg] = i;
  }
}


intersectionList::~intersectionList() {
}


void
intersectionList::logIntersections(void) {

  for (uint32 ii=0; ii<isects.size(); ii++) {
    intersectionPoint  *isect = &isects[ii];

    writeLog("INTERSECT from unitig %d frag %d end %d TO unitig %d frag %d end %d\n",
            Unitig::fragIn(isect->isectFrg), isect->isectFrg, isect->isect3p ? 3 : 5,
            Unitig::fragIn(isect->invadFrg), isect->invadFrg, isect->invad3p ? 3 : 5);
  }
}





void
breakUnitigs(UnitigVector &unitigs,
             char         *output_prefix,
             bool          enableIntersectionBreaking) {

  writeLog("==> BREAKING UNITIGS.\n");

  intersectionList  *ilist = new intersectionList(unitigs);

  //  Stop when we've seen all current unitigs.  Replace tiMax
  //  in the for loop below with unitigs.size() to recursively
  //  split unitigs.

  uint32 tiMax = unitigs.size();

  for (uint32 ti=0; ti<tiMax; ti++) {
    Unitig             *tig = unitigs[ti];

    if (tig == NULL)
      continue;

    vector<breakPoint>   breaks;

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode             *frg   = &tig->ufpath[fi];
      intersectionPoint  *isect = ilist->getIntersection(frg->ident, 0);

      if (isect == NULL)
        continue;

      for (; isect->isectFrg == frg->ident; isect++) {
        assert(tig->id() == Unitig::fragIn(isect->isectFrg));

        //  Grab the invading unitig

        Unitig *inv = unitigs[Unitig::fragIn(isect->invadFrg)];
        assert(inv->id() == Unitig::fragIn(isect->invadFrg));

        //  Grab the best edges off the invading fragment.

        BestEdgeOverlap  *best5 = OG->getBestEdgeOverlap(isect->invadFrg, false);
        BestEdgeOverlap  *best3 = OG->getBestEdgeOverlap(isect->invadFrg, true);

        //  Check if the incoming tig is a spur, and we should just ignore it immediately

        if ((inv->ufpath.size() == 1) &&
            ((best5->fragId() == 0) ||
             (best3->fragId() == 0))) {
          if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
            writeLog("unitig %d frag %d end %c' into unitig %d frag %d end %c' -- IS A SPUR, skip it\n",
                    inv->id(), isect->invadFrg, isect->invad3p ? '3' : '5',
                    tig->id(), isect->isectFrg, isect->isect3p ? '3' : '5');
          continue;
        }

        //  Keep only significant intersections

        if ((inv->getLength()   > MIN_BREAK_LENGTH) &&
            (inv->ufpath.size() > MIN_BREAK_FRAGS)) {
          if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
            writeLog("unitig %d frag %d end %c' into unitig %d frag %d end %c'\n",
                    inv->id(), isect->invadFrg, isect->invad3p ? '3' : '5',
                    tig->id(), isect->isectFrg, isect->isect3p ? '3' : '5');
          breaks.push_back(breakPoint(isect->isectFrg, isect->isect3p, true, false));
        }
      }  //  Over all incoming fragments

      //  If this is the last fragment, terminate the break point list with a 'fakeEnd' (in AS_BAT_Breaking.cc) break point
      //  at the end of the unitig.

      if ((fi+1 == tig->ufpath.size()) &&
          (breaks.size() > 0)) {
        breaks.push_back(breakPoint(frg->ident, (frg->position.bgn < frg->position.end), true, false));
      }
    }  //  Over all fragments in the unitig


    if (breaks.size() == 0)
      continue;

    //  Report where breaks occur.  'breaks' is a list, not a vector.
#if 0
    //  We've lost the fields in breaks[i] -- but the reports above aren't updated yet.
    if (logFileFlagSet(LOG_INTERSECTION_BREAKING) ||
        logFileFlagSet(LOG_MATE_SPLIT_COVERAGE_PLOT))
      for (uint32 i=0; i<breaks.size(); i++)
        writeLog("BREAK unitig %d at position %d,%d from inSize %d inFrags %d.\n",
                tig->id(),
                breaks[i].fragPos.bgn,
                breaks[i].fragPos.end,
                breaks[i].inSize,
                breaks[i].inFrags);
#endif

    //  Actually do the breaking.
    if (enableIntersectionBreaking)
      breakUnitigAt(unitigs, tig, breaks, true);

    breaks.clear();
  }  //  Over all unitigs
}

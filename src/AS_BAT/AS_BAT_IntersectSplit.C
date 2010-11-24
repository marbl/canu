
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

static const char *rcsid = "$Id: AS_BAT_IntersectSplit.C,v 1.1 2010-11-24 01:03:31 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_UnitigGraph.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "MultiAlignStore.h"



class intersectionEvidence {
public:
  intersectionEvidence() {
    frag5tig = 0;
    frag3tig = 0;

    frag5confirmed = false;
    frag3confirmed = false;

    frag5self = false;
    frag3self = false;
  };
  ~intersectionEvidence() {
  };


  BestEdgeOverlap   edge5;    //  fragID of the frag on our 5' best overlap
  BestEdgeOverlap   edge3;    //

  uint32   frag5tig;          //  tigID the frag on our 5' best overlap is in
  uint32   frag3tig;          //

  uint32   frag5confirmed:1;  //  true if our 5' end is confirmed by a best overlap in the same unitig
  uint32   frag3confirmed:1;  //

  uint32   frag5self:1;       //  true if our 5' end is intersecting the same unitig
  uint32   frag3self:1;
};



class intersectionPoint {
public:
  intersectionPoint() {
    isectFrg  = 0;
    isect3p   = false;

    sourceFrg = 0;
    source3p  = false;

    isSelf    = false;
  };
  intersectionPoint(BestEdgeOverlap edge, uint32 sId, bool s3p, bool self) {
    isectFrg  = edge.fragId();
    isect3p   = edge.frag3p();

    sourceFrg = sId;
    source3p  = s3p;

    isSelf    = self;
  };
  ~intersectionPoint() {
  };

  bool operator<(const intersectionPoint that) const {
    return(isectFrg < that.isectFrg);
  };

  uint32   isectFrg;  //  Fragment that is intersected into, we split on this.
  bool     isect3p;   //  True if we intersected onto the 3' end of the fragment.

  uint32   sourceFrg;
  bool     source3p;

  bool     isSelf;
};



static
void
writeBreakOVL(FILE             *breakFile,
              uint32            aFrag,
              uint32            bFrag,
              bool              best3p,
              BestEdgeOverlap  *bestEdge) {
  GenericMesg  pmesg;
  OverlapMesg  omesg;

  if (breakFile == NULL)
    return;

  omesg.aifrag          = aFrag;  //inFrag;
  omesg.bifrag          = bFrag;  ///f->ident;
  omesg.ahg             = bestEdge->ahang();
  omesg.bhg             = bestEdge->bhang();
  omesg.orientation.setIsUnknown();
  omesg.overlap_type    = AS_DOVETAIL;
  omesg.quality         = 0.0;
  omesg.min_offset      = 0;
  omesg.max_offset      = 0;
  omesg.polymorph_ct    = 0;
  omesg.alignment_trace = NULL;
#ifdef AS_MSG_USE_OVL_DELTA
  omesg.alignment_delta = NULL;
#endif

  if ((best3p == false) && (bestEdge->frag3p() == false))
    omesg.orientation.setIsOuttie();
  if ((best3p == false) && (bestEdge->frag3p() == true))
    omesg.orientation.setIsAnti();
  if ((best3p == true) && (bestEdge->frag3p() == false))
    omesg.orientation.setIsNormal();
  if ((best3p == true) && (bestEdge->frag3p() == true))
    omesg.orientation.setIsInnie();

  pmesg.t = MESG_OVL;
  pmesg.m = &omesg;

  WriteProtoMesg_AS(breakFile, &pmesg);
}



void
UnitigGraph::breakUnitigs(char *output_prefix,
                          bool enableIntersectionBreaking) {

  fprintf(logFile, "==> BREAKING UNITIGS.\n");

  //
  //  Analyze tigs for intersections
  //

  vector<intersectionPoint>   isects;

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
          //fprintf(logFile, "BUG1 F: %d,%d T %d,%d\n", minf, maxf, mint, maxt);
          //fprintf(logFile, "INTERSECT from unitig %d frag %d end %d TO unitig %d frag %d end %d (SELF)\n",
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
          //fprintf(logFile, "BUG2 F: %d,%d T %d,%d\n", minf, maxf, mint, maxt);
          //fprintf(logFile, "INTERSECT from unitig %d frag %d end %d TO unitig %d frag %d end %d (SELF)\n",
          //        tig->id(), frg->ident, 3, evidence[fi].frag3tig, evidence[fi].edge3.fragId(), evidence[fi].edge3.frag3p() ? 3 : 5);
        }
      }
    }

    //
    //  Dump.
    //

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode             *frg = &tig->ufpath[fi];

      if ((evidence[fi].frag5tig != 0) &&
          (evidence[fi].frag5tig != tig->id()) &&
          (evidence[fi].frag5confirmed == false)) {
        if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
          fprintf(logFile, "INTERSECT from unitig %d frag %d end %d TO unitig %d frag %d end %d\n",
                  tig->id(), frg->ident, 5, evidence[fi].frag5tig, evidence[fi].edge5.fragId(), evidence[fi].edge5.frag3p() ? 3 : 5);
        isects.push_back(intersectionPoint(evidence[fi].edge5, frg->ident, false, false));
      }

      if ((evidence[fi].frag5tig == tig->id()) &&
          (evidence[fi].frag5self == true) &&
          (evidence[fi].frag5confirmed == false)) {
        if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
          fprintf(logFile, "INTERSECT from unitig %d frag %d end %d TO unitig %d frag %d end %d (SELF)\n",
                  tig->id(), frg->ident, 5, evidence[fi].frag5tig, evidence[fi].edge5.fragId(), evidence[fi].edge5.frag3p() ? 3 : 5);
        isects.push_back(intersectionPoint(evidence[fi].edge5, frg->ident, false, true));
      }

      if ((evidence[fi].frag3tig != 0) &&
          (evidence[fi].frag3tig != tig->id()) &&
          (evidence[fi].frag3confirmed == false)) {
        if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
          fprintf(logFile, "INTERSECT from unitig %d frag %d end %d TO unitig %d frag %d end %d\n",
                  tig->id(), frg->ident, 3, evidence[fi].frag3tig, evidence[fi].edge3.fragId(), evidence[fi].edge3.frag3p() ? 3 : 5);
        isects.push_back(intersectionPoint(evidence[fi].edge3, frg->ident, true, false));
      }

      if ((evidence[fi].frag3tig == tig->id()) &&
          (evidence[fi].frag3self == true) &&
          (evidence[fi].frag3confirmed == false)) {
        if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
          fprintf(logFile, "INTERSECT from unitig %d frag %d end %d TO unitig %d frag %d end %d (SELF)\n",
                  tig->id(), frg->ident, 3, evidence[fi].frag3tig, evidence[fi].edge3.fragId(), evidence[fi].edge3.frag3p() ? 3 : 5);
        isects.push_back(intersectionPoint(evidence[fi].edge3, frg->ident, true, true));
      }
    }

    delete [] evidence;
  }


  //  Sort the intersections by the ID of the intersected fragment, then build an index into the array.

  std::sort(isects.begin(), isects.end());

  map<uint32,uint32>  isectsMap;

  for (uint32 i=0; i<isects.size(); i++)
    if (isectsMap.find(isects[i].isectFrg) == isectsMap.end())
      isectsMap[isects[i].isectFrg] = i;

  //  Stop when we've seen all current unitigs.  Replace tiMax
  //  in the for loop below with unitigs.size() to recursively
  //  split unitigs.
  uint32 tiMax = unitigs.size();

  for (uint32 ti=0; ti<tiMax; ti++) {
    Unitig             *tig = unitigs[ti];

    if (tig == NULL)
      continue;

    UnitigBreakPoints   breaks;

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode  *frg = &tig->ufpath[fi];

      //  Find the first piece

      if (isectsMap.find(frg->ident) == isectsMap.end())
        continue;

      uint32 ii = isectsMap[frg->ident];
      assert(isects[ii].isectFrg == frg->ident);
      assert(tig->id() == Unitig::fragIn(isects[ii].isectFrg));

      for (; isects[ii].isectFrg == frg->ident; ii++) {

        //  Grab the invading unitig

        Unitig *inv = unitigs[Unitig::fragIn(isects[ii].sourceFrg)];
        assert(inv->id() == Unitig::fragIn(isects[ii].sourceFrg));

        //  Grab the best edges off the invading fragment.

        BestEdgeOverlap  *best5 = OG->getBestEdgeOverlap(isects[ii].sourceFrg, false);
        BestEdgeOverlap  *best3 = OG->getBestEdgeOverlap(isects[ii].sourceFrg, true);

        //  Check if the incoming tig is a spur, and we should just ignore it immediately

        if ((inv->ufpath.size() == 1) &&
            ((best5->fragId() == 0) ||
             (best3->fragId() == 0))) {
          if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
            fprintf(logFile, "unitig %d frag %d end %c' into unitig %d frag %d end %c' -- IS A SPUR, skip it\n",
                    inv->id(), isects[ii].sourceFrg, isects[ii].source3p ? '3' : '5',
                    tig->id(), isects[ii].isectFrg,  isects[ii].isect3p  ? '3' : '5');
          continue;
        }

        //  
            
        if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
          fprintf(logFile, "unitig %d frag %d end %c' into unitig %d frag %d end %c'\n",
                  inv->id(), isects[ii].sourceFrg, isects[ii].source3p ? '3' : '5',
                  tig->id(), isects[ii].isectFrg,  isects[ii].isect3p  ? '3' : '5');

        breaks.push_back(UnitigBreakPoint(isects[ii].isectFrg, isects[ii].isect3p, frg->position,
                                          fi, tig->ufpath.size() - fi - 1,
                                          isects[ii].sourceFrg, isects[ii].source3p,
                                          inv->getLength(),
                                          inv->ufpath.size()));
      }  //  Over all incoming fragments

      //  If this is the last fragment, terminate the break point list with a 'fakeEnd' (in AS_BAT_Breaking.cc) break point
      //  at the end of the unitig.

      if ((fi+1 == tig->ufpath.size()) &&
          (breaks.size() > 0)) {
        breaks.push_back(UnitigBreakPoint(frg->ident, (frg->position.bgn < frg->position.end), frg->position,
                                          fi, 0,
                                          0, false,
                                          std::numeric_limits<int>::max(),
                                          0));
      }
    }  //  Over all fragments in the unitig

    //  If there are break points in the list, filter and break.

    if (breaks.size() > 0) {
      filterBreakPoints(tig, breaks);

      //  Report where breaks occur.  'breaks' is a list, not a vector.
      if (logFileFlagSet(LOG_INTERSECTION_BREAKING) ||
          logFileFlagSet(LOG_MATE_SPLIT_COVERAGE_PLOT))
        for (uint32 i=0; i<breaks.size(); i++)
          fprintf(logFile, "BREAK unitig %d at position %d,%d from inSize %d inFrags %d.\n",
                  tig->id(),
                  breaks[i].fragPos.bgn,
                  breaks[i].fragPos.end,
                  breaks[i].inSize,
                  breaks[i].inFrags);

      //  Actually do the breaking.
      if (enableIntersectionBreaking) {
        UnitigVector* newUs = breakUnitigAt(tig, breaks);

        if (newUs != NULL) {
          delete tig;
          unitigs[ti] = NULL;
          unitigs.insert(unitigs.end(), newUs->begin(), newUs->end());
        }

        delete newUs;
      }

      breaks.clear();
    }
  }  //  Over all unitigs
}

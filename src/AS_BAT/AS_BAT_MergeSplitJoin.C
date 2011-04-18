
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

static const char *rcsid = "$Id: AS_BAT_MergeSplitJoin.C,v 1.6 2011-04-18 02:00:20 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_IntersectSplit.H"

#include "AS_BAT_MergeSplitJoin.H"

#include "AS_BAT_Breaking.H"

#include "AS_UTL_intervalList.H"


//  Over every frag in 'target' see if the unitig from any invading fragment can
//  be merged in.

#define UNCOVERED_NOISE_FILTER 20
#define COVERED_BORDER         10

#define SPURIOUS_REPEAT_THRESHOLD 5

#define STUPID_REPEAT_SIZE 100

//  Define this to waste effort and return all placements for fragments.  This will just fail
//  asserts here, but they're easy to rearrange and avoid.
//
//  Note that this returns different data, due to an unstable sort.  In particular, the
//  overlapPlacement_byCluster sort only cares to group placements by cluster, the ordering within
//  the cluster is not specified.
//
#define RETURN_ALL_PLACEMENTS 0

class repeatUniqueJunction {
public:
  repeatUniqueJunction() {
    bgn      = 0;
    end      = 0;
    point    = 0;
    overlaps = 0;
    rptLeft  = false;
  };

  repeatUniqueJunction(uint32 bgn_, uint32 end_, uint32 point_, uint32 overlaps_, bool rptLeft_) {
    bgn      = bgn_;
    end      = end_;
    point    = point_;
    overlaps = overlaps_;
    rptLeft  = rptLeft_;
  };

  bool operator<(repeatUniqueJunction const that) const {
    return(point < that.point);
  };

  uint32  bgn;
  uint32  end;
  uint32  point;
  uint32  overlaps;
  bool    rptLeft;  //  Repeat is to the left of the point
};




bool
mergeBubbles_findEnds(UnitigVector &unitigs,
                      Unitig *bubble,
                      ufNode &fFrg,
                      ufNode &lFrg,
                      Unitig *target) {

  //  Search for edges.  For a bubble to exist, either the first or last non-contained fragment
  //  must have an edge to the 'merge' unitig it is a bubble of.  Ideally, both the first and
  //  last will have edges to the same unitig, but we'll test and allow only a single edge.

  uint32  zIdx = ~(uint32)0;
  uint32  fIdx = zIdx;
  uint32  lIdx = zIdx;

  //  We'd like to claim that all unitigs begin with a non-contained fragment, but zombie fragments
  //  (contained fragments that are in a circular containment relationship) violate this.  So, we
  //  could then claim that unitigs with more than one fragment begin with a non-contained fragment.
  //  But any zombie that has a bubble popped into it violate this.
  //
  //  Since nothing is really using the best edges from non-contained fragmnet, in cases where the
  //  first fragment is contained, we'll use the first and last fragments, instead of the first and
  //  last non-contained fragments.

  //  Search for the first and last non-contained fragments.

  for (uint32 ii=0; ((fIdx == zIdx) && (ii < bubble->ufpath.size())); ii++)
    if (OG->isContained(bubble->ufpath[ii].ident) == false)
      fIdx = ii;

  for (uint32 ii=bubble->ufpath.size(); ((lIdx == zIdx) && (ii-- > 0)); )
    if (OG->isContained(bubble->ufpath[ii].ident) == false)
      lIdx = ii;

  //  Didn't find a non-contained fragment!  Reset to the first/last fragments.

  if (fIdx == zIdx) {
    fprintf(logFile, "popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments STARTS WITH A CONTAINED FRAGMENT %d\n",
            bubble->id(), bubble->getLength(), bubble->ufpath.size(),
            bubble->ufpath[0].ident);
    fIdx = 0;
    lIdx = bubble->ufpath.size() - 1;
  }

  assert(fIdx <= lIdx);
  assert(fIdx != zIdx);
  assert(lIdx != zIdx);

  fFrg = bubble->ufpath[fIdx];  //  NOTE:  A COPY, not a pointer or reference.
  lFrg = bubble->ufpath[lIdx];  //         These get modified.

  //  Grab the best edges outside the unitig.  If the first fragment is reversed, we want
  //  to grab the edge off of the 3' end; opposite for the last fragment.
  //
  //  There is ALWAYS a best edge, even for contained fragments.  The edge might be empty, with fragId == 0.

  bool             f3p   = (isReverse(fFrg.position) == true);
  BestEdgeOverlap *fEdge = OG->getBestEdgeOverlap(fFrg.ident, f3p);
  uint32           fUtg  = Unitig::fragIn(fEdge->fragId());

  bool             l3p   = (isReverse(lFrg.position) == false);
  BestEdgeOverlap *lEdge = OG->getBestEdgeOverlap(lFrg.ident, l3p);
  uint32           lUtg  = Unitig::fragIn(lEdge->fragId());

  //  Just make sure...those edges should NOT to be to ourself.  But if they are, we'll just ignore
  //  them -- these can be from circular unitigs (in which case we can't really merge ourself into
  //  ourself at the correct spot) OR from a bubble that was already merged and just happened to tie
  //  for a fragment at the end.
  //
  //  aaaaaaaaaa
  //      aaaaaaaaaaa
  //    bbbbbbbbb
  //       bbbbbbbbbb
  //
  //  The 'b' unitig was merged into a.  The second b fragment now becomes the last non-contained
  //  fragment in the merged unitig, and it has a best edge back to the original 'a' fragment,
  //  ourself.

  if (fUtg == bubble->id()) {
    fEdge = NULL;
    fUtg  = 0;
  }
  if (lUtg == bubble->id()) {
    lEdge = NULL;
    lUtg  = 0;
  }

#if 0
  if ((fUtg != 0) && (lUtg != 0))
    fprintf(logFile, "popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  Edges (%d/%d') from frag %d/%d' and (%d/%d') from frag %d/%d'\n",
            bubble->id(), bubble->getLength(), bubble->ufpath.size(),
            fEdge->fragId(), (fEdge->frag3p() ? 3 : 5), fFrg.ident, (f3p ? 3 : 5),
            lEdge->fragId(), (lEdge->frag3p() ? 3 : 5), lFrg.ident, (l3p ? 3 : 5));
  else if (fUtg != 0)
    fprintf(logFile, "popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  Edge (%d/%d') from frag %d/%d'\n",
            bubble->id(), bubble->getLength(), bubble->ufpath.size(),
            fEdge->fragId(), (fEdge->frag3p() ? 3 : 5), fFrg.ident, (f3p ? 3 : 5));
  else if (lUtg != 0)
    fprintf(logFile, "popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  Edge (%d/%d') from frag %d/%d'\n",
            bubble->id(), bubble->getLength(), bubble->ufpath.size(),
            lEdge->fragId(), (lEdge->frag3p() ? 3 : 5), lFrg.ident, (l3p ? 3 : 5));
  else
    //  But then how do we get an intersection?!?!!  Intersections from a bubble that was
    //  already popped.  We pop A into B, and while iterating through fragments in B we find
    //  the -- now obsolete -- intersections we originally used and try to pop it again.
    //
    fprintf(logFile, "popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  NO EDGES, no bubble.\n",
            bubble->id(), bubble->getLength(), bubble->ufpath.size());
#endif

  if ((fUtg == 0) && (lUtg == 0))
    return(false);

#if 0
  //  The only interesting case here is if we have both edges and they point to different unitigs.
  //  We might as well place it aggressively.

  if ((fUtg != 0) && (lUtg != 0) && (fUtg != lUtg)) {
    fprintf(logFile, "popBubbles()--   bubble unitig %d has edges to both unitig %d and unitig %d, cannot place (yet)\n",
            bubble->id(), fUtg, lUtg);
    return(true);
  }
#endif

  return(true);
}


bool
mergeBubbles_checkEnds(UnitigVector &unitigs,
                       Unitig *bubble,
                       ufNode &fFrg,
                       ufNode &lFrg,
                       Unitig *target) {

  //  Compute placement of the two fragments.  Compare the size against the bubble.

  ufNode fFrgN = fFrg;
  ufNode lFrgN = lFrg;

  overlapPlacement    fFrgPlacement;
  overlapPlacement    lFrgPlacement;

  fFrgPlacement.errors  = 4.0e9;
  fFrgPlacement.aligned = 1;

  lFrgPlacement.errors  = 4.0e9;
  lFrgPlacement.aligned = 1;

  vector<overlapPlacement>   placements;

  placements.clear();

  placeFragUsingOverlaps(unitigs, (RETURN_ALL_PLACEMENTS) ? NULL : target, fFrg.ident, placements);
  for (uint32 i=0; i<placements.size(); i++) {
    assert(placements[i].tigID == target->id());
    if (placements[i].tigID != target->id()) continue;

    if (placements[i].fCoverage < 0.99)
      continue;

    if (placements[i].errors / placements[i].aligned < fFrgPlacement.errors / fFrgPlacement.aligned)
      fFrgPlacement = placements[i];
  }

  fFrgN.ident             = fFrgPlacement.frgID;
  fFrgN.contained         = 0;
  fFrgN.parent            = 0;
  fFrgN.ahang             = 0;
  fFrgN.bhang             = 0;
  fFrgN.position          = fFrgPlacement.position;
  fFrgN.containment_depth = 0;

  if ((fFrgN.position.bgn == 0) &&
      (fFrgN.position.end == 0)) {
    //fprintf(logFile, "popBubbles()--   failed to place fFrg.\n");
    return(false);
  }


  placements.clear();

  placeFragUsingOverlaps(unitigs, (RETURN_ALL_PLACEMENTS) ? NULL : target, lFrg.ident, placements);
  for (uint32 i=0; i<placements.size(); i++) {
    assert(placements[i].tigID == target->id());
    if (placements[i].tigID != target->id()) continue;

    if (placements[i].fCoverage < 0.99)
      continue;

    if (placements[i].errors / placements[i].aligned < lFrgPlacement.errors / lFrgPlacement.aligned)
      lFrgPlacement = placements[i];
  }

  lFrgN.ident             = lFrgPlacement.frgID;
  lFrgN.contained         = 0;
  lFrgN.parent            = 0;
  lFrgN.ahang             = 0;
  lFrgN.bhang             = 0;
  lFrgN.position          = lFrgPlacement.position;
  lFrgN.containment_depth = 0;

  if ((lFrgN.position.bgn == 0) &&
      (lFrgN.position.end == 0)) {
    //fprintf(logFile, "popBubbles()--   failed to place lFrg.\n");
    return(false);
  }


  int32 minL = MIN(fFrg.position.bgn, fFrg.position.end);
  int32 maxL = MAX(fFrg.position.bgn, fFrg.position.end);

  int32 minR = MIN(lFrg.position.bgn, lFrg.position.end);
  int32 maxR = MAX(lFrg.position.bgn, lFrg.position.end);

  int32 placedLen = MAX(maxL, maxR) - MIN(minL, minR);

  if (2 * placedLen < bubble->getLength()) {
    //  Too short.
    fprintf(logFile, "popBubbles()--   too short.  fFrg %d,%d lFrg %d,%d.  L %d,%d R %d,%d len %d\n",
            fFrg.position.bgn, fFrg.position.end,
            lFrg.position.bgn, lFrg.position.end,
            minL, maxL, minR, maxR, placedLen);
    return(false);
  }

  if (2 * bubble->getLength() < placedLen) {
    //  Too long.
    fprintf(logFile, "popBubbles()--   too long.  fFrg %d,%d lFrg %d,%d.  L %d,%d R %d,%d len %d\n",
            fFrg.position.bgn, fFrg.position.end,
            lFrg.position.bgn, lFrg.position.end,
            minL, maxL, minR, maxR, placedLen);
    return(false);
  }

  ////////////////////
  //
  //  Check orientations
  //
  ////////////////////

  //  If fFrg and lFrg are the same fragment (bubble is one uncontained fragment) then we're done.

  if (fFrg.ident == lFrg.ident) {
    fFrg = fFrgN;
    lFrg = lFrgN;
    return(true);
  }

  //  Otherwise, check that the orientation and positioning of the before and after fragments is the
  //  same.

  bool   bL    = (isReverse(fFrg.position));
  bool   bR    = (isReverse(lFrg.position));
  bool   bOrd  = (MIN(fFrg.position.bgn, fFrg.position.end) < MIN(lFrg.position.bgn, lFrg.position.end));

  bool   nL    = (isReverse(fFrgN.position));
  bool   nR    = (isReverse(lFrgN.position));
  bool   nOrd  = (MIN(fFrgN.position.bgn, fFrgN.position.end) < MIN(lFrgN.position.bgn, lFrgN.position.end));

  if (((bL == nL) && (bR == nR) && (bOrd == nOrd)) ||
      ((bL != nL) && (bR != nR) && (bOrd != nOrd))) {
    //  Yup, looks good!
    fFrg = fFrgN;
    lFrg = lFrgN;
    return(true);
  }

  //  Nope, something got screwed up in alignment.

#if 0
  fprintf(logFile, "popBubbles()--   Order/Orientation problem.  bL %d bR %d bOrd %d  nL %d nR %d nOrd %d\n",
          bL, bR, bOrd,
          nL, nR, nOrd);
#endif
  return(false);
}


//  False if any of the fragments in 'bubble' are not fully covered by overlaps to fragments in
//  'target'.  Such uncovered fragments would indicate a large bubble -- large enough that we failed
//  to find an overlap -- and would cause problems in consensus.
//
//  False if any of the fragments in 'bubble' cannot be placed between fFrg and lFrg.  This would
//  indicate the bubble contains a significant rearrangement and would cause problems in consensus.
//
//  If the above tests pass, 'bubble' is inserted into 'target' and 'bubble' is deleted.
//
static
bool
mergeBubbles_checkFrags(UnitigVector &unitigs,
                        Unitig *bubble,
                        ufNode &fFrg,
                        ufNode &lFrg,
                        Unitig *target) {

  //  Method:
  //
  //  * Call placeFragUsingOverlaps() for every fragment.  Save the placements returned.
  //  * Count the number of placements that are outside the fFrg/lFrg range.
  //  * Isolate down to one 'best' placement for each fragment.
  //    * Must be within fFrg/lFrg.
  //    * Resolve ties with
  //      * Placement in the original unitig       
  //      * Error rates on overlaps
  //      * Mate pairs

  bool success = false;

  vector<overlapPlacement>    *placements   = new vector<overlapPlacement> [bubble->ufpath.size()];
  overlapPlacement            *correctPlace = new        overlapPlacement  [bubble->ufpath.size()];

  for (uint32 fi=0; fi<bubble->ufpath.size(); fi++) {
    ufNode *frg = &bubble->ufpath[fi];

    placeFragUsingOverlaps(unitigs, (RETURN_ALL_PLACEMENTS) ? NULL : target, frg->ident, placements[fi]);

    //  Initialize the final placement to be bad, so we can pick the best.
    correctPlace[fi].fCoverage = 0.0;
    correctPlace[fi].errors    = 4.0e9;
    correctPlace[fi].aligned   = 1;
  }

  //  Some bizarre cases -- possibly only from bad data -- confound any logical attempt at finding the min/max extents.  Yes, even though
  //  this should work, it doesn't.  Or maybe it's just broken and I haven't seen how.
  //
  //int32  minE = (fFrg.position.bgn < lFrg.position.bgn) ? MIN(fFrg.position.bgn, fFrg.position.end) : MIN(lFrg.position.bgn, lFrg.position.end);
  //int32  maxE = (fFrg.position.bgn < lFrg.position.bgn) ? MAX(lFrg.position.bgn, lFrg.position.end) : MAX(fFrg.position.bgn, fFrg.position.end);
  //
  //  The one case that breaks it is a bubble unitig with a single chimeric fragment.
  //    fFrg ident = 367563, contained = 0, parent = 254673, ahang =  144,   bhang = 24,  bgn = 33406, end = 33238}
  //    lFrg ident = 367563, contained = 0, parent = 147697, ahang = -58,  bhang = -157,  bgn = 33406, end = 33574}
  //
  int32  minE = MIN(MIN(fFrg.position.bgn, fFrg.position.end), MIN(lFrg.position.bgn, lFrg.position.end));
  int32  maxE = MAX(MAX(fFrg.position.bgn, fFrg.position.end), MAX(lFrg.position.bgn, lFrg.position.end));
  int32  diff = maxE - minE;

  assert(minE < maxE);

  minE -= diff / 2;    if (minE < 0)  minE = 0;
  maxE += diff / 2;

  uint32  nCorrect = 0;

  for (uint32 fi=0; fi<bubble->ufpath.size(); fi++) {
    uint32  nNotPlaced = 0;
    uint32  nNotPlacedInCorrectPosition = 0;
    uint32  nNotPlacedFully = 0;
    uint32  nNotOriented = 0;

    if (placements[fi].size() == 0)
      nNotPlaced++;

    //  If we're contained, and our container is actually in the bubble, we can (or should be able
    //  to) safely allow almost any placement.

    bool requireFullAlignment = true;

    if ((OG->isContained(bubble->ufpath[fi].ident) == true) &&
        (bubble->fragIn(OG->getBestContainer(bubble->ufpath[fi].ident)->container) == bubble->id()))
      requireFullAlignment = false;

    for (uint32 pl=0; pl<placements[fi].size(); pl++) {
      assert(placements[fi][pl].tigID == target->id());
      if (placements[fi][pl].tigID != target->id()) continue;

      int32  minP = MIN(placements[fi][pl].position.bgn, placements[fi][pl].position.end);
      int32  maxP = MAX(placements[fi][pl].position.bgn, placements[fi][pl].position.end);

      if ((maxP < minE) || (maxE < minP)) {
        nNotPlacedInCorrectPosition++;
        continue;
      }

      if ((requireFullAlignment == true) && (placements[fi][pl].fCoverage < 0.99)) {
        nNotPlacedFully++;
        continue;
      }

      //if ((placements[fi][pl].nForward > 0) &&
      //    (placements[fi][pl].nReverse > 0)) {
      //  nNotOriented++;
      //  continue;
      //}

      //  The current placement seems like a good one.  Should we keep it?

      //  The length requirement was added to solve a problem during testing on hydra.  We tried to
      //  place a contained fragment -- so skipped the fCoverage test above.  This fragment had two
      //  placements in the correct location on the target unitig.  One plaement was fCoverage=1.00,
      //  the other was fCoverage=0.15.  Clearly the first was better, but the second had less
      //  error.  Without the length filter, we'd incorrectly pick the second placement.

      bool  keepIt = false;

      if (placements[fi][pl].fCoverage > correctPlace[fi].fCoverage)
        //  Yes!  The current placement has more coverage than the saved one.
        keepIt = true;

      if ((placements[fi][pl].fCoverage >= correctPlace[fi].fCoverage) &&
          (placements[fi][pl].errors / placements[fi][pl].aligned < correctPlace[fi].errors / correctPlace[fi].aligned))
        //  Yes!  The current placement is just as long, and lower error.
        keepIt = true;

      //  Yup, looks like a better placement.
      if (keepIt)
        correctPlace[fi] = placements[fi][pl];
    }  //  over all placements

    if (correctPlace[fi].fCoverage > 0) {
      nCorrect++;
    } else {
      //  We currently require ALL fragments to be well placed, so we can abort on the first fragment that
      //  fails.
      //fprintf(logFile, "popBubbles()--   Failed to place frag %d notPlaced %d notPlacedInCorrectPosition %d notPlacedFully %d notOriented %d\n",
      //        bubble->ufpath[fi].ident, nNotPlaced, nNotPlacedInCorrectPosition, nNotPlacedFully, nNotOriented);
      break;
    }
  }

  if (nCorrect != bubble->ufpath.size())
    goto finished;

  //  Now just move the fragments into the target unitig and delete the bubble unitig.
  //
  //  Explicitly DO NOT propagate the contained, parent, ahang or bhang from the bubble here.  We
  //  could figure all this stuff out, but it definitely is NOT just a simple copy from the bubble
  //  unitig (for example, we could add the bubble unitig reversed).
  //  
  //
  for (uint32 fi=0; fi<bubble->ufpath.size(); fi++) {
    ufNode  nFrg;

    nFrg.ident             = correctPlace[fi].frgID;
    nFrg.contained         = 0;
    nFrg.parent            = 0;
    nFrg.ahang             = 0;
    nFrg.bhang             = 0;
    nFrg.position          = correctPlace[fi].position;
    nFrg.containment_depth = 0;

    //target->addFrag(nFrg, 0, logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG));
    target->addFrag(nFrg, 0, false);
  }

  target->sort();

  success = true;

  fprintf(logFile, "popBubbles()--   merged bubble unitig %d with %ld frags into unitig %d now with %ld frags\n",
          bubble->id(), bubble->ufpath.size(), target->id(), target->ufpath.size());

 finished:
  delete [] placements;
  delete [] correctPlace;

  return(success);
}



void
mergeBubbles(UnitigVector &unitigs, Unitig *target, intersectionList *ilist) {

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode             *frg   = &target->ufpath[fi];
    intersectionPoint  *isect = ilist->getIntersection(frg->ident, 0);

    if (isect == NULL)
      continue;

    for (; isect->isectFrg == frg->ident; isect++) {
      assert(target->id() == Unitig::fragIn(isect->isectFrg));

      //  Grab the potential bubble unitig

      Unitig *bubble = unitigs[Unitig::fragIn(isect->invadFrg)];
      assert(bubble->id() == Unitig::fragIn(isect->invadFrg));

      if (bubble->id() == target->id())
        //  HEY!  We're not a bubble in ourself!
        continue;

      ufNode  fFrg;  //  First fragment in the bubble
      ufNode  lFrg;  //  Last fragment in the bubble

      if (mergeBubbles_findEnds(unitigs, bubble, fFrg, lFrg, target) == false)
        continue;

      if (mergeBubbles_checkEnds(unitigs, bubble, fFrg, lFrg, target) == false)
        continue;

      if (mergeBubbles_checkFrags(unitigs, bubble, fFrg, lFrg, target) == false)
        continue;

      //  Merged!
      //  o Delete the unitig we just merged in.
      //  o Skip the rest of the intersections for this fragment (because....)
      //  o Reset iteration over fragments in this unitig -- we'll try some failed merges
      //    again, but we might pick up a bunch more from the fragments we just added.  Plus,
      //    we changed the ufpath vector.

      unitigs[bubble->id()] = NULL;
      delete bubble;

      fi = 0;

      break;
    }
  }
}



void
stealBubbles(UnitigVector &unitigs, Unitig *target, intersectionList *ilist) {
}




double
computeAverageCoverage(Unitig *target, int32 bgn, int32 end) {
  double sum = 0;

  if (bgn < 0)                    bgn = 0;
  if (end > target->getLength())  end = target->getLength();

  assert(bgn < end);

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode     *frg   = &target->ufpath[fi];

    if (frg->position.end <= bgn)
      continue;
    if (frg->position.bgn >= end)
      continue;

    //  Intersecting.

    if (frg->position.bgn < frg->position.end)
      sum += MIN(end, frg->position.end) - MAX(bgn, frg->position.bgn);
    else
      sum += MIN(end, frg->position.bgn) - MAX(bgn, frg->position.end);
  }

  return(sum / (end - bgn));
}




//  Build a list of all the fragments that have overlaps to some fragment in this unitig.
//  Exclude from the list fragments that are already in this unitig.  We expect these fragments
//  to have multiple overlaps to the unitig, and we want to examine each one only once.
//
//  Use placeFragUsingOverlaps() to place each of these fragments.  Ideally, we'd restrict to just
//  this unitig, but for now we need to filter the results.  Three results:
//
//    o Fragment is fully contained in this unitig.  Why isn't it a bubble?  It could be contained
//      completely in a repeat, and the repeat is in two different unitigs.
//  
//    o Fragment is partially aligned.  This could be indicating a repeat or chimera that we should
//      be splitting.  We save the location of the break, and direction of the unaligned piece.
//  
//  After all fragments are 'placed' the list of breaks is examined.
//
//    o A chimer will induce about as many breaks as the local depth.  Breaks will be on both
//      sides of the chimeric point.
//
//    o A repeat will have many breaks on one side, and few to none on the other side.
//
//    o A spur will look like a repeat, but at the end of a unitig, with few to no fragments
//      following.
//
void
markRepeats(UnitigVector &unitigs,
            Unitig *target) {
  set<AS_IID>               ovlFrags;

  //            --------unitig----------
  //               ----         ----
  //              /         hangEnd \         (this is not a multiline comment!)
  //             / hangBgn           \        (neither is this - gcc complains about trailing \)
  //
  //intervalList              hang5;   //  Fragments hanging towards the bgn of the unitig.
  //intervalList              hang3;   //  Fragments hanging towards the end of the unitig.
  intervalList              align5;  //  Fragments with partial overlaps.
  intervalList              align3;  //  Fragments with partial overlaps.

  vector<overlapPlacement>  places;
  vector<breakPoint>        breaks;

  //  Build a list of all the fragments that have overlaps to this unitig.

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode     *frg   = &target->ufpath[fi];
    uint32      ovlLen = 0;
    BAToverlap *ovl    = OC->getOverlaps(frg->ident, ovlLen);

    for (uint32 i=0; i<ovlLen; i++) {
      if (Unitig::fragIn(ovl[i].b_iid) != target->id())
        ovlFrags.insert(ovl[i].b_iid);
    }
  }

  //  For each overlapping fragment, place it and process.

  for (set<AS_IID>::iterator it=ovlFrags.begin(); it!=ovlFrags.end(); it++) {
    AS_IID  iid = *it;

    vector<overlapPlacement>  op;

    placeFragUsingOverlaps(unitigs, (RETURN_ALL_PLACEMENTS) ? NULL : target, iid, op);

    //  placeFragUsingOverlaps() returns the expected placement for this fragment in 'position', and
    //  the amount of the fragment covered by evidence in 'covered'.
    //
    //  Below we'll try to decipher this into two intervals, covered by evidence and not covered by
    //  evidence.

    for (uint32 pl=0; pl<op.size(); pl++) {
      assert(op[pl].tigID == target->id());
      if (op[pl].tigID != target->id()) continue;

      fprintf(logFile, "markRepeats()-- op[%3d] tig %d fCoverage %f position %d %d verified %d %d\n",
              pl,
              op[pl].tigID,
              op[pl].fCoverage,
              op[pl].position.bgn, op[pl].position.end,
              op[pl].verified.bgn, op[pl].verified.end);

      if (op[pl].fCoverage > 0.99)
        //  No worries, fully placed.
        continue;

      if ((op[pl].position.bgn < 0) ||
          (op[pl].position.bgn > target->getLength()) ||
          (op[pl].position.end < 0) ||
          (op[pl].position.end > target->getLength()))
        //  Placed outside the range of the unitig
        continue;

      //  Otherwise, placed in the target unitig, at less than perfect coverage.  Compute
      //  the unitig coordinates that are covered by actual overlaps.

      int32   uncovered5bgn = 0;
      int32   uncovered5end = 0;

      int32   coveredbgn   = 0;
      int32   coveredend   = 0;

      int32   uncovered3bgn = 0;
      int32   uncovered3end = 0;

      FragmentEnd  *end3;
      FragmentEnd  *end5;

      if (op[pl].position.bgn < op[pl].position.end) {
        //  Fragment is placed forward.
        uncovered5bgn = op[pl].position.bgn;
        uncovered5end = op[pl].verified.bgn;

        coveredbgn   = op[pl].verified.bgn;
        coveredend   = op[pl].verified.end;

        uncovered3bgn = op[pl].verified.end;
        uncovered3end = op[pl].position.end;

        end5 = &op[pl].frag5p;
        end3 = &op[pl].frag3p;

      } else {
        //  Fragment is placed reverse.
        uncovered5bgn = op[pl].position.end;
        uncovered5end = op[pl].verified.end;

        coveredbgn   = op[pl].verified.end;
        coveredend   = op[pl].verified.bgn;

        uncovered3bgn = op[pl].verified.bgn;
        uncovered3end = op[pl].position.bgn;

        end5 = &op[pl].frag3p;
        end3 = &op[pl].frag5p;
      }

      if (uncovered5bgn < 0) {
        uncovered5bgn = 0;   //  Ignore ends that extend past the edge of the unitig
        uncovered5end = 0;
      }

      if (target->getLength() <= uncovered3end) {
        uncovered3bgn = 0;
        uncovered3end = 0;
      }

      bool  save5 = (uncovered5bgn + UNCOVERED_NOISE_FILTER < uncovered5end);
      bool  save3 = (uncovered3bgn + UNCOVERED_NOISE_FILTER < uncovered3end);

      if (save5 && save3) {
        //  Ignore this repeat-in-a-fragment alignment.  The only way (I can think of anyway) this could
        //  occur is from an alignment to a short contained fragment.
        *end5 = FragmentEnd();
        *end3 = FragmentEnd();

      } else if (save5) {
        //  Uncovered bases towards the start of the unitig.
        //hang5.add(uncovered5bgn, uncovered5end - uncovered5bgn + 1);
        places.push_back(op[pl]);
        align5.add(coveredbgn + COVERED_BORDER, coveredend - coveredbgn - 2 * COVERED_BORDER);
        *end3 = FragmentEnd();

      } else if (save3) {
        //  Uncovered bases towards the end of the unitig.
        //hang3.add(uncovered3bgn, uncovered3end - uncovered3bgn + 1);
        places.push_back(op[pl]);
        align3.add(coveredbgn + COVERED_BORDER, coveredend - coveredbgn - 2 * COVERED_BORDER);
        *end5 = FragmentEnd();

      } else {
        *end5 = FragmentEnd();
        *end3 = FragmentEnd();
      }

#if 1
      fprintf(logFile, "markRepeats()-- tig frag %8d ovl frag %8d %6d-%6d %5.2f%% tig pos %8d/%c' (%4d) %6d-%6d (%4d) %8d/%c'%s\n",
              op[pl].refID, op[pl].frgID,
              op[pl].covered.bgn,
              op[pl].covered.end,
              op[pl].fCoverage * 100.0,
              end5->fragId(), (end5->frag5p() ? '5' : '3'),
              uncovered5end - uncovered5bgn,
              coveredbgn, coveredend,
              uncovered3end - uncovered3bgn,
              end3->fragId(), (end3->frag5p() ? '5' : '3'),
              (save5 && save3) ? " ***" : "");
#endif
    }
  }

  align5.merge();
  align3.merge();

  //  Combine the begin-point and end-pointing repeat/unique junctions into one list.  This will let
  //  us ignore those that indicate a small repeat spanned by a single fragment, and those that
  //  are too close to the end of the unitig.

  vector<repeatUniqueJunction>  junctions;
  vector<repeatUniqueJunction>  junctionss;  //  Yuck.  'junctions' 'saved'

  //  FIRST FILTER -- discard any junctions with too few overlaps.

  //  align5 are those alignments that have an unaligned portion towards the start of the unitig.
  //  Thus, the repeat will be after the break point (rptLeft == false).  At this time, we don't
  //  know precisely what the break point will be.  The most conservative guess -- leaving too much
  //  of the repeat attached to the unique region -- is the hi() coord,

  for (uint32 i=0; i<align5.numberOfIntervals(); i++)
    if (align5.ct(i) > SPURIOUS_REPEAT_THRESHOLD)
      junctions.push_back(repeatUniqueJunction(align5.lo(i) - COVERED_BORDER, align5.hi(i) + COVERED_BORDER,
                                               align5.hi(i) + COVERED_BORDER,
                                               align5.ct(i),
                                               false));

  for (uint32 i=0; i<align3.numberOfIntervals(); i++)
    if (align3.ct(i) > SPURIOUS_REPEAT_THRESHOLD)
      junctions.push_back(repeatUniqueJunction(align3.lo(i) - COVERED_BORDER, align3.hi(i) + COVERED_BORDER,
                                               align3.lo(i) - COVERED_BORDER,
                                               align3.ct(i),
                                               true));

  if (junctions.size() == 0)
    //  Nothing interesting aligned.
    return;

  sort(junctions.begin(), junctions.end());

  fprintf(logFile, "markRepeats()--  unitig %d has %lu interesting intervals:\n", target->id(), junctions.size());

  //  Search for pairs that look like they are covered by a single fragment.
  //  This is kinda hard....

  for (uint32 bp=0; bp<junctions.size(); bp++)
    fprintf(logFile, "markRepeats()--  junction: %s %6d %s from %d fragments\n",
            junctions[bp].rptLeft ? "repeat" : "unique",
            junctions[bp].point,
            junctions[bp].rptLeft ? "unique" : "repeat",
            junctions[bp].overlaps);


  //  this is a stupid simple filter.  if there is another junction within STUPID_REPEAT_SIZE, throw
  //  out BOTH junctions.
  //
#if 1
  for (uint32 ba=0; ba<junctions.size(); ba++) {
    bool  keepIt = true;

    int32 lo = junctions[ba].point - STUPID_REPEAT_SIZE;
    int32 hi = junctions[ba].point + STUPID_REPEAT_SIZE;

    for (uint32 bb=0; bb<junctions.size(); bb++)
      if ((ba != bb) &&
          (lo <= junctions[bb].point) &&
          (junctions[bb].point <= hi))
        keepIt = false;

    if (keepIt) {
      fprintf(logFile, "SAVE junction: %s %6d %s from %d fragments\n",
            junctions[ba].rptLeft ? "repeat" : "unique",
            junctions[ba].point,
            junctions[ba].rptLeft ? "unique" : "repeat",
            junctions[ba].overlaps);
      junctionss.push_back(junctions[ba]);

    } else {
      fprintf(logFile, "JUNK junction: %s %6d %s from %d fragments\n",
            junctions[ba].rptLeft ? "repeat" : "unique",
              junctions[ba].point,
            junctions[ba].rptLeft ? "unique" : "repeat",
            junctions[ba].overlaps);
    }
  }

  junctions = junctionss;
#endif


  //  For each remaining break point, run through all the placements again.  For any placement that
  //  overlaps with the breakPoint, find the unitig fragment that we should break at to excise this
  //  repeat.  The picture is unfortunately hard to draw in ASCII.

  for (uint32 bp=0; bp<junctions.size(); bp++) {
    int32        breakPt = junctions[bp].point;
    bool         rptLeft = junctions[bp].rptLeft;
    double       avgCov  = computeAverageCoverage(target, breakPt - 10, breakPt + 10); 

    uint32       breakOrd = (rptLeft) ? UINT32_MAX : 0;
    uint32       breakPos = (rptLeft) ? UINT32_MAX : 0;
    FragmentEnd  breakEnd;

    for (uint32 pl=0; pl<places.size(); pl++) {
      assert(places[pl].tigID == target->id());
      assert(places[pl].fCoverage <= 0.99);

      //  Compute the min/max location of the verified portion, and remember the fragEnd of this
      //  placement.

      uint32          plMin,  plMax;
      FragmentEnd    *loEnd, *hiEnd;

      if (places[pl].position.bgn < places[pl].position.end) {
        plMin = places[pl].verified.bgn;
        plMax = places[pl].verified.end;
        loEnd = &places[pl].frag5p;
        hiEnd = &places[pl].frag3p;
      } else {
        plMin = places[pl].verified.end;
        plMax = places[pl].verified.bgn;
        loEnd = &places[pl].frag3p;
        hiEnd = &places[pl].frag5p;
      }

      //  (for rptLeft == false) Repeat to the right of the break point, so we should be looking for
      //  a placement to the left of the point we set at junctions.push_back() above.  We still need
      //  to ensure that this placement has a spur (indicated by a fragEnd) towards the start of the
      //  unitig, and isn't by coincidence.

      if (rptLeft == false) {
        if ((breakPt - 1 > plMax) ||  //  Break point doesn't overlap,
            (breakPt + 1 < plMax) ||  //  or no spur towards the start
            (loEnd->fragId() == 0))   //  of the unitig.
          continue;

      } else {
        if ((breakPt - 1 > plMin) ||
            (breakPt + 1 < plMin) ||
            (hiEnd->fragId() == 0))
          continue;
      }

      //  This placement is one that we care about.

      uint32  ord = target->pathPosition(places[pl].refID);
      ufNode &frg = target->ufpath[ord];
      uint32  pos = MIN(frg.position.bgn, frg.position.end);

      assert(frg.ident == places[pl].refID);

      fprintf(logFile, "markRepeats()-- breakPt=%d verified=%d,%d refID=%d frgID=%d (tig %d len "F_SIZE_T") frgPos=%d,%d\n",
              breakPt,
              places[pl].verified.bgn, places[pl].verified.end,
              places[pl].refID,
              places[pl].frgID, Unitig::fragIn(places[pl].frgID), unitigs[Unitig::fragIn(places[pl].frgID)]->ufpath.size(),
              frg.position.bgn, frg.position.end);

      //  Still thinking about rtpLeft == false, we want to find the last fragment in the unitig
      //  that supports this placement.  Then we want to break on the low end of that fragment -- if the fragment is
      //  forward, the low end is the 5' end, if reversed, the low end is 3'.
      //
      //  For rptLeft == true, we are looking for the earliest fragment that supports the placement,
      //  and we want to split on the high end of the fragment.

      if (rptLeft == false) {
        if (ord > breakOrd) {
          breakOrd = ord;
          breakPos = pos;
          breakEnd = FragmentEnd(frg.ident, (frg.position.bgn > frg.position.end));
        }

      } else {
        if (ord < breakOrd) {
          breakOrd = ord;
          breakPos = pos;
          breakEnd = FragmentEnd(frg.ident, (frg.position.bgn < frg.position.end));
        }
      }
    }

    if (breakEnd != FragmentEnd()) {
      fprintf(logFile, "markRepeats()--    %s-%6d-%s %.2f coverage -- BREAK %u/%c'\n",
              (rptLeft) ? "repeat" : "unique",
              breakPt,
              (rptLeft) ? "unique" : "repeat",
              avgCov,
              breakEnd.fragId(), breakEnd.frag3p() ? '3' : '5');
      breaks.push_back(breakPoint(breakEnd.fragId(), breakEnd.frag3p(), false, false));
    }
  }

  //  All break points added, now just break.

  //  Until we devise a method for filtering out spurious repeat/unique junctions, we instead
  //  analyze the break points and remove any that are near the end of the unitig.

#if 1
  vector<breakPoint>        fbreaks;

  uint32  lastStart = 0;

  for (uint32 fi=0; fi<target->ufpath.size(); fi++)
    if (lastStart < MIN(target->ufpath[fi].position.bgn, target->ufpath[fi].position.end))
      lastStart = MIN(target->ufpath[fi].position.bgn, target->ufpath[fi].position.end);

  for (uint32 i=0; i<breaks.size(); i++) {
    uint32  iid = breaks[i].fragEnd.fragId();
    uint32  i3p = breaks[i].fragEnd.frag3p();
    uint32  idx = target->pathPosition(iid);
    ufNode &frg = target->ufpath[idx];

    //  The 3' of the fragment is ALWAYS the 'end' coord:
    //    (bgn) 5' -------> 3' (end)
    //    (end) 3' <------- 5' (bgn)
    int32  pos = (i3p) ? frg.position.end : frg.position.bgn;


    if ((pos < 100) ||
        (idx < 2) ||
        (pos >= lastStart - 100) ||
        (idx >= target->ufpath.size() - 2))
      //fprintf(logFile, "BREAK pos idx=%u pos=%u tigLen idx="F_SIZE_T" pos=%u FILTERED\n",
      //        idx, pos, target->ufpath.size(), target->getLength());
      continue;

    fbreaks.push_back(breaks[i]);
  }

  breaks = fbreaks;
#endif

  UnitigVector *newTigs = breakUnitigAt(target, breaks);

  if (newTigs != NULL) {
    unitigs[target->id()] = NULL;
    delete target;

    unitigs.insert(unitigs.end(), newTigs->begin(), newTigs->end());
  }
}






void
markChimera(UnitigVector &unitigs,
            Unitig *target) {
}



void
mergeSplitJoin(UnitigVector &unitigs) {

  //logFileFlags |= LOG_PLACE_FRAG;
  //logFileFlags &= ~LOG_PLACE_FRAG;

  //  BUILD A LIST OF ALL INTERSECTIONS - build a reverse mapping of all BestEdges that are between
  //  unitigs.  For each fragment, we want to have a list of the incoming edges from other unitigs.

  intersectionList  *ilist = new intersectionList(unitigs);

  //ilist->logIntersections();

#if 0
  {
    Unitig        *target = unitigs[5];

    fprintf(logFile, "popBubbles()-- WORKING on unitig %d/"F_SIZE_T" with %ld fragments.\n",
            target->id(), unitigs.size(), target->ufpath.size());

    mergeBubbles(unitigs, target, ilist);
    stealBubbles(unitigs, target, ilist);
    markRepeats(unitigs, target);
    markChimera(unitigs, target);
    exit(1);
  }
#endif

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig        *target = unitigs[ti];

    if (target == NULL)
      continue;

    if (target->ufpath.size() < 5)
      continue;

    fprintf(logFile, "popBubbles()-- WORKING on unitig %d/"F_SIZE_T" with %ld fragments.\n",
            target->id(), unitigs.size(), target->ufpath.size());

    //  MERGE BUBBLES - for every unitig with at least one incoming edge (and possibly no edges
    //  to other unitigs) try to merge.  The incoming unitig must NOT extend the current unitig.

    mergeBubbles(unitigs, target, ilist);

    //  STEAL BUBBLES (PRE-JOIN) - for dovetail overlapping unitigs, try to merge one or the other
    //  end into the opposite unitig.  Requires a best overlap to force the intersection, then
    //  the split off piece can be attempted to be merged into the other unitig.

    stealBubbles(unitigs, target, ilist);

    //  MARK REPEATS and CHIMERA - using just overlaps, look for fragments that are partially
    //  covered by this unitig.

    markRepeats(unitigs, target);
    markChimera(unitigs, target);
  }

  //  JOIN EXPOSED BEST - after bubbles are stolen, this should leave some unitigs
  //  with exposed best edges that can now be connected.

  //  do we need to re-mark repeats after joining?

  //  SPLIT MARKED REPEATS - 

  //  SPLIT MARKED CHIMERA - 

  //  MERGE LEFTOVERS - these are the leftover pieces after repeats/chimera are split.  Hopefully
  //  they'll just be low coverage spurs

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig        *target = unitigs[ti];

    if (target == NULL)
      continue;

    mergeBubbles(unitigs, target, ilist);
  }

  delete ilist;

  logFileFlags &= ~LOG_PLACE_FRAG;
}

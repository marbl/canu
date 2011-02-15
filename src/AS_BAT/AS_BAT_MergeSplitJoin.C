
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

static const char *rcsid = "$Id: AS_BAT_MergeSplitJoin.C,v 1.1 2011-02-15 08:10:11 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_IntersectSplit.H"

#include "AS_BAT_MergeSplitJoin.H"


//  Over every frag in 'target' see if the unitig from any invading fragment can
//  be merged in.



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

  placeFragUsingOverlaps(unitigs, fFrg.ident, placements);
  for (uint32 i=0; i<placements.size(); i++) {
    if (placements[i].tigID != target->id())
      continue;

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

  placeFragUsingOverlaps(unitigs, lFrg.ident, placements);
  for (uint32 i=0; i<placements.size(); i++) {
    if (placements[i].tigID != target->id())
      continue;

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

  fprintf(logFile, "popBubbles()--   Order/Orientation problem.  bL %d bR %d bOrd %d  nL %d nR %d nOrd %d\n",
          bL, bR, bOrd,
          nL, nR, nOrd);

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

    placeFragUsingOverlaps(unitigs, frg->ident, placements[fi]);

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
    uint32  nNotPlacedInTarget = 0;
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
      if (placements[fi][pl].tigID != target->id()) {
        nNotPlacedInTarget++;
        continue;
      }

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
      //fprintf(logFile, "popBubbles()--   Failed to place frag %d notPlaced %d notPlacedInTarget %d notPlacedInCorrectPosition %d notPlacedFully %d notOriented %d\n",
      //        bubble->ufpath[fi].ident, nNotPlaced, nNotPlacedInTarget, nNotPlacedInCorrectPosition, nNotPlacedFully, nNotOriented);
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

    target->addFrag(nFrg, 0, logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG));
  }

  target->sort();

  success = true;

  fprintf(logFile, "popBubbles()--   merged bubble unitig %d into unitig %d\n",
          bubble->id(), target->id());

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
  vector<overlapPlacement>  placements;


  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode     *frg   = &target->ufpath[fi];
    uint32      ovlLen = 0;
    BAToverlap *ovl    = OC->getOverlaps(frg->ident, ovlLen);

    for (uint32 i=0; i<ovlLen; i++) {
      if (Unitig::fragIn(ovl[i].b_iid) != target->id())
        ovlFrags.insert(ovl[i].b_iid);
    }
  }

  //  As much as I dislike C++ iterators, sometimes they just gotta be used.

  for (set<AS_IID>::iterator it=ovlFrags.begin(); it!=ovlFrags.end(); it++) {
    AS_IID  iid = *it;

    placeFragUsingOverlaps(unitigs, iid, placements);

    for (uint32 pl=0; pl<placements.size(); pl++) {
      if (placements[pl].tigID != target->id())
        //  No worries, just not placed in the unitig we're looking at.
        continue;

      if (placements[pl].fCoverage > 0.99)
        //  No worries, fully placed.
        continue;

      //  Otherwise, placed in the target unitig, at less than perfect coverage.

      int32   unplaced5    = placements[pl].covered.bgn;
      int32   unplaced3    = FI->fragmentLength(iid) - placements[pl].covered.end;

      int32   unplaced5bgn = 0;
      int32   unplaced5end = 0;

      int32   unplaced3bgn = 0;
      int32   unplaced3end = 0;

      if (placements[pl].position.bgn < placements[pl].position.end) {
        //  Fragment is placed forward.
        unplaced5    = placements[pl].covered.bgn;
        unplaced3    = FI->fragmentLength(iid) - placements[pl].covered.end;

        if (unplaced5 > 0) {
          unplaced5bgn = placements[pl].position.bgn - unplaced5;
          unplaced5end = placements[pl].position.bgn;
        }

        if (unplaced3 > 0) {
          unplaced3bgn = placements[pl].position.end;
          unplaced3end = placements[pl].position.end + unplaced3;
        }

      } else {
        //  Fragment is placed reverse.
        unplaced5    = FI->fragmentLength(iid) - placements[pl].covered.end;
        unplaced3    = placements[pl].covered.bgn;

        if (unplaced5 > 0) {
          unplaced5bgn = placements[pl].position.end - unplaced5;
          unplaced5end = placements[pl].position.end;
        }

        if (unplaced3 > 0) {
          unplaced3bgn = placements[pl].position.bgn;
          unplaced3end = placements[pl].position.bgn + unplaced3;
        }
      }


      fprintf(logFile, "markRepeats()-- unitig %d fragment %d bases %d-%d coverage %5.2f%% at unitig position (%d-%d) %d-%d (%d-%d)\n",
              target->id(), iid,
              placements[pl].covered.bgn,
              placements[pl].covered.end,
              placements[pl].fCoverage * 100.0,
              unplaced5bgn, unplaced5end,
              placements[pl].position.bgn, placements[pl].position.end,
              unplaced3bgn, unplaced3end);
    }
  }
}



void
markChimera(UnitigVector &unitigs,
            Unitig *target) {
}



void
mergeSplitJoin(UnitigVector &unitigs) {

  logFileFlags |= LOG_PLACE_FRAG;

  //  BUILD A LIST OF ALL INTERSECTIONS - build a reverse mapping of all BestEdges that are between
  //  unitigs.  For each fragment, we want to have a list of the incoming edges from other unitigs.

  intersectionList  *ilist = new intersectionList(unitigs);

  ilist->logIntersections();

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig        *target = unitigs[ti];

    if (target == NULL)
      continue;

    fprintf(logFile, "popBubbles()-- WORKING on unitig %d with %ld fragments.\n", target->id(), target->ufpath.size());

    //  MERGE BUBBLES - for every unitig with at least one incoming edge (and possibly no edges
    //  to other unitigs) try to merge.  The incoming unitig must NOT extend the current unitig.

    mergeBubbles(unitigs, target, ilist);

    //  STEAL BUBBLES (PRE-JOIN) - for dovetail overlapping unitigs, try to merge one or the other
    //  end into the opposite unitig.  Requires a best overlap to force the intersection, then
    //  the split off piece can be attempted to be merged into the other unitig.

    //stealBubbles(unitigs, target, ilist);

    //  MARK REPEATS and CHIMERA - using just overlaps, look for fragments that are partially
    //  covered by this unitig.

    //markRepeats(unitigs, target);
    //markChimera(unitigs, target);
  }

  //  JOIN EXPOSED BEST - after bubbles are stolen, this should leave some unitigs
  //  with exposed best edges that can now be connected.

  //  do we need to re-mark repeats after joining?

  //  SPLIT MARKED REPEATS - 

  //  SPLIT MARKED CHIMERA - 

  //  MERGE LEFTOVERS - these are the leftover pieces after repeats/chimera are split.  Hopefully
  //  they'll just be low coverage spurs

  delete ilist;

  logFileFlags &= ~LOG_PLACE_FRAG;
}

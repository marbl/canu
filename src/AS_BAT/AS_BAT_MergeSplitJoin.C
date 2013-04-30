
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

static const char *rcsid = "$Id: AS_BAT_MergeSplitJoin.C,v 1.20 2012-08-06 23:36:44 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_IntersectSplit.H"

#include "AS_BAT_InsertSizes.H"

#include "AS_BAT_MergeSplitJoin.H"

#include "AS_BAT_Breaking.H"

#include "AS_BAT_Instrumentation.H"
#include "AS_BAT_EvaluateMates.H"

#include "AS_UTL_intervalList.H"


//  Over every frag in 'target' see if the unitig from any invading fragment can
//  be merged in.

#define FRAG_COVERS_REPEAT           10  //
#define REPEAT_COVERS_FRAG            3  //
#define SPURIOUS_REPEAT_THRESHOLD     5  //  repeat interval must have this many frags aligned to it
#define ISECT_NEEDED_TO_BREAK         5

#include "AS_BAT_RepeatJunctionEvidence.H"


#undef  LOG_BUBBLE_TESTS
#undef  LOG_BUBBLE_FAILURE
#define LOG_BUBBLE_SUCCESS


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
#ifdef LOG_BUBBLE_TESTS
    writeLog("popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments STARTS WITH A CONTAINED FRAGMENT %d\n",
            bubble->id(), bubble->getLength(), bubble->ufpath.size(),
            bubble->ufpath[0].ident);
#endif
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

#ifdef LOG_BUBBLE_TESTS
  if ((fUtg != 0) && (lUtg != 0))
    writeLog("popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  Edges (%d/%d') from frag %d/%d' and (%d/%d') from frag %d/%d'\n",
            bubble->id(), bubble->getLength(), bubble->ufpath.size(),
            fEdge->fragId(), (fEdge->frag3p() ? 3 : 5), fFrg.ident, (f3p ? 3 : 5),
            lEdge->fragId(), (lEdge->frag3p() ? 3 : 5), lFrg.ident, (l3p ? 3 : 5));
  else if (fUtg != 0)
    writeLog("popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  Edge (%d/%d') from frag %d/%d'\n",
            bubble->id(), bubble->getLength(), bubble->ufpath.size(),
            fEdge->fragId(), (fEdge->frag3p() ? 3 : 5), fFrg.ident, (f3p ? 3 : 5));
  else if (lUtg != 0)
    writeLog("popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  Edge (%d/%d') from frag %d/%d'\n",
            bubble->id(), bubble->getLength(), bubble->ufpath.size(),
            lEdge->fragId(), (lEdge->frag3p() ? 3 : 5), lFrg.ident, (l3p ? 3 : 5));
  else
    //  But then how do we get an intersection?!?!!  Intersections from a bubble that was
    //  already popped.  We pop A into B, and while iterating through fragments in B we find
    //  the -- now obsolete -- intersections we originally used and try to pop it again.
    //
    writeLog("popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  NO EDGES, no bubble.\n",
            bubble->id(), bubble->getLength(), bubble->ufpath.size());
#endif

  if ((fUtg == 0) && (lUtg == 0))
    return(false);

#ifdef LOG_BUBBLE_TESTS
  //  The only interesting case here is if we have both edges and they point to different unitigs.
  //  We might as well place it aggressively.

  if ((fUtg != 0) && (lUtg != 0) && (fUtg != lUtg)) {
    writeLog("popBubbles()--   bubble unitig %d has edges to both unitig %d and unitig %d\n",
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

  placeFragUsingOverlaps(unitigs, target, fFrg.ident, placements);

#ifdef LOG_BUBBLE_TESTS
  writeLog("popBubbles()-- fFrg %u has %u potential placements.\n", fFrg.ident, placements.size());
#endif

  for (uint32 i=0; i<placements.size(); i++) {
    assert(placements[i].tigID == target->id());

    if (placements[i].fCoverage < 0.99) {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- fFrg %u low coverage %f at %u,%u\n",
              fFrg.ident,
              placements[i].fCoverage,
              placements[i].position.bgn, placements[i].position.end);
#endif
      continue;
    } else {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- fFrg %u GOOD coverage %f at %u,%u\n",
              fFrg.ident,
              placements[i].fCoverage,
              placements[i].position.bgn, placements[i].position.end);
#endif
    }

    if (placements[i].errors / placements[i].aligned < fFrgPlacement.errors / fFrgPlacement.aligned) {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- fFrg %u GOOD identity %f at %u,%u\n",
              fFrg.ident,
              placements[i].errors / placements[i].aligned,
              placements[i].position.bgn, placements[i].position.end);
#endif
      fFrgPlacement = placements[i];
    } else {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- fFrg %u low identity %f at %u,%u\n",
              fFrg.ident,
              placements[i].errors / placements[i].aligned,
              placements[i].position.bgn, placements[i].position.end);
#endif
    }
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
#ifdef LOG_BUBBLE_FAILURE
    writeLog("popBubbles()--   failed to place fFrg.\n");
#endif
    return(false);
  }


  placements.clear();

  placeFragUsingOverlaps(unitigs, target, lFrg.ident, placements);

#ifdef LOG_BUBBLE_TESTS
  writeLog("popBubbles()-- lFrg %u has %u potential placements.\n", lFrg.ident, placements.size());
#endif

  for (uint32 i=0; i<placements.size(); i++) {
    assert(placements[i].tigID == target->id());

    if (placements[i].fCoverage < 0.99) {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- lFrg %u low coverage %f at %u,%u\n",
              lFrg.ident,
              placements[i].fCoverage,
              placements[i].position.bgn, placements[i].position.end);
#endif
      continue;
    } else {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- lFrg %u GOOD coverage %f at %u,%u\n",
              lFrg.ident,
              placements[i].fCoverage,
              placements[i].position.bgn, placements[i].position.end);
#endif
    }

    if (placements[i].errors / placements[i].aligned < lFrgPlacement.errors / lFrgPlacement.aligned) {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- lFrg %u GOOD identity %f at %u,%u\n",
              lFrg.ident,
              placements[i].errors / placements[i].aligned,
              placements[i].position.bgn, placements[i].position.end);
#endif
      lFrgPlacement = placements[i];
    } else {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- lFrg %u low identity %f at %u,%u\n",
              lFrg.ident,
              placements[i].errors / placements[i].aligned,
              placements[i].position.bgn, placements[i].position.end);
#endif
    }
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
#ifdef LOG_BUBBLE_FAILURE
    writeLog("popBubbles()--   failed to place lFrg.\n");
#endif
    return(false);
  }


  int32 minL = MIN(fFrg.position.bgn, fFrg.position.end);
  int32 maxL = MAX(fFrg.position.bgn, fFrg.position.end);

  int32 minR = MIN(lFrg.position.bgn, lFrg.position.end);
  int32 maxR = MAX(lFrg.position.bgn, lFrg.position.end);

  int32 placedLen = MAX(maxL, maxR) - MIN(minL, minR);

  if (2 * placedLen < bubble->getLength()) {
    //  Too short.
#ifdef LOG_BUBBLE_FAILURE
    writeLog("popBubbles()--   too short.  fFrg %d,%d lFrg %d,%d.  L %d,%d R %d,%d len %d\n",
            fFrg.position.bgn, fFrg.position.end,
            lFrg.position.bgn, lFrg.position.end,
            minL, maxL, minR, maxR, placedLen);
#endif
    return(false);
  }

  if (2 * bubble->getLength() < placedLen) {
    //  Too long.
#ifdef LOG_BUBBLE_FAILURE
    writeLog("popBubbles()--   too long.  fFrg %d,%d lFrg %d,%d.  L %d,%d R %d,%d len %d\n",
            fFrg.position.bgn, fFrg.position.end,
            lFrg.position.bgn, lFrg.position.end,
            minL, maxL, minR, maxR, placedLen);
#endif
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

#ifdef LOG_BUBBLE_FAILURE
  writeLog("popBubbles()--   Order/Orientation problem.  bL %d bR %d bOrd %d  nL %d nR %d nOrd %d\n",
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

    placeFragUsingOverlaps(unitigs, target, frg->ident, placements[fi]);

    //  Initialize the final placement to be bad, so we can pick the best.
    correctPlace[fi].fCoverage = 0.0;
    correctPlace[fi].errors    = 4.0e9;
    correctPlace[fi].aligned   = 1;
  }

  //  Some bizarre cases -- possibly only from bad data -- confound any logical attempt at finding
  //  the min/max extents.  Yes, even though this should work, it doesn't.  Or maybe it's just
  //  broken and I haven't seen how.
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
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()--   Failed to place frag %d notPlaced %d notPlacedInCorrectPosition %d notPlacedFully %d notOriented %d\n",
              bubble->ufpath[fi].ident, nNotPlaced, nNotPlacedInCorrectPosition, nNotPlacedFully, nNotOriented);
#endif
      break;
    }
  }

  if (nCorrect != bubble->ufpath.size()) {
    goto finished;
  }

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

#ifdef LOG_BUBBLE_SUCCESS
  writeLog("popBubbles()--   merged bubble unitig %d with %ld frags into unitig %d now with %ld frags\n",
          bubble->id(), bubble->ufpath.size(), target->id(), target->ufpath.size());
#endif

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

      if (bubble == NULL)
        //  Whoops!  Unitig was repeat/unique split and the repeats shattered
        continue;

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





void
markRepeats_buildOverlapList(Unitig *target, set<AS_IID> &ovlFrags) {

  ovlFrags.clear();

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode     *frg   = &target->ufpath[fi];
    uint32      ovlLen = 0;
    BAToverlap *ovl    = OC->getOverlaps(frg->ident, ovlLen);

    for (uint32 i=0; i<ovlLen; i++) {
      if (Unitig::fragIn(ovl[i].b_iid) != target->id())
        ovlFrags.insert(ovl[i].b_iid);
    }
  }
}




void
markRepeats_computeUnitigErrorRate(UnitigVector &unitigs,
                                   Unitig       *target,
                                   double       &meanError,
                                   double       &stddevError) {

  vector<overlapPlacement>  op;
  vector<double>            error;

  meanError   = 0;
  stddevError = 0;

#undef DUMPERROR
#ifdef DUMPERROR
  char  N[FILENAME_MAX];
  sprintf(N, "error.%08d.dat", target->id());
  FILE *F = fopen(N, "w");
#endif

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode     *frg   = &target->ufpath[fi];
    uint32      bgn   = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
    uint32      end   = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

    placeFragUsingOverlaps(unitigs, target, frg->ident, op);

    if (op.size() == 0)
      //  Huh?  Couldn't be placed in my own unitig?
      continue;
    assert(op.size() > 0);

    double  minError = 1.0;
    double  corError = 1.0;

    for (uint32 pl=0; pl<op.size(); pl++) {
      assert(op[pl].tigID == target->id());

      double e = op[pl].errors / op[pl].aligned;

      minError = MIN(minError, e);

      if (op[pl].position.bgn < op[pl].position.end) {
        if ((op[pl].position.end < bgn) ||
            (end < op[pl].position.bgn))
          continue;
      } else {
        if ((op[pl].position.bgn < bgn) ||
            (end < op[pl].position.end))
          continue;
      }

      corError = MIN(corError, e);
    }

#ifdef DUMPERROR
    fprintf(F, "%f\t%f\n", corError, minError);
#endif

    if (corError < 1.0)
      error.push_back(corError);
    else
      error.push_back(minError);
  }

#ifdef DUMPERROR
  fclose(F);
#endif

  for (uint32 i=0; i<error.size(); i++)
    meanError += error[i];

  meanError /= error.size();

  for (uint32 i=0; i<error.size(); i++)
    stddevError += (error[i] - meanError) * (error[i] - meanError);

  stddevError = sqrt(stddevError / error.size());

  writeLog("markRepeats_computeUnitigErrorRate()--  tig %d error %f +- %f\n",
          target->id(), meanError, stddevError);
}



void
markRepeats_placeAndProcessOverlaps(UnitigVector                     &unitigs,
                                    Unitig                           *target,
                                    double                           meanError,
                                    double                           stddevError,
                                    set<AS_IID>                      &ovlFrags,
                                    intervalList                     &aligned,
                                    vector<repeatJunctionEvidence>   &evidence) {

  aligned.clear();
  evidence.clear();

  for (set<AS_IID>::iterator it=ovlFrags.begin(); it!=ovlFrags.end(); it++) {
    AS_IID  iid = *it;

    vector<overlapPlacement>  op;

    placeFragUsingOverlaps(unitigs, target, iid, op);

    //  placeFragUsingOverlaps() returns the expected placement for this fragment in 'position', and
    //  the amount of the fragment covered by evidence in 'covered'.
    //
    //  Below we'll try to decipher this into two intervals, covered by evidence and not covered by
    //  evidence.

    for (uint32 pl=0; pl<op.size(); pl++) {
      assert(op[pl].tigID == target->id());

      double erate        = op[pl].errors / op[pl].aligned;
      bool   erateTooHigh = (meanError + 3 * stddevError < erate);

#if 0
      writeLog("markRepeats()-- op[%3d] tig %d frag %d fCoverage %f position %d %d verified %d %d erate %f%s\n",
              pl, op[pl].tigID, op[pl].frgID, op[pl].fCoverage,
              op[pl].position.bgn, op[pl].position.end,
              op[pl].verified.bgn, op[pl].verified.end,
              erate,
              erateTooHigh ? " - TOO HIGH" : "");
#endif

      if (erateTooHigh)
        continue;

      //  Save the aligned portion.  This will be used in conjunction with the partially aligned
      //  fragments to pick out where to split.
      //
      assert(op[pl].verified.bgn >= 0);
      assert(op[pl].verified.bgn <= target->getLength());  //  is bgn unsigned and underflowed?
      assert(op[pl].verified.end <= target->getLength());

      if (op[pl].verified.bgn < op[pl].verified.end)
        aligned.add(op[pl].verified.bgn, op[pl].verified.end - op[pl].verified.bgn);
      else
        aligned.add(op[pl].verified.end, op[pl].verified.bgn - op[pl].verified.end);

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

      repeatJunctionEvidence   ev(target, op[pl]);

      if (ev.tigFrag == FragmentEnd())
        //  Didn't pass muster; weak overhangs likely.
        continue;

      evidence.push_back(ev);
    }
  }
}





void
markRepeats_filterIntervalsSpannedByFragment(Unitig                    *target,
                                             intervalList              &aligned,
                                             vector<repeatRegion>      &regions) {
  
  regions.clear();
  aligned.merge();

  //  Discard the interval if there is a fragment that contains it with enough overhang to
  //  unambiguously place it in a unitig.
  
  for (uint32 i=0; i<aligned.numberOfIntervals(); i++) {
    bool     isWeak       = (aligned.ct(i) < SPURIOUS_REPEAT_THRESHOLD);
    ufNode  *containedIn  = NULL;

    if (isWeak == false) {
      for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
        ufNode     *frg   = &target->ufpath[fi];
        uint32      bgn   = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
        uint32      end   = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

        if (end < aligned.lo(i))
          //  Fragment is before the region.
          continue;

        if (aligned.hi(i) < bgn)
          //  Fragment is after the region, we're finished.
          break;

        if ((bgn + FRAG_COVERS_REPEAT < aligned.lo(i)) &&
            (aligned.hi(i)            < end - FRAG_COVERS_REPEAT)) {
          //  Fragment contains the region with no acceptable overhangs into the non-repeat area.
          containedIn  = frg;
          break;
        }
      }
    }

    if (containedIn != NULL) {
      writeLog("markRepeats()--  repeat alignment "F_U64","F_U64" count "F_U32" DISCARD - CONTAINED IN FRAGMENT "F_U32" "F_S32","F_S32"\n",
              aligned.lo(i), aligned.hi(i), aligned.ct(i),
              containedIn->ident,
              containedIn->position.bgn, containedIn->position.end);
      continue;
    }

    if (isWeak == true) {
      writeLog("markRepeats()--  repeat alignment "F_S64","F_S64" count %u DISCARD - TOO WEAK\n",
              aligned.lo(i), aligned.hi(i), aligned.ct(i));
      continue;
    }

    if ((regions.size() == 0) ||
        (regions.back().end + 100 < aligned.lo(i))) {
      regions.push_back(repeatRegion(aligned.lo(i), aligned.hi(i), aligned.ct(i)));

      writeLog("markRepeats()--  repeat alignment "F_U32","F_U32" count "F_U32"\n",
              regions.back().bgn, regions.back().end, regions.back().overlaps);

    } else {
      regions.back().bgn       = MIN(regions.back().bgn, aligned.lo(i));
      regions.back().end       = MAX(regions.back().end, aligned.hi(i));
      regions.back().overlaps += aligned.ct(i);

      writeLog("markRepeats()--  repeat alignment "F_U32","F_U32" count "F_U32" (merged from "F_U64","F_U64")\n",
              regions.back().bgn, regions.back().end, regions.back().overlaps,
              aligned.lo(i), aligned.hi(i));
    }

  }
}



void
markRepeats_filterIntervalsSpannedByMates(Unitig                    *target,
                                          vector<repeatRegion>      &regions) {

    //  Save the interval if it is small and has 'enough' mates spanning to convince us it is correct.  Spanning
    //  mates must be anchored in non-repeat marked areas.  If we choose to keep it, eject any fragment that
    //  is externally mated.  Also eject any non-mated fragments.

  for (uint32 i=0; i<regions.size(); i++) {
    assert(regions[i].bgn < regions[i].end);

    if (regions[i].end - regions[i].bgn > 500)
      //  Arbitrarily too large to be saved.
      continue;

    uint32        spanGood       = 0;  //  Mate spans the region
    uint32        spanGoodMaybe  = 0;  //  Mate spans the region, but is slightly off in position
    uint32        spanBad        = 0;  //  Mate should span, but is missing the other frag

    uint32        tigLen = target->getLength();

    for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
      ufNode     *frg      = &target->ufpath[fi];
      uint32      frgbgn   = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
      uint32      frgend   = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;
      bool        frg53    = (frg->position.bgn < frg->position.end) ? true : false;

      AS_IID      mid      = FI->mateIID(frg->ident);
      uint32      mrgtig   = target->fragIn(mid);
      ufNode     *mrg      = NULL;
      uint32      mrgbgn   = 0;
      uint32      mrgend   = 0;
      bool        mrg53    = false;

      if (mid == 0)
        //  Fragment is not mated.  Won't help us.
        continue;

      AS_IID      lib      = FI->libraryIID(frg->ident);

#if 1
      uint32      minD     = IS->mean(lib) - 3 * IS->stddev(lib);
      uint32      maxD     = IS->mean(lib) + 3 * IS->stddev(lib);
      uint32      stddev   = IS->stddev(lib);
#else
      uint32      minD     = FI->mean(lib) - 3 * FI->stddev(lib);
      uint32      maxD     = FI->mean(lib) + 3 * FI->stddev(lib);
      uint32      stddev   = FI->stddev(lib);
#endif

      bool        isGood   = false;
      bool        isBad    = false;

      bool        expout   = false;  //  Expected location is outside the unitig
      bool        expspan  = false;  //  Expected location spans the region
      uint32      expbgn   = 0;
      uint32      expend   = 0;

      if (frg53 == true) {
        expbgn = frgbgn + minD;
        expend = frgbgn + maxD;

        if ((frgend < regions[i].bgn) &&
            (regions[i].end < expbgn))
          expspan = true;

        if (tigLen < expend)
          expout = true;

      } else {
        expbgn = frgend - maxD;
        expend = frgend - minD;

        if ((expend < regions[i].bgn) &&
            (regions[i].end < frgbgn))
          expspan = true;

        if (frgend < maxD)
          expout = true;
      }

      if (expspan == false)
        //  Pair doesn't span the repeat.  Don't care.  Move along.
        continue;

      if (expout == true)
        //  Mate can be placed outside the unitig.  Don't care.  Move along.
        continue;

      if (target->id() != mrgtig) {
        //  Mate is not in this unitig, but it should be.  That's bad.
        spanBad++;

        //writeLog("frg "F_U32","F_U32",%d  region "F_U32","F_U32"  expected "F_U32","F_U32"  mate in tig "F_U32"  BAD\n",
        //        frgbgn, frgend, frg53, regions[i].bgn, regions[i].end, expbgn, expend, mrgtig);

        continue;
      }

      //  Otherwise, the mate is in this unitig, and it should be spanning the repeat.
      //  If it is in the expected location and oriented correctly, that's good!
      //  If it's almost in the expected location, that's maybe good.

      mrg      = &target->ufpath[ target->pathPosition(mid) ];
      mrgbgn   = (mrg->position.bgn < mrg->position.end) ? mrg->position.bgn : mrg->position.end;
      mrgend   = (mrg->position.bgn < mrg->position.end) ? mrg->position.end : mrg->position.bgn;
      mrg53    = (mrg->position.bgn < mrg->position.end) ? true : false;

      if ((expbgn <= mrgbgn) &&
          (mrgend <= expend) &&
          (frg53  != mrg53)) {
        spanGood++;
        //writeLog("frg "F_U32","F_U32",%d  region "F_U32","F_U32"   expected "F_U32","F_U32"  mate "F_U32","F_U32",%d  GOOD\n",
        //        frgbgn, frgend, frg53, regions[i].bgn, regions[i].end, expbgn, expend, mrgbgn, mrgend, mrg53);

      } else if ((expbgn              <= mrgbgn + 2 * stddev) &&
                 (mrgend + 2 * stddev <= expend) &&
                 (frg53  != mrg53)) {
        spanGoodMaybe++;

      } else {
        spanBad++;
        //writeLog("frg "F_U32","F_U32",%d  region "F_U32","F_U32"   expected "F_U32","F_U32"  mate "F_U32","F_U32",%d  BAD\n",
        //        frgbgn, frgend, frg53, regions[i].bgn, regions[i].end, expbgn, expend, mrgbgn, mrgend, mrg53);
      }
    }

    //  If we're mostly good/goodMaybe and very few bad, mark this is a repeat that we should
    //  kick out unmated and externally mated reads from.
    //
    //  We're allowing 

    uint32  totalSpan = spanBad + spanGood + spanGoodMaybe;

    double  sb = (double)spanBad       / totalSpan;
    double  sg = (double)spanGood      / totalSpan;
    double  sm = (double)spanGoodMaybe / totalSpan;

    if ((sb < 0.05) &&
        (sm < 0.10) &&
        (10 < spanGood))
      regions[i].ejectUnanchored = true;

    writeLog("markRepeats()-- region["F_U32"] "F_U32","F_U32" -- bad "F_U32" good "F_U32" goodmaybe "F_U32" -- sb=%.4f sg=%.4f sm-%.4f%s\n",
            i,
            regions[i].bgn, regions[i].end,
            spanBad, spanGood, spanGoodMaybe,
            sb, sg, sm,
            regions[i].ejectUnanchored ? " -- DON'T SPLIT" : "");
  }
}




//  Unitig fragment is completely within the repeat interval, or is close enough to the edge that
//  maybe it couldn't be placed uniquely.
//
void
markRepeats_findFragsInRegions(Unitig                    *target,
                               vector<repeatRegion>      &regions,
                               set<AS_IID>               &rptFrags,
                               set<AS_IID>               &ejtFrags) {

  for (uint32 i=0; i<regions.size(); i++) {
    for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
      ufNode     *frg   = &target->ufpath[fi];
      uint32      bgn   = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
      uint32      end   = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

      if ((bgn + REPEAT_COVERS_FRAG < regions[i].bgn) ||
          (regions[i].end + REPEAT_COVERS_FRAG < end))
        continue;

      if (regions[i].ejectUnanchored == false) {
        rptFrags.insert(frg->ident);
        continue;
      }

      AS_IID  mid = FI->mateIID(frg->ident);
      uint32  min = target->fragIn(mid);

      if ((mid == 0) || (min != target->id()))
        ejtFrags.insert(frg->ident);
    }
  }
}






//  Find any junctions in a region, and append them to list of fragments to split on.
//  This takes multiple lines of evidence (pointing to the same fragment end) and
//  combines them into a list of break points.
void
markRepeats_filterJunctions(Unitig                          *target,
                            vector<repeatRegion>            &regions,
                            vector<repeatJunctionEvidence>  &evidence,
                            vector<repeatUniqueBreakPoint>  &breakpoints) {

  map<FragmentEnd,uint32>   brkFrags;

  for (uint32 bi=0; bi<regions.size(); bi++)
    writeLog("markRepeats()--  region["F_U32"] "F_U32","F_U32" length "F_U32"\n",
            bi, regions[bi].bgn, regions[bi].end, regions[bi].end - regions[bi].bgn);

  sort(evidence.begin(), evidence.end());

  for (uint32 ai=0, bi=0; ai<evidence.size(); ai++) {
    repeatUniqueBreakPoint  ruj;

    //if (brkFrags[evidence[ai].tigFrag] >= ISECT_NEEDED_TO_BREAK)
    //  //  We've already logged a break on the fragment pointed to by the evidence.
    //  continue;

    if (evidence[ai].is3 == false) {
      assert(evidence[ai].uncovered5bgn < evidence[ai].uncovered5end);
      ruj = repeatUniqueBreakPoint(evidence[ai].point, evidence[ai].tigFrag, false);
    } else {
      assert(evidence[ai].uncovered3bgn < evidence[ai].uncovered3end);
      ruj = repeatUniqueBreakPoint(evidence[ai].point, evidence[ai].tigFrag, true);
    }

    //  Try to associate this junction with one of the repeat regions.  If there is no region,
    //  this is NOT a junction we care about.  The region must have been contained in a fragment.

    while ((bi < regions.size()) &&
           (regions[bi].end < ruj.point))
      //  Advance the region until it ends after the point.
      bi++;

    //  If this point is in the region, the region bgn will be lower (or equal) than the point.  We
    //  already ensured that the region end is after the point.

    if ((bi >= regions.size()) ||
        (ruj.point < regions[bi].bgn)) {
#if 1
      writeLog("markRepeats()--  junction: %s %6d %s at %d/%c' - DISCARD not in region\n",
              ruj.rptLeft ? "repeat" : "unique", ruj.point, ruj.rptLeft ? "unique" : "repeat",
              ruj.breakFrag.fragId(), ruj.breakFrag.frag3p() ? '3' : '5');
#endif
      continue;
    }

    //  If this region was marked as confirmed by mates, don't add the aplit points.

    if (regions[bi].ejectUnanchored == true) {
#if 1
      writeLog("markRepeats()--  junction: %s %6d %s at %d/%c' - DISCARD not in a breakable region\n",
              ruj.rptLeft ? "repeat" : "unique", ruj.point, ruj.rptLeft ? "unique" : "repeat",
              ruj.breakFrag.fragId(), ruj.breakFrag.frag3p() ? '3' : '5');
#endif
      continue;
    }

    assert(regions[bi].bgn <= ruj.point);
    assert(ruj.point       <= regions[bi].end);

    //  A new valid break point.

    brkFrags[evidence[ai].tigFrag]++;

    writeLog("markRepeats()-- junction: %s %6d %s at %d/%c' - with "F_U32" intersections\n",
            ruj.rptLeft ? "repeat" : "unique", ruj.point, ruj.rptLeft ? "unique" : "repeat",
            ruj.breakFrag.fragId(), ruj.breakFrag.frag3p() ? '3' : '5',
            brkFrags[evidence[ai].tigFrag]);

    if (brkFrags[evidence[ai].tigFrag] == 5)
      //if (brkFrags[evidence[ai].tigFrag] >= ISECT_NEEDED_TO_BREAK)
      breakpoints.push_back(ruj);
  }

  sort(breakpoints.begin(), breakpoints.end());

  writeLog("markRepeats()--  unitig %d has "F_SIZE_T" interesting junctions at the following regions:\n",
          target->id(), breakpoints.size());

  for (uint32 ji=0; ji<breakpoints.size(); ji++)
    writeLog("markRepeats()--  junction["F_U32"] at "F_IID"/%c' position "F_U32" repeat %s (count "F_U32"\n",
            ji,
            breakpoints[ji].breakFrag.fragId(), breakpoints[ji].breakFrag.frag3p() ? '3' : '5',
            breakpoints[ji].point,
            breakpoints[ji].rptLeft ? "<-" : "->",
            brkFrags[breakpoints[ji].breakFrag]);
  for (uint32 bi=0; bi<regions.size(); bi++)
    writeLog("markRepeats()--  region["F_U32"] "F_U32","F_U32" distToBreakPoint "F_U32","F_U32"\n",
            bi, regions[bi].bgn, regions[bi].end, regions[bi].bgnDist, regions[bi].endDist);
}




void
markRepeats_breakUnitigs(UnitigVector                    &unitigs,
                         Unitig                          *target,
                         vector<overlapPlacement>        &places,
                         vector<repeatUniqueBreakPoint>  &breakpoints,
                         set<AS_IID>                     &jctFrags,
                         set<AS_IID>                     &ejtFrags) {

  jctFrags.clear();

  if (breakpoints.size() == 0)
    return;

  uint32   *breakID = new uint32 [target->ufpath.size()];

  uint32    nextBreakPoint = 0;
  bool      currIsRepeat = (breakpoints[nextBreakPoint].rptLeft == true);
  uint32    curr = 1;
  uint32    next = 2;

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode     *frg   = &target->ufpath[fi];
    uint32      bgn   = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
    uint32      end   = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

    if (nextBreakPoint >= breakpoints.size()) {
      breakID[fi] = curr;

    } else if (breakpoints[nextBreakPoint].rptLeft == false) {
      //  Repeat to the right.  If the fragment starts at or after the junction, place this and
      //  future fragments into a new (repeat) unitig.

      if (breakpoints[nextBreakPoint].point <= bgn) {
        nextBreakPoint++;
        curr++;
        next++;
        jctFrags.insert(frg->ident);
      }

      //  If the fragment ends after the next junction, this fragment goes to the next
      //  unitig.  Otherwise, to the current one.
      //  
      if ((nextBreakPoint < breakpoints.size()) &&
          (breakpoints[nextBreakPoint].point < end) &&
          (breakpoints[nextBreakPoint].rptLeft == true))
        breakID[fi] = next;
      else
        breakID[fi] = curr;

    } else {
      //  Repeat to the left.  If the fragment ends before the junction, move to the current unitig

      if (end < breakpoints[nextBreakPoint].point) {
        breakID[fi] = curr;
        jctFrags.insert(frg->ident);
      } else {
        breakID[fi] = next;
      }

      //  Once we pass the junction, update pointers.  We are out of the repeat interval now.
      if (breakpoints[nextBreakPoint].point < bgn) {
        nextBreakPoint++;
        curr++;
        next++;
      }
    }
  }

  //  Append new unitigs.

  vector<Unitig *>      newTigs;
  Unitig              **uidToUnitig = new Unitig * [next + 1];
  uint32               *uidToOffset = new uint32   [next + 1];

  memset(uidToUnitig, 0, sizeof(Unitig *) * (next + 1));
  memset(uidToOffset, 0, sizeof(uint32)   * (next + 1));

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode  frg = target->ufpath[fi];
    uint32  bid = breakID[fi];

    if (ejtFrags.count(frg.ident) > 0) {
      target->removeFrag(frg.ident);
      continue;
    }

    if (uidToUnitig[bid] == NULL) {
      uidToUnitig[bid] = unitigs.newUnitig(false);
      uidToOffset[bid] = -MIN(frg.position.bgn, frg.position.end);

      newTigs.push_back(uidToUnitig[bid]);  //  For reporting below.
    }

    uidToUnitig[bid]->addFrag(frg, uidToOffset[bid], false);
  }

  delete [] breakID;
  delete [] uidToUnitig;
  delete [] uidToOffset;


  if (newTigs.size() > 0) {
    writeLog("markRepeats()-- SPLIT unitig %d of length %u with %ld fragments into "F_SIZE_T" unitigs:\n",
            target->id(), target->getLength(), target->ufpath.size(),
            newTigs.size());
    for (uint32 i=0; i<newTigs.size(); i++)
      writeLog("markRepeats()--   unitig %u of length %u with %ld fragments.\n",
              newTigs[i]->id(),
              newTigs[i]->getLength(),
              newTigs[i]->ufpath.size());

    unitigs[target->id()] = NULL;
    delete target;
  }

  //  Run back over the ejected frags, and place them with either their mate, or at their best location.

  for (set<AS_IID>::iterator it=ejtFrags.begin(); it!=ejtFrags.end(); it++)
    placeFragInBestLocation(unitigs, *it);
}



void
markRepeats_shatterRepeats(UnitigVector   &unitigs,
                           set<AS_IID>    &jctFrags,
                           set<AS_IID>    &covFrags) {

  for (set<AS_IID>::iterator it=jctFrags.begin(); it!=jctFrags.end(); it++) {
    AS_IID   iid = *it;
    uint32   ti  = Unitig::fragIn(iid);
    Unitig  *rpt = unitigs[ti];

    if ((ti == 0) || (rpt == NULL))
      //  Already shattered?
      continue;

    writeLog("markRepeats()--  shatter unitig %u with "F_SIZE_T" fragments from repeat frag %u\n",
            ti, rpt->ufpath.size(), iid);

    for (uint32 fi=0; fi<rpt->ufpath.size(); fi++)
      rpt->removeFrag(rpt->ufpath[fi].ident);

    unitigs[ti] = NULL;
    delete rpt;
  }

  for (set<AS_IID>::iterator it=covFrags.begin(); it!=covFrags.end(); it++) {
    AS_IID   iid = *it;
    uint32   ti  = Unitig::fragIn(iid);
    Unitig  *rpt = unitigs[ti];

    if ((ti == 0) || (rpt == NULL))
      //  Already shattered?
      continue;

    writeLog("markRepeats()--  frag "F_IID" covered by repeats, but still in unitig %d\n",
            iid, ti);
  }
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
            Unitig *target,
            bool shatterRepeats) {

  set<AS_IID>                     ovlFrags;
  set<AS_IID>                     covFrags;  //  Frag IIDs of fragments covered by repeat alignments
  set<AS_IID>                     jctFrags;  //  Frag IIDs of the first/last fragment in a repeat unitig
  set<AS_IID>                     ejtFrags;  //  Frag IIDs of frags we should eject instead of split

  double                          meanError = 0;
  double                          stddevError = 0;

  intervalList                    aligned;

  vector<overlapPlacement>        places;

  vector<repeatRegion>            regions;
  vector<repeatJunctionEvidence>  evidence;

  vector<repeatUniqueBreakPoint>  breakpoints;

  //  Build a list of all the fragments that have overlaps to this unitig.
  markRepeats_buildOverlapList(target, ovlFrags);

  //  Decide what a decent alignment should look like.
  markRepeats_computeUnitigErrorRate(unitigs, target, meanError, stddevError);

  //  For each overlapping fragment, place it and process.
  markRepeats_placeAndProcessOverlaps(unitigs, target, meanError, stddevError, ovlFrags, aligned, evidence);

  //  Convert 'aligned' into regions, throwing out weak ones and those contained in a fragment.
  markRepeats_filterIntervalsSpannedByFragment(target, aligned, regions);
  markRepeats_filterIntervalsSpannedByMates(target, regions);

  markRepeats_findFragsInRegions(target, regions, covFrags, ejtFrags);

  //  Discard junctions that are not in a remaining region.
  markRepeats_filterJunctions(target, regions, evidence, breakpoints);

  //  Split at whatever junctions remain.
  markRepeats_breakUnitigs(unitigs, target, places, breakpoints, jctFrags, ejtFrags);

  //  For each repeat unitig, shatter into fragments (not singleton unitigs) so we can later re-BOG.
  if (shatterRepeats)
    markRepeats_shatterRepeats(unitigs, jctFrags, covFrags);
}




void
markChimera(UnitigVector &unitigs,
            Unitig *target) {
}



void
mergeSplitJoin(UnitigVector &unitigs, const char *prefix, bool shatterRepeats) {

  //logFileFlags |= LOG_PLACE_FRAG;
  //logFileFlags &= ~LOG_PLACE_FRAG;

  if (IS == NULL)
    IS = new InsertSizes(unitigs);

  //  BUILD A LIST OF ALL INTERSECTIONS - build a reverse mapping of all BestEdges that are between
  //  unitigs.  For each fragment, we want to have a list of the incoming edges from other unitigs.

  intersectionList  *ilist = new intersectionList(unitigs);

  //ilist->logIntersections();

#if 0
  {
    Unitig        *target = unitigs[5];

    writeLog("popBubbles()-- WORKING on unitig %d/"F_SIZE_T" with %ld fragments.\n",
            target->id(), unitigs.size(), target->ufpath.size());

    mergeBubbles(unitigs, target, ilist);
    stealBubbles(unitigs, target, ilist);
    markRepeats(unitigs, target);
    markChimera(unitigs, target);
    exit(1);
  }
#endif

  //  Bubble popping
  //
  //  This used to be done right before each unitig was examined for repeats.
  //  It cannot be done in parallel -- there is a race condition when both unitigs
  //  A and B are considering merging in unitig C.

  setLogFile(prefix, "popBubbles");
  writeLog("popBubbles()-- working on "F_U64" unitigs.\n", unitigs.size());

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig        *target = unitigs[ti];

    if ((target == NULL) ||
        (target->ufpath.size() < 15) ||
        (target->getLength() < 300))
      continue;

    writeLog("popBubbles()-- WORKING on unitig %d/"F_SIZE_T" of length %u with %ld fragments.\n",
            target->id(), unitigs.size(), target->getLength(), target->ufpath.size());

    mergeBubbles(unitigs, target, ilist);
    stealBubbles(unitigs, target, ilist);
  }

  reportOverlapsUsed(unitigs, prefix, "popBubbles");
  reportUnitigs(unitigs, prefix, "popBubbles");
  evaluateMates(unitigs, prefix, "popBubbles");

  //  Since we create new unitigs for any of the splits, we need to remember
  //  where to stop.  We don't want to re-examine any of the split unitigs.
  //  Reevaluating seems to just trim off a few fragments at the end of the unitig.

  uint32  tiLimit = unitigs.size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  setLogFile(prefix, "mergeSplitJoin");
  writeLog("repeatDetect()-- working on "F_U32" unitigs, with "F_U32" threads.\n", tiLimit, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig        *target = unitigs[ti];

    if ((target == NULL) ||
        (target->ufpath.size() < 15) ||
        (target->getLength() < 300))
      continue;

    writeLog("repeatDetect()-- WORKING on unitig %d/"F_SIZE_T" of length %u with %ld fragments.\n",
            target->id(), unitigs.size(), target->getLength(), target->ufpath.size());

    markRepeats(unitigs, target, shatterRepeats);
    markChimera(unitigs, target);
  }

  reportOverlapsUsed(unitigs, prefix, "mergeSplitJoin");
  reportUnitigs(unitigs, prefix, "mergeSplitJoin");
  evaluateMates(unitigs, prefix, "mergeSplitJoin");

  //  JOIN EXPOSED BEST - after bubbles are stolen, this should leave some unitigs
  //  with exposed best edges that can now be connected.

  //  do we need to re-mark repeats after joining?

  //  SPLIT MARKED REPEATS - 

  //  SPLIT MARKED CHIMERA - 

  //  MERGE LEFTOVERS - these are the leftover pieces after repeats/chimera are split.  Hopefully
  //  they'll just be low coverage spurs

  delete ilist;

  delete IS;
  IS = NULL;

  logFileFlags &= ~LOG_PLACE_FRAG;
}

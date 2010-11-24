
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

static const char *rcsid = "$Id: AS_BAT_IntersectBubble.C,v 1.1 2010-11-24 01:03:31 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_UnitigGraph.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "MultiAlignStore.h"

#define MAX_OVERLAPS_PER_FRAG   (16 * 1024 * 1024)

//  IN THE FIRST PART:  Test the size of the bubble in the larger unitig using edges at the end.
//
//  Function is immediately true if the is only one edge from the bubble to the large unitig.
//  Otherwise, continue.
//
//  Function is immediately false if the edges are to different larger unitigs.  Function
//  fails if the edges are to the same unitig, but that unitig is not 'larger'.
//
//                ---bubble---
//               |            |
//               v            v
//     -----------------------------------
//
//  Function is immediately false if the sizes are wildly different, if the edges are to different
//  unitigs, etc.  Otherwise, continue.
//
//  IN THE SECOND PART: Test the orientation of the fragments at the ends of the bubble.  The FIRST
//  PART is unlikely to succeed if the SECOND PART is false.  The order and orientation of these two
//  fragments should be the same when placed in the larger unitig.  The bubble can be merged into the
//  larger unitig in the same order&orientation it is currently represented as, or the entire bubble
//  unitig can be reverse complemented - in which case the order reversed and the orientations flip.
//
//  lFrg and rFrg are UPDATED with the position of those fragments in the larger unitig.
//
bool
UnitigGraph::validateBubbleWithEdges(Unitig *bubble,
                                     ufNode &lFrg, BestEdgeOverlap *lEnd,
                                     ufNode &rFrg, BestEdgeOverlap *rEnd,
                                     Unitig *larger,
                                     OverlapStore *ovlStoreUniq,
                                     OverlapStore *ovlStoreRept) {

  //  Compute placement of the two fragments.  Compare the size against the bubble.

  ufNode lFrgN = lFrg;
  ufNode rFrgN = rFrg;

#if 0
  int32 lFrg5idx = -1;
  int32 lFrg3idx = -1;

  bool lPlaced = larger->placeFrag(lFrgN, lFrg5idx, (isReverse(lFrg.position) == true)  ? NULL : lEnd,
                                   lFrgN, lFrg3idx, (isReverse(lFrg.position) == false) ? lEnd : NULL);

  int32 rFrg5idx = -1;
  int32 rFrg3idx = -1;

  bool rPlaced = larger->placeFrag(rFrgN, rFrg5idx, (isReverse(rFrg.position) == true)  ? NULL : rEnd,
                                   rFrgN, rFrg3idx, (isReverse(rFrg.position) == false) ? rEnd : NULL);

  if (lPlaced == false) {
    //  Huh?  Didn't place?  Emit diagnostics.
    fprintf(logFile, "popBubbles()-- Failed to place lFrg.\n");
    return(false);
  }

  if (rPlaced == false) {
    //  Huh?  Didn't place?  Emit diagnostics.
    fprintf(logFile, "popBubbles()-- Failed to place rFrg.\n");
    return(false);
  }

#else

  overlapPlacement    lFrgPlacement;
  overlapPlacement    rFrgPlacement;

  lFrgPlacement.errors  = 4.0e9;
  lFrgPlacement.aligned = 1;

  rFrgPlacement.errors  = 4.0e9;
  rFrgPlacement.aligned = 1;

  vector<overlapPlacement>   placements;

  placements.clear();

  UG->placeFragUsingOverlaps(lFrg.ident, ovlStoreUniq, ovlStoreRept, placements);
  for (uint32 i=0; i<placements.size(); i++) {
    if (placements[i].tigID != larger->id())
      continue;

    if (placements[i].fCoverage < 0.99)
      continue;

    if ((placements[i].nForward > 0) &&
        (placements[i].nReverse > 0))
      continue;

    if (placements[i].errors / placements[i].aligned < lFrgPlacement.errors / lFrgPlacement.aligned)
      lFrgPlacement = placements[i];
  }

  placements.clear();

  UG->placeFragUsingOverlaps(rFrg.ident, ovlStoreUniq, ovlStoreRept, placements);
  for (uint32 i=0; i<placements.size(); i++) {
    if (placements[i].tigID != larger->id())
      continue;

    if (placements[i].fCoverage < 0.99)
      continue;

    if ((placements[i].nForward > 0) &&
        (placements[i].nReverse > 0))
      continue;

    if (placements[i].errors / placements[i].aligned < rFrgPlacement.errors / rFrgPlacement.aligned)
      rFrgPlacement = placements[i];
  }

  lFrgN.ident             = lFrgPlacement.frgID;
  lFrgN.contained         = 0;
  lFrgN.parent            = 0;
  lFrgN.ahang             = 0;
  lFrgN.bhang             = 0;
  lFrgN.position          = lFrgPlacement.position;
  lFrgN.containment_depth = 0;

  rFrgN.ident             = rFrgPlacement.frgID;
  rFrgN.contained         = 0;
  rFrgN.parent            = 0;
  rFrgN.ahang             = 0;
  rFrgN.bhang             = 0;
  rFrgN.position          = rFrgPlacement.position;
  rFrgN.containment_depth = 0;

  //fprintf(logFile, "lFrgN %d,%d rFrgN %d,%d\n",
  //        lFrgN.position.bgn, lFrgN.position.end,
  //        rFrgN.position.bgn, rFrgN.position.end);

  if ((lFrgN.position.bgn == 0) &&
      (lFrgN.position.end == 0)) {
    fprintf(logFile, "popBubbles()-- Failed to place lFrg.\n");
    return(false);
  }

  if ((rFrgN.position.bgn == 0) &&
      (rFrgN.position.end == 0)) {
    fprintf(logFile, "popBubbles()-- Failed to place rFrg.\n");
    return(false);
  }
#endif


  int32 minL = MIN(lFrg.position.bgn, lFrg.position.end);
  int32 maxL = MAX(lFrg.position.bgn, lFrg.position.end);

  int32 minR = MIN(rFrg.position.bgn, rFrg.position.end);
  int32 maxR = MAX(rFrg.position.bgn, rFrg.position.end);

  int32 placedLen = MAX(maxL, maxR) - MIN(minL, minR);

  if (2 * placedLen < bubble->getLength()) {
    //  Too short.
    fprintf(logFile, "popBubbles()-- Too short.  lFrg %d,%d rFrg %d,%d.  L %d,%d R %d,%d len %d\n",
            lFrg.position.bgn, lFrg.position.end,
            rFrg.position.bgn, rFrg.position.end,
            minL, maxL, minR, maxR, placedLen);
    return(false);
  }

  if (2 * bubble->getLength() < placedLen) {
    //  Too long.
    fprintf(logFile, "popBubbles()-- Too long.  lFrg %d,%d rFrg %d,%d.  L %d,%d R %d,%d len %d\n",
            lFrg.position.bgn, lFrg.position.end,
            rFrg.position.bgn, rFrg.position.end,
            minL, maxL, minR, maxR, placedLen);
    return(false);
  }

  ////////////////////
  //
  //  Check orientations
  //
  ////////////////////

  //  If lFrg and rFrg are the same fragment (bubble is one uncontained fragment) then we're done.

  if (lFrg.ident == rFrg.ident) {
    lFrg = lFrgN;
    rFrg = rFrgN;
    return(true);
  }

  //  Otherwise, check that the orientation and positioning of the before and after fragments is the
  //  same.

  bool   bL    = (isReverse(lFrg.position));
  bool   bR    = (isReverse(rFrg.position));
  bool   bOrd  = (MIN(lFrg.position.bgn, lFrg.position.end) < MIN(rFrg.position.bgn, rFrg.position.end));

  bool   nL    = (isReverse(lFrgN.position));
  bool   nR    = (isReverse(rFrgN.position));
  bool   nOrd  = (MIN(lFrgN.position.bgn, lFrgN.position.end) < MIN(rFrgN.position.bgn, rFrgN.position.end));

  if (((bL == nL) && (bR == nR) && (bOrd == nOrd)) ||
      ((bL != nL) && (bR != nR) && (bOrd != nOrd))) {
    //  Yup, looks good!
    lFrg = lFrgN;
    rFrg = rFrgN;
    return(true);
  }

  //  Nope, something got screwed up in alignment.

  fprintf(logFile, "popBubbles()-- Order/Orientation problem.  bL %d bR %d bOrd %d  nL %d nR %d nOrd %d\n",
          bL, bR, bOrd,
          nL, nR, nOrd);

  return(false);
}




//  False if any of the fragments in 'bubble' are not fully covered by overlaps to fragments in
//  'larger'.  Such uncovered fragments would indicate a large bubble -- large enough that we failed
//  to find an overlap -- and would cause problems in consensus.
//
//  False if any of the fragments in 'bubble' cannot be placed between lFrg and rFrg.  This would
//  indicate the bubble contains a significant rearrangement and would cause problems in consensus.
//
//  If the above tests pass, 'bubble' is inserted into 'larger' and 'bubble' is deleted.
//
bool
UnitigGraph::validateBubbleFragmentsWithOverlaps(Unitig *bubble,
                                                 ufNode &lFrg,
                                                 ufNode &rFrg,
                                                 Unitig *larger,
                                                 OverlapStore *ovlStoreUniq,
                                                 OverlapStore *ovlStoreRept) {

  //  Method:
  //
  //  * Call placeFragUsingOverlaps() for every fragment.  Save the placements returned.
  //  * Count the number of placements that are outside the lFrg/rFrg range.
  //  * Isolate down to one 'best' placement for each fragment.
  //    * Must be within lFrg/rFrg.
  //    * Resolve ties with
  //      * Placement in the original unitig       
  //      * Error rates on overlaps
  //      * Mate pairs

  bool success = false;

  vector<overlapPlacement>    *placements   = new vector<overlapPlacement> [bubble->ufpath.size()];
  overlapPlacement            *correctPlace = new        overlapPlacement  [bubble->ufpath.size()];

  for (uint32 fi=0; fi<bubble->ufpath.size(); fi++) {
    ufNode *frg = &bubble->ufpath[fi];

    UG->placeFragUsingOverlaps(frg->ident, ovlStoreUniq, ovlStoreRept, placements[fi]);

    //  Initialize the final placement to be bad, so we can pick the best.
    correctPlace[fi].fCoverage = 0.0;
    correctPlace[fi].errors    = 4.0e9;
    correctPlace[fi].aligned   = 1;
  }

  //  Some bizarre cases -- possibly only from bad data -- confound any logical attempt at finding the min/max extents.  Yes, even though
  //  this should work, it doesn't.  Or maybe it's just broken and I haven't seen how.
  //
  //int32  minE = (lFrg.position.bgn < rFrg.position.bgn) ? MIN(lFrg.position.bgn, lFrg.position.end) : MIN(rFrg.position.bgn, rFrg.position.end);
  //int32  maxE = (lFrg.position.bgn < rFrg.position.bgn) ? MAX(rFrg.position.bgn, rFrg.position.end) : MAX(lFrg.position.bgn, lFrg.position.end);
  //
  //  The one case that breaks it is a bubble unitig with a single chimeric fragment.
  //    lFrg ident = 367563, contained = 0, parent = 254673, ahang =  144,   bhang = 24,  bgn = 33406, end = 33238}
  //    rFrg ident = 367563, contained = 0, parent = 147697, ahang = -58,  bhang = -157,  bgn = 33406, end = 33574}
  //
  int32  minE = MIN(MIN(lFrg.position.bgn, lFrg.position.end), MIN(rFrg.position.bgn, rFrg.position.end));
  int32  maxE = MAX(MAX(lFrg.position.bgn, lFrg.position.end), MAX(rFrg.position.bgn, rFrg.position.end));
  int32  diff = maxE - minE;

  assert(minE < maxE);

  minE -= diff / 2;    if (minE < 0)  minE = 0;
  maxE += diff / 2;

  uint32  nCorrect = 0;

  for (uint32 fi=0; fi<bubble->ufpath.size(); fi++) {
    uint32  nNotPlaced = 0;
    uint32  nNotPlacedInLarger = 0;
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
      if (placements[fi][pl].tigID != larger->id()) {
        nNotPlacedInLarger++;
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

      if ((placements[fi][pl].nForward > 0) &&
          (placements[fi][pl].nReverse > 0)) {
        nNotOriented++;
        continue;
      }

      //  The current placement seems like a good one.  Is it the best one?
      //
      //  TODO:  We should also use mate pair information here!
      //
      if (placements[fi][pl].errors / correctPlace[fi].aligned < correctPlace[fi].errors / correctPlace[fi].aligned) {
        correctPlace[fi] = placements[fi][pl];
      }
    }  //  over all placements

    if (correctPlace[fi].fCoverage > 0)
      nCorrect++;
    else
      fprintf(logFile, "popBubbles()-- Failed to place frag %d notPlaced %d notPlacedInLarger %d notPlacedInCorrectPosition %d notPlacedFully %d notOriented %d\n",
              bubble->ufpath[fi].ident, nNotPlaced, nNotPlacedInLarger, nNotPlacedInCorrectPosition, nNotPlacedFully, nNotOriented);
  }

  if (nCorrect != bubble->ufpath.size())
    goto finished;

  //  Now just move the fragments into the larger unitig and delete the bubble unitig.
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

    larger->addFrag(nFrg, 0, logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG));
  }

  larger->sort();

  success = true;

  fprintf(logFile, "popBubbles()--  merged bubble unitig %d into unitig %d\n",
          bubble->id(), larger->id());

 finished:
  delete [] placements;
  delete [] correctPlace;

  return(success);
}






//  True if overlaps for every fragment in the bubble are colocated on the large unitig.  There
//  isn't a large insertion in the large unitig.
//
//bool
//UnitigGraph::validateBubblePlacementWithOverlaps(Unitig *bubble,
//                                                 Unitig *larger) {
//}


//  After all the above tests pass (the tests ALSO compute the placement) this
//  does the mechanics of inserting the short unitig into the large unitig.
//void
//UnitigGraph::insertBubble() {
//}




void
UnitigGraph::popIntersectionBubbles(OverlapStore *ovlStoreUniq, OverlapStore *ovlStoreRept) {
  uint32      ovlMax = MAX_OVERLAPS_PER_FRAG;
  uint32      ovlLen = 0;
  OVSoverlap *ovl    = new OVSoverlap [ovlMax];
  uint32     *ovlCnt = new uint32     [AS_READ_MAX_NORMAL_LEN];

  uint32      nBubblePopped    = 0;
  uint32      nBubbleTooBig    = 0;
  uint32      nBubbleConflict  = 0;
  uint32      nBubbleNoEdge    = 0;

  fprintf(logFile, "==> SEARCHING FOR BUBBLES\n");

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig        *shortTig = unitigs[ti];
    Unitig        *mergeTig = NULL;

    if (shortTig == NULL)
      //  Ain't no tig here!
      continue;

    //  Search for edges.  For a bubble to exist, either the first or last non-contained fragment
    //  must have an edge to the 'merge' unitig it is a bubble of.  Ideally, both the first and
    //  last will have edges to the same unitig, but we'll test and allow only a single edge.

    uint32  fIdx = 0;
    uint32  lIdx = shortTig->ufpath.size() - 1;

    //  Zombie fragments (circular containments) violate this constraint.  So we just ignore it for
    //  singleton unitigs, and completely skip searching for the last non-contained in that case.
    if (shortTig->ufpath.size() > 1) {
      assert(OG->isContained(shortTig->ufpath[fIdx].ident) == false);

      while (OG->isContained(shortTig->ufpath[lIdx].ident) == true)
        lIdx--;

      assert(lIdx >= 0);
    }

    ufNode  fFrg = shortTig->ufpath[fIdx];  //  NOTE:  A COPY, not a pointer or reference.
    ufNode  lFrg = shortTig->ufpath[lIdx];  //         These get modified.

    //  Grab the best edges outside the unitig.  If the first fragment is reversed, we want
    //  to grab the edge off of the 3' end; opposite for the last fragment.

    bool             f3p   = (isReverse(fFrg.position) == true);
    BestEdgeOverlap *fEdge = OG->getBestEdgeOverlap(fFrg.ident, f3p);

    bool             l3p   = (isReverse(lFrg.position) == false);
    BestEdgeOverlap *lEdge = OG->getBestEdgeOverlap(lFrg.ident, l3p);

    //  Just make sure...those edges should NOT to be to ourself.

    uint32  fUtg = Unitig::fragIn(fEdge->fragId());
    uint32  lUtg = Unitig::fragIn(lEdge->fragId());

    if ((fUtg == shortTig->id()) ||
        (lUtg == shortTig->id()))
      //  Must have been a circular unitig.  We'll split it later.
      continue;

    if ((fUtg == 0) &&
        (lUtg == 0))
      //  No edges, no bubble.
      continue;

    if ((fUtg != 0) && (lUtg != 0))
      fprintf(logFile, "popBubbles()-- Potential bubble unitig %d of length %d with %d fragments.  Edges (%d/%d') from frag %d/%d' and (%d/%d') from frag %d/%d'\n",
              shortTig->id(), shortTig->getLength(), shortTig->ufpath.size(),
              fEdge->fragId(), (fEdge->frag3p() ? 3 : 5), fFrg.ident, (f3p ? 3 : 5),
              lEdge->fragId(), (lEdge->frag3p() ? 3 : 5), lFrg.ident, (l3p ? 3 : 5));
    else if (fUtg != 0)
      fprintf(logFile, "popBubbles()-- Potential bubble unitig %d of length %d with %d fragments.  Edge (%d/%d') from frag %d/%d'\n",
              shortTig->id(), shortTig->getLength(), shortTig->ufpath.size(),
              fEdge->fragId(), (fEdge->frag3p() ? 3 : 5), fFrg.ident, (f3p ? 3 : 5));
    else if (lUtg != 0)
      fprintf(logFile, "popBubbles()-- Potential bubble unitig %d of length %d with %d fragments.  Edge (%d/%d') from frag %d/%d'\n",
              shortTig->id(), shortTig->getLength(), shortTig->ufpath.size(),
              lEdge->fragId(), (lEdge->frag3p() ? 3 : 5), lFrg.ident, (l3p ? 3 : 5));
    else
      assert(0);

    //  The only interesting case here is if we have both edges and they point to different unitigs.  We could try
    //  placing in both, but for now, we just give up.

    if ((fUtg != 0) && (lUtg != 0) && (fUtg != lUtg)) {
      fprintf(logFile, "popBubbles()-- bubble unitig %d has edges to both unitig %d and unitig %d, cannot place (yet)\n",
              shortTig->id(), fUtg, lUtg);
      continue;
    }

    mergeTig = (fUtg == 0) ? unitigs[lUtg] : unitigs[fUtg];

    if (validateBubbleWithEdges(shortTig, fFrg, fEdge, lFrg, lEdge, mergeTig, ovlStoreUniq, ovlStoreRept) == false) {
      fprintf(logFile, "popBubbles()-- failed to validate edges for bubble unitig %d into larger unitig %d\n",
              shortTig->id(), mergeTig->id());
      continue;
    }

    if (validateBubbleFragmentsWithOverlaps(shortTig, fFrg, lFrg, mergeTig, ovlStoreUniq, ovlStoreRept) == false) {
      fprintf(logFile, "popBubbles()-- failed to validate fragments for bubble unitig %d into larger unitig %d\n",
              shortTig->id(), mergeTig->id());
      continue;
    }

    //  Merged successfully!

    delete shortTig;
    unitigs[ti] = shortTig = NULL;
  }
}

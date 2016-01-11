
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
 *    src/AS_BAT/AS_BAT_IntersectBubble.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2015-JUN-03
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "AS_BAT_OverlapCache.H"

#include "AS_BAT_IntersectBubble.H"

//#include "AS_BAT_BestOverlapGraph.H"
//#include "MultiAlignStore.H"

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


////////////////////////////////////////
//
//  This is suboptimal.  It is possible to have two bubbles, A and B, that both really belong in
//  larger unitig C, but A -> B -> C.  If A is merged into B first, then there is a chance that
//  we'll lose the edges at the ends of B that place it in C (via ties; two fragments end at exaclty
//  the end of unitig B).
//
////////////////////////////////////////



static
bool
validateBubbleWithEdges(UnitigVector &unitigs,
                        double erateBubble,
                        Unitig *bubble,
                        ufNode &lFrg, BestEdgeOverlap *lEnd,
                        ufNode &rFrg, BestEdgeOverlap *rEnd,
                        Unitig *larger) {

  assert(0);

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
    writeLog("popBubbles()-- Failed to place lFrg.\n");
    return(false);
  }

  if (rPlaced == false) {
    //  Huh?  Didn't place?  Emit diagnostics.
    writeLog("popBubbles()-- Failed to place rFrg.\n");
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

  placeFragUsingOverlaps(unitigs, erateBubble, larger, lFrg.ident, placements);
  for (uint32 i=0; i<placements.size(); i++) {
    assert(placements[i].tigID == larger->id());
    if (placements[i].tigID != larger->id()) continue;

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
    writeLog("popBubbles()-- Failed to place lFrg.\n");
    return(false);
  }


  placements.clear();

  placeFragUsingOverlaps(unitigs, erateBubble, larger, rFrg.ident, placements);
  for (uint32 i=0; i<placements.size(); i++) {
    assert(placements[i].tigID == larger->id());
    if (placements[i].tigID != larger->id()) continue;

    if (placements[i].fCoverage < 0.99)
      continue;

    if (placements[i].errors / placements[i].aligned < rFrgPlacement.errors / rFrgPlacement.aligned)
      rFrgPlacement = placements[i];
  }

  rFrgN.ident             = rFrgPlacement.frgID;
  rFrgN.contained         = 0;
  rFrgN.parent            = 0;
  rFrgN.ahang             = 0;
  rFrgN.bhang             = 0;
  rFrgN.position          = rFrgPlacement.position;
  rFrgN.containment_depth = 0;

  if ((rFrgN.position.bgn == 0) &&
      (rFrgN.position.end == 0)) {
    writeLog("popBubbles()-- Failed to place rFrg.\n");
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
    writeLog("popBubbles()-- Too short.  lFrg %d,%d rFrg %d,%d.  L %d,%d R %d,%d len %d\n",
            lFrg.position.bgn, lFrg.position.end,
            rFrg.position.bgn, rFrg.position.end,
            minL, maxL, minR, maxR, placedLen);
    return(false);
  }

  if (2 * bubble->getLength() < placedLen) {
    //  Too long.
    writeLog("popBubbles()-- Too long.  lFrg %d,%d rFrg %d,%d.  L %d,%d R %d,%d len %d\n",
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

  writeLog("popBubbles()-- Order/Orientation problem.  bL %d bR %d bOrd %d  nL %d nR %d nOrd %d\n",
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
static
bool
validateBubbleFragmentsWithOverlaps(UnitigVector &unitigs,
                                    double erateBubble,
                                    Unitig *bubble,
                                    ufNode &lFrg,
                                    ufNode &rFrg,
                                    Unitig *larger) {

  assert(0);

  //  Method:
  //
  //  * Call placeFragUsingOverlaps() for every fragment.  Save the placements returned.
  //  * Count the number of placements that are outside the lFrg/rFrg range.
  //  * Isolate down to one 'best' placement for each fragment.
  //    * Must be within lFrg/rFrg.
  //    * Resolve ties with
  //      * Placement in the original unitig
  //      * Error rates on overlaps

  bool success = false;

  vector<overlapPlacement>    *placements   = new vector<overlapPlacement> [bubble->ufpath.size()];
  overlapPlacement            *correctPlace = new        overlapPlacement  [bubble->ufpath.size()];

  for (uint32 fi=0; fi<bubble->ufpath.size(); fi++) {
    ufNode *frg = &bubble->ufpath[fi];

    placeFragUsingOverlaps(unitigs, erateBubble, larger, frg->ident, placements[fi]);

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
      assert(placements[fi][pl].tigID == larger->id());
      if (placements[fi][pl].tigID != larger->id()) continue;

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

    if (correctPlace[fi].fCoverage > 0)
      nCorrect++;
    else
      writeLog("popBubbles()-- Failed to place frag %d notPlaced %d notPlacedInCorrectPosition %d notPlacedFully %d notOriented %d\n",
              bubble->ufpath[fi].ident, nNotPlaced, nNotPlacedInCorrectPosition, nNotPlacedFully, nNotOriented);
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

  writeLog("popBubbles()--  merged bubble unitig %d into unitig %d\n",
          bubble->id(), larger->id());

 finished:
  delete [] placements;
  delete [] correctPlace;

  return(success);
}



static
bool
popIntersectionBubble(UnitigVector &unitigs,
                      double        erateBubble,
                      Unitig       *shortTig) {

  //  Search for edges.  For a bubble to exist, either the first or last non-contained fragment
  //  must have an edge to the 'merge' unitig it is a bubble of.  Ideally, both the first and
  //  last will have edges to the same unitig, but we'll test and allow only a single edge.

  uint32  fIdx = 0;
  uint32  lIdx = shortTig->ufpath.size() - 1;

  //  We'd like to claim that all unitigs begin with a non-contained fragment, but zombie fragments
  //  (contained fragments that are in a circular containment relationship) violate this.  So, we
  //  could then claim that unitigs with more than one fragment begin with a non-contained fragment.
  //  But any zombie that has a bubble popped into it violate this.
  //
  //  We hope that any unitig that doesn't start with a non-contained fragment won't be a bubble, in
  //  particular, that there won't be non-contained fragments somewhere in that unitig.

  if (OG->isContained(shortTig->ufpath[fIdx].ident) == true) {
    writeLog("popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments STARTS WITH A CONTAINED FRAGMENT %d\n",
            shortTig->id(), shortTig->getLength(), shortTig->ufpath.size(),
            shortTig->ufpath[fIdx].ident);
    return(false);
  }

  //  Now, find the last non-contained fragment.

  while ((lIdx > 0) && (OG->isContained(shortTig->ufpath[lIdx].ident) == true))
    lIdx--;

  ufNode  fFrg = shortTig->ufpath[fIdx];  //  NOTE:  A COPY, not a pointer or reference.
  ufNode  lFrg = shortTig->ufpath[lIdx];  //         These get modified.

  //  Grab the best edges outside the unitig.  If the first fragment is reversed, we want
  //  to grab the edge off of the 3' end; opposite for the last fragment.

  bool             f3p   = (isReverse(fFrg.position) == true);
  BestEdgeOverlap *fEdge = OG->getBestEdgeOverlap(fFrg.ident, f3p);

  bool             l3p   = (isReverse(lFrg.position) == false);
  BestEdgeOverlap *lEdge = OG->getBestEdgeOverlap(lFrg.ident, l3p);

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
  //  The second b fragment now becomes the last non-contained fragment in the merged unitig, and it has
  //  a best edge to ourself.
  //
  //  We are no longer using fEdge or rEdge to place these fragments in the larger unitig; we're
  //  using all overlaps.  Just to be sure, we'll get rid of them.
  //
  uint32  fUtg = Unitig::fragIn(fEdge->fragId());
  uint32  lUtg = Unitig::fragIn(lEdge->fragId());

  if (fUtg == shortTig->id()) {
    fEdge = NULL;
    fUtg  = 0;
  }
  if (lUtg == shortTig->id()) {
    lEdge = NULL;
    lUtg  = 0;
  }

  if ((fUtg != 0) && (lUtg != 0))
    writeLog("popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  Edges (%d/%d') from frag %d/%d' and (%d/%d') from frag %d/%d'\n",
            shortTig->id(), shortTig->getLength(), shortTig->ufpath.size(),
            fEdge->fragId(), (fEdge->frag3p() ? 3 : 5), fFrg.ident, (f3p ? 3 : 5),
            lEdge->fragId(), (lEdge->frag3p() ? 3 : 5), lFrg.ident, (l3p ? 3 : 5));
  else if (fUtg != 0)
    writeLog("popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  Edge (%d/%d') from frag %d/%d'\n",
            shortTig->id(), shortTig->getLength(), shortTig->ufpath.size(),
            fEdge->fragId(), (fEdge->frag3p() ? 3 : 5), fFrg.ident, (f3p ? 3 : 5));
  else if (lUtg != 0)
    writeLog("popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  Edge (%d/%d') from frag %d/%d'\n",
            shortTig->id(), shortTig->getLength(), shortTig->ufpath.size(),
            lEdge->fragId(), (lEdge->frag3p() ? 3 : 5), lFrg.ident, (l3p ? 3 : 5));
  else {
    writeLog("popBubbles()-- Potential bubble unitig %d of length %d with %lu fragments.  NO EDGES, no bubble.\n",
            shortTig->id(), shortTig->getLength(), shortTig->ufpath.size());
    return(false);
  }

  //  The only interesting case here is if we have both edges and they point to different unitigs.  We could try
  //  placing in both, but for now, we just give up.

  if ((fUtg != 0) && (lUtg != 0) && (fUtg != lUtg)) {
    writeLog("popBubbles()-- bubble unitig %d has edges to both unitig %d and unitig %d, cannot place (yet)\n",
            shortTig->id(), fUtg, lUtg);
    return(true);
  }

  Unitig *mergeTig = (fUtg == 0) ? unitigs[lUtg] : unitigs[fUtg];

  if (validateBubbleWithEdges(unitigs, erateBubble, shortTig, fFrg, fEdge, lFrg, lEdge, mergeTig) == false) {
    writeLog("popBubbles()-- failed to validate edges for bubble unitig %d into larger unitig %d\n",
            shortTig->id(), mergeTig->id());
    return(false);
  }

  if (validateBubbleFragmentsWithOverlaps(unitigs, erateBubble, shortTig, fFrg, lFrg, mergeTig) == false) {
    writeLog("popBubbles()-- failed to validate fragments for bubble unitig %d into larger unitig %d\n",
            shortTig->id(), mergeTig->id());
    return(false);
  }

  //  Merged successfully!

  unitigs[shortTig->id()] = NULL;
  delete shortTig;

  return(true);
}


void
popIntersectionBubbles(UnitigVector &unitigs, double erateBubble) {
  uint32          nFrgToMerge      = 1;
  uint32          nFrgToMergeMax   = 500;

  uint32          nBubblePopped    = 0;

  logFileFlags |= LOG_PLACE_FRAG;

  while (1) {
    bool            keepPopping      = false;
    uint32          nBubbleFixed     = 1;
    vector<uint32>  tryAgain;

    //  Step 1:  Iterate over all possible merge sizes, popping whatever.

    for (nFrgToMerge=1; nFrgToMerge < nFrgToMergeMax; nFrgToMerge++) {
      writeLog("==> SEARCHING FOR BUBBLES of size %u fragments.\n", nFrgToMerge);

      for (uint32 ti=0; ti<unitigs.size(); ti++) {
        Unitig        *shortTig = unitigs[ti];

        if (shortTig == NULL)
          //  Ain't no tig here!
          continue;

        if (shortTig->ufpath.size() != nFrgToMerge)
          //  Wrong size.  We've either done it, or will do it, or it's just too big.
          continue;

        //  popIntersectionBubble() returns false if the shortTig cannot be merged.  It returns true
        //  if the shortTig was merged, or might be merged after some other merge.
        //
        if (popIntersectionBubble(unitigs, erateBubble, shortTig)) {
          if (unitigs[ti]) {
            tryAgain.push_back(ti);
          } else {
            nBubblePopped++;
            keepPopping = true;
          }
        }
      }  //  Over all unitigs
    }  //  Over all merge sizes

    //  Step 2:  If nothing changed, get out of here.

    if (keepPopping == false)
      break;

    //  Step 3: Attempt to merge bubbles that were across two unitigs.  Regardless of the order we
    //  merge in, it is possible for some bubble unitig B to have edges to A (a large unitig) and C
    //  (another bubble).  If B is examined before C, B will not be merged.  This is noted above
    //  ('tryAgain' will store B), and here we see if C was merged into A, thus allowing B to merge.

    while (nBubbleFixed > 0) {
      nBubbleFixed = 0;

      writeLog("==> SEARCHING FOR BUBBLES that spanned unitigs.\n");

      for (uint32 ta=0; ta<tryAgain.size(); ta++) {
        Unitig        *shortTig = unitigs[tryAgain[ta]];

        if (shortTig == NULL)
          continue;

        if (popIntersectionBubble(unitigs, erateBubble, shortTig)) {
          if (unitigs[tryAgain[ta]] == NULL) {
            nBubblePopped++;
            nBubbleFixed++;
            keepPopping = true;
          }
        }
      }
    }

    //  Step 4:  If nothing changed, get out of here.

    if (keepPopping == false)
      break;
  }  //  Until we break.

  logFileFlags &= ~LOG_PLACE_FRAG;

  writeLog("Popped %u bubbles.\n", nBubblePopped);
}

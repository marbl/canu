
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
 *    src/AS_BAT/AS_BAT_MergeSplitJoin.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2011-FEB-15 to 2014-MAY-03
 *      are Copyright 2011-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2015-AUG-05
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
#include "AS_BAT_IntersectSplit.H"

#include "AS_BAT_Instrumentation.H"

#include "intervalList.H"


#undef  LOG_BUBBLE_TESTS
#undef  LOG_BUBBLE_FAILURE
#define LOG_BUBBLE_SUCCESS


bool
mergeBubbles_findEnds(UnitigVector &UNUSED(unitigs),
                      double UNUSED(erateBubble),
                      Unitig *bubble,
                      ufNode &fFrg,
                      ufNode &lFrg,
                      Unitig *UNUSED(target)) {

  //  Search for edges.  For a bubble to exist, at least one of the first or last non-contained
  //  fragment must have an edge to the 'target' unitig (by construction of the inputs to this
  //  routine).  Ideally, both the first and last will have edges to the same unitig, but we'll test
  //  and allow only a single edge.

  uint32  zIdx = UINT32_MAX;
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
                       double erateBubble,
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

  placeFragUsingOverlaps(unitigs, erateBubble, target, fFrg.ident, placements);

#ifdef LOG_BUBBLE_TESTS
  writeLog("popBubbles()-- fFrg %u has %u potential placements in unitig %u.\n",
           fFrg.ident, placements.size(), target->id());
#endif

  for (uint32 i=0; i<placements.size(); i++) {
    assert(placements[i].tigID == target->id());

    if (placements[i].fCoverage < 0.99) {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- fFrg %u low coverage %f at unitig %u %u,%u\n",
               fFrg.ident,
               placements[i].fCoverage,
               placements[i].tigID,
               placements[i].position.bgn, placements[i].position.end);
#endif
      continue;
    } else {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- fFrg %u GOOD coverage %f at unitig %u %u,%u\n",
               fFrg.ident,
               placements[i].fCoverage,
               placements[i].tigID,
               placements[i].position.bgn, placements[i].position.end);
#endif
    }

    if (placements[i].errors / placements[i].aligned < fFrgPlacement.errors / fFrgPlacement.aligned) {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- fFrg %u GOOD identity %f at unitig %u %u,%u\n",
               fFrg.ident,
               placements[i].errors / placements[i].aligned,
               placements[i].tigID,
               placements[i].position.bgn, placements[i].position.end);
#endif
      fFrgPlacement = placements[i];
    } else {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- fFrg %u low identity %f at unitig %u %u,%u\n",
               fFrg.ident,
               placements[i].errors / placements[i].aligned,
               placements[i].tigID,
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

  placeFragUsingOverlaps(unitigs, erateBubble, target, lFrg.ident, placements);

#ifdef LOG_BUBBLE_TESTS
  writeLog("popBubbles()-- lFrg %u has %u potential placements in unitig %u.\n",
           lFrg.ident, placements.size(), target->id());
#endif

  for (uint32 i=0; i<placements.size(); i++) {
    assert(placements[i].tigID == target->id());

    if (placements[i].fCoverage < 0.99) {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- lFrg %u low coverage %f at unitig %u %u,%u\n",
               lFrg.ident,
               placements[i].fCoverage,
               placements[i].tigID,
               placements[i].position.bgn, placements[i].position.end);
#endif
      continue;
    } else {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- lFrg %u GOOD coverage %f at unitig %u %u,%u\n",
               lFrg.ident,
               placements[i].fCoverage,
               placements[i].tigID,
               placements[i].position.bgn, placements[i].position.end);
#endif
    }

    if (placements[i].errors / placements[i].aligned < lFrgPlacement.errors / lFrgPlacement.aligned) {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- lFrg %u GOOD identity %f at unitig %u %u,%u\n",
               lFrg.ident,
               placements[i].errors / placements[i].aligned,
               placements[i].tigID,
               placements[i].position.bgn, placements[i].position.end);
#endif
      lFrgPlacement = placements[i];
    } else {
#ifdef LOG_BUBBLE_FAILURE
      writeLog("popBubbles()-- lFrg %u low identity %f at unitig %u %u,%u\n",
               lFrg.ident,
               placements[i].errors / placements[i].aligned,
               placements[i].tigID,
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
                        double erateBubble,
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

  bool success = false;

  vector<overlapPlacement>    *placements   = new vector<overlapPlacement> [bubble->ufpath.size()];
  overlapPlacement            *correctPlace = new        overlapPlacement  [bubble->ufpath.size()];

  for (uint32 fi=0; fi<bubble->ufpath.size(); fi++) {
    ufNode *frg = &bubble->ufpath[fi];

    placeFragUsingOverlaps(unitigs, erateBubble, target, frg->ident, placements[fi]);

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

#ifdef LOG_BUBBLE_SUCCESS
  writeLog("popBubbles()--   merged bubble unitig %d with %ld frags into unitig %d now with %ld frags\n",
           bubble->id(), bubble->ufpath.size(), target->id(), target->ufpath.size());
#endif

 finished:
  delete [] placements;
  delete [] correctPlace;

  //  If not successful, mark this unitig as a potential bubble.

  if (success == false) {
    writeLog("popBubbles()--   bubble unitig %d (reads %u length %u) has large differences, not popped into unitig %d\n",
             bubble->id(), bubble->ufpath.size(), bubble->getLength(), target->id());
    bubble->_isBubble = true;
  }

  return(success);
}



void
mergeBubbles(UnitigVector &unitigs, double erateBubble, Unitig *target, intersectionList *ilist) {

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

      //  I don't like a number of reads filter - for 50x Illumina, 500 reads is only 1k of unitig,
      //  but for 10x PacBio, this is over 250k of unitig.

      if ((bubble == NULL) ||
          (bubble->getLength() > 50000)) {
#if 0
        writeLog("popBubbles()-- Skip bubble %u length %u with "F_SIZE_T" frags - edge from %d/%c' to utg %d %d/%c'\n",
                 bubble->id(), bubble->getLength(), bubble->ufpath.size(),
                 isect->invadFrg, isect->invad3p ? '3' : '5',
                 target->id(),
                 isect->isectFrg, isect->isect3p ? '3' : '5');
#endif
        continue;
      }

      if (bubble->id() == target->id())
        //  HEY!  We're not a bubble in ourself!
        continue;

      ufNode  fFrg;  //  First fragment in the bubble
      ufNode  lFrg;  //  Last fragment in the bubble

      //  We have no way of deciding if we've tested this bubble unitig already.  Each bubble unitig
      //  should generate two intersection edges.  If those edges are to the same target unitig, and
      //  the bubble fails to pop, we'll test the bubble twice.
      //
      //  This is kind of by design.  The two intersections could be to two different locations, and
      //  maybe one will work while the other doesn't.  Though, I think we accept a placement only
      //  if the two end reads are consistent implying that we'd double test a bubble if the
      //  placements are different, and that we'd fail both times.

      if (mergeBubbles_findEnds(unitigs, erateBubble, bubble, fFrg, lFrg, target) == false)
        continue;

      if (mergeBubbles_checkEnds(unitigs, erateBubble, bubble, fFrg, lFrg, target) == false)
        continue;

      if (mergeBubbles_checkFrags(unitigs, erateBubble, bubble, fFrg, lFrg, target) == false)
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




//  Bubble popping cannot be done in parallel -- there is a race condition when both unitigs
//  A and B are considering merging in unitig C.

void
popBubbles(UnitigVector &unitigs,
           double UNUSED(erateGraph), double erateBubble, double UNUSED(erateMerge), double erateRepeat,
           const char *prefix,
           uint32 minOverlap,
           uint64 genomeSize) {

  //logFileFlags |= LOG_PLACE_FRAG;
  //logFileFlags &= ~LOG_PLACE_FRAG;

  //  BUILD A LIST OF ALL INTERSECTIONS - build a reverse mapping of all BestEdges that are between
  //  unitigs.  For each fragment, we want to have a list of the incoming edges from other unitigs.

  intersectionList  *ilist = new intersectionList(unitigs);

  //ilist->logIntersections();

  writeLog("popBubbles()-- working on "F_U64" unitigs.\n", unitigs.size());

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig        *target = unitigs[ti];

    if ((target == NULL) ||
        (target->ufpath.size() < 1) ||  //  was 15
        (target->getLength() < 300))
      continue;

    //writeLog("popBubbles()-- WORKING on unitig %d/"F_SIZE_T" of length %u with %ld fragments.\n",
    //         target->id(), unitigs.size(), target->getLength(), target->ufpath.size());

    mergeBubbles(unitigs, erateBubble, target, ilist);
    //stealBubbles(unitigs, erateBubble, target, ilist);
  }

  delete ilist;
}

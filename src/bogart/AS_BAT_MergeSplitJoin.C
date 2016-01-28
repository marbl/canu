
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

#include "AS_BAT_MergeSplitJoin.H"

#include "AS_BAT_Breaking.H"

#include "AS_BAT_Instrumentation.H"

#include "AS_BAT_RepeatJunctionEvidence.H"

#include "intervalList.H"


uint32 SPURIOUS_COVERAGE_THRESHOLD  = 6;   //  Need to have more than this coverage in non-unitig reads aligned to call it a repeat area
uint32 ISECT_NEEDED_TO_BREAK        = 15;  //  Need to have at least  this number of reads confirming a repeat junction
uint32 REGION_END_WEIGHT            = 15;  //  Pretend there are this many intersections at the end points of each repeat region

#undef  LOG_BUBBLE_TESTS
#undef  LOG_BUBBLE_FAILURE
#define LOG_BUBBLE_SUCCESS

omp_lock_t  markRepeat_breakUnitigs_Lock;

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



//void
//stealBubbles(UnitigVector &unitigs, double erateBubble, Unitig *target, intersectionList *ilist) {
//}





void
markRepeats_buildOverlapList(Unitig *target, double erateRepeat, set<uint32> &ovlFrags) {

  ovlFrags.clear();

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode     *frg   = &target->ufpath[fi];
    uint32      ovlLen = 0;
    BAToverlap *ovl    = OC->getOverlaps(frg->ident, erateRepeat, ovlLen);

    for (uint32 i=0; i<ovlLen; i++) {
      if (Unitig::fragIn(ovl[i].b_iid) != target->id())
        ovlFrags.insert(ovl[i].b_iid);
    }
  }
}




void
markRepeats_computeUnitigErrorRate(UnitigVector &unitigs,
                                   double        erateRepeat,
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

    placeFragUsingOverlaps(unitigs, erateRepeat, target, frg->ident, op);

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
                                    double                            erateRepeat,
                                    Unitig                           *target,
                                    double                           meanError,
                                    double                           stddevError,
                                    set<uint32>                      &ovlFrags,
                                    intervalList<int32>              &aligned,
                                    vector<repeatJunctionEvidence>   &evidence) {

  aligned.clear();
  evidence.clear();

  for (set<uint32>::iterator it=ovlFrags.begin(); it!=ovlFrags.end(); it++) {
    uint32  iid = *it;

    vector<overlapPlacement>  op;

    placeFragUsingOverlaps(unitigs, erateRepeat, target, iid, op);

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


#if 0
uint32
markRepeats_computeUnitigCoverage(Unitig *tig) {
  intervalList<int32>   coverage;

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode  frg         = tig->ufpath[fi];

    if (frg.position.bgn < frg.position.end)
      coverage.add(frg.position.bgn, frg.position.end - frg.position.bgn);
    else
      coverage.add(frg.position.end, frg.position.bgn - frg.position.end);
  }

  coverage.merge();

  intervalDepth  depth(coverage);

  uint64   minCov = 0;           //  Probably 1 or 2
  uint64   medCov = 0;           //
  uint64   aveCov = 0;           //  If repeats present, might be high
  uint64   modCov = 0;           //  Probably not useful, due to read bias
  uint64   maxCov = UINT32_MAX;  //

  for (uint32 dd=0; dd<depth.numberOfIntervals(); dd++)
    maxCov = MAX(maxCov, depth.de(dd));
  maxCov++;

  uint32   *histogram = new uint32 [maxCov];

  memset(histogram, 0, sizeof(uint32) * maxCov);

  for (uint32 dd=0; dd<depth.numberOfIntervals(); dd++)
    histogram[depth.de(dd)] += depth.hi(dd) - depth.lo(dd);

  for (uint32 dd=0; dd<maxCov; dd++) {
    if (depth.de(dd) == 0)
      continue;

    if (minCov < depth.de(dd))
      minCov = depth.de(dd);

    if (depth.de(dd) < maxCov)
      maxCov = depth.de(dd);

    //  blah.
  }

  return(medCov);
}
#endif



//  Decide on a spurious coverage level using coverage in this unitig, global coverage, and
//  statistics from the potential repeats.
//
//  This is the end consumer of the 'aligned' data.  It is used only to populate 'regions', and
//  'regions' only cares about bgn,end coords, no underlying data.
//
void
markRepeats_filterIntervalsSpannedByFragment(Unitig                    *target,
                                             intervalList<int32>       &aligned,
                                             vector<repeatRegion>      &regions,
                                             uint32                     minOverlap) {
  uint32   tiglen  = target->getLength();

  uint32   spuriousNoiseThreshold = SPURIOUS_COVERAGE_THRESHOLD;
  uint32   filteredBases          = 0;
  uint32   filteredCovered        = 0;

  intervalList<int32>   depth(aligned);

  aligned.merge();  //  Just for a stupid log message

  writeLog("markRepeats()--  filtering low coverage spurious with t=%u in %u repeat regions (%u depth regions)\n",
           spuriousNoiseThreshold, aligned.numberOfIntervals(), depth.numberOfIntervals());

  regions.clear();
  aligned.clear();

  //
  //  Remove low coverage areas by making a new map for the high coverage areas.
  //

  for (uint32 dd=0; dd<depth.numberOfIntervals(); dd++) {
    if (depth.depth(dd) == 0)
      continue;

    if (depth.depth(dd) <= spuriousNoiseThreshold) {
      filteredBases += depth.hi(dd) - depth.lo(dd);
      continue;
    }

    aligned.add(depth.lo(dd), depth.hi(dd) - depth.lo(dd));
  }

  aligned.merge();

  writeLog("markRepeats()--  filtered %u bases, now with %u repeat regions\n",
           filteredBases, aligned.numberOfIntervals());

  //
  //  Adjust region boundaries so they land on the first read end that makes sense.
  //
  //  For the begin, decide if we need to expand or contract the region.  We will expand if the
  //  start of the read before the region is not anchored.  We will contract otherwise.
  //
  //  For the end, it is more complicated because the end points are not sorted.
  //

  for (uint32 i=0; i<aligned.numberOfIntervals(); i++) {
    uint32  intbgn  = aligned.lo(i);
    uint32  intend  = aligned.hi(i);

    for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
      ufNode     *frg    = &target->ufpath[fi];
      uint32      frgbgn = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
      uint32      frgend = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

      if (frgbgn + minOverlap / 2 < intbgn)
        //  Anchored.
        continue;

      if (frgbgn == intbgn) {
        //  Perfect!  Don't change a thing (unless we already expanded to get an unanchored read).
        break;
      }

      if ((frgbgn <= intbgn) && (intbgn <= frgbgn + minOverlap / 2)) {
        //  Not anchored, expand the region to this location
#ifdef VERBOSE_REGION_FITTING
        writeLog("markRepeats()--  region["F_U32"].bgn expands from "F_U32" to "F_U32" at frag "F_U32"\n", i, aligned.lo(i), frgbgn, frg->ident);
#endif
        aligned.lo(i) = frgbgn;
        break;
      }

      if (intbgn <= frgbgn) {
        //  First read begin that is inside the repeat region
#ifdef VERBOSE_REGION_FITTING
        writeLog("markRepeats()--  region["F_U32"].bgn contracts from "F_U32" to "F_U32" at frag "F_U32"\n", i, aligned.lo(i), frgbgn, frg->ident);
#endif
        aligned.lo(i) = frgbgn;
        break;
      }
    }

    uint32  newexp = 0, newexpid = UINT32_MAX;
    uint32  newcnt = 0, newcntid = UINT32_MAX;

    for (uint32 fi=target->ufpath.size(); fi-- > 0; ) {
      ufNode     *frg    = &target->ufpath[fi];
      uint32      frgbgn = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
      uint32      frgend = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

      if (intend + minOverlap / 2 < frgend)
        //  Anchored.
        continue;

      if (frgend == intend) {
        //  Perfect!
        //newexpid = UINT32_MAX;
        newcntid = UINT32_MAX;
        break;
      }

      if ((intend < frgend) && (frgend <= intend + minOverlap / 2) && (newexp < frgend)) {
        //  Not anchored, expand the region to this location (this will pick the largest expansion)
        newexp   = frgend;
        newexpid = frg->ident;
        continue;
      }

      if ((frgend <= intend) && (newcnt < frgend)) {
        //  Pick the largest read end that is within the repeat region
        newcnt   = frgend;
        newcntid = frg->ident;
        continue;
      }

      if (frgbgn + AS_MAX_READLEN < intend)
        //  All done, no more intersections possible
        break;
    }

    //  Expand the region if one was found, otherwise contract if one was found.

    if      (newexpid != UINT32_MAX) {
#ifdef VERBOSE_REGION_FITTING
      writeLog("markRepeats()--  region["F_U32"].end expands from "F_U32" to "F_U32" at frag "F_U32"\n", i, aligned.hi(i), newexp, newexpid);
#endif
      aligned.hi(i) = newexp;
    }

    else if (newcntid != UINT32_MAX) {
#ifdef VERBOSE_REGION_FITTING
      writeLog("markRepeats()--  region["F_U32"].end contracts from "F_U32" to "F_U32" at frag "F_U32"\n", i, aligned.hi(i), newcnt, newcntid);
#endif
      aligned.hi(i) = newcnt;
    }
  }


  {
    uint32 nc = 0;

    for (uint32 i=0; i<aligned.numberOfIntervals(); i++)
      if (aligned.hi(i) < aligned.lo(i))
        nc++;

    writeLog("markRepeats()--  filtered "F_U32" repeat regions after picking read endpoints, now with "F_U32" repeat regions.\n",
             nc, aligned.numberOfIntervals() - nc);
  }

  //
  //  Discard the interval if there is a fragment that contains it with enough overhang to
  //  unambiguously place it in a unitig.
  //

  for (uint32 i=0; i<aligned.numberOfIntervals(); i++) {
    ufNode  *containedIn  = NULL;

    //  If the region is backwards, then the region is contained in a read.  The easiest case to
    //  argue is:
    //
    //           ------------
    //                  -------RRR------
    //                              --------------
    //
    //  The read with the repeat region is anchored on both sides, so it is contracted to the next
    //  begin (for the start) and the previous end for the end)
    //
    if (aligned.hi(i) < aligned.lo(i)) {
      //writeLog("markRepeats()--  repeat alignment "F_U64","F_U64" DISCARD - CONTRACTED TO NULL\n",
      //         aligned.lo(i), aligned.hi(i));
      continue;
    }

    //  Ensure that reads not near the end of the unitig have enough non-repeat sequence to anchor the read in this location.
    //  This is done by increasing the size of the repeat.

    uint32   rptbgn  = aligned.lo(i);
    uint32   rptend  = aligned.hi(i);

    uint32   unique  = minOverlap / 2;

    bool     bgnFull = true;
    bool     endFull = true;

    if (unique <= rptbgn) {
      bgnFull = false;
      rptbgn -= unique;
    }

    if (rptend + unique <= tiglen) {
      endFull = false;
      rptend += unique;
    }

    //  Search for a covering fragment.

    for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
      ufNode     *frg    = &target->ufpath[fi];
      uint32      frgbgn = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
      uint32      frgend = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

      if (frgend < rptbgn)
        //  Fragment is before the region, keep searching for a spanning fragment.
        continue;

      if (rptend < frgbgn)
        //  Fragment is after the region, we're finished.
        break;

      if ((frgbgn <= rptbgn) &&
          (rptend <= frgend)) {
        //  Fragment contains the region with acceptable overhangs into the non-repeat area.
        containedIn = frg;
        break;
      }
    }

    if (containedIn != NULL) {
      filteredCovered++;
      //writeLog("markRepeats()--  repeat alignment %s"F_U64","F_U64"%s DISCARD - CONTAINED IN FRAGMENT "F_U32" "F_S32","F_S32"\n",
      //         (bgnFull ? "(end) " : ""), aligned.lo(i), aligned.hi(i), (endFull ? " (end)" : ""),
      //         containedIn->ident, containedIn->position.bgn, containedIn->position.end);
      continue;
    }

    if ((regions.size() == 0) ||
        (regions.back().end + 100 < aligned.lo(i))) {
      regions.push_back(repeatRegion(aligned.lo(i), aligned.hi(i)));
      //writeLog("markRepeats()--  repeat alignment "F_U32","F_U32"\n",
      //        regions.back().bgn, regions.back().end);

    } else {
      regions.back().bgn       = MIN(regions.back().bgn, aligned.lo(i));
      regions.back().end       = MAX(regions.back().end, aligned.hi(i));
      //writeLog("markRepeats()--  repeat alignment "F_U32","F_U32" (merged from "F_U64","F_U64")\n",
      //        regions.back().bgn, regions.back().end,
      //        aligned.lo(i), aligned.hi(i));
    }
  }

  writeLog("markRepeats()--  filtered %u repeat regions contained in a read, now with %u repeat regions\n",
           filteredCovered, regions.size());
}




//  Unitig fragment is completely within the repeat interval, or is close enough to the edge that
//  maybe it couldn't be placed uniquely.
//
void
markRepeats_findFragsInRegions(Unitig                    *target,
                               vector<repeatRegion>      &regions,
                               set<uint32>               &rptFrags,
                               set<uint32>               &UNUSED(ejtFrags),
                               uint32                     minOverlap) {

  for (uint32 i=0; i<regions.size(); i++) {
    for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
      ufNode     *frg   = &target->ufpath[fi];
      uint32      bgn   = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
      uint32      end   = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

      if (regions[i].bgn == bgn)
        regions[i].rujBgn = repeatUniqueBreakPoint(regions[i].bgn,
                                                   FragmentEnd(frg->ident, end < bgn),
                                                   false);

      if (regions[i].end == end)
        regions[i].rujEnd = repeatUniqueBreakPoint(regions[i].end,
                                                   FragmentEnd(frg->ident, bgn < end),
                                                   true);

      //  If the read has at least minOverlap/2 bases outside the repeat, assume it
      //  is placed correctly.
      if ((bgn + minOverlap/2 < regions[i].bgn) ||
          (regions[i].end + minOverlap/2 < end))
        continue;

      //  Otherwise, the read is 'contained' in a repeat region.  Remember it for later processing.
      rptFrags.insert(frg->ident);

      //  Read is unanchored in a repeat region, toss it out, but place it with the mate.
      //ejtFrags.insert(frg->ident);
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

  map<FragmentEnd,uint32>                   brkFrags;
  map<FragmentEnd,repeatUniqueBreakPoint>   brkJunct;

  //  Add breakpoints for each of the end points of the repeat region.

  for (uint32 rr=0; rr<regions.size(); rr++) {
    FragmentEnd &rujBgnFrg = regions[rr].rujBgn.breakFrag;
    FragmentEnd &rujEndFrg = regions[rr].rujEnd.breakFrag;

    if (rujBgnFrg.fragId() > 0) {
      brkFrags[rujBgnFrg] = REGION_END_WEIGHT;
      brkJunct[rujBgnFrg] = regions[rr].rujBgn;
    }

    if (rujEndFrg.fragId() > 0) {
      brkFrags[rujEndFrg] = REGION_END_WEIGHT;
      brkJunct[rujEndFrg] = regions[rr].rujEnd;
    }
  }

  //

  sort(evidence.begin(), evidence.end());

  for (uint32 ai=0, bi=0; ai<evidence.size(); ai++) {
    repeatUniqueBreakPoint  ruj;

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
        (ruj.point < regions[bi].bgn))
      continue;

    assert(regions[bi].bgn <= ruj.point);
    assert(ruj.point       <= regions[bi].end);

    //  A new valid break point.

    //  NOTE!  ruj's seem to be different.  We used to save the 5th ruj, and switching to saving the
    //  last showed differences in position (ruj.point) of the break.  The point comes directly from
    //  the evidence[] above, so no surprise.

    brkFrags[evidence[ai].tigFrag]++;
    brkJunct[evidence[ai].tigFrag] = ruj;
  }

  for (map<FragmentEnd,uint32>::iterator it=brkFrags.begin(); it != brkFrags.end(); it++) {
    uint32                  cnt = brkFrags[it->first];
    repeatUniqueBreakPoint  ruj = brkJunct[it->first];

    if (cnt < ISECT_NEEDED_TO_BREAK)
      continue;

    breakpoints.push_back(ruj);
  }

  sort(breakpoints.begin(), breakpoints.end());

  writeLog("markRepeats()--  unitig %d has "F_SIZE_T" interesting junctions:\n",
           target->id(), breakpoints.size());

  for (uint32 ji=0; ji<breakpoints.size(); ji++)
    writeLog("markRepeats()--  junction["F_U32"] at "F_U32"/%c' position "F_U32" repeat %s count "F_U32"\n",
             ji,
             breakpoints[ji].breakFrag.fragId(), breakpoints[ji].breakFrag.frag3p() ? '3' : '5',
             breakpoints[ji].point,
             breakpoints[ji].rptLeft ? "<-" : "->",
             brkFrags[breakpoints[ji].breakFrag]);
}




void
markRepeats_breakUnitigs(UnitigVector                    &unitigs,
                         double                           erateRepeat,
                         Unitig                          *target,
                         vector<overlapPlacement>        &UNUSED(places),
                         vector<repeatUniqueBreakPoint>  &breakpoints,
                         set<uint32>                     &jctFrags,
                         set<uint32>                     &rptFrags,
                         set<uint32>                     &ejtFrags) {

  jctFrags.clear();

  if (breakpoints.size() == 0)
    return;

  uint32   *breakID  = new uint32 [target->ufpath.size()];

  uint32    nextBreakPoint = 0;
  uint32    curr = 1;
  uint32    next = 2;

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode     *frg   = &target->ufpath[fi];
    uint32      bgn   = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
    uint32      end   = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

    //  If out of breakpoints, put all the remaining reads into the current tig.
    if (nextBreakPoint >= breakpoints.size()) {
      breakID[fi]  = curr;

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
          (breakpoints[nextBreakPoint].rptLeft == true)) {
        breakID[fi] = next;
      } else {
        breakID[fi] = curr;
      }

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
    ufNode  &frg = target->ufpath[fi];
    uint32   bid = breakID[fi];

    if (ejtFrags.count(frg.ident) > 0) {
      writeLog("markRepeats()-- EJECT unanchored frag %u from unitig %u\n",
               frg.ident, target->id());
      target->removeFrag(frg.ident);
      continue;
    }

    if (uidToUnitig[bid] == NULL) {
      uidToUnitig[bid] = unitigs.newUnitig(false);  //  Add a new unitig to the unitigs list
      //uidToUnitig[bid]->_isRepeat = true;

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

    for (uint32 ti=0; ti<newTigs.size(); ti++) {
      Unitig *tig   = newTigs[ti];
      uint32  nRept = 0;
      uint32  nUniq = 0;

      for (uint32 fi=0; fi < tig->ufpath.size(); fi++) {
        ufNode  &frg = tig->ufpath[fi];

        if (rptFrags.count(frg.ident) == 1)
          nRept++;
        else
          nUniq++;
      }

      if (nRept > nUniq)
        tig->_isRepeat = true;

      writeLog("markRepeats()--   unitig %u of length %u with %ld fragments (%u %.4f repeat and %u %.4f non-repeat).\n",
               tig->id(),
               tig->getLength(),
               tig->ufpath.size(),
               nRept, (double)nRept / (nRept + nUniq),
               nUniq, (double)nUniq / (nRept + nUniq));
    }

    writeLog("markRepeats()-- DELETE unitig %d\n", target->id());
    unitigs[target->id()] = NULL;
    delete target;
  }

  //  Run back over the ejected frags, and place them at their best location.

  for (set<uint32>::iterator it=ejtFrags.begin(); it!=ejtFrags.end(); it++) {
    writeLog("markRepeats()-- EJECT frag "F_U32"\n", *it);
    placeFragInBestLocation(unitigs, erateRepeat, *it);
  }

  writeLog("markRepeats()-- FINISHED.\n");
}



void
markRepeats_shatterRepeats(UnitigVector   &unitigs,
                           set<uint32>    &jctFrags,
                           set<uint32>    &rptFrags) {

  //  For each junction read (defined to be the first/last read in a repeat unitig), shatter the unitig.

  for (set<uint32>::iterator it=jctFrags.begin(); it!=jctFrags.end(); it++) {
    uint32   iid = *it;
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

  //  For each repeat read (defined to be a read contained nearly entirely in a repeat region), count
  //  the number that are still in a unitig.

  for (set<uint32>::iterator it=rptFrags.begin(); it!=rptFrags.end(); it++) {
    uint32   iid = *it;
    uint32   ti  = Unitig::fragIn(iid);
    Unitig  *rpt = unitigs[ti];

    if ((ti == 0) || (rpt == NULL))
      //  Already shattered?
      continue;

    writeLog("markRepeats()--  frag "F_U32" covered by repeats, but still in unitig %d\n",
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
            double erateRepeat,
            Unitig *target,
            uint32 minOverlap,
            bool shatterRepeats) {

  set<uint32>                     ovlFrags;
  set<uint32>                     rptFrags;  //  Frag IIDs of fragments covered by repeat alignments
  set<uint32>                     jctFrags;  //  Frag IIDs of the first/last fragment in a repeat unitig
  set<uint32>                     ejtFrags;  //  Frag IIDs of frags we should eject instead of split

  double                          meanError = 0;
  double                          stddevError = 0;

  intervalList<int32>             aligned;

  vector<overlapPlacement>        places;

  vector<repeatRegion>            regions;
  vector<repeatJunctionEvidence>  evidence;

  vector<repeatUniqueBreakPoint>  breakpoints;

  //  Build a list of all the fragments that have overlaps to this unitig.
  markRepeats_buildOverlapList(target, erateRepeat, ovlFrags);

  //  Decide what a decent alignment should look like.
  markRepeats_computeUnitigErrorRate(unitigs, erateRepeat, target, meanError, stddevError);

  //  For each overlapping fragment, place it and process.
  markRepeats_placeAndProcessOverlaps(unitigs, erateRepeat, target, meanError, stddevError, ovlFrags, aligned, evidence);

  //  Convert 'aligned' into regions, throwing out weak ones and those contained in a fragment.
  markRepeats_filterIntervalsSpannedByFragment(target, aligned, regions, minOverlap);

  markRepeats_findFragsInRegions(target, regions, rptFrags, ejtFrags, minOverlap);

  //  Discard junctions that are not in a remaining region.
  markRepeats_filterJunctions(target, regions, evidence, breakpoints);

  //  Split at whatever junctions remain.

  //  You'd think declaring this a critical region would work, but it resulted in deadlock on
  //    Linux 2.6.32-279.22.1.el6.x86_64
  //    g++ (GCC) 4.7.1
  //#pragma omp critical

  omp_set_lock(&markRepeat_breakUnitigs_Lock);
  markRepeats_breakUnitigs(unitigs, erateRepeat, target, places, breakpoints, jctFrags, rptFrags, ejtFrags);
  omp_unset_lock(&markRepeat_breakUnitigs_Lock);

  //  For each repeat unitig, shatter into fragments (not singleton unitigs) so we can later re-BOG.
  if (shatterRepeats)
    markRepeats_shatterRepeats(unitigs, jctFrags, rptFrags);
}




//void
//markChimera(UnitigVector &unitigs,
//            double erateRepeat,
//            Unitig *target) {
//}



void
mergeSplitJoin(UnitigVector &unitigs,
               double UNUSED(erateGraph), double erateBubble, double UNUSED(erateMerge), double erateRepeat,
               const char *prefix,
               uint32 minOverlap,
               bool shatterRepeats,
               uint64 genomeSize) {

  //logFileFlags |= LOG_PLACE_FRAG;
  //logFileFlags &= ~LOG_PLACE_FRAG;

  //  BUILD A LIST OF ALL INTERSECTIONS - build a reverse mapping of all BestEdges that are between
  //  unitigs.  For each fragment, we want to have a list of the incoming edges from other unitigs.

  intersectionList  *ilist = new intersectionList(unitigs);

  //ilist->logIntersections();

#if 0
  {
    Unitig        *target = unitigs[5];

    writeLog("popBubbles()-- WORKING on unitig %d/"F_SIZE_T" with %ld fragments.\n",
             target->id(), unitigs.size(), target->ufpath.size());

    mergeBubbles(unitigs, erateBubble, target, ilist);
    //stealBubbles(unitigs, erateBubble, target, ilist);
    markRepeats(unitigs, erateRepeat, target, minOverlap, shatterRepeats);
    //markChimera(unitigs, erateRepeat, target);
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
        (target->ufpath.size() < 1) ||  //  was 15
        (target->getLength() < 300))
      continue;

    //writeLog("popBubbles()-- WORKING on unitig %d/"F_SIZE_T" of length %u with %ld fragments.\n",
    //         target->id(), unitigs.size(), target->getLength(), target->ufpath.size());

    mergeBubbles(unitigs, erateBubble, target, ilist);
    //stealBubbles(unitigs, erateBubble, target, ilist);
  }

  reportOverlapsUsed(unitigs, prefix, "popBubbles");
  reportUnitigs(unitigs, prefix, "popBubbles", genomeSize);

  //  Since we create new unitigs for any of the splits, we need to remember
  //  where to stop.  We don't want to re-examine any of the split unitigs.
  //  Reevaluating seems to just trim off a few fragments at the end of the unitig.

  uint32  tiLimit = unitigs.size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  setLogFile(prefix, "mergeSplitJoin");
  writeLog("repeatDetect()-- working on "F_U32" unitigs, with "F_U32" threads.\n", tiLimit, numThreads);

  omp_init_lock(&markRepeat_breakUnitigs_Lock);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig        *target = unitigs[ti];

    if ((target == NULL) ||
        (target->ufpath.size() < 15) ||
        (target->getLength() < 300))
      continue;

    //writeLog("repeatDetect()-- WORKING on unitig %d/"F_SIZE_T" of length %u with %ld fragments.\n",
    //         target->id(), unitigs.size(), target->getLength(), target->ufpath.size());

    markRepeats(unitigs, erateRepeat, target, minOverlap, shatterRepeats);
    //markChimera(unitigs, erateRepeat, target);
  }

  omp_destroy_lock(&markRepeat_breakUnitigs_Lock);

  reportOverlapsUsed(unitigs, prefix, "mergeSplitJoin");
  reportUnitigs(unitigs, prefix, "mergeSplitJoin", genomeSize);

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

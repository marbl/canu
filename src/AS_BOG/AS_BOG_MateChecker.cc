
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

static const char *rcsid = "$Id: AS_BOG_MateChecker.cc,v 1.89 2010-06-18 16:21:15 skoren Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_MateChecker.hh"

#include "AS_OVL_overlap.h"  //  For DEFAULT_MIN_OLAP_LEN




void MateChecker::checkUnitigGraph(UnitigGraph& tigGraph, int badMateBreakThreshold) {

  computeGlobalLibStats(tigGraph);

  fprintf(stderr, "==> STARTING MATE BASED SPLITTING.\n");

  tigGraph.checkUnitigMembership();
  evaluateMates(tigGraph);

  //  Move unhappy contained fragments before attempting mate based
  //  splitting.  We know they're broken already, and if scaffolder
  //  can put them back, we'll let it.
  //
  fprintf(stderr, "==> MOVE CONTAINS #1\n");
  moveContains(tigGraph);

  tigGraph.reportOverlapsUsed("overlaps.aftermatecheck1");
  tigGraph.checkUnitigMembership();
  evaluateMates(tigGraph);

  //  This should do absolutely nothing.  If it does, something is
  //  broken.  By just ejecting unhappy contains, nothing should be
  //  disconnected.
  //
  fprintf(stderr, "==> SPLIT DISCONTINUOUS #1\n");
  splitDiscontinuousUnitigs(tigGraph);

  tigGraph.reportOverlapsUsed("overlaps.aftermatecheck2");
  tigGraph.checkUnitigMembership();
  evaluateMates(tigGraph);

  //  DON'T DO MATE BASED SPLITTING
  //return;

  fprintf(stderr, "==> SPLIT BAD MATES\n");
  {
    //  Need to get rid of this cMap guy
    ContainerMap       cMap;

    for (int  ti=0; ti<tigGraph.unitigs->size(); ti++) {
      Unitig  *tig = (*tigGraph.unitigs)[ti];

      if ((tig == NULL) || (tig->getNumFrags() < 2))
        continue;

      UnitigBreakPoints* breaks = computeMateCoverage(tig, tigGraph.bog_ptr, badMateBreakThreshold);
      UnitigVector*      newUs  = tigGraph.breakUnitigAt(cMap, tig, *breaks);

      if (newUs != NULL) {
        delete tig;
        (*tigGraph.unitigs)[ti] = NULL;
        tigGraph.unitigs->insert(tigGraph.unitigs->end(), newUs->begin(), newUs->end());
      }

      delete newUs;
      delete breaks;
    }
  }

  tigGraph.reportOverlapsUsed("overlaps.aftermatecheck3");
  tigGraph.checkUnitigMembership();
  evaluateMates(tigGraph);

  //  DO MATE BASED SPLITTING BUT NOTHING ELSE
  //return;

  //  The splitting code above is not smart enough to move contained
  //  fragments with containers.  This leaves unitigs disconnected.
  //  We break those unitigs here.
  //
  fprintf(stderr, "==> SPLIT DISCONTINUOUS #2\n");
  splitDiscontinuousUnitigs(tigGraph);

  tigGraph.reportOverlapsUsed("overlaps.aftermatecheck4");
  tigGraph.checkUnitigMembership();
  evaluateMates(tigGraph);

  //  But now, all the splitting probably screwed up happiness of
  //  contained fragments, of left some unhappy fragments in a unitig
  //  that just lost the container.
  //
  fprintf(stderr, "==> MOVE CONTAINS #2\n");
  moveContains(tigGraph);

  tigGraph.reportOverlapsUsed("overlaps.aftermatecheck5");
  tigGraph.checkUnitigMembership();
  evaluateMates(tigGraph);

  //  Do one last check for disconnected unitigs.
  //
  fprintf(stderr, "==> SPLIT DISCONTINUOUS #3\n");
  splitDiscontinuousUnitigs(tigGraph);

  tigGraph.reportOverlapsUsed("overlaps.aftermatecheck6");
  tigGraph.checkUnitigMembership();
  evaluateMates(tigGraph);
}



void
MateChecker::accumulateLibraryStats(Unitig *utg) {

  //  Check each frag in unitig for it's mate
  for (uint32 fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
    DoveTailNode  *frag = &(*utg->dovetail_path_ptr)[fi];

    if (_fi->mateIID(frag->ident) == 0)
      //  Unmated fragment.
      continue;

    if (utg->id() != utg->fragIn(_fi->mateIID(frag->ident)))
      //  Mate in a different unitig.
      continue;

    uint32        mi   = utg->pathPosition(_fi->mateIID(frag->ident));
    DoveTailNode *mate = &(*utg->dovetail_path_ptr)[mi];

    if (frag->ident < mate->ident)
      //  Only do this once for each mate pair.
      continue;

#warning assumes innie mate pairs
    if (isReverse(frag->position) == isReverse(mate->position))
      //  Same orient mates, not a good mate pair.
      continue;

    //  Compute the insert size.

    uint32  insertSize = 0;

    if (isReverse(frag->position)) {
      if (frag->position.end <= mate->position.bgn)
        //  Misordered, not a good mate pair.  NOTE: this is specifically allowing mates that
        //  overlap, i.e., insert size is less than the sum of fragment lengths.
        continue;

      //  We must have, at least, the following picture, where the relationship on the left is
      //  strict.  'frag' is allowed to move to the right, and 'mate' is allowed to move to the
      //  left, but moving either in the other direction turns this into an outtie relationship.
      //
      //  (end)  <----- (bgn)  frag
      //  (bgn) ----->  (end) mate
      //
      assert(mate->position.bgn < frag->position.end);

      //  However, we aren't guaranteed that the right side is ordered properly (i.e., misordered by
      //  one fragment contained in the other fragment).  We hope this doesn't happen, but if it
      //  does, we'll catch it.  The insertSize defaults to zero, which is then ignored below.
      //
      if (mate->position.bgn < frag->position.bgn)
        insertSize = frag->position.bgn - mate->position.bgn;

    } else {
      //  (See comments above)
      if (mate->position.end <= frag->position.bgn)
        continue;

      assert(frag->position.bgn < mate->position.end);

      if (frag->position.bgn < mate->position.bgn)
        insertSize = mate->position.bgn - frag->position.bgn;
    }

    if (insertSize > 0) {
      uint32 di = _fi->libraryIID(frag->ident);

      if (_globalStats[di].distancesMax <= _globalStats[di].distancesLen) {
        _globalStats[di].distancesMax *= 2;
        uint32 *d = new uint32 [_globalStats[di].distancesMax];
        memcpy(d, _globalStats[di].distances, sizeof(uint32) * _globalStats[di].distancesLen);
        delete [] _globalStats[di].distances;
        _globalStats[di].distances = d;
      }

      _globalStats[di].distances[_globalStats[di].distancesLen++] = insertSize;
    }
  }
}



void MateChecker::computeGlobalLibStats(UnitigGraph& tigGraph) {

  assert(_globalStats == NULL);

  _globalStats = new DistanceCompute [_fi->numLibraries() + 1];

  for (uint32 i=1; i < _fi->numLibraries()+1; i++) {
    _globalStats[i].mean     = _fi->mean(i);
    _globalStats[i].stddev   = _fi->stddev(i);
    _globalStats[i].samples  = _fi->numMatesInLib(i);
  }

  for (int  ti=0; ti<tigGraph.unitigs->size(); ti++) {
    Unitig        *utg = (*tigGraph.unitigs)[ti];

    if ((utg == NULL) ||
        (utg->dovetail_path_ptr->size() < 2))
      continue;

    accumulateLibraryStats(utg);
  }

  fprintf(stderr, "LIB\tnDist\tmedian\t1/3rd\t2/3rd\tmaxDiff\tmin\tmax\tnumGood\tmean\tstddev\n");

  //  Disregard outliers (those outside 5 (estimated) stddevs) and recompute global stddev

  for (uint32 i=1; i<_fi->numLibraries()+1; i++) {
    sort(_globalStats[i].distances, _globalStats[i].distances + _globalStats[i].distancesLen + 1);

    int median   = _globalStats[i].distances[_globalStats[i].distancesLen * 1 / 2];
    int oneThird = _globalStats[i].distances[_globalStats[i].distancesLen * 1 / 3];
    int twoThird = _globalStats[i].distances[_globalStats[i].distancesLen * 2 / 3];

    int aproxStd = MAX(median - oneThird, twoThird - median);

    int biggest  = median + aproxStd * 5;
    int smallest = median - aproxStd * 5;

    uint32    numPairs   = 0;
    double    sumDists   = 0.0;
    double    sumSquares = 0.0;

    for (uint32 d=0; d<_globalStats[i].distancesLen; d++)
      if ((smallest <= _globalStats[i].distances[d]) &&
          (_globalStats[i].distances[d] <= biggest)) {
        numPairs++;
        sumDists += _globalStats[i].distances[d];
      }

    if (numPairs > 0)
      _globalStats[i].mean = sumDists / numPairs;
    else
      _globalStats[i].mean = 0;

    for (uint32 d=0; d<_globalStats[i].distancesLen; d++)
      if ((smallest <= _globalStats[i].distances[d]) &&
          (_globalStats[i].distances[d] <= biggest))
        sumSquares += ((_globalStats[i].distances[d] - _globalStats[i].mean) *
                       (_globalStats[i].distances[d] - _globalStats[i].mean));

    if (numPairs > 1)
      _globalStats[i].stddev = sqrt(sumSquares / (numPairs - 1));
    else
      _globalStats[i].stddev = 0.0;

    _globalStats[i].samples = numPairs;

    fprintf(stderr, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.1f\n",
            i, _globalStats[i].distancesLen, median, oneThird, twoThird, aproxStd, smallest, biggest,
            numPairs, _globalStats[i].mean, _globalStats[i].stddev);
  }
}



void
MateChecker::evaluateMates(UnitigGraph &tigGraph) {

  //  [0] -- BOTH frag and mate are dovetail
  //  [1] -- ONE frag dovetail, ONE frag contained
  //  [2] -- BOTH frag and mate are contained

  uint64   unmated[3]         = { 0, 0, 0 };
  uint64   mated[3]           = { 0, 0, 0 };
  uint64   different[3][3]    = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
  uint64   happy[3]           = { 0, 0, 0 };
  uint64   grumpy[3]          = { 0, 0, 0 };

  for (int  ti=0; ti<tigGraph.unitigs->size(); ti++) {
    Unitig  *thisUtg = (*tigGraph.unitigs)[ti];

    if ((thisUtg == NULL) ||
        (thisUtg->dovetail_path_ptr->empty()) ||
        (thisUtg->dovetail_path_ptr->size() == 1))
      continue;

    MateLocation          positions(_fi, thisUtg, _globalStats);

    for (uint32 fi=0; fi<thisUtg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode  *thisFrg = &(*thisUtg->dovetail_path_ptr)[fi];

      uint32  thisFrgID = thisFrg->ident;
      uint32  mateFrgID = _fi->mateIID(thisFrg->ident);

      BestContainment *thiscont = tigGraph.bog_ptr->getBestContainer(thisFrgID);
      BestContainment *matecont = tigGraph.bog_ptr->getBestContainer(mateFrgID);

      uint32  type = (thiscont != NULL) + (matecont != NULL);

      //  Trivial case, not a mated fragment.

      if (mateFrgID == 0) {
        unmated[type]++;
        continue;
      }

      uint32  thisUtgID = thisUtg->fragIn(thisFrgID);
      uint32  mateUtgID = thisUtg->fragIn(mateFrgID);

      MateLocationEntry  mloc     = positions.getById(thisFrg->ident);

      //  Skip this fragment, unless it is mleFrgID1.  Fragments with mates in other unitigs
      //  are always listed in ID1, and mates completely in this unitig should be counted once.
      if (mloc.mleFrgID1 != thisFrgID)
        continue;

      mated[type]++;

      //  Easy case, both fragments in the same unitig.

      if (thisUtgID == mateUtgID) {
        assert(mloc.mleUtgID1 == thisUtg->id());
        assert(mloc.mleUtgID2 == thisUtg->id());
        assert(mloc.mleFrgID1 == thisFrgID);
        assert(mloc.mleFrgID2 == mateFrgID);

        if (mloc.isGrumpy == false)
          happy[type]++;
        else
          grumpy[type]++;
        continue;
      }

      //  Hard case, fragments in different unitigs.  We want to distinguish between
      //  three cases:
      //    1) mates at the end that could potentially join unitigs across a gap
      //    2) mate at the end to an interior mate -- possibly a repeat
      //    3) both interior mates

      assert(mloc.mleUtgID1 == thisUtg->id());
      assert(mloc.mleUtgID2 == 0);
      assert(mloc.mleFrgID1 == thisFrgID);
      assert(mloc.mleFrgID2 == 0);

      //  Get the mate frag.

      Unitig        *mateUtg = (*tigGraph.unitigs)[mateUtgID];
      DoveTailNode  *mateFrg = &(*mateUtg->dovetail_path_ptr)[mateUtg->pathPosition(mateFrgID)];

      //differentSum[type]++;

      bool  fragIsInterior = false;
      bool  mateIsInterior = false;

      uint32  lib           = _fi->libraryIID(thisFrg->ident);
      uint32  minInsertSize = _globalStats[lib].mean - 3 * _globalStats[lib].stddev;
      uint32  maxInsertSize = _globalStats[lib].mean + 3 * _globalStats[lib].stddev;

      if (thisFrg->position.bgn < thisFrg->position.end) {
        //  Fragment is forward, so mate should be after it.
        if (thisUtg->getLength() - thisFrg->position.bgn > maxInsertSize)
          fragIsInterior = true;
      } else {
        //  Fragment is reverse, so mate should be before it.
        if (thisFrg->position.bgn > maxInsertSize)
          fragIsInterior = true;
      }

      if (mateFrg->position.bgn < mateFrg->position.end) {
        //  Fragment is forward, so mate should be after it.
        if (mateUtg->getLength() - mateFrg->position.bgn > maxInsertSize)
          mateIsInterior = true;
      } else {
        //  Fragment is reverse, so mate should be before it.
        if (mateFrg->position.bgn > maxInsertSize)
          mateIsInterior = true;
      }

      uint32  dtyp = (fragIsInterior == true) + (mateIsInterior == true);

      different[type][dtyp]++;
    }
  }

  fprintf(stderr, "MATE HAPPINESS (dove/dove):  unmated %11"F_U64P"  mated %11"F_U64P"  sameTig: happy %11"F_U64P" grumpy %11"F_U64P"  diffTig: end-end %11"F_U64P" end-int %11"F_U64P" int-int %11"F_U64P"\n",
          unmated[0], mated[0], happy[0], grumpy[0], different[0][0], different[0][1], different[0][2]);
  fprintf(stderr, "MATE HAPPINESS (dove/cont):  unmated %11"F_U64P"  mated %11"F_U64P"  sameTig: happy %11"F_U64P" grumpy %11"F_U64P"  diffTig: end-end %11"F_U64P" end-int %11"F_U64P" int-int %11"F_U64P"\n",
          unmated[1], mated[1], happy[1], grumpy[1], different[1][0], different[1][1], different[1][2]);
  fprintf(stderr, "MATE HAPPINESS (cont/cont):  unmated %11"F_U64P"  mated %11"F_U64P"  sameTig: happy %11"F_U64P" grumpy %11"F_U64P"  diffTig: end-end %11"F_U64P" end-int %11"F_U64P" int-int %11"F_U64P"\n",
          unmated[2], mated[2], happy[2], grumpy[2], different[2][0], different[2][1], different[2][2]);
}



//  Make sure that contained fragments are in the same unitig
//  as their container.  Due to sorting, contained fragments
//  can come much later in the unitig:
//
//  ------------1
//    -------------2
//       --------------3
//         ----4 (contained in 1, too much error keeps it out of 2 and 3)
//
//  So, our first pass is to move contained fragments around.
//
void MateChecker::moveContains(UnitigGraph& tigGraph) {

  for (int  ti=0; ti<tigGraph.unitigs->size(); ti++) {
    Unitig  *thisUnitig = (*tigGraph.unitigs)[ti];

    if ((thisUnitig == NULL) ||
        (thisUnitig->dovetail_path_ptr->empty()) ||
        (thisUnitig->dovetail_path_ptr->size() == 1))
      continue;

    MateLocation positions(_fi, thisUnitig, _globalStats);

    DoveTailNode         *frags         = new DoveTailNode [thisUnitig->dovetail_path_ptr->size()];
    int                   fragsLen      = 0;

    bool                  verbose       = false;

    if (verbose)
      fprintf(stderr, "moveContain unitig %d\n", thisUnitig->id());

    for (DoveTailIter fragIter = thisUnitig->dovetail_path_ptr->begin();
         fragIter != thisUnitig->dovetail_path_ptr->end();
         fragIter++) {

      BestContainment   *bestcont   = tigGraph.bog_ptr->getBestContainer(fragIter->ident);
      MateLocationEntry  mloc       = positions.getById(fragIter->ident);

      uint32  thisFrgID = fragIter->ident;
      uint32  contFrgID = (bestcont) ? bestcont->container : 0;
      uint32  mateFrgID = _fi->mateIID(fragIter->ident);

      uint32  thisUtgID = thisUnitig->fragIn(thisFrgID);
      uint32  contUtgID = thisUnitig->fragIn(contFrgID);
      uint32  mateUtgID = thisUnitig->fragIn(mateFrgID);

      //  id1 != 0 -> we found the fragment in the mate happiness table
      //  isBad -> and the mate is unhappy.
      //
      //  What's id1 vs id2 in MateLocationEntry?  Dunno.  All I
      //  know is that if there is no mate present, one of those
      //  will be 0.  (Similar test used above too.)
      //
      bool    isMated    = (mateFrgID > 0);
      bool    isGrumpy   = ((isMated) && (mloc.mleFrgID1 != 0) && (mloc.mleFrgID2 != 0) && (mloc.isGrumpy == true));

      //
      //  Figure out what to do.
      //

      bool    moveToContainer = false;
      bool    moveToSingleton = false;

      if        ((fragIter->contained == 0) && (bestcont == NULL)) {
        //  CASE 1:  Not contained.  Leave the fragment here.
        //fprintf(stderr, "case1 frag %d fragsLen %d\n", thisFrgID, fragsLen);

      } else if (isMated == false) {
        //  CASE 2: Contained but not mated.  Move to be with the
        //  container (if the container isn't here).
        //fprintf(stderr, "case2 frag %d contID %d fragsLen %d\n", thisFrgID, contUtgID, fragsLen);

        if (thisUtgID != contUtgID)
          moveToContainer = true;

      } else if ((isGrumpy == true) && (thisUtgID == mateUtgID)) {
        //  CASE 3: Not happy, and the frag and mate are together.
        //  Kick out to a singleton.

        //fprintf(stderr, "case3 frag %d utg %d mate %d utg %d cont %d utg %d fragsLen %d\n",
        //        thisFrgID, thisUtgID, mateFrgID, mateUtgID, contFrgID, contUtgID, fragsLen);

        if (thisUtgID == mateUtgID)
          moveToSingleton = true;

      } else {

        //  This makes for some ugly code (we break the nice if else
        //  if else structure we had going on) but the next two cases
        //  need to know if there is an overlap to the rest of the
        //  unitig.

        bool  hasOverlap   = (thisUtgID == contUtgID);
        bool  allContained = false;


        if (hasOverlap == false) {
          if (fragsLen == 0) {
            //  The first fragment.  Check fragments after to see if
            //  there is an overlap (note only frags with an overlap
            //  in the layout are tested).  In rare cases, we ejected
            //  the container, and left a containee with no overlap to
            //  fragments remaining.
            //
            //  Note that this checks if there is an overlap to the
            //  very first non-contained (aka dovetail) fragment ONLY.
            //  If there isn't an overlap to the first non-contained
            //  fragment, then that fragment will likely NOT align
            //  correctly.

            DoveTailIter  ft = fragIter + 1;

            //  Skip all the contains.
            while ((ft != thisUnitig->dovetail_path_ptr->end()) &&
                   (tigGraph.bog_ptr->isContained(ft->ident) == true) &&
                   (MAX(fragIter->position.bgn, fragIter->position.end) < MIN(ft->position.bgn, ft->position.end)))
              ft++;

            //  If the frag is not contained (we could be the
            //  container), and overlaps in the layout, see if there
            //  is a real overlap.
            if ((ft != thisUnitig->dovetail_path_ptr->end()) &&
                (tigGraph.bog_ptr->isContained(ft->ident) == false) &&
                (MAX(fragIter->position.bgn, fragIter->position.end) < MIN(ft->position.bgn, ft->position.end)))
              hasOverlap = tigGraph.bog_ptr->containHaveEdgeTo(thisFrgID, ft->ident);
          } else {
            //  Not the first fragment, search for an overlap to an
            //  already placed frag.

            DoveTailIter  ft = fragIter;

            do {
              ft--;

              //  OK to overlap to a contained frag; he could be our
              //  container.

              hasOverlap = tigGraph.bog_ptr->containHaveEdgeTo(thisFrgID, ft->ident);

              //  Stop if we found an overlap, or we just checked the
              //  first frag in the unitig, or we no longer overlap in
              //  the layout.
            } while ((hasOverlap == false) &&
                     (ft != thisUnitig->dovetail_path_ptr->begin()) &&
                     (MIN(fragIter->position.bgn, fragIter->position.end) < MAX(ft->position.bgn, ft->position.end)));
          }
        }  //  end of hasOverlap


        //  An unbelievabe special case.  When the unitig is just a
        //  single container fragment (and any contained frags under
        //  it) rule 4 breaks.  The first fragment has no overlap (all
        //  later reads are contained) and so we want to eject it to a
        //  new unitig.  Since there are multiple fragments in this
        //  unitig, the ejection occurs.  Later, all the contains get
        //  moved to the new unitig.  And we repeat.  To prevent, we
        //  abort the ejection if the unitig is all contained in one
        //  fragment.
        //
        if (fragsLen == 0) {
          allContained = true;

          for (DoveTailIter  ft = fragIter + 1;
               ((allContained == true) &&
                (ft != thisUnitig->dovetail_path_ptr->end()));
               ft++)
            allContained = tigGraph.bog_ptr->isContained(ft->ident);
        }



        if (isGrumpy == true) {
          //  CASE 4: Not happy and not with the mate.  This one is a
          //  bit of a decision.
          //
          //  If an overlap exists to the rest of the unitig, we'll
          //  leave it here.  We'll also leave it here if it is the
          //  rest of the unitig is all contained in this fragment.
          //
          //  If no overlap, and the mate and container are in the
          //  same unitig, we'll just eject.  That also implies the
          //  other unitig is somewhat large, at least as big as the
          //  insert size.
          //
          //  Otherwise, we'll move to the container and cross our
          //  fingers we place it correctly.  The alternative is to
          //  eject, and hope that we didn't also eject the mate to a
          //  singleton.

          //fprintf(stderr, "case4 frag %d utg %d mate %d utg %d cont %d utg %d fragsLen %d\n",
          //        thisFrgID, thisUtgID, mateFrgID, mateUtgID, contFrgID, contUtgID, fragsLen);

          if ((hasOverlap == false) && (allContained == false))
            if (mateUtgID == contUtgID)
              moveToSingleton = true;
            else
              moveToContainer = true;

        } else {
          //  CASE 5: Happy!  If with container, or an overlap exists to
          //  some earlier fragment, leave it here.  Otherwise, eject it
          //  to a singleton.  The fragment is ejected instead of moved
          //  to be with its container since we don't know which is
          //  correct - the mate or the overlap.
          //
          //  If not happy, we've already made sure that the mate is not
          //  here (that was case 3).

          //fprintf(stderr, "case5 frag %d utg %d mate %d utg %d cont %d utg %d fragsLen %d\n",
          //        thisFrgID, thisUtgID, mateFrgID, mateUtgID, contFrgID, contUtgID, fragsLen);

          //  If no overlap (so not with container or no overlap to
          //  other frags) eject.
          if ((hasOverlap == false) && (allContained == false))
            moveToSingleton = true;
        }
      }  //  End of cases

      //
      //  Do it.
      //

      if (moveToContainer == true) {
        //  Move the fragment to be with its container.

        Unitig         *thatUnitig = (*tigGraph.unitigs)[contUtgID];
        DoveTailNode    containee  = *fragIter;

        assert(thatUnitig->id() == contUtgID);

        //  Nuke the fragment in the current list
        fragIter->ident        = 999999999;
        fragIter->contained    = 999999999;
        fragIter->position.bgn = 0;
        fragIter->position.end = 0;

        assert(thatUnitig->id() == contUtgID);

        if (verbose)
          fprintf(stderr, "Moving contained fragment %d from unitig %d to be with its container %d in unitig %d\n",
                  thisFrgID, thisUtgID, contFrgID, contUtgID);

        assert(bestcont->container == contFrgID);

        thatUnitig->addContainedFrag(thisFrgID, bestcont, verbose);
        assert(thatUnitig->id() == Unitig::fragIn(thisFrgID));

      } else if ((moveToSingleton == true) && (thisUnitig->getNumFrags() != 1)) {
        //  Eject the fragment to a singleton (unless we ARE the singleton)
        Unitig        *singUnitig  = new Unitig(verbose);
        DoveTailNode    containee  = *fragIter;

        //  Nuke the fragment in the current list
        fragIter->ident        = 999999999;
        fragIter->contained    = 999999999;
        fragIter->position.bgn = 0;
        fragIter->position.end = 0;

        if (verbose)
          fprintf(stderr, "Ejecting unhappy contained fragment %d from unitig %d into new unitig %d\n",
                  thisFrgID, thisUtgID, singUnitig->id());

        containee.contained = 0;

        singUnitig->addFrag(containee, -MIN(containee.position.bgn, containee.position.end), verbose);

        tigGraph.unitigs->push_back(singUnitig);
        thisUnitig = (*tigGraph.unitigs)[ti];  //  Reset the pointer; unitigs might be reallocated

      } else {
        //  Leave fragment here.  Copy the fragment to the list -- if
        //  we need to rebuild the unitig (because fragments were
        //  removed), the list is used, otherwise, we have already
        //  made the changes needed.
        //
        //  Also, very important, update our containment mark.  If our
        //  container was moved, but we stayed put because of a happy
        //  mate, we're still marked as being contained.  Rather than
        //  put this check in all the places where we stay put in the
        //  above if-else-else-else, it's here.

        if ((fragIter->contained) && (thisUtgID != contUtgID))
          fragIter->contained = 0;

        frags[fragsLen] = *fragIter;
        fragsLen++;
      }

    }  //  over all frags

    //  Now, rebuild this unitig if we made changes.

    if (fragsLen != thisUnitig->dovetail_path_ptr->size()) {
      if (verbose)
        fprintf(stderr, "Rebuild unitig %d after removing contained fragments.\n", thisUnitig->id());

      delete thisUnitig->dovetail_path_ptr;

      thisUnitig->dovetail_path_ptr = new DoveTailPath;

      //  Occasionally, we move all fragments out of the original unitig.  Might be worth checking
      //  if that makes sense!!
      //
#warning EMPTIED OUT A UNITIG
      if (fragsLen > 0) {
        //  No need to resort.  Offsets only need adjustment if the first fragment is thrown out.
        //  If not, splitOffset will be zero.
        //
        int splitOffset = -MIN(frags[0].position.bgn, frags[0].position.end);

        //  This is where we clean up from the splitting not dealing with contained fragments -- we
        //  force the first frag to be uncontained.
        //
        frags[0].contained = 0;

        for (int i=0; i<fragsLen; i++)
          thisUnitig->addFrag(frags[i], splitOffset, verbose);
      }
    }

    delete [] frags;
    frags = NULL;

  }  //  Over all unitigs
}




//  After splitting and ejecting some contains, check for discontinuous unitigs.
//
void MateChecker::splitDiscontinuousUnitigs(UnitigGraph& tigGraph) {

  bool  verbose = false;

  for (int  ti=0; ti<tigGraph.unitigs->size(); ti++) {
    Unitig  *unitig = (*tigGraph.unitigs)[ti];

    if ((unitig == NULL) ||
        (unitig->dovetail_path_ptr->empty()) ||
        (unitig->dovetail_path_ptr->size() == 1))
      continue;

    //  Check for discontinuities

    DoveTailConstIter     fragIter = unitig->dovetail_path_ptr->begin();
    int                   maxEnd   = 0;

    DoveTailNode         *splitFrags    = new DoveTailNode [unitig->dovetail_path_ptr->size()];
    int                   splitFragsLen = 0;

    while (fragIter != unitig->dovetail_path_ptr->end()) {

      //  If this is the first frag in this block (we are at
      //  the start of a unitig, or just split off a new
      //  unitig), remember the end location.
      //
      if (splitFragsLen == 0) {
        maxEnd =  MAX(fragIter->position.bgn, fragIter->position.end);
      }

      //  We require at least (currently 40bp, was 10bp hardcoded
      //  here) of overlap between fragments.  If we don't have that,
      //  split off the fragments we've seen.
      //
      //  10bp was a bad choice.  It caught most of the breaks, but
      //  missed one class; when a container fragment is moved out of
      //  the unitig, fragments contained in there are marked as
      //  uncontained.  That container fragment could have been the
      //  one holding the unitig together:
      //
      //  -----------------   <- container (removed)
      //    --------
      //      ---------
      //              -----------------
      //
      //  Because the two small guys are marked as uncontained, they
      //  are assumed to have a good dovetail overlap.
      //
      if (maxEnd - AS_OVERLAP_MIN_LEN < MIN(fragIter->position.bgn, fragIter->position.end)) {

        //  If there is exactly one fragment, and it's contained, and
        //  it's not mated, move it to the container.  (This has a
        //  small positive benefit over just making every read a
        //  singleton).
        //
        if ((splitFragsLen == 1) &&
            (_fi->mateIID(splitFrags[0].ident) == 0) &&
            (splitFrags[0].contained != 0)) {

          Unitig           *dangler  = (*tigGraph.unitigs)[unitig->fragIn(splitFrags[0].contained)];
          BestContainment  *bestcont = tigGraph.bog_ptr->getBestContainer(splitFrags[0].ident);

          assert(dangler->id() == unitig->fragIn(splitFrags[0].contained));

          if (verbose)
            fprintf(stderr, "Dangling contained fragment %d in unitig %d -> move them to container unitig %d\n",
                    splitFrags[0].ident, unitig->id(), dangler->id());

          dangler->addContainedFrag(splitFrags[0].ident, bestcont, verbose);
          assert(dangler->id() == Unitig::fragIn(splitFrags[0].ident));

        } else {
          Unitig *dangler = new Unitig(verbose);

          if (verbose)
            fprintf(stderr, "Dangling fragments in unitig %d -> move them to unitig %d\n", unitig->id(), dangler->id());

          int splitOffset = -MIN(splitFrags[0].position.bgn, splitFrags[0].position.end);

          //  This should already be true, but we force it still
          splitFrags[0].contained = 0;

          for (int i=0; i<splitFragsLen; i++)
            dangler->addFrag(splitFrags[i], splitOffset, verbose);

          tigGraph.unitigs->push_back(dangler);
          unitig = (*tigGraph.unitigs)[ti];
        }

        //  We just split out these fragments.  Reset the list.
        splitFragsLen = 0;
      }  //  End break

      splitFrags[splitFragsLen++] = *fragIter;

      maxEnd = MAX(maxEnd, MAX(fragIter->position.bgn, fragIter->position.end));

      fragIter++;
    }  //  End of unitig fragment iteration

    //  If we split this unitig, the length of the
    //  frags in splitFrags will be less than the length of
    //  the path in this unitg.  If so, rebuild this unitig.
    //
    if (splitFragsLen != unitig->dovetail_path_ptr->size()) {

      if (verbose)
        fprintf(stderr, "Rebuild unitig %d\n", unitig->id());

      delete unitig->dovetail_path_ptr;

      unitig->dovetail_path_ptr = new DoveTailPath;

      int splitOffset = -MIN(splitFrags[0].position.bgn, splitFrags[0].position.end);

      //  This should already be true, but we force it still
      splitFrags[0].contained = 0;

      for (int i=0; i<splitFragsLen; i++)
        unitig->addFrag(splitFrags[i], splitOffset, verbose);
    }

    delete [] splitFrags;
    splitFrags    = NULL;
  }  //  End of discontinuity splitting
}



//  True if interval a contains interval b.
//
bool contains(SeqInterval a, SeqInterval b) {
  int aMin,aMax,bMin,bMax;
  if (isReverse(a)) { aMin = a.end; aMax = a.bgn; }
  else              { aMin = a.bgn; aMax = a.end; }
  if (isReverse(b)) { bMin = b.end; bMax = b.bgn; }
  else              { bMin = b.bgn; bMax = b.end; }
  if (aMin <= bMin && aMax >= bMax)
    return true;
  else
    return false;
}

//  Returns the intersection of intervals a and b.
//
SeqInterval intersection(SeqInterval a, SeqInterval b) {
  SeqInterval retVal = NULL_SEQ_LOC;
  int aMin,aMax,bMin,bMax;
  if (isReverse(a)) { aMin = a.end; aMax = a.bgn; }
  else              { aMin = a.bgn; aMax = a.end; }
  if (isReverse(b)) { bMin = b.end; bMax = b.bgn; }
  else              { bMin = b.bgn; bMax = b.end; }

  if (aMax < bMin || bMax < aMin)
    return retVal;

  // so now aMax > bMin && bMax > aMin, thus intersection
  retVal.bgn = MAX(aMin, bMin);
  retVal.end = MIN(aMax, bMax);
  return retVal;
}


vector<SeqInterval> *
findPeakBad(int32 *badGraph, int tigLen, int badMateBreakThreshold) {
  vector<SeqInterval> *peakBads = new vector<SeqInterval>;
  SeqInterval          peak     = {0, 0};
  int32                peakBad  = 0;

  for (int32 i=0; i<tigLen; i++) {
    if (badGraph[i] <= badMateBreakThreshold) {
      //  We are below the bad threshold, start a new bad region, or extend an existing one.

      if (badGraph[i] < peakBad) {
        //  Reset the bad region, we found one that is worse.
        peakBad  = badGraph[i];
        peak.bgn = i;
        peak.end = i;
      }

      if (badGraph[i] <= peakBad)
        //  Extend the bad region into this base.
        peak.end = i;

    } else {
      //  Else, we are above the bad threshold, save any existing bad region and reset.
      if (peakBad < 0) {
        peakBads->push_back(peak);

        peakBad  = 0;
        peak.bgn = 0;
        peak.end = 0;
      }
    }
  }

  //  If there is still a bad region on the stack, save it too.

  if (peakBad < 0)
    peakBads->push_back(peak);

  return peakBads;
}


// hold over from testing if we should use 5' or 3' for range generation, now must use 3'
UnitigBreakPoints* MateChecker::computeMateCoverage(Unitig* tig, BestOverlapGraph* bog_ptr, int badMateBreakThreshold) {
  int tigLen = tig->getLength();

  bool verbose = false;

  MateLocation         positions(_fi, tig, _globalStats);
  vector<SeqInterval> *fwdBads = findPeakBad(positions.badFwdGraph, tigLen, badMateBreakThreshold);
  vector<SeqInterval> *revBads = findPeakBad(positions.badRevGraph, tigLen, badMateBreakThreshold);

  vector<SeqInterval>::const_iterator fwdIter = fwdBads->begin();
  vector<SeqInterval>::const_iterator revIter = revBads->begin();

  UnitigBreakPoints* breaks = new UnitigBreakPoints();

  uint32 backBgn; // Start position of final backbone unitig
  DoveTailNode backbone = tig->getLastBackboneNode(); //tig->getLastBackboneNode(backBgn);
  backBgn = isReverse(backbone.position) ? backbone.position.end :
    backbone.position.bgn ;

  bool combine = false;
  int32 currBackbonePredecessorEnd = 0;
  int32 currBackboneEnd = 0;
  int32 lastBreakBBEnd = 0;

  DoveTailConstIter tigIter = tig->dovetail_path_ptr->begin();
  // Go through the peak bad ranges looking for reads to break on
  while(fwdIter != fwdBads->end() || revIter != revBads->end()) {
    bool isFwdBad = false;
    SeqInterval bad;
    if (revIter == revBads->end() ||
         fwdIter != fwdBads->end() &&  *fwdIter < *revIter) {
      // forward bad group, break at 1st frag
      isFwdBad = true;
      bad = *fwdIter;
      fwdIter++;
      if (lastBreakBBEnd >= bad.bgn) {
        // Skip, instead of combine trying to detect in combine case
        if (verbose)
          fprintf(stderr,"Skip fwd bad range %d %d due to backbone %d\n",
                  bad.bgn, bad.end, lastBreakBBEnd);
        continue;
      }
    } else {                     // reverse bad group, break at last frag
      bad = *revIter;
      if (lastBreakBBEnd >= bad.bgn) {
        // Skip, instead of combine trying to detect in combine case
        if (verbose)
          fprintf(stderr,"Skip rev bad range %d %d due to backbone %d\n",
                  bad.bgn, bad.end, lastBreakBBEnd);
        revIter++;
        continue;
      }
      if (fwdIter != fwdBads->end()) {
        if (fwdIter->bgn < bad.end && bad.end - fwdIter->bgn > 500) {
          // if fwd and reverse bad overlap
          // and end of reverse is far away, do fwd 1st
          isFwdBad = true;
          bad = *fwdIter;
          fwdIter++;
        } else {
          // check for containment relations and skip them
          if (fwdIter->bgn >= bad.bgn && fwdIter->end <= bad.end) {
            if (verbose)
              fprintf(stderr,"Skip fwd bad range %d %d due to contained in rev %d %d\n",
            		  fwdIter->bgn, fwdIter->end, bad.bgn, bad.end);

        	fwdIter++;
            if (fwdIter == fwdBads->end()) {
              continue;
            }
          } else if (bad.bgn >= fwdIter->bgn && bad.end <= fwdIter->end) {
            if (verbose)
              fprintf(stderr,"Skip rev bad range %d %d due to contained in fwd %d %d\n",
                      bad.bgn, bad.end, fwdIter->bgn, fwdIter->end);
      	    revIter++;
            continue;
          }

          if (fwdIter->bgn < bad.end &&
               fwdIter->end > bad.end &&
               bad.end - fwdIter->end < 200) {
            if (verbose)
              fprintf(stderr,"Combine bad ranges %d - %d with %d - %d\n",
                      bad.bgn, bad.end, fwdIter->bgn, fwdIter->end);
            if (bad.bgn == 0) { // ignore reverse at start of tig
              bad.bgn = fwdIter->bgn;
              bad.end = fwdIter->end;
            } else {
              bad.bgn = bad.end;
              bad.end = fwdIter->bgn;
            }
            fwdIter++;
            combine = true;
          }
          revIter++;
        }
      } else {
        revIter++;
      }
    }

    if (verbose)
      fprintf(stderr,"Bad peak from %d to %d\n",bad.bgn,bad.end);

    for(;tigIter != tig->dovetail_path_ptr->end(); tigIter++) {
      DoveTailNode frag = *tigIter;
      SeqInterval loc = frag.position;
      if (isReverse(loc)) { loc.bgn = frag.position.end; loc.end = frag.position.bgn; }

      // Don't want to go past range and break in wrong place
      assert(loc.bgn <= bad.end+1 || loc.end <= bad.end+1);

      // keep track of current and previous uncontained contig end
      // so that we can split apart contained reads that don't overlap each other
      if (!bog_ptr->isContained(frag.ident)) {
        currBackbonePredecessorEnd = currBackboneEnd;
        currBackboneEnd = MAX(loc.bgn, loc.end);
      }

      bool breakNow = false;
      MateLocationEntry mloc = positions.getById(frag.ident);

      if (mloc.mleFrgID1 != 0 && mloc.isGrumpy) { // only break on bad mates
    	  if (isFwdBad && bad.bgn <= loc.end) {
			breakNow = true;
        } else if (!isFwdBad && (loc.bgn >= /* bad.bgn*/ bad.end) ||
                    (combine && loc.end >  bad.bgn) ||
                    (combine && loc.end == bad.end)) {
        	breakNow = true;
        } else if (bad.bgn > backBgn) {
          // fun special case, keep contained frags at end of tig in container
          // instead of in their own new tig where they might not overlap
        	breakNow = true;
        }
      }

      bool noContainerToBreak = false;
      bool incrementToNextFragment = true;
      if (breakNow) {
        if (bog_ptr->isContained(frag.ident)) {
	      // try to find a subsequent uncontained
          while (tigIter != tig->dovetail_path_ptr->end()) {
	        if (bog_ptr->isContained(tigIter->ident) != 0) {
	          tigIter++;
	        } else {
	          break;
	        }
	      }
	      // we couldn't find a downstream frag take upstream one then
	      if (tigIter == tig->dovetail_path_ptr->end()) {
	        tigIter--;
	        while (tigIter != tig->dovetail_path_ptr->begin()) {
	          if (bog_ptr->isContained(tigIter->ident) != 0) {
	            tigIter--;
	          } else {
		       break;
	          }
	        }
	        if (tigIter == tig->dovetail_path_ptr->begin()) {
	          noContainerToBreak = true;
	        } else {
              currBackboneEnd = currBackbonePredecessorEnd;
            }
          }
        }
        if (noContainerToBreak == false) {
          combine = false;
          lastBreakBBEnd = currBackboneEnd;
          if (verbose)
            fprintf(stderr,"Frg to break in peak bad range is %d fwd %d pos (%d,%d) backbone %d\n",
                  frag.ident, isFwdBad, loc.bgn, loc.end, currBackboneEnd);
          uint32 fragEndInTig = THREE_PRIME;
          // If reverse mate is 1st and overlaps its mate break at 5'
          if (mloc.mleUtgID2 == tig->id() && isReverse(loc) &&
               !isReverse(mloc.mlePos2) && loc.bgn >= mloc.mlePos2.bgn)
            fragEndInTig = FIVE_PRIME;

          // either adjust the break position to be before the container so the container travels together with containees
          if (fragEndInTig ==  FIVE_PRIME && !isReverse(tigIter->position)||
              fragEndInTig== THREE_PRIME && isReverse(tigIter->position)) {
          	  // do nothing we are breaking before the current fragment which is our container
        	  incrementToNextFragment = false;
          } else {
        	  // move one back we are breaking after the current fragment
        	  tigIter--;
          }
          frag = *tigIter;
          loc = frag.position;

          UnitigBreakPoint bp(frag.ident, fragEndInTig);
          bp.fragPos = frag.position;
          bp.inSize = 100000;
          bp.inFrags = 10;
          breaks->push_back(bp);
	    }
      }

      if (lastBreakBBEnd != 0 && lastBreakBBEnd > MAX(loc.bgn,loc.end)) {

        DoveTailConstIter nextPos = tigIter;
        nextPos++;

        if (nextPos != tig->dovetail_path_ptr->end()) {

          if (contains(loc, nextPos->position)) {
            // Contains the next one, so skip it
          } else {
            SeqInterval overlap = intersection(loc, nextPos->position);
            int diff = abs(overlap.end - overlap.bgn);

            //  No overlap between this and the next
            //  frag, or the overlap is tiny, or this
            //  frag is contained, but not contained
            //  in the next frag; Break after this
            //  frag.
            //
            if ((NULL_SEQ_LOC == overlap) ||
                (diff < DEFAULT_MIN_OLAP_LEN) ||
                (bog_ptr->isContained(frag.ident) && !bog_ptr->containHaveEdgeTo(frag.ident, nextPos->ident))) {

              uint32 fragEndInTig = THREE_PRIME;
              if (isReverse(loc))
                fragEndInTig = FIVE_PRIME;

              UnitigBreakPoint bp(frag.ident, fragEndInTig);
              bp.fragPos = loc;
              bp.inSize = 100001;
              bp.inFrags = 11;
              breaks->push_back(bp);
              if (verbose)
                fprintf(stderr,"Might make frg %d singleton, end %d size %d pos %d,%d\n",
                        frag.ident, fragEndInTig, breaks->size(),loc.bgn,loc.end);
            }
          }
        }
      }
      if (breakNow) { // Move to next breakpoint
        if (incrementToNextFragment) {
          tigIter++;  // make sure to advance past curr frg
        }
        break;
      }
    }
  }
  delete fwdBads;
  delete revBads;
  return breaks;
}


void
MateLocation::buildTable(Unitig *utg) {

  for (uint32 fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
    DoveTailNode  *frag = &(*utg->dovetail_path_ptr)[fi];

    if (_fi->mateIID(frag->ident) == 0)
      //  Not mated.
      continue;

    uint32  mid = _fi->mateIID(frag->ident);

    if (_iidToTableEntry.find(mid) == _iidToTableEntry.end()) {
      //  We didn't find the mate in the _table, so we know that we haven't seen
      //  either pair, and can create a new entry.
      //
      MateLocationEntry  mle;

      mle.mleFrgID1 = frag->ident;
      mle.mlePos1   = frag->position;
      mle.mleUtgID1 = utg->id();

      mle.mleFrgID2 = 0;
      mle.mlePos2   = NULL_SEQ_LOC;
      mle.mleUtgID2 = 0;

      mle.isGrumpy  = false;

      _iidToTableEntry[frag->ident] = _table.size();
      _table.push_back(mle);

    } else {
      //  Found the mate in the table.  Use that entry.
      //
      uint32  tid = _iidToTableEntry[mid];

      _table[tid].mleFrgID2 = frag->ident;
      _table[tid].mlePos2   = frag->position;
      _table[tid].mleUtgID2 = utg->id();

      _iidToTableEntry[frag->ident] = tid;
    }
  }

  std::sort(_table.begin(), _table.end());

  for (uint32 i=0; i<_table.size(); i++) {
    _iidToTableEntry[_table[i].mleFrgID1] = i;
    _iidToTableEntry[_table[i].mleFrgID2] = i;
  }
}


void
MateLocation::buildHappinessGraphs(Unitig *utg, DistanceCompute *globalStats) {

  for (uint32 mleidx=0; mleidx<_table.size(); mleidx++) {
    MateLocationEntry &loc = _table[mleidx];

    uint32 lib =  _fi->libraryIID(loc.mleFrgID1);

    if (globalStats[lib].samples < 10)
      // Don't check libs that we didn't generate good stats for
      continue;

    int32 badMax = static_cast<int32>(globalStats[lib].mean + 5 * globalStats[lib].stddev);
    int32 badMin = static_cast<int32>(globalStats[lib].mean - 5 * globalStats[lib].stddev);

    int32 dist = 0;
    int32 bgn  = 0;
    int32 end  = 0;

    int32 matBgn = loc.mlePos2.bgn;
    int32 matEnd = loc.mlePos2.end;
    int32 matLen = (matBgn < matEnd) ? (matEnd - matBgn) : (matBgn - matEnd);

    int32 frgBgn = loc.mlePos1.bgn;
    int32 frgEnd = loc.mlePos1.end;
    int32 frgLen = (frgBgn < frgEnd) ? (frgEnd - frgBgn) : (frgBgn - frgEnd);

    if ((matLen >= badMax) ||
        (frgLen >= badMax))
      //  Yikes, fragment longer than insert size!
      continue;


    //  Until reset, assume this is a bad mate pair.
    loc.isGrumpy = true;


    //  If the mate is in another unitig, mark bad only if there is enough space to fit the mate in
    //  this unitig.
    if (loc.mleUtgID1 != loc.mleUtgID2) {
      if ((isReverse(loc.mlePos1) == true)  && (badMax < frgBgn))
        goto markBad;

      if ((isReverse(loc.mlePos1) == false) && (badMax < _tigLen - frgBgn))
        goto markBad;

      //  Not enough space.  Not a grumpy mate pair.
      loc.isGrumpy = false;
      continue;
    }


    //  Both mates are in this unitig.


    //  Same orientation?
    if ((isReverse(loc.mlePos1) == false) &&
        (isReverse(loc.mlePos2) == false))
      goto markBad;
    if ((isReverse(loc.mlePos1) == true) &&
        (isReverse(loc.mlePos2) == true))
      goto markBad;


    //  Check a special case for a circular unitig, outtie mates, but close enough to the end to
    //  plausibly be linking the ends together.
    //
    //   <---              --->
    //  ========unitig==========
    //
    if ((isReverse(loc.mlePos1) == true) &&
        (badMin <= frgBgn + _tigLen - matBgn) &&
        (frgBgn + _tigLen - matBgn <= badMax)) {
      loc.isGrumpy = false;  //  IT'S GOOD, kind of.
      continue;
    }


    //  Outties?  True if pos1.end < pos2.bgn.  (For the second case, swap pos1 and pos2)
    //
    //  (pos1.end) <------   (pos1.bgn)
    //  (pos2.bgn)  -------> (pos2.end)
    //
    if ((isReverse(loc.mlePos1) == true)  && (loc.mlePos1.end < loc.mlePos2.bgn))
      goto markBad;
    if ((isReverse(loc.mlePos1) == false) && (loc.mlePos2.end < loc.mlePos1.bgn))
      goto markBad;

    //  So, now not NORMAL or ANTI or OUTTIE.  We must be left with innies.

    if (isReverse(loc.mlePos1) == false)
      //  First fragment is on the left, second is on the right.
      dist = loc.mlePos2.bgn - loc.mlePos1.bgn;
    else
      //  First fragment is on the right, second is on the left.
      dist = loc.mlePos1.bgn - loc.mlePos2.bgn;

    assert(dist >= 0);

    if ((badMin <= dist) && (dist <= badMax)) {
      bgn = frgBgn;
      end = matBgn;
      incrRange(goodGraph, 2, bgn, end);
      loc.isGrumpy = false;  //  IT'S GOOD!
      continue;
    }

  markBad:

    //  Mark bad from the 3' end of the fragment till the upper limit where the mate should go.

    if (loc.mleUtgID1 == utg->id()) {
      assert(loc.mleFrgID1 != 0);
      if (isReverse(loc.mlePos1) == false) {
        //  Mark bad for forward fagment 1
        bgn = frgEnd;
        end = frgBgn + badMax;
        incrRange(badFwdGraph, -1, bgn, end);
      } else {
        //  Mark bad for reverse fragment 1
        bgn = frgEnd - badMax;
        end = frgEnd;
        incrRange(badRevGraph, -1, bgn, end);
      }
    }

    if (loc.mleUtgID2 == utg->id()) {
      assert(loc.mleFrgID2 != 0);
      if (isReverse(loc.mlePos2) == false) {
        //  Mark bad for forward fragment 2
        bgn = matEnd;
        end = matBgn + badMax;
        incrRange(badFwdGraph, -1, bgn, end);
      } else {
        //  Mark bad for reverse fragment 2
        bgn = matEnd - badMax;
        end = matEnd;
        incrRange(badRevGraph, -1, bgn, end);
      }
    }
  }  //  Over all MateLocationEntries in the table
}  //  buildHappinessGraph()

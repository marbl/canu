
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

#include <math.h>
#include "AS_BOG_MateChecker.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_OVL_overlap.h"
#include "AS_UTL_alloc.h"



MateChecker::MateChecker(FragmentInfo *fi) {

  for(int i=1; i<=fi->numLibraries(); i++) {
    DistanceCompute dc;

    dc.mean      = fi->mean(i);
    dc.stddev    = fi->stddev(i);
    dc.numPairs  = fi->numMatesInLib(i);

    _globalStats[ i ] = dc;
  }

  _fi = fi;
}

MateChecker::~MateChecker() {
}

    
IntervalList* findPeakBad( std::vector<short>* badGraph, int tigLen);



void MateChecker::checkUnitigGraph(UnitigGraph& tigGraph, int badMateBreakThreshold) {

  fprintf(stderr, "==> STARTING MATE BASED SPLITTING.\n");

  tigGraph.checkUnitigMembership();

  computeGlobalLibStats( tigGraph );

  //  Move unhappy contained fragments before attempting mate based
  //  splitting.  We know they're broken already, and if scaffolder
  //  can put them back, we'll let it.
  //
  fprintf(stderr, "==> MOVE CONTAINS #1\n");
  moveContains(tigGraph);

  tigGraph.reportOverlapsUsed("overlaps.aftermatecheck1");
  tigGraph.checkUnitigMembership();

  //  This should do absolutely nothing.  If it does, something is
  //  broken.  By just ejecting unhappy contains, nothing should be
  //  disconnected.
  //
  fprintf(stderr, "==> SPLIT DISCONTINUOUS #1\n");
  splitDiscontinuousUnitigs(tigGraph);

  tigGraph.reportOverlapsUsed("overlaps.aftermatecheck2");
  tigGraph.checkUnitigMembership();

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

  //  The splitting code above is not smart enough to move contained
  //  fragments with containers.  This leaves unitigs disconnected.
  //  We break those unitigs here.
  //
  fprintf(stderr, "==> SPLIT DISCONTINUOUS #2\n");
  splitDiscontinuousUnitigs(tigGraph);

  tigGraph.reportOverlapsUsed("overlaps.aftermatecheck4");
  tigGraph.checkUnitigMembership();

  //  But now, all the splitting probably screwed up happiness of
  //  contained fragments, of left some unhappy fragments in a unitig
  //  that just lost the container.
  //
  fprintf(stderr, "==> MOVE CONTAINS #2\n");
  moveContains(tigGraph);

  tigGraph.reportOverlapsUsed("overlaps.aftermatecheck5");
  tigGraph.checkUnitigMembership();

  //  Do one last check for disconnected unitigs.
  //
  fprintf(stderr, "==> SPLIT DISCONTINUOUS #3\n");
  splitDiscontinuousUnitigs(tigGraph);

  tigGraph.reportOverlapsUsed("overlaps.aftermatecheck6");
  tigGraph.checkUnitigMembership();
}



LibraryStats* MateChecker::computeLibraryStats(Unitig* tig) {
  IdMap goodMates;
  iuid max = std::numeric_limits<iuid>::max(); // Sentinel value
  int numInTig = 0;
  int numWithMate = 0;
  LibraryStats *libs = new LibraryStats();
  std::vector<iuid> otherUnitig; // mates to another unitig
  IdMap otherTigHist; // count of mates to the other unitigs
  // Check each frag in unitig for it's mate
  DoveTailConstIter tigIter = tig->dovetail_path_ptr->begin();
  for(;tigIter != tig->dovetail_path_ptr->end(); tigIter++) {
    DoveTailNode frag = *tigIter;

    if (_fi->mateIID(frag.ident) != 0) {
      numWithMate++;
      if ( tig->id() == Unitig::fragIn( _fi->mateIID(frag.ident) ) ) {
        numInTig++;

        if (goodMates.find( frag.ident ) != goodMates.end()) {
          // already seen it's mate
          if (isReverse(frag.position)) {
            if ( goodMates[frag.ident] == max )
              continue;
            // 2nd frag seen is reverse, so good orientation
            long mateDist = frag.position.bgn - goodMates[frag.ident];
            goodMates[frag.ident] = mateDist;
            // Record the good distance for later stddev recalculation
            iuid distId = _fi->libraryIID(frag.ident);
            if (_dists.find( distId ) == _dists.end()) {
              DistanceList newDL;
              _dists[ distId ] = newDL;
            } 
            _dists[ distId ].push_back( mateDist );
                            
            // Record the distance for unitig local stddev
            if (libs->find( distId ) == libs->end()) {
              DistanceCompute d;
              d.numPairs   = 1;
              d.sumDists   = mateDist;
              (*libs)[ distId ] = d;
            } else {
              (*libs)[distId].numPairs++;
              (*libs)[distId].sumDists+=mateDist;
            }
          } else {
            // 2nd frag seen is forward, so bad
            goodMates[frag.ident] = max;
          }
        } else { // 1st of pair
          if (isReverse(frag.position)) {
            // 1st reversed, so bad
            goodMates[_fi->mateIID(frag.ident)] = max;
          } else {
            // 1st forward, so good store begin
            goodMates[_fi->mateIID(frag.ident)] = frag.position.bgn;
          }
        }
      } else {
        // mate in other unitig, just create histogram currently
        otherUnitig.push_back( frag.ident );
        iuid otherTig = Unitig::fragIn( _fi->mateIID(frag.ident) );
        if (otherTigHist.find( otherTig ) == otherTigHist.end())
          otherTigHist[ otherTig ] = 1;
        else
          otherTigHist[ otherTig ]++;
      }
    }
  }

  // Calculate the unitig local mean
  LibraryStats::iterator dcIter = libs->begin();
  for(; dcIter != libs->end(); dcIter++) {
    iuid lib = dcIter->first;
    DistanceCompute *dc = &(dcIter->second);
    dc->mean = dc->sumDists / dc->numPairs;
    //fprintf(stderr,"Distance lib %ld has %ld pairs with mean dist %.1f\n",
    //        lib, dc->numPairs, dc->mean );
  }
  // Sum of (x-mean)^2, not yet full stddev
  IdMapConstIter goodIter = goodMates.begin();
  for(;goodIter != goodMates.end(); goodIter++) {
    iuid fragId = goodIter->first;
    iuid mateDist = goodIter->second;
    if (mateDist == max)
      continue;
    iuid distId = _fi->libraryIID(fragId);
    (*libs)[distId].sumSquares+=pow(mateDist-(*libs)[distId].mean,2);
  }
  // Calculate the real stddev
  dcIter = libs->begin();
  for(; dcIter != libs->end(); dcIter++) {
    iuid lib = dcIter->first;
    DistanceCompute *dc = &(dcIter->second);
    // really need to just collect all values and calculate in checkUnitigGraph()
    if (dc->numPairs == 1)
      dc->stddev = 0.0;
    else
      dc->stddev = sqrt( dc->sumSquares / (dc->numPairs-1) );
    //fprintf(stderr,"Distance lib %ld has %ld pairs with stddev %.1f\n",
    //        lib, dc->numPairs, dc->stddev );
  }
  return libs;
}

    
void MateChecker::computeGlobalLibStats( UnitigGraph& tigGraph ) {
  LibraryStats::iterator dcIter;
  _dists.clear(); // reset to seperate multiple Graphs
  UnitigsConstIter tigIter = tigGraph.unitigs->begin();
  for(; tigIter != tigGraph.unitigs->end(); tigIter++) {
    if (*tigIter == NULL ) 
      continue;
    LibraryStats* libs = computeLibraryStats(*tigIter);
    // Accumulate per unitig stats to compute global stddev's
    for(dcIter = libs->begin(); dcIter != libs->end(); dcIter++) {
      iuid lib = dcIter->first;
      DistanceCompute dc = dcIter->second;
      if (_globalStats.find(lib) == _globalStats.end() ) {
        _globalStats[ lib ] = dc;
      }
      else {
        DistanceCompute *gdc = &(_globalStats[ lib ]);
        gdc->numPairs   += dc.numPairs;
        gdc->sumDists   += dc.sumDists;
        gdc->sumSquares += dc.sumSquares;
      }
    }
    delete libs; // Created in computeLibraryStats
  }
  // Calculate and output overall global mean
  for(dcIter= _globalStats.begin(); dcIter != _globalStats.end(); dcIter++){
    iuid lib = dcIter->first;
    DistanceCompute *dc = &(dcIter->second);
    dc->mean = dc->sumDists / dc->numPairs;
    fprintf(stderr,"Distance lib %ld has global %ld pairs with mean dist %.1f\n",
            lib, dc->numPairs, dc->mean );
  }
  // Calculate and output overall global stddev
  for(dcIter= _globalStats.begin(); dcIter != _globalStats.end(); dcIter++) {
    iuid lib = dcIter->first;
    DistanceCompute *dc = &(dcIter->second);
    if (dc->numPairs == 1)
      dc->stddev = 0.0;
    else
      dc->stddev = sqrt( dc->sumSquares / (dc->numPairs-1) );
    fprintf(stderr,"Distance lib %ld has global %ld pairs with stddev %.1f\n",
            lib, dc->numPairs, dc->stddev );
    // Now reset the global stats to zero so the real calculation below works
    dc->numPairs   = 0;
    dc->mean       = 0.0;
    dc->stddev     = 0.0;
    dc->sumSquares = 0.0;
    dc->sumDists   = 0.0;
  }
  // Tab delimited table header
  fprintf(stderr,
          "DistLib\tnumDists\tmedian\t1/3rd\t2/3rd\tmaxDiff\tmin\tmax\tnumGood\tmean\tstddev\n");

  // Disregard outliers and recalculate global stddev
  LibDistsConstIter libDistIter = _dists.begin();
  for(; libDistIter != _dists.end(); libDistIter++) {
    iuid libId = libDistIter->first;
    DistanceList dl = libDistIter->second;
    sort(dl.begin(),dl.end());
    int size = dl.size();
    int median   = dl[ size / 2 ];
    int third    = dl[ size / 3 ];
    int twoThird = dl[ size * 2 / 3 ];
    int aproxStd = MAX( median - third, twoThird - median);
    int biggest  = median + aproxStd * 5;
    int smallest = median - aproxStd * 5;

    // now go through the distances and calculate the real stddev
    // including everything within 5 stddevs
    iuid numBad = 0;
    DistanceCompute *gdc = &(_globalStats[ libId ]);
    DistanceListCIter dIter = dl.begin();
    for(;dIter != dl.end(); dIter++) {
      if (*dIter >= smallest && *dIter <= biggest ) {
        gdc->numPairs++;
        gdc->sumDists += *dIter;
      } else {
        numBad++;
      } 
    }
    gdc->mean = gdc->sumDists / gdc->numPairs;
    // Compute sum of squares for stddev
    for(dIter = dl.begin(); dIter != dl.end(); dIter++) {
      if (*dIter >= smallest && *dIter <= biggest )
        gdc->sumSquares += pow( *dIter - gdc->mean, 2);
    }
    if (gdc->numPairs > 1)
      gdc->stddev = sqrt( gdc->sumSquares / (gdc->numPairs-1) );

    //  Should use AS_IID, but we don't know it!
    fprintf(stderr, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.1f\n",
            libId, size, median, third, twoThird, aproxStd, smallest, biggest,
            gdc->numPairs, gdc->mean, gdc->stddev );
  }
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

    MateLocation positions(_fi);

    positions.buildTable(thisUnitig);
    delete positions.buildHappinessGraphs(thisUnitig->getLength(), _globalStats);

    DoveTailNode         *frags         = new DoveTailNode [thisUnitig->dovetail_path_ptr->size()];
    int                   fragsLen      = 0;

    bool                  verbose       = true;

    if (verbose)
      fprintf(stderr, "moveContain unitig %d\n", thisUnitig->id());

    for (DoveTailIter fragIter = thisUnitig->dovetail_path_ptr->begin();
         fragIter != thisUnitig->dovetail_path_ptr->end();
         fragIter++) {

      BestContainment   *bestcont   = tigGraph.bog_ptr->getBestContainer(fragIter->ident);
      MateLocationEntry  mloc       = positions.getById(fragIter->ident);

      if        ((fragIter->contained == 0) || (bestcont == NULL)) {
        //  CASE 1:  Not contained.  Leave the fragment here.

#if 0
        fprintf(stderr, "moveContain1 f=%d c=%d(%d b=%d(%d)\n",
                fragIter->ident,
                fragIter->contained,
                thisUnitig->fragIn(fragIter->contained),
                (bestcont == NULL) ? 0 : bestcont->container,
                (bestcont == NULL) ? 0 : thisUnitig->fragIn(bestcont->container));
#endif

        if (fragIter->contained)
          fprintf(stderr, "WARNING: frag %d is contained with no bestcontainer?!\n", fragIter->ident);
        assert(fragIter->contained == 0);

        //  Update the contained field in the current unitig, then
        //  copy the fragment to the list -- if we need to rebuild the
        //  unitig, the list is used, otherwise, we have already made
        //  the change needed.
        //
        fragIter->contained = 0;

        frags[fragsLen] = *fragIter;
        fragsLen++;
      } else if ((thisUnitig->fragIn(fragIter->contained) == thisUnitig->id()) ||
                 (thisUnitig->fragIn(bestcont->container) == thisUnitig->id())) {
        //  CASE 2: Frag is contained, and his container is in this
        //  unitig.  Leave here, UNLESS he is an unhappy mate, in
        //  which case we eject him to a new unitig.

#if 0
        fprintf(stderr, "moveContain2 f=%d c=%d(%d b=%d(%d)\n",
                fragIter->ident,
                fragIter->contained,
                thisUnitig->fragIn(fragIter->contained),
                (bestcont == NULL) ? 0 : bestcont->container,
                (bestcont == NULL) ? 0 : thisUnitig->fragIn(bestcont->container));
#endif

        //  Is he mated?  If not, we'll always leave here.
        //
        //  Is the mate here?  If not, we can leave here.
        //  (We've since decided to eject both, but here's the test)
        //    (thisUnitig->fragIn(_fi->mateIID(fragIter->ident)) == thisUnitig->fragIn(fragIter->ident)) &&
        //
        //  id1 != 0 -> we found the fragment in the mate happiness table
        //  isBad -> and the mate is unhappy.
        //
        //  What's id1 vs id2 in MateLocationEntry?  Dunno.  All I
        //  know is that if there is no mate present, one of those
        //  will be 0.  (Similar test used above too.)
        //
        //  See comment on similar if test below.
        //
        if ((_fi->mateIID(fragIter->ident) > 0) &&
            (mloc.id1 != 0) &&
            (mloc.id2 != 0) &&
            (mloc.isBad == true)) {
          Unitig        *sing = new Unitig(verbose);

          if (verbose)
            fprintf(stderr, "Ejecting unhappy contained fragment %d from unitig %d into new unitig %d\n",
                    fragIter->ident, thisUnitig->id(), sing->id());

          fragIter->contained = 0;

          sing->addFrag(*fragIter, -MIN(fragIter->position.bgn, fragIter->position.end), verbose);
          tigGraph.unitigs->push_back(sing);
          thisUnitig = (*tigGraph.unitigs)[ti];  //  Reset the pointer; unitigs might be reallocated
        } else {
          frags[fragsLen] = *fragIter;
          fragsLen++;
        }

      } else {
        //  CASE 3: Frag is contained, and his container is not here.
        //  If happy, leave him here if there is an overlap to some
        //  other fragment, but mark as not contained.
        //
        //  If not happy, or no overlap found above, move to the
        //  unitig of his container.

#if 0
        fprintf(stderr, "moveContain3 f=%d c=%d(%d b=%d(%d) mloc %d %d\n",
                fragIter->ident,
                fragIter->contained,
                thisUnitig->fragIn(fragIter->contained),
                (bestcont == NULL) ? 0 : bestcont->container,
                (bestcont == NULL) ? 0 : thisUnitig->fragIn(bestcont->container),
                mloc.id1, mloc.isBad);
#endif

        bool  hasOverlap = false;

        //  See comment on similar if test above.
        //
        if ((_fi->mateIID(fragIter->ident) > 0) &&
            (mloc.id1 != 0) &&
            (mloc.id2 != 0) &&
            (mloc.isBad == false)) {

          //  Mate is happy.  Look for an overlap.

          if (fragsLen == 0) {
            //  Nothing before?  Check fragments after if there is an
            //  overlap.  In rare cases, we eject the container, and the
            //  containee had no overlap to fragments remaining.
            //
            //  This isn't bullet proof.  One of those overlaps might be
            //  to a contained fragment that we eject later on.  So, we
            //  only accept overlaps to non-contained fragments.
            //
            //  And we only need to check fragments that overlap in
            //  the layout.
            //
            for (DoveTailIter ft=fragIter;
                 ((hasOverlap == false) && 
                  (ft != thisUnitig->dovetail_path_ptr->end()) &&
                  (MAX(fragIter->position.bgn, fragIter->position.end) < MIN(ft->position.bgn, ft->position.end)));
                 ft++)
              if (tigGraph.bog_ptr->isContained(ft->ident) == false)
                hasOverlap = tigGraph.bog_ptr->containHaveEdgeTo(fragIter->ident, ft->ident);
           
          } else {
            //  Not the first fragment.  Make sure this guy overlaps
            //  something alreay placed.
            //
            for (int ff=0; (hasOverlap == false) && (ff<fragsLen); ff++)
              hasOverlap = tigGraph.bog_ptr->containHaveEdgeTo(fragIter->ident, frags[ff].ident);
          }
        }

        if (hasOverlap) {
          //  Overlap exists to another fragment.  Add it to this unitig.
          //
          if (verbose)
            fprintf(stderr, "Removing containment relationship for fragment %d to keep it in unitig %d\n",
                    fragIter->ident, thisUnitig->id());

          //  See the comment with above about what this line does.
          //  Search backwards for the same line.
          fragIter->contained = 0;

          frags[fragsLen] = *fragIter;
          fragsLen++;
        } else {
          //  Mate is not happy, or is happy but has no overlap
          //  anymore, or the fragment is not mated.  Eject to
          //  container unitig (if not mated) or singleton (if mated).

#warning DANGEROUS assume unitig is at id-1 in vector
          Unitig         *thatUnitig = (*tigGraph.unitigs)[thisUnitig->fragIn(bestcont->container) - 1];
          DoveTailNode    containee  = *fragIter;  //  To make it read/write

          //  If fragment is mated, move to a singleton, instead of with container.
          //
          if (_fi->mateIID(containee.ident) != 0) {
            Unitig        *sing = new Unitig(verbose);

            if (verbose)
              fprintf(stderr, "Ejecting unhappy contained fragment %d from unitig %d into new unitig %d\n",
                      containee.ident, thisUnitig->id(), sing->id());

            containee.contained = 0;

            sing->addFrag(containee, -MIN(containee.position.bgn, containee.position.end), verbose);
            tigGraph.unitigs->push_back(sing);
            thisUnitig = (*tigGraph.unitigs)[ti];  //  Reset the pointer; unitigs might be reallocated
          } else {
            assert(thatUnitig->id() == thisUnitig->fragIn(bestcont->container));

            if (verbose)
              fprintf(stderr, "Moving contained fragment %d from unitig %d to be with its container %d in unitig %d\n",
                      containee.ident, thisUnitig->id(), bestcont->container, thatUnitig->id());

            //  Make sure it's marked as contained.
            containee.contained = bestcont->container;

            thatUnitig->addContainedFrag(containee, bestcont, verbose);
          }
        }  // end of moving unhappy contained to container unitig
      }  //  end of if blocks deciding what to do with a containee
    }  //  over all frags

    //  Now, rebuild this unitig if we made changes.

    if (fragsLen != thisUnitig->dovetail_path_ptr->size()) {
      if (verbose)
        fprintf(stderr, "Rebuild unitig %d after removing contained fragments.\n", thisUnitig->id());

      delete thisUnitig->dovetail_path_ptr;

      thisUnitig->dovetail_path_ptr = new DoveTailPath;

      //  No need to resort.  Offsets only need adjustment if
      //  the first fragment is thrown out.  If not,
      //  splitOffset will be zero.
      //
      int splitOffset = -MIN(frags[0].position.bgn, frags[0].position.end);

      //  This is where we clean up from the splitting not
      //  dealing with contained fragments -- we force the
      //  first frag to be uncontained.
      //
      frags[0].contained = 0;

      for (int i=0; i<fragsLen; i++)
        thisUnitig->addFrag(frags[i], splitOffset, verbose);
    }

    delete [] frags;
    frags = NULL;
    
  }  //  Over all unitigs
}


//  After splitting and ejecting some contains, check for discontinuous unitigs.
//
void MateChecker::splitDiscontinuousUnitigs(UnitigGraph& tigGraph) {

  bool  verbose = true;

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

      //  We require at least 10bp of overlap between fragments.  If
      //  we don't have that, split off the fragments we've seen.
      //
      if (maxEnd - 10 < MIN(fragIter->position.bgn, fragIter->position.end)) {

        //  If there is exactly one fragment, and it's contained, and
        //  it's not mated, move it to the container.  (This has a
        //  small positive benefit over just making every read a
        //  singleton).
        //
        if ((splitFragsLen == 1) &&
            (_fi->mateIID(splitFrags[0].ident) == 0) &&
            (splitFrags[0].contained != 0)) {

#warning DANGEROUS assume unitig is at id-1 in vector
          Unitig           *dangler  = (*tigGraph.unitigs)[unitig->fragIn(splitFrags[0].contained) - 1];
          BestContainment  *bestcont = tigGraph.bog_ptr->getBestContainer(splitFrags[0].ident);

          if (verbose)
            fprintf(stderr, "Dangling contained fragment in unitig %d -> move them to container unitig %d\n", unitig->id(), dangler->id());

          dangler->addContainedFrag(splitFrags[0], bestcont, verbose);
        } else {
          Unitig *dangler = new Unitig(verbose);

          if (verbose)
            fprintf(stderr, "Dangling fragments in unitig %d -> move them to unitig %d\n", unitig->id(), dangler->id());

          int splitOffset = -MIN(splitFrags[0].position.bgn, splitFrags[0].position.end);

          //  This should already be true, but we force it still
          splitFrags[0].contained = 0;

          for (int i=0; i<splitFragsLen; i++)
            dangler->addFrag(splitFrags[i], splitOffset, verbose);

          splitFragsLen = 0;

          tigGraph.unitigs->push_back(dangler);
          unitig = (*tigGraph.unitigs)[ti];
        }
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


void incrRange( std::vector<short>* graph, short val, iuid n, iuid m ) {
  if (n == m)
    return;
  int sz = graph->size();
  assert( m > n );
  assert( n < sz );
  assert( m <= sz );
  if (n < 0) n = 0;
  if (m >= sz) m = sz-1;

  for(iuid i=n; i <=m ; i++)
    graph->at(i) += val;
}

//  True if interval a contains interval b.
//
bool contains( SeqInterval a, SeqInterval b) {
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
SeqInterval intersection( SeqInterval a, SeqInterval b) {
  SeqInterval retVal = NULL_SEQ_LOC;
  int aMin,aMax,bMin,bMax;
  if (isReverse(a)) { aMin = a.end; aMax = a.bgn; }
  else              { aMin = a.bgn; aMax = a.end; }
  if (isReverse(b)) { bMin = b.end; bMax = b.bgn; }
  else              { bMin = b.bgn; bMax = b.end; }

  if (aMax < bMin || bMax < aMin)
    return retVal;

  // so now aMax > bMin && bMax > aMin, thus intersection
  retVal.bgn = MAX( aMin, bMin );
  retVal.end = MIN( aMax, bMax );
  return retVal;
}

// Assumes list is already sorted
void combineOverlapping( IntervalList* list ) {
  IntervalList::iterator iter = list->begin();
  IntervalList::iterator a = iter++;
  for(; iter != list->end() && iter != static_cast<IntervalList::iterator>(NULL);
      iter++) {
    SeqInterval aIb = intersection( *a, *iter );
    if (!(aIb == NULL_SEQ_LOC) && aIb.end - aIb.bgn > 1000) {
      a->bgn = aIb.bgn;
      a->end = aIb.end;
      list->erase( iter ); 
    }
  }
}


IntervalList* findPeakBad(std::vector<short>* badGraph, int tigLen, int badMateBreakThreshold) {
  IntervalList* peakBads = new IntervalList();
  SeqInterval   peak = NULL_SEQ_LOC;
  int badBegin, peakBad, lastBad;
  peakBad = lastBad = badBegin = 0;
  for(int i=0; i < tigLen; i++) {
    if( badGraph->at( i ) <= badMateBreakThreshold ) {
      if (badBegin == 0)  // start bad region
        badBegin = i;
      if(badGraph->at(i) < peakBad) {
        peakBad   = badGraph->at(i);
        peak.bgn = peak.end = i;
      } else if (lastBad < 0 && lastBad == peakBad) {
        peak.end = i-1;
      }
      lastBad = badGraph->at(i);
    } else {
      if (badBegin > 0) {  // end bad region
        peakBads->push_back( peak );
        peakBad = lastBad = badBegin = 0;
        peak = NULL_SEQ_LOC;
      }
    }
  }
  return peakBads;
}


// hold over from testing if we should use 5' or 3' for range generation, now must use 3'
UnitigBreakPoints* MateChecker::computeMateCoverage(Unitig* tig, BestOverlapGraph* bog_ptr, int badMateBreakThreshold) {
  int tigLen = tig->getLength();

  MateLocation positions(_fi);
  positions.buildTable( tig );
  MateCounts *unused = positions.buildHappinessGraphs( tigLen, _globalStats );
  delete unused;

  IntervalList *fwdBads = findPeakBad( positions.badFwdGraph, tigLen, badMateBreakThreshold );
  IntervalList *revBads = findPeakBad( positions.badRevGraph, tigLen, badMateBreakThreshold );

  UnitigBreakPoints* breaks = new UnitigBreakPoints();

  iuid backBgn; // Start position of final backbone unitig
  DoveTailNode backbone = tig->getLastBackboneNode(backBgn);
  backBgn = isReverse( backbone.position ) ? backbone.position.end :
    backbone.position.bgn ;

  bool combine = false;
  CDS_COORD_t currBackboneEnd = 0;
  CDS_COORD_t lastBreakBBEnd = 0;
  IntervalList::const_iterator fwdIter = fwdBads->begin();
  IntervalList::const_iterator revIter = revBads->begin();
  DoveTailConstIter tigIter = tig->dovetail_path_ptr->begin();
  // Go through the peak bad ranges looking for reads to break on
  while( fwdIter != fwdBads->end() || revIter != revBads->end() ) {
    bool isFwdBad = false;
    SeqInterval bad;
    if ( revIter == revBads->end() ||
         fwdIter != fwdBads->end() &&  *fwdIter < *revIter ) {
      // forward bad group, break at 1st frag
      isFwdBad = true;
      bad = *fwdIter;
      fwdIter++;
      if (lastBreakBBEnd >= bad.bgn) {
        // Skip, instead of combine trying to detect in combine case
        fprintf(stderr,"Skip fwd bad range %d %d due to backbone %d\n",
                bad.bgn, bad.end, lastBreakBBEnd);
        continue;
      }
    } else {                     // reverse bad group, break at last frag
      bad = *revIter;
      if (lastBreakBBEnd >= bad.bgn) {
        // Skip, instead of combine trying to detect in combine case
        fprintf(stderr,"Skip rev bad range %d %d due to backbone %d\n",
                bad.bgn, bad.end, lastBreakBBEnd);
        revIter++;
        continue;
      }
      if (fwdIter != fwdBads->end()) {
        if ( fwdIter->bgn < bad.end && bad.end - fwdIter->bgn > 500 ) {
          // if fwd and reverse bad overlap 
          // and end of reverse is far away, do fwd 1st
          isFwdBad = true;
          bad = *fwdIter;
          fwdIter++;
        } else {
          if ( fwdIter->bgn < bad.end &&
               fwdIter->end > bad.end &&
               bad.end - fwdIter->end < 200) {
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

    fprintf(stderr,"Bad peak from %d to %d\n",bad.bgn,bad.end);

    for(;tigIter != tig->dovetail_path_ptr->end(); tigIter++) {
      DoveTailNode frag = *tigIter;
      SeqInterval loc = frag.position;

      // Don't want to go past range and break in wrong place
      assert( loc.bgn <= bad.end+1 || loc.end <= bad.end+1 );

      // keep track of current and previous uncontained contig end
      // so that we can split apart contained reads that don't overlap each other
      if ( !bog_ptr->isContained(frag.ident) )
        currBackboneEnd = MAX(loc.bgn, loc.end);

      bool breakNow = false;
      MateLocationEntry mloc = positions.getById( frag.ident );

      if (mloc.id1 != 0 && mloc.isBad) { // only break on bad mates
        if ( isFwdBad && bad.bgn <= loc.end ) {
          breakNow = true;
        } else if ( !isFwdBad && (loc.bgn >= bad.end) ||
                    (combine && loc.end >  bad.bgn) ||
                    (combine && loc.end == bad.end) ) {
          breakNow = true;
        } else if (bad.bgn > backBgn) {
          // fun special case, keep contained frags at end of tig in container 
          // instead of in their own new tig where they might not overlap
          breakNow = true;
        }
      }

      if (breakNow) {
        combine = false;
        lastBreakBBEnd = currBackboneEnd;
        fprintf(stderr,"Frg to break in peak bad range is %d fwd %d pos (%d,%d) backbone %d\n",
                frag.ident, isFwdBad, loc.bgn, loc.end, currBackboneEnd );
        fragment_end_type fragEndInTig = THREE_PRIME;
        // If reverse mate is 1st and overlaps its mate break at 5'
        if ( mloc.unitig2 == tig->id() && isReverse( loc ) &&
             !isReverse(mloc.pos2) && loc.bgn >= mloc.pos2.bgn )
          fragEndInTig = FIVE_PRIME;

        UnitigBreakPoint bp( frag.ident, fragEndInTig );
        bp.fragPos = frag.position;
        bp.inSize = 100000;
        bp.inFrags = 10;
        breaks->push_back( bp );
      }

      if ( lastBreakBBEnd != 0 && lastBreakBBEnd > MAX(loc.bgn,loc.end)) {

        DoveTailConstIter nextPos = tigIter;
        nextPos++;

        if (nextPos != tig->dovetail_path_ptr->end()) {

          if ( contains( loc, nextPos->position ) ) {
            // Contains the next one, so skip it
          } else {
            SeqInterval overlap = intersection(loc, nextPos->position);
            int diff = abs( overlap.end - overlap.bgn);

            //  No overlap between this and the next
            //  frag, or the overlap is tiny, or this
            //  frag is contained, but not contained
            //  in the next frag; Break after this
            //  frag.
            //
            if ((NULL_SEQ_LOC == overlap) || 
                (diff < DEFAULT_MIN_OLAP_LEN) || 
                (bog_ptr->isContained( frag.ident ) && !bog_ptr->containHaveEdgeTo( frag.ident, nextPos->ident))) {

              fragment_end_type fragEndInTig = THREE_PRIME;
              if (isReverse( loc ))
                fragEndInTig = FIVE_PRIME;

              UnitigBreakPoint bp( frag.ident, fragEndInTig );
              bp.fragPos = loc;
              bp.inSize = 100001;
              bp.inFrags = 11;
              breaks->push_back( bp );
              fprintf(stderr,"Might make frg %d singleton, end %d size %d pos %d,%d\n",
                      frag.ident, fragEndInTig, breaks->size(),loc.bgn,loc.end);
            }
          }
        }
      }
      if (breakNow) { // Move to next breakpoint
        tigIter++;  // make sure to advance past curr frg
        break;
      }
    }
  }
  delete fwdBads;
  delete revBads;
  return breaks;
}


void MateLocation::buildTable( Unitig *tig) {

  for(DoveTailConstIter frag = tig->dovetail_path_ptr->begin();
      frag != tig->dovetail_path_ptr->end();
      frag++) {

    if ( _fi->mateIID(frag->ident) != NULL_FRAG_ID ) {
      if (hasFrag( _fi->mateIID(frag->ident) ) )
        addMate( tig->id(), frag->ident, frag->position );
      else
        startEntry( tig->id(), frag->ident, frag->position );
    }
  }
  sort();
}


MateCounts* MateLocation::buildHappinessGraphs( int tigLen, LibraryStats& globalStats ) {
  goodGraph->resize( tigLen+1 );
  badFwdGraph->resize( tigLen+1 );
  badRevGraph->resize( tigLen+1 );

  MateCounts *cnts = new MateCounts();

  for(MateLocIter  posIter  = begin(); posIter != end(); posIter++) {
    MateLocationEntry loc = *posIter;
    iuid fragId         =  loc.id1;
    iuid mateId         =  _fi->mateIID(fragId);
    iuid lib            =  _fi->libraryIID(fragId);
    cnts->total++;
    DistanceCompute *gdc = &(globalStats[ lib ]);
    // Don't check libs that we didn't generate good stats for
    if (gdc->numPairs < 10)
      continue;
    int badMax = static_cast<int>(gdc->mean + 5 * gdc->stddev);
    int badMin = static_cast<int>(gdc->mean - 5 * gdc->stddev);
    int frgBgn = loc.pos1.bgn;
    int frgEnd = loc.pos1.end;
    int frgLen = abs(frgEnd - frgBgn);
    if (frgLen >= badMax) {
      fprintf(stderr,"Warning skipping read %d with length %d > mean + 5*stddev %d\n",
              fragId, frgLen, badMax );
      continue; // Could assert instead
    }
    frgLen = abs( loc.pos2.end - loc.pos2.bgn );
    if (frgLen >= badMax) {
      fprintf(stderr,"Warning skipping read %d with length %d > mean + 5*stddev %d\n",
              loc.id2, frgLen, badMax );
      continue; // Could assert instead
    }
    if ( loc.unitig1 != loc.unitig2) {
      // mate in another tig, mark bad only if max range exceeded
      if (isReverse(loc.pos1)) {
        if ( frgBgn > badMax ) {
          incrRange( badRevGraph, -1, frgBgn - badMax, frgEnd );
          posIter->isBad = true;
#if 0
          fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                  fragId, frgBgn, loc.pos1.end, mateId, lib);
#endif
          cnts->badOtherTig++;
        } else {
          cnts->otherTig++;
        }
      } else {
        if ( frgBgn + badMax < tigLen ) {
          incrRange( badFwdGraph, -1, frgEnd, frgBgn + badMax );
          posIter->isBad = true;
#if 0
          fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                  fragId, frgBgn, frgEnd, mateId, lib);
#endif
          cnts->badOtherTig++;
        } else {
          cnts->otherTig++;
        }
      }
    } else {
      // both mates in this unitig
      int mateBgn =  loc.pos2.bgn;
      int mateEnd =  loc.pos2.end;
      if (isReverse( loc.pos1 )) {
        if (!isReverse( loc.pos2 )) {
          // reverse and forward, check for circular unitig
          int dist = frgBgn + tigLen - mateBgn; 
          if ( dist <= badMax && dist >= badMin) {
            cnts->goodCircular++;
            continue; // Good circular mates
          } 
        }
        // 1st reversed, so bad 
        iuid beg = MAX( 0, frgBgn - badMax );
        incrRange( badRevGraph, -1, beg, frgEnd);
        posIter->isBad = true;
#if 0
        fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                fragId, frgBgn, frgEnd, mateId, lib);
#endif

        if (isReverse( loc.pos2 )) {
          // 2nd mate is reversed, so mark bad towards tig begin
          beg = MAX( 0, mateBgn - badMax );
          incrRange( badRevGraph, -1, beg, mateEnd);
          cnts->badAntiNormal++;
        } else {
          // 2nd mate is forward, so mark bad towards tig end
          iuid end = MIN( tigLen, mateBgn + badMax );
          incrRange( badFwdGraph, -1, mateEnd, end);
          cnts->badOuttie++;
        }
#if 0
        fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                mateId, mateBgn, mateEnd, fragId, lib);
#endif

      } else {
        // 1st forward
        if (isReverse( loc.pos2 )) {
          // 2nd reverse so good orient, check distance
          uint16 mateLen = mateBgn - mateEnd;
          int mateDist = mateBgn - frgBgn;  

          if (mateDist >= badMin && mateDist <= badMax) {
            // For good graph only we mark from 5' to 5'
            // so overlapping mates can still be good
            incrRange(goodGraph,2, frgBgn, mateBgn);
            cnts->good++;
          }
          else {
            // both are bad, mate points towards tig begin
            iuid beg = MAX(0, mateBgn - badMax); 
            iuid end = mateEnd;

            incrRange(badRevGraph, -1, beg, end);
            posIter->isBad = true;
#if 0
            fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                    mateId, mateBgn, mateEnd, fragId, lib);
#endif

            end = MIN( tigLen, frgBgn + badMax );
            beg = frgEnd;

            incrRange(badFwdGraph,-1, beg, end);
#if 0
            fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                    fragId, frgBgn, frgEnd, mateId, lib);
#endif
            cnts->badInnie++;
          }
        } else {
          // 1st and 2nd forward so both bad 
          iuid end = MIN( tigLen, frgBgn + badMax );
          iuid beg = frgEnd;

          incrRange(badFwdGraph,-1, beg, end);
          posIter->isBad = true;

#if 0
          fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                  fragId, frgBgn, frgEnd, mateId, lib);
#endif

          // 2nd mate is forward, so mark bad towards tig end
          end = MIN( tigLen, mateBgn + badMax );
          beg = mateEnd;

          incrRange( badFwdGraph, -1, beg, end);
#if 0
          fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                  mateId, mateBgn, mateEnd, fragId, lib);
#endif
          cnts->badNormal++;
        }
      }
    }
  }
  return cnts;
}


bool MateLocation::startEntry(iuid unitigID, iuid fragID, SeqInterval fragPos) {
  if ( _iidIndex.find( fragID) != _iidIndex.end() )
    return false; // Entry already exists, can't start new

  assert( fragID != NULL_FRAG_ID );
  MateLocationEntry entry;
  entry.id1     = fragID;
  entry.pos1    = fragPos;
  entry.unitig1 = unitigID; 
  entry.id2     = NULL_FRAG_ID;
  entry.pos2    = NULL_SEQ_LOC;
  entry.unitig2 = NULL_FRAG_ID;; 
  entry.isBad   = false;

  _table.push_back( entry );
  _iidIndex[ fragID ] = _table.size()-1;
  return true;
}


bool MateLocation::addMate(iuid unitigId, iuid fragId, SeqInterval fragPos) {
  iuid mateId = _fi->mateIID(fragId);
  IdMapConstIter entryIndex = _iidIndex.find( mateId );
  if ( _iidIndex.find( fragId ) != _iidIndex.end() ||
       entryIndex == _iidIndex.end() )
    return false; // Missing mate or already added

  iuid idx = entryIndex->second;
  _table[ idx ].id2      = fragId;
  _table[ idx ].pos2     = fragPos;
  _table[ idx ].unitig2  = unitigId;
  _iidIndex[ fragId ]    = idx;
  return true;
}


MateLocationEntry MateLocation::getById( iuid fragId ) {
  IdMapConstIter entryIndex = _iidIndex.find( fragId );
  if ( entryIndex == _iidIndex.end() )
    return NULL_MATE_ENTRY;
  else
    return _table[ entryIndex->second ];
}


bool MateLocation::hasFrag(iuid fragId) {
  if ( _iidIndex.find( fragId ) == _iidIndex.end() )
    return false;
  else
    return true;
}


void MateLocation::sort() {
  std::sort(begin(),end());
  MateLocCIter iter = begin();
  int i = 0;
  for(; iter != end(); iter++, i++) {
    MateLocationEntry entry = *iter;
    _iidIndex[ entry.id1 ] = i;
    _iidIndex[ entry.id2 ] = i;
  }
}

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

/* RCS info
 * $Id: AS_BOG_MateChecker.cc,v 1.45 2008-03-21 22:05:28 brianwalenz Exp $
 * $Revision: 1.45 $
 */

#include <math.h>
#include "AS_BOG_MateChecker.hh"
#include "AS_OVL_overlap.h"  // For DEFAULT_MIN_OLAP_LEN
#include "AS_UTL_alloc.h"    // For safe_calloc

namespace AS_BOG{
    MateChecker::~MateChecker() {
    }

    iuid MateChecker::readStore(const char* gkpStorePath) {
        GateKeeperStore *gkpStore = openGateKeeperStore(gkpStorePath, FALSE);

        iuid numDists = getNumGateKeeperLibraries(gkpStore);
        iuid i;
        for(i = 1; i <= numDists; i++){
            GateKeeperLibraryRecord  *gkpl = getGateKeeperLibrary(gkpStore, i);
            DistanceCompute dc;
            dc.mean   = gkpl->mean;
            dc.stddev = gkpl->stddev;
            _globalStats[ i ] = dc;
        }

        StreamStruct *frags = openStream(gkpStore->frg);
        resetStream(frags, STREAM_FROMSTART, STREAM_UNTILEND);

        GateKeeperFragmentRecord gkpf;
        iuid frgIID = 1;
        while(nextStream(frags, &gkpf, 0, NULL)){
            if (gkpf.mateIID != 0) {
                MateInfo mi;
                mi.mate = gkpf.mateIID;
                mi.lib  = gkpf.libraryIID;
                _mates[ gkpf.readIID ] = mi;

                _globalStats[ mi.lib ].numPairs++;
            }
            frgIID++;
        }
        for(i = 1; i <= numDists; i++){
            DistanceCompute *dc = &(_globalStats[i]);
            fprintf(stderr, "Lib %d mean %.1f stddev %.1f numReads %d\n",
                    i, dc->mean, dc->stddev, dc->numPairs );

            assert( dc->numPairs % 2 == 0);
            dc->numPairs /= 2;
        }

        iuid numFrgs = getNumGateKeeperFragments( gkpStore );

        std::cerr << "Frg count " << frgIID << " num in store " << numFrgs << std::endl;

        closeStream(frags);
        closeGateKeeperStore(gkpStore);
        return numFrgs;
    }

    ///////////////////////////////////////////////////////////////////////////
    
    IntervalList* findPeakBad( std::vector<short>* badGraph, int tigLen);

    ///////////////////////////////////////////////////////////////////////////

    MateInfo MateChecker::getMateInfo(iuid fragID) {
        if (_mates.find( fragID ) == _mates.end())
            return NULL_MATE_INFO;
        else
            return _mates[ fragID ];
    }

    ///////////////////////////////////////////////////////////////////////////

    LibraryStats* MateChecker::computeLibraryStats(Unitig* tig) {
        fprintf(stderr,"Check mates for tig %ld\n",tig->id());
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
            iuid fragId = frag.ident;
            MateInfo mateInfo = getMateInfo( fragId );
            if ( mateInfo.mate != NULL_FRAG_ID ) {
                numWithMate++;
                if ( tig->id() == Unitig::fragIn( mateInfo.mate ) ) {
                    numInTig++;
                    SeqInterval fragPos = frag.position;
                    if (goodMates.find( fragId ) != goodMates.end()) {
                        // already seen it's mate
                        if (isReverse(fragPos)) {
                            if ( goodMates[fragId] == max )
                                continue;
                            // 2nd frag seen is reverse, so good orientation
                            long mateDist = fragPos.bgn - goodMates[fragId];
                            goodMates[fragId] = mateDist;
                            // Record the good distance for later stddev recalculation
                            iuid distId = mateInfo.lib;
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
                            goodMates[fragId] = max;
                        }
                    } else { // 1st of pair
                        if (isReverse(fragPos)) {
                            // 1st reversed, so bad
                            goodMates[mateInfo.mate] = max;
                        } else {
                            // 1st forward, so good store begin
                            goodMates[mateInfo.mate] = fragPos.bgn;
                        }
                    }
                } else {
                    // mate in other unitig, just create histogram currently
                    otherUnitig.push_back( fragId );
                    iuid otherTig = Unitig::fragIn( mateInfo.mate );
                    if (otherTigHist.find( otherTig ) == otherTigHist.end())
                        otherTigHist[ otherTig ] = 1;
                    else
                        otherTigHist[ otherTig ]++;
                }
            }
        }
        if (tig->dovetail_path_ptr->size() > 1 ) {
            fprintf(stderr,"Num frags in tig %ld\n",tig->dovetail_path_ptr->size());
            fprintf(stderr,"Num frags with mate %d\n",numWithMate);
            fprintf(stderr,"Num with mate in unitig %d\n",numInTig);
            fprintf(stderr,"Num other unitig %d\n",otherUnitig.size());
            IdMapConstIter histIter = otherTigHist.begin();
            for(;histIter != otherTigHist.end(); histIter++) {
                iuid libId = histIter->first;
                iuid cnt   = histIter->second;
                fprintf(stderr,"Num mates to unitig %ld is %ld.\n",libId,cnt);
            }
        }

        // Calculate the unitig local mean
        LibraryStats::iterator dcIter = libs->begin();
        for(; dcIter != libs->end(); dcIter++) {
            iuid lib = dcIter->first;
            DistanceCompute *dc = &(dcIter->second);
            dc->mean = dc->sumDists / dc->numPairs;
            fprintf(stderr,"Distance lib %ld has %ld pairs with mean dist %.1f\n",
                    lib, dc->numPairs, dc->mean );
        }
        // Sum of (x-mean)^2, not yet full stddev
        IdMapConstIter goodIter = goodMates.begin();
        for(;goodIter != goodMates.end(); goodIter++) {
            iuid fragId = goodIter->first;
            iuid mateDist = goodIter->second;
            if (mateDist == max)
                continue;
            iuid distId = getMateInfo( fragId ).lib;
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
            fprintf(stderr,"Distance lib %ld has %ld pairs with stddev %.1f\n",
                    lib, dc->numPairs, dc->stddev );
        }
        return libs;
    }

    ///////////////////////////////////////////////////////////////////////////
    
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

            fprintf(stderr, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%.1f\t%.1f\n",
                    libId, size, median, third, twoThird, aproxStd, smallest, biggest,
                    gdc->numPairs, gdc->mean, gdc->stddev );
        }
    }
    ///////////////////////////////////////////////////////////////////////////
    // main entry point into mate checking code
    void MateChecker::checkUnitigGraph( UnitigGraph& tigGraph ) {
        if ( ! BogOptions::useGkpStoreLibStats )
            computeGlobalLibStats( tigGraph );

        int numSplits = 1;
        int iterNum = 1;
        int prevNumSplits = 0;
        while (numSplits > 0) {
            UnitigVector* splits = new UnitigVector(); // save split unitigs
            UnitigsIter tigEditIter = tigGraph.unitigs->begin();
            for(; tigEditIter != tigGraph.unitigs->end(); tigEditIter++) {
                if (*tigEditIter == NULL || (*tigEditIter)->getNumFrags() < 2 ) 
                    continue;

                FragmentEnds* breaks = computeMateCoverage( *tigEditIter, tigGraph.bog_ptr );
                tigGraph.accumulateSplitUnitigs( tigEditIter, breaks, splits );
                delete breaks;
            }
            numSplits = splits->size();
            fprintf(stderr,"Num mate based splits %d iteration %d\n",numSplits,iterNum++);
            tigGraph.unitigs->insert( tigGraph.unitigs->end(),splits->begin(),splits->end());
            delete splits;

            // bugs in splitting can prevent numSplits from reaching zero
            if (numSplits == prevNumSplits)
                break;
            else
                prevNumSplits = numSplits;
        }


        //  After splitting, check for discontinuous unitigs.  What
        //  occasionally happens:
        //
        //  ---------1
        //    ----------2
        //       ----3
        //            -------4
        //
        //  We split at fragment 2.  We're not sure if fragment 3
        //  should be contained in fragment 2 -- it's mate could be
        //  later on (say, fragment 6).  BPW has seen this example,
        //  and fragment 3 is also not labeled as contained in 2 (at
        //  least, tigGraph.bog_ptr->isContained(3) == false and/or
        //  tigGraph.bog_ptr->containHaveEdgeTo(3, 2) == false).
        //
        //  After splitting, we have #1 and #2 in a new unitig, and #3
        //  and #4 in another unitig.  #3 and #4 are not overlapping,
        //  resulting in consensus failures.
        //
        for (UnitigsIter unitigIter = tigGraph.unitigs->begin();
             unitigIter != tigGraph.unitigs->end();
             unitigIter++) {

            Unitig  *unitig = *unitigIter;

            if ((unitig == NULL) ||
                (unitig->dovetail_path_ptr->empty()) ||
                (unitig->dovetail_path_ptr->size() == 1))
                continue;


            //  The original flavor of this code (currently in the
            //  #if0 below) would eject the first fragment if it had
            //  containment issues.  Perhaps we should still do that.
            //
            //  This has all sorts of problems, like what to do with
            //  contained fragments, and was the inspiration for this
            //  current block of code (the example in the comment at
            //  the start.


            //  Check for discontinuities

            DoveTailConstIter     fragIter = unitig->dovetail_path_ptr->begin();
            int                   maxEnd   = 0;

            DoveTailNode         *splitFrags    = new DoveTailNode [unitig->dovetail_path_ptr->size()];
            int                   splitFragsLen = 0;
            int                   splitOffset   = 0;

            while (fragIter != unitig->dovetail_path_ptr->end()) {

                //  True only if this is the first frag, or we just split.
                if (splitFragsLen == 0) {
                    maxEnd      = MAX(fragIter->position.bgn, fragIter->position.end);
                    splitOffset = MIN(fragIter->position.bgn, fragIter->position.end);
                }

                //  We require at least 10bp of overlap between
                //  fragments.  If we don't have that, split off the
                //  fragments we've seen.
                //
                if (maxEnd - 10 < MIN(fragIter->position.bgn, fragIter->position.end)) {
                    Unitig *dangler = new Unitig;

                    fprintf(stderr, "Dangling Fragments in unitig %d -> move them to unitig %d\n", unitig->id(), dangler->id());

                    for (int i=0; i<splitFragsLen; i++)
                        dangler->addFrag(splitFrags[i], splitOffset);

                    splitFragsLen = 0;

                    tigGraph.unitigs->push_back(dangler);
                }  //  End break

                splitFrags[splitFragsLen++] = *fragIter;

                maxEnd = MAX(maxEnd, MAX(fragIter->position.bgn, fragIter->position.end));

                fragIter++;
            }  //  End of unitig fragment iteration

            //  If we split this unitig, the length of the frags in
            //  splitFrags will be less than the length of the path in
            //  this unitg.  If so, rebuild this unitig.
            //
            if (splitFragsLen != unitig->dovetail_path_ptr->size()) {
                delete unitig->dovetail_path_ptr;

                unitig->dovetail_path_ptr = new DoveTailPath;

                for (int i=0; i<splitFragsLen; i++)
                    unitig->addFrag(splitFrags[i], splitOffset);
            }

            delete [] splitFrags;
        }  //  End of discontinuity splitting


        // Now we'll chuck out the contained frags that have unhappy mates
        // Plus some important cleanup of uncovered reads
        ejectUnhappyContains( tigGraph );
    }

    ///////////////////////////////////////////////////////////////////////////

    void MateChecker::ejectUnhappyContains(UnitigGraph& tigGraph ) {
        MateCounts allMates;
        UnitigVector* singletons = new UnitigVector(); // save singleton unitigs
        UnitigsIter tigEditIter = tigGraph.unitigs->begin();
        for(; tigEditIter != tigGraph.unitigs->end(); tigEditIter++) {
            if (*tigEditIter == NULL ) 
                continue;

            Unitig* tig = *tigEditIter;
            if (tig->dovetail_path_ptr->empty() || tig->dovetail_path_ptr->size() == 1)
                continue;

            int tigLen  = tig->getLength();

            MateLocation positions(this);

            // Build mate position table
            positions.buildTable( tig );

            // Build good and bad mate graphs
            MateCounts* cnts = positions.buildHappinessGraphs( tigLen, _globalStats );
            allMates += *cnts;
            delete cnts;

            // save position of singleton tig for case of a contain contain
            IdMap singleIds;

            // save original path to iterate through, make new one to
            // add reads back into

            DoveTailPath *oldpath = tig->dovetail_path_ptr;
            tig->dovetail_path_ptr = new DoveTailPath;

            DoveTailConstIter tigIter = oldpath->begin();
            int bgnShift = 0;

            //  Check for discontinuities

#if 0
            //  Make sure 1st and 2nd frag overlap, otherwise make 1st
            //  singleton

            bool overlapFound = false;
            while (!overlapFound) {
                if (oldpath->size() > 1 && tigIter != oldpath->end()) {
                    DoveTailNode first  = *(tigIter);
                    DoveTailNode second = *(tigIter+1);

                    BestContainment* cntn = tigGraph.bog_ptr->getBestContainer(second.ident);

                    if ( cntn != NULL && cntn->container == first.ident ) {
                        //  Second is contained in the first.
                        //
                        overlapFound = true;

                    } else if ( tigGraph.bog_ptr->isContained( first.ident ) &&
                                !tigGraph.bog_ptr->containHaveEdgeTo( first.ident, second.ident)) {
                        //  First is contained, and it's not contained in the second.
                        //
                        fprintf(stderr, "Non overlap start %d is contained in %d.\n",
                                first.ident,
                                tigGraph.bog_ptr->getBestContainer(first.ident)->container);

                        Unitig* singleton = new Unitig;
                        singleton->addFrag( first );
                        singleIds[ first.ident ] = singletons->size();

                        bgnShift = -MIN( second.position.bgn, second.position.end );

                        second = *(++tigIter);  //  Nop for second, needed to advance to the next frag.

#if 0
                        //  For any frags contained in the new unitig,
                        //  move them there.  bgnShift (should) never
                        //  change since we're just adding in
                        //  containees.

                        while ((tigGraph.bog_ptr->isContained(second.ident) &&
                                (tig->fragIn(tigGraph.bog_ptr->getBestContainer(second.ident)->container) == singleton->id()))) {
                            singleton->addFrag(second );
                            singleIds[ second.ident ] = singletons->size();
                            second = *(++tigIter);
                        }
#endif
                        singletons->push_back( singleton );
                    } else {
                        //  OK by default, I guess.
                        overlapFound = true;
                    }
                } else {
                    overlapFound = true;
                }
            }
#endif

            //  eject unhappies

            for(; tigIter != oldpath->end(); tigIter++) {
                DoveTailNode frag = *tigIter;
                if (BogOptions::ejectUnhappyContained && frag.contained ) {
                    MateLocationEntry mloc = positions.getById( frag.ident );
                    if (mloc.id1 != 0 && mloc.isBad ) {
                        // make contain bad a singleton unitig
                        Unitig* singleton = new Unitig;
                        singleton->addFrag( frag, bgnShift );
                        singleIds[ frag.ident ] = singletons->size();
                        singletons->push_back( singleton );
                    } else {    
                        // contained, but not bad mate what to do with it?
                        IdMapConstIter found = singleIds.find( frag.contained );
                        if ( found != singleIds.end() ) { // contained in unhappy singleton
                            if (mloc.id1 == 0) { // unmated, so follow container
                                iuid offset = found->second;
                                (*singletons)[ offset ]->addFrag( frag );
                                singleIds[ frag.ident ] = offset;
                                fprintf(stderr,"Contain placed in singleton frg %d off %d tig %d\n",
                                        frag.ident, offset, (*singletons)[ offset ]->id());
                            } else { // happy mate leave in place
                                frag.contained = 0; // Might pick something else
                                tig->addFrag( frag, bgnShift );
                            }
                        } else { // mated, so good mate
                            tig->addFrag( frag, bgnShift );
                        }
                    }
                } else {
                    // non-contained, just add back to dovetail path
                    tig->addFrag( frag, bgnShift );
                }
            }

            delete oldpath;
        }
        tigGraph.unitigs->insert( tigGraph.unitigs->end(), singletons->begin(), singletons->end() );

        std::cerr << allMates;
    }

    ///////////////////////////////////////////////////////////////////////////

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
    ///////////////////////////////////////////////////////////////////////////
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
    ///////////////////////////////////////////////////////////////////////////

    // hold over from testing if we should use 5' or 3' for range generation, now must use 3'
    FragmentEnds* MateChecker::computeMateCoverage( Unitig* tig, BestOverlapGraph* bog_ptr ) {
        int tigLen = tig->getLength();

        MateLocation positions(this);
        // Build mate position table
        positions.buildTable( tig );

        // Build good and bad mate graphs
        MateCounts *unused = positions.buildHappinessGraphs( tigLen, _globalStats );
        delete unused;

        // For debugging purposes output the table
        MateLocIter posIter = positions.begin();
        for(posIter = positions.begin(); posIter != positions.end(); posIter++){
            MateLocationEntry loc = *posIter;
            std::cerr << loc << std::endl;
        }

        // do something with the good and bad graphs
        fprintf(stderr,"Per 300 bases good graph unitig %ld size %ld:\n",tig->id(),tigLen);
        long sum = 0;
        for(int i=0; i < tigLen; i++) {
            if (i > 1 && i % 300 == 0) {
                fprintf(stderr,"%d ", sum / 300);
                sum = 0;
            }
            sum += positions.goodGraph->at( i );
        }
        IntervalList *fwdBads,*revBads;
        fprintf(stderr,"\nPer 300 bases bad fwd graph:\n");
        fwdBads = findPeakBad( positions.badFwdGraph, tigLen );
        fprintf(stderr,"\nPer 300 bases bad rev graph:\n");
        revBads = findPeakBad( positions.badRevGraph, tigLen );

        fprintf(stderr,"Num fwdBads is %d\n",fwdBads->size());
        fprintf(stderr,"Num revBads is %d\n",revBads->size());
        FragmentEnds* breaks = new FragmentEnds(); // return value

        posIter  = positions.begin();

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
                    bp.position = frag.position;
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
                                bp.position = loc;
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

    ///////////////////////////////////////////////////////////////////////////

    IntervalList* findPeakBad( std::vector<short>* badGraph, int tigLen ) {
        long sum = 0;
        for(int i=0; i < tigLen; i++) {
            if (i > 1 && i % 300 == 0) {
                fprintf(stderr,"%d ", sum / 300);
                sum = 0;
            }
            sum += badGraph->at( i );
        }
        fprintf(stderr,"\n");
        IntervalList* peakBads = new IntervalList();
        SeqInterval   peak = NULL_SEQ_LOC;
        int badBegin, peakBad, lastBad;
        peakBad = lastBad = badBegin = 0;
        for(int i=0; i < tigLen; i++) {
            if( badGraph->at( i ) <= BogOptions::badMateBreakThreshold ) {
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
                    fprintf(stderr,"Bad mates >%d from %d to %d peak %d from %d to %d\n",
                            -BogOptions::badMateBreakThreshold,
                            badBegin,i-1,peakBad,peak.bgn,peak.end);
                    peakBads->push_back( peak );
                    peakBad = lastBad = badBegin = 0;
                    peak = NULL_SEQ_LOC;
                }
            }
        }
        return peakBads;
    }

    ///////////////////////////////////////////////////////////////////////////
    
    void MateLocation::buildTable( Unitig *tig) {
        DoveTailConstIter tigIter = tig->dovetail_path_ptr->begin();
        for(;tigIter != tig->dovetail_path_ptr->end(); tigIter++) {
            DoveTailNode frag = *tigIter;
            iuid fragId = frag.ident;
            MateInfo mateInfo = _checker->getMateInfo( fragId );
            iuid mateId = mateInfo.mate;
            if ( mateInfo.mate != NULL_FRAG_ID ) {
                if (hasFrag( mateId ) )
                    addMate( tig->id(), fragId, frag.position );
                else
                    startEntry( tig->id(), fragId, frag.position );
            }
        }
        sort();
    }

    ///////////////////////////////////////////////////////////////////////////
    
    std::ostream& operator<< (std::ostream& os, MateCounts c) {
        int sum = c.badOtherTig + c.otherTig + c.goodCircular + c.good + c.badOuttie +
            c.badInnie + c.badAntiNormal + c.badNormal;

        os << std::endl << "Total mates " << c.total << " should equal sum " << sum
           << std::endl
           << "Good innies " << c.good << " good circular mates " << c.goodCircular
           << std::endl
           << "Other unitig " << c.otherTig << " other tig dist exceeded " << c.badOtherTig
           << std::endl
           << "Bad Outtie " << c.badOuttie << " Bad Innie " << c.badInnie << std::endl
           << "Bad Normal " << c.badNormal << " Bad Anti-normal " << c.badAntiNormal
           << std::endl << std::endl;
    }
    
    ///////////////////////////////////////////////////////////////////////////
    
    MateCounts* MateLocation::buildHappinessGraphs( int tigLen, LibraryStats& globalStats ) {
        goodGraph->resize( tigLen+1 );
        badFwdGraph->resize( tigLen+1 );
        badRevGraph->resize( tigLen+1 );

        MateCounts *cnts = new MateCounts();
        MateLocIter  posIter  = begin();
        for(;        posIter != end(); posIter++) {
            MateLocationEntry loc = *posIter;
            iuid fragId         =  loc.id1;
            MateInfo mateInfo   =  _checker->getMateInfo( fragId );
            iuid mateId         =  mateInfo.mate;
            iuid lib            =  mateInfo.lib;
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
                        fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                fragId, frgBgn, loc.pos1.end, mateId, lib);
                        cnts->badOtherTig++;
                    } else {
                        cnts->otherTig++;
                    }
                } else {
                    if ( frgBgn + badMax < tigLen ) {
                        incrRange( badFwdGraph, -1, frgEnd, frgBgn + badMax );
                        posIter->isBad = true;
                        fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                fragId, frgBgn, frgEnd, mateId, lib);
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
                    fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                            fragId, frgBgn, frgEnd, mateId, lib);

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
                    fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                            mateId, mateBgn, mateEnd, fragId, lib);

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
                            fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                    mateId, mateBgn, mateEnd, fragId, lib);

                            end = MIN( tigLen, frgBgn + badMax );
                            beg = frgEnd;

                            incrRange(badFwdGraph,-1, beg, end);
                            fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                    fragId, frgBgn, frgEnd, mateId, lib);
                            cnts->badInnie++;
                        }
                    } else {
                        // 1st and 2nd forward so both bad 
                        iuid end = MIN( tigLen, frgBgn + badMax );
                        iuid beg = frgEnd;

                        incrRange(badFwdGraph,-1, beg, end);
                        posIter->isBad = true;

                        fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                fragId, frgBgn, frgEnd, mateId, lib);

                        // 2nd mate is forward, so mark bad towards tig end
                        end = MIN( tigLen, mateBgn + badMax );
                        beg = mateEnd;

                        incrRange( badFwdGraph, -1, beg, end);
                        fprintf(stderr,"Bad mate %ld pos %ld %ld mate %ld lib %d\n",
                                mateId, mateBgn, mateEnd, fragId, lib);
                        cnts->badNormal++;
                    }
                }
            }
        }
        return cnts;
    }

    ///////////////////////////////////////////////////////////////////////////

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

    ///////////////////////////////////////////////////////////////////////////
    MateLocation::~MateLocation() {
        delete goodGraph;
        delete badFwdGraph;
        delete badRevGraph;
    }

    ///////////////////////////////////////////////////////////////////////////

    bool MateLocation::addMate(iuid unitigId, iuid fragId, SeqInterval fragPos) {
        iuid mateId = _checker->getMateInfo( fragId ).mate;
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

    ///////////////////////////////////////////////////////////////////////////

    MateLocationEntry MateLocation::getById( iuid fragId ) {
        IdMapConstIter entryIndex = _iidIndex.find( fragId );
        if ( entryIndex == _iidIndex.end() )
            return NULL_MATE_ENTRY;
        else
            return _table[ entryIndex->second ];
    }

    ///////////////////////////////////////////////////////////////////////////

    bool MateLocation::hasFrag(iuid fragId) {
        if ( _iidIndex.find( fragId ) == _iidIndex.end() )
            return false;
        else
            return true;
    }

    ///////////////////////////////////////////////////////////////////////////

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

    ///////////////////////////////////////////////////////////////////////////

    std::ostream& operator<< (std::ostream& os, MateLocationEntry& e) {
        int dist = e.pos2.bgn - e.pos1.bgn;
        os << e.pos1.bgn <<","<< e.pos1.end <<" -> "<< e.pos2.bgn <<","<< e.pos2.end
           << " "<< dist << " Frg " << e.id1 << " mate " << e.id2 << " isBad " << e.isBad ;
        return os;
    }

    ///////////////////////////////////////////////////////////////////////////
}

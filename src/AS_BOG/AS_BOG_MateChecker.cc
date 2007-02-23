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
 * $Id: AS_BOG_MateChecker.cc,v 1.8 2007-02-23 15:44:42 eliv Exp $
 * $Revision: 1.8 $
*/

#include "AS_BOG_MateChecker.hh"

namespace AS_BOG{
    MateChecker::~MateChecker() {
    }

    void MateChecker::readStore(const char* gkpStorePath) {
        GateKeeperStore *gkpStore = openGateKeeperStore(gkpStorePath, FALSE);

        StreamHandle frags = openStream(gkpStore->frg, NULL, 0);
        resetStream(frags, STREAM_FROMSTART, STREAM_UNTILEND);

        GateKeeperFragmentRecord gkpf;
        iuid frgIID = 1;
        while(nextStream(frags, &gkpf)){
            MateInfo mi;
            mi.mate = gkpf.mateIID;
            mi.lib  = gkpf.libraryIID;
            _mates[ gkpf.readIID ] = mi;

            frgIID++;
        }
        closeStream(frags);
        closeGateKeeperStore(gkpStore);
    }

    ///////////////////////////////////////////////////////////////////////////

    MateInfo MateChecker::getMateInfo(iuid fragID)
    {
        if (_mates.find( fragID ) == _mates.end())
            return NULL_MATE_INFO;
        else
            return _mates[ fragID ];
    }

    ///////////////////////////////////////////////////////////////////////////

    LibraryStats* MateChecker::checkUnitig(Unitig* tig)
    {
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
        for(;tigIter != tig->dovetail_path_ptr->end(); tigIter++)
        {
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
                                d.stddev     = 0.0;
                                d.mean       = 0.0;
                                d.sumSquares = 0.0;
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
                            // 1st forward, so good
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
        for(;goodIter != goodMates.end(); goodIter++)
        {
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

    void MateChecker::checkUnitigGraph( UnitigGraph& tigGraph )
    {
        LibraryStats globalStats;
        LibraryStats::iterator dcIter;
        _dists.clear(); // reset to seperate multiple Graphs
        UnitigsConstIter tigIter = tigGraph.unitigs->begin();
        for(; tigIter != tigGraph.unitigs->end(); tigIter++)
        {
            if (*tigIter == NULL ) 
                continue;
            LibraryStats* libs = checkUnitig(*tigIter);
            // Accumulate per unitig stats to compute global stddev's
            for(dcIter = libs->begin(); dcIter != libs->end(); dcIter++) {
                iuid lib = dcIter->first;
                DistanceCompute dc = dcIter->second;
                if (globalStats.find(lib) == globalStats.end() ) {
                    globalStats[ lib ] = dc;
                }
                else {
                    DistanceCompute *gdc = &(globalStats[ lib ]);
                    gdc->numPairs   += dc.numPairs;
                    gdc->sumDists   += dc.sumDists;
                    gdc->sumSquares += dc.sumSquares;
                }
            }
            delete libs; // Created in checkUnitig
        }
        // Calculate and output overall global mean
        for(dcIter= globalStats.begin(); dcIter != globalStats.end(); dcIter++){
            iuid lib = dcIter->first;
            DistanceCompute *dc = &(dcIter->second);
            dc->mean = dc->sumDists / dc->numPairs;
            fprintf(stderr,"Distance lib %ld has global %ld pairs with mean dist %.1f\n",
                    lib, dc->numPairs, dc->mean );
        }
        // Calculate and output overall global stddev
        for(dcIter= globalStats.begin(); dcIter != globalStats.end(); dcIter++)
        {
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
        for(; libDistIter != _dists.end(); libDistIter++)
        {
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
            DistanceCompute *gdc = &(globalStats[ libId ]);
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
            for(dIter = dl.begin(); dIter != dl.end(); dIter++)
            {
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
}

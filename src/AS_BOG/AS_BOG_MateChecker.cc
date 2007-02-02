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
 * $Id: AS_BOG_MateChecker.cc,v 1.1 2007-02-02 21:22:15 eliv Exp $
 * $Revision: 1.1 $
*/

#include "AS_BOG_MateChecker.hh"

namespace AS_BOG{
    MateChecker::~MateChecker() {
    }

    void MateChecker::readStore(const char* gkpStorePath) {
        GateKeeperStore gkpStore;
        InitGateKeeperStore(&gkpStore, gkpStorePath);
        OpenReadOnlyGateKeeperStore(&gkpStore);

        StreamHandle frags = openStream(gkpStore.frgStore,NULL,0);
        resetStream(frags, STREAM_FROMSTART, STREAM_UNTILEND);

        GateKeeperFragmentRecord gkpf;
        iuid frgIID = 1;
        while(nextStream(frags, &gkpf)){
            if(gkpf.numLinks == 1){
                GateKeeperLinkRecordIterator iterator;
                GateKeeperLinkRecord link;
                CreateGateKeeperLinkRecordIterator(gkpStore.lnkStore, gkpf.linkHead,
                                                    frgIID, &iterator);
                while(NextGateKeeperLinkRecordIterator(&iterator, &link)){
//                    fprintf(stderr,"Frg %ld link %ld to %ld dist:%ld type:%c ori:%c\n",
//                    frgIID, link.frag1, link.frag2, link.distance, link.type,
//                        getLinkOrientation( &link ));
                    if (frgIID == link.frag1)
                        mates[ link.frag1 ] = link.frag2;
                    else if (frgIID == link.frag2)
                        mates[ link.frag2 ] = link.frag1;
                    else
                        assert(0);
                    IdMapConstIter dIter = dists.find( link.distance );
                    if (dIter == dists.end())
                        // if not found add end to range
                        dists[ link.distance ] = link.frag2;
                    else if (dIter->second < link.frag2)
                        dists[ link.distance ] = link.frag2;
                }
            } else if (gkpf.numLinks > 1)
                assert(0);//Code doesn't handle multiple mates for a single frag

            frgIID++;
        }
        closeStream(gkpStore.frgStore);
        CloseGateKeeperStore(&gkpStore);
        iuid beginID = 1;
        IdMapConstIter distIter = dists.begin();
        for (; distIter != dists.end(); distIter++) {
            fprintf(stderr,"Dist %ld has iid range %ld to %ld.\n",
                    distIter->first, beginID, distIter->second );
            beginID = distIter->second + 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////

    iuid MateChecker::getMate(iuid fragID)
    {
        if (mates.find( fragID ) == mates.end())
            return NULL_FRAG_ID;
        else
            return mates[ fragID ];
    }

    ///////////////////////////////////////////////////////////////////////////

    iuid MateChecker::getDist(iuid fragID)
    {
        IdMapConstIter iter = dists.begin();
        iuid dist = NULL_FRAG_ID;
        for(;iter != dists.end(); iter++)
        {
            if (fragID < iter->second) {
                dist = iter->first;
                break;
            }
        }
        return dist;
    }

    ///////////////////////////////////////////////////////////////////////////
    struct DistanceCompute {
        double sumDists;
        int numFrags;
    };

    void MateChecker::checkUnitig(Unitig* tig)
    {
        fprintf(stderr,"Check mates for tig %ld\n",tig->id());
        IdMap goodMates;
        iuid max = std::numeric_limits<iuid>::max();
        int numInTig = 0;
        int numWithMate = 0;
        std::vector<iuid> otherUnitig;
        DoveTailConstIter tigIter = tig->dovetail_path_ptr->begin();
        for(;tigIter != tig->dovetail_path_ptr->end(); tigIter++)
        {
            DoveTailNode frag = *tigIter;
            iuid fragId = frag.ident;
            iuid mateId = getMate( fragId );
            if ( mateId != NULL_FRAG_ID ) {
                numWithMate++;
                if ( tig->id() == Unitig::fragIn( mateId ) ) {
                    numInTig++;
                    SeqInterval fragPos = frag.position;
                    if (goodMates.find( fragId ) != goodMates.end()) {
                        // already seen it's mate
                        if (isReverse(fragPos)) {
                            // 2nd frag seen is reverse, so good orientation
                            goodMates[fragId] = fragPos.end - goodMates[fragId];
                        } else {
                            // 2nd frag seen is forward, so bad
                            goodMates[fragId] = max;
                        }
                    } else {
                        // 1st of pair
                        if (isReverse(fragPos)) {
                            // 1st reversed, so bad
                            fprintf(stderr,"Skip bad orient %ld\n",mateId);
                            goodMates[mateId] = max;
                        } else {
                            // 1st forward, so good
                            fprintf(stderr,"Begin good %ld,%ld at %ld\n",fragId,mateId,fragPos.bgn);
                            goodMates[mateId] = fragPos.bgn;
                        }
                    }

                } else {
                    otherUnitig.push_back( fragId );
                }
            }
        }
        fprintf(stderr,"Num frags in tig %ld\n",tig->dovetail_path_ptr->size());
        fprintf(stderr,"Num frags with mate %d\n",numWithMate);
        fprintf(stderr,"Num with mate in unitig %d\n",numInTig);
        fprintf(stderr,"Num other unitig %d\n",otherUnitig.size());
        std::map<iuid,DistanceCompute> libs;
        IdMapConstIter goodIter = goodMates.begin();
        for(;goodIter != goodMates.end(); goodIter++)
        {
            iuid fragId = goodIter->first;
            iuid mateDist = goodIter->second;
            if (mateDist == max)
                continue;
            iuid distId = getDist( fragId );
            if (libs.find( distId ) == libs.end()) {
                DistanceCompute d;
                d.numFrags = 1;
                d.sumDists = mateDist;
            } else {
                libs[distId].numFrags++;
                libs[distId].sumDists+=mateDist;
            }
        }
        std::map<iuid,DistanceCompute>::const_iterator dcIter = libs.begin();
        for(; dcIter != libs.end(); dcIter++) {
            iuid lib = dcIter->first;
            DistanceCompute dc = dcIter->second;
            fprintf(stderr,"Distance lib %ld has %ld frags with avg dist %.1f\n",
                    lib, dc.numFrags, dc.sumDists / dc.numFrags );
        }
    }
    void MateChecker::checkUnitigGraph( UnitigGraph& tigGraph )
    {
        UnitigsConstIter tigIter = tigGraph.unitigs->begin();
        for(; tigIter != tigGraph.unitigs->end(); tigIter++)
        {
            if (*tigIter != NULL )
                checkUnitig(*tigIter);
        }
    }
}

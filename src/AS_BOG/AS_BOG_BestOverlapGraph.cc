
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
/*************************************************
* Module:  AS_BOG_BestOverlapGraph.cc
* Description:
*	Data structure to contain the best overlaps and containments
*	based on a defined metric.
* 
*    Programmer:  K. Li
*       Started:  1 August 2005
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: AS_BOG_BestOverlapGraph.cc,v 1.12 2005-09-07 18:03:47 eliv Exp $
 * $Revision: 1.12 $
*/

static const char CM_ID[] = "$Id: AS_BOG_BestOverlapGraph.cc,v 1.12 2005-09-07 18:03:47 eliv Exp $";

//  System include files
#include<iostream>

#include "AS_BOG_BestOverlapGraph.hh"
//#include "AS_BOG_BestOverlapGraphVisitor.hh"

extern "C" {
#include "AS_PER_fragStore.h"
}
namespace AS_BOG{

    fragment_end_type BestOverlapGraph::AEnd(const Long_Olap_Data_t& olap) {
        if (olap.a_hang < 0 && olap.b_hang < 0)
            return FIVE_PRIME;
        if (olap.a_hang > 0 && olap.b_hang > 0)
            return THREE_PRIME;

        assert(0); // no contained
    }
    fragment_end_type BestOverlapGraph::BEnd(const Long_Olap_Data_t& olap) {
        if (olap.a_hang < 0 && olap.b_hang < 0)
            if ( olap.flipped )
                return FIVE_PRIME;
            else
                return THREE_PRIME;

        if (olap.a_hang > 0 && olap.b_hang > 0)
            if ( olap.flipped )
                return THREE_PRIME;
            else
                return FIVE_PRIME;

        assert(0); // no contained
    }

	// BestOverlapGraph
	// Constructor
	BestOverlapGraph::BestOverlapGraph(int max_fragments)
        : _num_fragments(max_fragments), curFrag(0)
    {
		_best_overlaps = new BestFragmentOverlap[max_fragments+1];

        memset(_best_overlaps, 0, sizeof(BestFragmentOverlap)*(max_fragments+1));
	}

	// Destructor
	BestOverlapGraph::~BestOverlapGraph(){
		delete[] _best_overlaps;
	}

	// Interface to graph visitor
//	void BestOverlapGraph::accept(BestOverlapGraphVisitor bog_vis){
//		bog_vis.visit(this);
//	}

	// Accessor Get Functions
	BestEdgeOverlap *BestOverlapGraph::getBestEdge(
		iuid frag_id, fragment_end_type which_end){

		if(which_end == FIVE_PRIME)
			return(&_best_overlaps[frag_id].five_prime);
		else if(which_end == THREE_PRIME){
			return(&_best_overlaps[frag_id].three_prime);
        }
	}

	BestContainment* BestOverlapGraph::getBestContainer(iuid containee)
    {
        std::map<iuid,BestContainment>::iterator i = _best_containments.find( containee ); 
        if ( i != _best_containments.end() ) 
		    return &_best_containments[containee];
        else
            return NULL;
	}

    void BestOverlapGraph::transitiveContainment() {
        for(std::map<CDS_IID_t,BestContainment>::const_iterator it = _best_containments.begin();
                it != _best_containments.end(); it++)
        {
            CDS_IID_t id = it->first;
            BestContainment bst = it->second;
            bool sameOrient = bst.sameOrientation;
            std::map<CDS_IID_t,BestContainment>::iterator i2 =
                                        _best_containments.find( bst.container);
            std::map<CDS_IID_t,BestContainment> found;
            found[bst.container] = bst;
            while ( i2 != _best_containments.end() ) {
                BestContainment nb = i2->second;
                std::cout << id <<" "<<bst.container<<" "<< nb.container<< std::endl;
                if ( nb.container == id ) {
                    _best_containments.erase( id );
                    std::cout << "Erase self" << std::endl;
                    break;
                }
                std::map<CDS_IID_t,BestContainment>::iterator seen= found.find( nb.container );
                if ( seen != found.end() ) { 
                    _best_containments[ id ] = seen->second;
                    std::cout << "Circled " << seen->second.container<< std::endl;
                    _best_containments.erase( seen->second.container );
                    break;
                }
                _best_containments[ id ] = nb;
                found[ nb.container ] = nb;
                if (!sameOrient)
                    sameOrient = _best_containments[id].sameOrientation = ! nb.sameOrientation;
                found[nb.container].sameOrientation = sameOrient;
                i2 = _best_containments.find( nb.container);
            }
        }
    }

    void BestOverlapGraph::setBestEdge(const Long_Olap_Data_t& olap, float newScore) {

        if (AEnd(olap) == THREE_PRIME) {
            _best_overlaps[ olap.a_iid ].three_prime.frag_b_id = olap.b_iid;
            _best_overlaps[ olap.a_iid ].three_prime.score     = newScore;
            _best_overlaps[ olap.a_iid ].three_prime.bend      = BEnd(olap);

        }
        if (AEnd(olap) == FIVE_PRIME) {
            _best_overlaps[ olap.a_iid ].five_prime.frag_b_id = olap.b_iid;
            _best_overlaps[ olap.a_iid ].five_prime.score     = newScore;
            _best_overlaps[ olap.a_iid ].five_prime.bend      = BEnd(olap);
        }
    }

    uint16         *BestOverlapGraph::fragLength;
    ReadStructp     BestOverlapGraph::fsread = new_ReadStruct();
    FragStoreHandle BestOverlapGraph::fragStoreHandle;
    uint16          BestOverlapGraph::fragLen( iuid iid )
    {
        if (BestOverlapGraph::fragLength[ iid ] == 0) {
            uint32 clrBgn, clrEnd;
            getFragStore( fragStoreHandle, iid, FRAG_S_SEQUENCE, fsread);
            getClearRegion_ReadStruct( fsread, &clrBgn, &clrEnd, READSTRUCT_LATEST);
            BestOverlapGraph::fragLength[ iid ] = clrEnd - clrBgn;
        }
        return BestOverlapGraph::fragLength[ iid ];
    }

    short BestOverlapGraph::olapLength(const Long_Olap_Data_t& olap) {
        uint16 alen = fragLen(olap.a_iid);
        if (olap.a_hang < 0)
            return alen - abs(olap.b_hang);
        else
            return alen - olap.a_hang;
    }

    bool BestOverlapGraph::checkForNextFrag(const Long_Olap_Data_t& olap)
    {
        if (curFrag != olap.a_iid) {
            iuid bid = _best_overlaps[ curFrag ].three_prime.frag_b_id;
            switch(_best_overlaps[ curFrag ].three_prime.bend){
                case THREE_PRIME:
                    _best_overlaps[ bid ].three_prime.in_degree++; break;
                case FIVE_PRIME:
                    _best_overlaps[ bid ].five_prime.in_degree++; break;
                default: assert(0);
            }
            bid = _best_overlaps[ curFrag ].five_prime.frag_b_id;
            switch(_best_overlaps[ curFrag ].five_prime.bend){
                case THREE_PRIME:
                    _best_overlaps[ bid ].three_prime.in_degree++; break;
                case FIVE_PRIME:
                    _best_overlaps[ bid ].five_prime.in_degree++; break;
                default: assert(0);
            }
            _best_overlaps[ olap.a_iid ].five_prime.score = 0;

            curFrag  = olap.a_iid;
            _best_overlaps[ olap.a_iid ].three_prime.score = 0;
            _best_overlaps[ olap.a_iid ].five_prime.score = 0;
            bestLength = 0;
            return true;
        }
        return false;
    }

    void BestOverlapGraph::scoreOverlap(const Long_Olap_Data_t& olap)
    {
        float newScr = score(olap);
        if ( newScr <= 0 )
            return;
/*        if ( olap.a_hang == 0 && olap.b_hang == 0 )
         {
             //multiContain[ olap.a_iid ].equal[ olap.b_iid ] = erate;
             //handle identical containment
         }
         else
*/
         if ( olap.a_hang >= 0 && olap.b_hang <= 0 )
         {
             //handle a contains b
             BestContainment *best = getBestContainer( olap.b_iid );
             if (NULL == best || newScr > best->score || newScr == best->score
                      && fragLen(best->container) < fragLen(olap.a_iid) )
             {
//                 std::cout << olap.a_iid << " contains " << olap.b_iid <<" "<< fragLen(olap.a_iid) << std::endl; 
                 BestContainment newBest;
                 newBest.container = olap.a_iid;
                 newBest.score     = newScr;
                 newBest.sameOrientation = olap.flipped ? false : true;
                 _best_containments[ olap.b_iid ] = newBest;
             }
         }
         else if ( olap.a_hang <= 0 && olap.b_hang >= 0 )
         {
             //handle b contains a
/*             BestContainment *best = getBestContainer( olap.a_iid );
             if (NULL == best || newScr > best->score) {
                 BestContainment newBest;
                 newBest.container = olap.b_iid;
                 newBest.score     = newScr;
                 newBest.sameOrientation = olap.flipped ? false : true;
                 _best_containments[ olap.a_iid ] = newBest;
             }
*/
         } else {
             // no containment, so score
             checkForNextFrag(olap);
             BestEdgeOverlap *best = getBestEdge( olap.a_iid, AEnd(olap));
             short olapLen = olapLength(olap);
             if (newScr > best->score || newScr == best->score &&
                olapLen > bestLength )
             {
                 setBestEdge( olap, newScr );
                 bestLength = olapLen;
             }
         }
    }

    float ErateScore::score(const Long_Olap_Data_t& olap) {

        return 100 - Expand_Quality(olap.corr_erate) * 100;
    }
    
    float LongestEdge::score(const Long_Olap_Data_t& olap) {

        return olapLength(olap);
    }

    float LongestHighIdent::score(const Long_Olap_Data_t& olap) {

        short olapLen = olapLength(olap);
        float erate = Expand_Quality(olap.corr_erate) * 100;
        if (erate > mismatchCutoff)
            return 0;
        return olapLen;
    }
} //AS_BOG namespace

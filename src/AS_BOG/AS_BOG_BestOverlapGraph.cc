
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
 * $Id: AS_BOG_BestOverlapGraph.cc,v 1.5 2005-08-10 14:46:19 eliv Exp $
 * $Revision: 1.5 $
*/

static const char CM_ID[] = "$Id: AS_BOG_BestOverlapGraph.cc,v 1.5 2005-08-10 14:46:19 eliv Exp $";

//  System include files

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
	BestOverlapGraph::BestOverlapGraph(int max_fragments) : _num_fragments(max_fragments){
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

//	BestContainment *BestOverlapGraph::getBestContainer(iuid containee){
//		return(&_best_containments[containee]);
//	}

    void BestOverlapGraph::setBestEdge(const Long_Olap_Data_t& olap, float newScore) {

        if (AEnd(olap) == THREE_PRIME) {
            _best_overlaps[ olap.a_iid ].three_prime.frag_b_id = olap.b_iid;
            _best_overlaps[ olap.a_iid ].three_prime.score       = newScore;

        }
        if (AEnd(olap) == FIVE_PRIME) {
            _best_overlaps[ olap.a_iid ].five_prime.frag_b_id = olap.b_iid;
            _best_overlaps[ olap.a_iid ].five_prime.score       = newScore;
        }
        if (BEnd(olap) == THREE_PRIME)
            _best_overlaps[ olap.b_iid ].three_prime.in_degree++;

        if (BEnd(olap) == FIVE_PRIME)
            _best_overlaps[ olap.b_iid ].five_prime.in_degree++;
    }

    uint16 *BestOverlapGraph::fragLength;
    ReadStructp BestOverlapGraph::fsread;
    FragStoreHandle BestOverlapGraph::fragStoreHandle;
    uint16 BestOverlapGraph::fragLen( iuid iid ) {
        if (BestOverlapGraph::fragLength[ iid ] == 0) {
            uint32 clrBgn, clrEnd;
            getFragStore( fragStoreHandle, iid, FRAG_S_SEQUENCE, fsread);
            getClearRegion_ReadStruct( fsread, &clrBgn, &clrEnd, READSTRUCT_LATEST);
            BestOverlapGraph::fragLength[ iid ] = clrEnd - clrBgn;
        }
        return BestOverlapGraph::fragLength[ iid ];
    }

    bool BestOverlapGraph::checkForNextFrag(const Long_Olap_Data_t& olap, float scoreReset) {
        if (curFrag != olap.a_iid) {
            curFrag  = olap.a_iid;
            _best_overlaps[ olap.a_iid ].three_prime.score = scoreReset;
            _best_overlaps[ olap.a_iid ].five_prime.score = scoreReset;
            return true;
        }
        return false;
    }

    bool ErateScore::checkForNextFrag(const Long_Olap_Data_t& olap, float scoreReset) {
        if (BestOverlapGraph::checkForNextFrag( olap, scoreReset)) 
            bestLength = 0;
    }

    float ErateScore::score(const Long_Olap_Data_t& olap) {

        float erate = Expand_Quality(olap.corr_erate) * 100;
        checkForNextFrag(olap,100);
        BestEdgeOverlap *best = getBestEdge( olap.a_iid, AEnd(olap));
        short olapLen = olapLength(olap);
        if (erate < best->score || erate == best->score && olapLen > bestLength ) {
            setBestEdge( olap, erate );
            bestLength = olapLen;
        }
    }
    
    float LongestEdge::score(const Long_Olap_Data_t& olap) {

        uint16 alen = BestOverlapGraph::fragLen(olap.a_iid);
        short olapLen = olapLength(olap);
        checkForNextFrag(olap,0);
        BestEdgeOverlap *best = getBestEdge( olap.a_iid, AEnd(olap));
        if (best->score < olapLen) 
            setBestEdge( olap, olapLen );
    }


    float LongestHighIdent::score(const Long_Olap_Data_t& olap) {

        uint16 alen = BestOverlapGraph::fragLen(olap.a_iid);
        short olapLen = olapLength(olap);
        float erate = Expand_Quality(olap.corr_erate) * 100;
        checkForNextFrag(olap,0);

        if (erate > mismatchCutoff)
            return 0;

        BestEdgeOverlap* best = getBestEdge( olap.a_iid, AEnd(olap));
        if (best->score < olapLen) {
            setBestEdge( olap, olapLen );
        }
    }
} //AS_BOG namespace

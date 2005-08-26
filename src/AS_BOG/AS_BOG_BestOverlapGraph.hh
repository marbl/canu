
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
* Module:  AS_BOG_BestOverlapGraph.hh
* Description:
*	Data structure to contain the best overlaps and containments
*	based on a defined metric.
* 
*    Programmer:  K. Li
*       Started:  20 July 2005
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: AS_BOG_BestOverlapGraph.hh,v 1.11 2005-08-26 20:38:02 eliv Exp $
 * $Revision: 1.11 $
*/

//  System include files

#ifndef INCLUDE_AS_BOG_BESTOVERLAPGRAPH
#define INCLUDE_AS_BOG_BESTOVERLAPGRAPH

#include <map>
#include "AS_BOG_Datatypes.hh"

extern "C" {
#include "OlapStoreOVL.h"
}

namespace AS_BOG{

	///////////////////////////////////////////////////////////////////////

	struct BestEdgeOverlap{
		iuid frag_b_id;
		fragment_end_type bend;		
		int in_degree;
		float score;
	};

	struct BestFragmentOverlap{

		// Contains information on what a known fragment overlaps.
		// It is assumed that an index into an array of BestOverlap
		// will tell us what fragment has this best overlap

		BestEdgeOverlap five_prime; 
		BestEdgeOverlap three_prime;
	};

	///////////////////////////////////////////////////////////////////////

	struct BestContainment{

		// Contains what kind of containment relationship exists between
		// fragment a and fragment b

		iuid container;
		float score;
		bool sameOrientation;
	};

	///////////////////////////////////////////////////////////////////////


	struct BestOverlapGraph {

            // Class methods
        static uint16 *fragLength;
        static uint16 fragLen( iuid );
        static ReadStructp fsread;
        static FragStoreHandle fragStoreHandle;
        static overlap_type getType(const Long_Olap_Data_t & olap);
        fragment_end_type AEnd(const Long_Olap_Data_t & olap);
        fragment_end_type BEnd(const Long_Olap_Data_t& olap);
        short        olapLength(const Long_Olap_Data_t& olap) {
            uint16 alen = fragLen(olap.a_iid);
            if (olap.a_hang < 0)
                return alen - abs(olap.b_hang);
            else
                return alen - olap.a_hang;
        }

        // Constructor, parametrizing maximum number of overlaps
        BestOverlapGraph(int max_fragment_count);

        // Destructor
        ~BestOverlapGraph();

        // Interface to graph visitor
        //			accept(BestOverlapGraphVisitor bog_vis){
        //				bog_vis.visit(this);
//			}

			// Accessor Get Functions
			BestEdgeOverlap *getBestEdge(
				iuid frag_id, fragment_end_type which_end);

			void setBestEdge(const Long_Olap_Data_t& olap, float newScore);
			iuid getNumFragments() { return _num_fragments; }

			BestContainment *getBestContainer(iuid frag_id);

            void scoreOverlap(const Long_Olap_Data_t& olap);
            bool checkForNextFrag(const Long_Olap_Data_t& olap);
            virtual float score( const Long_Olap_Data_t& olap) =0;

		protected:
			iuid _num_fragments;
            iuid curFrag;
            int bestLength;

			BestFragmentOverlap* _best_overlaps;
            std::map<iuid, BestContainment> _best_containments;

	}; //BestOverlapGraph

    struct ErateScore : public BestOverlapGraph {
        ErateScore(int num) : BestOverlapGraph(num) {}
        float score( const Long_Olap_Data_t& olap);
    };

    struct LongestEdge : public BestOverlapGraph {
        LongestEdge(int num) : BestOverlapGraph(num) {}
        float score( const Long_Olap_Data_t& olap);
    };

    struct LongestHighIdent : public BestOverlapGraph {
        float mismatchCutoff;
        LongestHighIdent(int num, float maxMismatch)
            : BestOverlapGraph(num), mismatchCutoff(maxMismatch) {}
        float score( const Long_Olap_Data_t& olap);
    };

} //AS_BOG namespace


#endif //INCLUDE_AS_BOG_BESTOVERLAPGRAPH


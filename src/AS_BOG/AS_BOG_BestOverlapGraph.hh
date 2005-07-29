
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
* Module:  AS_BOG_BestOverlapGraph.c
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
 * $Id: AS_BOG_BestOverlapGraph.hh,v 1.2 2005-07-29 21:12:14 kli1000 Exp $
 * $Revision: 1.2 $
*/

static char CM_ID[] = "$Id: AS_BOG_BestOverlapGraph.hh,v 1.2 2005-07-29 21:12:14 kli1000 Exp $";

//  System include files

#ifndef INCLUDE_AS_BOG_BESTOVERLAPGRAPH
#define INCLUDE_AS_BOG_BESTOVERLAPGRAPH

#include "AS_BOG_Datatypes.hh"

namespace AS_BOG{

	///////////////////////////////////////////////////////////////////////

	class BestOverlap{

		// Contains information on what a known fragment overlaps.
		// It is assumed that an index into an array of BestOverlap
		// will tell us what fragment has this best overlap

		public:
			iuid ovl_frag_id;
			overlap_type type;		
			int three_prime_in_degree;
			int five_prime_in_degree;
			int score;
	}

	///////////////////////////////////////////////////////////////////////

	class BestContainment{

		// Contains what kind of containment relationship exists between
		// fragment a and fragment b

		public:
			iuid frag_a_id;
			iuid frag_b_id;
			overlap_type type;
			int score;
	}

	///////////////////////////////////////////////////////////////////////

	class BestOverlapGraph {

		public:

			// Constructor, parametrizing maximum number of overlaps
			BestOverlapGraph(int max_best_overlap_size){
				bestOverlaps=new BestOverlap(max_best_overlap_size);
			}

			// Destructor
			~BestOverlapGraph(void){
				delete bestOverlaps[];
				bestContainments.clear();
			}

			// Interface to graph visitor
			accept(BestOverlapGraphVisitor bog_vis){
				bog_vis.visit(this);
			}

			// Accessor Get Functions
			BestOverlap *getBestOverlap(iuid frag_id){
				return(&(bestOverlaps[frag_id]));
			}
	
			BestContainment *getBestContainment(iuid frag_id);

			// Accessor Set Functions
			void setBestOverlap(
				iuid frag_a_id, 
				iuid frag_b_id,
				int score,
				overlap_type ovl_type
			);

		private:
			BestOverlap _bestOverlaps[];
			int _numBestOverlaps;
			vector<BestContainment> _bestContainments;

	} //BestOverlapGraph

} //AS_BOG namespace


#endif //INCLUDE_AS_BOG_BESTOVERLAPGRAPH


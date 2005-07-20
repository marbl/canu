
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
 * $Id: AS_BOG_BestOverlapGraph.hh,v 1.1 2005-07-20 21:00:25 kli1000 Exp $
 * $Revision: 1.1 $
*/

static char CM_ID[] = "$Id: AS_BOG_BestOverlapGraph.hh,v 1.1 2005-07-20 21:00:25 kli1000 Exp $";

//  System include files

#ifndef INCLUDE_AS_BOG_BESTOVERLAPGRAPH
#define INCLUDE_AS_BOG_BESTOVERLAPGRAPH

namespace AS_BOG{

	typedef enum { 
		DOVE_NORMAL,		// AB_AB
		DOVE_ANTI_NORMAL,	// BA_BA
		DOVE_INNIE,		// AB_BA
		DOVE_OUTTIE,		// BA_AB
		CONT_A_CONTAINS,
		CONT_B_CONTAINS,
		CONT_MUTUAL
	} overlap_type;

	typedef unsigned int iuid;

	///////////////////////////////////////////////////////////////////////

	class BestOverlap{
		public:
			iuid ovl_frag_id;
			overlap_type type;		
			int three_prime_in_degree;
			int five_prime_in_degree;
			int score;
	}

	///////////////////////////////////////////////////////////////////////

	class BestContainment{
		public:
			iuid frag_a_id;
			iuid frag_b_id;
			overlap_type type;
			int score;
	}

	///////////////////////////////////////////////////////////////////////

	class BestOverlapGraph {

		public:
			BestOverlapGraph(int max_best_overlap_size){
				bestOverlaps=new BestOverlap(max_best_overlap_size);
			}

			~BestOverlapGraph(void){
				delete bestOverlaps[];
				bestContainments.clear();
			}

			accept(BestOverlapGraphVisitor bog_vis){
				bog_vis.visit(this);
			}

			BestOverlap *getBestOverlap(iuid frag_id){
				return(&(bestOverlaps[frag_id]));
			}
	
			BestContainment *getBestContainment(iuid frag_id);

			void setBestOverlap(
				iuid frag_a_id, 
				iuid frag_b_id,
				int score,
				overlap_type ovl_type
			);

		private:
			BestOverlap bestOverlaps[];
			int numBestOverlaps;
			vector<BestContainment> bestContainments;

	} //BestOverlapGraph

} //AS_BOG namespace


#endif //INCLUDE_AS_BOG_BESTOVERLAPGRAPH


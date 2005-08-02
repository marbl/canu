
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
 * $Id: AS_BOG_BestOverlapGraph.cc,v 1.1 2005-08-02 15:49:12 kli1000 Exp $
 * $Revision: 1.1 $
*/

static char CM_ID[] = "$Id: AS_BOG_BestOverlapGraph.cc,v 1.1 2005-08-02 15:49:12 kli1000 Exp $";

//  System include files

#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_BestOverlapGraphVisitor.hh"

namespace AS_BOG{


	// BestOverlap::BestOverlap

	// BestContainment::BestContainment

	// BestOverlapGraph
	BestOverlapGraph::BestOverlapGraph(int max_fragments){
		_bestOverlaps = new BestOverlap[max_fragments];

		for(iuid i=0; i<max_fragments; i++){
			_bestOverlaps[i].ovl_frag_id=0;
			_bestOverlaps[i].type=AS_BOG::UNDEFINED;
			_bestOverlaps[i].three_prime_in_degree=0;
			_bestOverlaps[i].five_prime_in_degree=0;
			_bestOverlaps[i].score=-1;
		}

		_numBestOverlaps = 0;
	}

	// Destructor
	BestOverlapGraph::~BestOverlapGraph(void){
		delete[] BestOverlap;
	}

	// Interface to graph visitor
	void BestOverlapGraph::accept(BestOverlapGraphVisitor bog_vis){
		bog_vis.visit(this);
	}

	// Accessor Get Functions
	BestOverlap *BestOverlapGraph::getBestOverlap(iuid frag_id){
		return(&(_bestOverlaps[frag_id]));
	}

	BestContainment *BestOverlapGraph::getBestContainer(iuid containee){
		iuid i;
		iuid container;
		for(i=_bestContainments.begin; i<_bestContainments.end; i++){
			switch(_bestContainments[i].overlap_type){

				case AS_BOG::CONT_A_CONTAINS:
					if(containee==_bestContainments[i].frag_b_id){
						return(&_bestContainments[i].frag_a_id);
					}
					break;
				case AS_BOG::CONT_B_CONTAINS:
					if(containee==_bestContainments[i].frag_a_id){
						return(&_bestContainments[i].frag_b_id);
					}
					break;
				case AS_BOG::CONT_MUTUAL:
					if(containee==_bestContainments[i].frag_a_id){
						return(&_bestContainments[i].frag_b_id);
					}else if(containee==_bestContainments[i].frag_b_id){
						return(&_bestContainments[i].frag_a_id);
					}
			}
		}
		return(NULL);
	}

	// Accessor Set Functions
	void BestOverlapGraph::setBestOverlap(
		iuid frag_a_id, 
		iuid frag_b_id,
		int score,
		overlap_type ovl_type
	){
		switch(ovl_type){
			case 
			AS_BOG::DOVE_NORMAL,
			AS_BOG::DOVE_INNIE,
			AS_BOG::DOVE_OUTTIE,
			AS_BOG::DOVE_ANTI_NORMAL:
				_best_overlaps[frag_a_id].frag_b_id=frag_b_id;
				_best_overlaps[frag_a_id].type=ovl_type;
				_best_overlaps[frag_a_id].score=score;
				switch(ovl_type){
					case AS_BOG::DOVE_NORMAL:
						_best_overlaps[frag_b_id].five_prime_in_degree++;
						break;
					case AS_BOG::DOVE_INNIE:
						_best_overlaps[frag_b_id].three_prime_in_degree++;
						break;
					case AS_BOG::DOVE_OUTTIE:
						_best_overlaps[frag_b_id].five_prime_in_degree++;
						break;
					case AS_BOG::DOVE_ANTI_NORMAL:
						_best_overlaps[frag_b_id].three_prime_in_degree++;
						break;
					default:	
				}

			case 
			AS_BOG::CONT_NORMAL,
			AS_BOG::CONT_INNIE,
			AS_BOG::CONT_OUTTIE:
			AS_BOG::CONT_ANTI_NORMAL:
			AS_BOG::CONT_MUTUAL:
			AS_BOG::CONT_MUTUAL_INNIE:
				BestContainment best_containment;

				best_containment.frag_a_id=frag_a_id;
				best_containment.frag_b_id=frag_b_id;
				best_containment.type=ovl_type;
				best_containment.score=score;

				_best_containments.push_back(best_containment);
			break;

			default:
		}
		
	}

} //AS_BOG namespace


#endif //INCLUDE_AS_BOG_BESTOVERLAPGRAPH



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
 * $Id: AS_BOG_BestOverlapGraph.cc,v 1.3 2005-08-02 21:29:25 kli1000 Exp $
 * $Revision: 1.3 $
*/

static char CM_ID[] = "$Id: AS_BOG_BestOverlapGraph.cc,v 1.3 2005-08-02 21:29:25 kli1000 Exp $";

//  System include files

#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_BestOverlapGraphVisitor.hh"

namespace AS_BOG{

	// BestOverlapGraph
	// Constructor
	BestOverlapGraph::BestOverlapGraph(int max_fragments){
		_best_overlaps = new BestFragmentOverlap[max_fragments];

		for(iuid i=0; i<max_fragments; i++){
			_best_overlaps[i].five_prime.ovl_frag_id=0;
			_best_overlaps[i].five_prime.type=UNDEFINED;
			_best_overlaps[i].five_prime.in_degree=0;
			_best_overlaps[i].five_prime.score=-1;

			_best_overlaps[i].three_prime.ovl_frag_id=0;
			_best_overlaps[i].three_prime.type=UNDEFINED;
			_best_overlaps[i].three_prime.in_degree=0;
			_best_overlaps[i].three_prime.score=-1;
		}

		_numBestOverlaps = 0;
	}

	// Destructor
	BestOverlapGraph::~BestOverlapGraph(void){
		delete[] _best_overlaps;
	}

	// Interface to graph visitor
	void BestOverlapGraph::accept(BestOverlapGraphVisitor bog_vis){
		bog_vis.visit(this);
	}

	// Accessor Get Functions
	BestOverlap *BestOverlapGraph::getBestOverlap(iuid frag_id){
		return(&_best_overlaps[frag_id]);
	}

	BestContainment *BestOverlapGraph::getBestContainer(iuid containee){
		return(&_best_containments[containee]);
	}

	// Accessor Set Functions
	void BestOverlapGraph::setBestOverlap(
		iuid frag_a_id, 
		iuid frag_b_id,
		int score,
		overlap_type ovl_type
	){
		switch(ovl_type){
		case DOVE_NORMAL: 
		case DOVE_INNIE: 
		case DOVE_OUTTIE: 
		case DOVE_ANTI_NORMAL: 
		case DOVE_NORMAL_BACK:

			_best_overlaps[frag_a_id].frag_b_id=frag_b_id;
			_best_overlaps[frag_a_id].type=ovl_type;
			_best_overlaps[frag_a_id].score=score;

			switch(ovl_type){

			case DOVE_NORMAL:
			case DOVE_OUTTIE:
				_best_overlaps[frag_b_id].five_prime_in_degree++;
				break;
			case DOVE_INNIE:
			case DOVE_ANTI_NORMAL:
				_best_overlaps[frag_b_id].three_prime_in_degree++;
				break;
			default:	
			}
		break;

		case CONT_NORMAL: 
		case CONT_INNIE:
		case CONT_OUTTIE: 
		case CONT_ANTI_NORMAL: 

			BestContainment best_containment;
			best_containment.container=frag_a_id;
			best_containment.type=ovl_type;
			best_containment.score=score;

			_best_containments[frag_b_id];
		break;

		default:
		}
		
	}

} //AS_BOG namespace


#endif //INCLUDE_AS_BOG_BESTOVERLAPGRAPH


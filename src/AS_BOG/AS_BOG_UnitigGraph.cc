
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
* Module:  AS_BOG_UnitigGraph.cc
* Description:
*	Data structure to contain the unitig paths and how they connect to each
*	other
* 
*    Programmer:  K. Li
*       Started:  2 Aug 2005
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: AS_BOG_UnitigGraph.cc,v 1.2 2005-08-08 21:49:02 kli1000 Exp $
 * $Revision: 1.2 $
*/

static char CM_ID[] = "$Id: AS_BOG_UnitigGraph.cc,v 1.2 2005-08-08 21:49:02 kli1000 Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"

namespace AS_BOG{

	///////////////////////////////////////////////////////////////////////
	// DoveTailPath

	DoveTailPath::SetNextFragment(iuid frag_id, orientation_type ori){
		_dt_path_node dtp_node;
		dtp_node.frag_id=frag_id;
		dtp_node.ori=ori;

		_dt_path.push_back(dtp_node);
	}

	ostream& operator<<(ostream& os, DoveTailPath &dtp){
		vector<_dt_path_node>::iterator i=_dt_path.begin();	
		vector<_dt_path_node>::iterator last=_dt_path.end();	
		
		while(i!=last){
			os << i->frag_id << " [" << i->ori << "]\n";
			i++;
		}

		return(os);
	}

	


	///////////////////////////////////////////////////////////////////////

	// Unitig::
	// Constructor
	Unitig::Unitig(iuid utg_id){
		unitig_id=utg_id;
	}

	Unitig::set_next_dovetailing_fragment(iuid frag_id){
		dovetail_path.push_back(frag_id);
	}

	Unitig::set_contained_fragment(iuid container, iuid containee){
		contained_frags[container]=containee;
	}

	// Computable statistics on this unitig
	float Unitig::Compute_Astat(void){
		return(0);
	}
	float Unitig::Compute_Rho(void){
		return(0);
	}

	///////////////////////////////////////////////////////////////////////

	// UnitigOverlap::
	UnitigOverlap::UnitigOverlap(iuid utg_a, iuid utg_b, overlap_type type, int score){
		unitig_a=utg_a;
		unitig_b=utg_b;
		ovl_type=type;
		ovl_score=score;
	}

	///////////////////////////////////////////////////////////////////////

	// UnitigGraph::

	void UnitigGraph::include_unitig(Unitig &utg){
		unitigs.push_back(utg);
	}
	void UnitigGraph::include_overlap(UnitigOverlap &utg_ovl){
		unitig_overlaps(utg_ovl);
	}

}



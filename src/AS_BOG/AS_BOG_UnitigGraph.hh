
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
* Module:  AS_BOG_UnitigGraph.hh
* Description:
*	Data structure to contain the unitig paths and how they connect to each
*	other
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
 * $Id: AS_BOG_UnitigGraph.hh,v 1.3 2005-08-02 21:55:22 kli1000 Exp $
 * $Revision: 1.3 $
*/

static char CM_ID[] = "$Id: AS_BOG_UnitigGraph.hh,v 1.3 2005-08-02 21:55:22 kli1000 Exp $";

#ifndef INCLUDE_AS_BOG_UNITIGGRAPH
#define INCLUDE_AS_BOG_UNITIGGRAPH

#include <vector>
#include <map>
#include "AS_BOG_Datatypes.hh"

namespace AS_BOG{

	///////////////////////////////////////////////////////////////////////

	class Unitig{

		// Contains the "chunk", based on a list of fragment iuids
		// which is the dovetail path.  The contained fragments are
		// also stored here.

		public:
			// Constructor
			Unitig(iuid utg_id);
			set_next_dovetailing_fragment(iuid frag_id, orientation_type ori);
			set_contained_fragment(iuid container, iuid containee);

			// Computable statistics on this unitig
			float Compute_Astat(void);
			float Compute_Rho(void);

			iuid unitig_id;
			vector<iuid> dovetail_path;
			map<iuid, iuid> contained_frags;
	}

	///////////////////////////////////////////////////////////////////////

	class UnitigOverlap{

		// Contains the description of how unitigs are connected to
		// each other
			
		public:
			UnitigOverlap(iuid utg_a, iuid utg_b, overlap_type type, int score);

			iuid unitig_a;
			iuid unitig_b;
			overlap_type ovl_type;
			int ovl_score;
	}

	///////////////////////////////////////////////////////////////////////

	class UnitigGraph{

		// Contains the set of unitigs and how they are connected together
		// in one object for a given assembly.

		public:
			void include_unitig(Unitig &utg);
			void include_overlap(UnitigOverlap &utg_ovl);

			vector<&Unitig> unitigs;
			vector<&UnitigOverlap> unitig_overlaps;

	}

}

#endif


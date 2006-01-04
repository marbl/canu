
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
 * $Id: AS_BOG_UnitigGraph.hh,v 1.6 2006-01-04 15:17:39 eliv Exp $
 * $Revision: 1.6 $
*/


#ifndef INCLUDE_AS_BOG_UNITIGGRAPH
#define INCLUDE_AS_BOG_UNITIGGRAPH

static char AS_BOG_UNITIG_GRAPH_HH_CM_ID[] = "$Id: AS_BOG_UnitigGraph.hh,v 1.6 2006-01-04 15:17:39 eliv Exp $";

#include <vector>
#include <map>
#include <iostream>
#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_ChunkGraph.hh"

#include "AS_MSG_pmesg.h"

namespace AS_BOG{

	///////////////////////////////////////////////////////////////////////
	// For the sake of self-documenting code

	typedef iuid container_id;
	typedef iuid containee_id;
	typedef iuid fragment_id;
	typedef std::vector<fragment_id> FragmentList;

	///////////////////////////////////////////////////////////////////////

	struct DoveTailNode{
		iuid frag_id;
		orientation_type ori;
		long olap_len; //3' overlap, last node has no overlap
	};
	typedef std::vector<DoveTailNode> DoveTailPath;


	struct ContaineeNode{
		iuid frag_id;
		bool same_ori;
		long olap_offset;
	};
	typedef std::vector<ContaineeNode> ContaineeList;	
	typedef std::map<container_id, ContaineeList> ContainerMap;

	///////////////////////////////////////////////////////////////////////

	// Position of fragment in unitig
	struct interval{ 
		long begin;
		long end;
	};
	typedef std::map<fragment_id, interval> FragmentPositionMap;
	std::ostream& operator<< (std::ostream& os, FragmentPositionMap *fpm_ptr);

	///////////////////////////////////////////////////////////////////////

	struct Unitig{

		Unitig(void);		
		~Unitig(void);		

		// Compute unitig based on given dovetails and containments
		FragmentPositionMap *computeFragmentPositions(void);

		// Accessor methods
		float getAvgRho(void);
		void setGlobalArrivalRate(float global_arrival_rate);
		float getCovStat(void);
		long getLength(void);
		long getNumFrags(void);
		long getNumRandomFrags(void); // For now, same as numFrags, but should be randomly sampled frag count
		// For proto I/O messages
		IntUnitigMesg *getIUM_Mesg(void);
		void freeIUM_Mesg(IntUnitigMesg *ium_ptr);

		friend std::ostream& operator<< (std::ostream& os, Unitig& utg);

		// Public Member Variables
		iuid id;
	        DoveTailPath *dovetail_path_ptr;
		ContainerMap *contained_frags_ptr;
		FragmentPositionMap *frag_pos_map_ptr;

		private:
			// Do not access these private variables directly, they may not be
			//  computed yet, use accessors!
			float _avgRho;
			float _covStat;
			long _length;
			long _numFrags;
			long _numRandomFrags;
			static float _globalArrivalRate;
	};

	///////////////////////////////////////////////////////////////////////

	struct UnitigOverlap{
		// Contains the description of how unitigs are connected to
		// each other
			
		iuid unitig_a;
		orientation_type ori_a;

		iuid unitig_b;
		orientation_type ori_b;

		overlap_type ovl_type;
		int ovl_score;
	};

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	typedef std::vector<Unitig*> UnitigVector;

	struct UnitigGraph{
		// This will store the entire set of unitigs that are generated
		// It's just a unitig container.

		// Call this on a chunk graph pointer to build a unitig graph
		void build(ChunkGraph *cg_ptr, BestOverlapGraph *bog_ptr, long num_rand_frags, long genome_size);

		// Debugging output operator
		friend std::ostream& operator<< (std::ostream& os, UnitigGraph& utgrph);

		// For compatibility with the rest of the assembler
		void writeIUMtoFile(char *filename);

		float getGlobalArrivalRate(long total_random_frags_in_genome=0, long genome_size=0);

		///////////////////////////////////////////////////////////////
		// Member Variables

		// Unitigs are the dove tails and their contained fragments
		UnitigVector unitigs;

		// Overlaps are unitig overlaps
		std::vector<UnitigOverlap*> overlaps;


		private:
			// Given a fragment, it will follow it's overlaps until 
			//   the end, and return the path that it traversed.
			DoveTailPath *_extract_dovetail_path(
				iuid src_frag_id, 
				orientation_type direction,
				ChunkGraph *cg_ptr,
				BestOverlapGraph *bog_ptr);

			// Inverts the containment map to key by container, instead of containee
			ContainerMap *_build_container_map(BestOverlapGraph *bog_ptr);

			// Build containee list
			ContainerMap *_extract_containees(DoveTailPath *dtp_ptr, 
				ContainerMap *cntnrmap_ptr);

			// Given a fragment and end, it will find the fragment it
			// overlaps and it's end.
			void _follow_to_next_fragment(
				iuid src_frag_id, orientation_type src_end,
				iuid& dst_frag_id, orientation_type& dst_end,
				ChunkGraph *cg_ptr);

			// Compute the global arrival rate based on the unitig rho's.
			float _compute_global_arrival_rate(void);

	};
		
	///////////////////////////////////////////////////////////////////////

};

#endif


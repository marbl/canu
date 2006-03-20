
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
* Module:  AS_BOG_ChunkGraph.hh
* Description:
*	Data structure to keep track of if overlaps can be chunked together.
*	We should be able to get rid of the data structure, but the rules
*	for deciding whether a overlaps can be chunked together should be
*	kept here.
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
 * $Id: AS_BOG_ChunkGraph.hh,v 1.3 2006-03-20 18:51:19 eliv Exp $
 * $Revision: 1.3 $
*/

#ifndef INCLUDE_AS_BOG_CHUNKGRAPH
#define INCLUDE_AS_BOG_CHUNKGRAPH

static char AS_BOG_CHUNK_GRAPH_HH_CM_ID[] = "$Id: AS_BOG_ChunkGraph.hh,v 1.3 2006-03-20 18:51:19 eliv Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_BestOverlapGraph.hh"

namespace AS_BOG{

	struct ChunkGraph{

	    public:

		// Constructor
		ChunkGraph(void);

		// Destructor
		~ChunkGraph(void);

		// Build the ChunkGraph, based on a BOG
		void build(BestOverlapGraph *bovlg);		

		// Chunkability rule
		bool isChunkable(iuid frag_id, fragment_end_type whichEnd);

	        bool isChunkable(
			BestEdgeOverlap *beo,
			BestOverlapGraph *bovlg);

	        bool isChunkable(
			iuid frag_a_id, fragment_end_type which_end,
			BestOverlapGraph *bovlg);

		// Returns IUID of 5' or 3' end of specified frag_id
		//  Since there should only be one out/incoming connection
		iuid getChunking(iuid src_frag_id, fragment_end_type whichEnd);

		void getChunking(
			iuid src_frag_id, 
			iuid* five_prime_dst_frag_id, iuid* three_prime_dst_frag_id);

		void setChunking(
			iuid src_frag_id, 
			iuid five_prime_dst_frag_id, iuid three_prime_dst_frag_id);

		void printFrom(iuid begin, iuid end);
		long getNumFragments(void);
		long countSingletons(void);

		friend std::ostream& operator<< (std::ostream& os, ChunkGraph &cg);

		void checkInDegree(BestOverlapGraph *bovlg);

	    private:

		struct _chunk_unit_struct{
			iuid five_prime;
			iuid three_prime;
		};

		_chunk_unit_struct *_chunkable_array;
		iuid _max_fragments;

	};

} //AS_BOG namespace


#endif


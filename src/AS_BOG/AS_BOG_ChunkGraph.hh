
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
 * $Id: AS_BOG_ChunkGraph.hh,v 1.1 2005-08-08 21:51:03 kli1000 Exp $
 * $Revision: 1.1 $
*/

static char CM_ID[] = "$Id: AS_BOG_ChunkGraph.hh,v 1.1 2005-08-08 21:51:03 kli1000 Exp $";

//  System include files

#ifndef INCLUDE_AS_BOG_CHUNKGRAPH
#define INCLUDE_AS_BOG_CHUNKGRAPH

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_BestOverlapGraphVisitor.hh"

namespace AS_BOG{

	class ChunkGraph{

	    public:

		ChunkGraph(iuid max_fragments);
		~ChunkGraph(void);

		void accept(ChunkGraphVisitor cgv);
	
		iuid getFivePrimeChunking(iuid src_frag_id);
		iuid getThreePrimeChunking(iuid src_frag_id);

		void getChunking(
			iuid src_frag_id, 
			iuid& five_prime_dst_frag_id, iuid& three_prime_dst_frag_id);

		void setChunking(
			iuid src_frag_id, 
			iuid five_prime_dst_frag_id, iuid three_prime_dst_frag_id);

	    private:

		struct _chunk_unit_struct{
			iuid five_prime;
			iuid three_prime;
		}

		_chunk_unit_struct _chunkable_array[];
		iuid _num_fragments;
		iuid _max_fragments;

	}

	//////////////////////////////////////////////////////////////////////////////

	class ChunkGraphBuilder : public BestOverlapGraphVisitor{

	    public:

		void visit(BestOverlapGraph& bovlg);

		bool isChunkable(
			iuid frag_id, fragment_end_type which_end, BestOverlapGraph& bovlg);

		bool isChunkable(
			BestEdgeOverlap *beo, BestOverlapGraph& bovlg);

	    private:

		ChunkGraph *_chunk_graph_ptr=NULL;

	}

	//////////////////////////////////////////////////////////////////////////////

	class ChunkGraphVisitor{
	    public:
		virtual void visit(ChunkGraph &cg);
	}
		

} //AS_BOG namespace


#endif


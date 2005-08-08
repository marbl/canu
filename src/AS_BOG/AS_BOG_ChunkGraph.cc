
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
* Module:  AS_BOG_ChunkGraph.cc
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
 * $Id: AS_BOG_ChunkGraph.cc,v 1.1 2005-08-08 21:51:03 kli1000 Exp $
 * $Revision: 1.1 $
*/

static char CM_ID[] = "$Id: AS_BOG_ChunkGraph.cc,v 1.1 2005-08-08 21:51:03 kli1000 Exp $";

//  System include files

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_BestOverlapGraphVisitor.hh"
#include "AS_BOG_ChunkGraph.hh"

namespace AS_BOG{

	//////////////////////////////////////////////////////////////////////////////
	// ChunkGraph

	ChunkGraph::ChunkGraph(iuid max_fragments){
		_chunkable_array=new _chunk_unit_struct[max_fragments];
		_max_fragments=max_fragments;
	}
	
	ChunkGraph::~ChunkGraph(void){
		delete[] _chunkable_array;
	}

	void ChunkGraph::accept(ChunkGraphVisitor cgv){
		cgv.visit(this);
	}


	//////////////////////////////////////////////////////////////////////////////
	// Accessor Functions

	iuid ChunkGraph::getFivePrimeChunking(iuid src_frag_id){
		return(_chunkable_array[frag_id].five_prime);
	}
	iuid ChunkGraph::getThreePrimeChunking(iuid src_frag_id){
		return(_chunkable_array[frag_id].three_prime);
	}

	void ChunkGraph::getChunking(
		iuid src_frag_id,
		iuid& five_prime_dst_frag_id, iuid& three_prime_dst_frag_id){

		five_prime_dst_frag_id=_chunkable_array[src_frag_id].five_prime;
		three_prime_dst_frag_id=_chunkable_array[src_frag_id].three_prime;
		
	}

	void ChunkGraph::setChunking(
		iuid src_frag_id,
		iuid five_prime_dst_frag_id, iuid three_prime_dst_frag_id){

		_chunkable_array[src_frag_id].five_prime=five_prime_dst_frag_id;
		_chunkable_array[src_frag_id].three_prime=three_prime_dst_frag_id;
		
	}

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	// ChunkGraphBuilder

	void ChunkGraphBuilder::visit(BestOverlapGraph& bovlg){
		// This will go through all the nodes in the bovlg and
		//   produce the ChunkGraph

		iuid i;
		iuid num_frags=bovlg.getNumFragments();

		if(_chunk_graph_ptr!=NULL){
			delete _chunk_graph_ptr;
		}
		_chunk_graph_ptr=new _chunk_graph(num_frags);

		for(frag_id=0; frag_id<num_frag; frag_id++){

			BestEdgeOverlap *fp_beo=
				bovlg.getBestEdgeOverlap(frag_id, FIVE_PRIME);
			bool fp_chunkability=isChunkable(fp_beo, bovlg);
			
			BestEdgeOverlap *tp_beo=
				bovlg.getBestEdgeOverlap(frag_id, THREE_PRIME);
			bool tp_chunkability=isChunkable(tp_beo, bovlg);

			_chunk_graph_ptr->setChunking(frag_id, 
				(fp_chunkability)?fp_beo->frag_b_id:NULL_FRAG_ID,
				(tp_chunkability)?tp_beo->frag_b_id:NULL_FRAG_ID
			);
		}
	}
	

	bool ChunkGraphBuilder::isChunkable(
		BestEdgeOverlap *beo,
		BestOverlapGraph& bovlg){
		
		BestEdgeOverlap *b_beo;
		b_beo=bovlg.getBestEdgeOverlap(a_beo.frag_b_id, overlap_end(a_beo.type, 'b'));

		// Check for buddy
		if(b_beo->frag_b_id == frag_a_id){
			// Check for unambiguous in degrees
			if(b_beo->in_degree==1 && a_beo->in_degree==1){
				return(true);
			}else{
				return(false);
			}
		}
	}

	bool ChunkGraphBuilder::isChunkable(
		iuid frag_a_id, fragment_end_type which_end, 
		BestOverlapGraph& bovlg){

		// Translate from fragment id to best overlap
		BestEdgeOverlap *a_beo;
		a_beo=bovlg.getBestEdgeOverlap(frag_a_id, which_end);

		return(isChunkable(a_beo));

	}

} //AS_BOG namespace



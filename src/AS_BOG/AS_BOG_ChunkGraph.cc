
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
 * $Id: AS_BOG_ChunkGraph.cc,v 1.12 2007-03-14 13:01:06 eliv Exp $
 * $Revision: 1.12 $
*/

static char AS_BOG_CHUNK_GRAPH_CC_CM_ID[] = "$Id: AS_BOG_ChunkGraph.cc,v 1.12 2007-03-14 13:01:06 eliv Exp $";

//  System include files

#include<iostream>
#include<vector>
#include<set>
#include<cassert>

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_ChunkGraph.hh"

namespace AS_BOG{

	//////////////////////////////////////////////////////////////////////////////
	// ChunkGraph

	ChunkGraph::ChunkGraph(void){
		_chunkable_array=NULL;
		_edgePathLen=NULL;
		_chunk_lengths=NULL;
		_max_fragments=0;
	}

	//////////////////////////////////////////////////////////////////////////////
	
	ChunkGraph::~ChunkGraph(void){
		if(_edgePathLen != NULL){
			delete[] _edgePathLen;
		}
		if(_chunkable_array != NULL){
			delete[] _chunkable_array;
		}
		if(_chunk_lengths != NULL){
            delete[] _chunk_lengths;
        }
	}


	//////////////////////////////////////////////////////////////////////////////
	// Accessor Functions

	iuid ChunkGraph::getChunking(iuid src_frag_id, fragment_end_type whichEnd){
		if(FIVE_PRIME == whichEnd){
			return(_chunkable_array[src_frag_id].five_prime);
		}else if(THREE_PRIME == whichEnd){
			return(_chunkable_array[src_frag_id].three_prime);
		}else{
			fprintf(stderr, "Bad fragment_end_type.\n");
			exit(-1);
		}
	}

	//////////////////////////////////////////////////////////////////////////////

	void ChunkGraph::getChunking(
		iuid src_frag_id,
		iuid& five_prime_dst_frag_id, iuid& three_prime_dst_frag_id){

		five_prime_dst_frag_id=_chunkable_array[src_frag_id].five_prime;
		three_prime_dst_frag_id=_chunkable_array[src_frag_id].three_prime;
	}

	//////////////////////////////////////////////////////////////////////////////

	void ChunkGraph::setChunking(
		iuid src_frag_id,
		iuid five_prime_dst_frag_id, iuid three_prime_dst_frag_id){

		_chunkable_array[src_frag_id].five_prime=five_prime_dst_frag_id;
		_chunkable_array[src_frag_id].three_prime=three_prime_dst_frag_id;
		
	}

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	void ChunkGraph::build(BestOverlapGraph *inBovlg){
		// This will go through all the nodes in the bovlg and
		//   produce the ChunkGraph
        bovlg = inBovlg;

		// Initialize the chunk graph if necessary
		iuid num_frags=bovlg->getNumFragments();
		if(_edgePathLen != NULL){
			delete[] _edgePathLen;
		}
		_edgePathLen = new iuid[(num_frags+1)*2];
        memset(_edgePathLen, 0, sizeof(iuid)*(num_frags+1)*2);
		if(_chunkable_array != NULL){
			delete[] _chunkable_array;
		}
		_chunkable_array=new _chunk_unit_struct[num_frags+1];
		if(_chunk_lengths != NULL){
			delete[] _chunk_lengths;
		}
		_chunk_lengths = new _chunk_length[num_frags];
		_max_fragments=num_frags;

		// Go through every fragment in the BOG
		iuid frag_id;
		for(frag_id=1; frag_id<=num_frags; frag_id++){

			// Grab 5' end, and check for chunkability
			BestEdgeOverlap *fp_beo=
				bovlg->getBestEdgeOverlap(frag_id, FIVE_PRIME);
			bool fp_chunkability=isChunkable(frag_id, FIVE_PRIME);
			
			// Grab 3' end, and check for chunkability
			BestEdgeOverlap *tp_beo=
				bovlg->getBestEdgeOverlap(frag_id, THREE_PRIME);
			bool tp_chunkability=isChunkable(frag_id, THREE_PRIME);

			// Save chunkability
			setChunking(frag_id, 
				(fp_chunkability)?fp_beo->frag_b_id:NULL_FRAG_ID,
				(tp_chunkability)?tp_beo->frag_b_id:NULL_FRAG_ID
			);
		}
		for(frag_id=1; frag_id<=num_frags; frag_id++){
            iuid fpCnt,tpCnt;
//            fpCnt = countChunkWidth(frag_id, FIVE_PRIME);
            fpCnt  = countFullWidth( frag_id, FIVE_PRIME);
////            tpCnt = countChunkWidth(frag_id, THREE_PRIME);
            tpCnt  = countFullWidth( frag_id, THREE_PRIME);
            
            _chunk_lengths[frag_id-1].fragId = frag_id;
            _chunk_lengths[frag_id-1].cnt    = MIN(fpCnt,tpCnt);
		}
        qsort( _chunk_lengths, num_frags, sizeof(_chunk_length), &ChunkGraph::sortChunkLens);
        fprintf(stderr,"Top chunkLength frgs %d cnt %d and %d cnt %d\n",
                _chunk_lengths[0].fragId, _chunk_lengths[0].cnt,
                _chunk_lengths[1].fragId, _chunk_lengths[1].cnt
               );

        delete[] _edgePathLen;
        _edgePathLen = NULL;
	}

    short ChunkGraph::countChunkWidth(iuid frag, fragment_end_type end) {
        short cnt = 0;
        std::set<iuid> seen;
        FragmentEnd anEnd( frag, end);
        while (cnt < FRAG_WALK_NUM && frag != NULL_FRAG_ID &&
                seen.find(frag) == seen.end())
        {
            cnt++;
            seen.insert(frag);
            anEnd = followPath(anEnd);
            frag = anEnd.id;
        }
        return cnt;
    }
    FragmentEnd ChunkGraph::followPath( FragmentEnd anEnd )
    {
        iuid frag = anEnd.id;
        if (anEnd.end == FIVE_PRIME) 
            frag = _chunkable_array[frag].five_prime ;
        else 
            frag = _chunkable_array[frag].three_prime ;

        // advances to the next fragment, opposite end
        bovlg->followOverlap( &anEnd );
        assert( frag == NULL_FRAG_ID || frag == anEnd.id ); 
        anEnd.id = frag;
        return anEnd;
    }
    iuid ChunkGraph::countFullWidth(iuid firstFrag, fragment_end_type end)
    {
        iuid index = firstFrag * 2 + end;
        if ( _edgePathLen[index] != 0 )
            return _edgePathLen[index] != 0;
        iuid cnt = 0;
        std::set<FragmentEnd> seen;
        FragmentEnd fragEnd( firstFrag, end );
        do {
            seen.insert(fragEnd);
            _edgePathLen[index] = ++cnt;
            FragmentEnd nextEnd;
            nextEnd = followPath( fragEnd );
            fragEnd = nextEnd;
            index = fragEnd.id * 2 + fragEnd.end;
        }
        while (fragEnd.id != NULL_FRAG_ID && _edgePathLen[index] == 0);
        if (fragEnd.id == NULL_FRAG_ID) {
            return cnt;
        }
        // if we end because of a circle, mark points in circle same cnt
        if (seen.find( fragEnd ) != seen.end()) {
            iuid circleLen = cnt - _edgePathLen[index];
            fprintf(stderr,"Circle len %d frag %d end %d\n",
                    circleLen, fragEnd.id, fragEnd.end);
            FragmentEnd currEnd = fragEnd;
            do {
                seen.erase( currEnd );
                _edgePathLen[index] = circleLen;
                bovlg->followOverlap( &currEnd );
                index = currEnd.id * 2 + currEnd.end;
            }
            while (!(fragEnd == currEnd));
        } else {
            // we encountered an already counted frag, so use it's length
            cnt += _edgePathLen[index];
        }
        // if the whole thing wasn't a circle we still have a linear path left
        iuid max = cnt;
        if (!seen.empty()) {
            FragmentEnd currEnd(firstFrag, end);
            do {
                _edgePathLen[currEnd.id*2+currEnd.end] = cnt--;
                bovlg->followOverlap( &currEnd );
            } while (currEnd.id != fragEnd.id);
        }
        return max;
    }
    int ChunkGraph::sortChunkLens( const void *a, const void *b) {
        struct _chunk_length *frg1 = (struct _chunk_length*)a;
        struct _chunk_length *frg2 = (struct _chunk_length*)b;
/*        short cnt1 = min(frg1->fpCnt,frg1->tpCnt);
        short cnt2 = min(frg2->fpCnt,frg2->tpCnt);
        if ( cnt1 != cnt2 )
            return cnt2 - cnt1;
*/
        if ( frg1->cnt != frg2->cnt )
            return frg2->cnt - frg1->cnt;
        else 
            return frg1->fragId - frg2->fragId;
    }
	//////////////////////////////////////////////////////////////////////////////

    iuid ChunkGraph::nextFragByChunkLength() {
        static iuid pos = 0;
        if (pos < _max_fragments)
            return _chunk_lengths[pos++].fragId;
        else {
            pos = 0;
            return pos;
        }
    }
	
	//////////////////////////////////////////////////////////////////////////////

	bool ChunkGraph::isChunkable(
		iuid frag_a_id, fragment_end_type which_end ) {
		// Given an edge (based on fragment ID and end), determines by looking at
		//   what the edge overlaps, whether the overlap is unambiguous.

		// Translate from fragment id to best overlap
		BestEdgeOverlap *a_beo;
		a_beo=bovlg->getBestEdgeOverlap(frag_a_id, which_end);

		if(a_beo->in_degree==1){

			// Get B's overlap information
			BestEdgeOverlap *b_beo;
			b_beo=bovlg->getBestEdgeOverlap(a_beo->frag_b_id, a_beo->bend);

			// If b points to a, as well
			if((b_beo->frag_b_id == frag_a_id) &&
			   (b_beo->bend == which_end)
			){
				// Check for unambiguous in degrees
				if(b_beo->in_degree==1){
					return(true);
				}else{
					return(false);
				}
			}else{
				return(false);
			}
		}else{
			return(false);
		}

	}
	bool PromiscuousChunkGraph::isChunkable(
		iuid frag_a_id, fragment_end_type which_end) { 

		// Translate from fragment id to best overlap
		BestEdgeOverlap *a_beo=bovlg->getBestEdgeOverlap(frag_a_id, which_end);
        return true;
        BestEdgeOverlap *b_beo=bovlg->getBestEdgeOverlap(a_beo->frag_b_id, a_beo->bend);

        // If b points to a, as well
        if((b_beo->frag_b_id == frag_a_id) && (b_beo->bend == which_end)) {
            return(true);
		}else{
			return(false);
		}
	}
	//////////////////////////////////////////////////////////////////////////////

	void ChunkGraph::printFrom(iuid begin, iuid end){
		iuid i;
		for(i=begin; i<=end; i++){
			std::cout << 
				_chunkable_array[i].five_prime << " <= " <<
				i << " => " <<
				_chunkable_array[i].three_prime <<
				std::endl;
		}

	}

	//////////////////////////////////////////////////////////////////////////////

	long ChunkGraph::getNumFragments(void){
		return(_max_fragments);		
	}

	//////////////////////////////////////////////////////////////////////////////

	long ChunkGraph::countSingletons(void){
		iuid i;
		long num_singletons=0;
		for(i=1; i<=_max_fragments; i++){
			if((_chunkable_array[i].five_prime == NULL_FRAG_ID) &&
			   (_chunkable_array[i].three_prime == NULL_FRAG_ID)){
				num_singletons++;
			}
		}
		return(num_singletons);
	}

	//////////////////////////////////////////////////////////////////////////////

	std::ostream& operator<< (std::ostream& os, ChunkGraph &cg){
		
		iuid i;
		for(i=1; i<=cg._max_fragments; i++){
			os << 
				cg._chunkable_array[i].five_prime << " " <<
				"[" << i << "] " <<
				cg._chunkable_array[i].three_prime << " " <<
				std::endl;

		}

	}

	//////////////////////////////////////////////////////////////////////////////

	void ChunkGraph::checkInDegree(){
		
		int *fivep_indegree_arr;
		int *threep_indegree_arr;
		long frag_id;

		long num_frags=bovlg->getNumFragments();
		fivep_indegree_arr=new int[num_frags+1];
		threep_indegree_arr=new int[num_frags+1];

		for(frag_id=0; frag_id<=num_frags; frag_id++){
			fivep_indegree_arr[frag_id]=0;
			threep_indegree_arr[frag_id]=0;
		}

		if(0){
		BestContainmentMap::iterator itr;
		for(
		    itr=bovlg->_best_containments.begin();
		    itr!=bovlg->_best_containments.end();
		    itr++){

			std::cerr << itr->first << " Contained in " << itr->second.container << 
				std::endl;
		}
		}

		BestEdgeOverlap	*fivep_overlap_ptr, *threep_overlap_ptr;
		for(frag_id=1; frag_id<=num_frags; frag_id++){

			fivep_overlap_ptr=bovlg->getBestEdgeOverlap(
				frag_id, FIVE_PRIME);
			threep_overlap_ptr=bovlg->getBestEdgeOverlap(
				frag_id, THREE_PRIME);

			if(bovlg->_best_containments.find(frag_id)==
			    bovlg->_best_containments.end()){

				if(bovlg->_best_containments.find(fivep_overlap_ptr->frag_b_id)!=
				    bovlg->_best_containments.end()){
					std::cerr << frag_id << "(5') best olaps w/ containee:" << 
					fivep_overlap_ptr->frag_b_id << " contained in " << 
		(bovlg->_best_containments[threep_overlap_ptr->frag_b_id].container) <<
					std::endl;
				}
				if(bovlg->_best_containments.find(threep_overlap_ptr->frag_b_id)!=
				    bovlg->_best_containments.end()){
					std::cerr << frag_id << "(3') best olaps w/ containee:" << 
					threep_overlap_ptr->frag_b_id << " contained in " << 
		(bovlg->_best_containments[threep_overlap_ptr->frag_b_id].container) <<
					std::endl;
				}


				if(fivep_overlap_ptr->bend==FIVE_PRIME){
					fivep_indegree_arr[fivep_overlap_ptr->frag_b_id]++;
				}
				if(fivep_overlap_ptr->bend==THREE_PRIME){
					threep_indegree_arr[fivep_overlap_ptr->frag_b_id]++;
				}
				if(threep_overlap_ptr->bend==FIVE_PRIME){
					fivep_indegree_arr[threep_overlap_ptr->frag_b_id]++;
				}
				if(threep_overlap_ptr->bend==THREE_PRIME){
					threep_indegree_arr[threep_overlap_ptr->frag_b_id]++;
				}
			}
		}
		

		for(frag_id=1; frag_id<=num_frags; frag_id++){
			BestEdgeOverlap *fp=bovlg->getBestEdgeOverlap(frag_id, FIVE_PRIME);
			BestEdgeOverlap *tp=bovlg->getBestEdgeOverlap(frag_id, THREE_PRIME);

			std::cerr << frag_id << " " << fivep_indegree_arr[frag_id] << " " <<
				threep_indegree_arr[frag_id] << " / " << fp->in_degree << " " << tp->in_degree <<
				" [" << bovlg->getBestContainer(frag_id) << "]" <<
				std::endl;

			//fp->in_degree=fivep_indegree_arr[frag_id];
			//tp->in_degree=threep_indegree_arr[frag_id];
		}


		delete[] fivep_indegree_arr;
		delete[] threep_indegree_arr;

	}

	//////////////////////////////////////////////////////////////////////////////

} //AS_BOG namespace



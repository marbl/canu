
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
 * $Id: AS_BOG_UnitigGraph.cc,v 1.7 2005-12-16 19:45:47 kli1000 Exp $
 * $Revision: 1.7 $
*/

//static char AS_BOG_UNITIG_GRAPH_CC_CM_ID[] = "$Id: AS_BOG_UnitigGraph.cc,v 1.7 2005-12-16 19:45:47 kli1000 Exp $";
static char AS_BOG_UNITIG_GRAPH_CC_CM_ID[] = "gen> @@ [0,0]";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include <float.h>
#include <stdlib.h>

extern "C" {
	#include "AS_global.h"
}

namespace AS_BOG{

	//////////////////////////////////////////////////////////////////////////////

	std::ostream& operator<< (std::ostream& os, FragmentPositionMap *fpm_ptr){

		FragmentPositionMap::iterator itr;
		for(
		    itr=fpm_ptr->begin();
		    itr!=fpm_ptr->end();
		    itr++){

			os << "[" << itr->first << 
				"]\t" << itr->second.begin << " - " << itr->second.end << 
				std::endl;
		}
		return(os);
	}

	//////////////////////////////////////////////////////////////////////////////

	void UnitigGraph::build(ChunkGraph *cg_ptr, BestOverlapGraph *bog_ptr, long num_rand_frags, long genome_size){

		iuid num_frags=cg_ptr->getNumFragments();
		iuid unitig_id=1;

		iuid frag_idx;
		iuid fp_dst_frag_id, tp_dst_frag_id;

		// Initialize where we've been to nowhere; "Do not retraverse list"
		std::map<iuid, bool> visited_map;
		visited_map.clear();

		ContainerMap *cntnrmap_ptr;
		BestContainmentMap *bcmp_ptr=&(bog_ptr->_best_containments);

		cntnrmap_ptr=_build_container_map(bog_ptr);

		// Step through all the fragments 
		std::cerr << "Building Unitigs from " << num_frags << " fragments.\n"; 
		for(frag_idx=1; frag_idx<=num_frags; frag_idx++){

			//std::cerr << "Working on " << frag_idx << std::endl; 

			// Check the map to so we don't visit a unitig twice (once from
			//   both ends)
			if(visited_map.find(frag_idx)==visited_map.end() && 
				bog_ptr->_best_containments.find(frag_idx)==bog_ptr->_best_containments.end()){

				cg_ptr->getChunking(
					frag_idx, 
					fp_dst_frag_id, 
					tp_dst_frag_id);

				if(  // If the 5' end is NULL but 3' end isn't, or vice versa
				   ((fp_dst_frag_id == NULL_FRAG_ID &&
				     tp_dst_frag_id != NULL_FRAG_ID ) ||
				    (tp_dst_frag_id == NULL_FRAG_ID &&
				     fp_dst_frag_id != NULL_FRAG_ID ) ||
				    (tp_dst_frag_id == NULL_FRAG_ID &&
				     fp_dst_frag_id == NULL_FRAG_ID ) 
				)){

					if(!(unitig_id%1000)){
						std::cerr << ".";
					}

					// Allocated a new unitig node
					Unitig *utg=new Unitig;

					// Store the dovetails
					if(fp_dst_frag_id == NULL_FRAG_ID && tp_dst_frag_id != NULL_FRAG_ID){
						utg->dovetail_path_ptr=
							_extract_dovetail_path(
								frag_idx, FORWARD, cg_ptr, bog_ptr);
					}else if(tp_dst_frag_id == NULL_FRAG_ID && fp_dst_frag_id != NULL_FRAG_ID){
						utg->dovetail_path_ptr=
							_extract_dovetail_path(
								frag_idx, REVERSE, cg_ptr, bog_ptr);
					}else{
						utg->dovetail_path_ptr=
							_extract_dovetail_path(
								frag_idx, FORWARD, cg_ptr, bog_ptr);
					}

					// Get the fragment ID of the last fragment in the
					//   dovetail, then store it in the "do not retraverse
					//   list"
					iuid last_frag_id_in_dovetail=
						(*utg->dovetail_path_ptr)
						[(*utg->dovetail_path_ptr).size()-1].frag_id;
					visited_map[last_frag_id_in_dovetail]=true;

					// Store the containees
					utg->contained_frags_ptr=
						_extract_containees(
							utg->dovetail_path_ptr,
							cntnrmap_ptr);
	
					//std::cerr << "Containment Extracted. " << std::endl;

					utg->frag_pos_map_ptr=
						utg->computeFragmentPositions();
					//std::cerr << utg->frag_pos_map_ptr << std::endl;

					// Set id
					utg->id=unitig_id;
					unitig_id++;

					// Store unitig in unitig graph
					unitigs.push_back(utg);

				}else{
					// We are either in the middle of a dovetail sequence,
					// it is a singleton, or we are going from 3' to 5'
					// across a dovetail.

				}
			}
		}

		std::cerr << "Setting Global Arrival Rate.\n";
		Unitig static_proxy;
		static_proxy.setGlobalArrivalRate(getGlobalArrivalRate(num_rand_frags, genome_size));

		std::cerr << std::endl << "There were " << unitigs.size() << " unitigs generated.\n";
		
	}

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	DoveTailPath *UnitigGraph::_extract_dovetail_path(
		iuid src_frag_id, orientation_type ori, ChunkGraph *cg_ptr, BestOverlapGraph *bog_ptr){

	// Note:  I only need BestOverlapGraph for it's frag_len and olap_length

		DoveTailPath *dtp_ptr=new DoveTailPath;

		iuid fp_frag_id, tp_frag_id;
		iuid current_frag_id=src_frag_id;
		iuid next_frag_id;
		fragment_end_type overlap_end;
		orientation_type travel_dir;

		travel_dir=ori;	
		next_frag_id=
			bog_ptr->
			getBestEdgeOverlap(src_frag_id, (travel_dir == FORWARD)? THREE_PRIME : FIVE_PRIME)->frag_b_id;

		// Start walking the chunk path
		orientation_type next_trav_dir;

		//std::cerr<<"Working on: "<<src_frag_id<< " Dir: " << travel_dir <<std::endl;

		iuid last_frag_id;
		do{
			// Store the current fragment into dovetail path
			DoveTailNode dt_node;
			dt_node.frag_id=current_frag_id;
			dt_node.ori=travel_dir;
			dt_node.olap_len=
			    bog_ptr->getBestEdgeOverlap(
				current_frag_id, 
				(travel_dir == FORWARD)? THREE_PRIME : FIVE_PRIME)->
				olap_len;
			dt_node.frag_len=bog_ptr->fragLen(current_frag_id);

			dtp_ptr->push_back(dt_node);

			// Get next fragment
			_follow_to_next_fragment(
				current_frag_id, travel_dir,
				next_frag_id, next_trav_dir,
				cg_ptr);
			
			// Set current to next
			current_frag_id=next_frag_id;
			travel_dir=next_trav_dir;

		}while(current_frag_id != NULL_FRAG_ID);


		return(dtp_ptr);
			
	}

	//////////////////////////////////////////////////////////////////////////////

	void UnitigGraph::_follow_to_next_fragment(
	    iuid src_frag_id, orientation_type src_end,
	    iuid& dst_frag_id, orientation_type& dst_end,
	    ChunkGraph *cg_ptr
	){

		//  5  A  3        5  B  3
		//  O=====O -----> O=====O
		//
		//  Given A's ID and end, will return B's ID and the end A is not pointing to.

		iuid dst_fp_frag_id, dst_tp_frag_id;	//Wht the source points at

		cg_ptr->getChunking(src_frag_id, dst_fp_frag_id, dst_tp_frag_id);
		//std::cerr << "[" << src_frag_id << "] " << dst_fp_frag_id << " " <<dst_tp_frag_id << 
		//	" (" << src_end <<  ")" << std::endl;

		// Follow the next node based on the suggested direction
		if(src_end == FORWARD){
			dst_frag_id = dst_tp_frag_id;
		}else if(src_end == REVERSE){
			dst_frag_id = dst_fp_frag_id;
		}

		// See what the next node is pointing to
		iuid dst_dst_fp_frag_id, dst_dst_tp_frag_id;	//What the destination points at
		cg_ptr->getChunking(dst_frag_id, dst_dst_fp_frag_id, dst_dst_tp_frag_id);

		// If the next node points back to the source node, we want to go the 
		// opposite direction
		if(dst_dst_fp_frag_id == src_frag_id){
			dst_end = FORWARD;
		}else if(dst_dst_tp_frag_id == src_frag_id){
			dst_end = REVERSE;
		}else{
			if(dst_frag_id != NULL_FRAG_ID){
				std::cerr << 
				    "Error, dst node does not point back at src node." << std::endl;
				std::cerr << 
				    "Src: " <<  src_frag_id << " End: " << src_end << std::endl;
				std::cerr << 
				    "Dst: " << dst_dst_fp_frag_id << " " << dst_dst_tp_frag_id << std::endl;
			}else{
				//std::cerr << "End of unitig reached.\n";
			}
		}	
		
		return;
	}

	//////////////////////////////////////////////////////////////////////////////

	ContainerMap *UnitigGraph::_build_container_map(BestOverlapGraph *bog_ptr){

	// I don't need the whole overlap graph, just the BestContainmentMap and
	//  the fragment lengths

		ContainerMap *cmptr=new ContainerMap;
		iuid container_id, containee_id;
		long containees=0;

		// I'm getting the pointer here, so in the future you can just
		// pass in the pointer to a BestContainmentMap instead of the whole
		// BestOverlapGraph
		BestContainmentMap *cntrmap_ptr=&(bog_ptr->_best_containments);

		BestContainmentMap::iterator bstcnmap_itr;
		for(
		    bstcnmap_itr=cntrmap_ptr->begin(); 
		    bstcnmap_itr!=cntrmap_ptr->end(); 
		    bstcnmap_itr++){

			ContaineeNode ctnee;
			ctnee.frag_id=bstcnmap_itr->first;
			ctnee.same_ori=bstcnmap_itr->second.sameOrientation;
			ctnee.olap_offset=bstcnmap_itr->second.a_hang;
			ctnee.frag_len=bog_ptr->fragLen(ctnee.frag_id);

			container_id=bstcnmap_itr->second.container;

			((*cmptr)[container_id]).push_back(ctnee);
			containees++;
		}
		
		std::cerr << "BestContainments Inverted into ContainerMap." << std::endl;
		std::cerr << "Num containees: " << containees << std::endl;
		std::cerr << "Num containers: " << cmptr->size() << std::endl;
		
		return(cmptr);
	}

	//////////////////////////////////////////////////////////////////////////////

	ContainerMap *UnitigGraph::_extract_containees(
		DoveTailPath *dtp_ptr, ContainerMap *all_cntnrmap_ptr){

		ContainerMap *ctnrmap_ptr=new ContainerMap;

		DoveTailPath::iterator dtp_itr;
		for(
		    dtp_itr=dtp_ptr->begin();
		    dtp_itr!=dtp_ptr->end();
		    dtp_itr++){
			
			// Copy the container entry over if it is found in the dovetail
			if(all_cntnrmap_ptr->find(dtp_itr->frag_id)!=all_cntnrmap_ptr->end()){
				(*ctnrmap_ptr)[dtp_itr->frag_id]=
					(*all_cntnrmap_ptr)[dtp_itr->frag_id];
			}
	
		}
		
		return(ctnrmap_ptr);

	}

	//////////////////////////////////////////////////////////////////////////////

	std::ostream& operator << (std::ostream& os, Unitig& utg){
		
		os << "Dovetails:" << std::endl;

		DoveTailPath::iterator dt_itr;
		for(dt_itr=utg.dovetail_path_ptr->begin();
			dt_itr!=utg.dovetail_path_ptr->end(); 
			dt_itr++){
			
			os << "  " << dt_itr->frag_id << " [" << dt_itr->ori << "]" << std::endl;
		}

		os << "Containments:" << std::endl;
	
		ContainerMap::iterator ctmp_itr;
		for(ctmp_itr=utg.contained_frags_ptr->begin();
			ctmp_itr!=utg.contained_frags_ptr->end(); 
			ctmp_itr++){
			
			os << "  " << ctmp_itr->first << std::endl;

			ContaineeList::iterator frgl_itr;
			for(frgl_itr=ctmp_itr->second.begin();
			    frgl_itr!=ctmp_itr->second.end(); 
			    frgl_itr++){
			       os << "    " << frgl_itr->frag_id << std::endl;
			}
		}

		os << "avgRho: " << utg.getAvgRho() << std::endl;
		os << "covStat: " << utg.getCovStat() << std::endl;
		os << "numFrags: " << utg.getNumFrags() << std::endl;
		os << "numRandomFrags: " << utg.getNumRandomFrags() << std::endl;
		os << "length: " << utg.getLength() << std::endl;

		return(os);

	}

	//////////////////////////////////////////////////////////////////////////////

	std::ostream& operator << (std::ostream& os, UnitigGraph& utgrph){
		
		UnitigVector::iterator utg_itr;
		iuid num_utgs=0;
		for(utg_itr=utgrph.unitigs.begin();
			utg_itr!=utgrph.unitigs.end();
			utg_itr++){
		
			std::cerr << num_utgs << std::endl;
			std::cerr << (**utg_itr) << std::endl;
			num_utgs++;
		
		}

		return(os);
	}

	//////////////////////////////////////////////////////////////////////////////

	float UnitigGraph::getGlobalArrivalRate(long total_random_frags_in_genome, long genome_size){
		
		float _globalArrivalRate;

		//  If the genome size has not been specified, estimate the GAR.
		if(genome_size == 0){

			float total_rho=0;
			float total_arrival_frags=0;

			// Go through all the unitigs to sum rho and unitig arrival frags
			UnitigVector::iterator iter;
			for(
			    iter=unitigs.begin();
			    iter!=unitigs.end();
			    iter++){
				
				total_rho += (*iter)->getAvgRho();
				float unitig_random_frags = (*iter)->getNumRandomFrags();
				total_arrival_frags += (unitig_random_frags > 0)?
					unitig_random_frags - 1 :
					0;
			}

			// Estimate GAR
			_globalArrivalRate = (total_rho > 0)?
				(total_arrival_frags / total_rho):
				0;
		}else{
			// Compute actual GAR
			_globalArrivalRate = (float)total_random_frags_in_genome / (float)genome_size;
		}

		return(_globalArrivalRate);		

	}

	//////////////////////////////////////////////////////////////////////////////

	Unitig::Unitig(void){
		// Initialize values to unlikely values
		_globalArrivalRate=-1;
		_covStat=FLT_MAX;
		_length=-1;
		_numFrags=-1;
		_numRandomFrags=-1;
		dovetail_path_ptr=NULL;
		contained_frags_ptr=NULL;
		frag_pos_map_ptr=NULL;
	}

	//////////////////////////////////////////////////////////////////////////////

	Unitig::~Unitig(void){
		if(dovetail_path_ptr!=NULL) delete dovetail_path_ptr;
		if(contained_frags_ptr!=NULL) delete contained_frags_ptr;
		if(frag_pos_map_ptr!=NULL) delete frag_pos_map_ptr;
		
	}

	//////////////////////////////////////////////////////////////////////////////

	float Unitig::getAvgRho(void){

		if(_avgRho==-1){
			return(_avgRho);
		}
		
		// We will compute the average rho.
		//
		// Since rho is the length(unitig) - length(last fragment),
		//   and the length(last fragment) is ambiguous depending on which
		//   direction we are walking the unitig from.  We will take the average 
		//   of the rhos through both directions.

		DoveTailPath::iterator dtp_iter;

		// Get first fragment's length
		dtp_iter=dovetail_path_ptr->begin();
		long first_frag_len = dtp_iter->frag_len;

		// Get last fragment's length
		dtp_iter=dovetail_path_ptr->end();
		dtp_iter--;
		long last_frag_len = dtp_iter->frag_len;

		// Get average of first and last fragment lengths
		double avg_frag_len = (last_frag_len + first_frag_len)/2.0;
		
		// Compute average rho
		long unitig_length=getLength();
		_avgRho = unitig_length - avg_frag_len;
		
		return(_avgRho);
	}

	//////////////////////////////////////////////////////////////////////////////

	float Unitig::_globalArrivalRate;

	void Unitig::setGlobalArrivalRate(float global_arrival_rate){
		_globalArrivalRate=global_arrival_rate;
	}

	//////////////////////////////////////////////////////////////////////////////

	float Unitig::getCovStat(void){

		const float ln2=0.69314718055994530941723212145818;

		// Note that we are using numFrags in this calculation.
		//   If the fragments in the unitig are not randomly sampled
		//   from the genome, this calculation will be wrong.
		//   Ie. if fragments being assembled are from a locally
		//   sequenced batch, the region may look artificially repetitive.
		//   The value should really be "number of randomly sampled
		//   fragments in the unitig".

		if(_globalArrivalRate==-1){
			std::cerr << "You have not set the _globalArrivalRate variable." << std::endl;
		}

		if(_covStat == FLT_MAX){
			if(_globalArrivalRate > 0.0){
				_covStat = (getAvgRho() * _globalArrivalRate) - 
					(ln2 * (getNumFrags() -1));
			}else{
				_covStat = 0.0;
			}
		}

		return(_covStat);

	}

	//////////////////////////////////////////////////////////////////////////////

	long Unitig::getLength(void){
	// The length of the unitig is just the position of the fragment 
	//   with the large begin or end value 

		if(_length!=-1){
			return(_length);
		}
		
		long max_pos=-1;

		if(frag_pos_map_ptr->size()==0){
			std::cerr << "This Unitig has an empty FragmentPositionMap." << std::endl;	
		}

		FragmentPositionMap::iterator fpm_itr;

		for(
		    fpm_itr=frag_pos_map_ptr->begin();
		    fpm_itr!=frag_pos_map_ptr->end();
		    fpm_itr++){

			interval intrvl=fpm_itr->second;

			if(max_pos<intrvl.begin){
				max_pos=intrvl.begin;
			}
			if(max_pos<intrvl.end){
				max_pos=intrvl.end;
			}
		}

		_length=max_pos;
		return(_length);

	}

	//////////////////////////////////////////////////////////////////////////////

	long Unitig::getNumFrags(void){

		if(_numFrags!=-1){
			return(_numFrags);
		}

		long num_dovetailing_frags=0;
		long num_contained_frags=0;

		// Count dovetailing fragments
		DoveTailPath::iterator dt_itr;
		for(dt_itr=dovetail_path_ptr->begin();
			dt_itr!=dovetail_path_ptr->end(); 
			dt_itr++){
			num_dovetailing_frags++;	
		}

		// Count contained fragments
		ContainerMap::iterator ctmp_itr;
		for(ctmp_itr=contained_frags_ptr->begin();
			ctmp_itr!=contained_frags_ptr->end(); 
			ctmp_itr++){

			ContaineeList::iterator clist_itr;
			for(clist_itr=ctmp_itr->second.begin();
			    clist_itr!=ctmp_itr->second.end(); 
			    clist_itr++){

				num_contained_frags++;
			}
		}

		// Total frags are the sum
		_numFrags = num_dovetailing_frags + num_contained_frags;
		return(_numFrags);

	}

	//////////////////////////////////////////////////////////////////////////////

	long Unitig::getNumRandomFrags(void){
	// This is a placeholder, random frags should not contain guides, or other
	//   fragments that are not randomly sampled across the whole genome.
		
		if(_numRandomFrags!=-1){
			return(_numRandomFrags);
		}	

		_numRandomFrags=getNumFrags();
		return(_numRandomFrags);
	}

	//////////////////////////////////////////////////////////////////////////////

	FragmentPositionMap *Unitig::computeFragmentPositions(void){

		FragmentPositionMap *fpmp=new FragmentPositionMap;
		
		long frag_ins_begin;
		long frag_ins_end;

		//std::cerr << "Positioning dovetails." << std::endl;
		// place dovetails in a row
		frag_ins_begin=0;
		DoveTailPath::iterator dt_itr;
		for(dt_itr=dovetail_path_ptr->begin();
		    dt_itr!=dovetail_path_ptr->end(); 
		    dt_itr++){

			frag_ins_end=frag_ins_begin + dt_itr->frag_len;

			interval intvl;
			if(dt_itr->ori==REVERSE){
				intvl.begin=frag_ins_end;
				intvl.end=frag_ins_begin;
			}else if(dt_itr->ori==FORWARD){
				intvl.begin=frag_ins_begin;
				intvl.end=frag_ins_end;
			}else{
				std::cerr << "Unknown orientation type for frag id:"
				     << dt_itr->frag_id << std::endl;
			}

			// Store interval in the map
			(*fpmp)[dt_itr->frag_id]=intvl;
			
			// Prep the start position of the next fragment
			frag_ins_begin = frag_ins_end - dt_itr->olap_len;
			
		}

		//std::cerr << "Positioning containees." << std::endl;
		// Place containees into unitig
		ContainerMap::iterator ctmp_itr;
		// For each container
		for(
		    ctmp_itr=contained_frags_ptr->begin();
		    ctmp_itr!=contained_frags_ptr->end(); 
		    ctmp_itr++){

			iuid container_id=ctmp_itr->first;
			interval cntnr_intvl=(*fpmp)[container_id];

			// For each containee
			ContaineeList::iterator cntee_itr;
			for(cntee_itr=ctmp_itr->second.begin();
			    cntee_itr!=ctmp_itr->second.end(); 
			    cntee_itr++){

				interval cntee_intvl;
				// Compute assuming that containee is the same orientation as container
				if(cntnr_intvl.begin < cntnr_intvl.end){

				// Container is in forward direction
				//
				// |=========F0============|
				// 0          |=============CER=============>|
				//                    |=====CEE=====>|
				//
				// |---Cro----|
				//            |--Ceo--|
				// |-------Cep--------|
				//
				// Cro = container offset from beginning of unitig = cntnr_intvl.begin
				// Ceo = containee offset from 5' end of container = cntee->olap_offset
				// Cep = containee offset from beginning of unitig = cntee_intvl.begin
				// CEE fragment can be from either orientation since
				//   definition of olap_offset is based on 3' origin.

					cntee_intvl.begin = 
					    cntnr_intvl.begin + cntee_itr->olap_offset;
					cntee_intvl.end   = 
					    cntee_intvl.begin + cntee_itr->frag_len;

					// Make sure the containee doesn't extend past the container
					cntee_intvl.begin = (cntee_intvl.begin < cntnr_intvl.begin)?
						cntnr_intvl.begin:cntee_intvl.begin;
					cntee_intvl.end = (cntee_intvl.end > cntnr_intvl.end)?
						cntnr_intvl.end:cntee_intvl.end;

				}else if(cntnr_intvl.begin > cntnr_intvl.end){

				// Container is in reverse direction
				//
				// |=========F0============|
				// 0          |<============CER==============|
				//                    |<====CEE======|
				//
				// |---Cro----|
				//                                   |--Ceo--|
				// |-------Cep-----------------------|
				//
				// Cro = container offset from beginning of unitig = cntnr_intvl.end
				// Ceo = containee offset from 5' end of container = cntee->olap_offset
				// Cep = containee offset from beginning of unitig = cntee_intvl.end
				// CEE fragment can be from either orientation since
				//   definition of olap_offset is based on 3' origin.

					cntee_intvl.begin = 
					    cntnr_intvl.begin - cntee_itr->olap_offset;
					cntee_intvl.end = 
					    cntee_intvl.begin - cntee_itr->frag_len;

					// Make sure the containee doesn't extend past the container
					cntee_intvl.begin = (cntee_intvl.begin > cntnr_intvl.begin)?
						cntnr_intvl.begin:cntee_intvl.begin;
					cntee_intvl.end = (cntee_intvl.end < cntnr_intvl.end)?
						cntnr_intvl.end:cntee_intvl.end;

				}else{
					std::cerr << 
					    "Error, container size is zero." << std::endl;
				}

				// Swap begin/ends, if containee is not same orientation as container.
				if(!cntee_itr->same_ori){
					long tmp;
					tmp=cntee_intvl.begin;
					cntee_intvl.begin=cntee_intvl.end;
					cntee_intvl.end=tmp;
				}

				(*fpmp)[cntee_itr->frag_id]=cntee_intvl;	
			}
		}

		return(fpmp);

	}

	//////////////////////////////////////////////////////////////////////////////

	void Unitig::freeIUM_Mesg(IntUnitigMesg *ium_ptr){
		delete[] (ium_ptr->f_list);
		delete ium_ptr;
	}

	//////////////////////////////////////////////////////////////////////////////

	int IntMultiPosCmp(const void *a, const void *b){
		IntMultiPos *impa=(IntMultiPos*)a;
		IntMultiPos *impb=(IntMultiPos*)b;
		long aleft = (impa->position.bgn < impa->position.end) ?
			impa->position.bgn : impa->position.end;
		long bleft = (impb->position.bgn < impb->position.end) ?
			impb->position.bgn : impb->position.end;
		if(aleft!=bleft){
			return(aleft - bleft);
		}else{
			if(impa->contained!=0)
				return(1);
			if(impb->contained!=0)
				return(-1);
			return(0);
		}
	}

	typedef struct{
		long ahang;
		long bhang;
		iuid a_id;
	}overlap_info_struct;

	IntUnitigMesg *Unitig::getIUM_Mesg(void){
		
		IntUnitigMesg *ium_mesg_ptr=new IntUnitigMesg;
		IntMultiPos *imp_msg_arr=new IntMultiPos[getNumFrags()];

		// Determine which order to print out fragments
		// Consensus needs contained fragments listed after their
		// containers.  We will order fragments by how they
		// dovetail, filling in containees as their containers
		// get traversed.
		std::map<containee_id, container_id> cntnee_map;

		typedef std::map<iuid, overlap_info_struct> olap_info_map_type;
		olap_info_map_type olap_info_map;

		FragmentList frag_order;
		DoveTailPath::iterator dp_itr;
		iuid prior_frag_id=0;
		for(
			dp_itr=dovetail_path_ptr->begin();
			dp_itr!=dovetail_path_ptr->end();
			dp_itr++
		){
			// Push back the dovetailers
			frag_order.push_back(dp_itr->frag_id);

			olap_info_map[dp_itr->frag_id].a_id=prior_frag_id;

			// Push back the containees
			ContaineeList *ctees=&(*contained_frags_ptr)[dp_itr->frag_id];
			ContaineeList::iterator ctee_itr;
			for(ctee_itr=ctees->begin();
			    ctee_itr!=ctees->end();
			    ctee_itr++){
				frag_order.push_back(ctee_itr->frag_id);
				cntnee_map[ctee_itr->frag_id]=dp_itr->frag_id;
				olap_info_map[ctee_itr->frag_id].a_id=dp_itr->frag_id;
			}

			prior_frag_id=dp_itr->frag_id;
		}


		// Estimate ahangs/bhangs between A and B fragment overlaps
		olap_info_map_type::iterator olap_info_map_itr;
		for(
		    olap_info_map_itr=olap_info_map.begin();
		    olap_info_map_itr!=olap_info_map.end();
		    olap_info_map_itr++){
		
			iuid a_frag_id=olap_info_map_itr->second.a_id;
			iuid b_frag_id=olap_info_map_itr->first;

			if(a_frag_id!=0){
				long a_begin=(*frag_pos_map_ptr)[a_frag_id].begin;
				long a_end  =(*frag_pos_map_ptr)[a_frag_id].end;
				long b_begin=(*frag_pos_map_ptr)[b_frag_id].begin;
				long b_end  =(*frag_pos_map_ptr)[b_frag_id].end;

				if(a_begin>a_end){
					long tmp;
					tmp=a_begin;
					a_begin=a_end;
					a_end=tmp;
				}

				if(b_begin>b_end){
					long tmp;
					tmp=b_begin;
					b_begin=b_end;
					b_end=tmp;
				}

				olap_info_map[b_frag_id].ahang=b_begin-a_begin;
				olap_info_map[b_frag_id].bhang=b_end-a_end;
			}else{
				olap_info_map[b_frag_id].ahang=0;
				olap_info_map[b_frag_id].bhang=0;
			}
		}
	
		// Populate IMP Messages
		FragmentList::iterator fp_itr;
		long i=0;
		
		for(
			fp_itr=frag_order.begin();
			fp_itr!=frag_order.end();
			fp_itr++

		){
			/*FragType*/        imp_msg_arr[i].type=AS_READ;
			/*IntFragment_ID*/  imp_msg_arr[i].ident=*fp_itr;
			# ifdef AS_ENABLE_SOURCE
			/*char* */          imp_msg_arr[i].source="";
			#ifdef i386
			/*int32*/           imp_msg_arr[i].ptrPad1;
			#endif
			# endif
			/*SeqInterval*/     imp_msg_arr[i].position.bgn=(*frag_pos_map_ptr)[*fp_itr].begin;
			/*SeqInterval*/     imp_msg_arr[i].position.end=(*frag_pos_map_ptr)[*fp_itr].end;
			/*IntFragment_ID*/  imp_msg_arr[i].contained=cntnee_map[*fp_itr];
			/*int32*/           imp_msg_arr[i].delta_length=0;
			/*int32*/           imp_msg_arr[i].delta=NULL;

			#ifdef NEW_UNITIGGER_INTERFACE
			/*IntFragment_ID*/  imp_msg_arr[i].ident2=olap_info_map[*fp_itr].a_id;
			/*int32*/           imp_msg_arr[i].ahang=olap_info_map[*fp_itr].ahang;
			/*int32*/           imp_msg_arr[i].bhang=olap_info_map[*fp_itr].bhang;
			#endif

			#ifdef i386
			/*int32*/           imp_msg_arr[i].ptrPad2;
			#endif
			i++;
		}

		//void qsort(void *base, size_t nmemb, size_t size,
                //  int(*compar)(const void *, const void *));
		// Sort by 5' end of fragment position
		qsort(imp_msg_arr, getNumFrags(), sizeof(IntMultiPos), &IntMultiPosCmp);

		// Populate the IUM message with unitig info
		/*IntChunk_ID*/		ium_mesg_ptr->iaccession=id-1;
		#ifdef AS_ENABLE_SOURCE
		/*char*/		ium_mesg_ptr->source=AS_BOG_UNITIG_GRAPH_CC_CM_ID;
		#endif
		/*float32*/		ium_mesg_ptr->coverage_stat=getCovStat();
		/*UnitigStatus*/	ium_mesg_ptr->status=AS_UNASSIGNED;
		/*CDS_COORD_t*/		ium_mesg_ptr->a_branch_point=0;
		/*CDS_COORD_t*/		ium_mesg_ptr->b_branch_point=0;
		/*CDS_COORD_t*/		ium_mesg_ptr->length=_length;
		/*char* */		ium_mesg_ptr->consensus="";
		/*char* */		ium_mesg_ptr->quality="";
		/*int32*/		ium_mesg_ptr->forced=0;
		/*int32*/		ium_mesg_ptr->num_frags=getNumFrags();
		/*IntMultiPos* */	ium_mesg_ptr->f_list=imp_msg_arr;
		/*int32*/		ium_mesg_ptr->num_vars=0;
		/*IntMultiVar* */	ium_mesg_ptr->v_list=NULL;

		return(ium_mesg_ptr);
	}

	//////////////////////////////////////////////////////////////////////////////

	void UnitigGraph::writeIUMtoFile(char *filename){
		
		// Open up the output file
		FILE *fptr=fopen(filename, "w");
		if(fptr==NULL){
			std::cerr << "Could not open " << filename << 
				"for writing IUM (IntUnitigMesg) output." << std::endl;
			
			exit(-1);
		}

		// Step through all the unitigs
		UnitigVector::iterator utg_itr;
		for(utg_itr=unitigs.begin(); utg_itr!=unitigs.end(); utg_itr++){

			IntUnitigMesg *ium_mesg_ptr;
			ium_mesg_ptr=(*utg_itr)->getIUM_Mesg();

			GenericMesg mesg;
			mesg.m=ium_mesg_ptr;
			mesg.t=MESG_IUM;

			WriteProtoMesg_AS(fptr, &mesg);
			(*utg_itr)->freeIUM_Mesg(ium_mesg_ptr);
		}

		fclose(fptr);
	}

	//////////////////////////////////////////////////////////////////////////////


} //AS_BOG namespace



















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
 * $Id: AS_BOG_UnitigGraph.cc,v 1.14 2006-03-22 16:39:26 eliv Exp $
 * $Revision: 1.14 $
*/

//static char AS_BOG_UNITIG_GRAPH_CC_CM_ID[] = "$Id: AS_BOG_UnitigGraph.cc,v 1.14 2006-03-22 16:39:26 eliv Exp $";
static char AS_BOG_UNITIG_GRAPH_CC_CM_ID[] = "gen> @@ [0,0]";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include <float.h>
#include <stdlib.h>
#include <set>

extern "C" {
	#include "AS_global.h"
}

namespace AS_BOG{

	//////////////////////////////////////////////////////////////////////////////

    UnitigGraph::~UnitigGraph() {
		UnitigVector::iterator utg_itr;
		for(utg_itr=unitigs.begin(); utg_itr!=unitigs.end(); utg_itr++)
            delete *utg_itr;
    }

	void UnitigGraph::build(ChunkGraph *cg_ptr, BestOverlapGraph *bog_ptr, long num_rand_frags, long genome_size){

		iuid num_frags=cg_ptr->getNumFragments();
		iuid unitig_id=1;

		iuid frag_idx;
		iuid fp_dst_frag_id, tp_dst_frag_id;

		// Initialize where we've been to nowhere; "Do not retraverse list"
		std::map<iuid, bool> visited_map;
		visited_map.clear();

		BestContainmentMap *best_cntr = &(bog_ptr->_best_containments);
		ContainerMap *cntnrmap_ptr = _build_container_map(best_cntr);

		// Step through all the fragments 
		std::cerr << "Building Unitigs from " << num_frags << " fragments.\n"; 
		for(frag_idx=1; frag_idx<=num_frags; frag_idx++){

			//std::cerr << "Working on " << frag_idx << std::endl; 

			// Check the map to so we don't visit a unitig twice (once from
			//   both ends)
			if(visited_map.find(frag_idx) == visited_map.end() && 
				 best_cntr->find(frag_idx) == best_cntr->end() )
            { 
				cg_ptr->getChunking(
					frag_idx, 
					fp_dst_frag_id, 
					tp_dst_frag_id);

				if(  // If the 5' end is NULL but 3' end isn't, or vice versa
				    (fp_dst_frag_id == NULL_FRAG_ID &&
				     tp_dst_frag_id != NULL_FRAG_ID ) ||
				    (tp_dst_frag_id == NULL_FRAG_ID &&
				     fp_dst_frag_id != NULL_FRAG_ID ) ||
				    (tp_dst_frag_id == NULL_FRAG_ID &&
				     fp_dst_frag_id == NULL_FRAG_ID ) 
				){

					if(!(unitig_id%10000)){
						std::cerr << ".";
					}

					// Allocated a new unitig node
					Unitig *utg=new Unitig;

					// Store the dovetails
					if(fp_dst_frag_id == NULL_FRAG_ID && tp_dst_frag_id != NULL_FRAG_ID){
						utg->dovetail_path_ptr=
							_extract_dovetail_path(
								frag_idx, THREE_PRIME, cg_ptr, bog_ptr);
					}else if(tp_dst_frag_id == NULL_FRAG_ID && fp_dst_frag_id != NULL_FRAG_ID){
						utg->dovetail_path_ptr=
							_extract_dovetail_path(
								frag_idx, FIVE_PRIME, cg_ptr, bog_ptr);
					}else{
						utg->dovetail_path_ptr=
							_extract_dovetail_path(
								frag_idx, THREE_PRIME, cg_ptr, bog_ptr);
					}

					// Get the fragment ID of the last fragment in the
					//   dovetail, then store it in the "do not retraverse
					//   list"
					iuid last_frag_id_in_dovetail= utg->dovetail_path_ptr->back().ident;
					visited_map[last_frag_id_in_dovetail]=true;

					//std::cerr << "Containment Extracted. " << std::endl;

                    utg->computeFragmentPositions(cntnrmap_ptr, best_cntr);

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
        delete cntnrmap_ptr;

		std::cerr << "Setting Global Arrival Rate.\n";
		Unitig static_proxy;
		static_proxy.setGlobalArrivalRate(getGlobalArrivalRate(num_rand_frags, genome_size));

		std::cerr << std::endl << "There were " << unitigs.size() << " unitigs generated.\n";
	}

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	DoveTailPath *UnitigGraph::_extract_dovetail_path(
		iuid src_frag_id, fragment_end_type firstEnd, ChunkGraph *cg_ptr, BestOverlapGraph *bog_ptr){

	// Note:  I only need BestOverlapGraph for it's frag_len and olap_length

		DoveTailPath *dtp_ptr=new DoveTailPath;

		iuid fp_frag_id, tp_frag_id;
		iuid current_frag_id=src_frag_id;
		iuid next_frag_id;
		fragment_end_type whichEnd = firstEnd;
        int frag_end,frag_begin;
        frag_begin = 0;

		//std::cerr<<"Working on: "<<src_frag_id<< " Dir: " << travel_dir <<std::endl;
		iuid last_frag_id;
		while(current_frag_id != NULL_FRAG_ID) {
			// Store the current fragment into dovetail path
			BestEdgeOverlap* bestEdge = bog_ptr->getBestEdgeOverlap(
				current_frag_id, whichEnd
                );
            DoveTailNode dt_node;
            dt_node.type         = AS_READ;
            dt_node.ident        = current_frag_id;
            dt_node.sourceInt    = -1;
            dt_node.contained    = 0;
            dt_node.delta_length = 0;
            dt_node.delta        = NULL;
            next_frag_id         = bestEdge->frag_b_id;
#ifdef NEW_UNITIGGER_INTERFACE
            dt_node.ident2       = next_frag_id;
            // consensus wants positive hangs, so swap
            if (bestEdge->ahang < 0 && bestEdge->bhang < 0 ) {
                dt_node.ahang = -bestEdge->bhang;
                dt_node.bhang = -bestEdge->ahang;
            } else {
                dt_node.ahang = bestEdge->ahang;
                dt_node.bhang = bestEdge->bhang;
            }
#endif
            frag_end = frag_begin + BestOverlapGraph::fragLen(current_frag_id);

			if(whichEnd == FIVE_PRIME){
				dt_node.position.bgn = frag_end;
				dt_node.position.end = frag_begin;
			}else {
				dt_node.position.bgn = frag_begin;
				dt_node.position.end = frag_end;
			}
			dtp_ptr->push_back(dt_node);

            int chunkNextId = cg_ptr->getChunking(current_frag_id, whichEnd); 
            if ( chunkNextId != NULL_FRAG_ID )
                assert( chunkNextId == next_frag_id );

			// Prep the start position of the next fragment
			frag_begin = frag_end - BestOverlapGraph::olapLength( current_frag_id,
                    next_frag_id, bestEdge->ahang, bestEdge->bhang );

			// Set current to next
			current_frag_id = chunkNextId;
            if ( bestEdge->bend == FIVE_PRIME ) {
                whichEnd = THREE_PRIME;
            } else {
                whichEnd = FIVE_PRIME;
            }
        }
		return(dtp_ptr);
	}

	//////////////////////////////////////////////////////////////////////////////

	ContainerMap *UnitigGraph::_build_container_map(BestContainmentMap *best_cntr){

		ContainerMap *cmptr=new ContainerMap;
		int containees=0;

		BestContainmentMap::const_iterator bstcnmap_itr;
		for( bstcnmap_itr  = best_cntr->begin(); 
		     bstcnmap_itr != best_cntr->end(); bstcnmap_itr++)
        {
			iuid ctnee_id     = bstcnmap_itr->first;
			iuid container_id = bstcnmap_itr->second.container;

			(*cmptr)[container_id].push_back(ctnee_id);
			containees++;
		}
		
		std::cerr << "BestContainments Inverted into ContainerMap." << std::endl;
		std::cerr << "Num containees: " << containees << std::endl;
		std::cerr << "Num containers: " << cmptr->size() << std::endl;
		
		return(cmptr);
	}

	//////////////////////////////////////////////////////////////////////////////

	std::ostream& operator << (std::ostream& os, Unitig& utg){
		
		os << "Dovetails:" << std::endl;

		DoveTailPath::const_iterator dt_itr;
		for(dt_itr=utg.dovetail_path_ptr->begin();
			dt_itr!=utg.dovetail_path_ptr->end(); 
			dt_itr++){
			
			os << "  " << dt_itr->ident << std::endl;
		}

		os << "Containments:" << std::endl;
	
		os << "avgRho: " << utg.getAvgRho() << std::endl;
		os << "covStat: " << utg.getCovStat() << std::endl;
		os << "numFrags: " << utg.getNumFrags() << std::endl;
		os << "numRandomFrags: " << utg.getNumRandomFrags() << std::endl;
		os << "length: " << utg.getLength() << std::endl;

		return(os);

	}

	//////////////////////////////////////////////////////////////////////////////

	std::ostream& operator << (std::ostream& os, UnitigGraph& utgrph){
		
		UnitigVector::const_iterator utg_itr;
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
			UnitigVector::const_iterator iter;
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
        _avgRho = -1;
		dovetail_path_ptr=NULL;
	}

	//////////////////////////////////////////////////////////////////////////////

	Unitig::~Unitig(void){
		if(dovetail_path_ptr!=NULL) delete dovetail_path_ptr;
	}

	//////////////////////////////////////////////////////////////////////////////

	float Unitig::getAvgRho(void){

		if(_avgRho!=-1){
			return(_avgRho);
		}
		
		// We will compute the average rho.
		//
		// Since rho is the length(unitig) - length(last fragment),
		//   and the length(last fragment) is ambiguous depending on which
		//   direction we are walking the unitig from.  We will take the average 
		//   of the rhos through both directions.

		DoveTailPath::const_iterator dtp_iter;

		// Get first fragment's length
		dtp_iter=dovetail_path_ptr->begin();
		long first_frag_len = BestOverlapGraph::fragLen(dtp_iter->ident);

		// Get last fragment's length
		dtp_iter=dovetail_path_ptr->end();
		dtp_iter--;
		long last_frag_len = BestOverlapGraph::fragLen(dtp_iter->ident);

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

		if(dovetail_path_ptr->size()==0){
			std::cerr << "This Unitig has an empty fragPositions." << std::endl;	
		}

        DoveTailPath::const_iterator fpm_itr;
		for(
		    fpm_itr=dovetail_path_ptr->begin();
		    fpm_itr!=dovetail_path_ptr->end();
		    fpm_itr++){

			SeqInterval intrvl=fpm_itr->position;

			if(max_pos<intrvl.bgn){
				max_pos=intrvl.bgn;
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

//		if(_numFrags!=-1){
//			return(_numFrags);
//		}
        if (!dovetail_path_ptr->empty()) {
            _numFrags = dovetail_path_ptr->size();
            return _numFrags;
        }

		long num_dovetailing_frags = dovetail_path_ptr->size();
		long num_contained_frags=0;

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
    //Recursively place all contains under this one into the FragmentPositionMap

    void Unitig::placeContains( const ContainerMap* cntnrp, BestContainmentMap *bestCtn,
                                const iuid container, const SeqInterval intvl)
    {
            ContainerMap::const_iterator ctmp_itr = cntnrp->find( container );
            if (ctmp_itr != cntnrp->end() ) {

                ContaineeList::const_iterator cntee_itr;
                for(cntee_itr  = ctmp_itr->second.begin();
                    cntee_itr != ctmp_itr->second.end(); cntee_itr++)
                {
                    iuid cntee = *cntee_itr;
                    BestContainment &best = (*bestCtn)[ cntee ];

                    if (best.isPlaced)
                        continue;

                    assert( best.container == container );

                    int offset = best.a_hang;
                    int bhang  = best.b_hang;
                    DoveTailNode pos;
                    pos.type         = AS_READ;
                    pos.ident        = cntee;
                    pos.sourceInt    = -1;
                    pos.contained    = container;
                    pos.delta_length = 0;
                    pos.delta        = NULL;
#ifdef NEW_UNITIGGER_INTERFACE
                    pos.ident2       = container;
                    pos.ahang        = offset;
                    pos.bhang        = bhang;
#endif

                    if(intvl.bgn < intvl.end) {
                        pos.position.bgn = intvl.bgn + offset;
                        pos.position.end = intvl.end + bhang;

                    } else if (intvl.bgn > intvl.end) {
                        pos.position.bgn = intvl.bgn - offset;
                        pos.position.end = intvl.end - bhang;

                    }else{
                        std::cerr << "Error, container size is zero." << std::endl;
                        assert(0);
                    }
                    // Swap ends if containee is not same strand as container
                    if(!best.sameOrientation){
                        int tmp          = pos.position.bgn;
                        pos.position.bgn = pos.position.end;
                        pos.position.end = tmp;
                    }

                    dovetail_path_ptr->push_back( pos );
                    best.isPlaced = true;
                    placeContains( cntnrp, bestCtn, cntee, pos.position);
                }
            }
    }

    void Unitig::computeFragmentPositions(ContainerMap *allcntnr_ptr,
                 BestContainmentMap *bestContain)
    {
		long frag_ins_begin = 0;
		long frag_ins_end;
        iuid lastFrag = 0;
#ifdef NEW_UNITIGGER_INTERFACE
        iuid nextFrag = 0;
#endif
		//std::cerr << "Positioning dovetails." << std::endl;
		// place dovetails in a row
        long numDoveTail = dovetail_path_ptr->size();
		for(long i=0; i < numDoveTail-1; i++) {
            DoveTailNode *dt_itr = &(*dovetail_path_ptr)[i];

            iuid fragId  = dt_itr->ident;
#ifdef NEW_UNITIGGER_INTERFACE
            if ( nextFrag != 0 )
                assert( nextFrag == fragId);
            nextFrag = dt_itr->ident2;
#endif
            lastFrag = fragId;

            placeContains( allcntnr_ptr, bestContain, fragId, dt_itr->position);
		}
        // Compute assuming that containee is the same orientation as container
        //	if(cntnr_intvl.begin < cntnr_intvl.end)

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

        // else if(cntnr_intvl.begin > cntnr_intvl.end)

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

	}

	//////////////////////////////////////////////////////////////////////////////

	void Unitig::freeIUM_Mesg(IntUnitigMesg *ium_ptr){
		delete ium_ptr;
	}

	//////////////////////////////////////////////////////////////////////////////

	int IntMultiPosCmp(const void *a, const void *b){
		IntMultiPos *impa=(IntMultiPos*)a;
		IntMultiPos *impb=(IntMultiPos*)b;
        long aleft,aright,bleft,bright;
		if (impa->position.bgn < impa->position.end) {
			aleft  = impa->position.bgn;
            aright = impa->position.end;
        } else {
			aright = impa->position.bgn;
            aleft  = impa->position.end;
        }
		if (impb->position.bgn < impb->position.end) {
			bleft  = impb->position.bgn;
            bright = impb->position.end;
        } else {
			bright = impb->position.bgn;
            bleft  = impb->position.end;
        }
		if(aleft!=bleft)
        {
			return(aleft - bleft);
		}
        else if (aright != bright)
        {
            return(bright - aright);
        }
        else {
			if(impa->contained == impb->ident)
                return(1);
			if(impb->contained == impa->ident)
                return(-1);
			if(impa->contained!=0)
				return(1);
			if(impb->contained!=0)
				return(-1);
			return(0);
		}
	}

	IntUnitigMesg *Unitig::getIUM_Mesg(){
		
		IntUnitigMesg *ium_mesg_ptr=new IntUnitigMesg;
		IntMultiPos *imp_msg_arr   = &(dovetail_path_ptr->front());

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

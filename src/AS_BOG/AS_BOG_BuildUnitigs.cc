
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
* Module:  AS_BOG_BuildUnitigs.cc
* Description:
*	Class definition for build unitigs algorithm
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
 * $Id: AS_BOG_BuildUnitigs.cc,v 1.1 2005-08-08 21:51:03 kli1000 Exp $
 * $Revision: 1.1 $
*/

static char CM_ID[] = "$Id: AS_BOG_BuildUnitigs.cc,v 1.1 2005-08-08 21:51:03 kli1000 Exp $";

//  System include files
#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_UnitigGraph.hh"

namespace AS_BOG{

	
	void BuildUnitigs::visit(ChunkGraph &cg){

		iuid num_frags=cg->getNumFragments();

		_reported_array = new bool[num_frags];

		// Initialize 
		iuid frag_idx;
		for(frag_idx=0; frag_idx<num_frags; frag_idx++){
			_reported_array[frag_idx]=false;
		}

		iuid fp_dst_frag_id, tp_dst_frag_id;
		iuid unitig_id=0;
		for(frag_idx=0; frag_idx<num_frags; frag_idx++){

			if(!_reported_array[frag_idx]){

				cg.getChunking(
					src_frag_id, 
					fp_dst_frag_id, 
					tp_dst_frag_id);
				
				// Only walk the fragment if either end is empty
				if(
				    fp_dst_frag_id == NULL_FRAG_ID ||
				    tp_dst_frag_id == NULL_FRAG_ID )
				{

					// For the given fragment end, follow the chunk
					DoveTailPath *dtp_ptr=
						_extract_unitig(src_frag_id, cg);

					// Store dove tail path in unitig graph
					Unitig *unitig_ptr=new Unitig(unitig_id);
					unitig_ptr->setDoveTailPath(dtp_ptr);
					_unitig_graph.AddUnitig(unitig_ptr);

					// Make sure we don't rewalk the unitig path
					_reported_array[last_fragment_id]=true;

					// Make sure unitig_id are unique
					unitig_id++;
				}

			}
		}
		delete[] _reported_array;
		
	}

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	DoveTailPath *BuildUnitigs::_extract_unitig(iuid src_frag_id, ChunkGraph& cg){

		DoveTailPath *dtp_ptr=new DoveTailPath;

		iuid fp_frag_id, tp_frag_id;
		cg.getChunking(src_frag_id, fp_frag_id, tp_frag_id);

		orientation_type trav_dir;

		// Decide which direction to go
		if(fp_frag_id == NULL_FRAG_ID){
			// Go three prime
			trav_dir=THREE_PRIME;
			src_frag_id=tp_frag_id;
		}else if(tp_frag_id == NULL_FRAG_ID){
			// Go five prime
			trav_dir=FIVE_PRIME;
			src_frag_id=fp_frag_id;
		}


		// Start walking the chunk path
		iuid next_frag_id;
		orientation_type next_trav_dir;

		do{
			dtp_ptr->SetNextFragment(src_frag_id, trav_dir);

			if(src_frag_id != NULL_FRAG_ID){
				_find_next_fragment(
					src_frag_id, trav_dir,
					next_frag_id, next_trav_dir,
					cg);
			}


			src_frag_id=next_frag_id;
			trav_dir=next_trav_dir;

		}while(src_frag_id != NULL_FRAG_ID);

		return(dtp_ptr);
			
	}

	//////////////////////////////////////////////////////////////////////////////

	void BuildUnitigs::_find_next_fragment(
	    iuid src_frag_id, orientation_type src_end,
	    iuid& dst_frag_id, orientation_type& dst_end,
	    ChunkGraph& cg
	){

		iuid dst_fp_frag_id, dst_tp_frag_id;	//Wht the source points at
		iuid rflx_fp_frag_id, rflx_tp_frag_id; 	//What the destination points at

		cg.getChunking(src_frag_id, dst_fp_frag_id, dst_tp_frag_id);

		// Follow the next node based on the suggested direction
		if(which_end == FIVE_PRIME){
			cg.getChunking(dst_fp_frag_id, rflx_fp_frag_id, rflx_tp_frag_id);
			dst_frag_id=dst_fp_frag_id;
		}else if(which_end == THREE_PRIME){
			cg.getChunking(dst_tp_frag_id, rflx_fp_frag_id, rflx_tp_frag_id);
			dst_frag_id=dst_tp_frag_id;
		}

		// If the destination fragment's 3' end matches the source fragment,
		//  then we want to go the opposite direction, or else we will be goin
		//  backwards!

		if(src_frag_id==rflx_tp_frag_id){
			dst_end=FIVE_PRIME;
		}else{
			dst_end=THREE_PRIME;
		}
	}


} //AS_BOG namespace



















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
* Module:  AS_BOG_BuildUnitigs.hh
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
 * $Id: AS_BOG_BuildUnitigs.hh,v 1.2 2005-08-08 21:49:02 kli1000 Exp $
 * $Revision: 1.2 $
*/

static char CM_ID[] = "$Id: AS_BOG_BuildUnitigs.hh,v 1.2 2005-08-08 21:49:02 kli1000 Exp $";

//  System include files

#ifndef INCLUDE_AS_BOG_BUILDUNITIGS
#define INCLUDE_AS_BOG_BUILDUNITIGS

#include <stack>
#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_UnitigGraph.hh"

namespace AS_BOG{

	class BuildUnitigs : public ChunkGraphVisitor{

	    public:

		void visit(ChunkGraph &bovlg);

	    private:

		DoveTailPath *_extract_unitig(iuid src_frag_id, ChunkGraph& cg);


		void _find_next_fragment(
			iuid src_frag_id, orientation_type src_end,
			iuid& dst_frag_id, orientation_type& dst_end,
			ChunkGraph& cg)

		bool _reported_array[]=NULL;
		UnitigGraph _unitig_graph;
	}

} //AS_BOG namespace


#endif


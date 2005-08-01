
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
* Module:  AS_BOG_Datatypes.c
* Description:
*	Common datatypes
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
 * $Id: AS_BOG_Datatypes.hh,v 1.2 2005-08-01 15:19:38 kli1000 Exp $
 * $Revision: 1.2 $
*/

static char CM_ID[] = "$Id: AS_BOG_Datatypes.hh,v 1.2 2005-08-01 15:19:38 kli1000 Exp $";

//  System include files

#ifndef INCLUDE_AS_BOG_DATATYPES
#define INCLUDE_AS_BOG_DATATYPES

namespace AS_BOG{

	typedef enum { 
		// Dovetailing overlaps
		DOVE_NORMAL,		// AB_AB
		DOVE_ANTI_NORMAL,	// BA_BA
		DOVE_INNIE,		// AB_BA
		DOVE_OUTTIE,		// BA_AB

		// Containment overlaps
		CONT_A_CONTAINS,
		CONT_B_CONTAINS,
		CONT_MUTUAL
	} overlap_type;

	typedef enum {
		FIVE_PRIME, 	// 5' End of fragment
		THREE_PRIME 	// 3' End of Fragment
	} fragment_end_type;

	typedef unsigned int iuid;

} //AS_BOG namespace


#endif


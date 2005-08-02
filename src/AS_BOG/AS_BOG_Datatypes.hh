
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
 * $Id: AS_BOG_Datatypes.hh,v 1.3 2005-08-02 15:49:47 kli1000 Exp $
 * $Revision: 1.3 $
*/

static char CM_ID[] = "$Id: AS_BOG_Datatypes.hh,v 1.3 2005-08-02 15:49:47 kli1000 Exp $";

//  System include files

#ifndef INCLUDE_AS_BOG_DATATYPES
#define INCLUDE_AS_BOG_DATATYPES

namespace AS_BOG{

	typedef enum { 
		UNDEFINED,

		// Dovetailing overlaps
		DOVE_NORMAL,		// AB_AB
			// A ----->
			// B    ----->
		DOVE_INNIE,		// AB_BA
			// A ----->
			// B    <-----
		DOVE_OUTTIE,		// BA_AB
			// A    ----->
			// B <-----
		DOVE_ANTI_NORMAL,	// BA_BA (should be same as Normal if B is rc'd)
			// A <-----
			// B    <-----

		// Containment overlaps
		CONT_NORMAL,
			// A -------->
			// B   --->
		CONT_INNIE,
			// A -------->
			// B   <---
		CONT_OUTTIE,
			// A <--------
			// B   --->
		CONT_ANTI_NORMAL,
			// A <--------
			// B   <---
		CONT_MUTUAL,
			// A -------->
			// B -------->
		CONT_MUTUAL_INNIE,
			// A -------->
			// B <--------
	} overlap_type;

	typedef enum {
		FIVE_PRIME, 	// 5' End of fragment
		THREE_PRIME 	// 3' End of Fragment
	} fragment_end_type;

	typedef unsigned int iuid;

} //AS_BOG namespace


#endif


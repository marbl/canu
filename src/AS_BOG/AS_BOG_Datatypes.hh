
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
 * $Id: AS_BOG_Datatypes.hh,v 1.5 2005-08-08 21:49:02 kli1000 Exp $
 * $Revision: 1.5 $
*/

static char CM_ID[] = "$Id: AS_BOG_Datatypes.hh,v 1.5 2005-08-08 21:49:02 kli1000 Exp $";

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
			// 	A <-----
			// 	B    <-----
			// AKA 
			// 	A    ----->
			// 	B  ----->

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
	} overlap_type;

	typedef enum {
		FIVE_PRIME, 	// 5' End of fragment
		THREE_PRIME 	// 3' End of Fragment
	} fragment_end_type;

	///////////////////////////////////////////////////////////////////////

	fragment_end_type opposite_end(fragment_end_type fet){
		if(fet == FIVE_PRIME){
			return(THREE_PRIME);
		}else if(fet == THREE_PRIME){
			return(FIVE_PRIME);
		}else{
			//Assert
		}
	}

	fragment_end_type overlap_end(overlap_type ot, char a_or_b){

		switch(a_or_b){
			case 'a':
			case 'A':
				switch(ot){
					case DOVE_NORMAL:
					case DOVE_INNIE:
						which_end=THREE_PRIME;
						break;
					case DOVE_ANTI_NORMAL:
					case DOVE_OUTTIE:
						which_end=FIVE_PRIME;
						break;
				}
				return(which_end);
			case 'b':
			case 'B':
				switch(ot){
					case DOVE_NORMAL:
					case DOVE_OUTTIE:
						which_end=FIVE_PRIME;
						break;
					case DOVE_INNIE:
					case DOVE_ANTI_NORMAL:
						which_end=THREE_PRIME;
						break;
				}
				return(which_end);
		}
	}

	///////////////////////////////////////////////////////////////////////
	

	typedef enum {
		UNKNOWN,
		FORWARD,
		REVERSE
	} orientation_type;

	typedef unint32 iuid;

	const iuid NULL_FRAG_ID=-1;

} //AS_BOG namespace


#endif


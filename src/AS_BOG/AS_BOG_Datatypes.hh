
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
 * $Id: AS_BOG_Datatypes.hh,v 1.16 2007-06-11 20:59:42 eliv Exp $
 * $Revision: 1.16 $
*/

#ifndef INCLUDE_AS_BOG_DATATYPES
#define INCLUDE_AS_BOG_DATATYPES

//  System include files

#include <map>
#include <list>
#include <vector>

extern "C" {
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
}

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

	enum fragment_end_type {
		FIVE_PRIME, 	// 5' End of fragment
		THREE_PRIME 	// 3' End of Fragment
	};

	typedef enum {
		UNKNOWN,
		FORWARD,
		REVERSE
	} orientation_type;

	typedef CDS_IID_t iuid;
	const iuid NULL_FRAG_ID=0;

    typedef std::list<SeqInterval> IntervalList;

    struct FragmentEnd {
        iuid id;
        fragment_end_type end;

        FragmentEnd(iuid id=0, fragment_end_type end=FIVE_PRIME) :
            id(id), end(end) {}

    };
    inline bool operator==(FragmentEnd a, FragmentEnd b) {
        if (a.id == b.id && a.end == b.end)
            return true;
        else
            return false;
    };
    inline bool operator<(FragmentEnd a, FragmentEnd b) {
        if (a.id != b.id)
            return a.id < b.id;
        else
            return a.end < b.end;
    };
} //AS_BOG namespace


#endif


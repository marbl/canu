/**************************************************************************
 * This file is part of A2Amapper.
 * Copyright (c) 2004 Applera Corporation
 * Author: Clark Mobarry
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
**************************************************************************/


/*********************************************************************
 *		  Confidential -- Do Not Distribute                  *
 *	    Copyright (c) 2003 Celera Genomics Corporation.          *
 *			 All Rights Reserved.                        *
 *********************************************************************/

#ifndef SUBDELREZ_H
#define SUBDELREZ_H

#include <assert.h>

// Define which SUB or DEL you want to have.
// This file is included by AS_ALN_aligners.h
#define UNIT
//#define QUALITY

#ifdef UNIT
#define DEL(a,qa)       1
#define SUB(a,qa,b,qb) ((a)!=(b))
#endif

//#ifdef QUALITY
//#define DEL(a,qa)  MismatchTableREZ[(qa-'0'-1)*60+(qa-'0'-1)]
//#define SUB(a,qa,b,qb) ( (a)==(b) ? MatchTableREZ[(qa-'0'-1)*60+(qb-'0'-1)] : MismatchTableREZ[(qa-'0'-1)*60+(qb-'0'-1)] )
//#endif

#endif









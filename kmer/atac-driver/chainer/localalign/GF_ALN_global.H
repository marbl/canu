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

/* This is the global include file that all C files in the AS subsystem should
   include.
*/
#ifndef GF_GLOBAL_H
#define GF_GLOBAL_H

#include "GF_ALN_cds.h"
#include "GF_ALN_pmesg.h"

// Constants that SHOULD be included

#ifndef max
#define max(a,b)		( ((a) > (b)) ? (a) : (b) )
#endif
#ifndef min
#define min(a,b)		( ((a) < (b)) ? (a) : (b) )
#endif

// cgw and cns use NULLINDEX for a NULL index value
#define NULLINDEX (-1)

// A convenient assert for testing whether ptrs are null
// without bothering lint
#define AssertPtr(ptr) (assert((ptr) != NULL))


#endif

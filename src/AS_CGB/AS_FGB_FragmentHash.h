
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
/*********************************************************************
 * $Id: AS_FGB_FragmentHash.h,v 1.2 2004-09-23 20:25:01 mcschatz Exp $
 *
 * Module:
 * Description:
 * Assumptions: Too many to count.
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_FRAGMENT_HASH_INCLUDE
#define AS_CGB_FRAGMENT_HASH_INCLUDE

#include "AS_MSG_pmesg.h"
#if 0
#include "AS_UTL_Var.h"
#endif

#define AS_CGB_NOT_SEEN_YET  CDS_INT32_MAX

//typedef IntFragment_ID FragmentHashInt;

typedef void FragmentHashObject;
// Hide the implementation of the object.  I encapsulate the object
// (1) data with a void pointer and (2) methods by extern functions.

extern FragmentHashObject * create_FragmentHash(IntFragment_ID max_iid);
// Create an object. The return value is non-NULL when successful.
// The maxiid is a pre-allocation amount for the maximum number of
// fragments.  You can specify zero for the pre-allocation amounts.

extern int destroy_FragmentHash(FragmentHashObject *);
// Destroy an object. The return value is zero when successful.

extern void set_vid_FragmentHash
(FragmentHashObject *self, IntFragment_ID iid, IntFragment_ID vid);

extern IntFragment_ID get_vid_FragmentHash
(FragmentHashObject *self, IntFragment_ID iid);

#endif // AS_CGB_FRAGMENT_HASH_INCLUDE


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
#ifndef AS_UTL_HASHCOMMON_H
#define AS_UTL_HASHCOMMON_H

#include "AS_global.h"

/* This file declares common items used by both
   AS_UTL_Hash (general open hashing) and
   AS_UTL_PHash (open persistent hashing of the UID->IID mapping 
*/
#define hashsize(n) ((uint32)1<<(n))
#define hashmask(n) (hashsize(n)-1)

#define HASH_FAILURE 0
#define HASH_SUCCESS 1

#include "math_AS.h"

uint32 Hash_AS( register uint8 *k,        /* the key */
	  register uint32  length,   /* the length of the key */
	  register uint32  initval);   /* the previous hash, or an arbitrary value */

#endif

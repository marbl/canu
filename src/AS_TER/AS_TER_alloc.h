
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
   CVS_ID:  $Id: AS_TER_alloc.h,v 1.1.1.1 2004-04-14 13:53:43 catmandew Exp $
 *********************************************************************/
#ifndef AS_TER_ALLOC_H
#define AS_TER_ALLOC_H

#include "AS_global.h"
#include "AS_TER_utils.h"
#include "AS_MSG_pmesg.h"

#define UID_CODE_OK                     101

/*****************************************************************/
/* UID ALLOCATION */
/*****************************************************************/

int32 get_uids(uint64 blockSize, uint64 *interval, int32 real);
/*****************************************************************/
// Allocates blockSize many UIDs from the UID server if real
// is TRUE. Otherwise it allocates some dummy numbers.
/*****************************************************************/

int32 get_next_uid(uint64 *uid, int32 real);
/*****************************************************************/
// Returns the next available uid. (A real if real==TRUE, o.w. a fake UID
/*****************************************************************/

void check_environment(void);
/*****************************************************************/
// Checks whether the SYS_UID* variables are set.
// If not it exits proposing a standard value
/*****************************************************************/


void set_start_uid(uint64 s);

#endif










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
/* $Id: AS_CVT_uid_interface.h,v 1.1.1.1 2004-04-14 13:51:28 catmandew Exp $ */
#ifndef AS_CVT_UID_INTERFACE_H
#define AS_CVT_UID_INTERFACE_H

/*********************************************************************/
// headers
/*********************************************************************/
// project headers
#include "AS_global.h"

// structure to hold variables needed in intracting with UID server
typedef struct
{
  cds_uint64  currUID;
} UIDInteractor;
typedef UIDInteractor * UIDInteractorp;

// function that sets up UID server interaction
UIDInteractorp CreateUIDInteractor( cds_uint64 block_size );

void DestroyUIDInteractor( UIDInteractorp ui );

// function to get a UID
int GetUID( UIDInteractorp ui, cds_uint64 * uid );

#endif // AS_CVT_UID_INTERFACE_H

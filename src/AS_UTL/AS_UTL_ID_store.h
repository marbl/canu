
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
/* $Id: AS_UTL_ID_store.h,v 1.2 2004-09-23 20:25:29 mcschatz Exp $ */

/*
  This is a set of utility functions for managing a list of UIDS.
  Written by Ian Dew for AS_MSG/remove_fragments, moved to 
  AS_UTL for general use by Karin Remington.

*/
#ifndef AS_UTL_IDSTORE_INCLUDE
#define AS_UTL_IDSTORE_INCLUDE
// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// project headers
#include "AS_global.h"

typedef struct
{
  CDS_UID_t * ids;
  cds_uint32   num_ids;
  cds_uint32   size;
  CDS_UID_t   min_id;
  CDS_UID_t   max_id;
} ID_Array;
typedef ID_Array * ID_Arrayp;

// qsort sorter, with some checking - there should be no dupes
int uid_compare( const CDS_UID_t * a, const CDS_UID_t * b );

void FreeID_Array( ID_Arrayp array );

ID_Arrayp AllocateID_Array( cds_uint32 num_ids );

int AppendToID_Array( ID_Arrayp array, CDS_UID_t id, int sort );

cds_int64 FindID_ArrayID( ID_Array * ids, CDS_UID_t id );

#endif

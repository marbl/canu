
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
/* $Id: AS_UTL_ID_store.c,v 1.5 2006-11-14 17:52:18 eliv Exp $ */

/*
  This is a set of utility functions for managing a list of UIDS.
  Written by Ian Dew for AS_MSG/remove_fragments, moved to 
  AS_UTL for general use by Karin Remington.

*/

// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// project headers
#include "AS_global.h"
#include "AS_UTL_ID_store.h"


// qsort sorter, with some checking - there should be no dupes
int uid_compare( const CDS_UID_t * a, const CDS_UID_t * b )
{
  if( *a > *b ) return 1;
  if( *a == *b )
  {
    fprintf( stderr, "Duplicate IDs!: " F_UID "\n", *a );
    return 0;
  }
  return -1;
}


void FreeID_Array( ID_Arrayp array )
{
  if( array != NULL )
  {
    if( array->ids != NULL )
      free( array->ids );
    free( array );
  }
}


ID_Arrayp AllocateID_Array( cds_uint32 num_ids )
{
  ID_Arrayp ret_array;
  ret_array = (ID_Arrayp) calloc( 1, sizeof( ID_Array ) );
  if( ret_array == NULL )
    return NULL;
  if( num_ids != 0 )
  {
    ret_array->ids = (CDS_UID_t *) calloc( num_ids, sizeof( CDS_UID_t ) );
    if( ret_array->ids == NULL )
    {
      FreeID_Array( ret_array );
      return NULL;
    }
    ret_array->size = num_ids;
  }
  return ret_array; 
}

int AppendToID_Array( ID_Arrayp array, CDS_UID_t id, int sort )
{
  if( array->num_ids == array->size )
    return 1;
  
  array->ids[array->num_ids] = id;
  array->num_ids++;
  array->min_id = min( id, array->min_id );
  array->max_id = MAX( id, array->max_id );
  if( array->num_ids > 1 && sort != 0 )
    qsort( array->ids, array->num_ids, sizeof( CDS_UID_t ),
           (int (*)(const void *, const void *)) uid_compare );

  return 0;
}


cds_int64 FindID_ArrayID( ID_Array * ids, CDS_UID_t id )
{
  double delta = ids->max_id - id;
  CDS_UID_t i = (CDS_UID_t) ( (ids->num_ids - 1) *
    (1.0f - delta / (ids->max_id - ids->min_id)) );

  if( id > ids->ids[ids->num_ids - 1] || id < ids->ids[0] )
    return -1;
  
  i = min( ids->num_ids - 1, i );
  
  if( ids->ids[i] > id )
  {
    for( ; i > 0 && ids->ids[i] > id; i-- );
  }
  else
  {
    for( ; i < ids->num_ids && ids->ids[i] < id; i++ );
  }
  
  if( i < ids->num_ids && ids->ids[i] == id )
    return i;
  //else
    return -1;
}

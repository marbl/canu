
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
/* $Id: remove_fragment.c,v 1.2 2004-09-23 20:25:24 mcschatz Exp $ */

/*
  This program, bcp2frg, converts a set of BAC-end batch copy files
  into one DST, and several FRG and LKG messages written to stdout.

  NOTE: this is likely to be of very temporary use. In the (near?) future
        BAC end data should be queryable in a database...
*/

// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// project headers
#include "AS_global.h"

typedef struct
{
  CDS_UID_t * ids;
  cds_uint32   num_ids;
  cds_uint32   size;
  CDS_UID_t min_id;
  CDS_UID_t max_id;
} ID_Array;
typedef ID_Array * ID_Arrayp;


// qsort sorter, with some checking - there should be no dupes
static int uid_compare( const CDS_UID_t * a, const CDS_UID_t * b )
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
  array->max_id = max( id, array->max_id );
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

  return -1;
}

CDS_UID_t a2UID( char * string )
{
  CDS_UID_t ret_val = 0;
  sscanf(string, F_UID, &ret_val);
  return ret_val;
}

// main function
int main( int argc, char ** argv )
{
  FILE       * fp_in;
  FILE       * fp_out;
  MesgReader   reader;
  MesgWriter   writer;
  GenericMesg      * gen;
  ID_Arrayp  frag_uids;
  ID_Arrayp  frag_uids_found;
  ID_Arrayp  frag_iids;
  cds_int64  i;
  
  if( argc < 4 )
  {
    fprintf( stderr,
             "Usage: %s input_file output_file ((-f list_file) | UID_list))\n", argv[0] );
    return 1;
  }

  // open files
  fp_in = fopen( argv[1], "r" );
  if( fp_in == NULL )
  {
    fprintf( stderr, "Failed to open file %s for reading.\n", argv[1] );
    return 1;
  }
  reader = InputFileType_AS( fp_in );

  // set up frag_uids array
  if( argv[3][0] != '-' )
  {
    // read in from command line
    frag_uids = AllocateID_Array( argc - 3 );
    frag_uids_found = AllocateID_Array( argc - 3 );
    frag_iids = AllocateID_Array( argc - 3 );
    if( frag_uids == NULL || frag_uids_found == NULL || frag_iids == NULL )
      return 1;
    for( i = 0; i < argc - 4; i++ )
      AppendToID_Array( frag_uids, a2UID( argv[i+3] ), 0 );
    AppendToID_Array( frag_uids, a2UID( argv[i+3] ), 1 );
  }
  else
  {
    // read in from file
    FILE * fp_list;
    char   string[1000];
    int    num_uids;

    fp_list = fopen( argv[4], "r" );
    if( fp_list == NULL )
    {
      fprintf( stderr, "Failed to open list file %s for reading.\n", argv[4] );
      return 1;
    }

    num_uids = 0;
    while( fgets( string, 1000, fp_list ) )
    {
      num_uids++;
    }
    rewind( fp_list );

    frag_uids = AllocateID_Array( num_uids );
    frag_uids_found = AllocateID_Array( num_uids );
    frag_iids = AllocateID_Array( num_uids );
    if( frag_uids == NULL || frag_uids_found == NULL || frag_iids == NULL )
      return 1;
    for( i = 0; i < num_uids - 1; i++ )
    {
      fgets( string, 1000, fp_list );
      AppendToID_Array( frag_uids, a2UID( string ), 0 );
    }
    fgets( string, 1000, fp_list );
    AppendToID_Array( frag_uids, a2UID( string ), 1 );

    fclose( fp_list );
  }

  // open output file
  fp_out = fopen( argv[2], "w" );
  if( fp_out == NULL )
  {
    fprintf( stderr, "Failed to open file %s for writing.\n", argv[2] );
    return 1;
  }
  writer =
    (reader == ReadBinaryMesg_AS) ? WriteBinaryMesg_AS : WriteProtoMesg_AS;
  
  while( reader( fp_in, &gen ) != EOF )
  {
    switch( gen->t )
    {
      case MESG_FRG:
	{
	  FragMesg  * frg = (FragMesg *) gen->m;
	  if( (i = FindID_ArrayID( frag_uids, frg->eaccession)) > -1 )
	    {
	      AppendToID_Array( frag_uids_found, frg->eaccession, 1 );
	      WriteProtoMesg_AS( stdout, gen );
	    }
	  else
	    {
	      writer( fp_out, gen );
	    }
	}
        break;
      case MESG_IFG:
      case MESG_SFG:
      case MESG_OFG:
	{
	  FragMesg * frg = (FragMesg *) gen->m;
	  if( (i = FindID_ArrayID( frag_uids, frg->eaccession)) > -1 ||
	      (i = FindID_ArrayID( frag_uids,
				   (CDS_UID_t) frg->iaccession)) > -1 )
	    {
	      AppendToID_Array( frag_uids_found, frg->eaccession, 1 );
	      AppendToID_Array( frag_iids, (CDS_UID_t) frg->iaccession, 1 );
	      WriteProtoMesg_AS( stdout, gen );
	    }
	  else
	    {
	      writer( fp_out, gen );
	    }
	}
        break;
      case MESG_LKG:
	{ 
	  LinkMesg * lkg = (LinkMesg *) gen->m;
	  if( frag_uids_found->num_ids == 0 ||
	      (FindID_ArrayID( frag_uids_found, lkg->frag1) == -1 &&
	       FindID_ArrayID( frag_uids_found, lkg->frag2) == -1) )
	    {
	      writer( fp_out, gen );
	    }
	  else
	    {
	      WriteProtoMesg_AS( stdout, gen );
	    }
	}
        break;
      case MESG_ILK:
	{ 
	  InternalLinkMesg * ilk = (InternalLinkMesg *) gen->m;
	  if( frag_iids->num_ids == 0 ||
	      (FindID_ArrayID( frag_iids, (CDS_UID_t) ilk->ifrag1) == -1 &&
	       FindID_ArrayID( frag_iids, (CDS_UID_t) ilk->ifrag2) == -1) )
	    {
	      writer( fp_out, gen );
	    }
	  else
	    {
	      WriteProtoMesg_AS( stdout, gen );
	    }
	}
	break;
      case MESG_OVL:
	{ 
	  OverlapMesg  * ovl = (OverlapMesg *) gen->m;
	  if( frag_iids->num_ids == 0 ||
	      (FindID_ArrayID( frag_iids, (CDS_UID_t) ovl->aifrag) == -1 &&
	       FindID_ArrayID( frag_iids, (CDS_UID_t) ovl->bifrag) == -1) )
	    {
	      writer( fp_out, gen );
	    }
	  else
	    {
	      WriteProtoMesg_AS( stdout, gen );
	    }
	}
        break;
      default:
        writer( fp_out, gen );
        break;
    }
  }

  fclose( fp_in );
  fclose( fp_out );
  return 0;
}

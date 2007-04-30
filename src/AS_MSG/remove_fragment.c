
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
/* $Id: remove_fragment.c,v 1.11 2007-04-30 13:00:30 brianwalenz Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_ID_store.h"

CDS_UID_t a2UID( char * string )
{
  CDS_UID_t ret_val = 0;
  sscanf(string, F_UID, &ret_val);
  return ret_val;
}

// main function
int main( int argc, char ** argv )
{
  FILE       * fp_out;
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
  
  while( ReadProtoMesg_AS(stdin, &gen ) != EOF )
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
	      WriteProtoMesg_AS( fp_out, gen );
	    }
	}
        break;
      case MESG_IFG:
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
	      WriteProtoMesg_AS( fp_out, gen );
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
	      WriteProtoMesg_AS( fp_out, gen );
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
	      WriteProtoMesg_AS( fp_out, gen );
	    }
	  else
	    {
	      WriteProtoMesg_AS( stdout, gen );
	    }
	}
        break;
      default:
        WriteProtoMesg_AS( fp_out, gen );
        break;
    }
  }

  fclose( fp_out );
  return 0;
}

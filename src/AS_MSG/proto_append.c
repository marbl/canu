
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
/* $Id: proto_append.c,v 1.4 2005-03-22 19:49:01 jason_miller Exp $ */

// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// project headers
#include "AS_global.h"

int main( int argc, char ** argv )
{
  FILE        * fp_from;
  FILE        * fp_to;
  GenericMesg * gen;
  char          code;

  // check parameters: need to-file & from-file names
  if( argc != 3 )
  {
    fprintf( stderr, "Usage: %s to_file from_file\n", argv[0] );
    return 1;
  }

  // check sanity - that from & to files are not the same
  if( strcmp( argv[1], argv[2] ) == 0 )
  {
    fprintf( stderr, "Error: Cannot append file to itself!\n" );
    return 1;
  }
  
  // open files
  fp_to = fopen( argv[1], "a+" );
  if( fp_to == NULL )
  {
    fprintf( stderr, "Failed to open file %s for appending.\n", argv[1] );
    return 1;
  }
  fp_from = fopen( argv[2], "r" );
  if( fp_from == NULL )
  {
    fprintf( stderr, "Failed to open file %s for reading.\n", argv[2] );
    return 1;
  }

  // check that both files are ASCII
  rewind( fp_to );
  fread( &code, sizeof( code ), 1, fp_to );
  if( code != '{' )
  {
    fprintf( stderr, "Error: %s is not ASCII protoIO file.\n", argv[1] );
    return 1;
  }
  rewind( fp_from );
  fread( &code, sizeof( code ), 1, fp_from );
  if( code != '{' )
  {
    fprintf( stderr, "Error: %s is not ASCII protoIO file.\n", argv[2] );
    return 1;
  }

  // move to end of to file & beginning of from file
  CDS_FSEEK( fp_to, (off_t) 0, SEEK_END );
  rewind( fp_from );

  // append from file to to file
  while( ReadProtoMesg_AS( fp_from, &gen ) != EOF )
    WriteProtoMesg_AS( fp_to, gen );

  // close files
  fclose( fp_to );
  fclose( fp_from );
  return 0;
}

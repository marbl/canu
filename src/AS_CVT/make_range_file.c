
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
static char CM_ID[] = "$Id: make_range_file.c,v 1.2 2004-09-23 20:25:21 mcschatz Exp $";

/*********************************************************************/
// headers
/*********************************************************************/
// Standard header files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Project header files
#include "AS_global.h"

int main( int argc, char ** argv )
{
  FILE     * fp_in;
  FILE     * fp_out;
  char       line[1024];
  cds_uint32 lo;
  cds_uint32 hi;

  if( argc != 3 )
  {
    fprintf( stderr, "Usage: %s input.range output.range\n", argv[0] );
    return 1;
  }

  if( (fp_in = fopen( argv[1], "r" )) == NULL )
  {
    fprintf( stderr, "Failed to open file %s for reading\n", argv[1] );
    return 1;
  }

  if( (fp_out = fopen( argv[2], "w" )) == NULL )
  {
    fprintf( stderr, "Failed to open file %s for writing\n", argv[2] );
    return 1;
  }

  // first line is batch - read it & write it
  fgets( line, 1024, fp_in );
  fprintf( fp_out, "%s", line );

  // second line has lo:num__hi:num
  fgets( line, 1024, fp_in );
  sscanf( line, "lo:%u  hi:%u", &lo, &hi );
  fprintf( fp_out, "lo:1  hi:%u", hi );

  fclose( fp_in );
  fclose( fp_out );

  return 0;
}

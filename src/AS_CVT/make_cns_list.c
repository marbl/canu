
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
static char CM_ID[] = "$Id: make_cns_list.c,v 1.2 2004-09-23 20:25:21 mcschatz Exp $";

/*********************************************************************/
// headers
/*********************************************************************/
// Standard header files
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

int main( int argc, char ** argv )
{
  FILE * fp_in;
  FILE * fp_out;
  char   filename[10000];
  int    i = 0;
  int    num_cns_files = 0;
  int ch;

  if( argc != 3 )
  {
    fprintf( stderr, "Usage: %s input.cns output.txt\n", argv[0] );
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

  // one line file with space-delimited .cns filenames
  i = 0;
  while( ch = fgetc(fp_in), filename[i] = ch, ch != EOF )
  {
    if( isspace( filename[i] ) )
    {
      filename[i] = '\0';
      fprintf( fp_out, "%s\n", filename );
      num_cns_files++;
      i = 0;
    }
    else
      i++;
  }
  fprintf( stderr, "%s\n", filename );
  fprintf( stderr, "character %d is %c\n", i-1, filename[i-1] );
  
  if( i != 0 && filename[i-1] == 's' )
  {
    filename[i] = '\0';
    fprintf( fp_out, "%s\n", filename );
    num_cns_files++;
  }
  
  fclose( fp_in );
  fclose( fp_out );

  fprintf( stdout, "%d consensus files in component.\n", num_cns_files );
  if( num_cns_files == 0 )
    return 1;
  else
    return 0;
}

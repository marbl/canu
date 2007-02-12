
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
/**********************************************************************
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/nodups.c,v $
$Revision: 1.7 $
**********************************************************************/

/**********************************************************************
Module: get-olaps

Description: Reads the message file produced by the overlapper and
             extracts and outputs the overlap messages in a condensed
             format suitable for sorting.

             Adapted from Ian Dew's overlap regressor analyzer.

**********************************************************************/


/*********************************************************************/
// headers
// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

// project headers
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_distStore.h"
#include "AS_ORA_fragments.h"
#include "AS_ORA_overlaps.h"
#include "AS_ORA_statistics.h"
#include "AS_ORA_inlines.h"


// defines
#ifndef MAX_SEQUENCE_LENGTH
#define MAX_SEQUENCE_LENGTH 2048
#endif

#ifndef MAX_SOURCE_LENGTH
#define MAX_SOURCE_LENGTH 512
#endif


// function prototype

int  Insert
    (int x);
void  Remove_Dups
    ( char * input_filename, char * output_filename);

int  Frag_Size;
int  * Frag_Start;
int  Start_Ct = 0;


int main( int argc, char ** argv )

{
  char              * input_filename = NULL;
  char              * output_filename = NULL;

  // parse the command line parameters
  // use getopt(): see "man 3 getopt"
  {
    int ch, errflg = 0;
    optarg = NULL;
    while( !errflg && ((ch = getopt( argc, argv, "s:i:o:l:" )) != EOF) )
    {
      switch( ch )
      {
        case 'i':
          input_filename = optarg;
          break;
        case 'o':
          output_filename = optarg;
          break;
        case '?':
          fprintf( stderr, "Unrecognized option -%c\n", optopt );
        default:
          errflg++;
          break;
      }
    }

    // need fragstore_name & min_overlap and one or both of
    // input and output ovl filenames
    if( errflg != 0
          ||  input_filename == NULL
          ||  output_filename == NULL)
    {
     fprintf( stderr, "Usage: %s\n"
               " -i input-filename -o output-filename\n",
               argv[0] );
     exit (-1);
    }
  }

  Frag_Size = 10000;
  Frag_Start = (int *) malloc (Frag_Size * sizeof (int));
  assert (Frag_Start != NULL);
  Start_Ct = 0;

  Remove_Dups (input_filename, output_filename);

  fprintf( stderr, "Done.\n" );
  return 0;
}



int  Insert
    (int x)

//  Insert  x  into global sorted array  Frag_Start  if it's not there already.
//  Return  TRUE  if it was inserted;  FALSE  if it was there
//  already.  Global  Frag_Size  is the number of entries
//  allocated for  Frag_Start  and  Start_Ct  is the number of
//  entries currently in it.

  {
   int  i, lo, hi;

   lo = 0;
   hi = Start_Ct - 1;

   while  (lo <= hi)
     {
      int  mid = (lo + hi) / 2;

      if  (x == Frag_Start [mid])
          return  FALSE;
      else if  (x < Frag_Start [mid])
          hi = mid - 1;
        else
          lo = mid + 1;
     }

   if  (Start_Ct == Frag_Size)
       {
        Frag_Size *= 2;
        Frag_Start = (int *) realloc (Frag_Start, Frag_Size * sizeof (int));
        assert (Frag_Start != NULL);
       }

   for  (i = Start_Ct;  i > lo;  i --)
     Frag_Start [i] = Frag_Start [i - 1];
   Frag_Start [lo] = x;
   Start_Ct ++;

   return  TRUE;
  }



void  Remove_Dups
    ( char * input_filename, char * output_filename)

//  Put all fragment start coordinates, in ascending order,
//  into global array  Frag_Start  and set  Start_Ct  to how
//  many there are.

{
  FILE * infile = fopen (input_filename, "r" );
  FILE * outfile = fopen (output_filename, "w" );
  GenericMesg * gmesg = NULL;
  int64  next_accession = -1;
  char  * p;

  
  if( infile == NULL )
  {
    fprintf( stderr,
             "Failed to open message file %s for reading.\n",
             input_filename );
    exit (-1);
  }
  if( outfile == NULL )
  {
    fprintf( stderr,
             "Failed to open output file %s .\n",
             output_filename );
    exit (-1);
  }

  while (ReadProtoMesg_AS (infile, & gmesg) != EOF)
    if  (gmesg)
        {
         if  (gmesg->t == MESG_LKG)   // Strip link messages, too
             continue;
         if  (gmesg->t == MESG_FRG)
             {
              FragMesg * fmesg = NULL;
              int  a_end, b_end, celsim_start;

              // deal with the generic message as an overlap message
              fmesg = (FragMesg *) gmesg->m;
              p = strchr (fmesg -> source, '[');
              if  (sscanf (p, "[%d,%d],", & a_end, & b_end) != 2)
                  {
                   fprintf (stderr, "ERROR: couldn't read source\n");
                   exit (-1);
                  }

              if  (a_end < b_end)
                  celsim_start = a_end;
                else
                  celsim_start = b_end;

              if  (! Insert (celsim_start))
                  continue;

              if  (next_accession == -1)
                  next_accession = fmesg -> eaccession + 1;
                else
                  fmesg -> eaccession = next_accession ++;
             }
         WriteProtoMesg_AS (outfile, gmesg);
        }

  fclose (infile);
  fclose (outfile);

  return;
}



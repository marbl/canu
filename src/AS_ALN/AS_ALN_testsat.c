
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
/*********************************************************************/
// headers
// standard headers
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <stdlib.h>

// project headers
#include "AS_global.h"
#include "AS_ALN_aligners.h"
/*********************************************************************/


/*********************************************************************/
// defines
#define NUM_SATELLITES  10
#define MIN_SAT_SIZE     3
#define MAX_SAT_SIZE    12

#define NUM_FRAGMENTS  100
#define MIN_FRAG_SIZE  150
#define MAX_FRAG_SIZE  AS_READ_MAX_LEN

#define MIN_LENGTH      30
#define MAX_LENGTH      70
#define VARIATION        0.05

/*********************************************************************/

/*********************************************************************/
// globals
#define NUM_BASES        4
char Bases[NUM_BASES] = { 'a', 'c', 'g', 't' };
/*********************************************************************/

/*
  create a bunch of satellite sequences
  create a bunch of fragments
  insert sequences in some fragments, add subindels & remember where
    (put a couple distant stretches in some of them)
  call Initialize_Sataligner_AS
  loop over fragments & call Run_SatAligner_AS
  call Free_Sataligner_AS
  compare actual to found insertions
 */
int main( int argc, char ** argv )
{
  char             * sats[NUM_SATELLITES];
  InternalFragMesg   frags[NUM_FRAGMENTS];
  ScreenMatch      * actual_matches[NUM_FRAGMENTS];
  ScreenMatch      * reported_matches[NUM_FRAGMENTS];
  ScreenMatch      * match;
  SatAlignSetp       sat_set;
  uint32             i, j, k, l;
  uint32             seq_size;
  uint32             num_insert_bases;
  uint32             curr_sat;
  uint32             num_regions;
  uint32             subindel;

  // allocate & populate satellites sequences
  for( i = 0; i < NUM_SATELLITES; i++ )
  {
    seq_size = rand() % (MAX_SAT_SIZE - MIN_SAT_SIZE) + MIN_SAT_SIZE;
    sats[i] = (char *) calloc( seq_size + 1, 1 );
    if( sats[i] == NULL )
    {
      fprintf( stderr, "Failed to allocate satellite sequence size %u\n",
               seq_size );
      return 1;
    }
    for( j = 0; j < seq_size; j++ )
    {
      sats[i][j] = Bases[rand() % NUM_BASES];
    }
  }

  // allocate & populate fragment sequences
  for( i = 0; i < NUM_FRAGMENTS; i++ )
  {
    seq_size = rand() % (MAX_FRAG_SIZE - MIN_FRAG_SIZE) + MIN_FRAG_SIZE;
    frags[i].sequence = (char *) calloc( seq_size + 1, 1 );
    if( frags[i].sequence == NULL )
    {
      fprintf( stderr, "Failed to allocate fragment sequence size %u\n",
               seq_size );
      return 1;
    }
    for( j = 0; j < seq_size; j++ )
    {
      frags[i].sequence[j] = Bases[rand() % NUM_BASES];
    }
    // set some clear range
    frags[i].clear_rng.bgn = (int32) (seq_size * 0.05);
    frags[i].clear_rng.end = (int32) (seq_size * 0.95);

    // insert satellites into half of the fragments
    actual_matches[i] = NULL;
    if( rand() % 2 )
    {
      // insert one or two regions
      if( (frags[i].clear_rng.end - frags[i].clear_rng.bgn) >
          (MAX_FRAG_SIZE - MIN_FRAG_SIZE) / 2 )
        num_regions = 2;
      else
        num_regions = 1;
      for( j = 0; j < num_regions; j++ )
      {
        match = (ScreenMatch *) calloc( 1, sizeof( ScreenMatch ) );
        if( match == NULL )
        {
          fprintf( stderr, "Failed to allocate screen match structure.\n" );
          return 1;
        }
        match->next = actual_matches[i];
        actual_matches[i] = match;

        // overload use of other available storage
        match->portion_of.bgn = 0;  // substitutions
        match->portion_of.end = 0;  // insertions
        match->repeat_id      = 0;  // deletions
        
        // insert satellites
        num_insert_bases = rand() % (MAX_LENGTH - MIN_LENGTH) + MIN_LENGTH;
        match->where.bgn = (j == 0) ?
          (frags[i].clear_rng.bgn + num_insert_bases) :
          (frags[i].clear_rng.end - 2 * num_insert_bases);
        match->where.end = match->where.bgn;
        k = 0;

        while( k < num_insert_bases )
        {
          curr_sat = rand() % NUM_SATELLITES;
          for( l = 0; l < strlen( sats[curr_sat] ); l++ )
          {
            if( rand() % 100 > 100 - VARIATION * 100 )
            {
              subindel = rand() % 3;
              switch( subindel )
              {
                case 0:
                  // substitution - may substitute same base...
                  frags[i].sequence[match->where.end++] =
                    Bases[rand() % NUM_BASES];
                  k++;
                  match->portion_of.bgn++;
                  break;
                case 1:
                  // insertion
                  frags[i].sequence[match->where.end++] =
                    Bases[rand() % NUM_BASES];
                  l--;
                  k++;
                  match->portion_of.end++;
                  break;
                default:
                  // deletion
                  l--;
                  match->repeat_id++;
                  break;
              }
            }
            else
            {
              frags[i].sequence[match->where.end++] = sats[curr_sat][l];
              k++;
            }
          }
        }
      }
    }
  }

  // initialize the satellite detection system
  sat_set = Initialize_SatAligner_AS( NUM_SATELLITES,
                                      sats,
                                      0, 0, MIN_LENGTH + 10, VARIATION, 1 );
  if( sat_set == NULL )
  {
    fprintf( stderr, "Failed to initialize satellite aligner.\n" );
    return 1;
  }

  fflush( stderr );
  for( i = 0; i < NUM_FRAGMENTS; i++ )
  {
    printf( "Fragment %4u\n", i );
    
    // do the satellite detection
    reported_matches[i] = Run_SatAligner_AS( &frags[i], sat_set );
    
    printf( "Fragment ** %u **:\n", i );

    // print out the actual matches
    printf( "\tActual matches:\n" );
    if( actual_matches[i] == NULL )
    {
      printf( "\t\tnone.\n" );
    }
    else
    {
      match = actual_matches[i];
      while( match != NULL )
      {
        printf( "\t\tbgn: %4" F_COORDP ",\tend: %4" F_COORDP "\tsub/in/del: " F_COORD "/" F_COORD "/" F_UID "\n",
                 match->where.bgn, match->where.end,
                 match->portion_of.bgn, match->portion_of.end,
                 match->repeat_id );
        match = match->next;
      }
    }

    // print out the actual matches
    printf( "\tReported matches:\n" );
    if( reported_matches[i] == NULL )
    {
      printf( "\t\tnone.\n" );
    }
    else
    {
      match = reported_matches[i];
      while( match != NULL )
      {
        printf( "\t\tbgn: %4" F_COORDP ",\tend: %4" F_COORDP "\n",
                 match->where.bgn, match->where.end );
        match = match->next;
      }
    }

  }

  // free the satellite set - other stuff should be freed too.
  Free_SatAligner_AS( sat_set );
  
  return 0;
}  

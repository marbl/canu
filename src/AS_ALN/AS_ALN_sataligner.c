
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
#include <stdlib.h>
#include <math.h>

// project headers
#include "AS_global.h"
#include "AS_ALN_aligners.h"
/*********************************************************************/


/*********************************************************************/
// defines
/*********************************************************************/


/*********************************************************************/
// IUB/IUPAC - or anything - reverse characters
char ReverseBase[256] =
{
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
/* @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
   P   Q   R   S   T   U   V   W   X   Y   Z                */
  'n','T','n','G','n','n','n','C','n','n','n','n','n','n','n','n',
  'n','n','n','n','A','A','n','n','n','n','n','n','n','n','n','n',
  'n','t','n','g','n','n','n','c','n','n','n','n','n','n','n','n',
  'n','n','n','n','a','a','n','n','n','n','n','n','n','n','n','n',
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
  'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n'
};
/*********************************************************************/



/* Function:
     Initialize_SatAligner_AS
   Description:
     allocates & populates memory associated with SatAlignSet structure
   Return Value:
     valid pointer to poplulated SatAlignSet structure if ok, else NULL
   Parameters:
     uint32 num_sats   - number of satellite sequences
     char ** sats      - array of satellite sequences
     int32 repeat_id   - number to use as repeat & screen item ID
     int32 relevance   - number to use as relevance
     uint32 min_length - minimum length of an acceptable alignment
     float32 variation - variation allowed in alignment
     int do_reverse    - flag: 1 = also do satellite reverses, 0 = don't
*/
SatAlignSetp Initialize_SatAligner_AS( uint32 num_sats,
                                       char ** sats,
                                       int32 repeat_id,
                                       int32 relevance,
                                       uint32 min_length,
                                       float32 variation,
                                       int do_reverse )
{
  SatAlignSetp sat_set;
  uint32 i, k, s;
  uint32 sat_length;

  // allocate a SatAlignSet structure to populate
  sat_set = (SatAlignSetp) calloc( 1, sizeof( SatAlignSet ) );
  if( sat_set == NULL )
  {
    fprintf( stderr, "Error: Initialize_SatAligner_AS - " );
    fprintf( stderr, " Failed to allocate SatAlignSet structure.\n" );
    return NULL;
  }

  // easy stuff
  sat_set->repeat_id  = repeat_id;
  sat_set->relevance  = relevance;
  sat_set->min_length = min_length;
  sat_set->variation  = variation;

  // multiplier * path length is the minimum required alignment score
  // questionable, using a fudge factor here
  sat_set->multiplier =
    SAT_FUDGE_FACTOR * (SAT_MS + variation * (SAT_MMS - SAT_MS));

  // allocate the array of satellite sequence objects
  sat_set->sats =
    (SatelliteObject *) calloc( num_sats * ((do_reverse == 0) ? 1 : 2),
                                sizeof( SatelliteObject ) );
  if( sat_set->sats == NULL )
  {
    fprintf( stderr, "Error: Initialize_SatAligner_AS - " );
    fprintf( stderr, " Failed to allocate SatelliteObject array.\n" );
    Free_SatAligner_AS( sat_set );
    return NULL;
  }

  // loop over satellites & set up local copy of satellites
  for( i = 0; i < num_sats; i++ )
  {
    sat_length = strlen( sats[i] );
    sat_set->max_sat_length = MAX( sat_set->max_sat_length, sat_length );
    
    // do forward
    sat_set->sats[sat_set->num_sats].seq = (char *)malloc( sat_length );
    if( sat_set->sats[sat_set->num_sats].seq == NULL )
    {
      fprintf( stderr, "Error: Initialize_SatAligner_AS - " );
      fprintf( stderr, " Failed to allocate satellite sequence array.\n" );
      Free_SatAligner_AS( sat_set );
      return NULL;
    }
    sat_set->sats[sat_set->num_sats].length = sat_length;
    sat_set->sats[sat_set->num_sats].start = sat_set->num_sat_scores;
    strcpy( sat_set->sats[sat_set->num_sats].seq, sats[i] );
    sat_set->num_sats++;
    sat_set->num_sat_scores += sat_length + 1;

    if( do_reverse )
    {
      sat_set->sats[sat_set->num_sats].seq = (char *)malloc( sat_length );
      if( sat_set->sats[sat_set->num_sats].seq == NULL )
      {
        fprintf( stderr, "Error: Initialize_SatAligner_AS - " );
        fprintf( stderr, " Failed to allocate satellite sequence array.\n" );
        Free_SatAligner_AS( sat_set );
        return NULL;
      }
      sat_set->sats[sat_set->num_sats].length = sat_length;
      sat_set->sats[sat_set->num_sats].start = sat_set->num_sat_scores;
      // copy reverse of string
      for( k = 0; k < sat_length; k++ )
      {
        sat_set->sats[sat_set->num_sats].seq[k] =
          ReverseBase[(int)sats[i][sat_length - 1 - k]];
      }
      sat_set->num_sats++;
      sat_set->num_sat_scores += sat_length + 1;
    }
  }

  // score arrays
  for( i = 0; i < SAT_NUM_ARRAYS; i++ )
  {
    sat_set->scores[i] = (SatAlignScore *) calloc( sat_set->num_sat_scores,
                                                   sizeof( SatAlignScore ) );
    if( sat_set->scores[i] == NULL )
    {
      fprintf( stderr, "Error: Initialize_SatAligner_AS - " );
      fprintf( stderr, " Failed to allocate SatAlignScore array.\n" );
      Free_SatAligner_AS( sat_set );
      return NULL;
    }
  }

  // populate initializer array
  for( i = 0; i < sat_set->num_sats; i++ )
  {
    for( s = sat_set->sats[i].start, k = 0;
         k <= sat_set->sats[i].length;
         k++ )
    {
      sat_set->scores[SAT_NUM_ARRAYS - 1][s + k].score = k * SAT_MMS;
      sat_set->scores[SAT_NUM_ARRAYS - 1][s + k].start = 0;
      sat_set->scores[SAT_NUM_ARRAYS - 1][s + k].end = 0;
    }
  }

  return sat_set;
}
  
/* Function:
     Free_SatAligner_AS
   Description:
     frees memory associated with SatAlignSet structure
   Return Value:
     none
   Parameters:
     SatAlignSetp sat_set - pointer to structure to free
*/
void Free_SatAligner_AS( SatAlignSetp sat_set )
{
  uint32 i;

  if( sat_set != NULL )
  {
    // free satellite sequences
    if( sat_set->sats != NULL )
    {
      for( i = 0; i < sat_set->num_sats; i++ )
      {
        if( sat_set->sats[i].seq != NULL )
          free( sat_set->sats[i].seq );
      }
      free( sat_set->sats );
    }
    
    // free score arrays
    for( i = 0; i < SAT_NUM_ARRAYS; i++ )
    {
      if( sat_set->scores[i] != NULL )
        free( sat_set->scores[i] );
    }

    // free the structure itself
    free( sat_set );
  }
}
    

/* Function:
     Create_SatAligner_ScreenMatch_AS
   Description:
     Allocates & populates a ScreenMatch structure given a score
   Return Value:
     valid pointer to populated ScreenMatch structure if ok, else NULL
   Parameters:
     SatAlignSetp sat_set - pointer to structure to free
     SatAlignScore score  - score structure to use to populate ScreenMatch
*/
ScreenMatch * Create_SatAligner_ScreenMatch_AS( SatAlignSetp sat_set,
                                                SatAlignScore score )
{
  ScreenMatch * sma = (ScreenMatch *) malloc( sizeof( ScreenMatch ) );

  if( sma == NULL )
  {
    fprintf( stderr, "Error: Create_SatAligner_ScreenMatch_AS - " );
    fprintf( stderr, "Failed to allocate ScreenMatch\n" );
    return NULL;
  }

  // populate the screen match message
  sma->next           = NULL;
  sma->where.bgn      = score.start;
  sma->where.end      = score.end;
  sma->what           = sat_set->repeat_id;
  sma->repeat_id      = sat_set->repeat_id;
  sma->relevance      = sat_set->relevance;
  sma->portion_of.bgn = 0;
  sma->portion_of.end = score.end - score.start;
  sma->direction      = AS_FORWARD;
  
  return sma;
}

/* Function:
     Run_SatAligner_AS
   Description:
     Performs an alignment over the Kleene closure of a union of
     satellite sequences using dynamic programming. Zero, one or
     more local alignments may be returned in the form of a linked
     list of ScreenMatch structures. The code is based on pseudo-code
     written by Gene Myers.

     Call Initialize_SatAligner_AS to create and populate a
     SatAlignSet data structure before calling this routine. Then
     call it as many times as desired. When finished, call
     Free_SatAligner_AS to free the memory.
     
   Return Value:
     pointer to valid ScreenMatch structure (linked list) if ok and
     matches found, else NULL
   Parameters:
     InternalFragMesg * frag - pointer to internal fragment structure
     SatAlignSetp sat_set - pointer to initialized SatAlignSet structure
*/
ScreenMatch * Run_SatAligner_AS( InternalFragMesg * frag,
                                 SatAlignSetp sat_set )
{
  uint32 s, i, j, k;
  SatAlignScore * score_switcher = NULL;
  SatAlignScore null_score;
  SatAlignScore c_del;
  SatAlignScore c_ins;
  SatAlignScore c_sub;
  SatAlignScore temp_score;
  SatAlignScore max_score;
  ScreenMatch * new_sma = NULL;
  ScreenMatch * ret_sma = NULL;

  // set an initial starting score
  max_score.score = 0.0f;
  max_score.start = 0;
  max_score.end = 0;
  
  // populate first row of scores
  memcpy( sat_set->scores[0],
          sat_set->scores[SAT_NUM_ARRAYS - 1],
          sat_set->num_sat_scores * sizeof( SatAlignScore ) );
  
  // loop over fragment bases
  null_score.score = 0.0f;
  for( j = frag->clear_rng.bgn; j < frag->clear_rng.end; j++ )
  {
    null_score.start = null_score.end = j + 1;
    temp_score = null_score;
    
    // loop over satellites
    for( i = 0; i < sat_set->num_sats; i++ )
    {
      s = sat_set->sats[i].start;
      
      // give each satellite the best possible starting score
      sat_set->scores[1][s] = sat_set->scores[0][s];
      
      // loop over bases in satellite
      for( k = 1; k <= sat_set->sats[i].length; k++ )
      {
        c_del.score = sat_set->scores[0][s+k].score + SAT_MMS;
        c_del.start = sat_set->scores[0][s+k].start;
        c_ins.score = sat_set->scores[1][s+k-1].score + SAT_MMS;
        c_ins.start = sat_set->scores[1][s+k-1].start;
        c_sub.score = sat_set->scores[0][s+k-1].score +
          ((sat_set->sats[i].seq[k-1] == frag->sequence[j]) ?
           SAT_MS : SAT_MMS);
        c_sub.start = sat_set->scores[0][s+k-1].start;
        
        sat_set->scores[1][s+k] = MaxSatScore_AS( c_del,
                                                  MaxSatScore_AS( c_ins,
                                                                  c_sub ) );
      }
      
      // possibly update the max ending score over all satellites
      temp_score =
        MaxSatScore_AS( sat_set->scores[1][s+sat_set->sats[i].length],
                        temp_score );
    }
    
    // update the bestest best score as appropriate
    if( temp_score.score > max_score.score &&
        j - temp_score.start > sat_set->min_length )
    {
      max_score = temp_score;
      max_score.end = j + 1;
    }
    else if( max_score.score > 0.0f &&
             j > sat_set->max_sat_length + max_score.end )
    {
      if( max_score.end - max_score.start >= sat_set->min_length &&
          max_score.score >=
          (max_score.end - max_score.start) * sat_set->multiplier )
      {
        new_sma = Create_SatAligner_ScreenMatch_AS( sat_set, max_score );
        if( new_sma != NULL )
        {
          new_sma->next = ret_sma;
          ret_sma = new_sma;
        }
      }
      max_score = null_score;
      temp_score = null_score;
      memcpy( sat_set->scores[1],
              sat_set->scores[SAT_NUM_ARRAYS - 1],
              sat_set->num_sat_scores * sizeof( SatAlignScore ) );
      for( i = 0; i < sat_set->num_sat_scores; i++ )
        sat_set->scores[1][i].start = j + 1;
    }
    
#ifdef DEBUG      
    fprintf( stderr, "%6d %c   %f   %f   %4u   %4u\n",
             j, frag->sequence[j],
             temp_score.score, max_score.score,
             max_score.start, max_score.end );
#endif
    
    for( i = 0; i < sat_set->num_sats; i++ )
      sat_set->scores[1][sat_set->sats[i].start] =
        MaxSatScore_AS( temp_score, null_score );
    
    // do second round
    for( i = 0; i < sat_set->num_sats; i++ )
    {
      s = sat_set->sats[i].start;
      sat_set->scores[1][s] = MaxSatScore_AS( temp_score,
                                              sat_set->scores[1][s] );
      k = 1;
      while( k <= sat_set->sats[i].length &&
             sat_set->scores[1][s+k-1].score + SAT_MMS >
             sat_set->scores[1][s+k].score )
      {
        sat_set->scores[1][s+k].score =
          sat_set->scores[1][s+k-1].score + SAT_MMS;
        sat_set->scores[1][s+k].start = sat_set->scores[1][s+k-1].start;
        k++;
      }
    }
    
    // switch score arrays around
    score_switcher = sat_set->scores[0];
    sat_set->scores[0] = sat_set->scores[1];
    sat_set->scores[1] = score_switcher;
  }

  // record sufficient match at the end
/*
  if( max_score.score > SAT_MIN_SCORE * (j - max_score.start)
      && j - max_score.start > sat_set->min_length )
*/
  if( max_score.end - max_score.start >= sat_set->min_length &&
      max_score.score >=
      (max_score.end - max_score.start) * sat_set->multiplier )
  {
    new_sma = Create_SatAligner_ScreenMatch_AS( sat_set, max_score );
    if( new_sma != NULL )
    {
      new_sma->next = ret_sma;
      ret_sma = new_sma;
    }
  }
  
  return ret_sma;
}

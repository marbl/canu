
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/AS_ORA_statistics.c,v $
$Revision: 1.6 $
**********************************************************************/


/*********************************************************************/
// headers
// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// project headers
#include "AS_global.h"
#include "AS_ORA_statistics.h"
#include "AS_ORA_inlines.h"
/*********************************************************************/

#define FN_SUFFIX  "fn.cgm"
#define CFN_SUFFIX "cfn.cgm"
#define FP_SUFFIX  "fp.cgm"

/*********************************************************************/
// would-be inlines
void switch_hangs( int32 * a, int32 * b )
{
  int32 temp = *a;
  *a = *b;
  *b = temp;
  if( *a * *b < 0 )
  {
    *a = -*a;
    *b = -*b;
  }
}
/*********************************************************************/


/*********************************************************************/
/* Function:
     AllocateHistogram
   Description:
     Allocates a histogram array structure
   Return Value:
     valid Histogramp pointer if ok
   Parameters:
     uint32 num_bins: number of bins for histogram
*/
Histogramp AllocateHistogram( uint32 num_bins )
{
  Histogramp histo;

  histo = (Histogramp) calloc( 1, sizeof( Histogram ) );
  if( histo == NULL )
  {
    fprintf( stderr, "Failed to allocate histogram object.\n" );
    return NULL;
  }
  
  histo->bins = (uint32 *) calloc( num_bins, sizeof( uint32 ) );
  if( histo->bins == NULL )
  {
    fprintf( stderr,
             "Failed to allocate %u histogram bins.\n",
             num_bins );
    FreeHistogram( histo );
    return NULL;
  }
  histo->num_bins = num_bins;
  
  return histo;
}

/* Function:
     FreeHistogram
   Description:
     Frees a Histogram array structure
   Return Value:
     none
   Parameters:
     Histogramp histo:  pointer to Histogram structure
*/
void FreeHistogram( Histogramp histo )
{
  if( histo != NULL )
  {
    if( histo->bins != NULL )
      free( histo->bins );
    free( histo );
  }
}

/* Function:
     InitializeStatistics
   Description:
     Sets initial values in an OverlapStatistics structure
   Return Value:
     none
   Parameters:
     OverlapStatisticsp statistics: pointer to OverlapStatistics structure
*/
void InitializeStatistics( OverlapStatisticsp stats )
{
  stats->min_overlap = CDS_INT16_MAX;
  stats->max_overlap = 0;
  stats->num_seq_gaps = 0;
  stats->seq_gap_bases = 0;
  stats->num_ovl_gaps = 0;
  stats->ovl_gap_bases = 0;
  stats->num_chunks = 1;
  stats->num_true = 0;
  stats->num_found = 0;
  stats->num_true_found = 0;
  stats->num_true_found_2x = 0;
  stats->num_near_misses = 0;
  stats->num_false_negatives = 0;
  stats->histo_false_negatives = NULL;
  stats->num_fn_critical = 0;
  stats->histo_fn_critical = NULL;
  stats->num_false_positives = 0;
  stats->histo_false_positives = NULL;
  stats->max_a_hang_delta = 0;
  stats->max_b_hang_delta = 0;
  stats->mean_a_hang_delta = 0.0f;
  stats->mean_b_hang_delta = 0.0f;
  stats->std_a_hang_delta = 0.0f;
  stats->std_b_hang_delta = 0.0f;
  stats->repeats = NULL;
}

/* Function:
     UpdateRepeatStats
   Description:
     updates the repeat linked list of composite repeat id bit vectors
     based on a new 'found' overlap that is not a 'true' overlap
   Return Value:
     0 if ok
   Parameters:
     OverlapStatisticsp stats: pointer statistics structure
     FragmentArrayp fragments: pointer to fragment array structure
     OverlapMesg * osp: overlap message of false positive
*/
int UpdateRepeatStats( OverlapStatisticsp stats,
                       FragmentArrayp fragments,
                       OverlapMesg * osp )
{
  RepeatObjectp s_repeat = stats->repeats;
  RepeatObjectp new_repeat = NULL;
  RepeatID repeat_id =
    fragments->data[osp->aifrag - fragments->min_id].repeat_id &
    fragments->data[osp->bifrag - fragments->min_id].repeat_id;

  /* TBD - check for location of repeat in fragment vis. overlap
  // Determine which repeats may have influenced false positive
  RepeatObjectp a_repeat =
    fragments->data[osp->aifrag - fragments->min_id].repeats;
  RepeatObjectp b_repeat = NULL;
  while( a_repeat != NULL )
  {
    b_repeat = fragments->data[osp->bifrag - fragments->min_id].repeats;
    while( b_repeat != NULL )
    {
      if( a_repeat->repeat_id == b_repeat->repeat_id )
      {
        // if the two repeats sync up in the overlap, add th repeat_id
        if( osp->orientation == AS_NORMAL )
        {
          // both frags are oriented as in frag file
        }
        else if( osp->orientation == AS_INNIE )
        {
          // A is oriented as in frag file
          // B is reversed relative to expression in frag file
        }
        else
        {
          // A is reversed relative to expression in frag file
          // B is oriented as in frag file
        }
      }
      b_repeat = b_repeat->next;
    }
    a_repeat = a_repeat->next;
  }
  */
  
  if( repeat_id != 0 )
  {
    // find the repeat_id in the stats repeats
    while( s_repeat != NULL && s_repeat->repeat_id < repeat_id )
      s_repeat = s_repeat->next;

    // if not found
    if( s_repeat == NULL || s_repeat->repeat_id != repeat_id )
    {
      // establish new repeat object
      new_repeat = GetRepeatObject( fragments->repeat_heap );
      if( new_repeat == NULL )
      {
        fprintf( stderr, "Failed to get repeat object from heap.\n" );
        return 1;
      }
      
      // populate it
      new_repeat->repeat_id = repeat_id;
      new_repeat->count = 1;
      
      // insert it in the stats->repeats list
      InsertionSortRepeat( &(stats->repeats), new_repeat );
    }
    else
    {
      // found: increment the count
      s_repeat->count++;
    }
  }
  
  return 0;
}

/* Function:
     AddToStats
   Description:
     Adds details of a true and found overlap comparison to stats
   Return Value:
     none
   Parameters:
     OverlapStatistcsp stats: pointer to statistics structure
     OverlapObjectp oop1: pointer to first overlap
     OverlapObjectp oop2: pointer to second overlap
     uint64 * sum_a_hang_deltas: pointer to intermediate accumulation variable
     uint64 * sum_a_hang_deltas2: pointer to intermediate accumulation variable
     uint64 * sum_b_hang_deltas: pointer to intermediate accumulation variable
     uint64 * sum_b_hang_deltas2: pointer to intermediate accumulation variable
*/
void AddToStats( OverlapStatisticsp stats,
                 OverlapObjectp t_oop,
                 OverlapObjectp f_oop,
                 uint64 * sum_a_hang_deltas,
                 uint64 * sum_a_hang_deltas2,
                 uint64 * sum_b_hang_deltas,
                 uint64 * sum_b_hang_deltas2 )
{
  uint64 d1, d2, d3;
  uint64 ahgd, bhgd;

  // d1 is normal comparison
  d1 = abs( t_oop->a_hang - f_oop->a_hang ) +
    abs( t_oop->b_hang - f_oop->b_hang );
  
  // d2 is reversed dovetail comparison
  d2 = abs( t_oop->a_hang - f_oop->b_hang ) +
    abs( t_oop->b_hang - f_oop->a_hang );
  
  // d3 is reversed containment comparison
  d3 = abs( t_oop->a_hang + f_oop->b_hang ) +
    abs( t_oop->b_hang + f_oop->a_hang );

  if( d1 < MIN( d2, d3 ) )
  {
    ahgd = abs( t_oop->a_hang - f_oop->a_hang );
    bhgd = abs( t_oop->b_hang - f_oop->b_hang );
  }
  else if( d2 < d3 )
  {
    ahgd = abs( t_oop->b_hang - f_oop->a_hang );
    bhgd = abs( t_oop->a_hang - f_oop->b_hang );
  }
  else
  {
    ahgd = abs( t_oop->b_hang + f_oop->a_hang );
    bhgd = abs( t_oop->a_hang + f_oop->b_hang );
  }

  // if a_hangs are 0, then b_hangs may have opposite signs because
  // one of us claimed containment & the other claimed dovetail
  if( t_oop->a_hang == 0 && f_oop->b_hang == 0 &&
      t_oop->b_hang * f_oop->b_hang < 0 )
  {
    ahgd = 0;
    bhgd = abs( t_oop->b_hang + f_oop->b_hang );
  }
  
  // collect stats
  // first a_hang
  stats->max_a_hang_delta = MAX( ahgd, stats->max_a_hang_delta );
  *sum_a_hang_deltas += ahgd;
  *sum_a_hang_deltas2 += ahgd * ahgd;

  // now b_hang
  stats->max_b_hang_delta = MAX( bhgd, stats->max_b_hang_delta );
  *sum_b_hang_deltas += bhgd;
  *sum_b_hang_deltas2 += bhgd * bhgd;
}

/* Function:
     FinalizeStatistics
   Description:
     loops through overlaps to accumulate statistics, then computes them
   Return Value:
     none
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     OverlapStatisticsp stats: pointer statistics structure
*/
void FinalizeStatistics( FragmentArrayp fragments, OverlapStatisticsp stats )
{
  uint32 i;
  OverlapObjectp oop;
  uint32 *hfn = stats->histo_false_negatives->bins;
  uint32 *hfnc = stats->histo_fn_critical->bins;
  uint64 sum_a_hang_deltas = 0;
  uint64 sum_a_hang_deltas2 = 0;
  uint64 sum_b_hang_deltas = 0;
  uint64 sum_b_hang_deltas2 = 0;
  uint32 overlap_size;

  // loop over fragments to gather statistics
  for( i = 0; i < fragments->num_items; i++ )
  {
    oop = fragments->data[i].overlaps;
    while( oop )
    {
      switch( oop->found_count )
      {
        case 0:
          stats->num_false_negatives++;
          overlap_size =
            GetOverlapInterval( &(fragments->data[i]),
                    &(fragments->data[oop->frag_id - fragments->min_id])) -
            stats->min_overlap;
          hfn[overlap_size]++;
          if( IsCriticalOverlap( fragments,
                                 i,
                                 oop->frag_id - fragments->min_id ) )
          {
            stats->num_fn_critical++;
            hfnc[overlap_size]++;
          }
          
          break;
        default:
          stats->num_true_found++;
          stats->num_true_found_2x += oop->found_count - 1;
          AddToStats( stats, oop, oop->best_match,
                      &sum_a_hang_deltas, &sum_a_hang_deltas2,
                      &sum_b_hang_deltas, &sum_b_hang_deltas2 );
          break;
      }
      oop = oop->next;
    }
  }
  
  stats->mean_a_hang_delta =
    ((float) sum_a_hang_deltas) / stats->num_true_found;
  stats->mean_b_hang_delta =
    ((float) sum_b_hang_deltas) / stats->num_true_found;

  stats->std_a_hang_delta =
    sqrt( ((double) sum_a_hang_deltas2) / stats->num_true_found -
           stats->mean_a_hang_delta * stats->mean_a_hang_delta );
  stats->std_b_hang_delta =
    sqrt( ((double) sum_b_hang_deltas2) / stats->num_true_found -
           stats->mean_b_hang_delta * stats->mean_b_hang_delta );
}

/* Function:
     PrintStatistics
   Description:
     prints statistics to stdout
   Return Value:
     none
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     char * fragstore_name: name of fragstore
     char * overlap_filename: name of overlap message file
     OverlapStatisticsp stats: pointer to statistics structure
*/
void PrintStatistics( FragmentArrayp fragments,
                      char * fragstore_name,
                      char * overlap_filename,
                      OverlapStatisticsp stats )
{
  int i;
  int max_fn, max_fp, max_index;
  uint32 *hfn = stats->histo_false_negatives->bins;
  uint32 *hfnc = stats->histo_fn_critical->bins;
  uint32 *hfp = stats->histo_false_positives->bins;
  RepeatObjectp repeat;
  uint32 prifps;
  char repeat_string[256];
  char fn_filename[256];
  char cfn_filename[256];
  char fp_filename[256];
  FILE * fnfp;
  FILE * cfnfp;
  FILE * fpfp;
  char * chr;

  // set up celagram filenames
  strcpy( fn_filename, overlap_filename );
  chr = strrchr( fn_filename, '.' );
  if( chr != NULL )
    chr[1] = '\0';
  strcpy( cfn_filename, fn_filename );
  strcpy( fp_filename, fn_filename );
  strcat( fn_filename, FN_SUFFIX );
  strcat( cfn_filename, CFN_SUFFIX );
  strcat( fp_filename, FP_SUFFIX );
  fnfp = fopen( fn_filename, "w" );
  cfnfp = fopen( cfn_filename, "w" );
  fpfp = fopen( fp_filename, "w" );
  if( fnfp == NULL || cfnfp == NULL || fpfp == NULL )
  {
    fprintf( stderr, "Warning: failed to open one of %s, %s, or %s.\n",
             fn_filename, cfn_filename, fp_filename );
  }
  if( fnfp != NULL )
    fprintf( fnfp, "bin 0 = %u bases\n", stats->min_overlap );
  if( cfnfp != NULL )
    fprintf( cfnfp, "bin 0 = %u bases\n", stats->min_overlap );
  if( fpfp != NULL )
    fprintf( fpfp, "bin 0 = %u bases\n", stats->min_overlap );

  // title
  printf( "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n" );
  printf( "Comparison of 'true' overlaps computed from fragment indices\n" );
  printf( " with 'found' overlaps determined by the overlap detector.\n" );

  // info sources
  printf( "\n-- Data sources\n" );
  printf( "  Fragment Store: %s\n", fragstore_name );
  printf( "  Overlap file:   %s\n", overlap_filename );
  printf( "    Number of fragments: %10u\n", fragments->num_items );
  printf( "    Minimum overlap size of interest: %4u\n",
          stats->min_overlap );
  printf( "    Maximum abs( A-hang + B-hang ) difference for"
          " overlap match: %4d\n",
          MAX_OVERLAP_MISMATCH );

  // sequencing sample
  printf( "\n-- Sampling\n" );
  printf( "   Number of sequencing gaps:       %6u\n", stats->num_seq_gaps );
  printf( "   Bases of sequencing gaps:        %6u\n", stats->seq_gap_bases );
  printf( "   Number of non-overlapable gaps:  %6u\n", stats->num_ovl_gaps );
  printf( "   Bases of non-overlapable gaps:   %6u\n", stats->ovl_gap_bases );
  
  // basic numbers
  printf( "\n-- Simple counts\n" );
  printf( "   Found overlaps (based on overlapper output):  %10" F_U64P "\n",
           stats->num_found );
  printf( "   True overlaps (based on fragment indices):    %10" F_U64P "\n",
           stats->num_true );
  printf( "                                     ----------------------\n" );
  printf( "                                     Difference: %10lu\n",
           labs( stats->num_true - stats->num_found ) );

  // found overlaps
  printf( "\n-- Composition of found overlaps\n" );
  printf( "  True overlaps that were found at least once:   %10" F_U64P "\n",
          stats->num_true_found );
  printf( "  True overlaps that were found more than once:  %10" F_U64P "\n",
          stats->num_true_found_2x );
  printf( "  Near misses:                                   %10" F_U64P "\n",
          stats->num_near_misses );
  printf( "  False positives:                               %10" F_U64P "\n",
          stats->num_false_positives );
  printf( "                                          -----------------\n" );
  printf( "                                          Total: %10" F_U64P "\n",
          stats->num_true_found +
          stats->num_true_found_2x +
          stats->num_near_misses +
          stats->num_false_positives );
  
  // true/found overlaps
  printf( "\n-- True overlaps that were found\n" );
  printf( "  Statistical comparison of A-hang and B-hang values\n" );
  printf( "    between matching true & found overlaps\n" );
  printf( "             Mean          Stdev\n" );
  printf( "  A-hang:  %7.4e   %7.4e\n",
           stats->mean_a_hang_delta,
           stats->std_a_hang_delta );
  printf( "  B-hang:  %7.4e   %7.4e\n",
           stats->mean_b_hang_delta,
           stats->std_b_hang_delta );

  // false positives
  printf( "\n-- False positives (limited by annotations in src field)\n" );
  // count the possible repeat-induced false positives
  repeat = stats->repeats;
  prifps = 0;
  while( repeat )
  {
    prifps += repeat->count;
    repeat = repeat->next;
  }
  printf( "  Possibly repeat-induced:       %10u\n", prifps );
  printf( "  Not repeat-induced:            %10" F_U64P "\n",
          stats->num_false_positives - prifps );

  // repeats
  printf( "\n-- Repeats\n" );
  if( fragments->repeats == NULL )
  {
    printf( "  No repeats detected.\n" );
  }
  else
  {
    // print basic repeat/fragment data
    printf( "     ID       Fragments containing\n" );
    repeat = fragments->repeats;
    while( repeat )
    {
      printf( "     %c               %5" F_U64P "\n",
              GetRepeatID( repeat->repeat_id ),
              repeat->count );
      repeat = repeat->next;
    }

    // print repeat/overlap data
    if( stats->repeats != NULL )
    {
      printf( "\n" );
      printf( "  Repeat(s)                    False positives containing\n" );
      repeat = stats->repeats;
      while( repeat )
      {
        GetRepeatString( repeat->repeat_id, repeat_string );
        printf( "   %s", repeat_string );
        // allow 30 spaces for prettiness
        for( i = 0; i < 30 - strlen( repeat_string ); i++ )
          printf( " " );
        printf( "%12" F_U64P "\n", repeat->count );
        repeat = repeat->next;
      }
    }
  }
  
  // false negatives
  printf( "\n-- False negatives\n" );
  printf( "  False negatives:                                  %10" F_U64P "\n",
           stats->num_false_negatives );
  printf( "  Critical false negatives:                         %10" F_U64P "\n",
           stats->num_fn_critical );
  printf( "  Number of separate, non-contained overlap graphs: %10u\n",
          stats->num_chunks );
  
  // histogram
  // find the meaningful end of the histogram
  // to avoid printing unnecessary zeros
  max_fn = stats->histo_false_negatives->num_bins - 1;
  while( max_fn && hfn[max_fn] == 0 )
  {
    max_fn--;
  }
  max_fn++;

  max_fp = stats->histo_false_positives->num_bins - 1;
  while( max_fp && hfp[max_fp] == 0 )
  {
    max_fp--;
  }
  max_fp++;

  max_index = MAX( max_fn, max_fp );

  printf( "\n-- Histogram of overlap sizes:\n" );
  printf( "(False negatives are in celagram file %s)\n", fn_filename );
  printf( "(Critical false negatives are in celagram file %s)\n",
          cfn_filename );
  printf( "(False positives are in celagram file %s)\n", fp_filename );
  
  printf( "                          Critical\n" );
  printf( "  Overlap      False        False        False\n" );
  printf( "   sizes     Negatives    Negatives    Positives\n" );
  for( i = 0; i < max_index; i++ )
  {
    if( i < stats->histo_false_positives->num_bins &&
        i < stats->histo_false_negatives->num_bins )
    {
      printf( "   %4u       %7d       %4d        %7d\n",
              stats->min_overlap + i, hfn[i], hfnc[i], hfp[i] );
      if( fnfp != NULL )
        fprintf( fnfp, "%u\n", hfn[i] );
      if( cfnfp != NULL )
        fprintf( cfnfp, "%u\n", hfnc[i] );
      if( fpfp != NULL )
        fprintf( fpfp, "%u\n", hfp[i] );
    }
    else if( i < stats->histo_false_positives->num_bins )
    {
      printf( "   %4u             0          0        %7d\n",
              stats->min_overlap + i, hfp[i] );
      if( fpfp != NULL )
        fprintf( fpfp, "%u\n", hfp[i] );
    }
    else
    {
      printf( "   %4u       %7d       %4d              0\n",
              stats->min_overlap + i, hfn[i], hfnc[i] );
      if( fnfp != NULL )
        fprintf( fnfp, "%u\n", hfn[i] );
      if( cfnfp != NULL )
        fprintf( cfnfp, "%u\n", hfnc[i] );
    }
  }

  if( fnfp != NULL )
    fclose( fnfp );
  if( cfnfp != NULL )
    fclose( cfnfp );
  if( fpfp != NULL )
    fclose( fpfp );
}

/* Function:
     FreeStatistics
   Description:
     frees histogram in statistics structure
   Return Value:
     none
   Parameters:
     OverlapStatisticsp stats: pointer to statistics structure
*/
void FreeStatistics( OverlapStatisticsp stats )
{
  FreeHistogram( stats->histo_false_negatives );
  FreeHistogram( stats->histo_fn_critical );
  FreeHistogram( stats->histo_false_positives );
}



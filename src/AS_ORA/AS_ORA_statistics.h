
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
#ifndef AS_ORA_STATISTICS_H
#define AS_ORA_STATISTICS_H
/**********************************************************************
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/AS_ORA_statistics.h,v $
$Revision: 1.3 $
**********************************************************************/


/*********************************************************************/
// headers
// project headers
#include "AS_global.h"
#include "AS_ORA_fragments.h"
#include "AS_ORA_overlaps.h"
/*********************************************************************/


/*********************************************************************/
// Statistics-related structures
// structure for histograms
typedef struct
{
  uint32 * bins;       // array of bins for counts
  uint32   num_bins;   // number of bins
} Histogram;
typedef Histogram * Histogramp;

// structure to hold statistics
typedef struct
{
  uint32     min_overlap;            // minimum overlap interval of interest
  uint32     max_overlap;            // maximum overlap interval not found
  uint32     num_seq_gaps;           // number of sequencing gaps
  uint32     seq_gap_bases;          // number of bases in sequencing gaps
  uint32     num_ovl_gaps;           // number of non-overlapable gaps
  uint32     ovl_gap_bases;          // number of bases in non-overlapable gaps
  uint32     num_chunks;            // number of distinct, non-contained chunks
  uint64     num_true;               // number of true overlaps
  uint64     num_found;              // number of found overlaps
  uint64     num_true_found;         // number in intersection of sets
  uint64     num_true_found_2x;      // number of overlaps reported twice
  uint64     num_near_misses;        // number of nearly correct found overlaps
  uint64     num_false_negatives;    // number of critical false negatives
  Histogramp histo_false_negatives;  // histogram of false negatives
  uint64     num_fn_critical;        // number of critical false negatives
  Histogramp histo_fn_critical;      // histogram of critical false negatives
  uint64     num_false_positives;    // number of false positives
  Histogramp histo_false_positives;  // histogram false positives
  uint32     max_a_hang_delta;       // max abs value diff of a_hang deltas
  uint32     max_b_hang_delta;       // max abs value diff of b_hang deltas
  float      mean_a_hang_delta;      // mean abs value diff of a_hang deltas
  float      mean_b_hang_delta;      // mean abs value diff of b_hang detlas
  float      std_a_hang_delta;       // stddev of a_hang deltas
  float      std_b_hang_delta;       // stddev of b_hang deltas
  RepeatObjectp repeats;             // false positive-related repeats
} OverlapStatistics;
typedef OverlapStatistics * OverlapStatisticsp;
/*********************************************************************/


/*********************************************************************/
// function declarations
/* Function:
     AllocateHistogram
   Description:
     Allocates a histogram array structure
   Return Value:
     valid Histogramp pointer if ok
   Parameters:
     uint32 num_bins: number of bins for histogram
*/
Histogramp AllocateHistogram( uint32 num_bins );

/* Function:
     FreeHistogram
   Description:
     Frees a Histogram array structure
   Return Value:
     none
   Parameters:
     Histogramp histo:  pointer to Histogram structure
*/
void FreeHistogram( Histogramp histo );

/* Function:
     InitializeStatistics
   Description:
     Sets initial values in an OverlapStatistics structure
   Return Value:
     none
   Parameters:
     OverlapStatisticsp statistics: pointer to OverlapStatistics structure
*/
void InitializeStatistics( OverlapStatisticsp stats );

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
                       OverlapMesg * osp );

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
                 OverlapObjectp oop1,
                 OverlapObjectp oop2,
                 uint64 * sum_a_hang_deltas,
                 uint64 * sum_a_hang_deltas2,
                 uint64 * sum_b_hang_deltas,
                 uint64 * sum_b_hang_deltas2 );

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
void FinalizeStatistics( FragmentArrayp fragments, OverlapStatisticsp stats );

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
                      OverlapStatisticsp stats );

/* Function:
     FreeStatistics
   Description:
     frees histogram in statistics structure
   Return Value:
     none
   Parameters:
     OverlapStatisticsp stats: pointer to statistics structure
*/
void FreeStatistics( OverlapStatisticsp stats );
/*********************************************************************/

#endif // #ifndef AS_ORA_STATISTICS_H

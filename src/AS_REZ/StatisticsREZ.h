
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
/* 	$Id: StatisticsREZ.h,v 1.4 2005-03-22 19:49:25 jason_miller Exp $	 */

/**********************************************************************************
 *  StatisticsRez.h
 *  
 *  Knut Reinert 12/99
 *
 *  This is the interface from CGW to REZ
 *
 **********************************************************************************/

#ifndef STATISTICS_REZ_H
#define STATISTICS_REZ_H

#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "AS_UTL_Var.h"
#include "AS_global.h"
#include "InputDataTypes_CGW.h"

typedef enum {
  ALL_SCAFFOLDS = -1,
  NO_SCAFFOLD = -2
} ScaffoldID;


/* the following struct is the atomic object for the
   walking statistics. it contains the information
   about a single gap */

typedef struct
{
  int lCID, rCID;
  // the contig ids of the two flanking contigs
  union{
    int64 all;
    struct{
      unsigned int walked:1; 
      // a boolean flag that indicates whether the gap was succesfully walked or not
      unsigned int trivial:1;
      // a boolean flag that indicates whether the gap was trivial
      unsigned int maxedOut:1;
      // a boolean flag that indicates whether we gave up on the gap exceeding
      // the upper bound for edges allowed to explore
      unsigned int walkedTooLong:1;
      // a boolean flag that indicates whether the gap was walked but the walk was too long
      unsigned int walkedTooShort:1;
      // a boolean flag that indicates whether the gap was walked but the walk was too short
      unsigned int tooShort:1;
      // a boolean that indicates that we did not even try too walk the gap, because
      // its overlap is too short
    }bits;
  }
  flags;
    
  LengthT gapEstimate;
  LengthT gapLength;
  // two LengthTs the first contains the estimation of the gap length
  // the second the actual length of a walked gap. For unwalked gaps
  // this is set to {0.0,0.0} 
    
  float bestTooShort;
  float bestTooLong;

  int exploredEdges;
  // the number of explored edges in a gap
  int walkedChunks;
  // the number of chunks on the best walk

  int 
    uuChunks,
    urChunks,
    ruChunks,
    rrChunks;
  // counts of the walked chunks. These numbers
  // reflect the result of the functions Is_UU, Is_UR, Is_RU, and Is_RR
  
}
GapStatisticsT;


/* ------------------------------------------------ */
/* functions to manipulate the GapStatistics struct */
/* ------------------------------------------------ */

void init_gap_stat_struct(GapStatisticsT* g);
// initializes all counters to zero

void print_gap_stat_struct(GapStatisticsT* g, FILE* file);
// prints a readable form of a gap statistic struct


/* the next struct contains walking information on a scaffold
   basis. Any information in the struct can be recomputed
   using information in the GapStatisticsT struct */

VA_DEF(GapStatisticsT)


typedef  struct
{
  int scaffoldID;
  // contains the id of the scaffold
  // the statistics are based on
  // if -1 then the statistics are accumulativ for all scaffolds
  // if -2 Scaffold is dead or not allocated

  VA_TYPE(GapStatisticsT) *GapStats;
  // a variable array containing statistics structs for all gaps in the scaffold

  int 
    uuChunks,
    urChunks,
    ruChunks,
    rrChunks;
  // counts of the walked chunks. These numbers
  // reflect the result of the functions Is_UU, Is_UR, Is_RU, and Is_RR

  int negativExploredEdges;
  // the total number of explored edges in successful negativ walks in the scaffold
  int negativWalkedChunks;
  // the total number of chunks on the best path in a successful negativ walk
  int smallExploredEdges;
  // the total number of explored edges in successful small walks in the scaffold
  int smallWalkedChunks;
  // the total number of chunks on the best path in a successful small walk
  int bigExploredEdges;
  // the total number of explored edges in successful big walks in the scaffold
  int bigWalkedChunks;
  // the total number of chunks on the best path in a successful big walk
  
  int insertedChunks;
  // number of the chunk actually inserted in the scaffold for all gaps

  int 
    bpsMaxGap,
    /* the size in bps of the biggest gap walked */
    bpsTried,
    /* the sum of the scaffold lengths of walked scaffolds */
    bpsWalked,
    /* the total number of contiguous base pairs added 
       through the walks of this scaffold, i.e. not adding up chunklengths,
       but adding up the corrected gap estimate means */
    bpsNotWalked,
    /* the total number of the estimated gap lengths that we 
       were not able to walk */
    bigGapsWalked,
    /* a gap is considered big if the estimated length is bigger than 
       SWITCH_THRESHOLD (defined in GWdrivers.h) this counter counts 
       the number of walked big gaps */
    bigGapsNotWalked,
    /* number of unwalked big gaps */
    smallGapsWalked,
    /* a gap is considered small if the estimated length is less than 
       SWITCH_THRESHOLD and greater than 0. 
       This counter counts the number of walked small gaps */
    smallGapsNotWalked,
    /* number of unwalked small gaps */
    negativGapsWalked,
    /* a gap is considered negativ if the estimated length is less than 0 
       (big surprise). This counter counts the number of walked negative gaps */
    negativGapsNotWalked,
    /* number of unwalked negativ gaps */
    gapsTooShort,
    /* number of gaps that we did not even try to walk, because it was too short */
    walkedMaxedOut,
    /* the number of gaps in which we explored more than MAXWALKCALLS 
       (defined in GapWalkerREZ.h) and hence stopped exploring the graph */
    walkedTooShort,
    /* the number of gaps in which we actually reached the other side of the gap
       but had a too short walk compared to the lower bound computed for the gap
       (depends on the level of walking) */
    walkedTooLong,
    /* the number of gaps in which we actually reached the other side of the gap
       but had a too long walk compared to the upper bound computed for the gap
       (3*STDDEV) */
    walkedTrivial;
    /* the number of tvival gaps, that means which have an edge between the 
       left and the right chunk */
     
    /* BEWARE that a single gap can contribute to several of the counts above,
       e.g. in a tandem gap we could
       - have a trivial walk count,
       - walk too short,
       - walk too long,
       and finally walk the gap within the bounds.
    */
  
} ScaffoldWalkStatisticsT;



VA_DEF(ScaffoldWalkStatisticsT)


typedef struct
{
  ScaffoldWalkStatisticsT allScaffoldStats;
  VA_TYPE(ScaffoldWalkStatisticsT) *ScaffoldStats;
} WalkStatisticsT;



/* -------------------------------------------------- */
/* functions to manipulate the WalkStatisticsT struct */
/* -------------------------------------------------- */

void allocate_walk_statistics(WalkStatisticsT *s, int num_scaff);
// allocates ScaffoldWalkStatisticsT structs in a VA
void free_walk_statistics(WalkStatisticsT *s);
// frees the above VA

//void store_walk_statistics(WalkStatisticsT *s);

//void read_walk_statistics(WalkStatisticsT *s);
void compute_combined_walk_statistics(WalkStatisticsT* ws);
// combine statistics of all scaffolds

void output_combined_celagram_files(WalkStatisticsT* ws);
// print celagram files for all gaps

void print_all_scaffold_walk_stat_struct(char* mesg, WalkStatisticsT* ws, FILE* file, int printGaps);
// prints a readable form of a statistic struct if printGaps is true it prints
// all GapStatisticsT structs in the VA GapStats

/* ----------------------------------------------------------- */
/* functions to manipulate the ScaffoldWalkStatisticsT struct  */
/* ----------------------------------------------------------- */

void allocate_scaffold_walk_statistics(ScaffoldWalkStatisticsT *s);
// preallocates room for 100  GapStatisticsStruct in a VA

void free_scaffold_walk_statistics(ScaffoldWalkStatisticsT *s);
// frees the VA holding the GapStatisticsStruct

void init_scaffold_walk_stat_struct(ScaffoldWalkStatisticsT* ws);
// initializes all counters to zero

void print_scaffold_walk_stat_struct(ScaffoldWalkStatisticsT* ws, FILE* file, int printGaps);
// prints a readable form of a statistic struct if printGaps is true it prints
// all GapStatisticsT structs in the VA GapStats

void compute_scaffold_statistics(ScaffoldWalkStatisticsT* ws);
// combine all gap statistics of gaps in that scaffold


int read_scaffold_walk_statistics(ScaffoldWalkStatisticsT *s);
// reads the actual atomic units the GapStatisticsT of one scaffold
// returns TRUE if succesful, FALSE otherwise

void store_scaffold_walk_statistics(ScaffoldWalkStatisticsT *s);
// stores the actual atomic units the GapStatisticsT of one scaffold

/* ------------------------------------------------ */
/* misc. functions                                  */
/* ------------------------------------------------ */

void init_stat_dir(void);
// checks for a directory stats in the working directory. If not present it 
// creates it

void init_cam_dir(void);
// checks for a directory cam in the working directory. If not present it 
// creates it

void init_uoms_dir(void);
// checks for a directory uoms in the working directory. If not present it 
// creates it
#endif

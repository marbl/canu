
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
/*
  $Id: AS_CGB_breakers.h,v 1.1.1.1 2004-04-14 13:49:47 catmandew Exp $
*/
#ifndef AS_CGB_BREAKERS_H
#define AS_CGB_BREAKERS_H

#define CHECK_FOR_BAD_BREAKERS

#define STRING_LENGTH 1024

#define CHIMERA_CHUNKS    5
#define CHIMERA_OVERLAPS  2

#define SPUR_CHUNKS       3
#define SPUR_OVERLAPS     1

#define STANDARD_RANGE   40
#define STANDARD_ERATE    0.06
#define STANDARD_THRESH   1e-6

#define RELAXED_RANGE    30
#define RELAXED_ERATE     0.08
#define RELAXED_THRESH    1e-4


// for indexing into an array of chunks
typedef enum
{
  Chunk_s = 0,
  Chunk_c,
  Chunk_d,
  Chunk_a,
  Chunk_b
} ChunkIndex;

// for indexing into an array of chunk/fragment overlaps
typedef enum
{
  Ovl_sd = 0,
  Ovl_sa
} OverlapIndex;

// structure for chimera/spur
typedef struct
{
   // s = 0, d = 1, a = 2
  IntUnitigMesg     chunks[CHIMERA_CHUNKS];
  int               suffixes[CHIMERA_CHUNKS];
  // s_d = 0, s_a = 1
  UnitigOverlapMesg overlaps[CHIMERA_OVERLAPS];
} Breaker;
typedef Breaker * Breakerp;

typedef struct
{
  IntChunk_ID iaccession;
  int         chim_i;
  int         chunk_i;
} ChunkItem;
typedef ChunkItem * ChunkItemp;

typedef enum
{
  Chimera,
  Spur
} BreakerType;

typedef struct
{
  BreakerType type;
  int         num_breakers;
  Breakerp    breakers;
  ChunkItemp  indexes;
} BreakerSet;
typedef BreakerSet * BreakerSetp;


BreakerSetp ReadBreakersFile( char * filename, BreakerType type );
void FreeBreakerSet( BreakerSetp bs );
int GetUnitigData( BreakerSetp chims,
                   BreakerSetp spurs,
                   char ** cgb_files,
                   int num_cgb_files );

#endif // AS_CGB_BREAKERS_H

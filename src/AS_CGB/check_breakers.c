
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <unistd.h>
#include <getopt.h>

#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_ALN_aligners.h"

// commented out until good solution found for UIDs
//#define NEED_REAL_UID

#ifdef NEED_REAL_UID
#include "SYS_UIDcommon.h"
#include "SYS_UIDclient.h"
#endif

#define CHECK_FOR_BAD_BREAKERS

#define STRING_LENGTH 1024

#define CHIMERA_CHUNKS    5
#define CHIMERA_OVERLAPS  6

#define CRAPPY_CHUNKS     3
#define CRAPPY_OVERLAPS   3

#define STANDARD_RANGE   40
#define STANDARD_ERATE    0.06
#define STANDARD_THRESH   1e-6

#define RELAXED_RANGE    30
#define RELAXED_ERATE     0.08
#define RELAXED_THRESH    1e-4

typedef enum
{
  Silent = 0,
  Terse  = 1,
  Wordy  = 2
} Verbosity;

typedef enum
{
  Chimera,
  Crappy
} BreakerType;

cds_uint64 BaseCycle = 0;
cds_int8 BaseComplement[256] =
{
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
/* @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
   P   Q   R   S   T   U   V   W   X   Y   Z                */
  'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N',
  'N','N','N','N','A','A','N','N','N','N','N','N','N','N','N','N',
  'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N',
  'N','N','N','N','A','A','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
  'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'
};

char Bases[4] = { 'A', 'C', 'G', 'T' };

#define COMPLEMENT_BASE(a) \
      ((BaseComplement[(int)(a)] != 'N') ? \
       BaseComplement[(int)(a)] : Bases[(BaseCycle++) % 4])

/***********************************************************************/
typedef enum
{
  Chunk_s = 0,
  Chunk_c,
  Chunk_d,
  Chunk_a,
  Chunk_b
} ChunkIndex;

typedef enum
{
  Ovl_sc = 0,
  Ovl_cd,
  Ovl_sd,
  Ovl_sb,
  Ovl_ab,
  Ovl_sa
} OverlapIndex;

// structure for chimera
typedef struct
{
   // s = 0, c = 1, d = 2, a = 3, b = 4
  IntUnitigMesg     chunks[CHIMERA_CHUNKS];
  int               suffixes[CHIMERA_CHUNKS];
  // s_c = 0, c_d = 1, s_d = 2 s_b = 3, a_b = 4, s_a = 5
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

typedef struct
{
  BreakerType type;
  int         num_breakers;
  Breakerp    breakers;
  ChunkItemp  indexes;
} BreakerSet;
typedef BreakerSet * BreakerSetp;


// called
void FreeBreakerSet( BreakerSetp bs )
{
  if( bs )
  {
    if( bs->breakers )
    {
      int i;
      for( i = 0; i < bs->num_breakers; i++ )
      {
        int j;
        for( j = 0; j < CHIMERA_CHUNKS; j++ )
        {
          if( bs->breakers[i].chunks[j].f_list )
          {
            int k;
            for( k = 0; k < bs->breakers[i].chunks[j].num_frags; k++ )
            {
              if( bs->breakers[i].chunks[j].f_list[k].delta )
                free( bs->breakers[i].chunks[j].f_list[k].delta );
            }
            free( bs->breakers[i].chunks[j].f_list );
          }
        }
      }
      free( bs->breakers );
    }
    if( bs->indexes )
      free( bs->indexes );
    free( bs );
  }
}


// called
void PrintIUM( FILE * fp, IntUnitigMesg * ium )
{
  GenericMesg gen;
  gen.t = MESG_IUM;
  gen.m = ium;
  WriteProtoMesg_AS( fp, &gen );
}


// called
void PrintBreaker( FILE * fp, Breakerp c, BreakerType type, int degree )
{
  if( type == Chimera )
  {
    fprintf( fp, "       s              a              b              c              d\n" );
    fprintf( fp, "(%10" F_IIDP ") (%10" F_IIDP ",%1d) (%10" F_IIDP ",%1d) (%10" F_IIDP ",%1d) (%10" F_IIDP ",%1d)\n",
             c->chunks[Chunk_s].iaccession,
             c->chunks[Chunk_a].iaccession, c->suffixes[Chunk_a],
             c->chunks[Chunk_b].iaccession, c->suffixes[Chunk_b],
             c->chunks[Chunk_c].iaccession, c->suffixes[Chunk_c],
             c->chunks[Chunk_d].iaccession, c->suffixes[Chunk_d] );
    if( degree > 0 )
    {
      // chunks s, a, b, c, d
      PrintIUM( fp, &(c->chunks[Chunk_s]) );
      PrintIUM( fp, &(c->chunks[Chunk_a]) );
      PrintIUM( fp, &(c->chunks[Chunk_b]) );
      PrintIUM( fp, &(c->chunks[Chunk_c]) );
      PrintIUM( fp, &(c->chunks[Chunk_d]) );
    }
  }
  else
  {
    fprintf( fp, "       s              c              d\n" );
    fprintf( fp, "(%10" F_IIDP ") (%10" F_IIDP ",%1d) (%10" F_IIDP ",%1d)\n",
             c->chunks[Chunk_s].iaccession,
             c->chunks[Chunk_c].iaccession, c->suffixes[Chunk_c],
             c->chunks[Chunk_d].iaccession, c->suffixes[Chunk_d] );
    if( degree > 0 )
    {
      // chunks s, a, b, c, d
      PrintIUM( fp, &(c->chunks[Chunk_s]) );
      PrintIUM( fp, &(c->chunks[Chunk_c]) );
      PrintIUM( fp, &(c->chunks[Chunk_d]) );
    }
  }
}


// called
void PrintBreakers( FILE * fp, BreakerSetp bs, int degree )
{
  if( bs )
  {
    int num_chunks = (bs->type == Chimera) ? CHIMERA_CHUNKS : CRAPPY_CHUNKS;
    if( bs->breakers && bs->num_breakers > 0 )
    {
      int i;
      for( i = 0; i < bs->num_breakers; i++ )
      {
        fprintf( fp, "%4d: ", i );
        PrintBreaker( fp, &(bs->breakers[i]), bs->type, degree );
      }
      if( degree == 0 && bs->indexes )
      {
        for( i = 0; i < bs->num_breakers * num_chunks; i++ )
        {
          fprintf( fp, "%d: " F_IID "\n", i, bs->indexes[i].iaccession );
        }
      }
    }
    else
    {
      fprintf( fp, "No breakers to print\n" );
    }
  }
}


// called
BreakerSetp AllocateBreakerSet( int num_breakers, BreakerType type )
{
  int num_chunks = (type == Chimera) ? CHIMERA_CHUNKS : CRAPPY_CHUNKS;
  BreakerSetp bs;
  bs = (BreakerSetp) calloc( 1, sizeof( BreakerSet ) );
  if( bs == NULL )
    return NULL;

  bs->type = type;
  
  bs->breakers = (Breakerp) calloc( num_breakers, sizeof( Breaker ) );
  if( bs->breakers == NULL )
  {
    FreeBreakerSet( bs );
    return NULL;
  }
  bs->num_breakers = num_breakers;

  bs->indexes = (ChunkItemp) calloc( num_breakers * num_chunks,
                                     sizeof( ChunkItem ) );
  if( bs->indexes == NULL )
  {
    FreeBreakerSet( bs );
    return NULL;
  }
  return bs;
}


/***********************************************************************/
// called
int AllocateAndCopyString( char ** dest, char * src )
{
  if( src != NULL )
  {
    *dest = (char *) calloc( strlen( src ) + 1, sizeof( char ) );
    if( *dest == NULL )
      return 1;
    strcpy( *dest, src );
  }
  else
    *dest = NULL;
  
  return 0;
}


// called
void ReverseComplement( char * to, char * from, cds_int32 length )
{
  cds_int32 i;

  for( i = 0; i < length; i++ )
    to[length - i - 1] = COMPLEMENT_BASE( from[i] );
 
  to[length] = '\0';
}


// called
int GetNextLine( FILE * fp, char * string, int length )
{
  string[0] = '\0';
  while( fgets( string, length, fp ) )
  {
    if( string[0] != '#' )
      return 0;
  }
  return 1;
}


// called
int CountLines( FILE * fp )
{
  int num_lines = 0;
  char string[STRING_LENGTH];

  while( !GetNextLine( fp, string, STRING_LENGTH ) )
    num_lines++;

  rewind( fp );
  return num_lines;
}


// called
int GetNextBreaker( FILE * fp, Breakerp b, BreakerType type )
{
  char string[STRING_LENGTH];
  int i;
  if( type == Chimera )
  {
    if( ! GetNextLine( fp, string, STRING_LENGTH ) )
    {
      sscanf( string, F_IID " : (" F_IID ",%d) (" F_IID ",%d) (" F_IID ",%d) (" F_IID ",%d)",
              &(b->chunks[Chunk_s].iaccession),
              &(b->chunks[Chunk_a].iaccession), &(b->suffixes[Chunk_a]),
              &(b->chunks[Chunk_b].iaccession), &(b->suffixes[Chunk_b]),
              &(b->chunks[Chunk_c].iaccession), &(b->suffixes[Chunk_c]),
              &(b->chunks[Chunk_d].iaccession), &(b->suffixes[Chunk_d]) );
      string[0] = '\0';
      for( i = 0; i < CHIMERA_CHUNKS; i++ )
      {
        b->chunks[i].length = 0;
        AllocateAndCopyString( &(b->chunks[i].consensus), string );
        AllocateAndCopyString( &(b->chunks[i].quality), string );
      }
      return 0;
    }
  }
  else
  {
    if( ! GetNextLine( fp, string, STRING_LENGTH ) )
    {
      sscanf( string, "(" F_IID ",%d) (" F_IID ",%d) (" F_IID ",%d)",
              &(b->chunks[Chunk_s].iaccession), &(b->suffixes[Chunk_s]),
              &(b->chunks[Chunk_c].iaccession), &(b->suffixes[Chunk_c]),
              &(b->chunks[Chunk_d].iaccession), &(b->suffixes[Chunk_d]) );
      string[0] = '\0';
      for( i = 0; i < CRAPPY_CHUNKS; i++ )
      {
        b->chunks[i].length = 0;
        AllocateAndCopyString( &(b->chunks[i].consensus), string );
        AllocateAndCopyString( &(b->chunks[i].quality), string );
      }
      return 0;
    }
  }
  return 1;
}


// called
static int CompareChunkItems( const ChunkItem * ci1, const ChunkItem * ci2 )
{
  if( ci1->iaccession < ci2->iaccession )
    return -1;
  if( ci1->iaccession > ci2->iaccession )
    return 1;
  return 0;
}


// called
BreakerSetp ReadBreakersFile( char * filename, BreakerType type )
{
  FILE * fp;
  int num_breakers;
  BreakerSetp bs = NULL;
  int num_chunks = (type == Chimera) ? CHIMERA_CHUNKS : CRAPPY_CHUNKS;

  if( (fp = fopen( filename, "r" )) == NULL )
  {
    fprintf( stderr, "Failed to open file %s for reading.\n", filename );
    return NULL;
  }

  num_breakers = CountLines( fp );
  if( num_breakers > 0 )
  {
    int i;
    bs = AllocateBreakerSet( num_breakers, type );
    if( bs == NULL )
    {
      fprintf( stderr, "Failed to allocate set of %d breakers\n",
               num_breakers );
      fclose( fp );
      return NULL;
    }

    for( i = 0; i < bs->num_breakers; i++ )
    {
      int j, k;
      // read the basic breaker chunk IDs & suffixes
      if( GetNextBreaker( fp, &(bs->breakers[i]), bs->type ) )
      {
        fprintf( stderr, "Failed to read breaker %i\n", i + 1 );
        FreeBreakerSet( bs );
        return NULL;
      }
#ifdef CHECK_FOR_BAD_BREAKERS
      if( bs->type == Chimera &&
          (bs->breakers[i].chunks[Chunk_s].iaccession ==
           bs->breakers[i].chunks[Chunk_a].iaccession ||
           bs->breakers[i].chunks[Chunk_s].iaccession ==
           bs->breakers[i].chunks[Chunk_b].iaccession ||
           bs->breakers[i].chunks[Chunk_s].iaccession ==
           bs->breakers[i].chunks[Chunk_c].iaccession ||
           bs->breakers[i].chunks[Chunk_s].iaccession ==
           bs->breakers[i].chunks[Chunk_d].iaccession ||
           bs->breakers[i].chunks[Chunk_a].iaccession ==
           bs->breakers[i].chunks[Chunk_b].iaccession ||
           bs->breakers[i].chunks[Chunk_a].iaccession ==
           bs->breakers[i].chunks[Chunk_c].iaccession ||
           bs->breakers[i].chunks[Chunk_a].iaccession ==
           bs->breakers[i].chunks[Chunk_d].iaccession ||
           bs->breakers[i].chunks[Chunk_b].iaccession ==
           bs->breakers[i].chunks[Chunk_c].iaccession ||
           bs->breakers[i].chunks[Chunk_b].iaccession ==
           bs->breakers[i].chunks[Chunk_d].iaccession ||
           bs->breakers[i].chunks[Chunk_c].iaccession ==
           bs->breakers[i].chunks[Chunk_d].iaccession) )
      {
        fprintf( stderr, "Chimera %d is bogus: non-unique IIDs!\n", i + 1 );
        fprintf( stderr, "s = " F_IID ", a = " F_IID ", b = " F_IID ", c = " F_IID ", d = " F_IID "\n",
                 bs->breakers[i].chunks[Chunk_s].iaccession,
                 bs->breakers[i].chunks[Chunk_a].iaccession,
                 bs->breakers[i].chunks[Chunk_b].iaccession,
                 bs->breakers[i].chunks[Chunk_c].iaccession,
                 bs->breakers[i].chunks[Chunk_d].iaccession );
        fprintf( stderr,
                 "Chimera %d has been eliminated from further checking\n",
                 i + 1 );
        bs->num_breakers--;
        i--;
        continue;
      }
      else if( bs->type == Crappy &&
               (bs->breakers[i].chunks[Chunk_s].iaccession ==
                bs->breakers[i].chunks[Chunk_c].iaccession ||
                bs->breakers[i].chunks[Chunk_s].iaccession ==
                bs->breakers[i].chunks[Chunk_d].iaccession ||
                bs->breakers[i].chunks[Chunk_c].iaccession ==
                bs->breakers[i].chunks[Chunk_d].iaccession) )
      {
        fprintf( stderr, "Crappy %d is bogus: non-unique IIDs!\n", i + 1 );
        fprintf( stderr, "s = " F_IID ", c = " F_IID ", d = " F_IID "\n",
                 bs->breakers[i].chunks[Chunk_s].iaccession,
                 bs->breakers[i].chunks[Chunk_c].iaccession,
                 bs->breakers[i].chunks[Chunk_d].iaccession );
        fprintf( stderr,
                 "Crappy %d has been eliminated from further checking\n",
                 i + 1 );
        bs->num_breakers--;
        i--;
        continue;
      }
#endif
      // set up the ChunkItems for easy indexing
      for( j = i * num_chunks, k = 0; k < num_chunks; j++, k++ )
      {
        bs->indexes[j].iaccession = bs->breakers[i].chunks[k].iaccession;
        bs->indexes[j].chim_i = i;
        bs->indexes[j].chunk_i = k;
	// check_breakers.c: warning: conversion from `int' to `enum ChunkIndex'
      }
    }
    if( bs->num_breakers > 0 )
      qsort( bs->indexes,
             bs->num_breakers * num_chunks,
             sizeof( ChunkItem ),
             (int (*) (const void *, const void *)) CompareChunkItems );
  }
  return bs;
}


// called
int AllocateAndCopyIMPs( cds_int32 * dest_num, IntMultiPos ** dest_f_list,
                         cds_int32 src_num, IntMultiPos * src_f_list )
{
  int i;

  *dest_f_list = (IntMultiPos *) calloc( src_num, sizeof( IntMultiPos ) );
  if( *dest_f_list == NULL )
    return 1;
  
  for( i = 0; i < src_num; i++ )
  {
    memcpy( &((*dest_f_list)[i]), &src_f_list[i], sizeof( IntMultiPos ) );
#ifdef AS_ENABLE_SOURCE
    if( AllocateAndCopyString( &((*dest_f_list)[i].source),
                               src_f_list[i].source ) )
    {
      fprintf( stderr, "Failed to allocate IMP source memory\n" );
      return 1;
    }
#endif
    // NOTE: these will be used to hold sequence for alignments
    (*dest_f_list)[i].delta_length = 0;
    (*dest_f_list)[i].delta = NULL;
  }
  *dest_num = src_num;
  return 0;
}


// called
int PopulateIUMWithIUM( IntUnitigMesg * dest, IntUnitigMesg * src )
{
  memcpy( dest, src, sizeof( IntUnitigMesg ) );
  dest->num_frags = 0;
  dest->f_list = NULL;

#ifdef AS_ENABLE_SOURCE
  if( AllocateAndCopyString( &(dest->source), src->source ) )
  {
    fprintf( stderr, "Failed to allocate IUM source memory\n" );
    return 1;
  }
#endif
  if( AllocateAndCopyString( &(dest->consensus), src->consensus ) )
  {
    fprintf( stderr, "Failed to allocate IUM consensus memory\n" );
    return 1;
  }
  if( AllocateAndCopyString( &(dest->quality), src->quality ) )
  {
    fprintf( stderr, "Failed to allocate IUM quality memory\n" );
    return 1;
  }
  if( AllocateAndCopyIMPs( &(dest->num_frags), &(dest->f_list),
                           src->num_frags, src->f_list ) )
  {
    fprintf( stderr, "Failed to allocate IUM IMP memory\n" );
    return 1;
  }
  return 0;
}


// called
void UpdateChunkLength( IntUnitigMesg * ium )
{
  int i;
  ium->length = 0;
  for( i = 0; i < ium->num_frags; i++ )
  {
    ium->length = max( ium->f_list[i].position.bgn,
                     max( ium->f_list[i].position.end,
                          ium->length ) );
  }

  // make sure the last fragment in the f_list is at the
  // end of the chunk
  // this is needed to compute chunk overlap sizes
  if( ium->length > max( ium->f_list[ium->num_frags - 1].position.bgn,
                         ium->f_list[ium->num_frags - 1].position.end ) )
  {
    // last fragment is not at end of unitig - switch
    IntMultiPos frag1 = ium->f_list[ium->num_frags - 1];
    for( i = ium->num_frags - 1;
         ium->length > max( ium->f_list[i].position.bgn,
                            ium->f_list[i].position.end );
         i-- );
    ium->f_list[ium->num_frags - 1] = ium->f_list[i];
    ium->f_list[i] = frag1;
  }
}


// called
int GetUnitigData( BreakerSetp chims, BreakerSetp craps, char * cgb_filename )
{
  FILE * fp;
  GenericMesg       * gen;
  IntUnitigMesg     * ium;
  MesgReader reader;
  int chim_i, crap_i;

  if( (fp = fopen( cgb_filename, "r" )) == NULL )
  {
    fprintf( stderr, "Failed to open file %s for reading\n",
             cgb_filename );
    return 1;
  }
  reader = InputFileType_AS( fp );

  chim_i = 0;
  crap_i = 0;
  while( reader( fp, &gen ) != EOF )
  {
    switch( gen->t )
    {
      case MESG_IUM:
        ium = (IntUnitigMesg *) gen->m;
        // NOTE: assuming the IUMs are in order in input file
        // for chunks belonging to multiple breaker patterns, populate all
        if( chims )
        {
          for( ;
               chim_i < chims->indexes[CHIMERA_CHUNKS *
                                      chims->num_breakers - 1].iaccession &&
                 chims->indexes[chim_i].iaccession == ium->iaccession;
               chim_i++ )
          {
            if( PopulateIUMWithIUM( &(chims->breakers[chims->indexes[chim_i].chim_i].chunks[chims->indexes[chim_i].chunk_i]), ium ) )
            {
              fprintf( stderr, "Failed to populate IUM with IUM\n" );
              fclose( fp );
              return 1;
            }
            else
            {
              UpdateChunkLength( &(chims->breakers[chims->indexes[chim_i].chim_i].chunks[chims->indexes[chim_i].chunk_i]) );
            }
          }
        }
        if( craps )
        {
          for( ;
               crap_i < craps->indexes[CRAPPY_CHUNKS *
                                      craps->num_breakers - 1].iaccession &&
                 craps->indexes[crap_i].iaccession == ium->iaccession;
               crap_i++ )
          {
            if( PopulateIUMWithIUM( &(craps->breakers[craps->indexes[crap_i].chim_i].chunks[craps->indexes[crap_i].chunk_i]), ium ) )
            {
              fprintf( stderr, "Failed to populate IUM with IUM\n" );
              fclose( fp );
              return 1;
            }
            else
            {
              UpdateChunkLength( &(craps->breakers[craps->indexes[crap_i].chim_i].chunks[craps->indexes[crap_i].chunk_i]) );
            }
          }
        }
        break;
      default:
        break;
    }
  }
  fclose( fp );
#ifndef CHECK_FOR_BAD_BREAKERS
  // here, two chunks in one breaker may have the same IID
  // and only one will have been populated, so populate others that match
  if( chims )
  {
    int l, m, n;
    for( chim_i = 0; chim_i < chims->num_breakers; chim_i++ )
    {
      // loop over all chunks in the current breaker (not n - 1)
      for( l = 0; l < CHIMERA_CHUNKS; l++ )
      {
        // loop over all other chunks in the current breaker
        for( m = l + 1; m < CHIMERA_CHUNKS; m++ )
        {
          if( chims->breakers[chim_i].chunks[l].iaccession ==
              chims->breakers[chim_i].chunks[m].iaccession )
          {
            // which one got populated? -
            // it's possible that neither are populated
            if( chims->breakers[chim_i].chunks[m].consensus == NULL )
            {
              if( chims->breakers[chim_i].chunks[l].consensus != NULL )
                PopulateIUMWithIUM( &(chims->breakers[chim_i].chunks[m]),
                                    &(chims->breakers[chim_i].chunks[l]) );
            }
            else
            {
              PopulateIUMWithIUM( &(chims->breakers[chim_i].chunks[l]),
                                  &(chims->breakers[chim_i].chunks[m]) );
            }
          }
        }
        // loop over all other breakers
        for( n = chim_i + 1; n < chims->num_breakers; n++ )
        {
          // loop over all chunks in the breaker
          for( m = 0; m < CHIMERA_CHUNKS; m++ )
          {
            if( chims->breakers[chim_i].chunks[l].iaccession ==
                chims->breakers[n].chunks[m].iaccession )
            {
              // which one got populated? -
              // it's possible that neither are populated
              if( chims->breakers[n].chunks[m].consensus == NULL )
              {
                if( chims->breakers[chim_i].chunks[l].consensus != NULL )
                  PopulateIUMWithIUM( &(chims->breakers[n].chunks[m]),
                                      &(chims->breakers[chim_i].chunks[l]) );
              }
              else
              {
                PopulateIUMWithIUM( &(chims->breakers[chim_i].chunks[l]),
                                    &(chims->breakers[n].chunks[m]) );
              }
            }
          }
        }
      }
    }
  }
#endif
  return 0;
}


// called
int PopulateFragmentSequence( IntMultiPos * f, FragStoreHandle fs )
{
  ReadStructp rs;
  char temp_seq[AS_READ_MAX_LEN];
  char temp_qvs[AS_READ_MAX_LEN];
  cds_uint32 bgn, end;

  if( (rs = new_ReadStruct()) == NULL )
  {
    fprintf( stderr, "Failed to get new ReadStruct\n" );
    return 1;
  }
    
  if( getFragStore( fs, f->ident, FRAG_S_ALL, rs ) )
  {
    fprintf( stderr, "Failed to get fragment " F_IID " sequence\n", f->ident );
    return 1;
  }
  
  getClearRegion_ReadStruct( rs, &bgn, &end, READSTRUCT_LATEST );
  getSequence_ReadStruct( rs, temp_seq, temp_qvs, AS_READ_MAX_LEN );
  f->delta_length = end - bgn;
  
  if( (f->delta = (cds_int32 *) calloc( f->delta_length + 1,
                                        sizeof( char ) )) == NULL )
  {
    fprintf( stderr, "Failed to allocate frag " F_IID " sequence of %d bases\n",
             f->ident, f->delta_length );
    return 1;
  }
  strncpy( (char *) f->delta, &(temp_seq[bgn]), f->delta_length );
  
  delete_ReadStruct( rs );
  return 0;
}


// called
#define BASES_PER_LINE  60
void PrintFragmentLine( FILE * fp, char * seq, CDS_COORD_t length,
                        CDS_COORD_t bgn, CDS_COORD_t end, int rev )
{
  CDS_COORD_t i;
  
  if( rev )
    fprintf( fp, "%5" F_COORDP ": ", length - bgn );
  else
    fprintf( fp, "%5" F_COORDP ": ", bgn );

  for( i = bgn; i < min( bgn + BASES_PER_LINE, end ); i++ )
  {
    fprintf( fp, "%c", seq[i] );
  }
  fprintf( fp, "\n" );
}


// called
OverlapMesg * OverlapFragments( FILE * fp, Verbosity verbose,
                                IntMultiPos * f1,
                                IntMultiPos * f2,
                                ChunkOrientationType orient )
{
  InternalFragMesg ifg1;
  InternalFragMesg ifg2;
  OverlapMesg    * ovl;
  char str1[AS_READ_MAX_LEN];
  char str2[AS_READ_MAX_LEN];
  int where;
  int dp_orient;

  // populate the ifg messages for DP_Compare_AS()
  ifg1.iaccession = f1->ident;
  ifg1.clear_rng.bgn = 0;
  ifg1.clear_rng.end = f1->delta_length;

  strcpy( str1, (char *) f1->delta );
  ifg1.sequence = str1;
  ifg1.quality = NULL;

  ifg2.iaccession = f2->ident;
  ifg2.clear_rng.bgn = 0;
  ifg2.clear_rng.end = f2->delta_length;

  strcpy( str2, (char *) f2->delta );
  ifg2.sequence = str2;
  ifg2.quality = NULL;

  if( orient == AB_AB || orient == BA_BA )
    dp_orient = 0;
  else
    dp_orient = 1;
  
  ovl = DP_Compare_AS( &ifg1, &ifg2,
                       -f2->delta_length, f1->delta_length,
                       dp_orient,
                       STANDARD_ERATE, STANDARD_THRESH, STANDARD_RANGE,
                       AS_FIND_OVERLAP, &where );
  if( ovl && verbose == Wordy )
  {
    fprintf( fp, "Overlap found with standard parameters:\n"
             "\trange = %d, erate = %f, thresh = %f\n",
             STANDARD_RANGE, STANDARD_ERATE, STANDARD_THRESH );
  }
  else
  {
    ovl = DP_Compare_AS( &ifg1, &ifg2,
                         -f2->delta_length, f1->delta_length,
                         dp_orient,
                         RELAXED_ERATE, RELAXED_THRESH, RELAXED_RANGE,
                         AS_FIND_OVERLAP, &where );
    if( ovl && verbose == Wordy )
      fprintf( fp, "Overlap found with relaxed parameters:\n"
               "\trange = %d, erate = %f, thresh = %f\n",
               RELAXED_RANGE, RELAXED_ERATE, RELAXED_THRESH );
  }

  if( ovl )
  {
    ovl->delta = (signed char *) "\0";
    if( verbose == Wordy )
    {
      Print_Overlap_AS( fp, &ifg1, &ifg2, ovl );
      fprintf( fp, "\n" );
    }
  }
  return ovl;
}


// called
int WriteSequence( FILE * fp, char * seq, int length, int bases_per_line )
{
  int i;

  fprintf( fp, "\t" );
  for( i = 0; i < length && seq[i] != '\0'; i++ )
  {
    if( i > 0 && i % bases_per_line == 0 )
      fprintf( fp, "\n\t" );
    fprintf( fp, "%c", seq[i] );
  }
  if( i % bases_per_line )
    fprintf( fp, "\n" );
  return 0;
}


// called
cds_int32 CheckFragmentOverlap( FILE * fp_log, Verbosity verbose,
                                FILE * fp_ovl, MesgWriter writer,
                                IntMultiPos * f1,
                                IntMultiPos * f2,
                                ChunkOrientationType orient,
                                FragStoreHandle fs )
{
  ChunkOrientationType frag_orient = XX_XX;
  
  // get the fragment sequences, if necessary
  if( f1->delta == NULL )
  {
    if( PopulateFragmentSequence( f1, fs ) )
    {
      fprintf( stderr,
               "Failed to populate fragment " F_IID " sequence\n",
               f1->ident );
      return 0;
    }
  }
  if( f2->delta == NULL )
  {
    if( PopulateFragmentSequence( f2, fs ) )
    {
      fprintf( stderr,
               "Failed to populate fragment " F_IID " sequence\n",
               f2->ident );
      return 0;
    }
  }

  if( f1->ident == f2->ident && verbose != Silent )
  {
    fprintf( fp_log, "Overlap of fragment " F_IID " & itself!\n\n", f1->ident );
    return f1->delta_length;
  }
  
  // do an alignment
  {
    OverlapMesg * ovl = NULL;

    // handle orientation of fragments within oriented chunks
    switch( orient )
    {
      case AB_AB:
        if( f1->position.end > f1->position.bgn )
        {
          if( f2->position.end > f2->position.bgn )
          {
            frag_orient = AB_AB;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
          else
          {
            frag_orient = AB_BA;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
        }
        else
        {
          if( f2->position.end > f2->position.bgn )
          {
            frag_orient = BA_AB;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
          else
          {
            frag_orient = BA_BA;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
        }
        break;
      case BA_BA:
        if( f1->position.end > f1->position.bgn )
        {
          if( f2->position.end > f2->position.bgn )
          {
            frag_orient = BA_BA;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
          else
          {
            frag_orient = BA_AB;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
        }
        else
        {
          if( f2->position.end > f2->position.bgn )
          {
            frag_orient = AB_BA;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
          else
          {
            frag_orient = AB_AB;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
        }
        break;
      case BA_AB:
        if( f1->position.end > f1->position.bgn )
        {
          if( f2->position.end > f2->position.bgn )
          {
            frag_orient = BA_AB;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
          else
          {
            frag_orient = BA_BA;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
        }
        else
        {
          if( f2->position.end > f2->position.bgn )
          {
            frag_orient = AB_AB;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
          else
          {
            frag_orient = AB_BA;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
        }
        break;
      case AB_BA:
        if( f1->position.end > f1->position.bgn )
        {
          if( f2->position.end > f2->position.bgn )
          {
            frag_orient = AB_BA;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
          else
          {
            frag_orient = AB_AB;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
        }
        else
        {
          if( f2->position.end > f2->position.bgn )
          {
            frag_orient = BA_BA;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
          else
          {
            frag_orient = BA_AB;
            ovl = OverlapFragments( fp_log, verbose, f1, f2, frag_orient );
          }
        }
        break;
      default:
        fprintf( stderr, "Undefined fragment overlap orientation: %d\n",
                 orient );
        break;
    }
    if( ovl )
    {
      GenericMesg gen;
      gen.t = MESG_OVL;
      gen.m = ovl;
      
      if( verbose == Wordy )
        WriteProtoMesg_AS( fp_log, &gen );
      
      // need to switch around the overlap orientation
      switch( frag_orient )
      {
        case AB_AB:
          ovl->orientation = AS_NORMAL;
          break;
        case AB_BA:
          ovl->orientation = AS_INNIE;
          break;
        case BA_AB:
          ovl->orientation = AS_OUTTIE;
          break;
        case BA_BA:
          ovl->orientation = AS_ANTI;
          break;
        default:
          break;
      }
      if( fp_ovl )
        writer( fp_ovl, &gen );

      return (f1->delta_length + f2->delta_length - ovl->ahg - ovl->bhg) / 2;
    }
    else if( verbose != Silent )
    {
      char seq1[AS_READ_MAX_LEN];
      char seq2[AS_READ_MAX_LEN];
      
      // no overlap!
      switch( frag_orient )
      {
        case AB_AB:
          strcpy( seq1, (char *) f1->delta );
          strcpy( seq2, (char *) f2->delta );
          fprintf( fp_log, "AB_AB " );
          break;
        case AB_BA:
          strcpy( seq1, (char *) f1->delta );
          ReverseComplement( seq2, (char *) f2->delta, f2->delta_length );
          fprintf( fp_log, "AB_BA " );
          break;
        case BA_AB:
          ReverseComplement( seq1, (char *) f1->delta, f1->delta_length );
          strcpy( seq2, (char *) f2->delta );
          fprintf( fp_log, "BA_AB " );
          break;
        case BA_BA:
          ReverseComplement( seq1, (char *) f1->delta, f1->delta_length );
          ReverseComplement( seq2, (char *) f2->delta, f2->delta_length );
          fprintf( fp_log, "BA_BA " );
          break;
        default:
          return 0;
          break;
      }

      fprintf( fp_log, "overlap between fragments " F_IID " & " F_IID " not found!\n",
               f1->ident, f2->ident );
      fprintf( fp_log, "sequences complemented as needed to 'expose' possible overlap\n\n" );
      fprintf( fp_log, "%10" F_IIDP ":\n", f1->ident );
      WriteSequence( fp_log, seq1, f1->delta_length, 60 );
      fprintf( fp_log, "%10" F_IIDP ":\n", f2->ident );
      WriteSequence( fp_log, seq2, f2->delta_length, 60 );
    }
  }
  return 0;
}


// called
int CheckBreakerChunkOverlap( FILE * fp_log, Verbosity verbose,
                              FILE * fp_ovl, MesgWriter writer,
                              Breakerp b, OverlapIndex oi,
                              ChunkIndex ci1, ChunkIndex ci2,
                              FragStoreHandle fs )
{
  int ci1_frag;
  int ci2_frag;
  
  // set the chunk IIDs
  b->overlaps[oi].chunk1 = b->chunks[ci1].iaccession;
  b->overlaps[oi].chunk2 = b->chunks[ci2].iaccession;
  
  // set the overlap orientation & indices of relevant fragments
  if( b->suffixes[ci1] )
  {
    ci1_frag = b->chunks[ci1].num_frags - 1;
    if( b->suffixes[ci2] )
    {
      b->overlaps[oi].orient = AB_BA;
      ci2_frag = b->chunks[ci2].num_frags - 1;
    }
    else
    {
      b->overlaps[oi].orient = AB_AB;
      ci2_frag = 0;
    }
  }
  else
  {
    ci1_frag = 0;
    if( b->suffixes[ci2] )
    {
      b->overlaps[oi].orient = BA_BA;
      ci2_frag = b->chunks[ci2].num_frags - 1;
    }
    else
    {
      b->overlaps[oi].orient = BA_AB;
      ci2_frag = 0;
    }
  }
  
  // look for an overlap
  b->overlaps[oi].best_overlap_length =
    CheckFragmentOverlap( fp_log, verbose, fp_ovl, writer,
                          &(b->chunks[ci1].f_list[ci1_frag]),
                          &(b->chunks[ci2].f_list[ci2_frag]),
                          b->overlaps[oi].orient,
                          fs );
  return 0;
}


// called
// since not all chunk overlaps were represented by UOMs
int CheckBreakerOverlaps( FILE * fp_log, Verbosity verbose,
                          FILE * fp_ovl, MesgWriter writer,
                          Breakerp b, BreakerType type,
                          FragStoreHandle fs )
{
  int s_suffix = b->suffixes[Chunk_s];
  // suffix of s is relevant in different ways for chimeras vs. crappies
  
  // say which breaker this is
  if( verbose != Silent )
  {
    fprintf( fp_log, "*************************************\n" );
    fprintf( fp_log, "Overlaps for %s:\n",
             (type == Chimera) ? "chimera" : "crappy" );
    PrintBreaker( fp_log, b, type, 0 );
  }
/*  
  // s-c overlap
  if( verbose != Silent )
    fprintf( fp_log, "\n########## s-c overlap #########\n" );
  if( type == Chimera )
    b->suffixes[Chunk_s] = 0;
  if( CheckBreakerChunkOverlap( fp_log, verbose, NULL, writer,
                                b, Ovl_sc, Chunk_s, Chunk_c, fs ) )
  {
    fprintf( stderr, "Failed to check breaker overlap\n" );
    return 1;
  }
    
  // c-d overlap
  if( verbose != Silent )
    fprintf( fp_log, "\n########## c-d overlap #########\n" );
  if( CheckBreakerChunkOverlap( fp_log, verbose, NULL, writer,
                                b, Ovl_cd, Chunk_c, Chunk_d, fs ) )
  {
    fprintf( stderr, "Failed to check breaker overlap\n" );
    return 1;
  }
*/  
  // s-d overlap
  if( verbose != Silent )
    fprintf( fp_log, "\n########## Possible s-d overlap #########\n" );
  if( type == Chimera )
    b->suffixes[Chunk_s] = 1;
  else
    b->suffixes[Chunk_s] = 1 - s_suffix;
  if( CheckBreakerChunkOverlap( fp_log, verbose, fp_ovl, writer,
                                b, Ovl_sd, Chunk_s, Chunk_d, fs ) )
  {
    fprintf( stderr, "Failed to check breaker overlap\n" );
    b->suffixes[Chunk_s] = s_suffix;
    return 1;
  }
  b->suffixes[Chunk_s] = s_suffix;

  if( type == Chimera )
  {
    // s-b overlap
    if( verbose != Silent )
      fprintf( fp_log, "\n########## s-b overlap #########\n" );
    b->suffixes[Chunk_s] = 1;
    if( CheckBreakerChunkOverlap( fp_log, verbose, NULL, writer,
                                  b, Ovl_sb, Chunk_s, Chunk_b, fs ) )
    {
      fprintf( stderr, "Failed to check breaker overlap\n" );
      return 1;
    }
    
    // a-b overlap
    if( verbose != Silent )
      fprintf( fp_log, "\n########## a-b overlap #########\n" );
    if( CheckBreakerChunkOverlap( fp_log, verbose, NULL, writer,
                                  b, Ovl_ab, Chunk_a, Chunk_b, fs ) )
    {
      fprintf( stderr, "Failed to check breaker overlap\n" );
      return 1;
    }
    
    // s-a overlap
    if( verbose != Silent )
      fprintf( fp_log, "\n########## Possible s-a overlap #########\n" );
    b->suffixes[Chunk_s] = 0;
    if( CheckBreakerChunkOverlap( fp_log, verbose, fp_ovl, writer,
                                  b, Ovl_sa, Chunk_s, Chunk_a, fs ) )
    {
      fprintf( stderr, "Failed to check breaker overlap\n" );
      return 1;
    }
  }
    
  return 0;
}


// called
int ShowAllBreakerOverlaps( FILE * fp_log, Verbosity verbose,
                            FILE * fp_ovl, MesgWriter writer,
                            BreakerSetp bs, char * frg_store_name )
{
  FragStoreHandle fstore;
  int bi;

  // open the fragment store to get sequences in lower functions
  fstore = openFragStore( frg_store_name, "r" );
  if( fstore == NULLSTOREHANDLE )
  {
    fprintf( stderr, "Failed to open frag store %s\n", frg_store_name );
    return 1;
  }

  // loop over all breakers to show the overlaps
  for( bi = 0; bi < bs->num_breakers; bi++ )
  {
    // some UOMs won't have been in input file. If so,
    // reverse engineer the overlaps
    if( CheckBreakerOverlaps( fp_log, verbose, fp_ovl, writer,
                              &(bs->breakers[bi]), bs->type, fstore ) )
    {
      fprintf( stderr, "Failed to check breaker overlaps\n" );
      closeFragStore( fstore );
      return 1;
    }
  }
  closeFragStore( fstore );
  return 0;
}

/* Function:
     Usage
   Description:
     prints a message, usage instructions, and exits with 1
   Return Value:
     none
   Parameters:
     char * program_name:  obvious
     char * message:       message to print out (w/o \n)
*/
void Usage( char * program_name, char * message )
{
  if( message != NULL )
    fprintf( stderr, "%s: %s\n", program_name, message );
  
  fprintf( stderr,
           "Usage: %s [-h chimeras_file] [-r crappies_file] [-c cgb_file]\n"
           "\t[-s frag_store] [-l log_file] [-v level] [-o ovl_file]\n"
           "where:\n"
           "-h chimeras_file   input chimeras file (optional)\n"
           "-r crappies_file   input crappies file (optional)\n"
           "                     one or both of -h & -r must be specified\n"
           "-c cgb_file        input cgb file (required)\n"
           "-s frag_store      input frag store path (required\n"
           "-l log_file        output log file (optional)\n"
           "                     default is stdout\n"
           "-v level           verbose level\n"
           "                     0 for none\n"
           "                     1 for only missing overlaps\n"
           "                     2 for all overlaps\n"
           "-o ovl_file        output ovl file for found overlaps(optional)\n"
           "-P                 specifies ASCII ovl_file\n",
           program_name );
  exit( 1 );
}

typedef struct
{
  char * program_name;
  char * version;
  char * chims_file;
  char * craps_file;
  char * cgb_file;
  char * fstore_path;
  char * log_file;
  FILE * fp_log;
  char * ovl_file;
  FILE * fp_ovl;
  Verbosity  verbose;
  MesgWriter writer;
  BreakerSetp chims;
  BreakerSetp craps;
  int realUID;
  char *euidServerName;
} CheckGlobals;
typedef CheckGlobals * CheckGlobalsp;

// called
void InitializeGlobals( CheckGlobalsp globals, char * program_name )
{
  globals->program_name = program_name;
  globals->version = "$Revision: 1.2 $";
  globals->chims_file = NULL;
  globals->craps_file = NULL;
  globals->cgb_file = NULL;
  globals->fstore_path = NULL;
  globals->log_file = NULL;
  globals->fp_log = stdout;
  globals->ovl_file = NULL;
  globals->fp_ovl = NULL;
  globals->verbose = Silent;
  globals->writer = WriteBinaryMesg_AS;
  globals->chims = NULL;
  globals->craps = NULL;
  globals->realUID = 0;
  globals->euidServerName = NULL;
}

// called
void ParseCommandLine( CheckGlobalsp globals, int argc, char ** argv )
{
  int ch, errflg = 0;
    
  optarg = NULL;
  /* check_breakers
           [-h chimeras]
           [-r crappies]
           [-c cgb]
           [-s store]
           [-l log]
           [-v level]
           [-o ovl]
           [-P]
  */
  while( !errflg && ((ch = getopt( argc, argv, "Ph:r:c:s:l:o:v:" )) != EOF) )
  {
    switch( ch )
    {
      case 'h':
        globals->chims_file = optarg;
        break;
      case 'r':
        globals->craps_file = optarg;
        break;
      case 'c':
        globals->cgb_file = optarg;
        break;
      case 's':
        globals->fstore_path = optarg;
        break;
      case 'l':
        globals->log_file = optarg;
        globals->fp_log = fopen( globals->log_file, "w" );
        if( globals->fp_log == NULL )
        {
          fprintf( stderr, "Failed to open log file %s for writing\n",
                   globals->log_file );
          fprintf( stderr, "  redirecting to stdout\n" );
          globals->fp_log = stdout;
        }
        break;
      case 'o':
        globals->ovl_file = optarg;
        globals->fp_ovl = fopen( globals->ovl_file, "w" );
        if( globals->fp_ovl == NULL )
        {
          fprintf( stderr, "Failed to open ovl file %s for writing\n",
                   globals->ovl_file );
          fprintf( stderr, "Aborting.\n" );
          exit( 1 );
        }
        break;
      case 'v':
        globals->verbose = (Verbosity) atoi( optarg );
        break;
      case 'P':
        globals->writer = WriteProtoMesg_AS;
        break;
      default:
        errflg++;
        break;
    }
  }

  if( errflg != 0 )
    Usage( globals->program_name, "bad option." );

  if( globals->cgb_file == NULL )
    Usage( globals->program_name, "need an input cgb file\n" );

  if( globals->fstore_path == NULL )
    Usage( globals->program_name, "need an input fragment store\n" );

  if( globals->chims_file == NULL &&
      globals->craps_file == NULL )
    Usage( globals->program_name,
           "need either a chimeras or crappies file." );
}


// called
void CleanUp( CheckGlobalsp globals )
{
  if( globals->fp_log )
    fclose( globals->fp_log );
  if( globals->fp_ovl )
    fclose( globals->fp_ovl );
  
  FreeBreakerSet( globals->chims );
  FreeBreakerSet( globals->craps );
}


// called
void InitializeOVLFile( CheckGlobalsp globals )
{
  GenericMesg gen;
  BatchMesg   bat;
  AuditMesg   adt;
  AuditLine   adl;
  
#ifdef NEED_REAL_UID
  cds_uint64  max_block_size;
  cds_int32   uid_status;
  cds_uint64  uid_interval[4];

  // set up the UID system to get 1 uid
  uid_status = SYS_UIDgetMaxUIDSize( &max_block_size );
  SYS_UIDsetUIDSize( 1 );
  uid_status = SYS_UIDgetNewUIDInterval( uid_interval );
  SYS_UIDgetNextUID( &bat.eaccession );
#endif

  // BAT messages
  bat.name       = globals->program_name;
  bat.created    = time(0);
  fprintf( stderr, "Batch UID is " F_U64 "\n", bat.eaccession );
  bat.comment    = globals->cgb_file;
  gen.t = MESG_BAT;
  gen.m = &bat;
  globals->writer( globals->fp_ovl, &gen );

  // ADT messages
  adt.list = NULL;
  AppendAuditLine_AS( &adt, &adl, time(0),
                      globals->program_name,
                      globals->version, NULL );
  gen.t = MESG_ADT;
  gen.m = &adt;
  globals->writer( globals->fp_ovl, &gen );
}


int main( int argc, char ** argv )
{
  CheckGlobals globals;

  InitializeGlobals( &globals, argv[0] );

  ParseCommandLine( &globals, argc, argv );

  if( globals.chims_file != NULL )
  {
    globals.chims = ReadBreakersFile( globals.chims_file, Chimera );
    if( globals.chims == NULL )
    {
      fprintf( stderr, "Failed to get chimera data from %s\n",
               globals.chims_file );
    }
    // PrintBreakers( globals.fp_log, globals.chims, 0 );
  }

  if( globals.craps_file != NULL )
  {
    globals.craps = ReadBreakersFile( globals.craps_file, Crappy );
    if( globals.craps == NULL )
    {
      fprintf( stderr, "Failed to get crappy data from %s\n",
               globals.craps_file );
    }
    // PrintCrappies( globals.fp_log, globals.craps, 0 );
  }

  if( globals.chims == NULL && globals.craps == NULL )
    return 1;
  
  if( GetUnitigData( globals.chims, globals.craps, globals.cgb_file ) )
  {
    fprintf( stderr, "Failed to get chimera & crappy unitig data\n" );
    return 1;
  }
  // PrintBreakers( globals.fp_log, globals.chims, 1 );
  // PrintBreakers( globals.fp_log, globals.craps, 1 );

  // write the start of the ovl file, if necessary
  if( globals.fp_ovl )
  {
    InitializeOVLFile( &globals );
  }
  
  if( globals.chims )
    if( ShowAllBreakerOverlaps( globals.fp_log, globals.verbose,
                                globals.fp_ovl, globals.writer,
                                globals.chims, globals.fstore_path ) )
      fprintf( stderr, "Failed to show chimeric overlaps\n" );

  if( globals.craps )
  {
    if( ShowAllBreakerOverlaps( globals.fp_log, globals.verbose,
                                globals.fp_ovl, globals.writer,
                                globals.craps, globals.fstore_path ) )
      fprintf( stderr, "Failed to show crappy overlaps\n" );
  }
  else
    fprintf( globals.fp_log, "No crappy fragments to show\n" );

  CleanUp( &globals );
  return 0;
}

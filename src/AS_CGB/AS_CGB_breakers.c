
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
static char CM_ID[] = "$Id: AS_CGB_breakers.c,v 1.4 2005-03-22 19:48:24 jason_miller Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>

#include "AS_global.h"
#include "AS_CGB_breakers.h"


static int GetNextLine( FILE * fp, char * string, int length )
{
  string[0] = '\0';
  while( fgets( string, length, fp ) )
  {
    if( string[0] != '#' )
      return 0;
  }
  return 1;
}


static int CountLines( FILE * fp )
{
  int num_lines = 0;
  char string[STRING_LENGTH];

  while( !GetNextLine( fp, string, STRING_LENGTH ) )
    num_lines++;

  rewind( fp );
  return num_lines;
}


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
# ifdef AS_ENABLE_SOURCE
            if( bs->breakers[i].chunks[j].f_list[k].source );
              free( bs->breakers[i].chunks[j].f_list[k].source );
#endif
              if( bs->breakers[i].chunks[j].f_list[k].delta )
                free( bs->breakers[i].chunks[j].f_list[k].delta );
            }
# ifdef AS_ENABLE_SOURCE
            if( bs->breakers[i].chunks[j].source )
              free( bs->breakers[i].chunks[j].source );
#endif
            if( bs->breakers[i].chunks[j].consensus )
              free( bs->breakers[i].chunks[j].consensus );
            if( bs->breakers[i].chunks[j].quality )
              free( bs->breakers[i].chunks[j].quality );
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


BreakerSetp AllocateBreakerSet( int num_breakers, BreakerType type )
{
  int num_chunks = (type == Chimera) ? CHIMERA_CHUNKS : SPUR_CHUNKS;
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


static int AllocateAndCopyString( char ** dest, char * src )
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


static int GetNextBreaker( FILE * fp, Breakerp b, BreakerType type )
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
      for( i = 0; i < SPUR_CHUNKS; i++ )
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


static int CompareChunkItems( const ChunkItem * ci1, const ChunkItem * ci2 )
{
  if( ci1->iaccession < ci2->iaccession )
    return -1;
  if( ci1->iaccession > ci2->iaccession )
    return 1;
  return 0;
}


BreakerSetp ReadBreakersFile( char * filename, BreakerType type )
{
  FILE * fp;
  int num_breakers;
  BreakerSetp bs = NULL;
  int num_chunks = (type == Chimera) ? CHIMERA_CHUNKS : SPUR_CHUNKS;

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
      else if( bs->type == Spur &&
               (bs->breakers[i].chunks[Chunk_s].iaccession ==
                bs->breakers[i].chunks[Chunk_c].iaccession ||
                bs->breakers[i].chunks[Chunk_s].iaccession ==
                bs->breakers[i].chunks[Chunk_d].iaccession ||
                bs->breakers[i].chunks[Chunk_c].iaccession ==
                bs->breakers[i].chunks[Chunk_d].iaccession) )
      {
        fprintf( stderr, "Spur %d is bogus: non-unique IIDs!\n", i + 1 );
        fprintf( stderr, "s = " F_IID ", c = " F_IID ", d = " F_IID "\n",
                 bs->breakers[i].chunks[Chunk_s].iaccession,
                 bs->breakers[i].chunks[Chunk_c].iaccession,
                 bs->breakers[i].chunks[Chunk_d].iaccession );
        fprintf( stderr,
                 "Spur %d has been eliminated from further checking\n",
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
	// AS_CGB_breakers.c: warning: conversion from `int' to `enum ChunkIndex'
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


static int AllocateAndCopyIMPs( cds_int32 * dest_num, IntMultiPos ** dest_f_list,
                                cds_int32 src_num, IntMultiPos * src_f_list )
{
  int i;

  *dest_f_list = (IntMultiPos *) calloc( src_num, sizeof( IntMultiPos ) );
  if( *dest_f_list == NULL )
    return 1;
  
  for( i = 0; i < src_num; i++ )
  {
    memset( &((*dest_f_list)[i]), 0, sizeof( IntMultiPos ) );
    (*dest_f_list)[i].type = src_f_list[i].type;
    (*dest_f_list)[i].ident = src_f_list[i].ident;
    (*dest_f_list)[i].contained = src_f_list[i].contained;
    (*dest_f_list)[i].position = src_f_list[i].position;
    (*dest_f_list)[i].delta_length = 0;
    (*dest_f_list)[i].delta = NULL;
    
#ifdef AS_ENABLE_SOURCE
    (*dest_f_list)[i].source = NULL;
    if( AllocateAndCopyString( &((*dest_f_list)[i].source),
                               src_f_list[i].source ) )
    {
      fprintf( stderr, "Failed to allocate IMP source memory\n" );
      return 1;
    }
#endif
    // NOTE: these will be used to hold sequence for alignments
  }
  *dest_num = src_num;
  return 0;
}


static int PopulateIUMWithIUM( IntUnitigMesg * dest, IntUnitigMesg * src )
{
  memset( dest, 0, sizeof( IntUnitigMesg ) );
  dest->iaccession = src->iaccession;
  dest->coverage_stat = src->coverage_stat;
  dest->status = src->status;
  dest->a_branch_point = src->a_branch_point;
  dest->b_branch_point = src->b_branch_point;
  dest->length = src->length;
  dest->forced = src->forced;
  dest->num_frags = 0;
  dest->f_list = NULL;

#ifdef AS_ENABLE_SOURCE
  dest->source = NULL;
  if( AllocateAndCopyString( &(dest->source), src->source ) )
  {
    fprintf( stderr, "Failed to allocate IUM source memory\n" );
    return 1;
  }
#endif
  if( dest->consensus != NULL )
    fprintf( stderr, "The impossible has happened.\n" );
  dest->consensus = NULL;
  if( AllocateAndCopyString( &(dest->consensus), src->consensus ) )
  {
    fprintf( stderr, "Failed to allocate IUM consensus memory\n" );
    return 1;
  }
  if( dest->quality != NULL )
    fprintf( stderr, "The impossible has happened.\n" );
  dest->quality = NULL;
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


int GetUnitigData( BreakerSetp chims,
                   BreakerSetp spurs,
                   char ** cgb_files,
                   int num_cgb_files )
{
  FILE * fp;
  GenericMesg       * gen;
  IntUnitigMesg     * ium;
  MesgReader reader;
  int chim_i, spur_i;
  int file_i;
  cds_uint32 num_chim_chunks_read = 0;
  cds_uint32 num_spur_chunks_read = 0;

  if( chims && spurs )
    fprintf( stderr, "Getting chimera- & spur-related unitigs\n" );
  else if( chims )
    fprintf( stderr, "Getting chimera-related unitigs\n" );
  else if( spurs )
    fprintf( stderr, "Getting spur-related unitigs\n" );
  else
    return 0;

  chim_i = 0;
  spur_i = 0;
  for( file_i = 0; file_i < num_cgb_files; file_i++ )
  {
    if( (fp = fopen( cgb_files[file_i], "r" )) == NULL )
    {
      fprintf( stderr, "Failed to open file %s for reading\n",
               cgb_files[file_i] );
      return 1;
    }
    reader = InputFileType_AS( fp );
    fprintf( stderr, "Reading cgb file %s\n", cgb_files[file_i] );
    
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
                 chim_i < CHIMERA_CHUNKS * chims->num_breakers &&
                   chims->indexes[chim_i].iaccession == ium->iaccession;
                 chim_i++ )
            {
              if( PopulateIUMWithIUM( &(chims->breakers[chims->indexes[chim_i].chim_i].chunks[chims->indexes[chim_i].chunk_i]), ium ) )
              {
                fprintf( stderr, "Failed to populate IUM with IUM\n" );
                fclose( fp );
                return 1;
              }
              num_chim_chunks_read++;
            }
          }
          if( spurs )
          {
            if( ium->iaccession > spurs->indexes[spur_i].iaccession &&
                spur_i < SPUR_CHUNKS * spurs->num_breakers )
            {
              fprintf( stderr, "Skipped an IUM?\n" );
              fprintf( stderr, "%10u\t%10u\n", ium->iaccession,
                       spurs->indexes[spur_i].iaccession );
              fclose( fp );
              return 1;
            }
            for( ;
                 spur_i < SPUR_CHUNKS * spurs->num_breakers &&
                   spurs->indexes[spur_i].iaccession == ium->iaccession;
                 spur_i++ )
            {
              if( PopulateIUMWithIUM( &(spurs->breakers[spurs->indexes[spur_i].chim_i].chunks[spurs->indexes[spur_i].chunk_i]), ium ) )
              {
                fprintf( stderr, "Failed to populate IUM with IUM\n" );
                fclose( fp );
                return 1;
              }
              num_spur_chunks_read++;
            }
          }
          break;
        default:
          break;
      }
    }
    fclose( fp );
  }
  
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

  if( chims && num_chim_chunks_read != chims->num_breakers * CHIMERA_CHUNKS )
  {
    fprintf( stderr, "Expected %d relevant IUMs, but found only " F_U32 "\n",
             chims->num_breakers * CHIMERA_CHUNKS, num_chim_chunks_read );
    return 1;
  }
  if( spurs && num_spur_chunks_read != spurs->num_breakers * SPUR_CHUNKS )
  {
    fprintf( stderr, "Expected %d relevant IUMs, but found only " F_U32 "\n",
             spurs->num_breakers * SPUR_CHUNKS, num_spur_chunks_read );
    return 1;
  }
    
  return 0;
}





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
static char *rcsid = "$Id:";

#include "AS_CGB_all.H"
#define UUNIQUE

/*
  Count instances of the following pattern, where (s) is the chimeric fragment:
     \                                /
    -(O2O(a),a)--------------(b,O2O(b))-
     /                       /        \
                            /
                           /
                          /
                 (s,O2O(s))
                 /
                /
               /
    \         /                       /
    -(O2O(c),c)--------------(d,O2O(d))-
    /                                 \

  pattern matching rules:
  s must consist of exactly one fragment
  s must have exactly one b-mate (b)
  s must have exactly one a-mate (c)
  b must have exactly two mates off s-side port (s,a)
  c must have exactly two mates off s-side port (s,d)
  b != c (chunk within a chunk?)
  a must have exactly one mate off b-side port (b)
  d must have exactly one mate off c-side port (c)
  a != d (quadrilateral, not a Z)
  a != c not a diamond
  b != d not a diamond
  a,b,c,d are u-unitigs

  when found, write out internal accession #'s: s b,a c,d

  To do: Mark the fragment s, and the overlaps (s,b) and (s,c) as
  chimeric.

*/

#define O2CI(a)    ((a) >> 1)
// From the chunk-end, return the chunk index.
#define O2S(a)     ((a) & 1)
// From the chunk-end, return the chunk suffix.
#define CIS2O(a,b) (((a) << 1) + ((b) & 1))
// From the chunk index and suffix, return the chunk-end.
#define O2O(a)     ((a) ^ 1)
// From the chunk-end, return the mate chunk-end.

static IntChunk_ID get_second_mate_index
(
  TChunkMesg chunks[],
  const Tfragment frags[],
  const Tedge edges[],
  const IntChunk_ID sentinel,
  const IntChunk_ID b_index,
  const IntChunk_ID other_index,
  const float cgb_unique_cutoff
)
{
  int i;
  const AChunkMesg * chunk = GetVA_AChunkMesg(chunks,O2CI(b_index));
  IntChunk_ID mate_index = sentinel;
  IntChunk_ID temp_index = sentinel;
  const int b_side = O2S(b_index);
  const int32 raw_degree =
    (b_side == 1) ? chunk->b_degree_raw : chunk->a_degree_raw;
  const IntEdge_ID edge_id = (b_side == 1) ? chunk->b_list_raw : chunk->a_list_raw;

  // loop over the edges
  for( i = 0; i < raw_degree; i++ )
  {
    const Tnes edge = get_nes_edge( edges, edge_id + i );

    // if it's the right kind of edge
    if( AS_CGB_INTERCHUNK_EDGE == edge )
    {
      // get the chunk index
      IntChunk_ID temp_chunk = get_chunk_index( chunks, frags,
					  edges, edge_id + i );
      int temp_suffix = get_chunk_suffix( chunks, frags,
					  edges, edge_id + i );
      temp_index = CIS2O( temp_chunk, temp_suffix);

      if( (temp_index != other_index )
#ifdef UUNIQUE
	  &&
	  (GetVA_AChunkMesg( chunks, temp_chunk )->coverage_stat >= cgb_unique_cutoff)
#endif
	  )
      {
        if( mate_index == sentinel )
        {
          // if it's the first non-other index, set the mate_index
          mate_index = temp_index;
        }
        else if( temp_index != mate_index )
        {
          // it's a third mate, return sentinel
          return sentinel;
        }
      }
    }
  }
  // return either unique second mate or sentinel
  return mate_index;
}

static IntChunk_ID get_lone_mate_index
(
 TChunkMesg chunks[],
 const Tfragment frags[],
 const Tedge edges[],
 const IntChunk_ID sentinel,
 const IntChunk_ID b_index,
 const float cgb_unique_cutoff
)
{
  int i;
  IntChunk_ID mate_index = sentinel;
  const AChunkMesg * chunk = GetVA_AChunkMesg( chunks, O2CI( b_index ) );
  const int b_side = O2S( b_index );
  const int32 raw_degree =
    (b_side == 1) ? chunk->b_degree_raw : chunk->a_degree_raw;
  const IntEdge_ID edge_id = (b_side == 1)
    ? chunk->b_list_raw : chunk->a_list_raw;

  // loop over the edges
  for( i = 0; i < raw_degree; i++ )
  {
    const Tnes edge = get_nes_edge( edges, edge_id + i );

    // if it's the right kind of edge
    if( AS_CGB_INTERCHUNK_EDGE == edge )
    {
      const IntChunk_ID temp_index = get_chunk_index( chunks, frags, edges, edge_id + i );

      if(
	 (mate_index == sentinel )
#ifdef UUNIQUE
	 &&
	 (GetVA_AChunkMesg( chunks, temp_index )->coverage_stat >= cgb_unique_cutoff )
#endif
	 ) {
	// first mate found
	mate_index =
	  CIS2O( get_chunk_index( chunks, frags, edges, edge_id + i ),
		 get_chunk_suffix( chunks, frags, edges, edge_id + i ) );
      }
      else if( O2CI( mate_index ) !=
	       get_chunk_index( chunks, frags, edges, edge_id + i ) ) {
	// it's a second mate, return sentinel
	return sentinel;
      }
    }
  }
  // return either lone mate or sentinel
  return mate_index;
}

static int is_hanging_chunk_end
(
 TChunkMesg chunks[],
 const Tfragment frags[],
 const Tedge edges[],
 const IntChunk_ID b_index
)
{
  int i;
  const AChunkMesg * chunk = GetVA_AChunkMesg( chunks, O2CI( b_index ) );
  const int b_side = O2S( b_index );
  const int32 raw_degree =
    (b_side == 1) ? chunk->b_degree_raw : chunk->a_degree_raw;
  const IntEdge_ID edge_id = (b_side == 1)
    ? chunk->b_list_raw : chunk->a_list_raw;

  // loop over the edges
  for( i = 0; i < raw_degree; i++ )
  {
    const Tnes edge = get_nes_edge( edges, edge_id + i );

    switch(edge) {
    // if it's any dovetail edge
    case AS_CGB_INTERCHUNK_EDGE:
    case AS_CGB_INTRACHUNK_EDGE:
    case AS_CGB_TOUCHES_CONTAINED_EDGE:
    case AS_CGB_BETWEEN_CONTAINED_EDGE:
      //case AS_CGB_TOUCHES_CRAPPY_EDGE:
    case AS_CGB_MARKED_BY_BRANCH_DVT:
    case AS_CGB_REMOVED_BY_THRESHOLD_DVT:
    case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
      return FALSE;
    default:
      break;
    }
  }
  return TRUE;
}


uint32 count_chimeras
(
 const char * const chimeras_report_filename,
 const float cgb_unique_cutoff,
 const Tfragment  frags[],
 const Tedge      edges[],
 TChunkFrag       chunk_frags[],
 TChunkMesg       chunks[]
)
{
  uint32         num_chimeras = 0;
  const IntChunk_ID  num_chunks = GetNumVA_AChunkMesg( chunks );
  const IntChunk_ID  sentinel = CIS2O( num_chunks, 0 );
  IntChunk_ID        schunk;
  FILE             * fp_chimeras;

  fp_chimeras = fopen( chimeras_report_filename, "w" );
  if( fp_chimeras == NULL )
  {
    fprintf( stderr, "Failed to open file %s for writing\n",
             chimeras_report_filename );
    return 0;
  }

  // loop over chunks
  for( schunk = 0; schunk < num_chunks; schunk++ )
  {
    const IntChunk_ID  s_index = CIS2O(schunk,0);
    IntChunk_ID  a_index = sentinel;
    IntChunk_ID  b_index = sentinel;
    IntChunk_ID  c_index = sentinel;
    IntChunk_ID  d_index = sentinel;

    assert( schunk == O2CI(s_index) );
    // s must consist of exactly one fragment
    if( (GetVA_AChunkMesg( chunks, O2CI(s_index)))->num_frags != 1 )
      continue;

    // s must have exactly one b-mate (b)
    b_index =
      get_lone_mate_index( chunks, frags, edges, sentinel, O2O(s_index), cgb_unique_cutoff );
    if( b_index == sentinel )
      continue;

    // s must have exactly one a-mate (c)
    c_index =
      get_lone_mate_index( chunks, frags, edges, sentinel, s_index, cgb_unique_cutoff );
    if( c_index == sentinel )
      continue;

    // b != c (chunk within a chunk?)
    if( O2CI( c_index ) == O2CI( b_index ) )
      continue;

    // b must have exactly two a-mates (s,a)
    a_index =
      get_second_mate_index( chunks, frags, edges,
			     sentinel, b_index, O2O(s_index), cgb_unique_cutoff );
    if( a_index == sentinel )
      continue;

    // c must have exactly two b-mates (s,d)
    d_index =
      get_second_mate_index( chunks, frags, edges,
			     sentinel, c_index, s_index, cgb_unique_cutoff );
    if( d_index == sentinel )
      continue;

    // a !=d (quadrilateral, not a Z)
    if( O2CI( a_index ) == O2CI( d_index ) )
      continue;

    // a != c (quadrilateral, not a Z)
    if( O2CI( a_index ) == O2CI( c_index ) )
      continue;

    // b != d (quadrilateral, not a Z)
    if( O2CI( b_index ) == O2CI( d_index ) )
      continue;

    // a must have exactly one b-mate (b)
    if( get_lone_mate_index( chunks, frags, edges,
                             sentinel, a_index, cgb_unique_cutoff ) != b_index )
      continue;

    // d must have exactly one a-mate (c)
    if( get_lone_mate_index( chunks, frags, edges,
                             sentinel, d_index, cgb_unique_cutoff ) != c_index )
      continue;

    // if here, it's a chimeric fragment
    if( num_chimeras == 0 )
    {
      fprintf( fp_chimeras,
             "# s = index of chunk consisting of single chimeric fragment\n" );
      fprintf( fp_chimeras,
        "# a & b are split by suffix of s, c & d are split by prefix of s\n" );
      fprintf( fp_chimeras,
               "# a2b is 1 if suffix of a has edge with b, and so on.\n" );
      fprintf( fp_chimeras, "#\n" );
      fprintf( fp_chimeras, "#    \\                     /\n" );
      fprintf( fp_chimeras, "#    -(a)---------------(b)-\n" );
      fprintf( fp_chimeras, "#    /                  /  \\\n" );
      fprintf( fp_chimeras, "#                     /\n" );
      fprintf( fp_chimeras, "#                   /\n" );
      fprintf( fp_chimeras, "#                 /\n" );
      fprintf( fp_chimeras, "#              (s)\n" );
      fprintf( fp_chimeras, "#             /\n" );
      fprintf( fp_chimeras, "#           /\n" );
      fprintf( fp_chimeras, "#         /\n" );
      fprintf( fp_chimeras, "#    \\  /                  /\n" );
      fprintf( fp_chimeras, "#    -(c)---------------(d)-\n" );
      fprintf( fp_chimeras, "#    /                     \\\n" );
      fprintf( fp_chimeras, "#\n" );
      fprintf( fp_chimeras,
               "# "
	       "si : (ai,as) (bi,bs) (ci,cs) (di,ds)\n");
      fprintf( fp_chimeras,
               "#-------------------------------------------------"
               "----------------------------------------\n" );
    }
    num_chimeras++;
    fprintf( fp_chimeras,
	     "%10" F_IIDP " :  (%10" F_IIDP ",%1d) (%10" F_IIDP ",%1d) (%10" F_IIDP ",%1d) (%10" F_IIDP ",%1d)\n",
             O2CI(s_index ),  //s_chunk->iaccession,
	     O2CI(a_index ),  //a_chunk->iaccession,
             O2S( a_index ),
             O2CI(b_index ),  //b_chunk->iaccession,
             O2S( b_index ),
             O2CI(c_index ),  //c_chunk->iaccession,
             O2S( c_index ),
             O2CI(d_index ),  //d_chunk->iaccession
             O2S( d_index )
	     );
  }

  fclose( fp_chimeras );
  return num_chimeras;
}





/*
  Count instances of the following pattern, where (s) is the crappy fragment:

                 (s,O2O(s))
                 /
                /
               /
    \         /                       /
    -(O2O(c),c)--------------(d,O2O(d))-
    /                                 \

  pattern matching rules:
  s must consist of exactly one fragment
  s must have exactly one a-mate (c)
  s must have exactly zero b-mates
  c must have exactly two mates off s-side port (s,d)
  d must have exactly one mate off c-side port (c)

  c,d are u-unitigs

  when found, write out internal accession #'s: s c,d

  To do: Mark the fragment s, and the overlap (s,c) as chimeric.

*/


uint32 count_crappies
(
 const char * const crappies_report_filename,
 const float cgb_unique_cutoff,
 const Tfragment  frags[],
 const Tedge      edges[],
 TChunkFrag       chunk_frags[],
 TChunkMesg       chunks[]
)
{
  uint32         num_crappies = 0;
  const IntChunk_ID  num_chunks = GetNumVA_AChunkMesg( chunks );
  const IntChunk_ID  sentinel = CIS2O( num_chunks, 0 );
  IntChunk_ID        schunk;
  int                ssuffix;
  FILE             * fp_crappies;

  fp_crappies = fopen( crappies_report_filename, "w" );
  if( fp_crappies == NULL )
  {
    fprintf( stderr, "Failed to open file %s for writing\n",
             crappies_report_filename );
    return 0;
  }

  // loop over chunks
  for( schunk = 0; schunk < num_chunks; schunk++ )
    for( ssuffix = 0; ssuffix < 2; ssuffix++ )
  {
    const IntChunk_ID  s_index = CIS2O(schunk,ssuffix);
    IntChunk_ID  c_index = sentinel;
    IntChunk_ID  d_index = sentinel;

    assert( schunk == O2CI(s_index) );
    // s must consist of exactly one fragment
    if( (GetVA_AChunkMesg( chunks, O2CI(s_index)))->num_frags != 1 )
      continue;

    // s must have exactly zero b-mates
    if( FALSE == is_hanging_chunk_end( chunks, frags, edges, O2O(s_index)) )
      continue;

    // s must have exactly one a-mate (c)
    c_index =
      get_lone_mate_index( chunks, frags, edges, sentinel, s_index, cgb_unique_cutoff );
    if( c_index == sentinel )
      continue;

    // c must have exactly two b-mates (s,d)
    d_index =
      get_second_mate_index( chunks, frags, edges,
			     sentinel, c_index, s_index, cgb_unique_cutoff );
    if( d_index == sentinel )
      continue;

    // d must have exactly one a-mate (c)
    if( get_lone_mate_index( chunks, frags, edges,
                             sentinel, d_index, cgb_unique_cutoff ) != c_index )
      continue;

    // if here, it's a crappy fragment.
    if( num_crappies == 0 )
    {
      fprintf( fp_crappies,
	       "# si is the chunk consisting of single crappy fragment\n" );
      fprintf( fp_crappies,
	       "# ci & di are split by prefix/suffix of s\n" );
      fprintf( fp_crappies, "#\n" );
      fprintf( fp_crappies, "#           (s)\n" );
      fprintf( fp_crappies, "#           /\n" );
      fprintf( fp_crappies, "#          /\n" );
      fprintf( fp_crappies, "#         /\n" );
      fprintf( fp_crappies, "#    \\   /                /\n" );
      fprintf( fp_crappies, "#    -(c)---------------(d)-\n" );
      fprintf( fp_crappies, "#    /                     \\\n" );
      fprintf( fp_crappies, "#\n" );
      fprintf( fp_crappies,
               "# "
	       "(si,ss) (ci,cs) (di,ds)\n");
      fprintf( fp_crappies,
               "#-------------------------------------------------"
               "----------------------------------------\n" );
    }
    num_crappies++;
    fprintf( fp_crappies,
	     "(%10" F_IIDP ",%1d) (%10" F_IIDP ",%1d) (%10" F_IIDP ",%1d)\n",
             O2CI(s_index ),  //s_chunk->iaccession,
	     O2S( s_index ),
             O2CI(c_index ),  //c_chunk->iaccession,
             O2S( c_index ),
             O2CI(d_index ),  //d_chunk->iaccession
             O2S( d_index )
	     );
  }

  fclose( fp_crappies );
  return num_crappies;
}


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

//  $Id: AS_CGB_unitigger_globals.h,v 1.12 2007-07-20 08:41:43 brianwalenz Exp $

#ifndef AS_CGB_UNITIGGER_GLOBALS_INCLUDE
#define AS_CGB_UNITIGGER_GLOBALS_INCLUDE

#include "AS_CGB_all.h"

typedef struct {
  IntFragment_ID  avx,bvx;
  int16       ahg,bhg; 

  uint32    quality : 16;  // zero is a perfect overlap
  uint32    nes : 8;      // The edge labeling.
  uint32    asx : 1;
  uint32    bsx : 1;
  uint32    reflected : 1; // Is this a direction mate edge.
  uint32    grangered : 1;
  uint32    invalid : 1; // Is this edge known to be invalid?
  uint32    blessed : 1;
  uint32    bit6 : 1;
  uint32    bit7 : 1;
} Aedge;


typedef struct { 
  //uint64  src;
  /* An index into the character array
     "the_source_text" holding the simulator source
     info. */

  IntFragment_ID  iid; /* The Celera Genomics Assembler Gatekeeper IID 
		      for this fragment read. */
  IntChunk_ID cid; /* The linkage to chunks. Needed due to the way
                      that interchunk edges are made. */
  IntEdge_ID  segbgn_prefix;/* The starting location of a segment of 
			   raw fragment overlaps to the prefix. */
  IntEdge_ID  segbgn_suffix;/* The starting location of a segment of 
			   fragment overlaps to the suffix. */
  IntEdge_ID  segend_prefix;/* The ending location of a segment of 
			   raw fragment overlaps to the prefix. */
  IntEdge_ID  segend_suffix;/* The ending location of a segment of 
			   fragment overlaps to the suffix. */

  int32   nprefix_all;  /* The segment length of the prefix-only 
			   fragment overlaps. */
  int32   nsuffix_all;  /* The segment length of the suffix-only 
			   fragment overlaps. */

  int32   nprefix_dvt; /* Recomputable */
  // The segment length of the prefix-only fragment dovetail edges.
  int32   nsuffix_dvt; /* Recomputable */
  // The segment length of the suffix-only fragment dovetail edges.
  int32   nprefix_frc; /* Recomputable */
  // The segment length of the prefix-only fragment containment edges.
  int32   nsuffix_frc; /* Recomputable */
  // The segment length of the suffix-only fragment containment edges.

  int32   raw_prefix_dvt_count;
  // The segment length of the raw prefix-only fragment dovetail overlaps.
  int32   raw_suffix_dvt_count;
  // The segment length of the raw suffix-only fragment dovetail overlaps.

  int32   raw_frc_count;
  int32   raw_toc_count;
  
  IntFragment_ID  container;  /* Recomputable */
  // A zero value means that the fragment is not contained by any
  // other fragment.  A positive value is the fragment IID of the
  // container fragment.  This scheme depends on the Assembler IO
  // convention that Fragment IIDs are positive integers.

  int64       offset3p;
  int64       offset5p;

  FragType    frag_type;
  // The laboratory designation of the fragment type.

  Tlab        label; /* Recomputable */
  // The chunk labeling of fragments.
  // AS_CGB_THRU_FRAG, AS_INTERCHUNK_FRAG,
  // AS_CGB_INTRACHUNK_FRAG, AS_CGB_CONTAINED_FRAG,
  // AS_CGB_DELETED_FRAG, etc ?

  int16   bp_length;
  // The length of the fragment read in bp.
  
  /* FGB bit flags: */
  unsigned int deleted : 1;
  // A flag indicating if this fragment is deleted in the fragment
  // graph.
  unsigned int contained : 1;  /* Recomputable */
  // A flag indicating if this fragment is contained in the fragment
  // graph.
  unsigned int spur : 1;

  unsigned int prefix_blessed : 1;
  unsigned int suffix_blessed : 1;
} Afragment;

typedef struct {
  IntFragment_ID vid;
} AChunkFrag;

typedef struct {
  int64         bp_length;  
  // the length in base pairs of the chunk.
  int64         rho;        
  IntChunk_ID   iaccession;
  // An arbitrary, but dense and non-negative enumeration of the
  // unitigs.
  float       coverage_stat;
  // Gene^s coverage statistic.

  IntFragment_ID    chunk_avx; /* The A and B vertices of the chunk. */
  IntFragment_ID    chunk_bvx;
  int32     chunk_asx;
  int32     chunk_bsx;
  IntFragment_ID    num_frags;
  IntFragment_ID    f_list; /* The index into a TChunkFrag array. */

  // The following is not necessary. The info available from
  // fragment-end info, chunk_avx, chunk_asx, chunk_bvx, and chunk_bsx.
  int32         a_degree_raw; 
  int32         b_degree_raw;

  IntEdge_ID        a_list_raw;
  IntEdge_ID        b_list_raw;

  size_t        seq_loc; /* An index into a character array that
			 stores a consensus sequence of the chunk. */
  size_t        seq_len; // The length of the consensus sequence of the chunk.
  size_t        qua_loc; /* An index into a character array that
			 stores a quality sequence of the chunk. */
  size_t        qua_len; // The length of the quality sequence of the chunk.

  int8 asl,bsl; 
  // Flags indicating whether the A and B chunk-ends are probably
  // touching a tandem repeat region.
} AChunkMesg;



VA_DEF(Aedge)
VA_DEF(Afragment)
VA_DEF(AChunkFrag)
VA_DEF(AChunkMesg)
VA_DEF(IntEdge_ID)

typedef VA_TYPE(Aedge)      Tedge;
typedef VA_TYPE(Afragment)  Tfragment;
typedef VA_TYPE(AChunkFrag) TChunkFrag;
typedef VA_TYPE(AChunkMesg) TChunkMesg;
typedef VA_TYPE(IntEdge_ID) TIntEdge_ID;


typedef struct { 
  Tfragment      *frags;
  Tedge          *edges;
  TIntEdge_ID    *next_edge_obj;
  TChunkFrag     *chunkfrags;
  TChunkMesg     *thechunks;
  int64           nbase_in_genome;
  float           global_fragment_arrival_rate;
} THeapGlobals;


typedef struct {
  char * frag_store;
  char * chimeras_file;
  char * spurs_file;
  char   bubble_overlaps_filename[FILENAME_MAX];
  char * blessed_overlaps_input_filename;
  char * blessed_overlaps_output_filename;
  
  /* Reaper/FGB input */
  int     num_ovl_files;
  char ** the_ovl_files;
  
  char * OVL_Store_Path;
  // The directory containing the overlap store.

  char * Output_Graph_Store_Prefix;
  // The directory containing the output graph store files.

  char * Dump_File_Name;
  // The file name of the output fragments and overlaps in OVL format.

  char * ovl_files_list_fname;

  int            recalibrate_global_arrival_rate;
  int            dechord_the_graph;
  int            create_dump_file;
  int            walk_depth;

  float          cgb_unique_cutoff;

  IntFragment_ID maxfrags;
  IntEdge_ID     maxedges;
  size_t         maxtext;

  int            work_limit_per_candidate_edge;
  int            dvt_double_sided_threshold_fragment_end_degree;
  int            con_double_sided_threshold_fragment_end_degree;
  int            intrude_with_non_blessed_overlaps_flag;
  int            cutoff_fragment_end_degree;
  uint32         overlap_error_threshold;

  int64          genome_length;

  int            work_limit_placing_contained_fragments;
  int            output_iterations_flag;
  int            aggressive_spur_fragment_marking;
  int            bubble_smoothing_flag;

  int            use_consensus;
  int            dont_count_chimeras;
  int            fragment_count_target;
} UnitiggerGlobals;

int main_fgb (THeapGlobals  * heapva,
              UnitiggerGlobals * rg);

int main_cgb(THeapGlobals  * heapva,
             UnitiggerGlobals * rg);

#endif // AS_CGB_UNITIGGER_GLOBALS_INCLUDE

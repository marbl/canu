
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
/*********************************************************************
 * $Id: AS_FGB_io.h,v 1.5 2007-01-29 20:40:58 brianwalenz Exp $
 *
 * Module: AS_FGB_io.h
 * Description: Header file for the code that reads and writes the 
 * check point data.
 * Assumptions:
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_FGB_IO_INCLUDE
#define AS_FGB_IO_INCLUDE

void InsertFragmentsIntoGraph
(
  THeapGlobals * heapva,
  VA_TYPE(OFGMesg)  * the_ofg_messages, // Additional fragments in a VA.
  VA_TYPE(char)     * the_ofg_source,
  FragmentHashObject * afr_to_avx,        // Part of a hash table replacement
  IntFragment_ID * Pnofg,
  IntFragment_ID * Pmin_frag_iid,
  IntFragment_ID * Pmax_frag_iid,
  VA_TYPE(char)  * frag_annotations
  );

void InsertOverlapsIntoGraph
(
  THeapGlobals * heapva,
  VA_TYPE(OverlapMesg) * the_ovl_messages, // Additional overlaps in a VA.
  VA_TYPE(char)        * the_ovl_source,
  FragmentHashObject   * afr_to_avx,        // Part of a hash table replacement
  const int dvt_double_sided_threshold_fragment_end_degree,
  const int con_double_sided_threshold_fragment_end_degree,
  const int intrude_with_non_blessed_overlaps_flag,
  const CGB_ERATE_TYPE overlap_error_threshold
  );

void input_messages_from_a_file
(int        argc, 
 char       *argv[],
 FILE       *fovl,
 FILE       *filk,
 Tfragment  frags[],
 // The internal representation of the fragment reads. 
 Tedge      edges[],
 // The internal representation of the overlaps.
 VA_TYPE(char) frag_annotations[],
 FragmentHashObject *afr_to_avx,
 IntFragment_ID    *min_frag_iid,
 IntFragment_ID    *max_frag_iid,
 TIntEdge_ID       *next_edge_obj,
 BPTYPE *nbase_in_genome,
 const int dvt_double_sided_threshold_fragment_end_degree,
 const int con_double_sided_threshold_fragment_end_degree,
 const int intrude_with_non_blessed_overlaps_flag,
 const CGB_ERATE_TYPE overlap_error_threshold
 );


void process_ovl_store
(
 char * OVL_Store_Path,
 Tfragment  frags[],
 // The internal representation of the fragment reads. 
 Tedge      edges[],
 // The internal representation of the overlaps.
 VA_TYPE(char) frag_annotations[],
 FragmentHashObject *afr_to_avx,
 IntFragment_ID    *min_frag_iid,
 IntFragment_ID    *max_frag_iid,
 TIntEdge_ID         next_edge_obj[],
 BPTYPE *nbase_in_genome,
 const int dvt_double_sided_threshold_fragment_end_degree,
 const int con_double_sided_threshold_fragment_end_degree,
 const int intrude_with_non_blessed_overlaps_flag,
 const CGB_ERATE_TYPE overlap_error_threshold
 );
  
#endif // AS_FGB_IO_INCLUDE

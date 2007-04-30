
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
 * $Id: AS_CGB_unitigger_globals.h,v 1.7 2007-04-30 13:00:29 brianwalenz Exp $
 *
 * Module: AS_CGB_unitigger_globals.h
 * Description: A subroutine interface for the unitigger.
 * Assumptions: Too many to count.
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_UNITIGGER_GLOBALS_INCLUDE
#define AS_CGB_UNITIGGER_GLOBALS_INCLUDE

#include "AS_CGB_all.h"

// Global parameters structure for this application
typedef struct
{
  char * frag_store;
  char * chimeras_file;
  char * ovl_file;
  char * spurs_file;
  char * iid_file;
  char * bubble_overlaps_filename;
  char * bubble_boundaries_filename;
  char * blessed_overlaps_input_filename;
  char * blessed_overlaps_output_filename;
  
  /* Reaper/FGB input */
  int     num_ovl_files;
  char ** the_ovl_files;
  
  char * program_name;
  char * OVL_Store_Path;
  // The directory containing the overlap store.

  char * Output_Graph_Store_Prefix;
  // The directory containing the output graph store files.

  char * Dump_File_Name;
  // The file name of the output fragments and overlaps in OVL format.

  char * ovl_files_list_fname;
  char * Fragment_Subset_IIDs_File_Name;
  // A file that subsets the active fragments by intersection.

  int            as_proto_output;

  int            dechord_the_graph;
  int            create_dump_file;
  int            check_point_level;
  int            compress_the_graph;
  int            walk_depth;
  int            reaper_pass;
  IntFragment_ID iv_start;
  int            debug_level;
  int            verbosity_level;
  int            analysis_level;
  int            analysis_flag;
  float          cgb_unique_cutoff;
  IntFragment_ID as_cgb_max_frag_iid;
  IntFragment_ID maxfrags;
  IntEdge_ID     maxedges;
  size_t         maxtext;
  int            work_limit_per_candidate_edge;
  int            remove_blizzard_overlaps;
  int            use_all_overlaps_in_reaper_pass;
  int            dvt_double_sided_threshold_fragment_end_degree;
  int            con_double_sided_threshold_fragment_end_degree;
  int            intrude_with_non_blessed_overlaps_flag;
  int            cutoff_fragment_end_degree;
  uint32         overlap_error_threshold;

  BPTYPE nbase_in_genome;

  int work_limit_placing_contained_fragments;
  int output_iterations_flag;
  int aggressive_spur_fragment_marking;
  int bubble_smoothing_flag;

  int use_consensus;
  int dont_find_branch_points;
  int dont_count_chimeras;
  int fragment_count_target;
  int num_cgb_passes;

} UnitiggerGlobals;

int main_fgb (TStateGlobals * gstate,
              THeapGlobals  * heapva,
              UnitiggerGlobals * rg);

int main_cgb(TStateGlobals * gstate,
             THeapGlobals  * heapva,
             UnitiggerGlobals * rg);

#endif // AS_CGB_UNITIGGER_GLOBALS_INCLUDE

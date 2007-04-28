
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
 * $Id: AS_CGB_walk.h,v 1.5 2007-04-28 08:46:21 brianwalenz Exp $
 *
 * Module: AS_CGB_walk.h
 * Description: 
 * Assumptions:
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_WALK_INCLUDE
#define AS_CGB_WALK_INCLUDE

#undef WALK_DEPTH_DIAGNOSTICS
#define MATCH_TARGET_EDGE

int is_there_an_overlap_path
(
 const Tfragment * const frags,
 const Tedge     * const edges,
 const IntFragment_ID target_avx,
 const int        target_asx,
#ifdef MATCH_TARGET_EDGE
 const IntFragment_ID target_bvx,
 const int        target_bsx,
 const int        target_ahg,
 const int        target_bhg,
 const Tnes       target_nes,
 const int        target_is_dovetail,
 const int        target_is_from_contained,
 const int        target_is_to_contained,
 const int        target_is_dgn_contained,
#endif // MATCH_TARGET_EDGE
 /* recursion variables: */
 const int        search_depth, // Use zero at the top level.
 const IntFragment_ID current_avx,  // Use target_avx at the top level.
 const int        current_asx,  // Use target_asx at the top level.
 const int        current_ahg,  // Use zero at the top level.
 const int        current_bhg,  // Use zero at the top level.
 const int        last_edge_was_contained, // From target_nes at the top level.
 /* search path limiting: */
 const int        walk_depth,   // The maximum depth of the stack.
 const int        tolerance,    // For the overlaps
 IntFragment_ID visited_a[],
#ifdef MATCH_TARGET_EDGE
 IntFragment_ID visited_b[],
#endif // MATCH_TARGET_EDGE
 // Was this fragment visited from the target overlap before? That is
 // by the same (target_avx, target_asx, target_bvx, target_bsx).
 const int        work_limit_per_candidate_edge,
 // Maximum number of edges to explore per candidate edge.
 int *            work_tally_per_candidate_edge,
 // Current number of edges explored per candidate edge.
 /* diagnostic variables: */
#ifdef WALK_DEPTH_DIAGNOSTICS
 int64      search_depth_histogram[], // How many times this depth is visited.
 int64      search_path_histogram[],  // How many paths were this length.
#endif // WALK_DEPTH_DIAGNOSTICS
 int64      *ntrans_test_fail // How many paths failed the last test.
);

#endif // AS_CGB_WALK_INCLUDE


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
 * $Id: AS_CGB_fgb.h,v 1.2 2004-09-23 20:25:01 mcschatz Exp $
 *
 * Module: AS_CGB_fgb.h
 * Description: 
 * Assumptions:
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_FGB_INCLUDE
#define AS_CGB_FGB_INCLUDE

void transitive_edge_marking
(
 TStateGlobals * gstate, // For time interval check pointing
 THeapGlobals  * heapva, // For time interval check pointing
 Tfragment     * frags,
 Tedge         * edges,
 TIntEdge_ID   * next_edge_obj,
 const int walk_depth,
 const int cutoff_fragment_end_degree,
 const int work_limit_per_candidate_edge,
 const IntFragment_ID iv_start,
 const int analysis_flag,
 const char Output_Graph_Store[]
 );

void reorder_edges
( Tfragment *frags,
  Tedge *edges,
  TIntEdge_ID *next_edge_obj);

#endif // AS_CGB_FGB_INCLUDE



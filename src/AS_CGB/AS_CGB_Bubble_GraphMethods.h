
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
#ifndef _AS_CGB_BUBBLE_GRAPH_METHODS_H_
#define _AS_CGB_BUBBLE_GRAPH_METHODS_H_

/* This method uses a depth first search to mark all(*) fragments with
   the AS_CGB_BUBBLE_V_VALID flag.  As part of the process, the position
   and orientation of each fragment is set, and some edges become marked
   as invalid if they represent back edges in the DFS tree.

   (*) Not all fragments are marked.  Searches are only begun with
   non-contained fragments, and only non-containment edges are used.
   Contained fragments which are not touched by dovetail edges can
   therefore be omitted.  */

/* Implemented in "AS_CGB_Bubble_dfs.c" */
void AS_CGB_Bubble_dfs(BubGraph_t bg);

/* This method removes the AS_CGB_BUBBLE_V_VALID flag from certain contained
   fragments.  The fragments "removed" either do not have the valid flag
   to begin with (and hence are not reachable by dovetails), or are
   "consistent".  A consistent contained fragment C with parent fragment
   P is one such that:

   (1) All fragments contained by C are contained by P
   (2) All fragments containing 
   (2) All fragments with dovetails to C have either a dovetail or containment
   relationship with P.

   Containments with multiple parents must satisfy the above criteria
   for all parents.  NOTE: This does leave open the odd possibility
   that a contained fragment has two parents P1 and P2 satisfying the
   conditions, but there is somehow no relationship between P1 and P2.
   Should this also be considered? */

/* Implemented in "AS_CGB_Bubble_rcc.c" */
/* void AS_CGB_Bubble_remove_consistent_contains(BubGraph_t bg); */



/* Performs a topological sort of the fragment graph.  The output is an array 
   of fragment IDs in topological order, pointed to by the out parameter (which
   is assumed to point to sufficient pre-allocated space).  Only those 
   fragments and edges which are marked as VALID are used. 

   The return value is the number of elements placed in topological order.
   If the graph is found to be cyclic, 0 is returned instead. */

/* Implemented in "AS_CGB_Bubble_top.c" */
IntFragment_ID
AS_CGB_Bubble_topo_sort(BubGraph_t bg, IntFragment_ID *out);

#endif


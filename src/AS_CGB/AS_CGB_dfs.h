
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
 * $Id: AS_CGB_dfs.h,v 1.4 2005-03-22 19:48:27 jason_miller Exp $
 *
 * Module: AS_CGB_dfs.h
 * Description: 
 * 
 * 1. Overview
 * 
 * This functional unit transverses the current state of the marked graph
 * and returns a best ranking of the verticies.
 * 
 * 2. Memory Usage
 * 

 * 3. Interface
 * 
 * 4. Design
 * 
 * The graph is transversed in a depth first manner.
 * The actual DFS tree is dependent on the sorting of the adjacency lists.
 * 
 * 5. Limitations
 * 
 * 
 * 6. Status
 * 
 * Assumptions:
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_DFS_INCLUDE
#define AS_CGB_DFS_INCLUDE

void as_dfs_graph_traversal
(
 FILE *fout,
 Tfragment frags[],
 Tedge edges[], 
 IntRank fragment_visited[]
 );

#endif /* AS_CGB_DFS_INCLUDE */

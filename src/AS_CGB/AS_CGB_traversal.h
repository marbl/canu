
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
 * $Id: AS_CGB_traversal.h,v 1.3 2005-03-22 19:02:32 jason_miller Exp $
 *
 * Module: AS_CGB_traversl.h
 * Description: 
 * 
 * 1. Overview
 * 
 * This functional unit transverses the current state of the marked graph
 * and returns a best ranking of the verticies.
 * 
 * 2. Memory Usage
 * 
 * There is one statically alloacted first-in first-out queue of length
 * MAXQUEUELEN.  There is no dynamicly alloacted memory.
 * 
 * 3. Interface
 * 
 * 4. Design
 * 
 * The graph is transversed in a breadth first manner, except over chunks
 * which are visisted in a depth first manner.
 * 
 * 5. Limitations
 * 
 * Currently it is a assumed that a chunk has two ends.  Thus a circular
 * chunk is not visited.
 * 
 * 6. Status
 * 
 * Assumptions:
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_TRAVERSAL_INCLUDE
#define AS_CGB_TRAVERSAL_INCLUDE

void as_graph_traversal(
			FILE *fout,
			Tfragment frags[],
			Tedge edges[], 
			IntRank fragment_visited[]
			);
#endif /* AS_CGB_TRAVERSAL_INCLUDE */

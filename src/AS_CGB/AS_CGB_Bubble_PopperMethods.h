
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
#ifndef _AS_CGB_BUBBLE_POPPER_METHODS_H_
#define _AS_CGB_BUBBLE_POPPER_METHODS_H_

/* Finds all of the fragments associated with the bubble (NOTE:
   currently defined to be those fragments on a dovetail path between the
   start and end, inclusively).  Fills in the following fields of the
   BubblePopper (which should not be used prior to the call): curBubSize,
   vidToBid, bubFrags, adj.  Uses dfsStack.

   Returns TRUE if successful, or FALSE if the operation failed.  This
   is most likely because more than POPPER_MAX_BUBBLE_SIZE fragments were
   found, indicating that the bubble may not be closed. */
int
BP_find_bubble_dfs(BubblePopper_t bp, IntFragment_ID start,
		   IntFragment_ID end);

/* Computes the transitive closure of the adjacency array, which must
   represent a DAG. */
void
BP_transitive_closure(BubblePopper_t bp);

/* Computes the length of the longest path in the bubble, based on the
   adjacency array.  Uses a DAG algorithm, so if a cycle is found then
   -1 is returned (and the bubble probably shouldn't be popped anyway). */
int
BP_DAG_longest_path(BubblePopper_t bp);

/* Computes an approximation of the discriminator statistic for the bubble,
   assuming that the bubble collapses properly. */
float
BP_discriminator(BubblePopper_t bp);

#endif

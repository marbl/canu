
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
 $Id: AS_FGB_hanging_fragment.h,v 1.1.1.1 2004-04-14 13:50:07 catmandew Exp $
 Module: 
 Description: 
 Assumptions:
 Author: Clark M. Mobarry
 *********************************************************************/

#ifndef AS_CGB_HANGING_INCLUDE
#define AS_CGB_HANGING_INCLUDE

void separate_fragments_as_solo_hanging_thru
(
 Tfragment frags[],
 Tedge edges[]);


void identify_early_spur_fragments
(
 Tfragment frags[],
 Tedge edges[]);
// The philosophy of early spur removal is to assume that a
// hanging fragment should be ignored while the chunk backbones are
// determined.  We must make exceptions when there is evidence to
// consider the hanging fragment important.  In particular we need
// to include hanging fragments at sequencing gaps in order to get a
// tentative DNA sequence near the gap.  We also need to include a
// hanging fragment as important if its dovetail overlap is used to
// assign a fragment its "thru" status.  Also doubleton unitigs
// should be assembled rather than output as two hanging fragments.

// The chosen heuristic is if a hanging fragment "competes" in the
// dovetail adjacency list of each neighboring fragment with a thru
// fragment (in every location it can be placed), then aggresively
// ignore the hanging fragment.  This is a generalization of the
// original spur definition.  We do not assume that we know which
// fragments are contained nor assume that we have a transitively
// reduced overlap graph yet. The implementation uses a restricted
// incoming dovetail degree computed from only thru fragments.

// Determine the incoming dovetail degree at all fragment-ends
// restricted so that only fragments-ends contribute to the degree.
// This can be implemented by a scatter-add from the thru
// fragment-ends using their outgoing dovetail adjacency list to the
// edge's distal fragment-end.

// The fragment re-labelling can be implemented by (1) first
// re-assigning all hanging fragments to AS_CGB_HANGING_CRAPPY_FRAG,
// and (2) for all active edges if the proximal fragment is hanging
// (AS_CGB_HANGING_CRAPPY_FRAG or AS_CGB_HANGING_FRAG or ...) and
// the distal fragment active (and non-contained?) and has its
// restricted dovetail degree zero, then re-assign the proximal
// fragment as AS_CGB_HANGING_FRAG.

#endif // AS_CGB_HANGING_INCLUDE

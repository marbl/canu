
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
 * $Id: AS_CGB_edgetrimmer.h,v 1.4 2005-03-22 19:48:27 jason_miller Exp $
 *
 * Module: AS_CGB_edgetrimmer.h
 *
 * Description: The Chunk Graph Builder header file.
 *
 * The basic chunk data structures are:
 * 
 *   (1) a segmented array of type "TChunkFrag" called chunkfrags[] used
 *   to store the fragments of each chunk,
 *
 *   (2) a segmented array of type "TChunkOverlap" called chunkedges[]
 *   used to store the overlaps between chunks,
 *
 *   (3) an array of type "char" called chunksrc[] used to store
 *   annotation strings, and
 *
 *   (4) an array of type "TChunkMesg" thechunks[].
 *
 *   The first three arrays are segmented to store variable length
 *   information for the chunks.  Each member of the last array has a
 *   segment header index into each of the three previous arrays. 
 *
 *
 * Assumptions:
 * Author: Clark Mobarry
 *******************************************************************/

#ifndef AS_CGB_EDGETRIMMER_INCLUDE
#define AS_CGB_EDGETRIMMER_INCLUDE

void chunk_end_edge_trimmer
(
 /* input only */
 const float       cgb_unique_cutoff,
 const float       global_fragment_arrival_rate,
 const Tfragment   frags[],
 const TChunkFrag  chunkfrags[],
 /* modify */
 Tedge             edges[],
 TChunkMesg        thechunks[]
);

void fragment_end_edge_trimmer
(
 /* input only */
 const Tfragment     frags[],
 /* modify */
 Tedge               edges[]
);

#endif // AS_CGB_EDGETRIMMER_INCLUDE

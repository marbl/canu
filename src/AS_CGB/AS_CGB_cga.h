
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
 * $Id: AS_CGB_cga.h,v 1.2 2004-09-23 20:25:01 mcschatz Exp $
 *
 * Module: AS_CGB_cga.h
 *
 * Description: Header file for the chunk graph analyser.
 *
 * Assumptions:
 *
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_CGA_INCLUDE
#define AS_CGB_CGA_INCLUDE

void chunk_graph_analysis
(/* Input Only */
 const int        analysis_flag,
 const IntFragment_ID max_frag_iid,
 // The maximum fragment IID assigned by the Celera assembler gatekeeper.
 Tfragment        frags[],
 // The internal representation of the fragment reads.
 Tedge            edges[],
 // The internal representation of the overlaps.
 VA_TYPE(char)    frag_annotations[], 
 // A store of simulator fragment annotations.
 const BPTYPE     nbase_in_genome,
 // The estimated length of the genome.
 const float cgb_unique_cutoff,
 // The threshold A-statistic score for classifying a unitig as a
 // u-unitig.
 const float      global_fragment_arrival_rate,
 //
 const char      *bubble_boundaries_filename,
 //
 TChunkFrag chunkfrags[],
 // An segmented array of LIDs for the chunks.

 /* Modify the annotation,  by setting the "isrc" index into 
    chunksrc array. */
 TChunkMesg thechunks[],
 /* Output only */
 VA_TYPE(char) chunksrcs[], 
 // The character array used to store the annotations.
 FILE    *fcga,
 FILE    *fcam,
 FILE    *fp_unitig_statistics,
 FILE    *fwrn
 );

#endif /*AS_CGB_CGA_INCLUDE*/

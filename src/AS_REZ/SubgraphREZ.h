
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
/**********************************************************************

        Module:  SubgraphREZ.h

   Description:  Declaration of common data types used in SubgraphREZ.c

    Programmer:  S. Lonardi (stelo@cs.purdue.edu)

       Written:  8 July 99
 **********************************************************************/

/*********************************************************************
   CVS_ID: $Id: SubgraphREZ.h,v 1.1.1.1 2004-04-14 13:53:33 catmandew Exp $
 *********************************************************************/

#ifndef SUBGRAPH_REZ_H
#define SUBGRAPH_REZ_H

#define CHECK_NEIGHBORS 0

//
// filters (can be used in Build_Subgraph)
// 
int All_Chunk_Graph_Nodes(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int All_Chunk_Subgraph_Nodes(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Has_Coordinates_And_Path_Bit(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Has_Path_Bit(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Has_Path_Or_Visited_Bit(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Has_Coordinates(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Is_Not_Unique(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Neighborhood(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Is_UU(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Is_UR(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Is_RU(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Is_RR(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Has_Min_Coverage_Stat_Positive(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Has_Min_Coverage_Stat_Minus_Two(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Has_Min_Coverage_Stat_Minus_Five(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Has_Min_Coverage_Stat_Minus_Ten(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Has_An_Edge_To_A_Unique_In_The_Scaffold(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Has_An_Edge_To_A_Unique(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Has_No_Edges_To_Uniques_In_Another_Scaffold(ChunkInstanceT *, chunk_subgraph *, int32, int32);

int Contains_Locale_Frags(ChunkInstanceT *, chunk_subgraph *, int32, int32);

//
// CIEdge filters (can be used in Build_Subgraph)
//

int Is_Not_Bogus_Is_Overlap(CIEdgeT * edge);

int Is_Overlap(CIEdgeT * edge);

int Is_Not_Bogus_Not_Contained_Is_Overlap(CIEdgeT * edge);

int All_Edges(CIEdgeT *);

int Is_Not_Guide(CIEdgeT *);

int Is_Not_Bogus(CIEdgeT *);

int Is_Not_Bogus_Not_Contained(CIEdgeT *);

int Is_Not_Contained_Is_Overlap(CIEdgeT * edge);

int Is_Not_Bogus_Not_Contained_Not_TransChunk(CIEdgeT *);

int Is_Not_Bogus_Not_Guide(CIEdgeT *);

int Is_Not_Guide_Overlap(CIEdgeT *);

int Is_Not_Guide_Not_Overlap(CIEdgeT *);

int Is_Not_Guide_Not_Pure_Overlap(CIEdgeT * edge);

int Is_Not_Pure_Overlap(CIEdgeT * edge);

int Is_Not_Confirmed_Not_Guide(CIEdgeT *);

int Is_Not_Bogus_Not_Guide_Not_Pure_Overlap(CIEdgeT * edge);

int Is_Not_Bogus_Not_Guide_Not_Pure_Overlap_And_PathConfirmed(CIEdgeT * edge);

int Is_Not_Bogus_Not_Guide_Has_Weight_2_Or_More(CIEdgeT * edge);

int Is_Not_Bogus_Not_Guide_Is_Overlap_Not_Tandem(CIEdgeT * edge);

// subgraph
//

//
// Build_Subgraph() returns a chunk_subgraph of the ChunkGraph where
// the filtering condition is specified by filter_nodes().  Only edges
// such that filter_edge(edge)==TRUE are considered
//
chunk_subgraph * Build_Subgraph(chunk_subgraph *,
				int32,
				int32,
				int32,
				int (*)(ChunkInstanceT *, chunk_subgraph *, int32, int32),
				int (*)(CIEdgeT *));

//
// Build_Subgraph_Gap() returns a chunk_subgraph representation for
// the chunks in the <gapid>-th gap of the <sid>-th scaffold. The size
// of the array is given by
// <fill_chunks[sid].gap[gapid].num_chunks>. Only edges such that
// filter(edge)==TRUE are considered
//
// NOTE: not useful anymore
//
chunk_subgraph * Build_Subgraph_Gap(Scaffold_Fill_t *,
				    int32,
				    int32,
				    int (*)(CIEdgeT *));


//
// Build_Subgraph_Bcc() creates a subgraph of the nodes in the
// biconnected component <bcc_id>.  Only edges such that
// filter(edge)==TRUE are considered
//
chunk_subgraph * Build_Subgraph_Bcc(bcc_array *,
				    int,
				    int (*)(CIEdgeT *));


//
// Build_Subgraph_Path() returns a chunk_subgraph representation for the
// chunks that lay on a path from CIs.curr to CIs.next plus the chunks
// that lay on a path from CIs.next and CIs.curr
//
// The path is searched in the subgraph f
//
// Only edges such that filter(edge)==TRUE are considered
//
chunk_subgraph * Build_Subgraph_Path(chunk_subgraph *,
				     CIScaffoldTIterator *,
				     LengthT *,
				     float (*)(CIEdgeT *, int32),
				     float qt,
				     int (*)(CIEdgeT *),
				     int,
				     int* hops,
				     long* calls,
				     float* tooShort,
				     float* tooLong);

//
// Subgraph_Union() does the union of <s1> with <s2>. Only edges such
// that filter(edge) == TRUE are added
//
// NOTE: it assumes that there are NO duplicates (i.e., the
// intersection of <s1> and <s2> is empty). If nodes with the same
// chunk id appears in both sets they will be probably merged but the
// results are unpredictable.
//
chunk_subgraph * Subgraph_Union(chunk_subgraph *,
				chunk_subgraph *,
				int (*)(CIEdgeT *));

//
// Add_Node() adds an element to the subgraph. Only edges such that
// filter(edge) == TRUE are added
//
void Add_Node(chunk_subgraph *,
	      int32,
	      int (*)(CIEdgeT *));


//
// Delete_Node() deletes an element from the subgraph
//
void Delete_Node(chunk_subgraph *,
			int32);


//
// output all the chunk in <f> that are in the
// range [coordinate of left_cid, coordinate or right_cid]
//
// highlight the chunks that are in the subgraph s
//
// draw the edges in the subgraph s
//
void Print_Subgraph_Cam(chunk_subgraph *,
			chunk_subgraph *,
			int32,
			int32,
			int);

// The same as Print_Subgraph_Cam except with calculated coords
void Print_Subgraph_Calc_Cam(chunk_subgraph *,
			chunk_subgraph *,
			int32,
			int32,
			int,
			LengthT *);
//
// Print_Dot_Subgraph () prints a .dot representation of subgraph. It
// can be printed on a postscript file using "dot -Tps <filename.dot>
// > <filename.ps>". It consider only edge such that filter(edge) ==
// TRUE
//
void Print_Dot_Subgraph(chunk_subgraph *,
			FILE *, char *,
			int (*)(CIEdgeT *));


//
// Print_Subgraph() prints the contents of the subgraph <s> in the
// GlobalData->gwlogfp
//
void Print_Subgraph(chunk_subgraph *);


//
// Belong_To() returns a pointer to the node if the cid is in the
// graph, NULL otherwise
//
chunk_subgraph_node * Belong_To(chunk_subgraph *,
				int32);


//
// Free_Subgraph() frees all the memory associated with the chunk
// subgraph
//
void Free_Subgraph(chunk_subgraph *);

#endif

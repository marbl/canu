
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

        Module:  GapWalkerREZ.h

   Description:  Declaration of common data types used in GapWalkerREZ.c

    Programmer:  S. Lonardi (stelo@cs.purdue.edu)

       Written:  20 May 99

 **********************************************************************/

/*********************************************************************
   CVS_ID: $Id: GapWalkerREZ.h,v 1.1.1.1 2004-04-14 13:53:17 catmandew Exp $
 ********************************************************************/

#ifndef GAPWALKER_REZ_H
#define GAPWALKER_REZ_H

#define DEBUG_GAP_WALKER              3 // debugging info (level 1,2,3,4 gives debug info, 0 nothing)
#define CREATE_CAM_FILE               0 // 0 will not produce cam files of various sorts
//#define USE_EDGE_QUALITY              1 // 1 will use the quality of the edges in the Find_Greedy_Path
#define DUMP_PATHS_NOT_CROSSED        0 // 1 will create a cam file of the path not crossed
#define ALLOW_RECOMPUTE_SHORTER_PATHS 0 // 1 will re-visit nodes if they happen to be on a shorter path (bug?)




#define MAXCALLS     1000000000
// the maximum number of edges we inspect in any call to Find_Greedy_Path
#define MAXWALKCALLS 30000
// the maximum number of edges we inspect totally in a walk
#define MAXOUTDEGREE 30
// The maximal number of edges we inspect going out from a node


#define MY_UNIQUE_COLOUR       0
#define NOT_VISITED_COLOUR     1
#define NOT_PATH_COLOUR        2
#define FIRST_PATH_COLOUR      3
#define PATH_NUM_COLOURS      12
#define SP_NUM_COLOURS        32

#define FIRST_CHUNK_ID         0 // once upon a time ...  the first chunk id was 1 ...

//#define ADMISSIBLE_STDDEV    3.0 // how many stdDev we allow in Find_Greedy_Path

#define MIN_SCAFF_EDGE    1000.0 // minimum length of Scaffold edges (for InterScaffold walking)
#define MIN_GAP_SIZE       -CGW_DP_MINLEN      
//#define MAX_GAP_SIZE       1000  
// minimum gap size for considering gap walking
#define SLOP               1000  // slop we allow in checking upper and lower bound.

#if DEBUG_GAP_WALKER > 3
#define ITERATOR_VERBOSE    GRAPH_EDGE_VERBOSE
#else
#define ITERATOR_VERBOSE   FALSE
#endif

#define STR_LEN              256 
// for statically allocated strings
#define NO_QUALITY_THRESH    10000.0 
// A dummy value. Should be higher than any used quality value


// ------
// protos
// ------

// ------
// filter
// ------

//
// termination condition used in Find_Greedy_Path
//
int Stop_At_Unique(ChunkInstanceT *,
		   CDS_CID_t,
		   chunk_subgraph *);

//
// termination condition used in Find_Greedy_Path
//
int Stop_At_The_Other_End(ChunkInstanceT *,
			  CDS_CID_t,
			  chunk_subgraph *);

//
// termination condition used in Find_Greedy_Path
//
int Stop_At_Path_Bit_Set(ChunkInstanceT *,
			 CDS_CID_t,
			 chunk_subgraph *);

// ------------
// edge quality
// ------------

//
// Edge_Quality_GW() returns a real value that indicates the reliability
// of the edge mate going out from <cid>. better edges
// have lower score
//
float Edge_Quality_GW(CIEdgeT *,
		      CDS_CID_t);

//
// No_Quality() returns always 1.0
//
float No_Quality(CIEdgeT *,
		 CDS_CID_t);


// 
// The bayesian quality approach.
// if theres is only a matelink it returns a better (= lower) value than for any overlap
//
float Bayesian_Quality(CIEdgeT *edge, CDS_CID_t cid);

//
// A quality function for walking shredded BACs.
// It relies on the overlap length and the edge->quality field
float Bac_Walking_Quality(CIEdgeT *edge, CDS_CID_t cid);


// -------
// general
// -------


//
// Visit_Subgraph() finds and prints all the paths from <begin_cid> to
// <end_cid> into a "rezwalker.cam" file.
// it uses Find_Greedy_Path() as a subroutine
// 
// NOTE: this version is directional, that is if it enters a chunk
// from the A end it goes out from the B end and viceversa
//
void Visit_Subgraph(chunk_subgraph *,
		    CDS_CID_t,
		    CDS_CID_t,
		    int,
		    LengthT *,
		    float (*)(CIEdgeT *, CDS_CID_t),
		    int (*)(ChunkInstanceT *, CDS_CID_t, chunk_subgraph *));


//
// Find_Greedy_Path() searches for a path in a subgraph <s> of chunks
// from the chunk <from_cid> to the chunk <to_cid> (or to a unique,
// depending on the function terminate() that you are passing) which has
// maximal length <max_distance> (if the CHECK_DISTANCE is defined).
// The search goes in the direction specified by <end>.  At each node,
// Find_Greedy_Path() evaluates the quality of each outgoing edge, and
// pick greedly the best one that has not been visited already
//
// Notes:
//   * <level> should be set to 0 (correspond to the depth level)
//
//   * the <quality()> function takes an edge and should evaluate to a float score
//     (see Edge_Quality2() as example)
//
//   * the <terminate()> function tells the procedure when to stop
//     (see Stop_At_Unique() or Stop_At_The_Other_End() as examples)
//
//   * the path is directional, that is if enters a chunk
//     from the A end it goes out from the B end and viceversa
//
// Returns:
//   * TRUE if a path has been found
//
// Modifies:
//   * the s->table[]->path* fields to record the path(s) found (see
//     the defn in GapWalkerREZ.h)
//
int Find_Greedy_Path(chunk_subgraph *,
		     CDS_CID_t,
		     CDS_CID_t,
		     CDS_CID_t,
		     LengthT*,
		     LengthT,
		     int,
		     int,
		     float (*)(CIEdgeT *, CDS_CID_t),
		     float qualityThresh,
		     int (*)(ChunkInstanceT *, CDS_CID_t, chunk_subgraph *),
		     int,
		     int,
		     long*,
		     long,
		     float*,
		     float*);


//
// Inter_Scaffold_Gap_walker() takes a scaffold id <A>, a scaffold
// id <B>
//
// * it finds the gap length between the last chunk of <A> and the
// first chunk in <B> by looking for some scaffold edge
//
// * if none, then return FALSE
//
// * it builds the subgraph of all the chunks in the gap
//
// * it looks for a path from <A> to <B> in the subgraph such that the
//   length of the path is consistent with the gap length

int Inter_Scaffold_Gap_Walker(CDS_CID_t, CDS_CID_t,
			      float (*)(CIEdgeT *, CDS_CID_t));


//
// Intra_Scaffold_Gap_walker() takes a chunk_subgraph of the graph we
// want to analyze, a scaff_id, and two adjacents chunks in the
// scaffold, then:
//
// * it estimates the gap length (if negative return NULL) between
//   the two chunks
//
// * it isolates the chunks in the gap calling Build_Subgraph_Path()
//
// * it returns the subgraph
//
// * it checks if the chunks created by the previous step actually
//   belongs to the gap (print debug info/a .dot file of that subgraph)
//
// the two counters hops and calls are passed along to count the number of explored
// edges in gap walking and the length of the path
// the two pointers tooShort and tooLong are passed down in order
// record the best walks that are too short and too long
//
chunk_subgraph * Intra_Scaffold_Gap_Walker(chunk_subgraph *,
					   CDS_CID_t,
					   CIScaffoldTIterator * CIs,
					   float (*)(CIEdgeT *, CDS_CID_t),
					   float qualityThresh,int* hops, long* calls,
					   float *tooShort, float *tooLong);


//
// Shortest_Path() finds the shortest path(s) from <begin_cid> to any
// node in <s> by calling Dijkstra() and it prints the informations
// into a "*.sp.cam" file
//
void Shortest_Path(chunk_subgraph *,
		   CDS_CID_t,
		   CDS_CID_t,
		   float (*)(CIEdgeT *, CDS_CID_t));


//
// Dijkstra's single source shortest path (see, e.g., CLR page 527)
//
// The substantial difference is the definition of path lenght that
// we assume as the maximum edge cost on the edges of the path
// (bottleneck distance); we break ties by using the standard
// distance (sum of edges weights on the path)
//
// NOTE: this version is directional, that is if enters a chunk
// from the A end it looks for edges from the B end and viceversa
//
// The sink is used only to avoid to pick the trivial path source->sink
//
void Dijkstra(chunk_subgraph *,
	      CDS_CID_t,
	      CDS_CID_t,
	      int,
	      float (*)(CIEdgeT *, CDS_CID_t));

//
// Compute_Outdegree() computes the outdegree of each side of each
// chunk stores into a celagram file
//
void Compute_Outdegree(chunk_subgraph *);


//
// Check_Edge_Distance() computes the error of the CI edges by
// comparing them with the simulator coordinates. It will consider
// only edges such that filter(e) == TRUE
//
void Check_Edge_Distance(chunk_subgraph *,
			 int (*)(CIEdgeT *));


//
// set the isUnique flag for all chunks in the struct <fill_chunks>
// (not used anymore)
//
void Mark_Unique(Scaffold_Fill_t *);


//
// Compute_Path_Outdegree() computes the number of uniques reachable
// from each non-unique in the subgraph <s>
//
void Compute_Path_Outdegree(chunk_subgraph *);


//
// FindGapLength() finds the gap size between two chunks: handles
// the case of intra scaffolds gaps; replaces Compute_Gap_Length
//
LengthT FindGapLength( ChunkInstanceT * lchunk,
					   ChunkInstanceT * rchunk,
					   int verbose);

//
// Compute_Gap_Length() finds the gap size between two chunks: handles
// the case of intra scaffolds gaps
//
LengthT Compute_Gap_Length(ChunkInstanceT *,
			   ChunkInstanceT *,
			   int verbose);

//
//
//
int Annotate_Edges_For_Consistent_Path(chunk_subgraph * s, chunk_subgraph * f);

// ------
// output
// ------

int Count_TransChunk_Edges(chunk_subgraph *);

int Count_Edges(chunk_subgraph *);

void  Force_Increasing_Variances_One_Scaffold(int scaff_id);
#endif

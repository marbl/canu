
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

        Module:  SubGraphREZ.c

   Description:  contains routines to handle/display/create/modify
                 subgraphs of the scaffold graph

		 log info is sent to <inputFile>.gwlog

    Programmer:  S. Lonardi (stelo@cs.purdue.edu)

       Written:  7 July 99

 **********************************************************************/


static char fileID[] = "$Id: SubgraphREZ.c,v 1.1.1.1 2004-04-14 13:53:33 catmandew Exp $";


#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
//
// AS_CGW
//
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
//
// AS_REZ
//
#include "DataTypesREZ.h"
#include "UtilsREZ.h"
#include "CommonREZ.h"
#include "GapWalkerREZ.h"
#include "BccREZ.h"
#include "SubgraphREZ.h"
#include "FbacREZ.h"
//
// externs
//
extern char
* GW_Filename_Prefix;

extern int
start_index[4],
  end_index[4],
  new_end[NUM_ORIENTATIONS],
  orient[NUM_ORIENTATIONS];

//
// node filters
//


int All_Chunk_Graph_Nodes(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes in the chunk graph
  //
  return TRUE;
}



int All_Chunk_Subgraph_Nodes(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes in the chunk subgraph s
  //
  return (Belong_To(s, chunk->id) != NULL);
}



int Has_Coordinates_And_Path_Bit(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes in the range [a,b] that has the path bit set
  //
  if (s)
    if (s->table[chunk->id])
      return (((chunk->aEndCoord >= a) && (chunk->aEndCoord <= a)) ||
			  ((chunk->bEndCoord >= b) && (chunk->bEndCoord <= b)) ||
			  s->table[chunk->id]->path_bit);
    else
      return FALSE;
  else
    return FALSE;
}



int Has_Path_Bit(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that has the path bit set
  //
  if (s)
    if (s->table[chunk->id])
      return (s->table[chunk->id]->path_bit);
    else
      return FALSE;
  else
    return FALSE;
}



int Has_Path_Or_Visited_Bit(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that has the path bit set
  //
  if (s)
    if (s->table[chunk->id])
      return (s->table[chunk->id]->path_bit || s->table[chunk->id]->visited);
    else
      return FALSE;
  else
    return FALSE;
}



int Has_Coordinates(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes in the range [a,b]
  //
  return (((chunk->aEndCoord >= a) && (chunk->aEndCoord <= a)) ||
		  ((chunk->bEndCoord >= b) && (chunk->bEndCoord <= b)));
}



int Is_Not_Unique(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that are not unique
  //
  return ((chunk->type != DISCRIMINATORUNIQUECHUNK_CGW) &&
		  (chunk->type != UNIQUECHUNK_CGW));
}



int Is_UU(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that are UU
  //
  return (chunk->flags.bits.cgbType == (unsigned int)UU_CGBTYPE);
}



int Is_UR(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that are UR
  //
  return (chunk->flags.bits.cgbType == (unsigned int)UR_CGBTYPE);
}



int Is_RU(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that are RU
  //
  return (chunk->flags.bits.cgbType == (unsigned int)RU_CGBTYPE);
}



int Is_RR(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that are RR
  //
  return (chunk->flags.bits.cgbType == (unsigned int)RR_CGBTYPE);
}



int Has_Min_Coverage_Stat_Positive(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that have coverage stat at least 0.0
  //
  return (GetCoverageStat(chunk) > 0.0);
}




int Has_Min_Coverage_Stat_Minus_Two(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that have coverage stat at least -2.0
  //
  return (GetCoverageStat(chunk) > -2.0);
}



int Has_Min_Coverage_Stat_Minus_Five(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that have coverage stat at least -5.0
  //
  return (GetCoverageStat(chunk) > -5.0);
}



int Has_Min_Coverage_Stat_Minus_Ten(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that have coverage stat at least -10.0
  //
  return (GetCoverageStat(chunk) > -10.0);
}



int Has_An_Edge_To_A_Unique_In_The_Scaffold(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that have a link to another unique in the same scaffold
  //
  CIEdgeT
    * e;
  GraphEdgeIterator
    edges;
  ChunkInstanceT
    * next_chunk;
  int32
    next_chunk_cid;
  InitGraphEdgeIterator(ScaffoldGraph->RezGraph, chunk->id, 
						ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
  while((e = NextGraphEdgeIterator(&edges)) != NULL){
    assert(e != NULL);
    if (chunk->id == e->idA)
      next_chunk_cid = e->idB;
    else
      next_chunk_cid = e->idA;
    next_chunk = GetGraphNode(ScaffoldGraph->RezGraph, next_chunk_cid);
    assert(next_chunk != NULL);
    if ((next_chunk->scaffoldID == chunk->scaffoldID)  &&
		(Is_Unique(next_chunk)))
      break;
  }
  return (e != NULL);
}

int Has_No_Edges_To_Uniques_In_Another_Scaffold(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that do not have links to a unique in another scaffold
  //
  CIEdgeT
    * e;
  GraphEdgeIterator
    edges;
  ChunkInstanceT
    * next_chunk;
  int32
    next_chunk_cid;
  InitGraphEdgeIterator(ScaffoldGraph->RezGraph, chunk->id, 
						ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
  while((e = NextGraphEdgeIterator(&edges))!= NULL) {
    assert(e != NULL);
    if (chunk->id == e->idA)
      next_chunk_cid = e->idB;
    else
      next_chunk_cid = e->idA;
    next_chunk = GetGraphNode(ScaffoldGraph->RezGraph, next_chunk_cid);
    assert(next_chunk != NULL);
    if ((next_chunk->scaffoldID != chunk->scaffoldID)  &&
		(Is_Unique(next_chunk)))
      break;
  }
  return (e == NULL);
}



int Has_An_Edge_To_A_Unique(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // get all the nodes that have a link to another unique
  //
  CIEdgeT
    * e;
  GraphEdgeIterator
    edges;
  ChunkInstanceT
    * next_chunk;
  int32
    next_chunk_cid;
  InitGraphEdgeIterator(ScaffoldGraph->RezGraph, chunk->id,
						ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
  while((e = NextGraphEdgeIterator(&edges)) != NULL){
    assert(e != NULL);
    if (chunk->id == e->idA)
      next_chunk_cid = e->idB;
    else
      next_chunk_cid = e->idA;
    next_chunk = GetGraphNode(ScaffoldGraph->RezGraph, next_chunk_cid);
    assert(next_chunk != NULL);
    if ((next_chunk->scaffoldID == chunk->scaffoldID)  &&
		(Is_Unique(next_chunk)))
      break;
  }
  return (e != NULL);
}



int Neighborhood(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // returns TRUE for all chunks is in the "neighborhood" of radius one
  // the subgraph s
  //
  CIEdgeT
    * e;
  GraphEdgeIterator
    edges;
  int32
    next_chunk_cid = NULLINDEX;
  InitGraphEdgeIterator(ScaffoldGraph->RezGraph, chunk->id,
                        ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
  while((e = NextGraphEdgeIterator(&edges)) != NULL){
    assert(e != NULL);
    if (chunk->id == e->idA)
      next_chunk_cid = e->idB;
    else
      next_chunk_cid = e->idA;
    if (Belong_To(s, next_chunk_cid) != NULL)
      break;
  }
  return (Belong_To(s, next_chunk_cid) != NULL);
}

// used to build the graph of nodes from the locale set in buildLocale
// should be replaced by calls to the BAC store, which will tell us what contigs
// contain frags from the locale of interest
int Contains_Locale_Frags_new(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // return only those nodes that include frags from BACs
  //
  localeInfoT *localeInfo, *localeInfoHead;
  int i, numLocales, foundLocale = 0;
  
  numLocales = getLocalesInNode(chunk, &localeInfo, A_END, 0);
  localeInfoHead = localeInfo;
  for (i = 0; i < numLocales; i++)
  {
	if (localeInfo->localeNumber == buildLocale)
	  foundLocale = 1;
	if (i < numLocales - 1)
	  localeInfo = localeInfo->next;
  }

  // now free
  if (numLocales > 0)
  {
	localeInfo = localeInfoHead;
	for (i = 0; i < numLocales - 1; i++)
	{
	  localeInfoHead = localeInfo->next;
	  free(localeInfo);
	  localeInfo = localeInfoHead;
	}
	free (localeInfoHead);  
  }  

  numLocales = getLocalesInNode(chunk, &localeInfo, B_END, 0);
  localeInfoHead = localeInfo;
  for (i = 0; i < numLocales; i++)
  {
	if (localeInfo->localeNumber == buildLocale)
	  foundLocale = 1;
	if (i < numLocales - 1)
	  localeInfo = localeInfo->next;
  }
  
  // now free
  if (numLocales > 0)
  {
	localeInfo = localeInfoHead;
	for (i = 0; i < numLocales - 1; i++)
	{
	  localeInfoHead = localeInfo->next;
	  free(localeInfo);
	  localeInfo = localeInfoHead;
	}
	free (localeInfoHead);  
  }

  if (foundLocale)
	return 1;
  else
	return 0;

  // return (chunk->flags.bits.includesFinishedBacFragments);
}

// used to build the graph of nodes from the locale set in buildLocale
// should be replaced by calls to the BAC store, which will tell us what contigs
// contain frags from the locale of interest
int Contains_Locale_Frags(ChunkInstanceT * chunk, chunk_subgraph * s, int32 a, int32 b) {
  //
  // filter on nodes (can be used in Build_Subgraph)
  //
  // return only those nodes that include frags from BACs
  //

  return (chunk->flags.bits.includesFinishedBacFragments);
}



//
// CIEdge filters
//



int All_Edges(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return TRUE;
}



int Is_Not_Guide(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return (! edge->flags.bits.hasGuide);
}



int Is_Not_Guide_Overlap(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return ((! edge->flags.bits.hasGuide) &&
		  (isOverlapEdge(edge)));
}



int Is_Not_Guide_Not_Overlap(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return ((! edge->flags.bits.hasGuide) &&
		  (! isOverlapEdge(edge)));
}



int Is_Not_Pure_Overlap(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return (! (isOverlapEdge(edge) &&
			 (edge->edgesContributing == 1)));
}



int Is_Not_Guide_Not_Pure_Overlap(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return (! ((edge->flags.bits.hasGuide) ||
			 (isOverlapEdge(edge) && 
			  (edge->edgesContributing == 1))));
}



int Is_Not_Bogus(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return (! isProbablyBogusEdge(edge));
}

int Is_Not_Bogus_Is_Overlap(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return (! isProbablyBogusEdge(edge) && isOverlapEdge(edge) );
}

int Is_Overlap(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return ( isOverlapEdge(edge) );
}

int Is_Not_Bogus_Not_Contained(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  if  (edge->flags.bits.hasContainmentOverlap != (edge->flags.bits.aContainsB || edge->flags.bits.bContainsA)){
    fprintf (stderr, "mismatch between hasContainmentOverlap and flags.bits.xContainsY\n");
    PrintGraphEdge(stderr,ScaffoldGraph->RezGraph," ", edge, edge->idA);
  }
  return (!isProbablyBogusEdge(edge) && !(edge->flags.bits.aContainsB || edge->flags.bits.bContainsA));
  // return (!isProbablyBogusEdge(edge) && (!edge->flags.bits.hasContainmentOverlap));
}



int Is_Not_Contained_Is_Overlap(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  if  (edge->flags.bits.hasContainmentOverlap != (edge->flags.bits.aContainsB || edge->flags.bits.bContainsA)){
    fprintf (stderr, "mismatch between hasContainmentOverlap and flags.bits.xContainsY\n");
    PrintGraphEdge(stderr,ScaffoldGraph->RezGraph," ", edge, edge->idA);
  }
  return ( isOverlapEdge(edge) && !(edge->flags.bits.aContainsB || edge->flags.bits.bContainsA) );
}





int Is_Not_Bogus_Not_Contained_Is_Overlap(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  if  (edge->flags.bits.hasContainmentOverlap != (edge->flags.bits.aContainsB || edge->flags.bits.bContainsA))
    fprintf (stderr, "mismatch between hasContainmentOverlap and flags.bits.xContainsY\n");
  return ( !isProbablyBogusEdge(edge) && isOverlapEdge(edge) && !(edge->flags.bits.aContainsB || edge->flags.bits.bContainsA) );
}


int Is_Not_Bogus_Not_Contained_Not_TransChunk(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return (! isProbablyBogusEdge(edge) && (!edge->flags.bits.hasContainmentOverlap) 
		  && !isTransChunkEdge(edge));
}



int Is_Not_Bogus_Not_Guide(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return (! (edge->flags.bits.hasGuide ||
			 isProbablyBogusEdge(edge)));
}



int Is_Not_Confirmed_Not_Guide(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return (! ((edge->edgesContributing > 1) ||
			 (edge->flags.bits.hasGuide)));
}



int Is_Not_Bogus_Not_Guide_Not_Pure_Overlap(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return (! (isProbablyBogusEdge(edge) ||
			 edge->flags.bits.hasGuide ||
			 (isOverlapEdge(edge) && 
			  (edge->edgesContributing == 1))));
}


int Is_Not_Bogus_Not_Guide_Not_Pure_Overlap_And_PathConfirmed(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  // did I break a record for the longest identifier?
  //
  return (! (isProbablyBogusEdge(edge) ||
			 edge->flags.bits.hasGuide ||
			 (isOverlapEdge(edge) && 
			  (edge->edgesContributing == 1)) ||
			 (! edge->flags.bits.hasConfirmingPath &&
			  (edge->edgesContributing == 1))));
}



int Is_Not_Bogus_Not_Guide_Has_Weight_2_Or_More(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return ((edge->edgesContributing >= 2) &&
		  (! edge->flags.bits.hasGuide) &&
		  (! isProbablyBogusEdge(edge)));
}


int Is_Not_Bogus_Not_Guide_Is_Overlap_Not_Tandem(CIEdgeT * edge) {
  //
  // filter on CIEdge (can be used in Build_Subgraph)
  //
  return (! (isProbablyBogusEdge(edge) || edge->flags.bits.hasGuide) 
          && ((edge->flags.bits.hasContributingOverlap || edge->flags.bits.hasRepeatOverlap) 
			  && !edge->flags.bits.hasTandemOverlap));
  //	  && (edge->edgesContributing == 1));
}


// ------------------
// <subgraph> methods
// ------------------


#define STEPS 10

void Print_Dot_Subgraph(chunk_subgraph * s,
						FILE * wfile,
						char * label,
						int (* filter)(CIEdgeT *)) {
  //
  // Print_Dot_Subgraph () prints a .dot representation of
  // subgraph. It can be printed on a postscript file using
  // "/home/stelo/bin/dot -Tps <filename.dot> > <filename.ps>". It
  // consider only edge such that filter(edge) == TRUE
  //
  int32
    i,
    k,
    o,
    cnt,
    cid,
    eid,
    next;
  CIEdgeT
    * e;
  char
    edge_style[STR_LEN],
    node_shape[STR_LEN];
  int
    min = INT_MAX,
    max = 0,
    delta;
  ChunkInstanceT 
    * chunk;

  fprintf(wfile,
		  "digraph G {\n");
  fprintf(wfile,
		  "size = \"7.5,10\";\n");
  fprintf(wfile,
		  "ratio = fill;\n");
  fprintf(wfile,
		  "center = true;\n");
  fprintf(wfile,
		  "fontname = \"Helvetica\";\n");
  fprintf(wfile, // page = \"8.5,11\";\n
		  "label = \"%s\";\n", label);

  //
  // find the min/max coordinate
  //
  for (i = 0; i < s->size; i++) {
    cid = s->node[i].cid;
    chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);
    assert(chunk != NULL);
    if (chunk->aEndCoord > max)
      max = chunk->aEndCoord;
    if (chunk->aEndCoord < min)
      min = chunk->aEndCoord;
  }

  //
  // print the coordinates
  //
  delta = (max - min) / STEPS;

  fprintf(wfile,
		  "{\n");
  fprintf(wfile,
		  "node [shape = plaintext];\n");
  for (i = min; i < max; i += delta)
    fprintf(wfile,
			"%d -> %d;\n", i , i + delta);
  fprintf(wfile,
		  "\n} /* delta = %d */\n", delta);
  //
  // print nodes and edges
  //
  for (i = 0; i < s->size; i++) {
    cid = s->node[i].cid;
    for (cnt = 0, k = 0; k < NUM_ORIENTATIONS; k++) {
      o = orient[k];
      assert(s->table[cid] != NULL);
      for (eid = 0; eid < s->table[cid]->num_edges[o]; eid++) {
		//
		// get the edge
		//
		e = s->table[cid]->edge[o][eid];
		assert(e != NULL);
		if (! filter(e))
		  continue;
		//
		// get other end
		//
		if (cid == e->idA)
		  next = e->idB;
		else
		  next = e->idA;
		//
		// we want to draw an undirected graph, then we
		// discard half of the edges
		//
		if (e->edgesContributing > 1)
		  sprintf (edge_style, ", style = bold");
		else
                  edge_style[0] = '\0';

		if (cid < next)
		  fprintf(wfile,
				  "%d -> %d [dir = none, label =\"(%5.2f,%5.2f)\"];\n",
				  cid,
				  next,
				  e->distance.mean,
				  e->distance.variance
			);
		cnt++;
      }
    }
    //
    // output the vertex if there are outgoing edges
    //
    if (cnt) {
      ChunkInstanceT 
		* chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);
      assert(chunk != NULL);

      if (Is_Unique(chunk))
		sprintf (node_shape, "box");
      else 
		sprintf (node_shape, "circle");

      fprintf(wfile,
			  "{ rank = same; %d; %d [shape = %s];} \n",
			  min + ((chunk->aEndCoord - min) / delta) * delta,
			  cid,
			  node_shape);
    }
  }

  fprintf(wfile,
		  "}\n");
}



void Print_Subgraph(chunk_subgraph * s) {
  //
  // print the contents of the subgraph <s> in the GlobalData->gwlogfp
  //
  int32
    i;

  assert(s != NULL);

# if DEBUG_GAP_WALKER > 0
  fprintf(GlobalData->gwlogfp,
		  "Subgraph size %d, max entry %d\n",
		  s->size,
		  s->max);
  fprintf(GlobalData->gwlogfp,
		  "--- Table ---\n");
  for (i = 0; i < s->max; i++)
    fprintf(GlobalData->gwlogfp,
			"id %d ptr %d \n",
			i,
			(int)s->table[i]);
# endif

  for (i = 0; i < s->max; i++)
    if (s->table[i]) {
      assert (i == s->table[i]->cid);
#     if DEBUG_GAP_WALKER > 0
      fprintf(GlobalData->gwlogfp,
			  "chunk %4d (uid %2d), dfs %d, fin %d, vis %d, pbit %d, pid %d, par %d, dist (%5.2f, %5.2f, %5.2f)\n",
			  s->table[i]->cid,
			  s->table[i]->union_id,
			  s->table[i]->dfs_time,
			  s->table[i]->fin_time,
			  s->table[i]->visited,
			  s->table[i]->path_bit,
			  s->table[i]->path_id,
			  s->table[i]->path_parent,
			  s->table[i]->distance.mean,
			  sqrt(s->table[i]->distance.variance),
			  s->table[i]->distance.variance);
	  fprintf(GlobalData->gwlogfp,
			  "           has edges (AB_AB) %d, (AB_BA) %d, (BA_AB) %d, (BA_BA) %d\n",
			  s->table[i]->num_edges[or2num(AB_AB)],
			  s->table[i]->num_edges[or2num(AB_BA)],
			  s->table[i]->num_edges[or2num(BA_AB)],
			  s->table[i]->num_edges[or2num(BA_BA)]);     
#     endif
    }
}


#define NUM_COLORS 8

void Print_Subgraph_Cam(chunk_subgraph * f,
						chunk_subgraph * s,
						int32 left_cid,
						int32 right_cid,
						int check_coordinates_on_s) {
  //
  // output all the chunk in <f> that are in the
  // range [coordinate of left_cid, coordinate or right_cid]
  //
  // highlight the chunks that are in the subgraph s
  //
  // draw the edges in the subgraph s
  // 
  // if check_coordinates_on_s == TRUE then we draw the edges
  // only for the chunks internal to the interval
  //
  int
    color,
    k,
    o,
    id,
    i;
  CIEdgeT
    * edge;
  int32
    low,
    high,
    a_end,
    b_end,
    cid,
    other_cid;
  ChunkInstanceT
    * lchunk = GetGraphNode(ScaffoldGraph->RezGraph, left_cid),
    * rchunk = GetGraphNode(ScaffoldGraph->RezGraph, right_cid),
    * other_chunk,
    * chunk;
  char
    unique[STR_LEN],
    orientation[STR_LEN],
    filename[STR_LEN],
    * Colour[NUM_COLORS] = {
      "CFFFF00 T2 S # unique_selected",
      "C0000F0 T2 S # unique_visited",
      "CAAAA00 T2 S # unique_not_visited",
      "C00FF00 T2 S # non_unique_selected",
      "CFF8000 T2 S # non_unique_visited",
      "CA0A0FF T2 S # non_unique_non_visited",
      "C0AAAA0 T2 S # gap",
      "C00AAAA T1 S # CIEdge"};
  FILE
    * walker_cam_file;
  
  //
  // get the interval
  //
  assert(lchunk != NULL);
  assert(rchunk != NULL);
  low = min(min(lchunk->aEndCoord, lchunk->bEndCoord),
			min(rchunk->aEndCoord, rchunk->bEndCoord));
  high = max(max(lchunk->aEndCoord, lchunk->bEndCoord),
			 max(rchunk->aEndCoord, rchunk->bEndCoord));
  
  //
  // open the cam file
  //
  if (check_coordinates_on_s)
    sprintf(filename, "./cam/_%s.%d.%d.ps.cam",
			GW_Filename_Prefix, left_cid, right_cid);
  else
    sprintf(filename, "./cam/%s.%d.%d.ps.cam",
			GW_Filename_Prefix, left_cid, right_cid);
  walker_cam_file = file_open (filename, "w");
  assert(walker_cam_file != NULL);
  
  //
  // output the colors
  //
  for (i = 0; i < NUM_COLORS; i++)
    fprintf(walker_cam_file, "%d: %s\n",
			i,
			Colour[i]);

  for (i = 0; i < f->size; i++) {
    cid = f->node[i].cid;
    chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);
    assert(chunk != NULL);
    
    if (((chunk->aEndCoord >= low) && (chunk->aEndCoord <= high)) ||
		((chunk->bEndCoord >= low) && (chunk->bEndCoord <= high)) ||
		f->node[i].path_bit || f->node[i].visited) {
      
      if (Is_Unique(chunk)) {
		strcpy(unique, "unique");
		color = 0;
      } else {
		strcpy(unique, "");
		color = 3;
      }
      
      if (f->node[i].path_bit)
		color += 0;
      else
	  {
		if (f->node[i].visited)
		  color += 1;
		else
		  color += 2;
	  }

      if (chunk->aEndCoord <= chunk->bEndCoord) {
		a_end = chunk->aEndCoord;
		b_end = chunk->bEndCoord;
		strcpy(orientation, "direct");
      } else {
		a_end = chunk->bEndCoord;
		b_end = chunk->aEndCoord;
		strcpy(orientation, "reverse");
      }
      
      fprintf(walker_cam_file,
			  "%d: %d A%d %d # %s %s chunk %d scaff_id " F_CID "\n",
			  cid,
			  a_end,
			  color,
			  b_end,
			  unique,
			  orientation,
			  cid,
			  chunk->scaffoldID);
    }
  }
  
  //
  // draw the edges
  //
  for (i = 0; i < s->size; i++) {
    cid = s->node[i].cid;
    chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);
    assert(chunk != NULL);
    
    for (k = 0; k < NUM_ORIENTATIONS; k++) {
      o = orient[k];
      assert(s->table[cid] != NULL);
      for (id = 0; id < s->table[cid]->num_edges[o]; id++) {
		//
		// get the edge
		//
		edge = s->table[cid]->edge[o][id];
		assert(edge != NULL);
		//
		// get the other end
		//
		if (cid == edge->idA)
		  other_cid = edge->idB;
		else
		  other_cid = edge->idA;
	
		other_chunk = GetGraphNode(ScaffoldGraph->RezGraph, other_cid);
		assert(other_chunk != NULL);

		if (check_coordinates_on_s) {
		  //
		  // we want both end to be inside the interval
		  //
		  if ((((chunk->aEndCoord >= low) && (chunk->aEndCoord <= high)) ||
			   ((chunk->bEndCoord >= low) && (chunk->bEndCoord <= high))) &&
			  (((chunk->aEndCoord >= low) && (chunk->aEndCoord <= high)) ||
			   ((chunk->bEndCoord >= low) && (chunk->bEndCoord <= high))) &&
			  (((other_chunk->aEndCoord >= low) && (other_chunk->aEndCoord <= high)) ||
			   ((other_chunk->bEndCoord >= low) && (other_chunk->bEndCoord <= high))) &&
			  (((other_chunk->aEndCoord >= low) && (other_chunk->aEndCoord <= high)) ||
			   ((other_chunk->bEndCoord >= low) && (other_chunk->bEndCoord <= high))) &&
			  (cid < other_cid))
			fprintf(walker_cam_file, "LNK: %d %d A6 # edge distance %f, contrib %d\n",
					cid,
					other_cid,
					edge->distance.mean,
					edge->edgesContributing);
		} else {
		  if (cid < other_cid)
			fprintf(walker_cam_file, "LNK: %d %d A6 # edge distance %f, contrib %d\n",
					cid,
					other_cid,
					edge->distance.mean,
					edge->edgesContributing);
		}
      }
    }
  }
  fclose (walker_cam_file);
}

void Print_Subgraph_Calc_Cam(chunk_subgraph * f,
							 chunk_subgraph * s,
							 int32 left_cid,
							 int32 right_cid,
							 int check_coordinates_on_s,
							 LengthT *original_gap) {
  //
  // output all the chunk in <f> that are in the
  // range [coordinate of left_cid, coordinate or right_cid]
  //
  // highlight the chunks that are in the subgraph s
  //
  // draw the edges in the subgraph s
  // 
  // if check_coordinates_on_s == TRUE then we draw the edges
  // only for the chunks internal to the interval
  //
  int
    color,
    i;
  int32
    low,
    min,
    max,
    cid;
#if 0
  int k, o, id;
  CIEdgeT * edge;
  ChunkInstanceT * other_chunk;
  int32 high, other_cid;
#endif
  ChunkInstanceT
    * lchunk = GetGraphNode(ScaffoldGraph->RezGraph, left_cid),
    * rchunk = GetGraphNode(ScaffoldGraph->RezGraph, right_cid),
    * chunk;
  char
    unique[STR_LEN],
    orientation[STR_LEN],
    filename[STR_LEN],
    * Colour[NUM_COLORS] = {
      "CFFFF00 T2 S # unique_selected",
      "C0000F0 T2 S # unique_visited",
      "CAAAA00 T2 S # unique_not_visited",
      "C00FF00 T2 S # non_unique_selected",
      "CFF8000 T2 S # non_unique_visited",
      "CA0A0FF T2 S # non_unique_non_visited",
      "C0AAAA0 T2 S # gap",
      "C00AAAA T1 S # CIEdge"};
  FILE
    * walker_cam_file;
  
  assert(s != NULL);

  // get the interval
  //
  assert(lchunk != NULL);
  assert(rchunk != NULL);

  low = min( min(lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean),
			 min(rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean));
  //
  // open the cam file
  //
  if (check_coordinates_on_s)
    sprintf(filename, "./cam/_%s.%d.%d.calc.ps.cam",
			GW_Filename_Prefix, left_cid, right_cid);
  else
    sprintf(filename, "./cam/%s.%d.%d.calc.ps.cam",
			GW_Filename_Prefix, left_cid, right_cid);
  walker_cam_file = file_open (filename, "w");
  assert(walker_cam_file != NULL);
  
  //
  // output the colors
  //
  for (i = 0; i < NUM_COLORS; i++)
    fprintf(walker_cam_file, "%d: %s\n",
			i,
			Colour[i]);

  // show the original gap 
  cid = 2147483647;
  min = max(lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean);
  color = 6;
  max = min + original_gap->mean;
  strcpy(unique, "gap");
  strcpy(orientation, "none");

  fprintf(walker_cam_file,
		  "%d: %d A%d %d # %s %s chunk %d scaff_id %d\n",
		  cid,
		  min,
		  color,
		  max,
		  unique,
		  orientation,
		  cid,
		  0);

  cid--;
  fprintf(walker_cam_file,
		  "%d: %d A%d %d # %s %s chunk %d scaff_id %d\n",
		  cid,
		  (int) max(0.0, min - 3.0 * sqrt(original_gap->variance)),
		  color,
		  (int) (max + 3.0 * sqrt(original_gap->variance)),
		  unique,
		  orientation,
		  cid,
		  0);

  for (i = 0; i < s->size; i++) {
    cid = s->node[i].cid;
    chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);
    fprintf(stderr, "s->size: %d, asserting on chunk with cid: %d\n", s->size, cid);
    assert(chunk != NULL);
    fprintf(stderr, "asserted on chunk with cid: %d\n", cid);
    
    if (Is_Unique(chunk)) {
      strcpy(unique, "unique");
      color = 0;
    } else {
      strcpy(unique, "");
      color = 3;
    }
    

    if (s->node[i].path_bit)
    {
	  color += 0;
    }
    else
	{
	  if (s->node[i].visited)
		color += 1;
	  else
		color += 2;
	}
    
    assert( s->table[cid] != NULL);
    
    fprintf(stderr, "%d\n", cid);
    fprintf(stderr, "%d\n", (int) s->table[cid]->A_end.mean);
    fprintf(stderr, "A%d\n", color);
    fprintf(stderr, "%d\n", (int) s->table[cid]->B_end.mean);
    fprintf(stderr, "%s\n", unique);
    fprintf(stderr, "%s\n", orientation);
    fprintf(stderr, "chunk %d\n", cid);
    fprintf(stderr, "scaff_id " F_CID "\n", chunk->scaffoldID);
    
    if (cid != left_cid && cid != right_cid && !Is_Unique(chunk))
	{
	  if (s->table[cid]->A_end.mean < s->table[cid]->B_end.mean)
	  {
	    min = (int) s->table[cid]->A_end.mean;
	    max = (int) s->table[cid]->B_end.mean;
	    strcpy(orientation, "direct");
	  }
	  else
	  {
	    min = (int) s->table[cid]->B_end.mean;
	    max = (int) s->table[cid]->A_end.mean;
	    strcpy(orientation, "reverse");
	  }
    
	  fprintf(walker_cam_file,
			  "%d: %d A%d %d # %s %s chunk %d scaff_id %d\n",
			  cid,
			  min,
			  color,
			  max,
			  unique,
			  orientation,
			  cid,
			  0);
	}
    else
	{
	  if (chunk->offsetAEnd.mean < chunk->offsetBEnd.mean)
	  {
	    min = chunk->offsetAEnd.mean;
	    max = chunk->offsetBEnd.mean;
	    strcpy(orientation, "direct");
	  }
	  else
	  {
	    min = chunk->offsetBEnd.mean;
	    max = chunk->offsetAEnd.mean;
	    strcpy(orientation, "reverse");
	  }
    
	  fprintf(walker_cam_file,
			  "%d: %d A%d %d # %s %s chunk %d scaff_id %d\n",
			  cid,
			  min,
			  color,
			  max,
			  unique,
			  orientation,
			  cid,
			  0);
	}
  }
  
#if 0
  //
  // draw the edges
  //
  for (i = 0; i < s->size; i++) {
    cid = s->node[i].cid;
    chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);
    assert(chunk != NULL);
    
    for (k = 0; k < NUM_ORIENTATIONS; k++) {
      o = orient[k];
      assert(s->table[cid] != NULL);
      for (id = 0; id < s->table[cid]->num_edges[o]; id++) {
		//
		// get the edge
		//
		edge = s->table[cid]->edge[o][id];
		assert(edge != NULL);
		//
		// get the other end
		//
		if (cid == edge->idA)
		  other_cid = edge->idB;
		else
		  other_cid = edge->idA;
	
		other_chunk = GetGraphNode(ScaffoldGraph->RezGraph, other_cid);
		assert(other_chunk != NULL);

		if (check_coordinates_on_s) {
		  //
		  // we want both end to be inside the interval
		  //
		  if ((((chunk->aEndCoord >= low) && (chunk->aEndCoord <= high)) ||
			   ((chunk->bEndCoord >= low) && (chunk->bEndCoord <= high))) &&
			  (((chunk->aEndCoord >= low) && (chunk->aEndCoord <= high)) ||
			   ((chunk->bEndCoord >= low) && (chunk->bEndCoord <= high))) &&
			  (((other_chunk->aEndCoord >= low) && (other_chunk->aEndCoord <= high)) ||
			   ((other_chunk->bEndCoord >= low) && (other_chunk->bEndCoord <= high))) &&
			  (((other_chunk->aEndCoord >= low) && (other_chunk->aEndCoord <= high)) ||
			   ((other_chunk->bEndCoord >= low) && (other_chunk->bEndCoord <= high))) &&
			  (cid < other_cid))
			fprintf(walker_cam_file, "LNK: %d %d A6 # edge distance %f, contrib %d\n",
					cid,
					other_cid,
					edge->distance.mean,
					edge->edgesContributing);
		} else {
		  if (cid < other_cid)
			fprintf(walker_cam_file, "LNK: %d %d A6 # edge distance %f, contrib %d\n",
					cid,
					other_cid,
					edge->distance.mean,
					edge->edgesContributing);
		}
      }
    }
  }
#endif
  fclose (walker_cam_file);
}



//
// Build_Edge_Arrays() copies the edges from the chunk graph (only the
// ones in this subgraph) only if filter(edge) == TRUE
//
// returns the number of edges that have been "copied"
//
static int Build_Edge_Arrays(chunk_subgraph * s,
                             int32 cid,
                             int (* filter)(CIEdgeT *)) {
  //
  // copy the edges from the chunk graph (only the ones in this subgraph)
  // only if filter(edge) == TRUE
  //
  // returns the number of edges that have been "copied"
  //
  CIEdgeT
    * e;
  GraphEdgeIterator
    edges;
  int32
    next_chunk_cid;
  int
    k,
    o,
    outgoing_edges = 0,
    count[NUM_ORIENTATIONS];

  //
  // initialize the number of edges
  //
  for (k = 0; k < NUM_ORIENTATIONS; k++)
	s->table[cid]->num_edges[k] = 0;

  //
  // iterates all the outgoing edges
  //
  InitGraphEdgeIterator(ScaffoldGraph->RezGraph, cid,
						ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
  while((e = NextGraphEdgeIterator(&edges)) != NULL){
    assert(e != NULL);
    //
    // if bogus ... discard!
    //

# if DEBUG_GAP_WALKER > 2
	fprintf(GlobalData->gwlogfp,"raw edge between end " F_CID " and " F_CID "\n", e->idA, e->idB);
# endif

	if (e->idA == 263 || e->idB == 263)
	{
	  fprintf( GlobalData->gwlogfp, "about to filter edge between " F_CID " and " F_CID "\n", e->idA, e->idB);
	  if (! filter(e))
		fprintf( GlobalData->gwlogfp, "edge between " F_CID " and " F_CID " will be filtered out\n", e->idA, e->idB);
	}
	
    if (! filter(e))
      continue;

	if (e->idA == 263 || e->idB == 263)
	  fprintf( GlobalData->gwlogfp, "after filtering edge between " F_CID " and " F_CID "\n", e->idA, e->idB);
	
    //
    // get the other end
    //
    if (cid == e->idA)
      next_chunk_cid = e->idB;
    else
      next_chunk_cid = e->idA;

    //
    // check if the next chunk is in the subgraph
    //
    if (Belong_To(s, next_chunk_cid)) {
      assert(s->table[next_chunk_cid] != NULL);
      s->table[cid]->num_edges[or2num(GetEdgeOrientationWRT(e, cid))]++;
#     if DEBUG_GAP_WALKER > 3
      fprintf(GlobalData->gwlogfp,
			  "%s\n",
			  Orientation_As_String(GetEdgeOrientationWRT(e, cid)));
#     endif
    }
  }
  
  //
  // allocate the edges
  //
  for (k = 0; k < NUM_ORIENTATIONS; k++)
    if (s->table[cid]->num_edges[k]) {
      outgoing_edges += s->table[cid]->num_edges[k];
      s->table[cid]->edge[k] = (CIEdgeT * *)safe_calloc(s->table[cid]->num_edges[k],
														sizeof(CIEdgeT *));
    } else
      s->table[cid]->edge[k] = (CIEdgeT * *)NULL;

  
# if DEBUG_GAP_WALKER > 2
  fprintf(GlobalData->gwlogfp,"Chunk %d has %d edges, of which\n%d (AB_AB)\n%d (AB_BA)\n%d (BA_AB)\n%d (BA_BA)\n",
		  cid,
		  outgoing_edges,
		  s->table[cid]->num_edges[or2num(AB_AB)],
		  s->table[cid]->num_edges[or2num(AB_BA)],
		  s->table[cid]->num_edges[or2num(BA_AB)],
		  s->table[cid]->num_edges[or2num(BA_BA)]);
# endif

  //
  // zero the counts
  //
  for (k=0; k < NUM_ORIENTATIONS; k++)
    count[k] = 0;
  
  //
  // store the edges
  //
  InitGraphEdgeIterator(ScaffoldGraph->RezGraph, cid, 
						ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
  while((e = NextGraphEdgeIterator(&edges)) != NULL){
    assert(e != NULL);
    //
    // if bogus ... discard!
    //
# if DEBUG_GAP_WALKER > 5
	PrintGraphEdge(GlobalData->gwlogfp,ScaffoldGraph->RezGraph,"LOOKING at edge  ", e, e->idA);	
#endif
    if (! filter(e))
	{
# if DEBUG_GAP_WALKER > 5
	  PrintGraphEdge(GlobalData->gwlogfp,ScaffoldGraph->RezGraph,"filtered graph edge with respect to idA ", e, e->idA);	
#endif
	  continue;
	}

    //
    // get the other end
    //
    if (cid == e->idA)
      next_chunk_cid = e->idB;
    else
      next_chunk_cid = e->idA;

    //
    // check if the next chunk is in the subgraph
    //
    if (Belong_To(s, next_chunk_cid)) {
      assert(s->table[next_chunk_cid] != NULL);
      o = or2num(GetEdgeOrientationWRT(e,cid));
      //
      // we are keeping pointer to the edges in the chunk graph
      // WARNING: if the chunk graph is reallocated we will end
      // up with dangling pointers!
      //
      s->table[cid]->edge[o][count[o]] = e;
      count[o]++;
    }
  }
  return outgoing_edges;
}



chunk_subgraph * Subgraph_Union(chunk_subgraph * s1,
                                chunk_subgraph * s2,
                                int (* filter)(CIEdgeT *)) {
  //
  // Subgraph_Union() does the union of <s1> with <s2>. Only edges such
  // that filter(edge) == TRUE are added
  //
  // NOTE: it assumes that there are NO duplicates (i.e., the
  // intersection of <s1> and <s2> is empty). If nodes with the same
  // chunk id appears in both sets they will be probably merged but the
  // results are unpredictable.
  //
  int
    i;
  chunk_subgraph
    * s = (chunk_subgraph *)safe_calloc(1, sizeof(chunk_subgraph));

  assert(s1 != NULL);
  assert(s2 != NULL);

  //
  // compute the max entry in the table and
  // allocate the subgraph and the table
  //
  s->size = s1->size + s2->size;
  s->max = max(s1->max, s2->max);
  s->table = (chunk_subgraph_node * *)safe_calloc(s->max,
												  sizeof(chunk_subgraph_node *));
  s->node = (chunk_subgraph_node *)safe_calloc(s->size,
											   sizeof(chunk_subgraph_node));

  //
  // build the subgraph and the table
  //
  for (i = 0; i < s1->size; i++) {
    s->node[i] = s1->node[i]; 
    s->table[s->node[i].cid] = &(s->node[i]);
    s->node[i].union_id = s1->node[i].union_id;
  }
  for (; i < s1->size + s2->size; i++) {
    s->node[i] = s2->node[i - s1->size];
    s->table[s->node[i].cid] = &(s->node[i]);
    s->node[i].union_id = 2 * s2->node[i].union_id + 1;
  }

  //
  // now copy the edges from the chunk graph (only the ones in this
  // subgraph) merging the two old sets does not work because we will
  // miss the edges between s1 and s2
  //
  for (i = 0;  i < s->size; i++)
    Build_Edge_Arrays(s, s->node[i].cid, filter);

  return s;
}



void Free_Subgraph (chunk_subgraph * s) {
  //
  // Free all the memory associated with the chunk subgraph
  //
  int
    i, k;
  assert(s != NULL);
  for (i = FIRST_CHUNK_ID; i < s->size; i++)
    for (k = 0; k < NUM_ORIENTATIONS; k++)
      if (s->node[i].num_edges[k])
		free(s->node[i].edge[k]);
  free(s->node);
  free(s->table);
  free(s);
}



chunk_subgraph_node * Belong_To(chunk_subgraph * s,
                                int32 cid) {
  //
  // returns a pointer to the node if the cid is in the graph, NULL
  // otherwise
  //
  assert(s != NULL);
  if (cid < s->max)
    return s->table[cid];
  else
    return NULL;
}



void Add_Node(chunk_subgraph * s,
			  int32 cid,
			  int (* filter)(CIEdgeT *)) {
  //
  // add an element to the subgraph
  // Only edges such that filter(edge) == TRUE are added
  //
  int
    i,
    k,
    old_size;
  chunk_subgraph_node * *
    old_table = NULL;
  chunk_subgraph_node *
    old_node;

  assert(s != NULL);

  //
  // check if already there
  //
  if (Belong_To(s,cid)) {
#   if DEBUG_GAP_WALKER > 2
    fprintf(GlobalData->gwlogfp,
			"-> %d already in the graph\n",
			cid);
#   endif
    return;
  }

# if DEBUG_GAP_WALKER > 2
  fprintf(GlobalData->gwlogfp,
		  "-> %d need to be added\n",
		  cid);
# endif

  if (cid >= s->max) {
    s->max = cid + 1;
#   if DEBUG_GAP_WALKER > 2
    fprintf(GlobalData->gwlogfp,
			"-> realloc the table to %d entries\n",
			s->max);
#   endif

    //
    // realloc the table to the new size
    //
    old_table = s->table;
    s->table = (chunk_subgraph_node * *)safe_calloc(s->max,
													sizeof(chunk_subgraph_node *));
  }

  old_size = s->size;
  s->size++;
 
# if DEBUG_GAP_WALKER > 2
  fprintf(GlobalData->gwlogfp,
		  "-> realloc the nodes to %d entries\n",
		  s->size);
# endif

  //
  // realloc the subgraph to the new size
  //
  old_node = s->node;
  s->node = (chunk_subgraph_node *)safe_calloc(s->size,
											   sizeof(chunk_subgraph_node));
  for (i = 0; i < old_size; i++) {
    s->node[i] = old_node[i];
    s->table[s->node[i].cid] = &(s->node[i]);
  }

  //
  // store the new entry
  //
  s->node[old_size].cid = cid;
  s->table[cid] = &(s->node[old_size]);
  
  //
  // now add the edges
  //
  Build_Edge_Arrays(s, cid, filter);

  //
  // free the old memory
  //
  for (i = FIRST_CHUNK_ID; i < old_size; i++)
    for (k = 0; k < NUM_ORIENTATIONS; k++)
      if (old_node[i].num_edges[k])
		free(old_node[i].edge[k]);
  free(old_node);
  free(old_table);
}



void Delete_Node(chunk_subgraph * s, int32 cid) {
  //
  // delete an element from the subgraph
  //
  int k;

  assert(s != NULL);

  //
  // if not there, return
  //
  if (! Belong_To(s, cid))
    return;

  //
  // deallocate the space for the edges
  //
  for (k=0; k < NUM_ORIENTATIONS; k++)
    if (s->table[cid]->edge[k])
      free(s->table[cid]->edge[k]);

  //
  // remove the entry from the table
  //
  s->table[cid] = NULL;
}



chunk_subgraph * Build_Subgraph_Bcc(bcc_array * b,
                                    int bcc_id,
                                    int (* filter)(CIEdgeT *)) {
  //
  // create a subgraph of the nodes in the biconnected component
  // <bcc_id>
  // Only edges such that filter(edge)==TRUE are considered
  //
  int
    i,
    num_edges = 0;
  int32
    cid;
  chunk_subgraph
    * s = (chunk_subgraph *)safe_calloc(1, sizeof(chunk_subgraph));

  //
  // size of the subgraph
  //
  s->size = b[bcc_id].size;

  //
  // size of the table
  //
  s->max = 0;
  for (i = 0; i < s->size; i++) {
    cid = b[bcc_id].cid[i];
    if (cid > s->max)
      s->max = cid;
  }
  s->max++;

  //
  // no elements?
  //
  if ((s->size == 0)||(s->max == 0))
    return s;

  //
  // alloc the subgraph
  //
  s->node = (chunk_subgraph_node *)safe_calloc(s->size, sizeof(chunk_subgraph_node));

  //
  // alloc the table
  //
  s->table = (chunk_subgraph_node * *)safe_calloc(s->max, sizeof(chunk_subgraph_node *));

  //
  // build the subgraph and the table
  //
  for (i = 0; i < s->size; i++) {
    cid = b[bcc_id].cid[i];
    s->node[i].cid = cid;
    s->table[cid] = &(s->node[i]);
  }

  //
  // copy the edges from the chunk graph (only the ones in this subgraph)
  //
  for (i = 0; i < s->size; i++)
    num_edges += Build_Edge_Arrays(s, s->node[i].cid, filter);

# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
		  "* Subgraph vertices %d, edges %d\n",
		  s->size, num_edges);
# endif

  return s;
}



chunk_subgraph * Build_Subgraph_Gap(Scaffold_Fill_t * fill_chunks,
                                    int32 sid,
                                    int32 gapid,
                                    int (* filter)(CIEdgeT *)) {
  //
  // Build_Subgraph_Gap[() returns a chunk_subgraph representation for
  // the chunks in the <gapid>-th gap of the <sid>-th scaffold. The
  // size of the array is given by
  // <fill_chunks[sid].gap[gapid].num_chunks>. Only edges such that
  // filter(edge)==TRUE are considered
  //
  // NOTE: not useful anymore
  //
  int
    k,
    num_edges = 0;
  chunk_subgraph
    * s = (chunk_subgraph *)safe_calloc(1, sizeof(chunk_subgraph));
  Gap_Fill_t
    * fc = &(fill_chunks[sid].gap[gapid]);
  CIScaffoldT
    * scaff = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, sid);

  assert(fc != NULL);
  assert(fill_chunks != NULL);
  assert(sid < GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds));
  assert(gapid < fill_chunks[sid].num_gaps);
  assert(scaff != NULL);
  assert(! isDeadCIScaffoldT(scaff));
  assert(scaff->type == REAL_SCAFFOLD);
  
# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
		  " Gap %3d:  (%5.2f, stdDev %5.0f)  (%5.2f, stdDev %5.0f)\n",
		  gapid,
		  fc->start.mean,
		  sqrt(fc->start.variance),
		  fc->end.mean,
		  sqrt(fc->end.variance));
# endif

  //
  // size of the subgraph
  //
  s->size = fill_chunks[sid].gap[gapid].num_chunks;

  //
  // size of the table
  //
  s->max = 0;
  for (k = 0;  k < s->size; k++) 
    if (fc->chunk[k].chunk_id > s->max)
      s->max = fc->chunk[k].chunk_id;
  s->max++;

  //
  // no elements?
  //
  if ((s->size == 0)||(s->max == 0))
    return s;

  //
  // alloc the subgraph
  //
  s->node = (chunk_subgraph_node *)safe_calloc(s->size, sizeof(chunk_subgraph_node));

  //
  // alloc the table
  //
  s->table = (chunk_subgraph_node * *)safe_calloc(s->max, sizeof(chunk_subgraph_node *));

  //
  // build the subgraph and the table
  //
  for (k = 0; k < s->size; k++) { 
    s->node[k].cid = fc->chunk[k].chunk_id;
    s->node[k].visited = FALSE;
    s->node[k].path_id = 0;
    s->node[k].path_bit = FALSE;
    s->table[fc->chunk[k].chunk_id] = &(s->node[k]);
  }

  //
  // copy the edges from the chunk graph (only the ones in this subgraph)
  //
  for (k = 0;  k < s->size; k++)
    num_edges += Build_Edge_Arrays(s, s->node[k].cid, filter);

# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
		  "* Subgraph vertices %d, edges %d\n",
		  s->size, num_edges);
# endif

  return s;
}



chunk_subgraph * Build_Subgraph_Path(chunk_subgraph * f0,
                                     CIScaffoldTIterator * CIs,
                                     LengthT * max,
                                     float (* quality)(CIEdgeT *, int32),
                                     float qualityThresh,
                                     int (* filter)(CIEdgeT *),
                                     int verbose,
                                     int *hops,
                                     long* calls,
                                     float *tooShort,
                                     float *tooLong) {
  //
  // Build_Subgraph_Path() returns a chunk_subgraph representation for the
  // chunks that lie on a path from CIs->curr to CIs->next plus the chunks
  //
  // The path is searched in the subgraph f
  //
  // Only edges such that filter(edge)==TRUE are considered
  //
  //
  int32
    left_cid,
    right_cid;
  chunk_subgraph
    * s;
  LengthT
    dist = {0.0, 0.0};
  int
    i,
    c = 0,
    from = ALL_END;

  ChunkInstanceT
    * lc,
    * rc;


  assert(f0 != NULL);
  assert(CIs != NULL);

  left_cid = CIs->curr;
  right_cid = CIs->next;
  //
  // clear everything
  //
  Clear_All_Path_Bit(f0);

  //
  // initialize the distances
  //
  for (i = 0; i < f0->size; i++) {
    f0->node[i].distance.mean = FLT_MAX;
    f0->node[i].distance.variance = FLT_MAX;
  }

  //
  // find all paths from left_cid to right_cid
  //

  Set_Path_Bit(f0, right_cid);  // we look for this to know we've arrived

  rc = GetGraphNode(ScaffoldGraph->RezGraph,right_cid);
  lc = GetGraphNode(ScaffoldGraph->RezGraph,left_cid);

  if( (abs(rc->offsetAEnd.mean-lc->offsetAEnd.mean)) >  
      (abs(rc->offsetAEnd.mean-lc->offsetBEnd.mean)) )
    from = B_END;
  else
    from = A_END;

  if((*hops = Find_Greedy_Path(f0, left_cid, left_cid, right_cid, max, dist,
                               0,from, quality, qualityThresh,
                               Stop_At_Path_Bit_Set, c, verbose,
                               calls,MAXWALKCALLS,tooShort,tooLong)) != 0)
    c = 1;

  // if we found no other walk than the trivial walk we take it
  if( lc->flags.bits.walkedTrivial && c == 0 )
    //    c = 1;

# if DEBUG_GAP_WALKER > -1
    if( *calls > MAXWALKCALLS )
      fprintf(GlobalData->gwlogfp,
              "\nOOO Exceeded maximal number of explored edges = %ld (%d allowed)\n",
              *calls,MAXWALKCALLS);
  fprintf(GlobalData->gwlogfp,
          "\n*** Number of explored edges = %ld (%d allowed)\n",
          *calls,MAXWALKCALLS);
#endif


# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp, "PATHS found in standard walking from %d to %d: %d with %d hops\n", left_cid, right_cid, c,*hops);
#endif

  Set_Path_Bit(f0, left_cid);

  //
  // clear everything
  //
  dist.mean = 0.0;
  dist.variance = 0.0;

  //
  // initialize the distances
  //
  for (i = 0; i < f0->size; i++) {
    f0->node[i].distance.mean = FLT_MAX;
    f0->node[i].distance.variance = FLT_MAX;
  }

# if CREATE_CAM_FILE > 0
  //
  // build a subgraph of f of all the chunks that
  // have the path or visited bit set
  //
  if (c > 0) {
    s = Build_Subgraph(f0, -1, -1, 0,
					   // Has_Path_Bit,
					   Has_Path_Or_Visited_Bit,
					   filter);
    Print_Subgraph_Cam(f0, s, left_cid, right_cid, FALSE);
    Free_Subgraph(s);
  }
# endif
  
  // clean up for next time through
  Clear_All_Visited_Bit(f0);

  //
  // now return a subgraph of f of all the chunks that
  // have the path bit set
  //
  if (c > 0) {
    s = Build_Subgraph(f0, -1, -1, 1,
					   Has_Path_Bit,
					   filter);
  } else 
    return NULL;

  //
  // and return
  //
  return s;
}




chunk_subgraph * Build_Subgraph(chunk_subgraph * f,
                                int32 begin_pos,
                                int32 end_pos,
                                int32 copy_best_edges,
                                int (* filter_nodes)(ChunkInstanceT *, chunk_subgraph *, int32, int32),
                                int (* filter_edge)(CIEdgeT *)) {
  //
  // Build_Subgraph() returns a chunk_subgraph of the chunk_subgraph <f> where
  // the filtering condition is specified by filter_nodes(). If <f>==NULL, it will
  // consider all the chunks in the ChunkGraph Only edges
  // such that filter_edge(edge)==TRUE are considered
  //
  int
    i,
    j,
    min,
    max,
    num_edges = 0;
  chunk_subgraph
    * s = (chunk_subgraph *)safe_calloc(1, sizeof(chunk_subgraph));
  ChunkInstanceT  
    * chunk;

  if (f) {
    min = 0;
    max = f->max;
  } else {
    min = FIRST_CHUNK_ID;
    max = GetNumGraphNodes(ScaffoldGraph->RezGraph);
  }

  //
  // compute the size of the subgraph and the size of the table
  // 
  s->size = 0;
  s->max = 0;
  for (i = min; i < max; i++) {
	chunk = GetGraphNode(ScaffoldGraph->RezGraph, i);
	assert(chunk != NULL);
	if (filter_nodes(chunk, f, begin_pos, end_pos) && !chunk->flags.bits.isDead) {
	  s->size++;
	  if (i > s->max)
		s->max = i;
	}
  }
  s->max++;
	
  //
  // if there are no elements in [begin,end], return
  //
  if ((s->size == 0)||(s->max == 0))
    return s;

  //
  // alloc the subgraph
  //
  s->node = (chunk_subgraph_node *)safe_calloc(s->size, sizeof(chunk_subgraph_node));

  //
  // alloc the table
  //
  s->table = (chunk_subgraph_node * *)safe_calloc(s->max, sizeof(chunk_subgraph_node *));

  //
  // build the table and the subgraph
  //
  j = 0;
  for (i = min; i < max; i++) {
    chunk = GetGraphNode(ScaffoldGraph->RezGraph, i);
    assert(chunk != NULL);
    if (filter_nodes(chunk, f, begin_pos, end_pos) && !chunk->flags.bits.isDead) 
	{
      s->node[j].cid = i;
      s->node[j].visited = FALSE;
      s->node[j].path_id = 0;
      s->node[j].path_bit = FALSE;
      s->table[i] = &(s->node[j]);
      if (copy_best_edges)
      {
		s->table[i]->best_edge = f->table[i]->best_edge;
		s->table[i]->path_bit = f->table[i]->path_bit;
      }
      j++;
    } 
	else // if (i < s->max) 
	{
	  s->table[i] = (chunk_subgraph_node *)NULL;
	}
  }
  
  //
  // copy the edges from the chunk graph (only the ones in this subgraph)
  //
  for (i = 0; i < s->size; i++)
    num_edges += Build_Edge_Arrays(s, s->node[i].cid, filter_edge);

# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
		  "* Subgraph vertices %d, edges %d\n",
		  s->size, num_edges);
# endif

  return s;
}


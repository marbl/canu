
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

        Module:  BccREZ.c

   Description:  code for computation/handling/display of biconnected
                 components and DFS

		 log info is sent to <inputFile>.gwlog

    Programmer:  S. Lonardi (stelo@cs.purdue.edu)

       Written:  30 June 99

 **********************************************************************/

static char fileID[] = "$Id: BccREZ.c,v 1.4 2005-03-22 19:49:21 jason_miller Exp $";

#include <stdio.h>
#include <assert.h>
#include <math.h>
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
#include "ConsistencyChecksREZ.h"
#include "GapWalkerREZ.h"
#include "BccREZ.h"
#include "SubgraphREZ.h"
#include "GWDriversREZ.h"

//
// static vars
//
static int
  time_stamp = 0,
  bcc_time_stamp = 1;

//
// externs
//
extern char
  * GW_Filename_Prefix;

extern int
  orient[NUM_ORIENTATIONS];


//
// code
//


//
// The function that compares two bcc prop
//
static int Compare_Bcc_Prop(const void *b1,
			    const void *b2) {
  //
  // The function that compares two bcc prop
  //
  bcc_prop *bb1 = (bcc_prop *)b1;
  bcc_prop *bb2 = (bcc_prop *)b2;

  assert((bb1 != NULL) && (bb2 != NULL));

  if (bb1->unique_prop < bb2->unique_prop)
    return 1;
  else if (bb1->unique_prop > bb2->unique_prop)
    return -1;
  else
    return 0;
}



void Print_Cam_BCC(bcc_array * b,
		   chunk_subgraph * s) {
  //
  // print all the chunks with their BCC component
  //
  int32
    a_end,
    b_end,
    c_id,
    id; 
  int
    i,
    bcc,
    color;
  ChunkInstanceT  
    * chunk;
  FILE
    * walker_cam_file;
  char
    * Colour[BCC_NUM_COLOURS] = {
      "CAAA434 T1 S # bcc links",
      "C5555AA T2 S # too_small_bcc",  //1
      "C0000FF T2 S # bcc_mod00",
      "C00FF00 T2 S # bcc_mod01", //3
      "CFF0000 T2 S # bcc_mod02",
      "C55AAAA T2 S # bcc_mod03", //5
      "C00AA00 T2 S # bcc_mod04",
      "CAA0000 T2 S # bcc_mod05", //7
      "CFF00FF T2 S # bcc_mod06",
      "CFFFF00 T2 S # bcc_mod07", //9
      "C00FFFF T2 S # bcc_mod08",
      "CAA00AA T2 S # bcc_mod09", //11
      "CAAAA00 T2 S # bcc_mod10",
      "C00AAAA T2 S # bcc_mod11", //13
      "CAAAAAA T4 S # unq_too_small_bcc",
      "C0000FF T4 S # unq_bcc_mod00",
      "C00FF00 T4 S # unq_bcc_mod01",
      "CFF0000 T4 S # unq_bcc_mod02",
      "C55AAAA T4 S # unq_bcc_mod03",
      "C00AA00 T4 S # unq_bcc_mod04",
      "CAA0000 T4 S # unq_bcc_mod05",
      "CFF00FF T4 S # unq_bcc_mod06",
      "CFFFF00 T4 S # unq_bcc_mod07",
      "C00FFFF T4 S # unq_bcc_mod08",
      "CAA00AA T4 S # unq_bcc_mod09",
      "CAAAA00 T4 S # unq_bcc_mod10",
      "C00AAAA T4 S # unq_bcc_mod11"};
  char
    orient[STR_LEN],
    unique[STR_LEN],
    filename[256];

  //
  // open the graphic file
  //
  sprintf(filename, "%s.bcc.cam", GW_Filename_Prefix);
  walker_cam_file = file_open (filename, "w");
  assert(walker_cam_file != NULL);

  // 
  // we use the "done" bit to avoid printing two times
  // the same chunk (celamy does not like it :-)
  //
  Clear_All_Done_Bit(s);

  //
  // output the colors
  //
  for (i = 0; i < BCC_NUM_COLOURS; i ++)
    fprintf(walker_cam_file, "%d: %s\n",
	    i,
	    Colour[i]);

  //
  // scan all the bcc
  //
  for (bcc = 0; bcc < bcc_time_stamp; bcc++) {
    
    //
    // scan the chunks in the bcc
    //
    for (c_id = 0; c_id < b[bcc].size; c_id++) {
      //
      // get the id number
      //
      id = b[bcc].cid[c_id];
      //
      // get this node
      //
      chunk = GetGraphNode(ScaffoldGraph->RezGraph, id);
      assert(chunk != NULL);

      //
      // assign the color (if the bcc is too small
      // give it a grey color)
      // 
      if (b[bcc].size < BCC_MIN_SIZE)
	color = 1;
      else
	color = (bcc % 12) + 2;
      
      if (Is_Unique(chunk)) {
	strcpy(unique, "unique");
	color += 13;
      } else
	strcpy(unique, "");
      
      if (chunk->aEndCoord <= chunk->bEndCoord) {
	a_end = chunk->aEndCoord;
	b_end = chunk->bEndCoord;
	strcpy(orient, "direct");
      } else {
	a_end = chunk->bEndCoord;
        b_end = chunk->aEndCoord;
	strcpy(orient, "reverse");
      }

      assert(s->table[id] != NULL);
      if (! s->table[id]->done) {
	fprintf (walker_cam_file,
		 "%d: %d A%d %d # %s %s chunk %d bcc (",
		 id,
		 a_end,
		 color,
		 b_end,
		 unique,
		 orient,
		 id);
      for (i = 0; i < s->table[id]->colors; i++)
        fprintf (walker_cam_file, " %d",
		 s->table[id]->bcc_id[i]);
      fprintf (walker_cam_file,")\n");
      }
      Set_Done_Bit(s, id);
    }

#   if 0
    //
    // print the link informations (we print
    // 30 LNK at the time due to celamy limits)
    //
    if (b[bcc].size >= BCC_MIN_SIZE) {
      int
	s = b[bcc].size / MAX_LNK_CHUNKS,
	i,
	j;
      c_id = 0;
      for (j = 0; j < s; j++) {
	fprintf(walker_cam_file,
		"LNK: ");
	for (i = 0; i < MAX_LNK_CHUNKS; i++) {
	  fprintf(walker_cam_file,
		  "%d ",
		  b[bcc].cid[c_id]);
	  c_id++;
	}
	fprintf(walker_cam_file,
		"A0 # bcc %d\n",
		bcc);
      }
      if (c_id < b[bcc].size) {
	fprintf(walker_cam_file,
		"LNK: ");
	for (; c_id < b[bcc].size; c_id++)
	  fprintf(walker_cam_file,
		  "%d ",
		  b[bcc].cid[c_id]);
	fprintf(walker_cam_file,
		"# bcc %d\n",
		bcc);
      }
    } 
#   endif
  }
  //
  // good boys close the files
  //
  fclose (walker_cam_file);
}



void Print_Cam_BCC_Links(bcc_array * b) {
  //
  // print the links within their BCC component
  //
  int
    bcc,
    c_id;

  for (bcc = 0; bcc < bcc_time_stamp; bcc++) {
    if (b[bcc].size >= BCC_MIN_SIZE) {
      //
      // print the link informations (we print
      // 30 LNK at the time due to celamy limits)
      //
      int
	i,
	j,
	s = b[bcc].size / MAX_LNK_CHUNKS;
      c_id = 0;
      for (j = 0; j < s; j++) {
	fprintf(Cam_File,
		"LNK: ");
	for (i = 0; i < MAX_LNK_CHUNKS; i++) {
	  fprintf(Cam_File,
		  "%dCHUNKREZ ",
		  b[bcc].cid[c_id]);
	  c_id++;
	}
	fprintf(Cam_File,
		"# bcc %d\n",
		bcc);
      }
      if (c_id < b[bcc].size) {
	fprintf(Cam_File,
		"LNK: ");
	for (; c_id < b[bcc].size; c_id++)
	  fprintf(Cam_File,
		  "%dCHUNKREZ ",
		  b[bcc].cid[c_id]);
	fprintf(Cam_File,
		"# bcc %d\n",
		bcc);
      }
    }
  }
}



//
// do a DFS of the subgraph from the <this> node
//
// assign start_time dfs, end_time, bcc_id (see CLR)
//
// avoid edges such that filter(e) == FALSE
//
void DFS_Visit(chunk_subgraph * s,
	       int32 this,
	       int32 parent,
	       nodes_stack * ns,
	       int level,
	       int (* filter)(CIEdgeT *)) {
  //
  // do a DFS of the subgraph
  //
  // assign start_time dfs, end_time, bcc_id (see CLR)
  // avoid edges such that filter(e) == FALSE
  //
  CIEdgeT
    * e;
  ChunkInstanceT
    * this_chunk = GetGraphNode(ScaffoldGraph->RezGraph, this),
    * to;
  int32
    id,
    tid;
  int
#   if DEBUG_GAP_WALKER > 1
    i,
#   endif
    k,
    o,
    min_dfs,                // keeps the minimum of the dfs of this node and the
                            // dfs-es of all the nodes that have a back edges from this
    popped,
    lowest_time = INT_MAX;  // keeps the minimum of the low value for all the descendents
                            // of this node

  assert(s != NULL);
  assert(this_chunk != NULL);

  //
  // mark this node as visited and assign the DFS number furthermore,
  // push the cid of this node on the stack
  //
  Set_Visited_Bit(s, this);
  s->table[this]->dfs_time = time_stamp;
  min_dfs = time_stamp;
  time_stamp++;
  Push_Node(ns, this);

# if DEBUG_GAP_WALKER > 2
  fprintf(GlobalData->gwlogfp,
	  "\n");
  for (i = 0; i < level; i++) 
    fprintf(GlobalData->gwlogfp,
	    " ");
  fprintf(GlobalData->gwlogfp,
	  "* node %d (dfs %d, stack %d, top %d)\n",
	  this,
	  s->table[this]->dfs_time,
	  ns->top,
	  ns->nodes[ns->top - 1]);
# endif

  //
  // scan the outgoing edges
  //
  for (k = 0; k < NUM_ORIENTATIONS; k++) {
    o = orient[k];
    for (id = 0; id < s->table[this]->num_edges[o]; id++) {
      //
      // get the edge
      //
      e = s->table[this]->edge[o][id];
      assert(e != NULL);

      //
      // check the edge
      //
      if (! filter(e))
	continue;

      //
      // get the "other" node
      //
      if (this == e->idA)
	tid = e->idB;
      else
	tid = e->idA;
      to = GetGraphNode(ScaffoldGraph->RezGraph, tid);
      assert(to != NULL);

#     if DEBUG_GAP_WALKER > 2
      for (i = 0; i < level; i++) 
	fprintf(GlobalData->gwlogfp," ");
      fprintf(GlobalData->gwlogfp,
	      "} edge from %d to %d orientation %s\n",
	      this,
	      tid,
	      Orientation_As_String(GetEdgeOrientationWRT(e,this)));
#     endif

      //
      // discard direct cycles with the parent
      //
      if (tid == parent)
	continue;

      //
      // if the node has not been visited (i.e., white in the
      // terminology of CLR)
      //
      if (! s->table[tid]->visited) {
	//
	// assign the this in the dfs tree
	//
	s->table[tid]->dfs_parent = this;

	//
	// visit the descendents
	//
	DFS_Visit(s, tid, this, ns, level + 1, filter);

	//
	// compute the lowest low_time of the children
	//
	if (s->table[tid]->low_time < lowest_time)
	  lowest_time = s->table[tid]->low_time;

#       if DEBUG_GAP_WALKER > 2
	for (i = 0; i < level; i++) 
	  fprintf(GlobalData->gwlogfp," ");
	fprintf(GlobalData->gwlogfp,
		"] child_low %d, dfs %d\n",
		s->table[tid]->low_time,
		s->table[this]->dfs_time);
#       endif
	//
	// check if this node is an articulation point
	// for the subtree we have visited
	//
	if (s->table[tid]->low_time >= s->table[this]->dfs_time) {
#         if DEBUG_GAP_WALKER > 2
	  for (i = 0; i < level; i++) 
	    fprintf(GlobalData->gwlogfp," ");
	  fprintf(GlobalData->gwlogfp,
		  "> articulation point (this %d, dfs %d, child_low %d)\n",
		  this,
		  s->table[this]->dfs_time,
		  s->table[tid]->low_time);
#         endif
	  //
	  // color all the nodes in the subtree with
	  // the <bcc_timestamp> id
	  //
	  do {
	    popped = Pop_Node(ns);
#           if DEBUG_GAP_WALKER > 2
	    for (i = 0; i < level; i++) 
	      fprintf(GlobalData->gwlogfp," ");
	    fprintf(GlobalData->gwlogfp, "> popping %d - bcc %d\n",
		    popped,
		    bcc_time_stamp);
#           endif
	    assert(s->table[popped] != NULL);
	    s->table[popped]->bcc_id[s->table[popped]->colors] = bcc_time_stamp;
	    s->table[popped]->colors++;
	    //
	    // if this assertion should fail, increase the number
	    // of MAX_COLORS
	    //
	    assert(s->table[popped]->colors < MAX_COLORS);
	  } while (popped != tid);
	  //
	  // if ==, the assign a bcc to <this> node too
	  //
	  if (s->table[tid]->low_time == s->table[this]->dfs_time) {
#           if DEBUG_GAP_WALKER > 2
	    for (i = 0; i < level; i++) 
	      fprintf(GlobalData->gwlogfp," ");
	    fprintf(GlobalData->gwlogfp,
		    "> extra topping %d - bcc %d\n",
		    this,
		    bcc_time_stamp);
#           endif
	    s->table[this]->bcc_id[s->table[this]->colors] = bcc_time_stamp;
	    s->table[this]->colors++;
	    //
	    // if this assertion should fail, increase the number
	    // of MAX_COLORS
	    //
	    assert(s->table[popped]->colors < MAX_COLORS);
	  } 
	  bcc_time_stamp++;
	}
      } else
	//
	// if the node has not been completed (i.e., it is grey)
	//
	if (! s->table[tid]->done) {
#         if DEBUG_GAP_WALKER > 2
	  for (i = 0; i < level; i++) 
	    fprintf(GlobalData->gwlogfp," ");
	  fprintf(GlobalData->gwlogfp,
		  "> there is a back edge (dfs %d, fin %d)\n",
		  s->table[tid]->dfs_time,
		  s->table[tid]->fin_time);
#         endif
	  if (s->table[tid]->dfs_time < min_dfs)
	    min_dfs = s->table[tid]->dfs_time;
	}
    }
  }
  //
  // compute the low of this node
  //
  s->table[this]->low_time = min(min_dfs, lowest_time);

  //
  // if this is a root then output all the nodes on the stack
  // and assign them a bcc number
  //
  if (parent == -1) {
    while (ns->top) {
      popped = Pop_Node(ns);
#     if DEBUG_GAP_WALKER > 2
      for (i = 0; i < level; i++) 
	fprintf(GlobalData->gwlogfp," ");
      fprintf(GlobalData->gwlogfp,
	      "> megapopping %d - bcc %d\n",
	      popped,
	      bcc_time_stamp);
#     endif
      assert(s->table[popped] != NULL);
      s->table[popped]->bcc_id[s->table[popped]->colors] = bcc_time_stamp;
      s->table[popped]->colors++;
    }
    bcc_time_stamp++;
  }

  // we are done with this node
  // assign the finishing number
  //
  Set_Done_Bit(s, this);
  s->table[this]->fin_time = time_stamp;
  time_stamp++;
# if DEBUG_GAP_WALKER > 2
  for (i = 0; i < level; i++) 
    fprintf(GlobalData->gwlogfp," ");
  fprintf(GlobalData->gwlogfp,
	  "* done with the node %d (dfs %d, fin %d, low %d, stack %d, top %d)\n",
	  this,
	  s->table[this]->dfs_time,
	  s->table[this]->fin_time,
	  s->table[this]->low_time,
	  ns->top,
	  ns->nodes[ns->top - 1]);
# endif
}



void DFS(chunk_subgraph * s,
	 int (* filter)(CIEdgeT *)) {
  //
  // do a DFS of the graph s
  // avoid edges such that filter(e) == FALSE
  // 
  int
    i;
  nodes_stack
    * t = Create_Stack(s->size);

  Clear_All_Done_Bit(s);
  Clear_All_Visited_Bit(s);
  time_stamp = 0;
  for (i = 0; i < s->size; i++)
    if (! s->node[i].visited) {
      s->node[i].dfs_parent = -1;
      DFS_Visit(s, s->node[i].cid, -1, t, 0, filter);
    }
# if DEBUG_GAP_WALKER > 2
  fprintf(GlobalData->gwlogfp,
	  "%d nodes still on the stack: ",
	  t->top);
  for (i = 0; i < t->top; i++)
    fprintf(GlobalData->gwlogfp,
	    "%d ",
	    t->nodes[i]);
  fprintf(GlobalData->gwlogfp, "\n");
# endif
  Free_Stack(t);
}



bcc_array * Compute_Bcc_Array(chunk_subgraph * s) {
  //
  // return a bcc_array that can be indexed by the component
  // number and get all the chunks that re in that component
  //
  // note that one chunk could be in more than one component
  // (articulation point)
  //
  // (the size of the array is <bcc_time_stamp>)
  //
  int32
    i;
  int
    j,
    bcc,
    * cnt = (int *)safe_calloc(bcc_time_stamp, sizeof(int));
  bcc_array
    * b = (bcc_array *)safe_calloc(bcc_time_stamp, sizeof(bcc_array));
  ChunkInstanceT  
    * chunk;

  //
  // count the size of each component
  //
  for (i = 0; i < s->size; i++)
    for (j = 0; j < s->node[i].colors; j++) {
      bcc = s->node[i].bcc_id[j];
      b[bcc].size++;
    }

  //
  // allocate the chunk-id arrays for each component
  //
  for (bcc = 0; bcc < bcc_time_stamp; bcc++) {
    if (b[bcc].size)
      b[bcc].cid = (int *)safe_calloc(b[bcc].size, sizeof(int));
    cnt[bcc] = 0;
  }

  //
  // fill in the arrays, count the unique unitigs
  //
  for (i = 0; i < s->size; i++) 
    for (j = 0; j < s->node[i].colors; j++) {
      bcc = s->node[i].bcc_id[j];
      b[bcc].cid[cnt[bcc]] = s->node[i].cid;
      cnt[bcc]++;
      chunk = GetGraphNode(ScaffoldGraph->RezGraph, s->node[i].cid);
      assert(chunk != NULL);
      if (Is_Unique(chunk))
	b[bcc].uniques++;
    }

# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
	   "GapWalker has found %d biconnected components\n", 
	   bcc_time_stamp);
# endif

# if DEBUG_GAP_WALKER > 2
  for (bcc = 0; bcc < bcc_time_stamp; bcc++) {
    fprintf(GlobalData->gwlogfp,
	    "bcc %4d, size %4d - (",
	    bcc,
	    b[bcc].size); 
    for (j = 0; j < b[bcc].size; j++)
      fprintf(GlobalData->gwlogfp,
	      "%d ",
	      b[bcc].cid[j]);
    fprintf(GlobalData->gwlogfp,
	    ")\n");
  }
# endif

  return b;
}



static void Check_Distance(chunk_subgraph * s,
			   bcc_array * b,
			   int32 source,
			   int32 bcc,
                           float (* quality)(CIEdgeT *, int32)) {
  //
  // compute the shortest distance from all the non-uniques in the
  // bcc and print it
  //
  int
    cid,
    id;
  ChunkInstanceT
    * from = GetGraphNode(ScaffoldGraph->RezGraph, source),
    * to;
  //
  // build a subgraph on the <bcc> component
  //
  chunk_subgraph
    * bcc_s = Build_Subgraph_Bcc(b, bcc, Is_Not_Bogus);

  assert(from != NULL);
  assert(bcc_s != NULL);

# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
	  "-> Check_Distance on %d bcc %d\n",source,bcc);
# endif

  //
  // run shortest path from <source>
  // 
  Dijkstra(bcc_s, source, source, ALL_END, quality);

  //
  // scan each chunk in the bcc
  //
  for (id = 0; id < b[bcc].size; id++) {
    cid = b[bcc].cid[id];
    assert(bcc_s->table[cid] != NULL);
    to = GetGraphNode(ScaffoldGraph->RezGraph, cid);
#   if DEBUG_GAP_WALKER > 1
    fprintf(GlobalData->gwlogfp,
            "shortest distance from %d:[%d,%d] to %d:[%d,%d] is %5.2f, var %5.2f, d %5.2f, d_neck %5.2f\n",
            source,
            from->aEndCoord - from->aEndCoord,
            from->bEndCoord - from->aEndCoord,
            cid,
            to->aEndCoord - from->aEndCoord,
            to->bEndCoord - from->aEndCoord,
            bcc_s->table[cid]->distance.mean,
            bcc_s->table[cid]->distance.variance,
            bcc_s->table[cid]->d,
            bcc_s->table[cid]->d_neck);
#   endif
  }

  //
  // we don't need the subgraph anymore
  //
  Free_Subgraph(bcc_s);
}



//
// compute a tentative position for the non unique
// chunks for the component <bcc> (TO BE DONE)
//
// write the position in the A_end/B_end fields of the
// chunk_subgraph
//
static void Compute_Tentative_Position_Bcc(chunk_subgraph * s,
					   bcc_array * b,
					   int bcc) {
  //
  // compute a tentative position for the non unique
  // chunks for the component <bcc> (TO BE DONE)
  //
  // write the position in the A_end/B_end fields of the
  // chunk_subgraph
  //
  int
    id,
    cid;
  ChunkInstanceT
    * this;

  for (id = 0; id < b[bcc].size; id++) {
    //
    // get the chunk id number
    //
    cid = b[bcc].cid[id];
 
    //
    // get this node
    //
    this = GetGraphNode(ScaffoldGraph->RezGraph, cid);
    assert(this != NULL);

    //
    // if it is not unique, compute its position relative to the
    // unique that are present in this bcc, but how?
    //
    if (! Is_Unique(this)) {
      //
      // ?????
      //
      s->table[cid]->A_end.mean = 0.0;
      s->table[cid]->A_end.variance = 0.0;
      s->table[cid]->B_end.mean = 0.0;
      s->table[cid]->B_end.variance = 0.0;
    }
  }
}




int Scan_Components(chunk_subgraph * s,
		     bcc_array * b,
		     float (* quality)(CIEdgeT *, int32),
		     int (* filter)(CIEdgeT *)) {
  //
  // scan the components and check for internal consistency:
  // out conjecture is that internal edges should be ok, while
  // bad edges should occur w.h.p. as bridges
  //
  // the procedure scans the components in order of proportion of
  // unique contents: it will consider first components with high
  // percentage of uniques, because they should represent the most
  // conservative moves (but we do not want components completely
  // composed by unique, not work for us there).
  //
  // The function does not check edges such that filter(e) == FALSE
  // note: never tested with CHECK_COMPONENTS == 0
  //
  int32
    c_id,
    cid,
    tid,
    eid;
  int
    i,
    j,
    o,
    or,
    bcc,
    big_bcc = 0,
    tot_edges = 0,     // # edges checked in each component
    tot_edges_ = 0,    // grand total of edges checked in the bcc scan
    bcc_edges,         // edges inside the bcc
    bri_edges,         // bridge edges (connecting chunks with different "colors")
    bcc_ok_edges,      // count of "ok" bcc edges
    bcc_notok_edges,   // count of "nonok" bcc edges
    bri_ok_edges,      // you can figure out
    bri_notok_edges,   // as well
    max_colors = 0,
    tot_non_uniques = 0;
  CIEdgeT
    * e;
  ChunkInstanceT  
    * from,
    * to;
  double
    d,
    delta,
    delta_max = 0.0,     // max error for each component
    delta_max_max = 0.0, // max error for all the components
    delta_avg,           // avg error for each component
    delta_avg_ = 0.0,    // avg error for all components
    left,
    right;
  char
    id[STR_LEN];
  bcc_prop
    * bp;
  int
    inserted = 0;        // no of inserted chunks

  //
  // allocate the array that will store the
  // proportion of uniques in each component
  //
  bp = (bcc_prop *)safe_calloc(bcc_time_stamp, sizeof(bcc_prop));

  //
  // fill in the array
  //
  for (i = 0; i < bcc_time_stamp; i++) {
    bp[i].bcc = i;
    if (b[i].size > 0)
      bp[i].unique_prop = (double)(b[i].uniques) / (double)(b[i].size);
    else
      bp[i].unique_prop = 0.0;
  }

  //
  // sort w.r.t. the uniques proportion
  //
  qsort(bp, bcc_time_stamp, sizeof(bcc_prop), Compare_Bcc_Prop); 

  //
  // scan all the bcc in order of uniques prop
  //
  for (i = 0; i < bcc_time_stamp; i++) {
    bcc = bp[i].bcc;

    //
    // but filter out components too small
    // or composed entirely by uniques or
    // composed by few uniques
    //
    if ((b[bcc].size >= BCC_MIN_SIZE) &&
	(bp[i].unique_prop < 1.0) &&
	(bp[i].unique_prop > MIN_UNIQUES_PROP)) {

#     if CHECK_COMPONENTS > 0
      //
      // initialize counters
      //
      tot_edges = 0;
      bcc_edges = 0;
      bri_edges = 0;
      bcc_ok_edges = 0;
      bri_ok_edges = 0;
      bcc_notok_edges = 0;
      bri_notok_edges = 0;
      //
      // max/avg error in each component
      //
      delta_max = 0;
      delta_avg = 0;
#     endif

#     if DEBUG_GAP_WALKER > 1
      assert(GlobalData->gwlogfp != NULL);
      fprintf(GlobalData->gwlogfp,
	      "\n------------- bcc %d -------------\n",
	      bcc);
#     endif

      //
      // compute the tentative placement
      // of all the non-uniques in this component
      //
      // TO BE DONE
      //
      Compute_Tentative_Position_Bcc(s, b, bcc);

#     if CHECK_COMPONENTS > 0
      //
      // scan the chunk id
      //
      for (c_id = 0; c_id < b[bcc].size; c_id++) {
	//
	// get the id number
	//
	cid = b[bcc].cid[c_id];

	//
	// compute the max number of colors
	//
	assert(s->table[cid] != NULL);
	if (s->table[cid]->colors > max_colors)
	  max_colors = s->table[cid]->colors;

	//
	// get this node
	//
	from = GetGraphNode(ScaffoldGraph->RezGraph, cid);
	assert(from != NULL);
	
#       if DEBUG_GAP_WALKER > 2
        if (Is_Unique(from))
	  fprintf(GlobalData->gwlogfp,
		  "-*> unique chunk %d colors ",
		  cid);
	else
	  fprintf(GlobalData->gwlogfp,
		  "-*> chunk %d colors ",
		  cid);
	for (j = 0; j < s->table[cid]->colors; j++)
	  fprintf(GlobalData->gwlogfp, " %d",
		  s->table[cid]->bcc_id[j]);
	fprintf(GlobalData->gwlogfp, "\n");
#       endif

	//
	// scan the edges and check if they are ok
	//
	for (or = 0; or < NUM_ORIENTATIONS; or++) {
	  o = orient[or];
	  for (eid = 0; eid < s->table[cid]->num_edges[o]; eid++) {
	    e = s->table[cid]->edge[o][eid];
	    assert(e != NULL);

	    //
	    // check the edge
	    //
	    if (! filter(e))
	      continue;

	    //
	    // get the "other" node
	    //
	    if (cid == e->idA)
	      tid = e->idB;
	    else
	      tid = e->idA;
	    to = GetGraphNode(ScaffoldGraph->RezGraph, tid);
	    assert(to != NULL);
	    
	    //
	    // disregard repeat chunks
	    //
	    if ((from->flags.bits.cgbType != (unsigned int)UU_CGBTYPE)||
		(to->flags.bits.cgbType != (unsigned int)UU_CGBTYPE)) {
#             if DEBUG_GAP_WALKER > 1
	      fprintf(GlobalData->gwlogfp,
		      "> Found a repeat chunk %d(%d,%d) or %d(%d,%d)\n",
		      cid,
		      from->aEndCoord,
		      from->bEndCoord,
		      tid,
		      to->aEndCoord,
		      to->bEndCoord);
#             endif
	      continue;
	    }

	    //
	    // compute the actual distance between the
	    // chunk <from> and the chunk <to> from
	    // the simulator coordinates
	    //
	    switch (GetEdgeOrientationWRT(e,cid)) {
	    case AB_AB :
	      d = (double)(to->aEndCoord - from->bEndCoord);
	      break;
	    case AB_BA :
	      d = (double)(to->bEndCoord - from->bEndCoord);
	      break;
	    case BA_AB :
	      d = (double)(to->aEndCoord - from->aEndCoord);
	      break;
	    case BA_BA :
	      d = (double)(to->bEndCoord - from->aEndCoord);
	      break;
	    default :
	      assert(0);
	    }
	    
	    //
	    // adjust the sign
	    //
	    if (min(from->aEndCoord, from->bEndCoord) > min(to->aEndCoord, to->bEndCoord))
	      d = -d;

	    //
	    // compute the interval
	    //
	    assert(e->distance.variance >= 0.0);
	    left = (double)e->distance.mean - 3.0 * sqrt(e->distance.variance);
	    right = (double)e->distance.mean + 3.0 * sqrt(e->distance.variance);

	    //
	    // compute delta (error)
	    //
	    if (d < left) 
	      delta = left - d;
	    else if (d > right)
	      delta = d - right;
	    else 
	      delta = 0.0;

	    //
	    // distinguish internal bcc and bridges
	    //
	    for (j = 0; j < s->table[tid]->colors; j++)
	      if (s->table[tid]->bcc_id[j] == bcc)
		break;

	    if (j == s->table[tid]->colors) {
	      strcpy(id, "Bridge edge");
	      bri_edges++;
	    } else {
	      strcpy(id, "BCC edge");
	      bcc_edges++;

	      //
	      // update stats
	      //
	      delta_avg += delta;
	      if (delta > delta_max)
		delta_max = delta;

	      //
	      // count the edges
	      //
	      tot_edges++;
	      tot_edges_++;
	    }

#           if DEBUG_GAP_WALKER > 2
	    fprintf(GlobalData->gwlogfp,
		    "%s (%4d->%4d), correct %4.2f, edge length %5.2f, std %4.2f, contributions %2d, orientation %4s\n",
		    id,
		    cid,
		    tid,
		    d,
		    e->distance.mean,
		    sqrt(e->distance.variance),
		    e->edgesContributing,
	            Orientation_As_String(GetEdgeOrientationWRT(e,cid)));
#           endif

	    //
 	    // count and print the the bad bcc/bri edges
	    //
	    if (delta > DELTA_THRE) {
	      if (j == s->table[tid]->colors)
		bcc_notok_edges++;
	      else
		bri_notok_edges++;

#           if DEBUG_GAP_WALKER > 1
	      
	      fprintf(GlobalData->gwlogfp,
		      "* %s out of the interval [%f,%f], correct %f, error %f, from %d:[%d,%d], to %d:[%d,%d]\n",
		      id,
		      left,
		      right,
		      d,
		      delta,
		      cid,
		      from->aEndCoord,
		      from->bEndCoord,
		      tid,
		      to->aEndCoord,
		      to->bEndCoord);
	      PrintCIEdgeT(GlobalData->gwlogfp, ScaffoldGraph, "\t", e, e->idA);
#           endif
	    } else if (j == s->table[tid]->colors)
	      bcc_ok_edges++;
	    else
	      bri_ok_edges++;
	  }
	}
#       if 0
	if (from->type == UNRESOLVEDCHUNK_CGW)
	  Check_Distance(s, b, cid, bcc, quality);
#       endif
      }
#     endif

      tot_non_uniques += (b[bcc].size - b[bcc].uniques);
#     if DEBUG_GAP_WALKER > 1
      fprintf(GlobalData->gwlogfp,
	      "bcc %d, nodes unique/tot (%d/%d) %4.2f, edges bcc/bri/tot (%d/%d/%d)\n",
	      bcc,
	      b[bcc].uniques,
	      b[bcc].size,
	      bp[i].unique_prop,
	      bcc_edges,
	      bri_edges,
	      tot_edges);
      fprintf(GlobalData->gwlogfp,
	      "consistent bcc/bri edges (%d/%d), not consistent bcc/bri edges (%d/%d)\n",
	      bcc_ok_edges,
	      bri_ok_edges,
	      bcc_notok_edges,
	      bri_notok_edges);
      assert(b[bcc].size != 0);
      fprintf(GlobalData->gwlogfp,
	      "delta_max %f, delta_avg %f\n\n",
	      delta_max,
	      delta_avg/(double)b[bcc].size);
#     endif
      
#     if CREATE_DOT_FILE > 0
      //
      // if a component is really BAD print the cam file
      //
      if (delta_max > 10000) {
	FILE *
	  dot_file;
	char
	  filename[256];
	chunk_subgraph
	  * s = Build_Subgraph_Bcc(b, bcc, Is_Not_Bogus);

	fprintf(stderr,
		"* printing the .dot file of component %d in the dot directory\n", bcc);
	sprintf(filename, "./dot/%s.%d.dot", GW_Filename_Prefix, bcc);
	dot_file = file_open(filename,"w");
	assert(dot_file != NULL);
	Print_Dot_Subgraph(s, dot_file, filename, filter);
	fclose(dot_file);
      }
#     endif
 
#     if CHECK_COMPONENTS > 0
      //
      // update stats for the "all components"
      //
      if (delta_max > delta_max_max)
	delta_max_max = delta_max;
      delta_avg_ += delta_avg/(double)b[bcc].size;
      big_bcc++;
#     endif
    }
  }
# if DEBUG_GAP_WALKER > 2
  fprintf(GlobalData->gwlogfp,
	  "max_colors %d\n",
	  max_colors);
# endif

# if DEBUG_GAP_WALKER > 0
  fprintf(GlobalData->gwlogfp,
	  "Total non uniques in these bcc %d chunks, max error %f\n",
	  tot_non_uniques,
	  delta_max);
  fprintf(GlobalData->gwlogfp,
	  "Total of %d bcc-edges considered\n",
	  tot_edges);
  if (big_bcc)
    fprintf(GlobalData->gwlogfp,
	    "Summary: minimum size %d, #bcc/tot %d/%d, edges %d, max_err %5.2f, avg_err %5.2f\n",
	    BCC_MIN_SIZE,
	    big_bcc,
	    bcc_time_stamp,
	    tot_edges_,
	    delta_max_max,
	    delta_avg_/(double)big_bcc);
# endif

  fprintf(stderr,
	  "                      Selected: %7d chunks\n",
	  tot_non_uniques);

  //
  // release the memory of bp
  //
  free(bp);

  return inserted;
}



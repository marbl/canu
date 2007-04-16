
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
/*************************************************
* Module:  DataTypesREZ.c
* Description:
*   Declaration of common data types used in repeat resolution
* 
*    Programmer:  A. Delcher
*                 S. Lonardi (stelo@cs.purdue.edu)
*       Written:  17 May 99
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: DataTypesREZ.h,v 1.5 2007-04-16 15:35:41 brianwalenz Exp $
 * $Revision: 1.5 $
*/

#ifndef DATA_TYPES_REZ_H
#define DATA_TYPES_REZ_H

#include "AS_global.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "ChunkOverlap_CGW.h"


#define  CHECK_CELSIM_COORDS    0
    // If  1  use celsim coordinate info for checking


//  A  Gap_Chunk_t  contains the information for a chunk that can
//  be inserted in a particular gap in a particular scaffold

typedef  enum
  {ROCKS, STONES, WALKS}  Kind_Of_Fill_t;

typedef  struct
  {
   int  chunk_id;
   int  scaff_id;
   int  gap;
   int  cover_stat, link_ct;
   int32  len;
   LengthT  start, end;      // Relative to scaffold coordinates
#if  CHECK_CELSIM_COORDS
   int  sim_start, sim_end;  // Celsim coordinates relative to beginning of scaffold
#endif
   float  avg_edge_quality;  // Average quality of edge mates used to determine
                             // the start and  end positions of this chunk
   double  reference_variance;  // Variance to add back to  start  and  end
                               // variances to make relative to scaffold start
   int  index;               // Relative position of this node for insertion
                             // purposes.  Set by position on confirming path.
   char  copy_letter;        // Indicates which surrogate copy this is in case
                             //   of multistones.
   unsigned int  keep : 1;
       // Set true if chunk belongs in this gap
   unsigned int  flipped : 1;
       // Set true iff a_coord of this chunk > b_coord of this chunk
   unsigned int  split : 1;
       // Set true if surrogate chunk should be split off when this chunk
       // is inserted
   unsigned int  best : 1;
       // Set true if chunk has tightest variance and could be in only this
       // gap
   unsigned int  candidate : 1;
       // Set true if this chunk is a possibility to be chosen best rock
       // in the gap.
   unsigned int  path_confirmed : 1;
       // Set true iff this chunk's position was confirmed by walking
       // an overlap path
   unsigned int  visited : 1;
       // Set true if chunk has been visited in DFS
   unsigned int  finished : 1;
       // Set true if chunk's subtree has been finished in DFS
 }  Gap_Chunk_t;


//  A  Gap_Fill_t  describes a set of chunks that can
//  be inserted in a one particular gap in a particular scaffold

typedef  struct
  {
   int  gap;
   LengthT  start, end;     // Of gap, relative to scaffold coordinates
   LengthT  adjustment;     // Amount by which to modify all scaffold
                            //   chunks after this gap.  Set by confirming
                            //   path across the gap.
   double  ref_variance;    // Reference variance value for all chunks
                            //   place in this gap
   double  len;             // Difference between start and end means.  Negative
                            //   if this is an "overlap" gap
   int  num_chunks;
   int  left_cid;           // ids of chunks that border this gap
   int  right_cid;          //   left_cid = -1 for gap 0 
                            //   right_cid = -1 for last gap in scaffold
   unsigned int  has_path : 1;  // True indicates chunks in this gap were
                            //   confirmed by an overlap path
   Gap_Chunk_t  * chunk;    // array of potential chunks to go in this gap
  }  Gap_Fill_t;


//  A  Scaffold_Fill_t  describes multiple sets of chunks that can
//  be inserted in all the gaps in a particular scaffold

typedef  struct
  {
   int  scaff_id;
   int  num_gaps;
   int  keep_ct;           // number of items with keep flag on
   int  added_to_id;       // scaff id to which this scaff is added
   int  added_to_ct;       // number of scaffs to which this scaff is added
   Gap_Fill_t  * gap;
  }  Scaffold_Fill_t;


//  Information about a  node  in the rez graph that is used
//  when looking for a path that passes through it, e.g., a "stone".

typedef  struct
  {
   int32  id;
   double  lo, hi;       // range of scaff coords of *high* end of target
   ChunkOrientationType
       orient;           // orientation of link that connects the start
                         //   of the search (i.e., "from") to this target
   double  where;        // scaff coord of *high* end of target on overlap
                         //   path if it was found
   double  total;        // sum of lengths of all fragments on overlap
                         //   path that determines  where .  Used to
                         //   estimate variance of  where .
   int  found;           // set true if found within above range
   int  next;            // subscript of next target on max path
  }  Target_Info_t;


typedef  enum
  {OVERLAP, LINK}  Kind_Of_Stone_Edge;

//  An edge between two stones in a gap.

typedef  struct
  {
   int32  from, to;
   int32  progress;      // number of bases "to" node extends past end of
                         //   "from" node on "right" side
   int32  a_hang;        // same on "left" side
   int32  length;        // of overlap region
   double  quality;      // number of errors (diffs) / length of overlap
   int32  next;
   Kind_Of_Stone_Edge  kind;
   unsigned int  bad : 1;
  }  Stone_Edge_t;


//  Information about best path in DAG for stones

typedef  struct
  {
   int  path_len;      // max number of predecessors on path
   int  from;          // subscript of predecessor on max path
   int  hi_position;   // dna bases of high position of this node
   int  lo_position;   // dna bases of low position of this node
                       //   determined from a_hangs
   int  a_hang;        // relative number of bases from low position
                       //   of previous node in path to low position of
                       //   this node
   int  total_olap;    // total bases used (including overlaps) to
                       //   get to  position;
   unsigned int  hit : 1;
   unsigned int  done : 1;
                       // set so we can eliminate cycles in olap graph
  }  Path_Info_t;


//  Information about overlapping chain of stones

typedef  struct
  {
   int  sub;           // subscript of end node in chain
   int  ct;            // number of nodes in chain
   unsigned int  has_start : 1;
                       // set true iff this component includes the
                       // scaffold contig at start of gap
   unsigned int  has_target : 1;
                       // set true iff this component includes the
                       // scaffold contig at end of gap
  }  Component_Info_t;

typedef enum { ASS_OVLP_TRUE, ASS_OVLP_FALSE, NO_ASS_OVLP } OverlapStatusREZ;

//
// the chunk_subgraph_node
//

#define NUM_ORIENTATIONS 4
#define MAX_COLORS 10

typedef struct {
  int
    cid,                         // chunk id
    num_edges[NUM_ORIENTATIONS], // # outgoing edges for AB_AB, AB_BA, BA_AB, and BA_BA
    dfs_time,                    // depth first search timestamp
    fin_time,                    // finishing timestamp
    low_time,                    // low timestamp (for biconnected components, see CLR page 496)
    dfs_parent,                  // parent in the dfs tree
    path_id,                     // if > 0 indicates the id_number of a path
    path_parent,                 // the parent id on the path
    union_id,                    // used for the union (each component will get a unique union_id, hopefully)
    colors,                      // number of colors (relative to bcc labeling)
    end,                         // see AS_CGW/dataTypes.h (can be NO_END, A_END, B_END, ALL_END) (used in Dijkstra)
    bcc_id[MAX_COLORS];          // biconnected component ID (see DFS_Visit())
  unsigned int
    already_used : 1,            // TRUE if used in multiple paths
    visited : 1,                 // TRUE if visited (correspond to GREY of CLR)
    done : 1,                    // TRUE if is done (correspond to BLACK of CLR)
    path_bit : 1;                // TRUE if this node is on some path
  float
    d,                           // shortest path distance from the source (sum of the weight on the edges)
    d_neck;                      // bootleneck distance (stores the "worst" edge on the path 'till here)
  LengthT
    A_end,                       // tentative position of the A end (relative to the scaffold coordinates)
    B_end,                       // tentative position of the B end (relative to the scaffold coordinates)
    distance;                    // distance from the source in the DFS visit
  CIEdgeT
    * * edge[NUM_ORIENTATIONS],  // a pointer (for all the possibile orientations)
                                 // to an array of pointers to the outgoing edges
    * best_edge;                 // the edge chosen by shortest path 
} chunk_subgraph_node;

/*
   the data structure for the subgraph is an array
   of <chunk_subgraph_nodes> and a table of entries
  
     2             123
   +---+---+------+---+---+
   ! * ! N ! .... ! * ! * !  chunk_table (some pointers could be NIL)
   +-!-+---+------+-!-+---+
      \__            \_
         v             v
   +---+---+------+---+---+
   ! 5 ! 2 ! .... !234!123! cid        } 
   +---+---+------+---+---+             >  chunk_subgraph_node(s)
   !   !   ! .... !   !   ! other info }
   +---+---+------+---+---+
  
   the following relation MUST be maintained for all 0 <= i < s->size
   (s is of type * chunk_subgraph)
  
             s->table[s->node[i].cid] == &(s->node[i])
  
   The edges are stored in four arrays, one for each
   orientation (AB_AB, AB_BA, BA_AB, BA_BA). The correspondence
   between orientation and the index in [0..3] is established by the
   function or2num()
*/

typedef struct {
  int
    max,               // entries in the chunk_table
    size;              // entries in the subgraph
  chunk_subgraph_node
    * * table; 
  chunk_subgraph_node
    * node;
} chunk_subgraph;

//
// the data structure to store the
// quality of each edge (see GapWalkerREZ.c)
//

typedef struct {
  CIEdgeT
    * edge;
  int
    orientation;     // orientation in 0..3 code (see or2num() )
  float
    quality;         // quality value
} edge_quality;

//
// the stack used in the DFS
//

typedef struct {
  int
    max_size, // the max # of elements
    top,      // the index of the top (# of components)
              // it is always the first empty position in the stack
    * nodes;  // the base pointer
} nodes_stack;

//
// biconnected components array
//

typedef struct {
  int
    size,     // the size of the array cid[]
    * cid,    // a sequence of chunk id for this component
    uniques;  // the number of unique chunks in this component
} bcc_array;

typedef struct {
  int
    bcc;         // bcc id
  float
    unique_prop; // unique proportion
} bcc_prop;

//
// heaps
//
typedef chunk_subgraph_node * chunk_subgraph_node_ptr;
typedef CIEdgeT * CIEdgeT_ptr;

#endif


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
#ifndef _AS_CGB_BUBBLE_GRAPH_H_
#define _AS_CGB_BUBBLE_GRAPH_H_

/*
 * FLAG VALUES
 */

#define AS_CGB_BUBBLE_E_ALL             0xFF /* "Any edge" flag for iterator */
#define AS_CGB_BUBBLE_E_DOVETAIL        0x01 /* Dovetail edge */
#define AS_CGB_BUBBLE_E_CONTAIN         0x02 /* Contain edge */
#define AS_CGB_BUBBLE_E_VALID           0x04 /* Tested by DFS and valid */
#define AS_CGB_BUBBLE_E_UNUSED          0x08 /* Not processed by DFS */
#define AS_CGB_BUBBLE_E_IN_BUBBLE       0x10 /* Used by Bubble Popper to mark
					        edges in current bubble */

#define AS_CGB_BUBBLE_V_REVERSED        0x01 /* Reverse ori. in component */
#define AS_CGB_BUBBLE_V_NEW             0x02 /* Not processed by DFS. */
#define AS_CGB_BUBBLE_V_STACKED         0x04 /* On DFS stack */
#define AS_CGB_BUBBLE_V_DONE            0x08 /* Processed by DFS */
#define AS_CGB_BUBBLE_V_VALID           0x08 /* Intentionally same as DONE. */
#define AS_CGB_BUBBLE_V_CONTAINED       0x10 /* Contained fragment */
#define AS_CGB_BUBBLE_V_INADMISSABLE    0x20 /* Can't be in bubble */
#define AS_CGB_BUBBLE_V_IN_BUBBLE       0x40 /* Used by Bubble Popper to mark
					        vertices in current bubble */


/*
 *
 * BUBGRAPH 
 *
 */

typedef struct BubGraph {
  Tedge *e;
  Tfragment *v;
  uint16 *vFlags;
  uint16 *eFlags;
  int64 *vPos;
} BubGraph;

typedef BubGraph * BubGraph_t;

/*
 * METHODS
 */

/* Initialize and destroy a bubble graph. */
void BG_initialize(BubGraph_t bg, Tfragment *frags, Tedge *edges);
void BG_destroy(BubGraph_t bg);

/* Get the vertex and edge sets.  To avoid writing lots of wrapper functions,
   most operations on the graph are performed by getting the vertex/edge set,
   then using one of Clark's operators on the set. */

Tfragment *BG_vertices(BubGraph_t bg);
Tedge *BG_edges(BubGraph_t bg);

/* Given an edge and an incident vertex, returns the ID of the other incident
   vertex. */
IntFragment_ID 
BG_getOppositeVertex(BubGraph_t bg, IntEdge_ID e, IntFragment_ID v); 

/* Returns TRUE if the vertex has forward orientation, or FALSE otherwise. */
int BG_vertexForward(BubGraph_t bg, IntFragment_ID v);

/* Returns the total degree of the given vertex, counting only those edges
   which have matching flags. */
int BG_degree(BubGraph_t bg, IntFragment_ID v, uint16 flags);

/* Returns the indegree of the given vertex, counting only those edges 
   which have the given flags set. */
int BG_inDegree(BubGraph_t bg, IntFragment_ID v, uint16 flags);

/* Returns the outdegree of the given vertex, counting only those edges 
   which have the given flags set. */
int BG_outDegree(BubGraph_t bg, IntFragment_ID v, uint16 flags);

/* Methods to manipulate the vertex and edge flags. */

/* Turns on a particular flag */
void BG_V_setFlag(BubGraph_t bg, IntFragment_ID v, uint16 flag);
void BG_E_setFlag(BubGraph_t bg, IntEdge_ID e, uint16 flag);
void BG_E_setFlagSymmetric(BubGraph_t bg, IntEdge_ID e, uint16 flag);

/* Turns off a particular flag. */
void BG_V_clearFlag(BubGraph_t bg, IntFragment_ID v, uint16 flag);
void BG_E_clearFlag(BubGraph_t bg, IntEdge_ID e, uint16 flag);
void BG_E_clearFlagSymmetric(BubGraph_t bg, IntEdge_ID e, uint16 flag);

/* Tests a particular flag. */
uint16 BG_V_isSetFlag(BubGraph_t bg, IntFragment_ID v, uint16 flag);
uint16 BG_E_isSetFlag(BubGraph_t bg, IntEdge_ID e, uint16 flag);

/* Returns the full flag vector for an edge/vertex. */
uint16 BG_V_getFlags(BubGraph_t bg, IntFragment_ID v);
uint16 BG_E_getFlags(BubGraph_t bg, IntEdge_ID e);

/* Sets and gets the position field of each vertex, which gives the
   position of the left endpoint relative to some path. */
void BG_V_setDistance(BubGraph_t bg, IntFragment_ID v, int64 p);
int64 BG_V_getDistance(BubGraph_t bg, IntFragment_ID v);

/*
 *
 * EDGE ITERATOR
 *
 */

typedef enum BG_E_IterType { bgeiIn, bgeiOut, bgeiBoth } BG_E_IterType;

typedef struct BG_E_Iter {
  BG_E_IterType type;
  IntEdge_ID cur;
  IntEdge_ID end;
  IntFragment_ID v;
} BG_E_Iter;

typedef BG_E_Iter * BG_E_Iter_t;

/*
 * METHODS
 */

/* Creates a new iterator and returns the first edge with respect to
   the given vertex matching both the type (either forward, back, or
   either) and the given flags.  If no such edge exists, a random value is
   returned and BGEI_end() returns TRUE. */
IntEdge_ID BGEI_bgn(BubGraph_t bg, BG_E_Iter_t it, IntFragment_ID v, 
		    BG_E_IterType t, uint16 flags);

/* Returns the ID of the next edge matching the flags and type.  Return
   value is undefined if there is no next edge. */
IntEdge_ID BGEI_next(BubGraph_t bg, BG_E_Iter_t it, uint16 flags);

/* Returns TRUE if the end of the edges satisfying the iterator have been
   reached, or false otherwise. */
int BGEI_end(BG_E_Iter_t it);

/* Returns the item currently pointed to by the iterator, or a random value if
   BGEI_end(). */
IntEdge_ID BGEI_cur(BG_E_Iter_t it);


#endif

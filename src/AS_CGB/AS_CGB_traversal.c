
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
static char CM_ID[] 
= "$Id: AS_CGB_traversal.c,v 1.1.1.1 2004-04-14 13:49:42 catmandew Exp $";
/*********************************************************************
 *
 * Module: AS_CGB_traversal.c
 *
 * Description: Graph traversal for the Chunk Graph Builder.
 * This module takes the current marked graph and returns
 * a best ranking of the verticies in fragment_rank[].
 *
 * Assumptions:
 *
 * 1. Overview
 * 
 * This functional unit transverses the current state of the marked graph
 * and returns a best ranking of the verticies.
 * 
 * 2. Memory Usage
 * 
 * There is one statically allocated first-in first-out queue of length
 * MAXQUEUELEN.  There is no dynamicly allocated memory.
 * 
 * 3. Interface
 * 
 * 4. Design
 * 
 * The graph is transversed in a breadth first manner, except over chunks
 * which are visisted in a depth first manner.
 * 
 * 5. Limitations
 * 
 * Currently it is a assumed that a chunk has two ends.  Thus a circular
 * chunk is not visited.
 * 
 * 6. Status
 * 
 * Author: Clark Mobarry
 ***********************************************************************/

/*********************************************************************/
/* System include files */
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>

/*************************************************************************/
/* Local include files */
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_CGB_all.h"
#include "AS_CGB_traversal.h"

/*************************************************************************/
/* Constant definitions; Macro definitions; type definitions */
#define MAXQUEUELEN   ((1 << 24)-1) // 0xffffff
//#define MAX_ADJACENCY_FOLLOWED 5
#undef MAX_ADJACENCY_FOLLOWED
#define NOT_RANKED           (CDS_UINT32_MAX)

/*************************************************************************/
/* Static Globals */

/*************************************************************************/
/* Function prototypes for internal static functions */

/*************************************************************************/
/* Basic Queue Operations */

static void as_addq
(IntFragment_ID item,
 IntFragment_ID thequeue[],
 int lenqueue, 
 IntFragment_ID *front, 
 IntFragment_ID *rear)
{
  /* Insert item into the circular queue stored in thequeue[0..lenqueue-1].
     rear points to one position counterclockwise from the last item and 
     front is the first item in thequeue.
  */
  thequeue[*rear] = item; /* insert the new item */
  *rear = (*rear + 1) & lenqueue; /* advance the rear clockwise */
  if( *front == *rear) { /* The queue is full */
    fprintf(stderr,"lenqueue=%d front=" F_IID " rear=" F_IID "\n",
            lenqueue, *front, *rear);
    assert( *front != *rear ); /* Increase MAXQUEUELEN ? */
    /* CALL queueFULL */
  }
}

static void as_delq
( IntFragment_ID *item,
  IntFragment_ID thequeue[],
  int lenqueue, 
  IntFragment_ID *front,
  IntFragment_ID *rear)
{
  /* Removes the rear element of the queue, and stores it
     in the item. This is a undo or stack pop operation. */
  if( *front == *rear) { /* The queue is empty. */
    assert( *front != *rear );
    /* CALL QUEUEempty */
  }
  *item = thequeue[*rear];
  *rear = (*rear - 1 + lenqueue) & lenqueue;
}

static void as_subq
( IntFragment_ID *item,
  IntFragment_ID thequeue[], 
  int lenqueue, 
  IntFragment_ID *front,
  IntFragment_ID *rear)
{
  /* Removes the front element of the queue, and stores it
     in the item. */
  if( *front == *rear) { /* The queue is empty. */
    assert( *front != *rear );
    /* CALL QUEUEempty */
  }
  *item = thequeue[*front];
  *front = (*front + 1) & lenqueue;
}

static void as_diagq
( FILE *fout,
  IntFragment_ID thequeue[],
  int lenqueue,
  IntFragment_ID front,
  IntFragment_ID rear)
{ 
  IntFragment_ID iq;
  fprintf(fout,"the queue front=" F_IID " rear=" F_IID "\n",front,rear);
  if( front != rear ) {
  if( front < rear ) {
    for(iq=front;iq<rear;iq++) {
      fprintf(fout,F_IID " ",thequeue[iq]);
    } 
  } else {
    for(iq=front;iq<lenqueue;iq++) {
      fprintf(fout,F_IID " ",thequeue[iq]);
    }
    for(iq=0;iq<rear;iq++) {
      fprintf(fout,F_IID " ",thequeue[iq]);
    }
  }}
  fprintf(fout,"\n");
}

/**************************************************************/

static void as_dfs_intrachunk
(
 IntFragment_ID iv0, 
 IntRank   *nfound,
 IntFragment_ID nfrag, 
 Tfragment  frags[],
 IntEdge_ID nedge, 
 Tedge      edges[], 
 IntRank    fragment_rank[])
{
  int is0;
  IntFragment_ID iv1;
  IntEdge_ID ie1,ie0;
  int ne0;

  fragment_rank[iv0] = *nfound;
  *nfound = *nfound + 1;
  for(is0=0;is0<2;is0++){
    /* Search all vertices "iv1" adjacent from "iv0" */
    ie0 = get_segstart_vertex(frags,iv0,is0);
    ne0 = get_seglen_vertex(frags,iv0,is0);
    for(ie1=ie0; ie1 < ie0+ne0; ie1++) {
      iv1 = get_bvx_edge(edges,ie1);
      if( (get_nes_edge(edges,ie1) == AS_CGB_INTRACHUNK_EDGE) && 
	  (fragment_rank[iv1] == NOT_RANKED) )
	{ /* iv1 is unexplored */
	  as_dfs_intrachunk(iv1, nfound,
			    nfrag, frags, nedge, edges, 
			    fragment_rank);
	}
    }
  }
}

static void as_dfs
(IntFragment_ID iv0, 
 IntRank   *nfound,
 IntFragment_ID nfrag, 
 Tfragment  frags[],
 IntEdge_ID nedge, 
 Tedge      edges[], 
 IntRank    fragment_rank[])
{
  int is0;
  IntFragment_ID iv1;
  IntEdge_ID ie1,ie0;
  int ne0;
  
  fragment_rank[iv0] = *nfound;
  *nfound = *nfound + 1;
  for(is0=0;is0<2;is0++) {
    /* Search all vertices "iv1" adjacent from "iv0" */
    ie0 = get_segstart_vertex(frags,iv0,is0);
    ne0 = get_seglen_vertex(frags,iv0,is0);
    for(ie1=ie0; ie1 < ie0+ne0; ie1++) {
      iv1 = get_bvx_edge(edges,ie1);
      if(
         (get_nes_edge(edges,ie1) != AS_CGB_REMOVED_BY_DUPLICATE_CON) 
	 && /* CMM be careful 
	       about the ordering of 
	       the nes_types. */
	 (fragment_rank[iv1] == NOT_RANKED) ) 
	{ /* iv1 is unexplored */
	  as_dfs(iv1, nfound,
		 nfrag, frags, nedge, edges, 
		 fragment_rank);
	}
    }
  }
}

static IntFragment_ID thequeue[MAXQUEUELEN];

static void as_bfs
(IntFragment_ID ivi, 
 IntRank   *nfound,
 IntFragment_ID nfrag, 
 Tfragment  frags[],
 IntEdge_ID nedge, 
 Tedge      edges[], 
 IntRank    fragment_rank[])
{ /* Breadth First Search of the graph. 
     A breadth first search of the graph is carried out
     beginning at fragment, ivi.
     All vertices that have been visited are marked as
     fragment_rank[.] != NOT_RANKED.
     The graph and the array fragment_rank[] are initialized.
  */

  IntFragment_ID iv0,iv1;
  IntFragment_ID front,rear;
  IntEdge_ID ie1,ie0;
  int ne0;
#ifdef DEBUG8
  FILE *fout;
  fout = stdout;
  { 
    IntFragment_ID iv2;
    fprintf(fout,"as_bfs: nfound=" F_SIZE_T "\n",*nfound);
    fprintf(fout,"as_bfs: iv2,fragment_rank[iv2]\n");
    for(iv2=0;iv2<nfrag;iv2++) {
      fprintf(fout,"%d, %d\n",iv2,fragment_rank[iv2]);
    }
  }
#endif

  iv0 = ivi;
  /* Initialize the queue to be an empty queue. 
     This queue is of unexplored vertices. */
  front = rear = 0;

  // A fragment that has been ranked and is in the queue is
  // BFS_GRAY.  fragment_color[iv0] = BFS_GRAY;

  fragment_rank[iv0] = (*nfound)++; /* Since the vertices indexes are zero based. */
  as_addq(iv0,thequeue,MAXQUEUELEN,&front,&rear);

#ifdef DEBUG8
  { 
    IntFragment_ID iv2;
    fprintf(fout,"as_bfs: nfound=" F_SIZE_T "\n",*nfound);
    fprintf(fout,"as_bfs: iv2,fragment_rank[iv2]\n");
    for(iv2=0;iv2<nfrag;iv2++) {
      fprintf(fout,"%d, %d\n",iv2,fragment_rank[iv2]);
    }
  }
  as_diagq(fout,thequeue,MAXQUEUELEN,front,rear);
#endif

  for(; front != rear /* The queue is not empty. */; ) {
    int is0;
    as_subq(&iv0,thequeue,MAXQUEUELEN,&front,&rear);

    /* Search all vertices "iv1" adjacent from "iv0" */
    for(is0=0;is0<2;is0++){
      ie0 = get_segstart_vertex(frags,iv0,is0);
      ne0 = get_seglen_vertex(frags,iv0,is0);
#ifdef DEBUG8
      fprintf(fout,"as_bfs: fragment %d\n",iv0,ie0,ne0);
#endif
#ifdef MAX_ADJACENCY_FOLLOWED
      ne0 = min(ne0,MAX_ADJACENCY_FOLLOWED);
#endif      
      /* Perform a Breadth First Search on the remaining edges. */
      for(ie1=ie0; ie1 < ie0+ne0; ie1++) {
	iv1 = get_bvx_edge(edges,ie1);
	if( 
           /* CMM be careful 
              about the ordering of 
              the nes_types. */
	   (fragment_rank[iv1] == NOT_RANKED) ) { /* iv1 is unexplored */
          // A fragment that has been ranked and is in the queue is
          // BFS_GRAY.  fragment_color[iv1] = BFS_GRAY;
	  fragment_rank[iv1] = (*nfound)++;
	  as_addq(iv1,thequeue,MAXQUEUELEN,&front,&rear);
	}
#ifdef DEBUG8
	{ 
	  IntFragment_ID iv2;
	  fprintf(fout,"as_bfs: nfound=" F_SIZE_T "\n",*nfound);
	  fprintf(fout,"as_bfs: iv2,fragment_rank[iv2]\n");
	  for(iv2=0;iv2<nfrag;iv2++) {
	    fprintf(fout,"%d, %d\n",iv2,fragment_rank[iv2]);
	  }
	}
#endif
      }
    }
    // A fragment that has been ranked and is no longer in the queue
    // is BFS_BLACK. fragment_color[iv0] = BFS_BLACK;
#ifdef DEBUG8
    as_diagq(fout,thequeue,MAXQUEUELEN,front,rear);
#endif    
  }
}

void as_graph_traversal
(
 FILE *    fout,
 Tfragment frags[],
 Tedge     edges[], 
 IntRank   fragment_rank[]
 )
{ /* A marked traversal of the graph. */
  IntRank nconnected=0,nfound=0,ncontained=0,ncircl=0;
  IntFragment_ID ivi;
  
  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);

  for(ivi=0; ivi<nfrag; ivi++) {
    fragment_rank[ivi] = NOT_RANKED;
  }
  for(ivi=0;ivi<nfrag;ivi++) {
    if(fragment_rank[ivi] == NOT_RANKED){
	nconnected++;
	as_bfs
	  (ivi, &nfound,
	   nfrag, frags,
	   nedge, edges,
	   fragment_rank
	   );
    }
  }

  for(ivi=0;ivi<nfrag;ivi++) {
    assert(fragment_rank[ivi] != NOT_RANKED);
    //ncircl++;
    //nfound++;
  }
  fprintf(fout,"FRAGMENT OVERLAP GRAPH TRAVERSAL STATISTICS\n\n"
	  "%15" F_IIDP " : number of fragments in overlap graph\n"
	  "%15" F_SIZE_TP " : number of fragments found by breadth first search\n"
	  "%15" F_SIZE_TP " : number of weakly connected graph components\n"
	  "%15" F_SIZE_TP " : number of contained fragments\n"
	  "%15" F_SIZE_TP " : number of circular chunks\n\n\n",
	  nfrag,nfound,nconnected,ncontained,ncircl);
  assert( nfound == nfrag );

  /* 
     Some of the search_method()s assume that circular chunks do not occur.
     If this assert ever fails than look for fragments with
     get_lab_fragment() == AS_CGB_INTRACHUNK_FRAG, but have not been visited.
  */
}

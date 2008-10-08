
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

 Module: UtilsREZ.c

 Description: Contains small utility functions

 Programmer: K. Reinert
             S. Lonardi (stelo@cs.purdue.edu)

 Assumptions: none

**********************************************************************/

static char *rcsid = "$Id: UtilsREZ.c,v 1.10 2008-10-08 22:03:00 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include "UtilsREZ.h"

// -------------------------------------
// bit manipulation (for chunk subgraph)
// -------------------------------------

void Clear_All_Path_Bit(chunk_subgraph * s) {
  int i;
  assert(s != NULL);
  for (i = 0; i < s->size; i++)
    s->node[i].path_bit = FALSE;
}



void Clear_Path_Bit(chunk_subgraph * s, int32 cid) {
  assert(s != NULL);
  assert(cid < s->max);
  assert(s->table[cid] != NULL);
  s->table[cid]->path_bit = FALSE;
}



void Set_Path_Bit(chunk_subgraph * s, int32 cid) {
  assert(s != NULL);
  assert(cid < s->max);
  assert(s->table[cid] != NULL);
  s->table[cid]->path_bit = TRUE;
}



void Clear_All_Visited_Bit(chunk_subgraph * s) {
  int i;
  assert(s != NULL);
  for (i = 0; i < s->size; i++)
    s->node[i].visited = FALSE;
}



void Clear_Visited_Bit(chunk_subgraph * s, int32 cid) {
  assert(s != NULL);
  assert(cid < s->max);
  assert(s->table[cid] != NULL);
  s->table[cid]->visited = FALSE;
}



void Set_Visited_Bit(chunk_subgraph * s, int32 cid) {
  assert(s != NULL);
  assert(cid < s->max);
  assert(s->table[cid] != NULL);
  s->table[cid]->visited = TRUE;
}



void Clear_All_Done_Bit(chunk_subgraph * s) {
  int i;
  assert(s != NULL);
  for (i = 0; i < s->size; i++)
    s->node[i].done = FALSE;
}



void Set_Done_Bit(chunk_subgraph * s, int32 cid) {
  assert(s != NULL);
  assert(cid < s->max);
  assert(s->table[cid] != NULL);
  s->table[cid]->done = TRUE;
}


// ---------------------
// <nodes_stack> methods
// ---------------------


void Push_Node(nodes_stack * s, int v) {
  assert(s != NULL);
  assert(s->top < s->max_size);
  s->nodes[s->top] = v;
  (s->top)++;
}



int Pop_Node(nodes_stack * s) {
  assert(s != NULL);
  assert(s->top);
  return s->nodes[--(s->top)];
}



// get the value of Top without popping
int Top(nodes_stack * s) {
  assert(s != NULL);
  assert(s->top);
  return s->nodes[(s->top - 1)];
}



// crate a stack of <no_element> size
nodes_stack * Create_Stack(int no_elements) {
  nodes_stack * s = (nodes_stack *)safe_calloc(1, sizeof(nodes_stack));
  assert(no_elements);
  s->top = 0;
  s->max_size = no_elements;
  s->nodes = (int *)safe_calloc(no_elements, sizeof(int));
  return s;
}



void Free_Stack(nodes_stack * s) {
  assert(s != NULL);
  safe_free(s->nodes);
  safe_free(s);
}

// ------------
// CIEdge stuff
// ------------

char * Orientation_As_String (ChunkOrientationType orient) {
  //
  //  Return string equivalent of orient
  //
  switch  (orient) {
  case  AB_AB :
    return  "AB_AB";
  case  AB_BA :
    return  "AB_BA";
  case  BA_BA :
    return  "BA_BA";
  case  BA_AB :
    return  "BA_AB";
  default :
    return  "*???*";
  }
}



int or2num(ChunkOrientationType o) {
  //
  // convert a <ChunkOrientationType> to a number
  // between 0 and 3
  //
  switch (o) {
  case AB_AB :
    return OR2NUM_AB_AB;
  case AB_BA :
    return OR2NUM_AB_BA;
  case BA_AB :
    return OR2NUM_BA_AB;
  case BA_BA :
    return OR2NUM_BA_BA;
  default :
    assert(0);
  }
  return -1;
}

/*--------------------------------------------------------------------*/
/* Interval Math */
/*--------------------------------------------------------------------*/


int Intersection (LengthT * a, LengthT * b) {
  //  Return the number of bases by which the closed interval  [mean - 3stdDev, mean + 3stdDev]
  //  intersects the closed interval  [mean - 3stdDev, mean + 3stdDev]
  return Interval_Intersection (a->mean - 3.0 * sqrt(a->variance), a->mean + 3.0 * sqrt(a->variance),
				b->mean - 3.0 * sqrt(b->variance), b->mean + 3.0 * sqrt(b->variance));
}



int Interval_Intersection (int a, int b, int c, int d) {
  //  Return the number of bases by which the closed interval [a, b]
  //  intersects the closed interval [c, d]
  if  (d < a || b < c)
    return  0;
  else
    return  1 + MIN (b, d) - MAX (a, c);
}


FILE*  file_open(const char* fileName, const char* mode)
{
   FILE* fp;
   errno = 0;
   fp = fopen (fileName,mode);
   if (errno) {
     fprintf(stderr, "Failed to open '%s' for '%s': %s\n", fileName, mode, strerror(errno));
     exit(1);
   }
   return  fp;
}

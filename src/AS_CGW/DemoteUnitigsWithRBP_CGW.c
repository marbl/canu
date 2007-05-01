
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "ScaffoldGraphIterator_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "GraphCGW_T.h"


void DemoteUnitigsWithRBP(FILE *stream, GraphCGW_T *graph){
  GraphNodeIterator nodes;
  GraphEdgeIterator edges;
  NodeCGW_T *node = NULL, *otherNode = NULL;
  CDS_CID_t otherNodeId;

  assert(graph->type == CI_GRAPH);
  
  
  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&nodes))){
    EdgeCGW_T *edge = NULL;
    int numAEndConfirmOverlap = 0, numBEndConfirmOverlap = 0;
    
    InitGraphEdgeIterator(graph, node->id, ALL_END, ALL_EDGES,
                          GRAPH_EDGE_DEFAULT , &edges);
    while(NULL != (edge = NextGraphEdgeIterator(&edges))){
      /* Check for correct end and confirmed statusCount for each end */
      if (edge->idA == node->id){
	otherNodeId = edge->idB;
      } else {
	assert(edge->idB == node->id);
	otherNodeId = edge->idA;
      }
      otherNode = GetNodeCGW_T(graph->nodes, otherNodeId);
      if((edge->edgesContributing >= 2) && isOverlapEdge(edge) && (otherNode->bpLength.mean > 1500.0)){
	switch(GetEdgeOrientationWRT(edge, node->id)){
	  /* EdgeMate from the A-End */
	case BA_BA:
	case BA_AB:
	  numAEndConfirmOverlap++;
	  break;
	case AB_BA:
	case AB_AB:
	  numBEndConfirmOverlap++;
	  break;
	default:
	  assert(0);
	}
      }
    }
    /* fprintf(stream, "Unitig %d: branchA %d branchB %d\n", node->id, numAEndConfirmOverlap,
	    numBEndConfirmOverlap); */
    if((numAEndConfirmOverlap > 1) && (numBEndConfirmOverlap > 1)){
      /* fprintf(stream, "Demoting unitig\n"); */
      node->flags.bits.isUnique = FALSE;
      node->type = UNRESOLVEDCHUNK_CGW;
    }
  }
}

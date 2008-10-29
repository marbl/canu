
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
static char *rcsid = "$Id: SplitScaffolds_CGW.c,v 1.13 2008-10-29 10:42:46 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "ChiSquareTest_CGW.h"
#include "ChunkOverlap_CGW.h"


void SplitScaffolds(ScaffoldGraphT *graph){

  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;
  CDS_CID_t firstSplitScaffold = NULLINDEX;

  fprintf(stderr, "* Start of SplitScaffolds\n");

  InitGraphNodeIterator(&scaffolds, graph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    int component;
    int numNodes, *nodes, *nodesPtr, *nodesEnd;
    NodeCGW_T *thisNode;
    CIScaffoldTIterator scaffoldNodes;

    if(scaffold->id == firstSplitScaffold){
      break;
    }
    if(isDeadCIScaffoldT(scaffold) || scaffold->type != REAL_SCAFFOLD){
      continue;
    }
    numNodes = scaffold->info.Scaffold.numElements;
    if(numNodes < 3){
      continue;
    }
    fprintf(stderr,"TESTING! Scaffold " F_CID " with %d nodes is splitting into 2 Scaffolds\n",
	    scaffold->id, numNodes);
    nodes = (int *)safe_malloc(numNodes * sizeof(int));
    AssertPtr(nodes);
    nodesEnd = nodes + numNodes;
    nodesPtr = nodes;
    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &scaffoldNodes);
    if((thisNode = NextCIScaffoldTIterator(&scaffoldNodes)) == NULL)
      assert(0);
    *nodesPtr++ = thisNode->id;
    thisNode->setID = 0;
    if((thisNode = NextCIScaffoldTIterator(&scaffoldNodes)) == NULL)
      assert(0);
    *nodesPtr++ = thisNode->id;
    thisNode->setID = 1;
    for(; (thisNode = NextCIScaffoldTIterator(&scaffoldNodes)) != NULL;
	nodesPtr++){
      assert(nodesPtr < nodesEnd);
      *nodesPtr = thisNode->id;
      thisNode->setID = rand() % 2;
    }
    scaffold->flags.bits.isDead = TRUE;  // Mark the old scaffold dead
    for(component = 0; component < 2; component++){
      LengthT NullLength = {0.0, 0.0};
      LengthT firstOffset = {0.0, 0.0};
      int seenFirstOffset;
      CIScaffoldT CIScaffold = {0};
      CDS_CID_t newScaffoldID;
      InitializeScaffold(&CIScaffold, REAL_SCAFFOLD);
      CIScaffold.info.Scaffold.AEndCI = NULLINDEX;
      CIScaffold.info.Scaffold.BEndCI = NULLINDEX;
      CIScaffold.info.Scaffold.numElements = 0;
      CIScaffold.edgeHead = NULLINDEX;
      CIScaffold.bpLength = NullLength;
      newScaffoldID = CIScaffold.id = GetNumGraphNodes(graph->ScaffoldGraph);
      CIScaffold.flags.bits.isDead = FALSE;
      CIScaffold.numEssentialA = CIScaffold.numEssentialB = 0;
      CIScaffold.essentialEdgeB = CIScaffold.essentialEdgeA = NULLINDEX;
      AppendGraphNode(graph->ScaffoldGraph, &CIScaffold);
      for(nodesPtr = nodes, seenFirstOffset = FALSE; nodesPtr < nodesEnd;
	  nodesPtr++){
	NodeCGW_T *thisNode = GetGraphNode(graph->RezGraph, *nodesPtr);
	if(thisNode->setID == component){
	  LengthT offsetAEnd, offsetBEnd;
	  if(!seenFirstOffset){
	    if(GetNodeOrient(thisNode) == A_B){
	      firstOffset = thisNode->offsetAEnd;
	    }else{
	      firstOffset = thisNode->offsetBEnd;
            }
	    seenFirstOffset = TRUE;
	  }
	  offsetAEnd.mean = thisNode->offsetAEnd.mean - firstOffset.mean;
	  offsetAEnd.variance = thisNode->offsetAEnd.variance -
	    firstOffset.variance;
	  offsetBEnd.mean = thisNode->offsetBEnd.mean - firstOffset.mean;
	  offsetBEnd.variance = thisNode->offsetBEnd.variance -
	    firstOffset.variance;
	  InsertCIInScaffold(graph, thisNode->id, newScaffoldID,
			     offsetAEnd, offsetBEnd, TRUE, FALSE);
	}
      }
      if(firstSplitScaffold == NULLINDEX){
	firstSplitScaffold = newScaffoldID;
      }
      assert((GetGraphNode(graph->ScaffoldGraph,
			   newScaffoldID))->info.Scaffold.numElements > 0);
      fprintf(stderr,"New Scaffold " F_CID " with %d nodes\n",
              newScaffoldID,
	      (GetGraphNode(graph->ScaffoldGraph,
			    newScaffoldID))->info.Scaffold.numElements);
    }
    safe_free(nodes);
  }
  return;
}

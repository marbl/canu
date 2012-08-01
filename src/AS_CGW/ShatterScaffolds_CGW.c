
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

static const char *rcsid = "$Id: ShatterScaffolds_CGW.c,v 1.2 2012-08-01 02:23:38 brianwalenz Exp $";

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

#if 0
void PrintContigContents(FILE *stream, ScaffoldGraphT *graph, char *message) {
	  GraphNodeIterator   contigs;
	  ContigT	     *contig;

	  InitGraphNodeIterator(&contigs, graph->ContigGraph, GRAPH_NODE_DEFAULT);
	  while ((contig = NextGraphNodeIterator(&contigs)) != NULL) {
	    assert(contig->id >= 0);
	    assert(contig->id < GetNumGraphNodes(graph->ContigGraph));

	    MultiAlignT *ma = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);

	    fprintf(stream, "%s contig %d placed status %c with %d unitigs: ", message, contig->id, graph->tigStore->getContigStatus(contig->id), GetNumIntMultiPoss(ma->u_list));


	    for (int i = 0; i < GetNumIntMultiPoss(ma->u_list); i++) {
	      IntUnitigPos *imp = GetIntUnitigPos(ma->u_list, i);

	      fprintf(stream, "%d ", imp->ident);
	    }
	    fprintf(stream, "\n");
	  }
}

void PrintScaffoldContents(FILE *stream, ScaffoldGraphT *graph, char *message) {
	// output initial scaffolds before breaking
		  GraphNodeIterator   scaffolds;
		  CIScaffoldT        *scaffold;

		  InitGraphNodeIterator(&scaffolds, graph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
		  while ((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL) {
		    if(scaffold->type != REAL_SCAFFOLD)
		      continue;

		    assert(scaffold->info.Scaffold.numElements > 0);

		    fprintf(stream, "%s scaffold %d contains %d elements: ", message, scaffold->id, scaffold->info.Scaffold.numElements);
		    CIScaffoldTIterator	    contigs;
		    ChunkInstanceT         *contig;


		    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &contigs);
		    while ((contig = NextCIScaffoldTIterator(&contigs)) != NULL) {
		    	// of course why would a sanity check work in cgw
		    	// line 2129 of TransitiveReduction_CGW.c sets the contig's scaffold ID to NULL but doesn't remove it from the scaffold!
		    	if (contig->scaffoldID != scaffold->id) {
		    		fprintf(stream, "Scaffold %d thinks it contains contig %d but contig thinks it belongs to %d\n", scaffold->id, contig->id, contig->scaffoldID);
		    	}
		    	assert(contig->scaffoldID == scaffold->id || contig->scaffoldID == NULLINDEX);

		    	fprintf(stream, "%d (%f) (%f) (%f) ", contig->id, contig->bpLength.mean, contig->offsetAEnd.mean, contig->offsetBEnd.mean);

          GraphEdgeIterator edges(ScaffoldGraph->ContigGraph, contig->id, ALL_END, ALL_EDGES);
		        CIEdgeT          *edge;

		        while ((edge = edges.nextMerged()) != NULL) {
		        	NodeCGW_T *otherNode = GetGraphNode(graph->ContigGraph, (edge->idA == contig->id ? edge->idB : edge->idA));
		        	//fprintf(stream, "The edge active %d is unique %d confirming %d ori:%c dist: %f connects %d and %d with weight %d and node is in scaffold %d and other end of edge is %d in scaffold %d\n", edge->flags.bits.isActive, edge->flags.bits.isUniquetoUnique, edge->flags.bits.isContigConfirming, edge->orient.toLetter(), edge->distance.mean, edge->idA, edge->idB, edge->edgesContributing, contig->scaffoldID, otherNode->id, otherNode->scaffoldID);
		        }
		    }
		    fprintf(stream, "\n");
		  }

		  PrintContigContents(stream, graph, message);
}
#endif

void ShatterScaffoldsConnectedByLowWeight(FILE *stream, ScaffoldGraphT *graph, uint32 minWeight, int verbose){
  GraphNodeIterator nodes;
  NodeCGW_T        *node;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while ((node = NextGraphNodeIterator(&nodes)) != NULL) {
    GraphEdgeIterator edges(ScaffoldGraph->ContigGraph, node->id, ALL_END, ALL_EDGES);
    CIEdgeT          *edge;

    int disconnected = (edges.nextMerged() == NULL ? FALSE : TRUE);	// don't disconnect a node if it has no edges

    while ((edge = edges.nextMerged()) != NULL) {
    	NodeCGW_T *otherNode = GetGraphNode(graph->ContigGraph, (edge->idA == node->id ? edge->idB : edge->idA));
    	if (verbose == TRUE)
    		fprintf(stream, "The edge ori:%c dist: %f connects %d and %d with weight %d and node is in scaffold %d and other end of edge is %d in scaffold %d\n", edge->orient.toLetter(), edge->distance.mean, edge->idA, edge->idB, edge->edgesContributing, node->scaffoldID, otherNode->id, otherNode->scaffoldID);

    	if (otherNode->scaffoldID == node->scaffoldID && otherNode->scaffoldID != NULLINDEX && edge->edgesContributing >= minWeight) {
    		disconnected = FALSE;
    		if (verbose == TRUE) {
    			fprintf(stream, "Node %d will not be disconnected from scaffold %d because it has edge %d higher than min %d\n", node->id, node->scaffoldID, edge->edgesContributing, minWeight);
    		}
    	}
    }

    if (disconnected == TRUE && node->scaffoldID != NULLINDEX) {
		if (verbose == TRUE)
			fprintf(stream, "Disconnecting contig with id %d from scaffold %d\n", node->id, node->scaffoldID);
		fprintf(stream, "Disconnecting contig with id %d from scaffold %d\n", node->id, node->scaffoldID);

		// is this all we need to do
		// don't set any of the flags for repeats, let it stay whatever it is now
        LengthT     offsetAEnd      = {0.0, 0.0};
        LengthT     offsetBEnd      = {0.0, 0.0};
        LengthT     firstOffset     = {0.0, 0.0};

        CIScaffoldT CIScaffold;
        InitializeScaffold(&CIScaffold, REAL_SCAFFOLD);
        CIScaffold.info.Scaffold.AEndCI = NULLINDEX;
        CIScaffold.info.Scaffold.BEndCI = NULLINDEX;
        CIScaffold.info.Scaffold.numElements = 0;
        CIScaffold.edgeHead = NULLINDEX;
        CIScaffold.bpLength = node->bpLength;
        CIScaffold.id = GetNumGraphNodes(graph->ScaffoldGraph);
        CIScaffold.flags.bits.isDead = FALSE;
        CIScaffold.numEssentialA = CIScaffold.numEssentialB = 0;
        CIScaffold.essentialEdgeB = CIScaffold.essentialEdgeA = NULLINDEX;
        AppendGraphNode(graph->ScaffoldGraph, &CIScaffold);

        node->numEssentialA = node->numEssentialB = 0;
        node->essentialEdgeA = node->essentialEdgeB = NULLINDEX;

        if(GetNodeOrient(node).isForward()){
          firstOffset = node->offsetAEnd;
        }else{
          firstOffset = node->offsetBEnd;
        }
        offsetAEnd.mean     = node->offsetAEnd.mean     - firstOffset.mean;
        offsetAEnd.variance = node->offsetAEnd.variance - firstOffset.variance;
        offsetBEnd.mean     = node->offsetBEnd.mean     - firstOffset.mean;
        offsetBEnd.variance = node->offsetBEnd.variance - firstOffset.variance;

        if (verbose == TRUE) {
			fprintf(stream, "Inserted node %d into scaffold %d at offsets (%f, %f) and (%f, %f) it used to be (%f, %f) and (%f, %f)\n",
					node->id, CIScaffold.id,
					offsetAEnd.mean, offsetAEnd.variance, offsetBEnd.mean, offsetBEnd.variance,
        			node->offsetAEnd.mean, node->offsetAEnd.variance, node->offsetBEnd.mean, node->offsetBEnd.variance);
        }

        CIScaffoldT *scaffold = GetGraphNode(graph->ScaffoldGraph, node->scaffoldID);
        RemoveCIFromScaffold(graph, scaffold, node, FALSE);
        if (scaffold->info.Scaffold.numElements == 0) {
        	scaffold->type = SCRATCH_SCAFFOLD;
			scaffold->flags.bits.isDead = 1;
        }
        InsertCIInScaffold(graph, node->id, CIScaffold.id, offsetAEnd, offsetBEnd, TRUE, FALSE);
    }
  }
}

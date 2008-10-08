
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
static char *rcsid = "$Id: DiagnosticsCGW.c,v 1.14 2008-10-08 22:02:55 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_interval.h"
#include "AS_UTL_heap.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "UnionFind_AS.h"
#include "Globals_CGW.h"
#include "AS_ALN_aligners.h"


void CheckSmallScaffoldGaps(ScaffoldGraphT *graph){
  CIScaffoldT *scaffold;
  GraphNodeIterator scaffolds;
  int overlapsConfirmed = 0;
  CDS_COORD_t overlapDistances = 0;
  int overlapsFound = 0;
  int overlapsChecked = 0;
  CDS_COORD_t gapDistances = 0;
  int gapsChecked = 0;
  int currEnd, nextEnd;
  char buffer[2048];
  FILE *fastaFile = NULL;
  LengthT gap;

  sprintf(buffer,"%s.fasta",graph->name);
  fastaFile = fopen(buffer,"w");

  fprintf(stderr,"* CheckSmallScaffoldGaps\n");
  InitGraphNodeIterator(&scaffolds, graph->ScaffoldGraph,  GRAPH_NODE_DEFAULT);


  while(NULL != (scaffold = NextGraphNodeIterator(&scaffolds))){
    CIScaffoldTIterator contigs;
    ContigT *curr, *next;
    ChunkOrientationType pairwiseOrient;

    fprintf(stderr,"* scaffold " F_CID "\n", scaffold->id);
    if(scaffold->info.Scaffold.numElements < 2)
      continue;

    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &contigs);

    curr = NextCIScaffoldTIterator(&contigs);
    while(NULL != (next = NextCIScaffoldTIterator(&contigs))){
      CDS_COORD_t actual = IntervalsOverlap(next->offsetAEnd.mean,
                                            next->offsetBEnd.mean,
                                            curr->offsetAEnd.mean,
                                            curr->offsetBEnd.mean, -15000);
      IntUnitigPos *nextUnitigPos, *currUnitigPos;
      MultiAlignT *currMA, *nextMA;

      /* If this is an overlap that will be contigged...skip it */
      if(actual >= CGW_DP_MINLEN)
        continue;

      currMA = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, curr->id, ScaffoldGraph->RezGraph->type == CI_GRAPH);
      nextMA = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, next->id, ScaffoldGraph->RezGraph->type == CI_GRAPH);

      gapsChecked++;
      if(actual > 0)
        overlapsChecked++;

      pairwiseOrient = AB_AB;
      /**** THIS WORKS ONLY IF THERE ARE NO CONTAINMENTS!!! ***/
      if(curr->offsetAEnd.mean < curr->offsetBEnd.mean){

        if(next->offsetAEnd.mean < next->offsetBEnd.mean){
          pairwiseOrient = AB_AB;
          currEnd = B_END;
          nextEnd = A_END;
          currUnitigPos = GetBendUnitigPos(currMA);
          nextUnitigPos = GetAendUnitigPos(nextMA);
          gap.mean = next->offsetAEnd.mean - curr->offsetBEnd.mean;
          gap.variance = next->offsetAEnd.variance - curr->offsetAEnd.variance;
        }else{
          pairwiseOrient = AB_BA;
          currEnd = B_END;
          nextEnd = B_END;
          currUnitigPos = GetBendUnitigPos(currMA);
          nextUnitigPos = GetBendUnitigPos(nextMA);
          gap.mean = next->offsetBEnd.mean - curr->offsetBEnd.mean;
          gap.variance = next->offsetBEnd.variance - curr->offsetBEnd.variance;
        }
      }else{
        if(next->offsetAEnd.mean < next->offsetBEnd.mean){
          pairwiseOrient = BA_AB;
          currEnd = A_END;
          nextEnd = A_END;
          currUnitigPos = GetAendUnitigPos(currMA);
          nextUnitigPos = GetAendUnitigPos(nextMA);
          gap.mean = next->offsetAEnd.mean - curr->offsetAEnd.mean;
          gap.variance = next->offsetAEnd.variance - curr->offsetAEnd.variance;
        }else{
          currEnd = A_END;
          nextEnd = B_END;
          currUnitigPos = GetAendUnitigPos(currMA);
          nextUnitigPos = GetBendUnitigPos(nextMA);
          pairwiseOrient = BA_BA;
          gap.mean = next->offsetBEnd.mean - curr->offsetAEnd.mean;
          gap.variance = next->offsetBEnd.variance - curr->offsetAEnd.variance;
        }
      }

      {
        NodeCGW_T *currUnitig = GetGraphNode(ScaffoldGraph->CIGraph, currUnitigPos->ident);
        NodeCGW_T *nextUnitig = GetGraphNode(ScaffoldGraph->CIGraph, nextUnitigPos->ident);
        if(currUnitig->type == RESOLVEDREPEATCHUNK_CGW &&
           nextUnitig->type == RESOLVEDREPEATCHUNK_CGW &&
           currUnitig->info.CI.baseID == nextUnitig->info.CI.baseID){
          fprintf(stderr,"**** >>> SURROGATE DUPLICATION!!! Unitigs on either side of this gap are surrogates (" F_CID " and " F_CID ") derived from unitig " F_CID "\n",
                  currUnitig->id, nextUnitig->id, currUnitig->info.CI.baseID);

        }
      }

      fprintf(stderr,"\n\n\n########### sc:" F_CID " (" F_CID "," F_CID ",%c) gap = (%g +/- %g)\n",
              scaffold->id,curr->id, next->id, pairwiseOrient, gap.mean, sqrt(gap.variance + 0.1));
      fprintf(stderr,"* Scaffold " F_CID " Contigs " F_CID " and " F_CID " oriented %c overlap in scaffold coordinates by " F_COORD "\n",
              scaffold->id,curr->id, next->id, pairwiseOrient,actual);

      {
        EdgeCGW_T *edge = FindGraphEdge(graph->RezGraph, curr->id, next->id, pairwiseOrient);
        CDS_COORD_t overlap;


        while(edge){
          PrintGraphEdge(stderr, graph->RezGraph, " Found associated edge ", edge, curr->id);
          if(edge->nextALink == edge->nextBLink &&
             edge->nextALink != NULLINDEX)
            edge = GetGraphEdge(graph->RezGraph, edge->nextALink);
          else
            break;
        }

        overlap = SmallOverlapExists(graph->RezGraph, curr->id, next->id, pairwiseOrient, CGW_DP_DESPERATION_MINLEN);

        if(overlap){
          fprintf(stderr,">>>>> YES SMALL Overlap of " F_COORD " found  gap mean: %g std: %g <<<<<<\n",
                  overlap, gap.mean, sqrt(gap.variance));
          curr = next;

          if(actual > 0){
            overlapsConfirmed++;
            overlapDistances += actual;
          }else{
            overlapsFound++;
            gapDistances += actual;
          }
        }else{
          fprintf(stderr,">>>>> No SMALL Overlap found <<<<<<\n");
        }

        overlap = LargeOverlapExists(graph->RezGraph, curr->id, next->id, pairwiseOrient, CGW_DP_MINLEN, 750);
        if(overlap){
          fprintf(stderr,">>>>> YES LARGE Overlap of " F_COORD " found  gap mean: %g std: %g <<<<<<\n",
                  overlap, gap.mean, sqrt(gap.variance));
          curr = next;

          if(actual > 0){
            overlapsConfirmed++;
            overlapDistances += actual;
          }else{
            overlapsFound++;
            gapDistances += actual;
          }
        }else{
          fprintf(stderr,">>>>> No LARGE Overlap found <<<<<<\n");
        }



	{
	  GraphEdgeIterator Edges;
	  EdgeCGW_T *edge;

	  fprintf(stderr,"* Node " F_CID " has the following overlaps on its %c end:\n",
		  curr->id, (currEnd == A_END?'A':'B'));
	  InitGraphEdgeIterator(graph->RezGraph, curr->id, currEnd, ALL_EDGES, GRAPH_EDGE_RAW_ONLY | GRAPH_EDGE_DEFAULT , &Edges);
	  while(NULL != (edge= NextGraphEdgeIterator(&Edges))){
	    if(!isOverlapEdge(edge) || isContainmentEdge(edge))
	      continue;
	    PrintGraphEdge(stderr, graph->RezGraph, "\t", edge, edge->idA);
	  }

	  fprintf(stderr,"* Node " F_CID " has the following overlaps on its %c end:\n",
		  next->id, (nextEnd == A_END?'A':'B'));
	  InitGraphEdgeIterator(graph->RezGraph, next->id, nextEnd, ALL_EDGES, GRAPH_EDGE_RAW_ONLY | GRAPH_EDGE_DEFAULT , &Edges);
	  while(NULL != (edge= NextGraphEdgeIterator(&Edges))){
	    if(!isOverlapEdge(edge) || isContainmentEdge(edge))
	      continue;
	    PrintGraphEdge(stderr, graph->RezGraph, "\t", edge, edge->idA);
	  }
	}
      }
      curr = next;
    }
  }
  fprintf(stderr,"* Gaps Checked %d Overlaps checked %d confirmed %d (avg " F_COORD ")    Overlaps found %d (avg " F_COORD ")\n",
          gapsChecked, overlapsChecked, overlapsConfirmed,  (overlapsConfirmed?overlapDistances/overlapsConfirmed:0), overlapsFound, (overlapsFound?gapDistances/overlapsFound:0));
}


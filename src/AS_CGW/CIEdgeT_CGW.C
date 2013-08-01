
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
static char *rcsid = "$Id: CIEdgeT_CGW.c,v 1.26 2012-08-28 21:09:39 brianwalenz Exp $";

//#define DEBUG 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.H"
#include "AS_UTL_Var.H"
#include "AS_UTL_interval.H"
#include "AS_CGW_dataTypes.H"
#include "Globals_CGW.H"
#include "ScaffoldGraph_CGW.H"
#include "ScaffoldGraphIterator_CGW.H"
#include "ChunkOverlap_CGW.H"




void PrintCIEdgeT(FILE *fp, ScaffoldGraphT *graph,
                  char *label, CIEdgeT *edge, CDS_CID_t cid){
  char actualOverlap[256];
  int32 delta = 0;
  char *flag = "  ", *flagTrans = "  ";
  ChunkInstanceT *ChunkInstanceA =
    GetChunkInstanceT(graph->ChunkInstances, edge->idA);
  ChunkInstanceT *ChunkInstanceB =
    GetChunkInstanceT(graph->ChunkInstances, edge->idB);

  assert(ChunkInstanceA && ChunkInstanceB);

  *actualOverlap = '\0';

  if(edge->flags.bits.isBogus){
    if(edge->flags.bits.isProbablyBogus)
      strcpy(actualOverlap," *Bogus and Prob Bogus*");
    else
      strcpy(actualOverlap," *Bogus*");
  }

  if(edge->flags.bits.hasContributingOverlap){
    if(edge->flags.bits.isPossibleChimera)
      flag = "$C";
    else
      flag = "$O";
  }
  if(edge->flags.bits.isEssential){
    if(edge->flags.bits.isInferred){
      flagTrans = "&I";
    }else{
      flagTrans = "&E";
    }
  }else if(edge->flags.bits.isTransitivelyRemoved){
    flagTrans = "&T";
  }else if(edge->flags.bits.isRedundantRemoved){
    flagTrans = "&R";
  }else if(edge->flags.bits.isInferredRemoved){
    flagTrans = "&t";
  }else if(edge->flags.bits.isActive){
    if(edge->flags.bits.isConfirmed){
      flagTrans = "&C";
    }else{
      flagTrans = "&A";
    }
  }else if(edge->flags.bits.isProbablyBogus){
    flagTrans = "&B";
  }

  fprintf(fp,"%s A:"F_CID" B:"F_CID" wgt:%d %s %s %s %s%s ori:%c con:%d dst:%d std:%g %s\n",
          label, edge->idA, edge->idB, edge->edgesContributing, flag,
          (edge->flags.bits.hasExtremalAFrag?"$A":"  "),
          (edge->flags.bits.hasExtremalBFrag?"$B":"  "),
          flagTrans, (edge->flags.bits.isLeastSquares?"L":" "),
          GetEdgeOrientationWRT(edge, cid).toLetter(),
          edge->flags.bits.hasContainmentOverlap,
          (int)edge->distance.mean, sqrt(edge->distance.variance),
          actualOverlap);

}

void PrintChunkInstanceHeader(FILE *stream, ScaffoldGraphT *graph,
                              ChunkInstanceT *chunk){
  fprintf(stream,"\n* CI "F_CID" cov:%d len:%d frags:%d\n",
          chunk->id,
          ScaffoldGraph->tigStore->getUnitigCoverageStat(chunk->id),
          (int)chunk->bpLength.mean,
          ScaffoldGraph->tigStore->getNumFrags(chunk->id, TRUE));
}




void DumpChunkInstance(FILE *stream, ScaffoldGraphT *graph,
                       ChunkInstanceT *chunk,
		       int confirmedOnly, int scaffoldedOnly,
		       int uniqueToUniqueOnly, int verbose){
  int aEndPrinted, bEndPrinted;


  PrintChunkInstanceHeader(stream, graph, chunk);

  if(scaffoldedOnly && chunk->scaffoldID == NULLINDEX)
    return;

  aEndPrinted = FALSE;
  bEndPrinted = FALSE;

  /* Iterate over A_END edges */
  {
    GraphEdgeIterator edgeMates(graph->CIGraph, chunk->id, A_END, ALL_EDGES);
    CIEdgeT          *edge;

    while((edge = edgeMates.nextMerged()) != NULL){
      if(uniqueToUniqueOnly && !edge->flags.bits.isUniquetoUnique){
        continue;
      }
      if(!aEndPrinted){
        PrintChunkInstanceHeader(stream,graph, chunk);
        fprintf(stream,"\n\tA_END\n");
        aEndPrinted = TRUE;
      }
      PrintGraphEdge(stream,graph->CIGraph,"\t", edge, chunk->id);
    }
  }

  /* Iterate over B_END edges */
  {
    GraphEdgeIterator edgeMates(graph->CIGraph, chunk->id, B_END, ALL_EDGES);
    CIEdgeT          *edge;

    while((edge = edgeMates.nextMerged()) != NULL){
      if(uniqueToUniqueOnly && !edge->flags.bits.isUniquetoUnique){
        continue;
      }
      if(!aEndPrinted){
        PrintChunkInstanceHeader(stream, graph, chunk);
        aEndPrinted = TRUE;
      }
      if(!bEndPrinted){
        fprintf(stream,"\n\tB_END\n");
        bEndPrinted = TRUE;
      }
      PrintGraphEdge(stream,graph->CIGraph,"\t", edge, chunk->id);
    }
  }
}
/****************************************************************************/
void DumpChunkInstances(FILE *stream, ScaffoldGraphT *graph, int confirmedOnly,
                        int scaffoldedOnly, int uniqueToUniqueOnly,
                        int verbose){
  CDS_CID_t cid;
  /* For each chunk */
  for(cid = 0; cid < GetNumGraphNodes(graph->CIGraph); cid++){
    ChunkInstanceT *chunk = GetGraphNode(graph->CIGraph,cid);
    DumpChunkInstance(stream, graph, chunk, confirmedOnly, scaffoldedOnly,
                      uniqueToUniqueOnly, verbose);
  }
}

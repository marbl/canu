
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
static char *rcsid = "$Id: CIEdgeT_CGW.c,v 1.17 2008-10-29 10:42:46 brianwalenz Exp $";

//#define DEBUG 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"





void PrintCIEdgeT(FILE *fp, ScaffoldGraphT *graph,
                  char *label, CIEdgeT *edge, CDS_CID_t cid){
  char actualOverlap[256];
  CDS_COORD_t delta = 0;
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

  fprintf(fp,"%s A:" F_CID " B:" F_CID " wgt:%d %s %s %s %s%s ori:%c con:%d dst:%d std:%g %s\n",
          label, edge->idA, edge->idB, edge->edgesContributing, flag,
          (edge->flags.bits.hasExtremalAFrag?"$A":"  "),
          (edge->flags.bits.hasExtremalBFrag?"$B":"  "),
          flagTrans, (edge->flags.bits.isLeastSquares?"L":" "),
          GetEdgeOrientationWRT(edge, cid),
          edge->flags.bits.hasContainmentOverlap,
          (int)edge->distance.mean, sqrt(edge->distance.variance),
          actualOverlap);

}

void PrintChunkInstanceHeader(FILE *stream, ScaffoldGraphT *graph,
                              ChunkInstanceT *chunk){
  fprintf(stream,"\n* CI " F_CID " cov:%d len:%d frags:%d\n",
          chunk->id,
          chunk->info.CI.coverageStat,
          (int)chunk->bpLength.mean,
          chunk->info.CI.numFragments);
}




void DumpChunkInstance(FILE *stream, ScaffoldGraphT *graph,
                       ChunkInstanceT *chunk,
		       int confirmedOnly, int scaffoldedOnly,
		       int uniqueToUniqueOnly, int verbose){
  int aEndPrinted, bEndPrinted;
  CIEdgeT *edge;
  CIEdgeTIterator edgeMates;

  if(chunk->edgeHead == NULLINDEX){
    PrintChunkInstanceHeader(stream, graph, chunk);
    return;
  }
  if(scaffoldedOnly && chunk->scaffoldID == NULLINDEX)
    return;

  aEndPrinted = FALSE;
  bEndPrinted = FALSE;
  /* Iterate over A_END edges */
  InitCIEdgeTIterator(graph, chunk->id, FALSE, confirmedOnly,
                      A_END, ALL_EDGES, verbose, &edgeMates);
  while((edge = NextCIEdgeTIterator(&edgeMates)) != NULL){
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
  /* Iterate over B_END edges */
  InitCIEdgeTIterator(graph, chunk->id, FALSE, confirmedOnly,
                      B_END, ALL_EDGES, verbose,  &edgeMates);
  while((edge = NextCIEdgeTIterator(&edgeMates)) != NULL){
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
  return;
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








/****************************************************************************/
/*
   Look for overlaps that are implicit in the extended chunk graph.
   If a unique 'seed' chunk has confirmed links with two chunks on the
   same side of the seed, then these links imply a relationship between
   the two chunks.  If this relationship implies a potential overlap,
   we want to check it directly
*/
int CheckImplicitOverlaps_(GraphCGW_T *graph, CDS_CID_t cid, int end){
  GraphEdgeIterator sources;
  EdgeCGW_T *source, *target;

  InitGraphEdgeIterator(graph,
                        cid,
                        end,
                        ALL_EDGES,
                        GRAPH_EDGE_CONFIRMED_ONLY,
                        &sources);

  while(NULL != (source = NextGraphEdgeIterator(&sources))){
    ChunkInstanceT *sourceChunk = GetGraphNode(graph, source->idB);
    GraphEdgeIterator targets = sources;
    while(NULL != (target = NextGraphEdgeIterator(&targets))){
      ChunkInstanceT *targetChunk = GetGraphNode(graph, target->idB);
      /* Check out the relationship between source and target */
      LengthT seedSource, seedTarget, sourceTarget;
      double sourceTargetDelta;
      int sourceCloserThanTarget;
      CDS_CID_t targetCid, sourceCid;

      /* Eliminate self references */
      if(source->idB == target->idB)
        continue;

      // Only look unique to unique
      if(!targetChunk->flags.bits.isUnique)
        continue;

      seedSource = source->distance;
      seedTarget = target->distance;

      // Always work with the chunk that is closest as the source
      sourceCloserThanTarget = seedSource.mean <= seedTarget.mean;

      if(sourceCloserThanTarget){
        targetCid = target->idB;
        sourceCid = source->idA;
      }else{
        targetCid = source->idA;
        sourceCid = target->idB;
        seedSource = target->distance;
        seedTarget = source->distance;
      }
      sourceChunk = GetGraphNode(graph, sourceCid);

#ifdef DEBUG_CGW
      fprintf(stderr,"* %c side source=" F_CID " target=" F_CID " seedSource = %d seedTarget = %d sourceLength=%d\n",
              (sources.end == 1?'A':'B'),
              sourceCid, targetCid,
              (int) seedSource.mean, (int) seedTarget.mean,
              (int) sourceChunk->bpLength.mean);
#endif

      sourceCloserThanTarget = seedSource.mean <= seedTarget.mean;

      // This is the way the edgeMates are sorted
      assert(sourceCloserThanTarget);

      sourceTarget.mean = seedTarget.mean -
        (seedSource.mean + sourceChunk->bpLength.mean);
      sourceTarget.variance = seedTarget.variance +
        seedSource.variance + sourceChunk->bpLength.variance;
      sourceTargetDelta = 3.0 * sqrt(sourceTarget.variance);

#ifdef DEBUG_CGW
      fprintf(stderr,"* sourceTarget.mean= %d variance=%g\n",
              (int) sourceTarget.mean, sourceTarget.variance);
#endif
      if(sourceTarget.mean - sourceTargetDelta <  CGW_DP_MINLEN){
        /* Figure out orientation of chunks overlap ranges, and
           CollectOverlap on the combination */
        ChunkOrient sourceOrient =
          GetRelativeChunkOrientation(GetEdgeOrientationWRT(source,
                                                            cid ));
        ChunkOrient targetOrient =
          GetRelativeChunkOrientation(GetEdgeOrientationWRT(target,
                                                            cid));
        ChunkOrientationType sourceTargetOrient =
          GetChunkPairOrientation(sourceOrient, targetOrient);



        assert(GetChunkSeedSide(GetEdgeOrientationWRT(target, cid)) == GetChunkSeedSide(GetEdgeOrientationWRT(source, cid)));
        {
          CDS_COORD_t delta = 3 * sqrt(sourceTarget.variance);
          CDS_COORD_t minOverlap =
            MAX(-sourceTarget.mean - delta, 0);
          CDS_COORD_t maxOverlap = -sourceTarget.mean + delta;
          ChunkOverlapCheckT olap = {0};

          int overlapCheckFound = LookupOverlap(graph,
                                                sourceCid, targetCid,
                                                sourceTargetOrient,
                                                &olap);
          if(overlapCheckFound)
            continue;

          fprintf(stderr,"** Found an implied potential overlap (" F_CID "," F_CID ") min:" F_COORD " max:" F_COORD "\n",
                  sourceCid, targetCid,
                  minOverlap, maxOverlap);

          assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));

          //handle suspicious
          olap = OverlapChunks(graph,
                               sourceCid, targetCid,
                               sourceTargetOrient,
                               minOverlap,
                               maxOverlap,
                               AS_CGW_ERROR_RATE,
                               TRUE);

          if(olap.overlap > 0){
            fprintf(stderr,"* Found Implied Overlap (" F_CID "," F_CID ",%c) with overlap " F_COORD "\n",
                    sourceCid, targetCid, sourceTargetOrient,
                    olap.overlap);
            if(olap.suspicious){
              fprintf(stderr,"* OVERLAP is SUSPICIOUS!!!!!\n");
            }
          }
        }
      }
    }
  }
  return TRUE;
}


/* AddImplicitOverlaps
   For all unique chunks, check for implicit overlaps on both
   sides of the unique seed.
*/

int AddImplicitOverlaps_(GraphCGW_T *graph){
  CDS_CID_t cid;
  GraphNodeIterator nodes;
  NodeCGW_T *chunk;
  InitGraphNodeIterator(&nodes, graph,GRAPH_NODE_UNIQUE_ONLY);
  while(NULL != (chunk = NextGraphNodeIterator(&nodes))){
    cid = chunk->id;

    /* For each pair of chunks to the A side of this 'seed'
       that don't have any links between them and potentially
       overlap, collect an overlap record for later evaluation */
    CheckImplicitOverlaps_(graph, cid, A_END);

    /* For each pair of chunks to the B side of this 'seed'
       that don't have any links between them and potentially
       overlap, collect an overlap record for later evaluation */

    CheckImplicitOverlaps_(graph, cid, B_END);

  }
  return(0);
}

/***************************************************************************/



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

        Module:  UpdateREZ.c

   Description:  Contains functions that take a Scaffold_Fill_t and
                 update the Scaffold Graph

		 log info is sent to <inputFile>.cgwlog

    Programmer:  S. Lonardi (stelo@cs.purdue.edu)
 
       Written:  17 May 99

 **********************************************************************/

static char fileID[] = "$Id: UpdateREZ.c,v 1.5 2006-06-14 19:57:23 brianwalenz Exp $";

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UnionFind_AS.h"
//
// AS_CGW
//
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#ifdef CREATE_CHUNK_GRAPH
#include "AS_CGW_Scaffold.h"
#include "AS_CGW_Scaffolds.h"
#endif
#include "ScaffoldGraph_CGW.h"
#include "ChiSquareTest_CGW.h"
//
// AS_REZ
//
#include "DataTypesREZ.h"
#include "UtilsREZ.h"
#include "CommonREZ.h"
#include "GapWalkerREZ.h"
#include "BccREZ.h"
#include "SubgraphREZ.h"
#include "UpdateREZ.h"
#include "ConsistencyChecksREZ.h"

int Trust_Edges(ScaffoldGraphT * sgraph,
				Scaffold_Fill_t * fill_chunks,
				Gap_Chunk_t * * table,
				int32 sid,
				int32 cid) {
  //
  // mark all the edges that the CI <cid> has to uniques and others CI
  // in the fill_chunks as tentatively trusted: returns the number of
  // RAW edge marked
  //
  GraphEdgeIterator
    edges;
  CIEdgeT
    * edge;
  int32
    next;
  int
    marked = 0;
  ChunkInstanceT 
    * next_chunk;
# if DEBUG_UPDATE > 1
  char
    buffer[256];
# endif

  InitGraphEdgeIterator(sgraph->RezGraph, cid,
						ALL_END, ALL_EDGES, FALSE, &edges);
  while((edge = NextGraphEdgeIterator(&edges)) != NULL){
    assert(edge != NULL);
    //
    // get the other end
    //
    if (cid == edge->idA)
      next = edge->idB;
    else
      next = edge->idA;

    next_chunk = GetGraphNode(sgraph->RezGraph, next);
    assert(next_chunk != NULL);

    if ((Is_Unique(next_chunk) && (next_chunk->scaffoldID == sid)) ||
        ((table[next] != NULL) && (table[next]->scaff_id == sid) && (table[next]->keep))) {
      //
      // mark the edge only if the distance is consistent
      //
      if (Is_Unique(next_chunk)) {
		Gap_Chunk_t
		  unique;
		unique.chunk_id = next_chunk->id;
		unique.scaff_id = next_chunk->scaffoldID;
		unique.start = next_chunk->offsetAEnd;
		unique.end = next_chunk->offsetBEnd;
	
		if (Is_Edge_Consistent(edge, table[cid], &unique)) {
		  if (Is_Edge_Orientation_Consistent(edge, table[cid], &unique)) {
			SetEdgeStatus(sgraph->RezGraph, edge, TENTATIVE_TRUSTED_EDGE_STATUS);
#           if DEBUG_UPDATE > 1
			sprintf(buffer, "pass Knut's test & orientation test * TRUSTED *\n");
#           endif
			marked += edge->edgesContributing;
		  } else {
			SetEdgeStatus(sgraph->RezGraph, edge, TENTATIVE_UNTRUSTED_EDGE_STATUS);
#           if DEBUG_UPDATE > 1
			sprintf(buffer, "pass Knut's test & but FAILED orientation test * UNTRUSTED *\n");
#           endif
		  }
		} else {
		  SetEdgeStatus(sgraph->RezGraph, edge, TENTATIVE_UNTRUSTED_EDGE_STATUS);
#         if DEBUG_UPDATE > 1
		  sprintf(buffer, "FAILED Knut's test * UNTRUSTED *\n");
#         endif
		}
      } else {
		if (Is_Edge_Consistent(edge, table[cid], table[next])) {
		  if (Is_Edge_Orientation_Consistent(edge, table[cid], table[next])) {
			SetEdgeStatus(sgraph->RezGraph, edge, TENTATIVE_TRUSTED_EDGE_STATUS);
#           if DEBUG_UPDATE > 1
			sprintf(buffer, "pass Knut's test & orientation test * TRUSTED *\n");
#           endif
			marked += edge->edgesContributing;
		  } else {
			SetEdgeStatus(sgraph->RezGraph, edge, TENTATIVE_UNTRUSTED_EDGE_STATUS);
#           if DEBUG_UPDATE > 1
			sprintf(buffer, "pass Knut's test & but FAILED orientation test * UNTRUSTED *\n");
#           endif
		  }
		} else {
		  SetEdgeStatus(sgraph->RezGraph, edge, TENTATIVE_UNTRUSTED_EDGE_STATUS);
#         if DEBUG_UPDATE > 1
		  sprintf(buffer, "FAILED Knut's test * UNTRUSTED *\n");
#         endif
		}
      }
    } else {
      //
	  // mark the edge as UNTRUSTED
      //
      SetEdgeStatus(sgraph->RezGraph, edge, TENTATIVE_UNTRUSTED_EDGE_STATUS);
#     if DEBUG_UPDATE > 1
      sprintf(buffer, "UNTRUSTED *\n");
#     endif
    }
#   if DEBUG_UPDATE > 1
    fprintf(stderr,"  -=> edge %d->%d (edge weight %d, new status %d, mean %f, var %f) (other end sid %d, %s) * %s", 
			cid,
			next,
			edge->edgesContributing,
			GetEdgeStatus(edge),
			edge->distance.mean,
			edge->distance.variance,
			next_chunk->scaffoldID,
			table[next] ? "in the set" : "not in the set",
			buffer);
#   endif
  }
  return marked;
}



int Update_Scaffold_Graph(ScaffoldGraphT *sgraph,
						  Scaffold_Fill_t *fill_chunks,
						  int ChiSquare,
						  int noRecompute,
						  int markTrusted,
						  int copyAllOverlaps,
						  int scaffId,
                                                  Kind_Of_Fill_t kind) {
  //
  // This function update the scaffolds by filling the appropriate
  // gaps as indicated by the infos in the gapAssignment (only the
  // ones which has a "keep" flag TRUE will be considered). if the flag
  // RebuildSEdges is TRUE then we will rebuild SEdges. if the flag
  // ChiSquare is TRUE then we will compute the ChiSquare test at
  // each scaffold. returns the number of inserted chunks
  // if scaffId is -1 we fill in all scaffolds, otherwise the scaffold 
  // with number scaffId
  // noRecompute is TRUE placing surrogates (walking/stones) and FALSE for placing rocks
  // copyAllOverlaps is FALSE for walking, TRUE for stones 

  int
    scaff_id,
    num_chunks = GetNumGraphNodes(sgraph->RezGraph),
    num_scaffs = GetNumCIScaffoldTs(sgraph->CIScaffolds),
    j,
    inserted = 0,
    trusted_edges,
    scaffold_modified,
    k;
  Gap_Fill_t
    * fc;
  CIScaffoldT
    *scaffold;
  
#if UPDATE_PARANOID_CHECKING > 0
  CIScaffoldTIterator  CIs;
#endif
  
#if 0
  RecomputeOffsetsStatus  res;
#endif
  
  Gap_Chunk_t
    * * table = (Gap_Chunk_t * *)safe_calloc(num_chunks, sizeof(Gap_Chunk_t *));
  ChunkInstanceT 
    * chunk;

# if DEBUG_UPDATE > 1
  fprintf(stderr,"-=> UPDATE started\n");
# endif

  //
  // we need a table so that given a chunk id we can go directly to
  // the relative Gap_Chunk_t (if any). Do some checking on the scaffolds
  // too ...
  //
  for (scaff_id = 0; scaff_id < num_scaffs;  scaff_id++) {
    //
    // get the scaffold
    //

    if( scaffId != -1 )
      if( scaffId != scaff_id )
		continue;

    scaffold = GetCIScaffoldT(sgraph->CIScaffolds, scaff_id);
    assert(scaffold != NULL);
    if ((isDeadCIScaffoldT(scaffold)) || (scaffold->type != REAL_SCAFFOLD))
      continue;
#   if DEBUG_UPDATE > 1
    fprintf(stderr,"\n-=> Scaffold %3d: (LSE %g, num %d)\n",
			scaff_id,
			scaffold->info.Scaffold.leastSquareError,
			scaffold->info.Scaffold.numLeastSquareClones);
#   endif

#   if UPDATE_PARANOID_CHECKING > 0
    //
    // make sure that we have ONE component to start with
    //
    if (IsScaffoldInternallyConnected(sgraph, scaffold) != 1)
      fprintf(stderr, "* error: the scaffold was initially not internally connected\n");

    //
    // make sure that trusted edges are internal to the scaffold
    //
    InitCIScaffoldTIterator(sgraph, scaffold, TRUE,
							FALSE, &CIs);
    while (chunk = NextCIScaffoldTIterator(&CIs))
      CheckTrustedEdges(sgraph, chunk->cid);
#   endif

    for (j = 0; j < fill_chunks[scaff_id].num_gaps;  j++) {
      fc = &(fill_chunks[scaff_id].gap[j]);
#     if DEBUG_UPDATE > 1
      fprintf(stderr,"\n -=> Gap %3d:  (%5.2f, %5.2f)-(%5.2f, %5.2f)\n",
			  j,
			  fc->start.mean,
			  sqrt(fc->start.variance),
			  fc->end.mean,
			  sqrt(fc->end.variance));
#     endif
      for  (k = 0;  k < fc->num_chunks;  k++) {
		table[fc->chunk[k].chunk_id] = &(fc->chunk[k]);
#       if DEBUG_UPDATE > 1
		fprintf(stderr,"  -=> Chunk %6d: <%5.2f, %5.2f, %5.2f>-<%5.2f, %5.2f, %5.2f>-<%s>-<%s>\n",
				fc->chunk[k].chunk_id,
				fc->chunk[k].start.mean,
				sqrt(fc->chunk[k].start.variance),
				fc->chunk[k].start.variance,
				fc->chunk[k].end.mean,
				sqrt(fc->chunk[k].end.variance),
				fc->chunk[k].end.variance,
				fc->chunk[k].keep ? "to keep" : "not to keep",
				fc->chunk[k].split ? "to split" : "not to split");
#       endif
      }
    }
  }

  //
  // now iterate over all scaffold, gaps, and chunks
  //
  for (scaff_id = 0; scaff_id < num_scaffs;  scaff_id++) {
    CIScaffoldT
      *internalScaffold;
    CIScaffoldTIterator
      internalChunks;
    ChunkInstanceT 
      *internalChunk;
    //
    // get the scaffold
    //
    scaffold = GetCIScaffoldT(sgraph->CIScaffolds, scaff_id);
    assert(scaffold != NULL);
    if ((isDeadCIScaffoldT(scaffold)) || (scaffold->type != REAL_SCAFFOLD))
      continue;
    // ScaffoldSanity( scaffold, sgraph);

#   if DEBUG_UPDATE > 1
    fprintf(stderr,"\n-=> --------------------------------------\n");
    fprintf(stderr,"-=> Scaffold %3d: (LSE %g, num %d, components %d)\n",
			scaff_id,
			scaffold->info.Scaffold.leastSquareError,
			scaffold->info.Scaffold.numLeastSquareClones,
			IsScaffoldInternallyConnected(sgraph, scaffold));
    fprintf(stderr,"-=> --------------------------------------\n");
#   endif

    //
    // now loop over all the gaps in "scaffold" (scaff_id)
    //
    scaffold_modified = FALSE;
    for (j = 0; j < fill_chunks[scaff_id].num_gaps;  j++) {
      //
      // a convenient pointer
      //
      fc = &(fill_chunks[scaff_id].gap[j]);

#     if DEBUG_UPDATE > 1
      fprintf(stderr,"\n -=> Gap %3d:  (%5.2f, %5.2f)-(%5.2f, %5.2f)\n",
			  j,
			  fc->start.mean,
			  sqrt(fc->start.variance),
			  fc->end.mean,
			  sqrt(fc->end.variance));
#     endif


      /* first we collect all chunks that are in DIFFERENT scaffolds than the current one.
		 Then we check whether we walk ALL chunks in the respective scaffold.
		 If so, we mark that scaffold as dead and all the chunks to keep */

      internalScaffold = GetCIScaffoldT(sgraph->CIScaffolds, scaff_id);
      // scaffold the gap is in. Has ID scaff_id

	  //
	  // loop over the chunks that are slated to go in this gap
	  //
      for  (k = 0;  k < fc->num_chunks;  k++) 
      {
        int keep = FALSE;
        int l;
        int chunkScaffID;
        CIScaffoldT  * chunkScaffold;
        
        if  (! fc -> chunk [k] . keep)
          continue;
        
        chunkScaffID
          = GetGraphNode (ScaffoldGraph -> RezGraph,
                          fc -> chunk [k] . chunk_id) -> scaffoldID;
        chunkScaffold = GetCIScaffoldT (sgraph -> CIScaffolds, chunkScaffID);
        // scaffold the current chunk is in. Has ID chunkScaffID
        
        /* scaffold ID must be different */
#if DEBUG_UPDATE > 1
        fprintf(stderr,"chunk %d, scaff_id = %d chunkScaffID = %d \n",fc->chunk[k].chunk_id,scaff_id,chunkScaffID);
#endif
        if( (chunkScaffID != scaff_id) && (chunkScaffID != NULLINDEX) )
        {
#if DEBUG_UPDATE > 1
          fprintf(stderr,"FOUND different Scaffold %d \n",chunkScaffID);
#endif
          //
          // iterate through the chunks that are in the same scaffold "chunkScaffold"
          // as the chunk slated for insertion in the current gap
          //
          InitCIScaffoldTIterator(sgraph, chunkScaffold, TRUE,
                                  FALSE,&internalChunks);
          /* for each chunk in the scaffold we look up whether he is in the fc */
          while((internalChunk = NextCIScaffoldTIterator(&internalChunks)) != NULL)
          {
#if DEBUG_UPDATE > 1
            fprintf(stderr,"chunk id = %d, internalChunkID %d \n",fc->chunk[k].chunk_id,internalChunk->id);  
#endif
            keep = FALSE;
            // see if the this chunk from "chunkScaffold" is among the chunks slated for this gap
            for  (l=0;  l<fc->num_chunks; l++) 
              if( fc->chunk[l].chunk_id == internalChunk->id )
                keep = TRUE;
            /* one chunk was not in the fc so we do not keep any chunk in the
               walk from this scaffold */
            if( keep == FALSE )
              break;
          }
          
          // now loop through the chunks slated for this gap and mark the ones from chunkScaffold
          // according to whether or not chunkScaffold is going to be absorbed
          for  (l=0;  l<fc->num_chunks; l++) 
          { 
            ChunkInstanceT *ci = GetGraphNode(ScaffoldGraph->RezGraph, fc->chunk[l].chunk_id);	
#if DEBUG_UPDATE > 1
            fprintf(stderr,"ci->scaffoldID = %d\n",ci->scaffoldID);
#endif
            if( ci->scaffoldID == chunkScaffID )
            {
              if( keep )
              {
#if DEBUG_UPDATE > 1
                fprintf(stderr,"SET to TRUE %d \n",fc->chunk[l].chunk_id);
#endif
                RemoveCIFromScaffold(ScaffoldGraph,chunkScaffold,ci, TRUE);
                fc->chunk[l].keep = TRUE;
                fc->chunk[l].split = FALSE;
                
              }
              else
              {
                fc->chunk[l].keep = FALSE;
#if DEBUG_UPDATE > 1
                fprintf(stderr,"SET to FALSE %d \n",fc->chunk[l].chunk_id);
#endif
              }
            }
          }
          if( keep )
          {
            chunkScaffold->flags.bits.isDead = TRUE;  
# if DEBUG_GAP_WALKER > -1
            fprintf (stderr,
                     "GW:  Killed Scaffold %d and inserted its contigs into Scaffold %d\n",
                     chunkScaffID, scaff_id);
#endif
          }
        }
      }

      for  (k = 0;  k < fc->num_chunks;  k++) {
		// here we have to filter out the scaffold ends
		// that means any chunk that is in the scaffold in which the gap is
#if DEBUG_UPDATE > 1
		fprintf(stderr,"Chunk Scaff ID = %d, scaff_id= %d \n",GetGraphNode(ScaffoldGraph->RezGraph, fc->chunk[k].chunk_id)->scaffoldID, scaff_id);
#endif
 
		if( GetGraphNode(ScaffoldGraph->RezGraph, fc->chunk[k].chunk_id)->scaffoldID == scaff_id ) 
		  fc->chunk[k].keep = FALSE;

#       if DEBUG_UPDATE > 1
		fprintf(stderr,"  -=> Chunk %6d: <%5.2f, %5.2f, %5.2f>-<%5.2f, %5.2f, %5.2f>-<%s><%s>\n",
				fc->chunk[k].chunk_id,
				fc->chunk[k].start.mean,
				sqrt(fc->chunk[k].start.variance),
				fc->chunk[k].start.variance,
				fc->chunk[k].end.mean,
				sqrt(fc->chunk[k].end.variance),
				fc->chunk[k].end.variance,
				fc->chunk[k].keep ? "to keep" : "not to keep",
				fc->chunk[k].split ? "to split" : "not to split");
#       endif
	

		//
		// if the bit keep is set AND the chunk is not the from or destination chunk  
		// then we try to insert the CI
		//
		//	if( GetGraphNode(ScaffoldGraph->RezGraph, fc->chunk[k].chunk_id)->scaffoldID != NULLINDEX ) 
		//	  fprintf(stderr,"*** Skipping chunk already in Scaffold\n");
		//	else
		if (fc->chunk[k].keep )
	    {
		  int cid = fc->chunk[k].chunk_id;

#         if DEBUG_UPDATE > 1
		  fprintf(stderr,
				  "\n  -=> Before insertion of CI %d\n",
				  fc->chunk[k].chunk_id);
#         if DEBUG_UPDATE > 2
		  DumpCIScaffold(stderr,sgraph, scaffold, FALSE);
#         endif
#         endif

	 
		  // mark the edges of this chunk: if no edges are trusted,
		  // avoid inserting this chunk
		  //
		  if( markTrusted )
		  {
			trusted_edges = Trust_Edges(sgraph, fill_chunks, table, scaff_id, cid);
	      
#         if DEBUG_UPDATE > 1
			fprintf(stderr,
					"-=> TOTAL marked edges for chunk %d is %d\n",
					cid,
					trusted_edges);
#         endif
			if (trusted_edges < 1) {
			  fprintf(stderr,
					  "* warning : attempted to insert CI %d with a total weight %d of trusted edges - skipped\n",
					  cid,
					  trusted_edges);
			  continue;
			}
		  }

		  //
		  // insert the CI in the scaffold
		  //
		  /*** mjf ***/ // fprintf(stderr, "in Update_Scaffold_Graph, inserting chunk %d\n", cid);
		  /* DON'T CONTIG NOW */

		  // Here we decide now, whether the chunk gets split or not

		  if( ! fc->chunk[k].split )
		  {
                        ContigT  * contig;
			//	      DumpCIScaffold(stderr,sgraph, scaffold, FALSE);
			// ScaffoldSanity(GetGraphNode(ScaffoldGraph->ScaffoldGraph, scaff_id), sgraph);
                        contig = GetGraphNode (sgraph -> RezGraph, cid);

                        if  (contig -> info . Contig . numCI == 1)
                            {
                             ChunkInstanceT  * unitig;

                             assert (contig -> info . Contig . AEndCI
                                       == contig -> info . Contig . BEndCI);
                             unitig
                                 = GetGraphNode (ScaffoldGraph -> CIGraph,
                                                 contig -> info . Contig . AEndCI);
                             switch  (kind)
                               {
                                case  ROCKS :
                                  unitig -> flags . bits . isRock = TRUE;
                                  break;
                                case  STONES :
                                  unitig -> flags . bits . isStone = TRUE;
                                  break;
                                case  WALKS :
                                  unitig -> flags . bits . isWalk = TRUE;
                                  break;
                                default :
                                  fprintf (stderr,
                                           "ERROR:  Unexpected insert type = %d\n",
                                           (int) kind);
                               }
                            }

                        switch  (kind)
                          {
                           case  ROCKS :
                             contig -> flags . bits . isRock = TRUE;
                             break;
                           case  STONES :
                             contig -> flags . bits . isStone = TRUE;
                             break;
                           case  WALKS :
                             contig -> flags . bits . isWalk = TRUE;
                             break;
                           default :
                             fprintf (stderr, "ERROR:  Unexpected insert type = %d\n",
                                      (int) kind);
                          }
			InsertCIInScaffold(sgraph, cid, scaff_id,
							   fc->chunk[k].start, fc->chunk[k].end, TRUE, NO_CONTIGGING);
		  }
		  else
		  {
                        ContigT  * contig;
			int splitCid;

                        contig = GetGraphNode (ScaffoldGraph -> RezGraph,
                                               cid);
                        if  (contig -> info . Contig . numCI == 1)
                            {
                             ChunkInstanceT  * unitig;

                             assert (contig -> info . Contig . AEndCI
                                       == contig -> info . Contig . BEndCI);
                             unitig
                                 = GetGraphNode (ScaffoldGraph -> CIGraph,
                                                 contig -> info . Contig . AEndCI);
                             switch  (kind)
                               {
                                case  ROCKS :
                                  unitig -> flags . bits . isRock = TRUE;
                                  break;
                                case  STONES :
                                  unitig -> flags . bits . isStone = TRUE;
                                  break;
                                case  WALKS :
                                  unitig -> flags . bits . isWalk = TRUE;
                                  break;
                                default :
                                  fprintf (stderr,
                                           "ERROR:  Unexpected insert type = %d\n",
                                           (int) kind);
                               }
                            }

			// first we split the chunk
			splitCid = SplitUnresolvedContig(sgraph->RezGraph, cid, NULL, copyAllOverlaps);
	      
			// now we insert the surrogate into the scaffold
			// remember that it has a different number
			//  fprintf(stderr, "BEFORE Inserting %d \n",splitCid);
			//     DumpCIScaffold(stderr,sgraph, scaffold, FALSE);
                        contig = GetGraphNode (sgraph -> RezGraph, splitCid);
                        switch  (kind)
                          {
                           case  ROCKS :
                             contig -> flags . bits . isRock = TRUE;
                             break;
                           case  STONES :
                             contig -> flags . bits . isStone = TRUE;
                             break;
                           case  WALKS :
                             contig -> flags . bits . isWalk = TRUE;
                             break;
                           default :
                             fprintf (stderr, "ERROR:  Unexpected insert type = %d\n",
                                      (int) kind);
                          }
			InsertCIInScaffold(sgraph, splitCid, scaff_id,
							   fc->chunk[k].start, fc->chunk[k].end, TRUE, NO_CONTIGGING);
		  }

		  chunk = GetGraphNode(sgraph->RezGraph, cid);    
		  assert(chunk != NULL);
	  
		  //
		  // remember that we touched this scaffold and the
		  // number of insertions
		  //
		  scaffold_modified = TRUE;
#if  SIMULATED_DATA
		  if  (chunk -> flags . bits . cgbType != UU_CGBTYPE)
			fprintf (stderr,
					 "### Inserted invalid chunk #%d  type = %s\n",
					 chunk -> id,
					 CGB_Type_As_String (chunk -> flags . bits . cgbType));
#endif
		  inserted++;
	  


#         if DEBUG_UPDATE > 1
		  fprintf(stderr,
				  "\n  -=> After insertion of CI %d in scaffold %d, ends:%5.2f,%5.2f components: %d\n",
				  fc->chunk[k].chunk_id,
				  scaff_id,
				  fc->chunk[k].start.mean,
				  fc->chunk[k].end.mean,
				  IsScaffoldInternallyConnected(sgraph, scaffold));
#         if DEBUG_UPDATE > 2
		  DumpCIScaffold(stderr,sgraph, scaffold, FALSE);
#         endif
#         endif
		}
      }

      // Beware
      //      fprintf(stderr,"**** FORCE variances after filling gap %d  ****\n",j);
      //      Force_Increasing_Variances();

    }
    //
    // if the scaffold has been modified ...
    //
    if (scaffold_modified) {
#     if UPDATE_PARANOID_CHECKING > 0
      //
      // make sure that trusted edges are internal to the scaffold
      //
      InitCIScaffoldTIterator(sgraph, scaffold, TRUE,
							  FALSE, &CIs);
      while (chunk = NextCIScaffoldTIterator(&CIs))
		CheckTrustedEdges(sgraph, chunk->cid);
#     endif   

#     if DEBUG_UPDATE > 1
      fprintf (stderr,
			   "  -=> Recompute scaffold %d, number of components %d\n",
			   scaff_id,
			   IsScaffoldInternallyConnected(sgraph, scaffold));
#     endif
     
#     if UPDATE_PARANOID_CHECKING > 0
      //
      // check the damage
      //
      if (IsScaffoldInternallyConnected(sgraph, scaffold) != 1)
		fprintf(stderr, "* error: the scaffold has been disconnected\n");
#     endif   

#if  0   // Set to 0 for Granger's new way
      //
      // recompute positions (and cross the fingers :-)
      //
      // If we split we do not Recompute the Offsets
      if( ! noRecompute )
	  {
#if  0
if  (scaff_id == 209)
    {
     fprintf (stderr, "\n### BEFORE  RecomputeOffsetsInScaffold  for scaff %d\n",
              scaff_id);
     DumpCIScaffold(stderr,sgraph, scaffold, FALSE);
    }
#endif
		res = RecomputeOffsetsInScaffold(sgraph, scaffold, TRUE, TRUE, FALSE);
//		res = RecomputeOffsetsInScaffold(sgraph, scaffold, TRUE, TRUE, TRUE);
#if  0
if  (scaff_id == 209)
    {
     fprintf (stderr, "\n### AFTER  RecomputeOffsetsInScaffold  for scaff %d\n",
              scaff_id);
     DumpCIScaffold(stderr,sgraph, scaffold, FALSE);
    }
#endif
		switch (res) {
		  case RECOMPUTE_OK:
#       if DEBUG_UPDATE > 0
			fprintf(stderr,
					"-=> recompute: ok\n");
#       endif
			CheckCIScaffoldT(sgraph, scaffold);  // NEW
			break;
		  case RECOMPUTE_SINGULAR:
#       if DEBUG_UPDATE > 0
			fprintf(stderr,
					"-=> recompute failed: the matrix is singular\n");
#       endif
			fprintf(stderr, "* error failed: the matrix is singular for scaffold %d\n", scaff_id);
			break;
		  case RECOMPUTE_LAPACK:
#       if DEBUG_UPDATE > 0
			fprintf(stderr,
					"-=> recompute failed: lapack failed\n");
#       endif
			fprintf(stderr, "* error failed: lapack failed on scaffold %d\n", scaff_id);
			break;
		  case RECOMPUTE_NO_GAPS:
#       if DEBUG_UPDATE > 0
			fprintf(stderr,
					"-=> recompute failed: the scaffold has only one CI\n");
#       endif
			fprintf(stderr, "* warning: scaffold %d has only one CI\n", scaff_id);
			break;
		  case RECOMPUTE_FAILED_REORDER_NEEDED:
#       if DEBUG_UPDATE > 0
			fprintf(stderr,
					"-=> recompute failed: reordering needed in scaffold %d\n", scaff_id);
#       endif
			fprintf(stderr, "* error: reordering needed in scaffold %d\n", scaff_id);
			break;
		  case RECOMPUTE_NOT_ENOUGH_CLONES:
#       if DEBUG_UPDATE > 0
			fprintf(stderr,
					"-=> recompute failed: not enough clones in scaffold %d\n", scaff_id);
#       endif
			fprintf(stderr, "* error: not enough clones in scaffold %d\n", scaff_id);
			break;
		  default:
#       if DEBUG_UPDATE > 0
			fprintf(stderr,
					"-=> recompute failed: unknow reasons in scaffold %d\n", scaff_id);
#       endif
			fprintf(stderr, "* error: recompute failed for unknown reasons for scaffold %d\n", scaff_id);
			break;
		}
	  
	  }
#endif

      //
      // do we want to change the TRUSTED/UNTRUSTED status?
      //
      if (ChiSquare) {
		MarkInternalEdgeStatus(sgraph, scaffold,
							   PAIRWISECHI2THRESHOLD_CGW,
							   1000000.0, TRUE, TRUE, 0, FALSE); // use 100000000000.0 if using guides
		//
		// check the damage again
		//
		if (IsScaffoldInternallyConnected(sgraph, scaffold, ALL_TRUSTED_EDGES) != 1)
		  fprintf(stderr, "* <REZ> * error: the scaffold has been disconnected after the MarkInternalEdgeStatus()\n");
      }

#     if UPDATE_PARANOID_CHECKING > 0
      //
      // make sure that trusted edges are internal to the scaffold
      //
      InitCIScaffoldTIterator(sgraph, scaffold, TRUE,
							  FALSE, &CIs);
      while (chunk = NextCIScaffoldTIterator(&CIs))
		CheckTrustedEdges(sgraph, chunk->cid);
#     endif   
   
#     if DEBUG_UPDATE > 1
      fprintf(stderr,
			  "\n-=> After adjustement\n");
#     if DEBUG_UPDATE > 1
      DumpCIScaffold(stderr,sgraph, scaffold, FALSE);
#     endif
      fprintf(stderr,
			  "\n--- End processing scaffold %3d: (LSE %g, num %d, components %d) ---\n",
			  scaff_id,
			  scaffold->info.Scaffold.leastSquareError,
			  scaffold->info.Scaffold.numLeastSquareClones,
			  IsScaffoldInternallyConnected(sgraph, scaffold));
#     endif
    }
  }
  fprintf(stderr,"* Inserted %d elements ... rebuilding scaffolds\n", inserted);
  // ****** NEW ******
  // Now rebuild the scaffolds from scratch (if we do not split)

#if  0   // Set to 0 for Granger's new way
  if( ! noRecompute )
    RebuildScaffolds(ScaffoldGraph, FALSE);
#endif

# if DEBUG_UPDATE > 1
  fprintf(stderr,
		  "-=> UPDATE ended\n");
# endif

  free(table);
  return inserted;
}

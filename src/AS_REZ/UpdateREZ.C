
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

static const char *rcsid = "$Id: UpdateREZ.c,v 1.25 2012-11-15 03:44:36 brianwalenz Exp $";

/**********************************************************************

        Module:  UpdateREZ.c

   Description:  Contains functions that take a Scaffold_Fill_t and
                 update the Scaffold Graph

		 log info is sent to <inputFile>.cgwlog

    Programmer:  S. Lonardi (stelo@cs.purdue.edu)

       Written:  17 May 99

 **********************************************************************/

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "AS_global.H"
#include "AS_UTL_Var.H"
#include "UnionFind_AS.H"

#include "AS_CGW_dataTypes.H"
#include "Globals_CGW.H"
#include "ScaffoldGraph_CGW.H"
#include "ChiSquareTest_CGW.H"

#include "DataTypesREZ.H"
#include "UtilsREZ.H"
#include "CommonREZ.H"
#include "GapWalkerREZ.H"
#include "UpdateREZ.H"
#include "ConsistencyChecksREZ.H"

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
  int32
    next;
  int
    marked = 0;
  ChunkInstanceT
    * next_chunk;

  GraphEdgeIterator  edges(sgraph->ContigGraph, cid, ALL_END, ALL_EDGES);
  CIEdgeT           *edge;

  while((edge = edges.nextMerged()) != NULL){
    assert(edge != NULL);
    //
    // get the other end
    //
    if (cid == edge->idA)
      next = edge->idB;
    else
      next = edge->idA;

    next_chunk = GetGraphNode(sgraph->ContigGraph, next);
    assert(next_chunk != NULL);

    if ((IsUnique(next_chunk) && (next_chunk->scaffoldID == sid)) ||
        ((table[next] != NULL) && (table[next]->scaff_id == sid) && (table[next]->keep))) {
      //
      // mark the edge only if the distance is consistent
      //
      if (IsUnique(next_chunk)) {
        Gap_Chunk_t
          unique;
        unique.chunk_id = next_chunk->id;
        unique.scaff_id = next_chunk->scaffoldID;
        unique.start = next_chunk->offsetAEnd;
        unique.end = next_chunk->offsetBEnd;

        if (Is_Edge_Consistent(edge, table[cid], &unique)) {
          if (Is_Edge_Orientation_Consistent(edge, table[cid], &unique)) {
            SetEdgeStatus(sgraph->ContigGraph, edge, TENTATIVE_TRUSTED_EDGE_STATUS);
            marked += edge->edgesContributing;
          } else {
            SetEdgeStatus(sgraph->ContigGraph, edge, TENTATIVE_UNTRUSTED_EDGE_STATUS);
          }
        } else {
          SetEdgeStatus(sgraph->ContigGraph, edge, TENTATIVE_UNTRUSTED_EDGE_STATUS);
        }
      } else {
        if (Is_Edge_Consistent(edge, table[cid], table[next])) {
          if (Is_Edge_Orientation_Consistent(edge, table[cid], table[next])) {
            SetEdgeStatus(sgraph->ContigGraph, edge, TENTATIVE_TRUSTED_EDGE_STATUS);
            marked += edge->edgesContributing;
          } else {
            SetEdgeStatus(sgraph->ContigGraph, edge, TENTATIVE_UNTRUSTED_EDGE_STATUS);
          }
        } else {
          SetEdgeStatus(sgraph->ContigGraph, edge, TENTATIVE_UNTRUSTED_EDGE_STATUS);
        }
      }
    } else {
      //
      // mark the edge as UNTRUSTED
      //
      SetEdgeStatus(sgraph->ContigGraph, edge, TENTATIVE_UNTRUSTED_EDGE_STATUS);
    }
#if DEBUG_UPDATE > 1
    fprintf(stderr,"  -=> edge %d->%d (edge weight %d, new status %d, mean %f, var %f) (other end sid %d, %s)",
            cid,
            next,
            edge->edgesContributing,
            GetEdgeStatus(edge),
            edge->distance.mean,
            edge->distance.variance,
            next_chunk->scaffoldID,
            table[next] ? "in the set" : "not in the set");
#endif
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
    num_chunks = GetNumGraphNodes(sgraph->ContigGraph),
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

  Gap_Chunk_t
    * * table = (Gap_Chunk_t * *)safe_calloc(num_chunks, sizeof(Gap_Chunk_t *));
  ChunkInstanceT
    * chunk;

  //  Compute how much work we need to do

  uint64  numScaffolds = 0, curScaffolds = 0;
  uint64  numGaps      = 0, curGaps      = 0;
  uint64  numChunks    = 0, curChunks    = 0;
  uint64                    curInserted  = 0;

  for (scaff_id = 0; scaff_id < num_scaffs;  scaff_id++) {
    if( scaffId != -1 )
      if( scaffId != scaff_id )
        continue;

    scaffold = GetCIScaffoldT(sgraph->CIScaffolds, scaff_id);
    assert(scaffold != NULL);
    if ((isDeadCIScaffoldT(scaffold)) || (scaffold->type != REAL_SCAFFOLD))
      continue;

    numScaffolds++;

    numGaps += fill_chunks[scaff_id].num_gaps;

    for (j = 0; j < fill_chunks[scaff_id].num_gaps;  j++)
      numChunks += fill_chunks[scaff_id].gap[j].num_chunks;
  }

  fprintf(stderr, "Update_Scaffold_Graph()-- place "F_U64" chunks into "F_U64" gaps in "F_U64" scaffolds.\n",
          numChunks, numGaps, numScaffolds);

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

    for (j = 0; j < fill_chunks[scaff_id].num_gaps;  j++) {
      fc = &(fill_chunks[scaff_id].gap[j]);
#if DEBUG_UPDATE > 1
      fprintf(stderr,"\n -=> Gap %3d:  (%5.2f, %5.2f)-(%5.2f, %5.2f)\n",
              j,
              fc->start.mean,
              sqrt(fc->start.variance),
              fc->end.mean,
              sqrt(fc->end.variance));
#endif
      for  (k = 0;  k < fc->num_chunks;  k++) {
        // SK - 10/28/08 - if we have just one chunk keep it
        // otherwise, we pick the first one that is kept for that specific cid
        // this was added for closure reads where we cannot tell their orientations (vs regular rocks where we can use mates to figure it out)
        // therefore, we add both orientations of a rock and let the overlaps decide which we keep
        if (table[fc->chunk[k].chunk_id] == NULL || (fc->chunk[k].isClosure && fc->chunk[k].keep) || fc->chunk[k].copy_letter != ' ') {
          if (table[fc->chunk[k].chunk_id] != NULL && table[fc->chunk[k].chunk_id]->keep == 1 && fc->chunk[k].copy_letter == ' ') {
            fprintf(stderr, "ERROR: chunk %d (%f, %f) copy_letter=%c keep=%d closure=%d path=%d has keep flag on but there is already a chunk %d (%f, %f) copy_letter=%c keep=%d closure=%d scf=%d curr_scf=%d gap=%d curr_gap=%d\n", 
                    fc->chunk[k].chunk_id, fc->chunk[k].start.mean, fc->chunk[k].end.mean,
                    fc->chunk[k].copy_letter, fc->chunk[k].keep, fc->chunk[k].isClosure, fc->chunk[k].path_confirmed,
                    table[fc->chunk[k].chunk_id]->chunk_id, table[fc->chunk[k].chunk_id]->start.mean, table[fc->chunk[k].chunk_id]->end.mean,
                    table[fc->chunk[k].chunk_id]->copy_letter, table[fc->chunk[k].chunk_id]->keep, table[fc->chunk[k].chunk_id]->isClosure, 
                    table[fc->chunk[k].chunk_id]->scaff_id, fc->chunk[k].scaff_id,
                    table[fc->chunk[k].chunk_id]->gap, fc->chunk[k].gap);
            continue;
          }
          table[fc->chunk[k].chunk_id] = &(fc->chunk[k]);
        }
#if DEBUG_UPDATE > 1
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
#endif
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

    curScaffolds++;

    //
    // now loop over all the gaps in "scaffold" (scaff_id)
    //
    scaffold_modified = FALSE;
    for (j = 0; j < fill_chunks[scaff_id].num_gaps;  j++) {
      curGaps++;

      //
      // a convenient pointer
      //
      fc = &(fill_chunks[scaff_id].gap[j]);

#if DEBUG_UPDATE > 1
      fprintf(stderr,"\n -=> Gap %3d:  (%5.2f, %5.2f)-(%5.2f, %5.2f)\n",
              j,
              fc->start.mean,
              sqrt(fc->start.variance),
              fc->end.mean,
              sqrt(fc->end.variance));
#endif


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

          curChunks++;
          curInserted++;

          if  (! fc -> chunk [k] . keep)
            continue;

          chunkScaffID
            = GetGraphNode (ScaffoldGraph -> ContigGraph,
                            fc -> chunk [k] . chunk_id) -> scaffoldID;
          chunkScaffold = GetCIScaffoldT (sgraph -> CIScaffolds, chunkScaffID);
          // scaffold the current chunk is in. Has ID chunkScaffID

          /* scaffold ID must be different */
          //fprintf(stderr,"chunk %d, scaff_id = %d chunkScaffID = %d \n",fc->chunk[k].chunk_id,scaff_id,chunkScaffID);

          if( (chunkScaffID != scaff_id) && (chunkScaffID != NULLINDEX) )
            {
              //
              // iterate through the chunks that are in the same scaffold "chunkScaffold"
              // as the chunk slated for insertion in the current gap
              //
              InitCIScaffoldTIterator(sgraph, chunkScaffold, TRUE,
                                      FALSE,&internalChunks);
              /* for each chunk in the scaffold we look up whether he is in the fc */
              while((internalChunk = NextCIScaffoldTIterator(&internalChunks)) != NULL)
                {
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
                  ChunkInstanceT *ci = GetGraphNode(ScaffoldGraph->ContigGraph, fc->chunk[l].chunk_id);
                  if( ci->scaffoldID == chunkScaffID )
                    {
                      if( keep )
                        {
                          RemoveCIFromScaffold(ScaffoldGraph,chunkScaffold,ci, TRUE);
                          fc->chunk[l].keep = TRUE;
                          fc->chunk[l].split = FALSE;

                        }
                      else
                        {
                          fc->chunk[l].keep = FALSE;
                        }
                    }
                }
              if( keep )
                {
                  chunkScaffold->flags.bits.isDead = TRUE;
                }
            }
        }

      for  (k = 0;  k < fc->num_chunks;  k++) {
        // here we have to filter out the scaffold ends
        // that means any chunk that is in the scaffold in which the gap is
        //fprintf(stderr,"Chunk Scaff ID = %d, scaff_id= %d \n",GetGraphNode(ScaffoldGraph->ContigGraph, fc->chunk[k].chunk_id)->scaffoldID, scaff_id);

        if( GetGraphNode(ScaffoldGraph->ContigGraph, fc->chunk[k].chunk_id)->scaffoldID == scaff_id )
          fc->chunk[k].keep = FALSE;

#if DEBUG_UPDATE > 1
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
#endif


        //
        // if the bit keep is set AND the chunk is not the from or destination chunk
        // then we try to insert the CI
        //
        //	if( GetGraphNode(ScaffoldGraph->ContigGraph, fc->chunk[k].chunk_id)->scaffoldID != NULLINDEX )
        //	  fprintf(stderr,"*** Skipping chunk already in Scaffold\n");
        //	else
        if (fc->chunk[k].keep )
          {
            int cid = fc->chunk[k].chunk_id;

            // mark the edges of this chunk: if no edges are trusted,
            // avoid inserting this chunk
            //
            if( markTrusted )
              {
                trusted_edges = Trust_Edges(sgraph, fill_chunks, table, scaff_id, cid);

                if (trusted_edges < 1) {
                  if (fc->chunk[k].isClosure) {
                    fprintf(stderr,
                            "Update_Scaffold_Graph()-- warning: closure CI %d with a total weight %d of trusted edges not being inserted - should it be?\n",
                            cid,
                            trusted_edges);
                    continue;
                  } else {
                    fprintf(stderr,
                            "Update_Scaffold_Graph()-- warning : attempted to insert CI %d with a total weight %d of trusted edges - skipped\n",
                            cid,
                            trusted_edges);
                    continue;
                  }
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
                contig = GetGraphNode (sgraph -> ContigGraph, cid);

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
                                   "Update_Scaffold_Graph()-- ERROR:  Unexpected insert type = %d\n",
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
                      fprintf (stderr, "Update_Scaffold_Graph()-- ERROR:  Unexpected insert type = %d\n",
                               (int) kind);
                  }
                InsertCIInScaffold(sgraph, cid, scaff_id,
                                   fc->chunk[k].start, fc->chunk[k].end, TRUE, FALSE);
              }
            else
              {
                ContigT  * contig;
                int splitCid;

                contig = GetGraphNode (ScaffoldGraph -> ContigGraph,
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
                          fprintf (stderr, "Update_Scaffold_Graph()-- ERROR:  Unexpected insert type = %d\n",
                                   (int) kind);
                      }
                  }

                // first we split the chunk
                splitCid = SplitUnresolvedContig(sgraph->ContigGraph, cid, NULL, copyAllOverlaps);

                // now we insert the surrogate into the scaffold
                // remember that it has a different number
                //  fprintf(stderr, "BEFORE Inserting %d \n",splitCid);
                //     DumpCIScaffold(stderr,sgraph, scaffold, FALSE);
                contig = GetGraphNode (sgraph -> ContigGraph, splitCid);
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
                      fprintf (stderr, "Update_Scaffold_Graph()-- ERROR:  Unexpected insert type = %d\n",
                               (int) kind);
                  }
                InsertCIInScaffold(sgraph, splitCid, scaff_id, fc->chunk[k].start, fc->chunk[k].end, TRUE, FALSE);
              }

            chunk = GetGraphNode(sgraph->ContigGraph, cid);
            assert(chunk != NULL);        

            //
            // remember that we touched this scaffold and the
            // number of insertions
            //
            scaffold_modified = TRUE;
            inserted++;
          }
      }

      // Beware
      //      fprintf(stderr,"**** FORCE variances after filling gap %d  ****\n",j);
      //      Force_Increasing_Variances();

    }

    if ((scaffold_modified) &&
        (curInserted > 100000)) {
      curInserted = 0;

      fprintf(stderr, "Update_Scaffold_Graph()-- placed "F_U64"/"F_U64" chunks into "F_U64"/"F_U64" gaps in "F_U64"/"F_U64" scaffolds.\n",
              curChunks,    numChunks,
              curGaps,      numGaps,
              curScaffolds, numScaffolds);

      //  This wasn't as useful as expected.
      //ScaffoldGraph->tigStore->flushCache();
    }

    if (scaffold_modified) {
      if (ChiSquare) {
        MarkInternalEdgeStatus(sgraph, scaffold, 0, FALSE);        //  OPERATES ON RAW

        //
        // check the damage again
        //
        //  true = merged; true = trusted
        if (IsScaffoldInternallyConnected(sgraph, scaffold, true, true) != 1)
          fprintf(stderr, "Update_Scaffold_Graph()-- error: the scaffold has been disconnected after the MarkInternalEdgeStatus()\n");
      }
    }
  }

  fprintf(stderr,"Update_Scaffold_Graph()-- Inserted %d elements ... rebuilding scaffolds\n", inserted);

  safe_free(table);
  return inserted;
}

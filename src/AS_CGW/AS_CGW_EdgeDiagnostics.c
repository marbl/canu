
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
static char CM_ID[] = "$Id: AS_CGW_EdgeDiagnostics.c,v 1.2 2004-09-23 20:25:02 mcschatz Exp $";


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "AS_CGW_EdgeDiagnostics.h"

/*
  1st set of functions:
    Check that contig orientation is preserved before & after
    a. create structure to hold IDs & orientations
    b. delete structure
    c. populate before
    d. compare after
 */


#define INITIAL_NUM_ORIENTS 1000

void GetFragment5pAndOrientationInChunk(ScaffoldGraphT * graph,
                                        CIFragT * frag,
                                        LengthT * offset5p,
                                        FragOrient * orient,
                                        ChunkInstanceT * chunk)
{
  LengthT * localFragOffset5p = NULL;
  LengthT * localFragOffset3p = NULL;
  FragOrient localOrient = X_X;
  // chunk could be scaffold, contig, or chunk
  /*
    contigID, contigOffset5p & 3p are in CIFragT
    offsetAEnd & offsetBEnd & bpLength are in NodeCGW_T (contig)

    Four possible configurations
    Scaffold:    ------------------------------------->
    Contig:             ------------------->
    Unitig:
    Fragment:                -------->
                 mean = contigAEnd + frag5p
                 var = contigAEndVar + frag5pVar
                 orientation = A_B

    Scaffold:    ------------------------------------->
    Contig:             ------------------->
    Fragment:                <--------
                 mean = contigAEnd + frag5p
                 var = contigAEndVar + frag5pVar
                 orientation = B_A

    Scaffold:    ------------------------------------->
    Contig:             <-------------------
    Fragment:                -------->
                 mean = contigAEnd - frag5p
                 var = contigAEndVar - frag5pVar
                 orientation = A_B

    Scaffold:    ------------------------------------->
    Contig:             <-------------------
    Fragment:                <--------
                 mean = contigAEnd - frag5p
                 var = contigAEndVar - frag5pVar
                 orientation = B_A
    
   */
  switch(chunk->type)
  {
    case DISCRIMINATORUNIQUECHUNK_CGW:
    case UNRESOLVEDCHUNK_CGW:
    case UNIQUECHUNK_CGW:
    case RESOLVEDREPEATCHUNK_CGW:
    {
      offset5p->mean = 0.0;
      offset5p->variance = 0.0;
      localFragOffset5p = &(frag->offset5p);
      localFragOffset3p = &(frag->offset3p);
      localOrient = A_B;
    }
    break;
    case CONTIG_CGW:
    case UNIQUECONTIG_CGW:
    case RESOLVEDCONTIG_CGW:
    case UNRESOLVEDCONTIG_CGW:
    {
      offset5p->mean = 0.0;
      offset5p->variance = 0.0;
      localFragOffset5p = &(frag->contigOffset5p);
      localFragOffset3p = &(frag->contigOffset3p);
      localOrient = A_B;
    }
    break;
    case REAL_SCAFFOLD:
    case OUTPUT_SCAFFOLD:
    case SCRATCH_SCAFFOLD:
    {
      ChunkInstanceT * contig = GetGraphNode(graph->RezGraph, frag->contigID);
      offset5p->mean = contig->offsetAEnd.mean;
      // use variances to end of contig, not to end of scaffold!
      offset5p->variance = 0.0;
      localFragOffset5p = &(frag->contigOffset5p);
      localFragOffset3p = &(frag->contigOffset3p);
      localOrient =
        (contig->offsetAEnd.mean < contig->offsetBEnd.mean) ? A_B : B_A;
    }
    break;
  }

  offset5p->variance = localFragOffset5p->variance;
  if(localOrient == A_B)
  {
    offset5p->mean += localFragOffset5p->mean;
    *orient =
      (localFragOffset5p->mean < localFragOffset3p->mean) ? A_B : B_A;
  }
  else
  {
    offset5p->mean -= localFragOffset5p->mean;
    *orient =
      (localFragOffset5p->mean < localFragOffset3p->mean) ? B_A : A_B;
  }
}


void ComputeFrag5pToChunkEndFromOffset(LengthT * fragOffset5p,
                                       ChunkOrientType whichEnd,
                                       LengthT * distFromEnd,
                                       ChunkInstanceT * chunk)
{
  if(whichEnd == A_END)
  {
    distFromEnd->mean = fragOffset5p->mean;
    distFromEnd->variance = fragOffset5p->variance;
  }
  else
  {
    distFromEnd->mean = chunk->bpLength.mean - fragOffset5p->mean;
    distFromEnd->variance =
      chunk->bpLength.variance - fragOffset5p->variance;
  }
  distFromEnd->variance =
    (distFromEnd->variance < 1.0) ? 1.0 : distFromEnd->variance;
}


void ComputeFrag5pToChunkEnd(ScaffoldGraphT * graph,
                             CIFragT * frag,
                             ChunkOrientType whichEnd,
                             LengthT * distFromEnd,
                             ChunkInstanceT * chunk)
{
  LengthT offset5p;
  FragOrient fragOrient;

  GetFragment5pAndOrientationInChunk(graph, frag, &offset5p, &fragOrient,
                                     chunk);

  ComputeFrag5pToChunkEndFromOffset(&offset5p,
                                    whichEnd,
                                    distFromEnd,
                                    chunk);
}


void ComputeFragToChunkEndForEdge(ScaffoldGraphT * graph,
                                  CIFragT * frag,
                                  ChunkOrientType * endFromFrag,
                                  LengthT * distFromEnd,
                                  ChunkInstanceT * chunk)
{
  LengthT offset5p;
  FragOrient fragOrient;

  GetFragment5pAndOrientationInChunk(graph, frag, &offset5p, &fragOrient,
                                     chunk);
  /*
    Which end of chunk is relevant depends on frag orientation
    and innie or outtie mate
    1. Chunk:     ------------------->
       Frag:            ----->          (A_B)
       i) Innie:        |------------|  = B_END
      ii) Outtie: |-----|               = A_END
    2. Chunk:     ------------------->
       Frag:            <-----          (B_A)
       i) Innie:  |----------|          = A_END
      ii) Outtie:            |-------|  = B_END
   */
  if(frag->flags.bits.innieMate)
    *endFromFrag = (fragOrient == A_B) ? B_END : A_END;
  else
    *endFromFrag = (fragOrient == A_B) ? A_END : B_END;

  ComputeFrag5pToChunkEndFromOffset(&offset5p,
                                    *endFromFrag,
                                    distFromEnd,
                                    chunk);
}


void PopulateChunkEdgeBasics(ScaffoldGraphT * graph,
                             CIFragT * fragA,
                             ChunkInstanceT * chunkA,
                             CIFragT * fragB,
                             ChunkInstanceT * chunkB,
                             DistT * dist,
                             EdgeCGW_T * edge)
{
  ChunkOrientType chunkEndFromFragA;
  ChunkOrientType chunkEndFromFragB;
  LengthT distFromAChunkEnd;
  LengthT distFromBChunkEnd;
  LengthT distBetweenChunks;
  cds_float32 distVariance = dist->stddev * dist->stddev;
  
  // determine distance from 5p of fragment to relevant end of CI
  ComputeFragToChunkEndForEdge(graph, fragA,
                               &chunkEndFromFragA,
                               &distFromAChunkEnd,
                               chunkA);
  ComputeFragToChunkEndForEdge(graph, fragB,
                               &chunkEndFromFragB,
                               &distFromBChunkEnd,
                               chunkB);
  
  distBetweenChunks.mean = dist->mean -
    (distFromAChunkEnd.mean + distFromBChunkEnd.mean);

  // need to be careful with containment edges
  if(distBetweenChunks.mean < 0.0)
  {
    /*
      There are two issues:
      1. must return the smallest of two possible edges
      2. must return a positive variance
     */
    /*
      May need to re-orient containment edge
      Consider the following example:
                        lowIID
                      <--------
                highIID
        <------------------------

      Since edges are stored with lowIID = idA, highIID = idB,
      This orientation had better be AS_NORMAL so that entities will be
      treated as
      
        lowIID            highIID
       -------->  ------------------------>

      instead of
                  
             highIID               lowIID
       ------------------------>  -------->


      What are the cases? (l = lowIID, h = highIID)
       They all simplify to requiring edge to have the
         smaller (in magnitude) of the two possible
         edge distance
         
      AS_OUTTIE - h must be off A end of l:
                    l              l
                   --->          <---
               h                       h
        <---------------        --------------->

                    h              h
                   --->          <---
               l                       l
        <---------------        --------------->

      AS_INNIE - h must be off B end of l:
          l                                  l
         --->                              <---
               h                       h
        <---------------        --------------->

          h                                  h
         --->                              <---
               l                       l
        <---------------        --------------->

      AS_NORMAL - h entity must be off B end of l:
                     l            l
                   <---          --->
               h                       h
        <---------------        --------------->

           h                                h
         <---                              --->
               l                       l
        <---------------        --------------->

      AS_ANTI - h entity must be off A end of l:
                     h            h
                   <---          --->
               l                       l
        <---------------        --------------->

           l                                l
         <---                              --->
               h                       h
        <---------------        --------------->
       
    */

    /*
      Don't bother to compute if magnitude of the edge is
      smaller of two possible edges
    */
    if(-distBetweenChunks.mean >
       (chunkA->bpLength.mean + chunkB->bpLength.mean) / 2)
    {
      chunkEndFromFragA =
        (chunkEndFromFragA == A_END) ? B_END : A_END;
      ComputeFrag5pToChunkEnd(graph, fragA,
                              chunkEndFromFragA,
                              &distFromAChunkEnd,
                              chunkA);
      
      chunkEndFromFragB =
        (chunkEndFromFragB == A_END) ? B_END : A_END;
      ComputeFrag5pToChunkEnd(graph, fragB,
                              chunkEndFromFragB,
                              &distFromBChunkEnd,
                              chunkB);
      
      distBetweenChunks.mean = -(dist->mean +
                                 distFromAChunkEnd.mean +
                                 distFromBChunkEnd.mean);
    }
  }

  // variances ALWAYS add
  distBetweenChunks.variance =
    distVariance + distFromAChunkEnd.variance + distFromBChunkEnd.variance;
  
  if(chunkEndFromFragA == A_END)
    edge->orient = (chunkEndFromFragB == A_END) ? BA_AB : BA_BA;
  else
    edge->orient = (chunkEndFromFragB == A_END) ? AB_AB : AB_BA;

  distBetweenChunks.variance = max(1.0, distBetweenChunks.variance);
  edge->distance = distBetweenChunks;
  edge->flags.bits.hasContainmentOverlap = (edge->distance.mean < 0.0) ? TRUE : FALSE;
}



void DeleteOrientProcessor(OrientProcessor * op)
{
  if(op)
  {
    if(op->array)
      DeleteVA_OrientHolder(op->array);
    if(op->ht)
      DeleteHashTable_AS(op->ht);
    free(op);
  }
}


int OrientHashFunction(const void * item, int length)
{
  return Hash_AS((uint8 *) item, length, 37);
}


int OrientCompareFunction(const void * item1, const void * item2)
{
  const OrientHolder * oh1 = (const OrientHolder *) item1;
  const OrientHolder * oh2 = (const OrientHolder *) item2;

  return(oh2->keyID - oh1->keyID);
}


OrientProcessor * CreateOrientProcessor(void)
{
  OrientProcessor * op = safe_calloc(sizeof(OrientProcessor), 1);
  if(op)
  {
    op->array = CreateVA_OrientHolder(INITIAL_NUM_ORIENTS);
    op->ht = CreateHashTable_AS(INITIAL_NUM_ORIENTS,
                                OrientHashFunction,
                                OrientCompareFunction);
  }
  if(!op || !op->array || !op->ht)
  {
    fprintf(stderr, "Failed to create orientation processor!\n");
    free(op);
    return NULL;
  }

  return op;
}


void ResetOrientProcessor(OrientProcessor * op)
{
  if(op)
  {
    if(op->array)
      ResetVA_OrientHolder(op->array);
    if(op->ht)
      ResetHashTable_AS(op->ht);
  }
}


void DestroyContigOrientChecker(ContigOrientChecker * coc)
{
  if(coc)
  {
    DeleteOrientProcessor(coc->contigs);
    DeleteOrientProcessor(coc->scaffolds);
    free(coc);
  }
}


ContigOrientChecker * CreateContigOrientChecker(void)
{
  ContigOrientChecker * coc = safe_calloc(sizeof(ContigOrientChecker), 1);

  if(coc == NULL)
  {
    fprintf(stderr, "Failed to allocate a ContigOrientChecker!\n");
    return NULL;
  }
  coc->contigs = CreateOrientProcessor();
  coc->scaffolds = CreateOrientProcessor();
  if(coc->contigs == NULL || coc->scaffolds == NULL)
  {
    fprintf(stderr, "Failed to allocate ContigOrientProcessor(s)!\n");
    DestroyContigOrientChecker(coc);
    return NULL;
  }
  return coc;
}


void ResetContigOrientChecker(ContigOrientChecker * coc)
{
  if(coc)
  {
    ResetOrientProcessor(coc->contigs);
    ResetOrientProcessor(coc->scaffolds);
  }
}


int AddScaffoldToContigOrientChecker(ScaffoldGraphT * graph,
                                     CIScaffoldT * scaffold,
                                     ContigOrientChecker * coc)
{
  OrientHolder oh;
  CIScaffoldTIterator ciTemp;

  // add scaffold to array - get orientation later based on first contig
  oh.keyID = scaffold->id;
  oh.orient = X_X;
  AppendVA_OrientHolder(coc->scaffolds->array, &oh);

  // iterate over scaffold & add contigs with orientations
  oh.secondID = scaffold->id;
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &ciTemp);
  while(NextCIScaffoldTIterator(&ciTemp))
  {
    ContigT * contig;

    // get the left contig
    if((contig = GetGraphNode(graph->RezGraph, ciTemp.curr)) == NULL)
    {
      fprintf(stderr, "Failed to get contig " F_CID " from graph!\n", ciTemp.curr);
      return 1;
    }

    // append it to the list, with its orientation
    oh.keyID = contig->id;
    oh.orient = (contig->offsetAEnd.mean < contig->offsetBEnd.mean) ?
      A_B : B_A;
    AppendVA_OrientHolder(coc->contigs->array, &oh);
  }
  return 0;
}


int PopulateOrientHTFromArray(HashTable_AS * ht, VA_TYPE(OrientHolder) * array)
{
  int i;
  
  ResetHashTable_AS(ht);
  for(i = 0; i < GetNumVA_OrientHolder(array); i++)
  {
    OrientHolder * oh = GetVA_OrientHolder(array, i);
    InsertInHashTable_AS(ht, &(oh->keyID), sizeof(oh->keyID), oh);
  }
  return 0;
}


int CompareNewOrientationsForScaffold(ScaffoldGraphT * graph,
                                      CIScaffoldT * scaffold,
                                      ContigOrientChecker * coc)
{
  int status = 0;
  CIScaffoldTIterator ciTemp;
  CDS_COORD_t lastBegin = -1;
  CDS_COORD_t lastEnd = lastBegin;
  CDS_CID_t lastContigID = NULLINDEX;
  
  // loop over contigs in new scaffold
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &ciTemp);
  while(NextCIScaffoldTIterator(&ciTemp))
  {
    ContigT * contig;
    OrientHolder * contigOH;
    OrientHolder * scaffoldOH;

    // get the left contig
    if((contig = GetGraphNode(graph->RezGraph, ciTemp.curr)) == NULL)
    {
      fprintf(stderr, "Failed to get contig " F_CID " from graph!\n", ciTemp.curr);
      return -1;
    }

    // check strange begin/end coordinates
    if(lastBegin != -1 && lastEnd != -1)
    {
      if((CDS_COORD_t) (contig->offsetAEnd.mean + 0.5) == lastBegin ||
         (CDS_COORD_t) (contig->offsetBEnd.mean + 0.5) == lastEnd ||
         (CDS_COORD_t) (contig->offsetBEnd.mean + 0.5) == lastBegin ||
         (CDS_COORD_t) (contig->offsetAEnd.mean + 0.5) == lastEnd)
      {
        fprintf(stderr, "!!!! Contig " F_CID " (%d, %d) has same coordinate as"
                " prior contig " F_CID " (" F_COORD "," F_COORD ")\n",
                contig->id,
                (int) (contig->offsetAEnd.mean + 0.5),
                (int) (contig->offsetBEnd.mean + 0.5),
                lastContigID,
                lastBegin,
                lastEnd);
        fprintf(stderr, "Scaffold length is %d\n",
                (int) (scaffold->bpLength.mean + 0.5));
      }
    }

    lastBegin = (CDS_COORD_t) (contig->offsetAEnd.mean + 0.5);
    lastEnd = (CDS_COORD_t) (contig->offsetBEnd.mean + 0.5);
    lastContigID = contig->id;
    
    // see which scaffold this contig was in
    if((contigOH = LookupInHashTable_AS(coc->contigs->ht,
                                        &(contig->id),
                                        sizeof(contig->id))) == NULL)
    {
      // don't want to barf if a contig was inserted but not in our set
      continue;
    }

    // see if any contigs from this contig's scaffold have been seen yet
    if((scaffoldOH = LookupInHashTable_AS(coc->scaffolds->ht,
                                          &(contigOH->secondID),
                                          sizeof(contigOH->secondID))) == NULL)
    {
      fprintf(stderr, "Failed to lookup scaffold " F_CID " in hashtable\n",
              contigOH->secondID);
      return -1;
    }
    else
    {
      // Set the scaffold orientation based on this contig, if appropriate
      switch(scaffoldOH->orient)
      {
        case X_X:
          // no contigs from this scaffold have been seen yet
          if((contigOH->orient == A_B &&
              contig->offsetAEnd.mean < contig->offsetBEnd.mean) ||
             (contigOH->orient == B_A &&
              contig->offsetAEnd.mean > contig->offsetBEnd.mean))
            scaffoldOH->orient = A_B;
          else
            scaffoldOH->orient = B_A;
          break;
        case A_B:
          // here, contig should have same orientation as before merging
          if((contigOH->orient == A_B &&
              contig->offsetAEnd.mean > contig->offsetBEnd.mean) ||
             (contigOH->orient == B_A &&
              contig->offsetAEnd.mean < contig->offsetBEnd.mean))
          {
            fprintf(stderr, "!!!!! Contig " F_CID " from scaffold " F_CID " has been "
                    "reversed in scaffold " F_CID "\n",
                    contig->id, contigOH->secondID, contig->scaffoldID);
            status = 1;
          }
          break;
        case B_A:
          // here, contig should be reversed
          if((contigOH->orient == A_B &&
              contig->offsetAEnd.mean < contig->offsetBEnd.mean) ||
             (contigOH->orient == B_A &&
              contig->offsetAEnd.mean > contig->offsetBEnd.mean))
          {
            fprintf(stderr, "!!!!! Contig " F_CID " from scaffold " F_CID " has been "
                    "reversed in scaffold " F_CID "\n",
                    contig->id, contigOH->secondID, contig->scaffoldID);
            status = 1;
          }
          break;
      }
    }
  }
  return status;
}


// returns 0 if all is okay
int CompareNewOrientations(ScaffoldGraphT * graph,
                           CIScaffoldT * scaffold,
                           ContigOrientChecker * coc)
{
  // populate hashtables - one for contigs, one for scaffolds
  PopulateOrientHTFromArray(coc->contigs->ht, coc->contigs->array);
  PopulateOrientHTFromArray(coc->scaffolds->ht, coc->scaffolds->array);

  return(CompareNewOrientationsForScaffold(graph, scaffold, coc));
}


int AddAllScaffoldsToContigOrientChecker(ScaffoldGraphT * graph,
                                         ContigOrientChecker * coc)
{
  CDS_CID_t sID;

  // loop over all scaffolds
  for(sID = 0; sID < GetNumGraphNodes(graph->ScaffoldGraph); sID++)
  {
    // add to ContigOrientChecker
    if(AddScaffoldToContigOrientChecker(graph,
                                        GetGraphNode(graph->ScaffoldGraph,
                                                     sID),
                                        coc))
    {
      fprintf(stderr, "Failed to add scaffold " F_CID " to ContigOrientChecker!\n",
              sID);
      return 1;
    }
  }
  return 0;
}



int CheckAllContigOrientationsInAllScaffolds(ScaffoldGraphT * graph,
                                             ContigOrientChecker * coc,
                                             int populateHashTables)
{
  CDS_CID_t sID;
  int status = 0;

  // populate hashtables - one for contigs, one for scaffolds
  if(populateHashTables)
  {
    PopulateOrientHTFromArray(coc->contigs->ht, coc->contigs->array);
    PopulateOrientHTFromArray(coc->scaffolds->ht, coc->scaffolds->array);
  }

  // loop over all scaffolds
  for(sID = 0; sID < GetNumGraphNodes(graph->ScaffoldGraph); sID++)
  {
    int thisStatus;
    
    // compare to ContigOrientChecker
    thisStatus =
      CompareNewOrientationsForScaffold(graph,
                                        GetGraphNode(graph->ScaffoldGraph,
                                                     sID),
                                        coc);
    if(thisStatus != 0)
      fprintf(stderr, "Contig re-orientation problem in scaffold " F_CID "!\n", sID);

    status |= thisStatus;
  }
  
  return status;
}

#define DIFFERENT_EDGE_FACTOR    0.1f

int CompareEdgeFloats(float32 edgeFloat,
                      float32 fragsFloat,
                      char * message)
{
  float32 delta = fabs(edgeFloat - fragsFloat);
  if(delta > 1.f &&
     (delta > fabs(edgeFloat) * DIFFERENT_EDGE_FACTOR ||
     delta > fabs(fragsFloat) * DIFFERENT_EDGE_FACTOR))
  {
    fprintf(stderr, "!!!!!!! Edge %s not supported by fragment matepairs: %.2f (edge) vs. %.2f (frags)\n",
            message, edgeFloat, fragsFloat);
    return 1;
  }
  return 0;
}


int CompareEdgeMeans(EdgeCGW_T * edge1,
                     EdgeCGW_T * edge2)
{
  return(CompareEdgeFloats(edge1->distance.mean, edge2->distance.mean, "means"));
}


int CompareEdgeVariances(EdgeCGW_T * edge1,
                         EdgeCGW_T * edge2)
{
  return(CompareEdgeFloats(edge1->distance.variance, edge2->distance.variance, "variances"));
}

void PrintOrientation(FILE * fp, ChunkOrientationType orient)
{
  switch(orient)
  {
    case AB_AB:
      fprintf(fp, "AB_AB");
      break;
    case AB_BA:
      fprintf(fp, "AB_BA");
      break;
    case BA_AB:
      fprintf(fp, "BA_AB");
      break;
    case BA_BA:
      fprintf(fp, "BA_BA");
      break;
    case XX_XX:
      fprintf(fp, "XX_XX");
      break;
  }
}

int CompareEdgeOrientations(EdgeCGW_T * edge1,
                            EdgeCGW_T * edge2)
{
  if(edge1->orient != edge2->orient)
  {
    fprintf(stderr, "!!!!!!! Edge orientations not supported by fragment matepairs:");
    PrintOrientation(stderr, edge1->orient);
    fprintf(stderr, " (edge) vs. ");
    PrintOrientation(stderr, edge2->orient);
    fprintf(stderr, " (frags\n");
    /* switch(edge1->orient)
    {
      
    } */
    return 1;
  }
  return 0;
}


int CompareEdges(EdgeCGW_T * edge1,
                 EdgeCGW_T * edge2)
{
  int retVal = 0;
  retVal |= CompareEdgeMeans(edge1, edge2);
  // retVal |= CompareEdgeVariances(edge1, edge2);
  retVal |= CompareEdgeOrientations(edge1, edge2);
  return retVal;
}


/*
  Function to verify that all inter-contig edges are correct
 */
int CheckAllEdgesForChunk(ScaffoldGraphT * graph,
                          ChunkInstanceT * chunk,
                          int recurse,
                          int doCanonical,
                          int fixBadOnes)
{
  GraphEdgeIterator edges;
  CIEdgeT * edge;
  EdgeCGW_T myEdge;
  int retVal = 0;
  
  switch(chunk->type)
  {
    case DISCRIMINATORUNIQUECHUNK_CGW:
    case UNRESOLVEDCHUNK_CGW:
    case UNIQUECHUNK_CGW:
    case RESOLVEDREPEATCHUNK_CGW:
      break;
    case CONTIG_CGW:
    case UNIQUECONTIG_CGW:
    case RESOLVEDCONTIG_CGW:
    case UNRESOLVEDCONTIG_CGW:
    {
      /* loop over all contig edges
         if(recurse)
           loop over all unitigs
           recurse for each unitig
      */
      InitGraphEdgeIterator(graph->RezGraph, chunk->id, ALL_END, ALL_EDGES, GRAPH_EDGE_RAW_ONLY, &edges);
      while((edge = NextGraphEdgeIterator(&edges)) != NULL)
      {
        int isA = (edge->idA == chunk->id);
        CDS_CID_t thiscid = (isA? edge->idA: edge->idB);
        CDS_CID_t othercid = (isA? edge->idB: edge->idA);
        ChunkInstanceT * otherChunk = GetGraphNode(graph->RezGraph, othercid);
        
        // RAW EDGES ONLY
        assert(edge->flags.bits.isRaw);
        
        // canonical edges only
        if(thiscid == othercid || (doCanonical && thiscid > othercid))
          continue;

        if(edge->edgesContributing != 1)
          continue;
        
        if(edge->fragA == NULLINDEX || edge->fragB == NULLINDEX)
          continue;

        PopulateChunkEdgeBasics(graph,
                                GetCIFragT(graph->CIFrags,
                                           (isA? edge->fragA: edge->fragB)),
                                chunk,
                                GetCIFragT(graph->CIFrags,
                                           (isA? edge->fragB: edge->fragA)),
                                otherChunk,
                                GetDistT(graph->Dists, edge->distIndex),
                                &myEdge);
        if(CompareEdges(edge, &myEdge))
        {
          retVal = 1;
          if(fixBadOnes)
          {
            edge->distance.mean = myEdge.distance.mean;
            edge->distance.variance = myEdge.distance.variance;
            edge->orient = myEdge.orient;
          }
        }
      }
    }
    break;
    case REAL_SCAFFOLD:
    case OUTPUT_SCAFFOLD:
    case SCRATCH_SCAFFOLD:
      break;
  }
  return retVal;
}


int ValidateAllContigEdges(ScaffoldGraphT * graph, int fixBadOnes)
{
  CDS_CID_t sID;
  int retVal = 0;
  
  for(sID = 0; sID < GetNumCIScaffoldTs(graph->CIScaffolds); sID++)
  {
    CIScaffoldT * scaffold = GetCIScaffoldT(graph->CIScaffolds, sID);
    CIScaffoldTIterator contigs;
    ChunkInstanceT * contig;

    if(isDeadCIScaffoldT(scaffold) || scaffold->type != REAL_SCAFFOLD)
      continue;

    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &contigs);
    while((contig = NextCIScaffoldTIterator(&contigs)) != NULL)
    {
      assert(contig->scaffoldID == sID);
      retVal |= CheckAllEdgesForChunk(graph, contig, 0, 1, fixBadOnes);
    }
  }
  
  return retVal;
}


typedef struct
{
  CDS_CID_t   scaffoldID;
  cds_float32 minOffsetScaffoldB;
  cds_float32 maxOffsetScaffoldB;
  cds_float32 minOffset;
  cds_float32 maxOffset;
  ChunkOrientationType orient;
  int         orientationsConsistent;
  int         weight;
} ScfLink;

VA_DEF(ScfLink)

int ScfLinkHashFn(const void * item, int length)
{
  return Hash_AS((uint8 *) item, length, 37);
}

int ScfLinkCompareFn(const void * item1, const void * item2)
{
  const ScfLink * scfl1 = (const ScfLink *) item1;
  const ScfLink * scfl2 = (const ScfLink *) item2;

  return(scfl1->scaffoldID - scfl2->scaffoldID);
}

void PrintFragPairAndEdge(CIFragT * fragA,
                          LengthT * offsetA,
                          FragOrient orientA,
                          CIFragT * fragB,
                          LengthT * offsetB,
                          FragOrient orientB,
                          CIEdgeT * edge)
{
  fprintf(stdout, "%8" F_CIDP " %7" F_CIDP " %.f %s\t%8" F_CIDP " %7" F_CIDP " %.f %s\t%.f     ",
          fragA->iid, fragA->contigID,
          offsetA->mean, (orientA == A_B) ? "A_B" : "B_A",
          fragB->iid, fragB->contigID,
          offsetB->mean, (orientB == A_B) ? "A_B" : "B_A",
          edge->distance.mean);
  switch(edge->orient)
  {
    case AB_AB:
      fprintf(stdout, "AB_AB");
      break;
    case AB_BA:
      fprintf(stdout, "AB_BA");
      break;
    case BA_AB:
      fprintf(stdout, "BA_AB");
      break;
    case BA_BA:
      fprintf(stdout, "BA_BA");
      break;
    case XX_XX:
      fprintf(stdout, "XX_XX");
      break;
  }
  fprintf(stdout, "\n");
}


void PrintScaffoldConnectivity(ScaffoldGraphT * graph,
                               CIScaffoldT * scaffold,
                               CDS_CID_t otherScaffoldID)
{
  /*
    Loop over all contigs in scaffold
      Loop over all raw edges of contig
        For inter-scaffold edges
          Interested in:
            interval in scaffold
            orientation & distance to other scaffold
            interval in other scaffold

    How?
      hashtable of structures of scaffold interval & other scaffold id
      pre-allocate total number of scaffolds for array & HT
   */
  VA_TYPE(ScfLink) * links;
  HashTable_AS * linkHT;
  CIScaffoldTIterator ciTemp;

  if(scaffold->id == NULLINDEX)
    return;
  
  linkHT = CreateHashTable_AS(GetNumGraphNodes(graph->ScaffoldGraph),
                              ScfLinkHashFn, ScfLinkCompareFn);
  if(linkHT == NULL)
  {
    fprintf(stderr, "Failed to allocate scaffold hashtable!\n");
    return;
  }
  links = CreateVA_ScfLink(GetNumGraphNodes(graph->ScaffoldGraph));
  if(links == NULL)
  {
    fprintf(stderr, "Failed to allocate scaffold array!\n");
    return;
  }

  if(otherScaffoldID != NULLINDEX)
  {
    fprintf(stdout, "Edges between scaffolds " F_CID " and " F_CID ":\n",
            scaffold->id, otherScaffoldID);
    fprintf(stdout, "FragID  ContigID Scf5p ScfOrient "
            "FragID  ContigID Scf5p ScfOrient\t"
            "ScfEdgeLength ScfEdgeOrient\n");
  }
  
  // iterate over contigs in scaffold
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &ciTemp);
  while(NextCIScaffoldTIterator(&ciTemp))
  {
    ContigT * chunkA;
    GraphEdgeIterator edges;
    CIEdgeT * edge;
    CIEdgeT myEdge;

    // get the left chunk
    if((chunkA = GetGraphNode(graph->RezGraph, ciTemp.curr)) == NULL)
    {
      fprintf(stderr, "Failed to get contig " F_CID " from graph!\n", ciTemp.curr);
      return;
    }

    // loop over edges of this chunk
    InitGraphEdgeIterator(graph->RezGraph,
                          chunkA->id,
                          ALL_END,
                          ALL_EDGES,
                          GRAPH_EDGE_RAW_ONLY,
                          &edges);
    while((edge = NextGraphEdgeIterator(&edges)) != NULL)
    {
      int isA = (edge->idA == chunkA->id);
      CDS_CID_t idB = (isA? edge->idB: edge->idA);
      ChunkInstanceT * chunkB = GetGraphNode(graph->RezGraph, idB);
      CIFragT * fragA = GetCIFragT(graph->CIFrags,
                                   (isA? edge->fragA: edge->fragB));
      CIFragT * fragB = GetCIFragT(graph->CIFrags,
                                   (isA? edge->fragB: edge->fragA));
      DistT * dist = GetDistT(graph->Dists, edge->distIndex);
      CIScaffoldT * scaffoldB = GetCIScaffoldT(graph->CIScaffolds,
                                               chunkB->scaffoldID);
      LengthT offsetA;
      LengthT offsetB;
      FragOrient orientA;
      FragOrient orientB;
      ScfLink * link;
      ScfLink newLink;

      if(chunkA == NULL || chunkB == NULL ||
         chunkB->scaffoldID == chunkA->scaffoldID ||
         chunkB->scaffoldID == NULLINDEX ||
         fragA == NULL || fragB == NULL ||
         (otherScaffoldID != NULLINDEX &&
          otherScaffoldID != chunkB->scaffoldID))
        continue;

      GetFragment5pAndOrientationInChunk(graph,
                                         fragA, &offsetA, &orientA,
                                         scaffold);

      GetFragment5pAndOrientationInChunk(graph,
                                         fragB, &offsetB, &orientB,
                                         scaffoldB);

      PopulateChunkEdgeBasics(graph,
                              fragA,
                              scaffold,
                              fragB,
                              scaffoldB,
                              dist,
                              &myEdge);

      if(otherScaffoldID != NULLINDEX)
      {
        PrintFragPairAndEdge(fragA, &offsetA, orientA,
                             fragB, &offsetB, orientB,
                             &myEdge);
      }
      // see link to scaffoldB is already present
      if((link = LookupInHashTable_AS(linkHT,
                                      (void *) &(scaffoldB->id),
                                      sizeof(scaffoldB->id))) == NULL)
      {
        // start a new link
        newLink.scaffoldID = scaffoldB->id;
        newLink.minOffset = newLink.maxOffset = offsetA.mean;
        newLink.minOffsetScaffoldB = newLink.maxOffsetScaffoldB = offsetB.mean;
        newLink.orient = myEdge.orient;
        newLink.orientationsConsistent = 1;
        newLink.weight = 1;
        AppendVA_ScfLink(links, &newLink);
        link = GetVA_ScfLink(links, GetNumVA_ScfLink(links) - 1);
        InsertInHashTable_AS(linkHT,
                             (void *) &(scaffoldB->id),
                             sizeof(scaffoldB->id),
                             (void *) link);
      }
      else
      {
        // update existing link
        link->minOffset = min(offsetA.mean, link->minOffset);
        link->maxOffset = max(offsetA.mean, link->minOffset);
        link->minOffsetScaffoldB = min(offsetB.mean, link->minOffsetScaffoldB);
        link->maxOffsetScaffoldB = max(offsetB.mean, link->minOffsetScaffoldB);
        link->weight++;
        if(link->orient != myEdge.orient)
          link->orientationsConsistent = 0;
      }
    }
  }

  // now print out
  if(otherScaffoldID == NULLINDEX)
  {
    HashTable_Iterator_AS iterator;
    void * key;
    void * value;
    
    fprintf(stdout, "*** Links from scaffold " F_CID ":\n", scaffold->id);
    InitializeHashTable_Iterator_AS(linkHT, &iterator);
    while(NextHashTable_Iterator_AS(&iterator, &key, &value) == HASH_SUCCESS)
    {
      ScfLink * link;
      if(key == NULL || value == NULL)
        continue;
      link = (ScfLink *) value;

      fprintf(stdout,
              "\tweight %d from (%.f,%.f) to scaffold " F_CID " (%.f,%.f) orient: ",
              link->weight, link->minOffset, link->maxOffset,
              link->scaffoldID,
              link->minOffsetScaffoldB, link->maxOffsetScaffoldB);
      switch(link->orient)
      {
        case AB_AB:
          fprintf(stdout, "AB_AB");
          break;
        case AB_BA:
          fprintf(stdout, "AB_BA");
          break;
        case BA_AB:
          fprintf(stdout, "BA_AB");
          break;
        case BA_BA:
          fprintf(stdout, "BA_BA");
          break;
        case XX_XX:
          fprintf(stdout, "XX_XX");
          break;
      }
      if(!link->orientationsConsistent)
        fprintf(stdout, ", (not consistent)");
      fprintf(stdout, "\n");
    }
  }
  fprintf(stdout, "\n");
  
  DeleteVA_ScfLink(links);
  DeleteHashTable_AS(linkHT);
  
}



void PrintScaffoldConnectivityForIID(ScaffoldGraphT * graph,
                                     CDS_CID_t scaffoldID)
{
  PrintScaffoldConnectivity(graph,
                            GetCIScaffoldT(graph->CIScaffolds, scaffoldID),
                            NULLINDEX);
}

void PrintConnectivityBetweenScaffolds(ScaffoldGraphT * graph,
                                       CDS_CID_t scaffoldID1,
                                       CDS_CID_t scaffoldID2)
{
  PrintScaffoldConnectivity(graph,
                            GetCIScaffoldT(graph->CIScaffolds, scaffoldID1),
                            scaffoldID2);
}


static int CDS_CID_compare(const CDS_CID_t * a, const CDS_CID_t * b)
{
  return (int) *a - *b;
}

int CountUniqueIDs(VA_TYPE(CDS_CID_t) * ids)
{
  int i, count;
  CDS_CID_t * id;
  CDS_CID_t * lastID = NULL;
  qsort(GetVA_CDS_CID_t(ids, 0),
        GetNumVA_CDS_CID_t(ids),
        sizeof(CDS_CID_t),
        (int (*) (const void *, const void *)) CDS_CID_compare);
  for(count = 0, i = 0; i < GetNumVA_CDS_CID_t(ids); i++)
  {
    id = GetVA_CDS_CID_t(ids, i);
    if(lastID == NULL || *id != *lastID)
      count++;
    lastID = id;
  }
  return count;
}

/* loop over contigs
       assert that each contig has one unitig with same id
       get unitig
       check Astat
       count number of contigs linked to this one off each end
       print
       if high Astat & links to multiple contigs off each end - repetitive?
*/
void DetectRepetitiveContigs(ScaffoldGraphT * graph)
{
  ChunkInstanceT * contig;
  GraphNodeIterator contigIterator;
  ChunkInstanceT * unitig;
  ContigTIterator unitigIterator;
  int numContigs = 0;
  int numRepetitiveContigs = 0;
  VA_TYPE(CDS_CID_t) * aEndIDs;
  VA_TYPE(CDS_CID_t) * bEndIDs;
  int numOffAEnd;
  int numOffBEnd;
  int maxOffEnd;
  CDS_COORD_t length = CDS_COORD_MIN;

  fprintf(stdout,
          "Identifying potentially repetitive contigs with Astat > %.f\n",
          CGB_UNIQUE_CUTOFF);
  fprintf(stdout, "  ContigID   Astat   #AEndContigs   #BendContigs   length\n");

  aEndIDs = CreateVA_CDS_CID_t(10000);
  bEndIDs = CreateVA_CDS_CID_t(10000);
  
  InitGraphNodeIterator(&contigIterator,
                        graph->ContigGraph,
                        GRAPH_NODE_UNIQUE_ONLY);
  while((contig = NextGraphNodeIterator(&contigIterator)) != NULL)
  {
    int numUnitigsInContig = 0;
    int coverageStat = 0;

    numContigs++;
    ResetVA_CDS_CID_t(aEndIDs);
    ResetVA_CDS_CID_t(bEndIDs);

    // Iterate over unitigs in contig & add data to contig instrumenter
    InitContigTIterator(graph, contig->id, TRUE, FALSE, &unitigIterator);
    while((unitig = NextContigTIterator(&unitigIterator)) != NULL)
    {
      numUnitigsInContig++;
      coverageStat = unitig->info.CI.coverageStat;
      length = unitig->bpLength.mean;
    }
    if(numUnitigsInContig > 1)
    {
      fprintf(stderr,
              "Only call DetectRepetitiveContigs with unitigs = contigs!\n");
      return;
    }
    if(coverageStat > CGB_UNIQUE_CUTOFF)
    {
      // loop over contig edges
      GraphEdgeIterator edges;
      CIEdgeT * edge;

      InitGraphEdgeIterator(graph->ContigGraph,
                            contig->id,
                            ALL_END,
                            ALL_EDGES,
                            GRAPH_EDGE_RAW_ONLY,
                            &edges);
      while((edge = NextGraphEdgeIterator(&edges)) != NULL)
      {
        int isA = (edge->idA == contig->id);
        CDS_CID_t othercid = (isA? edge->idB: edge->idA);

        if((isA && (edge->orient == AB_AB || edge->orient == AB_BA)) ||
           (!isA && (edge->orient == AB_BA || edge->orient == BA_BA)))
          AppendVA_CDS_CID_t(bEndIDs, &othercid);
        else
          AppendVA_CDS_CID_t(aEndIDs, &othercid);
      }

      // now sort & examine edges
      {
        if(GetNumVA_CDS_CID_t(aEndIDs) > 1)
          numOffAEnd = CountUniqueIDs(aEndIDs);
        else
          numOffAEnd = GetNumVA_CDS_CID_t(aEndIDs);

        if(GetNumVA_CDS_CID_t(bEndIDs) > 1)
          numOffBEnd = CountUniqueIDs(bEndIDs);
        else
          numOffBEnd = GetNumVA_CDS_CID_t(bEndIDs);

        maxOffEnd = (int) (15.04 + 3 * 16.99 + .5);
        if(numOffBEnd >= maxOffEnd || numOffAEnd >= maxOffEnd)
        {
          numRepetitiveContigs++;
          fprintf(stdout,"%10" F_CIDP "   %5d     %10d     %10d     %10" F_COORDP "\n",
                  contig->id, coverageStat, numOffAEnd, numOffBEnd, length);
        }
      }
    }
  }

  fprintf(stdout,
          "Graph contains %d contigs, %d of which may be repetitive.\n",
          numContigs, numRepetitiveContigs);
  
  DeleteVA_CDS_CID_t(aEndIDs);
  DeleteVA_CDS_CID_t(bEndIDs);
}


void DoSomethingWithUnitigsInScaffolds(ScaffoldGraphT * graph)
{
  GraphNodeIterator scaffolds;
  CIScaffoldT * scaffold;

   // loop over all scaffolds in the graph
  InitGraphNodeIterator(&scaffolds,
                        graph->ScaffoldGraph,
                        GRAPH_NODE_DEFAULT);
  while(NULL != (scaffold = NextGraphNodeIterator(&scaffolds)))
  {
    if(scaffold->flags.bits.isDead == FALSE &&
       scaffold->type == REAL_SCAFFOLD)
    {
      CIScaffoldTIterator cisTemp;
      
      // loop over all contigs in the scaffold
      InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &cisTemp);
      while( NextCIScaffoldTIterator(&cisTemp) && cisTemp.next != NULLINDEX)
      {
        ChunkInstanceT * contig;
        ContigTIterator unitigs;
        ChunkInstanceT * unitig;
        FragOrient contigOrient;

        if((contig = GetGraphNode(graph->RezGraph, cisTemp.curr)) == NULL)
        {
          fprintf(stderr, "Failed to get contig " F_CID " from graph!\n",
                  cisTemp.curr);
          return;
        }
        
        contigOrient =
          (contig->offsetAEnd.mean < contig->offsetBEnd.mean) ? A_B : B_A;
        // Iterate over unitigs in contig & add data to contig instrumenter
        InitContigTIterator(graph, cisTemp.curr, TRUE, FALSE, &unitigs);
        while((unitig = NextContigTIterator(&unitigs)) != NULL)
        {
          // if unitig is not a surrogate
          if(!unitig->flags.bits.isStoneSurrogate &&
             !unitig->flags.bits.isWalkSurrogate)
          {
            FragOrient unitigOrient;
            unitigOrient =
              (unitig->offsetAEnd.mean < unitig->offsetBEnd.mean) ? A_B : B_A;
            
            // do something
          }
        }
      }
    }
  }
 
}

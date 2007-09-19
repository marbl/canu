
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
static char CM_ID[] = "$Id: ContigT_CGW.c,v 1.16 2007-09-19 21:54:24 skoren Exp $";

//#define DEBUG 1
//#define TRY_IANS_EDGES

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
#include "MultiAlignment_CNS.h"  //  What a waste, only for SequenceComplement()

#ifdef TRY_IANS_SEDGES
#include "AS_CGW_EdgeDiagnostics.h"
#endif


void CheckContigs()
{
  GraphNodeIterator CIs;
  NodeCGW_T * contig;
  
  InitGraphNodeIterator(&CIs, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((contig = NextGraphNodeIterator(&CIs)) != NULL)
    {
      double meanDelta = fabs(contig->offsetAEnd.mean - contig->offsetBEnd.mean);
      double varianceDelta = fabs(contig->offsetAEnd.variance -
                                   contig->offsetBEnd.variance);
      if(meanDelta + .5 < contig->bpLength.mean ||
         meanDelta - .5 > contig->bpLength.mean ||
         varianceDelta + .5 < contig->bpLength.variance ||
         varianceDelta - .5 > contig->bpLength.variance)
        {
          fprintf(GlobalData->stderrc,
                  "Contig " F_CID " length (%f,%f) doesn't match offset difference (%f,%f)\n",
                  contig->id, contig->bpLength.mean, contig->bpLength.variance,
                  meanDelta, varianceDelta);
        }
    }
  
}

void
dumpContigInfo(ChunkInstanceT *contig) {
  int           contigOrientation;
  MultiAlignT  *ma;
  char         *seq1;
  int           len1;

  VA_TYPE(char) *consensus = CreateVA_char(2048);
  VA_TYPE(char) *quality   = CreateVA_char(2048);

  fprintf( stderr, "*********************** contig analysis **************************\n");
  fprintf( stderr, "analyzing contig: %d\n", contig->id);

  if (contig->offsetAEnd.mean < contig->offsetBEnd.mean)
    contigOrientation = 0;
  else
    contigOrientation = 1;

  fprintf(stderr, "contig orientation: %d\t length: %d  contig offsetAEnd: %d\t offsetBEnd: %d\n",
          contigOrientation,
          (int)contig->bpLength.mean,
          (int)contig->offsetAEnd.mean,
          (int)contig->offsetBEnd.mean);

  ma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, ScaffoldGraph->RezGraph->type == CI_GRAPH); 

  // Get the consensus sequences for the contig from the Store
  GetConsensus(ScaffoldGraph->ContigGraph, contig->id, consensus, quality);

  seq1 = Getchar(consensus, 0);
  len1 = strlen(seq1);

  if (contigOrientation == 1)
    SequenceComplement(seq1, NULL);
  
  if (len1 < 5000) {
    fprintf( stderr, ">contig%d consensus seq (flipped to reflect scaff orientation)\n", contig->id);
    fprintf( stderr, "%s\n", seq1);
  } else {
    char tmpchar = seq1[2500];
    seq1[2500] = '\0';

    fprintf( stderr, ">contig%d left end\n", contig->id);
    fprintf( stderr, "%s\n", seq1);

    seq1[2500] = tmpchar;

    fprintf( stderr, ">contig%d right end\n", contig->id);
    fprintf( stderr, "%s\n", seq1 + len1 - 2501);
  }

#if 1
  int numUnitigs = GetNumIntUnitigPoss(ma->u_list);
  fprintf( stderr, "number unitigs: %d\n", numUnitigs);

  int i;
  for (i = 0; i < numUnitigs; i++) {
    IntUnitigPos *upos = GetIntUnitigPos( ma->u_list, i);
    ChunkInstanceT *unitig = GetGraphNode( ScaffoldGraph->CIGraph, upos->ident);
    MultiAlignT *uma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, ScaffoldGraph->CIGraph->type == CI_GRAPH);
    IntMultiPos *ump;
    int icntfrag;

    fprintf( stderr, "  unitig: %d\t num frags: %ld surrogate: %d\n", unitig->id, GetNumIntMultiPoss(uma->f_list),
             (unitig->flags.bits.isStoneSurrogate || unitig->flags.bits.isWalkSurrogate));

    if (unitig->flags.bits.isStoneSurrogate ||
        unitig->flags.bits.isWalkSurrogate) {
      fprintf (stderr, "  surrogate unitig offsetAEnd: %f, offsetBEnd: %f\n", unitig->offsetAEnd.mean, unitig->offsetBEnd.mean);

      unitig = GetGraphNode( ScaffoldGraph->CIGraph, unitig->info.CI.baseID);
      fprintf ( stderr, "  using original unitig: %d\n", unitig->id);
      uma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, 
                                          ScaffoldGraph->CIGraph->type == CI_GRAPH); 
    }

    // now print out info on the frags in the unitig
    for (icntfrag = 0; icntfrag < GetNumIntMultiPoss(uma->f_list); icntfrag++) {
      IntMultiPos *imp = GetIntMultiPos(uma->f_list, icntfrag);
      InfoByIID   *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, imp->ident);
      CIFragT     *frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

      //assert(info->set);

      fprintf(stderr, "    frag: %6d\t contig pos (5p, 3p): %6d, %6d  type: %c\n",
              imp->ident, (int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean, frag->type);
    }
  }
#endif


#if 1
  CIEdgeT * e;
  GraphEdgeIterator edges;

  //  FALSE == ITERATOR_VERBOSE

  InitGraphEdgeIterator(ScaffoldGraph->RezGraph, contig->id, ALL_END, ALL_EDGES, FALSE, &edges);
  while((e = NextGraphEdgeIterator(&edges)) != NULL)
    PrintGraphEdge( stderr, ScaffoldGraph->RezGraph, "Analyzing edge", e, 0);
#endif

  DeleteVA_char(consensus);
  DeleteVA_char(quality);
}


void
GetContigPositionInScaffold(ChunkInstanceT *contig, int *left_end, int *right_end, 
                            int *contigScaffoldOrientation) {  
  if (contig->offsetAEnd.mean <= contig->offsetBEnd.mean) {
    *left_end = contig->offsetAEnd.mean;
    *right_end = contig->offsetBEnd.mean;
    *contigScaffoldOrientation = 0;
  } else {
    *left_end = contig->offsetBEnd.mean;
    *right_end = contig->offsetAEnd.mean;
    *contigScaffoldOrientation = 1;
  }
}


static
void
GetFragmentPositionInScaffoldFromContig(CIFragT *frag, int *left_end, int *right_end, 
                                        int *fragmentScaffoldOrientation, 
                                        int contigLeftEnd, int contigRightEnd, int contigScaffoldOrientation) {
  if (contigScaffoldOrientation == 0) {
    // contig is direct in scaffold
    if (frag->contigOffset5p.mean < frag->contigOffset3p.mean) {
      // frag is direct in contig
      *left_end = contigLeftEnd + frag->contigOffset5p.mean;
      *right_end = contigLeftEnd + frag->contigOffset3p.mean;
      *fragmentScaffoldOrientation = 0;
    } else {
      // frag is reversed in contig
      *left_end = contigLeftEnd + frag->contigOffset3p.mean;
      *right_end = contigLeftEnd + frag->contigOffset5p.mean;
      *fragmentScaffoldOrientation = 1;
    }
  } else {
    // contig is reversed in scaffold
    if (frag->contigOffset5p.mean < frag->contigOffset3p.mean) {
      // frag is direct in contig
      *left_end = contigRightEnd - frag->contigOffset3p.mean;
      *right_end = contigRightEnd - frag->contigOffset5p.mean;
      *fragmentScaffoldOrientation = 1;
    } else {
      // frag is reversed in contig
      *left_end = contigRightEnd - frag->contigOffset5p.mean;
      *right_end = contigRightEnd - frag->contigOffset3p.mean;
      *fragmentScaffoldOrientation = 0;
    }
  }
}


void
GetFragmentPositionInScaffold(CIFragT *frag, int *left_end, int *right_end, 
                                   int *fragmentScaffoldOrientation) {
  ContigT *containingContig = GetGraphNode(ScaffoldGraph->ContigGraph, frag->contigID);
  int contigLeftEnd, contigRightEnd, contigScaffoldOrientation;
  
  GetContigPositionInScaffold( containingContig, &contigLeftEnd, &contigRightEnd, &contigScaffoldOrientation);

  GetFragmentPositionInScaffoldFromContig( frag, left_end, right_end, fragmentScaffoldOrientation,
                                           contigLeftEnd, contigRightEnd, contigScaffoldOrientation);
}


void DumpContig(FILE *stream, ScaffoldGraphT *graph, ContigT *contig, int raw){
  int numCI;
  GraphEdgeIterator edges;
  ContigTIterator CIs;
  ChunkInstanceT *CI;
  CIEdgeT *edge;
  int flags = GRAPH_EDGE_DEFAULT;

  assert(contig->type == CONTIG_CGW);

  if(raw)
    flags |= GRAPH_EDGE_RAW_ONLY;

  fprintf(stream, "* Contig " F_CID " sc:" F_CID " [" F_COORD "," F_COORD "] aoff:%d boff:%d ANext:" F_CID " BNext:" F_CID " len:%d AEnd:" F_CID " BEnd:" F_CID " numCI:%d\nCIs:\n",
          contig->id, 
          contig->scaffoldID,
          contig->aEndCoord,
          contig->bEndCoord,
          (int)contig->offsetAEnd.mean,
          (int)contig->offsetBEnd.mean,
          contig->AEndNext, contig->BEndNext, 
          (int)contig->bpLength.mean,
          contig->info.Contig.AEndCI,
          contig->info.Contig.BEndCI,
          contig->info.Contig.numCI);

#ifdef DEBUG_CONTIG
  fprintf(stderr, "* Contig " F_CID " [" F_COORD "," F_COORD "] ANext:" F_CID " BNext:" F_CID " length:%d AEnd:" F_CID " BEnd:" F_CID " numCI:%d\nCIs:\n",
          contig->id, 
          contig->aEndCoord,
          contig->bEndCoord,
          contig->AEndNext, contig->BEndNext, (int)contig->bpLength.mean,
          contig->info.Contig.AEndCI,
          contig->info.Contig.BEndCI,
          contig->info.Contig.numCI);
  fflush(stderr);
#endif
  InitContigTIterator(graph,contig->id, TRUE, FALSE, &CIs);
  numCI = 0;
  while((CI = NextContigTIterator(&CIs)) != NULL){
    int isSurrogate = IsSurrogateNode(CI);
    numCI++;
    fprintf(stream,"\t%c:" F_CID " cov:%d [" F_COORD "," F_COORD "] aoff:%d boff:%d anxt:" F_CID " bnxt:" F_CID " ctg:" F_CID " scf:" F_CID "\n",
            (isSurrogate?'s':'u'),
            CI->id,
            CI->info.CI.coverageStat,
            CI->aEndCoord,
            CI->bEndCoord,
            (int)CI->offsetAEnd.mean,
            (int)CI->offsetBEnd.mean,
            CI->AEndNext,
            CI->BEndNext,
            CI->info.CI.contigID,
            CI->scaffoldID);
#ifdef DEBUG_CONTIG
    fprintf(stderr,"\t%c:" F_CID " cov:%d [" F_COORD "," F_COORD "] aoff:%d boff:%d anxt:" F_CID " bnxt:" F_CID " ctg:" F_CID " scf:" F_CID "\n",
            (isSurrogate?'s':'u'),
            CI->id,
            CI->info.CI.coverageStat,
            CI->aEndCoord,
            CI->bEndCoord,
            (int)CI->offsetAEnd.mean,
            (int)CI->offsetBEnd.mean,
            CI->AEndNext,
            CI->BEndNext,
            CI->info.CI.contigID,
            CI->scaffoldID);
    fflush(stderr);
#endif
  }

  fprintf(stream,"\t%s Edges from A End:\n", (raw?" Raw ":" Merged "));

  InitGraphEdgeIterator(graph->ContigGraph, contig->id, A_END, ALL_EDGES, flags, &edges);
  while((edge = NextGraphEdgeIterator(&edges)) != NULL){
    PrintGraphEdge(stream,graph->ContigGraph,"\t",edge, contig->id);
  }
  fprintf(stream,"\t%s Edges from B End:\n", (raw?" Raw ":" Merged "));
  InitGraphEdgeIterator(graph->ContigGraph, contig->id, B_END, ALL_EDGES, flags, &edges);
  while((edge = NextGraphEdgeIterator(&edges)) != NULL){
    PrintGraphEdge(stream,graph->ContigGraph,"\t",edge, contig->id);
  }
}

void DumpContigInScfContext(FILE *stream, ScaffoldGraphT *graph,
                            ContigT *contig, int raw){
  int numCI;
  GraphEdgeIterator edges;
  ContigTIterator CIs;
  ChunkInstanceT *CI;
  CIEdgeT *edge;
  int flags = GRAPH_EDGE_DEFAULT;

  assert(contig->type == CONTIG_CGW);

  if(raw)
    flags |= GRAPH_EDGE_RAW_ONLY;

  fprintf(stream, "* Contig " F_CID " sc:" F_CID " [" F_COORD "," F_COORD "] aoff:(%d,%e) boff:(%d,%e) ANext:" F_CID " BNext:" F_CID " len:%d AEnd:" F_CID " BEnd:" F_CID " numCI:%d\nCIs:\n",
          contig->id, 
          contig->scaffoldID,
          contig->aEndCoord,
          contig->bEndCoord,
          (int32)contig->offsetAEnd.mean,
          contig->offsetAEnd.variance,
          (int32)contig->offsetBEnd.mean,
          contig->offsetBEnd.variance,
          contig->AEndNext, contig->BEndNext, 
          (int)contig->bpLength.mean,
          contig->info.Contig.AEndCI,
          contig->info.Contig.BEndCI,
          contig->info.Contig.numCI);

#ifdef DEBUG_CONTIG
  fprintf(stderr, "* Contig " F_CID " [" F_COORD "," F_COORD "] ANext:" F_CID " BNext:" F_CID " length:%d AEnd:" F_CID " BEnd:" F_CID " numCI:%d\nCIs:\n",
          contig->id, 
          contig->aEndCoord,
          contig->bEndCoord,
          contig->AEndNext, contig->BEndNext, (int)contig->bpLength.mean,
          contig->info.Contig.AEndCI,
          contig->info.Contig.BEndCI,
          contig->info.Contig.numCI);
  fflush(stderr);
#endif
  InitContigTIterator(graph,contig->id, TRUE, FALSE, &CIs);
  numCI = 0;
  while((CI = NextContigTIterator(&CIs))!=NULL){
    int isSurrogate = IsSurrogateNode(CI);
    numCI++;
    fprintf(stream,"\t%c:" F_CID " cov:%d [" F_COORD "," F_COORD "] aoff:(%d,%e) boff:(%d,%e) anxt:" F_CID " bnxt:" F_CID " ctg:" F_CID " scf:" F_CID "\n",
            (isSurrogate?'s':'u'),
            CI->id,
            CI->info.CI.coverageStat,
            CI->aEndCoord,
            CI->bEndCoord,
            (int)CI->offsetAEnd.mean,
            CI->offsetAEnd.variance,
            (int)CI->offsetBEnd.mean,
            CI->offsetBEnd.variance,
            CI->AEndNext,
            CI->BEndNext,
            CI->info.CI.contigID,
            CI->scaffoldID);
#ifdef DEBUG_CONTIG
    fprintf(stderr,"\t%c:" F_CID " cov:%d [" F_COORD "," F_COORD "] aoff:%d boff:%d anxt:" F_CID " bnxt:" F_CID " ctg:" F_CID " scf:" F_CID "\n",
            (isSurrogate?'s':'u'),
            CI->id,
            CI->info.CI.coverageStat,
            CI->aEndCoord,
            CI->bEndCoord,
            (int)CI->offsetAEnd.mean,
            (int)CI->offsetBEnd.mean,
            CI->AEndNext,
            CI->BEndNext,
            CI->info.CI.contigID,
            CI->scaffoldID);
    fflush(stderr);
#endif
  }
  fprintf(stream,"\t%s Edges from A End:\n", (raw?" Raw ":" Merged "));

  InitGraphEdgeIterator(graph->ContigGraph, contig->id, A_END, ALL_EDGES, flags, &edges);
  while((edge = NextGraphEdgeIterator(&edges))!=NULL){
    PrintContigEdgeInScfContext(stream,graph->ContigGraph,"\t",edge, contig->id);
  }
  fprintf(stream,"\t%s Edges from B End:\n", (raw?" Raw ":" Merged "));
  InitGraphEdgeIterator(graph->ContigGraph, contig->id, B_END, ALL_EDGES, flags, &edges);
  while((edge = NextGraphEdgeIterator(&edges))!=NULL){
    PrintContigEdgeInScfContext(stream,graph->ContigGraph,"\t",edge, contig->id);
  }
}

void DumpContigs(FILE *stream, ScaffoldGraphT *graph, int raw){

  GraphNodeIterator CIs;
  NodeCGW_T *contig;

  fprintf(stream,"************Dumping Contigs ***********\n");

  InitGraphNodeIterator(&CIs, graph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((contig = NextGraphNodeIterator(&CIs)) != NULL){
    DumpContig(stream, graph, contig, raw);
  }
}


// When we can generate multi-alignments fro consensus sequences, this should die
#if 1 

/* Compute the offset and orientation of a fragment in its chunk and contig
   Offset is from 5p end of fragment to the end of the chunk/scaffold end in the
   direction of the 3p end of the fragment.

*/

int CIinContigOffsetAndOrientation(ScaffoldGraphT *graph, 
                                   CIFragT *frag,
                                   CDS_CID_t  cid, // input
                                   CIOrient   chunkOrient,    // input orientation of fragment in CI
                                   LengthT    *ciOffset,    // CI's offset within Contig 
                                   CIOrient   *ciOrient){    // CI's orientation within Contig

  ChunkInstanceT *CI = GetChunkInstanceT(graph->ChunkInstances, cid);
  ContigT *contig;
  CIOrient CIInContigOrient;

  if(CI->info.CI.contigID == NULLINDEX)
    return FALSE;

  contig = GetGraphNode(graph->ContigGraph, CI->info.CI.contigID);
  
  CIInContigOrient = GetNodeOrient(CI);


  /* Find the offset of the B end of the chunk within its chunkOfScaffolds.  The fragments
     offset is additive with this */

  switch(chunkOrient){
    case A_B:  // Fragment is A_B in CI
      switch(CIInContigOrient){
        case A_B:// Chunk is A_B in Contig
          /*    Contig   ------------------------------------------------------->
                CI      A_B     ---------------------->  
                Frag       A_B            ------>
                Offset                    |=====================================|
          */ 
	 
          *ciOrient = A_B;
          ciOffset->mean = 
            (contig->bpLength.mean -  CI->offsetBEnd.mean);
		 
          /*
            ciOffset->variance = 
            (contig->bpLength.variance - CI->offsetBEnd.variance);
          */

          // Take the proportinal part of the contig variance as the variance
          ciOffset->variance = (contig->bpLength.variance ) * (contig->bpLength.mean - ciOffset->mean)  / contig->bpLength.mean ;

          if(ciOffset->variance < 0.0){
            fprintf(stderr,"* 1) CIInContig ciOffset->variance = %g => 0.0\n",
                    ciOffset->variance);
            fprintf(stderr,"* 1) contig " F_CID " bpLength (%g,%g) CI " F_CID "  CI->offsetBEnd (%g,%g)\n",
                    CI->info.CI.contigID,
                    contig->bpLength.mean,
                    contig->bpLength.variance,
                    cid,
                    CI->offsetBEnd.mean,
                    CI->offsetBEnd.variance);

            ciOffset->variance = 0.0;
          }

          break;
        case B_A: // Chunk is B_A in Contig
          /*    Contig ------------------------------------------------------->
                Chunk     B_A  <----------------------   
                Frag      A_B         <------    
                Offset |====================|
          */ 
          *ciOrient = B_A;
          ciOffset->mean = CI->offsetBEnd.mean;
          ciOffset->variance = CI->offsetBEnd.variance;

          break;
        default:
          assert(0);
      }

      break;

    case B_A:  // Fragment is B_A in Chunk
      switch(CIInContigOrient){
        case A_B:// Chunk is A_B in Contig
          /*    Contig    ------------------------------------------------------->
                Chunk      A_B     ---------------------->  
                Frag       B_A            <------
                Offset |===========|
          */ 
	 
          *ciOrient = B_A;
          ciOffset->mean = CI->offsetAEnd.mean;
          ciOffset->variance = CI->offsetAEnd.variance;


          break;
        case B_A: // Chunk is B_A in Contig
          /*    Contig    ------------------------------------------------------->
                Chunk           <----------------------   B_A
                Frag                   ------>    B_A
                Offset                                 |=========================|
          */ 
          *ciOrient = A_B;
          ciOffset->mean =  
            (contig->bpLength.mean - CI->offsetAEnd.mean);

          ciOffset->variance =   
            (contig->bpLength.variance - CI->offsetAEnd.variance);

          if(ciOffset->variance < 0.0){
            fprintf(stderr,"* 2) CIInContig ciOffset->variance = %g => 0.0\n",
                    ciOffset->variance);
            fprintf(stderr,"* 2) contig " F_CID " bpLength (%g,%g) CI " F_CID "  CI->offsetAEnd (%g,%g)\n",
                    CI->info.CI.contigID,
                    contig->bpLength.mean,
                    contig->bpLength.variance,
                    cid,
                    CI->offsetAEnd.mean,
                    CI->offsetAEnd.variance);

            ciOffset->variance = 0.0;
          }
          break;

        default:
          assert(0);
      }


      break;
    default:
      assert(0);
  }

  return TRUE;

}


/* Using the CIEdges as a base, construct raw edges between contigs, and merge them */
int BuildContigEdges(ScaffoldGraphT *graph){
  CDS_CID_t sid;
  int oldSize = GetNumCIEdgeTs(graph->ContigEdges);
  CIOrient ciOrient, mciOrient;
  //CIOrient contigEdgeOrient;
  ChunkOrientationType contigEdgeOrient;
  LengthT ciOffset, mciOffset;
  LengthT distance;
  int edgesSucceeded = 0;
  CIEdgeT contigEdge;

  /* Recycle the SEdge VA */
  ResetCIEdgeT(graph->ContigGraph->edges);
  EnableRangeVA_CIEdgeT(graph->ContigGraph->edges, oldSize);
  
  /* Iterating over scaffolds rather than the Contigs should give us some
     locality of reference
  */
  for(sid = 0; sid < GetNumCIScaffoldTs(graph->CIScaffolds); sid++){
    CIScaffoldT *scaffold = GetCIScaffoldT(graph->CIScaffolds,sid);
    CIScaffoldTIterator Contigs;
    /* Iterate over all of the CIEdges */
    CIEdgeTIterator edges;
    CIEdgeT *edge;
    ChunkInstanceT *thisCI;
    ChunkInstanceT *thisContig;

    if(isDeadCIScaffoldT(scaffold))
      continue;

#if 0
    fprintf(stderr,"* BuildContigEdges Scaffold " F_CID "\n",
	    sid);
#endif
    assert(scaffold->flags.bits.containsCIs == 0); // iterate over scaffold containing CONTIGS

    //    fprintf(stderr,"* Building ContigEdges incident on Contigs of scaffold " F_CID "\n", sid);
    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &Contigs);

    
    while((thisContig = NextCIScaffoldTIterator(&Contigs)) != NULL){ /* Iterate over all CONTIGs in a Scaffold */

      assert(thisContig->scaffoldID == sid);

      /* Iterate over all CIs in a Contig */
      {
	ContigTIterator CIs;
#if 0
	CDS_CID_t nextcid = thisContig->info.Contig.AEndCI;
#endif
	InitContigTIterator(graph, thisContig->id, TRUE, FALSE, &CIs);
	while((thisCI = NextContigTIterator(&CIs)) != NULL){
#if 0
          CDS_CID_t thisCID = thisCI->id;
          fprintf(stderr,
                  "* Contig " F_CID " cid:" F_CID " nextcid:" F_CID "\n",
                  thisContig->id, thisCID, nextcid); 
#endif

	  InitCIEdgeTIterator(graph,  thisCI->id,   TRUE /* raw */, FALSE,  ALL_END,   ALL_EDGES, FALSE, &edges);// ONLY RAW
	  while((edge = NextCIEdgeTIterator(&edges)) != NULL){
	    int isA = (edge->idA == thisCI->id);
	    ChunkOrientationType edgeOrient;
	    CDS_CID_t otherCID = (isA? edge->idB: edge->idA);
	    ChunkInstanceT *otherCI = GetChunkInstanceT(graph->ChunkInstances, otherCID);
	    CIFragT *frag = GetCIFragT(graph->CIFrags, (isA? edge->fragA: edge->fragB));
	    CIFragT *otherFrag = GetCIFragT(graph->CIFrags, (isA? edge->fragB: edge->fragA));
	    FragOrient orient;
	    int CIok, mCIok;
	    CDS_CID_t contig, mcontig;

	    contig = thisCI->info.CI.contigID;
	    mcontig = otherCI->info.CI.contigID;

	    // RAW EDGES ONLY
	    assert(edge->flags.bits.isRaw);

	    // Only non-overlap  edges (i,j) s.t. i < j are interesting
	    if(/* (isOverlapCIEdgeT(edge)) || */
	       otherCI->info.CI.contigID == NULLINDEX ||              // not scaffolded/contiged
	       otherCI->info.CI.contigID <= thisCI->info.CI.contigID)       // not canonical
	      continue;

#ifndef TRY_IANS_EDGES
	    edgeOrient = GetEdgeOrientationWRT(edge, otherCI->id);
	    orient = ((edgeOrient == AB_BA || edgeOrient == AB_AB)?A_B:B_A);

#if 0
            fprintf(stderr,"* Edge (" F_CID "," F_CID ") %c dist: %g in scaffolds (" F_CID "," F_CID ") orient = %c\n",
		    thisCID, otherCID, edgeOrient, edge->distance.mean, 
		    thisCI->scaffoldID, otherCI->scaffoldID, orient);
#endif	    

	    /* Contigs are ALWAYS oriented A->B in scaffold, so we can reuse this code */

	    mCIok = CIinContigOffsetAndOrientation(graph,
                                                   otherFrag, // input
                                                   otherCI->id, // input
                                                   orient,    // input orientation of fragment in CI
                                                   &mciOffset,    // CI's offset within Scaffold 
                                                   &mciOrient);
	    if(!mCIok)
	      continue;

	    edgeOrient = GetEdgeOrientationWRT(edge, thisCI->id);
	    orient = ((edgeOrient == AB_BA || edgeOrient == AB_AB)?A_B:B_A);
	    CIok = CIinContigOffsetAndOrientation(graph,
                                                  frag,
                                                  thisCI->id, // input
                                                  orient,    // input orientation of fragment in CI
                                                  &ciOffset,    // CI's offset within Scaffold 
                                                  &ciOrient);
	    assert(CIok);

	    

	    /* Mate pairs must be oriented in opposite directions.  So, if they are oriented the same wrt
	       their own chunk, the chunks must be oriented opposite one another */
	    switch(ciOrient){
	      //
              case A_B:

                switch(mciOrient){
                  case A_B:
                    //           length - 5'             gap            length - 5'
                    //      |------------------------||---------------||-----------|
                    //  A --------------------------- B               B --------------------------- A
                    //    5'----->                                           <------5'
                    //      |-------------------------------------------------------|
                    //                             mate distance
                    //
                    contigEdgeOrient = AB_BA;  
                    break;
                  case B_A:
                    //           length - 5'             gap                5'
                    //      |------------------------||---------------||-----------|
                    //  A --------------------------- B               A --------------------------- B
                    //    5'----->                                           <------5'
                    //      |-------------------------------------------------------|
                    //                             mate distance
                    //
                    contigEdgeOrient = AB_AB; 
                    break;
                  default:
                    assert(0);
                    break;
                }
                break;
              case B_A:

                switch(mciOrient){
                  case A_B:
                    //                     5'             gap            length - 5'
                    //      |------------------------||---------------||-----------|
                    //  B --------------------------- A               B --------------------------- A
                    //    5'----->                                           <------5'
                    //      |-------------------------------------------------------|
                    //                             mate distance
                    //
                    contigEdgeOrient = BA_BA; 
                    break;
                  case B_A:
                    //                     5'             gap                5'
                    //      |------------------------||---------------||-----------|
                    //  B --------------------------- A               A --------------------------- B
                    //    5'----->                                           <------5'
                    //      |-------------------------------------------------------|
                    //                             mate/guide distance
                    //
                    contigEdgeOrient = BA_AB; 
                    break;
                  default:
                    assert(0);
                    break;
                }
                break;
              default:
                assert(0);
                break;
	    }
#if 0
            fprintf(stderr,"* CIEdge (" F_CID "," F_CID ") %c induced\n",
                    thisCID, otherCID, edge->orient);
            fprintf(stderr,"* (" F_CID "," F_CID ") %c ciedge:" F_CID " cifrag:" F_CID " otherFrag:" F_CID " mciOffset = %g mciOrient = %c  ciOffset = %g ciOrient = %c\n",
                    thisCI->contigID, otherCI->contigID, contigEdgeOrient,
                    (CDS_CID_t)GetVAIndex_CIEdgeT(graph->CIEdges, edge),
                    (frag?frag->iid:NULLINDEX),
                    (otherFrag?otherFrag->iid:NULLINDEX),
                    mciOffset.mean, mciOrient, ciOffset.mean, ciOrient);
#endif
	    distance.mean = edge->distance.mean - ciOffset.mean - mciOffset.mean;
	    // Since the two offsets and the dist are independent we SUM their variances
	    distance.variance = edge->distance.variance  + ciOffset.variance + mciOffset.variance;
	    contigEdge.orient = contigEdgeOrient;
	    contigEdge.distance = distance;
#else
            /* Ian's attempt to simplify, goes straight to fragments
               get edge between fragments, innie or outtie, distance
            */
            PopulateChunkEdgeBasics(graph,
                                    frag,
                                    thisCI,
                                    otherFrag,
                                    otherCI,
                                    GetDistT(graph->Dists, edge->distIndex),
                                    &contigEdge);
#endif
	    edgesSucceeded++;
	
	    contigEdge.referenceEdge = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->CIEdges, edge);
	    contigEdge.idA = thisCI->info.CI.contigID;
	    contigEdge.idB = otherCI->info.CI.contigID;
	    contigEdge.flags = edge->flags;
	    contigEdge.edgesContributing = edge->edgesContributing;
	    contigEdge.nextALink = contigEdge.nextBLink = NULLINDEX;
	    contigEdge.prevALink = contigEdge.prevBLink = NULLINDEX;
	    contigEdge.nextRawEdge = NULLINDEX;
	    contigEdge.minDistance = 0;
#if 0
	    contigEdge.maxDistance = contigEdge.fudgeDistance = 0;
#endif
	    if(isA){
	      contigEdge.fragA = edge->fragA;
	      contigEdge.fragB = edge->fragB;
	    }else{
	      contigEdge.fragA = edge->fragB;
	      contigEdge.fragB = edge->fragA;
	    }
	    {
	      CDS_CID_t edgeID = GetNumGraphEdges(graph->ContigGraph);
	      contigEdge.topLevelEdge = edgeID;
	      AppendGraphEdge(graph->ContigGraph, &contigEdge);
	      InsertGraphEdge(graph->ContigGraph, edgeID, FALSE);
	    }
	  }
	}
      }
    }
  }
  return TRUE;
}


#endif

/***************************************************************************/


void CreateInitialContig(ScaffoldGraphT *graph, CDS_CID_t cid){
  ChunkInstanceT *CI = GetGraphNode (graph->CIGraph, cid);
  ContigT contig = {0};

  CI->AEndNext = CI->BEndNext = NULLINDEX;
  contig = *CI;
  contig.type = CONTIG_CGW;
  CI->info.CI.contigID = contig.id = cid; // GetNumGraphNodes(graph->ContigGraph);
  contig.scaffoldID = NULLINDEX;
  contig.smoothExpectedCID = NULLINDEX;
  contig.numEssentialA = contig.numEssentialB = 0;
  contig.essentialEdgeA = contig.essentialEdgeB = NULLINDEX;
  contig.info.Contig.AEndCI = contig.info.Contig.BEndCI = cid;
  contig.info.Contig.numCI = 1;
  contig.indexInScaffold = NULLINDEX;
  contig.flags.bits.isCI = FALSE;
  contig.flags.bits.isContig = TRUE;
  contig.flags.bits.isChaff = CI->flags.bits.isChaff; // this property is inherited
  if(GlobalData->debugLevel > 0 &&
     contig.flags.bits.isChaff){
    fprintf(stderr,"* Contig " F_CID " is CHAFF\n", contig.id);
  }
  contig.microhetScore = NULLINDEX;
  contig.edgeHead = NULLINDEX;
  
  SetNodeCGW_T(graph->ContigGraph->nodes, cid, &contig);
}


//#define DEBUG_CONTIG
void CreateInitialContigEdges(ScaffoldGraphT *graph){
  CIEdgeT *edge, *redge;
  int i;
  int deletedSkipped = 0;
  int shouldHaveSkipped = 0;
  int actuallySkipped = 0;

  for(i = 0; i < GetNumGraphEdges(graph->CIGraph); i++){
    CIEdgeT newEdge, rawEdge;
    CDS_CID_t newCIEdge = GetNumGraphEdges(graph->ContigGraph);
    CDS_CID_t nextCIEdge = newCIEdge;
    edge = GetGraphEdge(graph->CIGraph, i);

    // Skip deleted edges
    if(edge->flags.bits.isDeleted){
      deletedSkipped++;
      continue;
    }
    newEdge = *edge;
    newEdge.prevALink = newEdge.prevBLink = NULLINDEX;
    newEdge.nextALink = newEdge.nextBLink = NULLINDEX;

    if(!newEdge.flags.bits.isRaw){
      newEdge.topLevelEdge = newCIEdge; // selfreference
      newEdge.nextRawEdge = ++nextCIEdge;
#ifdef DEBUG_CONTIG
      fprintf(stderr,"* Converting %s edge (" F_CID "," F_CID ") to edge (" F_CID "," F_CID ") edgeID " F_CID "\n",
              (edge->flags.bits.isRaw?" Raw":"Merged"),
              edge->idA, edge->idB,
              newEdge.idA, newEdge.idB, newCIEdge);
      fprintf(stderr,"* linked to raw edge " F_CID "\n",
	      newEdge.nextRawEdge);
#endif
      AppendGraphEdge(graph->ContigGraph, &newEdge);
      assert(newCIEdge == GetNumGraphEdges(graph->ContigGraph) -1);
      InsertGraphEdgeInList(graph->ContigGraph, newCIEdge, newEdge.idA, FALSE);
      InsertGraphEdgeInList(graph->ContigGraph, newCIEdge, newEdge.idB, FALSE);
      redge = edge;
#ifdef DEBUG_CONTIG
      fprintf(stderr,"* Converting edge (" F_CID "," F_CID ") to edge (" F_CID "," F_CID ") edgeID " F_CID "\n",
              edge->idA, edge->idB,
              newEdge.idA, newEdge.idB, newCIEdge);
#endif
      while(redge->nextRawEdge != NULLINDEX){
	actuallySkipped++;
	redge = GetGraphEdge(graph->CIGraph, redge->nextRawEdge);
	rawEdge = *redge;
	rawEdge.idA = newEdge.idA;
	rawEdge.idB = newEdge.idB;
	rawEdge.prevALink = rawEdge.prevBLink = NULLINDEX;
	rawEdge.nextALink = rawEdge.nextBLink = NULLINDEX;
	rawEdge.topLevelEdge = newCIEdge;
	rawEdge.referenceEdge = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->CIGraph->edges, edge);
	assert(GetNumGraphEdges(graph->ContigGraph) == nextCIEdge);
	if(rawEdge.nextRawEdge != NULLINDEX)
	  rawEdge.nextRawEdge = ++nextCIEdge;
#ifdef DEBUG_CONTIG
	fprintf(stderr,"* Adding raw edge " F_CID " linked  to raw edge " F_CID "\n",
                (CDS_CID_t) GetNumGraphEdges(graph->ContigGraph),
                rawEdge.nextRawEdge);
#endif
	AppendGraphEdge(graph->ContigGraph, &rawEdge);
      }
    }else{
      // If this isn't a top-level edge, skip it
      if(edge->topLevelEdge != (CDS_CID_t)GetVAIndex_CIEdgeT(graph->CIGraph->edges, edge)){
#ifdef DEBUG_CONTIG
        EdgeCGW_T *tle = GetGraphEdge(graph->CIGraph, edge->topLevelEdge);
	fprintf(stderr,"* !!!!! Skipping raw edge that is not topLevel (" F_CID "," F_CID ") edgeID " F_CID " \n",
		edge->idA, edge->idB,
		(CDS_CID_t) GetVAIndex_CIEdgeT(graph->CIGraph->edges, edge));
        PrintGraphEdge(stream,graph->CIGraph,"raw edge ",edge, edge->idA);
        PrintGraphEdge(stream,graph->CIGraph,"top edge ",tle, tle->idA);
#endif
	shouldHaveSkipped++;
	continue;
      }
      newEdge.referenceEdge = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->CIGraph->edges, edge);
      newEdge.topLevelEdge = newCIEdge; // selfreference
      AppendGraphEdge(graph->ContigGraph, &newEdge);
      assert(newCIEdge == GetNumGraphEdges(graph->ContigGraph) -1);
	

      InsertGraphEdgeInList(graph->ContigGraph, newCIEdge, newEdge.idA, FALSE);
      InsertGraphEdgeInList(graph->ContigGraph, newCIEdge, newEdge.idB, FALSE);
    }
  }
  fprintf(stderr,"* Skipped %d raw edges. %d deleted..should have been %d (%d <==> %d)\n",
	  actuallySkipped, deletedSkipped, shouldHaveSkipped, 
	  (int) GetNumGraphEdges(graph->ContigGraph),
          (int) GetNumGraphEdges(graph->CIGraph)  );

}

int BuildInitialContigs(ScaffoldGraphT *graph){
  CDS_CID_t cid;
  GraphNodeIterator CIs;
  NodeCGW_T *CI;
  
  /* Resize the ContigGraph to the same size as the CI Graph */

  fprintf(stderr,"* Allocating Contig Graph with %d nodes and %d edges\n",
	  (int) GetNumGraphNodes(graph->CIGraph),
          (int) GetNumGraphEdges(graph->CIGraph));

  fflush(NULL);

  // Would be nice if there was a Var_Array function to realloc a Var_Array, without
  // adjusting the number of active elements.
  DeleteVA_EdgeCGW_T(graph->ContigGraph->edges);
  DeleteVA_NodeCGW_T(graph->ContigGraph->nodes);
  graph->ContigGraph->edges = CreateVA_EdgeCGW_T(GetNumGraphEdges(graph->CIGraph));
  graph->ContigGraph->nodes = CreateVA_NodeCGW_T(GetNumGraphNodes(graph->CIGraph));

  EnableRange_VA(graph->ContigGraph->nodes, GetNumGraphNodes(graph->CIGraph));

  for(cid = 0; cid < GetNumGraphNodes(graph->ContigGraph); cid++){
    CI = GetGraphNode(graph->ContigGraph, cid);
    CI->flags.all = 0;
    CI->flags.bits.isContig = TRUE;
    CI->flags.bits.isDead = TRUE;
  }

  InitGraphNodeIterator(&CIs, graph->CIGraph, GRAPH_NODE_DEFAULT);
  while((CI = NextGraphNodeIterator(&CIs)) != NULL){
    /* For each chunk */
    //  graph->firstContig = GetNumChunkInstanceTs(graph->ChunkInstances);

    CreateInitialContig(graph, CI->id);
  }
  graph->numContigs = GetNumGraphNodes(graph->ContigGraph);

  // Create the Contig edges from the CIEdges
  CreateInitialContigEdges(graph);

  return TRUE;
}

int GetConsensus(GraphCGW_T *graph, CDS_CID_t CIindex,
                 VA_TYPE(char) *consensusVA, VA_TYPE(char) *qualityVA){
  // Return value is length of unitig or contig  sequence/quality (-1 if failure)
  ChunkInstanceT *CI = GetGraphNode(graph, CIindex);
  MultiAlignT *MA = NULL;

  ResetVA_char(consensusVA);
  ResetVA_char(qualityVA);
  if(CI->flags.bits.isCI){
    // Get it from the store of Unitig multi alignments
    MA = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, CIindex, TRUE);
  }else if(CI->flags.bits.isContig){// Get it from the store of Contig multi alignments
    assert(graph->type == CONTIG_GRAPH);
    MA = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, CIindex, FALSE);
  }else assert(0);

  GetMultiAlignUngappedConsensus(MA, consensusVA, qualityVA);

  return GetNumchars(consensusVA);
}


void SetCIScaffoldIds(ChunkInstanceT *CI, CDS_CID_t scaffoldID){
  // Set the scaffold ID of this CI
  CI->scaffoldID = scaffoldID;
  if(CI->flags.bits.isChaff){ // This can only happen once
    MultiAlignT *ma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, CI->id, TRUE);
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags,  (int)GetIntMultiPos(ma->f_list,0)->sourceInt);
    assert(frag->flags.bits.isSingleton);
    frag->flags.bits.isChaff = FALSE;
    CI->flags.bits.isChaff = FALSE;
    if(GlobalData->debugLevel > 0)
      fprintf(stderr,"* SetCIScaffoldIDs ci " F_CID " and frag " F_CID " are NOT chaff\n",
	      CI->id, frag->iid);
  }
}

void SetCIContigIds(ChunkInstanceT *CI, CDS_CID_t contigID){
  // Set the contig ID of this CI
  CI->info.CI.contigID = contigID;
}

void SetContigScaffoldIds(ContigT *contig, CDS_CID_t scaffoldID){
  ContigTIterator CIs;
  ChunkInstanceT *CI;

  // Set the scaffold ID of this CI
  contig->scaffoldID = scaffoldID;

  // Set the isUnique bit
  contig->flags.bits.isUnique = TRUE;
  contig->flags.bits.isChaff = FALSE;
  if(GlobalData->debugLevel > 0)
    fprintf(stderr,"* SetContigScaffoldIDs contig " F_CID " is NOT chaff\n",
            contig->id);

  // Set the scaffold ID of all of the contained CIs
  InitContigTIterator(ScaffoldGraph, contig->id, TRUE, FALSE, &CIs);
  while((CI = NextContigTIterator(&CIs)) != NULL){
    assert(CI->flags.bits.isCI);
    SetCIScaffoldIds(CI,scaffoldID);
  }
}


#if 0
/***************************************************************************/
void UpdateContigSimCoordinates(ContigT *contig){
  ContigTIterator CIs;
  ChunkInstanceT *CI;
  CDS_CID_t maxID = NULLINDEX, minID = NULLINDEX;
  int minEnd, maxEnd;
  CDS_COORD_t min = CDS_COORD_MAX, max = CDS_COORD_MIN;
  int invalid = TRUE;

  contig->aEndCoord = contig->bEndCoord = -1;
  contig->flags.bits.cgbType = RR_CGBTYPE;


  InitContigTIterator(ScaffoldGraph, contig->id, TRUE, FALSE, &CIs);
  while((CI = NextContigTIterator(&CIs)) != NULL){
    if(CI->flags.bits.cgbType == UU_CGBTYPE){
      if(CI->offsetAEnd.mean < CI->offsetBEnd.mean){
	if(min > CI->offsetAEnd.mean){
	  minID = CI->id;
	  minEnd = A_END;
	  min = CI->offsetAEnd.mean;
	}
	if(max < CI->offsetBEnd.mean){
	  maxID = CI->id;
	  maxEnd = B_END;
	  max = CI->offsetBEnd.mean;
	}
      }else{
	if(max < CI->offsetAEnd.mean){
	  maxID = CI->id;
	  maxEnd = A_END;
	  max = CI->offsetAEnd.mean;
	}
	if(min > CI->offsetBEnd.mean){
	  minID = CI->id;
	  min = CI->offsetBEnd.mean;
	  minEnd = B_END;
	}
      }
    }else{
      invalid = TRUE;
    }
  }
  
  contig->flags.bits.cgbType = (invalid? RR_CGBTYPE: UU_CGBTYPE);

  if(minID != NULLINDEX){
    NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, minID);
    if(minEnd == A_END)
      contig->aEndCoord = unitig->aEndCoord;
    else
      contig->aEndCoord = unitig->bEndCoord;
  }

  if(maxID != NULLINDEX){
    NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, maxID);
    if(maxEnd == A_END)
      contig->bEndCoord = unitig->aEndCoord;
    else
      contig->bEndCoord = unitig->bEndCoord;
  }
}
#endif
/****************************************************************************/
void UpdateContigSimCoordinates(NodeCGW_T *contig){
  ContigTIterator CIs;
  ChunkInstanceT *CI;
  CDS_CID_t maxID = NULLINDEX, minID = NULLINDEX;
  int minEnd = NO_END, maxEnd = NO_END;
  CDS_COORD_t min = CDS_COORD_MAX, max = CDS_COORD_MIN;
  int invalid = FALSE;
  int seenValid = FALSE;

  contig->aEndCoord = contig->bEndCoord = -1;
  contig->flags.bits.cgbType = RR_CGBTYPE;

  //  fprintf(stderr,"* Computing Sim coordinates for contig " F_CID "\n",
  //	  contig->id);

  InitContigTIterator(ScaffoldGraph, contig->id, TRUE, FALSE, &CIs);
  while((CI = NextContigTIterator(&CIs)) != NULL){
    if(CI->flags.bits.cgbType != UU_CGBTYPE){
      //      fprintf(stderr,"* CI " F_CID " is invalid ==> contig is invalid\n",
      //	      CI->id);
      invalid = TRUE;
      if(contig->info.Contig.numCI > 1)
	continue;
    }
    seenValid = TRUE;
    /*      fprintf(stderr,"* CI " F_CID " offset[%f,%f] sim[" F_COORD "," F_COORD "]\n",
            CI->id, 
            CI->offsetAEnd.mean, CI->offsetBEnd.mean,
            CI->aEndCoord, CI->bEndCoord); */
    assert(CI->offsetAEnd.mean != CI->offsetBEnd.mean);

    if(CI->offsetAEnd.mean < CI->offsetBEnd.mean){
      if(min > CI->offsetAEnd.mean){
        minID = CI->id;
        minEnd = A_END;
        min = CI->offsetAEnd.mean;
      }
      if(max < CI->offsetBEnd.mean){
        maxID = CI->id;
        maxEnd = B_END;
        max = CI->offsetBEnd.mean;
      }
    }else{
      if(max < CI->offsetAEnd.mean){
        maxID = CI->id;
        maxEnd = A_END;
        max = CI->offsetAEnd.mean;
      }
      if(min > CI->offsetBEnd.mean){
        minID = CI->id;
        min = CI->offsetBEnd.mean;
        minEnd = B_END;
      }
    }
  }
  
  if(seenValid == FALSE){ // This contig has NO unique unitigs...pick the first unitig and use its coordinates
    invalid = TRUE;
    InitContigTIterator(ScaffoldGraph, contig->id, TRUE, FALSE, &CIs);
    while((CI = NextContigTIterator(&CIs)) != NULL){
      assert(CI->offsetAEnd.mean != CI->offsetBEnd.mean);

      if(CI->offsetAEnd.mean < CI->offsetBEnd.mean){
	if(min > CI->offsetAEnd.mean){
	  minID = CI->id;
	  minEnd = A_END;
	  min = CI->offsetAEnd.mean;
	}
	if(max < CI->offsetBEnd.mean){
	  maxID = CI->id;
	  maxEnd = B_END;
	  max = CI->offsetBEnd.mean;
	}
      }else{
	if(max < CI->offsetAEnd.mean){
	  maxID = CI->id;
	  maxEnd = A_END;
	  max = CI->offsetAEnd.mean;
	}
	if(min > CI->offsetBEnd.mean){
	  minID = CI->id;
	  min = CI->offsetBEnd.mean;
	  minEnd = B_END;
	}
      }
      break;
    }

  }
  /*  fprintf(stderr,"* MinID = " F_CID " MinEnd = %c   MaxID = " F_CID " MaxEnd = %c\n",
      minID, (minEnd == A_END? 'A':'B'),
      maxID, (maxEnd == A_END? 'A':'B')); */

  contig->flags.bits.cgbType = (invalid? RR_CGBTYPE: UU_CGBTYPE);
  assert(maxID != NULLINDEX && minID != NULLINDEX);

  {
    NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, minID);
    if(minEnd == A_END)
      contig->aEndCoord = unitig->aEndCoord;
    else
      contig->aEndCoord = unitig->bEndCoord;
  }

  {
    NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, maxID);
    if(maxEnd == A_END)
      contig->bEndCoord = unitig->aEndCoord;
    else
      contig->bEndCoord = unitig->bEndCoord;
  }


  /*  fprintf(stderr,"* MinID = " F_CID " MinEnd = %c   MaxID = " F_CID " MaxEnd = %c [" F_COORD "," F_COORD "]\n",
      minID, (minEnd == A_END? 'A':'B'),
      maxID, (maxEnd == A_END? 'A':'B'),
      contig->aEndCoord, contig->bEndCoord); 
      fprintf(stderr,"* Contig " F_CID " computed sim coordinates [" F_COORD "," F_COORD "]\n",
      contig->id,
      contig->aEndCoord, contig->bEndCoord);
  */
}


/***************************************************************************/
void UpdateScaffoldSimCoordinates(NodeCGW_T *scaffold){
  CIScaffoldTIterator Contigs;
  ContigT *CI;
  CDS_CID_t maxID = NULLINDEX, minID = NULLINDEX;
  int minEnd = NO_END, maxEnd = NO_END;
  CDS_COORD_t min = CDS_COORD_MAX, max = CDS_COORD_MIN;
  int invalid = TRUE;

  scaffold->aEndCoord = scaffold->bEndCoord = -1;
  scaffold->flags.bits.cgbType = RR_CGBTYPE;


  InitCIScaffoldTIterator(ScaffoldGraph,
                          scaffold,
                          TRUE,
                          FALSE,
                          &Contigs);
  while((CI = NextCIScaffoldTIterator(&Contigs)) != NULL){
    if(CI->flags.bits.cgbType == UU_CGBTYPE){
      if(CI->offsetAEnd.mean < CI->offsetBEnd.mean){
	if(min > CI->offsetAEnd.mean){
	  minID = CI->id;
	  minEnd = A_END;
	  min = CI->offsetAEnd.mean;
	}
	if(max < CI->offsetBEnd.mean){
	  maxID = CI->id;
	  maxEnd = B_END;
	  max = CI->offsetBEnd.mean;
	}
      }else{
	if(max < CI->offsetAEnd.mean){
	  maxID = CI->id;
	  maxEnd = A_END;
	  max = CI->offsetAEnd.mean;
	}
	if(min > CI->offsetBEnd.mean){
	  minID = CI->id;
	  min = CI->offsetBEnd.mean;
	  minEnd = B_END;
	}
      }
    }else{
      invalid = TRUE;
    }
  }
  
  scaffold->flags.bits.cgbType = (invalid? RR_CGBTYPE: UU_CGBTYPE);

  if(minID != NULLINDEX){
    NodeCGW_T *contig = GetGraphNode(ScaffoldGraph->RezGraph, minID);
    if(minEnd == A_END)
      scaffold->aEndCoord = contig->aEndCoord;
    else
      scaffold->aEndCoord = contig->bEndCoord;
  }

  if(maxID != NULLINDEX){
    NodeCGW_T *contig = GetGraphNode(ScaffoldGraph->RezGraph, maxID);
    if(maxEnd == A_END)
      contig->bEndCoord = contig->aEndCoord;
    else
      contig->bEndCoord = contig->bEndCoord;
  }
}


void CheckAllContigFragments(void){
  GraphNodeIterator contigs;
  NodeCGW_T *contig;

  InitGraphNodeIterator(&contigs, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((contig = NextGraphNodeIterator(&contigs)) != NULL){
    MultiAlignT *ma  = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);
    int i;
    if(!ma){
      fprintf(stderr,"*CheckAllContigFragments -- Contig " F_CID " is missing\n", contig->id);
      continue;
    }
    for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++){
      IntMultiPos *mp = GetIntMultiPos(ma->f_list,i);
      CDS_CID_t fragID = (CDS_CID_t)mp->sourceInt;
      CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags,fragID);
      assert(frag->contigID == contig->id);
    }
  }
}

CDS_CID_t GetOriginalContigID(CDS_CID_t contigID){
  NodeCGW_T *contig = GetGraphNode(ScaffoldGraph->ContigGraph, contigID);
  NodeCGW_T *ci;

  if(contig->info.Contig.numCI > 1)
    return contig->id;

  // Get the CI in this contig
  ci = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);

  // If this is not a surrogate, we're done
  if(ci->type != RESOLVEDREPEATCHUNK_CGW)
    return contig->id;

  // Get the parent of the surrogate
  ci = GetGraphNode(ScaffoldGraph->CIGraph, ci->info.CI.baseID);

  // Return its contigID
  return ci->info.CI.contigID;

}


int IsDefinitelyUniqueContig(ContigT *contig){
  NodeCGW_T *ci;

  /*
    fprintf(stderr,"* IsDefinitelyUniqueContig " F_CID " numCI:%d\n",
    contig->id, contig->info.Contig.numCI);
  */
  // If it has been contigged already, it is unique
  if(contig->info.Contig.numCI > 1)
    return TRUE;

  // Get the CI in this contig
  ci = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);
  /*
    fprintf(stderr,"* type = %d  coverage = %d  (%d)\n",
    ci->type, ci->info.CI.coverageStat,
    GlobalData->cgbDefinitelyUniqueCutoff);
  */

  // when the flag says we are unique, we always return true
  if (ci->unique_rept == AS_FORCED_UNIQUE) {
    return TRUE;
  }

  // when the flag says we are repeat, we always return false
  if (ci->unique_rept == AS_FORCED_REPEAT) {
    return FALSE;
  }
  
  // If this is not a surrogate, we're done
  return( ci->info.CI.coverageStat > GlobalData->cgbDefinitelyUniqueCutoff);

}


int IsShakyContigAtScaffoldEnd(ContigT *contig){
  
  if(contig->numEssentialB <= 1 &&
     contig->numEssentialA <= 1){
    return FALSE;
  }

#if 0
  fprintf(stderr,"* IsShakyContigAtScaffoldEnd " F_CID " numA:%d numB:%d %s\n",
          contig->id, 
	  contig->numEssentialA,
	  contig->numEssentialB,
	  (IsDefinitelyUniqueContig(contig)?"unique":"shaky"));
#endif

  return(!IsDefinitelyUniqueContig(contig));
}

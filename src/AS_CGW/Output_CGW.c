
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
static char CM_ID[] = "$Id: Output_CGW.c,v 1.7 2005-08-02 02:40:51 gdenisov Exp $";

#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_PER_fragStore.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "ChiSquareTest_CGW.h"

VA_DEF(IntMate_Pairs)
static VA_TYPE(IntMate_Pairs) *JumpList = NULL;

     /*********/

/* This Routine outputs the MateDist messages and also set the 
   outMateStat field for the IntAugFrag messages */

void OutputMateDists(ScaffoldGraphT *graph){
  int                   i; 
  GenericMesg		pmesg;
  IntMateDistMesg	imd;
  DistT			*dptr;

  pmesg.m = &imd;
  pmesg.t = MESG_IMD;

  assert(graph->doRezOnContigs);

  for(i = 1; i < GetNumDistTs(graph->Dists); i++){
    dptr = GetDistT(graph->Dists, i);

#define ORIGINALFLAVOR 0

#if ORIGINALFLAVOR != 0
    if (dptr->numSamples == 0)
      continue;
    imd.num_buckets = dptr->bnum;
    imd.min = dptr->min;
    imd.max = dptr->max;
    imd.stddev = dptr->sigma;
    imd.mean = dptr->mu;
#else
    if (dptr->numSamples > 30 ){
      imd.stddev = dptr->sigma;
      imd.mean = dptr->mu;
      imd.num_buckets = dptr->bnum;
      imd.min = dptr->min;
      imd.max = dptr->max;
    } else {
      imd.num_buckets = 0;
      imd.min = CDS_COORD_MIN;
      imd.max = CDS_COORD_MAX;
      imd.stddev = dptr->stddev;
      imd.mean = dptr->mean;
    }
#endif
    imd.refines = i;
    imd.histogram = dptr->histogram;
    (GlobalData->writer)(GlobalData->outfp,&pmesg);
    free(dptr->histogram);
    dptr->histogram = NULL;
  }
  fflush(NULL);
}


/* Must be called after OutputMateDists */
void OutputFrags(ScaffoldGraphT *graph){
  CDS_CID_t		i;
  int numFrags = GetNumInfoByIIDs(graph->iidToFragIndex);
  GenericMesg		pmesg;
  IntAugFragMesg af_mesg;
  int goodMates = 0;
  //MateStatusType status;
  MateStatType status;
  InfoByIID *info = GetInfoByIID(graph->iidToFragIndex,0);
#if 0
  ReadStructp fsread = new_ReadStruct();
  uint32 clr_bgn_orig, clr_end_orig, clr_bgn_latest, clr_end_latest;
#endif
  
  pmesg.m = &af_mesg;
  pmesg.t = MESG_IAF;
  
  // Output fragments in iid order
  //
  for(i = 0; i < numFrags; i++, info++){
    CDS_CID_t fragID;
    CIFragT *cifrag;

    if(!info->set)
      continue;

    if(((i+1) % 1000000) == 0){
      fprintf(stderr,"* Outputing fragment " F_CID "\n",i);
      fflush(stderr);
    }
    fragID = info->fragIndex;
    status = BAD_MATE;

    cifrag = GetCIFragT(graph->CIFrags,fragID);

    assert(cifrag->iid == i);

    af_mesg.iaccession = cifrag->iid;
    af_mesg.type = (FragType)cifrag->type;
    af_mesg.chaff = cifrag->flags.bits.isChaff;
    switch(cifrag->flags.bits.edgeStatus){
    case INVALID_EDGE_STATUS:
      status = NO_MATE;
      break;
    case TRUSTED_EDGE_STATUS:
    case TENTATIVE_TRUSTED_EDGE_STATUS:
      status = GOOD_MATE;
      goodMates++;
      break;
    case  UNTRUSTED_EDGE_STATUS:
    case  TENTATIVE_UNTRUSTED_EDGE_STATUS:
      status = BAD_MATE;
      break;;
    case LARGE_VARIANCE_EDGE_STATUS:
    case INTER_SCAFFOLD_EDGE_STATUS:
    case UNKNOWN_EDGE_STATUS:
      status = UNRESOLVED_MATE;
      break;
    default:
      assert(0);
    }
    af_mesg.mate_status = status;
    af_mesg.chimeric = 0;

#if 0
	getFragStore( ScaffoldGraph->fragStore, i, FRAG_S_ALL, fsread);
	getClearRegion_ReadStruct( fsread, &clr_bgn_orig, &clr_end_orig, READSTRUCT_ORIGINAL);
	getClearRegion_ReadStruct( fsread, &clr_bgn_latest, &clr_end_latest, READSTRUCT_LATEST);

	if ( (clr_bgn_orig == clr_bgn_latest) && (clr_end_orig == clr_end_latest)) // clr range is unchanged
	{
	  af_mesg.clear_rng.bgn = -1;
	  af_mesg.clear_rng.end = -1;
	}
	else
	{
	  af_mesg.clear_rng.bgn = clr_bgn_latest;
	  af_mesg.clear_rng.end = clr_end_latest;
	}
#endif

	af_mesg.clear_rng.bgn = -1;
	af_mesg.clear_rng.end = -1;
	  
    (GlobalData->writer)(GlobalData->outfp,&pmesg);
  }
  fprintf(GlobalData->stderrc,"* Saw %d good mates\n", goodMates);
  fflush(NULL);

}


/* This routine not only outputs the PCM messages, but also sets up
   some values used for ICL & IMD messages.  Thus it must be called
   before OutputConigLinks and OutputMateDists */

void MarkContigEdges(void){
  assert(ScaffoldGraph->doRezOnContigs);
  
  fprintf(GlobalData->stderrc,"* MarkContigEdges\n");
  /* This block initializes some values for the mate distance distribution */
  {
    /* hopefully this is the max positive value for int32 */
    DistT *dptr;
    int i, dbound;
    
    dbound = (int) GetNumDistTs(ScaffoldGraph->Dists);
    for (i=0; i < dbound; ++i) {
      dptr = GetDistT(ScaffoldGraph->Dists,i);
      dptr->min = CDS_COORD_MAX;
      dptr->max = CDS_COORD_MIN;
      dptr->numSamples = 0;
      dptr->mu = dptr->sigma = 0.0;
      dptr->lower = dptr->mean - CGW_CUTOFF*dptr->stddev;
      dptr->upper = dptr->mean + CGW_CUTOFF*dptr->stddev;
    }
  } 
  
  fprintf(GlobalData->stderrc,"* Intra-scaffold inter-contig edges\n");
  // Mark the trustedness of the intra-scaffold, inter-contig edges
  {
    CIScaffoldT *scaffold;
    GraphNodeIterator scaffolds;
    int numContigs = 0, numEdges = 0, numScaffolds = 0;
    
    InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
    while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
      if(scaffold->type != REAL_SCAFFOLD)
	continue;
      MarkInternalEdgeStatus(ScaffoldGraph, scaffold, PAIRWISECHI2THRESHOLD_CGW,
			     100000000000.0, TRUE, TRUE, 0, FALSE);
    }
    
    InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
    while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
      ContigT *contig;
      CIScaffoldTIterator Contigs;
      
      numScaffolds++;
      if(numScaffolds % 100 == 0)
	fprintf(GlobalData->stderrc,"* Scaffold %d (" F_CID ") \n", numScaffolds, scaffold->id);
      
      /* Iterate over all contigs in scaffolds */
      
      InitCIScaffoldTIterator(ScaffoldGraph, scaffold,TRUE, FALSE, &Contigs);
      while((contig = NextCIScaffoldTIterator(&Contigs)) != NULL){
	GraphEdgeIterator edges;
	EdgeCGW_T *edge;
        
	numContigs++;
	InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, contig->id, ALL_END, ALL_EDGES, GRAPH_EDGE_RAW_ONLY , &edges);
	while((edge = NextGraphEdgeIterator(&edges)) != NULL){
          ContigT *mcontig;
          
	  assert(edge->flags.bits.isRaw);
	  if((edge->idA != contig->id) || isSingletonOverlapEdge(edge))
	    continue;
          
	  numEdges++;
          
	  if(numEdges % 1000 == 0)
	    fprintf(GlobalData->stderrc,"* Marked %d edges for %d contigs (contig " F_CID ")\n",
		    numEdges, numContigs, contig->id);
          
	  mcontig = GetGraphNode(ScaffoldGraph->ContigGraph, edge->idB);
          
	  if(contig->scaffoldID != mcontig->scaffoldID)
	    SetEdgeStatus(ScaffoldGraph->ContigGraph, edge, INTER_SCAFFOLD_EDGE_STATUS);
	  PropagateEdgeStatusToFrag(ScaffoldGraph->ContigGraph, edge);
	}
      }
    }
  }
  fprintf(GlobalData->stderrc,"* Calculating Intra-contig mate link stats *\n");
  fflush(NULL);
  
  // Now collect mate link stats for mate link 
  
  /* Accumulate Mate Link Statistics */
  
  fprintf(GlobalData->stderrc,"* Mate pair statistics\n");
  fflush(NULL);
  
  ComputeMatePairStatisticsRestricted( CONTIG_OPERATIONS, CDS_CID_MAX /* minSamplesForOverride */, "MarkContigEdges");
}



/****************************************************************************/
void OutputContigsFromMultiAligns(){
  GenericMesg		pmesg;
  IntConConMesg		icm_mesg;
  IntUnitigPos		*uptr;
  GraphCGW_T *graph = ScaffoldGraph->ContigGraph;
  GraphNodeIterator     nodes;
  ContigT		*ctg;
  MultiAlignT *ma = CreateEmptyMultiAlignT();
  int32 ubufSize = 100;
  
  pmesg.m = &icm_mesg;
  pmesg.t = MESG_ICM;
  
  icm_mesg.unitigs = (IntUnitigPos *) safe_malloc(ubufSize*sizeof(IntUnitigPos));
  
  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);
  /* 1st get min and max values */
  while((ctg = NextGraphNodeIterator(&nodes)) != NULL){
    CDS_IID_t i;
    
    if(ctg->flags.bits.isChaff){
      //      fprintf(GlobalData->stderrc,"* # Contig " F_CID " is CHAFF\n", ctg->id);
      continue;
    }
    
    {
      CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, ctg->scaffoldID);
      CDS_IID_t numFrag;
      CDS_IID_t numUnitig;
      CDS_IID_t * tmpSource;
      IntMultiPos *mp;
      IntUnitigPos *up;
      //    MultiAlignT *ma = GetMultiAlignInStore(graph->maStore, ctg->id);
      ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ma, ctg->id, FALSE);
      numFrag = GetNumIntMultiPoss(ma->f_list);
      mp = GetIntMultiPos(ma->f_list,0);
      numUnitig = GetNumIntUnitigPoss(ma->u_list);
      up = GetIntUnitigPos(ma->u_list,0);
      
      tmpSource = safe_malloc((GetNumIntMultiPoss(ma->f_list) + 1) * sizeof(CDS_IID_t));
      
      if(numUnitig >= ubufSize){
        ubufSize = numUnitig * 2;
        icm_mesg.unitigs = (IntUnitigPos *) safe_realloc(icm_mesg.unitigs, ubufSize*sizeof(IntUnitigPos));
      }
      uptr = icm_mesg.unitigs;
      for(i = 0; i < numUnitig; i++){
        IntUnitigPos *iup = up + i;
        NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, iup->ident);
        if(unitig->type == DISCRIMINATORUNIQUECHUNK_CGW){
          uptr[i].type = AS_UNIQUE_UNITIG;
        }else{
          if(unitig->scaffoldID != NULLINDEX){
            if(!unitig->flags.bits.isSurrogate){
              uptr[i].type = AS_ROCK_UNITIG;
            }else  if(unitig->flags.bits.isStoneSurrogate){
              uptr[i].type = AS_STONE_UNITIG;
            }else{
              uptr[i].type = AS_PEBBLE_UNITIG;
            }
          }else{
            uptr[i].type = AS_SINGLE_UNITIG;
          }
        }
        uptr[i].position = iup->position;
        uptr[i].delta_length = iup->delta_length;
        uptr[i].delta = iup->delta;
        if(unitig->type == RESOLVEDREPEATCHUNK_CGW){
          iup->ident = unitig->info.CI.baseID; // map back to the parent of this instance
        }
        uptr[i].ident = iup->ident;
      }
      // Null out the source field
      for(i = 0; i < numFrag; i++){
        IntMultiPos *mp_i = GetIntMultiPos(ma->f_list,i);
        CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags,
                                   (CDS_CID_t)mp_i->source);
        tmpSource[i] = (CDS_IID_t) mp_i->source;
#ifdef DEBUG_DATA
	// This is only turned on if GlobalData->debugLevel > 0  see Input_CGW.c
	mp_i->source = Getchar(ScaffoldGraph->SourceFields, frag->source);
	//      fprintf(GlobalData->stderrc,"* " F_IID "  Frag " F_IID " has source " F_IID " %s\n", i,tmpSource[i] ,(CDS_CID_t) frag->source, (CDS_CID_t) mp_i->source);
#else
	mp_i->source = NULL;
#endif
      }
      ctg->outputID = ctg->id ;  // cid++;
      icm_mesg.placed = (scaffold && (scaffold->type == REAL_SCAFFOLD)?AS_PLACED:AS_UNPLACED);
      icm_mesg.iaccession = ctg->outputID;
      icm_mesg.forced = 0;
      icm_mesg.num_pieces = numFrag;
      icm_mesg.pieces = mp;
      icm_mesg.num_unitigs = numUnitig;
      icm_mesg.length = GetMultiAlignLength(ma);
      icm_mesg.num_vars = 0;
      icm_mesg.v_list   = NULL;
      if(icm_mesg.num_unitigs > 1){
        icm_mesg.consensus = ""; // Getchar(ma->consensus,0);
        icm_mesg.quality = ""; // Getchar(ma->quality,0);
      }else{
        icm_mesg.consensus = Getchar(ma->consensus,0);
        icm_mesg.quality = Getchar(ma->quality,0);
      }
      
      if(icm_mesg.num_unitigs > 1){
        assert(ctg->scaffoldID != NULLINDEX);
        (GlobalData->writer)(GlobalData->outfp1,&pmesg);
      }else{
        if(ctg->scaffoldID == NULLINDEX) {// contig is not placed
          GenericMesg		mesg;
          IntDegenerateScaffoldMesg dsc_mesg;
          NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, ctg->info.Contig.AEndCI);
          
          assert(unitig != NULL);
          if(unitig->info.CI.numInstances == 0){ // If this unitig has been placed as a surrogate, don't output contig
            dsc_mesg.icontig = ctg->id;
            mesg.m = &dsc_mesg;
            mesg.t = MESG_IDS;
            
            (GlobalData->writer)(GlobalData->outfp,&pmesg); // write the contig
            (GlobalData->writer)(GlobalData->outfp,&mesg);  // write the associated degenerate scaffold
          }else{
            // do nothing. The unitig in this contig appears as a surrogate elsewhere in the assembly
          }
        }else{ // Contig is placed
	  (GlobalData->writer)(GlobalData->outfp,&pmesg); // write the contig
        }     
      }
      
      // Restore the source values
      for(i = 0; i < numFrag; i++){
        IntMultiPos *mp_i = GetIntMultiPos(ma->f_list,i);
        mp_i->source = (char *)tmpSource[i];
      }
      //    UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ctg->id, FALSE);
      free(tmpSource);
    }
  }
  free(icm_mesg.unitigs);
  DeleteMultiAlignT(ma);
  fflush(NULL);
}

int SurrogatedSingleUnitigContig( NodeCGW_T* contig)
{
  if(contig->info.Contig.numCI > 1)
  {
	return 0;
  }
  else
  {
	if(contig->scaffoldID == NULLINDEX) // contig is not placed
	{
	  NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);

	  assert(unitig != NULL);
	  
	  if(unitig->info.CI.numInstances == 0) // this unitig has not been placed as a surrogate
	  {
		return 0;
	  }
	  else
	  {
		return 1;  // The unitig in this contig appears as a surrogate elsewhere in the assembly
	  }
	}
	else
	{ // Contig is placed
	  return 0;
	}     
  }
}



void OutputContigLinks(ScaffoldGraphT *graph, int outputOverlapOnlyContigEdges)
{
  IntContigLinkMesg		clm;
  GenericMesg			pmesg;
  GraphNodeIterator nodes;
  ContigT *ctg;
  pmesg.m = &clm;
  pmesg.t = MESG_ICL;

  fprintf(GlobalData->stderrc,"* OutputContigLinks *\n");
    if(JumpList == NULL){
      fprintf(GlobalData->stderrc,"* Creating JumpList *\n");
      JumpList = CreateVA_IntMate_Pairs(256);
    }

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((ctg = NextGraphNodeIterator(&nodes)) != NULL){
    ContigT		*mate;
    GraphEdgeIterator	edges;
    CIEdgeT		*edge, *redge;
    CIFragT		*frag;
    int 		edgeTotal;
    int 		edgeCount;	// This var used for sanity checks
    IntMate_Pairs	imp;
    // MateStatType	mstat;

    if(ctg->flags.bits.isChaff)
      continue;

	if (SurrogatedSingleUnitigContig( ctg ))
	  continue;

    clm.contig1 = ctg->outputID;
    InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, ctg->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
    while((edge = NextGraphEdgeIterator(&edges)) != NULL){

      if (edge->idA != ctg->id)
	continue;

      ResetVA_IntMate_Pairs(JumpList);

      mate = GetGraphNode(ScaffoldGraph->ContigGraph, edge->idB);
      if(mate->flags.bits.isChaff)
	continue;

	  if (SurrogatedSingleUnitigContig( mate ))
		continue;

      clm.contig2 = mate->outputID;

      /* Don't need to map orientation, always using canonical orientation*/
      clm.orientation = edge->orient;
      if(!isOverlapEdge(edge)){
	clm.overlap_type = AS_NO_OVERLAP;
      }else if (edge->flags.bits.hasTandemOverlap){
	clm.overlap_type = AS_TANDEM_OVERLAP;
      }else {
	clm.overlap_type = AS_OVERLAP;
      }


      switch(GetEdgeStatus(edge)){
      case LARGE_VARIANCE_EDGE_STATUS:
      case UNKNOWN_EDGE_STATUS:
      case INTER_SCAFFOLD_EDGE_STATUS:
	clm.status = AS_UNKNOWN_IN_ASSEMBLY;
        break;
      case TENTATIVE_TRUSTED_EDGE_STATUS:
      case TRUSTED_EDGE_STATUS:
	clm.status = AS_IN_ASSEMBLY;
        break;
      case TENTATIVE_UNTRUSTED_EDGE_STATUS:
      case UNTRUSTED_EDGE_STATUS:
	clm.status = AS_BAD;
        break;
      default:
	assert(0 /* Invalid edge status */);
      }

      clm.is_possible_chimera = edge->flags.bits.isPossibleChimera;
      clm.includes_guide = edge->flags.bits.hasGuide;
      clm.mean_distance = edge->distance.mean;
      clm.std_deviation = sqrt(edge->distance.variance);
      edgeTotal = clm.num_contributing = edge->edgesContributing;
      if (clm.overlap_type != AS_NO_OVERLAP)
	--edgeTotal;
      if (!edgeTotal && !outputOverlapOnlyContigEdges)
	continue;	// don't output pure overlap edges


      if (edge->flags.bits.isRaw) {
	assert(edgeTotal <= 1);		// sanity check
	if(edgeTotal == 1){
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragA);
	  //	  frag->outMateStat = mstat;
	  frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
	  imp.in1 = frag->iid;
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragB);
	  frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
	  //	  frag->outMateStat = mstat;
	  imp.in2 = frag->iid;
	}else{
	  imp.in1 = imp.in2 = 0;
	}
	if(isOverlapEdge(edge)){
	  assert(outputOverlapOnlyContigEdges);
	  imp.type = 'X';
	}else{
	  if(edge->flags.bits.hasGuide)
	    imp.type = AS_BAC_GUIDE;
	  else if(edge->flags.bits.hasSTSGuide)
	    imp.type = AS_STS_GUIDE;
	  else if(edge->flags.bits.hasMayJoin)
	    imp.type = AS_MAY_JOIN;
	  else if(edge->flags.bits.hasMustJoin)
	    imp.type = AS_MUST_JOIN;
	  else
	    imp.type = AS_MATE;
	  AppendIntMate_Pairs(JumpList, &imp);
	}
      }
      else {
	redge = edge;

	assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge
	edgeCount = 0;

	while (redge->nextRawEdge != NULLINDEX) {
	  redge = GetGraphEdge(ScaffoldGraph->ContigGraph,redge->nextRawEdge);
	  if (isOverlapEdge(redge))
	    continue;		// overlap edges don't count
	  ++edgeCount;
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragA);
	  imp.in1 = frag->iid;
	  frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragB);
	  frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
	  imp.in2 = frag->iid;
	assert(!isOverlapEdge(redge));
	if(redge->flags.bits.hasGuide)
	  imp.type = AS_BAC_GUIDE;
	else if(redge->flags.bits.hasSTSGuide)
	  imp.type = AS_STS_GUIDE;
	else if(redge->flags.bits.hasMayJoin)
	  imp.type = AS_MAY_JOIN;
	else if(redge->flags.bits.hasMustJoin)
	  imp.type = AS_MUST_JOIN;
	else
	  imp.type = AS_MATE;
	AppendIntMate_Pairs(JumpList, &imp);
	}
	assert(GetNumIntMate_Pairss(JumpList) == edgeTotal);
	assert(edgeCount == edgeTotal);
      }		// if (edge . . . 
      clm.jump_list = GetIntMate_Pairs(JumpList,0);
      (GlobalData->writer)(GlobalData->outfp2,&pmesg);
    }	// while (edge . . .
  }	// for (i . . .
  fflush(NULL);

}


void OutputScaffoldLink(ScaffoldGraphT * graph,
                        CIScaffoldT * scaffold,
                        CIEdgeT * edge)
{
  InternalScaffoldLinkMesg slm;
  GenericMesg pmesg;
  CIScaffoldT *mate;
  CIEdgeT *redge;
  CIFragT *frag;
  IntMate_Pairs	imp;
  int edgeTotal = 0;
  int edgeCount = 0; // This var used for sanity checks
  
  pmesg.m = &slm;
  pmesg.t = MESG_ISL;
  slm.iscaffold1 = scaffold->id;
  
  if(JumpList == NULL)
    JumpList = CreateVA_IntMate_Pairs(256);
  else
    ResetVA_IntMate_Pairs(JumpList);
    
  mate = GetGraphNode(ScaffoldGraph->ScaffoldGraph, edge->idB);
  
  slm.iscaffold2 = mate->id;
  
  /* Don't need to map orientation, always using canonical orientation*/
  slm.orientation = edge->orient;
  assert(!isOverlapEdge(edge));
  
  slm.includes_guide = edge->flags.bits.hasGuide;
  slm.mean_distance = edge->distance.mean;
  slm.std_deviation = sqrt(edge->distance.variance);
  edgeTotal = slm.num_contributing = edge->edgesContributing;
  
  redge = edge;
#if 1
  if(edgeTotal < 2)
    return;
#endif
  if (edge->flags.bits.isRaw) {
    assert(edgeTotal <= 1);		// sanity check
    if(edgeTotal == 1){
      frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragA);
      //	  frag->outMateStat = mstat;
      frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
      imp.in1 = frag->iid;
      frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragB);
      frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
      //	  frag->outMateStat = mstat;
      imp.in2 = frag->iid;
    }else{
      imp.in1 = imp.in2 = 0;
    }
    if(edge->flags.bits.hasGuide)
      imp.type = AS_BAC_GUIDE;
    else if(edge->flags.bits.hasSTSGuide)
      imp.type = AS_STS_GUIDE;
    else if(edge->flags.bits.hasMayJoin)
      imp.type = AS_MAY_JOIN;
    else if(edge->flags.bits.hasMustJoin)
      imp.type = AS_MUST_JOIN;
    else
      imp.type = AS_MATE;
    AppendIntMate_Pairs(JumpList, &imp);
    edgeCount = 1;
  }else{
    
    assert(!edge->flags.bits.isRaw);
    
    assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge
    
    edgeCount = 0;
    
    while (redge->nextRawEdge != NULLINDEX) {
      redge = GetGraphEdge(ScaffoldGraph->ScaffoldGraph,redge->nextRawEdge);
      assert(!isOverlapEdge(redge));
      ++edgeCount;
      frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragA);
      imp.in1 = frag->iid;
      frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
      frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragB);
      frag->flags.bits.edgeStatus = GetEdgeStatus(edge);
      imp.in2 = frag->iid;
      if(redge->flags.bits.hasGuide)
        imp.type = AS_BAC_GUIDE;
      else if(redge->flags.bits.hasSTSGuide)
        imp.type = AS_STS_GUIDE;
      else if(redge->flags.bits.hasMayJoin)
        imp.type = AS_MAY_JOIN;
      else if(redge->flags.bits.hasMustJoin)
        imp.type = AS_MUST_JOIN;
      else
        imp.type = AS_MATE;
      AppendIntMate_Pairs(JumpList, &imp);
    }
  }
  assert(GetNumIntMate_Pairss(JumpList) == edgeTotal);
  assert(edgeCount == edgeTotal);
  slm.jump_list = GetIntMate_Pairs(JumpList,0);
  (GlobalData->writer)(GlobalData->outfp2,&pmesg);
}


void OutputScaffoldLinksForScaffold(ScaffoldGraphT * graph,
                                    CIScaffoldT * scaffold)
{
  GraphEdgeIterator	edges;
  CIEdgeT		*edge;
  
  InitGraphEdgeIterator(graph->ScaffoldGraph, scaffold->id,
                        ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
  while((edge = NextGraphEdgeIterator(&edges)) != NULL)
  {
    if (edge->idA != scaffold->id)
      continue;
    OutputScaffoldLink(graph, scaffold, edge);
  }
}


void OutputScaffoldLinks(ScaffoldGraphT *graph)
{
  GraphNodeIterator nodes;
  CIScaffoldT *scaffold;

  fprintf(GlobalData->stderrc,"* OutputScaffoldLinks *\n");
  if(JumpList == NULL){
    fprintf(GlobalData->stderrc,"* Creating JumpList *\n");
    JumpList = CreateVA_IntMate_Pairs(256);
  }

  InitGraphNodeIterator(&nodes, graph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&nodes)) != NULL)
  {
    OutputScaffoldLinksForScaffold(graph, scaffold);
  }

  fflush(NULL);
}

/********************************************************************************/
void OutputUnitigsFromMultiAligns(void){
  GenericMesg			pmesg;
  ContigT			*ci;
  IntUnitigMesg			ium_mesg;
  GraphNodeIterator nodes;
  int numCIs = (int) GetNumGraphNodes(ScaffoldGraph->CIGraph);
  CDS_CID_t cid = 0;
  MultiAlignT *ma = CreateEmptyMultiAlignT();

  pmesg.m = &ium_mesg;
  pmesg.t = MESG_IUM;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while((ci = NextGraphNodeIterator(&nodes)) != NULL){
    UnitigStatus   status;
    //    MultiAlignT *ma = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, ci->id);

    assert(ci->id>=0 && ci->id< numCIs);

    if(ci->flags.bits.isChaff){
      //      fprintf(GlobalData->stderrc,"* # Unitig " F_CID " is CHAFF\n", ci->id);
      continue;
    }
    switch(ci->type){
    case DISCRIMINATORUNIQUECHUNK_CGW:
      status = AS_UNIQUE;
      //	fprintf(GlobalData->stderrc,"* Unitig " F_CID " is UNIQUE: DISCRIMINATOR output " F_CID " \n",ci->id, cid);
      break;
    case UNIQUECHUNK_CGW:
      status = AS_UNIQUE;
      //	fprintf(GlobalData->stderrc,"* Unitig " F_CID " is UNIQUE: output " F_CID " \n",ci->id, cid);
      break;
    case UNRESOLVEDCHUNK_CGW:
      if(ci->info.CI.numInstances > 0){
	assert(!ci->flags.bits.isUnique);
	status = AS_SEP;
	//	fprintf(GlobalData->stderrc,"* Unitig " F_CID " has %d instances--- output " F_CID " SEP\n",ci->id, ci->info.CI.numInstances,cid);
      }else{
	if(ci->scaffoldID != NULLINDEX){
	  //	  fprintf(GlobalData->stderrc,"* Unitig " F_CID " has %d instances--- output " F_CID " UNIQUE\n",ci->id, ci->info.CI.numInstances,cid);
	  status = AS_UNIQUE;
	}else{
	  //	  fprintf(GlobalData->stderrc,"* Unitig " F_CID " has %d instances--- output " F_CID " NOTREZ\n",ci->id, ci->info.CI.numInstances,cid);
	  status = AS_NOTREZ;
	}
      }
      break;
    case RESOLVEDREPEATCHUNK_CGW:
      /* SKIP THESE */
      //      fprintf(GlobalData->stderrc,"* Skipping unitig " F_CID " --- RESOLVEDREPEAT\n",ci->id);
      continue;
    default:
      assert(0);
    }
    ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ma, ci->id, TRUE);
    {
    CDS_IID_t numFrag = GetNumIntMultiPoss(ma->f_list);
    CDS_IID_t * tmpSource;
    assert (ci->type != CONTIG_CGW);

    tmpSource = safe_malloc(sizeof(CDS_IID_t) * (GetNumIntMultiPoss(ma->f_list) + 1));
    {
      CDS_IID_t i;
      // Null out the source field
      for(i = 0; i < numFrag; i++){
	IntMultiPos *mp_i = GetIntMultiPos(ma->f_list,i);
	tmpSource[i] = (CDS_IID_t)mp_i->source;
	mp_i->source = NULL;
	assert(mp_i->ident);
      }
    }
    ci->outputID = cid++;
    //assert(ci->outputID == ci->id); // TRUE FOR UNITIGS UNTIL WE SPLIT
    
    ium_mesg.iaccession = ci->id;
#ifdef DEBUG_DATA
    ium_mesg.source = Getchar(ScaffoldGraph->SourceFields, ci->info.CI.source);
#else
    ium_mesg.source = NULL;
#endif

    ium_mesg.coverage_stat = ci->info.CI.coverageStat;
    ium_mesg.status = status;
    ium_mesg.a_branch_point = ci->info.CI.branchPointA;
    ium_mesg.b_branch_point = ci->info.CI.branchPointB;
    ium_mesg.length = GetMultiAlignLength(ma);
    ium_mesg.consensus = Getchar(ma->consensus,0);
    ium_mesg.quality = Getchar(ma->quality,0);
    ium_mesg.forced = 0;
    ium_mesg.num_frags = GetNumIntMultiPoss(ma->f_list);
    ium_mesg.f_list = GetIntMultiPos(ma->f_list,0);

    (GlobalData->writer)(GlobalData->outfp,&pmesg);
    {
      CDS_IID_t i;
      // Restore the source field
      for(i = 0; i < numFrag; i++){
	IntMultiPos *mp_i = GetIntMultiPos(ma->f_list,i);
	mp_i->source = (char *)tmpSource[i];
      }
    }
    //    UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ci->id, TRUE);
    free(tmpSource);
    }
  }	// while NextGraphNode
  DeleteMultiAlignT(ma);
  fflush(NULL);

}



void OutputUnitigLinksFromMultiAligns(void){
    IntUnitigLinkMesg		ulm;
    GenericMesg			pmesg;
    GraphNodeIterator nodes;
    ChunkInstanceT *ci;
    pmesg.m = &ulm;
    pmesg.t = MESG_IUL;

    fprintf(GlobalData->stderrc,"* OutputUnitigLinksFromMultiAligns *\n");

    if(JumpList == NULL){
      fprintf(GlobalData->stderrc,"* Creating JumpList *\n");
      JumpList = CreateVA_IntMate_Pairs(256);
      AssertPtr(JumpList);
    }
    fflush(NULL);

    InitGraphNodeIterator(&nodes, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
    while((ci = NextGraphNodeIterator(&nodes)) != NULL){
      ContigT		*mate;
      GraphEdgeIterator	edges;
      CIEdgeT		*edge, *redge;
      CIFragT		*frag;
      int 		edgeTotal;
      int 		edgeCount;	// This var used for sanity checks
      IntMate_Pairs     imp;

      AssertPtr(ci);
      assert (ci->type != CONTIG_CGW);
    
      // We skip these...
      if(ci->type == RESOLVEDREPEATCHUNK_CGW)
	continue;
      if(ci->flags.bits.isChaff)
	continue;

      if(ci->id % 50000 == 0)
	fprintf(GlobalData->stderrc,"* Outputing links incident on unitig " F_CID "\n", ci->id);

      ulm.unitig1 = ci->id;
      InitGraphEdgeIterator(ScaffoldGraph->CIGraph, ci->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
      while((edge = NextGraphEdgeIterator(&edges)) != NULL){
	AssertPtr(edge);
	ResetVA_IntMate_Pairs(JumpList);
	if (edge->idA != ci->id ||
	    edge->flags.bits.isInferred ||
	    edge->flags.bits.isInferredRemoved ||
	    edge->flags.bits.isMarkedForDeletion)
	  continue;
	ulm.unitig2 = edge->idB;
	/* Don't need to map orientation, always using canonical orientation*/
	ulm.orientation = edge->orient;
	if(!isOverlapEdge(edge)){
	  ulm.overlap_type = AS_NO_OVERLAP;
	}else if (edge->flags.bits.hasTandemOverlap){
	  ulm.overlap_type = AS_TANDEM_OVERLAP;
	}else {
	  ulm.overlap_type = AS_OVERLAP;
	}

	ulm.is_possible_chimera = edge->flags.bits.isPossibleChimera;
	ulm.includes_guide = edge->flags.bits.hasGuide;
	ulm.mean_distance = edge->distance.mean;
	ulm.std_deviation = sqrt(edge->distance.variance);
	edgeTotal = ulm.num_contributing = edge->edgesContributing;

	if (ulm.overlap_type != AS_NO_OVERLAP)
	  --edgeTotal;
	if (!edgeTotal)
	  continue;	// don't output pure overlap edges

	mate = GetGraphNode(ScaffoldGraph->CIGraph, edge->idB);
	if(mate->flags.bits.isChaff)
	  continue;

	  {
	  int numBad = 0;
	  int numGood = 0;
	  int numUnknown = 0;
	  CIFragT *fragA, *fragB;
	  // Look through the fragment pairs in this edge.  If any of the fragments are
	  // marked BAD ==> bad
	  // Else, if any are marked good ==> good
	  // Else, mark it unknown
	  if(edge->flags.bits.isRaw){
	    redge = edge;
	  }else{
	    redge = GetGraphEdge(ScaffoldGraph->CIGraph, edge->nextRawEdge);
	  }
	  assert(redge && edge);

	  for(; redge != NULL; redge = GetGraphEdge(ScaffoldGraph->CIGraph, redge->nextRawEdge)){
	    AssertPtr(redge);
	    if(isOverlapEdge(redge))
	      continue;
	    fragA = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragA);
	    fragB = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragB);

	    assert(fragA && fragB);
	    assert(fragA->flags.bits.edgeStatus== fragB->flags.bits.edgeStatus);
		
	    if(fragA->flags.bits.edgeStatus == UNTRUSTED_EDGE_STATUS ||
	       fragA->flags.bits.edgeStatus == TENTATIVE_UNTRUSTED_EDGE_STATUS){
	      numBad++;
	    }else if(fragA->flags.bits.edgeStatus == TRUSTED_EDGE_STATUS ||
		    fragA->flags.bits.edgeStatus == TENTATIVE_TRUSTED_EDGE_STATUS){
	      numGood++;
	    } else{
	      numUnknown++;
	    }
	  }
	  
	  if(numBad > 0){
	    	  ulm.status = AS_BAD;
	  }else if(numGood > 0){
	    	  ulm.status = AS_IN_ASSEMBLY;
	  }else ulm.status = AS_UNKNOWN_IN_ASSEMBLY;

	}
	      
	  ResetVA_IntMate_Pairs(JumpList);

	  if (edge->flags.bits.isRaw) {
	    assert(edgeTotal == 1);		// sanity check
	    edgeCount = 1;
	    frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragA);
	    AssertPtr(frag);
	    imp.in1 = frag->iid;
	    frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragB);
	    imp.in2 = frag->iid;
	    AssertPtr(frag);
	    assert(!isOverlapEdge(edge));
	    if(edge->flags.bits.hasGuide)
	      imp.type = AS_BAC_GUIDE;
	    else if(edge->flags.bits.hasSTSGuide)
	      imp.type = AS_STS_GUIDE;
	    else if(edge->flags.bits.hasMayJoin)
	      imp.type = AS_MAY_JOIN;
	    else if(edge->flags.bits.hasMustJoin)
	      imp.type = AS_MUST_JOIN;
	    else
	      imp.type = AS_MATE;
	    AppendIntMate_Pairs(JumpList,&imp);
	  } else { // not raw
	    redge = edge;

	    assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge
	    edgeCount = 0;
	    assert(edgeTotal > 0);

	    while(redge->nextRawEdge != NULLINDEX) {
	      redge = GetGraphEdge(ScaffoldGraph->CIGraph,redge->nextRawEdge);
	      if (isOverlapEdge(redge))
		continue;		// overlap edges don't count
	      ++edgeCount;
	      frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragA);
	      AssertPtr(frag);
	      imp.in1 = frag->iid;
	      frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragB);
	      AssertPtr(frag);
	      imp.in2 = frag->iid;
	      assert(!isOverlapEdge(redge));
	      if(redge->flags.bits.hasGuide)
		imp.type = AS_BAC_GUIDE;
	      else if(redge->flags.bits.hasSTSGuide)
		imp.type = AS_STS_GUIDE;
	      else if(redge->flags.bits.hasMayJoin)
		imp.type = AS_MAY_JOIN;
	      else if(redge->flags.bits.hasMustJoin)
		imp.type = AS_MUST_JOIN;
	      else
		imp.type = AS_MATE;

	      AppendIntMate_Pairs(JumpList,&imp);
	    }
	  }
	  assert(GetNumIntMate_Pairss(JumpList) == edgeTotal);

	  if(edgeCount != edgeTotal){
	    fprintf(GlobalData->stderrc,"* edgeCount = %d edgeTotal = %d\n",
		    edgeCount, edgeTotal);
	    PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->CIGraph," ", edge, edge->idA);
	    fflush(GlobalData->stderrc);

	    redge = edge;
	    assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge

	    while(redge->nextRawEdge != NULLINDEX) {
	      redge = GetGraphEdge(ScaffoldGraph->CIGraph,redge->nextRawEdge);
	      PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->CIGraph," ", redge, redge->idA);
	    }
	    assert(edgeCount == edgeTotal);
	  }
	ulm.jump_list = GetIntMate_Pairs(JumpList,0);
	(GlobalData->writer)(GlobalData->outfp,&pmesg);
      }
    }
  fflush(NULL);

}



void OutputContigLinksFromMultiAligns(void){

}


#if 0
void OutputUnitigLinks(ScaffoldGraphT *graph)
{
  IntUnitigLinkMesg		ulm;
  GenericMesg			pmesg;
  GraphNodeIterator nodes;
  ChunkInstanceT *ci;
  pmesg.m = &ulm;
  pmesg.t = MESG_IUL;


    if(JumpList == NULL){
      JumpList = CreateVA_IntMate_Pairs(256);
    }

  InitGraphNodeIterator(&nodes, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while((ci = NextGraphNodeIterator(&nodes)) != NULL){
    ContigT		*mate;
    GraphEdgeIterator	edges;
    CIEdgeT		*edge, *redge;
    CIFragT		*frag;
    int 		edgeTotal;
    int 		edgeCount;	// This var used for sanity checks
    IntMate_Pairs	imp;

    assert (ci->type != CONTIG_CGW);
    
    ulm.unitig1 = ci->id;
    InitGraphEdgeIterator(ScaffoldGraph->CIGraph, ci->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
    while((edge = NextGraphEdgeIterator(&edges)) != NULL){
      if (edge->idA != ci->id ||
	  edge->flags.bits.isInferred ||
	  edge->flags.bits.isInferredRemoved ||
	  edge->flags.bits.isMarkedForDeletion)
	continue;
      ulm.unitig2 = edge->idB;
      /* Don't need to map orientation, always using canonical orientation*/
      ulm.orientation = edge->orient;
      if(!isOverlapEdge(edge)){
	ulm.overlap_type = AS_NO_OVERLAP;
      }else if (edge->flags.bits.hasTandemOverlap){
	ulm.overlap_type = AS_TANDEM_OVERLAP;
      }else {
	ulm.overlap_type = AS_OVERLAP;
      }

      ulm.is_possible_chimera = edge->flags.bits.isPossibleChimera;
      ulm.includes_guide = edge->flags.bits.hasGuide;
      ulm.mean_distance = edge->distance.mean;
      ulm.std_deviation = sqrt(edge->distance.variance);
      edgeTotal = ulm.num_contributing = edge->edgesContributing;
      if (ulm.overlap_type != AS_NO_OVERLAP)
	--edgeTotal;
      if (!edgeTotal)
	continue;	// don't output pure overlap edges

      mate = GetGraphNode(ScaffoldGraph->CIGraph, edge->idB);

      // THIS IS TEMPORARY!!!!
      if (ci->scaffoldID != mate->scaffoldID || ci->scaffoldID < 0){
	ulm.status = AS_UNKNOWN_IN_ASSEMBLY;
      }else {
	float	m1,m2,v1,v2,rv;
	LengthT	*ciLeft, *ciRight, *mateLeft, *mateRight;
	int	oriFlag;	// to record relative orientation

	/* set the mean & variance according to the edge */
	m1 = edge->distance.mean;
	v1 = edge->distance.variance;

	/* now set the mean and varience according to the assembly */
	/* start by determining orientation of the unitigs in the assembly */
	if (ci->offsetAEnd.mean <= ci->offsetBEnd.mean) {
	  ciLeft = &ci->offsetAEnd;	// ci oriented AB
	  ciRight = &ci->offsetBEnd;
	  oriFlag = 0;
	} else { //
	  ciLeft = &ci->offsetBEnd;	// ci oriented BA
	  ciRight = &ci->offsetAEnd;
	  oriFlag = 1;
	}
	if (mate->offsetAEnd.mean <= mate->offsetBEnd.mean) {
	  mateLeft = &(mate->offsetAEnd);	// mate oriented AB
	  mateRight = &(mate->offsetBEnd);
	} else {
	  mateLeft = &(mate->offsetBEnd);	// mate oriented BA
	  mateRight = &(mate->offsetAEnd);
	  ++oriFlag;
	}

	/* If oriFlag is even, unitigs have same relative orientation */
	ulm.status = AS_IN_ASSEMBLY;
	if (oriFlag == 1) {
	  if (edge->orient == AB_AB || edge->orient == BA_BA)
	    ulm.status = AS_BAD;
	} else if (edge->orient == AB_BA || edge->orient == BA_AB)
	  ulm.status = AS_BAD;

	if (ulm.status == AS_IN_ASSEMBLY) {

	  /* now look at the distance from the right side of the left
	     unitig to the left side of the right unitig */
	  if (ciLeft->mean <= mateLeft->mean) {
	    m2 = mateLeft->mean - ciRight->mean;
	    v2 = mateLeft->variance - ciRight->variance;
	  } else {
	    m2 = ciLeft->mean - mateRight->mean;
	    v2 = ciLeft->variance - mateRight->variance;
	  }

	  //	  assert(v2 >= 0.0);
	  if(v2 <= 0.0){
	    fprintf(GlobalData->stderrc,"* Edge between CIs " F_CID " [%d+/-%g,%d+/-%g] and " F_CID " [%d+/-%g,%d+/-%g] has mean:%g var:%g =>1.0\n",
		    edge->idA, 
		    (int)ci->offsetAEnd.mean, ci->offsetAEnd.variance,
		    (int)ci->offsetBEnd.mean, ci->offsetBEnd.variance,
		    edge->idB,
		    (int)mate->offsetAEnd.mean, mate->offsetAEnd.variance,
		    (int)mate->offsetBEnd.mean, mate->offsetBEnd.variance,
		    m2, v2);
	    v2 = 1.0;
	  }
	  if (!PairwiseChiSquare(m1,v1,m2,v2,NULL,&rv,
				 PAIRWISECHI2THRESHOLD_CGW))
	    ulm.status = AS_BAD;
	}
      }

      ResetVA_IntMate_Pairs(JumpList);

      if (edge->flags.bits.isRaw) {
	assert(edgeTotal == 1);		// sanity check
	frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragA);
	imp.in1 = frag->iid;
	frag = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragB);
	imp.in2 = frag->iid;
	assert(!isOverlapEdge(edge));
	if(edge->flags.bits.hasGuide)
	  imp.type = AS_BAC_GUIDE;
	else if(edge->flags.bits.hasSTSGuide)
	  imp.type = AS_STS_GUIDE;
	else if(edge->flags.bits.hasMayJoin)
	  imp.type = AS_MAY_JOIN;
	else if(edge->flags.bits.hasMustJoin)
	  imp.type = AS_MUST_JOIN;
	else
	  imp.type = AS_MATE;
	AppendIntMate_Pairs(JumpList, &imp);
      } else {
	redge = edge;

	assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge
	edgeCount = 0;

	while(redge->nextRawEdge != NULLINDEX) {
	  redge = GetGraphEdge(ScaffoldGraph->CIGraph,redge->nextRawEdge);
	  if (isOverlapEdge(redge))
	    continue;		// overlap edges don't count
	  ++edgeCount;
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragA);
	  imp.in1 = frag->iid;
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragB);
	  imp.in2 = frag->iid;
	  assert(!isOverlapEdge(redge));
	  if(redge->flags.bits.hasGuide)
	    imp.type = AS_BAC_GUIDE;
	  else if(redge->flags.bits.hasSTSGuide)
	    imp.type = AS_STS_GUIDE;
	  else if(redge->flags.bits.hasMayJoin)
	    imp.type = AS_MAY_JOIN;
	  else if(redge->flags.bits.hasMustJoin)
	    imp.type = AS_MUST_JOIN;
	  else
	    imp.type = AS_MATE;
	  AppendIntMate_Pairs(JumpList, &imp);
	}
	if(edgeCount != edgeTotal){
	  fprintf(GlobalData->stderrc,"* edgeCount = %d edgeTotal = %d\n",
		  edgeCount, edgeTotal);
	  PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->CIGraph," ", edge, edge->idA);

	  redge = edge;
	  assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge

	while(redge->nextRawEdge != NULLINDEX) {
	  redge = GetGraphEdge(ScaffoldGraph->CIGraph,redge->nextRawEdge);
	  PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->CIGraph," ", redge, redge->idA);
	}
	assert(edgeCount == edgeTotal);
	}
      }
      ulm.jump_list = GetIntMate_Pairs(JumpList,0);
      (GlobalData->writer)(GlobalData->outfp,&pmesg);
    }	
  }	// while (edge . . .
}
#endif

void OutputScaffolds(ScaffoldGraphT *graph)
{
  CDS_IID_t			sid, pairCount;
  IntScaffoldMesg		ism;
  int				buffSize=2048;
  GenericMesg			pmesg;
  IntContigPairs		*cptr;
  CIScaffoldT			*scaf;
  CIScaffoldTIterator		Contigs;
  ChunkInstanceT		*curr, *last;
  GraphNodeIterator             scaffolds;
  CDS_IID_t cnt = 0;

  pmesg.m = &ism;
  pmesg.t = MESG_ISF;

  ism.contig_pairs = (IntContigPairs *) safe_malloc(sizeof(IntContigPairs)*buffSize);
  assert(ism.contig_pairs != NULL);

  InitGraphNodeIterator(&scaffolds, graph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaf = NextGraphNodeIterator(&scaffolds)) != NULL){
    ChunkOrient orientLast, orientCurr;
    sid = scaf->outputID = cnt++;

    ism.iaccession = scaf->id;
    ism.num_contig_pairs = scaf->info.Scaffold.numElements-1;
    if(scaf->type != REAL_SCAFFOLD)
      continue;
    assert(ism.num_contig_pairs >= 0);

    if (ism.num_contig_pairs > buffSize) {
      buffSize = ism.num_contig_pairs * 2;
      ism.contig_pairs = (IntContigPairs *) safe_realloc(ism.contig_pairs, sizeof(IntContigPairs)*buffSize);
      assert(ism.contig_pairs != NULL);
    }
    cptr = ism.contig_pairs;

    InitCIScaffoldTIterator(graph, scaf, TRUE, FALSE, &Contigs);
    last = NextCIScaffoldTIterator(&Contigs);
    orientLast = (last->offsetAEnd.mean < last->offsetBEnd.mean) ? A_B : B_A;

    assert(last->scaffoldID == scaf->id);

    if(ism.num_contig_pairs == 0){
	cptr->contig1 = last->outputID;
	cptr->contig2 = last->outputID;
	cptr->mean = 0.0;
	cptr->stddev = 0.0;
	cptr->orient = AB_AB; // got to put something
    }else{
    pairCount = 0;	// only used for sanity check
    while((curr = NextCIScaffoldTIterator(&Contigs)) != NULL){
      assert(pairCount < ism.num_contig_pairs);

      assert(curr->scaffoldID == scaf->id);

      cptr->contig1 = last->outputID;
      cptr->contig2 = curr->outputID;
      orientCurr = (curr->offsetAEnd.mean < curr->offsetBEnd.mean) ? A_B : B_A;
      if(orientLast == A_B){
	if(orientCurr == A_B){
	  cptr->mean = curr->offsetAEnd.mean - last->offsetBEnd.mean;
	  cptr->stddev = sqrt(curr->offsetAEnd.variance -
			      last->offsetBEnd.variance);
	  cptr->orient = AB_AB;
	}else{//orientCurr == B_A
	  cptr->mean = curr->offsetBEnd.mean - last->offsetBEnd.mean;
	  cptr->stddev = sqrt(curr->offsetBEnd.variance -
			      last->offsetBEnd.variance);
	  cptr->orient = AB_BA;
	}
      }else{//orientLast == B_A
	if(orientCurr == A_B){
	  cptr->mean = curr->offsetAEnd.mean - last->offsetAEnd.mean;
	  cptr->stddev = sqrt(curr->offsetAEnd.variance -
			      last->offsetAEnd.variance);
	  cptr->orient = BA_AB;
	}else{//orientCurr == B_A
	  cptr->mean = curr->offsetBEnd.mean - last->offsetAEnd.mean;
	  cptr->stddev = sqrt(curr->offsetBEnd.variance -
			      last->offsetAEnd.variance);
	  cptr->orient = BA_BA;
	}
      }
      last = curr;
      orientLast = orientCurr;
      ++cptr;
      ++pairCount;
    }		// while (curr . . .
    }
    (GlobalData->writer)(GlobalData->outfp2,&pmesg);
  }		// for (sid=0; . . .
  free(ism.contig_pairs);
  fflush(NULL);
  return;
}



#if 0
	/* set the mean & variance according to the edge */
	m1 = edge->distance.mean;
	v1 = edge->distance.variance;

	/* now set the mean and varience according to the assembly */
	/* start by determining orientation of the unitigs in the assembly */
	if (ci->offsetAEnd.mean <= ci->offsetBEnd.mean) {
	  ciLeft = &ci->offsetAEnd;	// ci oriented AB
	  ciRight = &ci->offsetBEnd;
	  oriFlag = 0;
	} else {
	  ciLeft = &ci->offsetBEnd;	// ci oriented BA
	  ciRight = &ci->offsetAEnd;
	  oriFlag = 1;
	}
	if (mate->offsetAEnd.mean <= mate->offsetBEnd.mean) {
	  mateLeft = &(mate->offsetAEnd);	// mate oriented AB
	  mateRight = &(mate->offsetBEnd);
	} else {
	  mateLeft = &(mate->offsetBEnd);	// mate oriented BA
	  mateRight = &(mate->offsetAEnd);
	  ++oriFlag;
	}

	/* If oriFlag is even, unitigs have same relative orientation */
	ulm.status = AS_IN_ASSEMBLY;
	if (oriFlag == 1) {
	  if (edge->orient == AB_AB || edge->orient == BA_BA)
	    ulm.status = AS_BAD;
	} else if (edge->orient == AB_BA || edge->orient == BA_AB)
	  ulm.status = AS_BAD;

	if (ulm.status == AS_IN_ASSEMBLY) {

	  /* now look at the distance from the right side of the left
	     unitig to the left side of the right unitig */
	  if (ciLeft->mean <= mateLeft->mean) {
	    m2 = mateLeft->mean - ciRight->mean;
	    v2 = mateLeft->variance - ciRight->variance;
	  } else {
	    m2 = ciLeft->mean - mateRight->mean;
	    v2 = ciLeft->variance - mateRight->variance;
	  }

	  //	  assert(v2 >= 0.0);
	  if(v2 <= 0.0){
	    fprintf(GlobalData->stderrc,"* Edge between CIs " F_CID " [%d+/-%g,%d+/-%g] and " F_CID " [%d+/-%g,%d+/-%g] has mean:%g var:%g =>1.0\n",
		    edge->idA, 
		    (int)ci->offsetAEnd.mean, ci->offsetAEnd.variance,
		    (int)ci->offsetBEnd.mean, ci->offsetBEnd.variance,
		    edge->idB,
		    (int)mate->offsetAEnd.mean, mate->offsetAEnd.variance,
		    (int)mate->offsetBEnd.mean, mate->offsetBEnd.variance,
		    m2, v2);
	    v2 = 1.0;
	  }
	  if (!PairwiseChiSquare(m1,v1,m2,v2,NULL,&rv,
				 PAIRWISECHI2THRESHOLD_CGW))
	    ulm.status = AS_BAD;
	}
      }
#endif

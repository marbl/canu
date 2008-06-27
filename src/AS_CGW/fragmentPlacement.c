
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_timer.h"
#include "GraphCGW_T.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "fragmentPlacement.h"


#define MAX_SIGMA_SLOP 3.0

#define EXTERNAL_MATES_ONLY 1
#define ALL_MATES 0


// new iterator for looping over the fragments in a chunk -- note need to call cleanup routine when done

  /* we need to have:
     an ma for the  top chunk
     if (type == contig)
     an iterator over chunks in top ma
     an ma for sub chunks
     an iterator over fragments in the (sub) ma

     the fragment iterator should be set to go after init

     each call to iterator should
     compute next in fragment iterator
     if non null, return
     else
     iterate to next sub chunk
     if null, return failure
     else
     restart fragment iterator
     compute next


     return the current
  */

typedef struct FragIterator{
  uint32 isUtg:1;
  uint32 isCtg:1;
  uint32 includeSurrogateFrgs:1;
  CDS_CID_t id;
  MultiAlignT *fragiterma;
  int32 fOrder; /* which fragment within multialign is next */
  ContigTIterator subchunks;
  struct FragIterator *subchunkIterator;
} CGWFragIterator;

// New iterator for finding the mates of a fragment
typedef struct  {
  CDS_CID_t  thisFragIID;
  CDS_CID_t  nextLink;
  NodeCGW_T *node;
  uint32     external_only:1;
} CGWMateIterator;






static
void GetRelationBetweenFragsOnChunks(CIFragT *frag,ChunkInstanceT*fragChunk,
				     CIFragT *mate,ChunkInstanceT*mateChunk,
				     LengthT*separation,OrientType* edgeOri){


  // we pass in chunks because if a fragment in question is in a surrogate, we can't
  // get the info we need by identifying the parent chunk containing the frag

  NodeOrient fOri,mOri,fCIOri,mCIOri,fCtgOri,mCtgOri,fOriOnScf,mOriOnScf;
  double f5pOnCI,m5pOnCI,f5pOnCtg,m5pOnCtg,f5pOnScf,m5pOnScf;
  NodeCGW_T *fragCtg, *mateCtg;

  assert(fragChunk->flags.bits.isCI);
  assert(mateChunk->flags.bits.isCI);


  f5pOnCI =  frag->offset5p.mean;
  fOri = (frag->offset5p.mean < frag->offset3p.mean) ? A_B : B_A;

  fCIOri = (fragChunk->offsetAEnd.mean < fragChunk->offsetBEnd.mean)? A_B : B_A;
  f5pOnCtg =  fragChunk->offsetAEnd.mean + ((fCIOri==A_B) ? f5pOnCI : -f5pOnCI);

  fragCtg = GetGraphNode(ScaffoldGraph->ContigGraph,fragChunk->info.CI.contigID);
  fCtgOri = (fragCtg->offsetAEnd.mean < fragCtg->offsetBEnd.mean) ? A_B : B_A;
  f5pOnScf = fragCtg->offsetAEnd.mean + ((fCtgOri == A_B) ? f5pOnCtg : -f5pOnCtg);


  m5pOnCI =  mate->offset5p.mean;
  mOri = (mate->offset5p.mean < mate->offset3p.mean) ? A_B : B_A;

  mCIOri = (mateChunk->offsetAEnd.mean < mateChunk->offsetBEnd.mean) ? A_B : B_A;
  m5pOnCtg =  mateChunk->offsetAEnd.mean + ((mCIOri==A_B) ? m5pOnCI : -m5pOnCI);

  mateCtg = GetGraphNode(ScaffoldGraph->ContigGraph,mateChunk->info.CI.contigID);
  mCtgOri = (mateCtg->offsetAEnd.mean < mateCtg->offsetBEnd.mean) ? A_B : B_A;
  m5pOnScf = mateCtg->offsetAEnd.mean + ((mCtgOri == A_B) ? m5pOnCtg : -m5pOnCtg);


  separation->mean = abs(f5pOnScf - m5pOnScf);

  separation->variance =
    abs(fragCtg->offsetAEnd.variance - mateCtg->offsetAEnd.variance) + ComputeFudgeVariance(f5pOnCtg) + ComputeFudgeVariance(m5pOnCtg);

  if( ((fOri == A_B)+(fCIOri == A_B)+(fCtgOri == A_B)) % 2 == 0 ){
    fOriOnScf = B_A;
  } else {
    fOriOnScf = A_B;
  }

  if( ((mOri == A_B)+(mCIOri == A_B)+(mCtgOri == A_B)) % 2 == 0 ){
    mOriOnScf = B_A;
  } else {
    mOriOnScf = A_B;
  }

  if(fOriOnScf == A_B){
    if(mOriOnScf == B_A){
      if(f5pOnScf < m5pOnScf){
	*edgeOri = AB_BA;
      } else {
	*edgeOri = BA_AB;
      }
    } else {
      *edgeOri = AB_AB;
    }
  } else { // fOriOnScf == B_A
    if(mOriOnScf == A_B){
      if(f5pOnScf < m5pOnScf){
	*edgeOri = BA_AB;
      } else {
	*edgeOri = AB_BA;
      }
    } else {
      *edgeOri = BA_BA;
    }
  }

}


// FragAndMateAreCompatible() is based partly on PlaceFragmentBasedOnMate()
// from AS_REZ/PlaceFragments.c

static
int FragAndMateAreCompatible(CIFragT *frag, ChunkInstanceT *fragChunk,
			     CIFragT *mate, ChunkInstanceT *mateChunk,
			     OrientType expectedOri){

  DistT *fragDist;
  LengthT separation;
  double jointstddev;
  OrientType edgeOri;

  GetRelationBetweenFragsOnChunks(frag,fragChunk,mate,mateChunk,&separation,&edgeOri);

  if(edgeOri != expectedOri) return FALSE;

  fragDist = GetDistT(ScaffoldGraph->Dists, frag->dist);

  jointstddev = sqrt ( separation.variance + fragDist->sigma * fragDist->sigma);

  if (abs (separation.mean - fragDist->mu) < MAX_SIGMA_SLOP * jointstddev) {
    return TRUE;
  } else {
    return FALSE;
  }

}



static
void InitCIFragTInChunkIterator(CGWFragIterator* frags,NodeCGW_T *chunk, int includeSurrogates){

  if(chunk->flags.bits.isScaffold){
    fprintf(stderr,"We haven't coded fragment iteration over scaffolds!\n");
    assert(0);
    exit(1);
    frags->isUtg = frags->isCtg = 0;
  }
  else if(chunk->flags.bits.isContig){
    frags->isUtg = 0;
    frags->isCtg = 1;
    frags->fOrder = 0;
  } else {
    frags->isCtg = 0;
    frags->isUtg = 1;
    frags->fOrder = 0;
  }

  frags->subchunkIterator = NULL;
  frags->fragiterma = NULL;

  frags->includeSurrogateFrgs = includeSurrogates;


  if(frags->isCtg && frags->includeSurrogateFrgs){

    // if a contig and we want surrogate frags, then we need to get
    // all frags out of all constituent unitigs,

    NodeCGW_T *ci;
    InitContigTIterator(ScaffoldGraph, chunk->id, TRUE, FALSE, &(frags->subchunks));
    ci = NextContigTIterator(&(frags->subchunks));
    assert(ci != NULL);
    if(frags->subchunkIterator == NULL){
      frags->subchunkIterator = (CGWFragIterator*) safe_malloc(sizeof(CGWFragIterator));
      frags->subchunkIterator->fragiterma=NULL;
      frags->subchunkIterator->subchunkIterator=NULL;
    }
    assert(frags->subchunkIterator != NULL);
    InitCIFragTInChunkIterator(frags->subchunkIterator,GetGraphNode(ScaffoldGraph->CIGraph,ci->id),frags->includeSurrogateFrgs);

  } else {

    // otherwise (either we have a unitig or we want only nonsurrogate
    // fragments), we can use the fragments already in the
    // multialignment

    frags->fragiterma = loadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, chunk->id, frags->isUtg);
  }

  frags->id = chunk->id;
  return;
}

static
int NextCIFragTInChunkIterator(CGWFragIterator* frags, CIFragT**nextfrg){

  if(frags->fragiterma!=NULL){

    assert(frags->isUtg || ! frags->includeSurrogateFrgs);

    if(frags->fOrder >= GetNumIntMultiPoss(frags->fragiterma->f_list)){
      *nextfrg=NULL;
      return FALSE;
    } else {
      CDS_CID_t fiid = GetIntMultiPos(frags->fragiterma->f_list,(frags->fOrder)++)->ident;
      *nextfrg = GetCIFragT(ScaffoldGraph->CIFrags,GetInfoByIID(ScaffoldGraph->iidToFragIndex,fiid)->fragIndex);
      return TRUE;
    }

  } else {
    int rv;
    CIFragT *retfrg;
    NodeCGW_T *ci;
    assert(frags->isCtg && frags->includeSurrogateFrgs);
    assert(frags->subchunkIterator!=NULL);
    while( (rv = NextCIFragTInChunkIterator(frags->subchunkIterator,&retfrg)) == FALSE){
      ci = NextContigTIterator(&(frags->subchunks));
      if(ci==NULL){
	*nextfrg=NULL;
	return FALSE;
      }
      assert(frags->subchunkIterator!=NULL);
      InitCIFragTInChunkIterator(frags->subchunkIterator,GetGraphNode(ScaffoldGraph->CIGraph,ci->id),frags->includeSurrogateFrgs);
    }
    *nextfrg = retfrg;
    return rv;

  }

}

static
void CleanupCIFragTInChunkIterator(CGWFragIterator* frags){
  if(frags->isCtg && frags->includeSurrogateFrgs)
    CleanupCIFragTInChunkIterator(frags->subchunkIterator);

  frags->subchunkIterator=NULL;
  frags->fragiterma=NULL;

  return;
}




static
void InitCGWMateIterator(CGWMateIterator* mates,CDS_CID_t fragIID, int external_only,ChunkInstanceT *node){

  CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags,
			     GetInfoByIID( ScaffoldGraph->iidToFragIndex,
					   fragIID)->fragIndex);

  mates->thisFragIID = fragIID;

  /* set the iterator up to fail next time around ... this will be revised
     if we find evidence of relevant links below */
  mates->nextLink = NULLINDEX;

  /* If this fragment has no constraints... continue */
  if(frag->flags.bits.hasMate == 0){
    assert(frag->mateOf == NULLINDEX);
    return;
  }

  //  An early OBT bug left in mates to reads that are deleted, so
  //  it's possibly to have links but no fragment.  Check and warn
  //  when this happens.
  //
  if (frag->mateOf == NULLINDEX) {
    if (frag->flags.bits.hasMate > 0)
      fprintf(stderr, "InitCGWMateIterator()-- WARNING!  Fragment "F_IID" has no mate, but still has a matelink!\n",
              fragIID);
    return;
  }


  /* Determine whether we want to know about all mates or only those outside the node
     of interest (asserted below to be one containing the fragment) */
  mates->external_only = external_only;

  if(mates->external_only){

    assert(node!=NULL);
    mates->node = node;

    // If this fragment only has links to fragments within the node of interesst ...continue
    if(node->flags.bits.isContig){
      if(frag->contigID != node->id){
	assert(0);
      }

      if( frag->flags.bits.hasInternalOnlyContigLinks){
	// return in state where iterator will fail ...
	return;
      }

    }else if(node->flags.bits.isCI){
      // If this fragment only has links to fragments within this CI...continue
      if(frag->cid != node->id){
	assert(0);
      }

      if( frag->flags.bits.hasInternalOnlyCILinks){
	// return in state where iterator will fail ...
	return;
      }
    }

  } else {
    assert(node==NULL);
    mates->node = NULL;
  }

  // at this point, we know there is at least one mate we need to
  // check ...  we will check external_only satisfaction in the
  // iterator itself, so we just want to set mates->nextLink and, as
  // necessary, the GateKeeper iterator

  //  The above comment ("at least one mate", "iterator") refers to
  //  old code where a read could have more than one mate.  Code now
  //  has at most one mate.

  mates->nextLink = GetCIFragT(ScaffoldGraph->CIFrags,frag->mateOf)->iid;

  return;

}

// determine whether a fragment is internal to a node
static
int check_internal(NodeCGW_T *node,CDS_CID_t frgIID){

  CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags,
			     GetInfoByIID(ScaffoldGraph->iidToFragIndex,
					  frgIID)->fragIndex);
  AssertPtr(frag);
  if(node->flags.bits.isScaffold){
    if(GetGraphNode(ScaffoldGraph->ContigGraph,frag->contigID)->scaffoldID ==
       node->id)
      return TRUE;
  } else if(node->flags.bits.isContig){
    if(node->id == frag->contigID)
      return TRUE;
  } else {
    assert(node->flags.bits.isCI);
    if(node->id == frag->cid || node->id == frag->CIid)
      return TRUE;
  }
  return FALSE;
}

// return the next link (external only as appropriate
static
int NextCGWMateIterator(CGWMateIterator* mates,CDS_CID_t * linkIID){

  int isInternal;

  if(mates->nextLink == NULLINDEX){
    return FALSE;
  }

  *linkIID = mates->nextLink;

  // there aren't any more mates, so set up to fail next time around ...
  mates->nextLink = NULLINDEX;


  if(mates->external_only)
    isInternal = check_internal(mates->node,*linkIID);
  else
    isInternal = FALSE;

  // if internal and we only want external, skip this one by recursing ...
  if(isInternal)
    return NextCGWMateIterator(mates,linkIID);

  else
    return TRUE;
}


static
int scaffoldOf(CDS_CID_t fiid){
  CIFragT *frg = GetCIFragT(ScaffoldGraph->CIFrags,
			    GetInfoByIID(ScaffoldGraph->iidToFragIndex,fiid)->fragIndex);
  NodeCGW_T *ci;
  if(frg->contigID != NULLINDEX){
    ci = GetGraphNode(ScaffoldGraph->ContigGraph,frg->contigID);
    assert(ci!=NULL);
    return ci->scaffoldID;
  }
  assert(frg->CIid!=NULLINDEX);
  ci = GetGraphNode(ScaffoldGraph->CIGraph,frg->CIid);
  assert(ci!=NULL);
  assert(ci->flags.bits.isSurrogate != TRUE);
  return ci->scaffoldID;
}

int matePlacedIn(CIFragT *frg, CDS_CID_t sid){
  CGWMateIterator mates;
  CDS_CID_t linkIID;
  InitCGWMateIterator(&mates,frg->iid,ALL_MATES,NULL);
  while(NextCGWMateIterator(&mates,&linkIID)){
    if(sid == scaffoldOf(linkIID)) return TRUE;
  }
  return FALSE;
}

static
int matePlacedOnlyIn(CIFragT *frg, CDS_CID_t sid, CIFragT **mate, ChunkInstanceT **mateChunk){
  CGWMateIterator mates;
  CDS_CID_t linkIID, mateiid;
  CDS_CID_t placedIn = NULLINDEX;
  InitCGWMateIterator(&mates,frg->iid,ALL_MATES,NULL);
  while(NextCGWMateIterator(&mates,&linkIID)){
    CDS_CID_t place = scaffoldOf(linkIID);
    if(place!=NULLINDEX){
      if(placedIn!=NULLINDEX){
	if(place!=placedIn){
	  *mate = NULL;
	  *mateChunk = NULL;
	  return FALSE;
	}
      } else {
	placedIn=place;
	mateiid = linkIID;
      }
    }
  }
  if(sid == placedIn){
    *mate = GetCIFragT(ScaffoldGraph->CIFrags,
		       GetInfoByIID(ScaffoldGraph->iidToFragIndex,mateiid)->fragIndex);
    assert((*mate)->cid==(*mate)->CIid);
    *mateChunk = GetGraphNode(ScaffoldGraph->CIGraph,(*mate)->cid);
    return TRUE;
  } else {
    *mate=NULL;
    *mateChunk=NULL;
    return FALSE;
  }
}




/* given an existing chunk and a list of IntMultiPoss, update the SeqDB entry for the chunk to add the IMPs, and update the consensus sequence */
static
void PlaceFragmentsInMultiAlignT(CDS_CID_t toID, int isUnitig,
				 VA_TYPE(IntMultiPos) *f_list){

  MultiAlignT *ma;

  //     1. get the old multialign from the seqDB
  ma =  loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, toID, isUnitig);
  //     2. add fragments to the f_list

  if (ma == NULL)
      return; // probably not good

  ConcatVA_IntMultiPos(ma->f_list,f_list);

  /* It might be a good idea to recompute consensus! */
  if(! isUnitig){
    //     2'. construct an ICM or IUM containing the new fragments
    //     2''. run consensus on it
    //     2'''. convert the returned ICM or IUM back to a multialignment
  }

  //     3. update the multialign
  updateMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB,toID,isUnitig,ma,TRUE);
}



/*
   Assign a subset of the fragments in an unresolved CI to one of its surrogates.
   The fragments listed are marked for membership (via their CIid field) int he new element
*/

static VA_TYPE(IntMultiPos) *f_list_CI = NULL;
static VA_TYPE(IntMultiPos) *f_list_Contig = NULL;

static
void ReallyAssignFragsToResolvedCI(GraphCGW_T *graph,
				   CDS_CID_t fromID, CDS_CID_t toID,
				   VA_TYPE(CDS_CID_t) *fragments){
  int i;
  int32 numFrags = GetNumCDS_CID_ts(fragments);
  NodeCGW_T *fromCI = GetGraphNode(graph, fromID);
  NodeCGW_T *toCI = GetGraphNode(graph, toID);
  ContigT *toContig = GetGraphNode(ScaffoldGraph->ContigGraph, toCI->info.CI.contigID);

  CDS_COORD_t surrogateAOffset = toCI->offsetAEnd.mean;
  CDS_COORD_t surrogateBOffset = toCI->offsetBEnd.mean;
  int flipped = (surrogateAOffset > surrogateBOffset);
  IntMultiPos fragPos;
  assert(graph->type == CI_GRAPH);
  assert(fromCI->type == UNRESOLVEDCHUNK_CGW);
  assert(toCI->type == RESOLVEDREPEATCHUNK_CGW);
  assert(toCI->scaffoldID != NULLINDEX);

  if(f_list_CI){
    ResetVA_IntMultiPos(f_list_CI);
    ResetVA_IntMultiPos(f_list_Contig);
  }else{
    f_list_CI = CreateVA_IntMultiPos(GetNumCDS_CID_ts(fragments));
    f_list_Contig = CreateVA_IntMultiPos(GetNumCDS_CID_ts(fragments));
  }

  fragPos.delta = NULL;
  fragPos.delta_length = 0;


  /* Check that the fragments CIid == node->id ! */
  /* Then, assign it to the new node and create a degenerate MultiAlignT */
  for(i = 0; i < numFrags; i++){
    CDS_CID_t fragID = *GetCDS_CID_t(fragments,i);
    int32 frgIdx = GetInfoByIID(ScaffoldGraph->iidToFragIndex,fragID)->fragIndex;
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, frgIdx);

    assert(frag->CIid == fromID);
    assert(frag->cid == fromID);
    frag->CIid = toID;              // Assign the fragment to the surrogate
    frag->contigID = toContig->id;  // Assign the fragment to the contig

    fragPos.type  = frag->type;
    fragPos.ident = fragID;
    fragPos.sourceInt = frgIdx;
    fragPos.position.bgn = frag->offset5p.mean;
    fragPos.position.end = frag->offset3p.mean;
    AppendIntMultiPos(f_list_CI, &fragPos);

    // Now figure out fragment position in target contig
    if(flipped){
      fragPos.position.bgn = surrogateAOffset - frag->offset5p.mean;
      fragPos.position.end = surrogateAOffset - frag->offset3p.mean;
      if(fragPos.position.end<0){
	fprintf(stderr,"DIRE WARNING: fragment " F_CID " left end rescued from < 0 coord on ctg " F_CID "\n",
		fragPos.ident,toContig->id);
	fragPos.position.end=0;
      }
      if(fragPos.position.bgn<0){
	fprintf(stderr,"DIRE WARNING: fragment " F_CID " left end rescued from < 0 coord on ctg " F_CID "\n",
		fragPos.ident,toContig->id);
	fragPos.position.bgn=0;
      }
    }else{
      fragPos.position.bgn = surrogateAOffset + frag->offset5p.mean;
      fragPos.position.end = surrogateAOffset + frag->offset3p.mean;
      if(fragPos.position.bgn<0){
	fprintf(stderr,"DIRE WARNING: fragment " F_CID " left end rescued from < 0 coord on ctg " F_CID "\n",
		fragPos.ident,toContig->id);
	fragPos.position.bgn=0;
      }
      if(fragPos.position.end<0){
	fprintf(stderr,"DIRE WARNING: fragment " F_CID " left end rescued from < 0 coord on ctg " F_CID "\n",
		fragPos.ident,toContig->id);
	fragPos.position.end=0;
      }
    }
    // We shouldn't need this!
    frag->contigOffset5p.variance = ComputeFudgeVariance(fragPos.position.bgn);
    frag->contigOffset5p.mean = fragPos.position.bgn;
    frag->contigOffset3p.variance = ComputeFudgeVariance(fragPos.position.end);
    frag->contigOffset3p.mean = fragPos.position.end;

    AppendIntMultiPos(f_list_Contig, &fragPos);
  }


  /* Copy IntMultiPos records from the source to destination CI, adjusting consensus sequence */
  PlaceFragmentsInMultiAlignT(toID, TRUE, f_list_CI);
  UpdateNodeFragments(ScaffoldGraph->CIGraph, toID, FALSE,FALSE);

  /* Copy IntMultiPos records to destination Contig, adjusting consensus sequence */
  PlaceFragmentsInMultiAlignT(toContig->id, FALSE, f_list_Contig);
  UpdateNodeFragments(ScaffoldGraph->ContigGraph, toContig->id,FALSE,FALSE);

  UpdateNodeUnitigs(loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, toContig->id, FALSE),
                    toContig);

  /* Do not Rebuild the Mate Edges of the target CI to reflect the changes in fragment membership */

  /* Do NOT Rebuild the Mate Edges of the target CI's Contig to reflect the changes in fragment membership.
     We will rebuild all mate edges when all fragments have been placed */


}






static
int
getChunkInstanceID(ChunkInstanceT *chunk, int index) {

  // chunk is not a surrogate -- just return chunk's id
  if (chunk->info.CI.numInstances == 0) {
    assert(index == 0);
    return(chunk->id);
  }

  // chunk is a surrogate

  if (chunk->info.CI.numInstances <= 2  && (index == 0))
    return(chunk->info.CI.instances.in_line.instance1);

  if (chunk->info.CI.numInstances == 2 && (index == 1))
    return(chunk->info.CI.instances.in_line.instance2);

  if (index < chunk->info.CI.numInstances)
    return(*(int32 *) Getint32(chunk->info.CI.instances.va, index));

  assert(index < chunk->info.CI.numInstances);

  return(-1);
}



void
resolveSurrogates(int    placeAllFragsInSinglePlacedSurros,
                  double cutoffToInferSingleCopyStatus) {
  int                i;
  GraphNodeIterator  CIGraphIterator;
  ChunkInstanceT    *parentChunk;
  int totalNumParentFrags=0;
  int numReallyPlaced=0;
  int                    allocedImpLists = 128;
  VA_TYPE(IntMultiPos) **impLists = NULL;

  assert((cutoffToInferSingleCopyStatus >= 0.0) &&
         (cutoffToInferSingleCopyStatus <= 1.0));

  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while ((parentChunk = NextGraphNodeIterator(&CIGraphIterator)) != NULL)
    if (allocedImpLists < parentChunk->info.CI.numInstances)
      allocedImpLists = parentChunk->info.CI.numInstances + 128;

  impLists = (VA_TYPE(IntMultiPos)**) safe_malloc(allocedImpLists*sizeof(VA_TYPE(IntMultiPos)*));
  for (i=0;i<allocedImpLists;i++)
    impLists[i] = CreateVA_IntMultiPos(20);

  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while ((parentChunk = NextGraphNodeIterator(&CIGraphIterator)) != NULL) {
    int numFrgsToPlace=0;
    HashTable_AS *fHash;
    int numInstances = parentChunk->info.CI.numInstances;
    int i,index, numFragmentsInParent;

    assert(numInstances <= allocedImpLists);

    // if numInstances >= 1 then it has a surrogate
    if(numInstances==0)
      continue;

    // count fragments and positions
    {
      MultiAlignT *maParent = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, parentChunk->id, TRUE);
      numFragmentsInParent = GetNumIntMultiPoss(maParent->f_list);
    }

    totalNumParentFrags += numFragmentsInParent;

    for(i=0;i<numInstances;i++){
      ChunkInstanceT *candidateChunk;
      CDS_CID_t sid,ctgiid;
      CGWFragIterator frags;
      CIFragT *nextfrg;

      index = getChunkInstanceID(parentChunk, i);
      candidateChunk = GetGraphNode(ScaffoldGraph->CIGraph, index);
      AssertPtr (candidateChunk);

      //  These were historically problems that were not asserts, but
      //  would just skip this instance.
      assert(parentChunk->type == UNRESOLVEDCHUNK_CGW);
      assert(candidateChunk->type == RESOLVEDREPEATCHUNK_CGW);
      assert(parentChunk != candidateChunk);

      if(candidateChunk->info.CI.baseID != parentChunk->id){
	if(candidateChunk==parentChunk){
	  fprintf(stderr,"resolveSurrogates: instance == parent for " F_CID " instance %d\n",
		  parentChunk->id,i);
	  continue;
	} else {
	  assert(candidateChunk->info.CI.baseID == parentChunk->id);
	}
      }

      sid=candidateChunk->scaffoldID;
      ctgiid = candidateChunk->info.CI.contigID;

      if(sid==NULLINDEX)
        continue;

      // now loop over fragments looking for those with mates in the scaffold
      InitCIFragTInChunkIterator(&frags,parentChunk,FALSE);
      while(NextCIFragTInChunkIterator(&frags, &nextfrg)){
        CIFragT *mate;
        ChunkInstanceT *mateChunk;
        int fragIsGood = 0;

        if (placeAllFragsInSinglePlacedSurros && numInstances == 1) {
          fragIsGood=1;
        } else {
          if(matePlacedOnlyIn(nextfrg,sid,&mate,&mateChunk)){
            assert(nextfrg->flags.bits.innieMate);
            if(FragAndMateAreCompatible(nextfrg,candidateChunk,mate,mateChunk,AS_INNIE)){
              fragIsGood= 1;
            }
          }
        }
        if(fragIsGood){
          // we're hot to trot ... now do something!
          IntMultiPos imp;
          imp.type = nextfrg->type;
          imp.ident = nextfrg->iid;
          imp.position.bgn = nextfrg->offset5p.mean;
          imp.position.end = nextfrg->offset3p.mean;
          imp.contained = 0; /* this might be wrong! */
          imp.delta_length=0;
          imp.delta=NULL;
          AppendVA_IntMultiPos(impLists[i],&imp);

          numFrgsToPlace++;
        }
      }
      CleanupCIFragTInChunkIterator(&frags);
    }  //  Over all instances


    if(numFrgsToPlace==0)
      continue;

    fHash = CreateScalarHashTable_AS(numFrgsToPlace);

    for(i=0;i<numInstances;i++){
      int j, numToPlace = GetNumIntMultiPoss(impLists[i]);
      for(j=0;j<numToPlace;j++){
	CDS_CID_t iid = GetIntMultiPos(impLists[i],j)->ident;
        ReplaceInHashTable_AS(fHash, iid, 0,
                              LookupValueInHashTable_AS(fHash,iid,0) + 1, 0);
      }
    }

    for(i=0;i<numInstances;i++){
      int j, numToPlace = GetNumIntMultiPoss(impLists[i]);
      VA_TYPE(CDS_CID_t) *toplace;

      if(numToPlace==0)
        continue;

      toplace = CreateVA_CDS_CID_t(numToPlace);

      //  Build the list of fragments to place
      for(j=0;j<numToPlace;j++){
	CDS_CID_t iid = GetIntMultiPos(impLists[i],j)->ident;

	int32 count_so_far = (int32)LookupValueInHashTable_AS(fHash,(uint64)iid,0);
	assert(count_so_far>0);
	if(count_so_far>1)
	  continue;

	AppendVA_CDS_CID_t(toplace,&iid);
      }


      // now, second-guess ourselves: if a sufficient fraction of
      // reads can be placed with mates, then place all fragments
      //
      if ((numInstances == 1) &&
          (GetNumCDS_CID_ts(toplace) > cutoffToInferSingleCopyStatus * numFragmentsInParent)) {

        CGWFragIterator frags;
        CIFragT *nextfrg;

        InitCIFragTInChunkIterator(&frags,parentChunk,FALSE);

        ResetVA_CDS_CID_t(toplace);

        while(NextCIFragTInChunkIterator(&frags, &nextfrg))
          AppendVA_CDS_CID_t(toplace,&(nextfrg->iid));

        CleanupCIFragTInChunkIterator(&frags);
      }

      // now really do the placement
      ReallyAssignFragsToResolvedCI(ScaffoldGraph->CIGraph,
                                    parentChunk->id,
                                    getChunkInstanceID(parentChunk,i),
                                    toplace);

      numReallyPlaced+=GetNumCDS_CID_ts(toplace);

      ResetVA_CDS_CID_t(toplace);

      ResetVA_IntMultiPos(impLists[i]);
    }

    DeleteHashTable_AS(fHash);
  }

  for (i=0;i<allocedImpLists;i++)
    DeleteVA_IntMultiPos(impLists[i]);

  fprintf(GlobalData->stderrc, "Placed %d surrogate fragments out of %d surrogate fragments\n", numReallyPlaced, totalNumParentFrags);
}



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

static const char *rcsid = "$Id: fragmentPlacement.c,v 1.37 2010-08-24 15:02:38 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
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

#undef  DEBUG_RS


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
				     LengthT*separation,PairOrient* edgeOri){


  // we pass in chunks because if a fragment in question is in a surrogate, we can't
  // get the info we need by identifying the parent chunk containing the frag

  SequenceOrient fOri,mOri,fCIOri,mCIOri,fCtgOri,mCtgOri,fOriOnScf,mOriOnScf;
  double f5pOnCI,m5pOnCI,f5pOnCtg,m5pOnCtg,f5pOnScf,m5pOnScf;
  NodeCGW_T *fragCtg, *mateCtg;

  assert(fragChunk->flags.bits.isCI);
  assert(mateChunk->flags.bits.isCI);


  f5pOnCI =  frag->offset5p.mean;
  fOri.setIsForward(frag->offset5p.mean < frag->offset3p.mean);

  fCIOri.setIsForward(fragChunk->offsetAEnd.mean < fragChunk->offsetBEnd.mean);
  f5pOnCtg =  fragChunk->offsetAEnd.mean + ((fCIOri.isForward()) ? f5pOnCI : -f5pOnCI);

  fragCtg = GetGraphNode(ScaffoldGraph->ContigGraph,fragChunk->info.CI.contigID);
  fCtgOri.setIsForward(fragCtg->offsetAEnd.mean < fragCtg->offsetBEnd.mean);
  f5pOnScf = fragCtg->offsetAEnd.mean + ((fCtgOri.isForward()) ? f5pOnCtg : -f5pOnCtg);


  m5pOnCI =  mate->offset5p.mean;
  mOri.setIsForward(mate->offset5p.mean < mate->offset3p.mean);

  mCIOri.setIsForward(mateChunk->offsetAEnd.mean < mateChunk->offsetBEnd.mean);
  m5pOnCtg =  mateChunk->offsetAEnd.mean + ((mCIOri.isForward()) ? m5pOnCI : -m5pOnCI);

  mateCtg = GetGraphNode(ScaffoldGraph->ContigGraph,mateChunk->info.CI.contigID);
  mCtgOri.setIsForward(mateCtg->offsetAEnd.mean < mateCtg->offsetBEnd.mean);
  m5pOnScf = mateCtg->offsetAEnd.mean + ((mCtgOri.isForward()) ? m5pOnCtg : -m5pOnCtg);


  separation->mean = abs(f5pOnScf - m5pOnScf);

  separation->variance =
    abs(fragCtg->offsetAEnd.variance - mateCtg->offsetAEnd.variance) + ComputeFudgeVariance(f5pOnCtg) + ComputeFudgeVariance(m5pOnCtg);

  if( ((fOri.isForward())+(fCIOri.isForward())+(fCtgOri.isForward())) % 2 == 0 ){
    fOriOnScf.setIsReverse();
  } else {
    fOriOnScf.setIsForward();
  }

  if( ((mOri.isForward())+(mCIOri.isForward())+(mCtgOri.isForward())) % 2 == 0 ){
    mOriOnScf.setIsReverse();
  } else {
    mOriOnScf.setIsForward();
  }

  if(fOriOnScf.isForward()){
    if(mOriOnScf.isReverse()){
      if(f5pOnScf < m5pOnScf){
	edgeOri->setIsAB_BA();
      } else {
	edgeOri->setIsBA_AB();
      }
    } else {
      edgeOri->setIsAB_AB();
    }
  } else { // fOriOnScf == B_A
    if(mOriOnScf.isForward()){
      if(f5pOnScf < m5pOnScf){
	edgeOri->setIsBA_AB();
      } else {
	edgeOri->setIsAB_BA();
      }
    } else {
      edgeOri->setIsBA_BA();
    }
  }

}


// FragAndMateAreCompatible() is based partly on PlaceFragmentBasedOnMate()
// from AS_REZ/PlaceFragments.c

static
int FragAndMateAreCompatible(CIFragT *frag, ChunkInstanceT *fragChunk,
			     CIFragT *mate, ChunkInstanceT *mateChunk){

  DistT *fragDist;
  LengthT separation;
  double jointstddev;
  PairOrient edgeOri;

  GetRelationBetweenFragsOnChunks(frag,fragChunk,mate,mateChunk,&separation,&edgeOri);

  //  This should really come from the library.  For now, we'll just assume innie.
#warning INNIE mate pairs assumed
  if(edgeOri.isInnie() == false)
    return FALSE;

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

    frags->fragiterma = ScaffoldGraph->tigStore->loadMultiAlign(chunk->id, frags->isUtg);
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
      *nextfrg = GetCIFragT(ScaffoldGraph->CIFrags, fiid);
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

  CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, fragIID);

  mates->thisFragIID = fragIID;

  /* set the iterator up to fail next time around ... this will be revised
     if we find evidence of relevant links below */
  mates->nextLink = NULLINDEX;

  /* If this fragment has no constraints... continue */
  if ((frag->flags.bits.hasMate == 0) ||
      (frag->mate_iid           == 0))
    return;

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

  mates->nextLink = frag->mate_iid;

  return;

}

// determine whether a fragment is internal to a node
static
int check_internal(NodeCGW_T *node,CDS_CID_t frgIID){

  CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, frgIID);
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
  CIFragT *frg = GetCIFragT(ScaffoldGraph->CIFrags, fiid);
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
  InitCGWMateIterator(&mates,frg->read_iid,ALL_MATES,NULL);
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
  InitCGWMateIterator(&mates,frg->read_iid,ALL_MATES,NULL);
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
    *mate = GetCIFragT(ScaffoldGraph->CIFrags, mateiid);
    assert((*mate)->cid==(*mate)->CIid);
    *mateChunk = GetGraphNode(ScaffoldGraph->CIGraph,(*mate)->cid);
    return TRUE;
  } else {
    *mate=NULL;
    *mateChunk=NULL;
    return FALSE;
  }
}





/*
   Assign a subset of the fragments in an unresolved CI to one of its surrogates.
   The fragments listed are marked for membership (via their CIid field) int he new element
*/

static
void ReallyAssignFragsToResolvedCI(CDS_CID_t frCIid,
                                   CDS_CID_t toCIid,
				   VA_TYPE(CDS_CID_t) *fragments){

  NodeCGW_T *frCI = GetGraphNode(ScaffoldGraph->CIGraph, frCIid);
  NodeCGW_T *toCI = GetGraphNode(ScaffoldGraph->CIGraph, toCIid);

  CDS_CID_t  toContigID = toCI->info.CI.contigID;
  ContigT   *toContig   = GetGraphNode(ScaffoldGraph->ContigGraph, toContigID);

  MultiAlignT  *utgMA = ScaffoldGraph->tigStore->loadMultiAlign(toCIid, TRUE);
  MultiAlignT  *ctgMA = ScaffoldGraph->tigStore->loadMultiAlign(toContigID, FALSE);

  assert(frCI->type == UNRESOLVEDCHUNK_CGW);
  assert(toCI->type == RESOLVEDREPEATCHUNK_CGW);

  assert(toCI->scaffoldID != NULLINDEX);

  for (int32 i=0; i<GetNumCDS_CID_ts(fragments); i++) {
    CDS_CID_t fragID  = *GetCDS_CID_t(fragments,i);
    CIFragT  *frag    = GetCIFragT(ScaffoldGraph->CIFrags, fragID);

    IntMultiPos fragPos;

#ifdef DEBUG_RS
    fprintf(stderr, "PLACE frag %d into instance %d\n", fragID, toCIid);
#endif

    //  Calling resolveSurrogates() multiple times (for closure reads, see AS_CGW_main.c) can result
    //  in placed surrogate reads placing other surrogate reads.
    //
    //  ----------surrogate----surrogate-------------
    //       A-------A              C-------C
    //                  B--------B
    //
    //  On the first pass, the right surrogate is examined first.  We cannot place fragment B, since
    //  the other end isn't placed.  C can be placed.  The left surrogate is examined, and placing A
    //  triggers the "second guess" to place all fragments in this surrogate.
    //
    //  On the second pass, B can now be placed in the right surrogate, since it is anchored in the
    //  left surrogate.  If that triggers the placement of all reads, we'll place C again.
    //
    //  So, we explicitly check that the fragment is already placed here.
    //
    if (frag->CIid != frCIid) {
      fprintf(stderr, "WARNING:  Fragment %d in surrogate CI %d, previously assigned to CI %d, now being (re)assigned to CI %d.\n",
              fragID, frCIid, frag->CIid, toCIid);

      IntMultiPos  *imp = GetIntMultiPos(utgMA->f_list, 0);

      for (int32 i=0; (imp != NULL) && (i<GetNumIntMultiPoss(utgMA->f_list)); i++) {
        if (imp[i].ident == fragID)
          imp = NULL;
      }

      //  If the imp is now NULL, we found the fragment already here, and can skip the rest.

      if (imp == NULL)
        continue;

      //  Otherwise, we did not find the fragment already here, thus, it is somewhere else, yet we
      //  want to put it here.  That's an error.

      fprintf(stderr, "ERROR:    Fragment %d in surrogate CI %d, previously assigned to CI %d, now being (re)assigned to CI %d.\n",
              fragID, frCIid, frag->CIid, toCIid);
      assert(frag->CIid == frCIid);
    }

    assert(frag->cid  == frCIid);

    frag->CIid            = toCIid;
    frag->contigID        = toContigID;

    //  Create a new fragment for the unitig.

    fragPos.type          = AS_READ;
    fragPos.ident         = fragID;
    fragPos.contained     = 0;
    fragPos.parent        = 0;
    fragPos.ahang         = 0;
    fragPos.bhang         = 0;
    fragPos.position.bgn  = frag->offset5p.mean;
    fragPos.position.end  = frag->offset3p.mean;
    fragPos.delta_length  = 0;
    fragPos.delta         = NULL;

    AppendIntMultiPos(utgMA->f_list, &fragPos);

    // Create a new fragment for the contig (by updating the placement of the fragment we just created)

    int32      surrogateAOffset = toCI->offsetAEnd.mean;
    int32      surrogateBOffset = toCI->offsetBEnd.mean;

    if (surrogateAOffset > surrogateBOffset) {
      fragPos.position.bgn = surrogateAOffset - frag->offset5p.mean;
      fragPos.position.end = surrogateAOffset - frag->offset3p.mean;
    } else {
      fragPos.position.bgn = surrogateAOffset + frag->offset5p.mean;
      fragPos.position.end = surrogateAOffset + frag->offset3p.mean;
    }

    if (fragPos.position.bgn < 0)
      fragPos.position.bgn = 0;
    if (fragPos.position.end < 0)
      fragPos.position.end = 0;

    // We shouldn't need this!
    frag->contigOffset5p.variance = ComputeFudgeVariance(fragPos.position.bgn);
    frag->contigOffset5p.mean     = fragPos.position.bgn;

    frag->contigOffset3p.variance = ComputeFudgeVariance(fragPos.position.end);
    frag->contigOffset3p.mean     = fragPos.position.end;

    AppendIntMultiPos(ctgMA->f_list, &fragPos);
  }

  //  Update.
  //
  //  Do NOT rebuild consensus sequences.
  //
  //  Do NOT Rebuild the Mate Edges of the target CI to reflect the changes in fragment membership
  //
  //  Do NOT Rebuild the Mate Edges of the target CI's Contig to reflect the changes in fragment
  //  membership.  We will rebuild all mate edges when all fragments have been placed

  ScaffoldGraph->tigStore->insertMultiAlign(utgMA, TRUE,  TRUE);
  ScaffoldGraph->tigStore->insertMultiAlign(ctgMA, FALSE, TRUE);

  UpdateNodeFragments(ScaffoldGraph->CIGraph,     toCIid,     FALSE, FALSE);
  UpdateNodeFragments(ScaffoldGraph->ContigGraph, toContigID, FALSE, FALSE);

  UpdateNodeUnitigs(ctgMA, toContig);
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

  assert(chunk->info.CI.numInstances == GetNumCDS_CID_ts(chunk->info.CI.instances.va));

  if (index < chunk->info.CI.numInstances)
    return(*(int32 *) Getint32(chunk->info.CI.instances.va, index));

  assert(index < chunk->info.CI.numInstances);

  return(-1);
}

int placedByClosureIn(int index, CDS_CID_t iid, CDS_CID_t sid, CDS_CID_t ctgiid, HashTable_AS *instanceList) {
   int result = FALSE;

   gkPlacement *gkpl = ScaffoldGraph->gkpStore->gkStore_getReadPlacement(iid);
   assert(gkpl);
   assert(gkpl->bound1);
   assert(gkpl->bound2);
   
   // get the reads indicated by the input line
   CIFragT *leftMate = GetCIFragT(ScaffoldGraph->CIFrags, gkpl->bound1); 
   CIFragT *rightMate = GetCIFragT(ScaffoldGraph->CIFrags, gkpl->bound2);

   // the reads aren't in contigs so there can't be gaps to fill
   if (leftMate->contigID == NULLINDEX || rightMate->contigID == NULLINDEX) {
      return result;
   }
   ChunkInstanceT * begin_chunk = GetGraphNode(ScaffoldGraph->CIGraph, leftMate->CIid);
   ChunkInstanceT * end_chunk   = GetGraphNode(ScaffoldGraph->CIGraph, rightMate->CIid);

   int instanceCount = 0;
   if (leftMate->contigID == ctgiid || rightMate->contigID == ctgiid) {
      if (begin_chunk->scaffoldID == sid && end_chunk->scaffoldID == sid) {
         // if our parent is already in this surrogate instance, place us here as well
         if (begin_chunk->id == index || end_chunk->id == index) {
            result = TRUE;
         } else {
            // try going off the aend of the leftNode
            while (begin_chunk->BEndNext != NULLINDEX && begin_chunk->BEndNext != end_chunk->id) {
               if (begin_chunk->BEndNext == index && instanceCount == 0) {
                  result = TRUE;
                  instanceCount++;
               }
               else if (ExistsInHashTable_AS(instanceList, begin_chunk->BEndNext, 0)) {
                  result = FALSE;
                  instanceCount++;
               } 
               begin_chunk = GetGraphNode(ScaffoldGraph->CIGraph, begin_chunk->BEndNext);
            }
            
            // if we ran off the end from the left and didn't find the surrogate instance, try from the right
            if (result && begin_chunk->BEndNext != NULLINDEX) { 
               // we stopped because we got to the end chunk or we found the instance we were looking for, no need to search other direction
            } else {
               result = FALSE;
               instanceCount = 0;
               begin_chunk = GetGraphNode(ScaffoldGraph->CIGraph, leftMate->CIid);
               while (end_chunk->BEndNext != NULLINDEX && end_chunk->BEndNext != begin_chunk->id) {
                  if (end_chunk->BEndNext == index && instanceCount == 0) {
                     result = TRUE;
                     instanceCount++;
                  }
                  else if (ExistsInHashTable_AS(instanceList, end_chunk->BEndNext, 0)) {
                     result = FALSE;
                     instanceCount++;
                  }
                  end_chunk = GetGraphNode(ScaffoldGraph->CIGraph, end_chunk->BEndNext);
               }
               if (end_chunk->BEndNext == NULLINDEX) {
                  result = FALSE;
               }
            }
         }
      }
   }

   return result;
}
      

void
resolveSurrogates(int    placeAllFragsInSinglePlacedSurros,
                  double cutoffToInferSingleCopyStatus) {
  int32 totalNumParentFrags  = 0;
  int32 numReallyPlaced      = 0;
  int32 allocatedPlacedLists   = 128;

  assert((cutoffToInferSingleCopyStatus >= 0.0) &&
         (cutoffToInferSingleCopyStatus <= 1.0));

  GraphNodeIterator  CIGraphIterator;
  ChunkInstanceT    *parentChunk;


  //  Allocate a list of IIDs for each possible surrogate instance.  We allocate enough lists for
  //  the largest number of instances in the whole assembly.

  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while ((parentChunk = NextGraphNodeIterator(&CIGraphIterator)) != NULL)
    if (allocatedPlacedLists < parentChunk->info.CI.numInstances)
      allocatedPlacedLists = parentChunk->info.CI.numInstances;

  VA_TYPE(int32) **placedList = (VA_TYPE(int32)**) safe_malloc(allocatedPlacedLists * sizeof(VA_TYPE(int32)*));

  for (int32 i=0; i<allocatedPlacedLists; i++)
    placedList[i] = CreateVA_int32(20);

  //  Over all unitigs, try to place surrogate fragments.

  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while ((parentChunk = NextGraphNodeIterator(&CIGraphIterator)) != NULL) {
    int32            numFragmentsInParent;
    int32            numFrgsToPlace = 0;
    int              numInstances = parentChunk->info.CI.numInstances;

    assert(numInstances <= allocatedPlacedLists);

    if (numInstances == 0)
      //  Not a surrogate
      continue;

    numFragmentsInParent  = ScaffoldGraph->tigStore->getNumFrags(parentChunk->id, TRUE);
    totalNumParentFrags  += numFragmentsInParent;

#ifdef DEBUG_RS
    fprintf(stderr, "----------------------------------------\n");
    fprintf(stderr, "RESOLVE for parent %d with %d fragments and %d instances.\n", parentChunk->id, numFragmentsInParent, numInstances);
#endif

    HashTable_AS *instanceList = CreateScalarHashTable_AS();

    //  Over all instances, do some sanity checking, clear the placedList, and remember the instance ID.

    for (int32 i=0; i<numInstances; i++) {
      int32           instanceID     = getChunkInstanceID(parentChunk, i);
      ChunkInstanceT *candidateChunk = GetGraphNode(ScaffoldGraph->CIGraph, instanceID);

      //  These were historically problems that were not asserts, but
      //  would just skip this instance.
      assert(parentChunk->type              == UNRESOLVEDCHUNK_CGW);
      assert(candidateChunk->type           == RESOLVEDREPEATCHUNK_CGW);
      assert(candidateChunk->info.CI.baseID == parentChunk->id);
      assert(parentChunk                    != candidateChunk);

      ResetVA_int32(placedList[i]);

      InsertInHashTable_AS(instanceList, (uint64)instanceID, 0, (uint64)1, 0);
    }

    //  Over all instances (again), build a list of the fragments we can place in each instance.

    HashTable_AS    *fragPlacedTimes = CreateScalarHashTable_AS();

    for (int32 i=0; i<numInstances; i++) {
      int32           instanceID     = getChunkInstanceID(parentChunk, i);
      ChunkInstanceT *candidateChunk = GetGraphNode(ScaffoldGraph->CIGraph, instanceID);

      int32       sid    = candidateChunk->scaffoldID;
      int32       ctgiid = candidateChunk->info.CI.contigID;

      if (sid == NULLINDEX)
        continue;

      CGWFragIterator frags;
      CIFragT        *frag;

      InitCIFragTInChunkIterator(&frags,parentChunk,FALSE);
      while(NextCIFragTInChunkIterator(&frags, &frag)){
        int fragIsGood = 0;

        if (frag->CIid != parentChunk->id)
          //  Placed in a previuos round
          continue;
        
        if ((placeAllFragsInSinglePlacedSurros) &&
            (numInstances == 1))
          fragIsGood = 1;

        CIFragT        *mate      = NULL;
        ChunkInstanceT *mateChunk = NULL;

        if ((matePlacedOnlyIn(frag,sid,&mate,&mateChunk)) &&
            (frag->flags.bits.innieMate) &&
            (FragAndMateAreCompatible(frag,candidateChunk,mate,mateChunk)))
          fragIsGood = 1;          

        // if this is closure read and we can place it in this location, do it
        if (ScaffoldGraph->gkpStore->gkStore_getFRGtoPLC(frag->read_iid) != 0 && 
            placedByClosureIn(instanceID, frag->read_iid, sid, ctgiid, instanceList))
          fragIsGood = 1;
        
        if (fragIsGood) {
#ifdef DEBUG_RS
          fprintf(stderr, "frag %d (CIid=%d cid=%d) from parent=%d to=%d\n", frag->read_iid, frag->CIid, frag->cid, parentChunk->id, instanceID);
#endif
          AppendVA_int32(placedList[i], &frag->read_iid);

          //  Count how many times we try to place this fragment
          ReplaceInHashTable_AS(fragPlacedTimes, frag->read_iid, 0,
                                LookupValueInHashTable_AS(fragPlacedTimes, frag->read_iid, 0) + 1, 0);

          numFrgsToPlace++;
        }
      }

      CleanupCIFragTInChunkIterator(&frags);
    }

    DeleteHashTable_AS(instanceList);

    if (numFrgsToPlace == 0) {
      //  Nothing to place.
      DeleteHashTable_AS(fragPlacedTimes);
      continue;
    }

    //  Over all instances, again, build the final list of fragments to place.  This discards any
    //  fragment we try to place more than once.

    for (int32 i=0; i<numInstances; i++) {
      int numToPlace = GetNumint32s(placedList[i]);

      if (numToPlace == 0)
        continue;

#ifdef DEBUG_RS
      fprintf(stderr, "PLACING for instance %d\n", getChunkInstanceID(parentChunk, i));
#endif

      VA_TYPE(CDS_CID_t) *toplace = CreateVA_CDS_CID_t(numToPlace);

      for (int32 j=0; j<numToPlace; j++){
	CDS_CID_t iid = *Getint32(placedList[i], j);

	int32 numTimes = (int32)LookupValueInHashTable_AS(fragPlacedTimes, (uint64)iid, 0);
	assert(numTimes > 0);

	if (numTimes > 1)
          //  Whoops!  Placed it more than once.
	  continue;

#ifdef DEBUG_RS
        fprintf(stderr, "frag %d is OK\n", iid);
#endif
	AppendVA_CDS_CID_t(toplace, &iid);
      }


      //  Second-guess ourselves: if a sufficient fraction of reads can be placed with mates, then
      //  place all fragments.
      //
      if ((numInstances == 1) &&
          (GetNumCDS_CID_ts(toplace) > cutoffToInferSingleCopyStatus * numFragmentsInParent)) {

        CGWFragIterator frags;
        CIFragT        *frag;

        InitCIFragTInChunkIterator(&frags,parentChunk,FALSE);

        ResetVA_CDS_CID_t(toplace);

        while (NextCIFragTInChunkIterator(&frags, &frag)) {
#ifdef DEBUG_RS
          fprintf(stderr, "frag %d (CIid=%d cid=%d) from parent=%d to=%d (placing all)\n", frag->read_iid, frag->CIid, frag->cid, parentChunk->id, getChunkInstanceID(parentChunk, i));
#endif

          AppendVA_CDS_CID_t(toplace, &frag->read_iid);
        }

        CleanupCIFragTInChunkIterator(&frags);
      }

      // now really do the placement
      ReallyAssignFragsToResolvedCI(parentChunk->id,
                                    getChunkInstanceID(parentChunk,i),
                                    toplace);

      numReallyPlaced+=GetNumCDS_CID_ts(toplace);

      DeleteVA_CDS_CID_t(toplace);

      ResetVA_int32(placedList[i]);
    }  //  Over all instances

    DeleteHashTable_AS(fragPlacedTimes);
  }  //  Over all unitigs

  for (int32 i=0; i<allocatedPlacedLists; i++)
    DeleteVA_int32(placedList[i]);
  safe_free(placedList);

  fprintf(stderr, "Placed %d surrogate fragments out of %d surrogate fragments\n", numReallyPlaced, totalNumParentFrags);
}


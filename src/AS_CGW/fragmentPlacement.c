
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
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
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

  jointstddev = sqrt ( separation.variance + fragDist->stddev * fragDist->stddev);
  
  if (abs (separation.mean - fragDist->mean) < MAX_SIGMA_SLOP * jointstddev) {
    return TRUE;
  } else {
    return FALSE;
  }

}




// stuff that was in localUnitigging.c initially



void InitCIFragTInChunkIterator(CGWFragIterator* frags,NodeCGW_T *chunk, int includeSurrogates){

  if(chunk->flags.bits.isScaffold){
    fprintf(stderr,"We haven't coded fragment iteration over scaffolds!\n");
    assert(0);
    exit(-99);
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
  frags->ma = NULL;

  frags->includeSurrogateFrgs = includeSurrogates;


  if(frags->isCtg && frags->includeSurrogateFrgs){

    // if a contig and we want surrogate frags, then we need to get all frags out of all constituent unitigs,

    NodeCGW_T *ci;
    InitContigTIterator(ScaffoldGraph, chunk->id, TRUE, FALSE, &(frags->subchunks));
    ci = NextContigTIterator(&(frags->subchunks));
    assert(ci != NULL);
    if(frags->subchunkIterator == NULL){
      frags->subchunkIterator = (CGWFragIterator*) malloc(sizeof(CGWFragIterator));
      frags->subchunkIterator->ma=NULL;
      frags->subchunkIterator->subchunkIterator=NULL;
    }
    assert(frags->subchunkIterator != NULL);
    InitCIFragTInChunkIterator(frags->subchunkIterator,GetGraphNode(ScaffoldGraph->CIGraph,ci->id),frags->includeSurrogateFrgs);

  } else {

    // otherwise (either we have a unitig or we want only nonsurrogate fragments), we can use the fragments already in the multialignment

    int rv;
    if(frags->ma==NULL){
      frags->ma = CreateEmptyMultiAlignT();
    } 
    rv = ReLoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, frags->ma, chunk->id, frags->isUtg);
    assert(rv==0);

  }

  frags->id = chunk->id;
  return;
}

int NextCIFragTInChunkIterator(CGWFragIterator* frags, CIFragT**nextfrg){

  if(frags->ma!=NULL){

    assert(frags->isUtg || ! frags->includeSurrogateFrgs);

    if(frags->fOrder >= GetNumIntMultiPoss(frags->ma->f_list)){
      *nextfrg=NULL;
      return FALSE;
    } else {
      CDS_CID_t fiid = GetIntMultiPos(frags->ma->f_list,(frags->fOrder)++)->ident;
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

void CleanupCIFragTInChunkIterator(CGWFragIterator* frags){
  if(frags->isCtg && frags->includeSurrogateFrgs){
    CleanupCIFragTInChunkIterator(frags->subchunkIterator);
    frags->subchunkIterator=NULL;
  } else {
    if(frags->ma!=NULL){
      DeleteMultiAlignT(frags->ma);
      frags->ma=NULL;
    }
  }
  return;
}




void InitCGWMateIterator(CGWMateIterator* mates,CDS_CID_t fragIID, int external_only,ChunkInstanceT *node){

  CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags,
			     GetInfoByIID( ScaffoldGraph->iidToFragIndex, 
					   fragIID)->fragIndex);
  int numLinks = frag->numLinks;

  mates->thisFragIID = fragIID;


  /* set the iterator up to fail next time around ... this will be revised
     if we find evidence of relevant links below */
  mates->nextLink = NULLINDEX;
  
  /* If this fragment has no constraints... continue */
  if(numLinks == 0){
    assert(frag->mateOf == NULLINDEX);
    return;
  }

  //  An early OBT bug left in mates to reads that are deleted, so
  //  it's possibly to have links but no fragment.  Check and warn
  //  when this happens.
  //
  if (frag->mateOf == NULLINDEX) {
    if (numLinks > 0)
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

  // at this point, we know there is at least one mate we need to check ...
  // we will check external_only satisfaction in the iterator itself, so
  // we just want to set mates->nextLink and, as necessary, the GateKeeper iterator

  mates->getLinksFromStore = frag->flags.bits.getLinksFromStore;

  if(! mates->getLinksFromStore){
    mates->nextLink = GetCIFragT(ScaffoldGraph->CIFrags,frag->mateOf)->iid;
  } else {
    assert(0);

#if 0
    GateKeeperLinkRecord GKPLink;
    int rv;
    assert(frag->linkHead != NULLINDEX);

    // set up the iterator
    CreateGateKeeperLinkRecordIterator(ScaffoldGraph->gkpStore.lnkStore,
				       frag->mateOf,
				       fragIID, &(mates->GKPLinks));
    // get the iid of the first linked fragment
    rv = NextGateKeeperLinkRecordIterator(&(mates->GKPLinks), &GKPLink);
    assert(rv==1);
    mates->nextLink = ( mates->thisFragIID == GKPLink.frag1 ? GKPLink.frag2 : GKPLink.frag1);
    assert(mates->nextLink>NULLINDEX);
#endif
  }

  return;
     
}

// determine whether a fragment is internal to a node
static int check_internal(NodeCGW_T *node,CDS_CID_t frgIID){

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
int NextCGWMateIterator(CGWMateIterator* mates,CDS_CID_t * linkIID){

  int isInternal;

  if(mates->nextLink == NULLINDEX){
    return FALSE;
  }

  *linkIID = mates->nextLink;

  if(!mates->getLinksFromStore){
    // there aren't any more mates, so set up to fail next time around ... 
    mates->nextLink = NULLINDEX;

  } else {
    assert(0);
#if 0
    // get the iid of the next link from the gatekeeper ... or set up to fail if no more
    GateKeeperLinkRecord GKPLink;
    int rv = NextGateKeeperLinkRecordIterator(&(mates->GKPLinks),&GKPLink);
    if(rv==FALSE){
      mates->nextLink = NULLINDEX;
    } else {
      mates->nextLink = ( mates->thisFragIID == GKPLink.frag1 ? GKPLink.frag2 : GKPLink.frag1);
    }
#endif
  }

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


int matePlacedIn(CIFragT *frg, CDS_CID_t sid){
  CGWMateIterator mates;
  CDS_CID_t linkIID;
  InitCGWMateIterator(&mates,frg->iid,ALL_MATES,NULL);
  while(NextCGWMateIterator(&mates,&linkIID)){
    if(sid == scaffoldOf(linkIID)) return TRUE;
  }
  return FALSE;
}

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


//VA_DEF(IntMultiPos);


/* given an existing chunk and a list of IntMultiPoss, update the SeqDB entry for the chunk to add the IMPs, and update the consensus sequence */
void PlaceFragmentsInMultiAlignT(CDS_CID_t toID, int isUnitig,
				 VA_TYPE(IntMultiPos) *f_list){
  
  MultiAlignT *ma;

  //     1. get the old multialign from the seqDB
  ma =  LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, toID, isUnitig);
  //     2. add fragments to the f_list

  ConcatVA_IntMultiPos(ma->f_list,f_list);

  /* It might be a good idea to recompute consensus! */
  if(! isUnitig){
    //     2'. construct an ICM or IUM containing the new fragments
    //     2''. run consensus on it
    //     2'''. convert the returned ICM or IUM back to a multialignment
  }  

  //     3. update the multialign
  UpdateMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB,toID,isUnitig,ma,FALSE);
  
  return;

}



/* 
   Assign a subset of the fragments in an unresolved CI to one of its surrogates.
   The fragments listed are marked for membership (via their CIid field) int he new element
*/

static VA_TYPE(IntMultiPos) *f_list_CI = NULL;
static VA_TYPE(IntMultiPos) *f_list_Contig = NULL;

void ReallyAssignFragsToResolvedCI(GraphCGW_T *graph,
				   CDS_CID_t fromID, CDS_CID_t toID,
				   VA_TYPE(CDS_CID_t) *fragments){
  int i;
  int32 numFrags = GetNumCDS_CID_ts(fragments);
  NodeCGW_T *fromCI = GetGraphNode(graph, fromID);
  NodeCGW_T *toCI = GetGraphNode(graph, toID);
  ContigT *toContig = GetGraphNode(ScaffoldGraph->ContigGraph, toCI->info.CI.contigID);

  MultiAlignT *toCIMA = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, toID, graph->type == CI_GRAPH);

  MultiAlignT *contigMA = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, toContig->id, FALSE);

  // GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, toContig->id);
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
    
    fragPos.type = AS_MSG_SafeConvert_charToFragType(frag->type,TRUE);
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
  //  PlaceFragmentsInMultiAlignT(toCIMA, f_list_CI);
  PlaceFragmentsInMultiAlignT(toID, TRUE, f_list_CI);
  UpdateNodeFragments(ScaffoldGraph->CIGraph, toID, FALSE,FALSE);
  
  /* Copy IntMultiPos records to destination Contig, adjusting consensus sequence */
  PlaceFragmentsInMultiAlignT(toContig->id, FALSE, f_list_Contig);
  UpdateNodeFragments(ScaffoldGraph->ContigGraph, toContig->id,FALSE,FALSE);
  UpdateNodeUnitigs(contigMA, toContig);
  
  
  /* Do not Rebuild the Mate Edges of the target CI to reflect the changes in fragment membership */
  
  /* Do NOT Rebuild the Mate Edges of the target CI's Contig to reflect the changes in fragment membership.
     We will rebuild all mate edges when all fragments have been placed */
  
  
}


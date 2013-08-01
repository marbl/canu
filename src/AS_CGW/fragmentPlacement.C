
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

static const char *rcsid = "$Id$";

#include "AS_global.H"

#include "fragmentPlacement.H"

#include "GraphCGW_T.H"
#include "AS_CGW_dataTypes.H"
#include "ScaffoldGraph_CGW.H"
#include "ScaffoldGraphIterator_CGW.H"
#include "Globals_CGW.H"
#include "ScaffoldGraph_CGW.H"
#include "Output_CGW.H"

#include <set>
#include <map>
#include <vector>

using namespace std;

#define MAX_SIGMA_SLOP 3.0

#define EXTERNAL_MATES_ONLY 1
#define ALL_MATES 0


typedef struct FragIterator{
  uint32               isUtg;
  uint32               isCtg;
  uint32               includeSurrogateFrgs;
  CDS_CID_t            id;
  MultiAlignT         *fragiterma;
  int32                 fOrder; /* which fragment within multialign is next */
  ContigTIterator       subchunks;
  struct FragIterator  *subchunkIterator;
} CGWFragIterator;


typedef struct  {
  CDS_CID_t  thisFragIID;
  CDS_CID_t  nextLink;
  NodeCGW_T *node;
  uint32     external_only;
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


  separation->mean = fabs(f5pOnScf - m5pOnScf);

  separation->variance =
    fabs(fragCtg->offsetAEnd.variance - mateCtg->offsetAEnd.variance) + ComputeFudgeVariance(f5pOnCtg) + ComputeFudgeVariance(m5pOnCtg);

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

  if (fabs (separation.mean - fragDist->mu) < MAX_SIGMA_SLOP * jointstddev) {
    return TRUE;
  } else {
    return FALSE;
  }

}



static
void
InitCIFragTInChunkIterator(CGWFragIterator* frags,NodeCGW_T *chunk, int includeSurrogates) {

  assert(chunk->flags.bits.isScaffold == 0);

  frags->isUtg  = (chunk->flags.bits.isContig == true) ? 0 : 1;
  frags->isCtg  = (chunk->flags.bits.isContig == true) ? 1 : 0;
  frags->fOrder = 0;

  frags->subchunkIterator = NULL;
  frags->fragiterma       = NULL;

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
    assert(frags->isCtg);
    assert(frags->includeSurrogateFrgs);
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
  CDS_CID_t linkIID=0, mateiid=0;
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





//  Assign a subset of the fragments in an unresolved CI to one of its surrogates.  The fragments
//  listed are marked for membership (via their CIid field) in the new element
//
static
void
ReallyAssignFragsToResolvedCI(CDS_CID_t           frCIid,
                              CDS_CID_t           toCIid,
                              vector<CDS_CID_t>  &fragments) {

  NodeCGW_T *frCI = GetGraphNode(ScaffoldGraph->CIGraph, frCIid);
  NodeCGW_T *toCI = GetGraphNode(ScaffoldGraph->CIGraph, toCIid);

  CDS_CID_t  toContigID = toCI->info.CI.contigID;
  ContigT   *toContig   = GetGraphNode(ScaffoldGraph->ContigGraph, toContigID);

  MultiAlignT  *utgMA = ScaffoldGraph->tigStore->loadMultiAlign(toCIid, TRUE);
  MultiAlignT  *ctgMA = ScaffoldGraph->tigStore->loadMultiAlign(toContigID, FALSE);

  assert(frCI->type == UNRESOLVEDCHUNK_CGW);
  assert(toCI->type == RESOLVEDREPEATCHUNK_CGW);

  assert(toCI->scaffoldID != NULLINDEX);

  for (uint32 i=0; i<fragments.size(); i++) {
    CDS_CID_t fragID  = fragments[i];
    CIFragT  *frag    = GetCIFragT(ScaffoldGraph->CIFrags, fragID);

    IntMultiPos fragPos;

    //fprintf(stderr, "PLACE frag %d into instance %d\n", fragID, toCIid);

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

    //  Create a new fragment for the contig (by updating the placement of the fragment we just created)

    int32      surrogateAOffset = toCI->offsetAEnd.mean;
    int32      surrogateBOffset = toCI->offsetBEnd.mean;
    int32      maxOffset;

    if (surrogateAOffset > surrogateBOffset) {
      //  AOffset is the upper coord.  Surrogate is backwards.
      maxOffset = surrogateAOffset;
      fragPos.position.bgn = surrogateAOffset - frag->offset5p.mean;
      fragPos.position.end = surrogateAOffset - frag->offset3p.mean;
    } else {
      //  AOffset is the lower coord.
      maxOffset = surrogateBOffset;
      fragPos.position.bgn = surrogateAOffset + frag->offset5p.mean;
      fragPos.position.end = surrogateAOffset + frag->offset3p.mean;
    }

    if (fragPos.position.bgn < 0)
      fragPos.position.bgn = 0;
    if (fragPos.position.end < 0)
      fragPos.position.end = 0;

    if (fragPos.position.bgn >= maxOffset)
      fragPos.position.bgn = maxOffset;
    if (fragPos.position.end >= maxOffset)
      fragPos.position.end = maxOffset;

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

static
int
placedByClosureIn(int index,
                  CDS_CID_t iid,
                  CDS_CID_t sid,
                  CDS_CID_t ctgiid,
                  set<CDS_CID_t> &instanceList) {
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
          else if (instanceList.count(begin_chunk->BEndNext) > 0) {
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
            else if (instanceList.count(end_chunk->BEndNext) > 0) {
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
  int32 totalNumParentFrags    = 0;
  int32 numReallyPlaced        = 0;
  int32 allocatedPlacedLists   = 128;

  assert(cutoffToInferSingleCopyStatus >= 0.0);
  assert(cutoffToInferSingleCopyStatus <= 1.0);

  GraphNodeIterator  CIs;
  ChunkInstanceT    *parentCI;

  //  Allocate a list of IIDs for each possible surrogate instance.  We allocate enough lists for
  //  the largest number of instances in the whole assembly.

  InitGraphNodeIterator(&CIs, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while ((parentCI = NextGraphNodeIterator(&CIs)) != NULL)
    if (allocatedPlacedLists < parentCI->info.CI.numInstances)
      allocatedPlacedLists = parentCI->info.CI.numInstances;

  fprintf(stderr, "resolveSurrogates()--  Allocating %d vectors.\n", allocatedPlacedLists);

  vector<CDS_CID_t>   *placedList = new vector<CDS_CID_t> [allocatedPlacedLists];

  //  Over all unitigs, try to place surrogate fragments.

  InitGraphNodeIterator(&CIs, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while ((parentCI = NextGraphNodeIterator(&CIs)) != NULL) {
    assert(parentCI->info.CI.numInstances <= allocatedPlacedLists);

    if (parentCI->info.CI.numInstances == 0)
      //  Not a surrogate
      continue;

    fprintf(stderr, "resolveSurrogates()--  Placing for surrogate CI %d with %d instances.\n",
            parentCI->id, parentCI->info.CI.numInstances);

    int32  numFrgToPlace = 0;
    int32  numFrgInParent  = ScaffoldGraph->tigStore->getNumFrags(parentCI->id, TRUE);

    totalNumParentFrags  += numFrgInParent;

    set<CDS_CID_t>        instanceList;
    map<CDS_CID_t,int32>  fragPlacedTimes;

    //  Over all instances, do some sanity checking, clear the placedList, and remember the instance ID.

    for (int32 i=0; i<parentCI->info.CI.numInstances; i++) {
      int32           instanceID     = getChunkInstanceID(parentCI, i);
      ChunkInstanceT *targetCI       = GetGraphNode(ScaffoldGraph->CIGraph, instanceID);

      //  These were historically problems that were not asserts, but
      //  would just skip this instance.
      assert(parentCI->type           == UNRESOLVEDCHUNK_CGW);
      assert(targetCI->type           == RESOLVEDREPEATCHUNK_CGW);
      assert(targetCI->info.CI.baseID == parentCI->id);
      assert(parentCI                 != targetCI);

      placedList[i].clear();

      instanceList.insert(instanceID);
    }

    //  Over all instances, build a list of the fragments we can place in each instance.

    for (int32 i=0; i<parentCI->info.CI.numInstances; i++) {
      int32           instanceID     = getChunkInstanceID(parentCI, i);
      ChunkInstanceT *targetCI = GetGraphNode(ScaffoldGraph->CIGraph, instanceID);

      if (targetCI->scaffoldID == NULLINDEX)
        continue;

      CGWFragIterator frags;
      CIFragT        *frag;

      InitCIFragTInChunkIterator(&frags, parentCI, FALSE);
      while (NextCIFragTInChunkIterator(&frags, &frag)) {
        bool            fragIsGood = false;
        CIFragT        *mate       = NULL;
        ChunkInstanceT *mateChunk  = NULL;

        if (frag->CIid != parentCI->id)
          //  Placed in a previuos round
          continue;
        
        if ((placeAllFragsInSinglePlacedSurros) &&
            (parentCI->info.CI.numInstances == 1))
          fragIsGood = true;

        if ((matePlacedOnlyIn(frag, targetCI->scaffoldID, &mate, &mateChunk)) &&
            (frag->flags.bits.innieMate) &&
            (FragAndMateAreCompatible(frag, targetCI, mate, mateChunk)))
          fragIsGood = true;

        // if this is closure read and we can place it in this location, do it
        if (ScaffoldGraph->gkpStore->gkStore_getFRGtoPLC(frag->read_iid) != 0 && 
            placedByClosureIn(instanceID, frag->read_iid, targetCI->scaffoldID, targetCI->info.CI.contigID, instanceList))
          fragIsGood = true;

        if (fragIsGood == false)
          continue;

        assert(frag->CIid == parentCI->id);
        assert(frag->cid  == parentCI->id);

        fprintf(stderr, "resolveSurrogates()--  frag %d from parent CI=%d to instance CI=%d in scaffold %d\n",
                frag->read_iid, parentCI->id, instanceID, targetCI->scaffoldID);

        placedList[i].push_back(frag->read_iid);
        fragPlacedTimes[frag->read_iid]++;

        numFrgToPlace++;
      }

      CleanupCIFragTInChunkIterator(&frags);
    }

    instanceList.clear();

    if (numFrgToPlace == 0)
      //  Nothing to place.
      continue;

    //  Over all instances, again, build the final list of fragments to place.  This discards any
    //  fragment we try to place more than once.

    for (int32 i=0; i<parentCI->info.CI.numInstances; i++) {
      if (placedList[i].size() == 0)
        continue;

      vector<CDS_CID_t>   toPlace;

      for (int32 j=0; j<placedList[i].size(); j++){
        CDS_CID_t iid = placedList[i][j];

        int32 numTimes = fragPlacedTimes[iid];
        assert(numTimes > 0);

        if (numTimes > 1)
          //  Whoops!  Placed it more than once.
          continue;

        toPlace.push_back(iid);
      }


      //  Second-guess ourselves: if a sufficient fraction of reads can be placed with mates, then
      //  place all fragments.
      //
      if ((parentCI->info.CI.numInstances == 1) &&
          (toPlace.size() > cutoffToInferSingleCopyStatus * numFrgInParent)) {
        CGWFragIterator frags;
        CIFragT        *frag;

        toPlace.clear();

        InitCIFragTInChunkIterator(&frags,parentCI,FALSE);

        while (NextCIFragTInChunkIterator(&frags, &frag))
          toPlace.push_back(frag->read_iid);

        CleanupCIFragTInChunkIterator(&frags);
      }

      ReallyAssignFragsToResolvedCI(parentCI->id, getChunkInstanceID(parentCI,i), toPlace);

      numReallyPlaced += toPlace.size();

      placedList[i].clear();
    }  //  Over all instances

    fragPlacedTimes.clear();
  }  //  Over all unitigs

  delete [] placedList;

  fprintf(stderr, "Placed %d surrogate fragments out of %d surrogate fragments\n",
          numReallyPlaced, totalNumParentFrags);
}


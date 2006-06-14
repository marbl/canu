
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

#include "ScaffoldGraphIterator_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "GraphCGW_T.h"
#include "MultiAlignment_CNS.h"
#include "fixZLFContigs.h"
#include "ChiSquareTest_CGW.h"
#include "Instrument_CGW.h"
#include "UtilsREZ.h"

void WriteIMPToFile(IntMultiPos * imp, FILE * fp)
{
  fwrite(&(imp->type), sizeof(imp->type), 1, fp);
  fwrite(&(imp->ident), sizeof(imp->ident), 1, fp);
  fwrite(&(imp->contained), sizeof(imp->contained), 1, fp);
  fwrite(&(imp->position), sizeof(imp->position), 1, fp);
  fwrite(&(imp->delta_length), sizeof(imp->delta_length), 1, fp);
  if(imp->delta_length > 0)
    fwrite(imp->delta, sizeof(int32), imp->delta_length, fp);
}
void ReadIMPFromFile(IntMultiPos * imp, FILE * fp)
{
  fread(&(imp->type), sizeof(imp->type), 1, fp);
  fread(&(imp->ident), sizeof(imp->ident), 1, fp);
  fread(&(imp->contained), sizeof(imp->contained), 1, fp);
#ifdef AS_ENABLE_SOURCE
  imp->sourceInt = -1;
#endif
  fread(&(imp->position), sizeof(imp->position), 1, fp);
  fread(&(imp->delta_length), sizeof(imp->delta_length), 1, fp);
  if(imp->delta_length > 0)
  {
    imp->delta = safe_malloc(imp->delta_length * sizeof(int32));
    fread(imp->delta, sizeof(int32), imp->delta_length, fp);
  }
  else
    imp->delta = NULL;
}


void WriteIMPsToFile(VA_TYPE(IntMultiPos) * imps, FILE * fp)
{
  int32 i;
  int32 numIMPs = (int32) GetNumVA_IntMultiPos(imps);
  fwrite(&(numIMPs), sizeof(numIMPs), 1, fp);
  for(i = 0; i < numIMPs; i++)
    WriteIMPToFile(GetVA_IntMultiPos(imps, i), fp);
}
VA_TYPE(IntMultiPos) * ReadIMPsFromFile(FILE * fp)
{
  VA_TYPE(IntMultiPos) * imps;
  int32 numIMPs;
  int32 i;
  
  fread(&numIMPs, sizeof(numIMPs), 1, fp);
  if(numIMPs > 0)
  {
    imps = CreateVA_IntMultiPos(numIMPs);
    for(i = 0; i < numIMPs; i++)
    {
      IntMultiPos imp;
      ReadIMPFromFile(&imp, fp);
      AppendVA_IntMultiPos(imps, &imp);
    }
  }
  else
    imps = CreateVA_IntMultiPos(1);
  return imps;
}



void WriteIUPToFile(IntUnitigPos * iup, FILE * fp)
{
  fwrite(&(iup->type), sizeof(iup->type), 1, fp);
  fwrite(&(iup->ident), sizeof(iup->ident), 1, fp);
  fwrite(&(iup->position), sizeof(iup->position), 1, fp);
  fwrite(&(iup->delta_length), sizeof(iup->delta_length), 1, fp);
  if(iup->delta_length > 0)
    fwrite(iup->delta, sizeof(int32), iup->delta_length, fp);
}
void ReadIUPFromFile(IntUnitigPos * iup, FILE * fp)
{
  fread(&(iup->type), sizeof(iup->type), 1, fp);
  fread(&(iup->ident), sizeof(iup->ident), 1, fp);
  fread(&(iup->position), sizeof(iup->position), 1, fp);
  fread(&(iup->delta_length), sizeof(iup->delta_length), 1, fp);
  if(iup->delta_length > 0)
  {
    iup->delta = safe_malloc(iup->delta_length * sizeof(int32));
    fread(iup->delta, sizeof(int32), iup->delta_length, fp);
  }
  else
    iup->delta = NULL;
}


void WriteIUPsToFile(VA_TYPE(IntUnitigPos) * iups, FILE * fp)
{
  int32 i;
  int32 numIUPs = (int32) GetNumVA_IntUnitigPos(iups);
  fwrite(&(numIUPs), sizeof(numIUPs), 1, fp);
  for(i = 0; i < numIUPs; i++)
    WriteIUPToFile(GetVA_IntUnitigPos(iups, i), fp);
}
VA_TYPE(IntUnitigPos) * ReadIUPsFromFile(FILE * fp)
{
  VA_TYPE(IntUnitigPos) * iups;
  int32 numIUPs;
  int32 i;
  
  fread(&numIUPs, sizeof(numIUPs), 1, fp);
  if(numIUPs > 0)
  {
    iups = CreateVA_IntUnitigPos(numIUPs);
    for(i = 0; i < numIUPs; i++)
    {
      IntUnitigPos iup;
      ReadIUPFromFile(&iup, fp);
      AppendVA_IntUnitigPos(iups, &iup);
    }
  }
  else
    iups = CreateVA_IntUnitigPos(1);
  
  return iups;
}


void WriteMAToFile(MultiAlignT * ma, FILE * fp)
{
  fwrite(&(ma->forced), sizeof(ma->forced), 1, fp);
  fwrite(&(ma->id), sizeof(ma->id), 1, fp);
  fwrite(&(ma->refCnt), sizeof(ma->refCnt), 1, fp);
  fwrite(&(ma->source_alloc), sizeof(ma->source_alloc), 1, fp);

  CopyToFileVA_char(ma->consensus, fp);
  CopyToFileVA_char(ma->quality, fp);
  CopyToFileVA_int32(ma->delta, fp);
  WriteIMPsToFile(ma->f_list, fp);
  CopyToFileVA_int32(ma->udelta, fp);
  WriteIUPsToFile(ma->u_list, fp);
}
int ReadMAFromFile(MultiAlignT * ma, FILE * fp)
{
  if(fread(&(ma->forced), sizeof(ma->forced), 1, fp) != 1)
    return 0;
  
  if(fread(&(ma->id), sizeof(ma->id), 1, fp) != 1)
    return 0;
  
  if(fread(&(ma->refCnt), sizeof(ma->refCnt), 1, fp) != 1)
    return 0;
  
  if(fread(&(ma->source_alloc), sizeof(ma->source_alloc), 1, fp) != 1)
    return 0;

  ma->consensus = CreateFromFileVA_char(fp, 0);
  if(ma->consensus == NULL)
    return 0;
  
  ma->quality = CreateFromFileVA_char(fp, 0);
  if(ma->quality == NULL)
    return 0;
  
  ma->delta = CreateFromFileVA_int32(fp, 0);
  if(ma->delta == NULL)
    return 0;
  
  ma->f_list = ReadIMPsFromFile(fp);
  if(ma->f_list == NULL)
    return 0;
  
  ma->udelta = CreateFromFileVA_int32(fp, 0);
  if(ma->udelta == NULL)
    return 0;
  
  ma->u_list = ReadIUPsFromFile(fp);
  if(ma->u_list == NULL)
    return 0;

  return 1;
}

/*
  Cannibalized from MergeScaffolds(), primarily
 */
void ReplaceOldScaffoldWithNew(CIScaffoldT * oldScaffold,
                               VA_TYPE(IEPish) * cCoords)
{
  CIScaffoldT newScaffoldO;
  CIScaffoldT * newScaffoldR;
  int32 i;

  oldScaffold->setID = 0;
  oldScaffold->flags.bits.isDead = TRUE;
  
  InitializeScaffold(&newScaffoldO, REAL_SCAFFOLD);
  newScaffoldO.info.Scaffold.AEndCI = NULLINDEX;
  newScaffoldO.info.Scaffold.BEndCI = NULLINDEX;
  newScaffoldO.info.Scaffold.numElements = 0;
  newScaffoldO.edgeHead = NULLINDEX;
  newScaffoldO.bpLength.mean = newScaffoldO.bpLength.variance = 0;
  newScaffoldO.id = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
  newScaffoldO.flags.bits.isDead = FALSE;
  newScaffoldO.aEndCoord = newScaffoldO.bEndCoord = -1;
  newScaffoldO.numEssentialA = newScaffoldO.numEssentialB = 0;
  newScaffoldO.essentialEdgeB = newScaffoldO.essentialEdgeA = NULLINDEX;
  newScaffoldO.microhetScore = 0.0;
  newScaffoldO.setID = NULLINDEX;
  AppendGraphNode(ScaffoldGraph->ScaffoldGraph, &newScaffoldO);
  newScaffoldR = GetGraphNode(ScaffoldGraph->ScaffoldGraph, newScaffoldO.id); 

  fprintf(stderr,
          "\tReplacing scaffold " F_CID " with new scaffold " F_CID "\n",
          oldScaffold->id, newScaffoldO.id);
  
  // insert CIs into new scaffold
  for(i = 0; i < GetNumVA_IEPish(cCoords); i++)
  {
    IEPish * iep = GetVA_IEPish(cCoords, i);
    
    InsertCIInScaffold(ScaffoldGraph,
                       iep->ident,
                       newScaffoldO.id,
                       iep->aEnd, iep->bEnd,
                       TRUE, ALL_CONTIGGING);
  }

  fprintf(stderr, "\t\tRemarking edges of scaffold " F_CID "\n",
          newScaffoldO.id);
  
  // try to make the scaffold clone-connected
  {
#ifdef NEVER
    int status = RECOMPUTE_SINGULAR;
    int recomputeIteration = 0;
    while(recomputeIteration < 3 &&
          (status == RECOMPUTE_SINGULAR ||
           status == RECOMPUTE_CONTIGGED_CONTAINMENTS))
    {
#endif
      // need to make sure scaffold is connected with trusted raw edges
      MarkInternalEdgeStatus(ScaffoldGraph,
                             GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                          newScaffoldO.id),
                             PAIRWISECHI2THRESHOLD_CGW,
                             1000.0 * SLOPPY_EDGE_VARIANCE_THRESHHOLD,
                             TRUE, TRUE, 0, TRUE);
#ifdef NEVER
      status =
        RecomputeOffsetsInScaffold(ScaffoldGraph,
                                   GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                                newScaffoldO.id),
                                   TRUE, TRUE, FALSE);
      recomputeIteration++;
    }
    if(status != RECOMPUTE_OK)
    {
      fprintf(GlobalData->stderrc,
              "\t\tReomputeOffsetsInScaffold failed (%d) "
              "for scaffold " F_CID " in MergeScaffolds\n",
              status, newScaffoldO.id);
    }
#endif
  }
}


static int positionCompare(const IEPish * a, const IEPish * b)
{
  return
    (int)(min(a->aEnd.mean, a->bEnd.mean) - min(b->aEnd.mean, b->bEnd.mean));
}


void PopulateContigUnitigCoords(ContigT * contig,
                                VA_TYPE(IEPish) * uCoords)
{
  ContigTIterator iterator;
  ChunkInstanceT * unitig;
  int32 flip = (contig->offsetAEnd.mean > contig->offsetBEnd.mean);

  fprintf(stderr, "\tComputing unitig coordinates in scaffold " F_CID "\n",
          contig->id);
  
  InitContigTIterator(ScaffoldGraph, contig->id, flip, FALSE, &iterator);
  while((unitig = NextContigTIterator(&iterator)) != NULL)
  {
    IEPish iep;

    iep.ident = unitig->id;

    if(!flip)
    {
      iep.aEnd.mean = contig->offsetAEnd.mean + unitig->offsetAEnd.mean;
      iep.aEnd.variance =
        contig->offsetAEnd.variance + unitig->offsetAEnd.variance;
      
      iep.bEnd.mean = contig->offsetAEnd.mean + unitig->offsetBEnd.mean;
      iep.bEnd.variance =
        contig->offsetAEnd.variance + unitig->offsetBEnd.variance;
    }
    else
    {
      iep.aEnd.mean = contig->offsetBEnd.mean +
        contig->bpLength.mean - unitig->offsetAEnd.mean;
      iep.aEnd.variance = contig->offsetBEnd.variance +
        contig->bpLength.variance - unitig->offsetAEnd.variance;
      
      iep.bEnd.mean = contig->offsetBEnd.mean +
        contig->bpLength.mean - unitig->offsetBEnd.mean;
      iep.bEnd.variance = contig->offsetBEnd.variance +
        contig->bpLength.variance - unitig->offsetBEnd.variance;
    }
    AppendVA_IEPish(uCoords, &iep);
  }
  qsort(GetVA_IEPish(uCoords, 0),
        GetNumVA_IEPish(uCoords),
        sizeof(IEPish),
        (int (*) (const void *, const void *)) positionCompare);
}


void PopulateContigScaffoldCoords(CIScaffoldT * scaffold,
                                  HashTable_AS * zlfCUOs,
                                  VA_TYPE(IEPish) * cCoords)
{
  CIScaffoldTIterator CIsTemp;
  ContigT * contig;
      
  fprintf(stderr, "\tSaving contig positions in scaffold " F_CID "\n",
          scaffold->id);

  // iterate over contigs in this scaffold
  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIsTemp);
  while((contig = NextCIScaffoldTIterator(&CIsTemp)) != NULL)
  {
    VA_TYPE(IEPish) * uCoords;
    if((uCoords = LookupInHashTable_AS(zlfCUOs,
                                       (void *) &(contig->id),
                                       sizeof(int32))) != NULL)
    {
      int i;
      for(i = 0; i < GetNumVA_IEPish(uCoords); i++)
      {
        AppendVA_IEPish(cCoords, GetVA_IEPish(uCoords, i));
      }
    }
    else
    {
      IEPish iep;
      iep.ident = contig->id;
      iep.aEnd = contig->offsetAEnd;
      iep.bEnd = contig->offsetBEnd;
      AppendVA_IEPish(cCoords, &iep);
    }
  }
}


// inefficient, but lists should be short
MultiAlignT * GetZLFUMA(ZLFContig * zlfContig, CDS_CID_t unitigID)
{
  int32 i;
  for(i = 0; i < GetNumVA_MultiAlignT(zlfContig->zlfUMAs); i++)
  {
    MultiAlignT * zlfUMA  = GetVA_MultiAlignT(zlfContig->zlfUMAs, i);
    if(unitigID == zlfUMA->id)
      return zlfUMA;
  }
  return NULL;
}

void UpdateFragLinkFlags(CIFragT * frag,
                         ChunkInstanceT * node,
                         int setMatesToo)
{
  CDS_CID_t fragIndex = GetVAIndex_CIFragT(ScaffoldGraph->CIFrags, frag);
  if(frag->numLinks == 0)
  {
    assert(frag->mateOf == NULLINDEX);
    return;
  }
  
  if(frag->flags.bits.getLinksFromStore)
  {
    GateKeeperLinkRecordIterator GKPLinks;
    GateKeeperLinkRecord GKPLink;
    int prevSetting;
    int newSetting;

    if(node->flags.bits.isContig)
    {
      prevSetting = frag->flags.bits.hasInternalOnlyContigLinks;
      frag->flags.bits.hasInternalOnlyContigLinks = TRUE;
    }
    else
    {
      prevSetting = frag->flags.bits.hasInternalOnlyCILinks;
      frag->flags.bits.hasInternalOnlyCILinks = TRUE;
    }
    newSetting = TRUE;
    
    assert(frag->linkHead != NULLINDEX);
    CreateGateKeeperLinkRecordIterator(ScaffoldGraph->gkpStore.lnkStore,
                                       frag->linkHead,
                                       fragIndex, &GKPLinks);
    while(NextGateKeeperLinkRecordIterator(&GKPLinks, &GKPLink))
    {
      CIFragT * mfrag;
      CDS_CID_t mfragIndex =
        (fragIndex == GKPLink.frag1) ? GKPLink.frag2 : GKPLink.frag1;
      
      assert(fragIndex == GKPLink.frag1 || fragIndex == GKPLink.frag2 );

      mfrag = GetCIFragT(ScaffoldGraph->CIFrags, mfragIndex);
      if(node->flags.bits.isContig)
      {
        if(frag->contigID != mfrag->contigID)
        {
          frag->flags.bits.hasInternalOnlyContigLinks = FALSE;
          newSetting = FALSE;
          break;
        }
      }
      else
      {
        if(frag->cid != mfrag->cid)
        {
          frag->flags.bits.hasInternalOnlyCILinks = FALSE;
          newSetting = FALSE;
          break;
        }
      }
    }

    if(prevSetting != newSetting && setMatesToo == TRUE)
    {
      // need to examine all mates & possibly reset their flags
      CreateGateKeeperLinkRecordIterator(ScaffoldGraph->gkpStore.lnkStore,
                                         frag->linkHead,
                                         fragIndex, &GKPLinks);
      while(NextGateKeeperLinkRecordIterator(&GKPLinks, &GKPLink))
      {
        CIFragT * mfrag = NULL;
        int32 mfragIndex =
          (fragIndex == GKPLink.frag1) ? GKPLink.frag2 : GKPLink.frag1;
        ChunkInstanceT * ci;

        if(mfragIndex != NULLINDEX)
        {
          mfrag = GetCIFragT(ScaffoldGraph->CIFrags, mfragIndex);
        }
        
        if(node->flags.bits.isContig)
          ci = GetGraphNode(ScaffoldGraph->ContigGraph, mfrag->contigID);
        else
          ci = GetGraphNode(ScaffoldGraph->CIGraph, mfrag->cid);
        UpdateFragLinkFlags(mfrag, ci, FALSE);
      }
    }
  }
  else
  {
    CDS_CID_t mfragIndex = frag->mateOf;
    CIFragT * mfrag = NULL;
    int prevSetting;
    int newSetting;
    
    if(mfragIndex != NULLINDEX)
    {
      mfrag = GetCIFragT(ScaffoldGraph->CIFrags, mfragIndex);
    }
    if(node->flags.bits.isContig)
    {
      prevSetting = frag->flags.bits.hasInternalOnlyContigLinks;
      if(frag->contigID == mfrag->contigID)
      {
        frag->flags.bits.hasInternalOnlyContigLinks = TRUE;
        newSetting = TRUE;
      }
      else
      {
        frag->flags.bits.hasInternalOnlyContigLinks = FALSE;
        newSetting = FALSE;
      }
    }
    else
    {
      prevSetting = frag->flags.bits.hasInternalOnlyCILinks;
      if(frag->cid == mfrag->cid)
      {
        frag->flags.bits.hasInternalOnlyContigLinks = TRUE;
        newSetting = TRUE;
      }
      else
      {
        frag->flags.bits.hasInternalOnlyContigLinks = FALSE;
        newSetting = FALSE;
      }
    }
    if(prevSetting != newSetting && setMatesToo == TRUE)
    {
      ChunkInstanceT * ci;
      
      if(node->flags.bits.isContig)
        ci = GetGraphNode(ScaffoldGraph->ContigGraph, mfrag->contigID);
      else
        ci = GetGraphNode(ScaffoldGraph->CIGraph, mfrag->cid);

      UpdateFragLinkFlags(mfrag, ci, FALSE);
    }
  }
}


void UpdateNodeFragmentLinkFlags(GraphCGW_T * cGraph, ChunkInstanceT * node)
{
  CDS_CID_t findex;
  MultiAlignT * localMA =
    LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,
                                  node->id, cGraph->type == CI_GRAPH);
  
  /* Determine extremal fragments so we can label the fragments */
  for(findex = 0; findex < GetNumIntMultiPoss(localMA->f_list); findex++)
  {
    IntMultiPos *mp = GetIntMultiPos(localMA->f_list, findex);
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32) mp->sourceInt);

    UpdateFragLinkFlags(frag, node, TRUE);
  }
}


/*
  cannibalized from CreateAContigInScaffold, primarily
 */
void CreateNewContigsFromUnitigs(ContigT * oldContig,
                                 ZLFContig * zlfContig,
                                 VA_TYPE(IEPish) * uCoords)
{
  int i;

  fprintf(stderr, "\tCreating new contigs for contig " F_CID " unitigs\n",
          oldContig->id);
  
  zlfContig->zlfUsFound = 0;
  for(i = 0; i < GetNumVA_IEPish(uCoords); i++)
  {
    IEPish * iep = GetVA_IEPish(uCoords, i);
    ChunkInstanceT * ci = GetGraphNode(ScaffoldGraph->CIGraph, iep->ident);
    ContigT * contig = CreateNewGraphNode(ScaffoldGraph->ContigGraph);
    MultiAlignT * zlfUMA;

    fprintf(stderr, "\t\tNew singleton contig " F_CID " corresponds to unitig " F_CID "\n",
            contig->id, ci->id);
    
    if((zlfUMA = GetZLFUMA(zlfContig, iep->ident)) != NULL)
    {
      if(ci->flags.bits.isStoneSurrogate || ci->flags.bits.isWalkSurrogate)
      {
        fprintf(stderr,
                "\t\tWARNING: Unitig " F_CID " in contig " F_CID " in scaffold " F_CID " "
                "is a surrogate!\n",
                ci->id, zlfContig->id, contig->scaffoldID);
        fprintf(stderr, "\tUnsafe to replace unitig multialignment.\n");
      }
      else
      {
        int j;
        ReadStructp fsread = new_ReadStruct();

        zlfContig->zlfUsFound++;
        fprintf(stderr,
                "\t\tChanging unitig " F_CID " multialignment in "
                "contig " F_CID " in scaffold " F_CID "\n",
                ci->id, zlfContig->id, oldContig->scaffoldID);

        fprintf(stderr, "\t\tReverting fragments to CNS clear ranges.\n");
        
        // update fragment clear ranges
        for(j = 0; j < GetNumVA_IntMultiPos(zlfUMA->f_list); j++)
        {
          IntMultiPos * imp = GetVA_IntMultiPos(zlfUMA->f_list, j);
          cds_uint32 preBgn, preEnd;
          cds_uint32 postBgn, postEnd;
          
          getFragStore(ScaffoldGraph->fragStore,
                       imp->ident, FRAG_S_ALL, fsread);
          getClearRegion_ReadStruct(fsread,
                                    &preBgn, &preEnd,
                                    READSTRUCT_CNS);
          getClearRegion_ReadStruct(fsread,
                                    &postBgn, &postEnd,
                                    READSTRUCT_CGW);
          if(preBgn != postBgn || preEnd != postEnd)
          {
            fprintf(stderr,
                    "\t\t\tFragment " F_IID " clear range reset "
                    "from %u %u to %u %u\n",
                    imp->ident, postBgn, postEnd, preBgn, preEnd);
            setClearRegion_ReadStruct(fsread, preBgn, preEnd, READSTRUCT_CGW);
            setFragStore(ScaffoldGraph->fragStore, imp->ident, fsread);
          }
        }
        delete_ReadStruct(fsread);
        
        fprintf(stderr, "\t\tReplacing unitig multialignment.\n");
        // revert to the old multialignment
        {
          UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,
                                          zlfUMA->id, TRUE);
          InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB,
                                        zlfUMA->id, TRUE, zlfUMA, FALSE);
        }
      }
    }

    // update the IEP for later use
    iep->ident = contig->id;
    
    // update the unitig structure
    ci->info.CI.contigID = contig->id;
    ci->offsetAEnd.mean = ci->offsetAEnd.variance = 0;
    ci->offsetBEnd = ci->bpLength;
    ci->prevScaffoldID = ci->scaffoldID;
    ci->scaffoldID = NULLINDEX;
    ci->indexInScaffold = NULLINDEX;
    ci->AEndNext = ci->BEndNext = NULLINDEX;

    // set up the new contig
    // several things are already initialized in CreateNewGraphNode()
    contig->bpLength = ci->bpLength;
    contig->offsetAEnd = ci->offsetAEnd;
    contig->offsetBEnd = ci->offsetBEnd;
    contig->info.Contig.AEndCI = contig->info.Contig.BEndCI = ci->id;
    contig->info.Contig.numCI = 1;
    contig->flags = oldContig->flags;
    contig->flags.bits.includesFinishedBacFragments =
      ci->flags.bits.includesFinishedBacFragments;
    
    // create contig multialign from the unitig's
    DuplicateEntryInSequenceDB(ScaffoldGraph->sequenceDB,
                               ci->id, TRUE, // from is unitig
                               contig->id, FALSE, // to is contig
                               FALSE);
    
    // update fragments
    // first FALSE means not to mark fragments as placed
    // second FALSE means not to mark fragment unitig coords
    UpdateNodeFragments(ScaffoldGraph->ContigGraph, contig->id, FALSE, FALSE);

    
  }

  // delete the original contig & its edges
  fprintf(stderr, "\t\tDeleting original contig " F_CID "\n", oldContig->id);
  DeleteGraphNode(ScaffoldGraph->ContigGraph, oldContig);
  
  // create edges for new contigs
  fprintf(stderr, "\t\tCreating edges for new contigs\n");
  for(i = 0; i < GetNumVA_IEPish(uCoords); i++)
  {
    GraphEdgeStatT stats;
    IEPish * iep = GetVA_IEPish(uCoords, i);
    ContigT * contig = GetGraphNode(ScaffoldGraph->ContigGraph, iep->ident);
    MultiAlignT * newMA =
      LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,
                                    iep->ident, FALSE);
    
    UpdateNodeUnitigs(newMA, contig);
    
    BuildGraphEdgesFromMultiAlign(ScaffoldGraph->ContigGraph,
                                  contig, newMA, &stats, TRUE);
    MergeNodeGraphEdges(ScaffoldGraph->ContigGraph, contig,
                        FALSE, TRUE, FALSE);
  }

  if(zlfContig->zlfUsFound != GetNumVA_MultiAlignT(zlfContig->zlfUMAs))
  {
    fprintf(stderr,
            "WARNING: Incorrect number of unitigs changed in contig " F_CID "\n",
            zlfContig->id);
    fprintf(stderr, "\t%d unitigs identified. %d unitigs changed.\n",
            (int) GetNumVA_MultiAlignT(zlfContig->zlfUMAs),
            zlfContig->zlfUsFound);
  }
}


int ShiftContigsInScaffoldToMin0(CIScaffoldT * scaffold)
{
  CIScaffoldTIterator CIsTemp;
  ContigT * contig;
  float minOffset = 0;
      
  // iterate over contigs in this scaffold
  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIsTemp);
  while((contig = NextCIScaffoldTIterator(&CIsTemp)) != NULL)
  {
    minOffset = min(minOffset, min(contig->offsetAEnd.mean,
                                   contig->offsetBEnd.mean));
  }
  if(minOffset != 0)
  {
    InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIsTemp);
    while((contig = NextCIScaffoldTIterator(&CIsTemp)) != NULL)
    {
      contig->offsetAEnd.mean -= minOffset;
      contig->offsetBEnd.mean -= minOffset;
    }
    return 1;
  }
  return 0;
}

int FixZLFContigs(VA_TYPE(ZLFScaffold) * zlfScaffolds,
                  int checkScaffolds,
                  int splitAsNeeded)
{
  /*
    For each zlfScaffold
      For each zlfContig in the scaffold
        break up the contig
          revert to older MA for identified unitigs in this contig
   */
  int scaffoldsFixed = 0;
  int i;
  static ScaffoldInstrumenter * si = NULL;

  fprintf(stderr, "Fixing scaffolds with messed up unitigs/contigs\n");
  for(i = 0; i < GetNumVA_ZLFScaffold(zlfScaffolds); i++)
  {
    int j;
    ZLFScaffold * zlfScaffold = GetVA_ZLFScaffold(zlfScaffolds, i);
    CIScaffoldT * scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                          zlfScaffold->id);
    VA_TYPE(IEPish) * cCoords = NULL;
    VA_TYPE(IEPish) ** cUOs = NULL;
    HashTable_AS * zlfCUOs = NULL;

    zlfScaffold->zlfContigsFound = 0;
    if(scaffold == NULL ||
       scaffold->flags.bits.isDead ||
       scaffold->type != REAL_SCAFFOLD)
    {
      fprintf(stderr, "\tScaffold " F_CID " is bad/dead/not real!\n",
              zlfScaffold->id);
      continue;
    }
    
    cCoords = CreateVA_IEPish(100);
    zlfCUOs =
      CreateHashTable_int32_AS(GetNumVA_ZLFContig(zlfScaffold->zlfContigs));

    cUOs =
      (VA_TYPE(IEPish) **) safe_malloc(GetNumVA_ZLFContig(zlfScaffold->zlfContigs) *
                                  sizeof(void *));

    fprintf(stderr, "\tFixing scaffold " F_CID "\n", scaffold->id);

#ifdef NEVER
    if(checkScaffolds)
    {
      fprintf(stderr,
              "****************** INCOMING SCAFFOLD *****************\n");
      DumpACIScaffold(stderr, ScaffoldGraph, scaffold, TRUE);

      if(si == NULL)
        si = CreateScaffoldInstrumenter(ScaffoldGraph, INST_OPT_ALL_MATES);
      InstrumentScaffold(ScaffoldGraph, scaffold, si,
                         InstrumenterVerbose3, stderr);
    }
#endif

    
    // need memory for unitig orders in each zlfContig
    for(j = 0; j < GetNumVA_ZLFContig(zlfScaffold->zlfContigs); j++)
    {
      ZLFContig * zlfContig = GetVA_ZLFContig(zlfScaffold->zlfContigs, j);
      ContigT * contig = GetGraphNode(ScaffoldGraph->ContigGraph,
                                      zlfContig->id);

      if(contig == NULL ||
         contig->flags.bits.isDead)
      {
        fprintf(stderr, "WARNING: Contig " F_CID " not found in scaffold " F_CID "!\n",
                zlfContig->id, zlfScaffold->id);
        continue;
      }
      zlfScaffold->zlfContigsFound++;
      
      cUOs[j] = CreateVA_IEPish(10);
      assert(cUOs[j] != NULL);
      
      // insert the contig ID & con/unitig order array in a hashtable
      InsertInHashTable_AS(zlfCUOs,
                           (void *) &(zlfContig->id),
                           sizeof(int32),
                           (void *) cUOs[j]);
      
      PopulateContigUnitigCoords(contig, cUOs[j]);
      CreateNewContigsFromUnitigs(contig, zlfContig, cUOs[j]);
    }
    PopulateContigScaffoldCoords(scaffold, zlfCUOs, cCoords);
    ReplaceOldScaffoldWithNew(scaffold, cCoords);

    if(zlfScaffold->zlfContigsFound !=
       GetNumVA_ZLFContig(zlfScaffold->zlfContigs))
    {
      fprintf(stderr,
              "WARNING: Incorrect number of contigs found in scaffold " F_CID "\n",
              zlfScaffold->id);
      fprintf(stderr, "\t%d contigs identified. %d contigs found.\n",
              (int) GetNumVA_MultiAlignT(zlfScaffold->zlfContigs),
              zlfScaffold->zlfContigsFound);
    }
    
    // NOTE: run scaffold checks here
    if(checkScaffolds)
    {
      scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                              GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph)-1);

#ifdef NEVER
      fprintf(stderr,
              "****************** OUTGOING SCAFFOLD *****************\n");
      DumpACIScaffold(stderr, ScaffoldGraph, scaffold, TRUE);
      
      InstrumentScaffold(ScaffoldGraph, scaffold, si,
                         InstrumenterVerbose3, stderr);
#endif
      
      fprintf(stderr, "****************** CheckCIScaffoldT:\n");
      CheckCIScaffoldT(ScaffoldGraph,scaffold);
      fprintf(stderr, "****************** CheckCIScaffoldTLength:\n");
      CheckCIScaffoldTLength(ScaffoldGraph, scaffold);
      fprintf(stderr, "****************** CheckLSScaffoldWierdnesses:\n");
      CheckLSScaffoldWierdnesses("CHECK", ScaffoldGraph, scaffold);
      fprintf(stderr, "****************** CheckScaffoldOrder:\n");
      CheckScaffoldOrder(scaffold, ScaffoldGraph);
      fprintf(stderr, "****************** CheckInternalEdgeStatus:\n");
      CheckInternalEdgeStatus(ScaffoldGraph, scaffold,
                              PAIRWISECHI2THRESHOLD_CGW,
                              100000000000.0, 0, FALSE);
      if(splitAsNeeded)
      {
        fprintf(stderr, "****************** RecomputeOffsetsInScaffold:\n");
        RecomputeOffsetsInScaffold(ScaffoldGraph, scaffold,
                                   TRUE, FALSE, FALSE);

        // iterate over contigs to make sure min offset = 0
        if(ShiftContigsInScaffoldToMin0(scaffold))
        {
          fprintf(stderr,
                  "RecomputeOffsetsInScaffold required contig shifting "
                  "to maintain 0 offset\n");
        }
        CheckCIScaffoldT(ScaffoldGraph,scaffold);
        RecomputeOffsetsInScaffold(ScaffoldGraph, scaffold,
                                   FALSE, FALSE, FALSE);
        InstrumentScaffold(ScaffoldGraph, scaffold, si,
                           InstrumenterVerbose3, stderr);
      

        fprintf(stderr,
                "****************** CheckScaffoldConnectivityAndSplit:\n");
        CheckScaffoldConnectivityAndSplit(ScaffoldGraph, scaffold,
                                          ALL_TRUSTED_EDGES, FALSE);
      }
    }
    
    DeleteVA_IEPish(cCoords);
    DeleteHashTable_AS(zlfCUOs);
    for(j = 0; j < GetNumVA_ZLFContig(zlfScaffold->zlfContigs); j++)
    {
      DeleteVA_IEPish(cUOs[j]);
    }
    free(cUOs);
    
    scaffoldsFixed++;
  }
  if(scaffoldsFixed != GetNumVA_ZLFScaffold(zlfScaffolds))
  {
    fprintf(stderr, "WARNING: Incorrect number of scaffolds fixed!\n");
    fprintf(stderr, "\t%d scaffolds identified. %d scaffolds fixed.\n",
            (int) GetNumVA_MultiAlignT(zlfScaffolds), scaffoldsFixed);
  }
  return scaffoldsFixed;
}

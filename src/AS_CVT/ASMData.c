
#warning ASMData.c not working.

#if 0
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
#include <memory.h>
#include <assert.h>

#include "AS_PER_gkpStore.h"
#include "AS_PER_asmStore.h"
#include "ASMData.h"

//#define DEBUG_ASMDATA

void PopulateInstance(ASM_InstanceRecord * ins,
                      uint32 index,
                      SeqInterval pos)
{
  ins->containerIndex = index;
  ins->pos = pos;
  ins->next = 0;
}


void PopulateDerivedInstance(ASM_InstanceRecord * childInGrand,
                             uint32 grandIndex,
                             SeqInterval parentInGrand,
                             SeqInterval childInParent)
{
  childInGrand->containerIndex = grandIndex;
  childInGrand->next = 0;
  
  if(parentInGrand.bgn < parentInGrand.end)
  {
    childInGrand->pos.bgn = parentInGrand.bgn + childInParent.bgn;
    childInGrand->pos.end = parentInGrand.bgn + childInParent.end;
  }
  else
  {
    childInGrand->pos.bgn = parentInGrand.bgn - childInParent.bgn;
    childInGrand->pos.end = parentInGrand.bgn - childInParent.end;
  }
}


void AppendInstance(ASM_InstanceStore insStore,
                    ASM_InstanceRecord ins,
                    int32 * prevIndex)
{
  ASM_InstanceRecord lastIns;
  int32 lastIndex;

  int32 newIndex = getNumASM_Instances(insStore) + 1;
  if(*prevIndex == 0)
  {
    *prevIndex = newIndex;
  }
  else
  {
    lastIndex = *prevIndex;
    getASM_InstanceStore(insStore, lastIndex, &lastIns);
    while(lastIns.next != 0)
    {
      lastIndex = lastIns.next;
      getASM_InstanceStore(insStore, lastIndex, &lastIns);
    }
    lastIns.next = newIndex;
    setASM_InstanceStore(insStore, lastIndex, &lastIns);
  }
  ins.next = 0;
  appendASM_InstanceStore(insStore, &ins);
}


void FixScaffoldInstances(AssemblyStore * asmStore)
{
  int32 i, j;
  int32 numUTGs = getNumASM_UTGs(asmStore->utgStore);
  ASM_UTGRecord utg;
  ASM_AFGRecord afg;
  ASM_CCORecord cco;
  ASM_InstanceRecord cIns;
  ASM_InstanceRecord sIns;
  ASM_IIDRecord fid;
  int32 numV= 0;

  for(i = 1; i <= numUTGs; i++)
  {
    getASM_UTGStore(asmStore->utgStore, i, &utg);
    if(utg.status != AS_UNIQUE)
      continue;

    for(j = utg.firstFrag; j < utg.firstFrag + utg.numFrags; j++, numV++)
    {
      if(numV != 0 && numV % 10000 == 0)
        fprintf(stderr, "\r%20d", numV);
      getASM_IIDStore(asmStore->utfStore, j, &fid);
      getASM_AFGStore(asmStore->afgStore, fid, &afg);
      if(afg.deleted || afg.unreferenced)
        continue;
      getASM_InstanceStore(asmStore->aciStore, afg.cInsIndex, &cIns);
      getASM_CCOStore(asmStore->ccoStore, cIns.containerIndex, &cco);
      getASM_InstanceStore(asmStore->asiStore, afg.sInsIndex, &sIns);
      assert(cIns.next == 0 && sIns.next == 0);
      PopulateDerivedInstance(&sIns, sIns.containerIndex,
                              cco.scaffoldPos, cIns.pos);
      setASM_InstanceStore(asmStore->asiStore, afg.sInsIndex, &sIns);
    }
  }
}


void AppendNewUnitigFragmentInstances(AssemblyStore * asmStore,
                                      int32 inDegenerate,
                                      int32 unitigStatus,
                                      int32 firstFrag, int32 numFrags,
                                      uint32 contigIndex,
                                      SeqInterval unitigInContig,
                                      uint32 scaffoldIndex,
                                      SeqInterval unitigInScaffold)
{
  int32 i;
  ASM_IIDRecord afgIndex;
  ASM_AFGRecord afg;
  ASM_CCORecord cco;
  ASM_InstanceRecord ins;
  ASM_InstanceRecord contigIns;
  
  for(i = firstFrag; i < firstFrag + numFrags; i++)
  {
    getASM_IIDStore(asmStore->utfStore, i, &afgIndex);
    getASM_AFGStore(asmStore->afgStore, afgIndex, &afg);

    if(afg.deleted) continue;  // shouldn't happen!
    
    afg.inDegenerate = inDegenerate;

    if(unitigStatus == AS_SEP)
    {
      // create contig instance
      PopulateDerivedInstance(&ins, contigIndex,
                              unitigInContig, afg.unitigPos);
      AppendInstance(asmStore->aciStore, ins, &(afg.cInsIndex));
      afg.unreferenced = FALSE;

      // create scaffold instance from UNITIG instance
      PopulateDerivedInstance(&ins, scaffoldIndex,
                              unitigInScaffold, afg.unitigPos);
      AppendInstance(asmStore->asiStore, ins, &(afg.sInsIndex));
    }
    else
    {
      if(afg.cInsIndex == 0)
      {
        // if at this point, the fragment has no contig instances,
        // it's unreferenced
        afg.unreferenced = TRUE;
        PopulateDerivedInstance(&ins, contigIndex,
                                unitigInContig, afg.unitigPos);
        AppendInstance(asmStore->aciStore, ins, &(afg.cInsIndex));
      }

      // create scaffold instance from CONTIG instance
      getASM_InstanceStore(asmStore->aciStore, afg.cInsIndex, &contigIns);
      getASM_CCOStore(asmStore->ccoStore, contigIns.containerIndex, &cco);
      PopulateDerivedInstance(&ins, cco.scaffoldIndex,
                              cco.scaffoldPos, contigIns.pos);
      AppendInstance(asmStore->asiStore, ins, &(afg.sInsIndex));
    }
    setASM_AFGStore(asmStore->afgStore, afgIndex, &afg);
  }
}


void PrintFragmentScaffoldCoordinates(AssemblyStore * asmStore,
                                      int doInstances,
                                      int doSingleSurrogates,
                                      int doDegenerates,
                                      int doChaff,
                                      int doLinks,
                                      int doUnreferenced,
                                      FILE * fo)
{
  uint32 i;
  ASM_AFGRecord afg;
  // ASM_UTGRecord utg;
  // ASM_IIDRecord itemIndex;
  // ASM_CCORecord cco;
  ASM_DSCRecord dsc;
  ASM_SCFRecord scf;
  ASM_InstanceRecord contigIns;
  ASM_InstanceRecord scaffIns;
  ASM_LKGRecord lkg;
  // int32 numUTGs = getNumASM_UTGs(asmStore->utgStore);
  // int32 numCCOs = getNumASM_CCOs(asmStore->ccoStore);
  // int32 numDSCs = getNumASM_DSCs(asmStore->dscStore);
  // int32 numSCFs = getNumASM_SCFs(asmStore->scfStore);
  int32 numAFGs = getNumASM_AFGs(asmStore->afgStore);
  int32 printed;
  
  // format is fragment, contig, scaffold, 5p in scaffold, 3p in scaffold
  // and if mated: mate, lib
  for(i = 1; i <= numAFGs; i++)
  {
    getASM_AFGStore(asmStore->afgStore, i, &afg);

    if(afg.deleted) continue;
    
    printed = FALSE;
    if(afg.chaff)
    {
      if(doChaff)
      {
        fprintf(fo, F_UID " %c C 0 (0,0) ", afg.uid, (char) afg.type);
        printed = TRUE;
      }
      else
        continue;
    }
    else
    {
      // getASM_UTGStore(asmStore->utgStore, afg.unitigIndex, &utg);
      if(afg.inDegenerate)
      {
        if(doDegenerates)
        {
          getASM_InstanceStore(asmStore->asiStore, afg.sInsIndex, &scaffIns);
          getASM_DSCStore(asmStore->dscStore, scaffIns.containerIndex, &dsc);
          fprintf(fo, F_UID " %c D " F_UID " (" F_COORD "," F_COORD ") ",
                  afg.uid, (char) afg.type,
                  dsc.uid, scaffIns.pos.bgn, scaffIns.pos.end);
          printed = TRUE;
        }
        else
          continue;
      }
      else
      {
        // in real scaffold(s) - or missing from contigs/scaffolds
        if(afg.cInsIndex != 0 && afg.sInsIndex != 0 &&
           (!afg.inSurrogate || doSingleSurrogates) &&
           (!afg.unreferenced || doUnreferenced))
        {
          getASM_InstanceStore(asmStore->asiStore, afg.sInsIndex, &scaffIns);
          getASM_SCFStore(asmStore->scfStore, scaffIns.containerIndex, &scf);
          if(scaffIns.next == 0 || doInstances)
          {
            fprintf(fo, F_UID " %c %c " F_UID " (" F_COORD "," F_COORD ") ",
                    afg.uid, (char) afg.type,
                    (afg.inSurrogate || scaffIns.next != 0) ? 'S' : 
                     (afg.unreferenced ? 'U' : 'R'),
                    scf.uid, scaffIns.pos.bgn, scaffIns.pos.end);
            printed = TRUE;
          }
          
          while(doInstances && scaffIns.next && contigIns.next)
          {
            getASM_InstanceStore(asmStore->asiStore, scaffIns.next, &scaffIns);
            getASM_SCFStore(asmStore->scfStore, scaffIns.containerIndex, &scf);
            fprintf(fo, "\n" F_UID " %c S " F_UID " (" F_COORD "," F_COORD ") ",
                    afg.uid, (char) afg.type,
                    scf.uid, scaffIns.pos.bgn, scaffIns.pos.end);
          }
        }
      }
    }

    // print mate & library
    if(afg.numLinks > 0 && doLinks && printed)
    {
      CreateGateKeeperLinkRecordIterator(asmStore->lkgStore, afg.linkHead,
                                         i, &iterator);
      while(NextGateKeeperLinkRecordIterator(&iterator, &lkg))
      {
        if(lkg.type == AS_MATE)
        {
          ASM_MDIRecord mdi;
          ASM_AFGRecord afg2;
          CDS_IID_t iid = (lkg.frag1 == i) ? lkg.frag2 : lkg.frag1;

          getASM_AFGStore(asmStore->afgStore, iid, &afg2);
          getASM_MDIStore(asmStore->mdiStore, lkg.distance, &mdi);
          fprintf(fo, "\t" F_UID "\t" F_UID, afg2.uid, mdi.uid);
          break;
        }
      }
    }
    if(printed)
      fprintf(fo, "\n");
  }
}


/*
  Format of reads.placed file:
    white-space delimited, one line per fragment
     # (a) NCBI ti number for read (or *, if none known)
     # (b) read name
     # (c) start of trimmed read on original read
     # (d) number of bases in trimmed read
     # (e) orientation on contig (0 = forward, 1 = reverse)
     # (f) contig name
     # (g) supercontig name
     # (h) approximate start of trimmed read on contig
     # (i) approximate start of trimmed read on supercontig.

     for (a), we write frag UID
         (b), concatenate frag type, parent type, and frag UID
 */
void PrintReadsPlaced(AssemblyStore * asmStore,
                      int doInstances,
                      int doSingleSurrogates,
                      int doDegenerates,
                      int doChaff,
                      int doUnreferenced,
                      FILE * fo)
{
  uint32 i;
  ASM_AFGRecord afg;
  ASM_DSCRecord dsc;
  ASM_CCORecord cco;
  ASM_SCFRecord scf;
  ASM_InstanceRecord contigIns;
  ASM_InstanceRecord scaffIns;
  ASM_LKGRecord lkg;
  GateKeeperLinkRecordIterator iterator;
  int32 numAFGs = getNumASM_AFGs(asmStore->afgStore);
  
  for(i = 1; i <= numAFGs; i++)
  {
    getASM_AFGStore(asmStore->afgStore, i, &afg);

    if(afg.deleted) continue;
    
    if(afg.chaff)
    {
      if(doChaff)
      {
        fprintf(fo, F_UID " %cC" F_UID " " F_COORD " " F_COORD " 0 - - - -\n",
                afg.uid, (char) afg.type, afg.uid,
                afg.asmClr.bgn + 1, afg.asmClr.end - afg.asmClr.bgn);
      }
      else
        continue;
    }
    else
    {
      if(afg.inDegenerate)
      {
        if(doDegenerates)
        {
          getASM_InstanceStore(asmStore->asiStore, afg.sInsIndex, &scaffIns);
          getASM_DSCStore(asmStore->dscStore, scaffIns.containerIndex, &dsc);
          fprintf(fo, F_UID " %cD" F_UID " " F_COORD " " F_COORD
                  " %d " F_UID " " F_UID " " F_COORD " " F_COORD "\n",
                  afg.uid, (char) afg.type, afg.uid,
                  afg.asmClr.bgn + 1, afg.asmClr.end - afg.asmClr.bgn,
                  ((scaffIns.pos.bgn < scaffIns.pos.end) ? 0 : 1),
                  dsc.uid, dsc.uid, scaffIns.pos.bgn, scaffIns.pos.bgn);
        }
        else
          continue;
      }
      else
      {
        // in real scaffold(s) - or missing from contigs/scaffolds
        if(afg.cInsIndex != 0 && afg.sInsIndex != 0 &&
           (!afg.inSurrogate || doSingleSurrogates) &&
           (!afg.unreferenced || doUnreferenced))
        {
          getASM_InstanceStore(asmStore->aciStore, afg.cInsIndex, &contigIns);
          getASM_CCOStore(asmStore->ccoStore, contigIns.containerIndex, &cco);
          getASM_InstanceStore(asmStore->asiStore, afg.sInsIndex, &scaffIns);
          getASM_SCFStore(asmStore->scfStore, scaffIns.containerIndex, &scf);
          if(scaffIns.next == 0 || doInstances)
          {
            fprintf(fo, F_UID " %c%c" F_UID " " F_COORD " " F_COORD
                    " %d " F_UID " " F_UID " " F_COORD " " F_COORD "\n",
                    afg.uid, (char) afg.type,
                    (afg.inSurrogate || scaffIns.next != 0) ? 'S' : 
                    (afg.unreferenced ? 'U' : 'R'), afg.uid,
                    afg.asmClr.bgn + 1, afg.asmClr.end - afg.asmClr.bgn,
                    ((contigIns.pos.bgn < contigIns.pos.end) ? 0 : 1),
                    cco.uid, scf.uid,
                    contigIns.pos.bgn, scaffIns.pos.bgn);
          }
          
          while(doInstances && scaffIns.next && contigIns.next)
          {
            getASM_InstanceStore(asmStore->aciStore, afg.cInsIndex, &contigIns);
            getASM_CCOStore(asmStore->ccoStore, contigIns.containerIndex, &cco);
            getASM_InstanceStore(asmStore->asiStore, scaffIns.next, &scaffIns);
            getASM_SCFStore(asmStore->scfStore, scaffIns.containerIndex, &scf);
            fprintf(fo, F_UID " %c%c" F_UID " " F_COORD " " F_COORD
                    " %d " F_UID " " F_UID " " F_COORD " " F_COORD "\n",
                    afg.uid, (char) afg.type,
                    (afg.inSurrogate || scaffIns.next != 0) ? 'S' : 
                    (afg.unreferenced ? 'U' : 'R'), afg.uid,
                    afg.asmClr.bgn + 1, afg.asmClr.end - afg.asmClr.bgn,
                    ((contigIns.pos.bgn < contigIns.pos.end) ? 0 : 1),
                    cco.uid, scf.uid,
                    contigIns.pos.bgn, scaffIns.pos.bgn);
          }
        }
      }
    }
  }
}


int AddMDI2Store(AssemblyStore * asmStore, SnapMateDistMesg * smdm)
{
  ASM_MDIRecord mdi;
  PHashValue_AS value;
  int i;

  // append buckets to bucketstore
  mdi.firstBucket = getNumASM_Buckets(asmStore->bktStore) + 1;
  mdi.numBuckets = smdm->num_buckets;
  for(i = 0; i < mdi.numBuckets; i++)
    appendASM_BucketStore(asmStore->bktStore, &smdm->histogram[i]);
  
  mdi.uid = smdm->erefines;
  mdi.asmMean = smdm->mean;
  mdi.asmStddev = smdm->stddev;
  mdi.inMean = mdi.inStddev = 0.0;

  if(asmStore->gkpStore != NULL)
  {
    if(HASH_FAILURE == LookupInPHashTable_AS(asmStore->gkpStore->hashTable,
                                             ASM_UID_NAMESPACE,
                                             mdi.uid,
                                             &value))
    {
      // look up original in gatekeeper store
      fprintf(stderr, "Failed to lookup " F_UID " in gatekeeper store hashtable\n",
              mdi.uid);
    }
    else
    {
      GateKeeperDistanceRecord gkdr;
      getGateKeeperDistanceStore(asmStore->gkpStore->dstStore,
                                 value.IID, &gkdr);
      mdi.inMean = gkdr.mean;
      mdi.inStddev = gkdr.stddev;
    }
  }
  appendASM_MDIStore(asmStore->mdiStore, &mdi);

  memset(&value, 0, sizeof(PHashValue_AS));
  value.type = AS_IID_MDI;
  InsertInPHashTable_AS(&(asmStore->hashTable), ASM_UID_NAMESPACE,
                        mdi.uid, &value, FALSE, TRUE);
  return 0;
}


int AddAFG2Store(AssemblyStore * asmStore, AugFragMesg * afg)
{
  ASM_AFGRecord myAFG;
  PHashValue_AS value;
  static ReadStructp rs;
  GateKeeperFragmentRecord gkfr;
  
  if(asmStore->gkpStore == NULL)
  { 
    memset(&myAFG, 0, sizeof(ASM_AFGRecord));
  
    myAFG.uid = afg->eaccession;
    myAFG.deleted = FALSE;
    appendASM_AFGStore(asmStore->afgStore, &myAFG);
    
    memset(&value, 0, sizeof(PHashValue_AS));
    value.type = AS_IID_AFG;
    InsertInPHashTable_AS(&(asmStore->hashTable), ASM_UID_NAMESPACE,
                          myAFG.uid, &value, FALSE, TRUE);
  }

  if(HASH_FAILURE == LookupInPHashTable_AS(asmStore->hashTable,
                                           ASM_UID_NAMESPACE,
                                           afg->eaccession,
                                           &value))
  {
    fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n",
            afg->eaccession);
    assert(0);
  }
  getASM_AFGStore(asmStore->afgStore, value.IID, &myAFG);
  
  myAFG.uid = afg->eaccession;
  myAFG.chaff = afg->chaff;
  myAFG.chimeric = afg->chimeric;
  myAFG.asmClr = afg->clear_rng;
  
  // get link fields - 1:1 correspondence between AFGs & gkfrs
  if(asmStore->gkpStore != NULL)
  {
    getGateKeeperFragmentStore(asmStore->gkpStore->frgStore,
                               value.IID, &gkfr);
    myAFG.type = gkfr.type;
    myAFG.numLinks = gkfr.numLinks;
    myAFG.linkHead = gkfr.linkHead;
  }
  
  // get input clear range
  if(asmStore->frgStore != NULLSTOREHANDLE)
  {
    rs = (rs == NULL) ? new_ReadStruct() : rs;
    getFragStore(asmStore->frgStore, afg->iaccession, FRAG_S_FIXED, rs);
    getClearRegion_ReadStruct(rs,
                              (uint32 *) &(myAFG.inClr.bgn),
                              (uint32 *) &(myAFG.inClr.end),
                              READSTRUCT_ORIGINAL);
  }
  setASM_AFGStore(asmStore->afgStore, value.IID, &myAFG);
  return 0;
}


int AddUTG2Store(AssemblyStore * asmStore, SnapUnitigMesg * sum)
{
  ASM_UTGRecord utg;
  PHashValue_AS value;
  int i;
  IntChunk_ID unitigIndex = getNumASM_UTGs(asmStore->utgStore) + 1;

  memset(&utg, 0, sizeof(ASM_UTGRecord));
  
  utg.uid = sum->eaccession;
  utg.coverageStat = sum->coverage_stat;
  utg.status = sum->status;
  utg.length = sum->length;
  utg.numFrags = sum->num_frags;
  utg.firstFrag = getNumASM_IIDs(asmStore->utfStore) + 1;

  for(i = 0; i < utg.numFrags; i++)
  {
    ASM_AFGRecord afg;

    if(HASH_FAILURE == LookupInPHashTable_AS(asmStore->hashTable,
                                             ASM_UID_NAMESPACE,
                                             sum->f_list[i].eident,
                                             &value))
    {
      fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n",
              sum->f_list[i].eident);
    }
    appendASM_IIDStore(asmStore->utfStore, &(value.IID));
    
    getASM_AFGStore(asmStore->afgStore, value.IID, &afg);
    afg.unitigIndex = unitigIndex;
    afg.type = sum->f_list[i].type;
    afg.inSurrogate = (utg.status == AS_SEP);
    afg.unitigPos = sum->f_list[i].position;
    setASM_AFGStore(asmStore->afgStore, value.IID, &afg);
  }
  appendASM_UTGStore(asmStore->utgStore, &utg);
  
  memset(&value, 0, sizeof(PHashValue_AS));
  value.type = AS_IID_UTG;
  InsertInPHashTable_AS(&(asmStore->hashTable), ASM_UID_NAMESPACE,
                        utg.uid, &value, FALSE, TRUE);
  return 0;
}


int AddCCO2Store(AssemblyStore * asmStore, SnapConConMesg * sccm)
{
  ASM_CCORecord cco;
  ASM_UTGRecord utg;
  ASM_AFGRecord afg;
  ASM_InstanceRecord ins;
  PHashValue_AS value;
  int i;
  IntContig_ID contigIndex = getNumASM_CCOs(asmStore->ccoStore) + 1;
  
#ifdef DEBUG_ASMDATA
  fprintf(stderr, "Adding CCO " F_UID "\n", sccm->eaccession);
#endif

  memset(&cco, 0, sizeof(ASM_CCORecord));
  
  cco.uid = sccm->eaccession;
  cco.length = sccm->length;
  cco.numFrags = sccm->num_pieces;
  cco.firstFrag = getNumASM_IIDs(asmStore->ccfStore) + 1;
  
  cco.numUnitigs = sccm->num_unitigs;
  cco.firstUnitig = getNumASM_IIDs(asmStore->ccuStore) + 1;
  
  for(i = 0; i < cco.numFrags; i++)
  {
    // look up the afg & set contig fields
    if(HASH_FAILURE == LookupInPHashTable_AS(asmStore->hashTable,
                                             ASM_UID_NAMESPACE,
                                             sccm->pieces[i].eident,
                                             &value))
    {
      fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n",
              sccm->pieces[i].eident);
    }
    appendASM_IIDStore(asmStore->ccfStore, &(value.IID));

    // each of these fragments should appear in only one contig
    getASM_AFGStore(asmStore->afgStore, value.IID, &afg);
    afg.unreferenced = FALSE;
    PopulateInstance(&ins, contigIndex, sccm->pieces[i].position);
    AppendInstance(asmStore->aciStore, ins, &(afg.cInsIndex));
    setASM_AFGStore(asmStore->afgStore, value.IID, &afg);
  }
  
  for(i = 0; i < cco.numUnitigs; i++)
  {
    // look up the afg & set contig fields
    if(HASH_FAILURE == LookupInPHashTable_AS(asmStore->hashTable,
                                             ASM_UID_NAMESPACE,
                                             sccm->unitigs[i].eident,
                                             &value))
    {
      fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n",
              sccm->unitigs[i].eident);
    }
    appendASM_IIDStore(asmStore->ccuStore, &(value.IID));

    // some of these unitigs may appear in multiple contigs
    getASM_UTGStore(asmStore->utgStore, value.IID, &utg);
    PopulateInstance(&ins, contigIndex, sccm->unitigs[i].position);
    AppendInstance(asmStore->uciStore, ins, &(utg.cInsIndex));
    utg.numInstances++;
    setASM_UTGStore(asmStore->utgStore, value.IID, &utg);
  }
  appendASM_CCOStore(asmStore->ccoStore, &cco);

  memset(&value, 0, sizeof(PHashValue_AS));
  value.type = AS_IID_CCO;
  InsertInPHashTable_AS(&(asmStore->hashTable), ASM_UID_NAMESPACE,
                        cco.uid, &value, FALSE, TRUE);
  return 0;
}


void PrintInstanceInterval(CDS_UID_t uid1, CDS_UID_t uid2,
                           CDS_COORD_t fiveP, CDS_COORD_t threeP,
                           int32 numInstances, FILE * fo)
{
  fprintf(fo, F_UID "\t" F_UID "\t" F_COORD " " F_COORD "\t%d\n",
          uid1, uid2, fiveP, threeP, numInstances);
}


void PrintSurrogateCoordinates(AssemblyStore * asmStore, FILE * fo)
{
  int i;
  int j;
  ASM_UTGRecord utg;
  ASM_SCFRecord scf;
  ASM_InstanceRecord scaffIns;
  int32 sInsIndex;
  int numUTGs = getNumASM_UTGs(asmStore->utgStore);
  
  for(j = 1, i = 1; i <= numUTGs; i++)
  {
    getASM_UTGStore(asmStore->utgStore, i, &utg);
    if(utg.status == AS_SEP)
    {
      for(sInsIndex = utg.sInsIndex;
          sInsIndex != 0;
          sInsIndex = scaffIns.next, j++)
      {
        getASM_InstanceStore(asmStore->usiStore, sInsIndex, &scaffIns);
        getASM_SCFStore(asmStore->scfStore, scaffIns.containerIndex, &scf);
        PrintInstanceInterval(utg.uid, scf.uid,
                              scaffIns.pos.bgn, scaffIns.pos.end,
                              utg.numInstances, fo);
      }
    }
  }
}


void PrintSurrogateSeqInterval(ASM_UTGStore utgStore,
                               int32 lastValue,
                               CDS_COORD_t firstCoordLastValue,
                               CDS_COORD_t k,
                               SeqInterval contigPos,
                               CDS_UID_t scfUID,
                               FILE * fo)
{
  ASM_UTGRecord utg;
  
  // print surrogate interval;
  if(contigPos.bgn < contigPos.end)
  {
    if(lastValue < 0)
    {
      getASM_UTGStore(utgStore, -lastValue, &utg);
      PrintInstanceInterval(utg.uid, scfUID,
                            contigPos.bgn + k,
                            contigPos.bgn + firstCoordLastValue,
                            utg.numInstances, fo);
    }
    else
    {
      getASM_UTGStore(utgStore, lastValue, &utg);
      PrintInstanceInterval(utg.uid, scfUID,
                            contigPos.bgn + firstCoordLastValue,
                            contigPos.bgn + k,
                            utg.numInstances, fo);
    }
  }
  else
  {
    if(lastValue < 0)
    {
      getASM_UTGStore(utgStore, -lastValue, &utg);
      PrintInstanceInterval(utg.uid, scfUID,
                            contigPos.bgn - k,
                            contigPos.bgn - firstCoordLastValue,
                            utg.numInstances, fo);
    }
    else
    {
      getASM_UTGStore(utgStore, lastValue, &utg);
      PrintInstanceInterval(utg.uid, scfUID,
                            contigPos.bgn - firstCoordLastValue,
                            contigPos.bgn - k,
                            utg.numInstances, fo);
    }
  }
}


void PrintSurrogateSequenceCoordinates(AssemblyStore * asmStore, FILE * fo)
{
  int i, j, k;
  ASM_CCORecord cco;
  ASM_IIDRecord iid;
  ASM_UTGRecord utg;
  ASM_SCFRecord scf;
  ASM_AFGRecord afg;
  ASM_InstanceRecord contigIns;
  int numCCOs = getNumASM_CCOs(asmStore->ccoStore);
  CDS_COORD_t maxLength = 0;
  CDS_IID_t * nonU;
  
  for(i = 1; i <= numCCOs; i++)
  {
    getASM_CCOStore(asmStore->ccoStore, i, &cco);
    maxLength = (maxLength < cco.length) ? cco.length : maxLength;
  }

  maxLength += 10000; // for safety...
  nonU = (CDS_IID_t *) calloc(maxLength, sizeof(CDS_IID_t));
  assert(nonU != NULL);
  
  for(i = 1; i <= numCCOs; i++)
  {
    getASM_CCOStore(asmStore->ccoStore, i, &cco);
    getASM_IIDStore(asmStore->ccfStore, cco.firstFrag, &iid);
    getASM_AFGStore(asmStore->afgStore, iid, &afg);

    // skip contigs in degenerate scaffolds
    if(afg.inDegenerate)
      continue;

    // populate surrogate intervals
    for(j = cco.firstUnitig; j < cco.firstUnitig + cco.numUnitigs; j++)
    {
      getASM_IIDStore(asmStore->ccuStore, j, &iid);
      getASM_UTGStore(asmStore->utgStore, iid, &utg);
      if(utg.status == AS_SEP || utg.numInstances > 1)
      {
        int32 ci;
        // unitig may be in contig multiple times...
        for(ci = utg.cInsIndex; ci != 0; ci = contigIns.next)
        {
          getASM_InstanceStore(asmStore->uciStore, ci, &contigIns);
          if(contigIns.containerIndex == i)
          {
            if(contigIns.pos.bgn < contigIns.pos.end)
            {
              // unitig is forward-oriented in contig
              for(k = contigIns.pos.bgn; k < contigIns.pos.end; k++)
                nonU[k] = iid;
            }
            else
            {
              // unitig is reverse-oriented in contig
              for(k = contigIns.pos.end; k < contigIns.pos.bgn; k++)
                nonU[k] = -iid;
            }
          }
        }
      }
    }

    // overwrite with non-surrogate intervals
    for(j = cco.firstUnitig; j < cco.firstUnitig + cco.numUnitigs; j++)
    {
      getASM_IIDStore(asmStore->ccuStore, j, &iid);
      getASM_UTGStore(asmStore->utgStore, iid, &utg);
      if(utg.status != AS_SEP && utg.numInstances < 2)
      {
        CDS_COORD_t start;
        CDS_COORD_t end;
        // unitig may be in contig multiple times...
        getASM_InstanceStore(asmStore->uciStore, utg.cInsIndex, &contigIns);
        assert(contigIns.containerIndex == i);
        
        start = MIN(contigIns.pos.bgn, contigIns.pos.end);
        end = MAX(contigIns.pos.bgn, contigIns.pos.end);
        for(k = start; k < end; k++)
              nonU[k] = 0;
      }
    }

    // iterate over nonU & print surrogate unitigs
    {
      int lastValue = 0;
      CDS_COORD_t firstCoordLastValue = 0;
      getASM_SCFStore(asmStore->scfStore, cco.scaffoldIndex, &scf);
      for(k = 0; k < cco.length; k++)
      {
        if(nonU[k] != lastValue)
        {
          if(lastValue != 0)
          {
            PrintSurrogateSeqInterval(asmStore->utgStore, lastValue,
                                      firstCoordLastValue, k,
                                      cco.scaffoldPos, scf.uid, fo);
          }
          lastValue = nonU[k];
          firstCoordLastValue = k;
        }
        nonU[k] = 0;
      }
      if(lastValue != 0)
        PrintSurrogateSeqInterval(asmStore->utgStore, lastValue,
                                  firstCoordLastValue, k,
                                  cco.scaffoldPos, scf.uid, fo);

    }
  }
  if(nonU != NULL)
    free(nonU);
}


int SetContigScaffoldFields(AssemblyStore * asmStore,
                            IntScaffold_ID scaffoldIndex,
                            IntContig_ID contigIndex,
                            int isDegenerate,
                            CDS_COORD_t * offset,
                            DirectionType direction)
{
  ASM_CCORecord cco;
  getASM_CCOStore(asmStore->ccoStore, contigIndex, &cco);

  cco.scaffoldIndex = scaffoldIndex;
  cco.inDegenerate = isDegenerate;
  if(direction == AS_FORWARD)
  {
    cco.scaffoldPos.bgn = *offset;
    cco.scaffoldPos.end = cco.length + *offset;
  }
  else
  {
    cco.scaffoldPos.bgn = cco.length + *offset;
    cco.scaffoldPos.end = *offset;
  }
  *offset += cco.length;
  setASM_CCOStore(asmStore->ccoStore, contigIndex, &cco);

  return 0;
}


int AddDSC2Store(AssemblyStore * asmStore, SnapDegenerateScaffoldMesg * sdsm)
{
  ASM_DSCRecord dsc;
  PHashValue_AS value;
  CDS_COORD_t offset = 0;
  
#ifdef DEBUG_ASMDATA
  fprintf(stderr, "Adding DSC " F_UID "\n", sdsm->eaccession);
#endif
  
  // look up the cco
  if(HASH_FAILURE == LookupInPHashTable_AS(asmStore->hashTable,
                                           ASM_UID_NAMESPACE,
                                           sdsm->econtig,
                                           &value))
  {
    fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n",
            sdsm->econtig);
  }
  dsc.uid = sdsm->eaccession;
  dsc.contigIndex = value.IID;

  // set scaffold fields in contig, unitigs, & afgs
  SetContigScaffoldFields(asmStore,
                          getNumASM_DSCs(asmStore->dscStore) + 1,
                          value.IID,
                          TRUE,
                          &offset, AS_FORWARD);
  
  appendASM_DSCStore(asmStore->dscStore, &dsc);
  
  memset(&value, 0, sizeof(PHashValue_AS));
  value.type = AS_IID_DSC;
  InsertInPHashTable_AS(&(asmStore->hashTable), ASM_UID_NAMESPACE,
                        dsc.uid, &value, FALSE, TRUE);
  return 0;
}


int AddSCF2Store(AssemblyStore * asmStore, SnapScaffoldMesg * ssm)
{
  ASM_SCFRecord scf;
  ASM_GapRecord gap;
  int i;
  CDS_COORD_t offset = 0;
  PHashValue_AS value;
  
#ifdef DEBUG_ASMDATA
  fprintf(stderr, "Adding SCF " F_UID "\n", ssm->eaccession);
#endif
  
  scf.uid = ssm->eaccession;
  scf.numContigs = ssm->num_contig_pairs + 1;
  scf.firstContig = getNumASM_IIDs(asmStore->sccStore) + 1;
  scf.firstGap = getNumASM_Gaps(asmStore->scgStore) + 1;

  for(i = 0; i < ssm->num_contig_pairs; i++)
  {
    // look up the cco
    if(HASH_FAILURE == LookupInPHashTable_AS(asmStore->hashTable,
                                             ASM_UID_NAMESPACE,
                                             ssm->contig_pairs[i].econtig1,
                                             &value))
    {
      fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n",
              ssm->contig_pairs[i].econtig1);
    }
    appendASM_IIDStore(asmStore->sccStore, &(value.IID));

    SetContigScaffoldFields(asmStore,
                            getNumASM_SCFs(asmStore->scfStore) + 1,
                            value.IID,
                            FALSE,
                            &offset,
                            (ssm->contig_pairs[i].orient == AS_NORMAL ||
                             ssm->contig_pairs[i].orient == AS_INNIE) ?
                            AS_FORWARD : AS_REVERSE);
    gap.asmMean = ssm->contig_pairs[i].mean;
    gap.asmStddev = ssm->contig_pairs[i].stddev;
    if(gap.asmMean < 20.0)
      gap.storeMean = 20;
    else
      gap.storeMean = (CDS_COORD_t) ssm->contig_pairs[i].mean;
    appendASM_GapStore(asmStore->scgStore, &gap);
    offset += gap.storeMean;
  }

  i -= (i == 0) ? 0 : 1;
  // look up the cco
  if(HASH_FAILURE == LookupInPHashTable_AS(asmStore->hashTable,
                                           ASM_UID_NAMESPACE,
                                           ssm->contig_pairs[i].econtig2,
                                           &value))
  {
    fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n",
            ssm->contig_pairs[i].econtig2);

  }
  appendASM_IIDStore(asmStore->sccStore, &(value.IID));
  
  SetContigScaffoldFields(asmStore,
                          getNumASM_SCFs(asmStore->scfStore) + 1,
                          value.IID,
                          FALSE,
                          &offset,
                          (ssm->contig_pairs[i].orient == AS_NORMAL ||
                           ssm->contig_pairs[i].orient == AS_OUTTIE) ?
                          AS_FORWARD : AS_REVERSE);
  appendASM_SCFStore(asmStore->scfStore, &scf);
  
  memset(&value, 0, sizeof(PHashValue_AS));
  value.type = AS_IID_SCF;
  InsertInPHashTable_AS(&(asmStore->hashTable), ASM_UID_NAMESPACE,
                        scf.uid, &value, FALSE, TRUE);

  {
    ASM_IIDRecord contigIID;
    ASM_CCORecord cco;
    
    // print out the length of the scaffold
    getASM_IIDStore(asmStore->sccStore, scf.firstContig + scf.numContigs - 1,
                    &contigIID);
    getASM_CCOStore(asmStore->ccoStore, contigIID, &cco);

    /*
    fprintf(stdout, F_UID " " F_COORD,
            scf.uid, MAX(cco.scaffoldPos.bgn, cco.scaffoldPos.end));
    */
  }
  return 0;
}


int AddGenericMesg2Store(AssemblyStore * asmStore, GenericMesg * gen)
{
  static int numMDI = 0;
  static int numAFG = 0;
  static int numUTG = 0;
  static int numCCO = 0;
  static int numDSC = 0;
  static int numSCF = 0;
  switch(gen->t)
  {
    case MESG_MDI:
      if(numMDI == 0)
        fprintf(stderr, "First MDI message.\n");
      AddMDI2Store(asmStore, gen->m);
      numMDI++;
      break;
    case MESG_AFG:
      if(numAFG == 0)
      {
        fprintf(stderr, "\n%d MDI messages.\n", numMDI);
        fprintf(stderr, "First AFG message.\n");
      }
      AddAFG2Store(asmStore, gen->m);
      numAFG++;
      break;
    case MESG_UTG:
      if(numUTG == 0)
      {
        fprintf(stderr, "\n%d AFG messages.\n", numAFG);
        fprintf(stderr, "First UTG message.\n");
      }
      AddUTG2Store(asmStore, gen->m);
      numUTG++;
      break;
    case MESG_CCO:
      if(numCCO == 0)
      {
        fprintf(stderr, "\n%d UTG messages.\n", numUTG);
        fprintf(stderr, "First CCO message.\n");
      }
      AddCCO2Store(asmStore, gen->m);
      numCCO++;
      break;
    case MESG_DSC:
      if(numDSC == 0)
        fprintf(stderr, "First DSC message.\n");
      AddDSC2Store(asmStore, gen->m);
      numDSC++;
      break;
    case MESG_SCF:
      if(numSCF == 0)
      {
        fprintf(stderr, "\n%d CCO messages.\n", numCCO);
        fprintf(stderr, "\n%d DSC messages.\n", numDSC);
        fprintf(stderr, "First SCF message.\n");
      }
      AddSCF2Store(asmStore, gen->m);
      numSCF++;
      break;
    default:
      break;
  }
  return 0;
}


void CopyGateKeeperLKGStore(AssemblyStore * asmStore)
{
  int32 numLNKs;
  int32 i;
  GateKeeperLinkRecord lnk;

  if(asmStore->gkpStore == NULL)
    return;

  numLNKs = getNumGateKeeperLinks(asmStore->gkpStore->lnkStore);
  assert(sizeof(GateKeeperLinkRecord) == sizeof(ASM_LKGRecord));
  assert(asmStore->gkpStore != NULL);
  for(i = 1; i <= numLNKs; i++)
  {
    getGateKeeperLinkStore(asmStore->gkpStore->lnkStore, i, &lnk);
    appendASM_LKGStore(asmStore->lkgStore, &lnk);
  }
}


void SetUnitigAndFragmentScaffoldCoordinates(AssemblyStore * asmStore)
{
  ASM_CCORecord cco;
  ASM_UTGRecord utg;
  ASM_InstanceRecord ucIns;
  ASM_InstanceRecord usIns;
  int32 numUTGs = getNumASM_UTGs(asmStore->utgStore);
  int32 i;
  int32 ucIndex;

  // all unitigs
  for(i = 1; i <= numUTGs; i++)
  {
    // get the unitig
    getASM_UTGStore(asmStore->utgStore, i, &utg);

    // all instances of the unitig in contigs
    for(ucIndex = utg.cInsIndex; ucIndex != 0; ucIndex = ucIns.next)
    {
      // get the unitig instance in the contig
      getASM_InstanceStore(asmStore->uciStore, ucIndex, &ucIns);
      // get the contig
      getASM_CCOStore(asmStore->ccoStore, ucIns.containerIndex, &cco);

      utg.inDegenerate = cco.inDegenerate;
      
      // populate the unitig's coords in the scaffold from
      // unitig in contig & contig in scaffold
      PopulateDerivedInstance(&usIns,
                              cco.scaffoldIndex,
                              cco.scaffoldPos, ucIns.pos);

      // add an instance for every fragment in this unitig
      AppendNewUnitigFragmentInstances(asmStore,
                                       cco.inDegenerate,
                                       utg.status,
                                       utg.firstFrag, utg.numFrags,
                                       ucIns.containerIndex,
                                       ucIns.pos,
                                       cco.scaffoldIndex,
                                       usIns.pos);

      // append the unitig instance in the scaffold
      AppendInstance(asmStore->usiStore, usIns, &(utg.sInsIndex));
    }
    // update the unitig
    setASM_UTGStore(asmStore->utgStore, i, &utg);
  }
}


void InitializeAFGStore(AssemblyStore * asmStore)
{
  int32 i;
  int32 numFRGs;
  GateKeeperFragmentRecord gkfr;
  ASM_AFGRecord afg;
  PHashValue_AS value;

  if(asmStore->gkpStore == NULL)
    return;
  numFRGs = getNumGateKeeperFragments(asmStore->gkpStore->frgStore);
  memset(&afg, 0, sizeof(ASM_AFGRecord));
  
  for(i = 1; i <= numFRGs; i++)
  {
    getGateKeeperFragmentStore(asmStore->gkpStore->frgStore, i, &gkfr);
    afg.uid = gkfr.readUID;
    afg.deleted = gkfr.deleted;
    appendASM_AFGStore(asmStore->afgStore, &afg);
    
    memset(&value, 0, sizeof(PHashValue_AS));
    value.type = AS_IID_AFG;
    InsertInPHashTable_AS(&(asmStore->hashTable), ASM_UID_NAMESPACE,
                          afg.uid, &value, FALSE, TRUE);
  }
}


AssemblyStore * CreateAssemblyStoreFromASMFile(FILE * fi,
                                               char * storePath,
                                               char * gkpStorePath,
                                               char * frgStorePath)
{
  AssemblyStore * asmStore;
  GenericMesg * gen;
  unsigned long mesgCount = 0;

  assert(fi != NULL &&
         storePath != NULL);
  
  asmStore = CreateAssemblyStore(storePath, gkpStorePath, frgStorePath);

  if(asmStore->gkpStore != NULL)
  {
    fprintf(stderr, "Sync'ing assembly store's AFGs with gatekeepers FRGs\n");
    InitializeAFGStore(asmStore);
  }

  fprintf(stderr, "Reading .asm file\n");
  while(ReadProtoMesg_AS(fi, &gen) != EOF)
  {
#ifndef DEBUG_ASMDATA
    if(++mesgCount % 10000 == 0)
      fprintf(stderr, "\r%20lu", mesgCount);
#endif
    if(AddGenericMesg2Store(asmStore, gen))
    {
      CloseAssemblyStore(asmStore);
      return NULL;
    }
  }

  if(asmStore->gkpStore != NULL)
  {
    fprintf(stderr, "Copying gatekeeeper store mate links\n");
    CopyGateKeeperLKGStore(asmStore);
  }

  fprintf(stderr, "Setting unitig and fragment scaffold coordinates\n");
  SetUnitigAndFragmentScaffoldCoordinates(asmStore);
  
  return asmStore;
}

int PopulateAndAppendMatePair(CloneData * cd,
                              ASM_MatePair * mp,
                              ASM_AFGRecord * afg1,
                              ASM_InstanceRecord * ins1,
                              ASM_AFGRecord * afg2,
                              ASM_InstanceRecord * ins2,
                              ASM_MDIRecord * mdi)
{
  ASM_AFGRecord * leftAFG;
  ASM_AFGRecord * rightAFG;
  ASM_InstanceRecord * leftIns;
  ASM_InstanceRecord * rightIns;

  if(ins1->pos.bgn < ins2->pos.bgn)
  {
    leftAFG = afg1;
    leftIns = ins1;
    rightAFG = afg2;
    rightIns = ins2;
  }
  else
  {
    leftAFG = afg2;
    leftIns = ins2;
    rightAFG = afg1;
    rightIns = ins1;
  }
  
  mp->leftUID = leftAFG->uid;
  mp->rightUID = rightAFG->uid;
  mp->fivePrimes.bgn = leftIns->pos.bgn;
  mp->fivePrimes.end = rightIns->pos.bgn;
  if(leftIns->pos.bgn < leftIns->pos.end)
  {
    if(rightIns->pos.bgn < rightIns->pos.end)
    {
      mp->orient = AS_NORMAL;
      AppendVA_ASM_MatePair(cd->normal, mp);
    }
    else
    {
      mp->orient = AS_INNIE;
      mp->basePairs = (rightIns->pos.bgn - leftIns->pos.bgn) - mdi->asmMean;
      mp->stddevs = (mp->basePairs) / mdi->asmStddev;
      AppendVA_ASM_MatePair(cd->innie, mp);
    }
  }
  else
  {
    if(rightIns->pos.bgn < rightIns->pos.end)
    {
      mp->orient = AS_OUTTIE;
      AppendVA_ASM_MatePair(cd->outtie, mp);
    }
    else
    {
      mp->orient = AS_ANTI;
      AppendVA_ASM_MatePair(cd->antinormal, mp);
    }
  }
  return 0;
}

void PopulateDCMatePair(ASM_MatePair * mp,
                        ASM_AFGRecord * afg1,
                        ASM_InstanceRecord * ins1,
                        ASM_AFGRecord * afg2,
                        ASM_InstanceRecord * ins2)
{
  mp->leftUID = afg1->uid;
  mp->rightUID = afg2->uid;
  mp->fivePrimes.bgn = ins1->pos.bgn;
  mp->fivePrimes.end = ins2->pos.bgn;
  if(ins1->pos.bgn < ins1->pos.end)
  {
    if(ins2->pos.bgn < ins2->pos.end)
      mp->orient = AS_NORMAL;
    else
      mp->orient = AS_INNIE;
  }
  else
  {
    if(ins2->pos.bgn < ins2->pos.end)
      mp->orient = AS_OUTTIE;
    else
      mp->orient = AS_ANTI;
  }
}


void AddContigFragmentsToCloneData(AssemblyStore * asmStore,
                                   CloneData * cd,
                                   ASM_IIDRecord contigIndex,
                                   int canonicalOnly)

{
  int i;
  ASM_AFGRecord afg1;
  ASM_LKGRecord lkg;
  ASM_CCORecord cco;
  ASM_SCFRecord scf;
  ASM_DSCRecord dsc;
  ASM_IIDRecord iid1;
  ASM_IIDRecord iid2;
  ASM_InstanceRecord ins1;
  ASM_InstanceRecord ins2;
  ASM_MatePair mp;
  ASM_MDIRecord mdi;
  ASM_AFGRecord afg2;
          
  getASM_CCOStore(asmStore->ccoStore, contigIndex, &cco);

  // iterate over fragments in this contig
  for(i = cco.firstFrag; i < cco.firstFrag + cco.numFrags; i++)
  {
    // get the fragment record
    getASM_IIDStore(asmStore->ccfStore, i, &iid1);
    getASM_AFGStore(asmStore->afgStore, iid1, &afg1);

    // only work with mated fragments
    if(afg1.numLinks > 0)
    {
      GateKeeperLinkRecordIterator iterator;
      CreateGateKeeperLinkRecordIterator(asmStore->lkgStore, afg1.linkHead,
                                         iid1, &iterator);

      // iterate over mates (should be only one)
      while(NextGateKeeperLinkRecordIterator(&iterator, &lkg))
      {
        // ignore non-clone type links
        if(lkg.type == AS_MATE)
        {
          // get the mate record
          iid2 = (lkg.frag1 == iid1) ? lkg.frag2 : lkg.frag1;
          getASM_AFGStore(asmStore->afgStore, iid2, &afg2);

          if(canonicalOnly && afg1.uid > afg2.uid)
            break;
          
          memset(&mp, 0, sizeof(ASM_MatePair));
          
          // get the fragment's scaffold instance & scaffold records
          getASM_InstanceStore(asmStore->asiStore, afg1.sInsIndex, &ins1);
          if(cco.inDegenerate)
          {
            getASM_DSCStore(asmStore->dscStore, ins1.containerIndex, &dsc);
            mp.containerUID = dsc.uid;
          }
          else
          {
            getASM_SCFStore(asmStore->scfStore, ins1.containerIndex, &scf);
            mp.containerUID = scf.uid;
          }
          
          // get the distance record
          getASM_MDIStore(asmStore->mdiStore, lkg.distance, &mdi);
          mp.distUID = mdi.uid;

          if(afg2.chaff || afg2.sInsIndex == 0)
          {
            // if the mate is chaff, register it as such & move on
            mp.leftUID = afg1.uid;
            mp.rightUID = afg2.uid;
            mp.fivePrimes.bgn = ins1.pos.bgn;
            mp.fivePrimes.end = 0;
            mp.orient = (ins1.pos.bgn < ins1.pos.end) ? AS_NORMAL : AS_ANTI;
            if(afg2.chaff)
              AppendVA_ASM_MatePair(cd->inChaff, &mp);
            else
              AppendVA_ASM_MatePair(cd->missing, &mp); // shouldn't happen!
            break;
          }
          else
          {
            // mate is uniquely placed
            getASM_InstanceStore(asmStore->asiStore, afg2.sInsIndex, &ins2);
          }

          // if the mate is multiply placed or is in a different scaffold
          if((!afg1.inDegenerate && afg2.inDegenerate) ||
             (afg1.inDegenerate && !afg2.inDegenerate) ||
             ins2.next != 0 ||
             ins1.containerIndex != ins2.containerIndex)
          {
            PopulateDCMatePair(&mp, &afg1, &ins1, &afg2, &ins2);
            if(afg2.inDegenerate)
              AppendVA_ASM_MatePair(cd->inDegenerate, &mp);
            else if(ins2.next != 0)
              AppendVA_ASM_MatePair(cd->inSurrogate, &mp);
            else
              AppendVA_ASM_MatePair(cd->elsewhere, &mp);
          }
          else if(afg1.uid < afg2.uid)
          {
            // this is the case of interest:
            // singly placed mate in the same scaffold
            PopulateAndAppendMatePair(cd, &mp,
                                      &afg1, &ins1, &afg2, &ins2, &mdi);
          } // else mate is in same scaffold & not in surrogate
          break;
        } // real link
      } // iteration over frag1's links
    } // if frag1 has links
  } // loop over contig fragment indexes
}


CloneData * CreateCloneData(void)
{
  CloneData * cd;
  cd = (CloneData *) calloc(1, sizeof(CloneData));
  cd->innie = CreateVA_ASM_MatePair(10);
  cd->normal = CreateVA_ASM_MatePair(10);
  cd->antinormal = CreateVA_ASM_MatePair(10);
  cd->outtie = CreateVA_ASM_MatePair(10);
  cd->elsewhere = CreateVA_ASM_MatePair(10);
  cd->inDegenerate = CreateVA_ASM_MatePair(10);
  cd->inSurrogate = CreateVA_ASM_MatePair(10);
  cd->inChaff = CreateVA_ASM_MatePair(10);
  cd->missing = CreateVA_ASM_MatePair(10);
  return cd;
}


CloneData * GetScaffoldCloneData(AssemblyStore * asmStore,
                                 int32 scaffoldIndex,
                                 int canonicalOnly,
                                 int isDegenerate)
{
  CloneData * cd;
  ASM_SCFRecord scf;
  ASM_DSCRecord dsc;
  ASM_IIDRecord contigIndex;
  int i;

  cd = CreateCloneData();
  if(isDegenerate)
  {
    getASM_DSCStore(asmStore->dscStore, scaffoldIndex, &dsc);
    AddContigFragmentsToCloneData(asmStore, cd, dsc.contigIndex,
                                  canonicalOnly);
  }
  else
  {
    getASM_SCFStore(asmStore->scfStore, scaffoldIndex, &scf);

    for(i = scf.firstContig; i < scf.firstContig + scf.numContigs; i++)
    {
      getASM_IIDStore(asmStore->sccStore, i, &contigIndex);
      AddContigFragmentsToCloneData(asmStore, cd, contigIndex,
                                    canonicalOnly);
    }
  }
  return cd;
}


void PrintInnieMatePairs(VA_TYPE(ASM_MatePair) * mps, FILE * fo)
{
  int i;
  ASM_MatePair * myMPs = GetVA_ASM_MatePair(mps, 0);
  int32 numMPs = GetNumVA_ASM_MatePair(mps);
  
  for(i = 0; i < numMPs; i++)
  {
    fprintf(fo,
            "%c\t" F_UID "\t" F_UID "\t" F_UID "\t" F_COORD "\t" F_COORD "\n",
            (char) AS_INNIE, myMPs[i].leftUID, myMPs[i].rightUID,
            myMPs[i].distUID,
            myMPs[i].fivePrimes.bgn, myMPs[i].fivePrimes.end);
  }
}


void PrintMatePairs(VA_TYPE(ASM_MatePair) * mps,
                    OrientType orient, FILE * fo)
{
  int i;
  ASM_MatePair * myMPs = GetVA_ASM_MatePair(mps, 0);
  int32 numMPs = GetNumVA_ASM_MatePair(mps);
  
  for(i = 0; i < numMPs; i++)
  {
    fprintf(fo, "%c\t" F_UID "\t" F_UID "\t" F_UID "\t" F_COORD "\t" F_COORD "\n",
            (char) orient,
            myMPs[i].leftUID, myMPs[i].rightUID, myMPs[i].distUID,
            myMPs[i].fivePrimes.bgn, myMPs[i].fivePrimes.end);
  }
}


void PrintPartialMatePairs(VA_TYPE(ASM_MatePair) * mps,
                           char type, FILE * fo)
{
  int i;
  ASM_MatePair * myMPs = GetVA_ASM_MatePair(mps, 0);
  int32 numMPs = GetNumVA_ASM_MatePair(mps);
  
  for(i = 0; i < numMPs; i++)
  {
    fprintf(fo, "%c\t" F_UID "\t" F_UID "\t" F_UID "\t" F_UID "\t" F_COORD "\n",
            type, myMPs[i].containerUID,
            myMPs[i].leftUID, myMPs[i].rightUID, myMPs[i].distUID,
            myMPs[i].fivePrimes.bgn);
  }
}


void PrintBadCloneData(Scaffold_ID containerUID,
                       CloneData * cd, FILE * fo)
{
  PrintMatePairs(cd->normal, AS_NORMAL, fo);
  PrintMatePairs(cd->antinormal, AS_ANTI, fo);
  PrintMatePairs(cd->outtie, AS_OUTTIE, fo);
}


void PrintMissingCloneData(Scaffold_ID containerUID,
                           CloneData * cd, FILE * fo)
{
  PrintPartialMatePairs(cd->elsewhere, 'E', fo);
  PrintPartialMatePairs(cd->inDegenerate, 'D', fo);
  PrintPartialMatePairs(cd->inSurrogate, 'S', fo);
  PrintPartialMatePairs(cd->inChaff, 'C', fo);
  PrintPartialMatePairs(cd->missing, 'M', fo);
}


void WriteBinaryMatePairs(char * fpref, char type, 
                          VA_TYPE(ASM_MatePair) * mps)
{
  char filename[4096];
  FILE * fp;
  ASM_MatePair * mp0 = GetVA_ASM_MatePair(mps, 0);
  int32 numMPs = GetNumVA_ASM_MatePair(mps);
  int numIW;
  
  sprintf(filename, "%s_%c.bdat", fpref, type);
  fp = fopen(filename, "w");
  assert(fp != NULL);

  numIW = fwrite(&numMPs, sizeof(numMPs), 1, fp);
  assert(numIW == 1);
  
  numIW = fwrite(mp0, sizeof(ASM_MatePair), numMPs, fp);
  assert(numIW == numMPs);
  
  fclose(fp);
}


void WriteBinaryCloneData(Scaffold_ID containerUID, CloneData * cd)
{
  char fpref[4096];
  sprintf(fpref, F_UID, containerUID);
  WriteBinaryMatePairs(fpref, 'I', cd->innie);
  WriteBinaryMatePairs(fpref, 'N', cd->normal);
  WriteBinaryMatePairs(fpref, 'A', cd->antinormal);
  WriteBinaryMatePairs(fpref, 'O', cd->outtie);
}


void PrintIntraCloneData(Scaffold_ID containerUID, CloneData * cd, FILE * fo)
{
  PrintInnieMatePairs(cd->innie, fo);
  PrintMatePairs(cd->normal, AS_NORMAL, fo);
  PrintMatePairs(cd->antinormal, AS_ANTI, fo);
  PrintMatePairs(cd->outtie, AS_OUTTIE, fo);
}


void PrintCloneData(Scaffold_ID containerUID,
                    CloneData * cd, char * which, FILE * fo)
{
  if(which[AS_INNIE])
    PrintInnieMatePairs(cd->innie, fo);
  if(which[AS_NORMAL])
    PrintMatePairs(cd->normal, AS_NORMAL, fo);
  if(which[AS_ANTI])
    PrintMatePairs(cd->antinormal, AS_ANTI, fo);
  if(which[AS_OUTTIE])
    PrintMatePairs(cd->outtie, AS_OUTTIE, fo);
  if(which['E'])
    PrintPartialMatePairs(cd->elsewhere, 'E', fo);
  if(which['D'])
    PrintPartialMatePairs(cd->inDegenerate, 'D', fo);
  if(which['S'])
    PrintPartialMatePairs(cd->inSurrogate, 'S', fo);
  if(which['C'])
    PrintPartialMatePairs(cd->inChaff, 'C', fo);
  if(which['M'])
    PrintPartialMatePairs(cd->missing, 'M', fo);
}


void DeleteCloneData(CloneData * cd)
{
  if(cd != NULL)
  {
    Delete_VA(cd->innie);
    Delete_VA(cd->normal);
    Delete_VA(cd->antinormal);
    Delete_VA(cd->outtie);
    Delete_VA(cd->elsewhere);
    Delete_VA(cd->inDegenerate);
    Delete_VA(cd->inSurrogate);
    Delete_VA(cd->inChaff);
    Delete_VA(cd->missing);
    free(cd);
  }
}


void PrintUnreferencedFrags(AssemblyStore * asmStore, FILE * fo)
{
  int32 i;
  int32 maxI;
  ASM_AFGRecord afg;
  ASM_UTGRecord utg;
  ASM_CCORecord cco;
  ASM_DSCRecord dsc;
  ASM_SCFRecord scf;
  ASM_InstanceRecord ins;

  maxI = getNumASM_AFGs(asmStore->afgStore);
  for(i = 1; i <= maxI; i++)
  {
    getASM_AFGStore(asmStore->afgStore, i, &afg);
    if(afg.deleted) continue;
    if(afg.unreferenced)
    {
      fprintf(fo, F_UID " %c", afg.uid, (char) afg.type);
      if(afg.unitigIndex != 0)
      {
         getASM_UTGStore(asmStore->utgStore, afg.unitigIndex, &utg);
         fprintf(fo, "\tU %c\t" F_UID, (char) utg.status, utg.uid);
      }
      if(afg.cInsIndex != 0)
      {
        getASM_InstanceStore(asmStore->aciStore, afg.cInsIndex, &ins);
        getASM_CCOStore(asmStore->ccoStore, ins.containerIndex, &cco);
        fprintf(fo, "\tC " F_UID, cco.uid);
      }
      if(afg.sInsIndex != 0)
      {
        getASM_InstanceStore(asmStore->asiStore, afg.sInsIndex, &ins);
        if(afg.inDegenerate)
        {
          getASM_DSCStore(asmStore->dscStore, ins.containerIndex, &dsc);
          fprintf(fo, "\tS D " F_UID, dsc.uid);
        }
        else
        {
          getASM_SCFStore(asmStore->scfStore, ins.containerIndex, &scf);
          fprintf(fo, "\tS R " F_UID, scf.uid);
        }
      }
      fprintf(fo, "\n");
    }
  }
}

   
void PrintMissingFrags(AssemblyStore * asmStore, FILE * fo)
{
  int32 i;
  int32 maxI;
  ASM_AFGRecord afg;
  ASM_UTGRecord utg;
  ASM_CCORecord cco;
  ASM_DSCRecord dsc;
  ASM_SCFRecord scf;
  ASM_InstanceRecord ins;

  maxI = getNumASM_AFGs(asmStore->afgStore);
  for(i = 1; i <= maxI; i++)
  {
    getASM_AFGStore(asmStore->afgStore, i, &afg);
    if(afg.deleted) continue;
    if(!afg.chaff && (afg.cInsIndex == 0 || afg.sInsIndex == 0))
    {
      fprintf(fo, "%d " F_UID " %c", i, afg.uid, (char) afg.type);
      if(afg.unitigIndex != 0)
      {
        getASM_UTGStore(asmStore->utgStore, afg.unitigIndex, &utg);
        fprintf(fo, "\tU %c\t" F_UID, (char) utg.status, utg.uid);
      }
      if(afg.cInsIndex != 0)
      {
        getASM_InstanceStore(asmStore->aciStore, afg.cInsIndex, &ins);
        getASM_CCOStore(asmStore->ccoStore, ins.containerIndex, &cco);
        fprintf(fo, "\tC " F_UID, cco.uid);
      }
      if(afg.sInsIndex != 0)
      {
        getASM_InstanceStore(asmStore->asiStore, afg.sInsIndex, &ins);
        if(afg.inDegenerate)
        {
          getASM_DSCStore(asmStore->dscStore, ins.containerIndex, &dsc);
          fprintf(fo, "\tS D " F_UID, dsc.uid);
        }
        else
        {
          getASM_SCFStore(asmStore->scfStore, ins.containerIndex, &scf);
          fprintf(fo, "\tS R " F_UID, scf.uid);
        }
      }
      fprintf(fo, "\n");
    }
  }
}


void PrintDeletedFrags(AssemblyStore * asmStore, FILE * fo)
{
  int32 i;
  int32 maxI;
  ASM_AFGRecord afg;

  maxI = getNumASM_AFGs(asmStore->afgStore);
  for(i = 1; i <= maxI; i++)
  {
    getASM_AFGStore(asmStore->afgStore, i, &afg);
    if(afg.deleted) fprintf(fo, F_UID "\n", afg.uid);
  }
}


void PrintDegenerateLength(AssemblyStore * asmStore, int32 index, FILE * fo)
{
  ASM_CCORecord cco;
  ASM_DSCRecord dsc;

  getASM_DSCStore(asmStore->dscStore, index, &dsc);
  getASM_CCOStore(asmStore->ccoStore, dsc.contigIndex, &cco);
  fprintf(fo, F_UID "\t" F_COORD "\n", dsc.uid, cco.length);
}


void getScaffoldLengths(AssemblyStore * asmStore, int32 index,
                        CDS_COORD_t * fastaLength, float32 * assemblyLength)
{
  ASM_SCFRecord scf;
  ASM_IIDRecord contigIndex;
  ASM_CCORecord cco;
  ASM_GapRecord gap;
  int i, j;

  *fastaLength = 0;
  *assemblyLength = 0.f;

  getASM_SCFStore(asmStore->scfStore, index, &scf);
  for(i = scf.firstContig, j = scf.firstGap;
      i < scf.firstContig + scf.numContigs - 1;
      i++, j++)
  {
    getASM_IIDStore(asmStore->sccStore, i, &contigIndex);
    getASM_CCOStore(asmStore->ccoStore, contigIndex, &cco);
    getASM_GapStore(asmStore->scgStore, j, &gap);
    *fastaLength += cco.length;
    *assemblyLength += cco.length;
    *fastaLength += gap.storeMean;
    *assemblyLength += gap.asmMean;
  }
  // get last contig
  getASM_IIDStore(asmStore->sccStore,
                   scf.firstContig + scf.numContigs - 1,
                   &contigIndex);
  getASM_CCOStore(asmStore->ccoStore, contigIndex, &cco);
  *fastaLength += cco.length;
  *assemblyLength += cco.length;
}


void PrintScaffoldLength(AssemblyStore * asmStore, int32 index, FILE * fo)
{
  ASM_SCFRecord scf;
  CDS_COORD_t fastaLength;
  float32 assemblyLength;

  getASM_SCFStore(asmStore->scfStore, index, &scf);
  getScaffoldLengths(asmStore, index, &fastaLength, &assemblyLength);
  fprintf(fo, F_UID "\t" F_COORD "\t%.2f\n",
          scf.uid, fastaLength, assemblyLength);
}




VA_TYPE(ASM_Quad) * IdentifyBadMateQuads(AssemblyStore * asmStore,
                                         CloneData * cd,
                                         char * fragTypes,
                                         BreakpointType bpType,
                                         float32 numStddevs)
{
  int32 i;
  PHashValue_AS value;
  ASM_MatePair * mp;
  ASM_AFGRecord afg;
  ASM_MDIRecord mdi;
  VA_TYPE(ASM_Quad) * quads = CreateVA_ASM_Quad(10);

  /*
    In the comments below,
    n = mean - numStddevs * stddev
    x = mean + numStddevs * stddev
   */
  switch(bpType)
  {
    case ASM_Stretched:
      /*
        STRETCHED INNIE MATES
        for each innie mate pair for which right 5' - left 5' > x
          
                  | |         2  3
        right bp  V |
                    |
                  n |   1
                    |
                  u |   0
                    |
                    --------------------------
                         ->   (  )
                        left bp
      */
      for(i = 0; i < GetNumVA_ASM_MatePair(cd->innie); i++)
      {
        mp = GetVA_ASM_MatePair(cd->innie, i);
        LookupInPHashTable_AS(asmStore->hashTable, ASM_UID_NAMESPACE,
                              mp->leftUID, &value);
        getASM_AFGStore(asmStore->afgStore, value.IID, &afg);
        if(fragTypes[afg.type] == 0)
          continue;
        
        LookupInPHashTable_AS(asmStore->hashTable, ASM_UID_NAMESPACE,
                              mp->distUID, &value);
        getASM_MDIStore(asmStore->mdiStore, value.IID, &mdi);

        if(mp->fivePrimes.end - mp->fivePrimes.bgn >
           mdi.asmMean + numStddevs * mdi.asmStddev)
        {
          ASM_Quad q;
          CDS_COORD_t n = mdi.asmMean - numStddevs * mdi.asmStddev;
          CDS_COORD_t x = mdi.asmMean + numStddevs * mdi.asmStddev;

          q.leftUID = afg.uid;
          // lower left
          q.pts[0].bgn = mp->fivePrimes.bgn;
          q.pts[0].end = mp->fivePrimes.end - x;

          // upper left
          q.pts[1].bgn = mp->fivePrimes.bgn;
          q.pts[1].end = mp->fivePrimes.end - n;

          // upper right
          q.pts[2].bgn = mp->fivePrimes.bgn + n;
          q.pts[2].end = mp->fivePrimes.end;

          // lower right
          q.pts[3].bgn = mp->fivePrimes.bgn + x;
          q.pts[3].end = mp->fivePrimes.end;

          AppendVA_ASM_Quad(quads, &q);
        }
      }
      break;
    case ASM_Compressed:
      /*
        COMPRESSED INNIE MATES
        for each innie mate pair for which right 5' - left 5' < n
          
                    |
       right end of |   1      2
        insertion   |
                    |
                    |   0      3
                    |
                    |
                    --------------------------
                      left bp
      */
      for(i = 0; i < GetNumVA_ASM_MatePair(cd->innie); i++)
      {
        mp = GetVA_ASM_MatePair(cd->innie, i);
        LookupInPHashTable_AS(asmStore->hashTable, ASM_UID_NAMESPACE,
                              mp->leftUID, &value);
        getASM_AFGStore(asmStore->afgStore, value.IID, &afg);
        if(fragTypes[afg.type] == 0)
          continue;
        
        LookupInPHashTable_AS(asmStore->hashTable, ASM_UID_NAMESPACE,
                              mp->distUID, &value);
        getASM_MDIStore(asmStore->mdiStore, value.IID, &mdi);

        if(mp->fivePrimes.end - mp->fivePrimes.bgn <
           mdi.asmMean - numStddevs * mdi.asmStddev)
        {
          ASM_Quad q;
          CDS_COORD_t n = mdi.asmMean - numStddevs * mdi.asmStddev;
          CDS_COORD_t x = mdi.asmMean + numStddevs * mdi.asmStddev;

          q.leftUID = afg.uid;
          // lower left
          q.pts[0].bgn = mp->fivePrimes.bgn;
          q.pts[0].end = mp->fivePrimes.bgn + n;

          // upper left
          q.pts[1].bgn = mp->fivePrimes.bgn;
          q.pts[1].end = mp->fivePrimes.bgn + x;

          // upper right
          q.pts[2].bgn = mp->fivePrimes.end;
          q.pts[2].end = mp->fivePrimes.bgn + x;

          // lower right
          q.pts[3].bgn = mp->fivePrimes.end;
          q.pts[3].end = mp->fivePrimes.bgn + n;

          AppendVA_ASM_Quad(quads, &q);
        }
      }
      break;
    case ASM_Normal:
      /*
        NORMAL MATES
          
                  n |  1
        right bp    |
                  u |  0
                    |
                  ^ |
                  | |         3  2
                    |
                    --------------------------
                         ->   (  )
                        left bp
      */
      for(i = 0; i < GetNumVA_ASM_MatePair(cd->normal); i++)
      {
        mp = GetVA_ASM_MatePair(cd->normal, i);
        LookupInPHashTable_AS(asmStore->hashTable, ASM_UID_NAMESPACE,
                              mp->leftUID, &value);
        getASM_AFGStore(asmStore->afgStore, value.IID, &afg);
        if(fragTypes[afg.type] == 0)
          continue;
        
        LookupInPHashTable_AS(asmStore->hashTable, ASM_UID_NAMESPACE,
                              mp->distUID, &value);
        getASM_MDIStore(asmStore->mdiStore, value.IID, &mdi);

        {
          ASM_Quad q;
          CDS_COORD_t n = mdi.asmMean - numStddevs * mdi.asmStddev;
          CDS_COORD_t x = mdi.asmMean + numStddevs * mdi.asmStddev;

          q.leftUID = afg.uid;
          // lower left
          q.pts[0].bgn = mp->fivePrimes.bgn;
          q.pts[0].end = mp->fivePrimes.end + n;

          // upper left
          q.pts[1].bgn = mp->fivePrimes.bgn;
          q.pts[1].end = mp->fivePrimes.end + x;

          // upper right
          q.pts[2].bgn = mp->fivePrimes.bgn + x;
          q.pts[2].end = mp->fivePrimes.end;

          // lower right
          q.pts[3].bgn = mp->fivePrimes.bgn + n;
          q.pts[3].end = mp->fivePrimes.end;

          AppendVA_ASM_Quad(quads, &q);
        }
      }
      break;
    case ASM_Antinormal:
      /*
        ANTINORMAL MATES
          
                  | |  0  1
        right bp  v |
                    |
                  n |           2
                    |
                  u |           3
                    |
                    --------------------------
                        (  )   <-
                         left bp
      */
      for(i = 0; i < GetNumVA_ASM_MatePair(cd->antinormal); i++)
      {
        mp = GetVA_ASM_MatePair(cd->antinormal, i);
        LookupInPHashTable_AS(asmStore->hashTable, ASM_UID_NAMESPACE,
                              mp->leftUID, &value);
        getASM_AFGStore(asmStore->afgStore, value.IID, &afg);
        if(fragTypes[afg.type] == 0)
          continue;
        
        LookupInPHashTable_AS(asmStore->hashTable, ASM_UID_NAMESPACE,
                              mp->distUID, &value);
        getASM_MDIStore(asmStore->mdiStore, value.IID, &mdi);

        {
          ASM_Quad q;
          CDS_COORD_t n = mdi.asmMean - numStddevs * mdi.asmStddev;
          CDS_COORD_t x = mdi.asmMean + numStddevs * mdi.asmStddev;

          q.leftUID = afg.uid;
          // lower left
          q.pts[0].bgn = mp->fivePrimes.bgn - x;
          q.pts[0].end = mp->fivePrimes.end;

          // upper left
          q.pts[1].bgn = mp->fivePrimes.bgn - n;
          q.pts[1].end = mp->fivePrimes.end;

          // upper right
          q.pts[2].bgn = mp->fivePrimes.bgn;
          q.pts[2].end = mp->fivePrimes.end - n;

          // lower right
          q.pts[3].bgn = mp->fivePrimes.bgn;
          q.pts[3].end = mp->fivePrimes.end - x;

          AppendVA_ASM_Quad(quads, &q);
        }
      }
      break;
    case ASM_OuttieLR:
      /*
        OUTTIE MATES
          
                  n |           1
        right bp    |
                  u |           2
                    |
                  ^ |
                  | |   0  3
                    |
                    --------------------------
                        (  )   <-
                         left bp
      */
      for(i = 0; i < GetNumVA_ASM_MatePair(cd->outtie); i++)
      {
        mp = GetVA_ASM_MatePair(cd->outtie, i);
        LookupInPHashTable_AS(asmStore->hashTable, ASM_UID_NAMESPACE,
                              mp->leftUID, &value);
        getASM_AFGStore(asmStore->afgStore, value.IID, &afg);
        if(fragTypes[afg.type] == 0)
          continue;
        
        LookupInPHashTable_AS(asmStore->hashTable, ASM_UID_NAMESPACE,
                              mp->distUID, &value);
        getASM_MDIStore(asmStore->mdiStore, value.IID, &mdi);

        {
          ASM_Quad q;
          CDS_COORD_t n = mdi.asmMean - numStddevs * mdi.asmStddev;
          CDS_COORD_t x = mdi.asmMean + numStddevs * mdi.asmStddev;

          q.leftUID = afg.uid;
          // lower left
          q.pts[0].bgn = mp->fivePrimes.bgn - x;
          q.pts[0].end = mp->fivePrimes.end;

          // upper left
          q.pts[1].bgn = mp->fivePrimes.bgn;
          q.pts[1].end = mp->fivePrimes.end + x;

          // upper right
          q.pts[2].bgn = mp->fivePrimes.bgn;
          q.pts[2].end = mp->fivePrimes.end + n;

          // lower right
          q.pts[3].bgn = mp->fivePrimes.bgn - n;
          q.pts[3].end = mp->fivePrimes.end;

          AppendVA_ASM_Quad(quads, &q);
        }
      }
      break;
    case ASM_OuttieM:
      /*
        Region between Outtie matepairs
          
                    |
       right frag   |   1      2
                    |
                    |
                    |   0      3
                    |
                    |
                    --------------------------
                      left frag
      */
      for(i = 0; i < GetNumVA_ASM_MatePair(cd->outtie); i++)
      {
        mp = GetVA_ASM_MatePair(cd->outtie, i);
        LookupInPHashTable_AS(asmStore->hashTable, ASM_UID_NAMESPACE,
                              mp->leftUID, &value);
        getASM_AFGStore(asmStore->afgStore, value.IID, &afg);
        if(fragTypes[afg.type] == 0)
          continue;
        
        {
          ASM_Quad q;
          
          q.leftUID = afg.uid;
          // lower left
          q.pts[0].bgn = mp->fivePrimes.bgn;
          q.pts[0].end = mp->fivePrimes.bgn;

          // upper left
          q.pts[1].bgn = mp->fivePrimes.bgn;
          q.pts[1].end = mp->fivePrimes.end;

          // upper right
          q.pts[2].bgn = mp->fivePrimes.end;
          q.pts[2].end = mp->fivePrimes.end;

          // lower right
          q.pts[3].bgn = mp->fivePrimes.end;
          q.pts[3].end = mp->fivePrimes.bgn;

          AppendVA_ASM_Quad(quads, &q);
        }
      }
      break;
    default:
      fprintf(stderr, "Bad breakpoint type %d\n", bpType);
      break;
  }
  return quads;
}
  

void PrintQuadrilaterals(VA_TYPE(ASM_Quad) * quads,
                         ASM_OutputType ot,
                         FILE * fo)
{
  int32 i;
  ASM_Quad * q;

  for(i = 0; i < GetNumVA_ASM_Quad(quads); i++)
  {
    q = GetVA_ASM_Quad(quads, i);

    switch(ot)
    {
      case ASM_CliqueFinderOutput:
        fprintf(fo, F_COORD " " F_COORD "\t" F_COORD " " F_COORD "\t" F_COORD " " F_COORD "\t" F_COORD " " F_COORD "\t" F_UID "\n",
                q->pts[0].bgn, q->pts[0].end,
                q->pts[1].bgn, q->pts[1].end,
                q->pts[2].bgn, q->pts[2].end,
                q->pts[3].bgn, q->pts[3].end,
                q->leftUID);
        break;
      case ASM_GnuplotOutput:
      default:
        fprintf(fo, F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n\n",
                q->pts[0].bgn, q->pts[0].end,
                q->pts[1].bgn, q->pts[1].end,
                q->pts[2].bgn, q->pts[2].end,
                q->pts[3].bgn, q->pts[3].end,
                q->pts[0].bgn, q->pts[0].end);
        break;
    }
  }
}


void CreateBMPFilename(char * fname,
                       CDS_UID_t uid,
                       char * fragTypes,
                       BreakpointType bpType,
                       float32 numStddevs)
{
  char fts[257];
  int i, j;

  for(i = 0, j = 0; i < 256; i++)
    if(fragTypes[i])
      fts[j++] = fragTypes[i];

  if(j == 256)
    sprintf(fts, "A");
  else
    fts[j] = '\0';

  switch(bpType)
  {
    case ASM_Stretched:
      sprintf(fname, F_UID "_S%.1f_%s.dat", uid, numStddevs, fts);
      break;
    case ASM_Compressed:
      sprintf(fname, F_UID "_C%.1f_%s.dat", uid, numStddevs, fts);
      break;
    case ASM_Normal:
      sprintf(fname, F_UID "_N_%s.dat", uid, fts);
      break;
    case ASM_Antinormal:
      sprintf(fname, F_UID "_A_%s.dat", uid, fts);
      break;
    case ASM_OuttieLR:
      sprintf(fname, F_UID "_OLR_%s.dat", uid, fts);
      break;
    case ASM_OuttieM:
      sprintf(fname, F_UID "_OM_%s.dat", uid, fts);
      break;
    default:
      fprintf(stderr, "Unknown breakpoint type %d\n", bpType);
      fname[0] = '\0';
      break;
  }
}


FILE * CreateGnuPlotFile(CDS_UID_t uid,
                         char * fragTypes,
                         BreakpointType bpType,
                         float32 numStddevs)
{
  FILE * fp;
  char fname[1024];
  char bpString[256];
  char xLabel[256];
  char yLabel[256];

  CreateBMPFilename(fname, uid, fragTypes, bpType, numStddevs);
  fp = fopen(fname, "w");
  assert(fp != NULL);

  switch(bpType)
  {
    case ASM_Stretched:
      sprintf(bpString, "stretched > %f stddevs", numStddevs);
      sprintf(xLabel, "Left-most end of required deletion");
      sprintf(yLabel, "Length of required deletion");
      break;
    case ASM_Compressed:
      sprintf(bpString, "compressed > %f stddevs", numStddevs);
      sprintf(xLabel, "Left-most end of required insertion");
      sprintf(yLabel, "Length of required insertion");
      break;
    case ASM_Normal:
      sprintf(bpString, "with normal orientation");
      xLabel[0] = '\0';
      yLabel[0] = '\0';
      break;
    case ASM_Antinormal:
      sprintf(bpString, "with antinormal orientation");
      xLabel[0] = '\0';
      yLabel[0] = '\0';
      break;
    case ASM_OuttieLR:
      sprintf(bpString, "left/right with outtie orientation");
      xLabel[0] = '\0';
      yLabel[0] = '\0';
      break;
    case ASM_OuttieM:
      sprintf(bpString, "middle with outtie orientation");
      xLabel[0] = '\0';
      yLabel[0] = '\0';
      break;
    default:
      assert(0);
      break;
  }

  fprintf(fp, "set title \"Mate pairs %s in scaffold " F_UID "\"\n",
          bpString, uid);
  fprintf(fp, "set xlabel \"%s\"\n", xLabel);
  fprintf(fp, "set ylabel \"%s\"\n", yLabel);

  return fp;
}


void PrintFastaFragmentCoordinates(AssemblyStore * asmStore, FILE * fo)
{
  int32 scfIndex;
  int32 numSCFs = getNumASM_SCFs(asmStore->scfStore);

  for(scfIndex = 1; scfIndex <= numSCFs; scfIndex++)
  {
    ASM_SCFRecord scf;
    int32 ccoIndex;
    CDS_COORD_t offset = 0;
    CDS_COORD_t lastCoord = 0;
    
    getASM_SCFStore(asmStore->scfStore, scfIndex, &scf);
    for(ccoIndex = scf.firstContig;
        ccoIndex < scf.firstContig + scf.numContigs;
        ccoIndex++)
    {
      ASM_IIDRecord contigIID;
      ASM_CCORecord cco;
      int32 afgIndex;
      CDS_COORD_t delta;
      
      getASM_IIDStore(asmStore->sccStore, ccoIndex, &contigIID);
      getASM_CCOStore(asmStore->ccoStore, contigIID, &cco);

      if(lastCoord != 0)
      {
        delta = MIN(cco.scaffoldPos.bgn, cco.scaffoldPos.end);
        delta = delta - lastCoord;
        delta = MAX(20, delta);
        offset += delta;
      }
      
      for(afgIndex = cco.firstFrag;
          afgIndex < cco.firstFrag + cco.numFrags;
          afgIndex++)
      {
        ASM_IIDRecord afgIID;
        ASM_AFGRecord afg;
        ASM_InstanceRecord cIns;
        
        getASM_IIDStore(asmStore->ccfStore, afgIndex, &afgIID);
        getASM_AFGStore(asmStore->afgStore, afgIID, &afg);

        if(afg.numLinks == 0 ||
           afg.chaff ||
           afg.inDegenerate ||
           afg.inSurrogate)
          continue;

        getASM_InstanceStore(asmStore->aciStore, afg.cInsIndex, &cIns);

        if(cco.scaffoldPos.bgn < cco.scaffoldPos.end)
        {
          fprintf(fo, F_UID " " F_UID " " F_COORD " " F_COORD "\n",
                  afg.uid, scf.uid,
                  offset + cIns.pos.bgn,
                  offset + cIns.pos.end);
        }
        else
        {
          fprintf(fo, F_UID " " F_UID " " F_COORD " " F_COORD "\n",
                  afg.uid, scf.uid,
                  offset + cco.length - cIns.pos.bgn,
                  offset + cco.length - cIns.pos.end);
        }
      }
      
    }
  }
}

void PrintScaffoldContigCoordinates(AssemblyStore * asmStore,
                                    uint32 index,
                                    int doGnuplotOutput,
                                    int doAssemblyCoords,
                                    FILE * fo)
{
  ASM_SCFRecord scf;
  ASM_CCORecord cco;
  ASM_IIDRecord contigIndex;
  ASM_GapRecord gap;
  CDS_COORD_t fastaLength;
  float32 assemblyLength;

  CDS_COORD_t offset = 0;
  int32 i, j;

  getASM_SCFStore(asmStore->scfStore, index, &scf);

  if(doGnuplotOutput)
    getScaffoldLengths(asmStore, index, &fastaLength, &assemblyLength);

  for(i = scf.firstContig, j = scf.firstGap;
      i < scf.firstContig + scf.numContigs - 1;
      i++, j++)
  {
    getASM_IIDStore(asmStore->sccStore, i, &contigIndex);
    getASM_CCOStore(asmStore->ccoStore, contigIndex, &cco);
    getASM_GapStore(asmStore->scgStore, j, &gap);

    if(doGnuplotOutput)
      fprintf(fo,
              F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n\n"
              F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n\n"
              F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n\n"
              F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n\n",
              offset, 0,
              offset,
              doAssemblyCoords ? (CDS_COORD_t) assemblyLength : fastaLength,
              offset + cco.length, 0,
              offset + cco.length,
              doAssemblyCoords ? (CDS_COORD_t) assemblyLength : fastaLength,
              0, offset,
              doAssemblyCoords ? (CDS_COORD_t) assemblyLength : fastaLength,
              offset, 0,
              offset + cco.length,
              doAssemblyCoords ? (CDS_COORD_t) assemblyLength : fastaLength,
              offset + cco.length);
    else
    {
      if(cco.scaffoldPos.bgn < cco.scaffoldPos.end)
        fprintf(fo, F_UID "\t" F_UID "\t" F_COORD "\t" F_COORD "\n",
                scf.uid, cco.uid, offset, offset + cco.length);
      else
        fprintf(fo, F_UID "\t" F_UID "\t" F_COORD "\t" F_COORD "\n",
                scf.uid, cco.uid, offset + cco.length, offset);
    }

    offset += cco.length;
    if(doAssemblyCoords)
      offset += gap.asmMean;
    else
      offset += gap.storeMean;
  }
  // get last contig
  getASM_IIDStore(asmStore->sccStore,
                   scf.firstContig + scf.numContigs - 1,
                   &contigIndex);
  getASM_CCOStore(asmStore->ccoStore, contigIndex, &cco);
  if(doGnuplotOutput)
      fprintf(fo,
              F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n\n"
              F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n\n"
              F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n\n"
              F_COORD " " F_COORD "\n" F_COORD " " F_COORD "\n\n",
              offset, 0,
              offset, doAssemblyCoords ? (CDS_COORD_t) assemblyLength : fastaLength,
              offset + cco.length, 0,
              offset + cco.length,
              doAssemblyCoords ? (CDS_COORD_t) assemblyLength : fastaLength,
              0, offset,
              doAssemblyCoords ? (CDS_COORD_t) assemblyLength : fastaLength,
              offset, 0,
              offset + cco.length,
              doAssemblyCoords ? (CDS_COORD_t) assemblyLength : fastaLength,
              offset + cco.length);
  else
    fprintf(fo, F_UID "\t" F_UID "\t" F_COORD "\t" F_COORD "\n",
            scf.uid, cco.uid,
            offset, offset + cco.length);
}


void PrintLibraries(AssemblyStore * asmStore, FILE * fo)
{
  ASM_MDIRecord mdi;
  int32 numMDIs = getNumASM_MDIs(asmStore->mdiStore);
  int i;
  for(i = 1; i <= numMDIs; i++)
  {
    getASM_MDIStore(asmStore->mdiStore, i, &mdi);
    fprintf(fo, F_UID " %.2f %.2f\n",
            mdi.uid, mdi.asmMean, mdi.asmStddev);
  }
}


MapStore * CreateMapStoreFromFiles(AssemblyStore * asmStore,
                                   FILE * chromFp,
                                   FILE * fragFp,
                                   char * mapStorePath)
{
  MapStore * mapStore = CreateMapStore(mapStorePath);
  char line[4096];

  // read in the chromosome file
  {
    ASM_CHRRecord chr;
    PHashValue_AS value;
    int numParsed;
    
    memset(&value, 0, sizeof(PHashValue_AS));
    
    while(fgets(line, 4095, chromFp))
    {
      if(line[0] == '#') continue;

      // format: chromosomeUID  chromosome[4]  arm[4]  description[256]
      numParsed = sscanf(line, F_UID " %s %s %s",
                         &chr.uid, chr.chromosome, chr.arm, chr.description);
      assert(numParsed == 4);
      chr.firstFrag = chr.lastFrag = 0;
      appendASM_CHRStore(mapStore->chrStore, &chr);
      
      value.type = AS_IID_CHR;
      InsertInPHashTable_AS(&(mapStore->hashTable), ASM_UID_NAMESPACE,
                            chr.uid, &value, FALSE, TRUE);
    }
  }

  // pre-populate the finStore to be in sync with asmStore's afgStore
  {
    int32 i;
    int32 numAFGs = getNumASM_AFGs(asmStore->afgStore);
    ASM_InstanceRecord ins;
    
    memset(&ins, 0, sizeof(ASM_InstanceRecord));
    for(i = 1; i <= numAFGs; i++)
      appendASM_InstanceStore(mapStore->finStore, &ins);
  }

  // read in fragment file
  {
    CDS_UID_t fragUID;
    CDS_UID_t chromUID;
    ASM_CHRRecord chr;
    ASM_InstanceRecord ins;
    PHashValue_AS value;
    int numParsed;
    int fragIndex;
    int chrIndex;
    CDS_COORD_t bgn, end;
    
    while(fgets(line, 4095, fragFp))
    {
      if(line[0] == '#') continue;

      // format: fragUID  chromosomeUID  5p  3p
      numParsed = sscanf(line, F_UID " " F_UID " " F_COORD " " F_COORD,
                         &fragUID, &chromUID, &bgn, &end);
      assert(numParsed == 4);

      // look up frag IID
      if(HASH_FAILURE == LookupInPHashTable_AS(asmStore->hashTable,
                                               ASM_UID_NAMESPACE,
                                               fragUID,
                                               &value))
      {
        fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n", fragUID);
        assert(0);
      }
      fragIndex = value.IID;
      
      // look up chr IID
      if(HASH_FAILURE == LookupInPHashTable_AS(mapStore->hashTable,
                                               ASM_UID_NAMESPACE,
                                               chromUID,
                                               &value))
      {
        fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n", chromUID);
        assert(0);
      }
      chrIndex = value.IID;
      
      // update the chromosome's list of fragments
      getASM_CHRStore(mapStore->chrStore, chrIndex, &chr);
      if(chr.lastFrag != 0)
      {
        assert(chr.firstFrag != 0);
        getASM_InstanceStore(mapStore->finStore, chr.lastFrag, &ins);
        ins.next = fragIndex;
        setASM_InstanceStore(mapStore->finStore, chr.lastFrag, &ins);
      }
      else
      {
        chr.firstFrag = fragIndex;
      }
      chr.lastFrag = fragIndex;
      setASM_CHRStore(mapStore->chrStore, chrIndex, &chr);

      // update this fragment instance
      getASM_InstanceStore(mapStore->finStore, fragIndex, &ins);
      ins.containerIndex = chrIndex;
      ins.pos.bgn = bgn;
      ins.pos.end = end;
      ins.next = 0;
      setASM_InstanceStore(mapStore->finStore, fragIndex, &ins);
    }
  }
  return mapStore;
}


void PrintMapStore(AssemblyStore * asmStore,
                   MapStore * mapStore,
                   FILE * fo)
{
  ASM_CHRRecord chr;
  int32 numIns = getNumASM_Instances(mapStore->finStore);
  ASM_InstanceRecord ins;
  ASM_AFGRecord afg;
  int32 i;

#ifdef NEVER
  int32 numChr = getNumASM_CHRs(mapStore->chrStore);
 for(i = 1; i <= numChr; i++)
  {
    getASM_CHRStore(mapStore->chrStore, i, &chr);
    fprintf(fo, F_UID " %s %s %s\n",
            chr.uid, chr.chromosome, chr.arm, chr.description);
  }
#endif

  for(i = 1; i <= numIns; i++)
  {
    getASM_InstanceStore(mapStore->finStore, i, &ins);
    if(ins.containerIndex == 0) continue;
    
    getASM_AFGStore(asmStore->afgStore, i, &afg);
    getASM_CHRStore(mapStore->chrStore, ins.containerIndex, &chr);
    fprintf(fo, F_UID " " F_UID " " F_COORD " " F_COORD "\n",
            afg.uid, chr.uid, ins.pos.bgn, ins.pos.end);
  }
}


void PrintChromosomeElsewheres(AssemblyStore * asmStore,
                               MapStore * mapStore,
                               Scaffold_ID containerUID,
                               FILE * fo)
{
  ASM_CHRRecord chr1;
  ASM_CHRRecord chr2;
  ASM_AFGRecord afg1;
  ASM_InstanceRecord ins1;
  PHashValue_AS value;
  ASM_LKGRecord lkg;
  GateKeeperLinkRecordIterator iterator;
  int nextIndex;
  int chrIndex;

  if(HASH_FAILURE == LookupInPHashTable_AS(mapStore->hashTable,
                                           ASM_UID_NAMESPACE,
                                           containerUID,
                                           &value))
  {
    fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n", containerUID);
    assert(0);
  }
  chrIndex = value.IID;
  getASM_CHRStore(mapStore->chrStore, chrIndex, &chr1);

  nextIndex = chr1.firstFrag;
  while(nextIndex != 0)
  {
    getASM_InstanceStore(mapStore->finStore, nextIndex, &ins1);
    getASM_AFGStore(asmStore->afgStore, nextIndex, &afg1);
    assert(ins1.containerIndex == chrIndex);
    
    // look up mate
    if(afg1.numLinks > 0)
    {
      CreateGateKeeperLinkRecordIterator(asmStore->lkgStore, afg1.linkHead,
                                         nextIndex, &iterator);
      while(NextGateKeeperLinkRecordIterator(&iterator, &lkg))
      {
        if(lkg.type == AS_MATE)
        {
          ASM_MDIRecord mdi;
          ASM_AFGRecord afg2;
          ASM_InstanceRecord ins2;
          
          CDS_IID_t iid = (lkg.frag1 == nextIndex) ? lkg.frag2 : lkg.frag1;

          getASM_AFGStore(asmStore->afgStore, iid, &afg2);
          getASM_InstanceStore(mapStore->finStore, iid, &ins2);
          getASM_MDIStore(asmStore->mdiStore, lkg.distance, &mdi);
          
          // are they in the same container?
          if(ins1.containerIndex != ins2.containerIndex &&
             ins2.containerIndex != 0)
          {
            // afg1 - uid, chrom, 5p, orient
            // afg2 - uid, chrom, 5p, orient
            // lib - uid
            getASM_CHRStore(mapStore->chrStore, ins2.containerIndex, &chr2);
            fprintf(fo, F_UID " " F_UID " " F_COORD " %s " F_UID " " F_UID " " F_COORD " %s " F_UID "\n",
                    afg1.uid, chr1.uid, ins1.pos.bgn,
                    (ins1.pos.bgn < ins1.pos.end) ? "A_B" : "B_A",
                    afg2.uid, chr2.uid, ins2.pos.bgn,
                    (ins2.pos.bgn < ins2.pos.end) ? "A_B" : "B_A",
                    mdi.uid);
          }
          break;
        }
      }
    }
    nextIndex = ins1.next;
  }
}


void PrintScaffoldElsewheres(AssemblyStore * asmStore,
                             int32 scaffoldIndex,
                             int canonicalOnly,
                             int isDegenerate,
                             FILE * fo)
{
  ASM_SCFRecord scf1;
  ASM_IIDRecord contigIndex;
  int s, i;

  if(isDegenerate)
  {
    // do nothing for now
  }
  else
  {
    getASM_SCFStore(asmStore->scfStore, scaffoldIndex, &scf1);

    for(s = scf1.firstContig; s < scf1.firstContig + scf1.numContigs; s++)
    {
      ASM_CCORecord cco;
      
      getASM_IIDStore(asmStore->sccStore, s, &contigIndex);
      getASM_CCOStore(asmStore->ccoStore, contigIndex, &cco);
      
      for(i = cco.firstFrag; i < cco.firstFrag + cco.numFrags; i++)
      {
        ASM_IIDRecord iid1;
        ASM_AFGRecord afg1;
        GateKeeperLinkRecordIterator iterator;
        ASM_LKGRecord lkg;
        
        getASM_IIDStore(asmStore->ccfStore, i, &iid1);
        getASM_AFGStore(asmStore->afgStore, iid1, &afg1);

        if(afg1.numLinks == 0 ||
           afg1.chaff ||
           afg1.inDegenerate ||
           afg1.inSurrogate)
          continue;

        CreateGateKeeperLinkRecordIterator(asmStore->lkgStore,
                                           afg1.linkHead,
                                           iid1, &iterator);

        // iterate over mates (should be only one)
        while(NextGateKeeperLinkRecordIterator(&iterator, &lkg))
        {
          // ignore non-clone type links
          if(lkg.type == AS_MATE)
          {
            // get the mate record
            ASM_IIDRecord iid2;
            ASM_AFGRecord afg2;
            ASM_InstanceRecord ins1;
            ASM_InstanceRecord ins2;
            ASM_MDIRecord mdi;
            
            iid2 = (lkg.frag1 == iid1) ? lkg.frag2 : lkg.frag1;
            getASM_AFGStore(asmStore->afgStore, iid2, &afg2);
            
            if(afg2.chaff ||
               afg2.sInsIndex == 0 ||
               afg2.inDegenerate ||
               afg2.inSurrogate)
              break;
            
            // get the distance record
            getASM_MDIStore(asmStore->mdiStore, lkg.distance, &mdi);
            
            // get the instances
            getASM_InstanceStore(asmStore->asiStore, afg1.sInsIndex, &ins1);
            getASM_InstanceStore(asmStore->asiStore, afg2.sInsIndex, &ins2);
            
            // if the mate is multiply placed or is in a different scaffold
            if(!afg2.inDegenerate && ins2.next == 0 &&
               ins1.containerIndex != ins2.containerIndex)
            {
              ASM_SCFRecord scf2;
              getASM_SCFStore(asmStore->scfStore,
                              ins2.containerIndex, &scf2);
              
              if(!canonicalOnly || scf1.uid > scf2.uid)
                fprintf(fo, F_UID " " F_UID " " F_COORD " %s " F_UID " " F_UID " " F_COORD " %s " F_UID "\n",
                        afg1.uid, scf1.uid, ins1.pos.bgn,
                        (ins1.pos.bgn < ins1.pos.end) ? "A_B" : "B_A",
                        afg2.uid, scf2.uid, ins2.pos.bgn,
                        (ins2.pos.bgn < ins2.pos.end) ? "A_B" : "B_A",
                        mdi.uid);
            }
            break;
          } // real link
        } // iteration over frag1's links
      } // loop over contig fragment indexes
    } // loop over scaffolds contigs
  } // if scaffold is not degenerate
}



CloneData * GetChromosomeCloneData(AssemblyStore * asmStore,
                                   MapStore * mapStore,
                                   Scaffold_ID containerUID)
{
  CloneData * cd;
  ASM_CHRRecord chr;
  ASM_AFGRecord afg1;
  ASM_InstanceRecord ins1;
  PHashValue_AS value;
  ASM_LKGRecord lkg;
  GateKeeperLinkRecordIterator iterator;
  int nextIndex;
  int chrIndex;

  cd = CreateCloneData();

  if(HASH_FAILURE == LookupInPHashTable_AS(mapStore->hashTable,
                                           ASM_UID_NAMESPACE,
                                           containerUID,
                                           &value))
  {
    fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n", containerUID);
    assert(0);
  }
  chrIndex = value.IID;
  getASM_CHRStore(mapStore->chrStore, chrIndex, &chr);

  nextIndex = chr.firstFrag;
  while(nextIndex != 0)
  {
    getASM_InstanceStore(mapStore->finStore, nextIndex, &ins1);
    getASM_AFGStore(asmStore->afgStore, nextIndex, &afg1);
    assert(ins1.containerIndex == chrIndex);
    
    // look up mate
    if(afg1.numLinks > 0)
    {
      CreateGateKeeperLinkRecordIterator(asmStore->lkgStore, afg1.linkHead,
                                         nextIndex, &iterator);
      while(NextGateKeeperLinkRecordIterator(&iterator, &lkg))
      {
        if(lkg.type == AS_MATE)
        {
          ASM_MDIRecord mdi;
          ASM_AFGRecord afg2;
          ASM_InstanceRecord ins2;
          ASM_MatePair mp;
          
          CDS_IID_t iid = (lkg.frag1 == nextIndex) ? lkg.frag2 : lkg.frag1;

          getASM_AFGStore(asmStore->afgStore, iid, &afg2);
          getASM_InstanceStore(mapStore->finStore, iid, &ins2);
          getASM_MDIStore(asmStore->mdiStore, lkg.distance, &mdi);
          mp.containerUID = containerUID;
          mp.distUID = mdi.uid;
          
          // are they in the same container?
          if(ins1.containerIndex != ins2.containerIndex)
          {
            PopulateDCMatePair(&mp, &afg1, &ins1, &afg2, &ins2);
            if(ins2.containerIndex == 0)
              AppendVA_ASM_MatePair(cd->missing, &mp);
            else
              AppendVA_ASM_MatePair(cd->elsewhere, &mp);
          }
          else if(afg1.uid < afg2.uid)
          {
            // this is the case of interest:
            // singly placed mate in the same container
            PopulateAndAppendMatePair(cd, &mp,
                                      &afg1, &ins1, &afg2, &ins2, &mdi);
          }
          break;
        }
      }
    }
    nextIndex = ins1.next;
  }
  return cd;
}


CDS_COORD_t GetScaffoldLength(AssemblyStore * asmStore, int32 index)
{
  ASM_SCFRecord scf;
  ASM_IIDRecord contigIndex;
  ASM_CCORecord cco;
  
  getASM_SCFStore(asmStore->scfStore, index, &scf);
  getASM_IIDStore(asmStore->sccStore,
                  scf.firstContig + scf.numContigs - 1,
                  &contigIndex);
  getASM_CCOStore(asmStore->ccoStore, contigIndex, &cco);
  return((cco.scaffoldPos.bgn < cco.scaffoldPos.end) ?
         cco.scaffoldPos.end : cco.scaffoldPos.bgn);
}


void PrintATACScaffoldGenomicAxis(AssemblyStore * asmStore,
                                  int32 index,
                                  char * parent,
                                  CDS_COORD_t * offset,
                                  FILE * fo)
{
  ASM_SCFRecord scf;
  CDS_COORD_t length;
  
  getASM_SCFStore(asmStore->scfStore, index, &scf);
  length = GetScaffoldLength(asmStore, index);
  fprintf(fo, "G s s" F_UID " %s " F_COORD " 0 . 0 " F_COORD " 0 0 0 0 " F_UID " " F_UID "\n",
          scf.uid, parent, *offset, length, scf.uid, scf.uid);
  *offset += length;
}


void PrintDeflineATACAxes(AssemblyStore * asmStore,
                          FILE * dfp,
                          char * parent,
                          FILE * fo)
{
  CDS_COORD_t offset = 0;
  char line[1024];
  PHashValue_AS value;
  
  while(fgets(line, 1023, dfp))
  {
    CDS_UID_t uid = STR_TO_UID(&line[1], NULL, 10);
    
    if(HASH_SUCCESS == LookupInPHashTable_AS(asmStore->hashTable,
                                             ASM_UID_NAMESPACE,
                                             uid, &value))
      PrintATACScaffoldGenomicAxis(asmStore, value.IID, parent, &offset, fo);
  }
}


void PrintATACSurrogates(AssemblyStore * asmStore,
                         char * parent,
                         FILE * fo)
{
  int i;
  int j;
  ASM_UTGRecord utg;
  ASM_SCFRecord scf;
  ASM_InstanceRecord scaffIns;
  int32 sInsIndex;
  int numUTGs = getNumASM_UTGs(asmStore->utgStore);

  for(i = 1; i <= numUTGs; i++)
  {
    getASM_UTGStore(asmStore->utgStore, i, &utg);
    if(utg.status == AS_SEP)
    {
      for(j = 0, sInsIndex = utg.sInsIndex;
          sInsIndex != 0;
          sInsIndex = scaffIns.next, j++)
      {
        getASM_InstanceStore(asmStore->usiStore, sInsIndex, &scaffIns);
        getASM_SCFStore(asmStore->scfStore, scaffIns.containerIndex, &scf);
        fprintf(fo, "F su su" F_UID "-%03d . %s " F_UID " " F_COORD " " F_COORD " %d 0\n",
                utg.uid, j,
                parent, scf.uid,
                MIN(scaffIns.pos.bgn, scaffIns.pos.end),
                abs(scaffIns.pos.end - scaffIns.pos.bgn),
                (scaffIns.pos.bgn < scaffIns.pos.end) ? 1 : -1);
      }
    }
  }
}


void InitializeATACFile(AssemblyStore * asmStore,
                        FILE * fo)
{
  // format version
  fprintf(fo, "# ATA format version\n");
  fprintf(fo, "! format ata 1.0\n");

#ifdef NEVER
  // types
  fprintf(fo, "# assembly data types\n");
  fprintf(fo, "T %d . chromosome .\n", AS_IID_CHR);
  fprintf(fo, "T %d . scaffold .\n", AS_IID_SCF);
  fprintf(fo, "T %d . degenerate .\n", AS_IID_DSC);
  fprintf(fo, "T %d . contig .\n", AS_IID_CCO);
  fprintf(fo, "T %d . unitig .\n", AS_IID_UTG);
  fprintf(fo, "T %d . fragment .\n", AS_IID_AFG);
  fprintf(fo, "T %d . distance .\n", AS_IID_MDI);
#endif
}

#endif

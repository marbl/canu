
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
/* 	$Id: MultiAlignStore_CNS.h,v 1.9 2007-02-18 14:04:48 brianwalenz Exp $	 */
#ifndef MULTIALIGNSTORE_H
#define MULTIALIGNSTORE_H

#include "AS_MSG_pmesg.h"
#include "AS_UTL_Var.h"
#include "AS_PER_gkpStore.h"
#include "PrimitiveVA_MSG.h"


#define MULTIALIGN_FORCED (-1)

typedef struct {
  int32 forced;// alignment is fishy if value == MULTIALIGN_FORCED
  int32 id;   /* this is the ID of the definition.  
		 Use this ONLY for persistence of references! */
  int32 refCnt;
  int source_alloc;
  VA_TYPE(char) *consensus;   // gapped consensus
  VA_TYPE(char) *quality;       // quality
  VA_TYPE(int32) *delta;       // deltas for all fragments in f_list
  VA_TYPE(IntMultiPos) *f_list;   // positions of fragments
  VA_TYPE(IntMultiVar) *v_list;   // variations                  
  VA_TYPE(int32) *udelta;       // deltas for all unitigs in u_list
  VA_TYPE(IntUnitigPos) *u_list;  // positions of unitigs
} MultiAlignT;

int32 GetMultiAlignLength(MultiAlignT *ma);
int32 GetMultiAlignUngappedLength(MultiAlignT *ma);

// Special accessors for the extremal unitigPos elements
IntUnitigPos *GetAendUnitigPos(MultiAlignT *ma);
IntUnitigPos *GetBendUnitigPos(MultiAlignT *ma);

void GetMultiAlignUngappedConsensus(MultiAlignT *ma, VA_TYPE(char) *ungappedSequence, VA_TYPE(char) *ungappedQuality);

void GetMultiAlignConsensusFromInterval(MultiAlignT *ma, SeqInterval range, VA_TYPE(char) *sequence, VA_TYPE(char) *quality);
void GetMultiAlignUngappedConsensusFromInterval(MultiAlignT *ma, SeqInterval range, VA_TYPE(char) *ungappedSequence, VA_TYPE(char) *ungappedQuality);


void GetMultiAlignUngappedOffsets(MultiAlignT *ma, VA_TYPE(int32) *ungappedOffsets);

// To maintain compatibility with old MultiAlignTs in checkpoints, we use
// a funky encoding of the forced condition.

static int  IsMultiAlignForced(MultiAlignT *ma){
  return ma->forced == MULTIALIGN_FORCED;
}

static void SetMultiAlignForced(MultiAlignT *ma, int32 trueFalse){
  if(trueFalse)
    ma->forced = MULTIALIGN_FORCED;
  else
    ma->forced = FALSE;
}

MultiAlignT *CreateMultiAlignT(void);
MultiAlignT *CreateEmptyMultiAlignT(void);
MultiAlignT *CreateMultiAlignTFromIUM(IntUnitigMesg *ium, int32 localFragID, int sequenceOnly);
MultiAlignT *CreateMultiAlignTFromICM(IntConConMesg *ium, int32 localFragID, int sequenceOnly);
MultiAlignT *CreateMultiAlignTFromCCO(SnapConConMesg *ium, int32 localFragID, int sequenceOnly);
MultiAlignT *CloneMultiAlignT(MultiAlignT *ma);
void CopyMultiAlignT(MultiAlignT *newma, MultiAlignT *oldma);

int CompareMultiAlignT(MultiAlignT *this_guy, MultiAlignT *other);

// Create Surrogate
MultiAlignT *CloneSurrogateOfMultiAlignT(MultiAlignT *oldMA, int32 newNodeID);

// Delete the multiAlign store and all of its referenced data
void DeleteMultiAlignT(MultiAlignT *multiAlign);

// Persistence
size_t SaveMultiAlignTToStream(MultiAlignT *ma, FILE *stream);
MultiAlignT *LoadMultiAlignTFromStream(FILE *stream, int32 *reference);
void ReLoadMultiAlignTFromStream(FILE *stream, MultiAlignT *ma, int32 *reference);

// Stats
size_t GetMemorySize(MultiAlignT *ma);


VA_DEF(MultiAlignT)

typedef struct {
    VA_TYPE(MultiAlignT *)  multiAligns;     
}MultiAlignStoreT;

typedef MultiAlignStoreT UnitigStoreT;

MultiAlignStoreT *CreateMultiAlignStoreT(int32 size);
// Delete the multiAlignStore and all of its referenced data
void DeleteMultiAlignStoreT(MultiAlignStoreT *multiAlignStore);

// Clear the multiAlignStore
size_t ClearMultiAlignStoreT(MultiAlignStoreT *multiAlignStore);

// Persistence
void SaveMultiAlignStoreTToStream(MultiAlignStoreT *ma, FILE *stream, int withReferences);
MultiAlignStoreT *LoadMultiAlignStoreTFromStream(FILE *stream);
MultiAlignStoreT *LoadMultiAlignStoreTFromStreamWithReferences(FILE *stream, MultiAlignStoreT *original);

// Reference Count Management
int32 AddReferenceMultiAlignT(MultiAlignT *ma);
int32 RemoveReferenceMultiAlignT(MultiAlignT *ma);
static int32 GetReferenceCountMultiAlignT(MultiAlignT *ma){
  return ma->refCnt;
}

// Clone
MultiAlignStoreT *CloneMultiAlignStoreT(MultiAlignStoreT *original);
int PrintMultiAlignT(FILE *out,
                     MultiAlignT *ma,
                     GateKeeperStore *gkp_store,
                     tFragStorePartition *pfrag_store,
                     int show_qv,
                     int dots,
                     uint32 clrrng_flag);

int GetCoverageInMultiAlignT(MultiAlignT *ma, SeqInterval range,
                VA_TYPE(int) *coverage, int includeExternal);

// Accessors
void SetMultiAlignInStore(MultiAlignStoreT *multiAlignStore, int index, MultiAlignT *multiAlign);
MultiAlignT *GetMultiAlignInStore(MultiAlignStoreT *multiAlignStore,  int index);
size_t RemoveMultiAlignFromStore(MultiAlignStoreT *mas, int index);
void CollectStats(MultiAlignT *ma,
                  GateKeeperStore *gkp_store, 
                  FILE *column_stats,
                  FILE *frag_stats,
                  uint32 clrrng_flag);

static int32 GetNumMultiAlignsInStore(MultiAlignStoreT *multiAlignStore){
  return GetNumPtrTs(multiAlignStore->multiAligns);
}
// Report Stats
int64 StatsMultiAlignStore(MultiAlignStoreT *maStore, FILE *fout, int owner);

// Sort Unitigs left to right
void MakeCanonicalMultiAlignT(MultiAlignT *ma);

MultiAlignT *RevcomplMultiAlignT(MultiAlignT *ma);
#endif

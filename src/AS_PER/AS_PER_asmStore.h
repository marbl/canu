
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
/* 	$Id: AS_PER_asmStore.h,v 1.8 2007-02-28 08:04:51 brianwalenz Exp $	 */
#ifndef AS_PER_ASMSTORE_H
#define AS_PER_ASMSTORE_H
/*************************************************************************
 Module:  AS_PER_asmStore
 Description:
    A thin layer on top of the IndexStore supporing the storage and
 retrieval of records from an assembly output file.
    The idea is to provide easier to use shortcuts for the common
 operations, and let the other operations be accessed through the
 generic Index Store API.

 Assumptions:
    Nothing special beyond genericStore.rtf

 Document:
      GenericStore.rtf

 *************************************************************************/


#include <time.h>

#include "AS_MSG_pmesg.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_PHash.h"

#define ASM_UID_NAMESPACE 'U'

// overload some iid types for assembly store
//
// BPW - why are we overlading and not using a different namespace?
//
#define AS_IID_AFG  0
#define AS_IID_MDI  1
#define AS_IID_UTG  2
#define AS_IID_CCO  3
#define AS_IID_DSC  4
#define AS_IID_SCF  5
#define AS_IID_CHR  6

typedef int32     ASM_BucketRecord;
typedef CDS_IID_t ASM_IIDRecord;

typedef struct
{
  uint32      containerIndex;
  SeqInterval pos;
  uint32      next;
} ASM_InstanceRecord;

typedef struct
{
  Dist_ID  uid;
  float32  inMean;
  float32  inStddev;
  float32  asmMean;
  float32  asmStddev;
  int32    numBuckets;
  int32    firstBucket;
} ASM_MDIRecord;


typedef struct
{
  Fragment_ID    uid;
  unsigned int   chaff:1;
  unsigned int   chimeric:1;
  unsigned int   inDegenerate:1;
  unsigned int   unreferenced:1;
  unsigned int   inSurrogate:1;
  unsigned int   deleted:1;
  unsigned int   spare1:2;
  unsigned int   type:8;
  unsigned int   unused:16;
  ASM_IIDRecord  mate;
  ASM_IIDRecord  library;
  SeqInterval    inClr;
  SeqInterval    asmClr;
  IntChunk_ID    unitigIndex;
  SeqInterval    unitigPos;
  int32          cInsIndex;
  int32          sInsIndex;
} ASM_AFGRecord;


typedef struct
{
  Chunk_ID       uid;
  float32        coverageStat;
  unsigned int   inDegenerate:1;
  unsigned int   spare1:7;
  unsigned int   status:8;
  unsigned int   spare2:16;
  uint32         numInstances;
  CDS_COORD_t    length;
  int32          numFrags;
  int32          firstFrag;
  int32          cInsIndex;
  int32          sInsIndex;
} ASM_UTGRecord;


typedef struct
{
  Contig_ID      uid;
  CDS_COORD_t    length;
  unsigned int   inDegenerate:1;
  unsigned int   numFrags:31;
  int32          firstFrag;
  int32          numUnitigs;
  int32          firstUnitig;
  IntScaffold_ID scaffoldIndex;
  SeqInterval    scaffoldPos;
} ASM_CCORecord;


typedef struct
{
  Scaffold_ID  uid;
  IntContig_ID contigIndex;
} ASM_DSCRecord;


typedef struct
{
  float32       asmMean;
  float32       asmStddev;
  CDS_COORD_t   storeMean;
} ASM_GapRecord;

typedef struct
{
  Scaffold_ID uid;
  int32       numContigs;
  int32       firstContig;
  int32       firstGap;
} ASM_SCFRecord;


#define ASMSTORE_DEF(type)\
typedef StoreHandle type ## Store;\
static int commit ## type ## Store(type ## Store sh){\
  return commitStore(sh);\
}\
static type ## Store reset ## type ## Store(type ## Store sh, int firstID){\
  return resetIndexStore(sh, firstID);\
}\
static int close ## type ## Store(type ## Store sh){\
  return closeStore(sh);\
}\
static int get ## type ## Store(type ## Store fs, int index, type ## Record *dr){\
  return getIndexStore(fs,index,dr); \
}\
static int set ## type ## Store(type ## Store fs, int index, type ## Record *dr){\
  return setIndexStore(fs,index,dr); \
}\
static type ## Store create ## type ## Store(char *StorePath, char *ext, int firstID){\
  type ## Store s = createIndexStore(StorePath,ext, sizeof(type ## Record), 1, firstID);\
  return s;\
}\
static type ## Store open ## type ## Store(char *StorePath, char *rw){\
  return openStore(StorePath, rw);\
}\
static int append ## type ## Store(type ## Store store, type ## Record *element){\
  return appendIndexStore(store,element);\
}\
static int32 getNum ## type ## s(type ## Store store){\
  StoreStat stat;\
  statsStore(store, &stat);\
  return(stat.lastElem);\
}

// Stores
#define NUM_ASM_FILES 17

// in addition to asm.phash, there are the following files

// for asm.mdi
ASMSTORE_DEF(ASM_MDI)

// for asm.bkt
ASMSTORE_DEF(ASM_Bucket)

// for asm.afg
ASMSTORE_DEF(ASM_AFG)

// for asm.aci, asm.asi, asm.uci, and asm.usi (instances)
ASMSTORE_DEF(ASM_Instance)

// for asm.utg
ASMSTORE_DEF(ASM_UTG)

// for asm.utf, asm.ccf, asm.ccu, asm.scc
ASMSTORE_DEF(ASM_IID)
  
// for asm.cco
ASMSTORE_DEF(ASM_CCO)

// for asm.dsc
ASMSTORE_DEF(ASM_DSC)

// for asm.scg
ASMSTORE_DEF(ASM_Gap)
  
// for asm.scf
ASMSTORE_DEF(ASM_SCF)

/* AssemblyStore */

typedef struct
{
  char gkpStorePath[FILENAME_MAX];
} ASM_Status;


typedef struct
{
  char storePath[FILENAME_MAX];

  // ASM_Status status;
  
  GateKeeperStore * gkpStore;
  
  ASM_MDIStore     mdiStore;
  ASM_BucketStore  bktStore;
  
  ASM_AFGStore     afgStore;
  ASM_InstanceStore aciStore;
  ASM_InstanceStore asiStore;
  
  ASM_UTGStore     utgStore;
  ASM_IIDStore     utfStore;
  ASM_InstanceStore uciStore;
  ASM_InstanceStore usiStore;
  
  ASM_CCOStore     ccoStore;
  ASM_IIDStore     ccfStore;
  ASM_IIDStore     ccuStore;
  
  ASM_DSCStore     dscStore;
  
  ASM_SCFStore     scfStore;
  ASM_GapStore     scgStore;
  ASM_IIDStore     sccStore;
  
  PHashTable_AS * hashTable;
}AssemblyStore;


int OpenGateKeeperStoreAssemblyStore(AssemblyStore * asmStore,
                                     char * gkpStorePath);
int CopyGateKeeperStoreLinksAssemblyStore(AssemblyStore * asmStore);

AssemblyStore * CreateAssemblyStore(char * path,
                                    char * gkpStorePath);

AssemblyStore * OpenAssemblyStore(char * path);
AssemblyStore * OpenReadOnlyAssemblyStore(char * path);
void CloseAssemblyStore(AssemblyStore *asmStore);

int CopyAssemblyStoreFiles(AssemblyStore *asmStore, char * destPath);
int RemoveAssemblyStoreFiles(AssemblyStore *asmStore);

//int TestOpenAssemblyStore(AssemblyStore *asmStore);

/*
  Map store - for additional mappings of scaffolds and/or fragments
 */
typedef struct
{
  CDS_UID_t  uid;
  char       chromosome[4];
  char       arm[4];
  char       description[256];
  int32      firstFrag;
  int32      lastFrag;
} ASM_CHRRecord;

#ifdef NEVER
#define NUM_MAP_FILES 4
ASMSTORE_DEF(ASM_Member)
#else
#define NUM_MAP_FILES 3
#endif
  
ASMSTORE_DEF(ASM_CHR)
  

typedef struct
{
  char storePath[FILENAME_MAX];
  ASM_CHRStore        chrStore; // chromosome store
#ifdef NEVER
  ASM_MemberStore     cfmStore; // frags in each chromosome;
#endif
  ASM_InstanceStore   finStore; // fragment instance store - ONE entry per frag
  PHashTable_AS *     hashTable; // for looking up chromosome UIDs
} MapStore;

MapStore * CreateMapStore(char * path);
MapStore * OpenMapStore(char * path);
MapStore * OpenReadOnlyMapStore(char * path);
void CloseMapStore(MapStore * mapStore);

int CopyMapStoreFiles(MapStore * mapStore, char * destPath);
int RemoveMapStoreFiles(MapStore * mapStore);


#endif

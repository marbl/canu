
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
static char CM_ID[] = "$Id: AS_UTL_PHash.c,v 1.7 2007-02-14 07:20:15 brianwalenz Exp $";
/*************************************************************************
 Module:  AS_UTL_PHash
 Description:
     This module defines the interface and implementation of the persistent
 hash table used to associated UIDs (aka accessions) with the assembler's
 IIDs (internal IDs).
     Persistent hash tables can either be created using a memory mapped file, or,
 transient tables can be created in memory.
     This module is not intended to be general purpose (for that see AS_UTL_Hash), but
 rather specifically focussed on the needs of the assembler.
     The hashing function used is from DDJ.  See AS_HashCommon.[hc].

 Assumptions:
      None.
 Document:
      PersistentHash.rtf

 *************************************************************************/
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "AS_global.h"
#include "AS_UTL_HashCommon.h"
#include "AS_UTL_PHash.h"

// to get around Tru64 UNIX v5 problem with memmapping > 500MB files
const unsigned long __sbrk_override = 1;

const char * const String_AS_IID[AS_IID_MAX + 2] = {
  "(Invalid)",
  "Fragment",
  "Distance",
  "Locale",
  "Sequence",
  "Bactig",
  "Plate",
  "Library",
  "Batch",
  "Donor",
  "(Invalid)"
};

//#define DEBUG_HASH

#define HASH_EMPTY -1
#define HASH_FN(key) (Hash_AS((unsigned char *)&(key),sizeof(CDS_UID_t),37))

/****************************************************************************/
/* Private functions used for hashtable management                                      */
/****************************************************************************/
/* Return to free list */
static void      FreePHashNode_AS(PHashTable_AS *table, PHashNode_AS *node);

/* Get from Free List */
/* Returns an index into the allocated array of the element */
/* COPIES value into the PHashNode */
static PHashNode_AS * AllocPHashNode_AS(PHashTable_AS **table, char nameSpace, CDS_UID_t key, PHashValue_AS *value, int refCount);

/***************************************************************************/
static int InsertNodeInPHashBucket(PHashTable_AS *table, PHashNode_AS *newnode);


/****************************************************************************/
/* Assign a type-specific IID  */
static void Assign_IID(PHashTable_AS *table, PHashNode_AS *node);


/****************************************************************************
 * Function: ConcatPHashTable_AS
 * Description:
 *     Inserts entries from source into target.
 *
 * Inputs:
 *     source         * to a PHashTable_AS.  
 * Input/Outputs
 *     target         ** to a PHashTable_AS. (May be reallocated as a result of operation )
 *
 * Return Value:
 *     HASH_SUCCESS if successful (value is valid)
 *     HASH_FAILURE if failure
 ***************************************************************************/
static int ConcatPHashTable_AS(PHashTable_AS **target, PHashTable_AS *source);


/***************************************************************************/

static int IndexOfPHashNode_AS(PHashTable_AS *table, PHashNode_AS *node){
  int nodeOffset = node - table->allocated;

  assert(nodeOffset >= 0 &&
	 nodeOffset <= table->numNodesAllocated);

  return(nodeOffset);

}
/*****************************************************************************/
static PHashNode_AS *GetPHashNode_AS(PHashTable_AS *table, int32 nodeOffset){

  if (!(nodeOffset >= 0 &&
	nodeOffset <= table->numNodesAllocated)){
#ifdef DEBUG_HASH
    fprintf(stderr," nodeOffset %d out of range (%d,%d)\n",
	    nodeOffset, 0, table->numNodesAllocated);
    assert(0);
#endif
    return NULL;
  }
  return(table->allocated + nodeOffset);

}

/***************************************************************************/
static void initializePHashTable(PHashTable_AS *table){

  table->lastKey = 0;
  table->collisions = 0;
  table->numNodes = 0;
  table->lastNodeAllocated = 0;
  {
    int i;

#ifdef DEBUG_HASH
    fprintf(stderr,"***Initializing %d buckets for %d elements, mask = 0x%x\n",
	    table->numBuckets, table->numNodesAllocated, table->hashmask);
#endif
    for(i = 0; i < table->numBuckets; i++){
      table->buckets[i] = HASH_EMPTY;
    }
  }	

  table->freeList = HASH_EMPTY;

#if COUNTS
  {
    int i;
#ifdef DEBUG_HASH
    fprintf(stderr,"***Initializing %d count buckets \n",
	    1<<LOG_NUM_TYPES);
#endif
    for(i = 0; i < 1<<LOG_NUM_TYPES; i++){
      table->counts[i] = 0;
    }
  }
#endif

}

/***************************************************************************/
static void  setupPHashTable(PHashTable_AS *table, char *allocate, int logBucketSize, int logAllocatedSize){

  table->buckets = (int *)(allocate + sizeof(PHashTable_AS));
  table->allocated = (PHashNode_AS *) (allocate +
				       sizeof(PHashTable_AS) +
				       hashsize(logBucketSize) * sizeof(uint32));

  table->numBuckets = hashsize(logBucketSize);
  table->hashmask = hashmask(logBucketSize);
  table->numNodesAllocated = hashsize(logAllocatedSize);

#ifdef DEBUG_HASH
  fprintf(stderr,"*setupPHashTable numBuckets %d numNodesAllocated = %d \n",
	  table->numBuckets, table->numNodesAllocated);
#endif
  initializePHashTable(table);
}


/***************************************************************************/
static int InsertNodeInPHashBucket(PHashTable_AS *table, PHashNode_AS *newnode){
  int32 hashkey = HASH_FN(newnode->key);
  int bucket = hashkey & table->hashmask ; 
  PHashNode_AS *node, *prevnode;
  int64 comparison;
  int nodeIndex;

#ifdef DEBUG_HASH
  fprintf(stderr,
	  "InsertNodeInHashBucket  value " F_IID " hashkey = " F_UID " bucket = %d\n",
	  newnode->value.IID, newnode->key, bucket);
#endif
  nodeIndex = table->buckets[bucket];

  /* ********* If bucket is empty, insert at head of bucket ********** */

  if(nodeIndex == HASH_EMPTY){   /* Should usually be this way */
    table->buckets[bucket] = IndexOfPHashNode_AS(table, newnode);
    return HASH_SUCCESS;
  }
  node = GetPHashNode_AS(table, nodeIndex);
#ifdef DEBUG_HASH
  fprintf(stderr,"*** Collision on insert in bucket %d for key " F_UID "\n",
	  bucket, newnode->key);
#endif

  table->collisions++;

  /* Use the nameSpaces first, and if they are the same, compare the keys */
  comparison = newnode->nameSpace - node->nameSpace; 
  if(comparison == 0){/* They are in the same namespace */
    comparison =  newnode->key - node->key;
  }
  if(comparison == 0){    /* element already present */
#ifdef DEBUG_HASH
    fprintf(stderr,"* Element Already Present!!!! ==> FAILURE\n");
#endif
    return(HASH_FAILURE_ALREADY_EXISTS);
  }

  /* *********            Insert at head of bucket         ********** */

  if(comparison > 0){ 
#ifdef DEBUG_HASH
    {
      fprintf(stderr,"*** Inserted before key " F_UID "\n",
	      node->key);
    }
#endif
    newnode->next = table->buckets[bucket];
    table->buckets[bucket] = IndexOfPHashNode_AS(table, newnode);
    return HASH_SUCCESS;
  }

  prevnode = node;
  nodeIndex = node->next;

  /* ********* Insert in the middle of the linked list  ********** */

  while(nodeIndex != HASH_EMPTY){
    node = GetPHashNode_AS(table,nodeIndex);
    /* Use the nameSpaces first, and if they are the same, compare the keys */
    comparison = newnode->nameSpace - node->nameSpace;
    if(comparison == 0){ /* They are in the same namespace */
      comparison = newnode->key - node->key;
    }
    if(comparison == 0){
#ifdef DEBUG_HASH
      fprintf(stderr,"* Element Already Present!!!! ==> FAILURE\n");
      /* element already present */
#endif
      /* element already present */
      return(HASH_FAILURE_ALREADY_EXISTS);
    }
    if(comparison > 0){
#ifdef DEBUG_HASH
      fprintf(stderr,"*** Inserted after key " F_UID " and before " F_UID "\n",
	      prevnode->key, node->key);
#endif
      newnode->next = IndexOfPHashNode_AS(table, node);
      prevnode->next = IndexOfPHashNode_AS(table, newnode);
      return(HASH_SUCCESS);
    }
    prevnode = node;
    nodeIndex = node->next;
  }
  /* ********* Append at the end of the linked list  ********** */

#ifdef DEBUG_HASH
  fprintf(stderr,"*** Inserted after key " F_UID "\n",    prevnode->key);
#endif
  newnode->next = HASH_EMPTY;
  prevnode->next = IndexOfPHashNode_AS(table,newnode);
  return(HASH_SUCCESS);


}

/***************************************************************************/
static int deleteFromPHashTable(PHashTable_AS *table, char nameSpace, CDS_UID_t key, int mark){
  int32 hashkey = HASH_FN(key);
  int bucket = hashkey & table->hashmask ; 
  PHashNode_AS  *previousNode;
  int32 nodeIndex;

  nodeIndex = table->buckets[bucket];
  /* If bucket is empty, return */
  if(nodeIndex == HASH_EMPTY)
    return HASH_FAILURE;  /* Failure, didn't find it */

  previousNode = NULL;

  while(nodeIndex != HASH_EMPTY){
    PHashNode_AS *node = GetPHashNode_AS(table, nodeIndex);
    int64 comparison;
    comparison = nameSpace - node->nameSpace;
    if(comparison == 0){ /* They are in the same namespace */
      comparison = key - node->key;
    }
    if(comparison == 0){  /* Found it */
      if(node->value.refCount != 1){ /* There are still references to this symbol, can't delete! */
#ifdef DEBUG_HASH
	fprintf(stderr,"* Outstanding Refs!!!! ==> FAILURE\n");
#endif
	return HASH_FAILURE_OUTSTANDING_REFS;
      }
      node->value.refCount = 0;
      if(mark){
	node->value.deleted = 1;
      }else{ /* Delete it */
	if(previousNode){
	  previousNode->next = node->next;
	}else{
	  table->buckets[bucket] = node->next;
	}
	FreePHashNode_AS(table, node);
      }
      table->isDirty = 1;
      return HASH_SUCCESS;
    }else if(comparison > 0){
      return HASH_FAILURE; /* The nodes are ordered by key */
    }

    previousNode = node;
    nodeIndex = node->next;
  }
  return HASH_FAILURE;

}
/***************************************************************************/
static void      FreePHashNode_AS(PHashTable_AS *table, PHashNode_AS *node){

  node->next = table->freeList;
  table->isDirty = 1;
  table->freeList = IndexOfPHashNode_AS(table, node);
  table->numNodes--;
#ifdef DEBUG_HASH
  fprintf(stderr,"* FreePHashNode_AS node %d IID " F_IID "\n",
	  table->freeList, node->value.IID);
  fflush(stderr);
#endif

}


/***************************************************************************/
static PHashTable_AS *ReallocToSize(PHashTable_AS *htable, int newsize){
  /* Reallocate the hashtable without writing it to disk */
  int oldNumNodes = htable->numNodesAllocated;
  // int oldNumBuckets = htable->numBuckets;
  // int logsize;
  int newNumNodes = oldNumNodes * 2;   /* Increase the size to the next power of two */
  PHashTable_AS *newTable;
  // char newBuffer[1024];
  // char *newFileName;
  // char oldBuffer[1024];
  // char *oldFileName;

  assert(newsize > 0);

  if(newsize > newNumNodes)
    newNumNodes =   hashsize(ceil_log2(newsize));


  /* Create a new PHashTable in a temporary file with a larger size */
#ifdef DEBUG_HASH
  fprintf(stderr,"** Reallocing PHashTable from %d to numNodes = %d\n", oldNumNodes, newNumNodes);
#endif
  newTable = CreatePHashTable_AS(newNumNodes,NULL);
  newTable->fp = htable->fp;                       /* Associate with the current open file */
  newTable->fileName = strdup(htable->fileName);   /* Associate with the current open file */
  assert(NULL != newTable); /* Rework to deal with new memory scheme and mmap */

  /* Concat the old hashtable to the new one */
  ConcatPHashTable_AS(&newTable,htable);

  /* Close the old hashtable and remove it */
  htable->fp = NULL;
  safe_free(htable->fileName);
  ClosePHashTable_AS(htable);  /* Free the memory associated with the old hashtable */

  newTable->isDirty = 1;
  return(newTable);

}


/***************************************************************************/
/* Assign a type-specific IID  */
static void Assign_IID(PHashTable_AS *table, PHashNode_AS *node){
  table->counts[node->value.type] = table->counts[node->value.type] + 1;
  node->value.IID = (table->counts[node->value.type]); /* NEVER recycle */
}


/***************************************************************************/
/* Get from Free List */
static PHashNode_AS *AllocPHashNode_AS(PHashTable_AS **table, 
				       char nameSpace, 
				       CDS_UID_t key, PHashValue_AS *value, int refCount){
  PHashNode_AS *node;
  int32 nodeIndex;


  if((*table)->numNodes + 1 >= (*table)->numNodesAllocated &&   // we've exhaused allocated space
     (*table)->freeList == HASH_EMPTY){                         // freelist is empty
    fprintf(stderr,"*#### Reallocating hashtable to size %d\n", (*table)->numNodesAllocated * 2);
    fflush(stderr);
    (*table) = ReallocToSize((*table), (*table)->numNodesAllocated * 2);
  }
  (*table)->numNodes++;
  if((*table)->freeList != HASH_EMPTY){
    nodeIndex = (*table)->freeList;
    node = GetPHashNode_AS((*table), nodeIndex);
    (*table)->freeList = node->next;
  }else{
    /* Allocate from the end of the file */
    nodeIndex = (*table)->lastNodeAllocated++;
    assert(nodeIndex <= (*table)->numNodesAllocated);
    // ####
#if 0
    if(nodeIndex > (*table)->numNodesAllocated - 10){
      fprintf(stderr,"* Allocated PHashNode %d/%d\n",
	      nodeIndex, (*table)->numNodesAllocated);
    }
#endif
    node = GetPHashNode_AS((*table), nodeIndex);
  }
  (*table)->isDirty = 1;
  node->next = HASH_EMPTY;
  node->key = key;
  node->nameSpace = nameSpace;
  node->value = *value;
  assert(refCount < PHASH_REFCOUNT_MAX);
  node->value.refCount = refCount;
  
  return node;
}

int updateRefCount(PHashTable_AS *table, char nameSpace, CDS_UID_t key, int increment);

/***************************************************************************/
/* Public functions for       hashtable management                                      */
/***************************************************************************/




/***************************************************************************/


PHashTable_AS *CreatePHashTable_AS(int numItemsToHash, char *pathToHashTable){
  PHashTable_AS *table;
  // PHashNode_AS *node;
  size_t totalHashSize;
  int logBucketSize ;
  int logAllocatedSize;
  char *allocate;

#ifdef DEBUG_HASH
  fprintf(stderr,"* numItemsToHash = %d hashsize(31) = %u\n", numItemsToHash, hashsize(31));
#endif

  assert(numItemsToHash > 0 &&
	 numItemsToHash < hashsize(31));

  /* Allocate 4 times as many buckets as items to hash to reduce collissions */
  logBucketSize = ceil_log2(numItemsToHash) + 2;

  /* Allocate the next power of two space for hash nodes */
  logAllocatedSize = ceil_log2(numItemsToHash);

  /* Compute the total size needed to hold header + buckets + nodes */
  totalHashSize = hashsize(logBucketSize) * sizeof(uint32) +            /* buckets */
    hashsize(logAllocatedSize) * sizeof(PHashNode_AS) +   /* nodes */
    sizeof(PHashTable_AS);                                /* header */

  allocate = (char *)safe_calloc(1, totalHashSize);

  table = (PHashTable_AS *)allocate;
  table->isDirty = 1; // Never been saved
  table->isReadWrite = 1;

  if(pathToHashTable == NULL){   /* In-memory "persistent" hash table */
    table->fp = NULL;
    table->fileName = NULL;
  }else{                         /* Memory saved to a file */

#ifdef DEBUG_HASH
    fprintf(stderr,"* pathToHashTable = %s totalHashSize = " F_SIZE_T "\n", pathToHashTable, totalHashSize);
#endif
    /* We need to create a file that is big enough for the entire memory
       mapped hash table.  Open it for write, and write one byte at the
       desired end of the file */

    table->fp = fopen(pathToHashTable,"w+b");

    /* We fail, if we can't create the file */
    if(table->fp == NULL)
      return NULL;

    table->fileName = strdup(pathToHashTable);
  }

#if 0
  fprintf(stderr,"* CreatePHashTable_AS logBucketSize:%d logAllocatedSize:%d numItemsToHash:%d\n",
	  logBucketSize, logAllocatedSize, numItemsToHash);
#endif
  setupPHashTable(table, allocate, logBucketSize, logAllocatedSize);

  return table;


}

/***************************************************************************/
static PHashTable_AS *OpenPHashTableCommon_AS(char *pathToHashTable, int32 readWrite){
  PHashTable_AS *table = NULL;
  FILE *fp;
  struct stat statbuf;
  char *allocate;
  PHashTable_AS header;


  fp = fopen(pathToHashTable,(readWrite?"r+":"r"));
#ifdef DEBUG_HASH
  fprintf(stderr,"* OpenHashTable pathToHashTable = %s\n" , pathToHashTable);
  assert(fp);
#endif
  if(fp == NULL)
    return NULL;

  /* fstat succeeds and filesize is reasonable */
  assert(0 == fstat(fileno(fp),&statbuf));

  if(statbuf.st_size <  (sizeof(PHashTable_AS)) ){
    fprintf(stderr,"* File size is screwed up " F_OFF_T "\n", statbuf.st_size);
    assert(0);
  }

  {
    /* Compute the total size needed to hold header + buckets + nodes */
    size_t totalHashSize, totalFileSize;
    size_t actualRead;

    CDS_FSEEK(fp,(off_t) 0,SEEK_SET);
    if(1 != fread(&header, sizeof(PHashTable_AS),1,fp))
      assert(0);

#ifdef DEBUG_HASH
    fprintf(stderr,"* Read file header numBuckets %d numNodesAllocated %d\n",
	    header.numBuckets, header.numNodesAllocated);
#endif

    totalFileSize = header.numBuckets * sizeof(uint32) +                /* buckets */
      header.lastNodeAllocated * sizeof(PHashNode_AS) +   /* nodes */
      sizeof(PHashTable_AS);                              /* header */

    totalHashSize = header.numBuckets * sizeof(uint32) +                /* buckets */
      header.numNodesAllocated * sizeof(PHashNode_AS) +   /* nodes */
      sizeof(PHashTable_AS);                              /* header */


#ifdef DEBUG_HASH
    fprintf(stderr,"* Open file read " F_SIZE_T " and allocate " F_SIZE_T "\n",
	    totalFileSize, totalHashSize);
#endif
    if(readWrite){
      allocate = (char *)safe_calloc(1, totalHashSize);

      CDS_FSEEK(fp,(off_t) 0,SEEK_SET);
      actualRead = fread(allocate, 1, totalFileSize, fp);
      if(actualRead < totalFileSize){
	fprintf(stderr,"* Read only " F_SIZE_T " bytes instead of " F_SIZE_T "\n",
		actualRead, totalFileSize);
	assert(0);
      }
      table = (PHashTable_AS *)allocate;
    }else{

      table = (PHashTable_AS *)safe_calloc(1, sizeof(PHashTable_AS));
      table->numBuckets = header.numBuckets;
      table->hashmask = header.hashmask;
      table->numNodesAllocated = header.numNodesAllocated;
#ifdef MAP_VARIABLE
      allocate = (char *) mmap(NULL,totalFileSize, PROT_READ,MAP_FILE|MAP_VARIABLE|MAP_PRIVATE,fileno(fp),0);
#else
      allocate = (char *) mmap(NULL,totalFileSize, PROT_READ,MAP_FILE|MAP_PRIVATE,fileno(fp),0);
#endif
      if(allocate == (caddr_t)-1){
	perror("mmap failed ");
	assert(0);
      }
    }
  }

  assert(table->numBuckets == header.numBuckets);
  table->isDirty = 0;
  table->fp = fp;
  table->fileName = strdup(pathToHashTable);
  table->buckets = (int *)(allocate + sizeof(PHashTable_AS));
  table->allocated = (PHashNode_AS *) (allocate +
				       sizeof(PHashTable_AS) +
				       table->numBuckets * sizeof(uint32));
  table->isReadWrite = readWrite;

  return table;
}


/***************************************************************************/
PHashTable_AS *OpenReadOnlyPHashTable_AS(char *pathToHashTable){

  return OpenPHashTableCommon_AS(pathToHashTable, FALSE);

}



/***************************************************************************/
PHashTable_AS *OpenPHashTable_AS(char *pathToHashTable){


  return OpenPHashTableCommon_AS(pathToHashTable, TRUE);
}

void PrintPHashValue_AS(FILE * fp, PHashValue_AS * value)
{
  fprintf(fp, "\t\tvalue->IID = " F_IID "\n", value->IID);
  fprintf(fp, "\t\tvalue->deleted = %d\n", value->deleted);
  fprintf(fp, "\t\tvalue->type = %d\n", value->type);
  fprintf(fp, "\t\tvalue->refCount = %d\n", value->refCount);
}

void PrintPHashNode_AS(FILE * fp, PHashNode_AS * node)
{
  fprintf(fp, "\tnode->key = " F_UID "\n", node->key);
  fprintf(fp, "\tnode->value:\n");
  PrintPHashValue_AS(fp, &(node->value));
  fprintf(fp, "\tnode->next = %d\n", node->next);
  fprintf(fp, "\tnode->nameSpace = %c\n", node->nameSpace);
  fprintf(fp, "\tnode->spare1 = %c\n", node->spare1);
  fprintf(fp, "\tnode->spare2 = %d\n", node->spare2);
}

void PrintPHashTableFields_AS(FILE * fp, PHashTable_AS * table)
{
  int i;
  
  fprintf(fp, "table->numBuckets = %d\n", table->numBuckets);
  fprintf(fp, "table->freeList = %d\n", table->freeList);
  fprintf(fp, "table->numNodes = %d\n", table->numNodes);
  fprintf(fp, "table->numNodesAllocated = %d\n", table->numNodesAllocated);
  fprintf(fp, "table->lastNodeAllocated = %d\n", table->lastNodeAllocated);
  fprintf(fp, "table->lastKey = " F_UID "\n", table->lastKey);
  fprintf(fp, "table->collisions = %d\n", table->collisions);
  fprintf(fp, "table->hashmask = %u\n", table->hashmask);
#if COUNTS
  for(i = 0; i < 1<<LOG_NUM_TYPES; i++)
  {
    fprintf(fp, "table->counts[%d] = " F_IID "\n", i, table->counts[i]);
  }
#endif
  /*
  fprintf(fp, "table->allocated:\n");
  for(i = 0; i < table->numNodes; i++)
  {
    PrintPHashNode_AS(fp, &(table->allocated[i]));
  }
  for(i = 0; i < table->numNodes; i++)
  {
    fprintf(fp, "table->buckets[%d] = %d\n", i, table->buckets[i]);
  }
  */
  fprintf(fp, "table->fileName = %s\n", table->fileName);
  fprintf(fp, "table->fp = %p\n", table->fp);
  fprintf(fp, "table->isDirty = %d\n", table->isDirty);
  fprintf(fp, "table->isReadWrite = %d\n", table->isReadWrite);
}

/***************************************************************************/
void ClosePHashTable_AS(PHashTable_AS *table){
  if(table->fp == NULL){
    //    fprintf(stderr,"* Closing memory-based table\n");
    safe_free(table);
  } else{
    size_t totalToWrite;

    if(table->isDirty){
#ifdef DEBUG_HASH
      fprintf(stderr,"* Closing dirty table\n");
#endif
      totalToWrite = table->numBuckets * sizeof(uint32) +   // buckets
	sizeof(PHashTable_AS) +                             // header
	table->lastNodeAllocated * sizeof(PHashNode_AS);    // nodes

      CDS_FSEEK(table->fp,(off_t) 0,SEEK_SET);
#if 0
      {
        const size_t numWritten = fwrite(table, 1, totalToWrite, table->fp);
        assert( numWritten == totalToWrite );
      }
#else // Implement the Digital UNIX fwrite bug work-around!
      // taken from AS_UTL_Var.c's CopyToFile_VA()
      {
        // The target size in bytes for a buffer
        const char * buff = (char *) table;
        const size_t maxWrite = (1 << 20);
        size_t increment;
        for( increment = 0; increment < totalToWrite; increment += maxWrite ){
          const size_t leftToWrite = (totalToWrite - increment);
          const size_t thisWrite = (leftToWrite < maxWrite ? leftToWrite : maxWrite);
          const size_t numWritten = fwrite((buff + increment), 1,
                                            thisWrite, table->fp);
          assert( numWritten == thisWrite );
        }
      }
#endif
    }else{
#ifdef DEBUG_HASH
      fprintf(stderr,"* Closing clean table\n");
#endif
    }

    safe_free(table->fileName);

    if(!table->isReadWrite)
    {
      PHashTable_AS header;
      size_t totalFileSize;
      char * allocate = ((char *) table->buckets) - sizeof(PHashTable_AS);
      
      CDS_FSEEK(table->fp,(off_t) 0,SEEK_SET);
      if(1 != fread(&header, sizeof(PHashTable_AS),1,table->fp))
        assert(0);
      totalFileSize = header.numBuckets * sizeof(uint32) + 
        header.lastNodeAllocated * sizeof(PHashNode_AS) +
        sizeof(PHashTable_AS);
      munmap(allocate, totalFileSize);
    }
    fclose(table->fp);
    safe_free(table);
  }
}

/***************************************************************************/
#if COUNTS
void ResetPHashTable_AS(PHashTable_AS *table, CDS_IID_t *counts){
#else
  void ResetPHashTable_AS(PHashTable_AS *table)  { 
#endif
    initializePHashTable(table);
#if COUNTS
    {
      int i;
      //fprintf(stderr,"***Initializing %d count buckets \n",    1<<LOG_NUM_TYPES);

      for(i = 0; i < 1<<LOG_NUM_TYPES; i++){
	table->counts[i] = counts[i];
      }
    }
#endif

  }

#if COUNTS

/***************************************************************************/
  void GetCountsPHashTable_AS(PHashTable_AS *table, CDS_IID_t *counts){
    int i;
    for(i = 0; i < 1<<LOG_NUM_TYPES; i++){
      counts[i] = table->counts[i];
    }

  }


/***************************************************************************/
  CDS_IID_t AllocateCountPHashTable_AS(PHashTable_AS *table, unsigned int type){

    assert(((int)type >= 0) && (type <= NUM_TYPES));
    table->counts[type]++;

    //   fprintf(stderr,"* returning " F_IID "\n", table->counts[type]);
    return   table->counts[type];

  }



#endif

/***************************************************************************/

  /* Insert in HashTable */
  /* Returns NULL if already in Hash Table */
  int InsertInPHashTable_AS
    (PHashTable_AS **table, char nameSpace, CDS_UID_t key,
     PHashValue_AS *value, int useRefCount, int assignIID){
    int refCount;
    PHashNode_AS *newnode;
    if(useRefCount)
      refCount = value->refCount;
    else
      refCount = 1;

    newnode = AllocPHashNode_AS(table, nameSpace, key, value, refCount);
    if(assignIID){
      Assign_IID(*table, newnode);
      value->IID = newnode->value.IID;
    }

    assert(NULL != newnode);
    newnode->value.deleted = 0;
    if(InsertNodeInPHashBucket(*table, newnode) == HASH_SUCCESS){
      (*table)->lastKey = key;
      (*table)->isDirty = 1;
      return HASH_SUCCESS;
    }
    // else
    {
      FreePHashNode_AS((*table), newnode);
      return HASH_FAILURE;
    }
  }
/***************************************************************************/

  int AddRefPHashTable_AS(PHashTable_AS *table, char nameSpace, CDS_UID_t key){

    table->isDirty = 1;
    return updateRefCount(table,nameSpace, key,1);
  }

/***************************************************************************/

  int UnRefPHashTable_AS(PHashTable_AS *table, char nameSpace, CDS_UID_t key){

    table->isDirty = 1;
    return updateRefCount(table,nameSpace, key,-1);
  }

/***************************************************************************/
  int updateRefCount(PHashTable_AS *table, char nameSpace, CDS_UID_t key, int increment){
    int32 hashkey = HASH_FN(key);
    int bucket = hashkey & table->hashmask ; 
    int32 nodeIndex;
#ifdef DEBUG_HASH
    fprintf(stderr,"UpdateRefCount of key " F_UID "\n", key);
#endif
    assert(bucket >= 0 && bucket <= table->numBuckets);

    nodeIndex = table->buckets[bucket];
    /* If bucket is empty, return */
    if(nodeIndex == HASH_EMPTY)
      return HASH_FAILURE;


    while(nodeIndex != HASH_EMPTY){
      PHashNode_AS *node = GetPHashNode_AS(table,nodeIndex);
      int64 comparison;

      comparison = nameSpace - node->nameSpace;
      if(comparison == 0){ /* They are in the same namespace */
	comparison = key - node->key;
      }

#ifdef DEBUG_HASH
      fprintf(stderr,"\tComparing key " F_UID " with node->key " F_UID ", comparison = %d\n",
	      key, node->key, comparison);
#endif
      if(comparison == 0){
#ifdef DEBUG_HASH
	fprintf(stderr,"\tFound it\n");
#endif
	if(node->value.deleted)
	  return HASH_FAILURE_FOUND_BUT_DELETED;

	assert(abs(node->value.refCount + increment) < PHASH_REFCOUNT_MAX);
	node->value.refCount += increment;
	/*      fprintf(stderr,"* Updated node with key " F_UID " to refcount %d\n",
		key, node->value.refCount); */

	return HASH_SUCCESS;
      }else if(comparison > 0){
#ifdef DEBUG_HASH
	fprintf(stderr,"\tGave Up\n");
#endif
	return HASH_FAILURE; /* The nodes are ordered by key */
      }

      nodeIndex = node->next;
    }
    return HASH_FAILURE;
  }


/***************************************************************************/
  int LookupTypeInPHashTable_AS(PHashTable_AS *phtbl, 
				char nameSpace,
				CDS_UID_t uid, 
				int type, 
				int reportFailure, 
				FILE *msgFile,
				PHashValue_AS *value){

    int lookup;

#ifdef DEBUG_HASH
    fprintf(stderr,"****LookupType key " F_UID " type %d lookupType %d\n",
	    uid, type, value->type);
#endif

    lookup = LookupInPHashTable_AS(phtbl, nameSpace, uid, value);
   

    if(lookup == HASH_FAILURE){
      if(reportFailure){
	fprintf(msgFile,"*LookupTypeInPHashTable_AS of uid "
		F_UID
		" (c%d) of type %d (%s) failed -- not found\n",
		uid, nameSpace, type,String_AS_IID[type]);
      }
      return HASH_FAILURE;
    }else if(value->deleted){
      if(reportFailure){
	fprintf(msgFile,"*LookupTypeInPHashTable_AS of uid "
		F_UID
		" (c%d) of type %d (%s) failed -- found but deleted\n",
		uid, nameSpace, type,String_AS_IID[type]);
      }
      return HASH_FAILURE_FOUND_BUT_DELETED;
    }else if(value->type != type){
      if(reportFailure){
	fprintf(msgFile,"*LookupTypeInPHashTable_AS of uid "
		F_UID
		" (c%d)of type %d (%s) failed"
		"-- found but wrong type %d (%s)\n",
		uid, nameSpace, type,String_AS_IID[type], value->type, String_AS_IID[value->type]);
      }
      return HASH_FAILURE_FOUND_BUT_WRONG_TYPE;
    }
    return HASH_SUCCESS;
  }
/***************************************************************************/
  int LookupInPHashTable_AS(PHashTable_AS *table, char nameSpace, CDS_UID_t key, PHashValue_AS *value){
    int32 hashkey = HASH_FN(key);
    int bucket = hashkey & table->hashmask ; 
    int32 nodeIndex;
#ifdef DEBUG_HASH
    fprintf(stderr,"Lookup of key " F_UID "\n", key);
#endif

    assert(bucket >= 0 && bucket <= table->numBuckets);

    value->deleted = 0;
    value->IID = 0;
    value->refCount = 0;

    nodeIndex = table->buckets[bucket];
    /* If bucket is empty, return */
    if(nodeIndex == HASH_EMPTY)
      return HASH_FAILURE;


    while(nodeIndex != HASH_EMPTY){
      PHashNode_AS *node = GetPHashNode_AS(table,nodeIndex);
      int64 comparison;


      comparison = nameSpace - node->nameSpace;
      if(comparison == 0){ /* They are in the same namespace */
	comparison = key - node->key;
      }

#ifdef DEBUG_HASH
      fprintf(stderr,"\tComparing key " F_UID " with node->key " F_UID ", comparison = %d\n",
	      key, node->key, comparison);
#endif

      if(comparison == 0){
#ifdef DEBUG_HASH
	fprintf(stderr,"\tFound it\n");
#endif
	*value = node->value;
	return HASH_SUCCESS;
      }else if(comparison > 0){
#ifdef DEBUG_HASH
	fprintf(stderr,"\tGave Up\n");
#endif
	return HASH_FAILURE; /* The nodes are ordered by key */
      }

      nodeIndex = node->next;
    }
    return HASH_FAILURE;
  }

/***************************************************************************/

  int DeleteFromPHashTable_AS(PHashTable_AS *table, char nameSpace, CDS_UID_t key){
    return deleteFromPHashTable(table,nameSpace, key,0);
  }
/***************************************************************************/
  int MarkAsDeletedPHashTable_AS(PHashTable_AS *table, char nameSpace, CDS_UID_t key){
    return deleteFromPHashTable(table,nameSpace, key,1);
  }


/***************************************************************************/

  int InitializePHashTable_Iterator_AS(PHashTable_AS *table, PHashTable_Iterator_AS *iterator){
    assert(table && iterator);
#ifdef DEBUG_HASH
    fprintf(stderr,"* initializeIterator nodes = (0, %d)\n",
	    table->numNodes);
#endif
    iterator->currentNodeIndex = 0;
    iterator->table = table;
    return HASH_SUCCESS;
  }


/***************************************************************************/

  int NextPHashTable_Iterator_AS(PHashTable_Iterator_AS *iterator, 
				 char *nameSpace, 
				 CDS_UID_t *key, 
				 PHashValue_AS *value){
    PHashNode_AS *node;

#ifdef DEBUG_HASH
    fprintf(stderr,"* NextIterator curr = %d numNodes = %d\n", 
	    iterator->currentNodeIndex, iterator->table->numNodes);
#endif
    if(iterator->currentNodeIndex >= iterator->table->lastNodeAllocated){
      return HASH_FAILURE;
    }

    node = GetPHashNode_AS(iterator->table, iterator->currentNodeIndex);
    iterator->currentNodeIndex++;
    *key = node->key;
    *value = node->value;
    *nameSpace = node->nameSpace;
#ifdef DEBUG_HASH
    fprintf(stderr,"* NextIterator key = " F_UID " value= %u nameSpace = %c\n", *key,*value, *nameSpace);
#endif
    return HASH_SUCCESS;
  }



/***************************************************************************/
  /* Concat two PHashTable_ASs...
     This really only works when one of them is EMPTY to start with.
     Utility function for resizing of Hashtables 
  */
static int ConcatPHashTable_AS(PHashTable_AS **target, PHashTable_AS *source){

    PHashTable_Iterator_AS iterator;
    PHashValue_AS value; 
    CDS_UID_t key;
    int result;
    char nameSpace;
    int i;

    /* If we will need to resize target during the concatenation, do it now! */

    result = MakeSpacePHashTable_AS(target, source->numNodes);
    assert(result == HASH_SUCCESS);
    if(result != HASH_SUCCESS)
      return HASH_FAILURE;

    // Copy the IID counters
    for(i = 0; i < 1<<LOG_NUM_TYPES; i++){
      (*target)->counts[i] = source->counts[i];
    }


#ifdef DEBUG_HASH
    fprintf(stderr,"** After MakeSpace target has nodes (%d/%d)\n",
	    (*target)->numNodes, (*target)->numNodesAllocated);
#endif
    /* ***********************************************************************/

    result = InitializePHashTable_Iterator_AS(source, &iterator);


    while(HASH_SUCCESS == NextPHashTable_Iterator_AS(&iterator, &nameSpace, &key, &value)){

    result = InsertInPHashTable_AS(target, nameSpace, key, &value, TRUE, FALSE /* Don't assign IIDs */);
#ifdef DEBUG_HASH
      fprintf(stderr,"* Concat:  Inserting Key " F_UID " ==> result = %d numNodes %d lastNode %d allocated %d\n", 
	      key, result,
	      (*target)->numNodes, (*target)->lastNodeAllocated, (*target)->numNodesAllocated);
#endif
      assert(result == HASH_SUCCESS);
      if(result != HASH_SUCCESS)
	return HASH_FAILURE;
    }
  
    return HASH_SUCCESS;
  }



/***************************************************************************/
  int MakeSpacePHashTable_AS(PHashTable_AS **table, int nodesToAdd){
   #ifdef DEBUG_HASH
    fprintf(stderr,"* MakeSpace for %d nodes (%d/%d already)\n",
	    nodesToAdd, (*table)->numNodes, (*table)->numNodesAllocated);
   #endif

    if((*table)->numNodes + nodesToAdd >= (*table)->numNodesAllocated ){
      *table = ReallocToSize(*table, (*table)->numNodes + nodesToAdd);
    }
    if(*table)
      return HASH_SUCCESS;
    //else
      return HASH_FAILURE;

  }






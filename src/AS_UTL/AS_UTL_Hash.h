
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
#ifndef AS_UTL_HASH_H
#define AS_UTL_HASH_H

#include "AS_UTL_heap.h"
#include "AS_UTL_HashCommon.h"
/***************** Hash Table Implementation ********************/


/* Type of function to compare hash table entries */
typedef int (*HashCmpFn_AS)(const void *, const void *);
typedef int (*HashFn_AS)(const void *, int length) ;

typedef struct HashNode_tag{
 void *key;
 uint32 keyLength;
  unsigned int isFree:1;
 void *value;
 struct HashNode_tag *next;
}HashNode_AS;

HEAP_DEF(HashNode_AS)

typedef struct{
  int32 numBuckets;
  HashNode_AS **buckets;
  HashNode_AS *freeList;
  HEAP_TYPE(HashNode_AS) *allocated;
  //HashNode_AS *allocated;
  int32 numNodes;
  int32 numNodesAllocated;
  int32 collisions;
  HashCmpFn_AS compare;
  HashFn_AS hash;
  uint32 hashmask;
} HashTable_AS;


int IntegerHashFn_AS(const void *pointerToInt, int length);
int StringHashFn_AS(const void *pointerToString, int length);

/********************* HashNodePool *********************/

HashTable_AS *CreateHashTable_AS(int numItemsToHash, HashFn_AS hash, HashCmpFn_AS cmp); /*);*/

HashTable_AS *CreateHashTable_int32_AS(int numItemsToHash);
HashTable_AS *CreateHashTable_uint64_AS(int numItemsToHash);

void ResetHashTable_AS(HashTable_AS *table);
void DeleteHashTable_AS(HashTable_AS *table);

static size_t ReportMemorySize_HT(HashTable_AS *ht,
			       const char *name, FILE *stream){
  size_t totalMemorySize = 0;
  fprintf(stream,"*\tStats for HashTable %s buckets = %d  hashnodes = %d\n",
	  name, ht->numBuckets, ht->numNodes);
  totalMemorySize += ht->numBuckets * sizeof(HashNode_AS);
  totalMemorySize += ReportMemorySize_HP(ht->allocated, "HashNodes", stream);
  return totalMemorySize;
}

/**************************************************************************************
 *  Insert in HashTable 
 *
 *   NOTE:  The hashtable stores only POINTERS to the key and value.
 *          Make sure you allocate space to store them that persists for the lifetime
 *	  of the hashtable!!!  Otherwise, all bets are off.
 *	  The AS_UTL_Heap.h functionality is strongly recommended
 **************************************************************************************/
int InsertInHashTable_AS(HashTable_AS *table, void *key, int keyLengthInBytes, void *value);


/* Insert in HashTable */
int DeleteFromHashTable_AS(HashTable_AS *table, void *key, int keyLengthInBytes);
void *LookupInHashTable_AS(HashTable_AS *table, void *key, int keyLengthInBytes);

/***********************************************************************************
 * HashTable_Iterator_AS
 *    Iterator for all allocated nodes in hashtable.
 *    In the absence of frees, values are returned in the order inserted
 *
 ***********************************************************************************/

typedef struct{
  HEAP_ITERATOR(HashNode_AS) iterator;
  HashTable_AS *table;
}HashTable_Iterator_AS;

/***********************************************************************************
 * Function: InitializeHashTable_Iterator_AS
 * Description:
 *     Initializes a PHashTable_Iterator_AS
 *
 * Inputs:
 *     table         * to a HashTable_AS.  
 * Input/Outputs
 *     iterator      * to a HashTable_Iterator_AS allocated by client
 *
 * Return Value:
 *     HASH_SUCCESS if successful (value is valid)
 *     HASH_FAILURE if failure
 ***********************************************************************************/
int InitializeHashTable_Iterator_AS(HashTable_AS *table, HashTable_Iterator_AS *iterator);

/***********************************************************************************
 * Function: NextHashTable_Iterator_AS
 * Description:
 *     Retrieves next element of HashTable using iterator
 *
 * Inputs:
 *     iterator      * to a HashTable_Iterator_AS
 * Input/Outputs
 *     value         * to a HashValue_AS *  (buffer is updated)
 *     key           * to a void *         (buffer is updated)
 *
 * Return Value:
 *     HASH_SUCCESS if successful (value is valid)
 *     HASH_FAILURE if failure
 ***********************************************************************************/
int NextHashTable_Iterator_AS(HashTable_Iterator_AS *iterator, 
			       void **key, 
			       void **value);


#endif




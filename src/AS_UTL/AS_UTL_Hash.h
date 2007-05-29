
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

#include "AS_global.h"
#include "AS_UTL_heap.h"
#include "math_AS.h"

#define HASH_FAILURE 0
#define HASH_SUCCESS 1

//  The hash function takes a 64-bit key (either a scalar or a pointer
//  to something) and associates it to a 64-bit value (again, scalar
//  or pointer).
//
//  The hash and compare functions need to be smart enough (as does
//  the client code) to know if the key/value is a pointer or a
//  scalar.

typedef int (*ASHashCompFn)(uint64 a, uint64 b);
typedef int (*ASHashHashFn)(uint64 k, uint32 l);

//  valueType is an annotation of the value present.  It is not used
//  by the hash table -- in particular, it will not distinguish
//  between identical keys.

typedef struct HashNode_AS{
  uint64               key;
  uint64               value;
  uint32               isFree:1;
  uint32               keyLength:23;
  uint32               valueType:8;
  struct HashNode_AS  *next;
}HashNode_AS;

HEAP_DEF(HashNode_AS)

typedef struct{
  uint32                   numBuckets;
  HashNode_AS            **buckets;
  HashNode_AS             *freeList;

  uint32                  numNodes;
  uint32                  numNodesAllocated;

  HEAP_TYPE(HashNode_AS) *allocated;

  uint32                  hashmask;

  uint32                  dirty;
  char                    filename[FILENAME_MAX];

  ASHashCompFn            compare;
  ASHashHashFn            hash;
} HashTable_AS;


typedef struct{
  HEAP_ITERATOR(HashNode_AS)   iterator;
  HashTable_AS                *table;
}HashTable_Iterator_AS;


uint32        Hash_AS(uint8 *k, uint32 length, uint32 initval);

HashTable_AS *CreateGenericHashTable_AS(uint32 numItemsToHash, ASHashHashFn hash, ASHashCompFn comp);
HashTable_AS *CreateScalarHashTable_AS(uint32 numItemsToHash);

void          ResetHashTable_AS (HashTable_AS *table);
void          DeleteHashTable_AS(HashTable_AS *table);

int           InsertInHashTable_AS     (HashTable_AS *table, uint64 key, uint32 keylen, uint64 value, uint32 valuetype);
int           DeleteFromHashTable_AS   (HashTable_AS *table, uint64 key, uint32 keylen);
uint64        ReplaceInHashTable_AS    (HashTable_AS *table, uint64 key, uint32 keylen, uint64 value, uint32 valuetype);

int           LookupInHashTable_AS     (HashTable_AS *table, uint64 key, uint32 keylen, uint64 *value, uint32 *valuetype);

int           ExistsInHashTable_AS     (HashTable_AS *table, uint64 key, uint32 keylen);
uint64        LookupValueInHashTable_AS(HashTable_AS *table, uint64 key, uint32 keylen);
uint32        LookupTypeInHashTable_AS (HashTable_AS *table, uint64 key, uint32 keylen);


void          SaveHashTable_AS(char *name, HashTable_AS *table);

HashTable_AS *LoadHashTable_AS(char *name,
                               ASHashHashFn  hashfn,
                               ASHashCompFn  compfn);
HashTable_AS *LoadUIDtoIIDHashTable_AS(char *name);

void          InitializeHashTable_Iterator_AS(HashTable_AS *table,
                                              HashTable_Iterator_AS *iterator);

int           NextHashTable_Iterator_AS(HashTable_Iterator_AS *iterator,
                                        uint64 *key,
                                        uint64 *value);

#endif




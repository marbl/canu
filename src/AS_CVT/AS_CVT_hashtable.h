
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
/* $Id: AS_CVT_hashtable.h,v 1.3 2005-03-22 19:04:49 jason_miller Exp $ */
#ifndef AS_CVT_HASHTABLE_H
#define AS_CVT_HASHTABLE_H

/*********************************************************************/
// headers
/*********************************************************************/
// project headers
#include "AS_global.h"
#include "AS_UTL_heap.h"

/* Bucket system */
// Buckets to store the lookup & return values
typedef struct BucketStruct
{
  cds_uint64   val;
  char       * ret;
  struct BucketStruct * next;
} Bucket;
typedef Bucket * Bucketp;

// other names for generics
typedef GenericArray  BucketArray;
typedef GenericArrayp BucketArrayp;
typedef GenericHeap   BucketHeap;
typedef GenericHeapp  BucketHeapp;

// bottomless stack
typedef struct
{
  cds_uint32  size;
  Bucketp     stack;
  BucketHeapp heap;
} BucketStack;
typedef BucketStack * BucketStackp;

typedef struct
{
  int            index;
  Bucketp        curr;
  BucketStackp   buckets;
  cds_uint32     num_nodes;
  Bucketp      * nodes;
} HashTable;
typedef HashTable * HashTablep;


/*********************************************************************/
// function declarations
/*********************************************************************/
/*
  Simple hash table with chaining to lookup uint64s & get something back
  If lookup value is not in hash table, a 0 is returned, so don't set
  the 'ret' value to 0
*/
/* Function:
     AllocateBucketHeap
   Description:
     allocates a heap of bucket objects
     wrapper around AllocateGenericHeap for bucket objects
   Return Value:
     valid BucketHeapp pointer if ok, else NULL
   Parameters:
     cds_uint32 num_buckets: initial number of bucket objects to allocate
*/
BucketHeapp AllocateBucketHeap( cds_uint32 num_buckets );

/* Function:
     FreeBucketHeap
   Description:
     frees a heap of bucket objects
     wrapper around FreeGenericHeap
   Return Value:
     none
   Parameters:
     BucketHeapp bucket_heap: heap of bucket objects
*/
void FreeBucketHeap( BucketHeapp bucket_heap );

/* Function:
     GetBucket
   Description:
     get a bucket from the bucket heap
     wrapper around GetGenericHeapItem
   Return Value:
     valid Bucketp if ok, else NULL
   Parameters:
     BucketHeapp bucket_heap: heap of bucket objects
     cds_uint32 size:         size of ret structure
*/
Bucketp GetBucket( BucketHeapp bucket_heap, cds_uint32 size );

/* Function:
     AllocateBucketArray
   Description:
     allocates an array of bucket objects
     wrapper around AllocateGenericArray
   Return Value:
     valid BucketArrayp if ok, else NULL
   Parameters:
     cds_uint32 num_buckets: number of bucket objects to allocate
*/
BucketArrayp AllocateBucketArray( cds_uint32 num_buckets );

/* Function:
     FreeBucketArray
   Description:
     frees an array of bucket objects
     wrapper around FreeGenericArray
   Return Value:
     none
   Parameters:
     BucketArrayp bucket_array: array of bucket objects
*/
void FreeBucketArray( BucketArrayp bucket_array );

/* Function:
     AllocateBucketStack
   Description:
     allocates bottomless stack of bucket objects using heap & stack
   Return Value:
     valid BucketStackp pointer if ok, else NULL
   Parameters:
     cds_uint32 num_buckets: initial number of buckets to allocate for heap
*/
BucketStackp AllocateBucketStack( cds_uint32 num_buckets, cds_uint32 size );

/* Function:
     FreeBucketStack
   Description:
     frees memory & structure of bucket stack
   Return Value:
     none
   Parameters:
     BucketStackp bs:  pointer to bucket stack structure
*/
void FreeBucketStack( BucketStackp bs );

/* Function:
     PushBucketOnStack
   Description:
     pushes one bucket object onto a bucket stack
   Return Value:
     none
   Parameters:
     BucketStackp bs:   pointer to bucket stack structure
     Bucketp bucket: pointer to bucket object
*/
void PushBucketOnStack( BucketStackp bs, Bucketp bucket );

/* Function:
     PopBucketOffStack
   Description:
     gets a bucket object from the bottomless stack of bucket objects
   Return Value:
     valid Bucketp if ok, else NULL
   Parameters:
     BucketStackp bs: pointer to bucket stack structure
*/
Bucketp PopBucketOffStack( BucketStackp bs );

void FreeHashTable( HashTablep ht );
HashTablep CreateHashTable( cds_uint32 num_nodes, cds_uint32 size );
char * LookupInHashTable( HashTablep ht, cds_uint64 val );
int InsertInHashTable( HashTablep ht, cds_uint64 val, char * ret );
int ChangeInHashTable( HashTablep ht, cds_uint64 val, char * ret );
void RewindHashTable( HashTablep ht );
char * GetNextHashTableRet( HashTablep ht );
int DeleteFromHashTable( HashTablep ht, cds_uint64 val );

#endif // AS_CVT_HASHTABLE_H


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
static char CM_ID[] = "$Id: AS_CVT_hashtable.c,v 1.2 2004-09-23 20:25:21 mcschatz Exp $";

/*********************************************************************/
// headers
/*********************************************************************/
// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

// project headers
#include "AS_UTL_heap.h"
#include "AS_CVT_hashtable.h"

/***********************************************************************/
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
BucketHeapp AllocateBucketHeap( cds_uint32 num_buckets )
{
  return( AllocateGenericHeap( num_buckets, sizeof( Bucket ) ) );
}

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
void FreeBucketHeap( BucketHeapp bucket_heap )
{
  FreeGenericHeap( bucket_heap );
}

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
Bucketp GetBucket( BucketHeapp bucket_heap, cds_uint32 size )
{
  Bucketp b = NULL;
  char * c;
  if( (c = (char *) calloc( 1, size )) != NULL )
  {
    if( (b = (Bucketp) GetGenericHeapItem( bucket_heap )) != NULL )
      b->ret = c;
    else
      free( c );
  }
  return b;
}

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
BucketArrayp AllocateBucketArray( cds_uint32 num_buckets )
{
  return( AllocateGenericArray( num_buckets, sizeof( Bucket ) ) );
}

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
void FreeBucketArray( BucketArrayp bucket_array )
{
  FreeGenericArray( bucket_array );
}

/* Function:
     AllocateBucketStack
   Description:
     allocates bottomless stack of bucket objects using heap & stack
   Return Value:
     valid BucketStackp pointer if ok, else NULL
   Parameters:
     cds_uint32 num_buckets: initial number of buckets to allocate for heap
*/
BucketStackp AllocateBucketStack( cds_uint32 num_buckets, cds_uint32 size )
{
  BucketStackp bs;
  
  if( (bs = (BucketStackp) calloc( 1, sizeof( BucketStack ) )) == NULL )
  {
    fprintf( stderr, "Failed to allocate bucket stack\n" );
    return NULL;
  }

  bs->size = size;
  bs->heap = AllocateBucketHeap( num_buckets );
  if( bs->heap == NULL )
  {
    fprintf( stderr,
             "Failed to allocate initial bucket heap for stack\n" );
    free( bs );
    return NULL;
  }
  return bs;
}

/* Function:
     FreeBucketStack
   Description:
     frees memory & structure of bucket stack
   Return Value:
     none
   Parameters:
     BucketStackp bs:  pointer to bucket stack structure
*/
void FreeBucketStack( BucketStackp bs )
{
  if( bs != NULL )
  {
    FreeBucketHeap( bs->heap );
    free( bs );
  }
}

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
void PushBucketOnStack( BucketStackp bs, Bucketp bucket )
{
  if( bucket->ret )
    free( bucket->ret );
  bucket->next = bs->stack;
  bs->stack = bucket;
}

/* Function:
     PopBucketOffStack
   Description:
     gets a bucket object from the bottomless stack of bucket objects
   Return Value:
     valid Bucketp if ok, else NULL
   Parameters:
     BucketStackp bs: pointer to bucket stack structure
*/
Bucketp PopBucketOffStack( BucketStackp bs )
{
  if( bs->stack != NULL )
  {
    char * c;
    if( (c = (char *) calloc( 1, bs->size )) != NULL )
    {
      Bucketp b = bs->stack;
      bs->stack = bs->stack->next;
      memset( b, 0, sizeof( Bucket ) );
      b->ret = c;
      return b;
    }
  }
  return( GetBucket( bs->heap, bs->size ) );
}

void FreeHashTable( HashTablep ht )
{
  if( ht )
  {
    int i;
  
    // need to free up all ret's
    for( i = 0; i < ht->num_nodes; i++ )
    {
      Bucketp b;
      for( b = ht->nodes[i]; b != NULL; b = b->next )
        free( b->ret );
    }
    FreeBucketStack( ht->buckets );
    if( ht->nodes )
      free( ht->nodes );
    free( ht );
  }
}


HashTablep CreateHashTable( cds_uint32 num_nodes, cds_uint32 size )
{
  HashTablep ht;
  
  ht = (HashTablep) calloc( 1, sizeof( HashTable ) );
  if( ht == NULL )
    return NULL;

  ht->buckets = AllocateBucketStack( num_nodes, size );
  if( ht->buckets == NULL )
  {
    FreeHashTable( ht );
    return NULL;
  }
  ht->curr = NULL;

  ht->nodes = (Bucketp *) calloc( num_nodes, sizeof( Bucketp ) );
  if( ht->nodes == NULL )
  {
    FreeHashTable( ht );
    return NULL;
  }
  ht->num_nodes = num_nodes;
  
  return ht;
}

char * LookupInHashTable( HashTablep ht, cds_uint64 val )
{
  Bucketp b;
  cds_uint32 key = val % ht->num_nodes;

  for( b = ht->nodes[key]; b != NULL; b = b->next )
  {
    if( b->val == val )
      return b->ret;
  }
  return NULL;

}

int InsertInHashTable( HashTablep ht, cds_uint64 val, char * ret )
{
  cds_uint32 key;
  Bucketp b;

  if( ret == NULL )
    return 1;
  if( LookupInHashTable( ht, val ) != NULL )
    return 1;
  if( (b = PopBucketOffStack( ht->buckets )) == NULL )
    return 1;
  key = val % ht->num_nodes;
  b->val = val;
  memcpy( b->ret, ret, ht->buckets->size );
  b->next = ht->nodes[key];
  ht->nodes[key] = b;
  return 0;
}

int ChangeInHashTable( HashTablep ht, cds_uint64 val, char * ret )
{
  Bucketp b;
  cds_uint32 key = val % ht->num_nodes;

  for( b = ht->nodes[key]; b != NULL; b = b->next )
  {
    if( b->val == val )
    {
      memcpy( b->ret, ret, ht->buckets->size );
      return 0;
    }
  }
  return 1;
}


void RewindHashTable( HashTablep ht )
{
  for( ht->index = 0; ht->index < ht->num_nodes; ht->index++ )
  {
    if( ht->nodes[ht->index] != NULL )
    {
      ht->curr = ht->nodes[ht->index];
      return;
    }
  }
}

char * GetNextHashTableRet( HashTablep ht )
{
  char * c = NULL;
  if( ht != NULL && ht->curr != NULL )
  {
    c = ht->curr->ret;
    if( ht->curr->next == NULL )
    {
      for( ++ht->index; ht->index < ht->num_nodes; ht->index++ )
      {
        if( ht->nodes[ht->index] != NULL )
        {
          ht->curr = ht->nodes[ht->index];
          return c;
        }
      }
      ht->curr = NULL;
      return c;
    }
    else
      ht->curr = ht->curr->next;
  }
  return c;
}


/*
int DeleteFromHashTable( HashTablep ht, cds_uint64 val )
{
  Bucketp b;
  cds_uint32 key = val % ht->num_nodes;

  if( ht->nodes[key]->val == val )
  {
    b = ht->nodes[key];
    ht->nodes[key] = ht->nodes[key]->next;
    PushBucketOnStack( ht->buckets, b );
    return 0;
  }
  for( b = ht->nodes[key]; b != NULL; b = b->next )
  {
    if( b->next && b->next->val == val )
    {
    }
    if( b->val == val )
    {
      memcpy( ht->buckets[ht->used_buckets].ret, ret, ht->buckets->size );
      return 0;
    }
  }
  return 1;
  
}
*/


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
static char CM_ID[] = "$Id: AS_UTL_heap.c,v 1.4 2005-03-22 19:49:29 jason_miller Exp $";

/*********************************************************************/
// headers
// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// project headers
#include "AS_global.h"
#include "AS_UTL_heap.h"

// constants
#define HEAP_FACTOR         1
#define HEAP_MIN            500

/*********************************************************************/

/* We probably don't need this...CM added it for Linux compatibility, but
   he manages with AS_UTL_Phash.c just fine...go figure. */
#define log2(x) (log(x)/log(2.))

/*********************************************************************/
/* Function:
     AllocateGenericHeap
   Description:
     allocates a generic heap structure & array
   Return Value:
     pointer to valid GenericHeap structure if ok
   Parameters:
     uint32 num_items: number of items to allocate
     size_t item_size: number of bytes per item to allocate
*/
GenericHeapp AllocateGenericHeap( uint32 num_items, size_t item_size )
{
  GenericHeapp heap = (GenericHeapp) calloc( 1, sizeof( GenericHeap ) );

  if( heap == NULL )
  {
    fprintf( stderr, "Failed to allocate generic heap.\n" );
    return( (GenericHeapp) NULL );
  }
  
  // allocate an initial heap array
  heap->first = AllocateGenericArray( num_items, item_size );
  if( heap->first == (GenericArrayp) NULL )
  {
    fprintf( stderr, "Failed to allocate first array in generic heap.\n" );
    FreeGenericHeap( heap );
    return( (GenericHeapp) NULL );
  }

  // set up other heap variables
  heap->current = heap->first;
  heap->first_num_items = heap->first->num_items;
  heap->item_size = item_size;
  
  return heap;
}


/* Function:
     GrowGenericHeap
   Description:
     increases the size of the heap by the original number of items divided
     by HEAP_FACTOR but no fewer than HEAP_MIN
   Return Value:
     0 if ok
   Parameters:
     GenericHeapp heap: pointer to heap structure
*/
int GrowGenericHeap( GenericHeapp heap )
{
   //  fprintf(stderr,"* Growing Heap %p\n", heap);
  return( GrowMinGenericHeap( heap, HEAP_MIN ) );
}

/* Function:
     GrowMinGenericHeap
   Description:
     increases the size of the heap by the original number of items divided
     by HEAP_FACTOR but no fewer than num_items
   Return Value:
     0 if ok
   Parameters:
     GenericHeapp heap: pointer to heap structure
     uint32 num_items: number by which to grow heap
*/
int GrowMinGenericHeap( GenericHeapp heap, uint32 num_items )
{
  // allocate an array & tack it onto the list of heap arrays
  heap->current->next =
    AllocateGenericArray( max( heap->first_num_items / HEAP_FACTOR,
                               num_items ), heap->item_size );
  if( heap->current->next == (GenericArrayp) NULL )
  {
    fprintf( stderr, "Failed to allocate generic array.\n" );
    return 1;
  }

  // set things to be ready to use the new array
  heap->current = heap->current->next;

  return 0;
}

/* Function:
     WipeGenericHeap
   Description:
     0s out the memory in a heap & makes it contiguous
   Return Value:
     0 if ok
   Parameters:
     GenericHeapp heap: pointer to heap structure
*/
int WipeGenericHeap( GenericHeapp heap )
{
  int num_items = 0;
  GenericArrayp arrayp1 = heap->first;
  GenericArrayp arrayp2 = heap->first;

  if( heap->current != heap->first )
  {
    // compute the size of heap to allocate
    while( arrayp1 != NULL )
    {
      num_items += arrayp1->num_items;
      arrayp1 = arrayp1->next;
      FreeGenericArray( arrayp2 );
      arrayp2 = arrayp1;
    }
    heap->current = 
       heap->first = AllocateGenericArray( num_items, heap->item_size );
    if( heap->first == NULL )
    {
      fprintf( stderr, "Failed to allocate array.\n" );
      return 1;
    }
    heap->first_num_items = heap->first->num_items;
  }
  else
  {
    assert(heap->current->num_items == heap->first_num_items);
    
    // Modification to make memset go faster - dewim 09/14/01
    assert(heap->first->curr_item <= heap->first_num_items &&
      (int) heap->first->array[heap->first->curr_item * heap->item_size] == 0);
    memset( (void *) heap->first->array, 0,
            heap->first->curr_item * heap->item_size );

    // Modification to keep heap from growing monotonically - dewim 09/14/01
    heap->first->curr_item = 0;
  }
  
  return 0;
}

/* Function:
     FreeGenericHeap
   Description:
     frees memory associated with a generic heap
   Return Value:
     none
   Parameters:
     GenericHeapp heap: pointer to heap structure
*/
void FreeGenericHeap( GenericHeapp heap )
{
  if( heap != NULL )
  {
    // walk through linked list of overlap arrays & free them
    GenericArrayp temp_array = heap->first;
    while( temp_array != (GenericArrayp) NULL )
    {
      heap->first = heap->first->next;
      FreeGenericArray( temp_array );
      temp_array = heap->first;
    }
    free( heap );
  }
}

/* Function:
     GetGenericHeapItem
   Description:
     returns a pointer to a GenericObject from the generic heap
     increases the size of the heap if necessary
     the number of writable bytes at that pointer is equal th item_size
   Return Value:
     valid pointer to a generic object if ok
   Parameters:
     GenericHeapp heap: pointer to heap structure
*/
GenericObjectp GetGenericHeapItem( GenericHeapp heap )
{
  // return something from the existing heap space, if possible
   if( heap->current->curr_item >= heap->current->num_items -1 ){
    // otherwise, try to increase the heap size
    if( GrowGenericHeap( heap ) )
    {
      fprintf( stderr, "Failed to grow generic heap.\n" );
      return( (GenericObjectp) NULL );
    }
  }
    return( heap->current->array + heap->item_size * heap->current->curr_item++ );
}

/* Function:
     GetGenericHeapItems
   Description:
     returns a pointer to an array of GenericObject from the generic heap
     increases the size of the heap if necessary
     the number of writable bytes at that pointer is equal to
     the item_size * the num_items specified
   Return Value:
     valid pointer to a generic object if ok
   Parameters:
     GenericHeapp heap: pointer to heap structure
     uint32 num_items: number of items to get
*/
GenericObjectp GetGenericHeapItems( GenericHeapp heap, uint32 num_items )
{
  uint32 ci = heap->current->curr_item;
  // return something from the existing heap space, if possible
  if( ci + num_items - 1 >= heap->current->num_items )
  {
    // otherwise, try to increase the heap size
    if( GrowMinGenericHeap( heap, num_items ) )
    {
      fprintf( stderr, "Failed to grow generic heap.\n" );
      return( (GenericObjectp) NULL );
    }
  }
  heap->current->curr_item += num_items;
  return heap->current->array + heap->item_size * ci;

}

/* Function:
     AllocateGenericArray
   Description:
     initial allocation of a generic array
   Return Value:
     valid GenericArrayp if ok
   Parameters:
     uint32 num_items: number of Generic objects to allocate
     size_t item_size: number of bytes per item
*/
GenericArrayp AllocateGenericArray( uint32 num_items, size_t item_size )
{
   int32 newSize, logNewSize, num_allocated;

  // allocate the generic array object
  GenericArrayp ret_ptr = (GenericArrayp) malloc( sizeof( GenericArray ) );
  if( ret_ptr == (GenericArrayp) NULL )
  {
    fprintf( stderr, "Failed to allocate generic array structure.\n" );
    return( (GenericArrayp) NULL );
  }

  newSize = num_items * item_size;
  logNewSize = (int) ceil(log2(newSize));  /* Find the next power of two that contains newSize */
  newSize = 1<<logNewSize; /* 2^logNewSize */
  num_allocated = newSize/item_size;

  // allocate the generic array itself
  ret_ptr->array =
    (GenericObjectp) calloc( num_allocated, item_size );
  if( ret_ptr->array == NULL )
  {
    fprintf( stderr,
             "Failed to allocate array of %d generic objects of size "
	     F_SIZE_T ".\n",
             num_items, item_size );
    free( ret_ptr );
    return( (GenericArrayp) NULL );
  }

  // set up other array variables
  ret_ptr->num_items = num_allocated;
  ret_ptr->curr_item = 0;
  ret_ptr->item_size = item_size;
  ret_ptr->next = (GenericArrayp) NULL;
  
  return ret_ptr;
}

/* Function:
     FreeGenericArray
   Description:
     Frees a GenericArray array & resets structure values
   Return Value:
     none
   Parameters:
     GenericArrayp generics:  pointer to GenericArray structure
*/
void FreeGenericArray( GenericArrayp generics )
{
  if( generics != NULL )
  {
    if( generics->array != NULL )
      free( generics->array );
    free( generics );
  }
}


/*** Heap Iterator Stuff ***/

/* Function:
     InitHeapIterator
   Description:
     Creates a heap iterator
   Return Value:
     none
   Parameters:
     GenericHeappp
*/

void InitHeapIterator(GenericHeap *heap, HeapIteratorT *iterator){
  iterator->heap = heap;
  iterator->current = heap->first;
  iterator->currentIndex = 0;
}

/* Function:
     NextHeapIterator
   Description:
     Returns next element ofheap, NULL if none
   Return Value:
     pointer to next element
   Parameters:
     GenericHeappp, HeapIteratorT*
*/

GenericObjectp NextHeapIterator(HeapIteratorT *iterator){
  GenericObjectp ret;

#ifdef DEBUG
  fprintf(stderr,"* NextHeapIterator currentIndex:%d current_items:%d current:%x\n",
	  iterator->currentIndex, iterator->current->curr_item, iterator->current);
#endif
  if(iterator->current->curr_item <= iterator->currentIndex){
    if(iterator->current->next == NULL)
      return NULL;
    assert(iterator->current != iterator->current->next);
    iterator->current = iterator->current->next;
    iterator->currentIndex = 0;
  }
  ret = (GenericObjectp)((char *)(iterator->current->array) + iterator->currentIndex * iterator->heap->item_size);
  iterator->currentIndex++;

  return ((GenericObjectp)ret);
}

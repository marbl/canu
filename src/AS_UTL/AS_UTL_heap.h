
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
#ifndef AS_UTL_HEAP_H
#define AS_UTL_HEAP_H

/**********************************************************************
$Id: AS_UTL_heap.h,v 1.4 2005-03-22 19:49:29 jason_miller Exp $
**********************************************************************/


/*********************************************************************/
// headers
// standard headers
#include <stdio.h>
#include <stdlib.h>

// project headers
#include "AS_global.h"
/*********************************************************************/

/*********************************************************************/


/*********************************************************************/
// structures
typedef char * GenericObjectp;

// structure for generic array
typedef struct GAStruct
{
  GenericObjectp    array;          // pointer for array of generic objects
  uint32            num_items;      // count of number of generic objects
  uint32            curr_item;      // current item, in issuing generic objects
  size_t            item_size;      // number of bytes per array item
  struct GAStruct * next;           // pointer to next array in linked list
} GenericArray;
typedef GenericArray * GenericArrayp;

// structure for heap of generics
typedef struct
{
  GenericArrayp first;           // pointer to first array in set
  GenericArrayp current;         // pointer to active array in set
  uint32        first_num_items; // number of items in first set
  size_t        item_size;       // number of bytes per item
} GenericHeap;
typedef GenericHeap * GenericHeapp;


static size_t ReportMemorySize_HP
( GenericHeap *heap,
  const char *name, FILE *stream)
{
  size_t totalMemorySize = 0;
  int totalItems = 0;
  if( heap != NULL )
  {
    // walk through linked list of overlap arrays & free them
    GenericArrayp temp_array = heap->first;
    while( temp_array != (GenericArrayp) NULL )
    {
      totalItems += temp_array->num_items;
      totalMemorySize += temp_array->num_items * temp_array->item_size;
      temp_array = temp_array->next;
    }
  }
  fprintf(stream,"*\tStats for Heap %s %d items " F_SIZE_T "bytes\n",  name,
	  totalItems, totalMemorySize);

  return totalMemorySize;
}



typedef struct{
  GenericHeap *heap;
  GenericArrayp current;
  int currentIndex;
}HeapIterator;

/*********************************************************************/
// function declarations
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
GenericHeapp AllocateGenericHeap( uint32 num_generics, size_t item_size );

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
int GrowGenericHeap( GenericHeapp heap );

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
int GrowMinGenericHeap( GenericHeapp heap, uint32 num_items );

/* Function:
     WipeGenericHeap
   Description:
     0s out the memory in a heap & makes it contiguous
   Return Value:
     0 if ok
   Parameters:
     GenericHeapp heap: pointer to heap structure
*/
int WipeGenericHeap( GenericHeapp heap );

/* Function:
     FreeGenericHeap
   Description:
     frees memory associated with a generic heap
   Return Value:
     none
   Parameters:
     GenericHeapp heap: pointer to heap structure
*/
void FreeGenericHeap( GenericHeapp heap );

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
GenericObjectp GetGenericHeapItem( GenericHeapp heap );

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
GenericObjectp GetGenericHeapItems( GenericHeapp heap, uint32 num_items );

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
GenericArrayp AllocateGenericArray( uint32 num_generics, size_t item_size );

/* Function:
     FreeGenericArray
   Description:
     Frees a GenericArray array & resets structure values
   Return Value:
     none
   Parameters:
     GenericArrayp generics:  pointer to GenericArray structure
*/
void FreeGenericArray( GenericArrayp generics );

/*** Heap Iterator Stuff ***/

typedef struct{
  GenericHeap *heap;
  GenericArrayp current;
  int currentIndex;
}HeapIteratorT;

/* Function:
     InitHeapIterator
   Description:
     Creates a heap iterator
   Return Value:
     none
   Parameters:
     GenericHeappp
*/

void InitHeapIterator(GenericHeap *heap, HeapIteratorT *iterator);

/* Function:
     NextHeapIterator
   Description:
     Returns next element ofheap, NULL if none
   Return Value:
     pointer to next element
   Parameters:
     GenericHeappp, HeapIteratorT*
*/

GenericObjectp NextHeapIterator(HeapIteratorT *iterator);

/********************************************************************/
/* TypeSafe Heap package
 * 
 *     Saul A. Kravitz
 *     March 1999
 *
 * This package is meant to simplify the coding and manipulation of
 * auto-resizing heaps for fixed size data.
 * It defines a basic set of operations, and provides a set of
 * macros that expand to support typesafe manipulation of the
 * heaps.
 *
 * Once a HEAP (variable array) has been defined, using the HEAP_DEF(<Type>)
 * macro, it may be manipulated as follows:
 *
 * A typedef is defined Heap<Type> for the VA, and is accessible via the
 * HEAP_TYPE(<Type>) macro.
 *
 * Allocate a new HEAP_Type(<Type>)
 * HEAP_TYPE(Type) * Allocate<Type>Heap (uint32 numElements)
 *
 * Free a HEAP_TYPE(<Type>)
 * void Free<Type>Heap (HEAP_TYPE(Type) *heap);
 *
 * Allocate an item fromt he heap
 * <Type> *Get<Type>HeapItem(HEAP_TYPE(Type) *heap);
 *
 * Allocate items from the heap
 * <Type> *Get<Type>HeapItems(HEAP_TYPE(Type) *heap);
 *
 */



#define HEAP_TYPE(Type) Heap ## Type
#define HEAP_ITERATOR(Type) Heap ## Type ## IteratorT

#define HEAP_DEF(Type)\
typedef GenericHeap Heap ## Type ;\
static HEAP_TYPE(Type) * Allocate ## Type ## Heap (uint32 numElements){\
     return ( (HEAP_TYPE(Type) *)AllocateGenericHeap(numElements, sizeof(Type)));\
}\
static void Free ## Type ## Heap (HEAP_TYPE(Type) *heap){\
     FreeGenericHeap(heap);\
}\
static Type *Get ## Type ## HeapItem(HEAP_TYPE(Type) *heap){\
     return ( (Type *)GetGenericHeapItem(heap));\
}\
static Type *Get ## Type ## HeapItems(HEAP_TYPE(Type) *heap, uint32 num_items){\
     return ( (Type *)GetGenericHeapItems(heap, num_items));\
}\
static int Wipe ## Type ## Heap(HEAP_TYPE(Type) *heap){\
      return WipeGenericHeap(heap);\
}\
typedef HeapIteratorT HEAP_ITERATOR(Type);\
static void Init ## Type ## HeapIterator(HEAP_TYPE(Type) *heap, HEAP_ITERATOR(Type) *iterator){\
  InitHeapIterator(heap, iterator);\
}\
static Type * Next ## Type ## HeapIterator (HEAP_ITERATOR(Type) *iterator){\
  return (( Type *)NextHeapIterator(iterator));\
}



/*********************************************************************/

#endif // #ifndef AS_UTL_HEAP_H

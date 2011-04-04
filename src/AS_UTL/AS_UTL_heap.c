
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

static char *rcsid = "$Id: AS_UTL_heap.c,v 1.12 2011-04-04 19:24:51 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_heap.h"

//  AllocateHeap_AS() returns an empty heap, with items_per_block
//  items allocated.  When those items are used, GetHeapItem_AS() will
//  allocate another block, but after doubling the items_per_block.
//  So, for small numbers of items, we keep our footprint small, and
//  for large numbers of items, we do not allocate thousands of
//  blocks.

Heap_AS *
AllocateHeap_AS(uint64 item_size, uint64 items_per_block) {
  Heap_AS *heap         = (Heap_AS *)safe_calloc(1, sizeof(Heap_AS));
  heap->items_per_block = items_per_block;
  heap->item_size       = item_size;
  heap->first           = (HeapArray_AS *)safe_calloc(1, sizeof(HeapArray_AS));
  heap->first->array    = safe_calloc(heap->items_per_block, item_size);
  heap->current         = heap->first;
  return(heap);
}

void
FreeHeap_AS(Heap_AS *heap) {
  if (heap == NULL)
    return;
  while (heap->first) {
    HeapArray_AS *t = heap->first->next;
    safe_free(heap->first->array);
    safe_free(heap->first);
    heap->first = t;
  }
  safe_free(heap);
}

void
ClearHeap_AS(Heap_AS *heap) {
  if (heap == NULL)
    return;
  //  While there is more than one block on the heap, remove the first
  //  block.  This saves the last (and largest) allocation.  VERY IMPORTANT!
  //  The last (and largest) allocation MUST be the one saved, because it
  //  is the current one -- heap->items_per_block is the size of the last
  //  block, not any other block.
  //
  //  One block is kept for performance reasons - in the usual use case
  //  of this function (reading messages) the heap has only one block
  //  allocated.
  //
  while (heap->first->next) {
    HeapArray_AS *t = heap->first->next;
    safe_free(heap->first->array);
    safe_free(heap->first);
    heap->first = t;
  }
  assert(heap->first        != NULL);
  assert(heap->first->array != NULL);
  assert(heap->first->next  == NULL);
  heap->first->nextAvail = 0;
  heap->current = heap->first;
}

void *
GetHeapItem_AS(Heap_AS *heap) {
  if (heap->current->nextAvail + 1 > heap->items_per_block) {
    heap->items_per_block     *= 2;
    heap->current->next        = (HeapArray_AS *)safe_calloc(1, sizeof(HeapArray_AS));
    heap->current->next->array = safe_calloc(heap->items_per_block, heap->item_size);
    heap->current              = heap->current->next;
  }
  uint64  na = heap->item_size * heap->current->nextAvail;
  heap->current->nextAvail += 1;
  return(((char *)heap->current->array) + na);
}

void *
GetHeapItems_AS(Heap_AS *heap, uint64 num_items) {
  if (num_items > heap->items_per_block) {
    //  We'd love to allocate exactly the size needed for just the next block, leaving the
    //  default block size the same, but cannot.  So, we'll just double the items
    //  per block until we fit.
    while (num_items > heap->items_per_block)
      heap->items_per_block *= 2;
    heap->current->next        = (HeapArray_AS *)safe_calloc(1, sizeof(HeapArray_AS));
    heap->current->next->array = safe_calloc(heap->items_per_block, heap->item_size);
    heap->current              = heap->current->next;
  }
  if (heap->current->nextAvail + num_items > heap->items_per_block) {
    //  This is the usual case; we just need another small block of memory,
    //  but the current block is too full.  Move along to the next block.
    heap->items_per_block     *= 2;
    heap->current->next        = (HeapArray_AS *)safe_calloc(1, sizeof(HeapArray_AS));
    heap->current->next->array = safe_calloc(heap->items_per_block, heap->item_size);
    heap->current              = heap->current->next;
  }
  uint64  na = heap->item_size * heap->current->nextAvail;
  heap->current->nextAvail += num_items;
  return(((char *)heap->current->array) + na);
}

void
InitHeapIterator_AS(Heap_AS *heap, HeapIterator_AS *iterator) {
  iterator->heap   = heap;
  iterator->array  = heap->first;
  iterator->item   = 0;
}

void *
NextHeapIterator_AS(HeapIterator_AS *iterator) {
  if (iterator->item >= iterator->array->nextAvail) {
    if (iterator->array->next == NULL)
      return(NULL);
    iterator->array = iterator->array->next;
    iterator->item  = 0;
  }
  return((char *)iterator->array->array + iterator->heap->item_size * iterator->item++);
}

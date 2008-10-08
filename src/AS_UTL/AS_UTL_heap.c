
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

static char *rcsid = "$Id: AS_UTL_heap.c,v 1.9 2008-10-08 22:03:00 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_heap.h"

Heap_AS *
AllocateHeap_AS(uint32 num_items, size_t item_size) {
  if (num_items < 16 * 1024)
    num_items = 16 * 1024;
  Heap_AS *heap = (Heap_AS *)safe_calloc(1, sizeof(Heap_AS));
  heap->first           = (HeapArray_AS *)safe_calloc(1, sizeof(HeapArray_AS));
  heap->first->array    = safe_calloc(num_items, item_size);
  heap->current         = heap->first;
  heap->num_items       = num_items;
  heap->item_size       = item_size;
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

void *
GetHeapItem_AS(Heap_AS *heap) {
  if (heap->current->nextAvail >= heap->num_items) {
    heap->current->next        = (HeapArray_AS *)safe_calloc(1, sizeof(HeapArray_AS));
    heap->current->next->array = safe_calloc(heap->num_items, heap->item_size);
    heap->current              = heap->current->next;
  }
  return((char *)heap->current->array + heap->item_size * heap->current->nextAvail++);
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


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

static const char *rcsid_AS_UTL_HEAP_H = "$Id: AS_UTL_heap.h,v 1.9 2008-12-05 19:06:12 brianwalenz Exp $";

#include "AS_global.h"

typedef struct HeapArray {
  void              *array;          // pointer for array of generic objects
  uint32             nextAvail;      // next available item, in issuing items
  struct HeapArray  *next;           // pointer to next array in linked list
} HeapArray_AS;

typedef struct {
  HeapArray_AS  *first;           // pointer to first array in set
  HeapArray_AS  *current;         // pointer to active array in set
  uint32         items_per_block; // number of items to allocate in the next set
  size_t         item_size;       // number of bytes per item
} Heap_AS;

typedef struct{
  Heap_AS      *heap;
  HeapArray_AS *array;
  int           item;
} HeapIterator_AS;

Heap_AS       *AllocateHeap_AS(size_t item_size);
void          *GetHeapItem_AS(Heap_AS *heap);
void           FreeHeap_AS(Heap_AS *heap);

void           InitHeapIterator_AS(Heap_AS *heap, HeapIterator_AS *iterator);
void          *NextHeapIterator_AS(HeapIterator_AS *iterator);

#endif // #ifndef AS_UTL_HEAP_H


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
#ifndef AS_ORA_REPEATS_H
#define AS_ORA_REPEATS_H
/**********************************************************************
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/AS_ORA_repeats.h,v $
$Revision: 1.4 $
**********************************************************************/


/*********************************************************************/
// headers
// project headers
#include "AS_global.h"
#include "AS_UTL_heap.h"
/*********************************************************************/


/*********************************************************************/
// types & structures
typedef uint64 RepeatID;

// RepeatObject
typedef struct ROStruct
{
  RepeatID          repeat_id;      // bit vector identifying the repeat(s)
  char              wc_complement;  // Watson-Crick complement sentinel
  uint64            count;          // number of repeats encountered
  uint64            rpt_begin;      // lower index of repeat in fragment
  uint64            rpt_end;        // upper index of repeat in fragment
  uint64            frag_begin;     // index in fragment where repeat begins
  uint64            frag_end;       // index in fragment where repeat ends
  struct ROStruct * next;           // pointer to next repeat in linked list
} RepeatObject;
typedef RepeatObject * RepeatObjectp;

// Generics
typedef GenericArray  RepeatArray;
typedef GenericArrayp RepeatArrayp;
typedef GenericHeap   RepeatHeap;
typedef GenericHeapp  RepeatHeapp;
/*********************************************************************/


/*********************************************************************/
// function declarations
// heap functions
/* Function:
     AllocateRepeatHeap
   Description:
     wrapper around AllocateGenericHeap to allocate a heap of RepeatObjects
   Return Value:
     valid RepeatHeap pointer if ok
   Parameters:
     uint32 num_repeats: number of repeats to allocate initially
*/
RepeatHeapp AllocateRepeatHeap( uint32 num_repeats );

/* Function:
     GrowRepeatHeap
   Description:
     wrapper around GrowGenericHeap
     increases size of existing repeat heap
   Return Value:
     0 if ok
   Parameters:
     RepeatHeapp repeat_heap: pointer to repeat heap
*/
int GrowRepeatHeap( RepeatHeapp repeat_heap );

/* Function:
     FreeRepeatHeap
   Description:
     wrapper around FreeGenericHeap
   Return Value:
     none
   Parameters:
     RepeatHeapp repeat_heap: pointer to repeat heap
*/
void FreeRepeatHeap( RepeatHeapp repeat_heap );

/* Function:
     GetRepeatObject
   Description:
     wrapper around GetGenericHeapItem
     gets a RepeatObject from a heap
   Return Value:
     valid RepeatObject pointer if ok
   Parameters:
     RepeatHeapp repeat_heap: pointer to repeat heap
*/
RepeatObjectp GetRepeatObject( RepeatHeapp repeat_heap );

/* Function:
     AllocateRepeatArray
   Description:
     wrapper around AllocateGenericArray
     initial allocation of an repeat array
   Return Value:
     valid RepeatArrayp if ok
   Parameters:
     uint32 num_repeats: number of Repeat objects to allocate
*/
RepeatArrayp AllocateRepeatArray( uint32 num_repeats );

/* Function:
     FreeRepeatArray
   Description:
     wrapper around FreeGenericArray
     frees a RepeatArray array & resets structure values
   Return Value:
     none
   Parameters:
     RepeatArrayp repeats:  pointer to RepeatArray structure
*/
void FreeRepeatArray( RepeatArrayp repeats );

/* Function:
     InsertionSortRepeat
   Description:
     inserts a new repeat object into a linked list of repeat objects in
     increasing order of repeat_id
   Return Value:
     none
   Parameters:
     RepeatObjectp * start: pointer to first pointer in linked list
     RepeatObjectp repeat: pointer to repeat object to insert
*/
void InsertionSortRepeat( RepeatObjectp * start, RepeatObjectp repeat );

/* Function:
     GetRepeatID
   Description:
     Returns the first or next annotation character from a repeat bit vector
   Return Value:
     char: the character representing the repeat
   Parameters:
     RepeatID id: the repeat id bit vector
*/
char GetRepeatID( RepeatID id );

/* Function:
     GetRepeatString
   Description:
     Produces a string listing all repeats in a repeat id bit vector
   Return Value:
     none
   Parameters:
     RepeatID id: the repeat id bit vector
     char * string: string to hold result
*/
void GetRepeatString( RepeatID id, char * string );
/*********************************************************************/
  
#endif // #ifndef AS_ORA_REPEATS_H
